from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
import Bio.PDB
import pandas as pd
import os
import requests
from bs4 import BeautifulSoup
from tabulate import tabulate

from app.residue_transformation import aa_conversion
from app.properties import ProteinFeatures, Interactions, Layer
from app.protein import Protein


class RandomForestModel:

######################################################################
  # Create a Random Forest classification model to predict if a residue
  # belongs or not to a binding site as well as assessing the model. 
  # Input: list of homologs from BLAST results
######################################################################


#  def __init__(self, directory):
#    self.directory =  directory

  def __init__(self):
    #self.directory =  directory
    pass
  
  def get_training_data(self, protein_object):
    
    # 0. AVOID REDUNDANT CALCULATIONS
    # Since we are working with chains separately, it is possible
    # that each chain has high similarity with others in terms of
    # homologous templates found after BLAST. So, we want to avoid
    # calculating features for a chain template, which is time expensive,
    # if these have already been calculated in a previous other chain

    # We start by creating a folder where the template datasets will be stored
    if not os.path.exists('./template_datasets'):
      os.makedirs('./template_datasets/')

    # If the folder already exists, we check if the csv with the features already exists  
    else:
      for filename in os.listdir('./template_datasets'):
      # Check if the search string is part of the filename
        if protein_object.protein_id[3:] in filename:

          # Read the dataframe, load it as a pandas dataframe and return it
          with open(os.path.join('./template_datasets', filename), 'r') as f:
            df = pd.read_csv(f, index_col=0)
            return df

    # 1. COMPUTE FEATURES OF EACH TEMPLATE
    # Calling data frame to add residues as first column
    #protein_object.dataframe_info()

    # Computing residue and atom features
    p_features = ProteinFeatures(protein_object, './app/data/atom_types.csv')
    p_features.residue_properties()
    p_features.atom_properties()
    p_features.is_cysteine()

    # Computing interactions:
    p_interactions = Interactions(protein_object, './app/data/atom_types.csv')
    p_interactions.calculate_interactions()

    # Compute layer features (atom and residue properties, and interactions)
    p_layer = Layer(protein_object, './app/data/atom_types.csv')
    p_layer.get_layer_properties()

    # 2. EXTRACT THE SITE LABELS OF THE TEMPLATE
    protein_object.dataframe["Label"] = pd.Series(dtype=int)

    protein_id = str(protein_object.structure.get_id())[3:7]
    chain_id = str(protein_object.structure.get_id())[7]

    list_binding_sites = ExtractBindingSites().extract_binding_sites(protein_id, chain_id)

    if list_binding_sites:

      for chain in protein_object.structure.get_chains():

        chain_id = chain.get_id()

        for residue in chain:

          if Bio.PDB.is_aa(residue, standard=True):  # ignore water molecules
            # Obtaining full residue name
            resname = residue.get_resname()
            resnum = residue.get_id()[1]
            res =  chain_id + '_' + resname + str(resnum)

            if res in list_binding_sites:
              protein_object.dataframe.loc[res,'Label'] = int(1)
        
            else:
              protein_object.dataframe.loc[res,'Label'] = int(0)
          
    else:
      return f"WARNING: {protein_object.protein_id} HAS NO BINDING SITES LABELS AND IT IS NOT CONSIDERED"
    
    # Saving dataframe into a folder
    protein_object.dataframe.to_csv('./template_datasets/' + protein_id + chain_id + '.csv')

    return protein_object.dataframe
  

  def concat_training_data(self, dataframes_list):
        
    # Create the template set
    train_df = pd.DataFrame()

    if len(dataframes_list) > 0:

      train_df = pd.concat(dataframes_list, axis = 0)

    return train_df
  

  def split_data(self, train_df):
    X = train_df.iloc[:, :-1]   # Select all the rows and remove first and last columns (residue_name and label)
    y = train_df.iloc[:, -1]     # Select as target the SITE label on the last column (Label)

    # Split the data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3 , random_state= 42, shuffle=True)

    return X_train, X_test, y_train, y_test

  def get_model(self, X_train, y_train):
    
    # We will set a threshold to which we consider class labels are unbalanced
    positive_data = y_train[y_train == 1]
    negative_data = y_train[y_train == 0]

    threshold = 3 # if negative data is 3 times more represented than positive data, we will
                  # consider that the dataset is highly unbalanced

    # if data is highly unbalanced, we add class_weight argument to give higher weight to
    # underrepresented class
    if len(negative_data)/len(positive_data) >= threshold:
      rf = RandomForestClassifier(n_estimators=100, max_depth = None, min_samples_split= 2, 
                                bootstrap = True, random_state=42, class_weight='balanced')
    
    # else, we will not consider data is too unbalanced and we will just run the classifier as is
    else:
      rf = RandomForestClassifier(n_estimators=100, max_depth = None, min_samples_split= 2, 
                                bootstrap = True, random_state=42)
    
    # fitting data into rf model
    rf.fit(X_train, y_train)

    return rf

  def get_predictions(self, dataset, cl_model):

    # Obtaining probability vector by applying classification model to query dataframe
    prob_vector = cl_model.predict_proba(dataset)[:,1]
    
    prob_threshold = 0.5

    pred_vector = []

    for prob in prob_vector:
        if prob < prob_threshold:
            pred_vector.append(0)
        else:
            pred_vector.append(1)

    dataset['prediction'] = pred_vector

    # This last condition is used so we can check the method worked when
    # calling the method from another script
    if len(pred_vector) > 0:
      return True
  

  def model_evaluation(self, y_prediction, y_test):
    accuracy = accuracy_score(y_test, y_prediction)
    precision = precision_score(y_test, y_prediction, average='macro', zero_division=0)
    recall = recall_score(y_test, y_prediction, average='macro')
    f1 = f1_score(y_test, y_prediction, average='macro')
    
    table = [['Metric', 'Score'],
                 ['Accuracy', round(accuracy, 2)],
                 ['Precision', round(precision, 2)],
                 ['Recall', round(recall, 2)],
                 ['F1', round(f1, 2)]]

    return tabulate(table)



class ExtractBindingSites:

######################################################################
  # Extract the residues that form binding sites to use them as 'labels' when 
  # creating the model to predict the binding sites, and allow us to evaluate 
  # the accuracy of the model
  # Input: Protein Object
  # Output: dictionary with binding sites
######################################################################

  def extract_binding_sites(self, protein_id, chain_id):

    # defining list where residues part of binding sites will appear
    final_residues = []

    # set up the URL of the BioLip site page
    for i in range(1,9):  

      url = "https://zhanggroup.org/BioLiP/pdb.cgi?pdb={}&chain={}&bs=BS0{}".format(protein_id, chain_id, i)

      # make a request to the page and retrieve the HTML content
      try:
        response = requests.get(url)
      except:
        print('Online BioLip database could not be accessed for current protein chain.')
        continue
      
      content = response.content

      # parse the HTML content using BeautifulSoup
      soup = BeautifulSoup(content, 'html.parser')

      # Find the tr tag that contains the specified text
      bindres_info = soup.find(string='(original residue number in PDB)')
      
      if bindres_info:
        bindres_line = bindres_info.parent.parent
      else:
        continue

      # Get the string containing the residues
      residue_str = bindres_line.text.strip()

      residues = residue_str.split(')')[1:]
      residues = residues[0]
      residues = residues.split(' ')

      # Converting residue names into data frame style
      residue_symbol = ''
      final_residues = []

      for residue in residues:
        res_num = residue[1:]
        one_let_res = residue[0]
        three_let_res = convert_residue(one_let_res)

        residue_symbol = str(chain_id + '_' + three_let_res + res_num)
        final_residues.append(residue_symbol)
    
    if len(final_residues) > 0:
      return final_residues
    
    # If no binding sites are found
    else:
      return None
    


### ASIDE FUNCTIONS

######################################################################
  # Converts 1-letter residues to 3-letter residues
  # Input: 1-letter residue
  # Output: 3-letter residue
######################################################################

def convert_residue(residue):

  converted_residue = aa_conversion[residue]
  return converted_residue