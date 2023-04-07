from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from protein import Protein
from residue_transformation import aa_conversion
from properties import ProteinFeatures, Interactions, Layer
import Bio.PDB
from Bio.SeqUtils.ProtParam import ProteinAnalysis, ProtParamData
import pandas as pd
import numpy as np
import sys
import os
import requests
from bs4 import BeautifulSoup


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

    # we start by creating a folder where the template datasets will be stored
    if not os.path.exists('./template_datasets'):
      os.makedirs('./template_datasets/')

    # if the folder already exists, we check if the csv with the features already exists  
    else:
      for filename in os.listdir('./template_datasets'):
      # Check if the search string is part of the filename
        if protein_object.protein_id[3:] in filename:
          # If it is, print the filename
          print(protein_object.protein_id[3:], "features have already been calculated")

          # and read the dataframe, load it as a pandas dataframe and return it
          with open(os.path.join('./template_datasets', filename), 'r') as f:
            df = pd.read_csv(f)
            return df

    # 1. COMPUTE FEATURES OF EACH TEMPLATE
    # Calling data frame to add residues as first column
    #protein_object.dataframe_info()

    # Computing residue and atom features
    p_features = ProteinFeatures(protein_object, './atom_types.csv')
    p_features.residue_properties()
    p_features.atom_properties()
    p_features.is_cysteine()

    # Computing interactions:
    p_interactions = Interactions(protein_object, './atom_types.csv')
    p_interactions.calculate_interactions()

    # Compute layer features (atom and residue properties, and interactions)
    p_layer = Layer(protein_object, './atom_types.csv')
    p_layer.get_layer_properties()

    # 2. EXTRACT THE SITE LABELS OF THE TEMPLATE
    protein_object.dataframe["Label"] = pd.Series(dtype=int)

    protein_id = str(protein_object.structure.get_id())[3:7]
    chain_id = str(protein_object.structure.get_id())[7]

    list_binding_sites = ExtractBindingSites().extract_binding_sites(protein_id, chain_id)
    #print(list_binding_sites)

    if list_binding_sites:

      for chain in protein_object.structure.get_chains():

        chain_id = chain.get_id()

        for residue in chain:

          if Bio.PDB.is_aa(residue):  # ignore water molecules
            
            # Obtaining full residue name
            resname = residue.get_resname()
            resnum = residue.get_id()[1]
            res =  chain_id + '_' + resname + str(resnum)

            if res in list_binding_sites:
              protein_object.dataframe.loc[res,'Label'] = int(1)
        
            else:
              protein_object.dataframe.loc[res,'Label'] = int(0)
          
    else:
      print("WARNING: THIS TEMPLATE HAS NO BINDING SITES LABELS AND SHOULD BE REMOVED")
    
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
    X = train_df.iloc[:, 1:]     # Select all the rows and all but the first column (residue_name)
    y = train_df.iloc[:, -1]     # Select as target the SITE label on the last column (Label)

    # Split the data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3 , random_state= 42, shuffle=True)

    return X_train, X_test, y_train, y_test

  def get_model(self, X_train, y_train):
    # Train a random forest classifier on the training data
    rf = RandomForestClassifier(n_estimators=100, max_depth = None, min_samples_split= 2, 
                                bootstrap = True, random_state=42)
    rf.fit(X_train, y_train)

    return rf

  def get_predictions(self, protein_object, cl_model):

    

    template_dataframe = self.get_training_data(protein_object)
    X_train, X_test, y_train, y_test = self.split_data(template_dataframe)
    class_model = self.get_model(X_train, y_train)

    y_prediction = class_model.predict(X_test)

    return y_prediction

  def model_evaluation(self, y_prediction, y_test):
    accuracy = accuracy_score(y_test, y_prediction)
    precision = precision_score(y_test, y_prediction, average='macro')
    recall = recall_score(y_test, y_prediction, average='macro')
    f1 = f1_score(y_test, y_prediction, average='macro')
    
    return accuracy, precision, recall, f1



'''
# Faz a predição a partir do vetor de probabilidades de um classificador
  def predict_prob(self,prob_vector, threshold):
      pred_vector = []
      for prob in prob_vector:
          if prob < threshold:
              pred_vector.append(0)
          else:
              pred_vector.append(1)
      return pred_vector

  # voting using calssification matrix
  def voting(self, binary_matrix, treshold):
      binary_matrix = binary_matrix.T
      confidence = []
      for i in range(binary_matrix.shape[0]):
          # computa a porcentagem de classificadores na votacao
          confidence.append(np.sum(binary_matrix[i])/float(binary_matrix[i].shape[0]))#
      return self.predict_prob(confidence, treshold)

  def balanced_prediction(self, protein, out_dir, naccess):
      # search and store templates matricies
      self.set_train_data(protein, naccess)
      pos_data, neg_data, N_PARTITIONS = self.split_data(self.train_set)
      # shuffle negative index
      permuted_indices = np.random.permutation(len(neg_data))
      # matrix with class values foa all predictions
      # each line is a prediction of one ensenble
      class_matrix = []
      for i in range(N_PARTITIONS):
          # Concat positive and a fragment os negative instances
          final_matrix = pd.concat([pos_data, neg_data.iloc[permuted_indices[i::N_PARTITIONS]]])
          class_model = self.get_model(final_matrix)
          # probability_vector
          probs = class_model.predict_proba(protein.matrix)[:,1]
          # voting probabilities
          class_matrix.append(self.predict_prob(probs, 0.5))
          # cleaning memory
          del class_model
      # Fzendo predicao final usando a combinacao dos ensembles
      vector_pred = self.voting(np.array(class_matrix), 0.5)
      protein.matrix['prediction'] = vector_pred
      protein.matrix.to_csv(out_dir + os.path.basename(protein.file_name).replace('.pdb', '.csv'), columns=['prediction'])


'''







class ExtractBindingSites:

######################################################################
  # Extract the residues that form binding sites to use them as 'labels' when 
  # creating the model to predict the binding sites, and allow us to evaluate 
  # the accuracy of the model
  # Input: Protein Object
  # Output: dictionary with binding sites
######################################################################

  def extract_binding_sites(self, protein_id, chain_id):

    # set up the URL of the BioLip site page
    url = "https://zhanggroup.org/BioLiP/pdb.cgi?pdb={}&chain={}&bs=BS01".format(protein_id, chain_id)

    # make a request to the page and retrieve the HTML content
    response = requests.get(url)
    content = response.content

    # parse the HTML content using BeautifulSoup
    soup = BeautifulSoup(content, 'html.parser')

    # Find the tr tag that contains the specified text
    bindres_line = soup.find(string='(original residue number in PDB)').parent.parent

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

    return final_residues
    


### ASIDE FUNCTIONS

######################################################################
  # Converts 1-letter residues to 3-letter residues
  # Input: 1-letter residue
  # Output: 3-letter residue
######################################################################

def convert_residue(residue):

  converted_residue = aa_conversion[residue]
  return converted_residue