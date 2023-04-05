from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from protein import Protein
from properties import ProteinFeatures, Interactions, Layer
import Bio.PDB
from Bio.SeqUtils.ProtParam import ProteinAnalysis, ProtParamData
import pandas as pd
import numpy as np
import sys
import os


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
          return protein_object.protein_id[3:], "has already been calculated"

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

    df_list = []


    # 2. EXTRACT THE SITE LABELS OF THE TEMPLATE
    protein_object.dataframe["Label"] = pd.Series(dtype=int)
    list_binding_sites = ExtractBindingSites().extract_binding_sites(protein_object)

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
              protein_object.dataframe.loc[res,'Label'] = 1
        
            else:
              protein_object.dataframe.loc[res,'Label'] = 0
        
      df_list.append((protein_object.protein_id, protein_object.dataframe))
    
    else:
      print("WE SHOULD REMOVE THIS TEMPLATE IT HAS NO BINDING SITES LABELS")

    # Create the template set
    train_df = pd.DataFrame()

    if df_list:
      
      # Create empty list to append only dataframes
      only_df_list = []

      for protein in df_list:

        file_name = protein[0]
        protein_id = file_name[3:]
        df = protein[1]

        df.to_csv('./template_datasets/' + protein_id + '.csv')

        # appending dataframe to only_df_list
        only_df_list.append(df)

      # now actually creating whole template set by concatenating all dfs
      train_df = pd.concat(only_df_list, axis = 0)

    return train_df


  def split_data(self, train_df):
    X = train_df.iloc[:, 1:]     # Select all the rows and all but the first column (residue_name)
    y = train_df.iloc[:,:-1]     # Select as target the SITE label on the last column (Label)

    # Split the data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3 , random_state= 42)

    return X_train, X_test, y_train, y_test

  def get_model(self, X_train, y_train):
    # Train a random forest classifier on the training data
    rf = RandomForestClassifier(n_estimators=100, max_depth = None, bootstrap = True, random_state=42)
    rf.fit(X_train, y_train)

    return rf

  def get_predictions(self, protein_object):
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

  def extract_binding_sites(self, protein):

    # Create a dictionary to store the binding sites
    binding_sites = {}

    # Use the PDB file to extract the binding sites
    with open(protein.file_name, 'r') as f:

      pdb_f = f.read()

    # Condition to avoid errors if the file does not contain 'SITE' lines
    line_site = 0


    for line in pdb_f.split('\n'):

      if line.startswith('SITE'):

        #Counter of SITE lines
        line_site += 1
        
        # Save the name of the binding site
        site_name = line.split()[2]
 
        # Generate a list of the elements (since the first amino acid name to the end of the line)
        site_residues = line.split()[4:]
        
        # Save the list of residues in site_residues in the correct format. If the line does not follow the format, continue
        try: 
          # Format: A_Ser1 (example)
          site_residues = [f"{site_residues[i+1]}_{site_residues[i]}{site_residues[i+2]}" for i in range(0, len(site_residues), 3)]
          
        except Exception as error:
          continue
    
        if site_name not in binding_sites:

          binding_sites[site_name] = site_residues

        else:
          
          binding_sites[site_name] += site_residues

      #else:
       #   print("No SITE section")
          
    if line_site > 0:

      # Return the binding sites dictionary
      return binding_sites

    else:

      return None