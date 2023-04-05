from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from protein import Protein
from properties import ProteinFeatures, Interactions, Layer
from Bio.PDB import PDBList
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


  def __init__(self, directory):
    self.directory =  directory

  
  def get_training_data(self, protein_object):
    
    # 1. COMPUTE FEATURES OF EACH TEMPLATE
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
    list_binding_sites = ExtractBindingSites(protein_object).extract_binding_sites()

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
              protein_object.dataframe.loc[res,'Label'] = '1'
        
            else:
              protein_object.dataframe.loc[res,'Label'] = '0'
        
      df_list.append(protein_object.dataframe)
    
    else:
      print("WE SHOULD REMOVE THIS TEMPLATE IT HAS NO BINDING SITES LABELS")

    # Create the template set
    train_df = pd.DataFrame()

    if df_list:
      train_df = pd.concat(df_list, axis = 0)

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

      else:
          print("No SITE section")
          
    if line_site > 0:

      # Return the binding sites dictionary
      return binding_sites

    else:

      return None