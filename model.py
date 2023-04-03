import pandas as pd
import numpy as np
import os
import sys
import Bio.PDB
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from Bio.SeqUtils.ProtParam import ProteinAnalysis, ProtParamData
from protein import Protein, ProteinFeatures
from interactions import Interactions

class RandomForestModel:
  def __init__(self, directory):
    self.directory =  directory

  
  def get_training_data(self, list_homologs_final):

    if list_homologs_final:
      # Create a folder to store all the templates PDBs
      path = "./templates"
      if not os.path.exists(path):
        os.makedirs(path)

      
      # Extract the PDBs of each template
      pdbl = PDBList()
      df_final = []

      for pdb_id in list_homologs_final:
        f = pdbl.retrieve_pdb_file(pdb_id, pdir = path, file_format="pdb") # en GRASP los comprimen
        p = Protein(f)

        ## Compute the features, interactions ... of each template
        p_features = ProteinFeatures(p, './atom_types.csv')
        p_features.residue_properties()
        p_features.atom_properties()
        p_features.is_cysteine()

        ## Extract the labels of the template
        p.dataframe["Label"] = pd.Series(dtype=int)
        list_binding_sites = ExtractBindingSites(p).extract_binding_sites()

        if not list_binding_sites == None:

          for chain in p.structure.get_chains():

            chain_id = chain.get_id()

            for residue in chain:

              if Bio.PDB.is_aa(residue):  # ignore water molecules
                
                # Obtaining full residue name
                resname = residue.get_resname()
                resnum = residue.get_id()[1]
                res =  chain_id + '_' + resname + str(resnum)

                if res in list_binding_sites:
                  p.dataframe.loc[res,'Label'] = '1'
            
                else:
                  p.dataframe.loc[res,'Label'] = '0'
        
        else:

          print("WE SHOULD REMOVE THIS TEMPLATE IT HAS NO BINDING SITES LABELS")

        #df_p = p.dataframe.style
      
        df_final.append(p.dataframe)
    
    return df_final
  
  def split_data(self, df_final):
    X_train, X_test, y_train, y_test = train_test_split(df_final, y, test_size = 0.3 , random_state= 42)
    # X: all the features in the matrix
    # Y: label class? LAST COLUMN
    return X_train, X_test, y_train, y_test

  def get_model(self, X_train, y_train):
    # Train a random forest classifier on the training data
    rf = RandomForestClassifier(n_estimators=100, max_depth = None, bootstrap = True, random_state=42)
    rf.fit(X_train, y_train)
    return rf

  def get_predictions(self, model, X_test):
    pass

  def model_evaluation(self):
    # Here we should use the SITE labels of each PDB?
    pass



class ExtractBindingSites:

######################################################################
  # Extract the residues that form binding sites to use them as 'labels' when 
  # creating the model to predict the binding sites, and allow us to evaluate 
  # the accuracy of the model
  # Input: Protein Object
  # Output: list with binding sites residues
######################################################################
  
  def __init__(self, protein):
    self.protein = protein
  
  def extract_binding_sites(self):

    # Create a list to store the binding sites residues
    binding_sites = []

    # Use the PDB file to extract the binding sites
    with open(self.protein.file_name, 'r') as f:

      pdb_f = f.read()

    # Condition to avoid errors if the file does not contain 'SITE' lines
    line_site = 0


    for line in pdb_f.split('\n'):

      if line.startswith('SITE'):

        #Counter of SITE lines
        line_site += 1
        
        # Save the name of the binding site
        #site_name = line.split()[2]
 
        # Generate a list of the elements (since the first amino acid name to the end of the line)
        site_residues = line.split()[4:]
        
        # Save the list of residues in site_residues in the correct format. If the line does not follow the format, continue
        try: 
          # Format: A_Ser1 (example)
          site_residues = [f"{site_residues[i+1]}_{site_residues[i]}{site_residues[i+2]}" for i in range(0, len(site_residues), 3)]
          
          for element in site_residues:

            binding_sites.append(element)

        except Exception as error:
          continue
          
    if line_site > 0:

      # Return the binding sites dictionary
      return binding_sites

    else:

      return None
    
