import pandas as pd
import os
import sys
import Bio.PDB
from Bio.SeqUtils.ProtParam import ProteinAnalysis, ProtParamData
import math

class ProteinFeatures:

  def __init__(self, protein, atom_props_dir):
    self.protein = protein
    self.atom_types_dict = read_atom_types(atom_props_dir)

######################################################################
  # Extract features of the residues
  # Input: list of chain objects
  # Output: dictionary with interesting features
######################################################################

  def residue_properties(self):

    # {'ChainID_ResNumResName' : [hydrophob, surf_acc], ...}
    #dict_properties = {}
    new_cols = ['hydrophob_value', 'surf_acc']

    # Create new columns for the dataframe
    for column in new_cols: 
      
      self.protein.dataframe[column] = pd.Series(dtype=float)

    for chain in self.protein.structure.get_chains():

      chain_id = chain.get_id()
      sequence = ''
      
      for residue in chain:

        if residue.get_resname() != 'HOH' and  Bio.PDB.is_aa(residue):  # ignore water molecules

          resname = residue.get_resname()
          resnum = residue.get_id()[1]
          res =  chain_id + '_' + resname + str(resnum)

          sequence += Bio.SeqUtils.seq1(residue.get_resname()) 

          for aa in sequence:
              
              hydrophobicity_value = ProtParamData.kd.get(aa, 0) 
              surface_acc_value = ProtParamData.ja.get(aa, 0)

              self.protein.dataframe.loc[res,'hydrophob_value'] = float( hydrophobicity_value )
              self.protein.dataframe.loc[res,'surf_acc'] = float( surface_acc_value )


    ### CALCULATING ATOM PROPERTIES
  
  def atom_properties(self):

    # Defining new columns for atom features
    atom_columns = ['don_count', 'acp_count', 'hpb_count', 'pos_count', 
                    'neg_count', 'arm_count']

    # Adding new columns to data frame
    for column in atom_columns: 
    
      self.protein.dataframe[column] = 0 #pd.Series(dtype=int)

    for chain in self.protein.structure.get_chains():

      chain_id = chain.get_id()
    
      # Iterating through residues of the protein
      for residue in chain:

        if Bio.PDB.is_aa(residue):  # ignore water molecules
        
          # Defining what we call 'res' for the dataframe rows
          resname = residue.get_resname()
          resnum = residue.get_id()[1]
          res =  chain_id + '_' + resname + str(resnum)

          # Now iterating through atoms of each residue to extract properties
          for atom in residue:
            
            # obtaining just the atom name
            at_id = atom.get_id() 

            # making sure the atom name is found in the atom_types_dict variable
            if at_id in self.atom_types_dict[resname]: 

              if 'DON' in self.atom_types_dict[resname][at_id]:
                self.protein.dataframe.loc[res,'don_count'] += 1
              if 'ACP' in self.atom_types_dict[resname][at_id]:
                self.protein.dataframe.loc[res,'acp_count'] += 1
              if 'HPB' in self.atom_types_dict[resname][at_id]:
                self.protein.dataframe.loc[res,'hpb_count'] += 1
              if 'POS' in self.atom_types_dict[resname][at_id]:
                self.protein.dataframe.loc[res,'pos_count'] += 1
              if 'NEG' in self.atom_types_dict[resname][at_id]:
                self.protein.dataframe.loc[res,'neg_count'] += 1
              if 'ARM' in self.atom_types_dict[resname][at_id]:
                self.protein.dataframe.loc[res,'arm_count'] += 1 

######################################################################
  # Is the amino acid a cysteine?
  # Input: -
  # Output: Append the column on the dataframe 
######################################################################

  def is_cysteine(self):

    # Create a new column for Cysteine
    self.protein.dataframe["Cysteine"] = pd.Series(dtype=str)
     
    for chain in self.protein.structure.get_chains():

      chain_id = chain.get_id()
    
    # Iterating through residues of the protein
      for residue in chain:

        resname = residue.get_resname()
        resnum = residue.get_id()[1]
        res =  chain_id + '_' + resname + str(resnum)


        if residue.get_resname() == 'CYS':
          self.protein.dataframe.loc[res,'Cysteine'] = 'Yes'
          
        elif Bio.PDB.is_aa(residue):
          self.protein.dataframe.loc[res,'Cysteine'] = 'No'

        else:
          continue


class Interactions:

  ######################################################################
    # Extract interactions between residues
    # Input: list of chain objects
    # Output: dictionary with interesting features
  ######################################################################
  
  def __init__(self, protein, atom_props_dir):

    self.protein = protein
    self.atom_types_dict = read_atom_types(atom_props_dir)
    self.calculate_interactions()
  
  def calculate_interactions(self):

    # Defining new columns for atom features

    interaction_columns = ['aromatic_stacking', 'hydrogen_bond', 'hydrophobic', 'repulsive', 
                    'attractive']

    # Adding new columns to data frame
    for column in interaction_columns: 
    
      self.protein.dataframe[column] = 0

    # Creating layer dictionary
    layer_dict = {}

    for chain in self.protein.structure.get_chains():

      redundancy_control = []
      chain_id = chain.get_id()

      for residue1 in chain:

        if Bio.PDB.is_aa(residue1):  # ignore water molecules

          # appending residue to redundancy control
          redundancy_control.append(residue1)   

          # Obtaining full residue name
          resname1 = residue1.get_resname()
          resnum1 = residue1.get_id()[1]
          res1 =  chain_id + '_' + resname1 + str(resnum1)

          # Creating empty set of first layer residues wrt current residue
          layer_dict[res1] = set()

          # Obtaining residue neighborhood to calculate distances only with 
          # residues close to current residue
          set_neigh = set()
          set_neigh = self.protein.get_neighborhood(residue1, 8)

          for residue2 in set_neigh:

            if Bio.PDB.is_aa(residue2) and residue2 not in redundancy_control:

              # Get the chain ID of the residue
              chain2_id = residue2.get_parent().id  
              
              # Obtaining full residue name
              resname2 = residue2.get_resname()
              resnum2 = residue2.get_id()[1]
              res2 =  chain2_id + '_' + resname2 + str(resnum2)

              for atom1 in residue1.get_atoms():

                # obtaining just the atom name
                at1_id = atom1.get_id()

                if at1_id not in self.atom_types_dict[resname1]:
                  continue

                for atom2 in residue2.get_atoms():

                  # obtaining just the atom name
                  at2_id = atom2.get_id()

                  distance = self.calculate_distance(atom1, atom2)

                  # making sure the atom name is found in the atom_types_dict variable
                  if at2_id not in self.atom_types_dict[resname2]:
                    continue

                  if ('ARM' in self.atom_types_dict[resname1][at1_id] and 'ARM' in self.atom_types_dict[resname2][at2_id]) and (1.5 <= distance <= 3.5):
                    self.protein.dataframe.loc[res1,'aromatic_stacking'] += 1
                    self.protein.dataframe.loc[res2,'aromatic_stacking'] += 1
                    layer_dict[res1].add(res2)

                  if ('ACP' in self.atom_types_dict[resname1][at1_id] and 'DON' in self.atom_types_dict[resname2][at2_id]) and (2.0 <= distance <= 3.0):
                    self.protein.dataframe.loc[res1,'hydrogen_bond'] += 1
                    self.protein.dataframe.loc[res2,'hydrogen_bond'] += 1
                    layer_dict[res1].add(res2)

                  if ('DON' in self.atom_types_dict[resname1][at1_id] and 'ACP' in self.atom_types_dict[resname2][at2_id]) and (2.0 <= distance <= 3.0):
                    self.protein.dataframe.loc[res1,'hydrogen_bond'] += 1
                    self.protein.dataframe.loc[res2,'hydrogen_bond'] += 1
                    layer_dict[res1].add(res2)

                  if ('HPB' in self.atom_types_dict[resname1][at1_id] and 'HPB' in self.atom_types_dict[resname2][at2_id]) and (2.0 <= distance <= 3.8):
                    self.protein.dataframe.loc[res1,'hydrophobic'] += 1
                    self.protein.dataframe.loc[res2,'hydrophobic'] += 1
                    layer_dict[res1].add(res2)

                  if ('POS' in self.atom_types_dict[resname1][at1_id] and 'POS' in self.atom_types_dict[resname2][at2_id]) and (2.0 <= distance <= 6.0):
                    self.protein.dataframe.loc[res1,'repulsive'] += 1
                    self.protein.dataframe.loc[res2,'repulsive'] += 1
                    layer_dict[res1].add(res2)

                  if ('NEG' in self.atom_types_dict[resname1][at1_id] and 'NEG' in self.atom_types_dict[resname2][at2_id]) and (2.0 <= distance <= 6.0):
                    self.protein.dataframe.loc[res1,'repulsive'] += 1  
                    self.protein.dataframe.loc[res2,'repulsive'] += 1  
                    layer_dict[res1].add(res2)  

                  if ('POS' in self.atom_types_dict[resname1][at1_id] and 'NEG' in self.atom_types_dict[resname2][at2_id]) and (2.0 <= distance <= 6.0):
                    self.protein.dataframe.loc[res1,'attractive'] += 1
                    self.protein.dataframe.loc[res2,'attractive'] += 1
                    layer_dict[res1].add(res2)

                  if ('NEG' in self.atom_types_dict[resname1][at1_id] and 'POS' in self.atom_types_dict[resname2][at2_id]) and (2.0 <= distance <= 6.0):
                    self.protein.dataframe.loc[res1,'attractive'] += 1
                    self.protein.dataframe.loc[res2,'attractive'] += 1  
                    layer_dict[res1].add(res2)

    # Obtaining a data frame from layer dictionary
    self.layer_df = pd.DataFrame.from_dict(layer_dict, orient='index')
    self.layer_df.index.name = 'residue_name'         

  ### CALCULATING DISTANCES BETWEEN ATOMS OF 2 RESIDUES
  def calculate_distance(self, atom1, atom2):

    coord1 = atom1.get_coord()
    coord2 = atom2.get_coord()
  
    distance = math.sqrt ( (coord1[0]- coord2[0])**2 +
                            (coord1[1]- coord2[1])**2 +
                            (coord1[2]- coord2[2])**2 )

    return distance # returning distance
  

######################################################################
  # Reading table of atom types
  # Input: directory to table file
  # Output: dictionary of table information 
######################################################################  
  
def read_atom_types(atom_props_dir):

  atom_types_dict = {}

  with open(atom_props_dir,"r") as atom_types:
    for line in atom_types:
      flist = line.strip().split(sep=",")

      if flist[0] not in atom_types_dict:
        atom_types_dict[flist[0]] = {}

      if len(flist) <= 2:
        atom_types_dict[flist[0]][flist[1]] = []
      else:
        atom_types_dict[flist[0]][flist[1]] = flist[2:]

  return atom_types_dict 