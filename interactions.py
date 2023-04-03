import pandas as pd
import os
import sys
import Bio.PDB
from Bio.SeqUtils.ProtParam import ProteinAnalysis, ProtParamData
from Bio.PDB import NeighborSearch
from protein import Protein, ProteinFeatures
import math

class Interactions:

  ######################################################################
    # Extract interactions between residues
    # Input: list of chain objects
    # Output: dictionary with interesting features
  ######################################################################
  
  def __init__(self, protein):

    self.protein = protein
    self.atom_types_dict = read_atom_types(atom_props_dir)
  
  def calculate_residue_interactions(self):

    # Defining new columns for atom features

    interaction_columns = ['aromatic_stacking', 'hydrogen_bond', 'hydrophobic', 'repulsive', 
                    'attractive']

    # Adding new columns to data frame
    for column in interaction_columns: 
    
      self.protein.dataframe[column] = 0

    for chain in self.protein.get_chains():

      redundancy_control = []
      chain_id = chain.get_id()

      for residue1 in chain:
        
        # Obtaining full residue name
        resname1 = residue1.get_resname()
        resnum1 = residue1.get_id()[1]
        res1 =  chain_id + '_' + resname1 + str(resnum1)

        redundancy_control.append(residue1)

        set_neigh = set()
        set_neigh = self.protein.get_neighborhood(residue, 8)

        for residue2 in set_neigh:

          # Get the chain ID of the residue
          chain2_id = residue2.get_parent().id  
          
          # Obtaining full residue name
          resname2 = residue2.get_resname()
          #resnum2 = residue2.get_id()[1]
          #res2 =  chain2_id + '_' + resname2 + str(resnum2)

          for atom1 in residue1.get_atoms():

            # obtaining just the atom name
            at1_id = atom.get_id()

            for atom2 in residue2.get_atoms():

              # obtaining just the atom name
              at2_id = atom.get_id()

              distance = self.calculate_distance(atom1, atom2)

              if ('ARM' in self.atom_types_dict[resname1][at1_id] and 'ARM' in self.atom_types_dict[resname2][at2_id]) and (1.5 <= distance <= 3.5):
                self.protein.dataframe.loc[res1,'aromatic_stacking'] += 1

              if ('ACP' in self.atom_types_dict[resname1][at1_id] and 'DON' in self.atom_types_dict[resname2][at2_id]) and (2.0 <= distance <= 3.0):
                self.protein.dataframe.loc[res1,'hydrogen_bond'] += 1

              if ('DON' in self.atom_types_dict[resname1][at1_id] and 'ACP' in self.atom_types_dict[resname2][at2_id]) and (2.0 <= distance <= 3.0):
                self.protein.dataframe.loc[res1,'hydrogen_bond'] += 1

              if ('HPB' in self.atom_types_dict[resname1][at1_id] and 'HPB' in self.atom_types_dict[resname2][at2_id]) and (2.0 <= distance <= 3.8):
                self.protein.dataframe.loc[res1,'hydrophobic'] += 1

              if ('POS' in self.atom_types_dict[resname1][at1_id] and 'POS' in self.atom_types_dict[resname2][at2_id]) and (2.0 <= distance <= 6.0):
                self.protein.dataframe.loc[res1,'repulsive'] += 1

              if ('NEG' in self.atom_types_dict[resname1][at1_id] and 'NEG' in self.atom_types_dict[resname2][at2_id]) and (2.0 <= distance <= 6.0):
                self.protein.dataframe.loc[res1,'repulsive'] += 1   

              if ('POS' in self.atom_types_dict[resname1][at1_id] and 'NEG' in self.atom_types_dict[resname2][at2_id]) and (2.0 <= distance <= 6.0):
                self.protein.dataframe.loc[res1,'attractive'] += 1

              if ('NEG' in self.atom_types_dict[resname1][at1_id] and 'POS' in self.atom_types_dict[resname2][at2_id]) and (2.0 <= distance <= 6.0):
                self.protein.dataframe.loc[res1,'attractive'] += 1           

  ### CALCULATING DISTANCES BETWEEN ATOMS OF 2 RESIDUES
  def calculate_distance(self, atom1, atom2):

    coord1 = atom1.get_coord()
    coord2 = atom2.get_coord()
  
    distance = math.sqrt ( (coord1[0]- coord2[0])**2 +
                            (coord1[1]- coord2[1])**2 +
                            (coord1[2]- coord2[2])**2 )

    return distance # returning distance