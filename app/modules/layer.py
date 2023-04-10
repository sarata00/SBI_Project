import pandas as pd
import Bio.PDB
from features_interactions import Interactions

class Layer:

  def __init__(self, protein, atom_props_dir):
    self.protein = protein
    self.layer_df = Interactions(protein, atom_props_dir).layer_df

  def get_layer_properties(self):

    # Generating new columns
    lresidue_cols = ['L_hydrophob_value', 'L_surf_acc', 'L_don_count', 'L_acp_count', 'L_hpb_count', 'L_pos_count', 
                    'L_neg_count', 'L_arm_count', 'L_aromatic_stacking', 'L_hydrogen_bond', 'L_hydrophobic', 'L_repulsive', 
                    'L_attractive']

    # Adding new columns to data frame
    for column in lresidue_cols: 
    
      self.protein.dataframe[column] = 0

    # Calculating properties
    # Iterating through residues in layer df
    for residue, row in self.layer_df.iterrows():

      # initializing variables
      # number of neighbors of current residue
      num_neighbors = 0
      # residue properties
      hydrophob_sum = 0
      # atom properties
      surf_acc_sum = 0
      don_count_sum = 0
      acp_count_sum = 0
      hpb_count_sum = 0
      pos_count_sum = 0
      neg_count_sum = 0
      arm_count_sum = 0
      # interactions
      aromatic_stacking_sum = 0
      hydrogen_bond_sum = 0
      hydrophobic_sum = 0
      repulsive_sum = 0
      attractive_sum = 0
      
      # Iterating through neighbors of each residue
      for colIndex, neighbor in row.items():

        if neighbor == None:
          continue

        num_neighbors += 1

        # Adding up values of neighbors in protein.dataframe for each property
        # residue properties
        hydrophob_sum += self.protein.dataframe.at[neighbor, 'hydrophob_value']
        surf_acc_sum += self.protein.dataframe.at[neighbor, 'surf_acc']
        # atom properties
        don_count_sum += self.protein.dataframe.at[neighbor, 'don_count']
        acp_count_sum += self.protein.dataframe.at[neighbor, 'acp_count']
        hpb_count_sum += self.protein.dataframe.at[neighbor, 'hpb_count']
        pos_count_sum += self.protein.dataframe.at[neighbor, 'pos_count']
        neg_count_sum += self.protein.dataframe.at[neighbor, 'neg_count']
        arm_count_sum += self.protein.dataframe.at[neighbor, 'arm_count']
        # interactions
        aromatic_stacking_sum += self.protein.dataframe.at[neighbor, 'aromatic_stacking']
        hydrogen_bond_sum += self.protein.dataframe.at[neighbor, 'hydrogen_bond']
        hydrophobic_sum += self.protein.dataframe.at[neighbor, 'hydrophobic']
        repulsive_sum += self.protein.dataframe.at[neighbor, 'repulsive']
        attractive_sum += self.protein.dataframe.at[neighbor, 'attractive']

      # Weighting values by number of neighbors, if any
      # if no neighbors, no values will be calculated and layer properties will 
      # be 0 as set above
      if num_neighbors > 0:
      
        # residue properties
        hydrophob_avg = hydrophob_sum/num_neighbors
        surf_acc_avg = surf_acc_sum/num_neighbors
        # atom properties
        don_count_avg = don_count_sum/num_neighbors
        acp_count_avg = acp_count_sum/num_neighbors
        hpb_count_avg = hpb_count_sum/num_neighbors
        pos_count_avg = pos_count_sum/num_neighbors
        neg_count_avg = neg_count_sum/num_neighbors
        arm_count_avg = arm_count_sum/num_neighbors
        # interactions
        aromatic_stacking_avg = aromatic_stacking_sum/num_neighbors
        hydrogen_bond_avg = hydrogen_bond_sum/num_neighbors
        hydrophobic_avg = hydrophobic_sum/num_neighbors
        repulsive_avg = repulsive_sum/num_neighbors
        attractive_avg = attractive_sum/num_neighbors

        # Adding values in protein.dataframe columns of layer
        # residue properties
        self.protein.dataframe.loc[residue,'L_hydrophob_value'] = hydrophob_avg
        self.protein.dataframe.loc[residue,'L_surf_acc'] = surf_acc_avg
        # atom properties
        self.protein.dataframe.loc[residue,'L_don_count'] = don_count_avg
        self.protein.dataframe.loc[residue,'L_acp_count'] = acp_count_avg
        self.protein.dataframe.loc[residue,'L_hpb_count'] = hpb_count_avg
        self.protein.dataframe.loc[residue,'L_pos_count'] = pos_count_avg
        self.protein.dataframe.loc[residue,'L_neg_count'] = neg_count_avg
        self.protein.dataframe.loc[residue,'L_arm_count'] = arm_count_avg
        # interactions
        self.protein.dataframe.loc[residue,'L_aromatic_stacking'] = aromatic_stacking_avg
        self.protein.dataframe.loc[residue,'L_hydrogen_bond'] = hydrogen_bond_avg
        self.protein.dataframe.loc[residue,'L_hydrophobic'] = hydrophobic_avg
        self.protein.dataframe.loc[residue,'L_repulsive'] = repulsive_avg
        self.protein.dataframe.loc[residue,'L_attractive'] = attractive_avg