import pandas as pd
import os
import sys
import Bio.PDB
from Bio.SeqUtils.ProtParam import ProteinAnalysis, ProtParamData
from Bio.PDB import NeighborSearch

# Class
class Protein:
  
  def __init__(self, pdbfile):
    self.file_name = pdbfile
    self.protein_id = pdbfile.split('/')[-1].split('.')[0]
    self.structure = self.read_pdb(self.protein_id, self.file_name) # read pdb
    self.dataframe_info()

################################################################################
  # Function to read the pdb file and obtain its Structure object.

  # Input: PDB code, PDB file
  # Output: Structure Object
################################################################################

  def read_pdb(self, pdbcode, pdbfilenm):
    """
    Read a PDB structure from a file.
    :param pdbcode: A PDB ID string
    :param pdbfilenm: The PDB file
    :return: a Bio.PDB.Structure object or None if something went wrong
    """
    try:
        pdbparser = Bio.PDB.PDBParser(QUIET=True)   # suppress PDBConstructionWarning
        struct = pdbparser.get_structure(pdbcode, pdbfilenm)
        return struct
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None 

######################################################################
  # Get the sequences for each chain
  # Input: list of chain objects
  # Output: dictionary with chain_ids and its sequences
######################################################################
  
  def get_sequence(self):
    
    dict_seq = {}
    for chain in self.structure.get_chains():
      chain_id = chain.get_id()
      sequence = ''
      
      if chain not in dict_seq:
        dict_seq[chain_id] = {}

      for residue in chain:
        if residue.get_resname() != 'HOH':  # ignore water molecules
          sequence += Bio.SeqUtils.seq1(residue.get_resname())

      dict_seq[chain_id] = sequence 
    
    return dict_seq
      
    

######################################################################
  # Gets the fasta sequence of each chain
  # Input: dictionary with chain ids and its sequences
  # Output: respective .fasta file for each chain in .pdb file
######################################################################
  
  def write_fasta_file(self):
    directory = os.getcwd() # no se si aquí añadir por si nos dan un directorio
    file_id = os.path.basename(self.file_name).split('.')[0]
    fasta_files_list = []
    
    for chain, seq in self.get_sequence().items():
      fasta_filename = os.path.join(directory, f"{file_id}_{chain}.fasta")
    
      with open(fasta_filename, 'w') as fasta_file:
      
        fasta_file.write(f">{file_id}_{chain}\n")
        fasta_file.write(f"{seq}\n")
      
      fasta_files_list.append(f'./{file_id}_{chain}.fasta')

    return fasta_files_list


################################################################################
  # Create a dataframe to store all the information retrieved for the analysis.
  # The features (columns) will be the info retrieved, the records (rows) will
  # be the residues of the protein.

  # Input: -
  # Output: Empty dataframe
################################################################################
 
  def dataframe_info(self):

    # List all residues
    residue_list = list(self.structure.get_residues())
    residue_names = []

    for residue in residue_list:

      # Iterate over the structure to obtain the Residue Name
      for chain in self.structure.get_chains():

        chain_id = chain.get_id()
      
        for resid in chain:
        
         # Only consider amino acids, and compute it when residue of the list = residue of the structure
         if Bio.PDB.is_aa(residue) and residue == resid: 

          resname = residue.get_resname()
          resnum = residue.get_id()[1]

          res =  chain_id + '_' + resname + str(resnum)
          residue_names.append(res)

    # Create empty dataframe
    self.dataframe = pd.DataFrame(index=residue_names)
    self.dataframe.index.name = 'residue_name'
  
################################################################################
  # Define neighborhood of residues.
  # Input: -
  # Output: list of neighbor residues
################################################################################

  def get_neighborhood(self, residue, radius):

    list_neighborhood = []
    neighborhood_set = set()

    # Using NeighborSearch from BioPython to define neighborhood given structure atoms
    ns = NeighborSearch(list(self.structure.get_atoms()))

    for atom in residue:
      neighborhood = ns.search(atom.get_coord(), radius, level='R')
      list_neighborhood += neighborhood
      #print(residue, atom, list_neighborhood)
    
    # Avoid repeated neighbors and the residue itself. Also avoid H20 molecules
    for neighbor in list_neighborhood:

      if (neighbor not in neighborhood_set and neighbor != residue) and (Bio.PDB.is_aa(residue) and Bio.PDB.is_aa(neighbor)):

            neighborhood_set.add(neighbor)
    
    return neighborhood_set




#### AUXILIARY FUNCTION ########
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

# Class
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

              #dict_properties[res] = [hydrophobicity_value, surface_acc_value]

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

        # Setting counters to 0 for each new residue
        don_count = 0
        acp_count = 0
        hpb_count = 0
        pos_count = 0
        neg_count = 0
        arm_count = 0

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
                don_count +=1
                self.protein.dataframe.loc[res,'don_count'] += 1
              if 'ACP' in self.atom_types_dict[resname][at_id]:
                acp_count +=1
                self.protein.dataframe.loc[res,'acp_count'] += 1
              if 'HPB' in self.atom_types_dict[resname][at_id]:
                hpb_count +=1
                self.protein.dataframe.loc[res,'hpb_count'] += 1
              if 'POS' in self.atom_types_dict[resname][at_id]:
                pos_count +=1
                self.protein.dataframe.loc[res,'pos_count'] += 1
              if 'NEG' in self.atom_types_dict[resname][at_id]:
                neg_count +=1 
                self.protein.dataframe.loc[res,'neg_count'] += 1
              if 'ARM' in self.atom_types_dict[resname][at_id]:
                arm_count +=1
                self.protein.dataframe.loc[res,'arm_count'] += 1  


              

######################################################################
  # The amino acid is a cysteine?
  # Input: -
  # Output: Append the column on the dataframe 
######################################################################

  def is_cysteine(self):

    # Create a new column for Cysteine
    self.protein.dataframe["Cysteine"] = pd.Series(dtype=int)
     
    for chain in self.protein.structure.get_chains():

      chain_id = chain.get_id()
    
    # Iterating through residues of the protein
      for residue in chain:

        resname = residue.get_resname()
        resnum = residue.get_id()[1]
        res =  chain_id + '_' + resname + str(resnum)


        if residue.get_resname() == 'CYS':
          self.protein.dataframe.loc[res,'Cysteine'] = '1'
          
        elif Bio.PDB.is_aa(residue):
          self.protein.dataframe.loc[res,'Cysteine'] = '0'

        else:
          continue