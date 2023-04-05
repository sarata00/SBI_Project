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

    if '/' in pdbfile:
      self.protein_id = pdbfile.split('/')[-1].split('.')[0]
    else:
      self.protein_id = pdbfile.split('.')[0]
      
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