import pandas as pd
import os
import sys
import argparse
import Bio.PDB
from Bio.SeqUtils.ProtParam import ProteinAnalysis, ProtParamData
from Bio.PDB import NeighborSearch
from protein import Protein
from properties import ProteinFeatures, Interactions, Layer
from BLAST import BLAST
from model import RandomForestModel, ExtractBindingSites


def main():
    # COMMAND FROM COMMAND-LINE
    ## name_program.py -p < protein_pdb > -db < biolip_db > (-a < atom_types_file -o < out_file > ??)

    parser = argparse.ArgumentParser(description= " This program does this")

    parser.add_argument('-p', dest = 'protein_file', default = None, required = True, help = 'Path to input PDB file')
    #parser.add_argument('-db', dest = 'BioLip_db', default = None, required = True, help = 'Path to the Biolip database, which is used to do the BLAST')
    #parser.add_argument('-a', dest = 'atom_types_file', default = None, required = True, help = 'Path to the atom_types.csv file, which stores the atom properties')
    # should we create one for the output ?
    parser.add_argument('-v', '--verbose', dest = 'verbose', action = 'store_true', default = False, help = 'Print log in stderr') # It saves the var verbose as True, so we can call it if we want to print to Stderr

    options = parser.parse_args()

    #print(options.verbose)
    #print(options.protein_file)  

    if os.path.isfile(options.protein_file) and options.protein_file.endswith('.pdb'):

        target_protein = Protein(options.protein_file)

        if options.verbose:
            sys.stderr.write('Starting prediction of ' + target_protein.file_name + '\n')
    
    else:
        print("PDB file needed")
    
    # Computing residue and atom properties:

    target_protein_features = ProteinFeatures(target_protein, './atom_types.csv')
    target_protein_features.residue_properties()
    target_protein_features.atom_properties()
    target_protein_features.is_cysteine()

    # Computing interactions:

    target_protein_interactions = Interactions(target_protein, './atom_types.csv')
    target_protein_interactions.calculate_interactions()
    
    print(target_protein.dataframe)

    # Compute layer features (atom and residue properties, and interactions)
    target_protein_layer = Layer(target_protein, './atom_types.csv')
    target_protein_layer.get_layer_properties()


    # Computing templates:
    templ = BLAST('./myproteindb')

    for fasta_file in target_protein.write_fasta_file(): 

        print("New fasta_file")
        list_homologs = templ.search_homologs(fasta_file)
        list_new = []

        for homolog in list_homologs:

            if target_protein.protein_id not in str(homolog):

                list_new.append(homolog)

        if len(list_new) >= 20:

            list_homologs_final = list_new[:20]

        print(list_homologs_final)


        # WE HAVE TO RETRIEVE THE PDB FROM PDB???
    


    





    
if __name__ == '__main__':
    
    main()


 