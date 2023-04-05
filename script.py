import pandas as pd
import os
import sys
import argparse
import subprocess
import Bio.PDB
from Bio.SeqUtils.ProtParam import ProteinAnalysis, ProtParamData
from Bio.PDB import NeighborSearch
from Bio.PDB import PDBList
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

        sys.stderr.write('Starting prediction of ' + target_protein.file_name + '\n') if options.verbose else None

    else:
        print("PDB file needed")

    # 1. COMPUTING TEMPLATES
    sys.stderr.write('Searching homologous templates to query protein...\n') if options.verbose else None

    # Create a folder to store everything related to the templates:
    path = "./templates"
    if not os.path.exists(path):
        os.makedirs(path)
    
    # 1.1 BLAST
    templ = BLAST('../BioLip_database/biolip_database')

    # Download our target fasta file
    for fasta_file in target_protein.write_fasta_file(): 

        print(f"New fasta_file: {fasta_file}")

        list_homologs = templ.search_homologs(fasta_file)

        if list_homologs != None:

            list_new = []

            for homolog in list_homologs:

                if target_protein.protein_id not in str(homolog):
                    list_new.append(homolog)

            if len(list_new) >= 20:
                list_homologs_final = list_new[:20]


            # 1.2 COMPUTING TEMPLATE FEATURES
            if list_homologs_final:
                
                print(list_homologs_final)            

            # Extract the PDBs of each template
                pdbl = PDBList()

                for pdb in list_homologs_final:
                    # Extract the first four characters of each PDB ID
                    pdb_id = pdb[:4]
                    print(pdb)
                
                    # Download the PDB file in .ent format
                    f = pdbl.retrieve_pdb_file(pdb_id, pdir = path, file_format="pdb")
                    
                    # Obtain chains from PDB file
                    template = Protein(f)
                    template.get_pdb_chains()

                    if pdb in 

                    template_dataset = RandomForestModel().get_training_data(template)
                    print(template_dataset)

                    # Convert the .ent file to .pdb format
                    #input_path = os.path.join("templates", f"pdb{pdb_id}.ent")
                    #output_path = os.path.join("templates", f"pdb{pdb_id}.pdb")
                    #subprocess.call(["pdb4amber", "-i", input_path, "-o", output_path])
            
                    # Create a Protein object for each template
                    #p = Protein(f)

        else:

            print(f"No homologs could be retrieved from {fasta_file}, stop the program")
            exit(3)

    # MAYBE PUT AN IF HERE IN TERMS OF IF IT HAS NOT FOUND ANY TEMPLATES TO STOP PROGRAM

    # 3. COMPUTING QUERY FEATURES
    sys.stderr.write('Calculating query protein features...\n') if options.verbose else None

    try:
        # Computing residue and atom features
        target_protein_features = ProteinFeatures(target_protein, './atom_types.csv')
        target_protein_features.residue_properties()
        target_protein_features.atom_properties()
        target_protein_features.is_cysteine()

        # Computing interactions
        target_protein_interactions = Interactions(target_protein, './atom_types.csv')
        target_protein_interactions.calculate_interactions()

        # Compute layer features (atom and residue properties, and interactions)
        target_protein_layer = Layer(target_protein, './atom_types.csv')
        target_protein_layer.get_layer_properties()

    except Exception as e:
        print("Error:", e, file=sys.stderr)
        raise
    
    sys.stderr.write('Query protein features calculated successfully.\n') if options.verbose else None

    print(target_protein.dataframe)
    
    #print(target_protein.dataframe)

    # Test the Random Forest model
    #rf = RandomForestModel("./templates")
    #template_dataframe = rf.get_training_data(p)

    #print(f"This is the template_df:\n {template_dataframe}")
   # print(f"Starting with predictions .......")
    #templ_predictions = rf.get_predictions(p)
    #y_test = rf.split_data(template_dataframe).y_test
    #accuracy = rf.model_evaluation(templ_predictions, y_test).accuracy

    #print(f"Model accuracy: {accuracy}")

    # Now depending on the accuracy of the model, use it to 
    # predict the binding sites of our target:

    # target_model = rf.get_training_data(target_protein)
    # target_prediction = rf.get_predictions(target_protein)

    


    
if __name__ == '__main__':
    main()


 