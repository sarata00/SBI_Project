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
import shutil # remove directory with files


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
        
        sys.stderr.write(' -----------------------------------\n') if options.verbose else None
        sys.stderr.write('                                                             \n') if options.verbose else None
        sys.stderr.write('  STARTING PREDICTION OF ' + target_protein.file_name + '    \n') if options.verbose else None
        sys.stderr.write('                                                             \n') if options.verbose else None
        sys.stderr.write(' -----------------------------------\n\n') if options.verbose else None

    else:
        print("PDB file needed")

    ##########################
    #                        #
    # 1. COMPUTING TEMPLATES #
    #                        #
    ##########################
    sys.stderr.write('--------------------------------------------------\n') if options.verbose else None
    sys.stderr.write('Searching homologous proteins to query sequence...\n') if options.verbose else None
    sys.stderr.write('--------------------------------------------------\n\n') if options.verbose else None

    # Create a folder to store everything related to the templates:
    path = "./templates"
    if not os.path.exists(path):
        os.makedirs(path)

    ###########
    # 1.1 BLAST
    ###########

    templ = BLAST('../BioLip_database/biolip_database')

    # Download our target fasta file
    for fasta_file in target_protein.write_fasta_file(): 

        sys.stderr.write(f"Working on following chain: {fasta_file}\n\n") if options.verbose else None

        list_homologs = templ.search_homologs(fasta_file)

        if list_homologs != None:

            list_new = []

            for homolog in list_homologs:

                if target_protein.protein_id not in str(homolog):
                    list_new.append(homolog)

            if len(list_new) >= 20:
                list_homologs_final = list_new[:20]

            #################################
            # 1.2 COMPUTING TEMPLATE FEATURES
            #################################

            if list_homologs_final:
                
                sys.stderr.write('Top 20 homologous sequences found:') if options.verbose else None
                sys.stderr.write(f"{' '.join(list_homologs_final)}\n") if options.verbose else None        

                # Extract the PDBs of each template
                pdbl = PDBList()

                # Creating empty list to store data frames
                dataframes_list = []

                for pdb in list_homologs_final:

                    # Extract the first four characters of each PDB ID
                    pdb_id = pdb[:4]
                
                    # Download the PDB file in .ent format
                    f = pdbl.retrieve_pdb_file(pdb_id, pdir = path, file_format="pdb")
                    
                    # Obtain chains from PDB file
                    template = Protein(f)
                    template.get_pdb_chains() # this produces the chain_directory mentioned below
                    pdb_id = template.structure.get_id()

                    chain_directory = str('./' + pdb_id + '_chains')
                    sys.stderr.write(f"We have created the folder {chain_directory} to store the chains of {pdb_id}\n")
                    
                    # now checking if current chain name is in the directory as a pdb chain file
                    template_chain = ''
                    for file_name in os.listdir(chain_directory):
                        if pdb in file_name:

                            sys.stderr.write(f"PDB file for {pdb} could be obtained. Now calculating chain features...\n")
                            template_chain = Protein(str(chain_directory + '/' + file_name))
                    
                    template_dataset = RandomForestModel().get_training_data(template_chain)

                    # we delete folder with chains because we want to avoid accumulation of folders
                    try:
                        shutil.rmtree(chain_directory)
                        #print(f"Directory {chain_directory} was successfully removed.")
                    except OSError as e:
                        print(f"Error: {chain_directory} : {e.strerror}")

                    # Appending each data frame
                    dataframes_list.append(template_dataset)

                # Concatenating data sets for each chain
                all_templates_dataset = RandomForestModel().concat_training_data(dataframes_list)
                print(all_templates_dataset)

                # Splitting data into train/test
                X_train, X_test, Y_train, Y_test = RandomForestModel().split_data(all_templates_dataset)

                # Generating classifier model using Random Forest
                rf_model = RandomForestModel().get_model(X_train, Y_train)
                print(rf_model)

                    #print(template_dataset)



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
    sys.stderr.write('\n-----------------------------------\n') if options.verbose else None
    sys.stderr.write('Calculating query protein features...\n') if options.verbose else None
    sys.stderr.write('-----------------------------------\n\n') if options.verbose else None

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


 