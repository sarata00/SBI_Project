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

    parser = argparse.ArgumentParser(description= ' This program does this')

    parser.add_argument('-p', dest = 'protein_file', default = None, required = True, help = 'Path to input PDB file')
    parser.add_argument('-v', '--verbose', dest = 'verbose', action = 'store_true', default = False, help = 'Print log in stderr') # It saves the var verbose as True, so we can call it if we want to print to Stderr

    options = parser.parse_args()

    if os.path.isfile(options.protein_file) and options.protein_file.endswith('.pdb'):

        target_protein = Protein(options.protein_file)
        
        sys.stderr.write(' -----------------------------------\n') if options.verbose else None
        sys.stderr.write('                                                             \n') if options.verbose else None
        sys.stderr.write('  STARTING PREDICTION OF ' + target_protein.file_name + '    \n') if options.verbose else None
        sys.stderr.write('                                                             \n') if options.verbose else None
        sys.stderr.write(' -----------------------------------\n\n') if options.verbose else None

    else:
        print("PDB file needed")
        exit(3)

    ##########################
    #                        #
    # 1. COMPUTING TEMPLATES #
    #                        #
    ##########################
    sys.stderr.write('--------------------------------------------------\n') if options.verbose else None
    sys.stderr.write('Searching homologous proteins to query sequence...\n') if options.verbose else None
    sys.stderr.write('--------------------------------------------------\n\n') if options.verbose else None

    # Create a folder to store everything related to the templates:
    path = './templates'
    if not os.path.exists(path):
        os.makedirs(path)

    ###########
    # 1.1 BLAST
    ###########

    templ = BLAST('../BioLip_database/biolip_database')

    # Generating a dictionary where we will store the random forest models obtained for each chain
    rf_models_dict = {}

    # Download the fasta files of our target protein (a fasta file per chain)
    for fasta_file in target_protein.write_fasta_file(): 

        sys.stderr.write(f"Working on following chain: {fasta_file}\n\n") if options.verbose else None

        # BLAST of each chain, the homologous proteins are stored in a list
        list_homologs = templ.search_homologs(fasta_file)

        if list_homologs != None:

            list_new = []

            for homolog in list_homologs:

                # Avoid storing the target protein as homolog
                if target_protein.protein_id not in str(homolog):
                    list_new.append(homolog)

            # We just need top 20 homologous proteins
            if len(list_new) >= 20:
                list_homologs_final = list_new[:20]

            ###################################
            # 1.2 COMPUTING TEMPLATE FEATURES #
            ###################################

            if list_homologs_final:
                
                sys.stderr.write('Top 20 homologous sequences found: ') if options.verbose else None
                sys.stderr.write(f"[{', '.join(list_homologs_final)}]\n") if options.verbose else None        

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
                    template.get_pdb_chains() # this method produces the chain_directory mentioned below
                    pdb_id = template.structure.get_id()

                    # Retrieving folder name for this template
                    chain_directory = template.chains_folder_name
                    #print(chain_directory)

                    if os.path.exists(chain_directory):
                        # now checking if current chain name is in the directory as a pdb chain file
                        template_chain = ''

                        for file_name in os.listdir(chain_directory):
                            if pdb in file_name:

                                sys.stderr.write(f"PDB file for {pdb} could be obtained. Now calculating chain features...\n") if options.verbose else None
                                template_chain = Protein(str(chain_directory + '/' + file_name))
                        
                        template_chain_dataset = RandomForestModel().get_training_data(template_chain)

                        if type(template_chain_dataset) == str:
                            sys.stderr.write(template_chain_dataset) if options.verbose else None

                        # we delete folder with chains because we want to avoid accumulation of folders
                        try:
                            shutil.rmtree(chain_directory)
                            #print(f"Directory {chain_directory} was successfully removed.")
                        except OSError as e:
                            print(f"Error: {chain_directory} : {e.strerror}")

                    else:
                        print(f"The directory {chain_directory} could not be created, current template will be dismissed.\n") if options.verbose else None
                        continue
                    
                    # Appending each data frame 
                    # (and avoid including templates without binding site labels)
                    if type(template_chain_dataset) != str:
                        dataframes_list.append(template_chain_dataset)

                # Concatenating data sets for each chain
                all_templates_dataset = RandomForestModel().concat_training_data(dataframes_list)
                sys.stderr.write(f"\nTraining set for current chain:\n") if options.verbose else None
                sys.stderr.write(str(all_templates_dataset)) if options.verbose else None
                sys.stderr.write('\n') if options.verbose else None

                # Splitting data into train/test
                X_train, X_test, Y_train, Y_test = RandomForestModel().split_data(all_templates_dataset)
                #print(X_train, X_test, Y_train, Y_test)
                sys.stderr.write(f"Data successfully split into train and test sets!\n") if options.verbose else None
                #sys.stderr.write(f"{X_train}\n{Y_train}")

                ###################################
                # 1.3 GENERATING CLASSIFIER MODEL #
                ###################################

                # Generating classifier model using Random Forest
                rf_model = RandomForestModel().get_model(X_train, Y_train.astype('int'))
                # Appending chain model to dictionary for later use
                rf_models_dict[fasta_file] = rf_model

                ###############################
                # 1.4 EVALUATION ON TEST DATA #
                ###############################

                RandomForestModel().get_predictions(X_test, rf_model)
                Y_pred = X_test.iloc[:, -1]

                # Evaluating model
                metrics_table = RandomForestModel().model_evaluation(Y_pred, Y_test)
                
                # Create a list of tuples with the metric names and values
                if metrics_table:
                    print("The following metrics were obtained from the test data:\n")
                    print(metrics_table, '\n')
                else:
                    print("WARNING: metrics could not be obtained.")
                    
            else:
                print(f"Final homologs list could not be retrieved from {fasta_file}, stop the program.")
                exit(3)
        else:

            print(f"No homologs could be retrieved from {fasta_file}, stop the program")
            exit(3)

    ##########################
    #                        #
    # 3. QUERY FEATURES      #
    #                        #
    ##########################

    sys.stderr.write('----------------------------------------------------------\n') if options.verbose else None
    sys.stderr.write('Calculating query protein features and making predictions\n') if options.verbose else None
    sys.stderr.write('------------------------------------------------------------\n\n') if options.verbose else None

     # Create a list to store the binding site residues (later on will be used for Chimera)
    binding_sites_chimera = []
    
    ##########################################
    # 3.1 FEATURES FOR EACH CHAIN SEPARATELY #
    ##########################################
    
    # Similar methodology applied in the templates section above
    target_protein.get_pdb_chains() # this produces the chain_directory mentioned below
    query_pdb_id = target_protein.structure.get_id()

    query_chain_directory = target_protein.chains_folder_name

    if os.path.exists(query_chain_directory):

        # (PARENTHESIS) - Creating file for chimera (see much below)
        try:
            chimera_cmd_file = open("chimera.cmd", "w")
        except FileExistsError:
            print("File already exists.") # realment aixo no cal pero no se q posar

        sys.stderr.write(f"We have created the folder {query_chain_directory} to store the chains of the query protein {query_pdb_id}\n")

        # Naming directory where we will store predictions
        out_dir = './query_predictions/'

        # now checking if current chain name is in the directory as a pdb chain file
        for file_name in os.listdir(query_chain_directory):
            sys.stderr.write(f"CALCULATING CHAIN FEATURES\n\n") if options.verbose else None
            sys.stderr.write(f"PDB file for {file_name[:-4]} could be obtained.\n\n") if options.verbose else None
            query_chain = Protein(str(query_chain_directory + '/' + file_name))  

            try:
                # Computing residue and atom features
                query_chain_features = ProteinFeatures(query_chain, './atom_types.csv')
                query_chain_features.residue_properties()
                query_chain_features.atom_properties()
                query_chain_features.is_cysteine()

                # Computing interactions
                query_chain_interactions = Interactions(query_chain, './atom_types.csv')
                query_chain_interactions.calculate_interactions()

                # Compute layer features (atom and residue properties, and interactions)
                query_chain_layer = Layer(query_chain, './atom_types.csv')
                query_chain_layer.get_layer_properties()

            except Exception as e:
                print("Error:", e, file=sys.stderr)
                raise
            
            #sys.stderr.write('Query protein features calculated successfully.\n') if options.verbose else None

            #print(query_chain.dataframe)

            ##########################
            #                        #
            #    4. PREDICTIONS      #
            #                        #
            ##########################
            
            sys.stderr.write(f'NOW PREDICTING BINDING SITES ON CURRENT QUERY PROTEIN CHAIN\n') if options.verbose else None

            # Checking if name of current chain in keys of rf_models_dict
            for q_chain in rf_models_dict:

                if file_name[:-4] in q_chain:

                    predictions = RandomForestModel().get_predictions(query_chain.dataframe, rf_models_dict[q_chain])
                    if predictions:

                        sys.stderr.write(f'Binding sites have been successfully predicted!\n\n') if options.verbose else None

                    if not os.path.exists('./query_predictions/'):
                        os.makedirs('./query_predictions/')

                    query_chain.dataframe.to_csv(out_dir + query_chain.structure.get_id() + '.csv', columns=['prediction'])

            ############################
            #                          #
            # 5. GENERATE CHIMERA FILE #
            #                          #
            ############################

            sys.stderr.write(str(query_chain.dataframe)) if options.verbose else None

            for index, row in query_chain.dataframe.iterrows():
            # Extract the residue name, number and chain
                
                res = index
                resnum = res[5:]
                chain = res[0]

                # Loop through the remaining columns and select the residues
                label = row['prediction']
                
                if int(label) == 1:
                    
                    # format must be for example 'select :10.A'
                    binding_sites_chimera.append(f'{resnum}.{chain}') 
                    

        selection = ','.join(binding_sites_chimera)            
        #chimera_cmd_file.write(f'select :{selection}\n')
        chimera_cmd_file.write('# Coloring binding sites\n')
        chimera_cmd_file.write(f'color green :{selection}')

        chimera_cmd_file.close()

        sys.stderr.write(f'Prediction tables can be found in {out_dir} \n') if options.verbose else None

if __name__ == '__main__':
    main()


 