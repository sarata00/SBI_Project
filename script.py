# Importing modules
import os
import sys
import argparse
import shutil
from Bio.PDB import PDBList

from app.protein import Protein
from app.properties import ProteinFeatures, Interactions, Layer
from app.BLAST import BLAST
from app.model import RandomForestModel, ExtractBindingSites

# Defining main() function for program to be executed only if name = main
def main():

    # COMMAND FROM COMMAND-LINE
    parser = argparse.ArgumentParser(description= 'Residue-centered machine learning algorithm that aims to predict binding sites of a given query protein using a training set of homologous proteins of known binding site information.')

    # adding options to command line program execution
    # -p option is for specifying query protein pdb file
    parser.add_argument('-p', dest = 'protein_file', default = None, required = True, help = 'Path to input PDB file')
    # -v option is for being verbose
    parser.add_argument('-v', '--verbose', dest = 'verbose', action = 'store_true', default = False, help = 'Print log in stderr') # It saves the var verbose as True, so we can call it if we want to print to Stderr

    options = parser.parse_args()

    # defining condition to obtain protein only if it ends with .pdb
    if os.path.isfile(options.protein_file) and options.protein_file.endswith('.pdb'):
        # defining input protein as protein() object
        target_protein = Protein(options.protein_file)
        
        sys.stdout.write(' -----------------------------------\n')
        sys.stderr.write('                                                             \n')
        sys.stderr.write('  STARTING PREDICTION OF ' + target_protein.file_name + '    \n')
        sys.stderr.write('                                                             \n')
        sys.stderr.write(' -----------------------------------\n\n')
    # exit program if input is not pdb file
    else:
        sys.stderr.write("PDB file needed")
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

    # defining variable as BLAST() class instance
    templ = BLAST('../BioLip_database/biolip_database')   

    # Generating a dictionary where we will store the random forest models obtained for each chain
    rf_models_dict = {}

    # Download the fasta files of our target protein (a fasta file per chain)
    for fasta_file in target_protein.write_fasta_file(): 

        chain_name = fasta_file[2:7]

        sys.stderr.write(f"Working on following chain: {fasta_file[2:7]}\n")

        # BLAST of each chain, the homologous proteins are stored in a list
        list_homologs = templ.search_homologs(fasta_file)
        sys.stderr.write(f"\nBLAST results can be found in {templ.out_file}\n\n") if options.verbose else None     

        # if homologs found iterate through each of the homologs (=templates)
        if list_homologs != None:

            list_new = []

            for homolog in list_homologs:

                # Avoid storing the target protein as homolog
                if target_protein.protein_id not in str(homolog):
                    list_new.append(homolog)

            # We will just use top 40 homologous proteins if there are so many
            if len(list_new) >= 40:
                list_homologs_final = list_new[:40]
            else:
                list_homologs_final = list_new

            ###################################
            # 1.2 COMPUTING TEMPLATE FEATURES #
            ###################################

            # if list of final homologs obtained, we proceed to
            # 3. 
            if list_homologs_final:
                
                sys.stderr.write(f'Top {len(list_homologs_final)} homologous sequences found: ') if options.verbose else None
                sys.stderr.write(f"[{', '.join(list_homologs_final)}]\n\n") if options.verbose else None        

                # Extract the PDBs of each template
                pdbl = PDBList(verbose=False)

                # Creating empty list to store data frames
                dataframes_list = []

                for pdb in list_homologs_final:

                    # Extract the first four characters of each PDB ID
                    pdb_id = pdb[:4]

                    # 1. retrieve template whole protein pdb chains files and assign them as Protein instances
                    # Download the PDB file in .ent format
                    f = pdbl.retrieve_pdb_file(pdb_id, pdir = path, file_format="pdb")

                    if not os.path.isfile(f):
                        sys.stderr.write(f'{pdb} template ignored. Jumping into next one...\n')
                        continue

                    # Obtain chains from PDB file
                    template = Protein(f)
                    template.get_pdb_chains() # this method produces the chain_directory mentioned below
                    pdb_id = template.structure.get_id()

                    # Retrieving folder name for this template
                    chain_directory = template.chains_folder_name
                    
                    if os.path.exists(chain_directory):
                        # now checking if current chain name is in the directory as a pdb chain file
                        template_chain = ''

                        # 2. select template chain that equals the homolog found by BLAST and assign it as Protein instance
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
                        except OSError as e:
                            sys.stderr.write(f"Error: {chain_directory} : {e.strerror}")

                    else:
                        sys.stderr.write(f"The directory {chain_directory} could not be created, current template will be dismissed.\n") if options.verbose else None
                        continue
                    
                    # Appending each data frame 
                    # (and avoid including templates without binding site labels)
                    if type(template_chain_dataset) != str:
                        dataframes_list.append(template_chain_dataset)

                # Concatenating data sets for each chain
                all_templates_dataset = RandomForestModel().concat_training_data(dataframes_list)
                all_templates_dataset.to_csv('./template_datasets/' + chain_name + '_dataset' + '.csv')

                sys.stderr.write(f"\nAll homologous templates data set for current chain:\n") if options.verbose else None
                sys.stderr.write(str(all_templates_dataset)) if options.verbose else None
                sys.stderr.write('\n') if options.verbose else None

                # Splitting data into train/test
                X_train, X_test, Y_train, Y_test = RandomForestModel().split_data(all_templates_dataset)
                sys.stderr.write(f"Data successfully split into train and test sets!\n") if options.verbose else None

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
                    sys.stderr.write("The following metrics were obtained from the test data:\n")
                    sys.stderr.write(f'{metrics_table}\n')
                else:
                    sys.stderr.write("WARNING: metrics could not be obtained.")
                    
            else:
                sys.stderr.write(f"Final homologs list could not be retrieved from {fasta_file}, stop the program.")
                exit(3)
        else:

            sys.stderr.write(f"No homologs could be retrieved from {fasta_file}, stop the program")
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
    binding_sites_predicted_chimera = []
    binding_sites_real_provisional_chimera = []
    binding_sites_real_chimera = []
    binding_sites_shared_chimera = []
    binding_sites_result = []
    
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
            chimera_cmd_file = open(f"chimera_{target_protein.protein_id}.cmd", "w")
        except FileExistsError:
            sys.stderr.write("File already exists.") # realment aixo no cal pero no se q posar

        sys.stderr.write(f"We have created the folder {query_chain_directory} to store the chains of the query protein {query_pdb_id}\n")

        # Naming directory where we will store predictions
        out_dir = './query_predictions/'

        # now checking if current chain name is in the directory as a pdb chain file
        for file_name in os.listdir(query_chain_directory):
            sys.stderr.write(f"\nCALCULATING CHAIN FEATURES\n\n") if options.verbose else None
            sys.stderr.write(f"PDB file for {file_name[:-4]} could be obtained.\n\n") if options.verbose else None
            query_chain = Protein(str(query_chain_directory + '/' + file_name))  

            try:
                # Computing residue and atom features
                query_chain_features = ProteinFeatures(query_chain, './app/data/atom_types.csv')
                query_chain_features.residue_properties()
                query_chain_features.atom_properties()
                query_chain_features.is_cysteine()

                # Computing interactions
                query_chain_interactions = Interactions(query_chain, './app/data/atom_types.csv')
                query_chain_interactions.calculate_interactions()

                # Compute layer features (atom and residue properties, and interactions)
                query_chain_layer = Layer(query_chain, './app/data/atom_types.csv')
                query_chain_layer.get_layer_properties()

            except Exception as e:
                sys.stderr.write(f"Error:, {e}, file={sys.stderr}")
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

                    if not os.path.exists(out_dir):
                        os.makedirs(out_dir)

                    # EXTRACTING REAL BINDING SITES IF QUERY PROTEIN HAS THEM
                    real_binding_site = ExtractBindingSites().extract_binding_sites(file_name[:-5], file_name[4])

                    if real_binding_site != None:

                        sys.stderr.write(f"Real binding sites for {target_protein.protein_id} could be retrieved from BioLip. Comparing real vs predicted...\n") if options.verbose else None
                        
                        # Creating column for real class labels
                        query_chain.dataframe["real_label"] = int(0)

                        for element in real_binding_site:

                            # binary classification of real binding site residues in real_label column
                            for index, row in query_chain.dataframe.iterrows():

                                if element == index:
                                    query_chain.dataframe.loc[index,'real_label'] = int(1)
                        
                            # Appending residues to list that will contain all real binding sites of whole protein
                            binding_sites_real_provisional_chimera.append(element)

                        # Storing dataframe with predictions and real labels
                        query_chain.dataframe.to_csv(out_dir + query_chain.structure.get_id() + '.csv', columns=['prediction', 'real_label'])

                        # Calculating metrics for query protein prediction for proteins with already known labels
                        query_chain_metrics_table = RandomForestModel().model_evaluation(query_chain.dataframe["prediction"], query_chain.dataframe["real_label"])
                            
                        if query_chain_metrics_table:
                            sys.stderr.write("Metrics (real vs predicted):\n")
                            sys.stderr.write(f'\n{query_chain_metrics_table}\n')
                        else:
                            sys.stderr.write("WARNING: metrics could not be obtained.")
                    
                    else:
                        # Storing dataframe with only predictions 
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
                    binding_sites_predicted_chimera.append(f'{resnum}.{chain}') 
                    # Generate the result to be printed
                    binding_sites_result.append(res)
                    

        selection_pred = ','.join(binding_sites_predicted_chimera)

        chimera_cmd_file.write('# Coloring predicted binding site residues\n')
        chimera_cmd_file.write(f'color green :{selection_pred}\n')

        if len(binding_sites_real_provisional_chimera) > 0:

            for residue1 in binding_sites_real_provisional_chimera:
                #for residue1 in list_ch:

                resnum = residue1[5:]
                chain = residue1[0]

                binding_sites_real_chimera.append(f'{resnum}.{chain}')

                if f'{resnum}.{chain}' in binding_sites_predicted_chimera:

                    binding_sites_shared_chimera.append(f'{resnum}.{chain}')
        
            selection_real = ','.join(binding_sites_real_chimera)
            selection_shared = ','.join(binding_sites_shared_chimera)

            chimera_cmd_file.write(f'# Coloring real binding site residues\n')
            chimera_cmd_file.write(f'color red :{selection_real}\n')

            chimera_cmd_file.write(f'# Coloring matching binding site residues\n')
            chimera_cmd_file.write(f'color blue :{selection_shared}\n')


        chimera_cmd_file.close()

        pred =  ','.join(binding_sites_result)

        if binding_sites_result:
            sys.stdout.write(f'\n\nPredicted residues:\n{pred} \n')
        else:
            sys.stdout.write(f'\n{target_protein.protein_id} does not present binding sites.\n')

        sys.stderr.write(f'\nPrediction tables can be found in {out_dir} \n') if options.verbose else None

if __name__ == '__main__':
    main()


 