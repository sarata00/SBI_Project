import csv
import chimera
import os
import pandas as pd

# Path to CSV file
predictions_dir = './query_predictions'

for csv_file in os.listdir(predictions_dir):
    # Read in the CSV file as a pandas dataframe
    df = pd.read_csv(csv_file)

    # Loop through the dataframe and select the residues
    for i, residue in df.iterrows():
        # Extract the residue number
        residue = row.iloc[0]
        res_num = int(residue[5:])
        print(residue)
        # Loop through the remaining columns and select the residues
        label = row.iloc[1]
        print(label)
            if int(label) == 1:
                # Select the residue in Chimera
                sel_str = f':{res_num}.{chimera.residues.one_letter_code(j)}'
                chimera.selection(sel_str).select(True)
