import csv
import chimera
import os

# Path to CSV file
predictions_dir = './query_predictions/'

for csv_file in os.listdir(predictions_dir):
    # Read in the CSV file as a csv object
    with open(str(predictions_dir + csv_file), 'r') as table:
        reader = csv.reader(table)
        # Skip the header row
        next(reader)
        # Loop through the rows and select the residues
        for row in reader:
            # Extract the residue number
            res = row[0]
            resnum = res[5:]
            res_name = res[2:5]
            chain = res[0]

            res_code = chimera.residues.one_letter_code(chimera.residues.three_letter_code(res_name))

            # Loop through the remaining columns and select the residues
            label = row[1]
            
            if int(label) == 1:
                # Select the residue in Chimera
                sel_str = ':' + resnum + '.' + chain
                chimera.selection(sel_str).select(True)
