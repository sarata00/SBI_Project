import os
import subprocess

from app.protein import Protein

# como input: un archivo fasta con diferentes cadenas
class BLAST:
    def __init__(self, db_path):
        self.db_path = db_path
        self.out_file = None
        
    def search_homologs(self, input_file):
        templates = 0
        evalue = 0.001

        while templates == 0:
            evalue *= 10

            if evalue > 50:
                break
            
            # Define the BLAST command
            self.out_file = input_file.replace('.fasta', '.out')
            blast_cmd = f"blastp -query {input_file} -db {self.db_path} -evalue {evalue} -outfmt 6 -out templates/{self.out_file}"
            
            # Run the BLAST command
            try:
                result = subprocess.run(blast_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error: BLAST command returned non-zero exit status ({e.returncode}):")
                print(e.stderr.decode())
                return None
        

            # Get homologs from BLAST result file
            if os.path.exists(f"templates/{self.out_file}"):
                homologs = BlastResult(f"templates/{self.out_file}").get_result()
                templates = len(homologs)
            else:
                print(f"Error: BLAST output file {self.out_file} not found")
                return None
          
        return homologs



class BlastResult:

    def __init__(self, file_name):
        self.file_name = file_name
        self.results = []


    def get_result(self):

        with open(self.file_name) as file:
            for line in file:
                self.results.append(line.split('\t')[1])

        return self.results




