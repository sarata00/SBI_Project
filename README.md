SBI Project 2022-23

Bind-Pred program

#### Authors
Alexandre Casadesús

Laia Fortuny

Sara Tolosa Alarcón 


# Introduction

Bind-Pred is a machine learning predictor of protein ligand binding site. It takes as input a protein structure ini PDB format and produce as output a list of aminoacids involved in the binding site and a `.cmd`file that allows us to visualize the predicted binding site with a molecular graphic software such as Chimera.

This program is tested in Linux and macOS using the `bash`command line to execute the program.

# Requeriments

It is reconmendable to create a new enviroment in conda to install all the dependencies and execute the program. This program has been developed using the Python 3.10.8 version. 


* [**Biopython**](https://anaconda.org/conda-forge/biopython), The Biopython Project is an international association of developers of freely available Python tools for computational molecular biology.

* [**UCSF Chimera**](https://pychimera.readthedocs.io/en/latest/install.html), UCSF Chimera is an extensible molecular visualization tool with a vast collection of modules gathered after years of development.

* [**BLAST from NCBI**](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html), the Basic Local Alignment Search Tool is used to find regions of similarity between biological sequences. The program compares nucleotide or protein sequences to sequence databases and calculates the statistical significance. 
  
* [**Scikit-learn**](https://scikit-learn.org/stable/install.html), a popular Python library for machine learning tasks such as classification, regression, and clustering. It provides a simple and efficient toolset for building predictive models and analyzing data.
  
* **Pandas**(), is an open-source data analysis library for the Python programming language. It provides data structures for efficiently storing and manipulating large datasets, as well as tools for analyzing and visualizing the data.
  
* **bs4**, another library of Python used for parsing HTML and XML documents, and extracting information from them. It allows us to search, navigate, and modify the parse tree using Python code, making it easier to extract specific information from web pages.
  
* **requests**, a Python library used for making HTTP requests, such as GET and POST requests, to web servers. It provides a simple API for sending HTTP requests and handling responses, making it easier to interact with web services and APIs.

* **tabulate**, a Pyhton library used for creating pretty, formatted tables from data in lists or other data structures


# Usage

In order to use the programme, you must first download all the modules and data of this GitHub repository. Then, in the same directory where the programme is downloaded, we can run it in the following way:

~~~bash
python3 script.py -p file.pdb                   # Run the program using a PDB file
python3 script.py -p path/to/pdb_files          # Run the program in a specific directory


python3 script.py -p file.pdb -v                # Use the argument *verbose* to print the standard error

~~~

Once the programme has been run, it will take some time to make predictions. As a consequence, different files and folders will be created:

* **FASTA files**: this files will be generated for each chain of the protein. If it consists of a single polypeptide chain, only one file will be created.
* **Chain folder**: a folder that contains separately the structure in PDB format each of the chains that make up the target protein.
* **Templates folder**: here we can find the BLAST results for each of the chains that are part of the target protein and the structures of these homologous proteins that have been found (in `.ent` format).
* **Template_dataset folder**:
* **Query_predictions**: for each chain of the protein a `.csv`file will be created with two columns: 'residue_name' and 'prediction'. The first one (residue_name) contains the residue information and the second one (prediction) the results of the prediction for each residue, represented as 1 if it is part of the binding site or 0 if not.
* **`chimera.cmd`file**: contains the necessary commands to select and colour the amino acids that are part of the binding site using Chimera.


## Example: 3bj1.pdb

Here is an example of how this program works using the structure of the hemoglobin protein (3bj1.pdb):

#### 1. Run the program

~~~bash
python3 script.py -p 3bj1.pdb -v
~~~

#### 2. Ouputs
Aquí deberíamos poner lo que aparece al poner verbose para que sea más intuitivo? o directamente las predicciones?

* `.csv`
* `chimera.cmd`  

#### 3. Visualize on Chimera
Once we have all the results, we can visualise the predicted amino acids in a visualiser such as Chimera. To do this we have to:
1. **Start chimera**
2. **Open the 3bj1.pdb file**
   
   -> *File...Open from the menu and file browse to locate the file*

   -> Command line: `open ./3bj1.pdb` (in case the PDB file is in the current directory)

3. **Open the chimera.cmd file**
   
   -> *File...Open from the menu and file browse to locate the file*

   -> Command line: `open ./chimera.cmd` (in case the script is in the current directory)


As a result, the predicted amino acids will appear green.

<p align="center">
    <img src="https://github.com/sarata00/SBI_Project/blob/main/3bj1_example.png?raw=true" width="600" alt="Bind-Pred image">
</p>