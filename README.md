# Bind-Pred program
*SBI Project 2022-23*

### **Authors**
Alexandre Casadesús Llimós

Laia Fortuny Galindo

Sara Tolosa Alarcón 


## **Introduction**

Bind-Pred is a machine learning predictor of protein ligand binding site. It takes as input a protein structure ini PDB format and produce as output a list of aminoacids involved in the binding site and a `.cmd`file that allows us to visualize the predicted binding site with a molecular graphic software such as Chimera.

This program is tested in Linux and macOS using the `bash`command line to execute the program.

## **Requeriments**

It is reconmendable to create a new enviroment in conda to install all the dependencies and execute the program. This program has been developed using the Python 3.10.8 version. 


* [**Biopython**](https://anaconda.org/conda-forge/biopython), The Biopython Project is an international association of developers of freely available Python tools for computational molecular biology.

* [**UCSF Chimera**](https://pychimera.readthedocs.io/en/latest/install.html), UCSF Chimera is an extensible molecular visualization tool with a vast collection of modules gathered after years of development.

* [**BLAST from NCBI**](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html), the Basic Local Alignment Search Tool is used to find regions of similarity between biological sequences. The program compares nucleotide or protein sequences to sequence databases and calculates the statistical significance. 
  
* [**Scikit-learn**](https://scikit-learn.org/stable/install.html), a popular Python library for machine learning tasks such as classification, regression, and clustering. It provides a simple and efficient toolset for building predictive models and analyzing data.
  
* [**Pandas**](https://pandas.pydata.org/), is an open-source data analysis library for the Python programming language. It provides data structures for efficiently storing and manipulating large datasets, as well as tools for analyzing and visualizing the data.
  
* [**bs4**](https://pypi.org/project/bs4/), another library of Python used for parsing HTML and XML documents, and extracting information from them. It allows us to search, navigate, and modify the parse tree using Python code, making it easier to extract specific information from web pages.
  
* [**requests**](https://pypi.org/project/requests/), a Python library used for making HTTP requests, such as GET and POST requests, to web servers. It provides a simple API for sending HTTP requests and handling responses, making it easier to interact with web services and APIs.

* [**tabulate**](https://pypi.org/project/tabulate/), a Pyhton library used for creating pretty, formatted tables from data in lists or other data structures


All this dependencies can be installed automatically by running the setup.py script as follows: 

~~~bash
python3 setup.py install
~~~


## **Usage**

Once we have all the requirements, we can start running our program.

### **1. Run the program**
In order to use the programme, you must first download all the modules and data of this GitHub repository. Then, in the same directory where the programme is downloaded, we can run it in the following way:

~~~bash
python3 script.py -p file.pdb                   # Run the program using a PDB file

python3 script.py -p file.pdb -v                # Use the argument *verbose* to print the standard error

~~~

### **2. Outputs**
Once the programme has been run, it will take some time to make predictions. As a consequence, different files and folders will be created:

* **FASTA files**: this files will be generated for each chain of the protein. If it consists of a single polypeptide chain, only one file will be created.
  
* **templates folder**: here we can find the BLAST results for each of the chains that are part of the target protein and the structures of these homologous proteins that have been found (in `.ent` format).
  
* **template_datasets folder**: it stores all for each homolog its dataframe with all the residues of the protein as records (rows) and 27 features as columns. Where the first 14 features are about the residue being considered, and the next 13 features are of the first shell of residues considered. The last column is the label that indicates if it is a residue that forms part of the binding site (1) or not (0).

* **query_predictions**: for each chain of the protein a `.csv` file will be created with two columns: 'residue_name' and 'prediction'. The first one (residue_name) contains the residue information and the second one (prediction) the results of the prediction for each residue, represented as 1 if it is part of the binding site or 0 if not.
  
* **target_chains**: it stores all pdb files for each protein chain.
  
* **`chimera_pdbID.cmd` file**: contains the necessary commands to select and colour the amino acids that are part of the binding site using Chimera.

* **Metrics (Accuracy, Precision, Recall, F1)**: in the command line for each chain will appear metrics to assess the algorithm predictions using the test data.
  

### **3. Visualize on Chimera**

Once we have all the results, we can visualise the predicted amino acids in a visualiser such as Chimera. To do this we have to:

1. **Start chimera**
   
2. **Open the target PDB file**
   
   -> *File...Open from the menu and file browse to locate the file*

   -> Command line: `open ./target.pdb` (in case the PDB file is in the current directory)

3. **Open the chimera.cmd file**
   
   -> *File...Open from the menu and file browse to locate the file*

   -> Command line: `open ./chimera.cmd`

## **Example: 3bj1.pdb**

Here is an example of how this program works using the structure of the hemoglobin protein with PDB ID 3bj1:

### 1. Run the program

~~~bash
python3 script.py -p 3bj1.pdb -v  # We use verbose option to follow all the process
~~~

### 2. Outputs

As a result, we get two main files:

* `3bj1_A.csv`, `3bj1_B.csv`, `3bj1_C.csv` and `3bj1_D.csv` files on **query_predictions folder**. These files store the predictions for each residue, represented as 1 if it is part of the binding site or 0 if not, as we have mentioned before.
   
* `chimera_3bj1.cmd`  script, which will be used to run on Chimera and visualize the predicted residues that could belong to the binding site.


In addition, in the command line we have some extra information such as some metrics that assess the predictions on the test set. In the current example, the results are:

|Protein_chain|Accuracy|Precision|Recall|F1|
|:----|:-----|:-----|:------|:-----|
|3bj1_A|0.99|1.00|0.97|0.98|
|3bj1_B|0.98|0.96|0.93|0.95|
|3bj1_C|0.99|0.97|1.00|0.98|
|3bj1_D|0.97|0.93|0.93|0.93|
|

Finally, the possible amino acids that could form part of the binding site are shown as follows:

~~~ bash
    Predicted residues:
    D_THR38,D_TYR41,D_PHE42,D_ASN44,D_PHE45,D_HIS63,D_THR66,D_ILE67,D_LEU71,D_LEU88,D_HIS92,D_LEU96,D_VAL98,D_ASN102,D_PHE103,D_LEU106,D_LEU141,B_TYR41,B_ASN44,B_PHE45,B_HIS63,B_THR66,B_ILE67,B_LEU71,B_LEU88,B_LEU91,B_HIS92,B_LEU96,B_VAL98,B_ASN102,B_LEU106,B_LEU141,A_MET32,A_THR39,A_TYR42,A_PHE43,A_HIS45,A_HIS59,A_THR62,A_ILE63,A_LEU67,A_LEU84,A_LEU87,A_HIS88,A_LEU92,A_VAL94,A_ASN98,A_PHE99,A_LEU102,A_LEU137,C_MET32,C_THR39,C_TYR42,C_PHE43,C_HIS45,C_HIS59,C_THR62,C_ILE63,C_LEU67,C_LEU84,C_LEU87,C_HIS88,C_LEU92,C_VAL94,C_ASN98,C_PHE99,C_LEU102,C_LEU137
~~~

Taking the first amino acid as an example, the format of the result means:

D_THR38 -> Amino acid Threonine (THR) at position 38 of the D chain.

### 3. Visualize on Chimera
Once we have our results, we can visualise the predicted amino acids using Chimera. To do this we have to:
1. **Start chimera**
2. **Open the 3bj1.pdb file**
   
   -> *File...Open from the menu and file browse to locate the file*

   -> Command line: `open ./3bj1.pdb` (in case the PDB file is in the current directory)

3. **Open the chimera.cmd file**
   
   -> *File...Open from the menu and file browse to locate the file*

   -> Command line: `open ./chimera_3bj1.cmd`


In the following image we can see the result obtained after running the previous commands:

<p align="center">
    <img src="https://github.com/sarata00/SBI_Project/blob/main/3bj1_example1.png?raw=true" width="600" alt="Bind-Pred image">
</p>

The different colours refer to:

* Green: amino acids predicted by Bind-Pred.
* Red: real amino acids that are part of the binding site according to information obtained from the BioLip database.
* Blue: predicted amino acids by Bind-Pred that match the real ones.

*In this case, we change manually the color of the ligand to pink*

## **Overall performance**

To test our Bind-Pred programme, we used several proteins and predicted the possible amino acids that are part of the binding site. The results are shown in the table below:












