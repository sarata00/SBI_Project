from distutils.core import setup
from setuptools import find_packages

setup(name="Bind-Pred",
      version="1.0",
      description="This program predicts the binding sites of a pdb query protein considering the physicochemical properties and interactions of its residues and their local environments using a machine learning approach.",
      long_description=open('README.md').read(),
      author="Alexandre Casadesús Llimós, Laia Fortuny Galindo and Sara Tolosa Alarcón",
      author_email="alexandre.casadesus01@estudiant.upf.edu, laia.fortuny01@estudiant.upf.edu, sara.tolosa02@estudiant.upf.edu",
      package_dir={"": "app"},
      packages=find_packages(where = "app"),
      install_requires =[
        "biopython >= 1.81",
        "pandas >= 2.0", 
        "scikit-learn >= 1.2.2", 
        "bs4 >= 0.0.1", 
        "numpy >= 1.23.4", 
        "requests >= 2.28.2", 
        "tabulate >= 0.9.0"],
      script_name= "BindPred.py",
      classifiers=["Programming Language :: Python :: 3.10.12"])
