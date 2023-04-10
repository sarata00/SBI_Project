from distutils.core import setup
from setuptools import find_packages

setup(name="Bind-Pred",
      version="1.0",
      description="Description of my project",
      long_description=open('README.md').read(),
      author="Alexandre Casadesús, Laia Fortuny and Sara Tolosa Alarcón",
      author_email="saratolosaalarcon@gmail.com",
      package_dir={"": "app"},
      packages=find_packages(where = "app"),
      install_requires =[
        "biopython >= 1.81",
        "pandas >= 2.0", 
        "scikit-learn", 
        "bs4", 
        "numpy", 
        "requests", 
        "tabulate"],
      script_name= "script.py")