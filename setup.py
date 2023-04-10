from distutils.core import setup, find_packages

setup(name="ProjectName",
      version="1.0",
      description="Description of my project",
      long_description=open('README.md').read()
      author="Alexandre Casadesús, Laia Fortuny and Sara Tolosa Alarcón",
      author_email="saratolosaalarcon@gmail.com",
      packages=find_packages(where = "app"),
      install_requires =[
        "bipython >= 1.81; python_version == 3.10.8",
        "pandas >= 2.0", 
        "scikit-learn", 
        "bs4", 
        "numpy", 
        "requests", 
        "tabulate"],
      py_modules=["Bio", "RandomForestClassifier"],
      script_name= "script.py")