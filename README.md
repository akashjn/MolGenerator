The molecular library generation with the synthesis-aware rule-based molecule generation algorithm (SRMGA) 

1. Run the script 'molecule_generator.py' in the terminal using the commond: 

Usage: python molecule_generator.py
  Library required to run this code are
  - rdkit (use: pip install rdkit-pypi)
  - tqdm (use: pip install tqdm)
  - numpy (use: pip install numpy)
  - matplotlib (use: pip install matplotlib)
  
After you hit the enter, provide the address of the input file in the terminal 
    - if you want to hard code the file address modify variable inp in the molecule_generator.py, for example, 
    inp = r'address of the input_file.csv'

2. The molecule_generator.py has a GenerateSmiles() class that will generate a library of molecules (canonical SMILES) for a given core molecule and chemical building blocks.

    It takes following arguments:             

    (1) content: address of the input file with information about the core molecule and chemical building blocks
        (Note: always end file with commas so that our code know where to stop)

    (2) complexity_max: maximum CS value allowed in the library

    (3) x0 (int): CS value below which all molecules with CS < x0 are allowed in the library

    (4) beta (int or float): β is the exponential factor that is analogous to the Boltzmann factor 1/kT in thermodynamics. 
    The higher is the “temperature” (the smaller is β), the weaker is the penalization of molecular complexity. 

    (5) nmax: number of molecules in the library

    The program will save new molecules in a csv file with name: {input_file}_mol_library.csv"

3. A sample input file (input_file.csv) is provided in this repo
   
   The first line 'core' is followed by the SMILES string of the core molecules with growth points  indicated by R1, R2 and so on. (User can change the SMILES string and the position of R1, R2). 

   The ',,,,' separtes different sections in the input file 


   The details of the casts are followed by the line which has keyword 'group'
   
   The details of the building blocks and their weights are followed by the line which has keyword 'types'
   
   The details of the position of growth points are followed by the line which has keyword 'position'
   
   Note: always end file and each section with the commas
