# Evaluating Molecules: A Python Tool for Fragment-to-Lead Optimization


This Python script automates the evaluation of large sets of molecules 
to identify and select candidates with improved binding properties
e.g. during a fragment-to-lead optimization process.
It was designed for a drug discovery process, where an initially identified small molecule 
has been varied and these variants have to be evaluated and ranked.
The script in its current form can handle any set of molecules 
like  variants of an original molecule, structurally diverse analogs, or new molecules 
generated from fragments by merging, growing, or linking strategies 
providing versatility across diverse workflows.
For the initially identified small molecule either a crystal structure or a docking prediction
is required. 

## Description

### Features

The script
- Converts 2D molecules to 3D conformers.
- Creates a receptor for docking based on an initial protein-molecule complex.  
- Performs docking calculations
- Re-scores the docked molecules using ScorePose and Smina.
- Predicts physicochemical properties xlogP, tPSA and ligand efficiencies.
- Creates a summary of the obtained results as a `.csv` file.
- Generates dotplots (`.jpg`) of scoring results for ranking of molecules.


### Workflow

This script follows a structured workflow:
1. **Input Molecule preparation**    
   All molecules of the collection (in 2D or 3D) are filtered, converted to 3D if necessary and relevant conformations are generated
   using OpenEye FILTER and OMEGA.
2. **Receptor Preparation**    
   A receptor model for docking is prepared based on the initial protein-molecule complex (crystal structure or docked structure)
   using OpenEye Spruce and Receptorindu.
3. **Docking**    
   All conformations of all molecules are docked into the receptor
   using OpenEye Posit. 
4. **File Conversion**             
   A PyMOL (`.pml`) script converts the resulting multimolecule.sdf file correctly into individual .pdb files.
5. **Rescoring**            
   All docked molecules (and the initial molecule) are rescored by OpenEye´s ScorePose and Autodock Smina.
6. **Property Calculation**     
   Additional physicochemical properties like molecular weight, XlogP, polar surface area (PSA) and ligand efficiencies are calculated.
7. **Result Summary**      
   The results are summarized in an spreadsheet in  `.csv` format.
8. **Plot Genaration**     
   A dotplot is genereted comparing the binding affinities of the two rescoring applications
   and with the initially identified molecule.
   
<br>

## Getting Started

### Dependencies

The following dependencies are required to run the script:

#### Python Packages
Packages not included in the standard Python library are labeled for installation. See also `requirements.txt`.
- **Python**: Version v3.11.5
   - `datetime` 
   - `subprocess`
   - `shutil`
   - `os`
   - `sys`
   - `pickle`
   - `pandas` as `pd` (Install via `pip install pandas`) == v2.1.4
   - `csv`
   - `matplotlib.pyplot` as `plt` (Install via `pip install matplotlib`) == v3.8.2
   - `importlib.util`
   - `time`
   - `glob`

#### Additional Tools

- **OpenBabel**: Version v2.4.1                         
   Used to convert `.pdb` files to `.pdbqt` files for Rescoring with Smina.
- **OpenEye Toolkit**: Version v2023.1.1               
   A valid OpenEye licence is required to use the toolkits.               
   Visit OpenEye, Cadence Molecular Sciences, Santa Fe, NM. [OpenEye](http://www.eyesopen.com) for more information.
- **PyMOL**: Version v2.5.0a0 or newer.                      
   Required for splitting and converting multifile `.sdf` into individual `.pdb` files.                      
- **Smina**: Version vOct 15 2019. Based on AutoDock Vina 1.1.2.                   
   Second rescoring application.

> Ensure that all required applications are installed and added to your PATH.

#### Installation Instructions

1. Clone the repository
2. Install the required Python packages using pip:                  
   `pip install -r requirements.txt`
3. Install Open Babel:                 
   Visit the *[Open Babel](http://openbabel.org/docs/Installation/install.html)* website for installation instructions.
4. Install PyMOL:            
   Visit *[PyMOL](http://www.pymol.org/pymol)*'s official page for installation:
   Schrödinger, L., & DeLano, W. (2020). PyMOL.
5. Install Smina:                  
   Please visit the original project homepage *[Smina](https://sourceforge.net/projects/smina/)* for installation.
6. Set up a conda environment and install OpenEye Toolkits:
   - Follow instructions for *[Conda Installation](https://educe-ubc.github.io/conda.html)* to install conda for OpenEye Toolkit installation
   - Follow OpenEye's virtual environment setup guide *[OpenEye_Toolkit_Installation](https://docs.eyesopen.com/toolkits/python/quickstart-python/linuxosx_anaconda.html)*
      - `conda create -n oepython -c openeye openeye-toolkits` creates a virtual environment named *oepython*
      - `source activate oepython` activates the newly created *oepython* environment
      - `export OE_LICENCE=/directory_to_licence/oe_licence.txt` localizes the OpenEye Licence
      - `conda install -c openeye openeye-toolkits` installs the OpenEye Python Toolkits into the new environment                  
      <br>
   > If the path to the OpenEye scripts differs from ./.conda/envs/oepython/bin/,
     make sure to change the definition of the path in the main script Script_fragment_evaluation.py
     by updating the definition of the variable openeyescript_paths [line 171, after =].
                     
                      
#### General informations and preparations

1. Directory structure

   Create a working directory named **Optimization_calculation**. This directory should contain three sub-directories:
      - ***Input_data*** contains also three subfolders with all important input data files:
         - *all_ligands*: contains the initial 3D molecule file (`.pdb`) of the ligand, for which optimized molecules were searched.             
           Name this ligand file by *`*_ligand.pdb`*.
         - *all_proteins*: stores the initial protein structure file (`.pdb`) in which the ligand was bound or docked, without a bound molecule.    
           Name this protein file by *`*_protein.pdb`*.
         - *all_proteins_with_ligand*: includes the initial protein file with the bound or docked molecule in `.pdb`.    
           Name this complex file by *`*.pdb`*.
         
         <br>
                       
         > Each folder may contain several structures.The structure you want to search for is defined by the name.
         >
         > Replace the placeholder `*` with your own name in all three folders (use the same name in each folder).
           When the script is started, this name is asked for as the `name of the initial molecule`.
         <br>
      - ***Results*** includes one subfolder *previous_runs*.       
        During the run all calculations are performed in this folder in an automatically created folder named by the date and the *`initial_molecule`*.          
        After finishing the calculations this folder is moved to the *previous_runs* folder.   
 
      - ***Scripts*** includes both Python scripts necessary for the workflow.
         - `Script_fragment_evaluation.py` (Main script)
         - `Pymol_script.pml`
               

2. Input and Output
   - **Input**: Molecular data files in `.sdf` format. Protein with bound or docked ligand, the protein without the ligand and the ligand alone.
   - **Output**: A summary report (`.csv`) detailing properties such as binding affinity, molecular weight, docking scores,  and ligand efficiencies.            
                 A dotplot (`.jpg`) comparing calculated affinities of two different Rescoring applications.
                 
#### Checks at the beginning of the script

The following checks are performed at the start of the script run to ensure that all requirements are met before the actual process begins.

1. **Python modules and external tools**   
   Verifies that all tools are installed and available in the path.
2. **OpenEye license**   
   Checks that the required Openeye licence file is present and active.
3. **Folders and files**  
   Ensures that all required input files and folders are available.
   
> The script is aborted with an error message if something is not found.
  

### Executing the program

1. **Prepare Input Molecules**              
   Generate molecules by individual optimization strategy and save them as `Optimization_*.sdf`
   in the **Optimization_calculation** folder.
2. **Activate Environment**              
   Activate the OpenEye environment:
   `conda activate oepython`
3. **Run the Script**                  
   `python Optimization_calculation/Scripts/Script_fragment_evaluation.py`

<br>

## Authors

**Anna-Maria Herzog**                  
Email: <anna-maria.herzog@uni-hohenheim.de>                       
LinkedIn: https://www.linkedin.com/in/anna-maria-herzog-35411b240/  

Institute of Biology, Department of Cellular Microbiology-190i
University of Hohenheim
70599 Stuttgart - Germany

<br>

## License

This project is licensed under the MIT License - see the LICENSE.md file for details

<br>

## References

1. **Python** v3.11.5 (https://www.python.org/)
   - pandas v2.1.4 (McKinney, W., & others. (2010). Data structures for statistical computing in python. In Proceedings of the 9th Python in Science Conference (Vol. 445, pp. 51–56))
   - matplotlib v3.8.2 (Hunter, J. D. (2007). Matplotlib: A 2D graphics environment. Computing in Science &amp; Engineering, 9(3), 90–95)                
2. **OpenBabel** v2.4.1 (http://openbabel.org/index.html)            
3. **PyMol** v2.5.0a0 (https://www.pymol.org/)
4. **Smina** vOct 15 2019. Based on AutoDock Vina 1.1.2 (http://pubs.acs.org/doi/abs/10.1021/ci300604z)             
5. **OpenEye Applications/Toolkits** v2023.1.1. Cadence Molecular Sciences, Santa Fe, NM. (https://docs.eyesopen.com)
   - OMEGA 4.2.1.1 (FILTER and OMEGA)  
     Hawkins, P.C.D.; Skillman, A.G.; Warren, G.L.; Ellingson, B.A.; Stahl, M.T.
     Conformer Generation with OMEGA: Algorithm and Validation Using High Quality Structures from the Protein Databank and the Cambridge Structural Database
     J. Chem. Inf. Model. 2010, 50, 572-584.    
   - SPRUCE 1.5.1.1
   - OEDOCKING 4.2.0.1 (ReceptorInDU, SorePose)
     Kelley, B.P.; Brown, S.P.; Warren, G.L.; Muchmore, S.W.
     POSIT: Flexible Shape-Guided Docking For Pose Prediction.
     J. Chem. Inf. Model., 2015, 55, 1771-1780. DOI: 10.1021/acs.jcim.5b00142
     McGann, M. FRED Pose Prediction and Virtual Screening Accuracy.
     J. Chem. Inf. Model., 2011, 51, 578-596. DOI: 10.1021/ci100436p
     McGann, M. FRED and HYBRID docking performance on standardized datasets.
     J. Comput. Aided Mol. Des., 2012, 26, 897-906. DOI: 10.1007/s10822-012-9584-8
   - POSIT 4.2.0.1
   - OEChem Toolkit 4.1.1.0 (convert.py, molpropcsv.py)

<br>

## Acknowledgements

- We gratefully acknowledge the free usage of OpenEye (https://docs.eyesopen.com ) software for academics, and funding from the Else Kröner-Fresenius Stiftung (to Prof. Dr. Julia Steuber and PD Dr. Günter Fritz).



