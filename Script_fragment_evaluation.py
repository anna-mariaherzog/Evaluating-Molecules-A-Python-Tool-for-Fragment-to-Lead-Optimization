###########################################
### General Information and Preparation ###
###########################################

# Working folder is named "Optimization_calculation".This contains three subfolders:
# - "Input_data" with "all_proteins", "all_proteins_with_ligand" and "all_ligands"
# - "Results", where a special folder with current date is created, in which all calculations are performed and saved
#    Results contains an additional folder "previous_runs", where the actual folder will be moved to after the calculations are done
# - "Scripts" contains all necessary scripts for running the calculations

# Run Optimization and export data as "optimization_*.sdf" file in the "Results" folder
# Activate environment for running oepython (oeb.gz conversion and molecule properties calculation) with: 'conda activate oepython'
# Run script with: 'python Optimization_calculation/Scripts/Script_fragment_evaluation.py' 


# help modules
import datetime
import subprocess
import shutil
import os
import sys
import pickle
import pandas as pd
import csv
import matplotlib.pyplot as plt
import importlib.util
import time
import glob

RED = "\033[91m"
GREEN = "\033[92m"
YELLOW = "\033[93m"
RESET = "\033[0m"

### Interactive query of variable parameters (number of cores and name of the initial molecule)

print(f"{YELLOW} Please specify the following parameters:{RESET}")
# Request number of cores (default 4)
while True:
    cores = input ("    Number of cores (default: 4. To confirm, press enter):"). strip()
    if cores == "":
        cores = 4 # default value
    try:
        cores = int(cores)
        if cores > 0:
            break
        else:
            print(f"{RED}      Number of cores must be at least 1!{RESET}")
    except ValueError:
        print(f"{RED}      Please enter a valid number!{RESET}")

# Definition of the initial molecule name
while True:
    molecule_file = input("    Name of the initial molecule:").strip()
    if molecule_file:
        break 
    print(f"{RED}      Please enter the initial name of your molecule!{RESET}")

# Assign values to existing variables
core_number = cores
initial_molecule = molecule_file

# Summary of inputs
print(f"\n{GREEN} The script will start with the following parameters:{RESET}")
print(f"   - Number of cores: {core_number}")
print(f"   - Molecule file: {initial_molecule}\n")


### Check for required tools, openeye license, folders and files

# List of used python modules
required_modules = ["datetime",
                    "subprocess",
                    "shutil",
                    "os",
                    "sys",
                    "pickle",
                    "pandas",
                    "csv",
                    "matplotlib",
                    "importlib",
                    "time",
                    "glob"
]  

# List of required external programs
required_tools = ["obabel",
                  "filter",
                  "oeomega",
                  "spruce",
                  "receptorindu",
                  "posit",
                  "scorepose",
                  "molpropcsv.py",
                  "pymol",
                  "smina"
] 

# Checking the python modules
missing_modules = [mod for mod in required_modules if importlib.util.find_spec(mod) is None]

# checking the external programs
missing_tools = [tool for tool in required_tools if shutil.which(tool) is None]

# Check results
if missing_modules:
    print(f"{RED} Missing Python modules: {RESET}", missing_modules)
    sys.exit(1)
else:
    print(f"{GREEN} All required Python modules are installed.{RESET}")

if missing_tools:
    print(f"{RED} Missing external programs:{RESET}", missing_tools)
    print(f"{YELLOW} Possible solutions{RESET}")
    print(f" -> Please make sure that the tools are installed and available in the path.")
    print(f" -> Make sure that the oepython environment is activated for molpropcsv.py. -> conda activate oepython{RESET}")
    sys.exit(1)
else:
    print(f"{GREEN} All programmes are installed and available.{RESET}")


# Checking OpenEye license file
default_license_path = os.path.expanduser("~/.oe_license.txt")  # Standard path for OpenEye licence
env_license_path = os.environ.get("OE_LICENSE", default_license_path)  # Alternatively from environment variable

if os.path.isfile(env_license_path):
    print(f"{GREEN} OpenEye licence file found: {env_license_path} {RESET}")
else:
    print(f"{RED} OpenEye licence file not found! Please check the path.{RESET}")
    sys.exit(1)

# Check Openeye license with Openeye module
def check_openeye_python_licence():
    try:
        oechem = importlib.import_module("openeye.oechem")
        if oechem.OEChemIsLicensed():
            print(f"{GREEN} OpenEye Python licence is active{RESET}")
        else:
            print(f"{RED} OpenEye Python licence not recognised{RESET}")
    except ModuleNotFoundError:
        print(f"{RED} OpenEye Python modules not found!{RESET}")
        print(f"{YELLOW} Possible solutions:{RESET}")
        print(" -> Make sure that the OpenEye toolkits are installed -> https://docs.eyesopen.com/toolkits/python/quickstart-python/linuxosx_anaconda.html)")
        print(" -> Make sure that the oepython environment is activated -> conda activate oepython")
        sys.exit(1)
    except Exception as e:
        print(f"{RED} Error when checking the OpenEye licence: {e}{RESET}")
        sys.exit(1)

check_openeye_python_licence()

# Define folders and files and check for their availability

date = datetime.datetime.now().strftime("%Y%m%d_")
target_file_path = "./Optimization_calculation/Results/"

subprocess.Popen("mkdir "+target_file_path+date+initial_molecule, shell=True)
time.sleep(0.1) # Short waiting time for the new folder to be recognized by the system

# subprocess.Popen() writes defined text in terminal commandline
# + links the variable with text, this text must be in quotation marks
# Folder in target_file_path is created named with actual date and initial_molecule name

target_folder = target_file_path+date+initial_molecule
protein_ligand_file = "./Optimization_calculation/Input_data/all_proteins_with_ligand/"+initial_molecule+".pdb"
protein_file = "./Optimization_calculation/Input_data/all_proteins/"+initial_molecule+"_protein.pdb"
ligand_file = "./Optimization_calculation/Input_data/all_ligands/"+initial_molecule+"_ligand.pdb"
optimization_file = "./Optimization_calculation/Results/optimization_*.sdf"

script_path = "./Optimization_calculation/Scripts/"
openeyescript_paths = "./.conda/envs/oepython/bin/"
plot_title = ""+initial_molecule+""

# Check foldersystem and files

required_directories = [
    "./Optimization_calculation/Input_data/all_proteins_with_ligand", # input_path_complex
    "./Optimization_calculation/Input_data/all_proteins", # input_path_protein
    "./Optimization_calculation/Input_data/all_ligands", # input_path_ligand
    "./Optimization_calculation/Scripts/", # script_path
    "./.conda/envs/oepython/bin/", # openeyescript_paths
     "./Optimization_calculation/Results", # target_file_path
    ""+target_file_path+date+initial_molecule+"", # target_folder
    "./Optimization_calculation/Results/previous_runs" # final results folder
]

required_files = [
    ""+protein_ligand_file+"",
    ""+protein_file+"",
    ""+ligand_file+"",
#    ""+optimization_file+"" matching file with * as placeholder is defined later
]

missing_directories = []
missing_files = []

for directory in required_directories:
    if not os.path.exists(directory):
        missing_directories.append(directory)
        print(f"{RED} Folder '{directory}' is missing.{RESET}")
        print(" -> Create with os.makedirs(directory)")
    else:
        print(f"{GREEN} Folder '{directory}' exists.{RESET}")

for file in required_files:
    if not os.path.exists(file):
        missing_files.append(file)
        print(f"{RED} Missing file: '{file}'{RESET}")
    else:
        print(f"{GREEN} File '{file}' exists. {RESET}")
        

# Skript abort if a folder or file misses
if missing_directories or missing_files:
    print(f"\n{RED} Error!{RESET}")
    
    if missing_directories:
        print(f"{RED} Missing folders:{RESET}")
        for directories in missing_directories:
            print(f"   - {directories}")
            
    if missing_files:
        print(f"{RED} Missing files: {RESET}")
        for file in missing_files:
            print(f"  - {file}")
    
    sys.exit(1)

matching_files = glob.glob(optimization_file)
if not matching_files:
    print(f"{RED} No file found that matches '{optimization_file}'{RESET}")
    sys.exit(1)
input_database_file = matching_files[0]
print(f"{GREEN} Found input file: {input_database_file}{RESET}")

print(f"\n{GREEN} All folders an files exist.{RESET}")

### Function to check if a file exists. Used several times during the run to check output generation.
def check_file_exists(filename, step_name=""):
    if not os.path.isfile(filename):
        print(f"{RED} Error: Required file '{filename}' after step '{step_name}' not found!{RESET}")
        sys.exit(1)
    else:
        print(f"{GREEN} File '{filename}' created successfully after '{step_name}'.{RESET}")

############################
### Database preparation ###
############################
# filter and conformation preparation of database optimization molecules

subprocess.Popen("mkdir -p "+target_folder+"/database", shell=True)
# mkdir -p creates folder with all parent directories (nested directories)
text  = "filter -in "+optimization_file+" -out \
        "+target_folder+"/database/optimization_filtered.oeb.gz -prefix filter"
subprocess.Popen(text,shell=True).wait()
subprocess.Popen("mv filter* "+target_folder+"/database",shell=True)

text = "oeomega pose -in "+target_folder+"/database/*_filtered.oeb.gz -out\
       "+target_folder+"/database/optimization_conformations.oeb.gz -prefix oeomega -mpi_np "+str(core_number)+""
subprocess.Popen(text,shell=True).wait()
subprocess.Popen("mv oeomega* "+target_folder+"/database",shell=True)

subprocess.Popen("mv "+target_file_path+"*.sdf "+target_folder+"/database",shell=True)

output_database_preparation = ""+target_folder+"/database/optimization_conformations.oeb.gz"
check_file_exists(output_database_preparation, "Database preparation")

# First runs OpenEye's FILTER to remove molecules that do not match the defined physicochemical properties (default physical property limits of OpenEye, preventing conformation generation of e.g. large polypetides and very flexible molecules)
# and then oeomega pose to generate optimal conformers for molecule alignment and pose prediction by docking.
# The generated output files and also the original input optimizationfile are moved to the folder for the generated database. 
 
###############
### Docking ###
###############

### Spruce
# Prepares experimental or modeled proteins for docking or simulation by adding missing protons/residues and optimizing hydrogen atoms
subprocess.Popen("mkdir -p "+target_folder+"/posit/Receptor", shell=True)
text = "spruce -in "+protein_ligand_file+ " -out \
        "+target_folder+"/posit/Receptor/"+date+"prepared_protein.oedu -prefix spruce"
subprocess.Popen(text,shell=True).wait()
subprocess.Popen("mv spruce* "+target_folder+"/posit/Receptor",shell=True)
 
print("Spruce is finished")


### Receptorindu
# creates a receptorfile for dockinginput containing localisation and shape of the active site
text = "receptorindu -in "+target_folder+"/posit/Receptor/"+date+"prepared_protein.oedu -out \
       "+target_folder+"/posit/Receptor/"+date+"receptor.oedu -prefix receptorindu"
subprocess.Popen(text,shell=True).wait()
subprocess.Popen("mv receptorindu* "+target_folder+"/posit/Receptor",shell=True)
  
print("Receptorindu is finished, OpenEye Receptor is prepared")

output_receptor_preparation = ""+target_folder+"/posit/Receptor/"+date+"receptor.oedu"
check_file_exists(output_receptor_preparation, "Docking-Receptor preparation")

### Posit
# docking of the multiconformer optimizationdatabase based on the inputmolecule to stay as close as possible to this initial molecule
# OpenEye's Posit is a pose-prediction tool choosing the best docking method based on molecule similarity and returns a probability that the docked pose is within 2 Â of the actual pose (minimum probability is set to 0 to output all poses)
text = "posit -receptor "+target_folder+"/posit/Receptor/"+date+"receptor.oedu -dbase "+target_folder+"/database/optimization_conformations.oeb.gz \
       -prefix posit -docked_molecule_file posit_"+initial_molecule+ "_optimization_docked.oeb.gz -mpi_np "+str(core_number)+""
subprocess.Popen(text,shell=True).wait()
subprocess.Popen("mv posit* "+target_folder+"/posit",shell=True)
# moves all posit outputs in the generated posit folder

# oeb.gz output from posit will be transformed to .sdf by openeyes convert.py

text = openeyescript_paths+"convert.py "+target_folder+"/posit/posit_"+initial_molecule+ "_optimization_docked.oeb.gz \
       "+target_folder+"/posit/posit_"+initial_molecule+ "_optimization_docked.sdf"
subprocess.Popen(text,shell=True).wait()
 
output_docking = ""+target_folder+"/posit/posit_"+initial_molecule+ "_optimization_docked.sdf"
check_file_exists(output_docking, "Docking")
print("Docking is finished, Outputfile is transformed to .sdf")


##############################
### Splitting states Pymol ###
##############################

text = "mkdir -p "+target_folder+"/pdb"
subprocess.Popen(text,shell=True)

transfer_file = target_file_path+"transfer"
file = open(transfer_file,"wb")
#write binary
pickle.dump(initial_molecule, file)
file.close()
text = "pymol "+script_path+"Pymol_script.pml"
subprocess.Popen(text,shell=True).wait()

os.remove(transfer_file)
# file with the initial_molecule name ("transfer") is created in the working folder to transfer the initial_molecule name from this script to the
# pymol script (Pymol_script.pml). After execution of the pymol script the transfer file is not required and deleted. 

shutil.copy(ligand_file, target_folder+"/pdb")
# shutil.copy(original, target)copies a file from one folder to another
# initial ligand from optimization search is also copied in the pdb folder, to also get rescored


#################
### ScorePose ###   
#################
 
text = "mkdir -p "+target_folder+"/scorepose"
subprocess.Popen(text,shell=True)
 
for file_name in os.listdir(target_folder+"/pdb"):
    if file_name.endswith(".pdb"):
        molecule_files = os.path.join(target_folder+"/pdb", file_name)
        command = "scorepose -receptor "+target_folder+"/posit/Receptor/"+date+"receptor.oedu \
                  -dbase "+molecule_files+" -prefix scorepose_"+file_name+" \
                  -no_extra_output_files true -score_file scorepose_"+os.path.splitext(file_name)[0]+"_score.txt -mpi_np "+str(core_number)+"" 
        os.system(command)
# rescores each .pdb file on its own. The output files are an oeb.gz file and a score.txt file with the score. os.path.splitext(file_name)[0] separates filename and pdb
# and keeps only the filename without the pdb ([0] > 1st list element) for naming the output files.
 
subprocess.Popen("mv scorepose* "+target_folder+"/scorepose",shell=True).wait()

# Function for writing out Affinities from each scorepose.txt file, by defining the second line not starting with T (first line), seperating the line by ' 'and taking the third part containing the affinity value
def readScore(path): 
    try:
        with open(path) as text:
            for line in text:
                if line[0]!="T":
                    a = str(line.split('	')[2])
                    return float(a) #float to get floating point number from Affinity text string
    except:
        return "Affinity Scorepose"
 
 
# write molecule names and ScorePose Affinities in one Scorepose.csv file
file_names = ["scorepose_Molecule_number_"] # needed for later split of scorepose_ at beginning an _conformation number_ at the end
file_names.extend([filename for filename in os.listdir(target_folder+"/scorepose") if os.path.isfile(os.path.join(target_folder+"/scorepose/", filename)) and filename[-1]=="t"]) # takes only the .txt files not the oeb.gz files (last letter "t")
scorepose_csv = target_folder+"/scorepose/scorepose_affinities.csv"
 
 
with open(scorepose_csv, mode='w', newline='') as file:
     writer = csv.writer(file)
     for filename in file_names:
         path = target_folder+"/scorepose/"+filename
         Affinity_Scorepose = readScore(path) # function defined above is called
         filename_short = filename.rsplit("_",2)[0].split("_", 1)[1] # removes scorepose at the begining and Conformationnumber_score at the end of the filename
         writer.writerow([filename_short,Affinity_Scorepose])

output_scorepose_affinities = ""+target_folder+"/scorepose/scorepose_affinities.csv"
check_file_exists(output_scorepose_affinities, "ScorePose")
print("ScorePose is finished")

#################
### Openbabel ###
#################

subprocess.Popen("mkdir "+target_folder+"/pdbqt",shell=True)
text = "obabel -i pdb "+target_folder+ "/pdb/*.pdb -o pdbqt -O \
       "+target_folder+ "/pdbqt/.pdbqt -m"
subprocess.Popen(text,shell=True).wait()
 
print("Openbabel is finished; pdbqt files for Smina are available")
# folder in destination file path is created and openbabel multiple inputfile command (-m) is executed


#############
### Smina ###
#############
# Fragments must be removed from the protein file before smina runs
# H atoms are added by default to the ligands
subprocess.Popen("mkdir "+target_folder+"/smina",shell=True)

for file_name in os.listdir(target_folder+ "/pdbqt"):
    if file_name.endswith(".pdbqt"):
        ligands_file = os.path.join(target_folder+ "/pdbqt", file_name)
        output_file = os.path.splitext(file_name)[0] + ".log"
        log_path = os.path.join(target_folder+ "/smina", output_file)
        command = "smina --score_only -r "+protein_file+" -l \
                  "+ligands_file+ " --log "+log_path+" -q --cpu 1"
        os.system(command)
# Search in pdb filepath for any file with .pdbqt extension
# The pdbqt filename with .log extension is used as outputfile
# q = suppress output message

# Function for writing out Affinities from each smina.log data, by defining the line starting with 'A', seperating the line by ' 'and taking the middle part containing the affinity value
# function will be executed later when it is called
def readSmina(path):
    try:
        with open(path) as text:
            for line in text:
                if line[0]=="A":
                    a = str(line.split(' ')[1])
                    return float(a) #float to get floating point number from Affinity text string
    except:
        return "Affinity Smina"

            
# write molecule names and ScorePose Affinities in one Smina.csv file
file_names = ["Molecule"]
file_names.extend([filename for filename in os.listdir(target_folder+"/smina") if os.path.isfile(os.path.join(target_folder+"/smina/", filename))]) 
smina_csv = target_folder+"/smina/smina_affinities.csv"


with open(smina_csv, mode='w', newline='') as file:
    writer = csv.writer(file)
    for filename in file_names:
        path = target_folder+"/smina/"+filename
        Affinity_Smina = readSmina(path) # function defined above is called
        filename_short = filename.rsplit("_",1)[0] # removes Conformationnumber.log at the end of the filename
        writer.writerow([filename_short,Affinity_Smina])
        
output_smina_affinities = target_folder+"/smina/smina_affinities.csv"
check_file_exists(output_smina_affinities, "Smina")        
print("Smina is finished")


#######################
### Results summary ###  
#######################
# Creates one .csv file containing all important informations and calculation results  
 
# defining paths of the outputfiles from all calculations to write out the necessary informations and of the final result outputfile
###path_oeomega = target_folder+"/database/oeomega_rpt.csv"
path_ligand_properties = target_folder+"/posit/"+initial_molecule+ "_ligand.csv"
path_posit = target_folder+"/posit/posit_"+initial_molecule+ "_optimization_docked.csv"
path_scorepose = target_folder+"/scorepose/scorepose_affinities.csv"
path_smina = target_folder+"/smina/smina_affinities.csv"
path_posit_folder = target_folder+"/posit/"
output_path = target_folder+"/"+initial_molecule+"_calculation_summary.csv"

# Apply openeyes "molpropcsv.py" first to the initial fragment, afterwards to the analogs to estimate molecule properties
# Then both lists are joined in a dataframe (dfMolPropFull)
text = openeyescript_paths+"molpropcsv.py "+ligand_file+ "\
       "+path_ligand_properties
subprocess.Popen(text,shell=True).wait()
dfLigandProp = pd.read_csv(path_ligand_properties, delimiter=",")
dfLigandProp['TITLE'] = initial_molecule
dfLigandProp.to_csv(path_ligand_properties, index=False)

text = openeyescript_paths+"molpropcsv.py "+target_folder+"/posit/posit_"+initial_molecule+ "_optimization_docked.oeb.gz \
       "+target_folder+"/posit/posit_"+initial_molecule+ "_optimization_docked.csv"
subprocess.Popen(text,shell=True).wait()

def shorten_Moleculename(text):
    return text.rsplit('_', 1)[0]
dfMolPropAnalogs = pd.read_csv(path_posit, delimiter=",")
dfMolPropAnalogs['TITLE'] = dfMolPropAnalogs['TITLE'].apply(shorten_Moleculename)
dfMolPropAnalogs.to_csv(path_posit, index=False)

files = os.listdir(path_posit_folder)
files_csv = [f for f in files if f[-3:] == 'csv']

df_list = []
for f in files_csv:
    data = pd.read_csv(os.path.join(path_posit_folder,f))
    df_list.append(data)
dfMolPropFull = pd.concat(df_list)
dfMolPropFull.to_csv(output_path, index=False)

# a new df(MolProp)is created containing only the columns with important informations.
# 2d PSA: 2-dimensional polar surface area (correlates with transport properties like intestinal absorption and blood-brain barrier)
# Posit method lists the choosen docking method for each molecule. Posit probability the probability that the docked pose is within 2 Â of the actual pose.  
dfMolProp = dfMolPropFull.loc[:,('TITLE','SMILES','atom count','molecular weight','XLogP','Solubility','2d PSA','POSIT::Probability','POSIT::Method')]

# 'Title' in dfMolProp corresponds to Molecule from Rescoringapplications + _conformation number.
# First Title is renamed to Molecule and the last part (conformation numer) is deleted so it now corresponds to 'Molecule' from rescorings for later sorting.
dfMolProp.rename(columns={'TITLE':'Molecule'}, inplace=True)
dfMolProp.rename(columns={'atom count':'Heavy Atoms'}, inplace=True)

# Merge the dfMolProp columns with the dfScorePose and dfSmina columns in dfOut. 
# The merging is done using the "Molecule" parameter of the "Smina" data frame (left),
# so that the original fragment is not removed from the list (not included in the MolProp File).
# merge smina
dfSmina = pd.read_csv(path_smina, delimiter=",")
dfOut = pd.merge(dfMolProp, dfSmina, on='Molecule', how="left")
# merge scorepose
dfScorepose = pd.read_csv(path_scorepose, delimiter=",")
dfOut = pd.merge(dfOut, dfScorepose, on='Molecule', how="left")
# adds three new columns for ligand efficiency, calculates them from the 'Smina affinities' and saves the whole table as .csv
dfOut['Ligand efficiency'] = -dfOut['Affinity Smina']/dfOut['Heavy Atoms']
dfOut['Size independent ligand efficiency'] = -dfOut['Affinity Smina']/(dfOut['Heavy Atoms']**0.3) # to the power of 0.3
dfOut['lipophilic efficiency dependent lipophilicity'] = dfOut['XLogP']/dfOut['Ligand efficiency']
dfOut.to_csv(output_path, index=False)

output_csv_summary = target_folder+"/"+initial_molecule+"_calculation_summary.csv"
check_file_exists(output_csv_summary, "Results summary - .csv")  

# plot all identified potentially optimized molecules in a dotplot (X=ScorePose_Affinity, Y=Smina_Affinity)
x_column = 'Affinity Scorepose'
y_column = 'Affinity Smina'
plt.scatter(dfOut[x_column], dfOut[y_column], marker='s', s=10, color='k') # black squares
plt.title(initial_molecule+plot_title)
plt.xlabel('ScorePose', fontsize=30)
plt.ylabel('Smina', fontsize=30)
plt.axis((-20, 10, -20, 50))
plt.subplots_adjust(left=0.15, bottom=0.15)
# Mark the initial fragment as red circle and name this point by initial_molecule
initial_fragment = dfOut[dfOut['Molecule']==initial_molecule]
plt.scatter(initial_fragment[x_column], initial_fragment[y_column], marker='o', color='r')
plt.annotate(initial_molecule, (initial_fragment[x_column].iloc[0], initial_fragment[y_column].iloc[0]), \
             textcoords='offset points',xytext=(5,-15))
plt.savefig(target_folder+"/ScorePose_Smina_Plot.jpg")

output_affinity_dotplot = target_folder+"/ScorePose_Smina_Plot.jpg"
check_file_exists(output_affinity_dotplot, "Results summary - dotplot")  

# Move working folder to "previous runs" and delete the old folder
shutil.move(target_folder, target_file_path+"previous_runs/")
# shutil module for high-level operations on files (also several files)
# shutil.move (x,y) recursively moves a file or directory x to another location y

print("Evaluation_Script_is_done")