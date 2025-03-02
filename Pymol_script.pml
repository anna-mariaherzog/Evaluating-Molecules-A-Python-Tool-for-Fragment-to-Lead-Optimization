python

# help modules
import os
import subprocess
import datetime
import sys
import pickle


date = datetime.datetime.now().strftime("%Y%m%d_")
target_file_path = "./Optimization_calculation/Results/"

# Import the Filename
transfer_file = target_file_path+"transfer"
file = open(transfer_file,"rb")
#read binary
initial_molecule = str(pickle.load(file))
target_folder = target_file_path+date+initial_molecule
input_folder = target_folder+"/posit/"
input_file = "posit_"+initial_molecule+ "_optimization_docked"


for dateiname in os.listdir(input_folder):  
	cmd.load(input_folder+input_file+".sdf")
	cmd.split_states(input_file)
	cmd.delete(input_file)
	cmd.multifilesave(target_folder+"/pdb/{name}.pdb")
	cmd.delete("all")
	cmd.quit()

# # file is opened in PyMol
	# individual states are split into individual PyMol objects, then the original multistates file is deleted
	# finally each object is saved individually with specific name {} as .pdb in determined target pdb folder
	# all files are deleted from Pymol and pymol is closed

print("pdb datas are available now")
python end

