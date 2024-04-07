import os


position = [1,2,3,4,5,6]
WT_PDB = "data/pten.pdb"
print(WT_PDB)
original_directory = os.getcwd()
folder_path = "tools"
os.chdir(folder_path)
from Feature_Cal.DSSP_RASA import DSSP_RASA
RASA = DSSP_RASA(position,WT_PDB)
os.chdir(original_directory)
print(original_directory)
print(RASA)

