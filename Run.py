from APMA import APMA
import os
import glob
import traceback
import warnings
def run_APMA(Position, Protein):
    warnings.filterwarnings("ignore")
    # pocess position file
    print("+++++++++++++++++++++++")
    print("++++++++ APMA +++++++++")
    print("+++++++++++++++++++++++")
    f = open(Position,"r")
    all = f.readlines()
    all_new = []
    for i in all:
        if i == '\n':
            pass
        else:
            string_split = i.split("\t")
            string_a = string_split[0]
            string_b = string_split[1]
            string_b = string_b[:1] + "A" + string_b[1:]

            FoldX_type_string = string_a + "\t" + string_b
            all_new.append(FoldX_type_string)
    f.close()

    f = open(Position,"w")
    for i in all_new:
        f.write(i)
    f.close()


    user_pdb = Protein
    print("Your Protein pdb file is",user_pdb)
    user_protein_name = user_pdb.split("/")[-1].rstrip(".pdb")
    print(f"Your protein name is {user_protein_name}")

    APMA(
        Protein_name = user_protein_name,
        file_path = "data/position.txt",
        WT_PDB = user_pdb
        )




import os

def delete_files_in_directory(directory):
    # 遍历目录中的所有文件和子目录
    for item in os.listdir(directory):
        item_path = os.path.join(directory, item)
        
        # 如果是文件，则删除
        if os.path.isfile(item_path):
            os.remove(item_path)
        
        # 如果是目录，则递归调用该函数
        elif os.path.isdir(item_path):
            delete_files_in_directory(item_path)
