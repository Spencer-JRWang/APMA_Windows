import os
import re
import pandas as pd
import subprocess
from sklearn.preprocessing import scale

def APMA(WT_PDB, Protein_name, file_path, MSA_data = "data/query_msa.fasta", FoldX = "tools/FoldX"):
    '''
    paramaters:
    Protein_name = input("Please provide your protein name")
    file_path = input("Please provide your description file")
    FoldX = input("Please provide your route to FoldX") default "/home/wangjingran/APMA/FoldX"
    WT_PDB = input("Please provide your route to Wild Type PDB")
    Mut_PDB = FoldX
    MSA_data = input("Please provide your MSA file") default home/wangjingran/APMA/data/query_msa.fasta
    '''
 
    from Feature_Cal.Blast_MSA import extract_sequence_from_pdb
    from Feature_Cal.Blast_MSA import blast_search
    from Feature_Cal.Blast_MSA import run_clustal
    phenotype_list = []
    site_list = []
    mutation_list = []
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            if len(re.findall(r'\d+', columns[1])) == 1:
                site_list.append(re.findall(r'\d+', columns[1])[0])
            phenotype_list.append(columns[0])
            mutation_list.append(columns[1][:1] + columns[1][2:])

    category = phenotype_list
    set_category = list(set(category))
    position = site_list
    position = [int(num) for num in position]
    pdb_sequences = extract_sequence_from_pdb(f'data/{Protein_name}.pdb')
    protein_sequence = ''.join(pdb_sequences)
    sequence = protein_sequence
##############################################################################################################################
    # Perform BLAST search and save results to a file
    output_file = "data/blast_results.fasta"
    blast_search(sequence, output_file)
    
    # 输入的FASTA文件，这里假设你已经有了一些同源序列的FASTA文件
    with open("data/blast_results.fasta", "r") as f:
        sequence_blast = []
        s_lines = f.readlines()
        for i in s_lines:
            if i.startswith(">") or i == "\n":
                pass
            else:
                sequence_blast.append(i)
        sequence_blast = list(set(sequence_blast))
        import random
        # 这里加上一个判断，如果少于了200个就把所有的都选上去
        if len(sequence_blast) > 200:
            random_numbers = random.sample(range(1, len(sequence_blast)), 200)
            sequence_blast = [sequence_blast[i] for i in random_numbers]
    
    with open("data/blast_results.fasta", "w") as f:
        f.write(">Input_Seq" + "\n")
        f.write(sequence + "\n")
        for i in range(len(sequence_blast)):
            f.write(">sequence" + str(i + 1) + "\n")
            f.write(sequence_blast[i])
##############################################################################################################################
    print("MSA started ...", end=" ")
    original_directory = os.getcwd()
    folder_path = "tools/clusteralo"
    os.chdir(folder_path)
    input_fasta = "../../data/blast_results.fasta"
    # 输出的FASTA文件，用于保存比对结果
    output_fasta = "../../data/query_msa.fasta"
    # 运行多序列比对
    run_clustal(input_fasta, output_fasta)
    with open("../../data/query_msa.fasta", 'r') as f:
        lines = f.readlines()
    lines[0] = '>Input_seq\n'
    
    with open("../../data/query_msa.fasta", 'w') as f:
        f.writelines(lines)
    os.chdir(original_directory)
    print("Success")
##############################################################################################################################
    # 计算蛋白质的相对可及表面积
    original_directory = os.getcwd()
    folder_path = "tools"
    os.chdir(folder_path)
    from Feature_Cal.DSSP_RASA import DSSP_RASA
    RASA = DSSP_RASA(position,WT_PDB)
    os.chdir(original_directory)
##############################################################################################################################
    # 计算熵和保守性
    from Feature_Cal.sequence import cal_entropy
    from Feature_Cal.sequence import cal_coevolution
    SI = cal_entropy(MSA_data,position)
    MI = cal_coevolution(MSA_data,position)
##############################################################################################################################
    # use rate4site to score each site of the protein
    import subprocess
    import time
    original_directory = os.getcwd()
    folder_path = "tools"
    print("rate4site started ...", end=" ")
    os.chdir(folder_path)
    process = subprocess.Popen("rate4site -s ../data/query_msa.fasta -o ../data/score.txt", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    #print("Output:", output.decode().strip())
    #print("Error:", error.decode().strip())
    os.chdir(original_directory)
    #print(original_directory)
    Consurf_Score = []
    f = open("data/score.txt","r")
    all = f.readlines()
    for i in range(len(all)):
        if i in [0,1,2,3,4,5,6,7,8,9,10,11,12,len(all)-2,len(all)-1,len(all)]:
            pass
        else:
            consurf = all[i].split()[2]
            Consurf_Score.append(float(consurf))
    f.close()
    Consurf_Scores = []
    for i in position:
        Consurf_Scores.append(Consurf_Score[i-1])
    print("Success")
##############################################################################################################################
    # 使用foldx构建突变体的pdb
    Mut_PDB = FoldX
    from mutation.FoldX import run_FoldX
    from mutation.FoldX import get_total_energy
    run_FoldX(FoldX, WT_PDB, file_path)
    tte = get_total_energy(FoldX, WT_PDB)
    # 氨基酸网络
    from Feature_Cal.AAWeb import AAWEB
    relative_path = "tools/dssp.exe"
    absolute_path = os.path.abspath(relative_path)
    absolute_path = absolute_path.replace("\\", "/")
    #print("Calculating Amino Acid Web Features...", end=" ")
    print("Amino Acid Contact Network Building...", end=" ")
    for i in set_category:
        AAWEB(absolute_path, i, category, Mut_PDB, WT_PDB, "data/AAWeb")
    # AAWEB(absolute_path,category,Protein_name,Mut_PDB,WT_PDB,"data",position)
    from Feature_Cal.AAWeb import data_AAW_gener
    AAWeb_data = data_AAW_gener(position, category)
    print("Success")
##############################################################################################################################
    # 计算弹性网络参数
    from Feature_Cal.prody_cal import dynamics_dat
    dynamics = dynamics_dat(Protein_name, position,WT_PDB)
##############################################################################################################################
    df_all = pd.DataFrame()

    df_all["Disease"] = category
    df_all["Site"] = position
    df_all["Mutation"] = mutation_list

    df_all["Co.evolution"] = MI
    df_all["Entropy"] = SI
    df_all["Consurf_Score"] = Consurf_Scores
    df_all["RASA"] = RASA
    df_all["ddG"] = tte

    df_all["Betweenness"] = [sublist[0] for sublist in AAWeb_data]
    df_all["Closeness"] = [sublist[1] for sublist in AAWeb_data]
    df_all["Degree"] = [sublist[2] for sublist in AAWeb_data]
    df_all["Eigenvector"] = [sublist[3] for sublist in AAWeb_data]
    df_all["Clustering.coefficient"] = [sublist[4] for sublist in AAWeb_data]
    
    df_all["Effectiveness"] = [sublist[0] for sublist in dynamics]
    df_all["Sensitivity"] = [sublist[1] for sublist in dynamics]
    df_all["MSF"] = [sublist[2] for sublist in dynamics]
    df_all["DFI"] = [sublist[3] for sublist in dynamics]
    df_all["Stiffness"] = [sublist[4] for sublist in dynamics]

    # 将结果保存到paras.txt文件中
    df_all.to_csv("data/paras.txt", sep='\t',index=False)
    # df_all.to_csv("/home/wangjingran/APMA/Outcome/paras.txt",sep = '\t', index=False)
############################################################################################################################## 
    print("..Machine Learning Starting...")
    from ML import ML_Build
    ML_Build(category)
