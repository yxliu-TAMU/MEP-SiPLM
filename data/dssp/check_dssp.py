import pandas as pd
import os
reference_files = "/scratch/user/yxliu/dataset/landscape/ProteinGym/reference_files/DMS_substitutions.csv"
dataset = pd.read_csv(reference_files)
structure_folder = "/scratch/user/yxliu/dataset/landscape/ProteinGym/AF2_structures/"

for i in range(len(dataset)):
    id = dataset["UniProt_ID"][i]
    print(id)
    #check if this id exists in the structure folder
    if not os.path.exists(structure_folder+id+".pdb"):
        print(f"ID: {id} does not exist in the structure folder")

#print(dataset["UniProt_ID"])
#with open('dssp.output', 'r') as f:
    #lines = f.readlines()
#for line in lines:
    #seq = None
    #dssp_seq = line.split(' ')[0]
    #id = line.split(' ')[1].split("/")[-1].split(".")[0]
    #print(id)
    #seq = dataset[dataset["UniProt_ID"]==id]["target_seq"].values[0]
    #if len(dssp_seq) != len(seq):
        #print(f"ID: {id}")
        #print(f"Original sequence: {len(seq)}")
        #print(f"DSSP sequence: {len(dssp_seq)}")
    #print(dssp_seq,id)
    
#print(f"Total number of sequences: {len(lines)}")
