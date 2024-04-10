import pickle
import pandas as pd
from tqdm import tqdm
import shutil
import os

sequence_list = pd.read_csv("/scratch/user/atharvagashe22/Bioinformatics/PST/ProteinGYM/reference_files/DMS_substitutions.csv")

err_seq_ids= {"A0A140D2T1_ZIKV_Sourisseau_2019", "BRCA2_HUMAN_Erwood_2022_HEK293T", "CAS9_STRP1_Spencer_2017_positive", "P53_HUMAN_Giacomelli_2018_Null_Etoposide", "P53_HUMAN_Giacomelli_2018_Null_Nutlin", "P53_HUMAN_Giacomelli_2018_WT_Nutlin",
"POLG_HCVJF_Qi_2014", "POLG_CXB3N_Mattenberger_2021"}
dms_to_uniprot_dict = {}
count = 0
for i in tqdm(range(len(sequence_list))):
    DMS_id = sequence_list.iloc[i]["DMS_id"]
    UniProt_ID = sequence_list.iloc[i]["UniProt_ID"]
    sequence_wt = sequence_list.iloc[i]["target_seq"]
    seq_len = sequence_list.iloc[i]["seq_len"]
    if DMS_id in err_seq_ids:
        print("Yes")
    if ~sequence_list.iloc[i]["includes_multiple_mutants"] and DMS_id not in err_seq_ids  :
        dms_to_uniprot_dict[DMS_id] = UniProt_ID
        count +=1
print(f"COUNT: {count}")
root="/scratch/user/atharvagashe22/Bioinformatics/PST/datasets/dms7"
dms_to_uniprot_dict_final = {}
directory_path = f"{root}"
if os.path.exists(directory_path):
    shutil.rmtree(directory_path)
os.makedirs(f"{directory_path}/raw")
for dataset, uniprot_id in tqdm(dms_to_uniprot_dict.items(), desc="Processing"):
            df = pd.read_csv(f"/scratch/user/atharvagashe22/Bioinformatics/PST/ProteinGYM/DMS_Substitutions/{dataset}.csv")
            df = df.rename(columns={"DMS_score": "y", "mutant": "mutations"})
            df["mutations"] = df["mutations"].map(lambda x: " ".join(x.split(":")))
            df = df[["mutations", "y"]]
            source_path = f"/scratch/user/atharvagashe22/Bioinformatics/PST/ProteinGYM/AF2_structures/{uniprot_id}.pdb"
            if os.path.exists(source_path):
                df.to_csv(f"{root}/raw/{dataset}.csv", index=False)
                destination_path = f"{root}/raw/{dataset}.pdb"
                shutil.copy(source_path, destination_path)
                dms_to_uniprot_dict_final[dataset] = uniprot_id
file_path="/scratch/user/atharvagashe22/Bioinformatics/PST/ProteinGYM/reference_files/dms_dict7.pkl"
print(dms_to_uniprot_dict_final)
with open(file_path, 'wb') as file:  # Note the 'wb' mode for writing in binary
    pickle.dump(dms_to_uniprot_dict_final, file)
