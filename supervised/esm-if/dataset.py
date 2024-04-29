import torch
from torch.utils.data import Dataset
import numpy as np
import pandas as pd
import os
import pickle
import re
from tqdm import tqdm

amino_acid_index = { 'L': 0, 'A': 1, 'G': 2, 'V': 3, 'S': 4, 'E': 5, 'R': 6, 'T': 7, 'I': 8, 'D': 9, 'P': 10, 'K': 11, 'Q': 12, 'N': 13, 'F': 14, 'Y': 15, 'M': 16, 'H': 17, 'W': 18, 'C': 19}

class ProteinGymDataset(Dataset):
    def __init__(self, data_path, embedding_path, split_method="random", test_fold=0,ignore_files=set(),split="train"):

        self.split = split  
        # self.data_path = data_path
        self.data_path = data_path
        self.embedding_path = embedding_path
        if os.path.exists("./data/dataset.pth"):
            data = torch.load("./data/dataset.pth")
            self.train_data = data["train_data"]
            self.test_data = data["test_data"]
            return
        sequence_list = pd.read_csv("/scratch/user/atharvagashe22/Bioinformatics/PST/ProteinGYM/reference_files/DMS_substitutions.csv")
        self.train_data=[]
        self.test_data=[]
        err_seq_ids= {"A0A140D2T1_ZIKV_Sourisseau_2019", "BRCA2_HUMAN_Erwood_2022_HEK293T", "CAS9_STRP1_Spencer_2017_positive", "P53_HUMAN_Giacomelli_2018_Null_Etoposide", "P53_HUMAN_Giacomelli_2018_Null_Nutlin", "P53_HUMAN_Giacomelli_2018_WT_Nutlin",
        "POLG_HCVJF_Qi_2014", "POLG_CXB3N_Mattenberger_2021"}
        ignore_files = err_seq_ids
        dms_to_uniprot_dict = {}
        for i in tqdm(range(len(sequence_list))):
            DMS_id = sequence_list.iloc[i]["DMS_id"]
            UniProt_ID = sequence_list.iloc[i]["UniProt_ID"]
            sequence_wt = sequence_list.iloc[i]["target_seq"]
            dms_to_uniprot_dict[DMS_id] = UniProt_ID

        #load all csv files in the data path
        self.csv_files = [f for f in os.listdir(self.data_path) if f.endswith('.csv') and f.split(".")[0] not in ignore_files]
        #load all the embeddings
        embedding_fliles = {f for f in os.listdir(embedding_path) if f.endswith('_wt.npy')}

        folds= {"random" : "fold_random_5", "modulo":"fold_modulo_5", "contiguous": "fold_contiguous_5"}
        selected_fold = folds[split_method]
        for f in tqdm(self.csv_files, desc="Processing files"):
            df = pd.read_csv(f"{self.data_path}/{f}")
             #Split the mutant  into wild AA, position of mutation and mutant AA
            df[['wildAA', 'position', 'mutantAA']] = df['mutant'].str.extract(r'([A-Za-z])(\d+)([A-Za-z])')

            test_df = df[df[selected_fold] == test_fold]
            train_df = df[df[selected_fold] != test_fold]

            #calculate mean and std for z-score normalization
            train_mean = train_df['DMS_score'].mean()
            train_std = train_df['DMS_score'].std()

            test_mean = test_df['DMS_score'].mean()
            test_std = test_df['DMS_score'].std()

            f = f.split(".")[0]

            print(f)

            #check if the embedding exists for the corresponding csv file
            if f"{dms_to_uniprot_dict[f]}_wt.npy" in embedding_fliles:
                embedding = np.load(f"{embedding_path}/{dms_to_uniprot_dict[f]}_wt.npy")
                # seq_len x 500
                position=int(df["position"][0])
                train_label = np.zeros(20)
                test_label = np.zeros(20)
                train_mask = np.zeros(20)
                test_mask = np.zeros(20)
                inp = embedding[position-1]
                for index, row in df.iterrows():   
                    if position != row["position"]:
                        if not np.all(train_label == 0):
                            self.train_data.append([inp,train_label,train_mask])
                        if not np.all(test_label == 0):
                            self.test_data.append([inp,test_label,test_mask])
                        train_label = np.zeros(20)
                        test_label = np.zeros(20)
                        train_mask = np.zeros(20)
                        test_mask = np.zeros(20)
                        #inp [1 x 512]
                        inp = embedding[int(row["position"])-1]

                    position = row["position"]
                    aa_index = amino_acid_index[row["mutantAA"]]
                    if test_fold == row[selected_fold]:
                        test_label[aa_index] = (row["DMS_score"]-test_mean)/test_std
                        test_mask[aa_index] = 1
                    else:
                        train_label[aa_index] = (row["DMS_score"]-train_mean)/train_std
                        train_mask[aa_index] = 1

                if not np.all(train_label == 0):
                    self.train_data.append([inp,train_label, train_mask])
                if not np.all(test_label == 0):
                    self.test_data.append([inp,test_label, test_mask])
            # print(self.train_data)
        torch.save({'train_data': self.train_data, 'test_data': self.test_data}, './data/dataset.pth')
        print(len(self.train_data),len(self.test_data))

    def __getitem__(self, idx):
        if self.split == "train":
            return self.train_data[idx]
        else:
            return self.test_data[idx]
            
    def __len__(self):
        return len(self.train_data)


data_path = "/scratch/user/atharvagashe22/Bioinformatics/P2/data/substitutions_singles"
embedding_path = "/scratch/user/atharvagashe22/Bioinformatics/P2/embeddings"
dataset = ProteinGymDataset(data_path,embedding_path)