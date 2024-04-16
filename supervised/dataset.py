#write a pytorch dataset for the ProteinGym dataset
from torch import Dataset
import pandas as pd
import os
from transformers import EsmTokenizer
from foldseek_util import get_struc_seq

class ProteinGymDataset(Dataset):
    def __init__(self, data_path, split_method="random", test_fold=0,model_path="",foldseek_path = "",ignore_files=[],split="train"):
        self.split = split
        self.data_path = data_path
        #load all csv files in the data path
        self.csv_files = [f for f in os.listdir(data_path) if f.endswith('.csv') and f not in ignore_files]
        # create 5 fold cross validation
        self.train_data = []
        self.test_data = []
        df_reference = pd.read_csv("../dataset/ProteinGym/reference_files/DMS_substitutions.csv")
        self.tokenizer = EsmTokenizer.from_pretrained(model_path)
        for f in self.csv_files:
            structure_file = df_reference[df_reference["DMS_filename"]==f]["pdb_file"]
            structure_file = os.path.join("../dataset/ProteinGym/AF2_structures",structure_file)
            mutations = pd.read_csv(os.path.join(data_path, f))
            sequences = mutations['mutated_sequence'].values
            labels = mutations['DMS_score'].values
            mutation_sites = mutations['mutant'].values
            if split_method == "random":
                folds = mutations["fold_random_5"].values
            elif split_method == "modulo":
                folds = mutations["fold_modulo_5"].values
            elif split_method == "contiguous":
                folds = mutations["fold_contiguous_5"].values
            else:
                raise ValueError("split_method must be one of 'random', 'modulo', 'contiguous'")
            #split mutations into 5 folds by fold number
            for i in range(len(folds)):
                if folds[i] == test_fold:
                    self.test_data.append((sequences[i], labels[i], mutation_sites[i],self.tokenize(structure_file,mutation_sites[i]),f))
                else:
                    self.train_data.append((sequences[i], labels[i], mutation_sites[i],self.tokenize(structure_file,mutation_sites[i]),f))

    def __len__(self):
        return len(self.train_data)
    
    def tokenize(self,structure_path,mutation_sites):
        #tokenize the sequence into a list of integers
        mutation_sites = mutation_sites.split(":")
        mutation_sites = [int(site[1:-1]) for site in mutation_sites]
        parsed_seqs = get_struc_seq("bin/foldseek",structure_path)
        _,_,combined_seq = parsed_seqs
        #apply mask to the mutation sites
        combined_seq = list(combined_seq)
        for site in mutation_sites:
            combined_seq[2*site+1] = "#"
        return self.tokenizer.tokenize(combined_seq)
    
    def __getitem__(self, idx):
        if self.split == "train":
            return self.train_data[idx]
        else:
            return self.test_data[idx]
