import os
from pathlib import Path
import pickle
import numpy as np
import pandas as pd
import requests
from proteinshake.datasets import Dataset
from proteinshake.utils import download_url, extract_tar
from tqdm import tqdm
import shutil
from itertools import islice


class MutationDataset(Dataset):
    exlude_args_from_signature = ["id"]

    def __init__(self, **kwargs):
        self.id = self.available_ids()
        kwargs["use_precomputed"] = False
        kwargs["root"] = str(kwargs["root"])
        super().__init__(**kwargs)
        self.mutations = [pd.read_csv(f"{self.root}/raw/{id}.csv") for id in self.id]
        #print(f"csv files:{self.mutations[0]}")
    @property
    def name(self):
        return f"{self.__class__.__name__}"

    def get_id_from_filename(self, filename):
        return filename.rstrip(".pdb")

    def download(self):
        """Implement me!
        This function downloads the raw data and preprocesses it in the following format, one for each protein:
        - {self.root}/raw/{self.id}.pdb # the reference pdb structure file
        - {self.root}/raw/{self.id}.csv # a table with columns 'mutations' (space separated), and 'y' (the measurements)
        """
        raise NotImplementedError


class DeepSequenceDataset(MutationDataset):
    file_path="/scratch/user/atharvagashe22/Bioinformatics/PST/ProteinGYM/reference_files/dms_dict7.pkl"
    with open(file_path, 'rb') as file:
    # Load the object from the file
        meta_data_raw = pickle.load(file)
    meta_data = meta_data_raw
    meta_data = dict(islice(meta_data_raw.items(), 86, 142))
    print(meta_data)

    # meta_data = {
    #     "BLAT_ECOLX_Tenaillon2013": "P62593",
    #     # "BLAT_ECOLX_Palzkill2012": "P62593",
    #     "DLG4_RAT_Ranganathan2012": "P31016",
    # }

    @classmethod
    def available_ids(cls):

        return list(cls.meta_data.keys())

    def get_raw_files(self):
        return [f"{self.root}/raw/{id}.pdb" for id in self.id]

    def download(self):
        pass
        # url = "https://dl.fbaipublicfiles.com/fair-esm/examples/variant-prediction/data/raw_df.csv"
        # download_url(url, f"{self.root}/raw/measurements.csv")
        # measurements = pd.read_csv(f"{self.root}/raw/measurements.csv")
        # for dataset, uniprot_id in tqdm(self.meta_data.items(), desc="Processing"):
        #     df = pd.read_csv(f"/scratch/user/atharvagashe22/Bioinformatics/PST/ProteinGYM/DMS_Substitutions/{dataset}.csv")
        #     # df = measurements[measurements["protein_name"] == dataset]
        #     df = df.rename(columns={"DMS_score": "y", "mutant": "mutations"})
        #     df["mutations"] = df["mutations"].map(lambda x: " ".join(x.split(":")))
        #     df = df[["mutations", "y"]]
        #     df.to_csv(f"{self.root}/raw/{dataset}.csv", index=False)

        #     source_path = f"/scratch/user/atharvagashe22/Bioinformatics/PST/ProteinGYM/AF2_structures/{uniprot_id}.pdb"

        #     # Specify the destination path where you want the file to be moved
        #     destination_path = f"{self.root}/raw/{dataset}.pdb"

        #     # Moving the file from source to destination
        #     shutil.move(source_path, destination_path)


            # url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
            # root_af2 = self.root.split("/")[0] + "/" + "non_af2_db"
            # root_af2 = f"{self.root}/non_af2_db"
            # if not Path(f"{root_af2}/{uniprot_id}.pdb").exists():
            #     try:
            #         download_url(url, f"{self.root}/raw/{dataset}.pdb")
            #     except requests.exceptions.HTTPError:
            #         print(f"Could not download {url}")

            # else:
            #     Path(f"{root_af2}/{uniprot_id}.pdb").rename(
            #         Path(f"{self.root}/raw/{dataset}.pdb")
            #     )


# class EnvisionDataset(MutationDataset):
#     def download(self):
#         pass


# class SkempiDataset(MutationDataset):
#     @classmethod
#     def available_ids(self):
#         return [
#             "1A22",
#             "1A4Y",
#             "1ACB",
#             "1AHW",
#             "1AK4",
#             "1BRS",
#             "1CBW",
#             "1CHO",
#             "1CSE",
#             "1DAN",
#             "1DQJ",
#             "1DVF",
#             "1E96",
#             "1EAW",
#             "1EMV",
#             "1F47",
#             "1FC2",
#             "1FCC",
#             "1FFW",
#             "1GC1",
#             "1GCQ",
#             "1H9D",
#             "1HE8",
#             "1IAR",
#             "1JCK",
#             "1JRH",
#             "1JTG",
#             "1KTZ",
#             "1LFD",
#             "1MAH",
#             "1MQ8",
#             "1NMB",
#             "1PPF",
#             "1R0R",
#             "1REW",
#             "1S1Q",
#             "1TM1",
#             "1UUZ",
#             "1VFB",
#             "1XD3",
#             "1Z7X",
#             "2FTL",
#             "2G2U",
#             "2I9B",
#             "2J0T",
#             "2JEL",
#             "2O3B",
#             "2PCB",
#             "2PCC",
#             "2SIC",
#             "2VLJ",
#             "2WPT",
#             "3BK3",
#             "3BN9",
#             "3HFM",
#             "3NPS",
#             "3SGB",
#             "4CPA",
#         ]

#     def get_raw_files(self):
#         return [f"{self.root}/raw/{self.id}.pdb"]

#     def download(self):
#         meta_url = "https://life.bsc.es/pid/mutation_database/SKEMPI_1.1.csv"
#         pdb_url = "https://life.bsc.es/pid/mutation_database/SKEMPI_pdbs.tar.gz"

#         data_path = os.path.join(self.root, "raw", "SKEMPI_1.1.csv")
#         pdb_path = os.path.join(self.root, "raw", "SKEMPI_pdbs.tar.gz")

#         label_keys = ["Affinity_mut (M)", "Affinity_wt (M)"]
#         mut_key = "Mutation(s)_cleaned"

#         download_url(meta_url, data_path)
#         download_url(pdb_url, pdb_path)
#         extract_tar(
#             os.path.join(self.root, "raw", "SKEMPI_pdbs.tar.gz"),
#             os.path.join(self.root, "raw"),
#         )

#         def mut_index(mut_num, mut_chain, protein_dict):
#             """Compute the index in the parsed sequence corresponding to the
#             mutated residue from SKEMPI"""
#             res_nums = protein_dict["residue"]["residue_number"]
#             chain_ids = protein_dict["residue"]["chain_id"]
#             try:
#                 hits = [
#                     i
#                     for i, (res, chain) in enumerate(zip(res_nums, chain_ids))
#                     if str(res) == mut_num and chain == mut_chain
#                 ]
#                 return hits[0] + 1
#             except IndexError:
#                 return None

#         mut_clean = (
#             lambda x, protein: x[0] + str(mut_index(x[2:-1], x[1], protein)) + x[-1]
#         )
#         df = pd.read_csv(data_path, sep=";")
#         rows = []
#         for protein, mutations in df.groupby("Protein"):
#             if len(mutations) < 5:
#                 continue
#             pdbid = protein.split("_")[0]
#             protein_dict = self.parse_pdb(
#                 os.path.join(self.root, "raw", pdbid + ".pdb")
#             )
#             if protein_dict is None:
#                 continue

#             mutations_col, y_col = [], []
#             for _, mutant in mutations.iterrows():
#                 mutations = [
#                     mut_clean(m, protein_dict) for m in mutant[mut_key].split(",")
#                 ]
#                 effect = np.log(mutant[label_keys[0]] / mutant[label_keys[1]])
#                 rows.append(
#                     {
#                         "pdbid": pdbid,
#                         "sequence_wt": protein_dict["protein"]["sequence"],
#                         "mutation": ":".join(mutations),
#                         "effect": effect,
#                     }
#                 )

#                 mutations_col.append(" ".join(mutations))
#                 y_col.append(effect)

#             pd.DataFrame(rows).to_csv(os.path.join(self.root, "measurements.csv"))
#             df = pd.DataFrame({"mutations": mutations_col, "y": y_col})
#             df.to_csv(f"{self.root}/raw/{pdbid}.csv", index=False)
