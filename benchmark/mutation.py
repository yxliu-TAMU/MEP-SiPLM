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
        meta_data_raw = pickle.load(file)
    meta_data = meta_data_raw
    meta_data = dict(islice(meta_data_raw.items(), 86, 142))

    @classmethod
    def available_ids(cls):

        return list(cls.meta_data.keys())

    def get_raw_files(self):
        return [f"{self.root}/raw/{id}.pdb" for id in self.id]

    def download(self):
        pass

