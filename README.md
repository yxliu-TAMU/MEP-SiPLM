# Protein Mutation Effect Prediction using structure information and protein language model

## Prerequisites:
A conda env with: Python, Pytorch, Pandas, Numpy, ESM

## Installation
> git clone https://github.com/yxliu-TAMU/MEP-SiPLM \
> Download dataset from Google Drive (https://drive.google.com/drive/folders/1zMvAAT9xHuAu6awTF8hnwa4zWDCyoXYe?usp=sharing)


## File tree

> --benchmark: scripts to evaluate the previous models performance\
> --data: scripts to preprocess the dataset\
> --dataset: ProteinGym dataset and related files.

## To Do:
```
1. 7 proteins' sequence and structure not match: seq_id: {A0A140D2T1_ZIKV_Sourisseau_2019, BRCA2_HUMAN_Erwood_2022_HEK293T, CAS9_STRP1_Spencer_2017_positive, P53_HUMAN_Giacomelli_2018_Null_Etoposide, P53_HUMAN_Giacomelli_2018_Null_Nutlin, P53_HUMAN_Giacomelli_2018_WT_Nutlin,
POLG_HCVJF_Qi_2014,}. skipped them for now.

2. Several sequence have multi-mutation sequences. Skipped them for now.
```