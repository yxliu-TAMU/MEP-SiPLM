# Protein Mutation Effect Prediction using structure information and protein language model

## Prerequisites:
A conda env with: Python, Pytorch, Pandas, Numpy, ESM

## Installation
> git clone https://github.com/yxliu-TAMU/MEP-SiPLM \
> Download dataset from Zenodo ([zenodo](https://zenodo.org/records/10951915?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjdmNDkzYjdjLWY3YzUtNGE1MC1hMGZhLWYyYmRkZWVkMDllMyIsImRhdGEiOnt9LCJyYW5kb20iOiJjMmM2MzVmZTY1YWYyY2JlYTE1YjBkMGI0NWJjNmQ3YSJ9.hx6zOm4OM-RnW4iMSUUlGulEhFbm5uCG3wT48V60nngr-a5dwEd7Z6sITZM7R2age66kDCQON3L3pXLZWccXgg))


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
