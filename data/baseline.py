import torch
from esm import Alphabet,pretrained
import pandas as pd
from tqdm import tqdm
from scipy.stats import spearmanr
import numpy as np
def label_row(row, sequence, token_probs, alphabet, offset_idx):
    wt, idx, mt = row[0], int(row[1:-1]) - offset_idx, row[-1]
    assert sequence[idx] == wt, "The listed wildtype does not match the provided sequence"
    wt_encoded, mt_encoded = alphabet.get_idx(wt), alphabet.get_idx(mt)
    score = token_probs[0, 1 + idx, mt_encoded] - token_probs[0, 1 + idx, wt_encoded]
    return score.item()

torch.hub.set_dir("/scratch/user/yxliu/.cache/torch/hub/checkpoints/")
sequence_list =  pd.read_csv("/scratch/user/yxliu/dataset/landscape/ProteinGym/reference_files/DMS_substitutions.csv")
#print(sequence_list.head())
spearman_list = []
for i in tqdm(range(len(sequence_list))):
    i = i+24
    DMS_id = sequence_list.iloc[i]["DMS_id"]
    sequence_wt = sequence_list.iloc[i]["target_seq"]
    seq_len = sequence_list.iloc[i]["seq_len"]
    if ~sequence_list.iloc[i]["includes_multiple_mutants"]:
        mutation_seqs = pd.read_csv("/scratch/user/yxliu/dataset/landscape/ProteinGym/substitution/"+DMS_id+".csv")
        torch.cuda.empty_cache()

        model,alphabet = pretrained.esm2_t33_650M_UR50D()
        repr_layers = [int("esm2_t33_650M_UR50D()".split('_')[1][1:])]
        device = torch.device('cuda')
        model = model.to(device)
        model.eval()
        batch_converter = alphabet.get_batch_converter()
        data = [("protein1",sequence_wt)]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        batch_tokens = batch_tokens.to(device)

        wt_marginals_list = []
        masked_marginals_list = []
        with torch.no_grad():
            result = model(batch_tokens,repr_layers, return_contacts=True)
            #print(result["logits"])
            wt_token_probs = torch.log_softmax(result["logits"],dim=-1) #wt-marginals
            all_token_probs = []
            #for i in tqdm(range(batch_tokens.size(1))):
                #batch_tokens_masked = batch_tokens.clone()
                #batch_tokens_masked[0, i] = alphabet.mask_idx
                #with torch.no_grad():
                    #token_probs = torch.log_softmax(model(batch_tokens_masked)["logits"], dim=-1)
                #all_token_probs.append(token_probs[:, i])
                #masked_token_probs = torch.cat(all_token_probs, dim=0).unsqueeze(0) #masked-marginals

            
            #print(token_probs.shape)
            #print(result["logits"].shape,result["representations"][6].shape,result["attentions"].shape,result["contacts"].shape)
            for i in tqdm(range(len(mutation_seqs))):
                mutant = mutation_seqs.iloc[i]["mutant"]
                wt,idx,mt = mutant[0],int(mutant[1:-1]),mutant[-1]
                assert sequence_wt[idx-1]==wt
                wt_encoded, mt_encoded = alphabet.get_idx(wt), alphabet.get_idx(mt)
                #print(wt_encoded,mt_encoded,wt_token_probs[0, idx, mt_encoded],wt_token_probs[0, idx, wt_encoded])
                wt_marginals_score = wt_token_probs[0, idx, mt_encoded] - wt_token_probs[0, idx, wt_encoded] #wt-marginals ESM-variants
                #masked_marginal_score = masked_token_probs[0, idx, mt_encoded] - wt_token_probs[0,idx,wt_encoded] #masked-marginals
                wt_marginals_list.append(wt_marginals_score.item())
                #masked_marginals_list.append(masked_marginal_score.item())
                """
                #pseudo-ppl is too time-consuming, skip it
                sequence_mt = sequence_wt[:idx-1] + mt + sequence_wt[(idx -1):]
                data = [("protein1", sequence_mt)]
                batch_converter = alphabet.get_batch_converter()
                atch_labels, batch_strs, batch_tokens = batch_converter(data)
                log_probs = []
                for i in range(1, len(sequence_mt) - 1):
                    batch_tokens_masked = batch_tokens.clone()
                    batch_tokens_masked[0, i] = alphabet.mask_idx
                    with torch.no_grad():
                        token_probs = torch.log_softmax(model(batch_tokens_masked.cuda())["logits"], dim=-1)
                    log_probs.append(token_probs[0, i, alphabet.get_idx(sequence_mt[i])].item())
                print(len(log_probs))
                pppl_score = sum(log_probs) #pseudo-ppl
                """
        #print(wt_marginals_list)
        #print(masked_marginals_list)
        wt_marginals_spearman = spearmanr(wt_marginals_list,mutation_seqs["DMS_score"])[0]
        #masked_marginals_spearman = spearmanr(masked_marginals_list,mutation_seqs["DMS_score"])[0]
        print(wt_marginals_spearman,seq_len)
    else:
        wt_marginals_spearman=np.nan
        masked_marginals_spearman=np.nan
        seq_len = sequence_list.iloc[i]["seq_len"]

    spearman_list.append([DMS_id,wt_marginals_spearman,seq_len])
spearman_df = pd.DataFrame(spearman_list,columns=["DMS_id","wt_marginals_spearman","masked_marginals_spearman","seq_len"])
spearman_df.to_csv("/scratch/user/yxliu/ecen766/course_project/data/esm-variant.csv",index=False)

