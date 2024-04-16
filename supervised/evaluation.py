# Evaluate the performance of the model on the test set by spearman correlation coefficient and mean squared error.

from model import SaProtMEP
from dataset import ProteinGymDataset
from torch.utils.data import DataLoader
import torch
from scipy.stats import spearmanr

if __name__ == "__main__":

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    data_path = "../dataset/ProteinGym/mutated_sequences"
    model_path = "../models/esm1_t34_670M_UR50S"
    test_fold = 0
    dataset = ProteinGymDataset(data_path,split_method="random",test_fold=test_fold,model_path=model_path,split="test")
    dataloader = DataLoader(dataset,batch_size=32,shuffle=False)
    #create a model object
    model = SaProtMEP(model_path)
    model.load_state_dict(torch.load(f"checkpoint/saprot_mep_100_fold{test_fold}.pt"))
    model.to(device)

    model.eval()
    val_loss = 0
    predictions = {}
    true_labels = {}
    for i, (sequences, labels, mutation_sites,tokens,file) in enumerate(dataloader):
        
        tokens = tokens.to(device)
        labels = labels.to(device)
        output = model(tokens)
        # if file is a key of predictions, append the output to the list, otherwise create a new list
        if file in predictions:
            predictions[file].append(output)
            true_labels[file].append(labels)
        else:
            predictions[file] = [output]
            true_labels[file] = [labels]

    #calculate spearman correlation
    for file in predictions:
        predictions[file] = torch.cat(predictions[file]).cpu().detach().numpy()
        true_labels[file] = torch.cat(true_labels[file]).cpu().detach().numpy()
        spearman_correlation = spearmanr(predictions[file],true_labels[file])
        MSE = ((predictions[file] - true_labels[file])**2).mean()
        print(f"File: {file}, Spearman Correlation: {spearman_correlation.correlation}, MSE: {MSE}")
    print("Mean Spearman Correlation: ",sum([spearmanr(predictions[file],true_labels[file]).correlation for file in predictions])/len(predictions))
    print("Mean MSE: ",sum([((predictions[file] - true_labels[file])**2).mean() for file in predictions])/len(predictions))
    