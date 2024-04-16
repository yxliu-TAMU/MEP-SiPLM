# Training script for the supervised learning model
from model import SaProtMEP
from dataset import ProteinGymDataset
from torch.utils.data import DataLoader
from torch import optim
import torch
import os
import time


if __name__ == "__main__":

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    data_path = "../dataset/ProteinGym/mutated_sequences"
    model_path = "../models/esm1_t34_670M_UR50S"
    test_fold = 0
    dataset = ProteinGymDataset(data_path,split_method="random",test_fold=test_fold,model_path=model_path)
    dataloader = DataLoader(dataset,batch_size=32,shuffle=True)
    #create a model object
    model = SaProtMEP(model_path)
    #Freeze the ESM model's parameters
    for param in model.esm.parameters():
        param.requires_grad = False
    model.to(device)

    optimizer = optim.Adam(model.parameters(),lr=1e-4)
    loss_fn = torch.nn.MSELoss()
    #train the model
    prev_val_loss = float("inf")
    val_loss_list = []
    for epoch in range(1000):
        model.train()
        #start time
        start = time.time()
        for i, (sequences, labels, mutation_sites,tokens,_) in enumerate(dataloader):
            #move the data to the device
            tokens = tokens.to(device)
            labels = labels.to(device)
            optimizer.zero_grad()
            #forward pass
            output = model(tokens)
            loss = loss_fn(output,labels)
            #backward pass
            loss.backward()
            optimizer.step()
            print(f"Epoch: {epoch}, Batch: {i}, Loss: {loss.item()}")
        # save model very 20 epochs
        if epoch % 20 == 0:
            torch.save(model.state_dict(),f"checkpoint/saprot_mep_{epoch}_fold{test_fold}.pt")
        #evaluate the model
        end = time.time()
        print(f"Epoch: {epoch}, Time: {end-start}")
        model.eval()
        val_loss = 0
        for i, (sequences, labels, mutation_sites,tokens,_) in enumerate(dataloader):
            tokens = tokens.to(device)
            labels = labels.to(device)
            output = model(tokens)
            val_loss += loss_fn(output,labels)
        if val_loss < prev_val_loss:
            prev_val_loss = val_loss
            torch.save(model.state_dict(),f"checkpoint/saprot_mep_best_fold{test_fold}.pt")
        print(f"Validation Loss: {val_loss}")
