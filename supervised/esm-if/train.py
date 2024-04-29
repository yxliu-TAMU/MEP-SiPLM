# Training script for the supervised learning model
from model import EsmIfMEP
from dataset import ProteinGymDataset
from torch.utils.data import DataLoader
from torch import optim
import torch
import os
import time


if __name__ == "__main__":

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    data_path = "/scratch/user/atharvagashe22/Bioinformatics/P2/data/substitutions_singles"
    embedding_path = "/scratch/user/atharvagashe22/Bioinformatics/P2/embeddings"    # model_path = "../models/esm1_t34_670M_UR50S"
    test_fold = 0
    dataset = ProteinGymDataset(data_path,embedding_path)
    dataloader = DataLoader(dataset,batch_size=32,shuffle=True)
    #create a model object
    model = EsmIfMEP()
    model.to(device)

    criterion = torch.nn.MSELoss(reduction='none')  # Compute per-element losses


    optimizer = optim.Adam(model.parameters(),lr=1e-4)
    # loss_fn = torch.nn.MSELoss()
    #train the model
    prev_val_loss = float("inf")
    val_loss_list = []
    for epoch in range(100):
        model.train()
        #start time
        start = time.time()
        for i, (embedding, labels, masks) in enumerate(dataloader):
            #move the data to the device
            embedding = embedding.float().to(device)
            labels = labels.float().to(device)
            masks = masks.to(device)

            optimizer.zero_grad()

            #forward pass
            outputs = model(embedding)

            # Compute per-element losses
            losses = criterion(outputs, labels)

            # Apply the mask - only keep losses for known labels
            masked_losses = losses * masks

            # Compute the final loss
            loss = masked_losses.sum() / masks.sum()
            #backward pass
            loss.backward()
            optimizer.step()
            # print(f"Epoch: {epoch}, Batch: {i}, Loss: {loss.item()}")
        # save model very 20 epochs
        if epoch % 20 == 0:
            torch.save(model.state_dict(),f"./checkpoint/train1h/esmif_mep_{epoch}_fold{test_fold}.pt")
        #evaluate the model
        end = time.time()
        print(f"Epoch: {epoch}, Time: {end-start}")
        model.eval()
        val_loss = 0
        for i, (embedding, labels, masks) in enumerate(dataloader):
            embedding = embedding.float().to(device)
            labels = labels.float().to(device)
            masks = masks.to(device)
            outputs = model(embedding)

            # Compute per-element losses
            losses = criterion(outputs, labels)

            # Apply the mask - only keep losses for known labels
            masked_losses = losses * masks

            # Compute the final loss
            val_loss += masked_losses.sum() / masks.sum()
        val_loss_list.append(val_loss)
        if val_loss < prev_val_loss:
            prev_val_loss = val_loss
            torch.save(model.state_dict(),f"./checkpoint/train1h/esmif_mep_best_fold{test_fold}.pt")
        torch.save(val_loss_list,f"./checkpoint/train1h/val_loss_list.pt")
        print(f"Validation Loss: {val_loss}")
    print(val_loss_list)