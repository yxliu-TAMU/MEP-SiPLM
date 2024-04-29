# Add a MLP to the SaProt's model class
# from transformers import EsmTokenizer, EsmForMaskedLM
from torch import nn
import torch
from torch.nn import functional as F

class EsmIfMEP(nn.Module):
    # SaProt model with a Multi-Layer Perceptron for regression task
    def __init__(self):
        super(EsmIfMEP, self).__init__()
        self.fc1 = nn.Linear(512, 128)
        self.fc2 = nn.Linear(128, 20)
    
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = self.fc2(x)
        return x 