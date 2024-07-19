import pickle
from PyFingerprint.fingerprint import get_fingerprint

with open("model.pkl", 'rb') as file:
    model = pickle.load(file)

with open("removed_cols.pkl", 'rb') as file:
    cols_2_remove = pickle.load(file)


import numpy as np

def featurize(smiles):
    r = []
    j = ['topological-torsion', 'hybridization', 'morgan', 'extended', 'fp2']
    for t in j:
        finger = list(get_fingerprint(smiles, t).to_numpy())
        if t == 'rdk-descriptor':
            if finger[42] is not None and finger[42] > np.finfo(np.float32).max:
                finger[42] = np.finfo(np.float32).max
        r += finger
    return r

import pandas as pd
import csv
import os

filename = os.path.join("..", "Molecule_Generation", "latest_csv", "smiles.csv")

def get_X(smiles):
    global filename
    data = {'SMILES': smiles}

    labels = labels = [f'fp_{i}' for i in range(7168)]
    
    df = pd.DataFrame(data)
    df['features'] = df['SMILES'].apply(featurize)
    df = df.dropna(subset=['features'])
    features_df = pd.DataFrame(df['features'].tolist(), columns=labels)
    print(features_df.shape)
    cols_2_remove_existing = [col for col in cols_2_remove if col in features_df.columns]
    features_df = features_df.drop(columns=cols_2_remove_existing)
    print(features_df.shape)
    return features_df

smiles = []
with open(filename, 'r') as file:
    csv_reader = csv.reader(file)
    for row in csv_reader:
        smiles.append(row[0])

IC50s = model.predict(get_X(smiles).to_numpy())
import json
with open("dump", 'w') as f:
    json.dump(IC50s.tolist(), f, indent=4)

print(IC50s)
print(min(IC50s))