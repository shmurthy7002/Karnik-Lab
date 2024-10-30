import pickle
from PyFingerprint.fingerprint import get_fingerprint

with open("new_model.pkl", 'rb') as file:
    model = pickle.load(file)

with open("new_removed_cols.pkl", 'rb') as file:
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

# smiles = []
# with open(filename, 'r') as file:
#     csv_reader = csv.reader(file)
#     for row in csv_reader:
#         smiles.append(row[0])

def remove_outliers(value):
    if(value < 5000):
        return value
    return None

df = pd.read_csv("ic50_mols.csv")
df = df.dropna(subset=["IC50 (nM)"])
df = df[pd.to_numeric(df['IC50 (nM)'], errors='coerce').notna()]
df['IC50 (nM)'] = df['IC50 (nM)'].astype(float)

df['IC50 (nM)'] = df['IC50 (nM)'].apply(remove_outliers)
df = df.dropna(subset=['IC50 (nM)'])
print(df.shape)
print(df.columns)
smiles = df['Ligand SMILES'].tolist()

df['IC50_nM_float'] = pd.to_numeric(df['IC50 (nM)'], errors='coerce')
df = df.dropna(subset=['IC50_nM_float'])
df['IC50 (M)'] = df['IC50_nM_float'] / 1e9
df['pIC50'] = -np.log10(df['IC50 (M)'])

y_pred = model.predict(get_X(smiles).to_numpy())
y_actual = df['pIC50'].tolist()

import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))

plt.scatter(y_actual, y_pred, label='Predicted vs True', color='blue', alpha=0.5)

# Line with slope 1 going through (0, 0)
min_val = min(min(y_actual), min(y_pred))
max_val = max(max(y_actual), max(y_pred))
plt.plot([min_val, max_val], [min_val, max_val], color='red', linestyle='--', label='Ideal Prediction')

plt.xlabel('Actual pIC50')
plt.ylabel('Predicted pIC50')
plt.title('Predicted vs Actual pIC50')
plt.legend()
plt.grid(True)
plt.show()

# import json
# with open("dump", 'w') as f:
#     json.dump(IC50s.tolist(), f, indent=4)

# print(IC50s)
# print(min(IC50s))