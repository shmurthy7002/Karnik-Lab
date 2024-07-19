import pandas as pd
import sys
import re

def clean_ic50(value):
    clean_value = re.sub(r'^[^\d]+', '', value)
    if float(clean_value) > 9999:
        clean_value = None
    return clean_value

def remove_outliers(value):
    if(value > int(sys.argv[1])):
        return None
    return value

df = pd.read_csv("ic50_mols.csv")

# removes rows that don't have ic50 values
df = df.dropna(subset=["IC50 (nM)"])

# Comment out below line to remove all rows with > and < signs
#df['IC50 (nM)'] = df['IC50 (nM)'].apply(clean_ic50)

df = df[pd.to_numeric(df['IC50 (nM)'], errors='coerce').notna()]
df['IC50 (nM)'] = df['IC50 (nM)'].astype(float)

df['IC50 (nM)'] = df['IC50 (nM)'].apply(remove_outliers)

rows_with_nan = df[df.isnull().any(axis=1)]
first = rows_with_nan['Ligand SMILES']
print(first.tolist())

df = df.dropna(subset=['IC50 (nM)'])

#df.to_csv("binding_affinities.csv", index=False)