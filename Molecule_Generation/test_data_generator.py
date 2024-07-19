import random
import os
import pandas as pd
from rdkit import Chem

def get_smiles_from_csv():
    df = pd.read_csv(os.path.join("..", "Binding_Affinity_Prediction", "binding_affinities.csv"))
    df['IC50 (nM)'] = df['IC50 (nM)'].astype(float)

    # This makes sure we only grab the best SMILES
    df = df[df['IC50 (nM)'] <= 1]
    return df['Ligand SMILES'].tolist()

def write_list_to_file(file_name, smiles):
    with open(os.path.join("INPUT", f"{file_name}.smiles"), "w") as file:
        for item in smiles:
            file.write(f"{item}\n")


def check_valid_smiles(smiles_list):
    valid_smiles_count = 0
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            valid_smiles_count += 1
    return valid_smiles_count

smiles_list = get_smiles_from_csv()

# smiles_list = ["CCOC1=NC2=CC=CC(=C2N1CC3=CC=C(C=C3)C4=CC=CC=C4C5=NOC(=O)N5)C(=O)",
#                "CCOc2nc1cccc(C(=O)O)c1n2Cc5ccc(c3ccccc3c4nn[nH]n4)cc5",
#                "CCCCC1=NC=C(/C=C(CC2=CC=CS2)/C(O)=O)N1CC3=CC=C(C(O)=O)C=C3.O=S(O)(C)=O",
#                "CCCCC1=NC2(CCCC2)C(=O)N1CC3=CC=C(C=C3)C4=CC=CC=C4C5=NNN=N5",
#                "CCCCC1=NC(=C(N1CC2=CC=C(C=C2)C3=CC=CC=C3C4=NNN=N4)CO)Cl",
#                "CCCC1=NC(=C(N1CC2=CC=C(C=C2)C3=CC=CC=C3C4=NNN=N4)C(=O)O)C(C)(C)O",
#                "CCCC1=NC2=C(N1CC3=CC=C(C=C3)C4=CC=CC=C4C(=O)O)C=C(C=C2C)C5=NC6=CC=CC=C6N5C",
#                "CCCCC(=O)N(CC1=CC=C(C=C1)C2=CC=CC=C2C3=NNN=N3)C(C(C)C)C(=O)O"]

random.shuffle(smiles_list)

train_split = int(0.8 * len(smiles_list))
valid_split = int(0.9 * len(smiles_list))

train_smiles = smiles_list[:train_split]
valid_smiles = smiles_list[train_split:valid_split]
test_smiles = smiles_list[valid_split:]

print(f"Number of molecules inputted: {len(smiles_list)}")
print(f"Number of valid smiles: {check_valid_smiles(smiles_list)}")
print(f"Total smiles {len(train_smiles + valid_smiles + test_smiles)}")

#print(train_smiles)
#print(valid_smiles)
#print(test_smiles)

write_list_to_file("test", test_smiles)
write_list_to_file("train", train_smiles)
write_list_to_file("valid", valid_smiles)

for smile in smiles_list:
    if len(smile) == 0:
        print("UH OH -1-312i31241209-40921-49-1294-01294-29104-")
        quit()
    print(len(smile))
