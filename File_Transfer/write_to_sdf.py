import os
from rdkit import Chem

f = open('ligands.txt', 'r')
ligands = f.readlines()


os.makedirs("Ligand_SDF", exist_ok=True)
for i, ligand in enumerate(ligands):
    mol = Chem.MolFromSmiles(ligand)
    writer = Chem.SDWriter(os.path.join("Ligand_SDF", f"molecule_{i+1}.sdf"))
    if mol is not None:
        writer.write(mol)
    else:
        print(f"Skipped molecule {i+1} : {ligand}, invalid molecule")