from rdkit import Chem
from rdkit.Chem import Draw
import sys
import os
import csv


def smiles_to_image(smiles, mol_number):
  mol = Chem.MolFromSmiles(smiles)
  img = Draw.MolToImage(mol)
  img.save(os.path.join("Images", f"molecule_{mol_number}.png"))

def delete_directory(directory_path):
  if os.path.exists(directory_path):
    for filename in os.listdir(directory_path):
      file_path = os.path.join(directory_path, filename)
      if os.path.isfile(file_path):
        os.remove(file_path)
      elif os.path.isdir(file_path):
        delete_directory(file_path)
    os.rmdir(directory_path)
  else:
    print(f"Directory '{directory_path}' does not exist.")


if(__name__ == "__main__"):
    if len(sys.argv) > 1:
        if(sys.argv[1] != "rm"):
            mol_num = sys.argv[1]
            smiles = ""
            with open(os.path.join("latest_csv", "smiles.csv"), newline="") as csvfile:
                reader = csv.reader(csvfile)
                smiles = list(reader)[int(mol_num)-1][0]
                print(f"Generating image for molecule {mol_num}: {smiles}")
            smiles_to_image(smiles, int(mol_num))
        else:
            delete_directory("Images")
            os.mkdir("Images")
    else:
        print("oops, this script takes a cmd argument")
        quit()
    print("Finished generating image! Check the 'Images' directory for your visual")