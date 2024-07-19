import rdkit.Chem as Chem
#from rdkit import DataStructs

from collections import deque
import subprocess
import os
import csv
import time
import sys



def capture_output_line_by_line(command, n):
  last_lines = deque(maxlen=n)

  with subprocess.Popen(command, stdout=subprocess.PIPE, universal_newlines=True) as process:
    for line in process.stdout:
      last_lines.append(line.rstrip())

  return list(last_lines)


def write_smiles_to_csv(smiles_list, output_dir="csv_output"):
  os.makedirs(output_dir, exist_ok=True)

  filename = f"smiles.csv"
  filepath = os.path.join(output_dir, filename)
  with open(filepath, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    for smiles in smiles_list:
      writer.writerow([smiles])

  print(f"Successfully wrote SMILES to CSV file in {output_dir}")


def write_smiles_to_sdf(data_lines, output_filename_base="sdf_molecule", output_dir="sdf_output"):
  os.makedirs(output_dir, exist_ok=True)
  for i, line in enumerate(data_lines):
    smiles = line.strip()
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
      filename = f"{output_filename_base}_{i+1}.sdf"
      with Chem.SDWriter(os.path.join(output_dir, filename)) as writer:
        writer.write(mol)
    else:
      print(f"Warning: Skipping invalid SMILES string: {line}")
  print(f"Successfully wrote SMILES to SDF file in {output_dir}")


if(__name__ == "__main__"):
    if(len(sys.argv) != 2):
      print("This script takes 1 cmd arg, num of molecules 2 generate!")
      quit()

    n = sys.argv[1]
    print(f"BEGINNING OUTPUT, attempting to generate {n} molecules")
    start_time = time.time()
    molecules = capture_output_line_by_line("molecule_generation sample now " + n, int(n))
    molecules = list(set(molecules))

    print("starting to write to sdf")

    sample = time.time()
    write_smiles_to_sdf(molecules, output_dir="100k_sdf")
    print(f"finished writing to sdf, took {time.time()-sample} seconds")

    sample = time.time()
    print("starting to write to csv")
    write_smiles_to_csv(molecules, output_dir="100k_csv")
    print(f"finished writing to csv, took {time.time()-sample} seconds")

    elapsed_time = time.time()-start_time
    print(f"Succesfully wrote data, took {elapsed_time} seconds to complete!")
 