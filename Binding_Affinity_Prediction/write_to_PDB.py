import os
from rdkit import Chem
from rdkit.Chem import PDBWriter
from rdkit.Chem import SDWriter

input_sdf_file = 'ligands.sdf'

output_directory = 'ligands/'

os.makedirs(output_directory, exist_ok=True)

# Read the SDF file
supplier = Chem.SDMolSupplier(input_sdf_file)

# Initialize a counter for valid molecules
valid_molecule_count = 0

frv = []

# Iterate over molecules in the SDF file
for idx, mol in enumerate(supplier):
    if mol is not None:  # Ensure the molecule is not None
        # Check if the molecule has an IC50 value
        ic50 = mol.GetProp('IC50 (nM)') if mol.HasProp('IC50 (nM)') else None
        if ic50 and ic50[0] != ">" and ic50[0] != "<" and float(ic50) <= 1000:
            # Create an output file name
            output_file = os.path.join(output_directory, f"molecule_{valid_molecule_count+1}.sdf")
            
            # Write the molecule to an individual SDF file
            writer = SDWriter(output_file)
            writer.write(mol)
            # if Chem.SDMolSupplier(output_file)[0] is not None:
            #     frv.append(valid_molecule_count+1)
                
            writer.close()
            
            # Increment the valid molecule counter
            valid_molecule_count += 1
print(f"Extracted {valid_molecule_count} molecules with IC50 values from the SDF file and saved them to {output_directory}")
print(f"frv: {frv}")
