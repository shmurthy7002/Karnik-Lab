import os
from rdkit import Chem

input_directory = 'sdf_output/'
output_sdf_file = 'combined.sdf'

sdf_files = [f for f in os.listdir(input_directory) if f.endswith('.sdf')]

writer = Chem.SDWriter(output_sdf_file)

for sdf_file in sdf_files:
    file_path = os.path.join(input_directory, sdf_file)
    supplier = Chem.SDMolSupplier(file_path)
    
    for mol in supplier:
        if mol is not None:
            writer.write(mol)

# Close the writer
writer.close()

print(f"Combined {len(sdf_files)} SDF files into {output_sdf_file}")
