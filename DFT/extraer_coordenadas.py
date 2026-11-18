import ase.db
from ase.io import write
import os

# Conectar a la base de datos
db = ase.db.connect('bio_molecules_smiles_esteroquimica.db')

# Crear carpeta para guardar archivos
output_dir = 'molecules_xyz'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    

print(f"Extrayendo {db.count()} moléculas...")

# Extraer cada molécula
for i in range(1, db.count() + 1):
    try:
        atoms = db.get_atoms(i)
        row = db.get(i)
        
        # Crear nombre del archivo
        chemical_name = row.Chemical if hasattr(row, 'Chemical') else f"molecule_{i}"
        filename = f"{output_dir}/{chemical_name}_{i}.xyz"
        
        # Guardar en formato XYZ
        write(filename, atoms)
        print(f"✓ Guardado: {filename}")
        
    except Exception as e:
        print(f"✗ Error en molécula {i}: {e}")

print(f"\nArchivos guardados en la carpeta: {output_dir}")