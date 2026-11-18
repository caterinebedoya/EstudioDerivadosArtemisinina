import pandas as pd
import ase.db
from tqdm import tqdm
from optimize_mace import optimize_molecule_mace
import numpy as np
import os

def create_database():
    input_csv = 'DHA_top50_con_estereoquimica.csv'
    db_name = 'bio_molecules_smiles_esteroquimica.db'
    
    if not os.path.exists(input_csv):
        print(f"Error: {input_csv} not found.")
        return
        

    # 1. Load and Sample Data
    print("Loading and sampling data...")
    df = pd.read_csv(input_csv)
    
    if len(df) < 10:
        print("Warning: Not enough rows to sample 10. Using all.")
        sample_df = df
    else:
        sample_df = df

    # 2. Initialize ASE Database
    # If it exists, we append to it (or you could delete it first to start fresh)
    if os.path.exists(db_name):
        os.remove(db_name)
        
    db = ase.db.connect(db_name)

    columns_to_save = [
        'page', 'ID', 'OID', 'Chemical', 
        'Rep_Skin_YF', 'Rep_Skin_M', 'Rep_Clothes_YF', 'Rep_Clothes_M', 
        'SMILES_Quiral'
    ]

    print(f"Starting optimization and database population for {len(sample_df)} molecules...")

    # 3. Process Loop with Progress Bar
    for i, row in tqdm(sample_df.iterrows(), total=len(sample_df), desc="Processing"):
        smiles = row['SMILES_Quiral']
        
        try:
            # Optimize using MACE (returns ASE Atoms object)
            # We don't need an intermediate file, so output_file=None
            atoms = optimize_molecule_mace(smiles, output_file=None, fmax=0.01)

            # Prepare metadata dictionary
            kv_data = {}
            for col in columns_to_save:
                if col in row:
                    val = row[col]
                    # Convert numpy types to native python types for JSON serialization safety
                    if pd.isna(val):
                        val = None
                    elif isinstance(val, (np.integer, int)):
                        val = int(val)
                    elif isinstance(val, (np.floating, float)):
                        val = float(val)
                    else:
                        val = str(val)
                    
                    kv_data[col] = val

            # Write to database
            db.write(atoms, **kv_data)
            
        except Exception as e:
            print(f"\nFailed to process {smiles}: {e}")

    print(f"\nDatabase creation complete: {db_name}")
    print(f"Total entries in DB: {db.count()}")

if __name__ == "__main__":
    create_database()
