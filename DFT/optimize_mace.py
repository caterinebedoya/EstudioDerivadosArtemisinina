import torch
from ase import Atoms
from ase.optimize import LBFGS
from rdkit import Chem
from rdkit.Chem import AllChem
from mace.calculators import mace_mp
import sys
import os

def optimize_molecule_mace(smiles, output_file=None, fmax=0.01):
    # 1. Generate 3D conformer with RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    
    mol = Chem.AddHs(mol)
    res = AllChem.EmbedMolecule(mol, randomSeed=42)
    if res == -1:
         res = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
         if res == -1:
             raise ValueError(f"Could not embed molecule: {smiles}")
    
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except:
        pass

    conf = mol.GetConformer()
    positions = conf.GetPositions()
    species = [atom.GetSymbol() for atom in mol.GetAtoms()]

    # 2. Setup ASE Atoms object
    atoms = Atoms(symbols=species, positions=positions)

    # 3. Setup MACE calculator
    # Using 'medium' model as a standard choice for MACE-MP
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    
    # mace_mp might download the model on first use
    calculator = mace_mp(model="medium", device=device, default_dtype="float64") 
    
    atoms.calc = calculator

    # 4. Optimize with LBFGS
    opt = LBFGS(atoms, logfile=None) #aquí poner True para ver el progreso
    opt.run(fmax=fmax)

    # 5. Save to file
    if output_file:
        atoms.write(output_file)
    return atoms

if __name__ == "__main__":
    if len(sys.argv) != 3 and len(sys.argv) != 4:
        print("Usage: python optimize_mace.py <smiles> <output_file> [fmax]")
        sys.exit(1)
    
    smiles = sys.argv[1]
    output_file = sys.argv[2]
    fmax = float(sys.argv[3]) if len(sys.argv) == 4 else 0.01
    
    try:
        optimize_molecule_mace(smiles, output_file, fmax=fmax)
        # print(f"Optimization complete. Saved to {output_file}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
