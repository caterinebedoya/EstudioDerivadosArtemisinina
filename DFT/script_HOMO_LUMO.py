import os
import time
import numpy as np

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

from pyscf import gto, dft
from tqdm import tqdm

# =====================================
# Parámetros de entrada
# =====================================
xyz_file = "artemisinina.xyz"   # archivo generado por MACE/ASE
basis_set = "6-31g**"
functional = "blyp"                 # puedes cambiar a "pbe0" o "b3lyp"
charge = 0                          # artemisinina neutra
spin = 0                            # singlete

# =====================================
# Función para leer formato extxyz
# =====================================
def read_extxyz(filename):
    """
    Lee un archivo tipo extxyz de ASE/MACE:
    - primera línea: número de átomos
    - segunda línea: propiedades/metadata
    - siguientes líneas: símbolo + posiciones + columnas extra
    Devuelve lista de átomos compatible con PySCF.
    """
    with open(filename, "r") as f:
        lines = f.readlines()

    n_atoms = int(lines[0].strip())
    comment = lines[1].strip()

    atoms = []
    for line in lines[2:2 + n_atoms]:
        parts = line.split()
        symbol = parts[0]
        x, y, z = map(float, parts[1:4])  # solo tomamos coordenadas
        atoms.append((symbol, (x, y, z)))

    return atoms, comment

# =====================================
# Inicio del cronómetro
# =====================================
start_time = time.time()

# =====================================
# Leer geometría
# =====================================
atoms, comment = read_extxyz(xyz_file)

print("=" * 60)
print(f"Archivo: {xyz_file}")
print(f"Comentario: {comment}")
print(f"Número de átomos: {len(atoms)}")
print("=" * 60)

# =====================================
# Definición de la molécula
# =====================================
mol = gto.M(
    atom=atoms,
    basis=basis_set,
    charge=charge,
    spin=spin,
    unit="Angstrom",
    verbose=4
)

# =====================================
# Configuración DFT
# =====================================
mf = dft.RKS(mol)
mf.xc = functional
mf.max_cycle = 100
mf.conv_tol = 1e-7

# =====================================
# Barra de progreso
# =====================================
pbar = tqdm(total=mf.max_cycle, desc="SCF Iterations", unit="iter")

def scf_callback(envs):
    pbar.update(1)

mf.callback = scf_callback

# =====================================
# Ejecutar cálculo
# =====================================
E_total = mf.kernel()
pbar.close()

if not mf.converged:
    print("\nADVERTENCIA: El SCF no convergió correctamente.")

# =====================================
# HOMO y LUMO
# =====================================
mo_energies = mf.mo_energy
mo_occ = mf.mo_occ

homo_idx = max(i for i, occ in enumerate(mo_occ) if occ > 0)
lumo_idx = homo_idx + 1

E_HOMO = mo_energies[homo_idx]
E_LUMO = mo_energies[lumo_idx]

# Conversión Hartree -> eV
Ha_to_eV = 27.211386245988

# =====================================
# Tiempo total
# =====================================
elapsed = time.time() - start_time

# =====================================
# Resultados
# =====================================
print("\n" + "=" * 60)
print("RESULTADOS PARA ARTEMISININA")
print("=" * 60)

print(f"Energía total : {E_total:.8f} Ha")
print(f"Energía total : {E_total * Ha_to_eV:.4f} eV")

print(f"\nHOMO : {E_HOMO:.8f} Ha")
print(f"HOMO : {E_HOMO * Ha_to_eV:.4f} eV")

print(f"\nLUMO : {E_LUMO:.8f} Ha")
print(f"LUMO : {E_LUMO * Ha_to_eV:.4f} eV")

gap_ha = E_LUMO - E_HOMO
gap_ev = gap_ha * Ha_to_eV

print(f"\nGap HOMO-LUMO : {gap_ha:.8f} Ha")
print(f"Gap HOMO-LUMO : {gap_ev:.4f} eV")

print(f"\nTiempo total : {elapsed:.2f} s")

print("\nOrbitales:")
for i, (energy, occ) in enumerate(zip(mo_energies, mo_occ)):
    label = ""
    if i == homo_idx:
        label = "  <-- HOMO"
    elif i == lumo_idx:
        label = "  <-- LUMO"

    print(
        f"MO {i+1:3d} | "
        f"E = {energy:12.6f} Ha | "
        f"{energy * Ha_to_eV:10.4f} eV | "
        f"occ = {occ}{label}"
    )