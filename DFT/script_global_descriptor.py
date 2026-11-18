import os
import time
import csv
import numpy as np

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

from pyscf import gto, dft
from tqdm import tqdm

# =====================================
# Parámetros de entrada
# =====================================
xyz_file = "artemisinina.xyz"
output_csv = "descriptores_artemisinina.csv"
molecule_name = "artemisinina"

basis_set = "6-31g**"
functional = "blyp"
charge = 0
spin = 0

# =====================================
# Función para leer formato extxyz
# =====================================
def read_extxyz(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
        parts = line.split()
        symbol = parts[0]
        x, y, z = map(float, parts[1:4])
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
    print("ADVERTENCIA: El SCF no convergió correctamente.")

# =====================================
# HOMO y LUMO
# =====================================
mo_energies = mf.mo_energy
mo_occ = mf.mo_occ

homo_idx = max(i for i, occ in enumerate(mo_occ) if occ > 0)
lumo_idx = homo_idx + 1

E_HOMO = mo_energies[homo_idx]
E_LUMO = mo_energies[lumo_idx]

# =====================================
# Descriptores globales de reactividad
# =====================================
potencial_quimico_ha = (E_HOMO + E_LUMO) / 2.0
dureza_quimica_ha = (E_LUMO - E_HOMO) / 2.0

# Electrofilicidad en Hartree
electrofilicidad_ha = (
    (potencial_quimico_ha ** 2) / (2.0 * dureza_quimica_ha)
    if dureza_quimica_ha != 0 else np.nan
)

# Conversión Hartree -> eV
Ha_to_eV = 27.211386245988

# Tiempo total
elapsed = time.time() - start_time

# =====================================
# Crear fila limpia para pandas/CSV
# =====================================
row = {
    "molecula": molecule_name,
    "archivo_entrada": xyz_file,
    "funcional": functional,
    "basis_set": basis_set,
    "carga": charge,
    "spin": spin,
    "energia_total_ha": E_total,
    "homo_ha": E_HOMO,
    "lumo_ha": E_LUMO,
    "homo_ev": E_HOMO * Ha_to_eV,
    "lumo_ev": E_LUMO * Ha_to_eV,
    "potencial_quimico_ha": potencial_quimico_ha,
    "potencial_quimico_ev": potencial_quimico_ha * Ha_to_eV,
    "dureza_quimica_ha": dureza_quimica_ha,
    "dureza_quimica_ev": dureza_quimica_ha * Ha_to_eV,
    "electrofilicidad_ha": electrofilicidad_ha,
    "electrofilicidad_ev": electrofilicidad_ha * Ha_to_eV if not np.isnan(electrofilicidad_ha) else np.nan,
    "convergido": mf.converged,
    "tiempo_total_s": elapsed
}

# =====================================
# Guardar en CSV limpio
# =====================================
with open(output_csv, mode="w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=row.keys())
    writer.writeheader()
    writer.writerow(row)

print("\n" + "=" * 60)
print("RESULTADOS GUARDADOS EN CSV")
print("=" * 60)
print(f"Archivo: {output_csv}")
print(f"HOMO: {E_HOMO:.8f} Ha")
print(f"LUMO: {E_LUMO:.8f} Ha")
print(f"Potencial químico: {potencial_quimico_ha:.8f} Ha")
print(f"Dureza química: {dureza_quimica_ha:.8f} Ha")
print(f"Electrofilicidad: {electrofilicidad_ha:.8f} Ha")