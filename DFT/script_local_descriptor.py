import os
import time
import csv

os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["OMP_NUM_THREADS"] = "16"
os.environ["MKL_NUM_THREADS"] = "16"


from pyscf import gto, dft
from tqdm import tqdm


# =====================================
# Parámetros de entrada
# =====================================

xyz_file = "artemisinina.xyz"
output_csv = "fukui_mulliken_artemisinina.csv"
molecule_name = "artemisinina"

basis_set = "6-31g**"
functional = "blyp"


# =====================================
# Función para leer extxyz
# =====================================

def read_extxyz(filename):

    atoms = []

    with open(filename, "r") as f:

        lines = f.readlines()

        for line in lines[2:]:
            parts = line.split()

            if len(parts) >= 4:

                symbol = parts[0]

                x, y, z = map(
                    float,
                    parts[1:4]
                )

                atoms.append(
                    (symbol, (x, y, z))
                )

    return atoms



# =====================================
# Función cálculo DFT + Mulliken
# =====================================

def calcular_mulliken(atoms, charge, spin):


    mol = gto.M(

        atom=atoms,
        basis=basis_set,
        charge=charge,
        spin=spin,
        unit="Angstrom",
        verbose=4
    )


    mf = dft.RKS(mol)

    mf.xc = functional

    mf.max_cycle = 100

    mf.conv_tol = 1e-7



    # Barra SCF

    pbar = tqdm(
        total=mf.max_cycle,
        desc=f"SCF carga {charge}",
        unit="iter"
    )


    def callback(envs):

        pbar.update(1)


    mf.callback = callback


    mf.kernel()


    pbar.close()



    if not mf.converged:

        print(
            f"ADVERTENCIA: SCF no convergió para carga {charge}"
        )



    # Cargas Mulliken

    _, cargas = mf.mulliken_pop(
        verbose=0
    )


    return mol, cargas



# =====================================
# Inicio cronómetro
# =====================================

start_time = time.time()



# =====================================
# Leer geometría
# =====================================

atoms = read_extxyz(xyz_file)



# =====================================
# Cálculos N, N+1, N-1
# =====================================


print("\nCalculando molécula neutra N")

mol, q_N = calcular_mulliken(
    atoms,
    charge=0,
    spin=0
)



print("\nCalculando anión N+1")

_, q_N1 = calcular_mulliken(
    atoms,
    charge=-1,
    spin=1
)



print("\nCalculando catión N-1")

_, q_N_1 = calcular_mulliken(
    atoms,
    charge=1,
    spin=1
)



# =====================================
# Funciones de Fukui
# ΔN = 1
# =====================================

f_plus = q_N1 - q_N

f_minus = q_N - q_N_1



# =====================================
# Tiempo total
# =====================================

elapsed = time.time() - start_time



# =====================================
# Guardar CSV
# =====================================


with open(
    output_csv,
    "w",
    newline="",
    encoding="utf-8"
) as f:


    writer = csv.DictWriter(

        f,

        fieldnames=[

            "molecula",
            "indice_atomo",
            "elemento",
            "f_plus",
            "f_minus"

        ]

    )


    writer.writeheader()



    for i in range(mol.natm):


        writer.writerow({

            "molecula":
                molecule_name,


            "indice_atomo":
                i,


            "elemento":
                mol.atom_symbol(i),


            "f_plus":
                f_plus[i],


            "f_minus":
                f_minus[i]

        })



# =====================================
# Resultados
# =====================================


print("\n" + "="*60)

print("FUNCIONES DE FUKUI GUARDADAS")

print("="*60)


print(f"Archivo: {output_csv}")

print(f"Tiempo total: {elapsed:.2f} s")


print("\nMayor f+")

idx_plus = f_plus.argmax()

print(
    mol.atom_symbol(idx_plus),
    idx_plus,
    f_plus[idx_plus]
)



print("\nMayor f-")

idx_minus = f_minus.argmax()

print(
    mol.atom_symbol(idx_minus),
    idx_minus,
    f_minus[idx_minus]
)