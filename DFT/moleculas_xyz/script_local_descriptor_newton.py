import os
import csv
# ==========================
# Número de hilos CPU
# ==========================

os.environ["OPENBLAS_NUM_THREADS"] = "32"
os.environ["OMP_NUM_THREADS"] = "32"
os.environ["MKL_NUM_THREADS"] = "32"


from pyscf import gto, dft



# ==========================
# Archivo xyz
# ==========================

xyz_file = "molecule_42_42.xyz"
output_csv = "newton_molecule_42_42.csv"



# ==========================
# Parámetros DFT
# ==========================

basis = "6-31G**"

functional = "blyp"



# ==========================
# Función cálculo DFT
# ==========================

# def calcular_cargas(charge):


#     mol = gto.M(

#         atom=xyz_file,

#         basis=basis,

#         charge=charge,

#         spin=0 if charge == 0 else 1,

#         unit="Angstrom",

#         verbose=4
#     )


#     mf = dft.RKS(mol)


#     mf.xc = functional


#     # criterios de convergencia

#     mf.max_cycle = 150

#     mf.conv_tol = 1e-7



#     # ejecutar SCF

#     mf.kernel()



#     # --------------------------
#     # Verificar convergencia
#     # --------------------------

#     if mf.converged:

#         print(
#             f"SCF convergido (carga {charge})"
#         )

#     else:

#         print(
#             f"ADVERTENCIA: SCF NO convergió (carga {charge})"
#         )



#     # --------------------------
#     # Cargas Mulliken
#     # --------------------------

#     _, cargas = mf.mulliken_pop(
#         verbose=0
#     )


#     return mol, cargas



def calcular_cargas(charge):

    mol = gto.M(
        atom=xyz_file,
        basis="6-31+G**",
        charge=charge,
        spin=0 if charge == 0 else 1,
        unit="Angstrom",
        verbose=4
    )


    # Neutro RKS, cargados UKS
    if charge == 0:
        mf = dft.RKS(mol)
    else:
        mf = dft.UKS(mol)


    mf.xc = "BLYP"


    # ==========================
    # Parámetros convergencia
    # ==========================

    mf.max_cycle = 300

    mf.conv_tol = 1e-7

    mf.conv_tol_grad = 1e-5


    # Evitar oscilaciones

    mf.diis_space = 8

    mf.diis_damp = 0.5


    # Separar orbitales para estabilizar

    #mf.level_shift = 0.3


    # Mejor densidad inicial

    mf.init_guess = "minao"



    print("\nSCF para carga:", charge)

    mf = mf.newton()
    mf.kernel()



    if mf.converged:

        print("SCF convergido")

    else:

        print("SCF NO convergido")



    _, cargas = mf.mulliken_pop(verbose=0)


    return mol, cargas


# ==========================
# Estados electrónicos
# ==========================


print("\nCalculando molécula neutra N")

mol, q_N = calcular_cargas(
    0
)



print("\nCalculando anión N+1")

_, q_N1 = calcular_cargas(
    -1
)



print("\nCalculando catión N-1")

_, q_N_1 = calcular_cargas(
    +1
)




# ==========================
# Funciones de Fukui
# ΔN = 1
# ==========================


# f_plus = q_N1 - q_N


# f_minus = q_N - q_N_1

f_plus = q_N - q_N1

f_minus = q_N_1 - q_N


# ==========================
# Resultados
# ==========================


print("\nFunciones de Fukui (Mulliken)")
print("--------------------------------")
print("Atomo   Elemento       f+        f-")
print("--------------------------------")


for i in range(mol.natm):

    print(

        f"{i:3d}      "
        f"{mol.atom_symbol(i):2s}      "
        f"{f_plus[i]: .6f}   "
        f"{f_minus[i]: .6f}"

    )


with open(
    output_csv,
    mode="w",
    newline="",
    encoding="utf-8"
) as f:


    writer = csv.DictWriter(

        f,

        fieldnames=[

            "indice_atomo",
            "elemento",
            "f_plus",
            "f_minus"

        ]

    )


    writer.writeheader()



    for i in range(mol.natm):

        writer.writerow({

            "indice_atomo":
                i,

            "elemento":
                mol.atom_symbol(i),

            "f_plus":
                f_plus[i],

            "f_minus":
                f_minus[i]

        })



print("\n" + "="*60)

print("FUNCIONES DE FUKUI GUARDADAS")

print("="*60)

print(f"Archivo: {output_csv}")