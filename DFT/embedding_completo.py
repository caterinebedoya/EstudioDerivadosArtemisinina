import csv
from rdkit import Chem
from rdkit.Chem import AllChem

# 1. Definir el "Molde" (Template) original con estereoquímica
smiles_template = "C[C@H]1CC[C@@H]2[C@@H]([C@@H](O[C@@H]3[C@]24[C@@H]1CC[C@@](O3)(OO4)C)O)C"
template_3d = Chem.AddHs(template)

# Generar la conformación 3D del molde y prepararlo
print("Preparando el andamio 3D de la Dihidroartemisinina...")
AllChem.EmbedMolecule(template_3d, AllChem.ETKDGv3())
AllChem.MMFFOptimizeMolecule(template_3d)
template_core = Chem.RemoveHs(template_3d) 

# Nombres de los archivos
archivo_entrada = "DHA_top50.csv"
archivo_salida = "DHA_top50_con_estereoquimica.csv"

datos_procesados = []
nombres_columnas = []

print(f"Leyendo derivados desde {archivo_entrada}...\n")

# 2. Leer el CSV original
with open(archivo_entrada, mode='r', encoding='utf-8') as csv_in:
    lector = csv.DictReader(csv_in)
    
    # Copiamos las columnas originales y añadimos nuestra nueva columna al final
    nombres_columnas = lector.fieldnames + ["SMILES_Quiral"]
    
    for i, fila in enumerate(lector):
        # Extraemos el SMILES plano de la primera columna
        smiles_plano = fila["SMILES"]
        mol_plana = Chem.MolFromSmiles(smiles_plano)
        
        if mol_plana is None:
            print(f"[{i+1}] Error leyendo el SMILES: {smiles_plano}")
            fila["SMILES_Quiral"] = "ERROR_LECTURA"
            datos_procesados.append(fila)
            continue
            
        try:
            # Añadir hidrógenos, encajar en el molde y calcular estereoquímica
            mol_derivado = Chem.AddHs(mol_plana)
            AllChem.ConstrainedEmbed(mol_derivado, template_core)
            Chem.AssignStereochemistryFrom3D(mol_derivado)
            
            # Limpiar y generar el nuevo SMILES
            mol_final = Chem.RemoveHs(mol_derivado)
            nuevo_smiles_quiral = Chem.MolToSmiles(mol_final, isomericSmiles=True)
            
            # Guardamos el éxito en la nueva columna de esta fila
            fila["SMILES_Quiral"] = nuevo_smiles_quiral
            print(f"[{i+1}] ¡Éxito! -> {nuevo_smiles_quiral}")
            
        except ValueError as e:
            print(f"[{i+1}] Fallo al encajar en 3D. Error: {e}")
            fila["SMILES_Quiral"] = "ERROR_EMBEDDING"
            
        # Añadimos la fila completa (con la nueva columna) a nuestra lista
        datos_procesados.append(fila)

# 3. Guardar el nuevo archivo CSV
print(f"\nGuardando resultados en {archivo_salida}...")

with open(archivo_salida, mode='w', newline='', encoding='utf-8') as csv_out:
    escritor = csv.DictWriter(csv_out, fieldnames=nombres_columnas)
    
    # Escribir cabeceras y todas las filas
    escritor.writeheader()
    escritor.writerows(datos_procesados)

print("¡Proceso terminado exitosamente!")