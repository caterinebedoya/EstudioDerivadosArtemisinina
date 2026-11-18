from rdkit import Chem
from rdkit.Chem import AllChem

# 1. Definir el "Molde" (Template) original con estereoquímica
smiles_template = "C[C@H]1CC[C@@H]2[C@@H]([C@@H](O[C@@H]3[C@]24[C@@H]1CC[C@@](O3)(OO4)C)O)C"
template = Chem.MolFromSmiles(smiles_template)
template_3d = Chem.AddHs(template)

# Generar la conformación 3D del molde
AllChem.EmbedMolecule(template_3d, AllChem.ETKDGv3())
AllChem.MMFFOptimizeMolecule(template_3d)
# Quitamos los hidrógenos al molde para facilitar el emparejamiento con los derivados
template_core = Chem.RemoveHs(template_3d) 

# 2. Tu lista de 40 SMILES planos (Aquí te pongo 2 ejemplos de derivados en el OH)
derivados_planos = [
    "COC1OC2OC3(C)CCC4C(C)CCC(C1C)C24OO3"
   
]

resultados_con_estereoquimica = []

# 3. Procesar cada derivado
for i, smiles_plano in enumerate(derivados_planos):
    mol_plana = Chem.MolFromSmiles(smiles_plano)
    
    if mol_plana is None:
        print(f"Error leyendo el SMILES {i}")
        continue
        
    try:
        # Añadir hidrógenos al derivado
        mol_derivado = Chem.AddHs(mol_plana)
        
        # EL PASO MÁGICO: Forzar al derivado a tomar la forma 3D del molde
        # Esto empareja el esqueleto y genera coordenadas 3D correctas para los centros quirales
        AllChem.ConstrainedEmbed(mol_derivado, template_core)
        
        # Leer la estructura 3D y reasignar las etiquetas quirales (R/S)
        Chem.AssignStereochemistryFrom3D(mol_derivado)
        
        # Quitar los hidrógenos de nuevo para que el SMILES quede limpio
        mol_final = Chem.RemoveHs(mol_derivado)
        
        # Generar el nuevo SMILES con la estereoquímica restaurada
        nuevo_smiles_quiral = Chem.MolToSmiles(mol_final, isomericSmiles=True)
        
        resultados_con_estereoquimica.append(nuevo_smiles_quiral)
        print(f"Derivado {i+1} procesado: {nuevo_smiles_quiral}")
        
    except ValueError as e:
        print(f"No se pudo encajar el derivado {i+1} en el molde. Error: {e}")

