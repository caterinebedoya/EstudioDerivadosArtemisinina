# Estudio Computacional de Derivados de Artemisinina

Trabajo de grado en Física — Universidad de Antioquia

Estrategia computacional integrada para diseñar y caracterizar derivados de artemisinina con potencial actividad antimalárica, combinando modelos generativos de inteligencia artificial, filtros de similitud a fármaco y descriptores de reactividad basados en la Teoría del Funcional de la Densidad (DFT).

## 🎯 Contexto

La malaria continúa representando una amenaza para la salud pública mundial, agravada por la resistencia creciente de *Plasmodium* a los antipalúdicos disponibles y por las limitaciones farmacocinéticas de la artemisinina, fármaco de referencia en su tratamiento. Este trabajo aborda el diseño racional de análogos que conserven el farmacóforo activo.

## 🔬 Metodología

1. **Generación de derivados.** Con REINVENT 4.0 se generaron derivados tipo éter a partir del *scaffold* de artemisinina.
2. **Filtrado.** Las moléculas se filtraron según los criterios de Lipinski y Veber, obteniendo un conjunto final de 40 candidatos.
3. **Optimización geométrica.** Las estructuras se optimizaron con el potencial MACE.
4. **Descriptores de reactividad.** Con PySCF se calcularon descriptores electrónicos globales (HOMO/LUMO, electrofilicidad ω) y locales (funciones de Fukui f⁻ₖ, filicidad ω⁻ₖ).

## 📊 Resultados

- El conjunto generado recuperó de manera no dirigida a **arteméter** y **arteéter**, lo que valida la metodología.
- El análisis de correlación reveló independencia entre los descriptores clásicos y los cuánticos globales.
- Los descriptores locales mostraron que los sustituyentes perturban de forma moderada, pero no predecible, la reactividad del farmacóforo.
- Se concluye que el diseño racional de análogos de artemisinina requiere evaluación computacional sistemática de cada candidato.

## 📂 Estructura del repositorio

| Carpeta | Contenido |
|---|---|
| [`FiltrosLipinskiVeber/`](FiltrosLipinskiVeber/) | Generación de derivados y filtrado por criterios de similitud a fármaco. |
| [`DFT/`](DFT/) | Optimización geométrica y cálculo de descriptores de reactividad. |
| [`Resultados/`](Resultados/) | Datos, gráficas y análisis de correlación entre descriptores. |

## 🛠️ Herramientas

- **REINVENT 4.0**: generación de moléculas con modelos generativos.
- **RDKit**: manejo de estructuras y cálculo de descriptores moleculares.
- **MACE**: potencial de aprendizaje automático para optimización geométrica.
- **PySCF**: cálculos de estructura electrónica (DFT).
- **Python**, con pandas, NumPy y Matplotlib para el análisis y la visualización.


