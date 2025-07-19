# AMR-Ecoli

# Evaluación de la precisión diagnóstica de la secuenciación genómica en la detección de resistencias a antibióticos en *Escherichia coli*

Repositorio asociado al Trabajo de Fin de Máster de Roser Ferragud Ferragud  
Máster en Investigación Biomédica, Universitat de València

## Descripción

Este repositorio contiene los scripts y plantillas utilizados para el análisis de datos en el TFM:  
*"Evaluación de la precisión diagnóstica de la secuenciación genómica en la detección de resistencias a antibióticos en Escherichia coli"*.  

El objetivo del estudio es evaluar la concordancia entre los perfiles genéticos de resistencia (detectados por WGS con AMRFinderPlus y ResFinder) y los fenotipos determinados por pruebas clásicas de sensibilidad antibiótica.  

Los scripts permiten:
- Recopilar y estandarizar datos genotípicos y fenotípicos de aislados de *E. coli* procedentes de publicaciones y bases de datos públicas (NCBI y BV-BRC).
- Filtrar genomas por calidad y especie.
- Determinar la concordancia entre determinantes genéticos de resistencia (detectados con AMRFinderPlus y ResFinder) y fenotipo.
- Calcular métricas diagnósticas (sensibilidad, especificidad, etc.) y generar las tablas y figuras principales del estudio.

## 📁 Estructura del repositorio

- `README.md`
- `scripts/`
  - `Script conffile def.R`
  - `eliminar duplicacions articles.R`
  - `NCBI_PATRIC – 2.R`
  - `estandarización MIC – 3.R`
  - `transformación MIC + análisis ResFinder.R`
  - `resultados amrf – 2.R`
  - `Filtros de calidad.R`
  - `corresp mic gen – 2.R`
  - `plot_MIC_distributions.R`
  - `análisis mutaciones.R`
  - `Dataset.R`
  - `Contexto muestras.R`
- `dataset/`
- `RFF-revised/`
- `conffile/`
  - `SOP/`
    - `example.configuration_file.choi2024.partA.v0.3.csv`
    - `example.configuration_file.choi2024.partB.v0.3.csv`
    - `pAST_SOP_RFF.docx`
  - `configuration_file.template.partA.v0.3.csv`
  - `configuration_file.template.partB.v0.3.csv`
- `genomes/`
- `Tables/`

**`scripts/`**: Código en R para las diferentes etapas del pipeline.  
**`dataset/`**: Listado con todos los bioprojects de E. coli.  
**`RFF-revised/`**: Estudios con las plantillas CSV completadas.  
**`conffile/`**: Plantillas CSV PartA y PartB, SOP con instrucciones y ejemplos completados.  
**`genomes/`**: Resultados de AMRFinderPlus, ResFinder y métricas de Sylph.  
**`Tables/`**: Resultados finales.


