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
  - `Filtros de calidad.R`
  - `transformación MIC + análisis ResFinder.R`
  - `resultados amrf – 2.R`
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

## Requisitos
- R ≥ 4.2
- Paquetes R:
  - tidyverse
  - data.table
  - readxl
  - httr
  - AMR
  - epiR

## Datos externos
Para reproducir los análisis, es necesario descargar:
Genomas de AllTheBacteria.
Datos fenotípicos y metadatos de artículos y bases públicas (NCBI, BV-BRC).
Tabla de correspondencia entre 

## Cómo ejecutar el pipeline
1️. Colocar los archivos de configuración (partA y partB) y las tablas suplementarias de cada estudio en una carpeta que tenga el nombre del estudio. Colocar todas las carpetas en RFF-revised.  
2️. Ejecutar el script de estandarización: scripts/Script conffile def.R`. Se generará un archivo log con errores o avisos y un csv con todas las muestras y metadatos estandraizados (df_cruzado.csv) en la carpeta RFF_revised.  
3. Ejecutar el script `eliminar duplicacions articles.R` para eliminar muestras duplicadas. Se generarán dos archivos csv: df_cruzado_completo_articles.csv con todas las muestras y metadatos sin duplicaciones y ENA_accessions.articles.csv con los identificadores de los genomas de las muestras. Se guardarán en Tables/.  
4. Descargar los datos de muestras que contengan información fenotípica y genotípica del NCBI (https://www.ncbi.nlm.nih.gov/pathogens/ast#escherichia%20coli) y del BV-BRC (https://www.bv-brc.org/view/Bacteria/2#view_tab=amr&filter=and(keyword(Escherichia),keyword(coli))) en tsv. Del BV-BRC se debe incluir tanto la tabla de AMR phenotypes como la de Genomes para poder cruzar el identificador del genoma con el fenotipo.  
5. Ejecutar el script `NCBI_PATRIC – 2.R` para obetener los identificadores de muestras que no se encuentran en la colección de los artículos. Se generarán 5 archivos csv en Tables/: patric_ENA_accessions.csv (identificadorees del PATRIC), NCBI_ENA_accessions.csv (identificadores del NCBI), patric_ENA_accessions.nuevas.temp.csv, muestras_NCBI.nuevas.coli.csv (identificadores + metadatos NCBI) y muestras_patric.nuevas.csv (identificadores `+ metadatos BV-BRC).  
6. Ejecutar el script `estandarización MIC – 3.R` para unir los metadatos de los artículos, el NCBI y el BV-BRC y seleccionarc las muestras que tienen MIC disponible. Se generarán 4 csv en Tables/: muestras_NCBI.nuevas.coli.st.csv (datos del NCBI estandarizados como la colección de artículos), patric_metadata_st.nuevas.csv (datos del BV-BRC estandarizados como la colección de artículos), broth_microdil.csv (muestras con MIC obtenido por broth microdilution de los artículos, NCBI y BV-BRC), muestras_mic_articles.csv (recuento de antibióticos con MIC por estudio).  
7. Descargar en AllTheBacteria los genomas asociados a los identificadores de Tables/patric_ENA_accessions.csv, Tables/NCBI_ENA_accessions.csv y Tables/ENA_accessions.articles.csv.  
8. Ejecutar AMRFinder y ResFinderPlus con los genomas que se han obtenido. Obtener las métricas de calidad de los genomas con Sylph y assembly-stats. Guardar todos estos resultados en genomes/.  
9. Ejecutar el script `Filtros de calidad.R` para obtener los identificadores de las muestras que tienen genomas de calidad. Se generarán 8 csv en Tables/: resumen_sylph.csv, atb_resf.csv (antibióticos que se han detectados con ResFinder), atb_amrf.csv (antibióticos que se han detectado con AMRFinderPlus), sylph_all.csv (todas las muestras con información de si pasan o no los filtros de calidad), all_summary_def_qc_clean.csv (muestras que pasan los filtros de calidad del genoma), amrf_results_accession_cruzado.completo.csv (resultados de AMRFinderPlus solo de las meustras que pasan los filtros), resf_results_accession_cruzado.csv  (resultados de ResFinder solo de las meustras que pasan los filtros).  
10. Ejecutar el script `transformación MIC + análisis ResFinder.R`.  Transforma el MIC a fenotipo de resistencia (sensible/no sensible; no tiene en cuenta los <, >, <=, >= no interpretables). Se conservan solo los fenotipos que pasan los filtros de calidad genómicos. Determina la concordancia entre determinantes genéticos de ResFinder y fenotipo de resistencia y también las métricas de precisión diagnóstica. Se generan X archivos csv en Tables/: broth_microdil_sir_complete3.new.csv, count_comparison.clean3.csv, count_comparison.completa.clean3.csv, atb_included_ns3, resultados_detallados_resfinder_eucast.ecoli_collapsed3.csv (correpondencia fenotipo-genotipo), resultados_detallados_resfinder_ecoff.ecoli_collapsed3.csv (correpondencia fenotipo-genotipo), resultados_detallados_resfinder_clsi.ecoli_collapsed3.csv (correpondencia fenotipo-genotipo), tabla_resultados_resf_eucast3.csv (precisón diagnóstica), tabla_resultados_resf_ecoff3.csv (precisón diagnóstica), tabla_resultados_resf_clsi3.csv (precisón diagnóstica).  


  

bash
Copiar
Editar
Rscript scripts/estandarizacion_datos.R
3️⃣ Descarga los genomas desde AllTheBacteria o ENA y colócalos en data/genomes/.

4️⃣ Ejecuta el filtrado y control de calidad:

bash
Copiar
Editar
Rscript scripts/filtrado_genomas.R
5️⃣ Corre el análisis de concordancia genotipo ↔ fenotipo:

bash
Copiar
Editar
Rscript scripts/analisis_concordancia.R
6️⃣ Revisa las salidas en results/.

📊 Resultados esperados
Tablas de métricas diagnósticas (sensibilidad, especificidad, etc.).

Gráficas de concordancia para cada antibiótico.

Reportes de calidad de los genomas procesados.

