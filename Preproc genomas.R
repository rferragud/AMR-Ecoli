library(tidyverse)
library(data.table)
library(dplyr)
library(readODS)
library(epiR)

rm(list=ls())

# 1. Preparar las tablas para que esten listas para el análisis: sylph, amrf, resf ####

# Arreglo del sylph para poner si las muestras han sido corridas en AMRfinder y Resfinder 

#resultados sylph sin yes/no 
all_summary_def_qc_sylph.raw <- read.csv("C:/Users/rferragud/Documents/Projecte/genomes/all_summary_def_qc_syplh_amrf_resf.new.csv") %>% select(-amrfplus, -resfinder)


  
accessions_articles <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/ENA_accessions.articles.csv", sep =",") 

accessions_articles_sam <- accessions_articles %>% pull(ENA_sample_accession)
  
accessions_NCBI <- read.csv("C:/Users/rferragud/Documents/Projecte/accessions/NCBI_ENA_accessions.csv", sep =",")

accessions_NCBI_sam <- accessions_NCBI %>% pull(ENA_sample_accession)

accessions_PATRIC <- read.csv("C:/Users/rferragud/Documents/Projecte/Scripts/scripts bons/patric_ENA_accessions.csv", sep =",")

accession_PATRIC_sam <- accessions_PATRIC %>% filter(!is.na(ENA_sample_accession)) %>% pull(ENA_sample_accession)
accession_PATRIC_run <- accessions_PATRIC %>% filter(!is.na(ENA_run_accession)) %>% pull(ENA_run_accession)

all_summary_def_qc_sylph.buenas <- all_summary_def_qc_sylph.raw %>% 
  filter(ENA_sample_accession %in% accessions_articles_sam |
           ENA_sample_accession %in% accessions_NCBI_sam |
           ENA_sample_accession %in% accession_PATRIC_sam |
           ENA_run_accession %in% accession_PATRIC_run) %>% distinct()

all_summary_def_qc_sylph.buenas %>% group_by(ENA_sample_accession) %>% count(ENA_sample_accession) %>% filter(n>1)
# resultados amrfplus

all_summary_def_qc_sylph.buenas %>% distinct(ENA_sample_accession)
amrfplus_results <- read.delim("C:/Users/rferragud/Documents/Projecte/genomes/amrfplus.concatenated.all.txt", sep = "\t")

#resultados de resistencias en resfinder (tabla ya arreglada)
resfinder_results_ecoli <- read.csv("C:/Users/rferragud/Documents/Projecte/genomes/all_samples.pheno_table_escherichia_coli.txt.gz", sep = "\t", header = F) %>%
  mutate(across(
    where(is.character),
    ~ na_if(., "")
  )) %>%
  rename(Biosample = "V1", Antibiotic = "V2", Class = "V3", WGS_predicted_phenotype = "V4", Match = "V5", genetic_background ="V6") %>%
  filter(!is.na(Antibiotic) & !Antibiotic == "parC;;1;;CP084529.1_108_v_AA" & !Antibiotic == "parC;;1;;CP084529.1_108_t_AA")

resfinder_results_general.raw <- read.csv("C:/Users/rferragud/Documents/Projecte/genomes/all_samples.pheno_table.txt.gz", sep = "\t", header = F) %>%
  mutate(across(
    where(is.character),
    ~ na_if(., "")
  )) %>%
  rename(Biosample = "V1", Antibiotic = "V2", Class = "V3", WGS_predicted_phenotype = "V4", Match = "V5", genetic_background ="V6") %>%
  filter(!is.na(Antibiotic) & !Antibiotic == "parC;;1;;CP084529.1_108_v_AA" & !Antibiotic == "parC;;1;;CP084529.1_108_t_AA")

# samples que se han corrido en amrfplus
amrfplus_tested <- 
  read.csv("C:/Users/rferragud/Documents/Projecte/genomes/amrfplus.tested.all.txt", sep = " ", header = F) %>% 
  mutate(Sample = gsub("_amrfplus.txt", "", V1)) %>%
  #have to edit assembly names to keep only accession
  mutate(Sample = ifelse(grepl("^GC[AF]", Sample), sub("^(GC[AF]_[^_]*).*", "\\1", Sample), Sample)) %>%
  select(Sample)

#samples que se han corrido en resfinder
resfinder_tested <- resfinder_results_ecoli %>%
  mutate(Sample = Biosample) %>% 
  distinct(Sample)

#sylph con yes/no en amrfplus y resfinder
all_summary_def_qc_sylph <- 
  all_summary_def_qc_sylph.buenas %>% # Este sería el archivo del summary que está en SACO
  mutate(
    amrfplus = ifelse(
      ENA_run_accession %in% amrfplus_tested$Sample |
        ENA_sample_accession %in% amrfplus_tested$Sample |
        assembly %in% amrfplus_tested$Sample,
      "yes",
      "no"
    )
  ) %>%
  mutate(
    resfinder = ifelse(
      ENA_run_accession %in% resfinder_tested$Sample |
        ENA_sample_accession %in% resfinder_tested$Sample |
        assembly %in% resfinder_tested$Sample,
      "yes",
      "no"
    )
  ) 


resumen_sylph <- all_summary_def_qc_sylph %>%
  count(source, assembly_pipeline, sylph, amrfplus, resfinder) 

write.csv(resumen_sylph, "C:/Users/rferragud/Documents/Projecte/Tables/resumen_sylph.csv", row.names = FALSE)

# antiboticos que se han detectado para e coli en resfinder
atb_resf_ecoli <- resfinder_results_ecoli%>% distinct(Antibiotic) %>% arrange(Antibiotic) %>% pull()

# antibioticos que se han detectado en general en resfinder
atb_resf_general <- resfinder_results_general.raw %>% distinct(Antibiotic) %>% arrange(Antibiotic) %>% pull(Antibiotic)

# antibioticos en amrfplus
atb_amrfplus <- amrfplus_results %>% select(Class, Subclass) %>% separate_rows(Subclass, sep = "/") %>% distinct() %>% arrange(Subclass)

# TODO tabla con todos los antibioticos para volver a hacer la correspondencia 
write.csv(atb_resf_general, "C:/Users/rferragud/Documents/Projecte/Tables/atb_resf.csv") 

write.csv(atb_amrfplus, "C:/Users/rferragud/Documents/Projecte/Tables/atb_amrf.csv")


# seleccionar los que no estan ya en e. coli
resfinder_results_general <- resfinder_results_general.raw %>% filter(!Antibiotic %in% atb_resf_ecoli) 

# arreglar resfinder para que quede la tabla de pointfinder unida y se pueda usar con la funcion que diseñamos
pointfinder.raw <- read.csv("C:/Users/rferragud/Documents/Projecte/genomes/all_samples.PointFinder_results.txt.gz", sep= "\t")

pointfinder <- pointfinder.raw %>% separate_rows(Resistance, sep = ", ")

point_mut_resf <- resfinder_results_ecoli %>% filter(WGS_predicted_phenotype == "Resistant" & is.na(genetic_background))


point_finder_renamed <- pointfinder %>%
  rename(
    Biosample  = Sample,
    Antibiotic = Resistance
  ) %>% 
  mutate(             
    Antibiotic = str_to_lower(str_trim(Antibiotic))
  ) 

point_finder_complete <- point_mut_resf %>%
  left_join(point_finder_renamed, by = c("Biosample", "Antibiotic"))

resf_no_resist <- resfinder_results_ecoli %>% filter(WGS_predicted_phenotype == "No resistance")

resf_gens <- resfinder_results_ecoli %>% filter(WGS_predicted_phenotype == "Resistant" & !is.na(genetic_background))

resfinder_results.complete <-bind_rows(point_finder_complete, resf_no_resist, resf_gens)

no_unidas <- point_finder_renamed %>% anti_join(point_mut_resf, by = c("Biosample", "Antibiotic") )

no_unidas %>% count(Antibiotic)

resfinder_results.complete %>% distinct(Biosample)
resfinder_results_ecoli %>% distinct(Biosample)

# 2. Pasar filtros de calidad a los genomas ####

# TODO comprobar que estas son las muestras que tocan (quitamos del patric algunas de campilobacter)

#Have a sylph call with at least 99 percent minimum abundance. 
#If a sample has more than one call (eg where it has more than one run), then require all species calls to be the same
#Minimum checkm2 completeness of 90%
#Maximum checkm2 contamination of 5%
#Total assembly length between 100kbp and 15Mbp
#Maximum number of contigs 2,000
#Minimum N50 2,000

#demasiado generales--> intentamos encontrar nuestros criterios de exclusión                                              

# 1. Filtrar las muestras que no tienen assembly pipeline ####
all_summary_def_qc_sylph.filtro1 <- 
  all_summary_def_qc_sylph %>% 
  mutate(
    filtro1 = if_else(
      is.na(assembly_pipeline),
      "no assembly pipeline",
      NA_character_
    )
  )

all_summary_def_qc_sylph.filtro2 <- all_summary_def_qc_sylph.filtro1 %>%
  mutate(
    filtro1 = if_else(
      is.na(filtro1) & str_detect(Comments, fixed("Biosample with more than 1 run acc")),
      "Biosample with more than 1 run acc",
      filtro1
    )
  )

all_summary_def_qc_sylph.filtro3 <- all_summary_def_qc_sylph.filtro2 %>%
  mutate(
    filtro1 = if_else(
      is.na(filtro1) & ((amrfplus == "no") | (resfinder == "no")),
      "no amrfplus or resfinder",
      filtro1
    )
  )


# 1. Ver de qué especies son los genomas que tenemos y filtrar los de E. coli ####

# Gráfica con las especie de los genomas 
all_summary_def_qc_sylph_species.temp <- 
  all_summary_def_qc_sylph.filtro3 %>% 
  mutate(
    species = word(Sylph_Contig_name, 2, 3))%>% 
  mutate(
    species = fct_infreq(species)
  )

all_summary_def_qc_sylph_species <- all_summary_def_qc_sylph_species.temp %>% 
  filter(is.na(filtro1)) %>%
  filter(!is.na(species))

especies_genomas <- ggplot(all_summary_def_qc_sylph_species, aes(x = species)) +
  geom_bar(fill  = "steelblue", 
           color = "black") +
  labs(
    x = "Especie",
    y = "Número de genomas",
    title = "Frecuencia de genomas por especie"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#getwd()

#ggsave("C:/Users/rferragud/Documents/Projecte/genomes/especies_genomas.png") 

#ggsave(
#filename = "C:/Users/rferragud/Documents/Projecte/genomes/especies_genomas.png",
#plot     = especies_genomas,    # tu objeto ggplot
#device   = "png")

# misma gráfica pero sin E. coli (mejor visualización de las otras especies)
all_summary_def_qc_sylph_species2 <- all_summary_def_qc_sylph %>% mutate(
  species = word(Sylph_Contig_name, 2, 3) 
)%>% 
  mutate(
    species = fct_infreq(species)
  )%>%
  filter(!species == "Escherichia coli")%>%
  filter(!is.na(species))

especies_genomas_scoli <- ggplot(all_summary_def_qc_sylph_species2, aes(x = species)) +
  geom_bar(fill  = "steelblue", 
           color = "black") +
  labs(
    x = "Especie",
    y = "Número de genomas",
    title = "Frecuencia de genomas por especie"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

especies_genomas_scoli

# df con solo genomas de E. coli
all_summary_def_qc_coli <- all_summary_def_qc_sylph_species %>% filter(species == "Escherichia coli")


all_summary_def_qc_sylph.filtro4 <- all_summary_def_qc_sylph_species.temp %>%
  mutate(
    filtro1 = if_else(
      is.na(filtro1) & is.na(Sylph_Taxonomic_abundance),
      "no sylph",
      filtro1
    ) 
  )


all_summary_def_qc_sylph.filtro5 <- all_summary_def_qc_sylph.filtro4 %>%
  mutate(
    filtro1 = if_else(
      is.na(filtro1) & !species == "Escherichia coli",
      "el genoma no pertenece a E. coli",
      filtro1
    )
  )



# 2. De los genomas que son E.coli miramos la abundancia taxonómica (%)

# gráfica que muestra el % de taxonomic abundance de E.coli
#1) Tus datos + rank + length_mb
all_summary_def_qc_arrange.sylph <- all_summary_def_qc_sylph.filtro5 %>%
  filter(is.na(filtro1)) %>%
  arrange(Sylph_Taxonomic_abundance) %>%
  mutate(
    rank_sylph = row_number(),
    sylph   = Sylph_Taxonomic_abundance
  )

# 2) Calcula media y sd; umbrales -2 SD y -3 SD
stats.sylph <- all_summary_def_qc_arrange.sylph %>%
  summarise(
    mean_sylph = mean(sylph, na.rm = TRUE),
    sd_sylph   = sd(sylph,   na.rm = TRUE)
  )
mean_sylph <- stats.sylph$mean_sylph
sd_sylph <- stats.sylph$sd_sylph

lower_2sd <- mean_sylph - 2 * sd_sylph
lower_3sd <- mean_sylph - 3 * sd_sylph

# 3) Cuenta outliers ±2 SD (opcional)
outlier_count <- all_summary_def_qc_arrange.sylph %>%
  filter(sylph < lower_2sd) %>%
  nrow()

# 4) Dibuja el gráfico con todas las líneas y etiquetas
g_sylph <-
  ggplot(all_summary_def_qc_arrange.sylph,
         aes(x = rank_sylph, y = sylph)) +
  geom_line() +
  geom_point(size = 1) +
  
  # líneas ±2 SD
  geom_hline(yintercept = lower_2sd, linetype = "dashed") +
  
  # líneas ±3 SD
  # geom_hline(yintercept = lower_3sd, linetype = "dotdash") +
  
  # etiquetas ±2 SD
  annotate("text",
           x     = max(all_summary_def_qc_arrange.sylph$rank_sylph) * 0.75,
           y     = lower_2sd,
           label = paste0("Media - 2 SD = ", round(lower_2sd, 2), "%"),
           vjust = 2) +
  
  # etiquetas ±3 SD
  #annotate("text",
  #         x     = max(all_summary_def_qc_arrange.sylph$rank_sylph) * 0.75,
  #         y     = lower_3sd,
  #         label = paste0("Mean - 3 SD = ", round(lower_3sd, 2)),
  #         vjust = 3) +
  
  # etiqueta de outliers ±2 SD (opcional)
  annotate("text",
           x     = max(all_summary_def_qc_arrange.sylph$rank_sylph),
           y     = max(all_summary_def_qc_arrange.sylph$sylph),
           hjust = 1, vjust = 4.5,
           label = paste0(outlier_count, " genomas fuera de rango")) +
  
  labs(x = "Genomas", y = "Abundancia taxonómica (%)") +
  theme_minimal()

g_sylph


# Vemos que siendo más conservadores si cortamos por -2SD, solo quitamos 137 genomas
# haremos eso para asegurarnos que no hay contaminaciones 

all_summary_def_qc_sylph.filtro6 <- all_summary_def_qc_sylph.filtro5 %>%
  mutate(
    filtro1 = if_else(
      is.na(filtro1) & Sylph_Taxonomic_abundance < lower_2sd,
      paste0("abundancia taxonómica menor a ", round(lower_2sd, 2), "%"),
      filtro1
    )
  )

#comprobación 
all_summary_def_qc_coli.abundance <- all_summary_def_qc_coli %>% filter(Sylph_Taxonomic_abundance > lower_2sd)


# 3. Mirar que el tamaño total de las lecturas sea el adecuado

# 1) Ordenar QC_total_length y establecer las dos variables a comparar
all_summary_def_qc_arrange.length <- all_summary_def_qc_sylph.filtro6 %>%
  filter(is.na(filtro1)) %>%
  arrange(QC_total_length) %>%
  mutate(
    rank      = row_number(),
    length_mb = QC_total_length / 1e6
  ) 

# 2) Calcula media y sd; umbrales ±2 SD y ±3 SD
stats <- all_summary_def_qc_arrange.length %>%
  summarise(
    mean_len = mean(length_mb, na.rm = TRUE),
    sd_len   = sd(length_mb,   na.rm = TRUE)
  )
mean_len <- stats$mean_len
sd_len   <- stats$sd_len

upper_2sd <- mean_len + 2 * sd_len
lower_2sd <- mean_len - 2 * sd_len
upper_3sd <- mean_len + 3 * sd_len
lower_3sd <- mean_len - 3 * sd_len

# 3) Cuenta outliers ±2 SD (opcional)
outlier_count <- all_summary_def_qc_arrange.length %>%
  filter(length_mb > upper_2sd | length_mb < lower_2sd) %>%
  nrow()

# 4) Dibuja el gráfico con todas las líneas y etiquetas
g_total_length <- 
  ggplot(all_summary_def_qc_arrange.length,
         aes(x = rank, y = length_mb)) +
  geom_line() +
  geom_point(size = 1) +
  
  # líneas ±2 SD
  geom_hline(yintercept = upper_2sd, linetype = "dashed") +
  geom_hline(yintercept = lower_2sd, linetype = "dashed") +
  
  # líneas ±3 SD
  #geom_hline(yintercept = upper_3sd, linetype = "dotdash") +
  #geom_hline(yintercept = lower_3sd, linetype = "dotdash") +
  
  # etiquetas ±2 SD
  annotate("text",
           x     = max(all_summary_def_qc_arrange.length$rank) * 0.75,
           y     = upper_2sd,
           label = paste0("Mean + 2 SD = ", round(upper_2sd, 2), " Mb"),
           vjust = -1) +
  annotate("text",
           x     = max(all_summary_def_qc_arrange.length$rank) * 0.75,
           y     = lower_2sd,
           label = paste0("Mean - 2 SD = ", round(lower_2sd, 2), " Mb"),
           vjust = 1.5) +
  
  # etiquetas ±3 SD
  #annotate("text",
  #         x     = max(all_summary_def_qc_arrange.length$rank) * 0.75,
  #         y     = upper_3sd,
  #         label = paste0("Mean + 3 SD = ", round(upper_3sd, 2), " Mb"),
  #         vjust = -2,
  #         fontface = "italic") +
  #annotate("text",
  #         x     = max(all_summary_def_qc_arrange.length$rank) * 0.75,
  #         y     = lower_3sd,
  #         label = paste0("Mean - 3 SD = ", round(lower_3sd, 2), " Mb"),
  #         vjust = 2.5,
  #         fontface = "italic") +
  
  # etiqueta de outliers ±2 SD (opcional)
  annotate("text",
           x     = max(all_summary_def_qc_arrange.length$rank),
           y     = max(all_summary_def_qc_arrange.length$length_mb),
           hjust = 1, vjust = 1,
           label = paste0(outlier_count, " genomas fuera de rango")) +
  
  labs(x = "Genomas", y = "Tamaño total (Mb)") +
  theme_minimal()

g_total_length

#Nos quedamos de momento con lo más conservador ±2 SD y eliminamos 575 genomas 

all_summary_def_qc_sylph.filtro7 <- all_summary_def_qc_sylph.filtro6 %>%
  mutate(
    filtro2 = if_else(
      is.na(filtro1) & ((QC_total_length / 1e6) < lower_2sd | (QC_total_length/ 1e6) > upper_2sd),
      paste0("tamaño total mayor que ", round(lower_2sd, 2), "Mb"),
      NA_character_
    )
  )

#comprobación
all_summary_def_qc_coli.abundance.length <- 
  all_summary_def_qc_coli.abundance %>% 
  filter((QC_total_length / 1e6) > lower_2sd | (QC_total_length / 1e6) < upper_2sd)



# 3. Mirar que el contig number sea el adecuado (arreglar todo, usar el filtrado)

# Gráfico de puntos  

# 1) Tus datos + rank + length_mb
all_summary_def_qc_sylph_number2 <- all_summary_def_qc_sylph.filtro7 %>%
  filter(is.na(filtro1)) %>%
  arrange(QC_number) %>%
  mutate(
    rank_number = row_number(),
    QC_number   = QC_number
  )

all_summary_def_qc_sylph_number2 %>% filter(QC_number == "0")

# 2) Calcula media y sd; umbrales ±2 SD y ±3 SD
stats.number <- all_summary_def_qc_sylph_number2 %>%
  summarise(
    mean_number = mean(QC_number, na.rm = TRUE),
    sd_number   = sd(QC_number,   na.rm = TRUE)
  )
mean_number <- stats.number$mean_number
sd_number   <- stats.number$sd_number
upper_2sd <- mean_number + 2 * sd_number
lower_2sd <- mean_number - 2 * sd_number
upper_3sd <- mean_number + 3 * sd_number
lower_3sd <- mean_number - 3 * sd_number

# 3) Cuenta outliers ±2 SD (opcional)
outlier_count <- all_summary_def_qc_sylph_number2 %>%
  filter(QC_number > upper_2sd) %>%
  nrow()

# 4) Dibuja el gráfico con todas las líneas y etiquetas
g_number <- ggplot(all_summary_def_qc_sylph_number2,
                   aes(x = rank_number, y = QC_number)) +
  geom_line() +
  geom_point(size = 1) +
  
  # líneas ±2 SD
  geom_hline(yintercept = upper_2sd, linetype = "dashed") +
  
  # líneas ±3 SD
  #geom_hline(yintercept = upper_3sd, linetype = "dotdash") +
  
  # etiquetas ±2 SD
  annotate("text",
           x     = max(all_summary_def_qc_sylph_number2$rank_number) * 0.75,
           y     = upper_2sd,
           label = paste0("Mean + 2 SD = ", round(upper_2sd, 2), "contigs"),
           hjust = 0.7, vjust = -1) +
  
  # etiquetas ±3 SD
  #annotate("text",
  #         x     = max(all_summary_def_qc_sylph_number2$rank_number) * 0.75,
  #         y     = upper_3sd,
  #         label = paste0("Mean + 3 SD = ", round(upper_3sd, 2), "contigs"),
  #         hjust =0.7, vjust = -1.5) +
  
  # etiqueta de outliers +2 SD (opcional)
  annotate("text",
           x     = max(all_summary_def_qc_sylph_number2$rank_number),
           y     = max(all_summary_def_qc_sylph_number2$QC_number),
           hjust = 1, vjust = 1,
           label = paste0(outlier_count, " genomas fuera de rango")) +
  
  labs(x = "Genomas", y = "Número de contigs") +
  theme_minimal()

g_number


all_summary_def_qc_sylph.filtro8 <- all_summary_def_qc_sylph.filtro7 %>%
  mutate(
    filtro3 = if_else(
      is.na(filtro1) & (QC_number > upper_2sd),
      paste0("número de contigs mayor que ", round(upper_2sd, 2)),
      NA_character_
    )
  )


# Gráfico qc_N50 (NO LO USAMOS)

# 1) Tus datos + rank + length_mb
all_summary_def_qc_arrange.n50 <- all_summary_def_qc_sylph.filtro8 %>%
  filter(is.na(filtro1)) %>%
  arrange(QC_N50) %>%
  mutate(
    rank_n50 = row_number(),
    n50_mb   = QC_N50 / 1e6
  )


# 2) Calcula media y sd; umbrales ±2 SD y ±3 SD
stats.n50 <- all_summary_def_qc_arrange.n50 %>%
  summarise(
    mean_n50 = mean(n50_mb, na.rm = TRUE),
    sd_n50   = sd(n50_mb,   na.rm = TRUE)
  )
mean_n50 <- stats.n50$mean_n50
sd_n50   <- stats.n50$sd_n50

upper_2sd <- mean_n50 + 2 * sd_n50
lower_2sd <- mean_n50 - 2 * sd_n50
upper_3sd <- mean_n50 + 3 * sd_n50
lower_3sd <- mean_n50 - 3 * sd_n50

# 3) Cuenta outliers ±2 SD (opcional)
outlier_count <- all_summary_def_qc_arrange.n50 %>%
  filter(n50_mb > upper_2sd) %>%
  nrow()

# 4) Dibuja el gráfico con todas las líneas y etiquetas
g_n50 <- ggplot(all_summary_def_qc_arrange.n50,
                aes(x = rank_n50, y = n50_mb)) +
  geom_line() +
  geom_point(size = 1) +
  
  # líneas ±2 SD
 #geom_hline(yintercept = upper_2sd, linetype = "dashed") +
 #geom_hline(yintercept = lower_2sd, linetype = "dashed") +
 #
 ## líneas ±3 SD
 #geom_hline(yintercept = upper_3sd, linetype = "dotdash") +
 #geom_hline(yintercept = lower_3sd, linetype = "dotdash") +
 #
 ## etiquetas ±2 SD
 #annotate("text",
 #         x     = max(all_summary_def_qc_arrange.n50$rank_n50) * 0.75,
 #         y     = upper_2sd,
 #         label = paste0("Mean + 2 SD = ", round(upper_2sd, 2), " Mb"),
 #         vjust = -0.5) +
 #annotate("text",
 #         x     = max(all_summary_def_qc_arrange.n50$rank) * 0.75,
 #         y     = lower_2sd,
 #         label = paste0("Mean - 2 SD = ", round(lower_2sd, 2), " Mb"),
 #         vjust = 1.5) +
  
  
  
  # etiqueta de outliers ±2 SD (opcional)
 # annotate("text",
 #          x     = max(all_summary_def_qc_arrange.n50$rank_n50),
 #          y     = max(all_summary_def_qc_arrange.n50$n50_mb),
 #          hjust = 1, vjust = 1,
 #          label = paste0(outlier_count, " genomas fuera de ±2 SD")) +
  
  labs(x = "Genomas", y = "N50 (Mb)") +
  theme_minimal(base_size = 16)

g_n50


# 3. Seleccionar solo las muestras que pasan los filtros ####

all_summary_def_qc_sylph.filtro8 %>% distinct(filtro1)

resumen <- all_summary_def_qc_sylph.filtro8 %>% count(filtro1, filtro2, filtro3)

todos_niveles <- c(
  "no assembly pipeline",
  "Biosample with more than 1 run acc",
  "no amrfplus or resfinder",
  "no sylph",
  "el genoma no pertenece a E. coli",
  "abundancia taxonómica menor a 98.47%"
)

resumen_ord <- resumen %>%
  mutate(
    filtro1 = fct_relevel(filtro1, !!!todos_niveles)
  )

resumen_ord <- resumen_ord %>%
  arrange(filtro1)

write.csv(resumen_ord, "C:/Users/rferragud/Documents/Projecte/genomes/resumen_sylph.csv", row.names = FALSE)
write.csv(all_summary_def_qc_sylph.filtro8, "C:/Users/rferragud/Documents/Projecte/genomes/sylph_all.csv", row.names = FALSE)

all_summary_def_qc_sylph.filtro8 <- read.csv("C:/Users/rferragud/Documents/Projecte/genomes/sylph_all.csv", check.names = FALSE)

# df con las muestras que pasan el QC
all_summary_def_qc_clean.temp <- all_summary_def_qc_sylph.filtro8 %>% filter(is.na(filtro1) & is.na(filtro2) & is.na(filtro3))

# eliminamos las filas repetidas 
all_summary_def_qc_clean <- all_summary_def_qc_clean.temp %>% distinct()

all_summary_def_qc_clean %>% distinct(ENA_sample_accession)

write.csv(all_summary_def_qc_clean, "C:/Users/rferragud/Documents/Projecte/Tables/all_summary_def_qc_clean", row.names = F)

# 4. Cruzar las tablas con los resultados de resf y amrf con la del sylph #### 
# para quedarnos con las muestras que pasan los filtros 

## 4.1. Cruzar AMRFinder con sylph ####

# cambiamos en el df de resfinder el assembly para que esté igual que en all_summary_def_qc_clean
amrf_results_assembly <- amrfplus_results %>% 
  mutate(Name = str_replace(Name, "^(([^_]+_[^_]+)).*", "\\1"))

# df con solo las accessions de all_summary_def_qc_clean_accession
all_summary_def_qc_clean_accession <- 
  all_summary_def_qc_clean %>% 
  select(ENA_sample_accession, ENA_run_accession, assembly)

# cruzamos resfinder con las accessions: ENA_sample_accession 
amrf_results_accession_cruzado <- 
  all_summary_def_qc_clean_accession %>% 
  filter(!is.na(ENA_sample_accession)) %>%
  left_join(amrf_results_assembly, by = c("ENA_sample_accession" = "Name"))

# samples que se han cruzado
sam_amrf <- amrf_results_accession_cruzado %>% pull(ENA_sample_accession) %>% unique()

# cruzamos resfinder con las accessions:ENA_run_accession 
#amrf_results_accession_cruzado1 <- 
#  all_summary_def_qc_clean_accession %>%
#  filter(!ENA_sample_accession %in% sam_amrf) %>%
#  filter(!is.na(ENA_run_accession))%>%
#  left_join(amrf_results_assembly, by = c("ENA_run_accession" = "Name"))
#
## run accession que se han cruzado
#run_amrf <- amrf_results_accession_cruzado1 %>% pull(ENA_run_accession) %>% unique()
#
## cruzamos resfinder con las accessions:assembly
#amrf_results_accession_cruzado2 <- 
#  all_summary_def_qc_clean_accession %>% 
#  filter(!ENA_sample_accession %in% sam_amrf) %>%
#  filter(!ENA_run_accession %in% run_amrf) %>%
#  filter(!is.na(assembly))%>%
#  left_join(amrf_results_assembly, by = c("assembly" = "Name"))
#
## assemblies que se han cruzado
#assembly_amrf <- amrf_results_accession_cruzado2 %>% pull(assembly) %>% unique()
#
## muestras que no se han cruzado -> no cumplen los criterios de inclusion
#sobrantes <- 
#  anti_join(amrf_results_assembly, all_summary_def_qc_clean_accession, by= c("Name" = "ENA_sample_accession")) %>% 
#  distinct(Name) %>% 
#  pull()
#
#all_summary_def_qc_remove <- all_summary_def_qc_sylph.filtro8 %>% filter(!is.na(filtro1) | !is.na(filtro2) | !is.na(filtro3)) %>% distinct()
#
#all_summary_def_qc_remove %>% filter(ENA_sample_accession %in% sobrantes)  
#all_summary_def_qc_remove %>% filter(ENA_run_accession %in% sobrantes)
#all_summary_def_qc_remove %>% filter(assembly %in% sobrantes)
#
#
#sobrantes <- 
#  anti_join(amrf_results_assembly, all_summary_def_qc_clean_accession, by= c("Name" = "ENA_sample_accession")) %>% 
#  distinct(Name)
#
#sam_sylph <- all_summary_def_qc_remove %>% distinct(ENA_sample_accession) %>% pull()
#run_sylph <- all_summary_def_qc_remove %>% distinct(ENA_run_accession) %>% pull()
#assembly_sylph <- all_summary_def_qc_remove %>% distinct(assembly) %>% pull()
#
#sobrantes %>% filter(Name %in% sam_sylph)
#sobrantes %>% filter(Name %in% run_sylph)
#sobrantes %>% filter(!Name %in% assembly_sylph & !Name %in% sam_sylph & !Name %in% run_sylph)
#

# df con todas las muestras que se han cruzado (solo tiene ya muestras que pasan el QC)
#amrf_results_accession_cruzado.completo <- 
#  bind_rows(amrf_results_accession_cruzado, amrf_results_accession_cruzado1, amrf_results_accession_cruzado2) %>%
#  rename(genetic_background = Element.symbol) %>%
#  mutate(
#    Subclass = str_replace_all(
#      Subclass,
#      regex("CLAVULANIC_ACID", ignore_case = TRUE),
#      "clavulanic acid"
#    ))

# # comprobación de que las muestras que se han unido y las que no son el total de las all_summary_def_qc_clean_accession 
# amrf_results_accession_cruzado.completo %>% filter(!is.na(Contig.id)) %>% pull(ENA_sample_accession) %>% unique()
# # 19020 (con resist)
# amrf_results_accession_cruzado.completo %>% filter(is.na(Contig.id)) %>% pull(ENA_sample_accession)%>% unique()
# # 3617 (no tienen genes de resistencia)
# all_summary_def_qc_clean_accession %>% pull(ENA_sample_accession)%>% unique()
# # 22637
# head(amrf_results_accession_cruzado.completo)

# comprobación de que las muestras que se han unido y las que no son el total de las all_summary_def_qc_clean_accession 
amrf_results_accession_cruzado %>% filter(!is.na(Contig.id)) %>% pull(ENA_sample_accession) %>% unique()
# 19020 (con resist)
amrf_results_accession_cruzado %>% filter(is.na(Contig.id)) %>% pull(ENA_sample_accession)%>% unique()
# 3617 (no tienen genes de resistencia)
all_summary_def_qc_clean_accession %>% pull(ENA_sample_accession)%>% unique()
# 22637
head(amrf_results_accession_cruzado)

write.csv(amrf_results_accession_cruzado, "C:/Users/rferragud/Documents/Projecte/Tables/amrf_results_accession_cruzado.completo.csv", row.names = F)

## 4.2 Cruzar resfinder con sylph ####

resf_results_accession_cruzado.temp <- 
  all_summary_def_qc_clean_accession %>% 
  filter(!is.na(ENA_sample_accession)) %>%
  left_join(resfinder_results.complete, by = c("ENA_sample_accession" = "Biosample"))

resf_results_accession_cruzado <- 
  resf_results_accession_cruzado.temp %>% 
  unite(genetic_background, genetic_background, Mutation, sep = ";", na.rm = TRUE) %>%
  rename(Subclass = Antibiotic) %>%
  mutate(Subclass = str_replace_all(
    Subclass,
    "(?<=\\w)\\+(?=\\w)",  # un ‘+’ precedido y seguido de un carácter de palabra
    "-"                    # lo cambiamos por un guión
  )
  ) %>%
  mutate(across(
    where(is.character),
    ~ na_if(., "")
  )) 
  

# comprobación de que las muestras que se han unido y las que no son el total de las all_summary_def_qc_clean_accession 
resf_results_accession_cruzado %>% pull(ENA_sample_accession) %>% unique()
# 22638 

sobrantes_resf <- 
  anti_join(resfinder_results.complete, all_summary_def_qc_clean_accession, by= c("Biosample" = "ENA_sample_accession")) %>% 
  distinct(Biosample)

sam_sylph <- all_summary_def_qc_remove %>% distinct(ENA_sample_accession) %>% pull()
run_sylph <- all_summary_def_qc_remove %>% distinct(ENA_run_accession) %>% pull()
assembly_sylph <- all_summary_def_qc_remove %>% distinct(assembly) %>% pull()

sobrantes_resf %>% filter(Biosample %in% sam_sylph)
sobrantes %>% filter(Name %in% run_sylph)
sobrantes %>% filter(!Name %in% assembly_sylph & !Name %in% sam_sylph & !Name %in% run_sylph)

write.csv(resf_results_accession_cruzado, "C:/Users/rferragud/Documents/Projecte/Tables/resf_results_accession_cruzado.csv", row.names =F)


# 5. Eliminar de las tablas de los fenotipos las muestras que no pasan el QC ####

# 5.1 EUCAST

EUCAST_sir <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/EUCAST_sir.csv", check.names = FALSE)

# primero obtener las ENA_sample_accession del EUCAST_sir con el all_summary_def_qc_clean_accession 

EUCAST_sir %>% filter(is.na(ENA_sample_accession) & is.na(ENA_run_accession) & !is.na(ENA_experiment_accession))

EUCAST_sir %>% filter(is.na(ENA_sample_accession) & !is.na(ENA_run_accession))

EUCAST_sir %>% filter(!is.na(ENA_sample_accession) & is.na(ENA_run_accession))

# las muestras que tienen ENA_experiment_accession van a tener tb sample o run accession así que se van a poder cruzar solo con sample y run 

# Obtener solo una columna que tenga una accession de cada muestra
EUCAST_coaslesce <- EUCAST_sir %>% mutate(ENA_accessions = coalesce(ENA_sample_accession, ENA_run_accession)) %>%
  select(-ENA_sample_accession, -ENA_run_accession)

EUCAST_coaslesce %>% distinct(ENA_accessions)

# cruzamos EUCAST_coalesce con las accessions: ENA_sample_accession 
eucast_results_accession_cruzado <- 
  all_summary_def_qc_clean_accession %>% 
  filter(!is.na(ENA_sample_accession)) %>%
  inner_join(EUCAST_coaslesce, by = c("ENA_sample_accession" = "ENA_accessions"))

# samples que se han cruzado
sam_eucast <- eucast_results_accession_cruzado %>% pull(ENA_sample_accession) %>% unique()

# cruzamos EUCAST_coalesce con las accessions:ENA_run_accession 
eucast_results_accession_cruzado1 <- 
  all_summary_def_qc_clean_accession %>%
  filter(!ENA_sample_accession %in% sam_eucast) %>%
  inner_join(EUCAST_coaslesce, by = c("ENA_run_accession" = "ENA_accessions"))

# run accession que se han cruzado
run_eucast <- eucast_results_accession_cruzado1 %>% pull(ENA_run_accession) %>% unique()

# muestras del EUCAST_sir que no se han cruzado -> no tienen genes de resitencia según resfinder
sobrantes_euc <- EUCAST_coaslesce %>% 
  filter(!ENA_accessions %in% sam_eucast) %>%
  filter(!ENA_accessions %in% run_eucast)

# df con todas las muestras que se han cruzado, solo tiene ya muestras que pasan el QC
EUCAST_sir_clean <- bind_rows(eucast_results_accession_cruzado, eucast_results_accession_cruzado1) %>%
  select(-assembly) %>%
  mutate(
    antibiotic = str_replace(
      antibiotic,
      "^amoxicillin\\.clavulanic\\.acid$",
      "amoxicillin-clavulanic acid"
    )
  )

sobrantes_acc_euc <- sobrantes_euc %>% distinct(ENA_accessions) %>% pull()

# comprobación de que las muestras que se han unido y las que no son el total de las all_summary_def_qc_clean_accession 
EUCAST_sir_clean %>% pull(ENA_sample_accession) %>% unique()
# 12458
sobrantes_acc_euc %>% unique()
# 2706 (se han eliminado por no pasar los criterios del QC)
EUCAST_coaslesce %>% pull(ENA_accessions)%>% unique()
# 15164

# comprobar que las muestras de sobrantes estan en el df de las muestras que no pasan el QC
remove_sam <- all_summary_def_qc_remove %>% pull(ENA_sample_accession)%>% unique()
remove_run <- all_summary_def_qc_remove %>% pull(ENA_run_accession)%>% unique()

sobrantes_euc %>% filter(ENA_accessions %in% remove_sam) %>% pull(ENA_accessions) %>% unique()
sobrantes_euc %>% filter(!ENA_accessions %in% remove_sam & ENA_accessions %in% remove_run) %>% 
  pull(ENA_accessions) %>% unique()


# ECOFF

ECOFF_sir <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/ECOFF_sir.csv", check.names = FALSE)

# Obtener solo una columna que tenga una accession de cada muestra
ECOFF_coaslesce <- ECOFF_sir %>% mutate(ENA_accessions = coalesce(ENA_sample_accession, ENA_run_accession)) %>%
  select(-ENA_sample_accession, -ENA_run_accession)

# cruzamos ECOFF_coalesce con las accessions: ENA_sample_accession 
ecoff_results_accession_cruzado <- 
  all_summary_def_qc_clean_accession %>% 
  filter(!is.na(ENA_sample_accession)) %>%
  inner_join(ECOFF_coaslesce, by = c("ENA_sample_accession" = "ENA_accessions"))

# samples que se han cruzado
sam_ecoff <- ecoff_results_accession_cruzado %>% pull(ENA_sample_accession) %>% unique()

# cruzamos ECOFF_coalesce con las accessions:ENA_run_accession 
ecoff_results_accession_cruzado1 <- 
  all_summary_def_qc_clean_accession %>%
  filter(!ENA_sample_accession %in% sam_ecoff) %>%
  inner_join(ECOFF_coaslesce, by = c("ENA_run_accession" = "ENA_accessions"))

# run accession que se han cruzado
run_ecoff<- ecoff_results_accession_cruzado1 %>% pull(ENA_run_accession) %>% unique()

# df con todas las muestras que se han cruzado, solo tiene ya muestras que pasan el QC
ECOFF_sir_clean <- bind_rows(ecoff_results_accession_cruzado, ecoff_results_accession_cruzado1) %>%
  select(-assembly)%>%
  mutate(
    antibiotic = str_replace(
      antibiotic,
      "^amoxicillin\\.clavulanic\\.acid$",
      "amoxicillin-clavulanic acid"
    )
  )

# comprobación de que las que no se han unido son las que no entran en el sylph 
sobrantes_ecoff <- ECOFF_coaslesce %>% 
  filter(!ENA_accessions %in% sam_ecoff) %>%
  filter(!ENA_accessions %in% run_ecoff)

sobrantes_acc_ecoff <- sobrantes_ecoff %>% distinct(ENA_accessions) %>% pull()

# comprobación de que las muestras que se han unido y las que no son el total de las all_summary_def_qc_clean_accession 
ECOFF_sir_clean %>% pull(ENA_sample_accession) %>% unique()
# 12253
sobrantes_acc_ecoff %>% unique()
# 2701 (se han eliminado por no pasar los criterios del QC)
ECOFF_coaslesce %>% pull(ENA_accessions)%>% unique()
# 14954

sobrantes_ecoff %>% filter(ENA_accessions %in% remove_sam) %>% pull(ENA_accessions) %>% unique()
sobrantes_ecoff %>% filter(!ENA_accessions %in% remove_sam & ENA_accessions %in% remove_run) %>% 
  pull(ENA_accessions) %>% unique()

# CLSI 

CLSI_sir <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/CLSI_sir.csv", check.names = FALSE)

# Obtener solo una columna que tenga una accession de cada muestra
CLSI_coaslesce <- CLSI_sir %>% mutate(ENA_accessions = coalesce(ENA_sample_accession, ENA_run_accession)) %>%
  select(-ENA_sample_accession, -ENA_run_accession)

# cruzamos CLSI_coalesce con las accessions: ENA_sample_accession 
clsi_results_accession_cruzado <- 
  all_summary_def_qc_clean_accession %>% 
  filter(!is.na(ENA_sample_accession)) %>%
  inner_join(CLSI_coaslesce, by = c("ENA_sample_accession" = "ENA_accessions"))

# samples que se han cruzado
sam_clsi <- clsi_results_accession_cruzado %>% pull(ENA_sample_accession) %>% unique()

# cruzamos CLSI_coalesce con las accessions:ENA_run_accession 
clsi_results_accession_cruzado1 <- 
  all_summary_def_qc_clean_accession %>%
  filter(!ENA_sample_accession %in% sam_clsi) %>%
  inner_join(CLSI_coaslesce, by = c("ENA_run_accession" = "ENA_accessions"))

# run accession que se han cruzado
run_clsi <- clsi_results_accession_cruzado1 %>% pull(ENA_run_accession) %>% unique()

# df con todas las muestras que se han cruzado, solo tiene ya muestras que pasan el QC
CLSI_sir_clean <- bind_rows(clsi_results_accession_cruzado, clsi_results_accession_cruzado1) %>%
  select(-assembly) %>% 
  mutate(
    antibiotic = str_replace(
      antibiotic,
      "amoxicillin.clavulanic.acid",
      "amoxicillin-clavulanic acid"
    )
  )

# comprobación de que las que no se han unido son las que no entran en el sylph 
sobrantes_clsi <- CLSI_coaslesce %>% 
  filter(!ENA_accessions %in% sam_clsi) %>%
  filter(!ENA_accessions %in% run_clsi)

sobrantes_acc_clsi <- sobrantes_clsi %>% distinct(ENA_accessions) %>% pull()

# comprobación de que las muestras que se han unido y las que no son el total de las all_summary_def_qc_clean_accession 
CLSI_sir_clean %>% pull(ENA_sample_accession) %>% unique()
# 12458
sobrantes_acc_clsi %>% unique()
# 2706 (se han eliminado por no pasar los criterios del QC)
CLSI_coaslesce %>% pull(ENA_accessions)%>% unique()
# 15164

sobrantes_clsi %>% filter(ENA_accessions %in% remove_sam) %>% pull(ENA_accessions) %>% unique()
sobrantes_clsi %>% filter(!ENA_accessions %in% remove_sam & ENA_accessions %in% remove_run) %>% 
  pull(ENA_accessions) %>% unique()

# renombrar con los atb estandarizados
EUCAST_sir_clean <- EUCAST_sir_clean %>% mutate(antibiotic = str_replace_all(antibiotic, "\\.", "-"))

ECOFF_sir_clean <- ECOFF_sir_clean %>% 
  mutate(antibiotic = str_replace_all(antibiotic, "\\.", "-")) %>%
  mutate(antibiotic = str_replace(antibiotic, "nalidixic-acid", "nalidixic acid")) %>%
  mutate(antibiotic = str_replace(antibiotic, "cefepime_taniborbactam", "cefepime-taniborbactam"))

CLSI_sir_clean <- CLSI_sir_clean %>% 
  mutate(antibiotic = str_replace_all(antibiotic, "\\.", "-")) %>%
  mutate(antibiotic = str_replace(antibiotic, "nalidixic-acid", "nalidixic acid"))
# contar fenotipos para elegir antibioticos 

count_eucast.temp.clean <- EUCAST_sir_clean %>% filter(!is.na(EUCAST_2025)) %>% count(antibiotic, EUCAST_2025, sort = TRUE) %>% rename(EUCAST = n)

count_eucast.clean <- count_eucast.temp.clean %>%
  group_by(antibiotic) %>%
  summarise(
    R_eucast = sum(EUCAST[EUCAST_2025 == "R"], na.rm = TRUE),
    S_eucast = sum(EUCAST[EUCAST_2025 == "S"], na.rm = TRUE),
    I_eucast = sum(EUCAST[EUCAST_2025 == "I"], na.rm = TRUE),
    total_eucast = sum(EUCAST, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  select(antibiotic, R_eucast, S_eucast, I_eucast, total_eucast)


count_ecoff.temp.clean <- ECOFF_sir_clean %>% filter(!is.na(ECOFF_2025)) %>% count(antibiotic, ECOFF_2025, sort = TRUE) %>% rename(ECOFF = n)

count_ecoff.clean <- count_ecoff.temp.clean %>%
  group_by(antibiotic) %>%
  summarise(
    R_ecoff = sum(ECOFF[ECOFF_2025 == "R"], na.rm = TRUE),
    S_ecoff = sum(ECOFF[ECOFF_2025 == "S"], na.rm = TRUE),
    I_ecoff = sum(ECOFF[ECOFF_2025 == "I"], na.rm = TRUE),
    total_ecoff = sum(ECOFF, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  select(antibiotic, R_ecoff, S_ecoff, I_ecoff, total_ecoff)


count_clsi.temp.clean <- CLSI_sir_clean %>% filter(!is.na(CLSI_2025)) %>% count(antibiotic, CLSI_2025, sort = TRUE) %>% rename(CLSI = n)

count_clsi.clean <- count_clsi.temp.clean %>%
  group_by(antibiotic) %>%
  summarise(
    R_clsi = sum(CLSI[CLSI_2025 == "R"], na.rm = TRUE),
    S_clsi = sum(CLSI[CLSI_2025 == "S"], na.rm = TRUE),
    I_clsi = sum(CLSI[CLSI_2025 == "I"], na.rm = TRUE),
    total_clsi = sum(CLSI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  select(antibiotic, R_clsi, S_clsi, I_clsi, total_clsi)

count_comparison.clean <- full_join(count_eucast.clean, count_ecoff.clean, by = "antibiotic") %>%
  full_join(count_clsi.clean, by = "antibiotic") %>% 
  select(antibiotic, R_eucast, S_eucast, I_eucast, total_eucast, R_ecoff, S_ecoff, I_ecoff, total_ecoff, R_clsi, S_clsi, I_clsi, total_clsi)

count_comparison.completa.clean <- count_comparison.clean %>%
  # 1) Llevamos sólo las columnas R_*, S_*, I_* a filas,
  #    dejando antibiotic, study_count y studies intactas
  pivot_longer(
    cols = matches("^(R|S|I)_.+"),      # todas las columnas que empiezan R_, S_ o I_
    names_to  = c("measure", "testing_standard"), # separa en medida (R/S/I) y método (eucast/ecoff/clsi)
    names_sep = "_",
    values_to = "count"
  ) %>%
  # 2) Reconstruimos columnas R, S e I
  pivot_wider(
    id_cols     = c(antibiotic, testing_standard),
    names_from  = measure,              # R, S, I
    values_from = count                 # conteos
  ) %>%
  # 3) (Opcional) ordenar
  mutate(
    total = R + S + I,
    testing_standard = factor(testing_standard, levels = c("eucast", "ecoff", "clsi"))
  ) %>%
  select(antibiotic, testing_standard, R, S, I, total)%>%
  arrange(desc(total))

write.csv(count_comparison.clean, "C:/Users/rferragud/Documents/Projecte/Tables/count_comparison.clean.csv", row.names = FALSE)
write.csv(count_comparison.completa.clean, "C:/Users/rferragud/Documents/Projecte/Tables/count_comparison.completa.clean.csv", row.names = FALSE)



# Criterios de inclusión de los antibióticos 

# detección de mic en más de 500 muestras 
# detección de mic en más de 100 muestras para cada tipo de fenotipo (R/S)

count_comparison.completa.filtro1 <- 
  count_comparison.completa.clean %>% 
  mutate(
    filtro = if_else(
      total < 500,
      "< 500 isolates with MIC",
      NA_character_
    )
  )


count_comparison.completa.filtro2 <- count_comparison.completa.filtro1 %>%
  mutate(
    filtro = if_else(
      is.na(filtro) & (R < 100 | S < 100 | is.na(R) | is.na(S)), 
      "< 100 R/S isolates", 
      filtro
    )
  )


count_comparison.completa.filtro3 <- count_comparison.completa.filtro2 %>%
  mutate(
    filtro = if_else(
      is.na(filtro), 
      "included", 
      filtro
    )
  )

write.csv(count_comparison.completa.filtro3, "C:/Users/rferragud/Documents/Projecte/Tables/atb_included.csv")

atb_ecoff <- count_comparison.completa.filtro3 %>% 
  filter(testing_standard == "ecoff" & filtro == "included")

atb_ecoff_names <- atb_ecoff %>% arrange(antibiotic) %>% pull(antibiotic)

atb_eucast <- count_comparison.completa.filtro3 %>% 
  filter(testing_standard == "eucast" & filtro == "included")

atb_eucast_names <- atb_eucast %>% arrange(antibiotic) %>% pull(antibiotic)

atb_clsi <- count_comparison.completa.filtro3 %>% 
  filter(testing_standard == "clsi" & filtro == "included")

atb_clsi_names <- atb_clsi %>% arrange(antibiotic) %>% pull(antibiotic)

atb_total <- count_comparison.completa.filtro3 %>% 
  filter(filtro == "included") %>% arrange(antibiotic) %>% distinct(antibiotic)


#-------------------------------------------------------------------------------

amrf_results_accession_cruzado.completo %>% distinct(Subclass)

resf_results_accession_cruzado %>% distinct(Subclass)

# TODO para resfinder cambiar pq los resistentes son los que tinen REsistant, no los del genetic_backgound
# función para todos los antibióticos 
analizar_antibiotico_resf <- function(antibiotic, program_class, sir_df, program) {
  
  # Filtrar EUCAST para cada antibiótico
  pheno_filtered <- sir_df %>%
    filter(antibiotic == !!antibiotic)
  # Construimos un patrón para que detecte antibiotic y resfinder_class (si no es NA)
  pattern <- if (!is.na(program_class)) {
    paste0(program_class, "|", antibiotic)
  } else {
    antibiotic
  }
  
  # 3) Filtrar ResFinder usando el patrón adecuado
  program_filtered <- program %>%
    filter(str_detect(
      Subclass,
      regex(pattern, ignore_case = TRUE)
    ))
  
  # Join usando muestras únicas
  cruzado <- pheno_filtered %>%
    left_join(
      program_filtered %>%
        select(ENA_sample_accession, WGS_predicted_phenotype, genetic_background),
      by = "ENA_sample_accession") %>%
    mutate(
      predition_type = case_when(
        sir == "R" & WGS_predicted_phenotype == "Resistant" ~ "TP",
        sir == "S" & WGS_predicted_phenotype == "Resistant"  ~ "FP",
        sir == "S" & WGS_predicted_phenotype == "No resistance" ~ "TN",
        sir == "R" & WGS_predicted_phenotype == "No resistance" ~ "FN"
      ),
      antibiotic = antibiotic  # asegura que quede como columna en cada fila
    ) 
  
  return(cruzado)
}



# Funcion para Resfinder 
corresp_atb_resfinder.ecoli <- read_ods("C:/Users/rferragud/Documents/Projecte/corresp_antibiotic_resf_amrf.new.ods", sheet = 2)

# eucast
corresp_atb_euc_resf.ecoli <- corresp_atb_resfinder.ecoli %>% filter(antibiotic %in% atb_eucast_names)
# ecoff
corresp_atb_ecoff_resf.ecoli <- corresp_atb_resfinder.ecoli %>% filter(antibiotic %in% atb_ecoff_names)
# clsi
corresp_atb_clsi_resf.ecoli <- corresp_atb_resfinder.ecoli %>% filter(antibiotic %in% atb_clsi_names)

EUCAST_sir_clean2 <- EUCAST_sir_clean %>% rename(sir = "EUCAST_2025") %>%
  mutate(sir = if_else(as.character(sir) == "I", "R", as.character(sir))) %>% 
  filter(!sir == "NI")

CLSI_sir_clean2 <- CLSI_sir_clean %>% rename(sir = "CLSI_2025") %>%
  mutate(sir = if_else(as.character(sir) == "I", "R", as.character(sir))) %>% 
  filter(!sir == "NI")

ECOFF_sir_clean2 <- ECOFF_sir_clean %>% rename(sir = "ECOFF_2025") %>%
  mutate(sir = if_else(as.character(sir) == "I", "R", as.character(sir))) %>% 
  filter(!sir == "NI")


# resultados resfinder
# a) eucast
resultados_detallados_resfinder_eucast.ecoli <- corresp_atb_euc_resf.ecoli %>%
  mutate(
    resultados = map2(
      antibiotic,
      program_class,
      analizar_antibiotico_resf,
      sir_df = EUCAST_sir_clean2, 
      program = resf_results_accession_cruzado
    )
  ) %>%
  unnest(resultados, names_sep = "_")


# b) ecoff
resultados_detallados_resfinder_ecoff.ecoli <- corresp_atb_ecoff_resf.ecoli %>%
  mutate(
    resultados = map2(
      antibiotic,
      program_class,
      analizar_antibiotico_resf,
      sir_df = ECOFF_sir_clean2, 
      program = resf_results_accession_cruzado
    )
  ) %>%
  unnest(resultados, names_sep = "_")

# c) clsi
resultados_detallados_resfinder_clsi.ecoli <- corresp_atb_clsi_resf.ecoli %>%
  mutate(
    resultados = map2(
      antibiotic,
      program_class,
      analizar_antibiotico_resf,
      sir_df = CLSI_sir_clean2, 
      program = resf_results_accession_cruzado
    )
  ) %>%
  unnest(resultados, names_sep = "_")

# --------------------------------------------------------------------------

# Resumen por antibiótico
resumen_por_antibiotico_resfinder_euc.ecoli <- resultados_detallados_resfinder_eucast.ecoli %>%
  count(antibiotic, resultados_predition_type) %>% 
  pivot_wider(names_from = antibiotic, values_from = n)























diagnostics_wide <-
  diagnosticos_por_ab %>% 
  group_by(resultados_antibiotic) %>% 
  pivot_wider(id_cols=c(resultados_antibiotic),  names_from = statistic, values_from = est) %>%
  select(resultados_antibiotic, se, sp, diag.ac, pv.pos, pv.neg, TP, TN, FP, FN) %>% 
  mutate(total_score = se + sp + diag.ac) %>%
  arrange(desc(total_score))


diagnostics_porcentaje <- diagnosticos_por_ab %>%
  mutate(
    est = case_when(
      statistic %in% c("ap", "tp", "se", "sp", "pv.pos", "pv.neg") ~
        sprintf(
          "%.1f (%.1f-%.1f)",
          est  * 100,    # convierte a porcentaje
          lower * 100,
          upper * 100
        ),
      TRUE ~ as.character(est)  # deja tal cual las demás estadísticas
    ),
    # Si quieres, puedes poner NA en lower/upper para esas filas:
    lower = if_else(statistic %in% c("ap", "tp", "se", "sp", "pv.pos", "pv.neg"),
                    NA_real_, lower),
    upper = if_else(statistic %in% c("ap", "tp", "se", "sp", "pv.pos", "pv.neg"),
                    NA_real_, upper)
  )


# 2. Filtramos todo el df para quitarles esos antibióticos
diagnostics_wide_porcentaje <-
  diagnostics_porcentaje %>% 
  group_by(resultados_antibiotic) %>% 
  pivot_wider(id_cols=c(resultados_antibiotic),  names_from = statistic, values_from = est) %>%
  select(resultados_antibiotic, se, sp, diag.ac, pv.pos, pv.neg, TP, TN, FP, FN) %>%
  rename(sensitivity = se, specificity = sp, PPV = pv.pos, NPV = pv.neg)


plot <- diagnostics_wide %>%
  select(resultados_antibiotic, se, sp) %>%
  pivot_longer(
    cols      = c(se, sp),
    names_to  = "métrica",
    values_to = "valor"
  ) %>%
  ggplot(aes(
    x    = resultados_antibiotic,
    y    = valor,
    fill = métrica
  )) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = "Antibiótico",
    y     = "Porcentaje",
    fill  = "Métrica",
    title = "Sensibilidad vs Especificidad por antibiótico"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )























# Amrfplus

analizar_antibiotico <- function(antibiotic, program_class, sir_df, program) {
  
  # Filtrar EUCAST para cada antibiótico
  pheno_filtered <- sir_df %>%
    filter(antibiotic == !!antibiotic)
  # Construimos un patrón para que detecte antibiotic y resfinder_class (si no es NA)
  pattern <- if (!is.na(program_class)) {
    paste0(program_class, "|", antibiotic)
  } else {
    antibiotic
  }
  
  # 3) Filtrar ResFinder usando el patrón adecuado
  program_filtered <- program %>%
    filter(str_detect(
      Subclass,
      regex(pattern, ignore_case = TRUE)
    ))
  
  # Join usando muestras únicas
  cruzado <- pheno_filtered %>%
    left_join(
      program_filtered %>%
        select(ENA_sample_accession, genetic_background),
      by = "ENA_sample_accession") %>%
    mutate(
      predition_type = case_when(
        sir == "R" & !is.na(genetic_background) ~ "TP",
        sir == "S" & is.na(genetic_background)  ~ "TN",
        sir == "I"                           ~ "unknown",
        sir == "S" & !is.na(genetic_background) ~ "FP",
        sir == "R" & is.na(genetic_background)  ~ "FN"
      ),
      antibiotic = antibiotic  # asegura que quede como columna en cada fila
    ) 
  
  return(cruzado)
}










corresp_atb_amrf <- read_ods("C:/Users/rferragud/Documents/Projecte/corresp_antibiotic_resf_amrf.ods", sheet = 2)

# eucast
corresp_atb_euc_amrf <- corresp_atb_amrfplus %>% filter(antibiotic %in% atb_eucast_names)

# ecoff
corresp_atb_ecoff_amrf <- corresp_atb_amrfplus %>% filter(antibiotic %in% atb_ecoff_names)

# clsi
corresp_atb_ecoff_amrf <- corresp_atb_amrfplus %>% filter(antibiotic %in% atb_clsi_names)












