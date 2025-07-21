
library(AMR)
library(tidyverse)
library(data.table)
library(dplyr)
library(stringr)

df_cruzado_completo <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/df_cruzado_completo_articles.csv", check.names = FALSE)

# 1. Saber cuantos atb tienen info del MIC (en la colección de articulos) #### 

MIC_columns <- 
  df_cruzado_completo %>% 
  select(ends_with("measurement"))%>%
  select(where(~ any(!is.na(.))))

dim(df_cruzado_completo)
dim(MIC_columns)

# Extraer nombres de antibióticos sin "_measurement"
antibiotics_mic <- gsub("_measurement", "", colnames(MIC_columns))

print(paste("Se encontraron valores MIC para los siguientes antibióticos" , antibiotics_mic))
cat("Se encontraron valores MIC para los siguientes antibióticos:\n")
cat(paste(antibiotics_mic, collapse = ", "))

# Contar valores no NA por antibiótico
mic_count <- colSums(!is.na(MIC_columns))

# Crear dataframe con antibiótico y conteo
mic_count.df <- data.frame(
  Antibiotic = antibiotics_mic,
  Muestras = mic_count
)

df_cruzado_completo <- df_cruzado_completo %>%
  mutate(across(everything(), as.character))

# Convertir all_summaries a formato largo para obtener internal_study_name
long_data <- df_cruzado_completo %>%
  select(internal_study_name, ends_with("measurement")) %>%
  pivot_longer(cols = -internal_study_name, 
               names_to = "Antibiotic", 
               values_to = "measurement") %>%
  filter(!is.na(measurement)) %>%  # Filtrar valores no NA
  mutate(Antibiotic = gsub("_measurement", "", Antibiotic))

# Obtener la lista de estudios por antibiótico
study_count.df <- long_data %>%
  group_by(Antibiotic) %>%
  summarize( #trata a los antibióticos cada uno por separado 
    Study_Count = n_distinct(internal_study_name),  # Número de estudios distintos
    Studies = paste(unique(internal_study_name), collapse = ", "),  # Lista de estudios
    .groups = "drop"
  )

# Unir los estudios por antibiótico y las muestras por antibiótico
study_mic_count.df <- mic_count.df %>%
  left_join(study_count.df, by = "Antibiotic")

# Ver el resultado final
print(study_mic_count.df)

# 2. crear df con las muestras que tienen información de MIC en los artículos ####
#extraer nombres de las columnas con el MIC
MIC_info_columns <- 
  df_cruzado_completo %>% 
  select(matches("(measurement$|measurement_sign$|measurement_units$|laboratory_typing_method$|laboratory_typing_platform$|testing_standard$)"))%>%
  select(starts_with(antibiotics_mic))

dim(MIC_info_columns)

#nombres de las columnas
MIC_info_columns_names <- MIC_info_columns %>% colnames()

#nombres de las columnas de mic, unidades y signo
MIC_columns_names <-
  df_cruzado_completo %>% 
  select(matches("(measurement$|measurement_sign$|measurement_units$)"), ) %>% colnames()

sort(antibiotics_mic)

MIC_columns_names

# Seleccionar columnas deseadas en el df
df_MIC <- df_cruzado_completo %>%
  select(internal_id, ENA_experiment_accession, ENA_run_accession, ENA_sample_accession, internal_study_name, ENA_scientific_name, all_of(MIC_info_columns_names)) %>%
  
  # Filtrar filas con al menos un valor no NA en columnas measurement
  filter(
    if_any(ends_with("_measurement"), ~ !is.na(.))
  )

tail(df_MIC)

#signs <- unique( na.omit(
#  unlist( df_MIC[ grep("measurement_sign$", names(df_MIC)) ] )
#) )
#
#ceftiofur_measurement_units
#
#sort(MIC_info_columns_names)
#
#NAS <- df_MIC %>%
#  rowwise() %>%
#  mutate(contains_NA_string = any(c_across(everything()) == "NA")) %>%
#  filter(contains_NA_string)
#
#
#df_MIC%>% str_detect("NA")
#
#
## convertir automaticamente todos los string NA a NA 
#
#platform <- unique( na.omit(
#  unlist( df_MIC[ grep("laboratory_typing_platform$", names(df_MIC)) ] )
#) )   
#
#methods <- unique( na.omit(
#  unlist( df_MIC[ grep("laboratory_typing_method$", names(df_MIC)) ] )
#) )   
#
#measurement <- unique( na.omit(
#  unlist( df_MIC[ grep("measurement$", names(df_MIC)) ] )
#) )
#
#
#df_MIC %>%
#  select(all_of(MIC_info_columns_names)) %>%
#  map(~ .x[!grepl("^\\s*(<=|>=|<|>)?\\s*[0-9.]+\\s*$", as.character(.x))] %>% unique()) %>%
#  discard(~ length(.x) == 0)


# quitamos los valores no numericos como "ND" o "NA" y lo convertimos a NA
df_MIC2 <- df_MIC %>%
  mutate(across(
    everything(),
    ~ replace(.x, .x %in% c("ND", "NA", "N/A"), NA)
  ))

# 3. Pasamos el df a formato long, con un atb por fila ####

df_long <- df_MIC2 %>%
  pivot_longer(
    cols     = all_of(MIC_info_columns_names),
    # primero grupo = antibiotic; segundo grupo = el sufijo que coincide con .value
    names_to = c("antibiotic", ".value"),
    names_pattern = "(.+)_(laboratory_typing_method|laboratory_typing_platform|measurement|measurement_sign|measurement_units|testing_standard)$"
  ) %>%
  # simplificamos nombres
  rename(
    sign  = measurement_sign,
    units = measurement_units
  ) %>%
  # dejamos solo y en orden las columnas finales
  select(
    internal_id,
    ENA_experiment_accession,
    ENA_run_accession,
    ENA_sample_accession,
    internal_study_name,
    ENA_scientific_name,
    antibiotic,
    sign, 
    measurement,
    units,
    laboratory_typing_method,
    laboratory_typing_platform, 
    testing_standard
  )

# Vemos un fragmento
print(head(df_long))

antibiotics_mic%>%sort()

length(unique(df_long$antibiotic)) 

# ¿Cuáles son?
sort(unique(df_long$antibiotic))
#hay más pq hay atb que tienen info de signo, unidades o testing standard

#conservamos solo los atb que tienen medida MIC
df_long_clean <- df_long %>%
  filter(!is.na(measurement))%>%
  mutate(
    across(
      where(is.character),  # solo columnas de texto
      ~ str_squish(.)  
    )
  )

unique(df_long_clean$antibiotic)

df_long_clean %>% distinct(ENA_sample_accession)

#------------------------------------------------------------------ #

# 4. Unión de las muestras del NCBI ####

muestras_NCBI.nuevas.coli <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/muestras_NCBI.nuevas.coli.csv", check.names = FALSE)

#cambiar los nombres de las columnas a nuestro formato 

#nombres de las columnas en el NCBI
muestras_NCBI.nuevas.coli %>% colnames()

#nombres de las columnas estandarizados en nuestros datos 
df_long_clean %>% colnames()

#correpondencia de los nombres de las columnas 
columns_st <- c(ENA_sample_accession = "#BioSample", 
                ENA_scientific_name = "Scientific name",
                bioproject_accession = "BioProject",
                antibiotic = "Antibiotic",
                measurement_disk = "Disk diffusion (mm)", 
                measurement_broth = "MIC (mg/L)",
                sign = "Measurement sign", 
                vendor = "Vendor", 
                testing_standard = "Testing standard",
                isolate_id_publication = "Isolate",  
                isolation_source = "Isolation type", 
                isolation_source_2 = "Isolation source", 
                location = "Location", 
                laboratory_typing_platform = "Laboratory typing platform", 
                laboratory_typing_method_version_or_reagent = "Laboratory typing method version or reagent",
                create_date = "Create date", 
                resitance_phenotype = "Resistance phenotype")


#df del NCBI con nuestros nombres estandarizados 
muestras_NCBI.nuevas.coli.st.temp1 <- 
  muestras_NCBI.nuevas.coli %>% 
  select(-`Organism group`)%>%
  rename(all_of(columns_st))

muestras_NCBI.nuevas.coli.st.temp2 <- 
  muestras_NCBI.nuevas.coli.st.temp1 %>% 
  unite(isolation_source, isolation_source, isolation_source_2, sep = "; ", remove = TRUE, na.rm = TRUE)

# eliminar muestras que tienen más de un resultado para el mismo atb 

muestras_NCBI.nuevas.coli.st.temp2 %>% distinct(ENA_sample_accession) # 9313

samples_repeated_ncbi <- muestras_NCBI.nuevas.coli.st.temp2 %>%  # 8 repetidas (las eliminamos)
  group_by(ENA_sample_accession, antibiotic) %>%
  filter(n() > 1) %>%
  ungroup() %>% 
  distinct(ENA_sample_accession) %>%
  pull()

muestras_NCBI.nuevas.coli.st <- muestras_NCBI.nuevas.coli.st.temp2 %>% 
  filter(!ENA_sample_accession %in% samples_repeated_ncbi)

muestras_NCBI.nuevas.coli.st %>% distinct(ENA_sample_accession)

write.csv(muestras_NCBI.nuevas.coli.st, "C:/Users/rferragud/Documents/Projecte/Tables/muestras_NCBI.nuevas.coli.st.csv", row.names = F)

#df del NCBI con mic i su info 
muestras_NCBI.nuevas.coli.st_mic.temp <- 
  muestras_NCBI.nuevas.coli.st %>% 
  select(ENA_sample_accession, 
         ENA_scientific_name, 
         antibiotic, 
         measurement_broth, 
         sign, 
         testing_standard, 
         laboratory_typing_platform, 
        )%>%
  filter(!is.na(measurement_broth))

muestras_NCBI.nuevas.coli.st_mic<- muestras_NCBI.nuevas.coli.st_mic.temp %>%
  mutate(
    laboratory_typing_method = "broth microdilution",
    units = "mg/L",
    internal_study_name = "NCBI"
  ) %>% 
  rename(measurement = "measurement_broth")%>%
  mutate(across(everything(), as.character))

muestras_NCBI.nuevas.coli.st_mic %>% distinct(ENA_sample_accession)

muestras_NCBI.nuevas.coli.st_disk.temp <- 
  muestras_NCBI.nuevas.coli.st %>% 
  select(ENA_sample_accession, 
         ENA_scientific_name, 
         antibiotic, 
         measurement_disk, 
         sign, 
         testing_standard, 
         laboratory_typing_platform, 
  )%>%
  filter(!is.na(measurement_disk))

muestras_NCBI.nuevas.coli.st_disk<- muestras_NCBI.nuevas.coli.st_disk.temp %>%
  mutate(
    laboratory_typing_method = "disk diffusion",
    units = "mm",
    internal_study_name = "NCBI"
  ) %>% 
  rename(measurement = "measurement_disk")%>%
  mutate(across(everything(), as.character))

muestras_NCBI.nuevas.coli.st_disk%>%pull(antibiotic)%>%unique()

muestras_NCBI.nuevas.coli.st.sep <- 
  bind_rows(muestras_NCBI.nuevas.coli.st_mic, muestras_NCBI.nuevas.coli.st_disk)

mic_articles_ncbi <- bind_rows(df_long_clean, muestras_NCBI.nuevas.coli.st.sep)%>%
  select(internal_id:ENA_run_accession, ENA_sample_accession, internal_study_name:testing_standard)

muestras_NCBI.nuevas.coli.st_mic %>% distinct(ENA_sample_accession)
mic_articles_ncbi %>% colnames()
mic_articles_ncbi %>% distinct(ENA_sample_accession)
df_long_clean %>% distinct(ENA_sample_accession)
mic_articles_ncbi %>% distinct(antibiotic) %>% print(n = Inf)
muestras_NCBI.nuevas.coli.st.sep %>% distinct(ENA_sample_accession)

#------------------------------------------------------------------------------#

# 5. unión de las muestras del PATRIC ####

getwd()

#muestras del PATRIC sin repeticiones con nuestra coleccion (pero sí con el NCBI)
muestras_patric.nuevas <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/muestras_patric.nuevas.csv", check.names = FALSE)

muestras_patric.nuevas %>% colnames()

df_long_clean %>% colnames()

#obtener los metadatos de interés 
patric_metadata2 <- muestras_patric.nuevas %>% 
  select(Antibiotic, Resistant.Phenotype:Testing.Standard.Year, Species, NCBI.Taxon.ID, 
         Culture.Collection, BioProject.Accession,
         Sequencing.Platform, Isolation.Source,
         Collection.Date, Isolation.Country, Geographic.Location, Host.Name, Host.Health, 
         Comments, match_sam.p:match_final_ers, -Measurement )

#correpondencia de las columnas con las que tenemos !(sin ENA_, no sé si deberiamos ponerlo) 
columns_metadata_st_patric <-  c(
  antibiotic = "Antibiotic",
  resistance_phenotype = "Resistant.Phenotype",
  measurement = "Measurement.Value",
  measurement_sign = "Measurement.Sign",
  measurement_units = "Measurement.Unit",
  laboratory_typing_method = "Laboratory.Typing.Method", 
  laboratory_typing_method_version_or_reagent = "Laboratory.Typing.Method.Version",
  laboratory_typing_platform = "Laboratory.Typing.Platform", 
  vendor = "Vendor", 
  testing_standard = "Testing.Standard", 
  testing_standard_2 = "Testing.Standard.Year", 
  tax_id = "NCBI.Taxon.ID", 
  ENA_scientific_name = "Species", 
  culture_collection = "Culture.Collection", 
  bioproject_accession = "BioProject.Accession",
  instrument_model = "Sequencing.Platform", 
  isolation_source = "Isolation.Source",
  collection_date = "Collection.Date", 
  geographic_location_country = "Isolation.Country", 
  geographic_location_region = "Geographic.Location", 
  host = "Host.Name", 
  host_disease = "Host.Health", 
  additional_information = "Comments", 
  ENA_sample_accession = "match_sam.p", 
  ENA_secondary_sample_accession = "match_ers.p", 
  ENA_experiment_accession = "match_erx.p", 
  ENA_run_accession = "match_err.p")

# computational_method, culture_collection no sé si la queremos 

#renombramiento del df
patric_metadata_st2 <- 
  patric_metadata2 %>% 
  rename(all_of(columns_metadata_st_patric)) %>% 
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~na_if(.x, ""))) %>%  
  unite("testing_standard", testing_standard, testing_standard_2, sep = " ", na.rm = TRUE)
  
dim(patric_metadata_st2)

#selección de las muestras nuevas del patric (respecto al NCBI)
patric_ENA_accessions.nuevas.temp <- read.csv( "C:/Users/rferragud/Documents/Projecte/Tables/patric_ENA_accessions.nuevas.temp.csv")

sample_acc_patric <- patric_ENA_accessions.nuevas.temp %>% filter(!is.na(ENA_sample_accession))%>% pull(ENA_sample_accession) %>% unique()

run_acc_patric <- patric_ENA_accessions.nuevas.temp %>% filter(!is.na(ENA_run_accession))%>% pull(ENA_run_accession) %>% unique()

sec_acc_patric <- patric_ENA_accessions.nuevas.temp %>% filter(!is.na(ENA_secondary_sample))%>% pull(ENA_secondary_sample) %>% unique()

  
# df con toda los metadatos de las muestras de PATRIC de interés 
patric_metadata_st.nuevas <- 
  patric_metadata_st2 %>% 
  filter(ENA_sample_accession %in% sample_acc_patric | ENA_run_accession %in% run_acc_patric | ENA_secondary_sample_accession %in% sec_acc_patric) 

write.csv(patric_metadata_st.nuevas, "C:/Users/rferragud/Documents/Projecte/Tables/patric_metadata_st.nuevas.csv", row.names = F)

# seleccion de las columnas de MIC
patric_metadata_st.nuevas_mic.temp <- 
  patric_metadata_st.nuevas %>% 
  select(ENA_sample_accession, 
         ENA_run_accession,
         ENA_scientific_name, 
         antibiotic, 
         measurement, 
         measurement_sign,
         measurement_units,
         testing_standard, 
         laboratory_typing_platform, 
         laboratory_typing_method
  )%>%
  rename(sign = "measurement_sign", units = "measurement_units") %>%
  mutate(internal_study_name = "PATRIC") %>%
  filter(!is.na(measurement)) %>%
  mutate(across(everything(), as.character))%>%
  mutate(antibiotic = str_replace_all(antibiotic, "/", "-"))%>%
  mutate(antibiotic = str_replace_all(antibiotic, "Â", ""))%>%
  mutate(across(everything(), str_trim))



#patric_metadata_st.nuevas_mic %>% filter(!is.na(sign))
#patric_metadata_st.nuevas %>% filter(!is.na(ENA_sample_accession)) %>% distinct(ENA_sample_accession) 
#patric_metadata_st.nuevas %>% filter(!is.na(ENA_run_accession)) %>% filter(is.na(ENA_sample_accession))%>%
#distinct(ENA_run_accession) 
#patric_metadata_st.nuevas_mic %>% distinct(laboratory_typing_method, ENA_sample_accession) %>% count(laboratory_typing_method)
#patric_metadata_st.nuevas_mic %>% filter(str_detect(antibiotic, "-")) %>% distinct(antibiotic)
#patric_metadata_st.nuevas_mic %>% filter(laboratory_typing_method == "MIC") %>% distinct(units)
#
#sample_patric <- patric_metadata_st.nuevas_mic %>% pull(ENA_sample_accession)%>% unique()
#patric_metadata_st.nuevas_mic %>% filter(!ENA_sample_accession %in% sample_patric) %>% pull(ENA_run_accession)


#laboratory_typing_method     n
#1            Agar dilution   610
#2           Broth dilution 11281
#3           Disk diffusion  3167
#4                      MIC 78617
#5                     <NA>    55

#laboratory_typing_method    n
#1            Agar dilution  252
#2           Broth dilution 1972
#3           Disk diffusion  299
#4                      MIC 2895
#5                     <NA>   37


# eliminar las muestras con más de un resultado para el mismo atb en PATRIC

samples_repeated_patric <- patric_metadata_st.nuevas_mic.temp %>%  # 320 repetidas (las eliminamos)
  select(ENA_sample_accession:units) %>%
  distinct()%>%
  group_by(ENA_sample_accession, antibiotic) %>%
  filter(n() > 1) %>%
  ungroup() 

samples_repeated_patric.names <- samples_repeated_patric %>% 
  distinct(ENA_sample_accession) %>%
  pull()

patric_metadata_st.nuevas_mic.temp %>% distinct(ENA_sample_accession)

patric_metadata_st.nuevas_mic <- patric_metadata_st.nuevas_mic.temp %>% 
  filter(!ENA_sample_accession %in% samples_repeated_patric.names)

patric_metadata_st.nuevas_mic %>% distinct(ENA_sample_accession)


# 6. df con toda la colección de muestras con mic ####

mic_articles_ncbi_patric <- bind_rows(mic_articles_ncbi, patric_metadata_st.nuevas_mic)

df_long_clean_sign <- mic_articles_ncbi_patric %>%
  unite(
    col      = "measurement_combined",   # nombre de la nueva columna
    sign, measurement,                   # columnas a unir
    sep      = "",                       # sin espacio: "<" + "0.25" → "<0.25"
    na.rm    = TRUE,                     # si sign es NA, solo deja el valor
    remove   = FALSE 
  )

df_long_clean_sign %>% distinct(ENA_sample_accession)
mic_articles_ncbi_patric %>% distinct(ENA_sample_accession)

broth_microdil <- df_long_clean_sign %>% filter(laboratory_typing_method %in% c("broth microdilution", "Broth dilution", "MIC"))

write.csv(broth_microdil, "C:/Users/rferragud/Documents/Projecte/Tables/broth_microdil.csv", row.names = FALSE)

disk_diffusion <- df_long_clean_sign %>% filter(laboratory_typing_method %in% c("disk diffusion", "Disk diffusion"))

disk <- disk_diffusion %>% distinct(ENA_sample_accession) %>% unique()

broth <- broth_microdil %>% distinct(ENA_sample_accession) %>% unique()

#------------------------------------------------------------#


write.csv(study_count.df, "C:/Users/rferragud/Documents/Projecte/Tables/muestras_mic_articles.csv")


