
#-------------------ELIMINACIÓN DE DUPLICADOS DE LOS ARTÍCULOS -----------------

library(tidyverse)
library(data.table)
library(dplyr)
library(stringr)

df_cruzado <- read.csv("C:/Users/rferragud/Documents/Projecte/RFF_revised/df_cruzado_2025-06-19_11-41-54.csv", check.names = FALSE)

#1.Comprobar si existe información diferente en las ENA_experiment_accession duplicadas####

#muestras con duplicaciones (el df incluye todos los duplicados/triplicados)
df_cruzado_duplicados <- df_cruzado %>%
  group_by(ENA_experiment_accession) %>%
  filter(n() > 1) %>%
  ungroup() %>% 
  select(internal_id:additional_info, match_sam:ENA_bioproject_accession, everything())

#duplicaciones colapsadas sin NA (run accession únicos)  
df_cruzado_colapsado <- df_cruzado_duplicados %>%
  group_by(ENA_experiment_accession) %>%
  summarise(across(everything(), ~paste(unique(na.omit(.)), collapse = "; ")), .groups = "drop") %>%
  mutate(
    across(
      where(is.character),               # todas las columnas de texto
      ~ str_remove(., "; $")             # quita uno o varios ';' al final
    )
  )

## 1.1. Mirar info contradictoria en los metadatos ####
#df con las columnas extraidas del conffile
conffile_collapse <- df_cruzado_colapsado %>%
  select(ENA_experiment_accession, internal_study_name, geographic_location_country:additional_info)

#muestras con diferencias en algunas de las columnas del conffile 
df_filtrado_conffile <- conffile_collapse %>%
  # 1. Selecciona la columna id y todas las columnas que en alguna celda contienen ";"
  select(
    ENA_experiment_accession,
    where(~ any(str_detect(as.character(.), ";"), na.rm = TRUE))
  ) %>%
  # 2. Filtra sólo las filas que en CUALQUIERA de esas columnas (excepto id) tienen ";"
  filter(
    if_any(
      -internal_study_name,
      ~ str_detect(as.character(.), ";")
    )
  )

conffile_contradictorio <- df_filtrado_conffile %>% pull(ENA_experiment_accession)%>% unique()

## 1.2. Mirar info contradictoria en los fenotipos ####
#df con las columnas con info de los antibióticos 
atb_collapse <- df_cruzado_colapsado %>%
  select(ENA_experiment_accession, internal_study_name, ciprofloxacin_laboratory_typing_method:last_col()) 

#diferencias en las columnas phenotype
df_filtrado_phenotype <- atb_collapse %>%
  select(
    ENA_experiment_accession,
    where(~ any(str_detect(as.character(.), ";"), na.rm = TRUE))
  ) %>%
  # 2. Filtra sólo las filas que en CUALQUIERA de esas columnas (excepto id) tienen ";"
  filter(
    if_any(
      ends_with("phenotype"),
      ~ str_detect(as.character(.), ";")
    )
  )

#muestras que tienen info contradictoria 
fenotipo_contradictorio <- df_filtrado_phenotype %>% pull(ENA_experiment_accession)%>% unique()

print(paste("Se deben eliminar las siguientes muestras por contradicciones en los metadatos:", conffile_contradictorio )) 
print(paste("Se deben eliminar las siguientes muestras por contradicciones en los fenotipos:", fenotipo_contradictorio))

## 1.3. Mirar info contradictoria en los AST ####
# df con muestras en info duplicada en cualquiera de las columnas con info de atb
df_filtrado_atb <- atb_collapse %>%
  select(
    ENA_experiment_accession,
    where(~ any(str_detect(as.character(.), ";"), na.rm = TRUE))
  ) %>%
  # 2. Filtra sólo las filas que en CUALQUIERA de esas columnas (excepto id) tienen ";"
  filter(
    if_any(
      -internal_study_name,
      ~ str_detect(as.character(.), ";")
    )
  ) 

# se confirma que hay muchas muestras con experimentos de AST duplicados y queremos conservar la info original 

#2. Conservar la info original de las duplicaciones ####

# Eliminamos las muestras con info contradictoria
df_cruzado_colapsado.scont <- df_cruzado_colapsado %>%
  filter(!is.na(ENA_experiment_accession)) %>%
  filter(!ENA_experiment_accession %in% c(
    "ERX3135001",
    "ERX3134925",
    "ERX451728"))

## 2.1. mirar en qué artículos hay muestras repetidas #### 
df_cruzado_colapsado.scont %>% select(internal_study_name)%>% pull() %>% unique()

#[1] "moradigaravand2018; lipworth2024"               "lipworth2024; moradigaravand2018"              
#[3] "moradigaravand2018; kallonen2017; lipworth2024" "gladstone2021; lipworth2024" 

# Se va a priorizar conservar la info original 

# "moradigaravand2018; lipworth2024"

# ENA_experiment_accession que estan en moradigaravand2018 y lipworth2024
dup_moradigaravand_lipworth <- 
  df_cruzado_colapsado.scont %>% 
  filter(internal_study_name == "moradigaravand2018; lipworth2024") %>%
  select(ENA_experiment_accession) %>%
  pull()

moradigaravand <- 
  df_cruzado %>%
  filter(ENA_experiment_accession %in% dup_moradigaravand_lipworth) %>% 
  filter(internal_study_name == "moradigaravand2018")

lipworth <- 
  df_cruzado %>%
  filter(ENA_experiment_accession %in% dup_moradigaravand_lipworth) %>% 
  filter(internal_study_name == "lipworth2024")

# 1. Unimos los data frames por la clave indicada
moradigaravand_lipworth <- full_join(
  moradigaravand,
  lipworth,
  by = "ENA_experiment_accession",
  suffix = c(".df1", ".df2")
)

# 2. Detectamos columnas en común (excepto la clave)
cols <- setdiff(intersect(names(moradigaravand), names(lipworth)), "ENA_experiment_accession")

# 3. Aplicamos coalesce para cada columna en común
for (col in cols) {
  moradigaravand_lipworth[[col]] <- coalesce(
    moradigaravand_lipworth[[paste0(col, ".df1")]],
    moradigaravand_lipworth[[paste0(col, ".df2")]]
  )
}

# 4. Mantenemos solo la clave + columnas finales coalescidas
moradigaravand_lipworth2 <- moradigaravand_lipworth %>%
  select(ENA_experiment_accession, all_of(cols))

unique(moradigaravand_lipworth2$internal_study_name)

moradigaravand_lipworth2$internal_study_name <- "moradigaravand2018; lipworth2024"

#rm(dup_moradigaravand_lipworth, moradigaravand, lipworth, moradigaravand_lipworth, cols)

# "lipworth2024; moradigaravand2018"

# ENA_experiment_accession que estan en moradigaravand2018 y lipworth2024
dup_lipworth_moradigaravand <- 
  df_cruzado_colapsado.scont %>% 
  filter(internal_study_name == "lipworth2024; moradigaravand2018") %>%
  filter(ENA_experiment_accession != "ERX3134925") %>% # eliminamos la muestra "ERX3134925" pq tiene info contradictoria en el fenotipo
  select(ENA_experiment_accession) %>%
  pull()

moradigaravand2 <- 
  df_cruzado %>%
  filter(ENA_experiment_accession %in% dup_lipworth_moradigaravand) %>% 
  filter(internal_study_name == "moradigaravand2018")

lipworth2 <- 
  df_cruzado %>%
  filter(ENA_experiment_accession %in% dup_lipworth_moradigaravand) %>% 
  filter(internal_study_name == "lipworth2024")

lipworth2 %>% select(ENA_experiment_accession) %>% unique()

# colapsar la informacion de muestras duplicadas en lipworth2024
lipworth2_colapsado <- lipworth2 %>%
  group_by(ENA_experiment_accession) %>%
  summarise(across(everything(), ~na.omit(.)[1]), .groups = "drop")

# 1. Unimos los data frames por la clave indicada
lipworth_moradigaravand <- full_join(
  moradigaravand2,
  lipworth2_colapsado,
  by = "ENA_experiment_accession",
  suffix = c(".df1", ".df2")
)

# 2. Detectamos columnas en común (excepto la clave)
cols2 <- setdiff(intersect(names(moradigaravand2), names(lipworth2_colapsado)), "ENA_experiment_accession")

# 3. Aplicamos coalesce para cada columna en común
for (col in cols2) {
  lipworth_moradigaravand[[col]] <- coalesce(
    lipworth_moradigaravand[[paste0(col, ".df1")]],
    lipworth_moradigaravand[[paste0(col, ".df2")]]
  )
}

# 4. Mantenemos solo la clave + columnas finales coalescidas
lipworth_moradigaravand2 <- lipworth_moradigaravand %>%
  select(ENA_experiment_accession, all_of(cols2)) 

unique(lipworth_moradigaravand2$internal_study_name)

lipworth_moradigaravand2$internal_study_name <- "moradigaravand2018; lipworth2024"

#rm(dup_lipworth_moradigaravand, moradigaravand2, lipworth2, lipworth2_colapsado, lipworth_moradigaravand2, cols2)

# "moradigaravand2018; kallonen2017; lipworth2024"

# ENA_experiment_accession que estan en moradigaravand2018 y kallonen2017
dup_moradigaravand_lipworth_kallonen <- 
  df_cruzado_colapsado.scont %>% 
  filter(internal_study_name == "moradigaravand2018; kallonen2017; lipworth2024") %>%
  select(ENA_experiment_accession) %>%
  pull()

moradigaravand3 <- 
  df_cruzado %>%
  filter(ENA_experiment_accession %in% dup_moradigaravand_lipworth_kallonen) %>% 
  filter(internal_study_name == "moradigaravand2018")

lipworth3 <- 
  df_cruzado %>%
  filter(ENA_experiment_accession %in% dup_moradigaravand_lipworth_kallonen) %>% 
  filter(internal_study_name == "lipworth2024")

kallonen <- 
  df_cruzado %>%
  filter(ENA_experiment_accession %in% dup_moradigaravand_lipworth_kallonen) %>% 
  filter(internal_study_name == "kallonen2017")

# 1. Unimos moradigaravand y kallonen
kallonen_moradigaravand <- full_join(
  kallonen, 
  moradigaravand3,
  by = "ENA_experiment_accession",
  suffix = c(".kal", ".mor")
)

# 2. Luego unimos el resultado con lipworth
moradigaravand_lipworth_kallonen <- full_join(
  kallonen_moradigaravand,
  lipworth3,
  by = "ENA_experiment_accession"
)

# 3. Detectar columnas comunes entre los tres (excepto la clave)
cols3 <- setdiff(
  Reduce(intersect, list(names(kallonen), names(moradigaravand3), names(lipworth3))),
  "ENA_experiment_accession"
)

# 4. Aplicar coalesce en orden de prioridad: mor > kal > lip
for (col in cols3) {
  moradigaravand_lipworth_kallonen[[col]] <- coalesce(
    moradigaravand_lipworth_kallonen[[paste0(col, ".kal")]],
    moradigaravand_lipworth_kallonen[[paste0(col, ".mor")]],
    moradigaravand_lipworth_kallonen[[col]]
  )
}

# 5. Conservar clave y columnas consolidadas
moradigaravand_lipworth_kallonen2 <- moradigaravand_lipworth_kallonen %>%
  select(ENA_experiment_accession, all_of(cols3))

unique(moradigaravand_lipworth_kallonen2$internal_study_name)

moradigaravand_lipworth_kallonen2$internal_study_name <- "kallonen2017; moradigaravand2018; lipworth2024"


# "gladstone2021; lipworth2024"

# ENA_experiment_accession que estan en moradigaravand2018 y lipworth2024
dup_gladstone_lipworth <- 
  df_cruzado_colapsado.scont %>% 
  filter(internal_study_name == "gladstone2021; lipworth2024") %>%
  select(ENA_experiment_accession) %>%
  pull()

gladstone <- 
  df_cruzado %>%
  filter(ENA_experiment_accession %in% dup_gladstone_lipworth) %>% 
  filter(internal_study_name == "gladstone2021")

lipworth4 <- 
  df_cruzado %>%
  filter(ENA_experiment_accession %in% dup_gladstone_lipworth) %>% 
  filter(internal_study_name == "lipworth2024")

# 1. Unimos los data frames por la clave indicada
gladstone_lipworth <- full_join(
  gladstone,
  lipworth4,
  by = "ENA_experiment_accession",
  suffix = c(".df1", ".df2")
)

# 2. Detectamos columnas en común (excepto la clave)
cols4 <- setdiff(intersect(names(gladstone), names(lipworth4)), "ENA_experiment_accession")

# 3. Aplicamos coalesce para cada columna en común
for (col in cols4) {
  gladstone_lipworth[[col]] <- coalesce(
    gladstone_lipworth[[paste0(col, ".df1")]],
    gladstone_lipworth[[paste0(col, ".df2")]]
  )
}

# 4. Mantenemos solo la clave + columnas finales coalescidas
gladstone_lipworth2 <- gladstone_lipworth %>%
  select(ENA_experiment_accession, all_of(cols4))

unique(gladstone_lipworth2$internal_study_name)

gladstone_lipworth2$internal_study_name <- "gladstone2021; lipworth2024"

# 3. Unir las muestras que tenían duplicaciones con las que eran únicas ####
# muestras que solo estaban en un estudio 
df_cruzado_unicos <- df_cruzado %>%
  group_by(ENA_experiment_accession) %>%
  filter(n() == 1) %>%
  ungroup()

# unir muestras unicas con las muestras duplicadas ya arregladas 
ENA_dup <- bind_rows(moradigaravand_lipworth2, lipworth_moradigaravand2, moradigaravand_lipworth_kallonen2, gladstone_lipworth2)

df_cruzado_completo.temp <- bind_rows(ENA_dup, df_cruzado_unicos)

# 4. Conservar en el df solo muestras de E. coli ####

# mirar que ENA_tax_id hay en el df
df_cruzado_completo.temp %>% pull(ENA_tax_id) %>% unique()
# [1] "562"     "573"     "624"     "625"     "564"     "354276"  "1344959" "548"     "300181"  "571"    
#[11] "615"     "588"     "2491876" "561"  

#tras  buscar estos tax id en el NCBI, corresponden a Escherichia 561 (n=315), 
#a E. coli 562 (n=12022) y 2491876 (n=1), a shigella 625 (n=13) y a shigella soneii 624 (n=84)

# filtrar solo las muestras con info de coli o Escherichia
df_cruzado_completo <- 
  df_cruzado_completo.temp %>%
  filter(ENA_tax_id %in% c(561, 562, 2491876, 625, 624))

# csv con toda la info de las muestras únicas sacadas de los artículos 
write.csv(df_cruzado_completo, "C:/Users/rferragud/Documents/Projecte/Tables/df_cruzado_completo_articles.csv", row.names = FALSE)

#comprobar que no ha quedado nada con ;
semicolon <-df_cruzado_completo %>% select(ENA_experiment_accession, -ends_with("testing_notes"), -internal_study_name)%>%
  mutate(across(everything(), as.character)) %>% 
  mutate(row = row_number()) %>%
  pivot_longer(-row, names_to = "col", values_to = "value") %>%
  filter(str_detect(as.character(value), ";"))

rm(df_cruzado_colapsado, df_cruzado_completo.temp, df_cruzado_colapsado.scont, 
   df_cruzado_duplicados, df_cruzado_unicos, df_cruzado)

# 5. Listado con todas las ENA_acessions de nuestras muestras ####
ENA_run_exp2 <- 
  df_cruzado_completo %>%
  select(ENA_experiment_accession, ENA_run_accession, ENA_sample_accession, ENA_secondary_sample_accession)

# df en csv
write.csv(ENA_run_exp2, "C:/Users/rferragud/Documents/Projecte/Tables/ENA_accessions.articles.csv", row.names = FALSE)




