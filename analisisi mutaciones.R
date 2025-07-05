

library(tidyverse)
library(data.table)
library(dplyr)

amr_resf_included_ns2 <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/amrf_gen_fen_ns4.csv", check.names = FALSE)


amr_resf_included_ns2_included <- amr_resf_included_ns2 %>% filter(QC == "included")

# mutaciones en ampicilina detectadas en amrf y resf
mut_resf <- 
  resultados_detallados_resfinder_eucast.ecoli_collapsed %>% 
  separate_rows(resultados_genetic_background, sep= ", ") %>% 
  mutate(resultados_genetic_background = str_remove(resultados_genetic_background, " \\(.*\\)")) %>%
  mutate(
  resf_simple = case_when(
    str_detect(resultados_genetic_background, regex("^blaOXA", ignore_case = TRUE)) ~ resultados_genetic_background,
    str_detect(resultados_genetic_background, regex("^(?i)bla.*-M(-\\d+)?$")) ~
      str_remove(resultados_genetic_background, "-\\d+$"),
    
    # 3) resto de bla… o mcr…: quitar todo a partir del primer guión
    str_detect(resultados_genetic_background, regex("^(?i)(bla|mcr)")) ~
      str_remove(resultados_genetic_background, "-.*$"),
    TRUE ~ resultados_genetic_background
  )
  ) %>%
  # 2) Si hay una palabra seguida de un espacio y la siguiente empieza por "p",
  #    sustituir esa última palabra por "mutpunt".
  mutate(
    resf_simple = str_replace(
      resf_simple,
      # capturamos todo hasta el espacio antes de la palabra que empieza por "p"
      "(.*\\s)(?i:p[^\\s]+)$",
      "\\1mutpunt"
    )
  ) 

mut_amrf <- 
  resultados_detallados_amrf_eucast.ecoli_collapsed %>% 
  separate_rows(resultados_genetic_background, sep= ", ") %>% 
  mutate(
    amrf_simple = case_when(
      str_detect(resultados_genetic_background, regex("^blaOXA", ignore_case = TRUE)) ~ resultados_genetic_background,
      str_detect(resultados_genetic_background, regex("^(?i)bla.*-M(-\\d+)?$")) ~
        str_remove(resultados_genetic_background, "-\\d+$"),
      
      # 3) resto de bla… o mcr…: quitar todo a partir del primer guión
      str_detect(resultados_genetic_background, regex("^(?i)(bla|mcr)")) ~
        str_remove(resultados_genetic_background, "-.*$"),
      TRUE ~ resultados_genetic_background, 
      
    )) %>% 
  mutate( 
    amrf_simple = str_replace(
      amrf_simple,
      "_.*$",            # desde el guión bajo hasta el final
      " mutpunt"         # lo reemplazamos con espacio + mutpunt
    )
  ) 
 



group_cols_resf <- setdiff(names(mut_resf), c("resf_simple", "resultados_genetic_background"))


# 2) Agrupa por esas columnas y colapsa los backgrounds únicos en una cadena separada por comas
mut_resf_collapsed <- mut_resf %>%
  group_by(across(all_of(group_cols_resf))) %>%
  summarise(
    resf_simple = str_c(
      unique(resf_simple),
      collapse = ","
    ),
    .groups = "drop"
  )

group_cols_amrf <- setdiff(names(mut_amrf), c("amrf_simple", "resultados_genetic_background"))

mut_amrf_collapsed <- 
  mut_amrf %>%
  group_by(across(all_of(group_cols_amrf))) %>%
  summarise(
    amrf_simple = str_flatten(
      unique(amrf_simple),
      collapse = ","
    ),
    .groups = "drop"
  )


rep <- resultados_detallados_amrf_eucast.ecoli_collapsed %>% group_by(across(all_of(group_cols_amrf))) %>% count() %>% filter(n>1)

cols_amrf <- colnames(resultados_detallados_amrf_eucast.ecoli_collapsed)

# Vamos a ver resultados generales

# genes que se detectaron 

mut_resf_total <- mut_resf %>% select(resultados_genetic_background:resf_simple) %>% group_by(resultados_genetic_background, resultados_predition_type) %>% count()  

mut_amrf_total <- mut_amrf %>% select(resultados_genetic_background:amrf_simple) %>% group_by(resultados_genetic_background, resultados_predition_type) %>% count()  

mut_resf_total_simple <- mut_resf %>% select(resultados_genetic_background:resf_simple) %>% group_by(resf_simple, resultados_predition_type) %>% count()  

mut_amrf_total_simple <- mut_amrf %>% select(resultados_genetic_background:amrf_simple) %>% group_by(amrf_simple, resultados_predition_type) %>% count()  


# resuktados para atb concretos 

#diferencias_mutaciones <- left_join(mut_amp_resf.generales, mut_amp_amrf, by = "resultados_genetic_background") %>%
#  rename(ResFinder = "n.x", AMRFinderPlus = "n.y") %>%
#  filter(!(AMRFinderPlus == ResFinder & !is.na(AMRFinderPlus) & !is.na(ResFinder)))
#
#anti_join(mut_amp_amrf, mut_amp_resf.generales)
#
#mut_resf_collapsed

amrf_resf_mut <- 
  full_join(mut_resf_collapsed, 
            mut_amrf_collapsed, 
            by = c("antibiotic", "resultados_ENA_sample_accession", "resultados_ENA_run_accession", "class", "subclass", "resultados_S_EUCAST", "resultados_R_EUCAST")) %>%
  select(antibiotic, class, resultados_ENA_sample_accession, resultados_ENA_run_accession, resultados_predition_type.x, resf_simple, resultados_predition_type.y, amrf_simple) %>%
  rename(predition_type_resf = "resultados_predition_type.x", predition_type_amrf = "resultados_predition_type.y") 
         

#amrf_resf_mut <- 
#  full_join(resultados_detallados_resfinder_eucast.ecoli_collapsed, 
#            resultados_detallados_amrf_eucast.ecoli_collapsed, 
#            by = c("antibiotic", "resultados_ENA_sample_accession", "resultados_ENA_run_accession")) %>%
#  select(antibiotic, resultados_ENA_sample_accession, resultados_ENA_run_accession, 
#         resultados_sir.x, resultados_WGS_predicted_phenotype, resultados_genetic_background.x, resultados_predition_type.x, 
#         resultados_sir.y, resultados_genetic_background.y, resultados_predition_type.y) %>%
#  rename(sir_ResF = "resultados_sir.x", predicted_pheno_ResF = "resultados_WGS_predicted_phenotype", 
#         predition_type_ResF = "resultados_predition_type.x", genetic_background_ResF = "resultados_genetic_background.x",
#         sir_AMRF = "resultados_sir.y", genetic_background_AMRF = "resultados_genetic_background.y", predition_type_AMRF = "resultados_predition_type.y"  )



amrf_resf_mut_nc <- 
  full_join(mut_resf, 
            mut_amrf, 
            by = c("antibiotic", "resultados_ENA_sample_accession", "resultados_ENA_run_accession", "class", "subclass", "resultados_S_EUCAST", "resultados_R_EUCAST")) %>%
  select(antibiotic, class, resultados_ENA_sample_accession, resultados_ENA_run_accession, resultados_predition_type.x, resf_simple, resultados_predition_type.y, amrf_simple) %>%
  rename(predition_type_resf = "resultados_predition_type.x", predition_type_amrf = "resultados_predition_type.y") 

amrf_resf_mut_nc2 <- 
  full_join(mut_resf, 
            mut_amrf, 
            by = c("antibiotic", "resultados_ENA_sample_accession", "resultados_ENA_run_accession", "class", "subclass", "resultados_S_EUCAST", "resultados_R_EUCAST")) %>%
  select(antibiotic, class, subclass, resultados_ENA_sample_accession, resultados_ENA_run_accession, resultados_predition_type.x, resf_simple, resultados_predition_type.y, amrf_simple) %>%
  rename(predition_type_resf = "resultados_predition_type.x", predition_type_amrf = "resultados_predition_type.y") 

amrf_resf_mut_nc2_clasific <- amrf_resf_mut_nc2 %>% group_by(antibiotic, class, subclass, resf_simple) %>% count() %>% filter(n>50)
amrf_resf_mut_nc2_clasific_amrf <- amrf_resf_mut_nc2 %>% group_by(antibiotic, class, subclass, amrf_simple) %>% count() %>% filter(n>50)

write.csv(amrf_resf_mut_nc2_clasific, "C:/Users/rferragud/Documents/Projecte/Tables/mut_por_clase_resf.csv", row.names = F)
write.csv(amrf_resf_mut_nc2_clasific_amrf, "C:/Users/rferragud/Documents/Projecte/Tables/mut_por_clase_amrf.csv", row.names = F)

# ampicillin

diferencias_mut_ampicillin1 <- amrf_resf_mut_nc %>% 
  filter(antibiotic == "ampicillin" & predition_type_resf == "TP" & predition_type_amrf == "FN" ) %>%
  group_by(resf_simple)%>%count
  
diferencias_mut_ampicillin2 <- amrf_resf_mut %>% 
  filter(antibiotic == "ampicillin" & predition_type_resf == "TP" & predition_type_amrf == "FN" ) %>%
  group_by(resf_simple)%>%count
  

# cefepime

diferencias_mut_cefepime1 <- amrf_resf_mut_nc %>% 
  filter(antibiotic == "cefepime" & predition_type_resf == "TN" & predition_type_amrf == "FP" ) %>%
  group_by(amrf_simple)%>%count

diferencias_mut_cefepime2 <- amrf_resf_mut %>% 
  filter(antibiotic == "cefepime" & predition_type_resf == "TN" & predition_type_amrf == "FP" ) %>%
  group_by(amrf_simple)%>%count

# piperacillin_tazobactam

# FN

diferencias_mut_piperacillintazobactam1_FN <- amrf_resf_mut_nc %>% 
  filter(antibiotic == "piperacillin-tazobactam" & predition_type_resf == "TP" & predition_type_amrf == "FN" ) %>%
  group_by(resf_simple)%>%count

diferencias_mut_piperacillintazobactam2_FN <- amrf_resf_mut %>% 
  filter(antibiotic == "piperacillin-tazobactam" & predition_type_resf == "TP" & predition_type_amrf == "FN" ) %>%
  group_by(resf_simple)%>%count

# FP

diferencias_mut_piperacillintazobactam1_FP <- amrf_resf_mut_nc %>% 
  filter(antibiotic == "piperacillin-tazobactam" & predition_type_resf == "FP" & predition_type_amrf == "TN" ) %>%
  group_by(resf_simple)%>%count

diferencias_mut_piperacillintazobactam2_FP <- amrf_resf_mut %>% 
  filter(antibiotic == "piperacillin-tazobactam" & predition_type_resf == "FP" & predition_type_amrf == "TN" ) %>%
  group_by(resf_simple)%>%count


# fosfomycin

# FN

diferencias_mut_fosfomycin1_FN <- amrf_resf_mut_nc %>% 
  filter(antibiotic == "fosfomycin" & predition_type_resf == "FN" & predition_type_amrf == "TP" ) %>%
  group_by(amrf_simple)%>%count

diferencias_mut_fosfomycin2_FN <- amrf_resf_mut %>% 
  filter(antibiotic == "fosfomycin" & predition_type_resf == "FN" & predition_type_amrf == "TP" ) %>%
  group_by(amrf_simple)%>%count

# FP

diferencias_mut_fosfomycin1_FP <- amrf_resf_mut_nc %>% 
  filter(antibiotic == "fosfomycin" & predition_type_amrf == "FP" & predition_type_resf == "TN" ) %>%
  group_by(amrf_simple)%>%count

diferencias_mut_fosfomycin2_FP <- amrf_resf_mut %>% 
  filter(antibiotic == "fosfomycin" & predition_type_amrf == "FP" & predition_type_resf == "TN" ) %>%
  group_by(amrf_simple)%>%count




  