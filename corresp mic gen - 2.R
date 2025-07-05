
library(tidyverse)
library(data.table)
library(dplyr)
library(readODS)
library(epiR)

broth_microdil_sir_complete.new

# # Función para sacar las predicciones 
# analizar_antibiotico <- function(antibiotic, program_class, sir_df, program) {
#   
#   # Filtrar EUCAST para cada antibiótico
#   pheno_filtered <- sir_df %>%
#     filter(antibiotic == !!antibiotic)
#   # Construimos un patrón para que detecte antibiotic y resfinder_class (si no es NA)
#   pattern <- if (!is.na(program_class)) {
#     paste0(program_class, "|", antibiotic)
#   } else {
#     antibiotic
#   }
#   
#   # 3) Filtrar ResFinder usando el patrón adecuado
#   program_filtered <- program %>%
#     filter(str_detect(
#       Subclass,
#       regex(pattern, ignore_case = TRUE)
#     ))
#   
#   # Join usando muestras únicas
#   cruzado <- pheno_filtered %>%
#     left_join(
#       program_filtered %>%
#         select(ENA_sample_accession, genetic_background),
#       by = "ENA_sample_accession")
#   
#   return(cruzado)
# }

# función resfinder
analizar_antibiotico_resf.sir <- function(antibiotic, program_class, sir_df, program) {
  
  # Filtrar fenotipos para cada antibiótico
  pheno_filtered <- sir_df %>%
    filter(antibiotic == !!antibiotic)
  # Construimos un patrón para que detecte antibiotic y resfinder_class (si no es NA)
  pattern <- if (!is.na(program_class)) {
    paste0("^", program_class,"$", "|", "^",antibiotic,"$")
  } else {
    paste0("^", antibiotic, "$")
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
      by = "ENA_sample_accession") 
  
  return(cruzado)
}

# función amrfinder

analizar_antibiotico_amrf.sir <- function(antibiotic, program_class, sir_df, program) {
  
  # Filtrar fenotipo para cada antibiótico
  pheno_filtered <- sir_df %>%
    filter(antibiotic == !!antibiotic)
  # Construimos un patrón para que detecte antibiotic y program_class (si no es NA)
  pattern <- if (!is.na(program_class)) {
    paste0("^", program_class,"$", "|", "^",antibiotic,"$")
  } else {
    paste0("^", antibiotic, "$")
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
      by = "ENA_sample_accession")
  
  return(cruzado)
}
# resultados amrf
corresp_atb_amrf <- read_ods("C:/Users/rferragud/Documents/Projecte/corresp_antibiotic_resf_amrf.new.ods", sheet = 3)

resultados_detallados_amrf_sir <- corresp_atb_amrf %>%
  mutate(
    resultados = map2(
      antibiotic,
      program_class,
      analizar_antibiotico_amrf.sir,
      sir_df = broth_microdil_sir_complete.new, 
      program = amrf_collapsed
    )
  ) %>%
  unnest(resultados, names_sep = "_") %>% 
  rename(amrf = resultados_genetic_background) %>%
  select(resultados_internal_id:amrf)

# resultados resfinder
corresp_atb_resfinder.ecoli <- read_ods("C:/Users/rferragud/Documents/Projecte/corresp_antibiotic_resf_amrf.new.ods", sheet = 2)

resultados_detallados_resf_sir <- corresp_atb_resfinder.ecoli %>%
  mutate(
    resultados = map2(
      antibiotic,
      program_class,
      analizar_antibiotico_resf.sir,
      sir_df = broth_microdil_sir_complete.new, 
      program = resf_collapsed
    )
  ) %>%
  unnest(resultados, names_sep = "_") %>% 
  rename(resf = resultados_genetic_background) %>%
  select(resultados_internal_id:resf)



amr_resf <- full_join(resultados_detallados_amrf_sir, resultados_detallados_resf_sir)


amr_resf_included <- amr_resf %>% left_join(
  all_summary_def_qc_clean_accession %>%
    select(ENA_sample_accession) %>%
    mutate(QC = "included"),
  by = c("resultados_ENA_sample_accession" = "ENA_sample_accession")
)

amr_resf_included_ns <- amr_resf_included %>%
  mutate(resultados_EUCAST_2025 = if_else(as.character(resultados_EUCAST_2025) == "I", "R", as.character(resultados_EUCAST_2025)),
         resultados_CLSI_2025 = if_else(as.character(resultados_CLSI_2025) == "I", "R", as.character(resultados_CLSI_2025)),
         resultados_ECOFF_2025 = if_else(as.character(resultados_ECOFF_2025) == "I", "R", as.character(resultados_ECOFF_2025)),
  ) %>%
  select(resultados_ENA_sample_accession:resultados_antibiotic, resultados_mic:resultados_mo, resultados_ECOFF_2025:QC, resultados_S_ECOFF:resultados_S_CLSI)

amr_resf_included %>% filter(QC == "included" & !is.na(resultados_EUCAST_2025))  %>% distinct(resultados_ENA_sample_accession)


#hola1 <- amr_resf %>% 
#  filter(!is.na(resultados_EUCAST_2025)) %>% filter(!is.na(resf)) %>% distinct(resultados_ENA_sample_accession)
#
#hola2 <- resultados_detallados_resfinder_eucast.ecoli %>% 
#  filter(!is.na(resultados_sir)) %>% filter(!is.na(resultados_genetic_background)) %>% distinct(resultados_ENA_sample_accession)
#
#anti_join(hola2, hola1)

# TODO posar en dos columnes més el fenotip de ecoff i clsi i passar els intermediate a resistant

#broth_microdil_sir_EUCAST.R <- broth_microdil_sir_EUCAST %>% 
#  rename(sir = "EUCAST_2025") %>%
#  mutate(EUCAST = if_else(as.character(sir) == "I", "R", as.character(sir)))
#
#broth_microdil_sir_ECOFF.R <- broth_microdil_sir_ECOFF %>%
#  rename(sir = "ECOFF_2025") %>%
#  mutate(ECOFF = if_else(as.character(sir) == "I", "R", as.character(sir)))
#
#broth_microdil_sir_CLSI.R <- broth_microdil_sir_CLSI %>%
#  rename(sir = "CLSI_2025") %>%
#  mutate(CLSI = if_else(as.character(sir) == "I", "R", as.character(sir)))
#
#broth_microdil_sir2 <- 
#  broth_microdil_sir_EUCAST.R %>%
#  left_join(
#    broth_microdil_sir_ECOFF.R %>% select(ENA_sample_accession, ENA_run_accession, antibiotic, ECOFF),
#  by = c("ENA_sample_accession", "ENA_run_accession", "antibiotic")
#) %>% 
#  left_join(
#    broth_microdil_sir_CLSI.R %>% select(ENA_sample_accession, ENA_run_accession, antibiotic, CLSI),
#    by = c("ENA_sample_accession", "ENA_run_accession", "antibiotic")) %>%
#  select(internal_id:antibiotic, laboratory_typing_method, units, ab, mo, mic, EUCAST:CLSI) 
#
#resultados_detallados_amrf_sir %>% distinct()

amr_resf_included2 <- amr_resf_included %>% select(resultados_ENA_sample_accession:QC) %>% 
  filter(!is.na(resultados_EUCAST_2025)|!is.na(resultados_CLSI_2025)|!is.na(resultados_ECOFF_2025))

write.csv(amr_resf_included2, "C:/Users/rferragud/Documents/Projecte/Tables/amrf_gen_fen4.csv", row.names = FALSE)

amr_resf_included_ns2 <- amr_resf_included_ns %>% select(resultados_ENA_sample_accession:resultados_S_CLSI) %>% 
  filter(!is.na(resultados_EUCAST_2025)|!is.na(resultados_CLSI_2025)|!is.na(resultados_ECOFF_2025))


write.csv(amr_resf_included_ns2, "C:/Users/rferragud/Documents/Projecte/Tables/amrf_gen_fen_ns4.csv", row.names = FALSE)


amr_resf_included2 %>% filter(!is.na(resultados_EUCAST_2025) & QC == "included") %>% distinct(resultados_ENA_sample_accession)

amr_resf_included_ns2 %>% filter(QC == "included") %>% distinct(resultados_ENA_sample_accession)

amr_resf_included_ns2_mutsi <- amr_resf_included_ns2 %>%
  mutate(
  resf_simple = case_when(
  str_detect(resf, regex("^blaOXA", ignore_case = TRUE)) ~ resultados_genetic_background,
  str_detect(resultados_genetic_background, regex("^(?i)bla.*-M(-\\d+)?$")) ~
    str_remove(resultados_genetic_background, "-\\d+$"),
  
  # 3) resto de bla… o mcr…: quitar todo a partir del primer guión
  str_detect(resultados_genetic_background, regex("^(?i)(bla|mcr)")) ~
    str_remove(resultados_genetic_background, "-.*$"),
  TRUE ~ resf
)) %>%
  # 2) Si hay una palabra seguida de un espacio y la siguiente empieza por "p",
  #    sustituir esa última palabra por "mutpunt".
  mutate(
    resf_simple = str_replace(
      resf_simple,
      # capturamos todo hasta el espacio antes de la palabra que empieza por "p"
      "(.*\\s)(?i:p[^\\s]+)$",
      "\\1mutpunt"
    )
  ) %>% 
  mutate(
    amrf_simple = case_when(
      str_detect(amrf, regex("^blaOXA", ignore_case = TRUE)) ~ amrf,
      str_detect(amrf, regex("^(?i)bla.*-M(-\\d+)?$")) ~
        str_remove(amrf, "-\\d+$"),
      
      # 3) resto de bla… o mcr…: quitar todo a partir del primer guión
      str_detect(amrf, regex("^(?i)(bla|mcr)")) ~
        str_remove(amrf, "-.*$"),
      TRUE ~ amrf, 
      
  )) %>% 
  mutate( 
      amrf = str_replace(
        amrf,
        "_.*$",            # desde el guión bajo hasta el final
        " mutpunt"         # lo reemplazamos con espacio + mutpunt
      )
    ) 
  
  
  
  
  