
library(tidyverse)
library(data.table)
library(dplyr)
library(readODS)
library(epiR)
library(AMR)

rm(list=ls())

rm(all_summary_def_qc_sylph.raw, amrf_results_accession_cruzado1, amrf_results_accession_cruzado, amrf_results_assembly, broth_microdil, EUCAST_sir, EUCAST_sir_clean, EUCAST_sir_clean2, CLSI_sir, CLSI_sir_clean, CLSI_sir_clean2, ECOFF_sir, ECOFF_sir_clean, ECOFF_sir_clean2)
AMR::set_AMR_locale("English")

amrf_results_accession_cruzado.completo <- read.csv( "C:/Users/rferragud/Documents/Projecte/Tables/amrf_results_accession_cruzado.completo.csv", check.names = F)%>%
  rename(genetic_background = Element.symbol) %>%
  mutate(
    Subclass = str_replace_all(
      Subclass,
      regex("CLAVULANIC_ACID", ignore_case = TRUE),
      "clavulanic acid"
    )
  )

resf_results_accession_cruzado <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/resf_results_accession_cruzado.csv", check.names =F)


# Nueva función para transformar mic a fenotipo

# tabla con los breakpoints de los atb que vamos a utilizar 
breakpoints_all_arreglada <- read_ods("C:/Users/rferragud/Documents/Projecte/Tables/breakpoints_all_arreglada.ods")

# misma tabla pero solo c0n los atb y breakpoints
breakpoints_R_S <- breakpoints_all_arreglada %>% select(ab_name, R_ECOFF, R_EUCAST, R_CLSI, S_ECOFF, S_EUCAST, S_CLSI)%>%
  mutate(across(everything(), ~ na_if(., "NA")))%>%
  mutate(R_ECOFF = as.double(R_ECOFF), 
         R_EUCAST = as.double(R_EUCAST), 
         R_CLSI = as.double(R_CLSI), 
         S_ECOFF = as.double(S_ECOFF), 
         S_EUCAST = as.double(S_EUCAST), 
         S_CLSI = as.double(S_CLSI))


# df con todos los mic que tenemos
broth_microdil <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/broth_microdil.csv", check.names = F)

# df con todas las muestras con mic y sus breakpoints
broth_microdil_BP <- left_join(broth_microdil, breakpoints_R_S, by = c("antibiotic" = "ab_name"))

# funcion para convertir mic a fenotipo eliminando mic que no se pueden interpretar

broth_microdil_sir_complete.new <- broth_microdil_BP %>%
  #filter(!internal_study_name == "rebelo2022") %>% 
  # 1) limpieza básica y unificación de formatos
  mutate(
    meas = measurement_combined %>%
      str_remove("^==|^=") %>%
      str_replace_all("≤", "<=") %>%
      str_replace(",", ".")
  ) %>%
  # 2) descartar combinaciones raras
  filter(!str_detect(meas, "[/+]")) %>%
  # 4) parsear mic y convertir nombres a códigos
  mutate(
    mic = as.mic(meas),
    ab  = as.ab(antibiotic), 
    mo = as.mo(ENA_scientific_name)
  ) %>%
  mutate( # quita los mic que tienen un < y numero mayor al R_breakpoint
    new_mic_ecoff = case_when(
      str_detect(mic, "^<")  & as.double(mic) > S_ECOFF ~ NA_mic_,
      str_detect(mic, "^<=")  & as.double(mic) > S_ECOFF ~ NA_mic_,
      str_detect(mic, "^>")  & as.double(mic) < S_ECOFF ~ NA_mic_,
      str_detect(mic, "^>=")  & as.double(mic) < S_ECOFF ~ NA_mic_,
      str_detect(mic, "^>")  & as.double(mic) > S_ECOFF ~ as.mic(as.double(mic)),
      str_detect(mic, "^>=") & as.double(mic) > S_ECOFF ~ as.mic(as.double(mic)),
      str_detect(mic, "^<=") & as.double(mic) == S_ECOFF ~ as.mic(as.double(mic)),
      #str_detect(mic, "^<")  & as.double(mic) > R_ECOFF ~ NA_mic_,
      #str_detect(mic, "^>")  & as.double(mic) < R_ECOFF ~ NA_mic_,
      #str_detect(mic, "^<=") & as.double(mic) > R_ECOFF ~ NA_mic_,
      #str_detect(mic, "^>=") & as.double(mic) < R_ECOFF ~ NA_mic_,
      #str_detect(mic, "^<=") & as.double(mic) > S_ECOFF ~ NA_mic_,
      # quita los valores que son <= al S_breakpoint
      TRUE ~ mic
    ), 
    new_mic_eucast = case_when(
      str_detect(mic, "^<")  & as.double(mic) > S_EUCAST ~ NA_mic_,
      str_detect(mic, "^<=")  & as.double(mic) > S_EUCAST ~ NA_mic_,
      str_detect(mic, "^>")  & as.double(mic) < S_EUCAST ~ NA_mic_,
      str_detect(mic, "^>=")  & as.double(mic) < S_EUCAST ~ NA_mic_,
      str_detect(mic, "^>")  & as.double(mic) > S_EUCAST ~ as.mic(as.double(mic)),
      str_detect(mic, "^>=") & as.double(mic) > S_EUCAST ~ as.mic(as.double(mic)),
      str_detect(mic, "^<=") & as.double(mic) == S_EUCAST ~ as.mic(as.double(mic)),
      #str_detect(mic, "^<")  & as.double(mic) > R_EUCAST ~ NA_mic_,
      #str_detect(mic, "^>")  & as.double(mic) < R_EUCAST ~ NA_mic_,
      #str_detect(mic, "^<=") & as.double(mic) > R_EUCAST ~ NA_mic_,
      #str_detect(mic, "^>=") & as.double(mic) < R_EUCAST ~ NA_mic_,
      #str_detect(mic, "^<=") & as.double(mic) > S_EUCAST ~ NA_mic_,
      # quita los valores que son <= al S_breakpoint
      TRUE ~ mic
    ), 
    new_mic_clsi = case_when(
      str_detect(mic, "^<")  & as.double(mic) > S_CLSI ~ NA_mic_,
      str_detect(mic, "^<=")  & as.double(mic) > S_CLSI ~ NA_mic_,
      str_detect(mic, "^>")  & as.double(mic) < S_CLSI ~ NA_mic_,
      str_detect(mic, "^>=")  & as.double(mic) < S_CLSI ~ NA_mic_,
      str_detect(mic, "^>")  & as.double(mic) > S_CLSI ~ as.mic(as.double(mic)),
      str_detect(mic, "^>=") & as.double(mic) > S_CLSI ~ as.mic(as.double(mic)),
      str_detect(mic, "^<=") & as.double(mic) == S_CLSI ~ as.mic(as.double(mic)),
      #str_detect(mic, "^<")  & as.double(mic) > R_CLSI ~ NA_mic_,
      #str_detect(mic, "^>")  & as.double(mic) < R_CLSI ~ NA_mic_,
      #str_detect(mic, "^<=") & as.double(mic) > R_CLSI ~ NA_mic_,
      #str_detect(mic, "^>=") & as.double(mic) < R_CLSI ~ NA_mic_,
      #str_detect(mic, "^<=") & as.double(mic) > S_CLSI ~ NA_mic_,
      # quita los valores que son <= al S_breakpoint
      TRUE ~ mic
    )
  )%>% 
  group_by(ab, mo) %>%
  mutate(
    ECOFF_2025 = as.sir(
      x         = new_mic_ecoff,
      ab        = ab[1],
      mo        = mo[1],
      guideline = "EUCAST", 
      breakpoint_type = "ECOFF",
      include_PKPD = FALSE
    ),
    EUCAST_2025 = as.sir(
      x         = new_mic_eucast,
      ab        = ab[1],
      mo        = mo[1],
      guideline = "EUCAST", 
      breakpoint_type = "human",
      include_PKPD = FALSE
    ),
    CLSI_2025 = as.sir(
      x         = new_mic_clsi,
      ab        = ab[1],
      mo        = mo[1],
      guideline = "CLSI", 
      breakpoint_type = "human",
      include_PKPD = FALSE
    )
  ) %>%
  ungroup() 


write.csv(broth_microdil_sir_complete.new, "C:/Users/rferragud/Documents/Projecte/Tables/broth_microdil_sir_complete3.new.csv", row.names = F)

broth_microdil_sir_complete.new <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/broth_microdil_sir_complete3.new.csv", check.names = F)


# Comprobar que están correctas las transformaciones de mic a fenotipo 

broth_microdil_sir_complete.new.long <- broth_microdil_sir_complete.new %>%
  select(antibiotic, R_ECOFF:CLSI_2025)%>%
  pivot_longer(
    cols = c(ECOFF_2025, EUCAST_2025, CLSI_2025),
    names_to = "guideline",
    values_to = "phenotype"
  ) %>%
  mutate(
    guideline = str_remove(guideline, "_2025"),
    breakpoint_R = case_when(
      guideline == "ECOFF"  ~ R_ECOFF,
      guideline == "EUCAST" ~ R_EUCAST,
      guideline == "CLSI"   ~ R_CLSI,
      TRUE ~ NA_real_
    ),
    breakpoint_S = case_when(
      guideline == "ECOFF"  ~ S_ECOFF,
      guideline == "EUCAST" ~ S_EUCAST,
      guideline == "CLSI"   ~ S_CLSI,
    ))

broth_microdil_sir_complete.new.summary <- broth_microdil_sir_complete.new.long %>%
  count(antibiotic, mic, guideline, phenotype, breakpoint_S, breakpoint_R, name = "n")

summary_ECOFF <- broth_microdil_sir_complete.new.summary %>% filter(guideline== "ECOFF" & !(is.na(breakpoint_R)))
summary_CLSI <- broth_microdil_sir_complete.new.summary %>% filter(guideline== "CLSI" & !(is.na(breakpoint_R)))
summary_EUCAST <- broth_microdil_sir_complete.new.summary %>% filter(guideline== "EUCAST" & !(is.na(breakpoint_R)))

# separar por guideline el df

EUCAST_sir <- 
  broth_microdil_sir_complete.new %>% 
  select(ENA_experiment_accession:ENA_sample_accession, antibiotic, EUCAST_2025, mic, new_mic_eucast, S_EUCAST, R_EUCAST) %>%
  filter(!is.na(EUCAST_2025))
 

ECOFF_sir <- 
  broth_microdil_sir_complete.new %>% 
  select(ENA_experiment_accession:ENA_sample_accession, antibiotic, ECOFF_2025, mic, new_mic_ecoff, S_ECOFF, R_ECOFF) %>%
  filter(!is.na(ECOFF_2025))

CLSI_sir <- 
  broth_microdil_sir_complete.new %>% 
  select(ENA_experiment_accession:ENA_sample_accession, antibiotic, CLSI_2025, mic, new_mic_clsi, S_CLSI, R_CLSI) %>%
  filter(!is.na(CLSI_2025)) 
  

# 5. Eliminar de las tablas de los fenotipos las muestras que no pasan el QC ####

# 5.1 EUCAST

# EUCAST_sir <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/EUCAST_sir.csv", check.names = FALSE)

# primero obtener las ENA_sample_accession del EUCAST_sir con el all_summary_def_qc_clean_accession 

EUCAST_sir %>% filter(is.na(ENA_sample_accession) & is.na(ENA_run_accession) & !is.na(ENA_experiment_accession))

EUCAST_sir %>% filter(is.na(ENA_sample_accession) & !is.na(ENA_run_accession))

EUCAST_sir %>% filter(!is.na(ENA_sample_accession) & is.na(ENA_run_accession))

# las muestras que tienen ENA_experiment_accession van a tener tb sample o run accession así que se van a poder cruzar solo con sample y run 
# todas tienen ENA_sample_accession

# Obtener solo una columna que tenga una accession de cada muestra
EUCAST_coaslesce <- EUCAST_sir %>% mutate(ENA_accessions = coalesce(ENA_sample_accession, ENA_run_accession)) %>%
  select(-ENA_sample_accession, -ENA_run_accession)

EUCAST_coaslesce %>% distinct(ENA_accessions)

all_summary_def_qc_clean <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/all_summary_def_qc_clean", check.names = F)

# columnas de identificadores de las muestras que pasan el qc
all_summary_def_qc_clean_accession <- 
  all_summary_def_qc_clean %>% 
  select(ENA_sample_accession, ENA_run_accession, assembly)

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
# 12458/12434/12625
sobrantes_acc_euc %>% unique()
# 2706//2710 (se han eliminado por no pasar los criterios del QC)
EUCAST_coaslesce %>% pull(ENA_accessions)%>% unique()
# 15164/15140/15335

# comprobar que las muestras de sobrantes estan en el df de las muestras que no pasan el QC
remove_sam <- all_summary_def_qc_remove %>% pull(ENA_sample_accession)%>% unique()
remove_run <- all_summary_def_qc_remove %>% pull(ENA_run_accession)%>% unique()

sobrantes_euc %>% filter(ENA_accessions %in% remove_sam) %>% pull(ENA_accessions) %>% unique()
sobrantes_euc %>% filter(!ENA_accessions %in% remove_sam & ENA_accessions %in% remove_run) %>% 
  pull(ENA_accessions) %>% unique()


# ECOFF

# ECOFF_sir <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/ECOFF_sir.csv", check.names = FALSE)

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
# 12253//12420
sobrantes_acc_ecoff %>% unique()
# 2701//2705 (se han eliminado por no pasar los criterios del QC)
ECOFF_coaslesce %>% pull(ENA_accessions)%>% unique()
# 14954//15125

sobrantes_ecoff %>% filter(ENA_accessions %in% remove_sam) %>% pull(ENA_accessions) %>% unique()
sobrantes_ecoff %>% filter(!ENA_accessions %in% remove_sam & ENA_accessions %in% remove_run) %>% 
  pull(ENA_accessions) %>% unique()

# CLSI 

# CLSI_sir <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/CLSI_sir.csv", check.names = FALSE)

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
# 12162/12138/12601
sobrantes_acc_clsi %>% unique()
# 2704//2710 (se han eliminado por no pasar los criterios del QC)
CLSI_coaslesce %>% pull(ENA_accessions)%>% unique()
# 14866/14842/15311

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

count_comparison.completa.clean_ns <- count_comparison.completa.clean %>% 
  mutate(R_I = R + I)

 write.csv(count_comparison.clean, "C:/Users/rferragud/Documents/Projecte/Tables/count_comparison.clean3.csv", row.names = FALSE)
 write.csv(count_comparison.completa.clean, "C:/Users/rferragud/Documents/Projecte/Tables/count_comparison.completa.clean3.csv", row.names = FALSE)



# Criterios de inclusión de los antibióticos 

# detección de mic en más de 500 muestras 
# detección de mic en más de 100 muestras para cada tipo de fenotipo (R/S)

count_comparison.completa.filtro1_ns <- 
  count_comparison.completa.clean_ns %>% 
  mutate(
    filtro = if_else(
      total < 500,
      "< 500 isolates with MIC",
      NA_character_
    )
  )


count_comparison.completa.filtro2_ns <- count_comparison.completa.filtro1_ns %>%
  mutate(
    filtro = if_else(
      is.na(filtro) & (R_I < 100 | S < 100 | is.na(R) | is.na(S)), 
      "< 100 R_I/S isolates", 
      filtro
    )
  )


count_comparison.completa.filtro3_ns <- count_comparison.completa.filtro2_ns %>%
  mutate(
    filtro = if_else(
      is.na(filtro), 
      "included", 
      filtro
    )
  )

 write.csv(count_comparison.completa.filtro3_ns, "C:/Users/rferragud/Documents/Projecte/Tables/atb_included_ns3.csv")

 #-----------------------------------------------------------------------
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
 
 #-----------------------------------------------------------------------
 
 
atb_ecoff <- count_comparison.completa.filtro3_ns %>% 
  filter(testing_standard == "ecoff" & filtro == "included")

atb_ecoff_names <- atb_ecoff %>% arrange(antibiotic) %>% pull(antibiotic)

atb_eucast <- count_comparison.completa.filtro3_ns %>% 
  filter(testing_standard == "eucast" & filtro == "included")

atb_eucast_names <- atb_eucast %>% arrange(antibiotic) %>% pull(antibiotic)

atb_clsi <- count_comparison.completa.filtro3_ns %>% 
  filter(testing_standard == "clsi" & filtro == "included")

atb_clsi_names <- atb_clsi %>% arrange(antibiotic) %>% pull(antibiotic)

atb_total <- count_comparison.completa.filtro3_ns %>% 
  filter(filtro == "included") %>% arrange(antibiotic) %>% distinct(antibiotic) %>% print(n = Inf)

pull(atb_total)

analizar_antibiotico_resf <- function(antibiotic, program_class, sir_df, program) {
  
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


resf_collapsed <- resf_results_accession_cruzado %>%
  group_by(
    ENA_sample_accession,
    ENA_run_accession,
    assembly,
    Subclass, 
    Class, 
    WGS_predicted_phenotype
  ) %>%
  summarise(
    genetic_background = na_if(
      paste(
        na.omit(unique(genetic_background)),  # quitamos NAs y duplicados
        collapse = ", "
      ),
      ""                                     # si no queda nada, "" → NA
    ), 
    .groups = "drop"
  )

# resultados resfinder
# a) eucast
resultados_detallados_resfinder_eucast.ecoli_collapsed <- corresp_atb_euc_resf.ecoli %>%
  mutate(
    resultados = map2(
      antibiotic,
      program_class,
      analizar_antibiotico_resf,
      sir_df = EUCAST_sir_clean2, 
      program = resf_collapsed
    )
  ) %>%
  unnest(resultados, names_sep = "_")


# b) ecoff
resultados_detallados_resfinder_ecoff.ecoli_collapsed <- corresp_atb_ecoff_resf.ecoli %>%
  mutate(
    resultados = map2(
      antibiotic,
      program_class,
      analizar_antibiotico_resf,
      sir_df = ECOFF_sir_clean2, 
      program = resf_collapsed
    )
  ) %>%
  unnest(resultados, names_sep = "_")

# c) clsi
resultados_detallados_resfinder_clsi.ecoli_collapsed <- corresp_atb_clsi_resf.ecoli %>%
  mutate(
    resultados = map2(
      antibiotic,
      program_class,
      analizar_antibiotico_resf,
      sir_df = CLSI_sir_clean2, 
      program = resf_collapsed
    )
  ) %>%
  unnest(resultados, names_sep = "_")

# csv con los resultados 

write.csv(resultados_detallados_resfinder_eucast.ecoli_collapsed, 
          "C:/Users/rferragud/Documents/Projecte/Tables/resultados_detallados_resfinder_eucast.ecoli_collapsed3.csv", 
          row.names = F)

write.csv(resultados_detallados_resfinder_ecoff.ecoli_collapsed, 
          "C:/Users/rferragud/Documents/Projecte/Tables/resultados_detallados_resfinder_ecoff.ecoli_collapsed3.csv", 
          row.names = F)

write.csv(resultados_detallados_resfinder_clsi.ecoli_collapsed, 
          "C:/Users/rferragud/Documents/Projecte/Tables/resultados_detallados_resfinder_clsi.ecoli_collapsed3.csv",
          row.names = F)

#---------------------------------- GRÁFICOS ----------------------------------#

# RESFINDER ECOLI 

# EUCAST
resultados_detallados_resfinder_eucast.ecoli_collapsed

graf_resf_euc_eco <- resultados_detallados_resfinder_eucast.ecoli_collapsed %>% 
  filter(!is.na(resultados_predition_type) & resultados_predition_type != "unknown") %>%
  mutate(
    pred_type_group = case_when(
      resultados_predition_type %in% c("TP", "TN") ~ "TP + TN",
      resultados_predition_type == "FP"           ~ "FP",
      resultados_predition_type == "FN"           ~ "FN"
    )
  ) %>%
  count(antibiotic, pred_type_group) %>%
  ggplot() +
  geom_col(aes(x = antibiotic, y = n, fill = pred_type_group)) +
  labs(
    x = "Antibiótico",
    y = "Número de predicciones",
    title = "Predicciones por antibiótico EUCAST",
    fill = "Tipo de predicción" 
  ) +
  scale_fill_manual(values = c(
    "TP + TN" = "lightgreen",
    "FN" = "red",
    "FP" = "darkred"
  )
  ) +
  #scale_y_continuous(
  #  breaks      = seq(0, by = 3000),
  #  minor_breaks= seq(0, by = 1500),
  #  limits      = c(0, 11250)
  #) +
  theme_minimal() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(r = 10), size = 12), 
    axis.title.x = element_text(vjust = 0.5, hjust = 0.5, size = 12),
    legend.title   = element_text(size = 10),    # tamaño del título de la leyenda
    legend.text    = element_text(size = 8),    # tamaño del texto de la leyenda
    legend.key.size= unit(0.8, "lines"),
    plot.margin     = margin(t = 5, r = 2, b = 5, l = 5, unit = "mm")
  )

# ECOFF

resultados_detallados_resfinder_ecoff.ecoli_collapsed

graf_resf_ecoff_eco <- resultados_detallados_resfinder_ecoff.ecoli_collapsed %>% 
  filter(!is.na(resultados_predition_type) & resultados_predition_type != "unknown") %>%
  mutate(
    pred_type_group = case_when(
      resultados_predition_type %in% c("TP", "TN") ~ "TP + TN",
      resultados_predition_type == "FP"           ~ "FP",
      resultados_predition_type == "FN"           ~ "FN"
    )
  ) %>%
  count(antibiotic, pred_type_group) %>%
  ggplot() +
  geom_col(aes(x = antibiotic, y = n, fill = pred_type_group)) +
  labs(
    x = "Antibiótico",
    y = "Número de predicciones",
    title = "Predicciones por antibiótico ECOFF",
    fill = "Tipo de predicción" 
  ) +
  scale_fill_manual(values = c(
    "TP + TN" = "lightgreen",
    "FN" = "red",
    "FP" = "darkred"
  )
  ) +
 #scale_y_continuous(
 #  breaks      = seq(0, by = 3000),
 #  minor_breaks= seq(0, by = 1500),
 #  limits      = c(0, 11250)
 # ) +
  theme_minimal() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(r = 10), size = 12), 
    axis.title.x = element_text(vjust = 0.5, hjust = 0.5, size = 12),
    legend.title   = element_text(size = 10),    # tamaño del título de la leyenda
    legend.text    = element_text(size = 8),    # tamaño del texto de la leyenda
    legend.key.size= unit(0.8, "lines"),
    plot.margin     = margin(t = 5, r = 2, b = 5, l = 5, unit = "mm")
  )


# CLSI

resultados_detallados_resfinder_clsi.ecoli_collapsed

graf_resf_clsi_eco <- resultados_detallados_resfinder_clsi.ecoli_collapsed %>% 
  filter(!is.na(resultados_predition_type) & resultados_predition_type != "unknown") %>%
  mutate(
    pred_type_group = case_when(
      resultados_predition_type %in% c("TP", "TN") ~ "TP + TN",
      resultados_predition_type == "FP"           ~ "FP",
      resultados_predition_type == "FN"           ~ "FN"
    )
  ) %>%
  count(antibiotic, pred_type_group) %>%
  ggplot() +
  geom_col(aes(x = antibiotic, y = n, fill = pred_type_group)) +
  labs(
    x = "Antibiótico",
    y = "Número de predicciones",
    title = "Predicciones por antibiótico CLSI",
    fill = "Tipo de predicción" 
  ) +
  scale_fill_manual(values = c(
    "TP + TN" = "lightgreen",
    "FN" = "red",
    "FP" = "darkred"
  )
  ) +
 #scale_y_continuous(
 #  breaks      = seq(0, by = 3000),
 #  minor_breaks= seq(0, by = 1500),
 #  limits      = c(0, 11250)
 #) +
  theme_minimal() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(r = 10), size = 12), 
    axis.title.x = element_text(vjust = 0.5, hjust = 0.5, size = 12),
    legend.title   = element_text(size = 10),    # tamaño del título de la leyenda
    legend.text    = element_text(size = 8),    # tamaño del texto de la leyenda
    legend.key.size= unit(0.8, "lines"),
    plot.margin     = margin(t = 5, r = 2, b = 5, l = 5, unit = "mm")
  )


# -------------------------------- TABLAS -------------------------------------#



compute_diagnostics_by_ab <- function(data) {
  data %>%
    # 1) Agrupamos por tu columna de antibiótico
    group_by(resultados_antibiotic) %>%
    # 2) Para cada antibiótico, aplicamos el mismo bloque de cálculo
    group_modify(~{
      df <- .x # df con los datos de cada antibiótico
      # 2a) Contar TP/FN/FP/TN (sin "unknown")
      confusion_counts <- df %>%
        filter(resultados_predition_type != "unknown") %>%
        count(resultados_predition_type) %>%
        pivot_wider(
          names_from  = resultados_predition_type,
          values_from = n,
          values_fill = list(n = 0)
        )
      
      TP <- ifelse("TP" %in% names(confusion_counts), confusion_counts$TP, 0)
      TN <- ifelse("TN" %in% names(confusion_counts), confusion_counts$TN, 0)
      FP <- ifelse("FP" %in% names(confusion_counts), confusion_counts$FP, 0)
      FN <- ifelse("FN" %in% names(confusion_counts), confusion_counts$FN, 0)
      
      # 2b) Montar matriz de confusión en double para evitar overflow
      cm <- matrix(
        as.numeric(c(TP, FN,
                     FP, TN)),
        nrow = 2,
        dimnames = list(
          "Phenotype"  = c("R", "S"),
          "Prediction" = c("R", "S")
        )
      )
      storage.mode(cm) <- "double"
      
      # 2c) Tibble con los conteos crudos
      counts <- tibble(
        statistic = c("TP", "TN", "FP", "FN"),
        est       = c(TP, TN, FP, FN),
        lower     = NA_real_,
        upper     = NA_real_
      )
      
      # 2d) Métricas diagnósticas
      diag_results <- summary(epi.tests(as.table(cm), digits = 2))
      
      # 2e) Unir conteos + métricas para este antibiótico
      bind_rows(counts, diag_results)
    }) %>%
    ungroup()
}

# Uso:
diagnosticos_por_ab_resf_eucast_collapsed <- compute_diagnostics_by_ab(resultados_detallados_resfinder_eucast.ecoli_collapsed)

diagnosticos_por_ab_resf_ecoff_collapsed <- compute_diagnostics_by_ab(resultados_detallados_resfinder_ecoff.ecoli_collapsed)

diagnosticos_por_ab_resf_clsi_collapsed <- compute_diagnostics_by_ab(resultados_detallados_resfinder_clsi.ecoli_collapsed)


# procesamiento de los df cruzados 

process_diagnostics <- function(diagnostics_df) {
  # 1. Wide numérico
  wide <- diagnostics_df %>%
    group_by(resultados_antibiotic) %>%
    pivot_wider(
      id_cols     = resultados_antibiotic,
      names_from  = statistic,
      values_from = est
    ) %>%
    select(resultados_antibiotic, se, sp, diag.ac, pv.pos, pv.neg, TP, TN, FP, FN) %>%
    mutate(
      total_score = se + sp + diag.ac
    ) %>%
    arrange(desc(total_score))
  
  # 2. Formatear porcentaje en est, limpiar bounds
  porcentaje <- diagnostics_df %>%
    mutate(
      est = case_when(
        statistic %in% c("ap","tp","se","sp","pv.pos","pv.neg") ~
          sprintf("%.1f (%.1f-%.1f)",
                  est  * 100,
                  lower * 100,
                  upper * 100),
        TRUE ~ as.character(est)
      ),
      lower = if_else(statistic %in% c("ap","tp","se","sp","pv.pos","pv.neg"),
                      NA_real_, lower),
      upper = if_else(statistic %in% c("ap","tp","se","sp","pv.pos","pv.neg"),
                      NA_real_, upper)
    )
  
  # 3. Wide con porcentaje
  wide_porcentaje <- porcentaje %>%
    group_by(resultados_antibiotic) %>%
    pivot_wider(
      id_cols     = resultados_antibiotic,
      names_from  = statistic,
      values_from = est
    ) %>%
    select(resultados_antibiotic, TP, TN, FP, FN, se, sp, pv.pos, pv.neg, ) %>%
    rename(
      sensibilidad = se,
      especificidad = sp,
      PPV = pv.pos,
      NPV = pv.neg
    )
  
  # 4. Gráfico se vs sp
  plot <- wide %>%
    select(resultados_antibiotic, se, sp) %>%
    pivot_longer(
      cols      = c(se, sp),
      names_to  = "métrica",
      values_to = "valor"
    ) %>%
    ggplot(aes(x = resultados_antibiotic, y = valor, fill = métrica)) +
    geom_col(position = "dodge") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1, suffix = "")) +
    labs(
      x     = "Antibiótico",
      y     = "Porcentaje (%)",
      fill  = "Métrica",
      title = "Sensibilidad vs Especificidad por antibiótico"
    ) +
    theme_minimal() +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 16),
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_text(margin = margin(r = 10), size = 12), 
      axis.title.x = element_text(vjust = 0.5, hjust = 0.5, size = 12),
      legend.title   = element_text(size = 10),    # tamaño del título de la leyenda
      legend.text    = element_text(size = 8),    # tamaño del texto de la leyenda
      legend.key.size= unit(0.8, "lines"),
      plot.margin     = margin(t = 5, r = 2, b = 5, l = 5, unit = "mm")
    )
  
  # Devolver todo en una lista
  list(
    wide               = wide,
    porcentaje         = porcentaje,
    wide_porcentaje    = wide_porcentaje,
    plot               = plot
  )
}

# eucast
resfinder_resultados_finales_euc <- process_diagnostics(diagnosticos_por_ab_resf_eucast_collapsed)

tabla_resf_euc.dbl <- resfinder_resultados_finales_euc$wide

tabla_resf_eucast.temp <- resfinder_resultados_finales_euc$wide_porcentaje

plot_resf_eucast <- resfinder_resultados_finales_euc$plot

head(resultados_detallados_resfinder_eucast.ecoli_collapsed)

recuento_resf_eucast <- resultados_detallados_resfinder_eucast.ecoli_collapsed %>% 
  count(resultados_antibiotic, resultados_sir) %>%
  rename(EUCAST = "n") %>%
  group_by(resultados_antibiotic) %>% 
  summarise(R_eucast = sum(EUCAST[resultados_sir == "R"], na.rm = TRUE),
            S_eucast = sum(EUCAST[resultados_sir == "S"], na.rm = TRUE),
            total_eucast = sum(EUCAST, na.rm = TRUE),
            .groups = "drop"
  ) %>%
  select(resultados_antibiotic, R_eucast, S_eucast, total_eucast) %>%
  mutate(
    aislados_resistentes = paste0(
      R_eucast,
      " (",
      scales::percent(R_eucast / total_eucast, accuracy = 0.1),
      ")"
    ),
    aislados_susceptibles = paste0(
      S_eucast,
      " (",
      scales::percent(S_eucast / total_eucast, accuracy = 0.1),
      ")"
    )
  )

tabla_resf_eucast <- full_join(tabla_resf_eucast.temp, recuento_resf_eucast, by = "resultados_antibiotic") %>%
  select(resultados_antibiotic, aislados_resistentes, aislados_susceptibles, TP:NPV)

write.csv(tabla_resf_eucast, "C:/Users/rferragud/Documents/Projecte/Tables/tabla_resultados_resf_eucast3.csv", row.names = FALSE)

# ecoff

diagnosticos_por_ab_resf_ecoff_collapsed


resfinder_resultados_finales_ecoff <- process_diagnostics(diagnosticos_por_ab_resf_ecoff_collapsed)

tabla_resf_ecoff.dbl <- resfinder_resultados_finales_ecoff$wide

tabla_resf_ecoff.temp <- resfinder_resultados_finales_ecoff$wide_porcentaje

plot_resf_ecoff <- resfinder_resultados_finales_ecoff$plot



recuento_resf_ecoff <- resultados_detallados_resfinder_ecoff.ecoli_collapsed %>% 
  count(resultados_antibiotic, resultados_sir) %>%
  rename(ECOFF = "n") %>%
  group_by(resultados_antibiotic) %>% 
  summarise(R_ecoff = sum(ECOFF[resultados_sir == "R"], na.rm = TRUE),
            S_ecoff = sum(ECOFF[resultados_sir == "S"], na.rm = TRUE),
            total_ecoff = sum(ECOFF, na.rm = TRUE),
            .groups = "drop"
  ) %>%
  select(resultados_antibiotic, R_ecoff, S_ecoff, total_ecoff) %>%
  mutate(
    aislados_resistentes = paste0(
      R_ecoff,
      " (",
      scales::percent(R_ecoff / total_ecoff, accuracy = 0.1),
      ")"
    ),
    aislados_susceptibles = paste0(
      S_ecoff,
      " (",
      scales::percent(S_ecoff / total_ecoff, accuracy = 0.1),
      ")"
    )
  )

tabla_resf_ecoff <- full_join(tabla_resf_ecoff.temp, recuento_resf_ecoff, by = "resultados_antibiotic") %>%
  select(resultados_antibiotic, aislados_resistentes, aislados_susceptibles, TP:NPV)

write.csv(tabla_resf_ecoff, "C:/Users/rferragud/Documents/Projecte/Tables/tabla_resultados_resf_ecoff3.csv", row.names = FALSE)


# clsi
diagnosticos_por_ab_resf_clsi_collapsed

resfinder_resultados_finales_clsi <- process_diagnostics(diagnosticos_por_ab_resf_clsi_collapsed)

tabla_resf_clsi.dbl <- resfinder_resultados_finales_clsi$wide

tabla_resf_clsi.temp <- resfinder_resultados_finales_clsi$wide_porcentaje

plot_resf_clsi <- resfinder_resultados_finales_clsi$plot



recuento_resf_clsi <- resultados_detallados_resfinder_clsi.ecoli_collapsed %>% 
  count(resultados_antibiotic, resultados_sir) %>%
  rename(CLSI = "n") %>%
  group_by(resultados_antibiotic) %>% 
  summarise(R_clsi = sum(CLSI[resultados_sir == "R"], na.rm = TRUE),
            S_clsi = sum(CLSI[resultados_sir == "S"], na.rm = TRUE),
            total_clsi = sum(CLSI, na.rm = TRUE),
            .groups = "drop"
  ) %>%
  select(resultados_antibiotic, R_clsi, S_clsi, total_clsi) %>%
  mutate(
    aislados_resistentes = paste0(
      R_clsi,
      " (",
      scales::percent(R_clsi / total_clsi, accuracy = 0.1),
      ")"
    ),
    aislados_susceptibles = paste0(
      S_clsi,
      " (",
      scales::percent(S_clsi / total_clsi, accuracy = 0.1),
      ")"
    )
  )

tabla_resf_clsi <- full_join(tabla_resf_clsi.temp, recuento_resf_clsi, by = "resultados_antibiotic") %>%
  select(resultados_antibiotic, aislados_resistentes, aislados_susceptibles, TP:NPV)

write.csv(tabla_resf_clsi, "C:/Users/rferragud/Documents/Projecte/Tables/tabla_resultados_resf_clsi3.csv", row.names = FALSE)


tabla_resf_eucast2 <- tabla_resf_euc.dbl %>%
  rename_with(~ paste0(.x, "_EUCAST")) %>% 
  select(resultados_antibiotic_EUCAST, se_EUCAST, sp_EUCAST) %>%
  rename(antibiotic = resultados_antibiotic_EUCAST)

tabla_resf_ecoff2 <- tabla_resf_ecoff.dbl %>%
  rename_with(~ paste0(.x, "_ECOFF")) %>% 
  select(resultados_antibiotic_ECOFF, se_ECOFF, sp_ECOFF) %>%
  rename(antibiotic = resultados_antibiotic_ECOFF)

tabla_resf_clsi2 <- tabla_resf_clsi.dbl %>%
  rename_with(~ paste0(.x, "_CLSI")) %>% 
  select(resultados_antibiotic_CLSI, se_CLSI, sp_CLSI) %>%
  rename(antibiotic = resultados_antibiotic_CLSI)

tabla_resf_completa <- tabla_resf_eucast2 %>% 
  full_join(tabla_resf_ecoff2, by = "antibiotic") %>%
  full_join(tabla_resf_clsi2, by = "antibiotic")

plot_se_sp_resf <- tabla_resf_completa %>%
  pivot_longer(
    cols = -antibiotic,
    names_to  = c("metric", "criterion"),
    names_sep = "_"
  ) %>% ggplot(aes(x = antibiotic, y = value, fill = criterion)) +
  geom_col(position = "dodge") +
  facet_wrap(
    ~ metric,
    nrow     = 2,
    labeller = labeller(
      metric = c(
        se = "Sensibilidad",
        sp = "Especificidad"
      )
    )
  ) +
  scale_fill_manual(
    values = c("se" = "#F8766D", "sp" = "#00BFC4"),
    labels = c(se = "Sensibilidad", sp = "Especificidad")
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1, suffix = "")) +
  labs(
    x     = "Antibiótico",
    y     = "Porcentaje (%)",
    fill  = "Criterio",
    title = "Sensibilidad vs Especificidad"
  ) +
  theme_minimal() +
  theme(
    plot.title   = element_text(hjust = 0.5),
    strip.text.x       = element_text(face = "bold", size = 12),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    legend.position    = "right", 
    legend.title   = element_text(size = 8),    # tamaño del título de la leyenda
    legend.text    = element_text(size = 6),    # tamaño del texto de la leyenda
    legend.key.size= unit(0.4, "lines"),    
  )




plot_resf_euc_clsi_ecoff <- tabla_resf_completa %>%
  pivot_longer(
    cols      = -antibiotic,
    names_to  = c("metric", "criterion"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(
    criterion = fct_relevel(criterion,
                            "EUCAST",  # primero EUCAST
                            "CLSI",    # luego CLSI
                            "ECOFF")   # luego ECOFF
  ) %>%
  ggplot(aes(x = antibiotic, y = value, fill = metric)) +
  geom_col(position = "dodge") +
  facet_wrap(
    ~ criterion,
    ncol     = 1,
    labeller = labeller(
      eucast = "EUCAST",
      ecoff  = "ECOFF",
      clsi   = "CLSI"
    )
  ) +
  scale_fill_manual(
    values = c("se" = "#F8766D", "sp" = "#00BFC4"),
    labels = c(se = "Sensibilidad", sp = "Especificidad")
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1, suffix = "")
  ) +
  labs(
    x     = "Antibiótico",
    y     = "Porcentaje (%)",
    fill  = "Métrica",
    title = "Sensibilidad vs Especificidad"
  ) +
  theme_minimal() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(r = 10), size = 12), 
    axis.title.x = element_text(vjust = 0.5, hjust = 0.5, size = 12),
    legend.title   = element_text(size = 10),    # tamaño del título de la leyenda
    legend.text    = element_text(size = 8),    # tamaño del texto de la leyenda
    legend.key.size= unit(0.8, "lines"),
    plot.margin     = margin(t = 5, r = 5, b = 5, l = 5, unit = "mm")
  )











