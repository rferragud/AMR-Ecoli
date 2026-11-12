
library(tidyverse)
library(data.table)
library(dplyr)
library(readODS)
library(epiR)
library(scales)

# AMRFINDERPLUS

analizar_antibiotico_amrf <- function(antibiotic, program_class, sir_df, program) {
  
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
      by = "ENA_sample_accession") %>%
    mutate(
      predition_type = case_when(
        sir == "R" & !is.na(genetic_background) ~ "TP",
        sir == "S" & is.na(genetic_background)  ~ "TN",
        sir == "S" & !is.na(genetic_background) ~ "FP",
        sir == "R" & is.na(genetic_background)  ~ "FN"
      ),
      antibiotic = antibiotic  # asegura que quede como columna en cada fila
    ) 
  
  return(cruzado)
}


corresp_atb_amrf <- read_ods("C:/Users/rferragud/Documents/Projecte/corresp_antibiotic_resf_amrf.ods", sheet = 2)

# eucast
corresp_atb_euc_amrf <- corresp_atb_amrf %>% filter(antibiotic %in% atb_eucast_names)

# ecoff
corresp_atb_ecoff_amrf <- corresp_atb_amrf %>% filter(antibiotic %in% atb_ecoff_names)

# clsi
corresp_atb_clsi_amrf <- corresp_atb_amrf %>% filter(antibiotic %in% atb_clsi_names)


EUCAST_sir_clean2
ECOFF_sir_clean2
CLSI_sir_clean2

amrf_results_accession_cruzado.completo

amrf_collapsed <- amrf_results_accession_cruzado.completo %>%
  separate_rows(Subclass, sep = "/") %>%      # separa en filas donde haya "/"
  mutate(subclass = str_trim(Subclass)) %>%      # quita espacios sobrantes
  group_by(
    ENA_sample_accession,
    ENA_run_accession,
    assembly,
    Subclass, 
    #Class, 
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

write.csv(amrf_collapsed, "C:/Users/rferragud/Documents/Projecte/Tables/amrf_collapsed3.csv", row.names = FALSE)

# resultados amrf
# a) eucast
resultados_detallados_amrf_eucast.ecoli_collapsed <- corresp_atb_euc_amrf %>%
  mutate(
    resultados = map2(
      antibiotic,
      program_class,
      analizar_antibiotico_amrf,
      sir_df = EUCAST_sir_clean2, 
      program = amrf_collapsed
    )
  ) %>%
  unnest(resultados, names_sep = "_")


# b) ecoff
resultados_detallados_amrf_ecoff.ecoli_collapsed <- corresp_atb_ecoff_amrf %>%
  mutate(
    resultados = map2(
      antibiotic,
      program_class,
      analizar_antibiotico_amrf,
      sir_df = ECOFF_sir_clean2, 
      program = amrf_collapsed
    )
  ) %>%
  unnest(resultados, names_sep = "_")

# c) clsi
resultados_detallados_amrf_clsi.ecoli_collapsed <- corresp_atb_clsi_amrf %>%
  mutate(
    resultados = map2(
      antibiotic,
      program_class,
      analizar_antibiotico_amrf,
      sir_df = CLSI_sir_clean2, 
      program = amrf_collapsed
    )
  ) %>%
  unnest(resultados, names_sep = "_")


#---------------------------------- GRÁFICOS ----------------------------------#

# AMRFINDERPLUS 

# EUCAST
resultados_detallados_amrf_eucast.ecoli_collapsed

graf_amrf_euc_eco <- resultados_detallados_amrf_eucast.ecoli_collapsed %>% 
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

resultados_detallados_amrf_ecoff.ecoli_collapsed

graf_amrf_ecoff_eco <- resultados_detallados_amrf_ecoff.ecoli_collapsed %>% 
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

resultados_detallados_amrf_clsi.ecoli_collapsed

graf_amrf_clsi_eco <- resultados_detallados_amrf_clsi.ecoli_collapsed %>% 
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
diagnosticos_por_ab_amrf_eucast_collapsed <- compute_diagnostics_by_ab(resultados_detallados_amrf_eucast.ecoli_collapsed)

diagnosticos_por_ab_amrf_ecoff_collapsed <- compute_diagnostics_by_ab(resultados_detallados_amrf_ecoff.ecoli_collapsed)

diagnosticos_por_ab_amrf_clsi_collapsed <- compute_diagnostics_by_ab(resultados_detallados_amrf_clsi.ecoli_collapsed)


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
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_text(margin = margin(r = 10), size = 12), 
      axis.title.x = element_text(vjust = 0.5, hjust = 0.5, size = 12),
      legend.title   = element_text(size = 10),    # tamaño del título de la leyenda
      legend.text    = element_text(size = 8),    # tamaño del texto de la leyenda
      legend.key.size= unit(0.8, "lines"),
      plot.margin     = margin(t = 5, r = 5, b = 5, l = 5, unit = "mm")
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
amrf_resultados_finales_euc <- process_diagnostics(diagnosticos_por_ab_amrf_eucast_collapsed)

tabla_amrf_euc.dbl <- amrf_resultados_finales_euc$wide

tabla_amrf_eucast.temp <- amrf_resultados_finales_euc$wide_porcentaje

plot_amrf_eucast <- amrf_resultados_finales_euc$plot

head(resultados_detallados_amrf_eucast.ecoli_collapsed)

recuento_amrf_eucast <- resultados_detallados_amrf_eucast.ecoli_collapsed %>% 
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

tabla_amrf_eucast <- full_join(tabla_amrf_eucast.temp, recuento_amrf_eucast, by = "resultados_antibiotic") %>%
  select(resultados_antibiotic, aislados_resistentes, aislados_susceptibles, TP:NPV)

write.csv(tabla_amrf_eucast, "C:/Users/rferragud/Documents/Projecte/Tables/tabla_resultados_amrf_eucast3.csv", row.names = FALSE)


# ecoff

diagnosticos_por_ab_amrf_ecoff_collapsed


amrf_resultados_finales_ecoff <- process_diagnostics(diagnosticos_por_ab_amrf_ecoff_collapsed)

tabla_amrf_ecoff.dbl <- amrf_resultados_finales_ecoff$wide

tabla_amrf_ecoff.temp <- amrf_resultados_finales_ecoff$wide_porcentaje

plot_amrf_ecoff <- amrf_resultados_finales_ecoff$plot



recuento_amrf_ecoff <- resultados_detallados_amrf_ecoff.ecoli_collapsed %>% 
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

tabla_amrf_ecoff <- full_join(tabla_amrf_ecoff.temp, recuento_amrf_ecoff, by = "resultados_antibiotic") %>%
  select(resultados_antibiotic, aislados_resistentes, aislados_susceptibles, TP:NPV)

write.csv(tabla_amrf_ecoff, "C:/Users/rferragud/Documents/Projecte/Tables/tabla_resultados_amrf_ecoff3.csv", row.names = FALSE)


# clsi
diagnosticos_por_ab_amrf_clsi_collapsed

amrf_resultados_finales_clsi <- process_diagnostics(diagnosticos_por_ab_amrf_clsi_collapsed)

tabla_amrf_clsi.dbl <- amrf_resultados_finales_clsi$wide 

tabla_amrf_clsi.dbl%>% print(n = Inf)

tabla_amrf_clsi.temp <- amrf_resultados_finales_clsi$wide_porcentaje

plot_amrf_clsi <- amrf_resultados_finales_clsi$plot



recuento_amrf_clsi <- resultados_detallados_amrf_clsi.ecoli_collapsed %>% 
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

tabla_amrf_clsi <- full_join(tabla_amrf_clsi.temp, recuento_amrf_clsi, by = "resultados_antibiotic") %>%
  select(resultados_antibiotic, aislados_resistentes, aislados_susceptibles, TP:NPV)

write.csv(tabla_amrf_clsi, "C:/Users/rferragud/Documents/Projecte/Tables/tabla_resultados_amrf_clsi3.csv", row.names = FALSE)



tabla_amrf_eucast2 <- tabla_amrf_euc.dbl %>%
  rename_with(~ paste0(.x, "_EUCAST")) %>% 
  select(resultados_antibiotic_EUCAST, se_EUCAST, sp_EUCAST) %>%
  rename(antibiotic = resultados_antibiotic_EUCAST)

tabla_amrf_ecoff2 <- tabla_amrf_ecoff.dbl %>%
  rename_with(~ paste0(.x, "_ECOFF")) %>% 
  select(resultados_antibiotic_ECOFF, se_ECOFF, sp_ECOFF) %>%
  rename(antibiotic = resultados_antibiotic_ECOFF)

tabla_amrf_clsi2 <- tabla_amrf_clsi.dbl %>%
  rename_with(~ paste0(.x, "_CLSI")) %>% 
  select(resultados_antibiotic_CLSI, se_CLSI, sp_CLSI) %>%
  rename(antibiotic = resultados_antibiotic_CLSI)

tabla_amrf_completa <- tabla_amrf_eucast2 %>% 
  full_join(tabla_amrf_ecoff2, by = "antibiotic") %>%
  full_join(tabla_amrf_clsi2, by = "antibiotic")

write.csv(tabla_amrf_completa, "C:/Users/rferragud/Documents/Projecte/Tables/tabla_amrf_completa3.csv", row.names = FALSE)


plot_se_sp_amrf <- tabla_amrf_completa %>%
  pivot_longer(
    cols = -antibiotic,
    names_to  = c("metric", "criterion"),
    names_sep = "_"
  ) %>% ggplot(aes(x = antibiotic, y = value, fill = criterion)) +
  geom_col(position = "dodge",) +
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
    legend.key.size= unit(0.4, "lines")    
  )




plot_amrf_euc_clsi_ecoff <- tabla_amrf_completa %>%
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




tabla_amrf_completa

tabla_resf_completa



library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# AMRF vs ResF (EUCAST)

# 1) Prepara y etiqueta cada tabla
amrf <- tabla_amrf_completa %>% filter(!is.na(se_EUCAST)|!is.na(sp_EUCAST)) %>%
  select(antibiotic, se = se_EUCAST, sp = sp_EUCAST) %>%
  mutate(method = "AMRFPlus") 

resf <- tabla_resf_completa %>% filter(!is.na(se_EUCAST)|!is.na(sp_EUCAST)) %>%
  select(antibiotic, se = se_EUCAST, sp = sp_EUCAST) %>%
  mutate(method = "ResF") 

atb_nuevos <- anti_join(amrf, resf, by = "antibiotic")%>%
  mutate(method = "ResF") %>%
  mutate(se = 0 , sp = 0)

# 2) Combina y pivota a largo
df_plot <- bind_rows(amrf, resf, atb_nuevos) %>%
  pivot_longer(
    cols      = c(se, sp),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric, levels = c("se", "sp"),
                    labels = c("Sensibilidad (%)", "Especificidad (%)")),
    method = factor(method, levels = c("AMRFPlus", "ResF")) 
  ) #%>%
 #group_by(antibiotic) %>%
 ## Sólo mantenemos los antibióticos con los dos métodos
 #filter(n_distinct(method) == 2) %>%
 #ungroup()%>%
 #filter(!is.na(value))



# 3) Dibuja el gráfico
ggplot(df_plot, aes(x = antibiotic, y = value, fill = method)) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_grid(
    metric ~ ., 
    switch = "y",     # mueve el strip al eje Y
    scales = "fixed",
    labeller = labeller(
      Sensibilidad = "Sensibilidad(%)",
      Especificidad = "Especificidad (%)"
    )
  )+
  scale_y_continuous(
    labels      = percent_format(accuracy = 1, suffix = ""),
    expand      = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(
    values = c(AMRFPlus = "steelblue", ResF = "orange"),
    name   = "Base de datos", 
    labels = c(AMRFPlus = "AMRFinderPlus", ResF = "ResFinder")
  ) +
  labs(
    x     = "Antibiótico",
    y = NULL,
    title = "Sensibilidad y Especificidad EUCAST: AMRFinderPlus vs ResFinder"
  ) +
  theme_minimal() +
  theme(
    strip.placement     = "outside",
    strip.text.y.left   = element_text(angle = 90, vjust = 0.5, size = 10),
    plot.title   = element_text(hjust = 0.7, size = 12),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(r = 10)), 
    axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
    legend.title   = element_text(size = 10),    # tamaño del título de la leyenda
    legend.text    = element_text(size = 8),    # tamaño del texto de la leyenda
    legend.key.size= unit(0.8, "lines"),
    plot.margin     = margin(t = 5, r = 10, b = 5, l = 10, unit = "mm")
   #strip.text.y.left   = element_text(angle = 90, face = "bold", vjust = 0.5),
   #plot.title      = element_text(hjust = 0.8, face = "bold", size = 11),
   #axis.text.x     = element_text(angle = 45, hjust = 1),
   #strip.text      = element_text(face = "bold"),
   #legend.position = "top", 
   #legend.title   = element_text(size = 10),    # tamaño del título de la leyenda
   #legend.text    = element_text(size = 8),    # tamaño del texto de la leyenda
   #legend.key.size= unit(0.6, "lines"), 
   #plot.margin     = margin(t = 5, r = 10, b = 5, l = 10, unit = "mm")
    
  )

# AMRF vs ResF (ECOFF)

# 1) Prepara y etiqueta cada tabla
amrf <- tabla_amrf_completa %>% filter(!is.na(se_ECOFF)|!is.na(sp_ECOFF)) %>%
  select(antibiotic, se = se_ECOFF, sp = sp_ECOFF) %>%
  mutate(method = "AMRFPlus") 

resf <- tabla_resf_completa %>% filter(!is.na(se_ECOFF)|!is.na(sp_ECOFF)) %>%
  select(antibiotic, se = se_ECOFF, sp = sp_ECOFF) %>%
  mutate(method = "ResF")

atb_nuevos <- anti_join(amrf, resf, by = "antibiotic")%>%
  mutate(method = "ResF") %>%
  mutate(se = 0 , sp = 0)


# 2) Combina y pivota a largo
df_plot <- bind_rows(amrf, resf, atb_nuevos) %>%
  pivot_longer(
    cols      = c(se, sp),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric, levels = c("se", "sp"),
                    labels = c("Sensibilidad (%)", "Especificidad (%)")),
    method = factor(method, levels = c("AMRFPlus", "ResF"))
  ) #%>%
  #group_by(antibiotic) %>%
  #  filter()
  ## Sólo mantenemos los antibióticos con los dos métodos
  #filter(n_distinct(method) == 2) %>%
  #ungroup()
  
# TODO llevar antes els NA de AMRF pq no tenen tampoc de resF 

# 3) Dibuja el gráfico
ggplot(df_plot, aes(x = antibiotic, y = value, fill = method)) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_grid(
    metric ~ ., 
    switch = "y",     # mueve el strip al eje Y
    scales = "fixed",
    labeller = labeller(
      se  = "Sensibilidad (%)",
      sp = "Especificidad (%)"
    )
  )+
  scale_y_continuous(
    labels      = percent_format(accuracy = 1, suffix = ""),
    expand      = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(
    values = c(AMRFPlus = "steelblue", ResF = "orange"),
    name   = "Base de datos", 
    labels = c(AMRFPlus = "AMRFinderPlus", ResF = "ResFinder")
  ) +
  labs(
    x     = "Antibiótico",
    y = NULL,
    title = "Sensibilidad y Especificidad ECOFF: AMRFinderPlus vs ResFinder"
  ) +
  theme_minimal() +
  theme(
    strip.placement     = "outside",
    strip.text.y.left   = element_text(angle = 90, vjust = 0.5, size = 10),
    plot.title   = element_text(hjust = 0.7, size = 12),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(r = 10)), 
    axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
    legend.title   = element_text(size = 10),    # tamaño del título de la leyenda
    legend.text    = element_text(size = 8),    # tamaño del texto de la leyenda
    legend.key.size= unit(0.8, "lines"),
    plot.margin     = margin(t = 5, r = 10, b = 5, l = 10, unit = "mm")
   #strip.placement     = "outside",
   #strip.text.y.left   = element_text(angle = 90, face = "bold", vjust = 0.5),
   #plot.title      = element_text(hjust = 0.8, face = "bold", size = 11),
   #axis.text.x     = element_text(angle = 45, hjust = 1),
   #strip.text      = element_text(face = "bold"),
   #legend.position = "top", 
   #legend.title   = element_text(size = 10),    # tamaño del título de la leyenda
   #legend.text    = element_text(size = 8),    # tamaño del texto de la leyenda
   #legend.key.size= unit(0.6, "lines"), 
   #plot.margin     = margin(t = 5, r = 10, b = 5, l = 10, unit = "mm")
    
  )

# amrf vs resf (CLSI)

# 1) Prepara y etiqueta cada tabla
amrf <- tabla_amrf_completa %>%  filter(!is.na(se_CLSI)|!is.na(sp_CLSI)) %>%
  select(antibiotic, se = se_CLSI, sp = sp_CLSI) %>%
  mutate(method = "AMRFPlus") 

resf <- tabla_resf_completa %>% filter(!is.na(se_CLSI)|!is.na(sp_CLSI)) %>%
  select(antibiotic, se = se_CLSI, sp = sp_CLSI) %>%
  mutate(method = "ResF")

atb_nuevos <- anti_join(amrf, resf, by = "antibiotic")%>%
  mutate(method = "ResF") %>%
  mutate(se = 0 , sp = 0)

# 2) Combina y pivota a largo
df_plot <- bind_rows(amrf, resf, atb_nuevos) %>%
  pivot_longer(
    cols      = c(se, sp),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric, levels = c("se", "sp"),
                    labels = c("Sensibilidad (%)", "Especificidad (%)")),
    method = factor(method, levels = c("AMRFPlus", "ResF"))
  ) # %>%
  #group_by(antibiotic) %>%
  ## Sólo mantenemos los antibióticos con los dos métodos
  #filter(n_distinct(method) == 2) %>%
  #ungroup()


# 3) Dibuja el gráfico
ggplot(df_plot, aes(x = antibiotic, y = value, fill = method)) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_grid(
    metric ~ ., 
    switch = "y",     # mueve el strip al eje Y
    scales = "fixed",
    labeller = labeller(
      se  = "Sensibilidad (%)",
      sp = "Especificidad (%)"
    )
  )+
  scale_y_continuous(
    labels      = percent_format(accuracy = 1, suffix = ""),
    expand      = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(
    values = c(AMRFPlus = "steelblue", ResF = "orange"),
    name   = "Base de datos", 
    labels = c(AMRFPlus = "AMRFinderPlus", ResF = "ResFinder")
  ) +
  labs(
    x     = "Antibiótico",
    y = NULL,
    title = "Sensibilidad y Especificidad CLSI: AMRFinderPlus vs ResFinder"
  ) +
  theme_minimal() +
  theme(
    strip.placement     = "outside",
    strip.text.y.left   = element_text(angle = 90, vjust = 0.5, size = 10),
    plot.title   = element_text(hjust = 0.7, size = 12),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(r = 10)), 
    axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
    legend.title   = element_text(size = 10),    # tamaño del título de la leyenda
    legend.text    = element_text(size = 8),    # tamaño del texto de la leyenda
    legend.key.size= unit(0.8, "lines"),
    plot.margin     = margin(t = 5, r = 10, b = 5, l = 10, unit = "mm")
   #strip.placement     = "outside",
   #strip.text.y.left   = element_text(angle = 90, face = "bold", vjust = 0.5),
   #plot.title      = element_text(hjust = 0.8, face = "bold", size = 11),
   #axis.text.x     = element_text(angle = 45, hjust = 1),
   #strip.text      = element_text(face = "bold"),
   #legend.position = "top", 
   #legend.title   = element_text(size = 10),    # tamaño del título de la leyenda
   #legend.text    = element_text(size = 8),    # tamaño del texto de la leyenda
   #legend.key.size= unit(0.6, "lines"), 
   #plot.margin     = margin(t = 5, r = 10, b = 5, l = 10, unit = "mm")
    
  )

