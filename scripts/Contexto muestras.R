
library(countrycode)

# Contextualización muestras 

df_cruzado_completo <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/df_cruzado_completo_articles.csv", check.names = FALSE)
patric_metadata_st.nuevas <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/patric_metadata_st.nuevas.csv", check.names = F)
muestras_NCBI.nuevas.coli.st <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/muestras_NCBI.nuevas.coli.st.csv", check.names = F)


NCBI <- muestras_NCBI.nuevas.coli.st %>% select(ENA_sample_accession:location)%>% distinct() 
articles <- df_cruzado_completo %>% select(ENA_experiment_accession:host_tissue_sampled, ENA_sample_accession)%>% distinct()
patric <- patric_metadata_st.nuevas%>% select(isolation_source:ENA_run_accession)%>% distinct()

muestras_total <- bind_rows(NCBI, articles, patric)%>%
  mutate(
    # extraemos antes de ":"  
    country_from_loc = sub(":.*$", "", location),
    # sustituimos solo donde geographic_location_country es NA o cadena vacía
    geographic_location_country = if_else(
      is.na(geographic_location_country) | geographic_location_country == "",
      country_from_loc,
      geographic_location_country
    )
  ) %>%
  select(-country_from_loc)

# seleccionar las que pasan el QC y tienen MIC

muestras_QC_MIC <- amr_resf_included_ns2 %>% filter(QC == "included") %>% distinct(resultados_ENA_sample_accession) %>% pull(resultados_ENA_sample_accession)

all_summary_def_qc_clean_accession_mic <- all_summary_def_qc_clean_accession %>% filter(ENA_sample_accession %in% muestras_QC_MIC)

ENA_sam_qc <- all_summary_def_qc_clean_accession_mic %>% filter(!is.na(ENA_sample_accession))%>% distinct(ENA_sample_accession) %>% pull()
ENA_run_qc <- all_summary_def_qc_clean_accession_mic %>% filter(!is.na(ENA_run_accession)) %>% distinct(ENA_run_accession) %>% pull()

muestras_buenas <- muestras_total %>% filter(ENA_sample_accession %in% ENA_sam_qc | ENA_run_accession %in% ENA_run_qc | ENA_secondary_sample_accession == "ERS1995029" | ENA_secondary_sample_accession == "ERS1995244" ) %>% distinct()

muestras_buenas

por_origen <- group_by(muestras_buenas, geographic_location_country, geographic_location_region)

origen_geografico <- summarise(por_origen, conteo =n()) 

por_pais <- group_by(muestras_buenas, geographic_location_country)

por_pais <- por_pais %>%
  mutate(
    geographic_location_country = str_replace(
      geographic_location_country, 
      "^England$", 
      "United Kingdom"
    )
  ) %>% group_by(geographic_location_country)

origen_pais <- summarise(por_pais, conteo =n()) 

origen_pais.st <- origen_pais %>%
  mutate(pais_estandarizado = countrycode(geographic_location_country, 
                                          origin = "country.name", 
                                          destination = "country.name"))

resumen_origen_geo <- origen_pais %>%
  mutate(
    continent = countrycode(
      sourcevar    = geographic_location_country,            # tu columna de país
      origin       = "country.name",     # formato de nombre de país
      destination  = "continent"         # queremos el continente
    )
  )

write.csv(resumen_origen_geo, "C:/Users/rferragud/Documents/Projecte/Tables/resumen_paises.csv", row.names = FALSE) 




df_cont <- resumen_origen_geo %>%
  mutate(
    continent = ifelse(is.na(continent), "Origen desconocido", continent),
    continent = recode(continent,
                       "Africa"   = "África",
                       "Europe"   = "Europa",
                       "Americas" = "América",
                       .default   = continent
    )
  ) %>%
  group_by(continent) %>%
  summarise(total = sum(conteo, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    prop  = total / sum(total) * 100,
    label = (paste0(round(prop, 1), "%"))
  
  )

# 3. Crea el pie chart
ggplot(df_cont, aes(x = "", y = total, fill = continent)) +
  geom_col(color = "white", width = 1) +
  coord_polar(theta = "y", start = 0) +
  
  # Sólo etiquetas para continentes distintos de Asia y origen desconocido
  geom_text(
    data = filter(df_cont, !continent %in% c("África", "Origen desconocido")),
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  
  # Títulos y leyenda
  labs(
    title = "Distribución de aislados por continente",
    fill  = "Continente"
  ) +
  
  # Paleta fría de RColorBrewer
  scale_fill_brewer(palette = "Blues") +
  
  # Tema más “limpio”
  theme_minimal(base_size = 14) +
  theme(
    plot.title      = element_text(face = "bold",
                                   hjust = 0,
                                   size = 15),        # tamaño del título
    legend.position = "right",
    legend.title    = element_text(size = 14),       # tamaño del título de leyenda
    legend.text     = element_text(size = 12),
    panel.grid       = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank()
  )


por_host <- group_by(muestras_buenas, host)

host <- summarise(por_host, conteo = n())





por_isolation_source <- group_by(muestras_buenas, isolation_source)


isolation_source <- summarise(por_isolation_source, conteo = n()) %>% print(n = Inf)


host_clasificado <- host %>%
  mutate(
    tipo = case_when(
      is.na(host)                                                   ~ "Desconocido",
      str_detect(host, regex("Homo sapiens|Human|human",    ignore_case = TRUE)) ~ "Humano",
      str_detect(host, regex("Bos taurus|Cow",               ignore_case = TRUE)) ~ "Mamífero",
      str_detect(host, regex("Sus scrofa|Pig",               ignore_case = TRUE)) ~ "Mamífero",
      str_detect(host, regex("Equus caballus|Horse",         ignore_case = TRUE)) ~ "Mamífero",
      str_detect(host, regex("Canis lupus|Dog",              ignore_case = TRUE)) ~ "Mamífero",
      str_detect(host, regex("Felis silvestris|Cat",         ignore_case = TRUE)) ~ "Mamífero",
      # Aves de corral
      str_detect(host, regex("Gallus gallus|Chicken",        ignore_case = TRUE)) ~ "Ave",
      str_detect(host, regex("Meleagris gallopavo|Turkey",   ignore_case = TRUE)) ~ "Ave",
      # Productos de pescado
      str_detect(host, regex("Siluriformes|Catfish",   ignore_case = TRUE)) ~ "Pez",
      # Por si queda algún otro
      TRUE                                                           ~ "Otros"
    )
  ) %>%
  group_by(tipo) %>%
  summarise(total_isolados = sum(conteo, na.rm = TRUE), .groups = "drop")


host_clasificado_animales <- host %>%
  mutate(
    tipo = case_when(
      is.na(host)                                                   ~ "Huésped desconocido",
      str_detect(host, regex("Homo sapiens|Human|human",    ignore_case = TRUE)) ~ "Humano",
      str_detect(host, regex("Bos taurus|Cow",               ignore_case = TRUE)) ~ "Vaca",
      str_detect(host, regex("Sus scrofa|Pig",               ignore_case = TRUE)) ~ "Cerdo",
      str_detect(host, regex("Equus caballus|Horse",         ignore_case = TRUE)) ~ "Caballo",
      str_detect(host, regex("Canis lupus|Dog",              ignore_case = TRUE)) ~ "Perro",
      str_detect(host, regex("Felis silvestris|Cat",         ignore_case = TRUE)) ~ "Gato",
      # Aves de corral
      str_detect(host, regex("Gallus gallus|Chicken",        ignore_case = TRUE)) ~ "Gallo",
      str_detect(host, regex("Meleagris gallopavo|Turkey",   ignore_case = TRUE)) ~ "Pavo",
      # Productos de pescado
      str_detect(host, regex("Siluriformes|Catfish",   ignore_case = TRUE)) ~ "Pez",
      # Por si queda algún otro
      TRUE                                                           ~ "Otros"
    )
  ) %>%
  group_by(tipo) %>%
  summarise(total_isolados = sum(conteo, na.rm = TRUE), .groups = "drop")


ggplot(host_clasificado_animales, aes(x = "", y = total_isolados, fill = tipo)) +
  geom_col(color = "white", width = 1) +
  coord_polar(theta = "y", start = 0) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Distribución de aislados por tipo de huésped",
    fill  = "Tipo de huésped"
  ) +
  theme_void() +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "right",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 10)
  )







library(dplyr)
library(stringr)

df_iso_cat <- isolation_source %>%
  mutate(
    source = ifelse(is.na(isolation_source) | str_trim(isolation_source)== "",
                    "Desconocido",
                    str_to_lower(isolation_source)),
    categoria_iso = case_when(
      str_detect(source, "clinical|culture|blood|urine|peritoneal|aspirat|wound|tissue") ~ "Clínico (Humano)",
      str_detect(source, "human|homo sapiens")                                                     ~ "Clínico (Humano)",
      str_detect(source, "caecum|bile|pond|animal|catfish|chicken|pig|cow|horse|dog|duck|turkey") ~ "Animal",
      str_detect(source, "meat|filet|filets|ground|breast|kebab|sandwich|nugget|dressed")        ~ "Alimentos",
      str_detect(source, "water|dust|environmental|fly")                                         ~ "Ambiental",
      TRUE                                                                                        ~ "Otros"
    )
  ) %>%
  group_by(categoria_iso) %>%
  summarise(total = sum(conteo, na.rm = TRUE), .groups = "drop")

df_iso_cat


por_host_isolation <- group_by(muestras_buenas, host, isolation_source)

isolation_host <- summarise(por_host_isolation, conteo = n()) %>% print(n = Inf)




df_clas <- isolation_host %>%
  mutate(
    # 1. Normalizo isolation_source y host
    iso_src   = if_else(is.na(isolation_source) | str_trim(isolation_source) == "",
                        "desconocido",
                        str_to_lower(isolation_source)),
    host_lc   = if_else(is.na(host), "", str_to_lower(host))
  ) %>%
  mutate(
    categoria = case_when(
      # 2. Primero: todo lo que tenga huésped animal
      str_detect(host_lc, "bos taurus|equus caballus|canis lupus|felis|sus scrofa|gallus|meleagris|ictaluridae|canine") ~ "Animales",
      
      # 3. Humano clínico: muestras clínicas de humanos
      str_detect(host_lc, "homo sapiens|human") &
        str_detect(iso_src, "clinical|blood|stool|urine|culture|swab|aspirat|fluid|wound|peritoneal|abscess") ~ "Humano (Clínico)",
      
      # 4. Alimentos: carne y derivados
      str_detect(iso_src, "meat|filet|filets|ground|nugget|kebab|sandwich|drumstick|wing|patties|retail|prok")  ~ "Alimentos",
      
      # 5. Ambiental: agua, polvo, insectos, superficies…
      str_detect(iso_src, "water|dust|environmental|fly|pond|farm|soil|surface")                             ~ "Ambiental",
      
      # 6. Humano (otros): muestras humanas no clínicas
      str_detect(host_lc, "homo sapiens|human")                                                           ~ "Humano (Clínico)",
      
      # 7. Desconocido
      iso_src == "desconocido"                                                                            ~ "Desconocido",
      
      # 8. Resto
      TRUE                                                                                  ~ "Otro"
    )
  ) %>%
  # 4. Agrupar y sumar
  group_by(categoria) %>%
  summarise(total = sum(conteo, na.rm = TRUE), .groups = "drop")

write.csv(df_clas, "C:/Users/rferragud/Documents/Projecte/Tables/origen_aislamiento.csv", row.names = FALSE) 


df_clas

df_plot <- df_clas %>%
  mutate(
    prop  = total / sum(total) * 100,
    label = if_else(categoria == "Desconocido", "", paste0(round(prop, 1), "%"))
  )



ggplot(df_plot, aes(x = "", y = total, fill = categoria)) +
  geom_col(color = "white", width = 1) +     # width = 1 para eliminar el hueco
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 3) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    title = "Distribución de aislados por huésped",
    fill  = "Huésped"
  ) +
  theme_void() +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "right",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 10)
  )

library(ggplot2)
library(RColorBrewer)

ggplot(df_cont, aes(x = "", y = total, fill = continent)) +
  geom_col(color = "white", width = 1) +
  coord_polar(theta = "y", start = 0) +
  
  # 1. Quitar espacios vacíos en X e Y
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  
  # 2. Etiquetas (sin África y Origen desconocido)
  geom_text(
    data = filter(df_cont, !continent %in% c("África", "Origen desconocido")),
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  
  # 3. Paleta y títulos
  scale_fill_brewer(palette = "Blues") +
  labs(
    title = "Distribución de aislados por continente",
    fill  = "Continente"
  ) +
  
  # 4. Tema limpio + aspecto cuadrado + márgenes cero
  theme_minimal(base_size = 14) +
  theme(
    aspect.ratio    = 1,                              # 1:1 para círculo perfecto
    plot.margin     = unit(rep(0, 4), "cm"),          # sin margen
    plot.title      = element_text(face = "bold",
                                   hjust = 0,
                                   size = 15),
    legend.position = "right",
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 12),
    panel.grid      = element_blank(),
    axis.text       = element_blank(),
    axis.title      = element_blank()
  )

