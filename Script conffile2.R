# 1.  SET UP ####
## 1.1 Cargar la librería tidyverse y readxl ####
library(tidyverse)
library(readxl)
library(data.table)

## 1.2 Definir la ruta base ####
#TODO! base_path Argumento de entrada al script. args[1]
base_path <- "//nasibv05/Datos/Grupos/UAMG/rferragud/5.eco_amr_pred/collections/RFF_revised"

# Validar si la ruta base existe
if (!dir.exists(base_path)) {
  stop("La ruta base no existe o no es accesible: ", base_path)
}

# 2. SCRIPT ####

all_summaries <- data.frame("isolate_id_publication" = as.character())

## 2.1 Crear un data.frame vacío con nombres de columnas únicos ####
summary_df <- 
  data.frame(
    isolate_id_publication = character(),
    internal_study_name = character(),
    bioproject_accession = character(),
    sample_accession = character(),
    secondary_sample_accession = character(),
    sample_alias = character(),
    run_accession = character(),
    geographic_location = character(),
    latitude_and_longitude = character(),
    collection_date = character(),
    collected_by = character(),
    host = character(),
    isolation_source = character(),
    host_disease = character(),
    host_tissue_sampled = character()
  )

## 2.2 Obtener la lista de estudios (subdirectorios en "collections") ####

studies <- 
  list.dirs(base_path, recursive = FALSE)

for (study in studies) {
  # Validar si el directorio del estudio existe
  if (!dir.exists(study)) {
    warning("El directorio del estudio no existe: ", study)
    next
  } else {
    
    #study = "//nasibv05/Datos/Grupos/UAMG/rferragud/5.eco_amr_pred/collections/RFF_revised/murphy2021"
    
    # set working directory
    setwd(study)
    
    getwd()
    
    study_name <- basename(study)  # Extrae solo el nombre del directorio
    
    print(head(study_name))
    
    ## 2.3 Buscar las rutas de los archivos de configuración (partA) ####
    
    # Buscar archivos PartA y partB usando Sys.glob para manejar regular expressions
    partA_files <- 
      Sys.glob(file.path(study, paste0("configuration_file.", study_name, "*.partA.v0.*.csv")))
    
    partB_files <- 
      Sys.glob(file.path(study, paste0("configuration_file.", study_name, "*.partB.v0.*.csv")))
    
    if (length(partA_files) == 0 || length(partB_files) == 0) { 
      warning("No se encontraron archivos PartA o PartB para el estudio: ", study_name)
      next } else {
        
    # 3. partA ####
    
    # Lista para almacenar todos los summary_df_partA que se generen en el bucle for 
    all_summary_df_partA <- list() 
    
      #la partB seguiria fent-se no?
      
      for (partA_file in partA_files) { # iterar sobre las partsA del estudio
        
        print("for partA_file in partA_files")
        
        partA <- read.csv(partA_file) %>%
          mutate(across(everything(), str_trim)) # leer partA file
        
        ## 3.2. Obtener las tablas suplementarias e iterar sobre estas ####
        
        # obtener número de tablas suplementarias
        supp_data_tables <- 
          partA %>% 
          filter(source_field != "missing" & source_field != "not applicable") %>% 
          select(source_field) %>% 
          distinct() %>% 
          pull(source_field)
        
        # iterar sobre las tablas suplementarias 
        for (supp_table in supp_data_tables) {
        
          # comprobar cuantas tabs hay 
          supp_table.tabs <-  
            partA %>% 
            filter(source_field == supp_table) %>% 
            select(tab) %>% 
            distinct() %>% 
            pull(tab)
          
          # comprobar tipo de archivo (csv o excel)
          file_extension <- tools::file_ext(supp_table)
          
          ### 3.2.1 Extraer y leer los campos de interés en csv ####
          
          if (file_extension == "csv") {  
            
            print("if (file_extension == csv)")
            
            # Correspondencia standard_field/value de los campos de la supp_table (sin espacios)
            supp_table.fields2extract.temp <-  
              partA %>% 
              filter(source_field == supp_table) %>% 
              select(standard_field, value) %>%
              mutate(value = str_squish(value)) %>% #quitar espacios antes y despues de las palabras
              mutate(value = str_replace_all(value, " ", "\\.")) #quitar los puntos que queden entre las palabras 
            # Si usamos make.names quita las ; 
            
            # Obtener los nombres de las columnas 
            col2extract <- 
              supp_table.fields2extract.temp %>% 
              pull(value) %>% 
              strsplit(";") %>% 
              unlist()   %>%
              map_chr(str_trim) %>%
              str_replace_all( "^\\.", "")
            
            print("col2extract")
            
            # Leer las columnas de interés de las supp_table
            suptab2 <- 
              fread(supp_table, check.names = T) %>%
              select(all_of(col2extract))
          
            ## 3.2.2 Extraer y leer los campos de interés en excel ####
            
          }  else if (file_extension %in% c("xls", "xlsx")) { 
            
            print("else if (file_extension %in% c(xls, xlsx)")
            
            isolate_id_publication.1 <- 
              partA %>% 
              filter(standard_field == "isolate_id_publication") %>% 
              pull(value)
            
            #crear un df vacío con una sola columna isolate_id_publication.1 para después poder unir la info de las tabs
            suptab2 <- setNames(data.frame(as.character()), isolate_id_publication.1)
            
            for (tab in supp_table.tabs) {
              tab.t = supp_table.tabs[as.integer(tab)]
              # Damos por hecho que en todas las tabs hay las mismas isolate_id_publication
              # para cada tab hay un isolate_id_publication column 
              
              print("for (tab in supp_table.tabs)")
              
              # 3.3 Obtener la tabla suplementaria donde se encuentran los sample ids ####
             
              # Extraemos primer los isolate_id_publication de todas las tablas y tabs
              isolate_id_publication.suptable.tab <- 
                partA %>% 
                filter(source_field == supp_table & tab == tab.t &
                         str_starts(standard_field, "isolate_id_publication")) %>% 
                pull(value)
              
              isolate_id_publication.suptable.tab
              
              # TODO comprobar que funciona ####
              if (isolate_id_publication.suptable.tab == "missing") {
                isolate_id_publication.suptable.tab <-  
                  partA %>% 
                  filter(standard_field == "sample_accession") %>% 
                  pull(source_field)
                
                print(head(isolate_id_publication.suptable.tab))
              }  
              
              if (isolate_id_publication.suptable.tab == "missing") {
                isolate_id_publication.suptable.tab <-  
                  partA %>% 
                  filter(standard_field == "secondary_sample_accession") %>% 
                  pull(source_field)
                
                print(head(isolate_id_publication.suptable.tab))
              }  
              
              # Correspondencia standard_field/value de los campos de la supp_table (con espacios)
              supp_table.fields2extract.temp <-  
                partA %>% 
                filter(source_field == supp_table & tab == tab.t) %>% 
                select(standard_field, value) 
              
              # Obtener los nombres de las columnas 
              col2extract <-
                supp_table.fields2extract.temp %>% 
                pull(value) %>% 
                strsplit(";") %>%
                unlist() %>% 
                map_chr(str_trim)
              
              print("col2extract")
              
              # Leer las columnas de interés de las supp_table
              suptab2.temp <- 
                read_excel(supp_table, sheet = as.integer(tab.t)) %>% 
                select(all_of(col2extract)) %>% 
                rename_with(~ isolate_id_publication.1, !!isolate_id_publication.suptable.tab)
              
              # Unir las info de las tabs al df vacío
              suptab2 <- 
                full_join(
                  suptab2, suptab2.temp, 
                  by= isolate_id_publication.1) 
              
            } #  for (tab in supp_table.tabs)
            
            print(head(suptab2))
            
           supp_table.fields2extract.temp <-  
             partA  %>%
             filter(source_field == supp_table) %>% 
             select(standard_field, value) %>% 
             mutate(value =str_replace_all(value, "^\\.", ""))
            
           print("supp_table.fields2extract.temp")
            
          } # else if (file_extension %in% c("xls", "xlsx"))
          
          ## 3.4 Obtener la correspondencia entre standard_field y value #### 
          
          # correspondencia entre standard_field y value con todos los elementos separados (sin ;)
          supp_table.fields2extract <- 
            supp_table.fields2extract.temp %>%
            separate_rows(value, sep = ";") %>%
            mutate(value = str_trim(value)) %>%
            mutate(value = str_replace_all(value, "^\\.", ""))
          
          # Crear identificadores únicos para la columna standard_field
          supp_table.fields2extract$standard_field <- #poner en el if??
            ave(
              supp_table.fields2extract$standard_field,
              supp_table.fields2extract$standard_field,
              FUN = function(x) {
                if (length(x) > 1) {
                  paste0(x, "_", seq_along(x))
                } else {
                  x
                }
              }
            )
          
          print("supp_table.fields2extract")
          
          print(head(supp_table.fields2extract))
          
         ## 3.5 Cambiar los nombres de las variables por los normalizados ####
          
          # Crear el vector de correspondencia 
          
          new_colnames <- setNames(supp_table.fields2extract$standard_field, 
                                   supp_table.fields2extract$value)
          
          # Renombrar las columnas de suptab2 usando la correspondencia
          colnames(suptab2) <- new_colnames[colnames(suptab2)]
          
          print("tabla normalizada")
          print(head(suptab2))
          
          ## 3.6 Unir las columnas que aportan info de la misma variable ####
          
          #encontrar el standard_field para el cual el value tiene ; 
          if (any(str_detect(supp_table.fields2extract.temp$value, ";"))) {
            standard_fields_with_semicolon <- supp_table.fields2extract.temp %>%
              filter(str_detect(value, ";")) %>%
              pull(standard_field)
            
            #hacer una lista que contiene los standard_fields con ; seguidos de _número
            matching_columns.lists <- map(standard_fields_with_semicolon, function(prefix) {
              names(suptab2)[str_detect(names(suptab2), paste0("^", prefix, "_\\d+$"))]
            })
            
            #poner como nombre a las listas el standard_field correspondiente
            names(matching_columns.lists) <- standard_fields_with_semicolon
            
            #bucle en el cual para cada nombre de las listas se unen las columnas 
            for (new_col_name in names(matching_columns.lists)) {
              # Extraer las columnas a unir
              cols_to_unite <- matching_columns.lists[[new_col_name]]
              # Aplicar unite
              suptab2 <- suptab2 %>%
                unite(!!new_col_name, all_of(cols_to_unite), sep = ";", remove = TRUE)
            }
            
            suptab_united <- suptab2 %>%
              mutate(across(everything(), ~ str_remove(., ";$")))
            
            tail(suptab_united)
            
          } else {
            
            suptab_united <- suptab2
          }
          
          print("suptab_united")
          print(head(suptab_united))
          
          ## 3.7 Unir la info obtenida al summary_df ####
          
          # Obtener los nombres de las columnas de ambos dataframes
          suptab_united.columns <- colnames(suptab_united)
          summary_df.columns <- colnames(summary_df)
          
          # Setdiff indica las columnas de summary_df_columns que no están en suptab2_columns
          missing_columns <- setdiff(summary_df.columns, suptab_united.columns)
          
          # Crear un dataframe con las columnas faltantes en suptab2, llenando con NA (numero de filas igual que suptab2)
          missing_data <- data.frame(matrix(NA, ncol = length(missing_columns), nrow = nrow(suptab_united)))
          colnames(missing_data) <- missing_columns
          
          # Unir los datos de suptab2 con los missing_data
          combined_df <- bind_cols(suptab_united, missing_data)
          
          # Poner las columnas ordenadas como en summary_df
          summary_df_suptab <- combined_df %>%
            select(all_of(summary_df.columns))
          
          # Imprimir el dataframe final
          
          print("summary_df_suptab")
          print(head(summary_df_suptab))
          
          ## 3.8 Unir valores del manuscript al summary_df ####
          
          manuscript_fields <-  
            partA %>% 
            filter(source_type == "manuscript") %>% 
            select(standard_field, value)
          
          for (i in 1:nrow(manuscript_fields)) {
            # 1. Extraemos el nombre de la columna del dataframe 'summary_df' que queremos actualizar.
            column_to_update <- manuscript_fields$standard_field[i]
            
            # 2. Extraemos el valor que queremos asignar a dicha columna.
            value_to_assign <- manuscript_fields$value[i]
            
            # 3. Actualizamos la columna correspondiente en 'summary_df'.
            summary_df_suptab[[column_to_update]] <- ifelse(
              is.na(summary_df_suptab[[column_to_update]]),  # Si el valor en la columna es NA...
              value_to_assign,                              # Asignamos el nuevo valor.
              summary_df_suptab[[column_to_update]]          # Si no, mantenemos el valor existente.
            )
          }
          
          print(head(summary_df_suptab))
          
          # Add study name
          summary_df_partA <- 
            summary_df_suptab %>% 
            mutate(internal_study_name = !!study_name)
          
          rm(summary_df_suptab)
          
          # Añadir a la lista los df con partA 
          all_summary_df_partA <- append(all_summary_df_partA, list(summary_df_partA))
          
        } # for table in supp_tables     
        
        #!TODO  Merge info supptables
        
      } # for (partA_file in partA_files)
    
    #!TODO añadir para que pase de este estudio al siguiente
    
    
    #for part in partA 
    # if lenght(partA_file) != 0
    
    
    # Unir todos los summary_df_partA que hay en la lista all_summary_df_partA en un solo data frame
    final_summary_df_partA <- bind_rows(all_summary_df_partA)
    
    # Imprimir la tabla combinada
    print(head(final_summary_df_partA))
    
    view(final_summary_df_partA)
    
    
    
    # 4. ENA file ===========================================================
    ## TODO ####
    # aqui incluir el codigo para sacar el PRJ con wget
    # Buscar el archivo CSV con los accessions del proyecto ENA usando Sys.glob para manejar comodines
    ena_accessions_files <- Sys.glob(file.path(study, "P*_ena_accessions.csv"))
    
    all_ena_accessions_files <- list()
    
    # Validar si los archivos existen y leerlos
    for (ena_accession in ena_accessions_files) {
      #for i in ena_accessions_files
      if (length(ena_accessions_files) == 0) {
        warning("No se encontró ningún archivo de accesiones ENA en: ", study)
      } else {
        # Leer y combinar todos los archivos en un solo data frame
        ena_accessions.df <- read.csv(ena_accession, sep="\t") # no sé si puc gastar tsv
        
        print(ena_accessions.df) # Ver los datos combinados
      }
      
      all_ena_accessions_files <- append(all_ena_accessions_files, list(ena_accessions.df))
    }
    
    final_ena_accessions_files <- bind_rows(all_ena_accessions_files)
    
    
    # 5. Part B ====================================================================
    
    ## 5.1 Iterar sobre todas las partsB del estudio ####
    
    #Lista para almacenar todos los summary_df_partB que se generen en el bucle for 
    
    all_summary_df_partB <- list()
    
    for (partB_file in partB_files) { # iterar sobre las partsB del estudio
      # Validar si los archivos existen
      if (length(partB_file) != 0 ) { # Validar si existen archivos partB 
        partB <- read.csv(partB_file)
        
        print(head(partB))
        
        ## 5.2 Limpiar el archivo de configuración ####
        
        #Eliminar las filas que no nos dan info ("missing") ####
        partB_tidy <- 
          partB %>% 
          filter((standard_field == "to all" & source_type != "missing") | (standard_field != "to all")) %>%
          mutate(across(everything(), str_trim))
        
        print(head(partB_tidy))
        
        ## 5.3 Duplicar las filas "to all" para todos los antibióticos ####
        
        # Filas base "to all"
        to_all_rows <- partB_tidy %>% filter(standard_field == "to all")
        
        # Filas con antibióticos
        antibiotic_rows <- partB_tidy %>% filter(standard_field != "to all")
        
        # Expandir las filas "to all" para cada antibiótico
        new_rows <- antibiotic_rows %>% 
          distinct(standard_field) %>% 
          cross_join(to_all_rows) %>% 
          mutate(standard_field = standard_field.x) %>% 
          select(-standard_field.x, -standard_field.y)
        
        # Combinar los "to all" de los antibióticos con el df original (tiene info propia del atb) 
        partB_tidy_atb <- bind_rows(partB_tidy, new_rows) %>% 
          filter(standard_field != "to all") %>% 
          arrange(standard_field)
        
        print(head(partB_tidy_atb))
        
        #unir standard_field con standard_subfield por _ 
        partB_tidy_atb.suptab <- partB_tidy_atb %>%
          unite(standard_subfield_atb, standard_field, standard_subfield, sep="_")
        
        ## 5.4 Extraer la información de la tabla suplementaria  ####
        
        # Obtener la tabla suplementaria de la partB 
        suptabB <- 
          partB_tidy_atb %>% 
          filter(source_field != "missing" & source_field != "not applicable") %>% 
          select(source_field) %>% 
          distinct() %>% 
          pull(source_field)
        
        # comprobar tipo de archivo (csv o excel)
        file_extension <- tools::file_ext(suptabB)
        
        ### 5.4.1 Extraer y leer los campos de interés en csv ####
        
        if (file_extension == "csv") {
          
          #Obtener los diferentes id de las muestras de la partA. Sirve para unir con la partA   
          isolate_id_publication.partA <-  
            partA %>% 
            filter(standard_field %in% c("isolate_id_publication", "sample_accession", "secondary_sample_accession")) %>% 
            select(standard_field, value) %>% 
            mutate(
              value = str_trim(value)) %>% 
            filter(value != "missing") %>%
            rename(standard_subfield_atb = standard_field)
          
          #TODO alomejor da problemas si hay más de un isolate_id_publication
          
          #Obtener las columnas donde está la información en la tabla suplementaria de la partB
          suptabB.fields2extract <-  
            partB_tidy_atb.suptab %>% 
            filter(source_field == suptabB) %>% 
            select(standard_subfield_atb, value) %>%
            mutate(value = str_replace_all(value, " ", "\\.")) 
          
          #unir la info de la partB con en isolate_id_publication de la partA
          suptabB.fields2extract_id <- 
            bind_rows(suptabB.fields2extract, isolate_id_publication.partA) %>%
            mutate(value = make.names(value)) #quita espacios y caracteres especailes (guiones)
          
          #info + isolate_id_publications value
          col2extract_partB <- suptabB.fields2extract_id %>% 
            pull(value)
          
          #Leer los campos de interés en la tabla suplementaria
          suptabB_info <- fread(suptabB, check.names = T) %>% 
            select(all_of(col2extract_partB))
          
          ### 5.4.2 Extraer y leer los campos de interés en excel #### 
          
        } else if (file_extension %in% c("xls", "xlsx")) {
          
          #Buscar las tabs de la suptabB
          tab.B <- 
            partB_tidy %>%
            filter(source_field == suptabB) %>% 
            select(tab) %>%
            distinct() %>%
            pull(tab)
          
          # Buscar el isolate_id_publication que se encuentra en la tab de los AST
          isolate_id_publication.partB <-
            partA %>% 
            filter(source_field == suptabB & tab == tab.B &
                     str_starts(standard_field, "isolate_id_publication")) %>%
            select(standard_field, value) %>% 
            rename(standard_subfield_atb = standard_field)
          
          isolate_id_publication.partB.value <-
            isolate_id_publication.partB %>%
            pull(value)
          
          # Correspondencia de los campos a extraer
          suptabB.fields2extract <-  
            partB_tidy_atb.suptab %>% 
            filter(source_field == suptabB) %>% 
            select(standard_subfield_atb, value)
          
          # Correspondencia de los campos a extraer con el isolate_id_publication
          suptabB.fields2extract_id.temp <- bind_rows(suptabB.fields2extract, isolate_id_publication.partB)
          
          # Nombres de las columnas a extraer
          col2extract_partB <- suptabB.fields2extract_id.temp %>% 
            pull(value) 
          
          # Obtener isolate_id_publication.partA para después poder unirlo con la partB
          isolate_id_publication.partA <- 
            partA %>% 
            filter(standard_field == "isolate_id_publication") %>%
            select(standard_field, value)%>%
            rename(standard_subfield_atb = standard_field)
          
          suptabB.fields2extract_id <- bind_rows(suptabB.fields2extract, isolate_id_publication.partA)
          
          # Leer los campos de interés de la tabla suplementaria
          suptabB_info <- read_excel(supp_table, sheet = as.integer(tab.B)) %>% 
            select(all_of(col2extract_partB)) %>% 
            rename_with(~ isolate_id_publication.1, !!isolate_id_publication.partB.value)
          # renombramos el isolate_id_publication.partB.value para poder unir todas las partB con el df vacío
          
        } else {
          
          stop("Formato de archivo no soportado. Solo se admiten .csv, .xls y .xlsx")
        }  
        
        ## 5.5 Cambiar los nombres de las variables de la tabla suplementaria por los normalizados ####
        
        # Crear el vector de correspondencia nuevamente
        new_colnames <- setNames(suptabB.fields2extract_id$standard_subfield_atb, 
                                 suptabB.fields2extract_id$value)
        
        # Renombrar las columnas de suptab2 usando la correspondencia
        colnames(suptabB_info) <- new_colnames[colnames(suptabB_info)]
        
        print(head(suptabB_info))
        
        ## 5.6 Unir valores del manuscript al summary_df_partA_suptabB ####
        
        # Obtener los values y standard_subfield_atb de los campos que vienen del manuscript
        manuscript_fields.partB <-  
          partB_tidy_atb.suptab %>% 
          filter(source_type == "manuscript") %>% 
          select(standard_subfield_atb, value)
        
        head(manuscript_fields.partB)
        
        # Pivotar para obtener el df con el standard_subfield_atb y los valores 
        # correspondinetes y poner los valores para todas las filas del summary_df_partA_suptabB 
        expanded_manuscript_fields.partB <- 
          manuscript_fields.partB %>%
          pivot_wider(names_from = standard_subfield_atb, values_from = value) %>%
          slice(rep(1, nrow(suptabB_info)))
        
        # Unir info extraída de la tabla suplementaria (suptabB_info) y info del manuscript 
        summary_df_partB <- bind_cols(suptabB_info, expanded_manuscript_fields.partB)
        
        # Ver el resultado
        print(head(summary_df_partB))
        
        # Unir todas las partsB
        all_summary_df_partB <- append(all_summary_df_partB, list(summary_df_partB))
        
        final_summary_df_partB <- bind_rows(all_summary_df_partB)
        
        view(final_summary_df_partB)
        
      } else {
        warning("No se encontró archivo PartB para el estudio: ", study_name)
      }
    }
    
    # 6. Unir info de partA y partB ####
    
    # Realizar la unión con merge
    
      summary_df_partA_partB <-
        merge(final_summary_df_partA, final_summary_df_partB,
               by = intersect(names(final_summary_df_partA), names(final_summary_df_partB)), 
               all = TRUE)
  
    # 7. Unir todos los estudios #### 
    all_summaries <-
      merge(all_summaries, summary_df_partA_partB,
             by = intersect(names(all_summaries), names(summary_df_partA_partB)),
             all = TRUE)
    
    view(all_summaries)
    
      }
    
  } #if (!dir.exists(study))
  
}#for study in studies



  

