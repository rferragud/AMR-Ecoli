
#--------------- LISTADO DE MUESTRAS DE OTRAS BASES DE DATOS  ------------------

setwd("C:/Users/rferragud/Documents/Projecte/RFF_revised")

# 1. Recopilar info del NCBI #### 

## 1.1. leer csv con todas las muestras de E. coli del NCBI con AST y genoma ####
muestras_NCBI <- read_tsv("C:/Users/rferragud/Documents/Projecte/Tables/asts.tsv")

# muestras (+ metadatos) recopiladas en los artículos 
df_cruzado_completo <- read.csv("C:/Users/rferragud/Documents/Projecte/Tables/df_cruzado_completo_articles.csv", check.names = FALSE)

## 1.2. df con info de muestras del NCBI que no coinciden con las muestras de los artículos ####
muestras_NCBI.nuevas <- anti_join(muestras_NCBI, df_cruzado_completo, by = c("#BioSample" = "ENA_sample_accession"))

#biosample de las muestras que no coinciden (10717)
biosample_NCBI.nuevas <- muestras_NCBI.nuevas %>% pull(`#BioSample`)%>% unique()

#df de las muestras coincidentes 
muestras_NCBI.dup <- inner_join(muestras_NCBI, df_cruzado_completo, by = c("#BioSample" = "ENA_sample_accession"))

#biosample de las muestras coincidentes (53)
muestras_NCBI.dup %>% pull(`#BioSample`) %>% unique()

# biosamples del NCBI
muestras_NCBI %>% pull(`#BioSample`) %>% unique()

# biosamples que vienen de los artículos
df_cruzado_completo %>% pull(ENA_sample_accession)%>% unique()

## 1.3. Comprobar que no he hecho un error y que no tenemos esos bioprojects 
#cuántas muestras hay de cada bioproject
muestras_NCBI.nuevas %>% select(`#BioSample`, BioProject) %>%  distinct() %>%  count(BioProject) %>% arrange(desc(n))

#1 PRJNA292663  2626
#2 PRJNA809394  2076
#3 PRJNA966974  1128
#4 PRJNA292664  1083
#5 PRJNA278886   51
#6 PRJNA324573   379
#7 PRJNA885502   337
#8 PRJNA318591   328

#biosamples del bioproject X (no deberían estar en nuestra colección)
muestras_NCBI.nuevas %>% filter(BioProject == "PRJNA292664") %>% pull(`#BioSample`) %>% unique()

#buscar en el df de los artículos qué muestras estan de un bioproject concreto (no debería haber)
df_cruzado_completo %>% filter(str_detect(bioproject_accession, "PRJNA809394")) 


## 1.4. Seleccionar solo las muestras de E. coli ####

#cuántas muestras hay de cada organismo
muestras_NCBI.nuevas %>% select(`#BioSample`, `Organism group`) %>%  distinct() %>%  count(`Organism group`) %>% arrange(desc(n))

#`Organism group`         n
#<chr>                <int>
#  1 E.coli and Shigella   9313
#2 Campylobacter jejuni  1404

#Nos quedamos solo con E.coli 
muestras_NCBI.nuevas.coli <-  muestras_NCBI.nuevas %>% filter(`Organism group` == "E.coli and Shigella")

## 1.5. df con solo las biosample 
NCBI_ENA_accessions.nuevas <- 
  muestras_NCBI.nuevas.coli %>% 
  select(`#BioSample`) %>% 
  rename(ENA_sample_accession = `#BioSample`) %>%
  distinct() 


#----------------------------------------------------------------------------####

# 2. Recopilar info del PATRIC #### 

## 2.1. Leer la tabla con la info de las muestras de E. coli del PATRIC con AST y genoma ####
muestras_patric <- read.csv("C:/Users/rferragud/Documents/Projecte/RFF_revised/BVBRC_genome_amr (1).csv", colClasses = "character")

# muestras del PATRIC
muestras_patric %>% pull(Genome.ID)%>% unique()

#Como esta tabla no se puede cruzar con la nuestra pq no tiene la biosample
#leemos otra tabla que tiene la correspondencia 
muestras_patric_biosample <- read.csv("C:/Users/rferragud/Documents/Projecte/RFF_revised/BVBRC_genome.csv", colClasses = "character")

## 2.2. df con las muestras del PATRIC y sus biosample ####
muestras_patric_cruzadas <- left_join(muestras_patric, muestras_patric_biosample, by = "Genome.ID")

#seleccionar solo las muestras de E. coli
muestras_patric_cruzadas %>% count(Species)
muestras_patric_cruzadas.Ecoli <- muestras_patric_cruzadas %>% filter(Species == "Escherichia coli")

## 2.3. comprobar qué muestras tenemos ya nosotros ####

#como algunas de las muestras no estan bien clasificadas, clasificamos por regex las err y sam 
muestras_patric_cruzadas_accessions <- muestras_patric_cruzadas.Ecoli %>%
  mutate(
    match_sam.p = apply(., 1, function(row) {
      str_extract(paste(row, collapse = " "), "SAM(E|D|N)[A-Z]?[0-9]+")
    }),
    match_ers.p = apply(., 1, function(row) {
      str_extract(paste(row, collapse = " "), "(E|D|S)RS[0-9]{6,}")
    }),
    match_erx.p = apply(., 1, function(row) {
      str_extract(paste(row, collapse = " "), "(E|D|S)RX[0-9]{6,}")
    }),
    match_err.p = apply(., 1, function(row) {
      str_extract(paste(row, collapse = " "), "(E|D|S)RR[0-9]{6,}")
    })
  )

dim_ERS.p <-  muestras_patric_cruzadas_accessions %>%
  filter(!is.na(match_ers.p))%>%
  dim()

dim_ERX.p <-  muestras_patric_cruzadas_accessions %>%
  filter(!is.na(match_erx.p))%>%
  dim()

dim_ERR.p <-  muestras_patric_cruzadas_accessions %>%
  filter(!is.na(match_err.p)) %>%
  distinct()%>%
  dim()

dim_SAM.p <- muestras_patric_cruzadas_accessions %>%
  filter(!is.na(match_sam.p)) %>%
  distinct()%>%
  dim()


#> dim_ERS.p
#[1] 38721   111
#> dim_ERX.p
#[1]   0 111
#> dim_ERR.p
#[1] 167364    111
#> dim_SAM.p
#[1] 128829    111

# empezamos cruzando las match_sam.p con ENA_sample_accession pq lo podemos cruzar con el NCBI 

data_SAM.p <- muestras_patric_cruzadas_accessions %>%
  filter(!is.na(match_sam.p)) %>% 
  distinct()%>%
  anti_join(df_cruzado_completo, by = c("match_sam.p" = "ENA_sample_accession"))
  
dim(data_SAM.p)

#[1] 124896    111

# obtener los match_sam.p cruzados 
sam.p <- data_SAM.p %>% pull(match_sam.p)

# cruzar match_err con ENA_run_accession (sin tener en cuenta las match_sam.p)
data_ERR.p <-     
  muestras_patric_cruzadas_accessions %>%
  filter(!is.na(match_err.p)) %>% 
  distinct()%>%
  filter(!match_sam.p %in% sam.p)%>%
  anti_join(df_cruzado_completo, by = c("match_err.p" = "ENA_run_accession"))

dim(data_ERR.p)

data_ERR.p %>% pull(SRA.Accession) %>% unique()

err.p <- data_ERR.p %>% pull(match_err.p)
# [1] 527 111

# obtener los match_err cruzados 
err.p <- data_ERR.p %>% pull(match_err.p)

data_ERR.p %>% filter(!match_err.p %in% err.p)

# cruzar match_ers con ENA_run_accession (sin tener en cuenta las muestras ya cruzadas)
data_ERS.p <-     
  muestras_patric_cruzadas_accessions %>%
  filter(!is.na(match_ers.p)) %>% 
  distinct()%>%
  filter(!match_err.p %in% err.p)%>%
  filter(!match_sam.p %in% sam.p)%>%
  anti_join(df_cruzado_completo, by = c("match_ers.p" = "ENA_secondary_sample_accession")) 

dim(data_ERS.p)
# [1]  27 111

# obtener los match_ers cruzados 
ers.p <- data_ERS.p %>% pull(match_ers.p)

# df con las muestras nuevas que no teníamos 
muestras_patric.nuevas.temp  <- 
  bind_rows(data_ERR.p, data_ERS.p, data_SAM.p)

muestras_patric.nuevas <- muestras_patric.nuevas.temp %>%
  mutate(
    match_final_sam = ifelse(!is.na(match_sam.p), match_sam.p, NA),
    match_final_err = ifelse(is.na(match_sam.p) & !is.na(match_err.p), match_err.p, NA),
    match_final_ers = ifelse(is.na(match_sam.p) & is.na(match_err.p) & !is.na(match_ers.p), match_ers.p, NA)
  )

muestras_patric.nuevas %>% distinct(match_final_err)%>% dim()
# [1] 33  1
muestras_patric.nuevas %>% distinct(match_final_ers)%>% dim()
# [1] 3 1
muestras_patric.nuevas %>% distinct(match_final_sam)%>%dim()
# [1] 6290    1

# seleccionar las columnas match_err.p, match_ers.p, match_sam.p
patric_ENA_accessions.temp <- muestras_patric.nuevas %>% select(match_final_err, match_final_ers, match_final_sam)

patric_ENA_accessions.temp%>%pull(match_final_sam)%>%unique()
#6290

#quitamos las muestras duplicadas 
patric_ENA_accessions <- patric_ENA_accessions.temp %>%
  filter(
    !duplicated(match_final_err) | is.na(match_final_err),
    !duplicated(match_final_ers) | is.na(match_final_ers),
    !duplicated(match_final_sam) | is.na(match_final_sam)
  )

dim(patric_ENA_accessions)
# [1] 6323    3


# seleccionamos las biosample que no estén ya en las muestras del NCBI
# las run_accession si no hay biosample no las podemos comprobar 
patric_ENA_accessions.nuevas.temp <- 
  patric_ENA_accessions %>%
  filter(!match_final_sam %in% biosample_NCBI.nuevas) %>% 
  select(match_final_sam, match_final_err, match_final_ers) %>%
  rename(
    ENA_sample_accession = match_final_sam, 
    ENA_run_accession = match_final_err, 
    ENA_secondary_sample = match_final_ers
  )

patric_ENA_accessions.nuevas.temp %>% dim()
#[1] 5775    3


patric_ENA_accessions.nuevas.temp %>% filter(is.na(ENA_sample_accession) & is.na(ENA_run_accession) & is.na(ENA_secondary_sample))

patric_ENA_accessions.nuevas <- patric_ENA_accessions.nuevas.temp %>%
  mutate(
    ENA_sample_accession = case_when(
      ENA_secondary_sample == "ERS1995029" ~ "SAMEA104370051",
      ENA_secondary_sample == "ERS1995244" ~ "SAMEA104370266",
      TRUE                                        ~ ENA_sample_accession
    )
  )%>% select(-ENA_secondary_sample)



patric_ENA_accessions.nuevas %>% filter(!is.na(ENA_run_accession)) %>% pull(ENA_run_accession)

# 3. Escribir los archivos csv con los nombres de las muestras  

write.csv(patric_ENA_accessions.nuevas, "C:/Users/rferragud/Documents/Projecte/Tables/patric_ENA_accessions.csv", row.names = FALSE)

write.csv(NCBI_ENA_accessions.nuevas, "C:/Users/rferragud/Documents/Projecte/Tables/NCBI_ENA_accessions.csv", row.names = FALSE)

write.csv(patric_ENA_accessions.nuevas.temp, "C:/Users/rferragud/Documents/Projecte/Tables/patric_ENA_accessions.nuevas.temp.csv", row.names = FALSE)

# 4. df con todas las biosamples y sus metadatos + datos fenotipos

write.csv(muestras_NCBI.nuevas.coli, "C:/Users/rferragud/Documents/Projecte/Tables/muestras_NCBI.nuevas.coli.csv", row.names = FALSE)

write.csv(muestras_patric.nuevas, "C:/Users/rferragud/Documents/Projecte/Tables/muestras_patric.nuevas.csv", row.names = FALSE)

muestras_patric.nuevas %>% distinct(Genome.ID)
getwd()
