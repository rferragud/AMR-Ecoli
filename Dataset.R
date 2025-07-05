
library(readxl)
dataset <- read_xlsx("C:/Users/rferragud/Documents/Projecte/dataset/RFF-Escherichia_coli.bioproject_accessions.170624.v1.xlsx")

warnings()

dataset_eco_run <- dataset %>% filter(num_runs_eco > 0)

dataset_eco_run_200 <- dataset %>% filter(num_runs_eco >= 200)


dataset %>% colnames()
