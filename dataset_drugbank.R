set.seed(9427)
setwd("D:/Bioestadística/gam-farmacovigilancia")
source("00_functions.R", local = TRUE)

install.packages("dbparser")


library(dbparser)
library(dplyr)

# 1. Cargar la base de datos de DrugBank (cambia el nombre por el de tu archivo XML)
drug_bank <- parseDrugBank("full_database.xml")

# 2. Extraer los códigos ATC de todas las drogas
atc_data <- subset_drugbank_dvobject(drug_bank, drugs_ids) 
# Esto te da un dataframe con "drugbank_id" y "atc_code"

# 3. Extraer los targets (los genes/proteínas sobre los que actúa la droga)
targets_data <- targets_polypeptides() 
# Esto te da el "drugbank_id" y el "gene_name" (ej. CYP3A4, EGFR, etc.)

# 4. Unir ambas tablas para tener tu diccionario final ATC -> Gen
diccionario_drogas_genes <- atc_data %>%
  inner_join(targets_data, by = "drugbank_id") %>%
  select(atc_code, gene_name) %>%
  distinct() # Eliminamos duplicados

