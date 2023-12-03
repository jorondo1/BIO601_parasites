library(pacman); p_load(magrittr)
spe_table <- readRDS("projet_BIO601/brainstorm/plant_data/spe_table.rds")
taxa_table <-  readRDS("projet_BIO601/brainstorm/plant_data/RDP_taxa.rds")

spe_table %>% rownames
