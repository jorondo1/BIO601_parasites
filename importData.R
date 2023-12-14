# Main: assemble all 3 datasets in 3 coherent ps objects
library(pacman); p_load(ape, tidyverse, magrittr, phyloseq, tools, genefilter)
source("fun.R")
taxLvls <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

##########################
### TENEBRIO MOLITOR ######
##########################

beetle.raw <- readxl::read_xlsx("data/beetle_data/pone.0227561.s004.xlsx") %>% 
  rename_with(tools::toTitleCase, .cols = domain:species) %>% 
  mutate(OTU = sprintf("OTU.%03d", row_number())) %>% 
  column_to_rownames('OTU')

beetle.tax <- beetle.raw %>% 
  transmute_at(taxLvls, rename.fun) %>% 
  mutate(Species = paste0(Genus, '_', Species)) %>% 
           as.matrix

beetle.otu <- beetle.raw %>% select(-taxLvls)

beetle.sam <- beetle.otu %>% colnames %>% 
  tibble(sampleID = .) %>% 
  mutate(
    status = case_when(sampleID %>% startsWith('I') ~ 'Parasités',
                       sampleID %>% startsWith('U') ~ 'Contrôles'),
    time = str_extract(sampleID, "(?<=^[IU])[0-9]"),
    beetleID = str_extract(sampleID, "(?<=\\.)\\d+"),
    cohort = paste0(status,time)
  ) %>% 
  column_to_rownames('sampleID')

beetle.ps <- phyloseq(otu_table(beetle.otu, taxa_are_rows = T),
                      sample_data(beetle.sam),
                      tax_table(beetle.tax)) %>% 
  depth.fun()

###################
### BUMBLEBEE ######
###################

bumble.raw <- read.delim("data/bumblebee_data/D2_OTUcounts.txt") %>%  t 
colnames(bumble.raw) <- bumble.raw[1,]
bumble.otu <- bumble.raw %>% .[-1,] %>% as.data.frame %>% 
  mutate(across(everything(), as.numeric))

bumble.tax <- readxl::read_xlsx("data/bumblebee_data/D4_OTUTaxonomy.xlsx") %>% 
  column_to_rownames("OTU_name") %>% 
  rename_with(~ str_replace_all(., "RDP_", "")) %>% 
  mutate(Domain = Rank1) %>% 
  select(all_of(taxLvls)) %>% 
  as.matrix 

bumble.sam <- readxl::read_xlsx("data/bumblebee_data/D1_Metadata_Samples.xlsx") %>% 
  column_to_rownames("unique.sample.id") %>% 
  filter(treatment.line != 'neg.contr' & 
           microbiota.sequenced == 'yes') %>% 
  transmute(type = case_when(sample.type == 'faeces' ~ 'Fèces',
                             sample.type == 'gut' ~ 'Intestin'),
            colony = as.factor(colony.id),
            treatment = case_when(treatment.line == 'sham' ~ 'Contrôles',
                                  treatment.line == 'cocktail' ~ 'Parasités'),
            bee.id = as.factor(bee.id),
            origin = as.factor(population.origin))

bumble.ps <- phyloseq(otu_table(bumble.otu, taxa_are_rows = T),
                      sample_data(bumble.sam),
                      tax_table(bumble.tax)) %>% 
  # Remove contaminated sample
  prune_samples(sample_names(.) != "996", .)  %>% 
  # Remove taxa without any OTU across all samples
  prune_taxa(taxa_sums(.) > 0,.) %>% 
  # Rarefy
  rarefy_even_depth(
    sample.size = 13820,
    verbose = TRUE,  
    trimOTUs = TRUE,
    rngseed = 4466)

################
### PLANTS ######
################

plant.otu.raw <- readRDS("data/plant_data/spe_table.rds") %>% t %>% as.data.frame
# There are 21865 in there... let's remove the very low abundance ones.
keepSpecies <- apply(plant.otu.raw, 1, function(row) max(row) > 100) 
plant.otu <- plant.otu.raw %>% 
  .[keepSpecies, ] %>% 
  mutate(row_id = row.names(.)) %>% 
  select(-c("mock", 'PAO1', 'Undetermined', 'water')) %>% 
  mutate(across(everything(), as.numeric))

plant.tree <- read.tree(file= "data/plant_data/PPM_bac.tre")

plant.tax <- readRDS("data/plant_data/RDP_taxa.rds")

plant.sam <- plant.otu %>% t %>% rownames %>% as_tibble %>% 
  separate(value, into = c("site","ID","type","replicate"), "-", remove = F) %>% 
  mutate(
    compart = case_when(type == "IR" ~ "Racines",
                        type == 'R' ~ 'Racines',
                        type == 'L' ~ 'Feuilles',
                        type == 'S' ~ 'Sol'),
    IDtype = interaction(ID,type),
    IDtypesite = interaction(IDtype, site),
    ID = case_when(ID == 'I' ~ 'Parasités',
                   ID == 'U' ~ 'Non-parasités')
    ) %>% 
  left_join(read_csv("data/plant_data/sample_yields.csv"),
            join_by(value == sample)) %>% 
  filter(ID!='P') %>%  # enlever les séquneces du parasite
  column_to_rownames('value')

plant.ps <- phyloseq(otu_table(plant.otu, taxa_are_rows = T), 
                sample_data(plant.sam), 
                tax_table(plant.tax),
                phy_tree(plant.tree)) %>% 
  prune_taxa(taxa_sums(.) > 0,.) %>% 
  subset_taxa(!is.na(Family) & !Class == 'Chloroplast') %>% 
  subset_samples(site != "Undetermined" & site != "mock" & 
                   site != "PAO1" & site != "water" & compart != 'Feuilles') %>% 
  filter_taxa(kOverA(5,A=25),TRUE) %>% 
# Il y a une grosse variabilité de profondeur de séquencage dans ces données,
# on l'ajoute aux métadonnées pour pouvoir la modéliser plus tard:
  depth.fun()


### EXPORT
write_rds(plant.ps, "data/ps/plant.ps.RDS")
write_rds(bumble.ps, "data/ps/bumble.ps.RDS")
write_rds(beetle.ps, "data/ps/beetle.ps.RDS")


