source("fun.R")
library(pacman)
p_load(tidyverse, magrittr, phyloseq, vegan, DESeq2, knitr, kableExtra,
       microViz, patchwork)
plant.ps <- read_rds("data/ps/plant.ps.RDS")
bumble.ps <- read_rds("data/ps/bumble.ps.RDS")
beetle.ps <- read_rds("data/ps/beetle.ps.RDS")
###############################################################################
### Visualisation des communautés #############################################
###############################################################################

# Groups
plant.ps@sam_data %>% 
  group_by(compart, ID) %>% 
  summarise(n = n())

bumble.ps@sam_data %>% 
  group_by(treatment, type) %>% 
  summarise(n = n())

beetle.ps@sam_data %>% 
  group_by(status, time) %>% 
  summarise(n = n())

# Community overview
comm_barplot(bumble.ps, 'treatment', 'Species', 8, 'Bombus terrestris', 'type',
             'Abondance relative des espèces microbiennes chez B. terrestris.')

comm_barplot(beetle.ps, 'status', 'Species', 14, 'Tenebrio molitor', 'time',
             'Abondance relative des espèces microbiennes chez T. molitor.') # ADD A FACET for time

comm_barplot(plant.ps, 'ID', 'Order',14, 'Hedera','compart', 
             'Abondance relative des ordres microbiens chez Hedera spp. par compartiment.')

###############################################################################
### PCoA sur les dissimilarités Bray Curtis ###################################
###############################################################################

######################
### Méthode du cours #
#######################
dbRDA_results <- list()

### **************** ###
### TENEBRIO MOLITOR ###
### **************** ###

# PLOTTING BY TIME 
ord.fun(beetle.ps, groupVar = "time", group = "2")
ord.fun(beetle.ps, groupVar = "time", group = "7")

beetle.p1 <- plot_PCoA(beetle_time2_bray, group = "status", title = 'Jour 2') 
beetle.p2 <- plot_PCoA(beetle_time7_bray, group = "status", title = 'Jour 7') 
beetle.col <- c("Parasités" = "orangered", "Contrôles" = "darkolivegreen3")
ord_Tmolitor <- customPatchwork(beetle.p1 + beetle.p2) &
  plot_annotation(title = "Microbiote des ténébrions 2 et 7 jours après infection par Hymenoleps diminuta.") &
  scale_colour_manual(values = beetle.col) &
  scale_fill_manual(values = beetle.col) &
  labs(colour = "Traitement") 

ggsave("out/ord_Tmolitor.png", plot = ord_Tmolitor, 
       width = 48, height = 27, units = 'cm')

# REDUNDANCY ANALYSIS
# Temps 2
(dbRDA_results[["Tmolitor_2"]] <- 
  adonis2(beetle_time2_bray$dist ~ status, 
        data = beetle_time2_bray$data,
        permutations=9999, method="bray") # p = 0.16
)
# Temps 7
(dbRDA_results[["Tmolitor_7"]] <- 
  adonis2(beetle_time7_bray$dist ~ status, 
        data = beetle_time2_bray$data,
        permutations=9999, method="bray") # p = 0.00001
)

## RCLR ??


### ***************** ###
### BOMBUS TERRESTRIS ###
### ***************** ###

# PLOTTING BY TIME
ord.fun(bumble.ps, groupVar = "type", group = "faeces")
ord.fun(bumble.ps, groupVar = "type", group = "gut")
ord.fun(bumble.ps)

bumble.p1 <- plot_PCoA(bumble_typefaeces_bray, group = "treatment", title = 'Jour 0') 
bumble.p2 <- plot_PCoA(bumble_typegut_bray, group = "treatment", title = 'Jour 8') 

bumble.col <- c("Parasités" = "darkorchid4", "Contrôles" = "yellow3")
ord_Bterrestris <- customPatchwork(bumble.p1 + bumble.p2) & 
  plot_annotation(title = "Microbiote du bourbon avant et après infection par Chritidia bombi.") &
  scale_colour_manual(values = bumble.col) &
  scale_fill_manual(values = bumble.col) &
  labs(colour = "Traitement") 

ggsave("out/ord_Bterrestris.png", plot = ord_Bterrestris, 
       width = 48, height = 27, units = 'cm')

# REDUNDANCY ANALYSIS
# Jour 0
(dbRDA_results[["Bterrestris_0"]] <- 
  adonis2(bumble_typefaeces_bray$dist ~ treatment, 
        data = bumble_typefaeces_bray$data,
        permutations=9999, method="bray") # p = 0.76
)
# Jour 8
(dbRDA_results[["Bterrestris_8"]] <-
  adonis2(bumble_typegut_bray$dist ~ treatment, 
        data = bumble_typegut_bray$data,
        permutations=9999, method="bray") # p = 0.98
)# Rien !

# PLOTTING BY COMPARTMENT
bumble.col2 <- c("faeces" = "darkkhaki", "gut" = "orangered")
bumble.p3 <- plot_PCoA(bumble__bray, group = "type") 
ord_Bterrestris_compart <- bumble.p3 +
  plot_annotation(title = "Comparaison des compartiments chez le bourdon.") &
  theme_light(base_family = "Baskerville", base_size = 26) &
  theme(plot.title = element_text(size = 32)) &
  theme(legend.position = 'bottom') &
  scale_colour_manual(values = bumble.col2,
                      labels = c("faeces" = "Fèces", "gut" = "Intestin")) &
  scale_fill_manual(values = bumble.col2) &
  labs(colour = "Compartiment") 

ggsave("out/ord_Bterrestris_compart.png", plot = ord_Bterrestris_compart, 
       width = 27, height = 27, units = 'cm')

# REDUNDANCY ANALYSIS
# Restreindre les permutations par colonie
(dbRDA_results[["Bterrestris_compart"]] <- 
  adonis2(bumble__bray$dist ~ type + treatment, 
        data = bumble__bray$data,
        permutations=9999, method="bray") # p = 0.0001 
)
# Se concentrer sur la différence entre gut et feces; si différence, 
# souligner l'impact potentiel sur les analyses réalisées par les auteurs
# Interesting!

### ********** ###
### HEDERA SPP ###
### ********** ###

# PLOTTING BY COMPARTMENT
ord.fun(plant.ps, groupVar = "compart", group = "L")
ord.fun(plant.ps, groupVar = "compart", group = "R")
ord.fun(plant.ps, groupVar = "compart", group = "S")

plant.p1 <- plot_PCoA(plant_compartL_bray, group = "ID", title = 'Feuilles') 
plant.p2 <- plot_PCoA(plant_compartR_bray, group = "ID", title = 'Racines') 
plant.p3 <- plot_PCoA(plant_compartS_bray, group = "ID", title = 'Sol') 

plant.col <- c("Parasités" = "orangered", "Non-parasités" = "darkolivegreen3")
ord_Hedera <- customPatchwork(plant.p1 + plant.p2 + plant.p3) &
  plot_annotation(title = "Microbiote de Hedera en présence ou absence de Orobranche hederae.") &
  scale_colour_manual(values = plant.col) &
  scale_fill_manual(values = plant.col) &
  labs(colour = "Statut") 

ggsave("out/ord_Hedera.png", plot = ord_Hedera, 
       width = 50, height = 23, units = 'cm')

# REDUNDANCY ANALYSIS
# Leaves
(dbRDA_results[["Hedera_L"]] <- 
  adonis2(plant_compartL_bray$dist ~ ID, 
          data = plant_compartL_bray$data,
          permutations=99999, method="bray") # p = 0.77
)

# Racines
(dbRDA_results[["Hedera_R"]] <- 
  adonis2(plant_compartR_bray$dist ~ ID, 
          data = plant_compartR_bray$data,
          permutations=99999, method="bray") # p = 0.001
)

# Sol
(dbRDA_results[["Hedera_S"]] <- 
  adonis2(plant_compartS_bray$dist ~ ID, 
          data = plant_compartS_bray$data,
          permutations=99999, method="bray") # p = 0.03
)

# Cependant, il y a une énorme variabilité d'un réplicata à l'autre:
depth.plot <- ggplot(plant.ps@sam_data, aes(x = compart, y = depth, fill = ID)) +
  geom_boxplot() + 
  theme_light(base_family = "Baskerville", base_size = 26) +
  labs(fill = "Compartiment", x = "",
       y = "Nombre de séquences totales par échantillon",
       title = "Profondeur de séquençage par compartiment.") +
  scale_fill_manual(values = plant.col) +
  scale_x_discrete(labels = c("L" = "Feuilles", "R" = "Racines", "S" = "Sol"))

ggsave("out/depth_boxplot.png", plot = depth.plot, 
       width = 50, height = 30, units = 'cm')

# Est-ce ce que l'ordination capture en fait la différence de profondeur de séquençage?
# Vérifions avec les échantillons Racines et Sol, qui semblent avoir une profondeur
# différente entre les parasités et contrôles?
plot_PCoA(plant_compartS_bray, group = "depth", ellipse = FALSE) 
plot_PCoA(plant_compartR_bray, group = "depth", ellipse = FALSE) 

(dbRDA_results[["Hedera_S_depth"]] <- 
  adonis2(plant_compartS_bray$dist ~ depth, 
        data = plant_compartS_bray$data,
        permutations=99999, method="bray") # p = 0.016
  )

#### dbRDA conditionnelle 

# RACINES
ord_R <- capscale(plant_compartR_bray$dist ~ ID + Condition(depth), 
                  data = plant_compartR_bray$data)

ord_R_anova <- anova(ord_R, permutations = how(nperm=99999))
ord_R_anova$R2 <- ord_R$CCA$tot.chi/ord_R$tot.chi
dbRDA_results[["Hedera_R_cond"]] <- ord_R_anova

# SOL
ord_S <- capscale(plant_compartS_bray$dist ~ ID + Condition(depth), 
                data = plant_compartS_bray$data)

ord_S_anova <- anova(ord_S, permutations = how(nperm=99999))
ord_S_anova$R2 <- ord_S$CCA$tot.chi/ord_S$tot.chi
dbRDA_results[["Hedera_S_cond"]] <- ord_S_anova

# Fooling around
ord.fun(plant.ps)

capscale(plant__bray$dist ~ ID + Condition(depth) + type + IDtypesite, 
         data = plant__bray$data) %>% 
  anova(permutations = how(nperm=9999))

adonis2(plant__bray$dist ~ ID + depth + type + IDtypesite,
        data = plant__bray$data, permutations = 9999)
############################
### COMPILE RESULTS ########
############################

results <- lapply(names(dbRDA_results), function(test_name) {
  test_result <- dbRDA_results[[test_name]]
  
  # Extract R2 and Pr(>F) for the 'treatment' row
  data.frame(
    Test = test_name,
    `R2` = test_result$R2[1],
    p = test_result$`Pr(>F)`[1]
  )
})


do.call(rbind, results) %>%
  mutate(` `= case_when(p<0.001 ~ '***',
                       p<0.01 ~ '**',
                       p<0.05 ~ '*',
                       TRUE ~''),
         R2 = round(R2, 2), 
         p = round(p, 4)) %>%  # Round values to two significant digits
  kable(caption = "Résultats des tests de permutations par analyses de redondance basée sur les distances.") %>% 
  kable_styling(bootstrap_options = "striped", full_width = F)


# Est-ce qu'une méthode mieux adaptée permettrait de comparer les échantillons
# d'un compartiment à l'autre?
plantRCLR.ps <- tax_transform(plant.ps, "rclr")

ord.fun(plantRCLR.ps, groupVar = "compart", group = "R",  d = "euclidean")
ord.fun(plantRCLR.ps, groupVar = "compart", group = "S",  d = "euclidean")

plot_PCoA(plantRCLR_compartR_euclidean, group = "ID") 
plot_PCoA(plantRCLR_compartS_euclidean, group = "ID") 

adonis2(plantRCLR_compartR_euclidean$dist ~ ID + depth, 
        data = plantRCLR_compartR_euclidean$data,
        permutations=99999, method="euclidean") # significant AF

# Can we use all datasets together to see if an effect can be detected?
ord.fun(plantRCLR.ps, d = "euclidean")

plot_PCoA(plantRCLR__euclidean, group = "ID") 

rest_perm2 <- with(plantRCLR__euclidean$data, how(nperm= 9999, blocks = compart)) 
adonis2(plantRCLR__euclidean$dist ~ ID*compart + depth, 
        data = plantRCLR__euclidean$data,
        permutations=rest_perm2, method="euclidean") # p = 0.016


####################################################################################




##########################################
### Méthode du lab ########################
### Variance stabilizing transformation ###
##########################################

vst.fun(beetle.ps, groupVar = "time", group = "2", group = "status")
vst.fun(beetle.ps, groupVar = "time", group = "7", group = "status")
vst.fun(beetle.ps, group = "cohort")
# flagrant

vst.fun(bumble.ps, groupVar = "type", group = "faeces", group = "treatment")
vst.fun(bumble.ps, groupVar = "type", group = "gut", group = "treatment")
vst.fun(bumble.ps, groupVar = "treatment", group = "sham",  group = "type") # faeces vs gut seulement
      # À reessayer avec CLR! 
# Énorme différence entre gut et faeces
vst.fun(bumble.ps, group = "treatment") # traitement seulement

vst.fun(plant.ps, groupVar = "type", group = "L", group = "ID")
vst.fun(plant.ps, groupVar = "type", group = "R", group = "ID")
vst.fun(plant.ps, groupVar = "type", group = "S", group = "ID")

# Très peu de samples par groupe, essayons tous ensemble:
vst.fun(plant.ps, group = "ID")
# Patron similaire à l'approche capscale


### TOUTES
# Sensibilité du test stat à un contrôle sur la profondeur de séquençage

#### Hedera 
# Ordination Roots (avec IR), une Leaf
# Une comparaison avec Unifrac?? Attendre
# Modèle à effets mixtes avec réplicats?

#### Tenebrio
# Une ordination pour chaque temps 

#### Bommus
# Une ordination par compartiment +
# Une comparaison Sham avec CLR (à moins que ce soit celui-ci qui )