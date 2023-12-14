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
             'Abondance relative des espèces microbiennes chez T. molitor.')

comm_barplot(plant.ps, 'ID', 'Order',14, 'Hedera','compart', 
             'Abondance relative des ordres microbiens chez Hedera spp. par compartiment.')

###############################################################################
### PCoA sur les dissimilarités Bray Curtis ###################################
###############################################################################
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

# Cette différence est-elle explicable par la profondeur de séquençage?
ord.fun(beetle.ps)
adonis2(beetle__bray$dist ~ depth + status, 
        data = beetle__bray$data,
        permutations=99999, method="bray") 
# Il y a une proportion importante de la variance expliquée par la profondeur!

# RDA conditionnelle
(ord_beetle <- capscale(beetle__bray$dist ~ status + Condition(depth), 
                  data = beetle__bray$data, distance = "bray"))

ord_beetle_anova <- anova(ord_beetle, permutations = how(nperm=99999))
ord_beetle_anova$R2 <- ord_beetle$CCA$tot.chi/ord_beetle$tot.chi
dbRDA_results[["beetle_cond"]] <- ord_beetle_anova

# transformation clr 
beetle_RCLR.ps <- tax_transform(beetle.ps, "rclr")
ord.fun(beetle_RCLR.ps, groupVar = "time", group = "2", d = "euclidean")
ord.fun(beetle_RCLR.ps, groupVar = "time", group = "7", d = "euclidean")
adonis2(beetle_RCLR_time2_euclidean$dist ~ status + depth, 
        data = beetle_RCLR_time2_euclidean$data,
        permutations=99999, method="euclidean")
adonis2(beetle_RCLR_time7_euclidean$dist ~ status + depth, 
        data = beetle_RCLR_time7_euclidean$data,
        permutations=99999, method="euclidean")

### ***************** ###
### BOMBUS TERRESTRIS ###
### ***************** ###

# PLOTTING BY TIME
ord.fun(bumble.ps, groupVar = "type", group = "Fèces")
ord.fun(bumble.ps, groupVar = "type", group = "Intestin")
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
bumble.col2 <- c("Fèces" = "darkkhaki", "Intestin" = "orangered")
bumble.p3 <- plot_PCoA(bumble__bray, group = "type") 
ord_Bterrestris_compart <- bumble.p3 +
  plot_annotation(title = "Comparaison des compartiments chez le bourdon.") &
  theme_light(base_family = "Baskerville", base_size = 26) &
  theme(plot.title = element_text(size = 32)) &
  theme(legend.position = 'bottom') &
  scale_colour_manual(values = bumble.col2) &
  scale_fill_manual(values = bumble.col2) &
  labs(colour = "Compartiment") 

ggsave("out/ord_Bterrestris_compart.png", plot = ord_Bterrestris_compart, 
       width = 27, height = 27, units = 'cm')

ord_beeID <- capscale(bumble__bray$dist ~ type,
                  data = bumble__bray$data)

anova_beeID <- anova(ord_beeID, permutations = how(nperm=9999))
anova_beeID$R2 <- ord_beeID$CCA$tot.chi/ord_beeID$tot.chi
(dbRDA_results[["Bterrestris_compart"]] <- anova_beeID)

### ********** ###
### HEDERA SPP ###
### ********** ###

# PLOTTING BY COMPARTMENT
# ord.fun(plant.ps, groupVar = "compart", group = "L")
ord.fun(plant.ps, groupVar = "compart", group = "Racines")
ord.fun(plant.ps, groupVar = "compart", group = "Sol")

#plant.p1 <- plot_PCoA(plant_compartL_bray, group = "ID", title = 'Feuilles') 
plant.p2 <- plot_PCoA(plant_compartR_bray, group = "ID", title = 'Racines') 
plant.p3 <- plot_PCoA(plant_compartS_bray, group = "ID", title = 'Sol') 

plant.col <- c("Parasités" = "orangered", "Non-parasités" = "darkolivegreen3")
ord_Hedera <- customPatchwork(plant.p2 + plant.p3) &
  plot_annotation(title = "Microbiote de Hedera en présence ou absence de Orobranche hederae.") &
  scale_colour_manual(values = plant.col) &
  scale_fill_manual(values = plant.col) &
  labs(colour = "Statut") 

ggsave("out/ord_Hedera.png", plot = ord_Hedera, 
       width = 48, height = 27, units = 'cm')

# REDUNDANCY ANALYSIS

# Racines
ord_R <- capscale(plant_compartR_bray$dist ~ ID, 
                       data = plant_compartR_bray$data)

anova_R <- anova(ord_R, permutations = how(nperm=99999))
anova_R$R2 <- ord_R$CCA$tot.chi/ord_R$tot.chi
dbRDA_results[["Hedera_R"]] <- anova_R

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
depthBox.plot <- ggplot(plant.ps@sam_data %>% data.frame %>% filter(compart != 'L'), 
                     aes(x = compart, y = depth, fill = ID)) +
  geom_boxplot() +
  labs(fill = "Compartiment", x = "",
       y = "Nombre de séquences",
       title = "A. Profondeur de séquençage des échantillons.") +
  scale_fill_manual(values = plant.col) +
  scale_x_discrete(labels = c("L" = "Feuilles", "R" = "Racines", "S" = "Sol")) +
  scale_y_continuous(limits = c(0, NA))

depthPCoA.plot <- plot_PCoA(plant_compartS_bray, group = "depth_log", ellipse = FALSE) +
  labs(colour = "Profondeur de \nséquençage (log)",
       title = 'B. Dissimilarité Bray-Curtis des échantillons du sol.') +
  theme(plot.title = element_text(size = 32, hjust = 0.5)) 

depth.plot <- depthBox.plot + depthPCoA.plot &
  theme_light(base_family = "Baskerville", base_size = 20) &
  theme(plot.title = element_text(size = 24)) &
  theme(legend.position = 'bottom')
  
ggsave("out/depth_plots.png", plot = depth.plot, 
       width = 48, height = 27, units = 'cm')

adonis2(plant_compartS_bray$dist ~ ID + depth, 
        data = plant_compartS_bray$data,
        permutations=99999, method="bray") 

#### dbRDA conditionnelle 
# RACINES
ord_R_cond <- capscale(plant_compartR_bray$dist ~ ID + Condition(depth), 
                  data = plant_compartR_bray$data)

anova_R_cond <- anova(ord_R_cond, permutations = how(nperm=99999))
anova_R_cond$R2 <- ord_R_cond$CCA$tot.chi/ord_R_cond$tot.chi
(dbRDA_results[["Hedera_R_cond"]] <- anova_R_cond)

# SOL
ord_S_cond <- capscale(plant_compartS_bray$dist ~ ID + Condition(depth), 
                data = plant_compartS_bray$data)

anova_S_cond <- anova(ord_S_cond, permutations = how(nperm=99999))
anova_S_cond$R2 <- ord_S_cond$CCA$tot.chi/ord_S_cond$tot.chi
(dbRDA_results[["Hedera_S_cond"]] <- anova_S_cond)

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
         R2 = round(R2, 4), 
         p = round(p, 4)) %>%  # Round values to two significant digits
  kable(caption = "Résultats des tests de permutations par analyses de redondance basée sur les distances.") %>% 
  kable_styling(bootstrap_options = "striped", full_width = F) %>% return


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

##### PROFONDEUR DE SÉQUENCAGe et T. molitor jour 7

ggplot(beetle.ps@sam_data %>% data.frame, 
                        aes(x = time, y = depth, fill = status)) +
  geom_boxplot() + theme_minimal() +
  labs(fill = "Compartiment", x = "",
       y = "Nombre de séquences totales.",
       title = "Profondeur de séquençage des échantillons T. molitor.") +
  scale_y_continuous(limits = c(0, NA))

plot_PCoA(beetle_time7_bray, group = "depth_log", ellipse = FALSE) +
  labs(colour = "Profondeur de \nséquençage (log)",
       title = 'Dissimilarité Bray-Curtis des échantillons du jour 7.') +
  theme(plot.title = element_text(size = 32, hjust = 0.5)) 

### RDA partielle
ord_7_cond <- capscale(beetle_time7_bray$dist ~ status + Condition(depth), 
                       data = beetle_time7_bray$data)

anova_7_cond <- anova(ord_7_cond, permutations = how(nperm=99999))
anova_7_cond$R2 <- ord_7_cond$CCA$tot.chi/ord_7_cond$tot.chi
(dbRDA_results[["Tmolitor_7_cond"]] <- anova_7_cond)
