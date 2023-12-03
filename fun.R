library(pacman); p_load(RColorBrewer, tidyverse, magrittr)
rename.fun <- function(x) {str_remove(x,".__")}

comm_barplot <- function(ps, byGroup, taxLvl, topN, host, facet = FALSE, title) {
  
  melt <- psmelt(ps)
  mycolors1 <- colorRampPalette(brewer.pal(8, "Set2"))(topN+1)
  
  # Choose top taxa
  topTaxa <- melt %>%  
    group_by(!!sym(taxLvl)) %>% # melt the table and group by tax level
    summarise(Abundance = sum(Abundance)) %>% # find most overall abundant taxa
    arrange(desc(Abundance)) %>%  # Species them 
    mutate(aggTaxo = as.factor(case_when( # aggTaxo will become the plot legend
      row_number()<=topN ~ !!sym(taxLvl), #+++ We'll need to manually order the species!
      row_number()>topN~'Others'))) %>%  # +1 to include the Others section!
    select(-Abundance)
  
  topTaxaLvls <- topTaxa$aggTaxo %>% 
    head(topN) %>% 
    as.character %>% sort %>% c(.,"Others")
  
  ### We create an object that will be fed into ggplot
  df_comm <- melt %>% 
    left_join(.,topTaxa, by = taxLvl) %>% # use topTaxa to join aggTaxo variable
    {
      if (facet != FALSE) {
        group_by(., aggTaxo, .data[[byGroup]], .data[[facet]])
      } else {
        group_by(., aggTaxo, .data[[byGroup]])
      }
    } %>%
    summarize(Abundance = sum(Abundance)) %>% # ...sums "Others" taxa by specified variables
    mutate(aggTaxo = factor(aggTaxo,# reorder 
                            levels = topTaxaLvls))
  
  ### Let's plot:
  p <- df_comm %>% ggplot(aes(x = .data[[byGroup]], y = Abundance, fill = aggTaxo)) +
    geom_bar(stat="identity", position="fill") + # change to "stack" to plot total counts
    labs(title=title, 
         fill=taxLvl, 
         x="", 
         y="Abondance relative") +
    scale_fill_manual(values = mycolors1,
                      breaks = topTaxaLvls) +
    theme_light(base_family="Baskerville", base_size = 20) +
    theme(axis.line = element_blank(), 
          axis.ticks = element_blank(),
          plot.margin = unit(c(1,0.5,1,2), "cm"), 
          plot.caption.position="plot",
          axis.title.y.left = element_text(vjust = 3), 
          legend.spacing.y = unit(0.7, 'cm')) +
    guides(fill = guide_legend(byrow = TRUE))
  
  if (facet!=FALSE) {
    p <- p + facet_grid(cols = vars(get(facet)),
               scales="free", # hide samples without counts
               space="free_x",
               switch="both") # constant bar width despite different # samples per group
  }
  return(p)
}

### Fonction pour plot la PCoA
plot_PCoA <- function(df, group, title = "", caption = "", ellipse = TRUE) {
  eig <- df[["eig"]]
  stEig1 <- eig[1]/sum(abs(eig))
  stEig2 <- eig[2]/sum(abs(eig))
  
  p <- ggplot(df[["data"]], 
         aes(PCOA1, PCOA2, colour = get(group))) +
    geom_point(size = 4) +
    guides(fill = 'none') + 
    labs(
      x = paste0("PCoA 1 [",round(stEig1*100,1),"% ]"),
      y = paste0("PCoA 2 [",round(stEig2*100,1),"% ]"),
      colour = group,
      title = title,
      caption = caption)
  if (ellipse == TRUE) {
    p + stat_ellipse(level=0.9, geom = "polygon", 
                   alpha = 0.18, aes(fill = get(group)))
  } else {p}
}

### Fonction de subsetting ####
subset.fun <- function(ps, groupVar, group) {
  if (groupVar != FALSE) {
    ps %>% 
      sample_data %$% get(groupVar) %in% group %>% 
      prune_samples(ps) %>% 
      prune_taxa(taxa_sums(.) > 0,.)
  } else {
    ps
  }
}

### Patchwork themes
customPatchwork <- function(plot) {
    plot + 
      plot_layout(guides = 'collect') & 
      theme_light(base_family = "Baskerville", base_size = 26) &
      theme(plot.title = element_text(size = 32, hjust = 0.5)) &
      theme(legend.position = 'bottom')
}
### PCoA selon la méthode vue en cours
ord.fun <- function(ps, groupVar = FALSE, group = FALSE, d = "bray") {
  # Subset the phyloseq object using this index
  subset.ps <- subset.fun(ps, groupVar, group)
  
  # Compute distance matrix
  dist.mx <- subset.ps %>%
    otu_table %>% t %>% # on veut une composante par échantillon...
    vegan::vegdist(d) 
  
  # Calculer les PCoA
  PCoA <- dist.mx %>%
    wcmdscale(k = 2, eig = T, add = "lingoes")
  
  # Plot data
  plot.df <- data.frame(PCOA1 = PCoA %>% scores %>% .[,1], 
             PCOA2 = PCoA %>% scores %>% .[,2]) %>%
    cbind(subset.ps %>% sample_data %>% data.frame)
  
  # Utiliser le nom de l'object PS pour extraire le nom du jeu de données
  set <- match.call()[["ps"]] %>% deparse %>% str_remove(".ps")
  
  # Exporter une liste contenant toutes les infos nécessaires
  exportList.fun(set, groupVar, group, dist.mx, PCoA$eig, plot.df, d)
}

### PCoA selon la méthode VST
vst.fun <- function(ps, groupVar = FALSE, group = FALSE) {
  subset.ps <- subset.fun(ps, groupVar, group)
  
  # Variance-stabilizing transformation
  vst.mx <- subset.ps %>% 
    phyloseq_to_deseq2(~1) %>% # DESeq2 object
    estimateSizeFactors(., geoMeans = apply(
      counts(.), 1, function(x) exp(sum(log(x[x>0]))/length(x)))) %>% 
    DESeq2::varianceStabilizingTransformation(blind=T) %>% # VST
    SummarizedExperiment::assay(.) %>% t %>% 
    { .[. < 0] <- 0; . } # Replace negative distances with 0s
  
  # Distance-based redundancy analysis
  dist.mx <- vegan::vegdist(vst.mx, distance = "bray")
  PCoA <- capscale(dist.mx~1)

  # Plot data
  plot.df <- data.frame(PCOA1 = PCoA %>% scores %$% sites %>% .[,1], 
             PCOA2 = PCoA %>% scores %$% sites %>% .[,2]) %>%
    cbind(subset.ps %>% sample_data %>% data.frame)

  set <- match.call()[["ps"]] %>% deparse %>% str_remove(".ps")
  exportList.fun(set, groupVar, group, dist.mx, PCoA$CA$eig, plot.df, "_vst")
}

# Générer une liste contenant la matrice de distances, les eigenvalues et 
# le jeu de données avec les composantes calculées (servira pour les plots)
exportList.fun <- function(set, groupVar, group, dist.mx, eig, plot.df, method){
  if(group==FALSE) {group=""}
  if(groupVar==FALSE) {groupVar=""}
  
  objectName <- paste0(set, "_", groupVar, group, "_", method)
  print(objectName) # sanity check
  
  assign(objectName, # name of list object
         list(dist = dist.mx,
              eig = eig,
              data = plot.df),
         envir = .GlobalEnv)
}
