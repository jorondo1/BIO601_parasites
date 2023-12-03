##### Parasitic plant microbiome #####

# libraries
library(ggplot2)
library(phyloseq)
library(dada2)
library(vegan)
library(lme4)
library(lmerTest)
library(DESeq2)
library(ALDEx2)
library(ggmap)
library(genefilter)
library(Biostrings)
library(ape)
library(ade4)
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(phytools)
library(dplyr)
library(reshape2)
library(broom)
library(seqinr)
library(psych)
library(picante)
library(ade4)
library(gplots)
library(corrplot)
library(SpiecEasi)
library(igraph)
library(decontam)

# functions #
st.err <- function(x) {
  sd(x)/sqrt(length(x))
}
is.finite.data.frame <- function(obj){
  sapply(obj,FUN = function(x) all(is.finite(x)))
} # Thanks joran @ stack exchange

rel.abund <- function(x){x / sum(x)}

# Define zero-tolerant function for calculating geometric means (courtesy of P.J. McMurdie) #
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# ggplot theme #
# *** Note: all figures made in R, arranged in Inkscape #
theme.gg<-theme(panel.background = element_rect(fill = "transparent", 
                colour = NA), plot.background = element_rect(colour = 'NA',                               fill = 'transparent'),panel.grid.major=element_line(color=NA),
                axis.line=element_line(size=1.2,color="black"),
                axis.ticks=element_line(color="black"),
                axis.text=element_text(color="black",size=26),
                axis.title=element_text(color="black",size=30),
                panel.grid.minor = element_line(colour = NA),
                legend.text=element_text(size=17),legend.key=element_rect(fill="white"),
                legend.title=element_text(size=18,face="bold"),
                axis.text.x=element_text(angle=60,size=26,vjust=0.5),
                panel.border = element_rect(colour = "black", fill=NA, size=3),
                legend.position = "none")

theme.gg.leg<-theme(panel.background = element_rect(fill = "transparent", 
                    colour = NA), plot.background = element_rect(colour = 'NA',                                                                                                 fill = 'transparent'),panel.grid.major=element_line(color=NA),
                    axis.line=element_line(size=1.2,color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(color="black",size=26),
                    axis.title=element_text(color="black",size=30),
                    panel.grid.minor = element_line(colour = NA),
                    legend.text=element_text(size=17),legend.key=element_rect(fill="white"),
                    legend.title=element_text(size=18,face="bold"),
                    axis.text.x=element_text(angle=60,size=26,vjust=0.5),
                    panel.border = element_rect(colour = "black", fill=NA, size=3))


# bacterial phyla colours
phylum.colors <- c("Acidobacteria" = "#fb8072", "Actinobacteria" = "#d95f02", "Bacteroidetes"="#1b9e77", 
                   "Armatimonadetes"="#fdbf6f", "Chlamydiae"="#CDB79E", "Chloroflexi"="#B22222", 
                   "Firmicutes"="#E69F00", "Parcubacteria"="#DCDCDC", "Planctomycetes"="#cab2d6",
                   "Proteobacteria"="#1f78b4", "Verrucomicrobia"="#ffff33", "Gemmatimonadetes"="#F0E68C", 
                   "Cyanobacteria/Chloroplast"="#548B54", "Spirochaetes"="#CC79A7","Ignavibacteriae"="#8B4726",
                   "Low Abundance"="#000000","candidate_division_WPS-2"="#1E90FF")

phylum.colors.grey <- c("Acidobacteria" = "#d7191c", "Actinobacteria" = "#abdda4", "Bacteroidetes"="#ffffbf", 
                   "Armatimonadetes"="#fdbf6f", "Chlamydiae"="#CDB79E", "Chloroflexi"="#B22222", 
                   "Firmicutes"="#E69F00", "Parcubacteria"="#DCDCDC", "Planctomycetes"="#cab2d6",
                   "Proteobacteria"="#2b83ba", "Verrucomicrobia"="#fdae61", "Gemmatimonadetes"="#F0E68C", 
                   "Cyanobacteria/Chloroplast"="#548B54", "Spirochaetes"="#CC79A7","Ignavibacteriae"="#8B4726",
                   "Low Abundance"="#000000","candidate_division_WPS-2"="#1E90FF")

sample.colors.alpha <- c("I.IR"="#00441b","I.L"="#006d2c","I.R"="#00441b","U.R"="#00441b","U.L"="#006d2c",
                   "P.R"="#3f007d","P.L"="#6a51a3","I.S"="#993404","U.S"="#662506")

sample.colors.beta <- c("I.IR"="#004529","I.L"="#004529","I.R"="#238443","U.R"="#addd8e","U.L"="#addd8e",
                         "P.R"="#7a0177","P.L"="#f768a1")

sample.colors.beta.grey <- c("I.IR"="#d7191c","I.L"="#d7191c","I.R"="#fdae61","U.R"="#ffffbf","U.L"="#ffffbf",
                        "P.R"="#2b83ba","P.L"="#abdda4")

sample.colors.alpha.grey <- c("I.IR"="#d7191c","I.L"="#d7191c","I.R"="#fdae61","U.R"="#ffffbf","U.L"="#ffffbf",
                             "P.R"="#2b83ba","P.L"="#abdda4")

### UC Berkeley map figure ###
map <- get_map(location = c(lon = -122.263725, lat = 37.871650),
               source = "google", maptype = "satellite", zoom = 16)

pdf(file="Map.pdf", width=10, height=10)
ggmap(map)
dev.off()

##### Build Phyloseq object #####

seqtab.total<-readRDS("spe_table.rds")
bac.tree<-read.tree(file= "PPM_bac.tre")

samples.out<-t(as.data.frame(sapply(strsplit(rownames(seqtab.total), "-"), `[`)))
meta.data<-as.data.frame(samples.out)
rownames(meta.data) <- rownames(seqtab.total)
colnames(meta.data)<-c("site","ID","type","replicate")
meta.data$IDtype<-interaction(meta.data$ID,meta.data$type)
meta.data$IDtypesite<-interaction(meta.data$IDtype,meta.data$site)
meta.data$species<-c(rep("Ha",9),rep("soil",2),rep("Oh",6),rep("Ha",6),rep("soil",2),
                     rep("Ha",9),rep("soil",2),rep("Oh",6),rep("Ha",6),rep("soil",2),
                     rep("Ha",9),rep("soil",2),rep("Oh",6),rep("Ha",6),rep("soil",2),
                     rep("Ha",9),rep("soil",2),rep("Oh",6),rep("Ha",6),rep("soil",1),
                     "mock","PAO1","undetermined","water")

# add PCR and DNA conc. #
meta.data$samp<-rownames(meta.data)
pcr.quant<-read.csv("sample_yields.csv")
pcr.quant$quant<-as.numeric(pcr.quant$quant)
pcr.quant$dna<-as.numeric(pcr.quant$dna)
meta.data<-merge(meta.data,pcr.quant,by.x="samp",by.y="sample")
rownames(meta.data)<-meta.data$samp

taxa<-readRDS("RDP_taxa.rds")

ps0 <- phyloseq(otu_table(seqtab.total, taxa_are_rows=FALSE), 
               sample_data(meta.data), 
               tax_table(taxa),
               phy_tree(bac.tree))

colnames(tax_table(ps0)) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus")

ps0<-subset_samples(ps0, site != "Undetermined") # remove undetermined reads


# remove Archaea, unknowns, chloroplasts #
ps = subset_taxa(ps0, Domain == "Bacteria") 
ps 

ps1 <- subset_taxa(ps, !is.na(Phylum)) 
ps1 

bad.class<-c("Chloroplast")
ps2 = subset_taxa(ps1, !G %in% bad.class)
ps2 

# caluclate # of usable reads #
sample_data(ps2)$UsableReads<-sample_sums(ps2)

readsumsdf = data.frame(nreads = sort(taxa_sums(ps2), TRUE), sorted = 1:ntaxa(ps2), 
                        type = "OTUs")
sum(readsumsdf$nreads) # total # of reads after removing unwanted samples/ASVs

# identifying contaminants #
head(sample_data(ps2))
contamdf.freq <- isContaminant(ps2, method="frequency", conc="quant")
table(contamdf.freq$contaminant) # 80 potential contaminants
head(which(contamdf.freq$contaminant)) # not particularly abundant

bad.asv<-contamdf.freq[contamdf.freq$contaminant=="TRUE",]
bad.asv<-bad.asv[bad.asv$p<0.05,]
bad.asv<-rownames(bad.asv)

good.asv <- setdiff(taxa_names(ps2), bad.asv)
ps3 <- prune_taxa(good.asv, ps2)

sum(taxa_sums(ps3))/sum(taxa_sums(ps2))

# obtain table of potential contaminants #
psCon <- prune_taxa(bad.asv, ps2)
psCon<-prune_taxa(taxa_sums(psCon)>0,psCon); ntaxa(psCon) 
tax.con<-as.data.frame(tax_table(psCon))
spe.con<-t(as.data.frame(otu_table(psCon)))
data.con<-merge(tax.con,spe.con,by=0)
data.con$sums<-rowSums(data.con[,c(8:109)])
write.csv(data.con,"contaminants.csv")

##### 1) Alpha diversity analyses #####
alpha.1<-as.data.frame(estimate_richness(ps3,measures=c("Observed","Shannon","InvSimpson")))
alpha.1.sd<-as(sample_data(ps3),'data.frame')
alpha.2<-cbind(alpha.1.sd,alpha.1)
alpha.2$E2<-alpha.2$InvSimpson/alpha.2$Observed
alpha.2$log.reads<-log(alpha.2$UsableReads)
alpha.2$IDtype<-factor(alpha.2$IDtype, 
                       levels = c("I.L","U.L","I.IR","I.R","U.R","P.L","P.R","I.S","U.S"),
                       ordered = T)

alpha.2$IDtypesite <- factor(alpha.2$IDtypesite,
                          levels = c("I.L.1","I.L.2","I.L.3","I.L.4",
                                     "U.L.1","U.L.2","U.L.3","U.L.4",
                                     "I.IR.1","I.IR.2","I.IR.3","I.IR.4",
                                     "I.R.1","I.R.2","I.R.3","I.R.4",
                                     "U.R.1","U.R.2","U.R.3","U.R.4",
                                     "P.L.1","P.L.2","P.L.3","P.L.4",
                                     "P.R.1","P.R.2","P.R.3","P.R.4",
                                     "I.S.1","I.S.2","I.S.3","I.S.4",
                                     "U.S.1","U.S.2","U.S.3","U.S.4"),
                          ordered = T)


# phylogenetic diversity #
samp<-as.data.frame(otu_table(ps3))
pd.total<-pd(samp,bac.tree)
alpha.2$pd<-pd.total$PD

# creating subsets of aplha diversity data #
a2<-alpha.2[alpha.2$type=="IR" | alpha.2$type=="R" | alpha.2$type=="S" | alpha.2$type=="L",]
a2$sample<-ifelse(a2$type=="IR","R",
                  ifelse(a2$type=="R","R",
                         ifelse(a2$type=="L","L",
                                ifelse(a2$type=="S","S","NA"))))
a2.leaf<-a2[a2$type=="L",]
a2.root<-a2[a2$type=="IR" | a2$type=="R",]
a2.root.soil<-a2[a2$type=="IR" | a2$type=="R" | a2$type=="S",]
a2.root.leaf<-a2[a2$type=="IR" | a2$type=="R" | a2$type=="L",]

# calculating mean and st.err alpha diversity values for each sample type
mean.alpha <- a2 %>% group_by(IDtypesite) %>% summarise_all(funs(mean))
se.alpha <- a2 %>% group_by(IDtypesite) %>% summarise_all(funs(st.err))
summ.alpha <- cbind(mean.alpha[,-c(2:7,13)],se.alpha[,-c(2:7,13)])
write.csv(summ.alpha,"alpha.summ.csv")

plot(log(Observed)~log.reads,data=a2)
plot(log(InvSimpson)~log.reads,data=a2)
plot(E2~log.reads,data=a2)

pdf(file="Alpha.obs.plantsamples.pdf",
    width=6,height=6,bg = "transparent")
ggplot(data=a2.root.leaf[a2.root.leaf$Observed<1300,], aes(x=IDtype, y=Observed, fill=IDtype)) +
  geom_jitter(alpha=.5,width=.3)+geom_boxplot(outlier.color = "transparent") +
  geom_boxplot() + scale_y_continuous(breaks=seq(0,1500,200)) +
  scale_fill_manual(values=sample.colors.beta) +
  theme.gg  
dev.off()

pdf(file="Alpha.InvS.plantsamples.pdf",
    width=6,height=6,bg = "transparent")
ggplot(data=a2.root.leaf, aes(x=IDtype, y=InvSimpson, fill=IDtype)) +
  geom_jitter(alpha=.5,width=.3)+geom_boxplot(outlier.color = "transparent") +
  geom_boxplot() + scale_y_continuous(breaks=seq(0,200,20)) +
  scale_fill_manual(values=sample.colors.beta) +
  theme.gg
dev.off()

pdf(file="Alpha.E2.plantsamples.pdf",
    width=6,height=6,bg = "transparent")
ggplot(data=a2.root.leaf, aes(x=IDtype, y=E2, fill=IDtype)) +
  geom_jitter(alpha=.5,width=.3)+geom_boxplot(outlier.color = "transparent") +
  geom_boxplot() + scale_y_continuous(breaks=seq(0,1,0.1)) +
  scale_fill_manual(values=sample.colors.beta) +
  theme.gg 
dev.off()

pdf(file="Alpha.PD.plantsamples.pdf",
    width=6,height=6,bg = "transparent")
ggplot(data=a2.root.leaf, aes(x=IDtype, y=pd, fill=IDtype)) +
  geom_jitter(alpha=.5,width=.3)+geom_boxplot(outlier.color = "transparent") +
  geom_boxplot() + scale_y_continuous(breaks=seq(0,300,25)) +
  scale_fill_manual(values=sample.colors.beta) +
  theme.gg 
dev.off()

# Leaf and root together - effect of species and tissue 
m1<-lmer(log(Observed) ~ species*sample + log.reads + dna + (1|sample:site) +(1|species:site) + (1|sample:species:site) + (1|site),data=a2.root.leaf,REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(log(InvSimpson) ~ species*sample + log.reads + dna + (1|sample:site) +(1|species:site) + (1|sample:species:site) + (1|site),data=a2.root.leaf,REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(log(InvSimpson) ~ sample + log.reads + (1|sample:site) +(1|species:site) + (1|sample:species:site) + (1|site),data=a2.root.leaf,REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(log(E2) ~ species*sample + log.reads + (1|sample:site) +(1|species:site) + (1|sample:species:site) + (1|site),data=a2.root.leaf,REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(log(E2) ~ log.reads + (1|sample:site) +(1|species:site) + (1|sample:species:site) + (1|site),data=a2.root.leaf,REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(log(pd) ~ species*sample + log.reads + (1|sample:site) +(1|species:site) + (1|sample:species:site) + (1|site),data=a2.root.leaf,REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
ranova(m1)
plot(m1)
hist(resid(m1))

# Leaf and root separately - effect of infection
m1<-lmer(log(Observed) ~ ID + log.reads + (1|ID:site) + (1|site), data=a2.root.leaf[a2.root.leaf$species == "Ha" & a2.root.leaf$sample == "L",],REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(log(Observed) ~ log.reads + (1|ID:site) + (1|site), data=a2.root.leaf[a2.root.leaf$species == "Ha" & a2.root.leaf$sample == "L",],REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(log(InvSimpson) ~ ID + log.reads + (1|ID:site) + (1|site), data=a2.root.leaf[a2.root.leaf$species == "Ha" & a2.root.leaf$sample == "L",],REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(log(InvSimpson) ~ log.reads + (1|ID:site) + (1|site), data=a2.root.leaf[a2.root.leaf$species == "Ha" & a2.root.leaf$sample == "L",],REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(log(E2) ~ ID + log.reads + (1|ID:site) + (1|site), data=a2.root.leaf[a2.root.leaf$species == "Ha" & a2.root.leaf$sample == "L",],REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(log(E2) ~ log.reads + (1|ID:site) + (1|site), data=a2.root.leaf[a2.root.leaf$species == "Ha" & a2.root.leaf$sample == "L",],REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(log(pd) ~ ID + log.reads + (1|ID:site) + (1|site), data=a2.root.leaf[a2.root.leaf$species == "Ha" & a2.root.leaf$sample == "L",],REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(log(pd) ~ log.reads + (1|ID:site) + (1|site), data=a2.root.leaf[a2.root.leaf$species == "Ha" & a2.root.leaf$sample == "L",],REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(log(Observed) ~ ID + log.reads + (1|ID:site) + (1|site), data=a2.root.leaf[a2.root.leaf$species == "Ha" & a2.root.leaf$sample == "R",],REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(log(Observed) ~ log.reads + (1|ID:site) + (1|site), data=a2.root.leaf[a2.root.leaf$species == "Ha" & a2.root.leaf$sample == "R",],REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(log(InvSimpson) ~ ID + log.reads + (1|ID:site) + (1|site), data=a2.root.leaf[a2.root.leaf$species == "Ha" & a2.root.leaf$sample == "R",],REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(log(InvSimpson) ~ 1 + (1|ID:site) + (1|site), data=a2.root.leaf[a2.root.leaf$species == "Ha" & a2.root.leaf$sample == "R",],REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(log(E2) ~ ID + log.reads + (1|ID:site) + (1|site), data=a2.root.leaf[a2.root.leaf$species == "Ha" & a2.root.leaf$sample == "R",],REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(log(E2) ~ 1 + (1|ID:site) + (1|site), data=a2.root.leaf[a2.root.leaf$species == "Ha" & a2.root.leaf$sample == "R",],REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(log(pd) ~ ID + log.reads + (1|ID:site) + (1|site), data=a2.root.leaf[a2.root.leaf$species == "Ha" & a2.root.leaf$sample == "R",],REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(log(pd) ~ log.reads + (1|ID:site) + (1|site), data=a2.root.leaf[a2.root.leaf$species == "Ha" & a2.root.leaf$sample == "R",],REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

# Are host and parasite leaf and root microbial diversity correlated?
ivy.iir<-a2[a2$IDtype == "I.IR", ]
colnames(ivy.iir)<-paste(colnames(ivy.iir),"IIR",sep = ".")
ivy.iur<-a2[a2$IDtype == "I.R", ]
colnames(ivy.iur)<-paste(colnames(ivy.iur),"IUR",sep = ".")
ivy.il<-a2[a2$IDtype == "I.L", ]
colnames(ivy.il)<-paste(colnames(ivy.il),"IL",sep = ".")
oro.r<-a2[a2$IDtype == "P.R", ]
colnames(oro.r)<-paste(colnames(oro.r),"PR",sep = ".")
oro.L<-a2[a2$IDtype == "P.L", ]
colnames(oro.L)<-paste(colnames(oro.L),"PL",sep = ".")
div<-cbind(ivy.iir[,c(8,10,11)],ivy.iur[,c(8,10,11)],ivy.il[,c(8,10,11)],oro.r[,c(8,10,11)],oro.L[,c(8,10,11)])
res <- rcorr(as.matrix(div))

corrplot(res$r, type="upper", 
         p.mat = res$P, sig.level = 0.05, insig = "blank")

plot(div$Observed.PR~div$InvSimpson.IUR)
plot(div$Observed.PR~div$InvSimpson.IIR)
plot(div$Observed.PL~div$InvSimpson.IL)

##### 2) Beta diversity analyses #####
# Standard thresholding: 1% occurence x 25 reads (throw out "non-reproducible" OTUs) #
ps4<-subset_samples(ps3, site != "Undetermined" & site != "mock" & site != "PAO1" & site != "water") # remove toothpick controls

threshold1<-kOverA(5,A=25) # set threshold values (require k samples with A reads)
ps5<-filter_taxa(ps4,threshold1,TRUE)
ntaxa(ps5) 
sum(taxa_sums(ps5))/sum(taxa_sums(ps4)) # common set of 2241 ASVs represents 78% of all observations
sample_data(ps5)$sample<-ifelse(sample_data(ps5)$type=="L" | sample_data(ps5)$type=="R" | sample_data(ps5)$type=="S",
                                as.character(sample_data(ps5)$type), "R")
sample_data(ps5)$site.rep<-interaction(sample_data(ps5)$site,sample_data(ps5)$replicate)
sample_data(ps5)$site.rep

sample.colors.site.rep <- c("1.1"="#bcbddc","1.2"="#6a51a3","1.3"="#3f007d",
                            "2.1"="#fcbba1","2.2"="#fb6a4a","2.3"="#a50f15",
                            "3.1"="#9ecae1","3.2"="#2171b5","3.3"="#08306b",
                            "4.1"="#a1d99b","4.2"="#238b45","4.3"="#00441b")
sample_data(ps5)

# subsets #
ps5.root<-subset_samples(ps5, type == "IR" | type == "R")
ps5.leaf<-subset_samples(ps5, type == "L")
ps5.leaf.root<-subset_samples(ps5, type == "IR" | type == "R" | type == "L")
ps5.parasite.root<-subset_samples(ps5, IDtype == "P.R")
ps5.parasite.leaf<-subset_samples(ps5, IDtype == "P.L")
ps5.host.root<-subset_samples(ps5, IDtype == "I.IR")
ps5.host.uroot<-subset_samples(ps5, IDtype == "I.R")
ps5.host.leaf<-subset_samples(ps5, IDtype == "I.L")
ps5.uninf.root<-subset_samples(ps5, IDtype == "U.R")
ps5.uninf.leaf<-subset_samples(ps5, IDtype == "U.L")
ps5.allsoil<-subset_samples(ps5, type == "S")
ps5.soil<-subset_samples(ps5, IDtype == "I.S")
ps5.usoil<-subset_samples(ps5, IDtype == "U.S")
ps6<-prune_taxa(taxa_sums(ps5)>0,ps5); ntaxa(ps6) 
ps6.root<-prune_taxa(taxa_sums(ps5.root)>0,ps5.root); ntaxa(ps6.root) 
ps6.leaf<-prune_taxa(taxa_sums(ps5.leaf)>0,ps5.leaf); ntaxa(ps6.leaf) 
ps6.leaf.root<-prune_taxa(taxa_sums(ps5.leaf.root)>0,ps5.leaf.root); ntaxa(ps6.leaf.root) 
ps6.parasite.root<-prune_taxa(taxa_sums(ps5.parasite.root)>0,ps5.parasite.root); ntaxa(ps6.parasite.root) 
ps6.parasite.leaf<-prune_taxa(taxa_sums(ps5.parasite.leaf)>0,ps5.parasite.leaf); ntaxa(ps6.parasite.leaf) 
ps6.host.root<-prune_taxa(taxa_sums(ps5.host.root)>0,ps5.host.root); ntaxa(ps6.host.root) 
ps6.host.uroot<-prune_taxa(taxa_sums(ps5.host.uroot)>0,ps5.host.uroot); ntaxa(ps6.host.uroot) 
ps6.host.leaf<-prune_taxa(taxa_sums(ps5.host.leaf)>0,ps5.host.leaf); ntaxa(ps6.host.leaf) 
ps6.uninf.root<-prune_taxa(taxa_sums(ps5.uninf.root)>0,ps5.uninf.root); ntaxa(ps6.uninf.root) 
ps6.uninf.leaf<-prune_taxa(taxa_sums(ps5.uninf.leaf)>0,ps5.uninf.leaf); ntaxa(ps6.uninf.leaf) 
ps6.soil<-prune_taxa(taxa_sums(ps5.soil)>0,ps5.soil); ntaxa(ps6.soil) 
ps6.usoil<-prune_taxa(taxa_sums(ps5.usoil)>0,ps5.usoil); ntaxa(ps6.usoil) 
ps6.allsoil<-prune_taxa(taxa_sums(ps5.allsoil)>0,ps5.allsoil); ntaxa(ps6.allsoil) 

## first change counts into relative abundance ##
ps6.ra = transform_sample_counts(ps6, rel.abund)
ps6.root.ra = transform_sample_counts(ps6.root, rel.abund)
ps6.leaf.ra = transform_sample_counts(ps6.leaf, rel.abund)
ps6.leaf.root.ra = transform_sample_counts(ps6.leaf.root, rel.abund)
ps6.parasite.root.ra = transform_sample_counts(ps6.parasite.root, rel.abund)
ps6.parasite.leaf.ra = transform_sample_counts(ps6.parasite.leaf, rel.abund)
ps6.host.root.ra = transform_sample_counts(ps6.host.root, rel.abund)
ps6.host.uroot.ra = transform_sample_counts(ps6.host.uroot, rel.abund)
ps6.host.leaf.ra = transform_sample_counts(ps6.host.leaf, rel.abund)
ps6.uninf.root.ra = transform_sample_counts(ps6.uninf.root, rel.abund)
ps6.uninf.leaf.ra = transform_sample_counts(ps6.uninf.leaf, rel.abund)
ps6.soil.ra = transform_sample_counts(ps6.soil, rel.abund)
ps6.usoil.ra = transform_sample_counts(ps6.usoil, rel.abund)
ps6.allsoil.ra = transform_sample_counts(ps6.allsoil, rel.abund)

##### 2.1) Procrustes tests ####
## need to get mean axis scores for each distance metric ##
out.wunifrac.ps6.soil <- ordinate(ps6.soil.ra, method = "MDS", distance = "wunifrac")
out.wunifrac.ps6.usoil <- ordinate(ps6.usoil.ra, method = "MDS", distance = "wunifrac")
out.wunifrac.ps6.pr <- ordinate(ps6.parasite.root.ra, method = "MDS", distance = "wunifrac")
out.wunifrac.ps6.pl <- ordinate(ps6.parasite.leaf.ra, method = "MDS", distance = "wunifrac")
out.wunifrac.ps6.iir <- ordinate(ps6.host.root.ra, method = "MDS", distance = "wunifrac")
out.wunifrac.ps6.iur <- ordinate(ps6.host.uroot.ra, method = "MDS", distance = "wunifrac")
out.wunifrac.ps6.il <- ordinate(ps6.host.leaf.ra, method = "MDS", distance = "wunifrac")
out.wunifrac.ps6.ur <- ordinate(ps6.uninf.root.ra, method = "MDS", distance = "wunifrac")
out.wunifrac.ps6.ul <- ordinate(ps6.uninf.leaf.ra, method = "MDS", distance = "wunifrac")

pr.iir.pro<-protest(out.wunifrac.ps6.iir$vectors[,1:2], 
        out.wunifrac.ps6.pr$vectors[,1:2], 
        permutations = 999) # 0.5572; 0.039

# Fig. S7 #
pdf(file="pro.errors.pr.iir.pdf",
    width=6,height=4,bg = "transparent")
plot(pr.iir.pro,kind=2,ylim=c(0,0.6))
dev.off()

dist.wu.pr <- phyloseq::distance(ps6.parasite.root.ra, method="wunifrac")
dend.pr<-as.dendrogram(hclust(dist.wu.pr, method = "average"))
labels(dend.pr) <- c("2-1", "2-3", "1-3", "1-1", "1-2", "3-2", "3-1", "3-3",
                     "4-3", "4-2", "2-2", "4-1")
dend.pr %>% set("labels_to_char") %>% labels

dist.wu.iir <- phyloseq::distance(ps6.host.root.ra, method="wunifrac")
dend.iir<-as.dendrogram(hclust(dist.wu.iir, method = "average"))
labels(dend.iir) <- c("1-1", "1-2", "1-3", "3-1", "2-1", "2-2", "4-2", "4-3",
                      "3-3", "4-1", "2-3", "3-2")
dend.iir %>% set("labels_to_char") %>% labels


pdf(file="tanglegram.pr.iir.pdf",
    width=6,height=4,bg = "transparent")
tanglegram(untangle(dend.pr,dend.iir,method = "step2side"),
           lwd=3,color_lines="#8c96c6",edge.lwd=3, lab.cex=1.5,
           common_subtrees_color_lines = FALSE, highlight_distinct_edges  = FALSE, 
           highlight_branches_lwd = FALSE)
dev.off()

iir.iur.pro<-protest(out.wunifrac.ps6.iir$vectors[,1:2], 
        out.wunifrac.ps6.iur$vectors[,1:2], 
        permutations = 999) # 0.739; 0.002

pdf(file="pro.errors.iur.iir.pdf",
    width=6,height=4,bg = "transparent")
plot(iir.iur.pro,kind=2,ylim=c(0,0.6))
dev.off()

dist.wu.iur <- phyloseq::distance(ps6.host.uroot.ra, method="wunifrac")
dend.iur<-as.dendrogram(hclust(dist.wu.iur, method = "average"))
labels(dend.iur) <- c("2-3","4-1","4-2","2-1","4-3","1-1","1-2","1-3",
                       "3-2","2-2","3-1","3-3")
dend.iur %>% set("labels_to_char") %>% labels


pdf(file="tanglegram.iur.iir.pdf",
    width=6,height=4,bg = "transparent")
tanglegram(untangle(dend.iur,dend.iir,method = "step2side"),
           lwd=3,color_lines="#8c96c6",edge.lwd=3, lab.cex=1.5,
           common_subtrees_color_lines = FALSE, highlight_distinct_edges  = FALSE, 
           highlight_branches_lwd = FALSE)
dev.off()


protest(out.wunifrac.ps6.iir$vectors[,1:2], 
        out.wunifrac.ps6.il$vectors[,1:2], 
        permutations = 999) # 0.739; 0.002

iir.pl.pro<-protest(out.wunifrac.ps6.iir$vectors[,1:2], 
        out.wunifrac.ps6.pl$vectors[,1:2], 
        permutations = 999) # 0.739; 0.002

il.pl.pro<-protest(out.wunifrac.ps6.il$vectors[,1:2], 
        out.wunifrac.ps6.pl$vectors[,1:2], 
        permutations = 999) # 0.739; 0.002

pdf(file="pro.errors.pl.il.pdf",
    width=6,height=4,bg = "transparent")
plot(il.pl.pro,kind=2,ylim=c(0,0.6))
dev.off()

dist.wu.pl <- phyloseq::distance(ps6.parasite.leaf.ra, method="wunifrac")
dend.pl<-as.dendrogram(hclust(dist.wu.pl, method = "average"))
labels(dend.pl) <- c("2-2","3-1","3-2","3-3","2-1","1-1","1-3","4-1",
                     "4-2","2-3","1-2","4-3")
dend.pl %>% set("labels_to_char") %>% labels

dist.wu.il <- phyloseq::distance(ps6.host.leaf.ra, method="wunifrac")
dend.il<-as.dendrogram(hclust(dist.wu.il, method = "average"))
labels(dend.il) <- c("3-3","2-1","4-1","3-2","1-1","1-3","1-2","4-3",
                     "2-2","4-2","2-3","3-1")
dend.il %>% set("labels_to_char") %>% labels

pdf(file="tanglegram.pl.il.pdf",
    width=6,height=4,bg = "transparent")
tanglegram(untangle(dend.pl,dend.il,method = "step2side"),
           lwd=3,color_lines="#8c96c6",edge.lwd=3, lab.cex=1.5,
           common_subtrees_color_lines = FALSE, highlight_distinct_edges  = FALSE, 
           highlight_branches_lwd = FALSE)
dev.off()

protest(out.wunifrac.ps6.il$vectors[,1:2], 
        out.wunifrac.ps6.pr$vectors[,1:2], 
        permutations = 999) # 0.739; 0.002

protest(out.wunifrac.ps6.il$vectors[,1:2], 
        out.wunifrac.ps6.iur$vectors[,1:2], 
        permutations = 999) # 0.739; 0.002

protest(out.wunifrac.ps6.pl$vectors[,1:2], 
        out.wunifrac.ps6.pr$vectors[,1:2], 
        permutations = 999) # 0.739; 0.002

protest(out.wunifrac.ps6.pl$vectors[,1:2], 
        out.wunifrac.ps6.iur$vectors[,1:2], 
        permutations = 999) # 0.739; 0.002

protest(out.wunifrac.ps6.iur$vectors[,1:2], 
        out.wunifrac.ps6.pr$vectors[,1:2], 
        permutations = 999) # 0.5604; 0.03

protest(out.wunifrac.ps6.pl$vectors[,1:2], 
        out.wunifrac.ps6.il$vectors[,1:2], 
        permutations = 999) # 0.739; 0.002

protest(out.wunifrac.ps6.ur$vectors[,1:2], 
        out.wunifrac.ps6.ul$vectors[,1:2], 
        permutations = 999) # 0.739; 0.002

# Infected Soil #
sum(out.wunifrac.ps6.soil$values$Eigenvalues[1:3])/sum(out.wunifrac.ps6.soil$values$Eigenvalues) # 0.86
a1.wu.ps6.soil<-as.data.frame(out.wunifrac.ps6.soil$vectors[,1])
a2.wu.ps6.soil<-as.data.frame(out.wunifrac.ps6.soil$vectors[,2])
a3.wu.ps6.soil<-as.data.frame(out.wunifrac.ps6.soil$vectors[,3])
pcoa.axes.wunifrac.ps6.soil<-cbind(a1.wu.ps6.soil,a2.wu.ps6.soil,a3.wu.ps6.soil)
colnames(pcoa.axes.wunifrac.ps6.soil)<-c("a1.wu","a2.wu","a3.wu")
pcoa.axes.wunifrac.ps6.soil$IDtypesite<-as.data.frame(sample_data(ps6.soil))$IDtypesite

mean.PCoA1.wu.ps6.soil<-aggregate(pcoa.axes.wunifrac.ps6.soil$a1.wu ~ pcoa.axes.wunifrac.ps6.soil$IDtypesite, pcoa.axes.wunifrac.ps6.soil, mean)
mean.PCoA2.wu.ps6.soil<-aggregate(pcoa.axes.wunifrac.ps6.soil$a2.wu ~ pcoa.axes.wunifrac.ps6.soil$IDtypesite, pcoa.axes.wunifrac.ps6.soil, mean)
mean.PCoA3.wu.ps6.soil<-aggregate(pcoa.axes.wunifrac.ps6.soil$a3.wu ~ pcoa.axes.wunifrac.ps6.soil$IDtypesite, pcoa.axes.wunifrac.ps6.soil, mean)
soil.pcoa<-cbind(mean.PCoA1.wu.ps6.soil,mean.PCoA2.wu.ps6.soil[,2],mean.PCoA1.wu.ps6.soil[,2])
colnames(soil.pcoa)<-cbind("site","a1.soil","a2.soil","a3.soil")
euc.dist.wu.soil<-as.matrix(dist(soil.pcoa[,c(2:4)]))

# Uninfected Soil #
sum(out.wunifrac.ps6.usoil$values$Eigenvalues[1:3])/sum(out.wunifrac.ps6.usoil$values$Eigenvalues) # 0.85
a1.wu.ps6.usoil<-as.data.frame(out.wunifrac.ps6.usoil$vectors[,1])
a2.wu.ps6.usoil<-as.data.frame(out.wunifrac.ps6.usoil$vectors[,2])
a3.wu.ps6.usoil<-as.data.frame(out.wunifrac.ps6.usoil$vectors[,3])
pcoa.axes.wunifrac.ps6.usoil<-cbind(a1.wu.ps6.usoil,a2.wu.ps6.usoil,a3.wu.ps6.usoil)
colnames(pcoa.axes.wunifrac.ps6.usoil)<-c("a1.wu","a2.wu","a3.wu")
pcoa.axes.wunifrac.ps6.usoil$IDtypesite<-as.data.frame(sample_data(ps6.usoil))$IDtypesite

mean.PCoA1.wu.ps6.usoil<-aggregate(pcoa.axes.wunifrac.ps6.usoil$a1.wu ~ pcoa.axes.wunifrac.ps6.usoil$IDtypesite, pcoa.axes.wunifrac.ps6.usoil, mean)
mean.PCoA2.wu.ps6.usoil<-aggregate(pcoa.axes.wunifrac.ps6.usoil$a2.wu ~ pcoa.axes.wunifrac.ps6.usoil$IDtypesite, pcoa.axes.wunifrac.ps6.usoil, mean)
mean.PCoA3.wu.ps6.usoil<-aggregate(pcoa.axes.wunifrac.ps6.usoil$a3.wu ~ pcoa.axes.wunifrac.ps6.usoil$IDtypesite, pcoa.axes.wunifrac.ps6.usoil, mean)
usoil.pcoa<-cbind(mean.PCoA1.wu.ps6.usoil,mean.PCoA2.wu.ps6.usoil[,2],mean.PCoA1.wu.ps6.usoil[,2])
colnames(usoil.pcoa)<-cbind("site","a1.usoil","a2.usoil","a3.usoil")
euc.dist.wu.usoil<-as.matrix(dist(usoil.pcoa[,c(2:4)]))

# Parasite root #
sum(out.wunifrac.ps6.pr$values$Eigenvalues[1:3])/sum(out.wunifrac.ps6.pr$values$Eigenvalues) # 0.74
a1.wu.ps6.pr<-as.data.frame(out.wunifrac.ps6.pr$vectors[,1])
a2.wu.ps6.pr<-as.data.frame(out.wunifrac.ps6.pr$vectors[,2])
a3.wu.ps6.pr<-as.data.frame(out.wunifrac.ps6.pr$vectors[,3])
pcoa.axes.wunifrac.ps6.pr<-cbind(a1.wu.ps6.pr,a2.wu.ps6.pr,a3.wu.ps6.pr)
colnames(pcoa.axes.wunifrac.ps6.pr)<-c("a1.wu","a2.wu","a3.wu")
pcoa.axes.wunifrac.ps6.pr$IDtypesite<-as.data.frame(sample_data(ps6.parasite.root.ra))$IDtypesite

mean.PCoA1.wu.ps6.pr<-aggregate(pcoa.axes.wunifrac.ps6.pr$a1.wu ~ pcoa.axes.wunifrac.ps6.pr$IDtypesite, pcoa.axes.wunifrac.ps6.pr, mean)
mean.PCoA2.wu.ps6.pr<-aggregate(pcoa.axes.wunifrac.ps6.pr$a2.wu ~ pcoa.axes.wunifrac.ps6.pr$IDtypesite, pcoa.axes.wunifrac.ps6.pr, mean)
mean.PCoA3.wu.ps6.pr<-aggregate(pcoa.axes.wunifrac.ps6.pr$a3.wu ~ pcoa.axes.wunifrac.ps6.pr$IDtypesite, pcoa.axes.wunifrac.ps6.pr, mean)
pr.pcoa<-cbind(mean.PCoA1.wu.ps6.pr,mean.PCoA2.wu.ps6.pr[,2],mean.PCoA1.wu.ps6.pr[,2])
colnames(pr.pcoa)<-cbind("site","a1.pr","a2.pr","a3.pr")
euc.dist.wu.pr<-as.matrix(dist(pr.pcoa[,c(2:4)]))

# Parasite leaf #
sum(out.wunifrac.ps6.pl$values$Eigenvalues[1:3])/sum(out.wunifrac.ps6.pl$values$Eigenvalues) # 0.74
a1.wu.ps6.pl<-as.data.frame(out.wunifrac.ps6.pl$vectors[,1])
a2.wu.ps6.pl<-as.data.frame(out.wunifrac.ps6.pl$vectors[,2])
a3.wu.ps6.pl<-as.data.frame(out.wunifrac.ps6.pl$vectors[,3])
pcoa.axes.wunifrac.ps6.pl<-cbind(a1.wu.ps6.pl,a2.wu.ps6.pl,a3.wu.ps6.pl)
colnames(pcoa.axes.wunifrac.ps6.pl)<-c("a1.wu","a2.wu","a3.wu")
pcoa.axes.wunifrac.ps6.pl$IDtypesite<-as.data.frame(sample_data(ps6.parasite.leaf.ra))$IDtypesite

mean.PCoA1.wu.ps6.pl<-aggregate(pcoa.axes.wunifrac.ps6.pl$a1.wu ~ pcoa.axes.wunifrac.ps6.pl$IDtypesite, pcoa.axes.wunifrac.ps6.pl, mean)
mean.PCoA2.wu.ps6.pl<-aggregate(pcoa.axes.wunifrac.ps6.pl$a2.wu ~ pcoa.axes.wunifrac.ps6.pl$IDtypesite, pcoa.axes.wunifrac.ps6.pl, mean)
mean.PCoA3.wu.ps6.pl<-aggregate(pcoa.axes.wunifrac.ps6.pl$a3.wu ~ pcoa.axes.wunifrac.ps6.pl$IDtypesite, pcoa.axes.wunifrac.ps6.pl, mean)
pl.pcoa<-cbind(mean.PCoA1.wu.ps6.pl,mean.PCoA2.wu.ps6.pl[,2],mean.PCoA1.wu.ps6.pl[,2])
colnames(pl.pcoa)<-cbind("site","a1.pl","a2.pl","a3.pl")
euc.dist.wu.pl<-as.matrix(dist(pl.pcoa[,c(2:4)]))

# Host root #
sum(out.wunifrac.ps6.iir$values$Eigenvalues[1:3])/sum(out.wunifrac.ps6.iir$values$Eigenvalues) # 0.78
a1.wu.ps6.iir<-as.data.frame(out.wunifrac.ps6.iir$vectors[,1])
a2.wu.ps6.iir<-as.data.frame(out.wunifrac.ps6.iir$vectors[,2])
a3.wu.ps6.iir<-as.data.frame(out.wunifrac.ps6.iir$vectors[,3])
pcoa.axes.wunifrac.ps6.iir<-cbind(a1.wu.ps6.iir,a2.wu.ps6.iir,a3.wu.ps6.iir)
colnames(pcoa.axes.wunifrac.ps6.iir)<-c("a1.wu","a2.wu","a3.wu")
pcoa.axes.wunifrac.ps6.iir$IDtypesite<-as.data.frame(sample_data(ps6.host.root.ra))$IDtypesite

mean.PCoA1.wu.ps6.iir<-aggregate(pcoa.axes.wunifrac.ps6.iir$a1.wu ~ pcoa.axes.wunifrac.ps6.iir$IDtypesite, pcoa.axes.wunifrac.ps6.iir, mean)
mean.PCoA2.wu.ps6.iir<-aggregate(pcoa.axes.wunifrac.ps6.iir$a2.wu ~ pcoa.axes.wunifrac.ps6.iir$IDtypesite, pcoa.axes.wunifrac.ps6.iir, mean)
mean.PCoA3.wu.ps6.iir<-aggregate(pcoa.axes.wunifrac.ps6.iir$a3.wu ~ pcoa.axes.wunifrac.ps6.iir$IDtypesite, pcoa.axes.wunifrac.ps6.iir, mean)
iir.pcoa<-cbind(mean.PCoA1.wu.ps6.iir,mean.PCoA2.wu.ps6.iir[,2],mean.PCoA1.wu.ps6.iir[,2])
colnames(iir.pcoa)<-cbind("site","a1.iir","a2.iir","a3.iir")
euc.dist.wu.iir<-as.matrix(dist(iir.pcoa[,c(2:4)]))

# Host uninfected root #
sum(out.wunifrac.ps6.iur$values$Eigenvalues[1:3])/sum(out.wunifrac.ps6.iur$values$Eigenvalues) # 0.75
a1.wu.ps6.iur<-as.data.frame(out.wunifrac.ps6.iur$vectors[,1])
a2.wu.ps6.iur<-as.data.frame(out.wunifrac.ps6.iur$vectors[,2])
a3.wu.ps6.iur<-as.data.frame(out.wunifrac.ps6.iur$vectors[,3])
pcoa.axes.wunifrac.ps6.iur<-cbind(a1.wu.ps6.iur,a2.wu.ps6.iur,a3.wu.ps6.iur)
colnames(pcoa.axes.wunifrac.ps6.iur)<-c("a1.wu","a2.wu","a3.wu")
pcoa.axes.wunifrac.ps6.iur$IDtypesite<-as.data.frame(sample_data(ps6.host.uroot.ra))$IDtypesite

mean.PCoA1.wu.ps6.iur<-aggregate(pcoa.axes.wunifrac.ps6.iur$a1.wu ~ pcoa.axes.wunifrac.ps6.iur$IDtypesite, pcoa.axes.wunifrac.ps6.iur, mean)
mean.PCoA2.wu.ps6.iur<-aggregate(pcoa.axes.wunifrac.ps6.iur$a2.wu ~ pcoa.axes.wunifrac.ps6.iur$IDtypesite, pcoa.axes.wunifrac.ps6.iur, mean)
mean.PCoA3.wu.ps6.iur<-aggregate(pcoa.axes.wunifrac.ps6.iur$a3.wu ~ pcoa.axes.wunifrac.ps6.iur$IDtypesite, pcoa.axes.wunifrac.ps6.iur, mean)
iur.pcoa<-cbind(mean.PCoA1.wu.ps6.iur,mean.PCoA2.wu.ps6.iur[,2],mean.PCoA1.wu.ps6.iur[,2])
colnames(iur.pcoa)<-cbind("site","a1.iur","a2.iur","a3.iur")
euc.dist.wu.iur<-as.matrix(dist(iur.pcoa[,c(2:4)]))

# Host leaf #
sum(out.wunifrac.ps6.il$values$Eigenvalues[1:3])/sum(out.wunifrac.ps6.il$values$Eigenvalues) # 0.77
a1.wu.ps6.il<-as.data.frame(out.wunifrac.ps6.il$vectors[,1])
a2.wu.ps6.il<-as.data.frame(out.wunifrac.ps6.il$vectors[,2])
a3.wu.ps6.il<-as.data.frame(out.wunifrac.ps6.il$vectors[,3])
pcoa.axes.wunifrac.ps6.il<-cbind(a1.wu.ps6.il,a2.wu.ps6.il,a3.wu.ps6.il)
colnames(pcoa.axes.wunifrac.ps6.il)<-c("a1.wu","a2.wu","a3.wu")
pcoa.axes.wunifrac.ps6.il$IDtypesite<-as.data.frame(sample_data(ps6.host.leaf.ra))$IDtypesite

mean.PCoA1.wu.ps6.il<-aggregate(pcoa.axes.wunifrac.ps6.il$a1.wu ~ pcoa.axes.wunifrac.ps6.il$IDtypesite, pcoa.axes.wunifrac.ps6.il, mean)
mean.PCoA2.wu.ps6.il<-aggregate(pcoa.axes.wunifrac.ps6.il$a2.wu ~ pcoa.axes.wunifrac.ps6.il$IDtypesite, pcoa.axes.wunifrac.ps6.il, mean)
mean.PCoA3.wu.ps6.il<-aggregate(pcoa.axes.wunifrac.ps6.il$a3.wu ~ pcoa.axes.wunifrac.ps6.il$IDtypesite, pcoa.axes.wunifrac.ps6.il, mean)
il.pcoa<-cbind(mean.PCoA1.wu.ps6.il,mean.PCoA2.wu.ps6.il[,2],mean.PCoA1.wu.ps6.il[,2])
colnames(il.pcoa)<-cbind("site","a1.il","a2.il","a3.il")
euc.dist.wu.il<-as.matrix(dist(il.pcoa[,c(2:4)]))

# Uninfected root #
sum(out.wunifrac.ps6.ur$values$Eigenvalues[1:3])/sum(out.wunifrac.ps6.ur$values$Eigenvalues) # 0.75
a1.wu.ps6.ur<-as.data.frame(out.wunifrac.ps6.ur$vectors[,1])
a2.wu.ps6.ur<-as.data.frame(out.wunifrac.ps6.ur$vectors[,2])
a3.wu.ps6.ur<-as.data.frame(out.wunifrac.ps6.ur$vectors[,3])
pcoa.axes.wunifrac.ps6.ur<-cbind(a1.wu.ps6.ur,a2.wu.ps6.ur,a3.wu.ps6.ur)
colnames(pcoa.axes.wunifrac.ps6.ur)<-c("a1.wu","a2.wu","a3.wu")
pcoa.axes.wunifrac.ps6.ur$IDtypesite<-as.data.frame(sample_data(ps6.uninf.root.ra))$IDtypesite

mean.PCoA1.wu.ps6.ur<-aggregate(pcoa.axes.wunifrac.ps6.ur$a1.wu ~ pcoa.axes.wunifrac.ps6.ur$IDtypesite, pcoa.axes.wunifrac.ps6.ur, mean)
mean.PCoA2.wu.ps6.ur<-aggregate(pcoa.axes.wunifrac.ps6.ur$a2.wu ~ pcoa.axes.wunifrac.ps6.ur$IDtypesite, pcoa.axes.wunifrac.ps6.ur, mean)
mean.PCoA3.wu.ps6.ur<-aggregate(pcoa.axes.wunifrac.ps6.ur$a3.wu ~ pcoa.axes.wunifrac.ps6.ur$IDtypesite, pcoa.axes.wunifrac.ps6.ur, mean)
ur.pcoa<-cbind(mean.PCoA1.wu.ps6.ur,mean.PCoA2.wu.ps6.ur[,2],mean.PCoA1.wu.ps6.ur[,2])
colnames(ur.pcoa)<-cbind("site","a1.ur","a2.ur","a3.ur")
euc.dist.wu.ur<-as.matrix(dist(ur.pcoa[,c(2:4)]))

# Uninfected leaf #
sum(out.wunifrac.ps6.ul$values$Eigenvalues[1:3])/sum(out.wunifrac.ps6.ul$values$Eigenvalues) # 0.81
a1.wu.ps6.ul<-as.data.frame(out.wunifrac.ps6.ul$vectors[,1])
a2.wu.ps6.ul<-as.data.frame(out.wunifrac.ps6.ul$vectors[,2])
a3.wu.ps6.ul<-as.data.frame(out.wunifrac.ps6.ul$vectors[,3])
pcoa.axes.wunifrac.ps6.ul<-cbind(a1.wu.ps6.ul,a2.wu.ps6.ul,a3.wu.ps6.ul)
colnames(pcoa.axes.wunifrac.ps6.ul)<-c("a1.wu","a2.wu","a3.wu")
pcoa.axes.wunifrac.ps6.ul$IDtypesite<-as.data.frame(sample_data(ps6.uninf.leaf.ra))$IDtypesite

mean.PCoA1.wu.ps6.ul<-aggregate(pcoa.axes.wunifrac.ps6.ul$a1.wu ~ pcoa.axes.wunifrac.ps6.ul$IDtypesite, pcoa.axes.wunifrac.ps6.ul, mean)
mean.PCoA2.wu.ps6.ul<-aggregate(pcoa.axes.wunifrac.ps6.ul$a2.wu ~ pcoa.axes.wunifrac.ps6.ul$IDtypesite, pcoa.axes.wunifrac.ps6.ul, mean)
mean.PCoA3.wu.ps6.ul<-aggregate(pcoa.axes.wunifrac.ps6.ul$a3.wu ~ pcoa.axes.wunifrac.ps6.ul$IDtypesite, pcoa.axes.wunifrac.ps6.ul, mean)
ul.pcoa<-cbind(mean.PCoA1.wu.ps6.ul,mean.PCoA2.wu.ps6.ul[,2],mean.PCoA1.wu.ps6.ul[,2])
colnames(ul.pcoa)<-cbind("site","a1.ul","a2.ul","a3.ul")
euc.dist.wu.ul<-as.matrix(dist(ul.pcoa[,c(2:4)]))

protest(iir.pcoa[,c(2:3)], 
        soil.pcoa[,c(2:3)], 
        permutations = 999) # 0.9118; 0.042

protest(iur.pcoa[,c(2:3)], 
        soil.pcoa[,c(2:3)], 
        permutations = 999) # 0.9118; 0.042

protest(il.pcoa[,c(2:3)], 
        soil.pcoa[,c(2:3)], 
        permutations = 999) # 0.6562; 0.58

protest(pr.pcoa[,c(2:3)], 
        soil.pcoa[,c(2:3)], 
        permutations = 999) # 0.952; 0.08

protest(pl.pcoa[,c(2:3)], 
        soil.pcoa[,c(2:3)], 
        permutations = 999) # 0.9118; 0.042

protest(ur.pcoa[,c(2:3)], 
        usoil.pcoa[,c(2:3)], 
        permutations = 999) # 0.9118; 0.042

protest(ul.pcoa[,c(2:3)], 
        usoil.pcoa[,c(2:3)], 
        permutations = 999) # 0.9118; 0.042

#### 2.1.1) Procrustes: leave taxa out analysis ####
# First generate phyloseq objects of parasite root and leaf microbiomes#
ps5.parasite.root<-subset_samples(ps5, IDtype == "P.R")
ps5.parasite.leaf<-subset_samples(ps5, IDtype == "P.L")
ps5.host.root<-subset_samples(ps5, IDtype == "I.IR")
ps5.host.uroot<-subset_samples(ps5, IDtype == "I.R")

ps6.parasite.root<-prune_taxa(taxa_sums(ps5.parasite.root)>0,ps5.parasite.root); ntaxa(ps6.parasite.root) 
ps6.parasite.leaf<-prune_taxa(taxa_sums(ps5.parasite.leaf)>0,ps5.parasite.leaf); ntaxa(ps6.parasite.leaf) 
ps6.host.root<-prune_taxa(taxa_sums(ps5.host.root)>0,ps5.host.root); ntaxa(ps6.host.root) 
ps6.host.uroot<-prune_taxa(taxa_sums(ps5.host.uroot)>0,ps5.host.uroot); ntaxa(ps6.host.uroot) 

ps6.parasite.root.ra = transform_sample_counts(ps6.parasite.root, rel.abund)
ps6.parasite.leaf.ra = transform_sample_counts(ps6.parasite.leaf, rel.abund)
ps6.host.root.ra = transform_sample_counts(ps6.host.root, rel.abund)
ps6.host.uroot.ra = transform_sample_counts(ps6.host.uroot, rel.abund)

out.wunifrac.pr <- phyloseq::ordinate(ps6.parasite.root.ra, method = "MDS", distance = "wunifrac")
out.wunifrac.pl <- phyloseq::ordinate(ps6.parasite.leaf.ra, method = "MDS", distance = "wunifrac")
out.wunifrac.iir <- phyloseq::ordinate(ps6.host.root.ra, method = "MDS", distance = "wunifrac")
out.wunifrac.iur <- phyloseq::ordinate(ps6.host.uroot.ra, method = "MDS", distance = "wunifrac")

# IIR versus PR: compare goodness of fit with taxa missing 
pro.pr.iir<-protest(out.wunifrac.iir$vectors[,1:2], 
                    out.wunifrac.pr$vectors[,1:2], 
                    permutations = 999) # 0.6895; 0.038
plot(pro.pr.iir)
plot(pro.pr.iir,kind=2)
taxa.pr<-as.data.frame(tax_table(ps6.parasite.root.ra))
phyla.pr<-as.character(unique(taxa.pr$Phylum))

phylum.results<-data.frame("t0"=numeric(),"pval"=numeric(),"numer.taxa.exluded"=numeric())
for (name in phyla.pr) 
{
  ps.phy<-subset_taxa(ps6.parasite.root, Phylum!= name)
  ps.phy.diff<-dim(otu_table(ps6.parasite.root))[2] - dim(otu_table(ps.phy))[2]
  ps.phy.ra<-transform_sample_counts(ps.phy, rel.abund)
  wu<-ordinate(ps.phy.ra, method = "MDS", distance = "wunifrac")
  pro<-protest(out.wunifrac.iir$vectors[,1:2], 
               wu$vectors[,1:2], permutations = 999)
  this.result<-c(pro$t0,pro$signif,ps.phy.diff)
  phylum.results<-rbind(phylum.results,this.result)
}

phylum.results$taxa<-phyla.pr
colnames(phylum.results)[1:3]<-c("t0","pval","ASVs.excluded")
phylum.results$deltat0<-phylum.results$t0-pro.pr.iir$t0
write.csv(phylum.results, "Manuscript/leave_one_out.phy.csv")

pdf(file="leave_one_out.phy.pdf",
    width=8,height=10,bg = "transparent")
ggplot(data=phylum.results, aes(x=taxa, y=t0)) +
  geom_point(aes(fill=taxa),size=12,color="black") +
  geom_hline(yintercept=pro.pr.iir$t0, linetype="dashed", color = "red")+
  scale_fill_manual(values=phylum.colors) +
  scale_y_continuous(breaks=seq(0.5,0.62,0.01))+
  theme.gg
dev.off()

# Order level #
prot.ord<-droplevels(unique(taxa.pr[which(taxa.pr$Phylum=="Proteobacteria"),4]))
prot.ord<-prot.ord[!is.na(prot.ord)]
actin.ord<-droplevels(unique(taxa.pr[which(taxa.pr$Phylum=="Actinobacteria"),4]))
actin.ord<-actin.ord[!is.na(actin.ord)]
planc.ord<-droplevels(unique(taxa.pr[which(taxa.pr$Phylum=="Planctomycetes"),4]))
planc.ord<-planc.ord[!is.na(planc.ord)]
arm.ord<-droplevels(unique(taxa.pr[which(taxa.pr$Phylum=="Armatimonadetes"),4]))
arm.ord<-arm.ord[!is.na(arm.ord)]
bac.ord<-droplevels(unique(taxa.pr[which(taxa.pr$Phylum=="Bacteroidetes"),4]))
bac.ord<-bac.ord[!is.na(bac.ord)]
acid.ord<-droplevels(unique(taxa.pr[which(taxa.pr$Phylum=="Acidobacteria"),4]))
acid.ord<-acid.ord[!is.na(acid.ord)]
order.list<-as.character(unlist(list(prot.ord,actin.ord,planc.ord,arm.ord,bac.ord,acid.ord)))

order.results<-data.frame("t0"=numeric(),"pval"=numeric(),"numer.taxa.exluded"=numeric())
for (name in order.list) 
{
  ps.ord.rem<-subset_taxa(ps6.parasite.root, Order == name)
  keep.asv<-setdiff(taxa_names(ps6.parasite.root),taxa_names(ps.ord.rem))
  ps.ord<-prune_taxa(keep.asv, ps6.parasite.root)
  ps.ord.diff<-ntaxa(ps6.parasite.root) - ntaxa(ps.ord)
  ps.ord.ra<-transform_sample_counts(ps.ord, rel.abund)
  wu<-ordinate(ps.ord.ra, method = "MDS", distance = "wunifrac")
  pro<-protest(out.wunifrac.iir$vectors[,1:2], 
               wu$vectors[,1:2], permutations = 999)
  this.result<-c(pro$t0,pro$signif,ps.ord.diff)
  order.results<-rbind(order.results,this.result)
}


order.results$taxa<-order.list
order.results$phylum<-c(rep("Proteobacteria",length(prot.ord)),rep("Actinobacteria",length(actin.ord)),
                        rep("Plancotmycetes",length(planc.ord)),rep("Armatimonadetes",length(arm.ord)),
                        rep("Bacteroidetes",length(bac.ord)),rep("Acidobacteria",length(acid.ord)))

order.results$taxa <- factor(order.results$taxa,
                             levels = order.list,ordered = TRUE)
colnames(order.results)[1:3]<-c("t0","pval","ASVs.excluded")
order.results$deltat0<-order.results$t0-pro.pr.iir$t0
write.csv(order.results, "Manuscript/leave_one_out.ord.csv")

pdf(file="leave_one_out.ord.pdf",
    width=28,height=20,bg = "transparent")
ggplot(data=order.results, aes(x=taxa, y=t0)) +
  geom_point(aes(fill=phylum),size=12,color="black") +
  geom_hline(yintercept=pro.pr.iir$t0, linetype="dashed", color = "red")+
  scale_fill_manual(values=phylum.colors) +
  scale_y_continuous(breaks=seq(0.2,0.9,0.01))+
  theme.gg
dev.off()

# sanity check #
ps.ord.rem<-subset_taxa(ps6.parasite.root, Order == "Flavobacteriales")
keep.asv<-setdiff(taxa_names(ps6.parasite.root),taxa_names(ps.ord.rem))
ps.ord<-prune_taxa(keep.asv, ps6.parasite.root)
ps.ord.ra<-transform_sample_counts(ps.ord, rel.abund)
wu<-ordinate(ps.ord.ra, method = "MDS", distance = "wunifrac")
pro<-protest(out.wunifrac.iir$vectors[,1:2], 
             wu$vectors[,1:2], permutations = 999) # matches with order.results
pro$t0

# IIR versus PL: compare goodness of fit with taxa missing 
pro.pl.iir<-protest(out.wunifrac.iir$vectors[,1:2], 
                    out.wunifrac.pl$vectors[,1:2], 
                    permutations = 999) # 0.6095; 0.011
plot(pro.pl.iir)
plot(pro.pl.iir,kind=2)
taxa.pl<-as.data.frame(tax_table(ps6.parasite.leaf.ra))
phyla.pl<-as.character(unique(taxa.pl$Phylum))

phylum.results<-data.frame("t0"=numeric(),"pval"=numeric(),"numer.taxa.exluded"=numeric())
for (name in phyla.pl) 
{
  ps.phy<-subset_taxa(ps6.parasite.leaf, Phylum!= name)
  ps.phy.diff<-ntaxa(ps6.parasite.leaf) - ntaxa(ps.phy)
  ps.phy.ra<-transform_sample_counts(ps.phy, rel.abund)
  wu<-ordinate(ps.phy.ra, method = "MDS", distance = "wunifrac")
  pro<-protest(out.wunifrac.iir$vectors[,1:2], 
               wu$vectors[,1:2], permutations = 999)
  this.result<-c(pro$t0,pro$signif,ps.phy.diff)
  phylum.results<-rbind(phylum.results,this.result)
}

phylum.results$taxa<-phyla.pl
colnames(phylum.results)[1:3]<-c("t0","pval","ASVs.excluded")
phylum.results$deltat0<-phylum.results$t0-pro.pl.iir$t0
write.csv(phylum.results, "Manuscript/leave_one_out.phy.leaf.csv")

pdf(file="leave_one_out.phy.leaf.pdf",
    width=6,height=10,bg = "transparent")
ggplot(data=phylum.results, aes(x=taxa, y=t0)) +
  geom_hline(yintercept=pro.pl.iir$t0, linetype="dashed", color = "red")+
  geom_point(aes(fill=taxa),size=12,color="black") +
  scale_fill_manual(values=phylum.colors) +
  scale_y_continuous(breaks=seq(0.5,0.65,0.01))+
  theme.gg
dev.off()

# Order level #
prot.ord<-droplevels(unique(taxa.pl[which(taxa.pl$Phylum=="Proteobacteria"),4]))
prot.ord<-prot.ord[!is.na(prot.ord)]
actin.ord<-droplevels(unique(taxa.pl[which(taxa.pl$Phylum=="Actinobacteria"),4]))
actin.ord<-actin.ord[!is.na(actin.ord)]
bac.ord<-droplevels(unique(taxa.pl[which(taxa.pl$Phylum=="Bacteroidetes"),4]))
bac.ord<-bac.ord[!is.na(bac.ord)]
order.list<-as.character(unlist(list(prot.ord,actin.ord,bac.ord)))

order.results<-data.frame("t0"=numeric(),"pval"=numeric(),"numer.taxa.exluded"=numeric())
for (name in order.list) 
{
  ps.ord.rem<-subset_taxa(ps6.parasite.leaf, Order == name)
  keep.asv<-setdiff(taxa_names(ps6.parasite.leaf),taxa_names(ps.ord.rem))
  ps.ord<-prune_taxa(keep.asv, ps6.parasite.leaf)
  ps.ord.diff<-ntaxa(ps6.parasite.leaf) - ntaxa(ps.ord)
  ps.ord.ra<-transform_sample_counts(ps.ord, rel.abund)
  wu<-ordinate(ps.ord.ra, method = "MDS", distance = "wunifrac")
  pro<-protest(out.wunifrac.iir$vectors[,1:2], 
               wu$vectors[,1:2], permutations = 999)
  this.result<-c(pro$t0,pro$signif,ps.ord.diff)
  order.results<-rbind(order.results,this.result)
}


order.results$taxa<-order.list
order.results$phylum<-c(rep("Proteobacteria",length(prot.ord)),rep("Actinobacteria",length(actin.ord)),
                        rep("Bacteroidetes",length(bac.ord)))

order.results$taxa <- factor(order.results$taxa,
                             levels = order.list,ordered = TRUE)
colnames(order.results)[1:3]<-c("t0","pval","ASVs.excluded")
order.results$deltat0<-order.results$t0-pro.pl.iir$t0
write.csv(order.results, "Manuscript/leave_one_out.ord.leaf.csv")

pdf(file="leave_one_out.ord.leaf.pdf",
    width=16,height=20,bg = "transparent")
ggplot(data=order.results, aes(x=taxa, y=t0)) +
  geom_hline(yintercept=pro.pr.iir$t0, linetype="dashed", color = "red")+
  geom_point(aes(fill=phylum),size=12,color="black") +
  scale_fill_manual(values=phylum.colors) +
  scale_y_continuous(breaks=seq(0.5,0.75,0.01))+
  theme.gg
dev.off()

# Relative abundance of taxa contributing to Procrustes results #
ps6.ord<-tax_glom(ps6,taxrank="Order")
taxa_names(ps6.ord)<-as.data.frame(tax_table(ps6.ord)[,4])$Order  # update taxa names
ps6.ord = transform_sample_counts(ps6.ord, rel.abund)

ps6.ord.pr<-subset_samples(ps6.ord, IDtype == "P.R")
ps6.ord.iir<-subset_samples(ps6.ord, IDtype == "I.IR")
ps6.ord.pl<-subset_samples(ps6.ord, IDtype == "P.L")

ps6.ord.pr.dat<-as.data.frame(otu_table(ps6.ord.pr))
ps6.ord.iir.dat<-as.data.frame(otu_table(ps6.ord.iir))
ps6.ord.pl.dat<-as.data.frame(otu_table(ps6.ord.pl))

cor.test(ps6.ord.pr.dat$Burkholderiales[-6],ps6.ord.iir.dat$Burkholderiales[-6],
         method="pearson") # r = 0.85, P = 0.000452, without outlier r = 0.54, p = 0.06

cor.test(ps6.ord.pr.dat$Actinomycetales,ps6.ord.iir.dat$Actinomycetales,
         method="pearson") # r = 0.44, P = 0.15

pro.corr.dat<-as.data.frame(cbind(ps6.ord.pr.dat$Burkholderiales,ps6.ord.iir.dat$Burkholderiales,
                    ps6.ord.pr.dat$Actinomycetales,ps6.ord.iir.dat$Actinomycetales))
colnames(pro.corr.dat)<-c("pr.burk","iir.burk","pr.act","iir.act")

pdf(file="Burk.corr.pdf",
    width=8,height=6,bg = "transparent")
ggplot(data=pro.corr.dat,aes(x=iir.burk,y=pr.burk))+
  geom_smooth(method='lm',formula=y~x,se=FALSE,fullrange=TRUE,size=6,colour="black") +
  geom_point(size=10,color="black") +
  scale_x_continuous(breaks=seq(0,0.5,0.02)) +
  scale_y_continuous(breaks=seq(0,0.5,0.1)) +
  theme.gg
dev.off()

pdf(file="act.corr.pdf",
    width=8,height=6,bg = "transparent")
ggplot(data=pro.corr.dat,aes(x=iir.act,y=pr.act))+
  geom_point(size=10,color="black") +
  scale_x_continuous(breaks=seq(0,0.5,0.05)) +
  scale_y_continuous(breaks=seq(0,0.3,0.05)) +
  theme.gg
dev.off()


#### 2.2) Linear mixed effects models #####
## Leaves roots and soil ##
out.wunifrac <- phyloseq::ordinate(ps6.ra, method = "MDS", distance = "wunifrac")
out.uunifrac <- ordinate(ps6.ra, method = "MDS", distance = "unifrac")
out.bray <- phyloseq::ordinate(ps6.ra, method = "MDS", distance = "bray")

# Leaf samples separately
out.wunifrac.l <- phyloseq::ordinate(ps6.leaf.ra, method = "MDS", distance = "wunifrac")
out.uunifrac.l <- ordinate(ps6.leaf.ra, method = "MDS", distance = "unifrac")
out.bray.l <- phyloseq::ordinate(ps6.leaf.ra, method = "MDS", distance = "bray")

# Root samples separately
out.wunifrac.r <- phyloseq::ordinate(ps6.root.ra, method = "MDS", distance = "wunifrac")
out.uunifrac.r <- ordinate(ps6.root.ra, method = "MDS", distance = "unifrac")
out.bray.r <- phyloseq::ordinate(ps6.root.ra, method = "MDS", distance = "bray")

# Leaves and roots together #
out.wunifrac.lr <- phyloseq::ordinate(ps6.leaf.root.ra, method = "MDS", distance = "wunifrac")
out.uunifrac.lr <- ordinate(ps6.leaf.root.ra, method = "MDS", distance = "unifrac")
out.bray.lr <- phyloseq::ordinate(ps6.leaf.root.ra, method = "MDS", distance = "bray")

# Leaf and root microbiome together
# wunifrac #
sum(out.wunifrac.lr$values$Eigenvalues[1:3])/sum(out.wunifrac.lr$values$Eigenvalues) # 0.62
a1.w.lr<-as.data.frame(out.wunifrac.lr$vectors[,1])
a2.w.lr<-as.data.frame(out.wunifrac.lr$vectors[,2])
a3.w.lr<-as.data.frame(out.wunifrac.lr$vectors[,3])
pcoa.axes.wunifrac.lr<-cbind(a1.w.lr,a2.w.lr,a3.w.lr)
colnames(pcoa.axes.wunifrac.lr)<-c("a1.w","a2.w","a3.w")

## uunifrac ##
sum(out.uunifrac.lr$values$Eigenvalues[1:3])/sum(out.uunifrac.lr$values$Eigenvalues) # 0.53
a1.u.lr<-as.data.frame(out.uunifrac.lr$vectors[,1])
a2.u.lr<-as.data.frame(out.uunifrac.lr$vectors[,2])
a3.u.lr<-as.data.frame(out.uunifrac.lr$vectors[,3])
pcoa.axes.uunifrac.lr<-cbind(a1.u.lr,a2.u.lr,a3.u.lr)
colnames(pcoa.axes.uunifrac.lr)<-c("a1.u","a2.u","a3.u")

## Bray-Curtis ##
sum(out.bray.lr$values$Eigenvalues[1:3])/sum(out.bray.lr$values$Eigenvalues) # 0.33
a1.b.lr<-as.data.frame(out.bray.lr$vectors[,1])
a2.b.lr<-as.data.frame(out.bray.lr$vectors[,2])
a3.b.lr<-as.data.frame(out.bray.lr$vectors[,3])
pcoa.axes.bray.lr<-cbind(a1.b.lr,a2.b.lr,a3.b.lr)
colnames(pcoa.axes.bray.lr)<-c("a1.b","a2.b","a3.b")

## Combining all measures 
meta.data.ps6.lr<-as(sample_data(ps6.leaf.root.ra),"data.frame")
ps6.pcoa.lr<-cbind(meta.data.ps6.lr,pcoa.axes.bray.lr,pcoa.axes.uunifrac.lr,pcoa.axes.wunifrac.lr)
ps6.pcoa.lr$log.reads<-log(ps6.pcoa.lr$UsableReads)
ps6.pcoa.lr$sample<-ifelse(ps6.pcoa.lr$type=="IR","R",
                  ifelse(ps6.pcoa.lr$type=="R","R",
                         ifelse(ps6.pcoa.lr$type=="L","L","NA")))


# LMMs testing for effect of species, tissue #
# wunifrac #
m1<-lmer(a1.w ~ log.reads + species*sample + (1|species:site) + (1|sample:site) + (1|species:sample:site) + (1|site), data=ps6.pcoa.lr, REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(a2.w ~ log.reads + species*sample + (1|species:site) + (1|sample:site) + (1|species:sample:site) + (1|site), data=ps6.pcoa.lr, REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a2.w~ species + (1|species:site) + (1|sample:site) + (1|species:sample:site) + (1|site), data=ps6.pcoa.lr, REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(a3.w ~ log.reads + species*sample + (1|species:site) + (1|sample:site) + (1|species:sample:site) + (1|site), data=ps6.pcoa.lr, REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a3.w~ species*sample + (1|species:site) + (1|sample:site) + (1|species:sample:site) + (1|site), data=ps6.pcoa.lr, REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

# Bray Curtis #
m1<-lmer(a1.b ~ log.reads + species*sample + (1|species:site) + (1|sample:site) + (1|species:sample:site) + (1|site), data=ps6.pcoa.lr, REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(a2.b ~ log.reads + species*sample + (1|species:site) + (1|sample:site) + (1|species:sample:site) + (1|site), data=ps6.pcoa.lr, REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a2.b~ species*sample + (1|species:site) + (1|sample:site) + (1|species:sample:site) + (1|site), data=ps6.pcoa.lr, REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(a3.b ~ log.reads + species*sample + (1|species:site) + (1|sample:site) + (1|species:sample:site) + (1|site), data=ps6.pcoa.lr, REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a3.b ~ species*sample + (1|species:site) + (1|sample:site) + (1|species:sample:site) + (1|site), data=ps6.pcoa.lr, REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

# Just Ivy leaf samples - effect of infection
# wunifrac #
m1<-lmer(a1.w ~ log.reads + ID + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="L",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a1.w ~ 1 + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="L",], REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(a2.w ~ log.reads + ID + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="L",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a2.w ~ 1 + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="L",], REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(a3.w ~ log.reads + ID + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="L",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a3.w ~ 1 + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="L",], REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

# Bray-Curtis
m1<-lmer(a1.b ~ log.reads + ID + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="L",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a1.b ~ log.reads + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="L",], REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(a2.b ~ log.reads + ID + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="L",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a2.b ~ 1 + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="L",], REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(a3.b ~ log.reads + ID + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="L",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a3.b ~ 1 + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="L",], REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

# Just Ivy root samples - effect of infection
# wunifrac #
m1<-lmer(a1.w ~ log.reads + ID + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="R",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a1.w ~ 1 + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="R",], REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(a2.w ~ log.reads + ID + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="R",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a2.w ~ 1 + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="R",], REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(a3.w ~ log.reads + ID + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="R",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a3.w ~ 1 + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="R",], REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

# Bray-Curtis
m1<-lmer(a1.b ~ log.reads + ID + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="R",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a1.b ~ 1 + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="R",], REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(a2.b ~ log.reads + ID + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="R",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a2.b ~ 1 + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="R",], REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

m1<-lmer(a3.b ~ log.reads + ID + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="R",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a3.b ~ 1 + (1|ID:site) + (1|site), data=ps6.pcoa.lr[ps6.pcoa.lr$species=="Ha" & ps6.pcoa.lr$sample=="R",], REML=T)
ranova(m1)
plot(m1)
hist(resid(m1))

# Leaf microbiome #
## axis scores for each distance metric ##
## wunifrac ##
sum(out.wunifrac.l$values$Eigenvalues[1:3])/sum(out.wunifrac.l$values$Eigenvalues) # 0.66
a1.w.l<-as.data.frame(out.wunifrac.l$vectors[,1])
a2.w.l<-as.data.frame(out.wunifrac.l$vectors[,2])
a3.w.l<-as.data.frame(out.wunifrac.l$vectors[,3])
pcoa.axes.wunifrac.l<-cbind(a1.w.l,a2.w.l,a3.w.l)
colnames(pcoa.axes.wunifrac.l)<-c("a1.w","a2.w","a3.w")

## uunifrac ##
sum(out.uunifrac.l$values$Eigenvalues[1:3])/sum(out.uunifrac.l$values$Eigenvalues) # 0.44
a1.u.l<-as.data.frame(out.uunifrac.l$vectors[,1])
a2.u.l<-as.data.frame(out.uunifrac.l$vectors[,2])
a3.u.l<-as.data.frame(out.uunifrac.l$vectors[,3])
pcoa.axes.uunifrac.l<-cbind(a1.u.l,a2.u.l,a3.u.l)
colnames(pcoa.axes.uunifrac.l)<-c("a1.u","a2.u","a3.u")

## Bray-Curtis ##
sum(out.bray.l$values$Eigenvalues[1:3])/sum(out.bray.l$values$Eigenvalues) # 0.43
a1.b.l<-as.data.frame(out.bray.l$vectors[,1])
a2.b.l<-as.data.frame(out.bray.l$vectors[,2])
a3.b.l<-as.data.frame(out.bray.l$vectors[,3])
pcoa.axes.bray.l<-cbind(a1.b.l,a2.b.l,a3.b.l)
colnames(pcoa.axes.bray.l)<-c("a1.b","a2.b","a3.b")

## Combining all measures ## Dataframe of PCoA scores of endo and rhizo together ##
meta.data.ps6.l<-as(sample_data(ps6.leaf.ra),"data.frame")
ps6.pcoa.l<-cbind(meta.data.ps6.l,pcoa.axes.bray.l,pcoa.axes.uunifrac.l,pcoa.axes.wunifrac.l)
ps6.pcoa.l$log.reads<-log(ps6.pcoa.l$UsableReads)

# Parasite and ivy leaf microbiome - Effect of plant species #
# weighted unifrac axis scores
m1<-lmer(a1.w~ log.reads + species + (1|species:site) + (1|site), data=ps6.pcoa.l, REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a1.w ~ species + (1|species:site) + (1|site), data=ps6.pcoa.l, REML=T)
ranova(m1)
plot(m1)

m1<-lmer(a2.w~ log.reads + species + (1|species:site) + (1|site), data=ps6.pcoa.l, REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a2.w ~ 1 + (1|species:site) + (1|site), data=ps6.pcoa.l, REML=T)
ranova(m1)
plot(m1)

m1<-lmer(a3.w~ log.reads + species + (1|species:site) + (1|site), data=ps6.pcoa.l, REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a3.w ~ 1 + (1|species:site) + (1|site), data=ps6.pcoa.l, REML=T)
ranova(m1)
plot(m1)

# Bray-Curtis axis scores
m1<-lmer(a1.b~ log.reads + species + (1|species:site) + (1|site), data=ps6.pcoa.l, REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a1.b ~ species + (1|species:site) + (1|site), data=ps6.pcoa.l, REML=T)
ranova(m1)
plot(m1)

m1<-lmer(a2.b~ log.reads + species + (1|species:site) + (1|site), data=ps6.pcoa.l, REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a2.b ~ 1 + (1|species:site) + (1|site), data=ps6.pcoa.l, REML=T)
ranova(m1)
plot(m1)

m1<-lmer(a3.b~ log.reads + species + (1|species:site) + (1|site), data=ps6.pcoa.l, REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a3.b ~ 1 + (1|species:site) + (1|site), data=ps6.pcoa.l, REML=T)
ranova(m1)
plot(m1)

# Ivy leaf microbiome - Effect of infection #
# weighted unifrac axis scores
m1<-lmer(a1.w~ log.reads + ID + (1|ID:site) + (1|site), data=ps6.pcoa.l[ps6.pcoa.l$species=="Ha",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a1.w ~ 1 + (1|ID:site) + (1|site), data=ps6.pcoa.l[ps6.pcoa.l$species=="Ha",], REML=T)
ranova(m1)
plot(m1)

m1<-lmer(a2.w~ log.reads + ID + (1|ID:site) + (1|site), data=ps6.pcoa.l[ps6.pcoa.l$species=="Ha",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a2.w ~ log.reads + (1|ID:site) + (1|site), data=ps6.pcoa.l[ps6.pcoa.l$species=="Ha",], REML=T)
ranova(m1)
plot(m1)

m1<-lmer(a3.w~ log.reads + ID + (1|ID:site) + (1|site), data=ps6.pcoa.l[ps6.pcoa.l$species=="Ha",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a3.w ~ 1 + (1|ID:site) + (1|site), data=ps6.pcoa.l[ps6.pcoa.l$species=="Ha",], REML=T)
ranova(m1)
plot(m1)

# Bray-Curtis axis scores
m1<-lmer(a1.b~ log.reads + ID + (1|ID:site) + (1|site), data=ps6.pcoa.l[ps6.pcoa.l$species=="Ha",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a1.b ~ 1 + (1|ID:site) + (1|site), data=ps6.pcoa.l[ps6.pcoa.l$species=="Ha",], REML=T)
ranova(m1)
plot(m1)

m1<-lmer(a2.b~ log.reads + ID + (1|ID:site) + (1|site), data=ps6.pcoa.l[ps6.pcoa.l$species=="Ha",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a2.b ~ 1 + (1|ID:site) + (1|site), data=ps6.pcoa.l[ps6.pcoa.l$species=="Ha",], REML=T)
ranova(m1)
plot(m1)

m1<-lmer(a3.b~ log.reads + ID + (1|ID:site) + (1|site), data=ps6.pcoa.l[ps6.pcoa.l$species=="Ha",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a3.b ~ 1 + (1|ID:site) + (1|site), data=ps6.pcoa.l[ps6.pcoa.l$species=="Ha",], REML=T)
ranova(m1)
plot(m1)

### Root samples ###
## axis scores for each distance metric ##
## wunifrac ##
# Leaf microbiome #
## axis scores for each distance metric ##
## wunifrac ##
sum(out.wunifrac.r$values$Eigenvalues[1:3])/sum(out.wunifrac.r$values$Eigenvalues) # 0.65
a1.w.r<-as.data.frame(out.wunifrac.r$vectors[,1])
a2.w.r<-as.data.frame(out.wunifrac.r$vectors[,2])
a3.w.r<-as.data.frame(out.wunifrac.r$vectors[,3])
pcoa.axes.wunifrac.r<-cbind(a1.w.r,a2.w.r,a3.w.r)
colnames(pcoa.axes.wunifrac.r)<-c("a1.w","a2.w","a3.w")

## uunifrac ##
sum(out.uunifrac.r$values$Eigenvalues[1:3])/sum(out.uunifrac.r$values$Eigenvalues) # 0.38
a1.u.r<-as.data.frame(out.uunifrac.r$vectors[,1])
a2.u.r<-as.data.frame(out.uunifrac.r$vectors[,2])
a3.u.r<-as.data.frame(out.uunifrac.r$vectors[,3])
pcoa.axes.uunifrac.r<-cbind(a1.u.r,a2.u.r,a3.u.r)
colnames(pcoa.axes.uunifrac.r)<-c("a1.u","a2.u","a3.u")

## Bray-Curtis ##
sum(out.bray.r$values$Eigenvalues[1:3])/sum(out.bray.r$values$Eigenvalues) # 0.35
a1.b.r<-as.data.frame(out.bray.r$vectors[,1])
a2.b.r<-as.data.frame(out.bray.r$vectors[,2])
a3.b.r<-as.data.frame(out.bray.r$vectors[,3])
pcoa.axes.bray.r<-cbind(a1.b.r,a2.b.r,a3.b.r)
colnames(pcoa.axes.bray.r)<-c("a1.b","a2.b","a3.b")

## Combining all measures ## Dataframe of PCoA scores of endo and rhizo together ##
meta.data.ps6.r<-as(sample_data(ps6.root.ra),"data.frame")
ps6.pcoa.r<-cbind(meta.data.ps6.r,pcoa.axes.bray.r,pcoa.axes.uunifrac.r,pcoa.axes.wunifrac.r)
ps6.pcoa.r$log.reads<-log(ps6.pcoa.r$UsableReads)

# Parasite and ivy leaf microbiome - Effect of plant species #
# weighted unifrac axis scores
m1<-lmer(a1.w~ log.reads + species + (1|species:site) + (1|site), data=ps6.pcoa.r, REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a1.w ~ species + (1|species:site) + (1|site), data=ps6.pcoa.r, REML=T)
ranova(m1)
plot(m1)

m1<-lmer(a2.w~ log.reads + species + (1|species:site) + (1|site), data=ps6.pcoa.r, REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a2.w ~ 1 + (1|species:site) + (1|site), data=ps6.pcoa.r, REML=T)
ranova(m1)
plot(m1)

m1<-lmer(a3.w~ log.reads + species + (1|species:site) + (1|site), data=ps6.pcoa.r, REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a3.w ~ 1 + (1|species:site) + (1|site), data=ps6.pcoa.r, REML=T)
ranova(m1)
plot(m1)

# Bray-Curtis axis scores
m1<-lmer(a1.b~ log.reads + species + (1|species:site) + (1|site), data=ps6.pcoa.r, REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a1.b ~ species + (1|species:site) + (1|site), data=ps6.pcoa.r, REML=T)
ranova(m1)
plot(m1)

m1<-lmer(a2.b~ log.reads + species + (1|species:site) + (1|site), data=ps6.pcoa.r, REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a2.b ~ species + (1|species:site) + (1|site), data=ps6.pcoa.r, REML=T)
ranova(m1)
plot(m1)

m1<-lmer(a3.b~ log.reads + species + (1|species:site) + (1|site), data=ps6.pcoa.r, REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a3.b ~ 1 + (1|species:site) + (1|site), data=ps6.pcoa.r, REML=T)
ranova(m1)
plot(m1)

# Ivy root microbiome - Effect of infection #
# weighted unifrac axis scores
m1<-lmer(a1.w~ log.reads + IDtype + (1|IDtype:site) + (1|site), data=ps6.pcoa.r[ps6.pcoa.r$species=="Ha",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a1.w ~ 1 + (1|IDtype:site) + (1|site), data=ps6.pcoa.r[ps6.pcoa.r$species=="Ha",], REML=T)
ranova(m1)
plot(m1)

m1<-lmer(a2.w~ log.reads + IDtype + (1|IDtype:site) + (1|site), data=ps6.pcoa.r[ps6.pcoa.r$species=="Ha",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a2.w ~ 1 + (1|IDtype:site) + (1|site), data=ps6.pcoa.r[ps6.pcoa.r$species=="Ha",], REML=T)
ranova(m1)
plot(m1)

m1<-lmer(a3.w~ log.reads + IDtype + (1|IDtype:site) + (1|site), data=ps6.pcoa.r[ps6.pcoa.r$species=="Ha",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a3.w ~ 1 + (1|IDtype:site) + (1|site), data=ps6.pcoa.r[ps6.pcoa.r$species=="Ha",], REML=T)
ranova(m1)
plot(m1)

# Bray-Curtis axis scores
m1<-lmer(a1.b~ log.reads + IDtype + (1|IDtype:site) + (1|site), data=ps6.pcoa.r[ps6.pcoa.r$species=="Ha",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a1.b ~ 1 + (1|IDtype:site) + (1|site), data=ps6.pcoa.r[ps6.pcoa.r$species=="Ha",], REML=T)
ranova(m1)
plot(m1)

m1<-lmer(a2.b~ log.reads + IDtype + (1|IDtype:site) + (1|site), data=ps6.pcoa.r[ps6.pcoa.r$species=="Ha",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a2.b ~ 1 + (1|IDtype:site) + (1|site), data=ps6.pcoa.r[ps6.pcoa.r$species=="Ha",], REML=T)
ranova(m1)
plot(m1)

m1<-lmer(a3.b~ log.reads + IDtype + (1|IDtype:site) + (1|site), data=ps6.pcoa.r[ps6.pcoa.r$species=="Ha",], REML=T)
car::Anova(m1,type="III",test.statistic="F",contrasts=list(topic=contr.sum, sys=contr.sum))
m1<-lmer(a3.b ~ IDtype + (1|IDtype:site) + (1|site), data=ps6.pcoa.r[ps6.pcoa.r$species=="Ha",], REML=T)
ranova(m1)
plot(m1)

# Beta diverity figures #
# ordination of leaf and root samples together
sample_data(ps6.leaf.root.ra)$sample<-ifelse(sample_data(ps6.leaf.root.ra)$type=="IR","R",
                                                                 ifelse(sample_data(ps6.leaf.root.ra)$type=="R","R",
                                                                        ifelse(sample_data(ps6.leaf.root.ra)$type=="L","L","NA")))

# weighted unifrac
pdf(file="PCoA.wunifrac.pdf",
    width=9,height=7,bg = "transparent")
plot_ordination(ps6.leaf.root.ra, out.wunifrac.lr, type="samples", color="IDtype", shape = "sample") +
  geom_point(aes(fill=IDtype,shape=sample),size=10,color="black") + 
  scale_shape_manual(values=c(21,24)) +
  scale_fill_manual(values=sample.colors.beta) +
  scale_x_continuous(breaks=seq(-0.4,0.5,0.1)) +
  scale_y_continuous(breaks=seq(-0.4,0.4,0.1)) +
  theme.gg.leg
dev.off()

# ordination of samples that exhibit high congruence from Procrustes #
sample.shapes1 <- c(I.IR=3,I.L=21,I.R=3,U.R=3,U.L=3,
                             P.R=3,P.L=24,U.S=3,I.S=3)

pdf(file="PCoA.wunifrac.sig.procrust.pl.il.pdf",
    width=9,height=7,bg = "transparent")
plot_ordination(ps6.ra, out.wunifrac, type="samples", color="site.rep", shape = "IDtype") +
  geom_point(aes(fill=site.rep,shape=IDtype),size=10,color="black") + 
  scale_shape_manual(values=sample.shapes1) +
  scale_fill_manual(values=sample.colors.site.rep) +
  scale_x_continuous(breaks=seq(-0.4,0.5,0.1)) +
  scale_y_continuous(breaks=seq(-0.4,0.4,0.1)) +
  theme.gg.leg
dev.off()

# unweighted unifrac
pdf(file="PCoA.uunifrac.pdf",
    width=9,height=7,bg = "transparent")
plot_ordination(ps6.leaf.root.ra, out.uunifrac.lr, type="samples", color="IDtype", shape = "sample") +
  geom_point(aes(fill=IDtype,shape=sample),size=10,color="black") + 
  scale_shape_manual(values=c(21,24)) +
  scale_fill_manual(values=sample.colors.beta) +
  scale_x_continuous(breaks=seq(-0.4,0.5,0.1)) +
  scale_y_continuous(breaks=seq(-0.4,0.4,0.1)) +
  theme.gg.leg
dev.off()

# Bray-Curtis
pdf(file="PCoA.bray.pdf",
    width=9,height=7,bg = "transparent")
plot_ordination(ps6.leaf.root.ra, out.bray.lr, type="samples", color="IDtype", shape = "sample") +
  geom_point(aes(fill=IDtype,shape=sample),size=10,color="black") + 
  scale_shape_manual(values=c(21,24)) +
  scale_fill_manual(values=sample.colors.beta) +
  scale_x_continuous(breaks=seq(-0.4,0.5,0.1)) +
  scale_y_continuous(breaks=seq(-0.4,0.5,0.1)) +
  theme.gg.leg
dev.off()

# Ordination of just soil samples #
out.wunifrac.soil <- phyloseq::ordinate(ps6.allsoil.ra, method = "MDS", distance = "wunifrac")
pdf(file="PCoA.soil.wunifrac.pdf",
    width=7,height=6,bg = "transparent")
plot_ordination(ps6.allsoil.ra, out.wunifrac.soil, type="samples", color="site", shape = "ID") +
  geom_point(aes(fill=site,shape=ID),size=10,color="black") + 
  scale_shape_manual(values=c(21,24)) +
  scale_fill_manual(values=c("#8dd3c7","#ffffb3","#bebada","#fb8072")) +
  scale_x_continuous(breaks=seq(-0.12,0.1,0.02)) +
  scale_y_continuous(breaks=seq(-0.1,0.1,0.02)) +
  theme.gg
dev.off()
dist.wu.soil <- phyloseq::distance(ps6.allsoil.ra, method="wunifrac")
md.soil<-data.frame(sample_data(ps6.allsoil.ra))
adonis(dist.wu.soil ~ ID+site, data = md.soil)

# scree plots #
wu.full<-as.data.frame(out.wunifrac.lr$values[1:20,])
wu.full$pc<-as.factor(rownames(wu.full))
uu.full<-as.data.frame(out.uunifrac.lr$values[1:20,])
uu.full$pc<-as.factor(rownames(uu.full))
b.full<-as.data.frame(out.bray.lr$values[1:20,])
b.full$pc<-as.factor(rownames(b.full))

pdf(file="scree.wu.pdf",
    width=9,height=6,bg = "transparent")
ggplot(wu.full,aes(x=reorder(pc, -Relative_eig),y=Relative_eig)) +
  geom_bar(stat="identity", fill="#2171b5", color="black") + theme.gg
dev.off()

pdf(file="scree.uu.pdf",
    width=9,height=6,bg = "transparent")
ggplot(uu.full,aes(x=reorder(pc, -Relative_eig),y=Relative_eig)) +
  geom_bar(stat="identity", fill="#2171b5", color="black") + theme.gg
dev.off()

pdf(file="scree.bray.pdf",
    width=9,height=6,bg = "transparent")
ggplot(b.full,aes(x=reorder(pc, -Relative_eig),y=Relative_eig)) +
  geom_bar(stat="identity", fill="#2171b5", color="black") + theme.gg
dev.off()


# Relative abudnance of phyla across sample types - stacked bar charts #
ps6.phy<-tax_glom(ps6,taxrank="Phylum")
ps6.cla<-tax_glom(ps6,taxrank="Class")
taxa_names(ps6.phy)<-as.data.frame(tax_table(ps6.phy)[,2])$Phylum
taxa_names(ps6.cla)<-as.data.frame(tax_table(ps6.cla)[,3])$Class
ps6.phy.ra = transform_sample_counts(ps6.phy, rel.abund)
ps6.cla.ra = transform_sample_counts(ps6.cla, rel.abund)
ps6.phy.ra <-cbind(as(sample_data(ps6.phy.ra),'data.frame'),as.data.frame(otu_table(ps6.phy.ra))) # data frame of phyla abundance and sample data
ps6.cla.ra <-cbind(as(sample_data(ps6.cla.ra),'data.frame'),as.data.frame(otu_table(ps6.cla.ra))) # data frame of phyla abundance and sample data
ps6.phy.ra$IDtypesite<-interaction(ps6.phy.ra$IDtype,ps6.phy.ra$site)
ps6.cla.ra$IDtypesite<-interaction(ps6.cla.ra$IDtype,ps6.cla.ra$site)
# ID X Type
df1.ra <- ps6.phy.ra %>% group_by(IDtype) %>% summarise_all(funs(mean))
df1.ra <- df1.ra[,-c(2:8)] 

df1.ra <- ps6.cla.ra %>% group_by(IDtype) %>% summarise_all(funs(mean))
df1.ra <- df1.ra[,-c(2:8)] 

# ID X Type X site
df1.ra <- ps6.phy.ra %>% group_by(IDtypesite) %>% summarise_all(funs(mean))
df1.ra <- df1.ra[,-c(2:8)] 

df1.ra <- ps6.cla.ra %>% group_by(IDtypesite) %>% summarise_all(funs(mean))
df1.ra <- df1.ra[,-c(2:8)]

mdf1 <- melt(df1.ra, id.vars = "IDtype")
mdf1$variable<-as.character(mdf1$variable)
mdf1$Taxon<-ifelse(mdf1$value>0.01,mdf1$variable,"Low Abundance")

mdf1$Taxon <- factor(mdf1$Taxon,
                            levels = c("Low Abundance","Ignavibacteriae","Latescibacteria",
                                       "Hydrogenedentes","Nitrospirae","Chlamydiae","Cyanobacteria/Chloroplast",
                                       "candidate_division_WPS-1","Candidatus_Saccharibacteria","Gemmatimonadetes",
                                       "candidate_division_WPS-2","Armatimonadetes","Spirochaetes",
                                       "Planctomycetes","Chloroflexi","Firmicutes","Acidobacteria",
                                       "Verrucomicrobia","Bacteroidetes","Actinobacteria","Proteobacteria"),ordered = TRUE)

mdf1$IDtype <- factor(mdf1$IDtype,
                      levels = c("I.L","U.L","I.IR","I.R","U.R","P.L","P.R","I.S","U.S"),
                      ordered = T)

mdf1$IDtypesite <- factor(mdf1$IDtypesite,
                      levels = c("I.L.1","I.L.2","I.L.3","I.L.4",
                                 "U.L.1","U.L.2","U.L.3","U.L.4",
                                 "I.IR.1","I.IR.2","I.IR.3","I.IR.4",
                                 "I.R.1","I.R.2","I.R.3","I.R.4",
                                 "U.R.1","U.R.2","U.R.3","U.R.4",
                                 "P.L.1","P.L.2","P.L.3","P.L.4",
                                 "P.R.1","P.R.2","P.R.3","P.R.4",
                                 "I.S.1","I.S.2","I.S.3","I.S.4",
                                 "U.S.1","U.S.2","U.S.3","U.S.4"),
                      ordered = T)

pdf(file="barchart.phyla.ra.pdf",
    width=8,height=6,bg = "transparent")
ggplot(mdf1, aes(x= IDtype, y= value, fill = Taxon)) + 
  geom_bar(stat = "identity",col="black") + theme.gg.leg +
  scale_fill_manual(values=phylum.colors)
dev.off()

# WUNIFRAC distance between leaves and roots for paraistes and hosts #
ps5.p<-subset_samples(ps5, species == "Oh")
ps5.h<-subset_samples(ps5, IDtype == "I.IR" | IDtype == "I.L")
ps6.p<-prune_taxa(taxa_sums(ps5.p)>0,ps5.p); ntaxa(ps6.p) 
ps6.h<-prune_taxa(taxa_sums(ps5.h)>0,ps5.h); ntaxa(ps6.h) 
ps6.h.ra = transform_sample_counts(ps6.h, rel.abund)
ps6.p.ra = transform_sample_counts(ps6.p, rel.abund)

dist.wu.p <- phyloseq::distance(ps6.p.ra, method="wunifrac")
dist.wu.h <- phyloseq::distance(ps6.h.ra, method="wunifrac")

dist.wu.p<-as.matrix(dist.wu.p)
dist.wu.h<-as.matrix(dist.wu.h)
# extract average dissimilarity between leaves and roots 
dwp1<-dist.wu.p[-c(1:3,7:9,13:15,19:21),c(1:3,7:9,13:15,19:21)]
mean(diag(dwp1)) # average dissimilarity between paired parasite root and leaf = 0.37

dwh1<-dist.wu.h[-c(1:3,7:9,13:15,19:21),c(1:3,7:9,13:15,19:21)]
mean(diag(dwh1)) # average dissimilarity between paired parasite root and leaf = 0.47

# are they sig. different?
t.test(diag(dwp1),diag(dwh1),paired=T,alternative="less") # t = -2.52 df = 11 p = 0.01

# Figure #
mp<-mean(diag(dwp1))
mh<-mean(diag(dwh1))
ep<-st.err(diag(dwp1))
eh<-st.err(diag(dwh1))
m<-c(mp,mh)
e<-c(ep,eh)
agg.dat<-data.frame(m,e)
agg.dat$host<-c("p","h")

pdf(file="mean.wunifrac.leafroot.pdf",
    width=3,height=3,bg = "transparent")
ggplot(data=agg.dat, aes(x=host, y=m)) +
  geom_errorbar(aes(ymin=m-e, ymax=m+e), width=.5) +
  geom_point(aes(fill=host),size=12,color="black") +
  scale_fill_manual(values=c("#004529","#7a0177")) +
  scale_y_continuous(breaks=seq(0.3,0.5,0.05))+
  theme.gg
dev.off()

#### 2.3) Are there unique taxa to parasites or infection status? ####
alltaxa<-taxa_names(ps6)
IStaxa<-taxa_names(ps6.soil)
UStaxa<-taxa_names(ps6.usoil)
PRtaxa<-taxa_names(ps6.parasite.root)
PLtaxa<-taxa_names(ps6.parasite.leaf)
IIRtaxa<-taxa_names(ps6.host.root)
IURtaxa<-taxa_names(ps6.host.uroot)
ILtaxa<-taxa_names(ps6.host.leaf)
URtaxa<-taxa_names(ps6.uninf.root)
ULtaxa<-taxa_names(ps6.uninf.leaf)

venn(list(PRtaxa, IIRtaxa, IURtaxa, URtaxa))
venn(list(PLtaxa, ILtaxa, ULtaxa))

pdf(file="leaf_venn.pdf",
    width=9,height=9,bg = "transparent")
draw.triple.venn(area1 = 200+20+51+13, area2 = 13+51+13+12, area3 = 51+13+21+20, n12 = 13+51, n23 = 51+13, n13 = 20+51, 
                 n123 = 51, category = c("PL", "IL", "UL"), 
                 fill = c("#f768a1", "#004529","#238443"),
                 alpha = rep(0.8, 3),scaled=T)
dev.off()

pdf(file="root_venn.pdf",
    width=9,height=9,bg = "transparent")
draw.quad.venn(area1 = 26+18+147+25+787+53+39+53, 
               area2 = 18+22+147+113+787+406+73+53, 
               area3 = 63+94+113+406+787+147+53+25,
               area4 = 94+131+406+73+787+53+53+39,
               n12 = 53+787+147+18, n13 = 25+53+147+787, n14 = 787+53+53+39, 
               n23 = 147+113+787+406, n24 = 787+406+53+73, n34 = 94+406+787+53,
               n123 = 147+787, n124 = 787+53, n134 = 787+53, n234 = 406+787, n1234 = 787,
               category = c("PR", "IIR", "IUR", "UR"), 
               fill = c("#7a0177", "#004529","#238443", "#addd8e"),
               alpha = rep(0.8, 4),scaled=T)
dev.off()

venn(list(PRtaxa, IIRtaxa, PLtaxa, ILtaxa))
pdf(file="leaf.root_venn.pdf",
    width=9,height=9,bg = "transparent")
draw.quad.venn(area1 = 117+818+159+16+21+8+7+2, 
               area2 = 571+818+32+159+7+21+4+7, 
               area3 = 13+28+32+7+159+21+16+8,
               area4 = 28+12+7+4+21+7+8+2,
               n12 = 818+159+21+7, n13 = 159+21+16+8, n14 = 21+7+8+2, 
               n23 = 32+7+159+21, n24 = 21+7+7+4, n34 = 21+7+7+4,
               n123 = 159+21, n124 = 21+7, n134 = 21+8, n234 = 21+7, n1234 = 21,
               category = c("PR", "IIR", "PL", "IL"), 
               fill = c("#7a0177", "#004529","#f768a1", "#004529"),
               alpha = rep(0.8, 4),scaled=T)
dev.off()


pronly<-(setdiff(PRtaxa,c(IIRtaxa,IURtaxa,URtaxa))); length(pronly) # N=26
plonly<-(setdiff(PLtaxa,c(ILtaxa,ULtaxa))); length(plonly) # N=200
plonly_root<-(setdiff(PLtaxa,c(ILtaxa,ULtaxa,IIRtaxa))); length(plonly_root) # N=25
iironly<-(setdiff(IIRtaxa,c(IURtaxa,URtaxa,PRtaxa))); length(iironly) # N=22
ilonly<-(setdiff(ILtaxa,c(ULtaxa,PLtaxa))); length(ilonly) # N=12

# How prevalent are these unique taxa?
# PR only ASVs
# obtain df of ASV table
pr.otu<-as.data.frame(otu_table(ps6.parasite.root))
pr.taxa<-as.data.frame(tax_table(ps6.parasite.root))
pr.taxa$taxon<-rownames(pr.taxa)
# subset to unique taxa
pr.otu.unique<-pr.otu[,colnames(pr.otu) %in% pronly]

#convert to presence/absence for each otu
pr.otu.unique1<- as.data.frame(apply(pr.otu.unique, 2, function(y) as.numeric(y > 0)))
pr.only.sums<-as.data.frame(colSums(pr.otu.unique1))
pr.only.sums$taxon<-rownames(pr.only.sums)
pr.only.full<-merge(pr.taxa,pr.only.sums,by="taxon")
pr.only.full$sample<-rep("parasite root",nrow(pr.only.full))
colnames(pr.only.full)[8]<-"prevalence"

# PL only ASVs
# obtain df of ASV table
pl.otu<-as.data.frame(otu_table(ps6.parasite.leaf))
pl.taxa<-as.data.frame(tax_table(ps6.parasite.leaf))
pl.taxa$taxon<-rownames(pl.taxa)
# subset to unique taxa
pl.otu.unique<-pl.otu[,colnames(pl.otu) %in% plonly]

#convert to presence/absence for each otu
pl.otu.unique1<- as.data.frame(apply(pl.otu.unique, 2, function(y) as.numeric(y > 0)))
pl.only.sums<-as.data.frame(colSums(pl.otu.unique1))
pl.only.sums$taxon<-rownames(pl.only.sums)
pl.only.full<-merge(pl.taxa,pl.only.sums,by="taxon")
pl.only.full$sample<-rep("parasite leaf",nrow(pl.only.full))
colnames(pl.only.full)[8]<-"prevalence"

# IIR only ASVs
# obtain df of ASV table
iir.otu<-as.data.frame(otu_table(ps6.host.root))
iir.taxa<-as.data.frame(tax_table(ps6.host.root))
iir.taxa$taxon<-rownames(iir.taxa)
# subset to unique taxa
iir.otu.unique<-iir.otu[,colnames(iir.otu) %in% iironly]

#convert to presence/absence for each otu
iir.otu.unique1<- as.data.frame(apply(iir.otu.unique, 2, function(y) as.numeric(y > 0)))
iir.only.sums<-as.data.frame(colSums(iir.otu.unique1))
iir.only.sums$taxon<-rownames(iir.only.sums)
iir.only.full<-merge(iir.taxa,iir.only.sums,by="taxon")
iir.only.full$sample<-rep("infected root",nrow(iir.only.full))
colnames(iir.only.full)[8]<-"prevalence"

# IL only ASVs
# obtain df of ASV table
il.otu<-as.data.frame(otu_table(ps6.host.leaf))
il.taxa<-as.data.frame(tax_table(ps6.host.leaf))
il.taxa$taxon<-rownames(il.taxa)
# subset to unique taxa
il.otu.unique<-il.otu[,colnames(il.otu) %in% ilonly]

#convert to presence/absence for each otu
il.otu.unique1<- as.data.frame(apply(il.otu.unique, 2, function(y) as.numeric(y > 0)))
il.only.sums<-as.data.frame(colSums(il.otu.unique1))
il.only.sums$taxon<-rownames(il.only.sums)
il.only.full<-merge(il.taxa,il.only.sums,by="taxon")
il.only.full$sample<-rep("infected leaf",nrow(il.only.full))
colnames(il.only.full)[8]<-"prevalence"

write.csv(rbind(pr.only.full,pl.only.full,iir.only.full,il.only.full),
          "Manuscript/unique.ASV.csv" )

#### 3) Differential abundance testing ####
# 5 comparisons #
# 1) Infected root versus Parasite root # IIR vs PR
# 2) Infected leaf versus parasite leaf # IL vs PL
# 3) Infected root versus root from infected plant # IIR vs IR
# 4) Infected root versus root from uninfected plant # IIR vs UR
# 5) Uninfected root from infected plant versus root from uninfected plant # IUR vs UR
# 6) Infected leaf versus uninfected leaf # IL vs UL

# remove soil #
ps5.root<-subset_samples(ps5, type == "IR" | type == "R")
ps5.leaf<-subset_samples(ps5, type == "L")
ps5.root<-prune_taxa(taxa_sums(ps5.root)>0,ps5.root); ntaxa(ps5.root) 
ps5.leaf<-prune_taxa(taxa_sums(ps5.leaf)>0,ps5.leaf); ntaxa(ps5.leaf) 

# need to create subsets of data for ALDEx2 #
# 1) Infected root versus Parasite root # IIR vs PR
ps6.IIR.PR<-subset_samples(ps5.root, IDtype == "I.IR" | IDtype == "P.R")
ps7.IIR.PR<-prune_taxa(taxa_sums(ps6.IIR.PR)>0,ps6.IIR.PR)
conds<-as.data.frame(sample_data(ps7.IIR.PR))
conds.IDtype.IIR.PR<-as.character(conds$IDtype)
# 2) Infected root versus root from infected plant # IIR vs IR
ps6.IIR.IR<-subset_samples(ps5.root, IDtype == "I.IR" | IDtype == "I.R")
ps7.IIR.IR<-prune_taxa(taxa_sums(ps6.IIR.IR)>0,ps6.IIR.IR)
conds<-as.data.frame(sample_data(ps7.IIR.IR))
conds.IDtype.IIR.IR<-as.character(conds$IDtype)
# 3) Infected root versus root from uninfected plant # IIR vs UR
ps6.IIR.UR<-subset_samples(ps5.root, IDtype == "I.IR" | IDtype == "U.R")
ps7.IIR.UR<-prune_taxa(taxa_sums(ps6.IIR.UR)>0,ps6.IIR.UR)
conds<-as.data.frame(sample_data(ps7.IIR.UR))
conds.IDtype.IIR.UR<-as.character(conds$IDtype)
# 4) Uninfected root from infected plant versus root from uninfected plant # IUR vs UR
ps6.IUR.UR<-subset_samples(ps5.root, IDtype == "I.R" | IDtype == "U.R")
ps7.IUR.UR<-prune_taxa(taxa_sums(ps6.IUR.UR)>0,ps6.IUR.UR)
conds<-as.data.frame(sample_data(ps7.IUR.UR))
conds.IDtype.IUR.UR<-as.character(conds$IDtype)
# 5) Infected leaf versus parasite leaf # IL vs PL
ps6.IL.PL<-subset_samples(ps5.leaf, IDtype == "I.L" | IDtype == "P.L")
ps7.IL.PL<-prune_taxa(taxa_sums(ps6.IL.PL)>0,ps6.IL.PL)
conds<-as.data.frame(sample_data(ps7.IL.PL))
conds.IDtype.IL.PL<-as.character(conds$IDtype)
# 6) Infected leaf versus uninfected leaf # IL vs UL
ps6.IL.UL<-subset_samples(ps5.leaf, IDtype == "I.L" | IDtype == "U.L")
ps7.IL.UL<-prune_taxa(taxa_sums(ps6.IL.UL)>0,ps6.IL.UL)
conds<-as.data.frame(sample_data(ps7.IL.UL))
conds.IDtype.IL.UL<-as.character(conds$IDtype)
# 7) Uninfected root versus Parasite root # IUR vs PR
ps6.IUR.PR<-subset_samples(ps5.root, IDtype == "I.R" | IDtype == "P.R")
ps7.IUR.PR<-prune_taxa(taxa_sums(ps6.IUR.PR)>0,ps6.IUR.PR)
conds<-as.data.frame(sample_data(ps7.IUR.PR))
conds.IDtype.IUR.PR<-as.character(conds$IDtype)
# 7) Uninfected root versus Parasite root # UR vs PR
ps6.UR.PR<-subset_samples(ps5.root, IDtype == "U.R" | IDtype == "P.R")
ps7.UR.PR<-prune_taxa(taxa_sums(ps6.UR.PR)>0,ps6.UR.PR)
conds<-as.data.frame(sample_data(ps7.UR.PR))
conds.IDtype.UR.PR<-as.character(conds$IDtype)

## Taxonomic agglomerations ##
# use raw count data from thresholded phyloseq onjects
# tax.glom by Phylum #
ps5.root.phy<-tax_glom(ps5.root,taxrank="Phylum")
taxa_names(ps5.root.phy)<-as.data.frame(tax_table(ps5.root.phy)[,2])$Phylum  # update taxa names
ps5.leaf.phy<-tax_glom(ps5.leaf,taxrank="Phylum")
taxa_names(ps5.leaf.phy)<-as.data.frame(tax_table(ps5.leaf.phy)[,2])$Phylum  # update taxa names
ps7.IIR.PR.phy<-tax_glom(ps7.IIR.PR,taxrank="Phylum")
taxa_names(ps7.IIR.PR.phy)<-as.data.frame(tax_table(ps7.IIR.PR.phy)[,2])$Phylum
ps7.IIR.IR.phy<-tax_glom(ps7.IIR.IR,taxrank="Phylum")
taxa_names(ps7.IIR.IR.phy)<-as.data.frame(tax_table(ps7.IIR.IR.phy)[,2])$Phylum
ps7.IIR.UR.phy<-tax_glom(ps7.IIR.UR,taxrank="Phylum")
taxa_names(ps7.IIR.UR.phy)<-as.data.frame(tax_table(ps7.IIR.UR.phy)[,2])$Phylum
ps7.IUR.UR.phy<-tax_glom(ps7.IUR.UR,taxrank="Phylum")
taxa_names(ps7.IUR.UR.phy)<-as.data.frame(tax_table(ps7.IUR.UR.phy)[,2])$Phylum
ps7.IL.PL.phy<-tax_glom(ps7.IL.PL,taxrank="Phylum")
taxa_names(ps7.IL.PL.phy)<-as.data.frame(tax_table(ps7.IL.PL.phy)[,2])$Phylum
ps7.IL.UL.phy<-tax_glom(ps7.IL.UL,taxrank="Phylum")
taxa_names(ps7.IL.UL.phy)<-as.data.frame(tax_table(ps7.IL.UL.phy)[,2])$Phylum
ps7.IUR.PR.phy<-tax_glom(ps7.IUR.PR,taxrank="Phylum")
taxa_names(ps7.IUR.PR.phy)<-as.data.frame(tax_table(ps7.IUR.PR.phy)[,2])$Phylum
ps7.UR.PR.phy<-tax_glom(ps7.UR.PR,taxrank="Phylum")
taxa_names(ps7.UR.PR.phy)<-as.data.frame(tax_table(ps7.UR.PR.phy)[,2])$Phylum

# tax.glom by Class #
ps5.root.cla<-tax_glom(ps5.root,taxrank="Class")
taxa_names(ps5.root.cla)<-as.data.frame(tax_table(ps5.root.cla)[,3])$Class  # update taxa names
ps5.leaf.cla<-tax_glom(ps5.leaf,taxrank="Class")
taxa_names(ps5.leaf.cla)<-as.data.frame(tax_table(ps5.leaf.cla)[,3])$Class  # update taxa names
ps7.IIR.PR.cla<-tax_glom(ps7.IIR.PR,taxrank="Class")
taxa_names(ps7.IIR.PR.cla)<-as.data.frame(tax_table(ps7.IIR.PR.cla)[,3])$Class
ps7.IIR.IR.cla<-tax_glom(ps7.IIR.IR,taxrank="Class")
taxa_names(ps7.IIR.IR.cla)<-as.data.frame(tax_table(ps7.IIR.IR.cla)[,3])$Class
ps7.IIR.UR.cla<-tax_glom(ps7.IIR.UR,taxrank="Class")
taxa_names(ps7.IIR.UR.cla)<-as.data.frame(tax_table(ps7.IIR.UR.cla)[,3])$Class
ps7.IUR.UR.cla<-tax_glom(ps7.IUR.UR,taxrank="Class")
taxa_names(ps7.IUR.UR.cla)<-as.data.frame(tax_table(ps7.IUR.UR.cla)[,3])$Class
ps7.IL.PL.cla<-tax_glom(ps7.IL.PL,taxrank="Class")
taxa_names(ps7.IL.PL.cla)<-as.data.frame(tax_table(ps7.IL.PL.cla)[,3])$Class
ps7.IL.UL.cla<-tax_glom(ps7.IL.UL,taxrank="Class")
taxa_names(ps7.IL.UL.cla)<-as.data.frame(tax_table(ps7.IL.UL.cla)[,3])$Class
ps7.IUR.PR.cla<-tax_glom(ps7.IUR.PR,taxrank="Class")
taxa_names(ps7.IUR.PR.cla)<-as.data.frame(tax_table(ps7.IUR.PR.cla)[,3])$Class
ps7.UR.PR.cla<-tax_glom(ps7.UR.PR,taxrank="Class")
taxa_names(ps7.UR.PR.cla)<-as.data.frame(tax_table(ps7.UR.PR.cla)[,3])$Class

# tax.glom by Order #
ps5.root.ord<-tax_glom(ps5.root,taxrank="Order")
taxa_names(ps5.root.ord)<-as.data.frame(tax_table(ps5.root.ord)[,4])$Order  # update taxa names
ps5.leaf.ord<-tax_glom(ps5.leaf,taxrank="Order")
taxa_names(ps5.leaf.ord)<-as.data.frame(tax_table(ps5.leaf.ord)[,4])$Order  # update taxa names
ps7.IIR.PR.ord<-tax_glom(ps7.IIR.PR,taxrank="Order")
taxa_names(ps7.IIR.PR.ord)<-as.data.frame(tax_table(ps7.IIR.PR.ord)[,4])$Order
ps7.IIR.IR.ord<-tax_glom(ps7.IIR.IR,taxrank="Order")
taxa_names(ps7.IIR.IR.ord)<-as.data.frame(tax_table(ps7.IIR.IR.ord)[,4])$Order
ps7.IIR.UR.ord<-tax_glom(ps7.IIR.UR,taxrank="Order")
taxa_names(ps7.IIR.UR.ord)<-as.data.frame(tax_table(ps7.IIR.UR.ord)[,4])$Order
ps7.IUR.UR.ord<-tax_glom(ps7.IUR.UR,taxrank="Order")
taxa_names(ps7.IUR.UR.ord)<-as.data.frame(tax_table(ps7.IUR.UR.ord)[,4])$Order
ps7.IL.PL.ord<-tax_glom(ps7.IL.PL,taxrank="Order")
taxa_names(ps7.IL.PL.ord)<-as.data.frame(tax_table(ps7.IL.PL.ord)[,4])$Order
ps7.IL.UL.ord<-tax_glom(ps7.IL.UL,taxrank="Order")
taxa_names(ps7.IL.UL.ord)<-as.data.frame(tax_table(ps7.IL.UL.ord)[,4])$Order
ps7.IUR.PR.ord<-tax_glom(ps7.IUR.PR,taxrank="Order")
taxa_names(ps7.IUR.PR.ord)<-as.data.frame(tax_table(ps7.IUR.PR.ord)[,4])$Order
ps7.UR.PR.ord<-tax_glom(ps7.UR.PR,taxrank="Order")
taxa_names(ps7.UR.PR.ord)<-as.data.frame(tax_table(ps7.UR.PR.ord)[,4])$Order

# tax.glom by Family #
ps5.root.fam<-tax_glom(ps5.root,taxrank="Family")
taxa_names(ps5.root.fam)<-as.data.frame(tax_table(ps5.root.fam)[,5])$Family  # update taxa names
ps5.leaf.fam<-tax_glom(ps5.leaf,taxrank="Family")
taxa_names(ps5.leaf.fam)<-as.data.frame(tax_table(ps5.leaf.fam)[,5])$Family  # update taxa names
ps7.IIR.PR.fam<-tax_glom(ps7.IIR.PR,taxrank="Family")
taxa_names(ps7.IIR.PR.fam)<-as.data.frame(tax_table(ps7.IIR.PR.fam)[,5])$Family
ps7.IIR.IR.fam<-tax_glom(ps7.IIR.IR,taxrank="Family")
taxa_names(ps7.IIR.IR.fam)<-as.data.frame(tax_table(ps7.IIR.IR.fam)[,5])$Family
ps7.IIR.UR.fam<-tax_glom(ps7.IIR.UR,taxrank="Family")
taxa_names(ps7.IIR.UR.fam)<-as.data.frame(tax_table(ps7.IIR.UR.fam)[,5])$Family
ps7.IUR.UR.fam<-tax_glom(ps7.IUR.UR,taxrank="Family")
taxa_names(ps7.IUR.UR.fam)<-as.data.frame(tax_table(ps7.IUR.UR.fam)[,5])$Family
ps7.IL.PL.fam<-tax_glom(ps7.IL.PL,taxrank="Family")
taxa_names(ps7.IL.PL.fam)<-as.data.frame(tax_table(ps7.IL.PL.fam)[,5])$Family
ps7.IL.UL.fam<-tax_glom(ps7.IL.UL,taxrank="Family")
taxa_names(ps7.IL.UL.fam)<-as.data.frame(tax_table(ps7.IL.UL.fam)[,5])$Family
ps7.IUR.PR.fam<-tax_glom(ps7.IUR.PR,taxrank="Family")
taxa_names(ps7.IUR.PR.fam)<-as.data.frame(tax_table(ps7.IUR.PR.fam)[,5])$Family
ps7.UR.PR.fam<-tax_glom(ps7.UR.PR,taxrank="Family")
taxa_names(ps7.UR.PR.fam)<-as.data.frame(tax_table(ps7.UR.PR.fam)[,5])$Family

# tax.glom by Genus #
ps5.root.gen<-tax_glom(ps5.root,taxrank="Genus")
taxa_names(ps5.root.gen)<-as.data.frame(tax_table(ps5.root.gen)[,6])$Genus  # update taxa names
ps5.leaf.gen<-tax_glom(ps5.leaf,taxrank="Genus")
taxa_names(ps5.leaf.gen)<-as.data.frame(tax_table(ps5.leaf.gen)[,6])$Genus  # update taxa names
ps7.IIR.PR.gen<-tax_glom(ps7.IIR.PR,taxrank="Genus")
taxa_names(ps7.IIR.PR.gen)<-as.data.frame(tax_table(ps7.IIR.PR.gen)[,6])$Genus
ps7.IIR.IR.gen<-tax_glom(ps7.IIR.IR,taxrank="Genus")
taxa_names(ps7.IIR.IR.gen)<-as.data.frame(tax_table(ps7.IIR.IR.gen)[,6])$Genus
ps7.IIR.UR.gen<-tax_glom(ps7.IIR.UR,taxrank="Genus")
taxa_names(ps7.IIR.UR.gen)<-as.data.frame(tax_table(ps7.IIR.UR.gen)[,6])$Genus
ps7.IUR.UR.gen<-tax_glom(ps7.IUR.UR,taxrank="Genus")
taxa_names(ps7.IUR.UR.gen)<-as.data.frame(tax_table(ps7.IUR.UR.gen)[,6])$Genus
ps7.IL.PL.gen<-tax_glom(ps7.IL.PL,taxrank="Genus")
taxa_names(ps7.IL.PL.gen)<-as.data.frame(tax_table(ps7.IL.PL.gen)[,6])$Genus
ps7.IL.UL.gen<-tax_glom(ps7.IL.UL,taxrank="Genus")
taxa_names(ps7.IL.UL.gen)<-as.data.frame(tax_table(ps7.IL.UL.gen)[,6])$Genus
ps7.IUR.PR.gen<-tax_glom(ps7.IUR.PR,taxrank="Genus")
taxa_names(ps7.IUR.PR.gen)<-as.data.frame(tax_table(ps7.IUR.PR.gen)[,6])$Genus
ps7.UR.PR.gen<-tax_glom(ps7.UR.PR,taxrank="Genus")
taxa_names(ps7.UR.PR.gen)<-as.data.frame(tax_table(ps7.UR.PR.gen)[,6])$Genus

# Phylum level differential abundance testing # 
dds.phy.root<-phyloseq_to_deseq2(ps5.root.phy, ~ IDtype) # create full model without interactions
dds.phy.root = estimateSizeFactors(dds.phy.root, geoMeans = apply(counts(dds.phy.root), 1, gm_mean))
dds.phy.lfce.root<-DESeq(dds.phy.root,test="Wald",fitType="parametric") # perform Wald test to obtain lfce #
dds.phy.leaf<-phyloseq_to_deseq2(ps5.leaf.phy, ~ IDtype) # create full model without interactions
dds.phy.leaf = estimateSizeFactors(dds.phy.leaf, geoMeans = apply(counts(dds.phy.leaf), 1, gm_mean))
dds.phy.lfce.leaf<-DESeq(dds.phy.leaf,test="Wald",fitType="parametric") # perform Wald test to obtain lfce #

res.dds.phy.lfce.1 = as(results(dds.phy.lfce.root,contrast=c("IDtype", "P.R", "I.IR")),'data.frame')
res.dds.phy.lfce.1$Phylum<-rownames(res.dds.phy.lfce.1)
res.dds.phy.lfce.2 = as(results(dds.phy.lfce.leaf,contrast=c("IDtype", "P.L", "I.L")),'data.frame')
res.dds.phy.lfce.2$Phylum<-rownames(res.dds.phy.lfce.2)
res.dds.phy.lfce.3 = as(results(dds.phy.lfce.root,contrast=c("IDtype", "I.R", "I.IR")),'data.frame')
res.dds.phy.lfce.3$Phylum<-rownames(res.dds.phy.lfce.3)
res.dds.phy.lfce.4 = as(results(dds.phy.lfce.root,contrast=c("IDtype", "U.R", "I.IR")),'data.frame')
res.dds.phy.lfce.4$Phylum<-rownames(res.dds.phy.lfce.4)
res.dds.phy.lfce.5 = as(results(dds.phy.lfce.root,contrast=c("IDtype", "U.R", "I.R")),'data.frame')
res.dds.phy.lfce.5$Phylum<-rownames(res.dds.phy.lfce.5)
res.dds.phy.lfce.6 = as(results(dds.phy.lfce.leaf,contrast=c("IDtype", "U.L", "I.L")),'data.frame')
res.dds.phy.lfce.6$Phylum<-rownames(res.dds.phy.lfce.6)
res.dds.phy.lfce.7 = as(results(dds.phy.lfce.root,contrast=c("IDtype", "P.R", "I.R")),'data.frame')
res.dds.phy.lfce.7$Phylum<-rownames(res.dds.phy.lfce.7)
res.dds.phy.lfce.8 = as(results(dds.phy.lfce.root,contrast=c("IDtype", "P.R", "U.R")),'data.frame')
res.dds.phy.lfce.8$Phylum<-rownames(res.dds.phy.lfce.8)

res.phy.root<-cbind(as(tax_table(ps5.root.phy)[rownames(dds.phy.root), ], "matrix"),
                    res.dds.phy.lfce.1[,c(2:6)],
                    res.dds.phy.lfce.7[,c(2:6)],
                    res.dds.phy.lfce.8[,c(2:6)],
                    res.dds.phy.lfce.3[,c(2:6)],
                    res.dds.phy.lfce.4[,c(2:6)],
                    res.dds.phy.lfce.5[,c(2:6)])
res.phy.leaf<-cbind(as(tax_table(ps5.leaf.phy)[rownames(dds.phy.leaf), ], "matrix"),
                    res.dds.phy.lfce.2[,c(2:6)],
                    res.dds.phy.lfce.6[,c(2:6)])
res.phy.total<-merge(res.phy.root,res.phy.leaf[,-c(1,3:6)],by="Phylum",all=T)
res.phy.total<-res.phy.total[,c(1:21,37:41,22:36,42:46)]
names(res.phy.total)[7:46]<-c("lfc.1","lfcse.1","stat.1","p.1","padj.1",
                                 "lfc.2","lfcse.2","stat.2","p.1","padj.2",
                                 "lfc.3","lfcse.3","stat.3","p.1","padj.3",
                                 "lfc.4","lfcse.4","stat.4","p.1","padj.4",
                                 "lfc.5","lfcse.5","stat.5","p.1","padj.5",
                                 "lfc.6","lfcse.6","stat.6","p.1","padj.6",
                                "lfc.7","lfcse.7","stat.7","p.1","padj.7",
                              "lfc.8","lfcse.8","stat.8","p.1","padj.8")
res.phy.total.p<-res.phy.total[,c(11,16,21,26,31,36)]  

# ALDEx2 #
# need different subsets of data for each comparison #
aldex2.phy.IIR.PR.phy<-as.data.frame(t(otu_table(ps7.IIR.PR.phy)))
aldex2.phy.IIR.PR <- aldex(aldex2.phy.IIR.PR.phy, conds.IDtype.IIR.PR, mc.samples=128, test="t", effect=TRUE,
                        include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.phy.IIR.PR$taxa<-rownames(aldex2.phy.IIR.PR)
colnames(aldex2.phy.IIR.PR) <- paste(colnames(aldex2.phy.IIR.PR), "IIR.PR", sep = "_")
aldex2.phy.IUR.PR.phy<-as.data.frame(t(otu_table(ps7.IUR.PR.phy)))
aldex2.phy.IUR.PR <- aldex(aldex2.phy.IUR.PR.phy, conds.IDtype.IUR.PR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.phy.IUR.PR$taxa<-rownames(aldex2.phy.IUR.PR)
colnames(aldex2.phy.IUR.PR) <- paste(colnames(aldex2.phy.IUR.PR), "IUR.PR", sep = "_")
aldex2.phy.UR.PR.phy<-as.data.frame(t(otu_table(ps7.UR.PR.phy)))
aldex2.phy.UR.PR <- aldex(aldex2.phy.UR.PR.phy, conds.IDtype.UR.PR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.phy.UR.PR$taxa<-rownames(aldex2.phy.UR.PR)
colnames(aldex2.phy.UR.PR) <- paste(colnames(aldex2.phy.UR.PR), "UR.PR", sep = "_")
aldex2.phy.IL.PL.phy<-as.data.frame(t(otu_table(ps7.IL.PL.phy)))
aldex2.phy.IL.PL <- aldex(aldex2.phy.IL.PL.phy, conds.IDtype.IL.PL, mc.samples=128, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.phy.IL.PL$taxa<-rownames(aldex2.phy.IL.PL)
colnames(aldex2.phy.IL.PL) <- paste(colnames(aldex2.phy.IL.PL), "IL.PL", sep = "_")
aldex2.phy.IIR.IR.phy<-as.data.frame(t(otu_table(ps7.IIR.IR.phy)))
aldex2.phy.IIR.IR <- aldex(aldex2.phy.IIR.IR.phy, conds.IDtype.IIR.IR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.phy.IIR.IR$taxa<-rownames(aldex2.phy.IIR.IR)
colnames(aldex2.phy.IIR.IR) <- paste(colnames(aldex2.phy.IIR.IR), "IIR.IR", sep = "_")
aldex2.phy.IIR.UR.phy<-as.data.frame(t(otu_table(ps7.IIR.UR.phy)))
aldex2.phy.IIR.UR <- aldex(aldex2.phy.IIR.UR.phy, conds.IDtype.IIR.UR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.phy.IIR.UR$taxa<-rownames(aldex2.phy.IIR.UR)
colnames(aldex2.phy.IIR.UR) <- paste(colnames(aldex2.phy.IIR.UR), "IIR.UR", sep = "_")
aldex2.phy.IUR.UR.phy<-as.data.frame(t(otu_table(ps7.IUR.UR.phy)))
aldex2.phy.IUR.UR <- aldex(aldex2.phy.IUR.UR.phy, conds.IDtype.IUR.UR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.phy.IUR.UR$taxa<-rownames(aldex2.phy.IUR.UR)
colnames(aldex2.phy.IUR.UR) <- paste(colnames(aldex2.phy.IUR.UR), "IUR.UR", sep = "_")
aldex2.phy.IL.UL.phy<-as.data.frame(t(otu_table(ps7.IL.UL.phy)))
aldex2.phy.IL.UL <- aldex(aldex2.phy.IL.UL.phy, conds.IDtype.IL.UL, mc.samples=128, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.phy.IL.UL$taxa<-rownames(aldex2.phy.IL.UL)
colnames(aldex2.phy.IL.UL) <- paste(colnames(aldex2.phy.IL.UL), "IL.UL", sep = "_")

ald.phy1<-merge(aldex2.phy.IIR.PR[,c(4:6,9,11:12)],aldex2.phy.IUR.PR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IUR.PR",all=T)
ald.phy2<-merge(ald.phy1,aldex2.phy.UR.PR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_UR.PR",all=T)
ald.phy3<-merge(ald.phy2,aldex2.phy.IL.PL[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IL.PL",all=T)
ald.phy4<-merge(ald.phy3,aldex2.phy.IIR.IR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IIR.IR",all=T)
ald.phy5<-merge(ald.phy4,aldex2.phy.IIR.UR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IIR.UR",all=T)
ald.phy6<-merge(ald.phy5,aldex2.phy.IUR.UR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IUR.UR",all=T)
ald.phy<-merge(ald.phy6,aldex2.phy.IL.UL[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IL.UL",all=T)

# merge DESeq2 and ALDex2 results #
res.phy.AD.total<-merge(res.phy.total,ald.phy,by.x="Phylum",by.y="taxa_IIR.PR",all=T)
res.phy.AD.total<-res.phy.AD.total[,c(1:11,47:51,12:16,52:56,17:21,57:61,22:26,62:66,27:31,67:71,32:36,72:76,37:41,77:81,42:46,82:86)]

# add relative abundance of taxa in roots and leaves #
root.phy.ra = transform_sample_counts(ps5.root.phy, rel.abund)
root.phy.ra <-cbind(as(sample_data(root.phy.ra),'data.frame'),as.data.frame(otu_table(root.phy.ra))) # data frame of phyla abundance and sample data
df1.phy.ra.root <- root.phy.ra[,-c(1:8)] %>% summarise_all(funs(mean))
mdf1.phy.ra.root <- melt(df1.phy.ra.root)
colnames(mdf1.phy.ra.root)<-c("Phylum","root.ra")

leaf.phy.ra = transform_sample_counts(ps5.leaf.phy, rel.abund)
leaf.phy.ra <-cbind(as(sample_data(leaf.phy.ra),'data.frame'),as.data.frame(otu_table(leaf.phy.ra))) # data frame of phyla abundance and sample data
df1.phy.ra.leaf <- leaf.phy.ra[,-c(1:8)] %>% summarise_all(funs(mean))
mdf1.phy.ra.leaf <- melt(df1.phy.ra.leaf)
colnames(mdf1.phy.ra.leaf)<-c("Phylum","leaf.ra")
phy.ra<-merge(mdf1.phy.ra.root,mdf1.phy.ra.leaf,by="Phylum",all=T)
res.phy.AD.total<-merge(res.phy.AD.total,phy.ra,by="Phylum",all=T)
res.phy.AD.total<-res.phy.AD.total[,c(1:6,87:88,7:86)]
write.csv(res.phy.AD.total,"Manuscript/DA/res.phy.AD.total1.csv")

# Class level differential abundance testing # 
dds.cla.root<-phyloseq_to_deseq2(ps5.root.cla, ~ IDtype) # create full model without interactions
dds.cla.root = estimateSizeFactors(dds.cla.root, geoMeans = apply(counts(dds.cla.root), 1, gm_mean))
dds.cla.lfce.root<-DESeq(dds.cla.root,test="Wald",fitType="parametric") # perform Wald test to obtain lfce #
dds.cla.leaf<-phyloseq_to_deseq2(ps5.leaf.cla, ~ IDtype) # create full model without interactions
dds.cla.leaf = estimateSizeFactors(dds.cla.leaf, geoMeans = apply(counts(dds.cla.leaf), 1, gm_mean))
dds.cla.lfce.leaf<-DESeq(dds.cla.leaf,test="Wald",fitType="parametric") # perform Wald test to obtain lfce #

res.dds.cla.lfce.1 = as(results(dds.cla.lfce.root,contrast=c("IDtype", "P.R", "I.IR")),'data.frame')
res.dds.cla.lfce.1$Class<-rownames(res.dds.cla.lfce.1)
res.dds.cla.lfce.2 = as(results(dds.cla.lfce.leaf,contrast=c("IDtype", "P.L", "I.L")),'data.frame')
res.dds.cla.lfce.2$Class<-rownames(res.dds.cla.lfce.2)
res.dds.cla.lfce.3 = as(results(dds.cla.lfce.root,contrast=c("IDtype", "I.R", "I.IR")),'data.frame')
res.dds.cla.lfce.3$Class<-rownames(res.dds.cla.lfce.3)
res.dds.cla.lfce.4 = as(results(dds.cla.lfce.root,contrast=c("IDtype", "U.R", "I.IR")),'data.frame')
res.dds.cla.lfce.4$Class<-rownames(res.dds.cla.lfce.4)
res.dds.cla.lfce.5 = as(results(dds.cla.lfce.root,contrast=c("IDtype", "U.R", "I.R")),'data.frame')
res.dds.cla.lfce.5$Class<-rownames(res.dds.cla.lfce.5)
res.dds.cla.lfce.6 = as(results(dds.cla.lfce.leaf,contrast=c("IDtype", "U.L", "I.L")),'data.frame')
res.dds.cla.lfce.6$Class<-rownames(res.dds.cla.lfce.6)
res.dds.cla.lfce.7 = as(results(dds.cla.lfce.root,contrast=c("IDtype", "P.R", "I.R")),'data.frame')
res.dds.cla.lfce.7$Class<-rownames(res.dds.cla.lfce.1)
res.dds.cla.lfce.8 = as(results(dds.cla.lfce.root,contrast=c("IDtype", "P.R", "U.R")),'data.frame')
res.dds.cla.lfce.8$Class<-rownames(res.dds.cla.lfce.1)

res.cla.root<-cbind(as(tax_table(ps5.root.cla)[rownames(dds.cla.root), ], "matrix"),
                    res.dds.cla.lfce.1[,c(2:6)],
                    res.dds.cla.lfce.7[,c(2:6)],
                    res.dds.cla.lfce.8[,c(2:6)],
                    res.dds.cla.lfce.3[,c(2:6)],
                    res.dds.cla.lfce.4[,c(2:6)],
                    res.dds.cla.lfce.5[,c(2:6)])
res.cla.leaf<-cbind(as(tax_table(ps5.leaf.cla)[rownames(dds.cla.leaf), ], "matrix"),
                    res.dds.cla.lfce.2[,c(2:6)],
                    res.dds.cla.lfce.6[,c(2:6)])
res.cla.total<-merge(res.cla.root,res.cla.leaf[,-c(1:2,4:6)],by="Class",all=T)
res.cla.total<-res.cla.total[,c(1:21,37:41,22:36,42:46)]
names(res.cla.total)[7:46]<-c("lfc.1","lfcse.1","stat.1","p.1","padj.1",
                              "lfc.2","lfcse.2","stat.2","p.1","padj.2",
                              "lfc.3","lfcse.3","stat.3","p.1","padj.3",
                              "lfc.4","lfcse.4","stat.4","p.1","padj.4",
                              "lfc.5","lfcse.5","stat.5","p.1","padj.5",
                              "lfc.6","lfcse.6","stat.6","p.1","padj.6",
                              "lfc.7","lfcse.7","stat.7","p.1","padj.7",
                              "lfc.8","lfcse.8","stat.8","p.1","padj.8")

# ALDEx2 #
# need different subsets of data for each comparison #
aldex2.cla.IIR.PR.cla<-as.data.frame(t(otu_table(ps7.IIR.PR.cla)))
aldex2.cla.IIR.PR <- aldex(aldex2.cla.IIR.PR.cla, conds.IDtype.IIR.PR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.cla.IIR.PR$taxa<-rownames(aldex2.cla.IIR.PR)
colnames(aldex2.cla.IIR.PR) <- paste(colnames(aldex2.cla.IIR.PR), "IIR.PR", sep = "_")
aldex2.cla.IUR.PR.cla<-as.data.frame(t(otu_table(ps7.IUR.PR.cla)))
aldex2.cla.IUR.PR <- aldex(aldex2.cla.IUR.PR.cla, conds.IDtype.IUR.PR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.cla.IUR.PR$taxa<-rownames(aldex2.cla.IUR.PR)
colnames(aldex2.cla.IUR.PR) <- paste(colnames(aldex2.cla.IUR.PR), "IUR.PR", sep = "_")
aldex2.cla.UR.PR.cla<-as.data.frame(t(otu_table(ps7.UR.PR.cla)))
aldex2.cla.UR.PR <- aldex(aldex2.cla.UR.PR.cla, conds.IDtype.UR.PR, mc.samples=128, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.cla.UR.PR$taxa<-rownames(aldex2.cla.UR.PR)
colnames(aldex2.cla.UR.PR) <- paste(colnames(aldex2.cla.UR.PR), "UR.PR", sep = "_")
aldex2.cla.IL.PL.cla<-as.data.frame(t(otu_table(ps7.IL.PL.cla)))
aldex2.cla.IL.PL <- aldex(aldex2.cla.IL.PL.cla, conds.IDtype.IL.PL, mc.samples=128, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.cla.IL.PL$taxa<-rownames(aldex2.cla.IL.PL)
colnames(aldex2.cla.IL.PL) <- paste(colnames(aldex2.cla.IL.PL), "IL.PL", sep = "_")
aldex2.cla.IIR.IR.cla<-as.data.frame(t(otu_table(ps7.IIR.IR.cla)))
aldex2.cla.IIR.IR <- aldex(aldex2.cla.IIR.IR.cla, conds.IDtype.IIR.IR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.cla.IIR.IR$taxa<-rownames(aldex2.cla.IIR.IR)
colnames(aldex2.cla.IIR.IR) <- paste(colnames(aldex2.cla.IIR.IR), "IIR.IR", sep = "_")
aldex2.cla.IIR.UR.cla<-as.data.frame(t(otu_table(ps7.IIR.UR.cla)))
aldex2.cla.IIR.UR <- aldex(aldex2.cla.IIR.UR.cla, conds.IDtype.IIR.UR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.cla.IIR.UR$taxa<-rownames(aldex2.cla.IIR.UR)
colnames(aldex2.cla.IIR.UR) <- paste(colnames(aldex2.cla.IIR.UR), "IIR.UR", sep = "_")
aldex2.cla.IUR.UR.cla<-as.data.frame(t(otu_table(ps7.IUR.UR.cla)))
aldex2.cla.IUR.UR <- aldex(aldex2.cla.IUR.UR.cla, conds.IDtype.IUR.UR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.cla.IUR.UR$taxa<-rownames(aldex2.cla.IUR.UR)
colnames(aldex2.cla.IUR.UR) <- paste(colnames(aldex2.cla.IUR.UR), "IUR.UR", sep = "_")
aldex2.cla.IL.UL.cla<-as.data.frame(t(otu_table(ps7.IL.UL.cla)))
aldex2.cla.IL.UL <- aldex(aldex2.cla.IL.UL.cla, conds.IDtype.IL.UL, mc.samples=128, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.cla.IL.UL$taxa<-rownames(aldex2.cla.IL.UL)
colnames(aldex2.cla.IL.UL) <- paste(colnames(aldex2.cla.IL.UL), "IL.UL", sep = "_")

ald.cla1<-merge(aldex2.cla.IIR.PR[,c(4:6,9,11:12)],aldex2.cla.IUR.PR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IUR.PR",all=T)
ald.cla2<-merge(ald.cla1,aldex2.cla.UR.PR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_UR.PR",all=T)
ald.cla3<-merge(ald.cla2,aldex2.cla.IL.PL[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IL.PL",all=T)
ald.cla4<-merge(ald.cla3,aldex2.cla.IIR.IR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IIR.IR",all=T)
ald.cla5<-merge(ald.cla4,aldex2.cla.IIR.UR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IIR.UR",all=T)
ald.cla6<-merge(ald.cla5,aldex2.cla.IUR.UR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IUR.UR",all=T)
ald.cla<-merge(ald.cla6,aldex2.cla.IL.UL[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IL.UL",all=T)

# merge DESeq2 and ALDex2 results #
res.cla.AD.total<-merge(res.cla.total,ald.cla,by.x="Class",by.y="taxa_IIR.PR",all=T)
res.cla.AD.total<-res.cla.AD.total[,c(1:11,47:51,12:16,52:56,17:21,57:61,22:26,62:66,27:31,67:71,32:36,72:76,37:41,77:81,42:46,82:86)]

# add relative abundance of taxa in roots and leaves #
root.cla.ra = transform_sample_counts(ps5.root.cla, rel.abund)
root.cla.ra <-cbind(as(sample_data(root.cla.ra),'data.frame'),as.data.frame(otu_table(root.cla.ra))) # data frame of clala abundance and sample data
df1.cla.ra.root <- root.cla.ra[,-c(1:8)] %>% summarise_all(funs(mean))
mdf1.cla.ra.root <- melt(df1.cla.ra.root)
colnames(mdf1.cla.ra.root)<-c("Class","root.ra")

leaf.cla.ra = transform_sample_counts(ps5.leaf.cla, rel.abund)
leaf.cla.ra <-cbind(as(sample_data(leaf.cla.ra),'data.frame'),as.data.frame(otu_table(leaf.cla.ra))) # data frame of clala abundance and sample data
df1.cla.ra.leaf <- leaf.cla.ra[,-c(1:8)] %>% summarise_all(funs(mean))
mdf1.cla.ra.leaf <- melt(df1.cla.ra.leaf)
colnames(mdf1.cla.ra.leaf)<-c("Class","leaf.ra")
cla.ra<-merge(mdf1.cla.ra.root,mdf1.cla.ra.leaf,by="Class",all=T)
res.cla.AD.total<-merge(res.cla.AD.total,cla.ra,by="Class",all=T)
res.cla.AD.total<-res.cla.AD.total[,c(1:6,87:88,7:86)]
write.csv(res.cla.AD.total,"Manuscript/DA/res.cla.AD.total1.csv")

# order level differential abundance testing # 
dds.ord.root<-phyloseq_to_deseq2(ps5.root.ord, ~ IDtype) # create full model without interactions
dds.ord.root = estimateSizeFactors(dds.ord.root, geoMeans = apply(counts(dds.ord.root), 1, gm_mean))
dds.ord.lfce.root<-DESeq(dds.ord.root,test="Wald",fitType="parametric") # perform Wald test to obtain lfce #
dds.ord.leaf<-phyloseq_to_deseq2(ps5.leaf.ord, ~ IDtype) # create full model without interactions
dds.ord.leaf = estimateSizeFactors(dds.ord.leaf, geoMeans = apply(counts(dds.ord.leaf), 1, gm_mean))
dds.ord.lfce.leaf<-DESeq(dds.ord.leaf,test="Wald",fitType="parametric") # perform Wald test to obtain lfce #

res.dds.ord.lfce.1 = as(results(dds.ord.lfce.root,contrast=c("IDtype", "P.R", "I.IR")),'data.frame')
res.dds.ord.lfce.1$Order<-rownames(res.dds.ord.lfce.1)
res.dds.ord.lfce.2 = as(results(dds.ord.lfce.leaf,contrast=c("IDtype", "P.L", "I.L")),'data.frame')
res.dds.ord.lfce.2$Order<-rownames(res.dds.ord.lfce.2)
res.dds.ord.lfce.3 = as(results(dds.ord.lfce.root,contrast=c("IDtype", "I.R", "I.IR")),'data.frame')
res.dds.ord.lfce.3$Order<-rownames(res.dds.ord.lfce.3)
res.dds.ord.lfce.4 = as(results(dds.ord.lfce.root,contrast=c("IDtype", "U.R", "I.IR")),'data.frame')
res.dds.ord.lfce.4$Order<-rownames(res.dds.ord.lfce.4)
res.dds.ord.lfce.5 = as(results(dds.ord.lfce.root,contrast=c("IDtype", "U.R", "I.R")),'data.frame')
res.dds.ord.lfce.5$Order<-rownames(res.dds.ord.lfce.5)
res.dds.ord.lfce.6 = as(results(dds.ord.lfce.leaf,contrast=c("IDtype", "U.L", "I.L")),'data.frame')
res.dds.ord.lfce.6$Order<-rownames(res.dds.ord.lfce.6)
res.dds.ord.lfce.7 = as(results(dds.ord.lfce.root,contrast=c("IDtype", "P.R", "I.R")),'data.frame')
res.dds.ord.lfce.7$Order<-rownames(res.dds.ord.lfce.1)
res.dds.ord.lfce.8 = as(results(dds.ord.lfce.root,contrast=c("IDtype", "P.R", "U.R")),'data.frame')
res.dds.ord.lfce.8$Order<-rownames(res.dds.ord.lfce.1)

res.ord.root<-cbind(as(tax_table(ps5.root.ord)[rownames(dds.ord.root), ], "matrix"),
                    res.dds.ord.lfce.1[,c(2:6)],
                    res.dds.ord.lfce.7[,c(2:6)],
                    res.dds.ord.lfce.8[,c(2:6)],
                    res.dds.ord.lfce.3[,c(2:6)],
                    res.dds.ord.lfce.4[,c(2:6)],
                    res.dds.ord.lfce.5[,c(2:6)])
res.ord.leaf<-cbind(as(tax_table(ps5.leaf.ord)[rownames(dds.ord.leaf), ], "matrix"),
                    res.dds.ord.lfce.2[,c(2:6)],
                    res.dds.ord.lfce.6[,c(2:6)])
res.ord.total<-merge(res.ord.root,res.ord.leaf[,-c(1:3,5:6)],by="Order",all=T)

res.ord.total<-res.ord.total[,c(1:21,37:41,22:36,42:46)]
names(res.ord.total)[7:46]<-c("lfc.1","lfcse.1","stat.1","p.1","padj.1",
                              "lfc.2","lfcse.2","stat.2","p.1","padj.2",
                              "lfc.3","lfcse.3","stat.3","p.1","padj.3",
                              "lfc.4","lfcse.4","stat.4","p.1","padj.4",
                              "lfc.5","lfcse.5","stat.5","p.1","padj.5",
                              "lfc.6","lfcse.6","stat.6","p.1","padj.6",
                              "lfc.7","lfcse.7","stat.7","p.1","padj.7",
                              "lfc.8","lfcse.8","stat.8","p.1","padj.8")

# ALDEx2 #
# need different subsets of data for each comparison #
aldex2.ord.IIR.PR.ord<-as.data.frame(t(otu_table(ps7.IIR.PR.ord)))
aldex2.ord.IIR.PR <- aldex(aldex2.ord.IIR.PR.ord, conds.IDtype.IIR.PR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.ord.IIR.PR$taxa<-rownames(aldex2.ord.IIR.PR)
colnames(aldex2.ord.IIR.PR) <- paste(colnames(aldex2.ord.IIR.PR), "IIR.PR", sep = "_")
aldex2.ord.IUR.PR.ord<-as.data.frame(t(otu_table(ps7.IUR.PR.ord)))
aldex2.ord.IUR.PR <- aldex(aldex2.ord.IUR.PR.ord, conds.IDtype.IUR.PR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.ord.IUR.PR$taxa<-rownames(aldex2.ord.IUR.PR)
colnames(aldex2.ord.IUR.PR) <- paste(colnames(aldex2.ord.IUR.PR), "IUR.PR", sep = "_")
aldex2.ord.UR.PR.ord<-as.data.frame(t(otu_table(ps7.UR.PR.ord)))
aldex2.ord.UR.PR <- aldex(aldex2.ord.UR.PR.ord, conds.IDtype.UR.PR, mc.samples=128, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.ord.UR.PR$taxa<-rownames(aldex2.ord.UR.PR)
colnames(aldex2.ord.UR.PR) <- paste(colnames(aldex2.ord.UR.PR), "UR.PR", sep = "_")
aldex2.ord.IL.PL.ord<-as.data.frame(t(otu_table(ps7.IL.PL.ord)))
aldex2.ord.IL.PL <- aldex(aldex2.ord.IL.PL.ord, conds.IDtype.IL.PL, mc.samples=128, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.ord.IL.PL$taxa<-rownames(aldex2.ord.IL.PL)
colnames(aldex2.ord.IL.PL) <- paste(colnames(aldex2.ord.IL.PL), "IL.PL", sep = "_")
aldex2.ord.IIR.IR.ord<-as.data.frame(t(otu_table(ps7.IIR.IR.ord)))
aldex2.ord.IIR.IR <- aldex(aldex2.ord.IIR.IR.ord, conds.IDtype.IIR.IR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.ord.IIR.IR$taxa<-rownames(aldex2.ord.IIR.IR)
colnames(aldex2.ord.IIR.IR) <- paste(colnames(aldex2.ord.IIR.IR), "IIR.IR", sep = "_")
aldex2.ord.IIR.UR.ord<-as.data.frame(t(otu_table(ps7.IIR.UR.ord)))
aldex2.ord.IIR.UR <- aldex(aldex2.ord.IIR.UR.ord, conds.IDtype.IIR.UR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.ord.IIR.UR$taxa<-rownames(aldex2.ord.IIR.UR)
colnames(aldex2.ord.IIR.UR) <- paste(colnames(aldex2.ord.IIR.UR), "IIR.UR", sep = "_")
aldex2.ord.IUR.UR.ord<-as.data.frame(t(otu_table(ps7.IUR.UR.ord)))
aldex2.ord.IUR.UR <- aldex(aldex2.ord.IUR.UR.ord, conds.IDtype.IUR.UR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.ord.IUR.UR$taxa<-rownames(aldex2.ord.IUR.UR)
colnames(aldex2.ord.IUR.UR) <- paste(colnames(aldex2.ord.IUR.UR), "IUR.UR", sep = "_")
aldex2.ord.IL.UL.ord<-as.data.frame(t(otu_table(ps7.IL.UL.ord)))
aldex2.ord.IL.UL <- aldex(aldex2.ord.IL.UL.ord, conds.IDtype.IL.UL, mc.samples=128, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.ord.IL.UL$taxa<-rownames(aldex2.ord.IL.UL)
colnames(aldex2.ord.IL.UL) <- paste(colnames(aldex2.ord.IL.UL), "IL.UL", sep = "_")

ald.ord1<-merge(aldex2.ord.IIR.PR[,c(4:6,9,11:12)],aldex2.ord.IUR.PR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IUR.PR",all=T)
ald.ord2<-merge(ald.ord1,aldex2.ord.UR.PR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_UR.PR",all=T)
ald.ord3<-merge(ald.ord2,aldex2.ord.IL.PL[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IL.PL",all=T)
ald.ord4<-merge(ald.ord3,aldex2.ord.IIR.IR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IIR.IR",all=T)
ald.ord5<-merge(ald.ord4,aldex2.ord.IIR.UR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IIR.UR",all=T)
ald.ord6<-merge(ald.ord5,aldex2.ord.IUR.UR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IUR.UR",all=T)
ald.ord<-merge(ald.ord6,aldex2.ord.IL.UL[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IL.UL",all=T)

# merge DESeq2 and ALDex2 results #
res.ord.AD.total<-merge(res.ord.total,ald.ord,by.x="Order",by.y="taxa_IIR.PR",all=T)
res.ord.AD.total<-res.ord.AD.total[,c(1:11,47:51,12:16,52:56,17:21,57:61,22:26,62:66,27:31,67:71,32:36,72:76,37:41,77:81,42:46,82:86)]

# add relative abundance of taxa in roots and leaves #
root.ord.ra = transform_sample_counts(ps5.root.ord, rel.abund)
root.ord.ra <-cbind(as(sample_data(root.ord.ra),'data.frame'),as.data.frame(otu_table(root.ord.ra))) # data frame of ordla abundance and sample data
df1.ord.ra.root <- root.ord.ra[,-c(1:8)] %>% summarise_all(funs(mean))
mdf1.ord.ra.root <- melt(df1.ord.ra.root)
colnames(mdf1.ord.ra.root)<-c("Order","root.ra")

leaf.ord.ra = transform_sample_counts(ps5.leaf.ord, rel.abund)
leaf.ord.ra <-cbind(as(sample_data(leaf.ord.ra),'data.frame'),as.data.frame(otu_table(leaf.ord.ra))) # data frame of ordla abundance and sample data
df1.ord.ra.leaf <- leaf.ord.ra[,-c(1:8)] %>% summarise_all(funs(mean))
mdf1.ord.ra.leaf <- melt(df1.ord.ra.leaf)
colnames(mdf1.ord.ra.leaf)<-c("Order","leaf.ra")
ord.ra<-merge(mdf1.ord.ra.root,mdf1.ord.ra.leaf,by="Order",all=T)
res.ord.AD.total<-merge(res.ord.AD.total,ord.ra,by="Order",all=T)
res.ord.AD.total<-res.ord.AD.total[,c(1:6,87:88,7:86)]
write.csv(res.ord.AD.total,"Manuscript/DA/res.ord.AD.total1.csv")

# family level differential abundance testing # 
dds.fam.root<-phyloseq_to_deseq2(ps5.root.fam, ~ IDtype) # create full model without interactions
dds.fam.root = estimateSizeFactors(dds.fam.root, geoMeans = apply(counts(dds.fam.root), 1, gm_mean))
dds.fam.lfce.root<-DESeq(dds.fam.root,test="Wald",fitType="parametric") # perform Wald test to obtain lfce #
dds.fam.leaf<-phyloseq_to_deseq2(ps5.leaf.fam, ~ IDtype) # create full model without interactions
dds.fam.leaf = estimateSizeFactors(dds.fam.leaf, geoMeans = apply(counts(dds.fam.leaf), 1, gm_mean))
dds.fam.lfce.leaf<-DESeq(dds.fam.leaf,test="Wald",fitType="parametric") # perform Wald test to obtain lfce #

res.dds.fam.lfce.1 = as(results(dds.fam.lfce.root,contrast=c("IDtype", "P.R", "I.IR")),'data.frame')
res.dds.fam.lfce.1$Family<-rownames(res.dds.fam.lfce.1)
res.dds.fam.lfce.2 = as(results(dds.fam.lfce.leaf,contrast=c("IDtype", "P.L", "I.L")),'data.frame')
res.dds.fam.lfce.2$Family<-rownames(res.dds.fam.lfce.2)
res.dds.fam.lfce.3 = as(results(dds.fam.lfce.root,contrast=c("IDtype", "I.R", "I.IR")),'data.frame')
res.dds.fam.lfce.3$Family<-rownames(res.dds.fam.lfce.3)
res.dds.fam.lfce.4 = as(results(dds.fam.lfce.root,contrast=c("IDtype", "U.R", "I.IR")),'data.frame')
res.dds.fam.lfce.4$Family<-rownames(res.dds.fam.lfce.4)
res.dds.fam.lfce.5 = as(results(dds.fam.lfce.root,contrast=c("IDtype", "U.R", "I.R")),'data.frame')
res.dds.fam.lfce.5$Family<-rownames(res.dds.fam.lfce.5)
res.dds.fam.lfce.6 = as(results(dds.fam.lfce.leaf,contrast=c("IDtype", "U.L", "I.L")),'data.frame')
res.dds.fam.lfce.6$Family<-rownames(res.dds.fam.lfce.6)
res.dds.fam.lfce.7 = as(results(dds.fam.lfce.root,contrast=c("IDtype", "P.R", "I.R")),'data.frame')
res.dds.fam.lfce.7$Family<-rownames(res.dds.fam.lfce.1)
res.dds.fam.lfce.8 = as(results(dds.fam.lfce.root,contrast=c("IDtype", "P.R", "U.R")),'data.frame')
res.dds.fam.lfce.8$Family<-rownames(res.dds.fam.lfce.1)

res.fam.root<-cbind(as(tax_table(ps5.root.fam)[rownames(dds.fam.root), ], "matrix"),
                    res.dds.fam.lfce.1[,c(2:6)],
                    res.dds.fam.lfce.7[,c(2:6)],
                    res.dds.fam.lfce.8[,c(2:6)],
                    res.dds.fam.lfce.3[,c(2:6)],
                    res.dds.fam.lfce.4[,c(2:6)],
                    res.dds.fam.lfce.5[,c(2:6)])
res.fam.leaf<-cbind(as(tax_table(ps5.leaf.fam)[rownames(dds.fam.leaf), ], "matrix"),
                    res.dds.fam.lfce.2[,c(2:6)],
                    res.dds.fam.lfce.6[,c(2:6)])
res.fam.total<-merge(res.fam.root,res.fam.leaf[,-c(1:4,6)],by="Family",all=T)

res.fam.total<-res.fam.total[,c(1:21,37:41,22:36,42:46)]
names(res.fam.total)[7:46]<-c("lfc.1","lfcse.1","stat.1","p.1","padj.1",
                              "lfc.2","lfcse.2","stat.2","p.1","padj.2",
                              "lfc.3","lfcse.3","stat.3","p.1","padj.3",
                              "lfc.4","lfcse.4","stat.4","p.1","padj.4",
                              "lfc.5","lfcse.5","stat.5","p.1","padj.5",
                              "lfc.6","lfcse.6","stat.6","p.1","padj.6",
                              "lfc.7","lfcse.7","stat.7","p.1","padj.7",
                              "lfc.8","lfcse.8","stat.8","p.1","padj.8")

# ALDEx2 #
# need different subsets of data for each comparison #
aldex2.fam.IIR.PR.fam<-as.data.frame(t(otu_table(ps7.IIR.PR.fam)))
aldex2.fam.IIR.PR <- aldex(aldex2.fam.IIR.PR.fam, conds.IDtype.IIR.PR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.fam.IIR.PR$taxa<-rownames(aldex2.fam.IIR.PR)
colnames(aldex2.fam.IIR.PR) <- paste(colnames(aldex2.fam.IIR.PR), "IIR.PR", sep = "_")
aldex2.fam.IUR.PR.fam<-as.data.frame(t(otu_table(ps7.IUR.PR.fam)))
aldex2.fam.IUR.PR <- aldex(aldex2.fam.IUR.PR.fam, conds.IDtype.IUR.PR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.fam.IUR.PR$taxa<-rownames(aldex2.fam.IUR.PR)
colnames(aldex2.fam.IUR.PR) <- paste(colnames(aldex2.fam.IUR.PR), "IUR.PR", sep = "_")
aldex2.fam.UR.PR.fam<-as.data.frame(t(otu_table(ps7.UR.PR.fam)))
aldex2.fam.UR.PR <- aldex(aldex2.fam.UR.PR.fam, conds.IDtype.UR.PR, mc.samples=128, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.fam.UR.PR$taxa<-rownames(aldex2.fam.UR.PR)
colnames(aldex2.fam.UR.PR) <- paste(colnames(aldex2.fam.UR.PR), "UR.PR", sep = "_")
aldex2.fam.IL.PL.fam<-as.data.frame(t(otu_table(ps7.IL.PL.fam)))
aldex2.fam.IL.PL <- aldex(aldex2.fam.IL.PL.fam, conds.IDtype.IL.PL, mc.samples=128, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.fam.IL.PL$taxa<-rownames(aldex2.fam.IL.PL)
colnames(aldex2.fam.IL.PL) <- paste(colnames(aldex2.fam.IL.PL), "IL.PL", sep = "_")
aldex2.fam.IIR.IR.fam<-as.data.frame(t(otu_table(ps7.IIR.IR.fam)))
aldex2.fam.IIR.IR <- aldex(aldex2.fam.IIR.IR.fam, conds.IDtype.IIR.IR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.fam.IIR.IR$taxa<-rownames(aldex2.fam.IIR.IR)
colnames(aldex2.fam.IIR.IR) <- paste(colnames(aldex2.fam.IIR.IR), "IIR.IR", sep = "_")
aldex2.fam.IIR.UR.fam<-as.data.frame(t(otu_table(ps7.IIR.UR.fam)))
aldex2.fam.IIR.UR <- aldex(aldex2.fam.IIR.UR.fam, conds.IDtype.IIR.UR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.fam.IIR.UR$taxa<-rownames(aldex2.fam.IIR.UR)
colnames(aldex2.fam.IIR.UR) <- paste(colnames(aldex2.fam.IIR.UR), "IIR.UR", sep = "_")
aldex2.fam.IUR.UR.fam<-as.data.frame(t(otu_table(ps7.IUR.UR.fam)))
aldex2.fam.IUR.UR <- aldex(aldex2.fam.IUR.UR.fam, conds.IDtype.IUR.UR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.fam.IUR.UR$taxa<-rownames(aldex2.fam.IUR.UR)
colnames(aldex2.fam.IUR.UR) <- paste(colnames(aldex2.fam.IUR.UR), "IUR.UR", sep = "_")
aldex2.fam.IL.UL.fam<-as.data.frame(t(otu_table(ps7.IL.UL.fam)))
aldex2.fam.IL.UL <- aldex(aldex2.fam.IL.UL.fam, conds.IDtype.IL.UL, mc.samples=128, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.fam.IL.UL$taxa<-rownames(aldex2.fam.IL.UL)
colnames(aldex2.fam.IL.UL) <- paste(colnames(aldex2.fam.IL.UL), "IL.UL", sep = "_")

ald.fam1<-merge(aldex2.fam.IIR.PR[,c(4:6,9,11:12)],aldex2.fam.IUR.PR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IUR.PR",all=T)
ald.fam2<-merge(ald.fam1,aldex2.fam.UR.PR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_UR.PR",all=T)
ald.fam3<-merge(ald.fam2,aldex2.fam.IL.PL[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IL.PL",all=T)
ald.fam4<-merge(ald.fam3,aldex2.fam.IIR.IR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IIR.IR",all=T)
ald.fam5<-merge(ald.fam4,aldex2.fam.IIR.UR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IIR.UR",all=T)
ald.fam6<-merge(ald.fam5,aldex2.fam.IUR.UR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IUR.UR",all=T)
ald.fam<-merge(ald.fam6,aldex2.fam.IL.UL[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IL.UL",all=T)

# merge DESeq2 and ALDex2 results #
res.fam.AD.total<-merge(res.fam.total,ald.fam,by.x="Family",by.y="taxa_IIR.PR",all=T)
res.fam.AD.total<-res.fam.AD.total[,c(1:11,47:51,12:16,52:56,17:21,57:61,22:26,62:66,27:31,67:71,32:36,72:76,37:41,77:81,42:46,82:86)]

# add relative abundance of taxa in roots and leaves #
root.fam.ra = transform_sample_counts(ps5.root.fam, rel.abund)
root.fam.ra <-cbind(as(sample_data(root.fam.ra),'data.frame'),as.data.frame(otu_table(root.fam.ra))) # data frame of famla abundance and sample data
df1.fam.ra.root <- root.fam.ra[,-c(1:8)] %>% summarise_all(funs(mean))
mdf1.fam.ra.root <- melt(df1.fam.ra.root)
colnames(mdf1.fam.ra.root)<-c("Family","root.ra")

leaf.fam.ra = transform_sample_counts(ps5.leaf.fam, rel.abund)
leaf.fam.ra <-cbind(as(sample_data(leaf.fam.ra),'data.frame'),as.data.frame(otu_table(leaf.fam.ra))) # data frame of famla abundance and sample data
df1.fam.ra.leaf <- leaf.fam.ra[,-c(1:8)] %>% summarise_all(funs(mean))
mdf1.fam.ra.leaf <- melt(df1.fam.ra.leaf)
colnames(mdf1.fam.ra.leaf)<-c("Family","leaf.ra")
fam.ra<-merge(mdf1.fam.ra.root,mdf1.fam.ra.leaf,by="Family",all=T)
res.fam.AD.total<-merge(res.fam.AD.total,fam.ra,by="Family",all=T)
res.fam.AD.total<-res.fam.AD.total[,c(1:6,87:88,7:86)]
write.csv(res.fam.AD.total,"Manuscript/DA/res.fam.AD.total1.csv")

# genus level differential abundance testing # 
dds.gen.root<-phyloseq_to_deseq2(ps5.root.gen, ~ IDtype) # create full model without interactions
dds.gen.root = estimateSizeFactors(dds.gen.root, geoMeans = apply(counts(dds.gen.root), 1, gm_mean))
dds.gen.lfce.root<-DESeq(dds.gen.root,test="Wald",fitType="parametric") # perform Wald test to obtain lfce #
dds.gen.leaf<-phyloseq_to_deseq2(ps5.leaf.gen, ~ IDtype) # create full model without interactions
dds.gen.leaf = estimateSizeFactors(dds.gen.leaf, geoMeans = apply(counts(dds.gen.leaf), 1, gm_mean))
dds.gen.lfce.leaf<-DESeq(dds.gen.leaf,test="Wald",fitType="parametric") # perform Wald test to obtain lfce #

res.dds.gen.lfce.1 = as(results(dds.gen.lfce.root,contrast=c("IDtype", "P.R", "I.IR")),'data.frame')
res.dds.gen.lfce.1$Genus<-rownames(res.dds.gen.lfce.1)
res.dds.gen.lfce.2 = as(results(dds.gen.lfce.leaf,contrast=c("IDtype", "P.L", "I.L")),'data.frame')
res.dds.gen.lfce.2$Genus<-rownames(res.dds.gen.lfce.2)
res.dds.gen.lfce.3 = as(results(dds.gen.lfce.root,contrast=c("IDtype", "I.R", "I.IR")),'data.frame')
res.dds.gen.lfce.3$Genus<-rownames(res.dds.gen.lfce.3)
res.dds.gen.lfce.4 = as(results(dds.gen.lfce.root,contrast=c("IDtype", "U.R", "I.IR")),'data.frame')
res.dds.gen.lfce.4$Genus<-rownames(res.dds.gen.lfce.4)
res.dds.gen.lfce.5 = as(results(dds.gen.lfce.root,contrast=c("IDtype", "U.R", "I.R")),'data.frame')
res.dds.gen.lfce.5$Genus<-rownames(res.dds.gen.lfce.5)
res.dds.gen.lfce.6 = as(results(dds.gen.lfce.leaf,contrast=c("IDtype", "U.L", "I.L")),'data.frame')
res.dds.gen.lfce.6$Genus<-rownames(res.dds.gen.lfce.6)
res.dds.gen.lfce.7 = as(results(dds.gen.lfce.root,contrast=c("IDtype", "P.R", "I.R")),'data.frame')
res.dds.gen.lfce.7$Genus<-rownames(res.dds.gen.lfce.1)
res.dds.gen.lfce.8 = as(results(dds.gen.lfce.root,contrast=c("IDtype", "P.R", "U.R")),'data.frame')
res.dds.gen.lfce.8$Genus<-rownames(res.dds.gen.lfce.1)

res.gen.root<-cbind(as(tax_table(ps5.root.gen)[rownames(dds.gen.root), ], "matrix"),
                    res.dds.gen.lfce.1[,c(2:6)],
                    res.dds.gen.lfce.7[,c(2:6)],
                    res.dds.gen.lfce.8[,c(2:6)],
                    res.dds.gen.lfce.3[,c(2:6)],
                    res.dds.gen.lfce.4[,c(2:6)],
                    res.dds.gen.lfce.5[,c(2:6)])
res.gen.leaf<-cbind(as(tax_table(ps5.leaf.gen)[rownames(dds.gen.leaf), ], "matrix"),
                    res.dds.gen.lfce.2[,c(2:6)],
                    res.dds.gen.lfce.6[,c(2:6)])
res.gen.total<-merge(res.gen.root,res.gen.leaf[,-c(1:5)],by="Genus",all=T)

res.gen.total<-res.gen.total[,c(1:21,37:41,22:36,42:46)]
names(res.gen.total)[7:46]<-c("lfc.1","lfcse.1","stat.1","p.1","padj.1",
                              "lfc.2","lfcse.2","stat.2","p.1","padj.2",
                              "lfc.3","lfcse.3","stat.3","p.1","padj.3",
                              "lfc.4","lfcse.4","stat.4","p.1","padj.4",
                              "lfc.5","lfcse.5","stat.5","p.1","padj.5",
                              "lfc.6","lfcse.6","stat.6","p.1","padj.6",
                              "lfc.7","lfcse.7","stat.7","p.1","padj.7",
                              "lfc.8","lfcse.8","stat.8","p.1","padj.8")

# ALDEx2 #
# need different subsets of data for each comparison #
aldex2.gen.IIR.PR.gen<-as.data.frame(t(otu_table(ps7.IIR.PR.gen)))
aldex2.gen.IIR.PR <- aldex(aldex2.gen.IIR.PR.gen, conds.IDtype.IIR.PR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.gen.IIR.PR$taxa<-rownames(aldex2.gen.IIR.PR)
colnames(aldex2.gen.IIR.PR) <- paste(colnames(aldex2.gen.IIR.PR), "IIR.PR", sep = "_")
aldex2.gen.IUR.PR.gen<-as.data.frame(t(otu_table(ps7.IUR.PR.gen)))
aldex2.gen.IUR.PR <- aldex(aldex2.gen.IUR.PR.gen, conds.IDtype.IUR.PR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.gen.IUR.PR$taxa<-rownames(aldex2.gen.IUR.PR)
colnames(aldex2.gen.IUR.PR) <- paste(colnames(aldex2.gen.IUR.PR), "IUR.PR", sep = "_")
aldex2.gen.UR.PR.gen<-as.data.frame(t(otu_table(ps7.UR.PR.gen)))
aldex2.gen.UR.PR <- aldex(aldex2.gen.UR.PR.gen, conds.IDtype.UR.PR, mc.samples=128, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.gen.UR.PR$taxa<-rownames(aldex2.gen.UR.PR)
colnames(aldex2.gen.UR.PR) <- paste(colnames(aldex2.gen.UR.PR), "UR.PR", sep = "_")
aldex2.gen.IL.PL.gen<-as.data.frame(t(otu_table(ps7.IL.PL.gen)))
aldex2.gen.IL.PL <- aldex(aldex2.gen.IL.PL.gen, conds.IDtype.IL.PL, mc.samples=128, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.gen.IL.PL$taxa<-rownames(aldex2.gen.IL.PL)
colnames(aldex2.gen.IL.PL) <- paste(colnames(aldex2.gen.IL.PL), "IL.PL", sep = "_")
aldex2.gen.IIR.IR.gen<-as.data.frame(t(otu_table(ps7.IIR.IR.gen)))
aldex2.gen.IIR.IR <- aldex(aldex2.gen.IIR.IR.gen, conds.IDtype.IIR.IR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.gen.IIR.IR$taxa<-rownames(aldex2.gen.IIR.IR)
colnames(aldex2.gen.IIR.IR) <- paste(colnames(aldex2.gen.IIR.IR), "IIR.IR", sep = "_")
aldex2.gen.IIR.UR.gen<-as.data.frame(t(otu_table(ps7.IIR.UR.gen)))
aldex2.gen.IIR.UR <- aldex(aldex2.gen.IIR.UR.gen, conds.IDtype.IIR.UR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.gen.IIR.UR$taxa<-rownames(aldex2.gen.IIR.UR)
colnames(aldex2.gen.IIR.UR) <- paste(colnames(aldex2.gen.IIR.UR), "IIR.UR", sep = "_")
aldex2.gen.IUR.UR.gen<-as.data.frame(t(otu_table(ps7.IUR.UR.gen)))
aldex2.gen.IUR.UR <- aldex(aldex2.gen.IUR.UR.gen, conds.IDtype.IUR.UR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.gen.IUR.UR$taxa<-rownames(aldex2.gen.IUR.UR)
colnames(aldex2.gen.IUR.UR) <- paste(colnames(aldex2.gen.IUR.UR), "IUR.UR", sep = "_")
aldex2.gen.IL.UL.gen<-as.data.frame(t(otu_table(ps7.IL.UL.gen)))
aldex2.gen.IL.UL <- aldex(aldex2.gen.IL.UL.gen, conds.IDtype.IL.UL, mc.samples=128, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.gen.IL.UL$taxa<-rownames(aldex2.gen.IL.UL)
colnames(aldex2.gen.IL.UL) <- paste(colnames(aldex2.gen.IL.UL), "IL.UL", sep = "_")

ald.gen1<-merge(aldex2.gen.IIR.PR[,c(4:6,9,11:12)],aldex2.gen.IUR.PR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IUR.PR",all=T)
ald.gen2<-merge(ald.gen1,aldex2.gen.UR.PR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_UR.PR",all=T)
ald.gen3<-merge(ald.gen2,aldex2.gen.IL.PL[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IL.PL",all=T)
ald.gen4<-merge(ald.gen3,aldex2.gen.IIR.IR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IIR.IR",all=T)
ald.gen5<-merge(ald.gen4,aldex2.gen.IIR.UR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IIR.UR",all=T)
ald.gen6<-merge(ald.gen5,aldex2.gen.IUR.UR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IUR.UR",all=T)
ald.gen<-merge(ald.gen6,aldex2.gen.IL.UL[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IL.UL",all=T)

# merge DESeq2 and ALDex2 results #
res.gen.AD.total<-merge(res.gen.total,ald.gen,by.x="Genus",by.y="taxa_IIR.PR",all=T)
res.gen.AD.total<-res.gen.AD.total[,c(1:11,47:51,12:16,52:56,17:21,57:61,22:26,62:66,27:31,67:71,32:36,72:76,37:41,77:81,42:46,82:86)]

# add relative abundance of taxa in roots and leaves #
root.gen.ra = transform_sample_counts(ps5.root.gen, rel.abund)
root.gen.ra <-cbind(as(sample_data(root.gen.ra),'data.frame'),as.data.frame(otu_table(root.gen.ra))) # data frame of genla abundance and sample data
df1.gen.ra.root <- root.gen.ra[,-c(1:8)] %>% summarise_all(funs(mean))
mdf1.gen.ra.root <- melt(df1.gen.ra.root)
colnames(mdf1.gen.ra.root)<-c("Genus","root.ra")

leaf.gen.ra = transform_sample_counts(ps5.leaf.gen, rel.abund)
leaf.gen.ra <-cbind(as(sample_data(leaf.gen.ra),'data.frame'),as.data.frame(otu_table(leaf.gen.ra))) # data frame of genla abundance and sample data
df1.gen.ra.leaf <- leaf.gen.ra[,-c(1:8)] %>% summarise_all(funs(mean))
mdf1.gen.ra.leaf <- melt(df1.gen.ra.leaf)
colnames(mdf1.gen.ra.leaf)<-c("Genus","leaf.ra")
gen.ra<-merge(mdf1.gen.ra.root,mdf1.gen.ra.leaf,by="Genus",all=T)
res.gen.AD.total<-merge(res.gen.AD.total,gen.ra,by="Genus",all=T)
res.gen.AD.total<-res.gen.AD.total[,c(1:6,87:88,7:86)]
write.csv(res.gen.AD.total,"Manuscript/DA/res.gen.AD.total1.csv")

# asv level differential abundance testing # 
dds.root<-phyloseq_to_deseq2(ps5.root, ~ IDtype) # create full model without interactions
dds.root = estimateSizeFactors(dds.root, geoMeans = apply(counts(dds.root), 1, gm_mean))
dds.lfce.root<-DESeq(dds.root,test="Wald",fitType="parametric") # perform Wald test to obtain lfce #
dds.leaf<-phyloseq_to_deseq2(ps5.leaf, ~ IDtype) # create full model without interactions
dds.leaf = estimateSizeFactors(dds.leaf, geoMeans = apply(counts(dds.leaf), 1, gm_mean))
dds.lfce.leaf<-DESeq(dds.leaf,test="Wald",fitType="parametric") # perform Wald test to obtain lfce #

res.dds.lfce.1 = as(results(dds.lfce.root,contrast=c("IDtype", "P.R", "I.IR")),'data.frame')
res.dds.lfce.1$ASV<-rownames(res.dds.lfce.1)
res.dds.lfce.2 = as(results(dds.lfce.leaf,contrast=c("IDtype", "P.L", "I.L")),'data.frame')
res.dds.lfce.2$ASV<-rownames(res.dds.lfce.2)
res.dds.lfce.3 = as(results(dds.lfce.root,contrast=c("IDtype", "I.R", "I.IR")),'data.frame')
res.dds.lfce.3$ASV<-rownames(res.dds.lfce.3)
res.dds.lfce.4 = as(results(dds.lfce.root,contrast=c("IDtype", "U.R", "I.IR")),'data.frame')
res.dds.lfce.4$ASV<-rownames(res.dds.lfce.4)
res.dds.lfce.5 = as(results(dds.lfce.root,contrast=c("IDtype", "U.R", "I.R")),'data.frame')
res.dds.lfce.5$ASV<-rownames(res.dds.lfce.5)
res.dds.lfce.6 = as(results(dds.lfce.leaf,contrast=c("IDtype", "U.L", "I.L")),'data.frame')
res.dds.lfce.6$ASV<-rownames(res.dds.lfce.6)
res.dds.lfce.7 = as(results(dds.lfce.root,contrast=c("IDtype", "P.R", "I.R")),'data.frame')
res.dds.lfce.7$ASV<-rownames(res.dds.lfce.1)
res.dds.lfce.8 = as(results(dds.lfce.root,contrast=c("IDtype", "P.R", "U.R")),'data.frame')
res.dds.lfce.8$ASV<-rownames(res.dds.lfce.1)

res.root<-cbind(as(tax_table(ps5.root)[rownames(dds.root), ], "matrix"),
                    res.dds.lfce.1[,c(2:6)],
                    res.dds.lfce.7[,c(2:6)],
                    res.dds.lfce.8[,c(2:6)],
                    res.dds.lfce.3[,c(2:6)],
                    res.dds.lfce.4[,c(2:6)],
                    res.dds.lfce.5[,c(2:6)])
res.root$ASV<-rownames(res.root)
res.leaf<-cbind(as(tax_table(ps5.leaf)[rownames(dds.leaf), ], "matrix"),
                    res.dds.lfce.2[,c(2:6)],
                    res.dds.lfce.6[,c(2:6)])
res.leaf$ASV<-rownames(res.leaf)
res.total<-merge(res.root,res.leaf[,-c(1:6)],by="ASV",all=T)

res.total<-res.total[,c(1:22,38:42,23:37,43:47)]

names(res.total)[8:47]<-c("lfc.1","lfcse.1","stat.1","p.1","padj.1",
                              "lfc.2","lfcse.2","stat.2","p.1","padj.2",
                              "lfc.3","lfcse.3","stat.3","p.1","padj.3",
                              "lfc.4","lfcse.4","stat.4","p.1","padj.4",
                              "lfc.5","lfcse.5","stat.5","p.1","padj.5",
                              "lfc.6","lfcse.6","stat.6","p.1","padj.6",
                              "lfc.7","lfcse.7","stat.7","p.1","padj.7",
                              "lfc.8","lfcse.8","stat.8","p.1","padj.8")

# ALDEx2 #
# need different subsets of data for each comparison #
aldex2.IIR.PR<-as.data.frame(t(otu_table(ps7.IIR.PR)))
aldex2.IIR.PR <- aldex(aldex2.IIR.PR, conds.IDtype.IIR.PR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.IIR.PR$taxa<-rownames(aldex2.IIR.PR)
colnames(aldex2.IIR.PR) <- paste(colnames(aldex2.IIR.PR), "IIR.PR", sep = "_")
aldex2.IUR.PR<-as.data.frame(t(otu_table(ps7.IUR.PR)))
aldex2.IUR.PR <- aldex(aldex2.IUR.PR, conds.IDtype.IUR.PR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.IUR.PR$taxa<-rownames(aldex2.IUR.PR)
colnames(aldex2.IUR.PR) <- paste(colnames(aldex2.IUR.PR), "IUR.PR", sep = "_")
aldex2.UR.PR<-as.data.frame(t(otu_table(ps7.UR.PR)))
aldex2.UR.PR <- aldex(aldex2.UR.PR, conds.IDtype.UR.PR, mc.samples=128, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.UR.PR$taxa<-rownames(aldex2.UR.PR)
colnames(aldex2.UR.PR) <- paste(colnames(aldex2.UR.PR), "UR.PR", sep = "_")
aldex2.IL.PL<-as.data.frame(t(otu_table(ps7.IL.PL)))
aldex2.IL.PL <- aldex(aldex2.IL.PL, conds.IDtype.IL.PL, mc.samples=128, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.IL.PL$taxa<-rownames(aldex2.IL.PL)
colnames(aldex2.IL.PL) <- paste(colnames(aldex2.IL.PL), "IL.PL", sep = "_")
aldex2.IIR.IR<-as.data.frame(t(otu_table(ps7.IIR.IR)))
aldex2.IIR.IR <- aldex(aldex2.IIR.IR, conds.IDtype.IIR.IR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.IIR.IR$taxa<-rownames(aldex2.IIR.IR)
colnames(aldex2.IIR.IR) <- paste(colnames(aldex2.IIR.IR), "IIR.IR", sep = "_")
aldex2.IIR.UR<-as.data.frame(t(otu_table(ps7.IIR.UR)))
aldex2.IIR.UR <- aldex(aldex2.IIR.UR, conds.IDtype.IIR.UR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.IIR.UR$taxa<-rownames(aldex2.IIR.UR)
colnames(aldex2.IIR.UR) <- paste(colnames(aldex2.IIR.UR), "IIR.UR", sep = "_")
aldex2.IUR.UR<-as.data.frame(t(otu_table(ps7.IUR.UR)))
aldex2.IUR.UR <- aldex(aldex2.IUR.UR, conds.IDtype.IUR.UR, mc.samples=128, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.IUR.UR$taxa<-rownames(aldex2.IUR.UR)
colnames(aldex2.IUR.UR) <- paste(colnames(aldex2.IUR.UR), "IUR.UR", sep = "_")
aldex2.IL.UL<-as.data.frame(t(otu_table(ps7.IL.UL)))
aldex2.IL.UL <- aldex(aldex2.IL.UL, conds.IDtype.IL.UL, mc.samples=128, test="t", effect=TRUE,
                          include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2.IL.UL$taxa<-rownames(aldex2.IL.UL)
colnames(aldex2.IL.UL) <- paste(colnames(aldex2.IL.UL), "IL.UL", sep = "_")

ald1<-merge(aldex2.IIR.PR[,c(4:6,9,11:12)],aldex2.IUR.PR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IUR.PR",all=T)
ald2<-merge(ald1,aldex2.UR.PR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_UR.PR",all=T)
ald3<-merge(ald2,aldex2.IL.PL[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IL.PL",all=T)
ald4<-merge(ald3,aldex2.IIR.IR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IIR.IR",all=T)
ald5<-merge(ald4,aldex2.IIR.UR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IIR.UR",all=T)
ald6<-merge(ald5,aldex2.IUR.UR[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IUR.UR",all=T)
ald<-merge(ald6,aldex2.IL.UL[,c(4:6,9,11:12)],by.x="taxa_IIR.PR",by.y="taxa_IL.UL",all=T)

# merge DESeq2 and ALDex2 results #
res.AD.total<-merge(res.total,ald,by.x="ASV",by.y="taxa_IIR.PR",all=T)
res.AD.total<-res.AD.total[,c(1:12,48:52,13:17,53:57,18:22,58:62,23:27,63:67,28:32,68:72,33:37,73:77,38:42,78:82,43:47,83:87)]

# add relative abundance of taxa in roots and leaves #
root.ra = transform_sample_counts(ps5.root, rel.abund)
root.ra <-cbind(as(sample_data(root.ra),'data.frame'),as.data.frame(otu_table(root.ra))) # data frame of genla abundance and sample data
df1.ra.root <- root.ra[,-c(1:8)] %>% summarise_all(funs(mean))
mdf1.ra.root <- melt(df1.ra.root)
colnames(mdf1.ra.root)<-c("ASV","root.ra")

leaf.ra = transform_sample_counts(ps5.leaf, rel.abund)
leaf.ra <-cbind(as(sample_data(leaf.ra),'data.frame'),as.data.frame(otu_table(leaf.ra))) # data frame of genla abundance and sample data
df1.ra.leaf <- leaf.ra[,-c(1:8)] %>% summarise_all(funs(mean))
mdf1.ra.leaf <- melt(df1.ra.leaf)
colnames(mdf1.ra.leaf)<-c("ASV","leaf.ra")
asv.ra<-merge(mdf1.ra.root,mdf1.ra.leaf,by="ASV",all=T)
res.AD.total<-merge(res.AD.total,asv.ra,by="ASV",all=T)
res.AD.total<-res.AD.total[,c(1:7,88:89,8:87)]
write.csv(res.AD.total,"Manuscript/DA/res.AD.total1.csv")

#### 3.1) DA figures ####
# Figure demonstrating distribution of LFC across bacterial taxa and specific comparisons
# import DESeq2 and ALDEx2 results
phy.res<-read.csv("Manuscript/DA/res.phy.AD.total.csv")
cla.res<-read.csv("Manuscript/DA/res.cla.AD.total.csv")
ord.res<-read.csv("Manuscript/DA/res.ord.AD.total.csv")
fam.res<-read.csv("Manuscript/DA/res.fam.AD.total.csv")
gen.res<-read.csv("Manuscript/DA/res.gen.AD.total.csv")
asv.res<-read.csv("Manuscript/DA/res.asv.AD.total.csv")

colnames(phy.res)[15:19]<-paste(colnames(phy.res)[15:19], "1", sep = ".")
colnames(phy.res)[25:29]<-paste(colnames(phy.res)[25:29], "2", sep = ".")
colnames(phy.res)[35:39]<-paste(colnames(phy.res)[35:39], "3", sep = ".")
colnames(phy.res)[45:49]<-paste(colnames(phy.res)[45:49], "4", sep = ".")
colnames(phy.res)[55:59]<-paste(colnames(phy.res)[55:59], "5", sep = ".")
colnames(phy.res)[65:69]<-paste(colnames(phy.res)[65:69], "6", sep = ".")
colnames(cla.res)[15:19]<-paste(colnames(cla.res)[15:19], "1", sep = ".")
colnames(cla.res)[25:29]<-paste(colnames(cla.res)[25:29], "2", sep = ".")
colnames(cla.res)[35:39]<-paste(colnames(cla.res)[35:39], "3", sep = ".")
colnames(cla.res)[45:49]<-paste(colnames(cla.res)[45:49], "4", sep = ".")
colnames(cla.res)[55:59]<-paste(colnames(cla.res)[55:59], "5", sep = ".")
colnames(cla.res)[65:69]<-paste(colnames(cla.res)[65:69], "6", sep = ".")
colnames(ord.res)[15:19]<-paste(colnames(ord.res)[15:19], "1", sep = ".")
colnames(ord.res)[25:29]<-paste(colnames(ord.res)[25:29], "2", sep = ".")
colnames(ord.res)[35:39]<-paste(colnames(ord.res)[35:39], "3", sep = ".")
colnames(ord.res)[45:49]<-paste(colnames(ord.res)[45:49], "4", sep = ".")
colnames(ord.res)[55:59]<-paste(colnames(ord.res)[55:59], "5", sep = ".")
colnames(ord.res)[65:69]<-paste(colnames(ord.res)[65:69], "6", sep = ".")
colnames(fam.res)[15:19]<-paste(colnames(fam.res)[15:19], "1", sep = ".")
colnames(fam.res)[25:29]<-paste(colnames(fam.res)[25:29], "2", sep = ".")
colnames(fam.res)[35:39]<-paste(colnames(fam.res)[35:39], "3", sep = ".")
colnames(fam.res)[45:49]<-paste(colnames(fam.res)[45:49], "4", sep = ".")
colnames(fam.res)[55:59]<-paste(colnames(fam.res)[55:59], "5", sep = ".")
colnames(fam.res)[65:69]<-paste(colnames(fam.res)[65:69], "6", sep = ".")
colnames(gen.res)[15:19]<-paste(colnames(gen.res)[15:19], "1", sep = ".")
colnames(gen.res)[25:29]<-paste(colnames(gen.res)[25:29], "2", sep = ".")
colnames(gen.res)[35:39]<-paste(colnames(gen.res)[35:39], "3", sep = ".")
colnames(gen.res)[45:49]<-paste(colnames(gen.res)[45:49], "4", sep = ".")
colnames(gen.res)[55:59]<-paste(colnames(gen.res)[55:59], "5", sep = ".")
colnames(gen.res)[65:69]<-paste(colnames(gen.res)[65:69], "6", sep = ".")
asv.res<-asv.res[,-1]
colnames(asv.res)[15:19]<-paste(colnames(asv.res)[15:19], "1", sep = ".")
colnames(asv.res)[25:29]<-paste(colnames(asv.res)[25:29], "2", sep = ".")
colnames(asv.res)[35:39]<-paste(colnames(asv.res)[35:39], "3", sep = ".")
colnames(asv.res)[45:49]<-paste(colnames(asv.res)[45:49], "4", sep = ".")
colnames(asv.res)[55:59]<-paste(colnames(asv.res)[55:59], "5", sep = ".")
colnames(asv.res)[65:69]<-paste(colnames(asv.res)[65:69], "6", sep = ".")

# Phylum
# Find all taxa that exhibit sig. differential abundance #
sig.lfc.phy1<-phy.res[ which(phy.res$padj.1<0.05 & phy.res$we.eBH_IIR.PR.1<0.05),]
sig.lfc.phy2<-phy.res[ which(phy.res$padj.2<0.05 & phy.res$we.eBH_IL.PL.2<0.05),]
sig.lfc.phy3<-phy.res[ which(phy.res$padj.3<0.05 & phy.res$we.eBH_IIR.IR.3<0.05),]
sig.lfc.phy4<-phy.res[ which(phy.res$padj.4<0.05 & phy.res$we.eBH_IIR.UR.4<0.05),]
sig.lfc.phy5<-phy.res[ which(phy.res$padj.5<0.05 & phy.res$we.eBH_IUR.UR.5<0.05),]
sig.lfc.phy6<-phy.res[ which(phy.res$padj.6<0.05 & phy.res$we.eBH_IL.UL.6<0.05),]

# Family
# Find all taxa that exhibit sig. differential abundance
sig.lfc.fam1<-fam.res[ which(fam.res$padj.1<0.05 & fam.res$wi.eBH_IIR.PR.1<0.05),]
sig.lfc.fam2<-fam.res[ which(fam.res$padj.2<0.05 & fam.res$wi.eBH_IL.PL.2<0.05),]
sig.lfc.fam3<-fam.res[ which(fam.res$padj.3<0.05 & fam.res$wi.eBH_IIR.IR.3<0.05),]
sig.lfc.fam4<-fam.res[ which(fam.res$padj.4<0.05 & fam.res$wi.eBH_IIR.UR.4<0.05),]
sig.lfc.fam5<-fam.res[ which(fam.res$padj.5<0.05 & fam.res$wi.eBH_IUR.UR.5<0.05),]
sig.lfc.fam6<-fam.res[ which(fam.res$padj.6<0.05 & fam.res$wi.eBH_IL.UL.6<0.05),]
sig.lfc.fam<-rbind(sig.lfc.fam1,sig.lfc.fam2)

sig.lfc.fam$lfc.1<-ifelse(sig.lfc.fam$we.eBH_IIR.PR.1<0.05,sig.lfc.fam$lfc.1,0)
sig.lfc.fam$lfc.2<-ifelse(sig.lfc.fam$we.eBH_IL.PL.2<0.05,sig.lfc.fam$lfc.2,0)
sig.lfc.fam$lfc.3<-ifelse(sig.lfc.fam$we.eBH_IIR.IR.3<0.05,sig.lfc.fam$lfc.3,0)
sig.lfc.fam$lfc.4<-ifelse(sig.lfc.fam$we.eBH_IIR.UR.4<0.05,sig.lfc.fam$lfc.4,0)
sig.lfc.fam$lfc.5<-ifelse(sig.lfc.fam$we.eBH_IUR.UR.5<0.05,sig.lfc.fam$lfc.5,0)
sig.lfc.fam$lfc.6<-ifelse(sig.lfc.fam$we.eBH_IL.UL.6<0.05,sig.lfc.fam$lfc.6,0)

sig.fam<-as.character(sig.lfc.fam$Family)
# Prune bacterial tree to just sig. taxa
sig.fam.tree <- prune_taxa(sig.fam, ps6.fam)
p <- ggtree(phy_tree(sig.fam.tree)) + geom_tiplab(size=3.5,align=T)
# get LFC into heatmap form
sig.lfc.fam.hm<-sig.lfc.fam[,c(10,20,30,40,50,60)]
sig.lfc.fam.hm.1<-as.data.frame(ifelse(sig.lfc.fam.hm>0,1,
                                       ifelse(sig.lfc.fam.hm<0,-1,NA)))
rownames(sig.lfc.fam.hm.1)<-sig.lfc.fam$Family

pdf(file="DA.fam.overlap.pdf",
    width=10,height=8,bg = "transparent")
gheatmap(p, sig.lfc.fam.hm.1, offset = 1.6, width=0.5, font.size=3, 
         colnames_angle=-45, hjust=0,low = "#b2182b", high = "#2166ac")
dev.off()

# very few taxa found to be sig differentially abundant by ALDEx2
# How about just DESeq2
sig.lfc.fam<-fam.res[ which(fam.res$padj.1<0.05 | fam.res$padj.2 <0.05 |
                              fam.res$padj.3<0.05 | fam.res$padj.4<0.05 |
                              fam.res$padj.5<0.05 | fam.res$padj.6<0.05),]

# set non-sig. LFC to 0 
sig.lfc.fam$lfc.1<-ifelse(sig.lfc.fam$padj.1<0.05,sig.lfc.fam$lfc.1,0)
sig.lfc.fam$lfc.2<-ifelse(sig.lfc.fam$padj.2<0.05,sig.lfc.fam$lfc.2,0)
sig.lfc.fam$lfc.3<-ifelse(sig.lfc.fam$padj.3<0.05,sig.lfc.fam$lfc.3,0)
sig.lfc.fam$lfc.4<-ifelse(sig.lfc.fam$padj.4<0.05,sig.lfc.fam$lfc.4,0)
sig.lfc.fam$lfc.5<-ifelse(sig.lfc.fam$padj.5<0.05,sig.lfc.fam$lfc.5,0)
sig.lfc.fam$lfc.6<-ifelse(sig.lfc.fam$padj.6<0.05,sig.lfc.fam$lfc.6,0)

sig.fam<-as.character(sig.lfc.fam$Family)
# Prune bacterial tree to just sig. taxa
sig.fam.tree <- prune_taxa(sig.fam, ps6.fam)
p <- ggtree(phy_tree(sig.fam.tree)) + geom_tiplab(size=3.5,align=T)
# get LFC into heatmap form
sig.lfc.fam.hm<-sig.lfc.fam[,c(10,20,30,40,50,60)]
sig.lfc.fam.hm.1<-as.data.frame(ifelse(sig.lfc.fam.hm>0,1,
                                       ifelse(sig.lfc.fam.hm<0,-1,NA)))
rownames(sig.lfc.fam.hm.1)<-sig.lfc.fam$Family

pdf(file="DA.fam.pdf",
    width=10,height=8,bg = "transparent")
gheatmap(p, sig.lfc.fam.hm.1, offset = 1.6, width=0.5, font.size=3, 
         colnames_angle=-45, hjust=0,low = "#b2182b", high = "#2166ac")
dev.off()

# Genus
# Find all taxa that exhibit sig. differential abundance
sig.lfc.gen1<-gen.res[ which(gen.res$padj.1<0.05 & gen.res$we.eBH_IIR.PR.1<0.05),]
sig.lfc.gen2<-gen.res[ which(gen.res$padj.2<0.05 & gen.res$we.eBH_IL.PL.2<0.05),]
sig.lfc.gen3<-gen.res[ which(gen.res$padj.3<0.05 & gen.res$we.eBH_IIR.IR.3<0.05),]
sig.lfc.gen4<-gen.res[ which(gen.res$padj.4<0.05 & gen.res$we.eBH_IIR.UR.4<0.05),]
sig.lfc.gen5<-gen.res[ which(gen.res$padj.5<0.05 & gen.res$we.eBH_IUR.UR.5<0.05),]
sig.lfc.gen6<-gen.res[ which(gen.res$padj.6<0.05 & gen.res$we.eBH_IL.UL.6<0.05),]

sig.lfc.gen<-rbind(sig.lfc.gen1,sig.lfc.gen2)

sig.lfc.gen$lfc.1<-ifelse(sig.lfc.gen$we.eBH_IIR.PR.1<0.05,sig.lfc.gen$lfc.1,0)
sig.lfc.gen$lfc.2<-ifelse(sig.lfc.gen$we.eBH_IL.PL.2<0.05,sig.lfc.gen$lfc.2,0)
sig.lfc.gen$lfc.3<-ifelse(sig.lfc.gen$we.eBH_IIR.IR.3<0.05,sig.lfc.gen$lfc.3,0)
sig.lfc.gen$lfc.4<-ifelse(sig.lfc.gen$we.eBH_IIR.UR.4<0.05,sig.lfc.gen$lfc.4,0)
sig.lfc.gen$lfc.5<-ifelse(sig.lfc.gen$we.eBH_IUR.UR.5<0.05,sig.lfc.gen$lfc.5,0)
sig.lfc.gen$lfc.6<-ifelse(sig.lfc.gen$we.eBH_IL.UL.6<0.05,sig.lfc.gen$lfc.6,0)

sig.gen<-as.character(sig.lfc.gen$Genus)
# Prune bacterial tree to just sig. taxa
sig.gen.tree <- prune_taxa(sig.gen, ps6.gen)
p <- ggtree(phy_tree(sig.gen.tree)) + geom_tiplab(size=3.5,align=T)
# get LFC into heatmap form
sig.lfc.gen.hm<-sig.lfc.gen[,c(10,20,30,40,50,60)]
sig.lfc.gen.hm.1<-as.data.frame(ifelse(sig.lfc.gen.hm>0,1,
                                       ifelse(sig.lfc.gen.hm<0,-1,NA)))
rownames(sig.lfc.gen.hm.1)<-sig.lfc.gen$Genus

pdf(file="DA.gen.overlap.pdf",
    width=10,height=8,bg = "transparent")
gheatmap(p, sig.lfc.gen.hm.1, offset = 1.6, width=0.5, font.size=3, 
         colnames_angle=-45, hjust=0,low = "#b2182b", high = "#2166ac")
dev.off()

# very few taxa found to be sig differentially abundant by ALDEx2
# How about just DESeq2
sig.lfc.gen<-gen.res[ which(gen.res$padj.1<0.05 | gen.res$padj.2 <0.05 |
                              gen.res$padj.3<0.05 | gen.res$padj.4<0.05 |
                              gen.res$padj.5<0.05 | gen.res$padj.6<0.05),]
# set non-sig. LFC to 0
sig.lfc.gen$lfc.1<-ifelse(sig.lfc.gen$padj.1<0.05,sig.lfc.gen$lfc.1,0)
sig.lfc.gen$lfc.2<-ifelse(sig.lfc.gen$padj.2<0.05,sig.lfc.gen$lfc.2,0)
sig.lfc.gen$lfc.3<-ifelse(sig.lfc.gen$padj.3<0.05,sig.lfc.gen$lfc.3,0)
sig.lfc.gen$lfc.4<-ifelse(sig.lfc.gen$padj.4<0.05,sig.lfc.gen$lfc.4,0)
sig.lfc.gen$lfc.5<-ifelse(sig.lfc.gen$padj.5<0.05,sig.lfc.gen$lfc.5,0)
sig.lfc.gen$lfc.6<-ifelse(sig.lfc.gen$padj.6<0.05,sig.lfc.gen$lfc.6,0)

sig.gen<-as.character(sig.lfc.gen$Genus)
# Prune bacterial tree to just sig. taxa
sig.gen.tree <- prune_taxa(sig.gen, ps6.gen)
p <- ggtree(phy_tree(sig.gen.tree)) + geom_tiplab(size=2, align = T)
# get LFC into heatmap form
sig.lfc.gen.hm<-sig.lfc.gen[,c(10,20,30,40,50,60)]
sig.lfc.gen.hm.1<-as.data.frame(ifelse(sig.lfc.gen.hm>0,1,
                                       ifelse(sig.lfc.gen.hm<0,-1,NA)))
rownames(sig.lfc.gen.hm.1)<-sig.lfc.gen$Genus

pdf(file="DA.gen.pdf",
    width=10,height=8,bg = "transparent")
gheatmap(p, sig.lfc.gen.hm.1, offset = 0.7, width=0.5, font.size=3, 
         colnames_angle=-45, hjust=0,low = "#b2182b", high = "#2166ac")
dev.off()

# ASV
# Find all taxa that exhibit sig. differential abundance
sig.lfc.asv1<-asv.res[ which(asv.res$padj.1<0.05 & asv.res$we.eBH_IIR.PR.1<0.05),]
sig.lfc.asv2<-asv.res[ which(asv.res$padj.2<0.05 & asv.res$we.eBH_IIR.IR.2<0.05),]
sig.lfc.asv3<-asv.res[ which(asv.res$padj.3<0.05 & asv.res$we.eBH_IIR.UR.3<0.05),]
sig.lfc.asv4<-asv.res[ which(asv.res$padj.4<0.05 & asv.res$we.eBH_IL.PL.4<0.05),]
sig.lfc.asv5<-asv.res[ which(asv.res$padj.5<0.05 & asv.res$we.eBH_IL.UL.5<0.05),]

# very few taxa found to be sig differentially abundant by ALDEx2
# How about just DESeq2
sig.lfc.asv<-res.asv.total[ which(res.asv.total$padj.1<0.05 | res.asv.total$padj.2 <0.05 |
                                    res.asv.total$padj.3<0.05 | res.asv.total$padj.4<0.05 |
                                    res.asv.total$padj.5<0.05),]
# set non-sig. LFC to 0
sig.lfc.asv$lfc.1<-ifelse(sig.lfc.asv$padj.1<0.05,sig.lfc.asv$lfc.1,0)
sig.lfc.asv$lfc.2<-ifelse(sig.lfc.asv$padj.2<0.05,sig.lfc.asv$lfc.2,0)
sig.lfc.asv$lfc.3<-ifelse(sig.lfc.asv$padj.3<0.05,sig.lfc.asv$lfc.3,0)
sig.lfc.asv$lfc.4<-ifelse(sig.lfc.asv$padj.4<0.05,sig.lfc.asv$lfc.4,0)
sig.lfc.asv$lfc.5<-ifelse(sig.lfc.asv$padj.5<0.05,sig.lfc.asv$lfc.5,0)

sig.asv<-rownames(sig.lfc.asv)
# Prune bacterial tree to just sig. taxa
sig.asv.tree <- prune_taxa(sig.asv, ps6)
p <- ggtree(phy_tree(sig.asv.tree))
# get LFC into heatmap form
sig.lfc.asv.hm<-sig.lfc.asv[,c(7,12,17,22,27)]
sig.lfc.asv.hm.1<-as.data.frame(ifelse(sig.lfc.asv.hm>0,1,
                                       ifelse(sig.lfc.asv.hm<0,-1,NA)))

pdf(file="DA.asv.pdf",
    width=10,height=8,bg = "transparent")
gheatmap(p, sig.lfc.asv.hm.1, offset = 0, width=0.5, font.size=3, 
         colnames_angle=-45, hjust=0,low = "#b2182b", high = "#2166ac")
dev.off()

#### Taxonomy proportions affected by factors ##
sig.lfc.phy1<-phy.res[ which(phy.res$padj.1<0.05 | phy.res$wi.eBH_IIR.PR.1<0.05),]
sig.lfc.phy2<-phy.res[ which(phy.res$padj.2<0.05 | phy.res$wi.eBH_IL.PL.2<0.05),]
sig.lfc.phy3<-phy.res[ which(phy.res$padj.3<0.05 | phy.res$wi.eBH_IIR.IR.3<0.05),]
sig.lfc.phy4<-phy.res[ which(phy.res$padj.4<0.05 | phy.res$wi.eBH_IIR.UR.4<0.05),]
sig.lfc.phy5<-phy.res[ which(phy.res$padj.5<0.05 | phy.res$wi.eBH_IUR.UR.5<0.05),]
sig.lfc.phy6<-phy.res[ which(phy.res$padj.6<0.05 | phy.res$wi.eBH_IL.UL.6<0.05),]

sig.lfc.cla1<-cla.res[ which(cla.res$padj.1<0.05 | cla.res$wi.eBH_IIR.PR.1<0.05),]
sig.lfc.cla2<-cla.res[ which(cla.res$padj.2<0.05 | cla.res$wi.eBH_IL.PL.2<0.05),]
sig.lfc.cla3<-cla.res[ which(cla.res$padj.3<0.05 | cla.res$wi.eBH_IIR.IR.3<0.05),]
sig.lfc.cla4<-cla.res[ which(cla.res$padj.4<0.05 | cla.res$wi.eBH_IIR.UR.4<0.05),]
sig.lfc.cla5<-cla.res[ which(cla.res$padj.5<0.05 | cla.res$wi.eBH_IUR.UR.5<0.05),]
sig.lfc.cla6<-cla.res[ which(cla.res$padj.6<0.05 | cla.res$wi.eBH_IL.UL.6<0.05),]

sig.lfc.ord1<-ord.res[ which(ord.res$padj.1<0.05 | ord.res$wi.eBH_IIR.PR.1<0.05),]
sig.lfc.ord2<-ord.res[ which(ord.res$padj.2<0.05 | ord.res$wi.eBH_IL.PL.2<0.05),]
sig.lfc.ord3<-ord.res[ which(ord.res$padj.3<0.05 | ord.res$wi.eBH_IIR.IR.3<0.05),]
sig.lfc.ord4<-ord.res[ which(ord.res$padj.4<0.05 | ord.res$wi.eBH_IIR.UR.4<0.05),]
sig.lfc.ord5<-ord.res[ which(ord.res$padj.5<0.05 | ord.res$wi.eBH_IUR.UR.5<0.05),]
sig.lfc.ord6<-ord.res[ which(ord.res$padj.6<0.05 | ord.res$wi.eBH_IL.UL.6<0.05),]

sig.lfc.fam1<-fam.res[ which(fam.res$padj.1<0.05 | fam.res$wi.eBH_IIR.PR.1<0.05),]
sig.lfc.fam2<-fam.res[ which(fam.res$padj.2<0.05 | fam.res$wi.eBH_IL.PL.2<0.05),]
sig.lfc.fam3<-fam.res[ which(fam.res$padj.3<0.05 | fam.res$wi.eBH_IIR.IR.3<0.05),]
sig.lfc.fam4<-fam.res[ which(fam.res$padj.4<0.05 | fam.res$wi.eBH_IIR.UR.4<0.05),]
sig.lfc.fam5<-fam.res[ which(fam.res$padj.5<0.05 | fam.res$wi.eBH_IUR.UR.5<0.05),]
sig.lfc.fam6<-fam.res[ which(fam.res$padj.6<0.05 | fam.res$wi.eBH_IL.UL.6<0.05),]

sig.lfc.gen1<-gen.res[ which(gen.res$padj.1<0.05 | gen.res$wi.eBH_IIR.PR.1<0.05),]
sig.lfc.gen2<-gen.res[ which(gen.res$padj.2<0.05 | gen.res$wi.eBH_IL.PL.2<0.05),]
sig.lfc.gen3<-gen.res[ which(gen.res$padj.3<0.05 | gen.res$wi.eBH_IIR.IR.3<0.05),]
sig.lfc.gen4<-gen.res[ which(gen.res$padj.4<0.05 | gen.res$wi.eBH_IIR.UR.4<0.05),]
sig.lfc.gen5<-gen.res[ which(gen.res$padj.5<0.05 | gen.res$wi.eBH_IUR.UR.5<0.05),]
sig.lfc.gen6<-gen.res[ which(gen.res$padj.6<0.05 | gen.res$wi.eBH_IL.UL.6<0.05),]

sig.lfc.asv1<-asv.res[ which(asv.res$padj.1<0.05 | asv.res$wi.eBH_IIR.PR.1<0.05),]
sig.lfc.asv2<-asv.res[ which(asv.res$padj.2<0.05 | asv.res$wi.eBH_IL.PL.2<0.05),]
sig.lfc.asv3<-asv.res[ which(asv.res$padj.3<0.05 | asv.res$wi.eBH_IIR.IR.3<0.05),]
sig.lfc.asv4<-asv.res[ which(asv.res$padj.4<0.05 | asv.res$wi.eBH_IIR.UR.4<0.05),]
sig.lfc.asv5<-asv.res[ which(asv.res$padj.5<0.05 | asv.res$wi.eBH_IUR.UR.5<0.05),]
sig.lfc.asv6<-asv.res[ which(asv.res$padj.6<0.05 | asv.res$wi.eBH_IL.UL.6<0.05),]

prop.phy.1<-nrow(sig.lfc.phy1)/nrow(phy.res[phy.res$root.ra>0,])
prop.phy.2<-nrow(sig.lfc.phy2)/nrow(phy.res[complete.cases(phy.res[ ,9]),])
prop.phy.3<-nrow(sig.lfc.phy3)/nrow(phy.res[phy.res$root.ra>0,])
prop.phy.4<-nrow(sig.lfc.phy4)/nrow(phy.res[phy.res$root.ra>0,])
prop.phy.5<-nrow(sig.lfc.phy5)/nrow(phy.res[phy.res$root.ra>0,])
prop.phy.6<-nrow(sig.lfc.phy6)/nrow(phy.res[complete.cases(phy.res[ ,9]),])

prop.cla.1<-nrow(sig.lfc.cla1)/nrow(cla.res[cla.res$root.ra>0,])
prop.cla.2<-nrow(sig.lfc.cla2)/nrow(cla.res[complete.cases(cla.res[ ,9]),])
prop.cla.3<-nrow(sig.lfc.cla3)/nrow(cla.res[cla.res$root.ra>0,])
prop.cla.4<-nrow(sig.lfc.cla4)/nrow(cla.res[cla.res$root.ra>0,])
prop.cla.5<-nrow(sig.lfc.cla5)/nrow(cla.res[cla.res$root.ra>0,])
prop.cla.6<-nrow(sig.lfc.cla6)/nrow(cla.res[complete.cases(cla.res[ ,9]),])

prop.ord.1<-nrow(sig.lfc.ord1)/nrow(ord.res[ord.res$root.ra>0,])
prop.ord.2<-nrow(sig.lfc.ord2)/nrow(ord.res[complete.cases(ord.res[ ,9]),])
prop.ord.3<-nrow(sig.lfc.ord3)/nrow(ord.res[ord.res$root.ra>0,])
prop.ord.4<-nrow(sig.lfc.ord4)/nrow(ord.res[ord.res$root.ra>0,])
prop.ord.5<-nrow(sig.lfc.ord5)/nrow(ord.res[ord.res$root.ra>0,])
prop.ord.6<-nrow(sig.lfc.ord6)/nrow(ord.res[complete.cases(ord.res[ ,9]),])

prop.fam.1<-nrow(sig.lfc.fam1)/nrow(fam.res[fam.res$root.ra>0,])
prop.fam.2<-nrow(sig.lfc.fam2)/nrow(fam.res[complete.cases(fam.res[ ,9]),])
prop.fam.3<-nrow(sig.lfc.fam3)/nrow(fam.res[fam.res$root.ra>0,])
prop.fam.4<-nrow(sig.lfc.fam4)/nrow(fam.res[fam.res$root.ra>0,])
prop.fam.5<-nrow(sig.lfc.fam5)/nrow(fam.res[fam.res$root.ra>0,])
prop.fam.6<-nrow(sig.lfc.fam6)/nrow(fam.res[complete.cases(fam.res[ ,9]),])

prop.gen.1<-nrow(sig.lfc.gen1)/nrow(gen.res[gen.res$root.ra>0,])
prop.gen.2<-nrow(sig.lfc.gen2)/nrow(gen.res[complete.cases(gen.res[ ,9]),])
prop.gen.3<-nrow(sig.lfc.gen3)/nrow(gen.res[gen.res$root.ra>0,])
prop.gen.4<-nrow(sig.lfc.gen4)/nrow(gen.res[gen.res$root.ra>0,])
prop.gen.5<-nrow(sig.lfc.gen5)/nrow(gen.res[gen.res$root.ra>0,])
prop.gen.6<-nrow(sig.lfc.gen6)/nrow(gen.res[complete.cases(gen.res[ ,9]),])

prop.asv.1<-nrow(sig.lfc.asv1)/nrow(asv.res[asv.res$root.ra>0,])
prop.asv.2<-nrow(sig.lfc.asv2)/nrow(asv.res[complete.cases(asv.res[ ,9]),])
prop.asv.3<-nrow(sig.lfc.asv3)/nrow(asv.res[asv.res$root.ra>0,])
prop.asv.4<-nrow(sig.lfc.asv4)/nrow(asv.res[asv.res$root.ra>0,])
prop.asv.5<-nrow(sig.lfc.asv5)/nrow(asv.res[asv.res$root.ra>0,])
prop.asv.6<-nrow(sig.lfc.asv6)/nrow(asv.res[complete.cases(asv.res[ ,9]),])

sig.prop<-c(prop.phy.1,prop.cla.1,prop.ord.1,prop.fam.1,prop.gen.1,prop.asv.1,
            prop.phy.2,prop.cla.2,prop.ord.2,prop.fam.2,prop.gen.2,prop.asv.2,
            prop.phy.3,prop.cla.3,prop.ord.3,prop.fam.3,prop.gen.3,prop.asv.3,
            prop.phy.4,prop.cla.4,prop.ord.4,prop.fam.4,prop.gen.4,prop.asv.4,
            prop.phy.5,prop.cla.5,prop.ord.5,prop.fam.5,prop.gen.5,prop.asv.5,
            prop.phy.6,prop.cla.6,prop.ord.6,prop.fam.6,prop.gen.6,prop.asv.6)
taxon.prop<-rep(c("Phylum","Class","Order","Family","Genus","ASV"),6)
factor.prop<-c(rep("PR_IIR",6),rep("PL_IL",6),rep("IUR_IIR",6),rep("UR_IIR",6),rep("UR_IUR",6),rep("UL_IL",6))
ppm.prop<-as.data.frame(cbind(factor.prop,taxon.prop,sig.prop))
ppm.prop$taxon.prop<-factor(ppm.prop$taxon.prop,
                                levels = c("Phylum","Class","Order","Family",
                                           "Genus","ASV", ordered=T))
ppm.prop$sig.prop<-as.numeric(as.character(ppm.prop$sig.prop))

pdf(file="DA.PR.IIR.pdf",
    width=6,height=11,bg = "transparent")
ggplot(data=ppm.prop[c(1:6),], aes(x=taxon.prop, y=sig.prop))+
  geom_bar(stat="identity",fill="black",colour="black",alpha=0.8,width=0.9,position = position_dodge(width = 0.85)) + 
  theme.gg + coord_cartesian(ylim=c(0,1))
dev.off()

pdf(file="DA.PL.IL.pdf",
    width=6,height=11,bg = "transparent")
ggplot(data=ppm.prop[c(7:12),], aes(x=taxon.prop, y=sig.prop))+
  geom_bar(stat="identity",fill="black",colour="black",alpha=0.8,width=0.9,position = position_dodge(width = 0.85)) + 
  theme.gg + coord_cartesian(ylim=c(0,1))
dev.off()

pdf(file="DA.IUR.IIR.pdf",
    width=6,height=11,bg = "transparent")
ggplot(data=ppm.prop[c(13:18),], aes(x=taxon.prop, y=sig.prop))+
  geom_bar(stat="identity",fill="black",colour="black",alpha=0.8,width=0.9,position = position_dodge(width = 0.85)) + 
  theme.gg + coord_cartesian(ylim=c(0,1))
dev.off()

pdf(file="DA.UR.IIR.pdf",
    width=6,height=11,bg = "transparent")
ggplot(data=ppm.prop[c(19:24),], aes(x=taxon.prop, y=sig.prop))+
  geom_bar(stat="identity",fill="black",colour="black",alpha=0.8,width=0.9,position = position_dodge(width = 0.85)) + 
  theme.gg + coord_cartesian(ylim=c(0,1))
dev.off()

pdf(file="DA.UR.IUR.pdf",
    width=6,height=11,bg = "transparent")
ggplot(data=ppm.prop[c(25:30),], aes(x=taxon.prop, y=sig.prop))+
  geom_bar(stat="identity",fill="black",colour="black",alpha=0.8,width=0.9,position = position_dodge(width = 0.85)) + 
  theme.gg + coord_cartesian(ylim=c(0,1))
dev.off()

pdf(file="DA.UL.IL.pdf",
    width=6,height=11,bg = "transparent")
ggplot(data=ppm.prop[c(26:31),], aes(x=taxon.prop, y=sig.prop))+
  geom_bar(stat="identity",fill="black",colour="black",alpha=0.8,width=0.9,position = position_dodge(width = 0.85)) + 
  theme.gg + coord_cartesian(ylim=c(0,1))
dev.off()

#### 4) Network Analysis ####
# begin with ps5 = control samples removed and threshold at 5 x 25
ps.nt<-ps5
taxa_names(ps.nt)<-c(1:2241)
# remove ASVs not found in at least 50% of samples with at least 10 reads for each community type
# IIR #
ps.nt.iir<-subset_samples(ps.nt,IDtype=="I.IR") 
threshold50<-kOverA(6,A=10)
threshold25<-kOverA(3,A=10)
ps.nt.iir.50<-filter_taxa(ps.nt.iir,threshold50,TRUE) # 229 taxa
sum(taxa_sums(ps.nt.iir.50))/sum(taxa_sums(ps.nt.iir)) # 63%
ps.nt.iir.25<-filter_taxa(ps.nt.iir,threshold25,TRUE) # 723 taxa
sum(taxa_sums(ps.nt.iir.25))/sum(taxa_sums(ps.nt.iir)) # 87%
ps.iir.nt.50.clr<-as.numeric(rowMeans(clr(otu_table(ps.nt.iir.50), 1)))
min(ps.iir.nt.50.clr)
size.iir.50<-ps.iir.nt.50.clr+1.300738
ps.iir.nt.25.clr<-as.numeric(rowMeans(clr(otu_table(ps.nt.iir.25), 1)))
min(ps.iir.nt.25.clr)
size.iir.25<-ps.iir.nt.25.clr+1.05311
taxa.ps.iir.nt.50<-as.data.frame(tax_table(ps.nt.iir.50))
taxa.ps.iir.nt.25<-as.data.frame(tax_table(ps.nt.iir.25))

# IR #
ps.nt.ir<-subset_samples(ps.nt,IDtype=="I.R") 
threshold50<-kOverA(6,A=10)
threshold25<-kOverA(3,A=10)
ps.nt.ir.50<-filter_taxa(ps.nt.ir,threshold50,TRUE) # 285 taxa
sum(taxa_sums(ps.nt.ir.50))/sum(taxa_sums(ps.nt.ir)) # 61%
ps.nt.ir.25<-filter_taxa(ps.nt.ir,threshold25,TRUE) # 854 taxa
sum(taxa_sums(ps.nt.ir.25))/sum(taxa_sums(ps.nt.ir)) # 88%
ps.ir.nt.50.clr<-as.numeric(rowMeans(clr(otu_table(ps.nt.ir.50), 1)))
min(ps.ir.nt.50.clr)
size.ir.50<-ps.ir.nt.50.clr+1.143919
ps.ir.nt.25.clr<-as.numeric(rowMeans(clr(otu_table(ps.nt.ir.25), 1)))
min(ps.ir.nt.25.clr)
size.ir.25<-ps.ir.nt.25.clr+1.24686
taxa.ps.ir.nt.50<-as.data.frame(tax_table(ps.nt.ir.50))
taxa.ps.ir.nt.25<-as.data.frame(tax_table(ps.nt.ir.25))

# UR #
ps.nt.ur<-subset_samples(ps.nt,IDtype=="U.R") 
threshold50<-kOverA(6,A=10)
threshold25<-kOverA(3,A=10)
ps.nt.ur.50<-filter_taxa(ps.nt.ur,threshold50,TRUE) # 289 taxa
sum(taxa_sums(ps.nt.ur.50))/sum(taxa_sums(ps.nt.ur)) # 60%
ps.nt.ur.25<-filter_taxa(ps.nt.ur,threshold25,TRUE) # 828 taxa
sum(taxa_sums(ps.nt.ur.25))/sum(taxa_sums(ps.nt.ur)) # 88%
ps.ur.nt.50.clr<-as.numeric(rowMeans(clr(otu_table(ps.nt.ur.50), 1)))
min(ps.ur.nt.50.clr)
size.ur.50<-ps.ur.nt.50.clr+1.086457
ps.ur.nt.25.clr<-as.numeric(rowMeans(clr(otu_table(ps.nt.ur.25), 1)))
min(ps.ur.nt.25.clr)
size.ur.25<-ps.ur.nt.25.clr+0.8836636
taxa.ps.ur.nt.50<-as.data.frame(tax_table(ps.nt.ur.50))
taxa.ps.ur.nt.25<-as.data.frame(tax_table(ps.nt.ur.25))

# PR #
ps.nt.pr<-subset_samples(ps.nt,IDtype=="P.R") 
threshold50<-kOverA(6,A=10)
threshold25<-kOverA(3,A=10)
ps.nt.pr.50<-filter_taxa(ps.nt.pr,threshold50,TRUE) # 113 taxa
sum(taxa_sums(ps.nt.pr.50))/sum(taxa_sums(ps.nt.pr)) # 51%
ps.nt.pr.25<-filter_taxa(ps.nt.pr,threshold25,TRUE) # 378 taxa
sum(taxa_sums(ps.nt.pr.25))/sum(taxa_sums(ps.nt.pr)) # 86%
ps.pr.nt.50.clr<-as.numeric(rowMeans(clr(otu_table(ps.nt.pr.50), 1)))
min(ps.pr.nt.50.clr)
size.pr.50<-ps.pr.nt.50.clr+1.11405
ps.pr.nt.25.clr<-as.numeric(rowMeans(clr(otu_table(ps.nt.pr.25), 1)))
min(ps.pr.nt.25.clr)
size.pr.25<-ps.pr.nt.25.clr+1.517776
taxa.ps.pr.nt.50<-as.data.frame(tax_table(ps.nt.pr.50))
taxa.ps.pr.nt.25<-as.data.frame(tax_table(ps.nt.pr.25))

# IL #
ps.nt.il<-subset_samples(ps.nt,IDtype=="I.L") 
threshold50<-kOverA(6,A=10)
threshold25<-kOverA(3,A=10)
ps.nt.il.50<-filter_taxa(ps.nt.il,threshold50,TRUE) # 11 taxa
sum(taxa_sums(ps.nt.il.50))/sum(taxa_sums(ps.nt.il)) # 50%
ps.nt.il.25<-filter_taxa(ps.nt.il,threshold25,TRUE) # 37 taxa
sum(taxa_sums(ps.nt.il.25))/sum(taxa_sums(ps.nt.il)) # 87%
ps.il.nt.50.clr<-as.numeric(rowMeans(clr(otu_table(ps.nt.il.50), 1)))
min(ps.il.nt.50.clr)
size.il.50<-ps.il.nt.50.clr+0.6632998
ps.il.nt.25.clr<-as.numeric(rowMeans(clr(otu_table(ps.nt.il.25), 1)))
min(ps.il.nt.25.clr)
size.il.25<-ps.il.nt.25.clr+0.553652
taxa.ps.il.nt.50<-as.data.frame(tax_table(ps.nt.il.50))
taxa.ps.il.nt.25<-as.data.frame(tax_table(ps.nt.il.25))

# UL #
ps.nt.ul<-subset_samples(ps.nt,IDtype=="U.L") 
threshold50<-kOverA(6,A=10)
threshold25<-kOverA(3,A=10)
ps.nt.ul.50<-filter_taxa(ps.nt.ul,threshold50,TRUE) # 10 taxa
sum(taxa_sums(ps.nt.ul.50))/sum(taxa_sums(ps.nt.ul)) # 54%
ps.nt.ul.25<-filter_taxa(ps.nt.ul,threshold25,TRUE) # 40 taxa
sum(taxa_sums(ps.nt.ul.25))/sum(taxa_sums(ps.nt.ul)) # 86%
ps.ul.nt.50.clr<-as.numeric(rowMeans(clr(otu_table(ps.nt.ul.50), 1)))
min(ps.ul.nt.50.clr)
size.ul.50<-ps.ul.nt.50.clr+0.5410275
ps.ul.nt.25.clr<-as.numeric(rowMeans(clr(otu_table(ps.nt.ul.25), 1)))
min(ps.ul.nt.25.clr)
size.ul.25<-ps.ul.nt.25.clr+0.3497377
taxa.ps.ul.nt.50<-as.data.frame(tax_table(ps.nt.ul.50))
taxa.ps.ul.nt.25<-as.data.frame(tax_table(ps.nt.ul.25))

# PL #
ps.nt.pl<-subset_samples(ps.nt,IDtype=="P.L") 
threshold50<-kOverA(6,A=10)
threshold25<-kOverA(3,A=10)
ps.nt.pl.50<-filter_taxa(ps.nt.pl,threshold50,TRUE) # 23 taxa
sum(taxa_sums(ps.nt.pl.50))/sum(taxa_sums(ps.nt.pl)) # 56%
ps.nt.pl.25<-filter_taxa(ps.nt.pl,threshold25,TRUE) # 72 taxa
sum(taxa_sums(ps.nt.pl.25))/sum(taxa_sums(ps.nt.pl)) # 94%
ps.pl.nt.50.clr<-as.numeric(rowMeans(clr(otu_table(ps.nt.pl.50), 1)))
min(ps.pl.nt.50.clr)
size.pl.50<-ps.pl.nt.50.clr+1.29306
ps.pl.nt.25.clr<-as.numeric(rowMeans(clr(otu_table(ps.nt.pl.25), 1)))
min(ps.pl.nt.25.clr)
size.pl.25<-ps.pl.nt.25.clr+1.030113
taxa.ps.pl.nt.50<-as.data.frame(tax_table(ps.nt.pl.50))
taxa.ps.pl.nt.25<-as.data.frame(tax_table(ps.nt.pl.25))

#### 4.1) Network inference ####
#### 4.1.1) SPIEC-EASI ####
#### Root ####
# IIR network #
se.iir <- spiec.easi(ps.nt.iir.50, method='mb',lambda.min.ratio=2e-3,
                     nlambda=100, icov.select.params=list(rep.num=30))
ig.iir <- adj2igraph(se.iir$refit,vertex.attr=list(name=taxa_names(ps.nt.iir.50)))

print(se.iir$variability)
betaMat.iir <- as.matrix(symBeta(getOptBeta(se.iir)))
positive.iir <- length(betaMat.iir[betaMat.iir>0])/2 
negative.iir <- length(betaMat.iir[betaMat.iir<0])/2 
total.iir <- length(betaMat.iir[betaMat.iir!=0])/2 
total.iir.prop<-total.iir/(length(betaMat.iir)/2)
positive.iir/total.iir
negative.iir/total.iir

otu.ids.iir<-colnames(se.iir$data)
edges.iir<-E(ig.iir)
edge.colors.iir<-c()
for(e.index in 1:length(edges.iir)){
  adj.nodes.iir=ends(ig.iir,edges.iir[e.index])
  xindex=which(otu.ids.iir==adj.nodes.iir[1])
  yindex=which(otu.ids.iir==adj.nodes.iir[2])
  beta=betaMat.iir[xindex,yindex]
  if(beta>0){
    edge.colors.iir=append(edge.colors.iir,"forestgreen")
  }else if(beta<0){
    edge.colors.iir=append(edge.colors.iir,"red")
  }
}

edge.values.iir<-c()
for(e.index in 1:length(edges.iir)){
  adj.nodes.iir=ends(ig.iir,edges.iir[e.index])
  xindex=which(otu.ids.iir==adj.nodes.iir[1])
  yindex=which(otu.ids.iir==adj.nodes.iir[2])
  beta=betaMat.iir[xindex,yindex]
  if(beta>0){
    edge.values.iir=append(edge.values.iir,beta)
  }else if(beta<0){
    edge.values.iir=append(edge.values.iir,beta)
  }
}

E(ig.iir)$color=edge.colors.iir
ig.iir1 <- ig.iir
V(ig.iir1)$name <- as.character(rownames(taxa.ps.iir.nt.50))
V(ig.iir1)$Phylum <- as.character(taxa.ps.iir.nt.50$Phylum)
V(ig.iir1)$Class <- as.character(taxa.ps.iir.nt.50$Class)
V(ig.iir1)$Order <- as.character(taxa.ps.iir.nt.50$Order)
V(ig.iir1)$Family <- as.character(taxa.ps.iir.nt.50$Family)
V(ig.iir1)$Genus <- as.character(taxa.ps.iir.nt.50$Genus)
E(ig.iir1)$arrow.size<-5
E(ig.iir1)$width <- abs(edge.values.iir) * 10
V(ig.iir1)$color <- phylum.colors[V(ig.iir1)$Phylum]
V(ig.iir1)$size <- size.iir.50 * 5
V(ig.iir1)$frame.color="black"
V(ig.iir1)$label <- taxa_names(ps.iir.nt)
vertex_attr(ig.iir1)
edge_attr(ig.iir1)

l.iir <- layout_with_fr(ig.iir1)

pdf(file="iir.network.se.50.l.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.iir1,layout=l.iir,vertex.label.font=2,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

pdf(file="iir.network.se.50.nl.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.iir1,layout=l.iir,vertex.label.font=2,vertex.label=NA,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

# IUR Network #
se.iur <- spiec.easi(ps.nt.ir.50, method='mb',lambda.min.ratio=2e-3,
                     nlambda=100, icov.select.params=list(rep.num=30))
ig.iur <- adj2igraph(se.iur$refit,vertex.attr=list(name=taxa_names(ps.nt.ir.50)))

print(se.iur$variability)
betaMat.iur <- as.matrix(symBeta(getOptBeta(se.iur)))
positive.iur <- length(betaMat.iur[betaMat.iur>0])/2 
negative.iur <- length(betaMat.iur[betaMat.iur<0])/2 
total.iur <- length(betaMat.iur[betaMat.iur!=0])/2 
total.iur.prop<-total.iur/(length(betaMat.iur)/2)
positive.iur/total.iur
negative.iur/total.iur

otu.ids.iur<-colnames(se.iur$data)
edges.iur<-E(ig.iur)
edge.colors.iur<-c()
for(e.index in 1:length(edges.iur)){
  adj.nodes.iur=ends(ig.iur,edges.iur[e.index])
  xindex=which(otu.ids.iur==adj.nodes.iur[1])
  yindex=which(otu.ids.iur==adj.nodes.iur[2])
  beta=betaMat.iur[xindex,yindex]
  if(beta>0){
    edge.colors.iur=append(edge.colors.iur,"forestgreen")
  }else if(beta<0){
    edge.colors.iur=append(edge.colors.iur,"red")
  }
}

edge.values.iur<-c()
for(e.index in 1:length(edges.iur)){
  adj.nodes.iur=ends(ig.iur,edges.iur[e.index])
  xindex=which(otu.ids.iur==adj.nodes.iur[1])
  yindex=which(otu.ids.iur==adj.nodes.iur[2])
  beta=betaMat.iur[xindex,yindex]
  if(beta>0){
    edge.values.iur=append(edge.values.iur,beta)
  }else if(beta<0){
    edge.values.iur=append(edge.values.iur,beta)
  }
}

E(ig.iur)$color=edge.colors.iur
ig.iur1 <- ig.iur
V(ig.iur1)$name <- as.character(rownames(taxa.ps.ir.nt.50))
V(ig.iur1)$Phylum <- as.character(taxa.ps.ir.nt.50$Phylum)
V(ig.iur1)$Class <- as.character(taxa.ps.ir.nt.50$Class)
V(ig.iur1)$Order <- as.character(taxa.ps.ir.nt.50$Order)
V(ig.iur1)$Family <- as.character(taxa.ps.ir.nt.50$Family)
V(ig.iur1)$Genus <- as.character(taxa.ps.ir.nt.50$Genus)
E(ig.iur1)$arrow.size<-5
E(ig.iur1)$width <- abs(edge.values.iur) * 10
V(ig.iur1)$color <- phylum.colors[V(ig.iur1)$Phylum]
V(ig.iur1)$size <- size.ir.50 * 5
V(ig.iur1)$frame.color="black"
V(ig.iur1)$label <- rownames(taxa.ps.ir.nt.50)
vertex_attr(ig.iur1)
edge_attr(ig.iur1)

l.ir <- layout_with_fr(ig.iur1)

pdf(file="iur.network.se.50.l.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.iur1,layout=l.ir,vertex.label.font=2,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

pdf(file="iur.network.se.50.nl.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.iur1,layout=l.ir,vertex.label.font=2,vertex.label=NA,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

# UR Network #
se.ur <- spiec.easi(ps.nt.ur.50, method='mb',lambda.min.ratio=2e-3,
                     nlambda=100, icov.select.params=list(rep.num=30))
ig.ur <- adj2igraph(se.ur$refit,vertex.attr=list(name=taxa_names(ps.nt.ur.50)))

print(se.ur$variability)
betaMat.ur <- as.matrix(symBeta(getOptBeta(se.ur)))
positive.ur <- length(betaMat.ur[betaMat.ur>0])/2 
negative.ur <- length(betaMat.ur[betaMat.ur<0])/2 
total.ur <- length(betaMat.ur[betaMat.ur!=0])/2 
total.ur.prop<-total.ur/(length(betaMat.ur)/2)
positive.ur/total.ur
negative.ur/total.ur

otu.ids.ur<-colnames(se.ur$data)
edges.ur<-E(ig.ur)
edge.colors.ur<-c()
for(e.index in 1:length(edges.ur)){
  adj.nodes.ur=ends(ig.ur,edges.ur[e.index])
  xindex=which(otu.ids.ur==adj.nodes.ur[1])
  yindex=which(otu.ids.ur==adj.nodes.ur[2])
  beta=betaMat.ur[xindex,yindex]
  if(beta>0){
    edge.colors.ur=append(edge.colors.ur,"forestgreen")
  }else if(beta<0){
    edge.colors.ur=append(edge.colors.ur,"red")
  }
}

edge.values.ur<-c()
for(e.index in 1:length(edges.ur)){
  adj.nodes.ur=ends(ig.ur,edges.ur[e.index])
  xindex=which(otu.ids.ur==adj.nodes.ur[1])
  yindex=which(otu.ids.ur==adj.nodes.ur[2])
  beta=betaMat.ur[xindex,yindex]
  if(beta>0){
    edge.values.ur=append(edge.values.ur,beta)
  }else if(beta<0){
    edge.values.ur=append(edge.values.ur,beta)
  }
}

E(ig.ur)$color=edge.colors.ur
ig.ur1 <- ig.ur
V(ig.ur1)$name <- as.character(rownames(taxa.ps.ur.nt.50))
V(ig.ur1)$Phylum <- as.character(taxa.ps.ur.nt.50$Phylum)
V(ig.ur1)$Class <- as.character(taxa.ps.ur.nt.50$Class)
V(ig.ur1)$Order <- as.character(taxa.ps.ur.nt.50$Order)
V(ig.ur1)$Family <- as.character(taxa.ps.ur.nt.50$Family)
V(ig.ur1)$Genus <- as.character(taxa.ps.ur.nt.50$Genus)
E(ig.ur1)$arrow.size<-5
E(ig.ur1)$width <- abs(edge.values.ur) * 10
V(ig.ur1)$color <- phylum.colors[V(ig.ur1)$Phylum]
V(ig.ur1)$size <- size.ur.50 * 5
V(ig.ur1)$frame.color="black"
V(ig.ur1)$label <- taxa_names(ps.ur.nt)
vertex_attr(ig.ur1)
edge_attr(ig.ur1)

l.ur<-layout_with_fr(ig.ur1)
  
pdf(file="ur.network.se.50.l.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.ur1,layout=l.ur,vertex.label.font=2,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

pdf(file="ur.network.se.50.nl.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.ur1,layout=l.ur,vertex.label.font=2,vertex.label=NA,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

# PR network #
se.pr <- spiec.easi(ps.nt.pr.50, method='mb',lambda.min.ratio=2e-3,
                     nlambda=100, icov.select.params=list(rep.num=30))
ig.pr <- adj2igraph(se.pr$refit,vertex.attr=list(name=taxa_names(ps.nt.pr.50)))

print(se.pr$variability)
betaMat.pr <- as.matrix(symBeta(getOptBeta(se.pr)))
positive.pr <- length(betaMat.pr[betaMat.pr>0])/2 
negative.pr <- length(betaMat.pr[betaMat.pr<0])/2 
total.pr <- length(betaMat.pr[betaMat.pr!=0])/2 
total.pr.prop<-total.pr/(length(betaMat.pr)/2)

otu.ids.pr<-colnames(se.pr$data)
edges.pr<-E(ig.pr)
edge.colors.pr<-c()
for(e.index in 1:length(edges.pr)){
  adj.nodes.pr=ends(ig.pr,edges.pr[e.index])
  xindex=which(otu.ids.pr==adj.nodes.pr[1])
  yindex=which(otu.ids.pr==adj.nodes.pr[2])
  beta=betaMat.pr[xindex,yindex]
  if(beta>0){
    edge.colors.pr=append(edge.colors.pr,"forestgreen")
  }else if(beta<0){
    edge.colors.pr=append(edge.colors.pr,"red")
  }
}

edge.values.pr<-c()
for(e.index in 1:length(edges.pr)){
  adj.nodes.pr=ends(ig.pr,edges.pr[e.index])
  xindex=which(otu.ids.pr==adj.nodes.pr[1])
  yindex=which(otu.ids.pr==adj.nodes.pr[2])
  beta=betaMat.pr[xindex,yindex]
  if(beta>0){
    edge.values.pr=append(edge.values.pr,beta)
  }else if(beta<0){
    edge.values.pr=append(edge.values.pr,beta)
  }
}

E(ig.pr)$color=edge.colors.pr
ig.pr1 <- ig.pr
V(ig.pr1)$name <- as.character(rownames(taxa.ps.pr.nt.50))
V(ig.pr1)$Phylum <- as.character(taxa.ps.pr.nt.50$Phylum)
V(ig.pr1)$Class <- as.character(taxa.ps.pr.nt.50$Class)
V(ig.pr1)$Order <- as.character(taxa.ps.pr.nt.50$Order)
V(ig.pr1)$Family <- as.character(taxa.ps.pr.nt.50$Family)
E(ig.pr1)$arrow.size<-5
E(ig.pr1)$width <- abs(edge.values.pr) * 25
V(ig.pr1)$color <- phylum.colors[V(ig.pr1)$Phylum]
V(ig.pr1)$size <- size.pr.50 * 5
V(ig.pr1)$frame.color="black"
V(ig.pr1)$label <- taxa_names(ps.pr.nt)
vertex_attr(ig.pr1)
edge_attr(ig.pr1)

l.pr<-layout_with_fr(ig.pr1)

pdf(file="pr.network.50.l.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.pr1,layout=l.pr,vertex.label.font=2,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

pdf(file="pr.network.se.50.nl.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.pr1,layout=l.pr,vertex.label.font=2,vertex.label=NA,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

#### Leaf ####
# IL network #
se.il <- spiec.easi(ps.nt.il.50, method='mb',lambda.min.ratio=2e-3,
                     nlambda=100, icov.select.params=list(rep.num=30))
ig.il <- adj2igraph(se.il$refit,vertex.attr=list(name=taxa_names(ps.nt.il.50)))

print(se.il$variability)
betaMat.il <- as.matrix(symBeta(getOptBeta(se.il)))
positive.il <- length(betaMat.il[betaMat.il>0])/2 
negative.il <- length(betaMat.il[betaMat.il<0])/2 
total.il <- length(betaMat.il[betaMat.il!=0])/2 
total.il.prop<-total.il/(length(betaMat.il)/2)

otu.ids.il<-colnames(se.il$data)
edges.il<-E(ig.il)
edge.colors.il<-c()
for(e.index in 1:length(edges.il)){
  adj.nodes.il=ends(ig.il,edges.il[e.index])
  xindex=which(otu.ids.il==adj.nodes.il[1])
  yindex=which(otu.ids.il==adj.nodes.il[2])
  beta=betaMat.il[xindex,yindex]
  if(beta>0){
    edge.colors.il=append(edge.colors.il,"forestgreen")
  }else if(beta<0){
    edge.colors.il=append(edge.colors.il,"red")
  }
}

edge.values.il<-c()
for(e.index in 1:length(edges.il)){
  adj.nodes.il=ends(ig.il,edges.il[e.index])
  xindex=which(otu.ids.il==adj.nodes.il[1])
  yindex=which(otu.ids.il==adj.nodes.il[2])
  beta=betaMat.il[xindex,yindex]
  if(beta>0){
    edge.values.il=append(edge.values.il,beta)
  }else if(beta<0){
    edge.values.il=append(edge.values.il,beta)
  }
}

E(ig.il)$color=edge.colors.il
ig.il1 <- ig.il
V(ig.il1)$name <- as.character(rownames(taxa.ps.il.nt.50))
V(ig.il1)$Phylum <- as.character(taxa.ps.il.nt.50$Phylum)
V(ig.il1)$Class <- as.character(taxa.ps.il.nt.50$Class)
V(ig.il1)$Order <- as.character(taxa.ps.il.nt.50$Order)
V(ig.il1)$Family <- as.character(taxa.ps.il.nt.50$Family)
V(ig.il1)$Genus <- as.character(taxa.ps.il.nt.50$Genus)
E(ig.il1)$arrow.size<-5
E(ig.il1)$width <- abs(edge.values.il) * 25
V(ig.il1)$color <- phylum.colors[V(ig.il1)$Phylum]
V(ig.il1)$size <- size.il.50 * 5
V(ig.il1)$frame.color="black"
V(ig.il1)$label <- taxa_names(ps.il.nt)
vertex_attr(ig.il1)
edge_attr(ig.il1)

l.il<-layout_with_fr(ig.il1)

pdf(file="il.network.se.50.l.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.il1,layout=l.il,vertex.label.font=2,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

pdf(file="il.network.se.50.nl.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.il1,layout=l.il,vertex.label.font=2,vertex.label=NA,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

# UL network #
se.ul <- spiec.easi(ps.nt.ul.50, method='mb',lambda.min.ratio=2e-3,
                    nlambda=100, icov.select.params=list(rep.num=30))
ig.ul <- adj2igraph(se.ul$refit,vertex.attr=list(name=taxa_names(ps.nt.ul.50)))

print(se.ul$variabulity)
betaMat.ul <- as.matrix(symBeta(getOptBeta(se.ul)))
positive.ul <- length(betaMat.ul[betaMat.ul>0])/2 
negative.ul <- length(betaMat.ul[betaMat.ul<0])/2 
total.ul <- length(betaMat.ul[betaMat.ul!=0])/2 
total.ul.prop<-total.ul/(length(betaMat.ul)/2)

otu.ids.ul<-colnames(se.ul$data)
edges.ul<-E(ig.ul)
edge.colors.ul<-c()
for(e.index in 1:length(edges.ul)){
  adj.nodes.ul=ends(ig.ul,edges.ul[e.index])
  xindex=which(otu.ids.ul==adj.nodes.ul[1])
  yindex=which(otu.ids.ul==adj.nodes.ul[2])
  beta=betaMat.ul[xindex,yindex]
  if(beta>0){
    edge.colors.ul=append(edge.colors.ul,"forestgreen")
  }else if(beta<0){
    edge.colors.ul=append(edge.colors.ul,"red")
  }
}

edge.values.ul<-c()
for(e.index in 1:length(edges.ul)){
  adj.nodes.ul=ends(ig.ul,edges.ul[e.index])
  xindex=which(otu.ids.ul==adj.nodes.ul[1])
  yindex=which(otu.ids.ul==adj.nodes.ul[2])
  beta=betaMat.ul[xindex,yindex]
  if(beta>0){
    edge.values.ul=append(edge.values.ul,beta)
  }else if(beta<0){
    edge.values.ul=append(edge.values.ul,beta)
  }
}

E(ig.ul)$color=edge.colors.ul
ig.ul1 <- ig.ul
V(ig.ul1)$name <- as.character(rownames(taxa.ps.ul.nt.50))
V(ig.ul1)$Phylum <- as.character(taxa.ps.ul.nt.50$Phylum)
V(ig.ul1)$Class <- as.character(taxa.ps.ul.nt.50$Class)
V(ig.ul1)$Order <- as.character(taxa.ps.ul.nt.50$Order)
V(ig.ul1)$Famuly <- as.character(taxa.ps.ul.nt.50$Famuly)
V(ig.ul1)$Genus <- as.character(taxa.ps.ul.nt.50$Genus)
E(ig.ul1)$arrow.size<-5
E(ig.ul1)$width <- abs(edge.values.ul) * 10
V(ig.ul1)$color <- phylum.colors[V(ig.ul1)$Phylum]
V(ig.ul1)$size <- size.ul.50 * 5
V(ig.ul1)$frame.color="black"
V(ig.ul1)$label <- taxa_names(ps.ul.nt)
vertex_attr(ig.ul1)
edge_attr(ig.ul1)

l.ul<-layout_with_fr(ig.ul1)

pdf(file="ul.network.50.l.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.ul1,layout=l.ul,vertex.label.font=2,
     vertex.label.famuly="sans",vertex.label.color="black")
dev.off()

pdf(file="ul.network.50.nl.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.ul1,layout=l.ul,vertex.label.font=2,vertex.label=NA,
     vertex.label.famuly="sans",vertex.label.color="black")
dev.off()

# PL network #
se.pl <- spiec.easi(ps.nt.pl.50, method='mb',lambda.min.ratio=2e-3,
                    nlambda=100, icov.select.params=list(rep.num=30))
ig.pl <- adj2igraph(se.pl$refit,vertex.attr=list(name=taxa_names(ps.nt.pl.50)))

print(se.pl$variabplity)
betaMat.pl <- as.matrix(symBeta(getOptBeta(se.pl)))
positive.pl <- length(betaMat.pl[betaMat.pl>0])/2 
negative.pl <- length(betaMat.pl[betaMat.pl<0])/2 
total.pl <- length(betaMat.pl[betaMat.pl!=0])/2 
total.pl.prop<-total.pl/(length(betaMat.pl)/2)

otu.ids.pl<-colnames(se.pl$data)
edges.pl<-E(ig.pl)
edge.colors.pl<-c()
for(e.index in 1:length(edges.pl)){
  adj.nodes.pl=ends(ig.pl,edges.pl[e.index])
  xindex=which(otu.ids.pl==adj.nodes.pl[1])
  yindex=which(otu.ids.pl==adj.nodes.pl[2])
  beta=betaMat.pl[xindex,yindex]
  if(beta>0){
    edge.colors.pl=append(edge.colors.pl,"forestgreen")
  }else if(beta<0){
    edge.colors.pl=append(edge.colors.pl,"red")
  }
}

edge.values.pl<-c()
for(e.index in 1:length(edges.pl)){
  adj.nodes.pl=ends(ig.pl,edges.pl[e.index])
  xindex=which(otu.ids.pl==adj.nodes.pl[1])
  yindex=which(otu.ids.pl==adj.nodes.pl[2])
  beta=betaMat.pl[xindex,yindex]
  if(beta>0){
    edge.values.pl=append(edge.values.pl,beta)
  }else if(beta<0){
    edge.values.pl=append(edge.values.pl,beta)
  }
}

E(ig.pl)$color=edge.colors.pl
ig.pl1 <- ig.pl
V(ig.pl1)$name <- as.character(rownames(taxa.ps.pl.nt.50))
V(ig.pl1)$Phylum <- as.character(taxa.ps.pl.nt.50$Phylum)
V(ig.pl1)$Class <- as.character(taxa.ps.pl.nt.50$Class)
V(ig.pl1)$Order <- as.character(taxa.ps.pl.nt.50$Order)
V(ig.pl1)$Famply <- as.character(taxa.ps.pl.nt.50$Family)
V(ig.pl1)$Genus <- as.character(taxa.ps.pl.nt.50$Genus)
E(ig.pl1)$arrow.size<-5
E(ig.pl1)$width <- abs(edge.values.pl) * 25
V(ig.pl1)$color <- phylum.colors[V(ig.pl1)$Phylum]
V(ig.pl1)$size <- size.pl.50 * 5
V(ig.pl1)$frame.color="black"
V(ig.pl1)$label <- taxa_names(ps.pl.nt)
vertex_attr(ig.pl1)
edge_attr(ig.pl1)

l.pl<-layout_with_fr(ig.pl1)

pdf(file="pl.network.se.50.l.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.pl1,layout=l.pl,vertex.label.font=2,
     vertex.label.famply="sans",vertex.label.color="black")
dev.off()

pdf(file="pl.network.se.50.nl.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.pl1,layout=l.pl,vertex.label.font=2,vertex.label=NA,
     vertex.label.famply="sans",vertex.label.color="black")
dev.off()

#### 4.1.2) SparCC ####
##### Root ####
#### IIR #
sp.iir<-sparcc(as.matrix(otu_table(ps.nt.iir.50)))
sp.iir.boot<-sparccboot(as.matrix(otu_table(ps.iir.nt)),R=100)
sp.iir.pval<-pval.sparccboot(sp.iir.boot)

# generate a symmetrix matrix of P values #
mat<-matrix(0,ntaxa(ps.iir.nt),ntaxa(ps.iir.nt))
mat[lower.tri(mat, diag=FALSE)] <- sp.iir.pval$pvals
mat <- mat + t(mat)
diag(mat) <- 1

# create correlation matrix thresholded at r >= 0.6, P < 0.05 #
sparcc.graph.iir <- sp.iir$Cor
diag(sparcc.graph.iir) <- 0
sparcc.graph.iir <- abs(sp.iir$Cor) >= 0.6 & mat <= 0.05
sparcc.graph.iir <- Matrix(sparcc.graph.iir, sparse=TRUE)
ig.sparcc.iir <- adj2igraph(sparcc.graph.iir)
ig.sparcc.iir<-readRDS("Network/iGraph/ig.sparcc.iir.50.rds")

otu.ids.iir<-rownames(taxa.ps.iir.nt.50)
V(ig.sparcc.iir)$label <- rownames(taxa.ps.iir.nt.50)
edges.iir<-E(ig.sparcc.iir)
iir.cor <- sp.iir$Cor
rownames(iir.cor)<-rownames(taxa.ps.iir.nt.50)
colnames(iir.cor)<-rownames(taxa.ps.iir.nt.50)

edge.colors.iir<-c()
for(e.index in 1:length(edges.iir)){
  adj.nodes.iir=ends(ig.sparcc.iir,edges.iir[e.index])
  xindex=adj.nodes.iir[1]
  yindex=adj.nodes.iir[2]
  vcor=iir.cor[xindex,yindex]
  if(vcor>0){
    edge.colors.iir=append(edge.colors.iir,"forestgreen")
  }else if(vcor<0){
    edge.colors.iir=append(edge.colors.iir,"red")
  }
}

edge.values.iir<-c()
for(e.index in 1:length(edges.iir)){
  adj.nodes.iir=ends(ig.sparcc.iir,edges.iir[e.index])
  xindex=adj.nodes.iir[1]
  yindex=adj.nodes.iir[2]
  vcor=iir.cor[xindex,yindex]
  if(vcor>0){
    edge.values.iir=append(edge.values.iir,vcor)
  }else if(vcor<0){
    edge.values.iir=append(edge.values.iir,vcor)
  }
}

E(ig.sparcc.iir)$color=edge.colors.iir
V(ig.sparcc.iir)$name <- as.character(rownames(taxa.ps.iir.nt.50))
V(ig.sparcc.iir)$Phylum <- as.character(taxa.ps.iir.nt.50$Phylum)
V(ig.sparcc.iir)$Class <- as.character(taxa.ps.iir.nt.50$Class)
V(ig.sparcc.iir)$Order <- as.character(taxa.ps.iir.nt.50$Order)
V(ig.sparcc.iir)$Family <- as.character(taxa.ps.iir.nt.50$Family)
V(ig.sparcc.iir)$Genus <- as.character(taxa.ps.iir.nt.50$Genus)
E(ig.sparcc.iir)$arrow.size<-5
E(ig.sparcc.iir)$width <- abs(edge.values.iir) * 2
V(ig.sparcc.iir)$color <- phylum.colors[V(ig.sparcc.iir)$Phylum]
V(ig.sparcc.iir)$size <- size.iir.50 * 5
V(ig.sparcc.iir)$frame.color="black"
vertex_attr(ig.sparcc.iir)
edge_attr(ig.sparcc.iir)


pdf(file="iir.network.sparcc.50.nl.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.sparcc.iir,layout=l.iir,vertex.label.font=2,vertex.label=NA,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

pdf(file="iir.network.sparcc.50.l.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.sparcc.iir,layout=l.iir,vertex.label.font=2,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

##### IUR Network #
sp.iur<-sparcc(as.matrix(otu_table(ps.nt.ir.50)))
sp.iur.boot<-sparccboot(as.matrix(otu_table(ps.iur.nt)),R=100)
sp.iur.pval<-pval.sparccboot(sp.iur.boot)

# generate a symmetrix matrix of P values #
mat<-matrix(0,ntaxa(ps.iur.nt),ntaxa(ps.iur.nt))
mat[lower.tri(mat, diag=FALSE)] <- sp.iur.pval$pvals
mat <- mat + t(mat)
diag(mat) <- 1

# create correlation matrix thresholded at r >= 0.6, P < 0.05 #
sparcc.graph.iur <- sp.iur$Cor
diag(sparcc.graph.iur) <- 0
sparcc.graph.iur <- abs(sp.irr$Cor) >= 0.6 & mat <= 0.05
sparcc.graph.iur <- Matrix(sparcc.graph.iur, sparse=TRUE)
ig.sparcc.iur <- adj2igraph(sparcc.graph.iur)
ig.sparcc.iur<-readRDS("Network/iGraph/ig.sparcc.ir.50.rds")

otu.ids.iur<-rownames(taxa.ps.ir.nt.50)
V(ig.sparcc.iur)$label <- rownames(taxa.ps.ir.nt.50)
edges.iur<-E(ig.sparcc.iur)
iur.cor <- sp.iur$Cor
rownames(iur.cor)<-rownames(taxa.ps.ir.nt.50)
colnames(iur.cor)<-rownames(taxa.ps.ir.nt.50)

edge.colors.iur<-c()
for(e.index in 1:length(edges.iur)){
  adj.nodes.iur=ends(ig.sparcc.iur,edges.iur[e.index])
  xindex=adj.nodes.iur[1]
  yindex=adj.nodes.iur[2]
  vcor=iur.cor[xindex,yindex]
  if(vcor>0){
    edge.colors.iur=append(edge.colors.iur,"forestgreen")
  }else if(vcor<0){
    edge.colors.iur=append(edge.colors.iur,"red")
  }
}

edge.values.iur<-c()
for(e.index in 1:length(edges.iur)){
  adj.nodes.iur=ends(ig.sparcc.iur,edges.iur[e.index])
  xindex=adj.nodes.iur[1]
  yindex=adj.nodes.iur[2]
  vcor=iur.cor[xindex,yindex]
  if(vcor>0){
    edge.values.iur=append(edge.values.iur,vcor)
  }else if(vcor<0){
    edge.values.iur=append(edge.values.iur,vcor)
  }
}

E(ig.sparcc.iur)$color=edge.colors.iur
V(ig.sparcc.iur)$name <- as.character(rownames(taxa.ps.ir.nt.50))
V(ig.sparcc.iur)$Phylum <- as.character(taxa.ps.ir.nt.50$Phylum)
V(ig.sparcc.iur)$Class <- as.character(taxa.ps.ir.nt.50$Class)
V(ig.sparcc.iur)$Order <- as.character(taxa.ps.ir.nt.50$Order)
V(ig.sparcc.iur)$Family <- as.character(taxa.ps.ir.nt.50$Family)
V(ig.sparcc.iur)$Genus <- as.character(taxa.ps.ir.nt.50$Genus)
E(ig.sparcc.iur)$arrow.size<-5
E(ig.sparcc.iur)$width <- abs(edge.values.iur) * 2
V(ig.sparcc.iur)$color <- phylum.colors[V(ig.sparcc.iur)$Phylum]
V(ig.sparcc.iur)$size <- size.ir.50 * 5
V(ig.sparcc.iur)$frame.color="black"
vertex_attr(ig.sparcc.iur)
edge_attr(ig.sparcc.iur)

pdf(file="iur.network.sparcc.50.l.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.sparcc.iur,layout=l.ir,vertex.label.font=0.2,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

pdf(file="iur.network.sparcc.50.nl.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.sparcc.iur,layout=l.ir,vertex.label.font=0.2,vertex.label=NA,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

##### UR Network #
sp.ur<-sparcc(as.matrix(otu_table(ps.nt.ur.50)))
sp.ur.boot<-sparccboot(as.matrix(otu_table(ps.ur.nt)),R=100)
sp.ur.pval<-pval.sparccboot(sp.ur.boot)

# generate a symmetrix matrix of P values #
mat<-matrix(0,ntaxa(ps.ur.nt),ntaxa(ps.ur.nt))
mat[lower.tri(mat, diag=FALSE)] <- sp.ur.pval$pvals
mat <- mat + t(mat)
diag(mat) <- 1

# create correlation matrix thresholded at r >= 0.6, P < 0.05 #
sparcc.graph.ur <- sp.ur$Cor
diag(sparcc.graph.ur) <- 0
sparcc.graph.ur <- abs(sp.irr$Cor) >= 0.6 & mat <= 0.05
sparcc.graph.ur <- Matrix(sparcc.graph.ur, sparse=TRUE)
ig.sparcc.ur <- adj2igraph(sparcc.graph.ur)
ig.sparcc.ur<-readRDS("Network/iGraph/ig.sparcc.ur.50.rds")

otu.ids.ur<-rownames(taxa.ps.ur.nt.50)
V(ig.sparcc.ur)$label <- rownames(taxa.ps.ur.nt.50)
edges.ur<-E(ig.sparcc.ur)
ur.cor <- sp.ur$Cor
rownames(ur.cor)<-rownames(taxa.ps.ur.nt.50)
colnames(ur.cor)<-rownames(taxa.ps.ur.nt.50)

edge.colors.ur<-c()
for(e.index in 1:length(edges.ur)){
  adj.nodes.ur=ends(ig.sparcc.ur,edges.ur[e.index])
  xindex=adj.nodes.ur[1]
  yindex=adj.nodes.ur[2]
  vcor=ur.cor[xindex,yindex]
  if(vcor>0){
    edge.colors.ur=append(edge.colors.ur,"forestgreen")
  }else if(vcor<0){
    edge.colors.ur=append(edge.colors.ur,"red")
  }
}

edge.values.ur<-c()
for(e.index in 1:length(edges.ur)){
  adj.nodes.ur=ends(ig.sparcc.ur,edges.ur[e.index])
  xindex=adj.nodes.ur[1]
  yindex=adj.nodes.ur[2]
  vcor=ur.cor[xindex,yindex]
  if(vcor>0){
    edge.values.ur=append(edge.values.ur,vcor)
  }else if(vcor<0){
    edge.values.ur=append(edge.values.ur,vcor)
  }
}

E(ig.sparcc.ur)$color=edge.colors.ur
V(ig.sparcc.ur)$name <- as.character(rownames(taxa.ps.ur.nt.50))
V(ig.sparcc.ur)$Phylum <- as.character(taxa.ps.ur.nt.50$Phylum)
V(ig.sparcc.ur)$Class <- as.character(taxa.ps.ur.nt.50$Class)
V(ig.sparcc.ur)$Order <- as.character(taxa.ps.ur.nt.50$Order)
V(ig.sparcc.ur)$Family <- as.character(taxa.ps.ur.nt.50$Family)
V(ig.sparcc.ur)$Genus <- as.character(taxa.ps.ur.nt.50$Genus)
E(ig.sparcc.ur)$arrow.size<-5
E(ig.sparcc.ur)$width <- abs(edge.values.ur) * 2
V(ig.sparcc.ur)$color <- phylum.colors[V(ig.sparcc.ur)$Phylum]
V(ig.sparcc.ur)$size <- size.ur.50 * 5
V(ig.sparcc.ur)$frame.color="black"
vertex_attr(ig.sparcc.ur)
edge_attr(ig.sparcc.ur)

pdf(file="ur.network.sparcc.50.l.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.sparcc.ur,layout=l.ur,vertex.label.font=0.2,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

pdf(file="ur.network.sparcc.50.nl.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.sparcc.ur,layout=l.ur,vertex.label.font=0.2,vertex.label=NA,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

##### PR network #
sp.pr<-sparcc(as.matrix(otu_table(ps.nt.pr.50)))
sp.pr.boot<-sparccboot(as.matrix(otu_table(ps.pr.nt)),R=100)
sp.pr.pval<-pval.sparccboot(sp.pr.boot)

# generate a symmetrix matrix of P values #
mat<-matrix(0,ntaxa(ps.pr.nt),ntaxa(ps.pr.nt))
mat[lower.tri(mat, diag=FALSE)] <- sp.pr.pval$pvals
mat <- mat + t(mat)
diag(mat) <- 1

# create correlation matrix thresholded at r >= 0.6, P < 0.05 #
sparcc.graph.pr <- sp.pr$Cor
diag(sparcc.graph.pr) <- 0
sparcc.graph.pr <- abs(sp.irr$Cor) >= 0.6 & mat <= 0.05
sparcc.graph.pr <- Matrix(sparcc.graph.pr, sparse=TRUE)
ig.sparcc.pr <- adj2igraph(sparcc.graph.pr)
ig.sparcc.pr<-readRDS("Network/iGraph/ig.sparcc.pr.50.rds")

otu.ids.pr<-taxa_names(ps.nt.pr.50)
V(ig.sparcc.pr)$label <- taxa_names(ps.nt.pr.50)
edges.pr<-E(ig.sparcc.pr)
pr.cor <- sp.pr$Cor
rownames(pr.cor)<-taxa_names(ps.nt.pr.50)
colnames(pr.cor)<-taxa_names(ps.nt.pr.50)

edge.colors.pr<-c()
for(e.index in 1:length(edges.pr)){
  adj.nodes.pr=ends(ig.sparcc.pr,edges.pr[e.index])
  xindex=adj.nodes.pr[1]
  yindex=adj.nodes.pr[2]
  vcor=pr.cor[xindex,yindex]
  if(vcor>0){
    edge.colors.pr=append(edge.colors.pr,"forestgreen")
  }else if(vcor<0){
    edge.colors.pr=append(edge.colors.pr,"red")
  }
}

edge.values.pr<-c()
for(e.index in 1:length(edges.pr)){
  adj.nodes.pr=ends(ig.sparcc.pr,edges.pr[e.index])
  xindex=adj.nodes.pr[1]
  yindex=adj.nodes.pr[2]
  vcor=pr.cor[xindex,yindex]
  if(vcor>0){
    edge.values.pr=append(edge.values.pr,vcor)
  }else if(vcor<0){
    edge.values.pr=append(edge.values.pr,vcor)
  }
}

E(ig.sparcc.pr)$color=edge.colors.pr
V(ig.sparcc.pr)$name <- as.character(rownames(taxa.ps.pr.nt.50))
V(ig.sparcc.pr)$Phylum <- as.character(taxa.ps.pr.nt.50$Phylum)
V(ig.sparcc.pr)$Class <- as.character(taxa.ps.pr.nt.50$Class)
V(ig.sparcc.pr)$Order <- as.character(taxa.ps.pr.nt.50$Order)
V(ig.sparcc.pr)$Family <- as.character(taxa.ps.pr.nt.50$Family)
V(ig.sparcc.pr)$Genus <- as.character(taxa.ps.pr.nt.50$Genus)
E(ig.sparcc.pr)$arrow.size<-5
E(ig.sparcc.pr)$width <- abs(edge.values.pr) * 2
V(ig.sparcc.pr)$color <- phylum.colors[V(ig.sparcc.pr)$Phylum]
V(ig.sparcc.pr)$size <- size.pr.50 * 5
V(ig.sparcc.pr)$frame.color="black"
vertex_attr(ig.sparcc.pr)
edge_attr(ig.sparcc.pr)

pdf(file="pr.network.sparcc.50.l.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.sparcc.pr,layout=l.pr,vertex.label.font=0.2,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

pdf(file="pr.network.sparcc.50.nl.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.sparcc.pr,layout=l.pr,vertex.label.font=0.2,vertex.label=NA,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

#### Leaf ####
### IL #
sp.il<-sparcc(as.matrix(otu_table(ps.nt.il.50)))
sp.il.boot<-sparccboot(as.matrix(otu_table(ps.il.nt)),R=100)
sp.il.pval<-pval.sparccboot(sp.il.boot)

# generate a symmetrix matrix of P values #
mat<-matrix(0,ntaxa(ps.il.nt),ntaxa(ps.il.nt))
mat[lower.tri(mat, diag=FALSE)] <- sp.il.pval$pvals
mat <- mat + t(mat)
diag(mat) <- 1

# create correlation matrix thresholded at r >= 0.6, P < 0.05 #
sparcc.graph.il <- sp.il$Cor
diag(sparcc.graph.il) <- 0
sparcc.graph.il <- abs(sp.il$Cor) >= 0.6 & mat <= 0.05
sparcc.graph.il <- Matrix(sparcc.graph.il, sparse=TRUE)
ig.sparcc.il <- adj2igraph(sparcc.graph.il)
ig.sparcc.il<-readRDS("Network/iGraph/ig.sparcc.il.50.rds")

otu.ids.il<-taxa_names(ps.nt.il.50)
V(ig.sparcc.il)$label <- taxa_names(ps.nt.il.50)
edges.il<-E(ig.sparcc.il)
il.cor <- sp.il$Cor
rownames(il.cor)<-taxa_names(ps.nt.il.50)
colnames(il.cor)<-taxa_names(ps.nt.il.50)

edge.colors.il<-c()
for(e.index in 1:length(edges.il)){
  adj.nodes.il=ends(ig.sparcc.il,edges.il[e.index])
  xindex=adj.nodes.il[1]
  yindex=adj.nodes.il[2]
  vcor=il.cor[xindex,yindex]
  if(vcor>0){
    edge.colors.il=append(edge.colors.il,"forestgreen")
  }else if(vcor<0){
    edge.colors.il=append(edge.colors.il,"red")
  }
}

edge.values.il<-c()
for(e.index in 1:length(edges.il)){
  adj.nodes.il=ends(ig.sparcc.il,edges.il[e.index])
  xindex=adj.nodes.il[1]
  yindex=adj.nodes.il[2]
  vcor=il.cor[xindex,yindex]
  if(vcor>0){
    edge.values.il=append(edge.values.il,vcor)
  }else if(vcor<0){
    edge.values.il=append(edge.values.il,vcor)
  }
}

E(ig.sparcc.il)$color=edge.colors.il
V(ig.sparcc.il)$name <- as.character(rownames(taxa.ps.il.nt.50))
V(ig.sparcc.il)$Phylum <- as.character(taxa.ps.il.nt.50$Phylum)
V(ig.sparcc.il)$Class <- as.character(taxa.ps.il.nt.50$Class)
V(ig.sparcc.il)$Order <- as.character(taxa.ps.il.nt.50$Order)
V(ig.sparcc.il)$Family <- as.character(taxa.ps.il.nt.50$Family)
V(ig.sparcc.il)$Genus <- as.character(taxa.ps.il.nt.50$Genus)
E(ig.sparcc.il)$arrow.size<-5
E(ig.sparcc.il)$width <- abs(edge.values.il) * 2
V(ig.sparcc.il)$color <- phylum.colors[V(ig.sparcc.il)$Phylum]
V(ig.sparcc.il)$size <- size.il.50 * 5
V(ig.sparcc.il)$frame.color="black"
vertex_attr(ig.sparcc.il)
edge_attr(ig.sparcc.il)

pdf(file="il.network.sparcc.50.l.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.sparcc.il,layout=l.il,vertex.label.font=2,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

### UL #
sp.ul<-sparcc(as.matrix(otu_table(ps.nt.ul.50)))
sp.ul.boot<-sparccboot(as.matrix(otu_table(ps.ul.nt)),R=100)
sp.ul.pval<-pval.sparccboot(sp.ul.boot)

# generate a symmetrix matrix of P values #
mat<-matrix(0,ntaxa(ps.ul.nt),ntaxa(ps.ul.nt))
mat[lower.tri(mat, diag=FALSE)] <- sp.ul.pval$pvals
mat <- mat + t(mat)
diag(mat) <- 1

# create correlation matrix thresholded at r >= 0.6, P < 0.05 #
sparcc.graph.ul <- sp.ul$Cor
diag(sparcc.graph.ul) <- 0
sparcc.graph.ul <- abs(sp.ul$Cor) >= 0.6 & mat <= 0.05
sparcc.graph.ul <- Matrix(sparcc.graph.ul, sparse=TRUE)
ig.sparcc.ul <- adj2igraph(sparcc.graph.ul)
ig.sparcc.ul<-readRDS("Network/iGraph/ig.sparcc.ul.50.rds")

otu.ids.ul<-taxa_names(ps.nt.ul.50)
V(ig.sparcc.ul)$label <- taxa_names(ps.nt.ul.50)
edges.ul<-E(ig.sparcc.ul)
ul.cor <- sp.ul$Cor
rownames(ul.cor)<-taxa_names(ps.nt.ul.50)
colnames(ul.cor)<-taxa_names(ps.nt.ul.50)

edge.colors.ul<-c()
for(e.index in 1:length(edges.ul)){
  adj.nodes.ul=ends(ig.sparcc.ul,edges.ul[e.index])
  xindex=adj.nodes.ul[1]
  yindex=adj.nodes.ul[2]
  vcor=ul.cor[xindex,yindex]
  if(vcor>0){
    edge.colors.ul=append(edge.colors.ul,"forestgreen")
  }else if(vcor<0){
    edge.colors.ul=append(edge.colors.ul,"red")
  }
}

edge.values.ul<-c()
for(e.index in 1:length(edges.ul)){
  adj.nodes.ul=ends(ig.sparcc.ul,edges.ul[e.index])
  xindex=adj.nodes.ul[1]
  yindex=adj.nodes.ul[2]
  vcor=ul.cor[xindex,yindex]
  if(vcor>0){
    edge.values.ul=append(edge.values.ul,vcor)
  }else if(vcor<0){
    edge.values.ul=append(edge.values.ul,vcor)
  }
}

E(ig.sparcc.ul)$color=edge.colors.ul
V(ig.sparcc.ul)$name <- as.character(rownames(taxa.ps.ul.nt.50))
V(ig.sparcc.ul)$Phylum <- as.character(taxa.ps.ul.nt.50$Phylum)
V(ig.sparcc.ul)$Class <- as.character(taxa.ps.ul.nt.50$Class)
V(ig.sparcc.ul)$Order <- as.character(taxa.ps.ul.nt.50$Order)
V(ig.sparcc.ul)$Family <- as.character(taxa.ps.ul.nt.50$Family)
V(ig.sparcc.ul)$Genus <- as.character(taxa.ps.ul.nt.50$Genus)
E(ig.sparcc.ul)$arrow.size<-5
E(ig.sparcc.ul)$width <- abs(edge.values.ul) * 2
V(ig.sparcc.ul)$color <- phylum.colors[V(ig.sparcc.ul)$Phylum]
V(ig.sparcc.ul)$size <- size.ul.50 * 5
V(ig.sparcc.ul)$frame.color="black"
vertex_attr(ig.sparcc.ul)
edge_attr(ig.sparcc.ul)

pdf(file="ul.network.sparcc.50.nl.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.sparcc.ul,layout=l.ul,vertex.label.font=2,vertex.label=NA,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

### PL #
sp.pl<-sparcc(as.matrix(otu_table(ps.nt.pl.50)))
sp.pl.boot<-sparccboot(as.matrix(otu_table(ps.pl.nt)),R=100)
sp.pl.pval<-pval.sparccboot(sp.pl.boot)

# generate a symmetrix matrix of P values #
mat<-matrix(0,ntaxa(ps.pl.nt),ntaxa(ps.pl.nt))
mat[lower.tri(mat, diag=FALSE)] <- sp.pl.pval$pvals
mat <- mat + t(mat)
diag(mat) <- 1

# create correlation matrix thresholded at r >= 0.6, P < 0.05 #
sparcc.graph.pl <- sp.pl$Cor
diag(sparcc.graph.pl) <- 0
sparcc.graph.pl <- abs(sp.pl$Cor) >= 0.6 & mat <= 0.05
sparcc.graph.pl <- Matrix(sparcc.graph.pl, sparse=TRUE)
ig.sparcc.pl <- adj2igraph(sparcc.graph.pl)
ig.sparcc.pl<-readRDS("Network/iGraph/ig.sparcc.pl.50.rds")

otu.ids.pl<-taxa_names(ps.nt.pl.50)
V(ig.sparcc.pl)$label <- taxa_names(ps.nt.pl.50)
edges.pl<-E(ig.sparcc.pl)
pl.cor <- sp.pl$Cor
rownames(pl.cor)<-taxa_names(ps.nt.pl.50)
colnames(pl.cor)<-taxa_names(ps.nt.pl.50)

edge.colors.pl<-c()
for(e.index in 1:length(edges.pl)){
  adj.nodes.pl=ends(ig.sparcc.pl,edges.pl[e.index])
  xindex=adj.nodes.pl[1]
  yindex=adj.nodes.pl[2]
  vcor=pl.cor[xindex,yindex]
  if(vcor>0){
    edge.colors.pl=append(edge.colors.pl,"forestgreen")
  }else if(vcor<0){
    edge.colors.pl=append(edge.colors.pl,"red")
  }
}

edge.values.pl<-c()
for(e.index in 1:length(edges.pl)){
  adj.nodes.pl=ends(ig.sparcc.pl,edges.pl[e.index])
  xindex=adj.nodes.pl[1]
  yindex=adj.nodes.pl[2]
  vcor=pl.cor[xindex,yindex]
  if(vcor>0){
    edge.values.pl=append(edge.values.pl,vcor)
  }else if(vcor<0){
    edge.values.pl=append(edge.values.pl,vcor)
  }
}

E(ig.sparcc.pl)$color=edge.colors.pl
V(ig.sparcc.pl)$name <- as.character(rownames(taxa.ps.pl.nt.50))
V(ig.sparcc.pl)$Phylum <- as.character(taxa.ps.pl.nt.50$Phylum)
V(ig.sparcc.pl)$Class <- as.character(taxa.ps.pl.nt.50$Class)
V(ig.sparcc.pl)$Order <- as.character(taxa.ps.pl.nt.50$Order)
V(ig.sparcc.pl)$Family <- as.character(taxa.ps.pl.nt.50$Family)
V(ig.sparcc.pl)$Genus <- as.character(taxa.ps.pl.nt.50$Genus)
E(ig.sparcc.pl)$arrow.size<-5
E(ig.sparcc.pl)$width <- abs(edge.values.pl) * 2
V(ig.sparcc.pl)$color <- phylum.colors[V(ig.sparcc.pl)$Phylum]
V(ig.sparcc.pl)$size <- size.pl.50 * 5
V(ig.sparcc.pl)$frame.color="black"
vertex_attr(ig.sparcc.pl)
edge_attr(ig.sparcc.pl)

pdf(file="pl.network.sparcc.50.l.pdf",
    width=15,height=15,bg = "transparent")
plot(ig.sparcc.pl,layout=l.pl,vertex.label.font=2,
     vertex.label.family="sans",vertex.label.color="black")
dev.off()

#### 4.2) Network attributes #####
# Degree distributions #
dd.iir <- degree.distribution(ig.iir)
dd.pr <- degree.distribution(ig.pr)
dd.ur <- degree.distribution(ig.ur)
dd.iur <- degree.distribution(ig.iur)

dd.iir<-c(dd.iir,rep("0",c(length(dd.iur)-length(dd.iir))))
dd.pr<-c(dd.pr,rep("0",c(length(dd.iur)-length(dd.pr))))
dd.ur<-c(dd.ur,rep("0",c(length(dd.iur)-length(dd.ur))))
degree.num<-as.factor(c(1:13))

dd.iur<-cbind(dd.iur,degree.num)
dd.pr<-cbind(dd.pr,degree.num)
dd.iir<-cbind(dd.iir,degree.num)
dd.ur<-cbind(dd.ur,degree.num)
dd<-as.data.frame(rbind(dd.iur,dd.pr,dd.iir,dd.ur))
dd$community<-c(rep("I.R",13),rep("P.R",13),rep("I.IR",13),rep("U.R",13))
colnames(dd)[1]<-c("freq")
dd$degree.num<-factor(dd$degree.num, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13"))
dd$freq<-as.numeric(as.character(dd$freq))

pdf(file="degree_distributions.se.50.pdf",
    width=10,height=8,bg = "transparent")
ggplot(dd, aes(x=degree.num, y=freq, group=community)) +
  geom_line(aes(color=community),size=2)+
  geom_point(aes(color=community),size=3)+
  scale_y_continuous(breaks=seq(0,0.5,0.1)) +
  scale_color_manual(values=sample.colors.beta) +
  theme.gg.leg
dev.off()

dd.iir <- degree.distribution(ig.sparcc.iir)
dd.pr <- degree.distribution(ig.sparcc.pr)
dd.ur <- degree.distribution(ig.sparcc.ur)
dd.iur <- degree.distribution(ig.sparcc.iur)

dd.iir<-c(dd.iir,rep("0",c(length(dd.ur)-length(dd.iir))))
dd.pr<-c(dd.pr,rep("0",c(length(dd.ur)-length(dd.pr))))
dd.iur<-c(dd.iur,rep("0",c(length(dd.ur)-length(dd.iur))))
degree.num<-as.factor(c(1:24))

dd.iur<-cbind(dd.iur,degree.num)
dd.pr<-cbind(dd.pr,degree.num)
dd.iir<-cbind(dd.iir,degree.num)
dd.ur<-cbind(dd.ur,degree.num)
dd<-as.data.frame(rbind(dd.iur,dd.pr,dd.iir,dd.ur))
dd$community<-c(rep("I.R",24),rep("P.R",24),rep("I.IR",24),rep("U.R",24))
colnames(dd)[1]<-c("freq")
dd$degree.num<-factor(dd$degree.num, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13",
                                              "14","15","16","17","18","19","20","21","22","23","24"))
dd$freq<-as.numeric(as.character(dd$freq))

pdf(file="degree_distributions.sparcc.50.pdf",
    width=10,height=8,bg = "transparent")
ggplot(dd, aes(x=degree.num, y=freq, group=community)) +
  geom_line(aes(color=community),size=2)+
  geom_point(aes(color=community),size=3)+
  scale_y_continuous(breaks=seq(0,0.5,0.1)) +
  scale_color_manual(values=sample.colors.beta) +
  theme.gg.leg
dev.off()

# Density #
edge_density(ig.iir, loops=F) 
edge_density(ig.iur, loops=F) 
edge_density(ig.ur, loops=F) 
edge_density(ig.pr, loops=F) 
edge_density(ig.il, loops=F) 
edge_density(ig.ul, loops=F) 
edge_density(ig.pl, loops=F) 

edge_density(ig.sparcc.iir, loops=F) 
edge_density(ig.sparcc.iur, loops=F) 
edge_density(ig.sparcc.ur, loops=F)
edge_density(ig.sparcc.pr, loops=F) 
edge_density(ig.sparcc.il, loops=F) 
edge_density(ig.sparcc.ul, loops=F) 
edge_density(ig.sparcc.pl, loops=F) 


centr_betw(ig.iir, directed = FALSE,normalized = T)$centralization
centr_betw(ig.iur, directed = FALSE,normalized = T)$centralization
centr_betw(ig.ur, directed = FALSE,normalized = T)$centralization
centr_betw(ig.pr, directed = FALSE,normalized = T)$centralization
centr_betw(ig.il, directed = FALSE,normalized = T)$centralization
centr_betw(ig.ul, directed = FALSE,normalized = T)$centralization
centr_betw(ig.pl, directed = FALSE,normalized = T)$centralization

centr_betw(ig.sparcc.iir, directed = FALSE,normalized = T)$centralization
centr_betw(ig.sparcc.iur, directed = FALSE,normalized = T)$centralization
centr_betw(ig.sparcc.ur, directed = FALSE,normalized = T)$centralization
centr_betw(ig.sparcc.pr, directed = FALSE,normalized = T)$centralization
centr_betw(ig.sparcc.il, directed = FALSE,normalized = T)$centralization
centr_betw(ig.sparcc.ul, directed = FALSE,normalized = T)$centralization
centr_betw(ig.sparcc.pl, directed = FALSE,normalized = T)$centralization


# Centrality indices
degree.iir<-degree(ig.iir, mode="in")
degree.iur<-degree(ig.iur, mode="in")
degree.ur<-degree(ig.ur, mode="in")
degree.pr<-degree(ig.pr, mode="in")

betweeness.iir<-betweenness(ig.iir, directed=F, weights=NA, normalized=T)
betweeness.iur<-betweenness(ig.iur, directed=F, weights=NA, normalized=T)
betweeness.ur<-betweenness(ig.ur, directed=F, weights=NA, normalized=T)
betweeness.pr<-betweenness(ig.pr, directed=F, weights=NA, normalized=T)

sample.bc<-function(x){
  df<- sample(x,50,replace=T)
  mean.bc<-mean(df)
}

bc.iir.samp <- as.data.frame(replicate(10000,sample.bc(betweeness.iir)))
bc.iir.samp$community<-rep("I.IR",10000)
colnames(bc.iir.samp)[1]<-"bc"
bc.iur.samp <- as.data.frame(replicate(10000,sample.bc(betweeness.iur)))
bc.iur.samp$community<-rep("I.R",10000)
colnames(bc.iur.samp)[1]<-"bc"
bc.ur.samp <- as.data.frame(replicate(10000,sample.bc(betweeness.ur)))
bc.ur.samp$community<-rep("U.R",10000)
colnames(bc.ur.samp)[1]<-"bc"
bc.pr.samp <- as.data.frame(replicate(10000,sample.bc(betweeness.pr)))
bc.pr.samp$community<-rep("P.R",10000)
colnames(bc.pr.samp)[1]<-"bc"

sampl.df<-rbind(bc.iir.samp,bc.iur.samp,bc.ur.samp,bc.pr.samp)

pdf(file="centrality.resamp.se.50.pdf",
    width=10,height=8,bg = "transparent")
ggplot(sampl.df,aes(x=bc, fill=community)) + 
  geom_density(alpha=0.5) +
  scale_fill_manual(values=sample.colors.beta) +
  scale_x_continuous(breaks=seq(0,0.05,0.002)) +
  theme.gg.leg
dev.off()

ggplot(sampl.df[sampl.df$community!="P.R",],aes(x=bc, fill=community)) + 
  geom_density(alpha=0.5) +
  scale_fill_manual(values=sample.colors.beta) +
  theme.gg.leg


betweeness.iir<-betweenness(ig.sparcc.iir, directed=F, weights=NA, normalized=T)
betweeness.iur<-betweenness(ig.sparcc.iur, directed=F, weights=NA, normalized=T)
betweeness.ur<-betweenness(ig.sparcc.ur, directed=F, weights=NA, normalized=T)
betweeness.pr<-betweenness(ig.sparcc.pr, directed=F, weights=NA, normalized=T)

bc.iir.samp <- as.data.frame(replicate(10000,sample.bc(betweeness.iir)))
bc.iir.samp$community<-rep("I.IR",10000)
colnames(bc.iir.samp)[1]<-"bc"
bc.iur.samp <- as.data.frame(replicate(10000,sample.bc(betweeness.iur)))
bc.iur.samp$community<-rep("I.R",10000)
colnames(bc.iur.samp)[1]<-"bc"
bc.ur.samp <- as.data.frame(replicate(10000,sample.bc(betweeness.ur)))
bc.ur.samp$community<-rep("U.R",10000)
colnames(bc.ur.samp)[1]<-"bc"
bc.pr.samp <- as.data.frame(replicate(10000,sample.bc(betweeness.pr)))
bc.pr.samp$community<-rep("P.R",10000)
colnames(bc.pr.samp)[1]<-"bc"

sampl.df<-rbind(bc.iir.samp,bc.iur.samp,bc.ur.samp,bc.pr.samp)

pdf(file="centrality.resamp.sparcc.50.pdf",
    width=10,height=8,bg = "transparent")
ggplot(sampl.df,aes(x=bc, fill=community)) + 
  geom_density(alpha=0.5) +
  scale_fill_manual(values=sample.colors.beta) +
  scale_y_continuous(breaks=seq(0,5000,1000)) +
  scale_x_continuous(breaks=seq(0,0.02,0.002)) +
  theme.gg.leg
dev.off()

##### DADA2 Pipeline ##### 
# completed using RStudio AMI courtesy of Louis Aslett, thank you! #
# used a M4.2XL instance #

# courtesy of Benji Callahan - DADA2 Tutorial #
### Filtering and Trimming ###
fastqs <- fns[grepl(".fastq$", fns)] # List names of all files ending with .fastq
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files
# Get sample names from the first part of the forward read filenames
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Fully specify the path for the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)
saveRDS(out, "out_Jun13.rds")

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errF, "errF_Jun13.rds")
saveRDS(errR, "errR_Jun13.rds")

plotErrors(errF, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
saveRDS(dadaFs, "dadaF_Jun13.rds")

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
saveRDS(dadaRs, "dadaR_Jun13.rds")

#merge#
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
saveRDS(mergers, "mergers_Jun13.rds")

seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "seqtab_Jun13.rds")

dim(seqtab) # 103 22507
table(nchar(colnames(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)
table(nchar(colnames(seqtab.nochim))) # frequency of read length
saveRDS(seqtab.nochim, "seqtab.nochim_Jun13.rds")

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
saveRDS(track, "track_Jun13.rds")

taxa <- assignTaxonomy(seqtab.nochim, "rdp_train_set_16.fa.gz", multithread=TRUE)
saveRDS(taxa, "taxa_Jun13.rds")
taxaSilva <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=TRUE)
saveRDS(taxaSilva, "taxaSilva_Jun13.rds")

##### Build bacterial phylogeny #####
seqs<-getSequences(seqtab.total)
seq.names<-colnames(seqtab.total)
write.fasta(as.list(seqs), names = seq.names, file.out="seqtab.total.fasta")

#### Rarefaction ####
set.seed(10021987)
ps3rare = rarefy_even_depth(ps3, sample.size = 1000)
sum(taxa_sums(ps3rare))/sum(taxa_sums(ps3)) # 0.014 
ntaxa(ps3rare)
# rarefying yields about 1.5% of total reads and 8291 ASVs
# from this point repeat analyses above with rarefied data #

## Alpha diversity analyses ##
alpha.1.rare<-as.data.frame(estimate_richness(ps3rare,measures=c("Observed","Shannon","InvSimpson")))
alpha.1.rare.sd<-as(sample_data(ps3rare),'data.frame')
alpha.2.rare<-cbind(alpha.1.rare.sd,alpha.1.rare)
alpha.2.rare$E2<-alpha.2.rare$InvSimpson/alpha.2.rare$Observed
samp<-as.data.frame(otu_table(ps3rare))
rare.tree<-phy_tree(ps3rare)
pd.total.rare<-pd(samp,rare.tree)
alpha.2.rare$pd<-pd.total.rare$PD

# are diversity estimates from rarefied versus non-rarefied data correlated?
alpha.2.rare$individual<-rownames(alpha.2.rare)
a2$individual<-rownames(a2)
alpha.test<-merge(a2,alpha.2.rare,by="individual")

cor.test(alpha.test$Observed.x,alpha.test$Observed.y,data=alpha.test) # r = 0.95
cor.test(alpha.test$InvSimpson.x,alpha.test$InvSimpson.y,data=alpha.test) # r = 0.99
cor.test(alpha.test$E2.x,alpha.test$E2.y,data=alpha.test) # r = 0.79
cor.test(alpha.test$pd.x,alpha.test$pd.y,data=alpha.test) # r = 0.97