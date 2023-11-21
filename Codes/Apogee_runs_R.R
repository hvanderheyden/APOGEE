

#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library("phyloseq")
library("ggplot2")     
library("readxl")      
library("dplyr")
library("tibble")
library("vegan")
library("DESeq2") # BiocManager::install("DESeq2")
library("speedyseq") # install with remotes::install_github("mikemc/speedyseq") 
library("ape")
library("ggstar")
library("forcats")
library(patchwork)
library(ggpubr)
library("cowplot")

setwd("/media/herve/10TB/Apogee")

## load the biom_table, taxonomy and metadata ####
biom_spaghetti<- read.csv("/media/herve/10TB/Apogee/3_Fastq/FASTQ_runs/fastq_merged/otu_table.csv", 
                          header=TRUE, 
                          sep=",")
head(biom_spaghetti)

taxo_spaghetti<- read.csv("/media/herve/10TB/Apogee/3_Fastq/FASTQ_runs/fastq_merged/phyloseq_taxonomy.csv", 
                          header=TRUE, 
                          sep=",")
head(taxo_spaghetti)

metaS <- read.csv("/media/herve/10TB/Apogee/7_pipeline_runs/Meta_full.tsv", 
                  header=TRUE, 
                  sep="\t")
head(metaS)

metaS$DOY<-as.factor(metaS$DOY)
metaS$Year<-as.factor(metaS$Year)

metaS <- metaS %>% 
  tibble::column_to_rownames("Sample_ID")

samplesS = sample_data(metaS)



# define the row names from the otu column ####

biom_spaghetti <- biom_spaghetti %>%
  tibble::column_to_rownames("OTU") 

taxo_spaghetti <- taxo_spaghetti %>%
  tibble::column_to_rownames("OTU") 

# Transform into matrixes otu and tax tables (sample table can be left as data frame) ####

biom_spaghetti <- as.matrix(biom_spaghetti)
taxo_spaghetti <- as.matrix(taxo_spaghetti)

class(biom_spaghetti)
class(biom_spaghetti)

# convert to phyloseq objects ####

OTU_spaghetti = otu_table(biom_spaghetti, taxa_are_rows = TRUE)
TAX_spaghetti = phyloseq::tax_table(taxo_spaghetti)


phylo_spaghetti <- phyloseq(OTU_spaghetti, TAX_spaghetti, samplesS)
phylo_spaghetti

sample_variables(phylo_spaghetti)


random_tree = rtree(ntaxa(phylo_spaghetti), rooted=TRUE, tip.label=taxa_names(phylo_spaghetti))
plot(random_tree)

phylo_spaghetti <- phyloseq(OTU_spaghetti, TAX_spaghetti, samplesS, random_tree)
phylo_spaghetti


#subset the buffer vs field samples ####

buffer<- subset_samples(phylo_spaghetti, SampleType =="Buffer")
buffer


plot_bar(buffer, fill="Phylum")

field_samples <- subset_samples(phylo_spaghetti, SampleType =="Field")

field_samples

saveRDS(field_samples,"/media/herve/HERVE_256/Appogee/APOGEE.RDS")

plot_bar(field_samples, fill="Phylum")+
  ylim(0, 10000)


# filter the data to remove low depth samples #####

library(microbiome) # BiocManager::install("microbiome")
library(microbiomeutilities) #remotes::install_github("microsud/microbiomeutilities")

summarize_phyloseq(field_samples)

Dep1<-plot_read_distribution(field_samples, groups = "Year", 
                             plot.type = "histogram")+
  theme_biome_utils()+
  scale_x_continuous(trans='log10', limits = c(100, 500000))+
  scale_fill_manual(values=c("#111111"))+ 
  geom_vline(xintercept = 3000, colour = "black", linetype="dashed")+
  theme(legend.position="none")+
  labs(x = "", y = "Count")

field_samples1 <- prune_samples(sample_sums(field_samples) >= 3000, field_samples)

summarize_phyloseq(field_samples1)

Dep2<-plot_read_distribution(field_samples1, groups = "Year", 
                             plot.type = "histogram")+
  theme_biome_utils()+
  scale_x_continuous(trans='log10', limits = c(100, 500000))+
  scale_fill_manual(values=c("#111111"))+ 
  geom_vline(xintercept = 3000, colour = "black", linetype="dashed")+
  theme(legend.position="none")+
  labs(x = "Reads per samples", y = "Count")


depth<-plot_grid(Dep1+theme(legend.position="none"),
                 Dep2+theme(legend.position="none"), 
                 align="vh",
                 labels = c("A", "B"),
                 hjust = -1,
                 vjust= 2,
                 nrow = 2)

depth_final<-plot_grid(depth, ncol = 1, rel_heights = c(0.8, .05))
depth_final

######### rarefaction 

field_samplesP <- filter_taxa(field_samples1, function(x) sum(x > 0) > (0.1*length(x)), TRUE)

summarize_phyloseq(field_samplesP)

set.seed(1024)
rarecurve  <- ggrarecurve(obj=field_samplesP,
                          factorNames="Year",
                          indexNames=c("Observe", 
                                       "Chao1", 
                                       "ACE")) 

rarecurve +
  theme(legend.spacing.y=unit(0.01,"cm"),
        legend.text=element_text(size=4))+
  xlim(0, 50000)+
  ylim(0, 1000)

field_samplesR <- rarefy_even_depth(field_samplesP, 
                                    rngseed=123,
                                    sample.size=3000,
                                    replace=F)

########################Subset data ######

Y2021 <- subset_samples(field_samplesR, Year =="2021")
Y2021

Y2022 <- subset_samples(field_samplesR, Year =="2022")
Y2022

DIEC <- subset_samples(field_samplesR, Site =="DIEC")
DIEC

rank_names(field_samplesR)

oom <- subset_taxa(field_samplesR, Phylum=="p__Oomycota")
oom

O2021 <- subset_samples(oom, Year =="2021")
O2021

O2022 <- subset_samples(oom, Year =="2022")
O2022


O2021_taxa <- get_taxadf(obj = O2021, taxlevel=7)

graph_O2021_taxa <- ggbartax(obj = O2021_taxa,
                           facetNames=factor("Site"),
                           topn=24)+
  theme(axis.text.x = element_text(angle = 90))
graph_O2021_taxa

O2022_taxa <- get_taxadf(obj = O2022, taxlevel=7)

graph_O2022_taxa <- ggbartax(obj = O2022_taxa,
                             facetNames=factor("Site"),
                             topn=12)+
  theme(axis.text.x = element_text(angle = 90))

graph_O2022_taxa

##### plot taxa ####


plot_bar(field_samplesR, fill="Phylum")

BiocManager::install("MicrobiotaProcess")
library("MicrobiotaProcess")

fY2021_taxa <- get_taxadf(obj = Y2021, taxlevel=7)

graph_fY2021_taxa <- ggbartax(obj = fY2021_taxa,
                                      facetNames=factor("DOY"),
                                      topn=24)+
  theme(axis.text.x = element_text(angle = 90))


graph_fY2021_taxa


fY2022_taxa <- get_taxadf(obj = Y2022, taxlevel=7)

graph_fY2022_taxa <- ggbartax(obj = fY2022_taxa,
                              facetNames=factor("DOY"),
                              topn=24)

graph_fY2022_taxa



graph_fY2022_taxa$data

sample_variables(field_samplesR)

DIEC_taxa <- get_taxadf(obj = DIEC, taxlevel=7)

graph_DIEC_taxa <- ggbartax(obj = DIEC_taxa,
                                      ,facetNames="Year"
                                      ,topn=24)
graph_DIEC_taxa

# distmethod
# "unifrac",  "wunifrac", "manhattan", "euclidean", "canberra", "bray", "kulczynski" ...(vegdist, dis



pcoares <- get_pcoa(obj=Y2022, distmethod="bray", method="hellinger")
# Visulizing the result
pcaplot1 <- ggordpoint(obj=pcoares, biplot=TRUE, speciesannot=FALSE,
                       factorNames=c("Site"), ellipse=TRUE) 
  #scale_color_manual(values=c("#00AED7", "#FD9347")) +
  #scale_fill_manual(values=c("#00AED7", "#FD9347"))
# pc = c(1, 3) to show the first and third principal components.
pcaplot2 <- ggordpoint(obj=pcoares, pc=c(1, 3), biplot=TRUE, speciesannot=TRUE,
                       factorNames=c("Site"), ellipse=TRUE)
  #scale_color_manual(values=c("#00AED7", "#FD9347")) +
  #scale_fill_manual(values=c("#00AED7", "#FD9347"))
pcaplot1 | pcaplot2


#################### 
set.seed(123)
ig <- make_network(Y2022, 
                   dist.fun="unifrac",
                   max.dist=0.5)

plot_network(ig, Y2022, 
             color="DOY", 
             shape="Site", 
             line_weight=0.4, 
             label=NULL)

