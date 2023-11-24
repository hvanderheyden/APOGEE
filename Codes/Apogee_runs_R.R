
    
library("readxl")      

library("tibble")
library("vegan")
library("DESeq2") # BiocManager::install("DESeq2")
library("speedyseq") # install with remotes::install_github("mikemc/speedyseq") 

library("ggstar")
library("forcats")
library("patchwork")
library("ggpubr")
library("cowplot")


## load the biom_table, taxonomy and metadata ####
Apogee_biom<- read.csv("Data/runs/otu_table.csv", 
                          header=TRUE, 
                          sep=",")
head(Apogee_biom)

Apogee_taxo<- read.csv("Data/runs/phyloseq_taxonomy.csv", 
                          header=TRUE, 
                          sep=",")
head(Apogee_taxo)

Apogee_meta <- read.csv("Data/runs/Meta_full.tsv", 
                  header=TRUE, 
                  sep="\t")
head(Apogee_meta)

library("dplyr")
Apogee_meta$DOY<-as.factor(Apogee_meta$DOY)
Apogee_meta$Year<-as.factor(Apogee_meta$Year)

Apogee_meta <- Apogee_meta %>% 
  tibble::column_to_rownames("Sample_ID")

library("phyloseq")
samplesS = sample_data(Apogee_meta)

# define the row names from the otu column ####

Apogee_biom <- Apogee_biom %>%
  tibble::column_to_rownames("OTU") 

Apogee_taxo <- Apogee_taxo %>%
  tibble::column_to_rownames("OTU") 

# Transform into matrixes otu and tax tables (sample table can be left as data frame) ####

Apogee_biom <- as.matrix(Apogee_biom)
Apogee_taxo <- as.matrix(Apogee_taxo)

class(Apogee_biom)
class(Apogee_taxo)

# convert to phyloseq objects ####

Apogee_OTU = otu_table(Apogee_biom, taxa_are_rows = TRUE)
Apogee_TAX = phyloseq::tax_table(Apogee_taxo)


Apogee_PS <- phyloseq(Apogee_OTU, 
                      Apogee_TAX, 
                      samplesS)
Apogee_PS

sample_variables(Apogee_PS)

library("ape")
random_tree = rtree(ntaxa(Apogee_PS), 
                    rooted=TRUE, 
                    tip.label=taxa_names(Apogee_PS))
#plot(random_tree)

Apogee_PS <- phyloseq(Apogee_OTU, 
                            Apogee_TAX, 
                            samplesS, 
                            random_tree)
Apogee_PS

# filter the data to remove low depth samples #####

library(microbiome) # BiocManager::install("microbiome")
library(microbiomeutilities) #remotes::install_github("microsud/microbiomeutilities")

summarize_phyloseq(Apogee_PS)

library("ggplot2") 
Dep1<-plot_read_distribution(Apogee_PS, 
                             groups = "Year", 
                             plot.type = "histogram")+
  theme_biome_utils()+
  scale_x_continuous(trans='log10', 
                     limits = c(100, 300000))+
  scale_fill_manual(values=c("#111111"))+ 
  geom_vline(xintercept = 3000, 
             colour = "black", 
             linetype="dashed")+
  theme(legend.position="none")+
  labs(x = "", y = "Count")

Apogee_PS1 <- prune_samples(sample_sums(Apogee_PS) >= 3000, Apogee_PS)

summarize_phyloseq(Apogee_PS1)

Dep2<-plot_read_distribution(Apogee_PS1, groups = "Year", 
                             plot.type = "histogram")+
  theme_biome_utils()+
  scale_x_continuous(trans='log10', limits = c(100, 300000))+
  scale_fill_manual(values=c("#111111"))+ 
  geom_vline(xintercept = 3000, colour = "black", linetype="dashed")+
  theme(legend.position="none")+
  labs(x = "Reads per samples", y = "Count")

library("cowplot")
depth<-plot_grid(Dep1+theme(legend.position="none"),
                 Dep2+theme(legend.position="none"), 
                 align="vh",
                 labels = c("A", "B"),
                 hjust = -1,
                 vjust= 2,
                 nrow = 2)

depth_final<-plot_grid(depth, ncol = 1, rel_heights = c(0.8, .05))
depth_final

#subset the buffer vs field samples ####

Apogee_buffer<- subset_samples(Apogee_PS1, 
                        SampleType =="Buffer")
Apogee_buffer


plot_bar(Apogee_buffer, fill="Phylum")

Apogee_field <- subset_samples(Apogee_PS1, 
                          SampleType =="Field")

Apogee_field

saveRDS(field_samples,"/media/herve/HERVE_256/Appogee/APOGEE.RDS")

plot_bar(Apogee_field, fill="Phylum")+
  ylim(0, 10000)

######### rarefaction #####

Apogee_fieldP <- filter_taxa(Apogee_field, 
                             function(x) sum(x > 10) > (0.01*length(x)), 
                             TRUE)

summarize_phyloseq(Apogee_fieldP)


library("MicrobiotaProcess")
detach(package:MicrobiotaProcess)


set.seed(1024)
rarecurve  <- ggrarecurve(obj=Apogee_field,
                          factorNames="Year",
                          indexNames=c("Observe", 
                                       "Chao1", 
                                       "ACE")) 

rarecurve +
  theme(legend.spacing.y=unit(0.01,"cm"),
        legend.text=element_text(size=4))+
  xlim(0, 50000)+
  ylim(0, 2000)



field_samplesR <- rarefy_even_depth(Apogee_fieldP, 
                                    rngseed=1024,
                                    sample.size=5000,
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
                                      facetNames="Year",
                                      topn=24)
graph_DIEC_taxa

# distmethod
# "unifrac",  "wunifrac", "manhattan", "euclidean", "canberra", "bray", "kulczynski" ...(vegdist, dis



pcoares <- get_pcoa(obj=field_samplesR, 
                    distmethod="wunifrac", 
                    method="hellinger")
# Visulizing the result
pcaplot1 <- ggordpoint(obj=pcoares, 
                       biplot=FALSE, 
                       speciesannot=FALSE,
                       factorNames=c("DOY"), 
                       ellipse=TRUE) 
  #scale_color_manual(values=c("#00AED7", "#FD9347")) +
  #scale_fill_manual(values=c("#00AED7", "#FD9347"))
# pc = c(1, 3) to show the first and third principal components.
pcaplot2 <- ggordpoint(obj=pcoares, 
                       pc=c(1, 3), 
                       biplot=FALSE, 
                       speciesannot=FALSE,
                       factorNames=c("DOY"), 
                       ellipse=TRUE)
  #scale_color_manual(values=c("#00AED7", "#FD9347")) +
  #scale_fill_manual(values=c("#00AED7", "#FD9347"))
pcaplot1 | pcaplot2


#################### 
set.seed(123)
ig <- make_network(field_samplesR, 
                   dist.fun="bray",
                   max.dist=0.6)

plot_network(ig, field_samplesR, 
             color="DOY", 
             shape="Site", 
             line_weight=0.4, 
             label=NULL)

