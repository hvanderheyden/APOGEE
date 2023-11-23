

library("qiime2R") # devtools::install_github("jbisanz/qiime2R")
library("phyloseq")
library("readxl")      
library("tibble")
library("vegan")
library("DESeq2") # BiocManager::install("DESeq2")
library("speedyseq") # remotes::install_github("mikemc/speedyseq") 
library("ape")
library("ggstar")
library("forcats")
library("patchwork")
library("ggpubr")
library("plotROC")
library("viridis")
library("cowplot")

library("ggtree") # BiocManager::install("ggtree")
library("ggtreeExtra") #install.packages("ggExtra")
library('MicrobiotaProcess') # BiocManager::install("MicrobiotaProcess")
library("tidytree") # install.packages("tidytree")

setwd("C:/Users/vanderheydenh/OneDrive - AGR-AGR/Projets/2023/Biovigilance/APOGEE/Apogee_mock")

mock <- readRDS("C:/Users/vanderheydenh/OneDrive - AGR-AGR/Projets/2023/Biovigilance/APOGEE/Apogee_mock/To_phyloseq1/mock.rds")

sample_variables(mock)

str(mock)


mock <- subset_samples(mock, Name != "mock")
mock <- subset_samples(mock, Name != "mock2")

mockR = rarefy_even_depth(mock, 
                           rngseed=123, 
                           sample.size=0.5*min(sample_sums(mock)), 
                           replace=F)


plot_bar(mockR, fill = "Kingdom")

classtaxa_mock <- get_taxadf(obj=mockR, taxlevel=7)

####### Mock3 #####

mock3 <- subset_samples(mockR, Name =="mock3")

plot_bar(mock3, fill = "Kingdom")

classtaxa3R <- get_taxadf(obj=mock3, taxlevel=7)

test3<- ggbartax(obj=classtaxa3R, 
                 facetNames=c("Description"), 
                 plotgroup=TRUE, 
                 topn=12)+
  xlab(NULL) +
  ylab("Relative abundance (%)") +
  guides(fill= guide_legend(keywidth = 0.5, 
                            keyheight = 0.8, 
                            ncol=5))

test3

test3$data

####### Mock4 #####

mock4 <- subset_samples(mockR, Name =="mock4")

plot_bar(mock4, fill = "Kingdom")


classtaxa4R <- get_taxadf(obj=mock4, taxlevel=7)

test4<- ggbartax(obj=classtaxa4R, 
                 facetNames=c("Description"), 
                 plotgroup=TRUE, 
                 topn=13)+
  xlab(NULL) +
  ylab("Relative abundance (%)") +
  guides(fill= guide_legend(keywidth = 0.5, 
                            keyheight = 0.8, 
                            ncol=4))

test4

test4$data

####### Mock5 #####

mock5 <- subset_samples(mockR, Name =="mock5")

plot_bar(mock5, fill = "Kingdom")

classtaxa5R <- get_taxadf(obj=mock5, taxlevel=7)

test5<- ggbartax(obj=classtaxa5R, 
                 topn=13, 
                 facetNames=c("Description"),
                 plotgroup=TRUE, 
                 )+
  xlab(NULL) +
  ylab("Relative abundance (%)") +
  guides(fill= guide_legend(keywidth = 0.5, 
                            keyheight = 0.8, 
                            ncol=4))

test5

test5$data