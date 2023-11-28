
library("qiime2R") # devtools::install_github("jbisanz/qiime2R")

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


colours = c("#a91919","#582630", "#e1e919", "#4DAA57",
                     "#D55E00","#0072B2", "#F9ECCC", "#679289", "#33658A",
                     "#F6AE2D","#00AFBB", "grey","#998540" ,"black", "purple")
                     
mock <- readRDS("R_objects/mock.rds")

library("phyloseq")
sample_variables(mock)
rank_names(mock)

mock <- subset_samples(mock, Name != "mock")
mock <- subset_samples(mock, Name != "mock2")

predefined_species = c("Botrytis_porri",
             "Alternaria_alternata",
             "Botrytis_cinerea",
             "Peronospora_variabilis",
             "Epicoccum_nigrum",
             "Stemphylium_vesicarium",
             "Peronospora_destructor",
             "Sclerotinia_sclerotiorum",
             "Alternaria_brassicicola",
             "Botrytis_squamosa",
             "Fusarium_oxysporum",
             "Hyaloperonospora_brassicae")

#################################### 
mock_dat <- psmelt(mock)

library(tidyverse)
library(ggtext)

library("dplyr")
mock_dat <- mock_dat %>%
  mutate(Species = case_when(
    Species %in% 
      predefined_species ~ Species,  # Keep valid species as is
    TRUE ~ "other"  # Replace other species with "other"
  ))

otu_rel_abund <- mock_dat %>%
  group_by(Name, Description) %>%
  mutate(rel_abund = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  select(-Abundance) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
               names_to="level",
               values_to="taxon") %>% 
  mutate(Description = factor(Description, 
                              levels=c("Theoretical", "Minimap2", "Bracken", "Kraken", "Qiime2")))

str(otu_rel_abund)

###################################

otu_rel_abund %>%
  filter(level=="Species") %>%
  group_by(Name, Description, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Name, Description, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  mutate(taxon=str_replace(taxon, "_", " ")) %>%
  mutate(taxon = factor(taxon, 
                              levels=c("Botrytis cinerea",
                                       "Botrytis porri",
                                       "Alternaria alternata",
                                       "Peronospora variabilis",
                                       "Epicoccum nigrum",
                                       "Stemphylium vesicarium",
                                       "Peronospora destructor",
                                       "Sclerotinia sclerotiorum",
                                       "Alternaria brassicicola",
                                       "Botrytis squamosa",
                                       "Fusarium oxysporum",
                                       "Hyaloperonospora brassicae",
                                       "Other"))
         ) %>% 
  ggplot(aes(x=Description, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  facet_wrap(vars(Name), nrow = 2)+
  theme(legend.title=element_blank())+
  labs(x=NULL,
       y="Mean Relative Abundance (%)") +
  theme(legend.text = element_text(face="italic"))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=colours)+
  theme(legend.position="right")+
  guides(fill= guide_legend(keywidth = 0.6, 
                            keyheight = 0.7, 
                            ncol=1))+
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0))+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.4))

                     