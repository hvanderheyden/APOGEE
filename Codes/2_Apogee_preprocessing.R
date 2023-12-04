
library("phyloseq")
library("microbiomeutilities")

###############################################################
########## preprocessing for the mock experiments #############
###############################################################

# Import the raw phyloseq mock_raw object 

mock_raw <- readRDS("R_objects/mock_raw.rds")

library("tidyverse")
taxM <- data.frame(phyloseq::tax_table(mock_raw))
tax.cleanM <- data.frame(row.names = row.names(taxM),
                        Kingdom = str_replace(taxM[,1], "k__",""),
                        Phylum = str_replace(taxM[,2], "p__",""),
                        Class = str_replace(taxM[,3], "c__",""),
                        Order = str_replace(taxM[,4], "o__",""),
                        Family = str_replace(taxM[,5], "f__",""),
                        Genus = str_replace(taxM[,6], "g__",""),
                        Species = str_replace(taxM[,7], "s__",""),
                        stringsAsFactors = FALSE)
for (i in 1:7){ tax.cleanM[,i] <- as.character(tax.cleanM[,i])}

tax.cleanM<-tax.cleanM %>% 
mutate(Species=str_replace(Species, "_", " "))

phyloseq::tax_table(mock_raw) <- as.matrix(tax.cleanM) ; phyloseq::tax_table(mock_raw)

library("MicrobiotaProcess")

# For the mock experiments, the data were first rarefied

mockR = rarefy_even_depth(mock_raw, 
                               rngseed=1024, 
                               replace=FALSE)

# Because singletons still remains after rarefaction, taxa not supported with 
# at least 10 ASVs, were removed from the final "mock" object. 

mockRF = prune_taxa(taxa_sums(mockR) > 09, mockR); mockRF

saveRDS(mockRF, "R_objects/mock_RF.rds")

###############################################################
########## preprocessing for the field experiments ############
###############################################################

Runs_raw <- readRDS("R_objects/Runs_raw.rds"); Runs_raw

library("tidyverse")
taxR <- data.frame(phyloseq::tax_table(Runs_raw))
tax.cleanR <- data.frame(row.names = row.names(taxR),
                         Kingdom = str_replace(taxR[,1], "k__",""),
                         Phylum = str_replace(taxR[,2], "p__",""),
                         Class = str_replace(taxR[,3], "c__",""),
                         Order = str_replace(taxR[,4], "o__",""),
                         Family = str_replace(taxR[,5], "f__",""),
                         Genus = str_replace(taxR[,6], "g__",""),
                         Species = str_replace(taxR[,7], "s__",""),
                         stringsAsFactors = FALSE)
for (i in 1:7){ tax.cleanR[,i] <- as.character(tax.cleanR[,i])}

tax.cleanR<-tax.cleanR %>% 
  mutate(Species=str_replace(Species, "_", " "))

phyloseq::tax_table(Runs_raw) <- as.matrix(tax.cleanR) ; phyloseq::tax_table(Runs_raw)


# plot depth 

library("microbiome")
library("microbiomeutilities")
library("ggplot2") 


Dep1<-plot_read_distribution(Runs_raw, 
                             groups = "Year", 
                             plot.type = "histogram")+
  theme_biome_utils()+
  scale_x_continuous(trans='log10', 
                     limits = c(200, 15000))+
  scale_fill_manual(values=c("#111111"))+ 
  geom_vline(xintercept = 3000, 
             colour = "black", 
             linetype="dashed")+
  theme(legend.position="none")+
  labs(x = "", y = "Count")

## remove samples with less than 3000 reads 

Runs_P <- prune_samples(sample_sums(Runs_raw) > 4999, Runs_raw); Runs_P

## re-plot depth  

Dep2<-plot_read_distribution(Runs_P, groups = "Year", 
                             plot.type = "histogram")+
  theme_biome_utils()+
  scale_x_continuous(trans='log10', limits = c(200, 150000))+
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

ggsave(file="Figures/supp1_depth_final.pdf", width=8, height=5, units="in", dpi=300)

### Rarefaction 

Runs_PR = rarefy_even_depth(Runs_P, 
                          rngseed=1024, 
                          replace=FALSE)

summarize_phyloseq(Runs_PR)

Runs_PRF = prune_taxa(taxa_sums(Runs_PR) > 09, Runs_PR); Runs_PRF

summarize_phyloseq(Runs_PRF)


saveRDS(Runs_PRF, "R_objects/Runs_PRF.rds")
