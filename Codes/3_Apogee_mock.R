

library("phyloseq")

# Load the .rds object 

mock <- readRDS("R_objects/mock_RF.rds")

################################################
mock_dat <- psmelt(mock)

str(mock_dat)

predefined_species = c("Botrytis porri",
             "Alternaria alternata",
             "Botrytis cinerea",
             "Peronospora variabilis",
             "Epicoccum nigrum",
             "Stemphylium vesicarium",
             "Peronospora destructor",
             "Sclerotinia sclerotiorum",
             "Alternaria brassicicola",
             "Botrytis squamosa",
             "Fusarium oxysporum",
             "Hyaloperonospora brassicae")

#################################### 

library(tidyverse)
library(ggtext)

library("dplyr")
mock_dat <- mock_dat %>%
  mutate(Species = case_when(
    Species %in% 
      predefined_species ~ Species,  # Keep valid species as is
    TRUE ~ "Other"  # Replace other species with "other"
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
                              levels=c("Theoretical", "Minimap2_F", "Minimap2", "Bracken", "Kraken", "Qiime2")))


###################################

mock_abund<-otu_rel_abund %>%
  filter(level=="Species") %>%
  group_by(Name, Description, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Name, Description, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  mutate(taxon = factor(taxon, 
                              levels=c("Botrytis cinerea",
                                       "Botrytis porri",
                                       "Botrytis squamosa",
                                       "Sclerotinia sclerotiorum",
                                       "Alternaria alternata",
                                       "Alternaria brassicicola",
                                       "Peronospora variabilis",
                                       "Peronospora destructor",
                                       "Hyaloperonospora brassicae",
                                       "Epicoccum nigrum",
                                       "Stemphylium vesicarium",
                                       "Fusarium oxysporum",
                                       "Other")))

summary <- mock_abund %>% 
  group_by(Name, Description, taxon) %>%
  summarize(mean_rel_abund) %>% 
pivot_wider(names_from = Description, values_from = mean_rel_abund)


#write.csv(summary, "Data/mocks/summary_mock.csv", row.names=FALSE)

my_cols<- c(  
"#8dd3c7",
"#ffffb3",
"#bebada",
"#fb8072",
"#80b1d3",
"#fdb462",
"#b3de69",
"#fccde5",
"#d9d9d9",
"#bc80bd",
"#ccebc5",
"#ffed6f",
"black")

# plot the stacked bar chart 
mock_stacked<-ggplot(data=mock_abund, 
         aes(x=Description, 
             y=mean_rel_abund, 
             fill=taxon)) +
  geom_col(colour = "black", width=0.8, linewidth=0.1) +
  facet_wrap(vars(Name), nrow = 2)+
  #theme(strip.text = element_blank())+
  theme(legend.title=element_blank())+
  labs(x=NULL,
       y="Relative Abundance (%)") +
  theme(legend.text = element_text(face="italic"))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=my_cols)+
  theme(legend.position="right")+
  guides(fill= guide_legend(keywidth = 0.6, 
                            keyheight = 0.7, 
                            ncol=1))+
  theme(legend.position = c(0.9, -0.05),
        legend.justification = c(0.9, -0.05))+
  theme(axis.text.x = element_text(angle = 60, 
                                   vjust = 0.5, 
                                   hjust=0.4))+
  theme(legend.background=element_rect(fill = alpha("white", 0.1)));mock_stacked


################################################

# estimate weighted unifrac distance (dissimilarity)

library("MicrobiotaProcess")
library("ggpubr")
library("ape")

taxaPS <- get_taxadf(mock, 
                     taxlevel=7,
                     type = "species")

taxaPS_tree = rtree(ntaxa(taxaPS), 
                    rooted=TRUE, 
                    tip.label=taxa_names(taxaPS))

taxaPS_tax  = phyloseq::tax_table(taxaPS)
taxaPS_otu = otu_table(taxaPS)

meta <- read.csv("Data/mocks/meta.csv", 
                 header=TRUE, sep=";")
head(meta)

meta <- meta %>% 
  tibble::column_to_rownames("Sample_ID")

samples = sample_data(meta)

taxaPS <- phyloseq(taxaPS_tax,
                   taxaPS_otu,
                   samples,
                   taxaPS_tree)



###############################################

mock <- subset_samples(taxaPS, Name  == "mock")
mock2 <- subset_samples(taxaPS, Name  == "mock2")
mock3 <- subset_samples(taxaPS, Name  == "mock3")

distm <- get_dist(mock, 
                    distmethod ="bray", 
                    method="hellinger");distm

distm2 <- get_dist(mock2, 
                  distmethod ="bray", 
                  method="hellinger");distm2

distm3 <- get_dist(mock3, 
                  distmethod ="bray",
                  method="hellinger");distm3
  
#############################
# Get the PCOA results 
pcoa_results <- get_pcoa(obj=taxaPS, 
                      distmethod="bray", 
                      method="hellinger")
  
# Visualizing the PCOA results
pcoaplot1 <- ggordpoint(obj=pcoa_results, 
                          biplot=FALSE,
                          factorNames=c("Description"), 
                          ellipse=FALSE,
                          poinsize = 1.5,
                          stroke = 0.8)+
  xlim(-0.45, 0.45)+
  ylim(-0.45, 0.45)+
  theme(legend.box = "horizontal", legend.position = c(0.2, 0.85))+
  theme(legend.text=element_text(size=8))+
  theme(legend.title = element_blank())+ 
  theme(plot.title = element_blank())+
  guides(fill = guide_legend(keywidth = 0.6, 
                             keyheight = 0.7,
                             ncol = 2))+
  theme(legend.background=element_rect(fill = alpha("white", 0.1)));pcoaplot1



ggarrange(ggarrange(mock_stacked,
                    pcoaplot1,
                    nrow = 2, labels = c("A", "B"),
                    heights = c(2.5, 1)))

ggsave(file="Figures/Fig2_mock.pdf", width=5, height=8, units="in", dpi=900)


        