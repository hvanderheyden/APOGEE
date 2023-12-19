
    
library("readxl")      

library("tibble")
library("vegan")
library("DESeq2")a # BiocManager::install("DESeq2")
library("speedyseq") # install with remotes::install_github("mikemc/speedyseq") 

library("ggstar")
library("forcats")
library("patchwork")
library("ggpubr")
library("cowplot")


my_cols<- c(  
           "#8dd3c7",
           "#ffed6f",         
           "#80b1d3",
           "#bebada",
           "#920000",
           "#fb8072",
           "#004949",
           "#fdb462",
           "#b3de69",
           "#fccde5",
           "#d9d9d9",
           "#bc80bd",
           "#ccebc5",
           "#ffffb3", 
           "#490092",
           "#009209",
           "#006ddb",
           "#000000"
)

library("phyloseq") #BiocManager::install("phyloseq")

Runs_PRF <- readRDS("R_objects/Runs_PRF.rds");Runs_PRF

##################################################################
# compute and plot stacked bar charts to sumarize by Week/Year####
##################################################################


Runs_PRF_dat <- psmelt(Runs_PRF)

library("forcats")
Runs_PRF_dat$Month <- as_factor(Runs_PRF_dat$Month)
Runs_PRF_dat$Week <- as_factor(Runs_PRF_dat$Week)


str(Runs_PRF_dat)

predefined_species = c("Alternaria alternata",
                       "Botrytis cinerea",
                       "Botrytis squamosa",
                       "Botrytis porrii",
                       "Blumeria graminis",
                       "Cladosporium aphidis",
                       "Cladosporium herbarum",
                       "Epicoccum nigrum",
                       "Ganoderma destructans",
                       "Hyaloperonospora camelinae",
                       "Neoascochyta europaea",
                       "Neosetophoma guiyangensis",
                       "Peniophora tamaricicola",
                       "Peronospora destructor",
                       "Pseudopithomyces rosae",
                       "Rhodotorula babjevae",
                       "Rhodotorula diobovata",
                       "Stemphylium solani",
                       "Stemphylium vesicarium"
                       );predefined_species

#################################### 

library("tidyverse")
library("ggtext")

library("dplyr")

Runs_PRF_dat2 <-  Runs_PRF_dat %>%
  mutate(Species = case_when(
    Species %in% 
      predefined_species ~ Species,  # Keep valid species as is
    TRUE ~ "other"  # Replace other species with "other"
  ));Runs_PRF_dat2

tail(Runs_PRF_dat2)
head(Runs_PRF_dat2)

otu_rel_abund <- Runs_PRF_dat2 %>%
  group_by(Year, Week, Site) %>%
  mutate(rel_abund = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  select(-Abundance) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
               names_to="level",
               values_to="taxon")

otu_rel_abund$taxon

###################################

 run_abund<-otu_rel_abund %>%
  filter(level=="Species") %>%
  group_by(Year, Week, Site, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Year, Week, Site, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  mutate(taxon = factor(taxon, 
                        levels=c("Stemphylium vesicarium",
                                 "Cladosporium herbarum",
                                 "Alternaria alternata",
                                 "Epicoccum nigrum",
                                 "Peronospora destructor",
                                 "Botrytis cinerea",
                                 "Botrytis squamosa",
                                 "Botrytis porrii",
                                 "Blumeria graminis",
                                 "Cladosporium aphidis",
                                 "Ganoderma destructans",
                                 "Hyaloperonospora camelinae",
                                 "Neoascochyta europaea",
                                 "Neosetophoma guiyangensis",
                                 "Peniophora tamaricicola",
                                 "Pseudopithomyces rosae",
                                 "Rhodotorula babjevae",
                                 "Rhodotorula diobovata",
                                 "Stemphylium solani",
                                 "Other")))%>%
           mutate(Week = factor(Week, 
                                 levels=c("23", "24", "25", "26",
                                          "27", "28", "29", "30",
                                          "31", "32", "33", "34")))


# plot the stacked bar chart 
run_stacked<-ggplot(data=run_abund, 
                     aes(x=Week, 
                         y=mean_rel_abund, 
                         fill=taxon)) +
  geom_col(colour = "black", width=0.8, linewidth=0.1) +
  facet_wrap(vars(Year, Site), nrow = 4)+
  theme(legend.title=element_blank())+
  labs(x="Week number",
       y="Relative Abundance (%)") +
  theme(legend.text = element_text(face="italic"))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=my_cols)+
  theme(legend.position="bottom")+
  guides(fill= guide_legend(keywidth = 0.6, 
                            keyheight = 0.7, 
                            ncol=4))+
  theme(axis.text.x = element_text(angle = 60, 
                                   vjust = 0.5, 
                                   hjust=0.4));run_stacked

ggsave(file="Figures/Fig3_Run_stacked.pdf", width=6, height=6, units="in", dpi=900)








