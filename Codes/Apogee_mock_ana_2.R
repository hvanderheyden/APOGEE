

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


write.csv(summary, "Data/mocks/summary_mock.csv", row.names=FALSE)


# plot the stacked bar chart 
mock_stacked<-ggplot(data=mock_abund, 
         aes(x=Description, 
             y=mean_rel_abund, 
             fill=taxon)) +
  geom_col() +
  facet_wrap(vars(Name), nrow = 2)+
  #theme(strip.text = element_blank())+
  theme(legend.title=element_blank())+
  labs(x=NULL,
       y="Relative Abundance (%)") +
  theme(legend.text = element_text(face="italic"))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=c(
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
    "#ffed6f"))+
  theme(legend.position="right")+
  guides(fill= guide_legend(keywidth = 0.6, 
                            keyheight = 0.7, 
                            ncol=1))+
  theme(legend.position = c(0.9, -0.05),
        legend.justification = c(0.9, -0.05))+
  theme(axis.text.x = element_text(angle = 60, 
                                   vjust = 0.5, 
                                   hjust=0.4));mock_stacked



################################################

# estimate weighted unifrac distance (dissimilarity)

library(MicrobiotaProcess)
  
  mock1<- subset_samples(mock, Name =="mock")
  
  distm <- get_dist(mock1, distmethod ="wunifrac")

  mock2<- subset_samples(mock, Name =="mock2")
  
  distm2 <- get_dist(mock2, distmethod ="wunifrac")
  
  mock3<- subset_samples(mock,Name =="mock3")
  
  distm3 <- get_dist(mock3, distmethod ="wunifrac")

  
  distm
  distm2
  distm3
    
                     