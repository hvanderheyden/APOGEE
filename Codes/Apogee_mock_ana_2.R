

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
  theme(strip.text = element_blank())+
  theme(legend.title=element_blank())+
  labs(x=NULL,
       y="Mean Relative Abundance (%)") +
  theme(legend.text = element_text(face="italic"))+
  scale_y_continuous(expand=c(0,10))+
  scale_fill_manual(values=c("#a91919","#582630","#998540","black",
                             "#D55E00","#0072B2","#F9ECCC","#679289", 
                             "#33658A","#F6AE2D","#00AFBB","grey"))+
  theme(legend.position="right")+
  guides(fill= guide_legend(keywidth = 0.6, 
                            keyheight = 0.7, 
                            ncol=1))+
  theme(legend.position = c(0.9, -0.05),
        legend.justification = c(0.9, -0.05))+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.4))


dat_text <- data.frame(Description=c("mock3", "mock4", "mock5"), 
                       label = c("A", "B", "C"))
                      
mock_stacked +
  geom_text(x = 6, y = 40, aes(label = label), data = dat_text)


  ################################################



library(MicrobiotaProcess)
  
  mock3<- subset_samples(mockRF, Name =="mock3")
  
  distm3 <- get_dist(mock3, distmethod ="wunifrac")

  mock4<- subset_samples(mockRF,Name =="mock4")
  
  distmem4 <- get_dist(mock4, distmethod ="wunifrac")
  
  mock5<- subset_samples(mockRF,Name =="mock5")
  
  distmem5 <- get_dist(mock5, distmethod ="wunifrac")

  
  distm3
  distmem4
  distmem5
    
                     