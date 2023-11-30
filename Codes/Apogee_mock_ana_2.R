
library("phyloseq")

mock <- readRDS("R_objects/mock.rds")

sample_variables(mock)
rank_names(mock)

mock <- subset_samples(mock, Name != "mock")
mock <- subset_samples(mock, Name != "mock2")

library("microbiome")
summarize_phyloseq(mock)

library("MicrobiotaProcess")

mockR = rarefy_even_depth(mock, 
                          rngseed=1024, 
                          replace=F)

summarize_phyloseq(mockR)

mockRF=prune_taxa(taxa_sums(mockR) > 19, mockR); mockRF

summarize_phyloseq(mockRF)

#
################################################
mock_dat <- psmelt(mockRF)

str(mock_dat)

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
                                       "Other")))

summary <- mock_abund %>% 
  group_by(Name, Description, taxon) %>%
  summarize(mean_rel_abund) %>% 
pivot_wider(names_from = Description, values_from = mean_rel_abund)

write.csv(summary, "Data/mocks/summary_mock.csv", row.names=FALSE)


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


dat_text <- data.frame(label = c("A", "B", "C"))

mock_stacked +
  geom_text(
  data    = dat_text,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -0.1,
  vjust   = -1
)



  ################################################

library(MicrobiotaProcess)
  
  mock3<- subset_samples(mockRF, 
                         Name =="mock3")
  
  distmeDF <- get_dist(mock3, 
                       distmethod ="wunifrac")
  distmeDF

  
  
  
    
                     