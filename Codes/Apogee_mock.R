



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


#######################################
########## Bracken ####################
#######################################

## load the biom_table, taxonomy and metadata ####
biom_bracken<- read.csv("C:/Users/vanderheydenh/OneDrive - AGR-AGR/Projets/2023/Biovigilance/APOGEE/Apogee_mock/To_phyloseq/bracken.biom.tsv", header=TRUE, sep="\t")
head(biom_bracken)

taxo_bracken<- read.csv("C:/Users/vanderheydenh/OneDrive - AGR-AGR/Projets/2023/Biovigilance/APOGEE/Apogee_mock/To_phyloseq/bracken.taxo.csv", header=TRUE, sep=";")
head(taxo_bracken)

# define the row names from the otu column ####

biom_bracken <- biom_bracken %>%
  tibble::column_to_rownames("otu") 

taxo_bracken <- taxo_bracken %>%
  tibble::column_to_rownames("otu") 

# Transform into matrixes otu and tax tables (sample table can be left as data frame) ####

biom_bracken <- as.matrix(biom_bracken)
taxo_bracken <- as.matrix(taxo_bracken)

class(taxo_bracken)
class(biom_bracken)

# convert to phyloseq objects ####
OTU_bracken = otu_table(biom_bracken, taxa_are_rows = TRUE)
TAX_bracken = phyloseq::tax_table(taxo_bracken)

phylo_bracken <- phyloseq(OTU_bracken, TAX_bracken)
phylo_bracken

#######################################
########## Kracken ####################
#######################################

## load the biom_table, taxonomy and metadata ####
biom_kracken<- read.csv("C:/Users/vanderheydenh/OneDrive - AGR-AGR/Projets/2023/Biovigilance/APOGEE/Apogee_mock/To_phyloseq/kracken.biom.tsv", header=TRUE, sep="\t")
head(biom_kracken)

taxo_kracken<- read.csv("C:/Users/vanderheydenh/OneDrive - AGR-AGR/Projets/2023/Biovigilance/APOGEE/Apogee_mock/To_phyloseq/kraken.taxo.tsv", header=TRUE, sep="\t")
head(taxo_kracken)

# define the row names from the otu column ####

biom_kracken <- biom_kracken %>%
  tibble::column_to_rownames("otu") 

taxo_kracken <- taxo_kracken %>%
  tibble::column_to_rownames("otu") 

# Transform into matrixes otu and tax tables (sample table can be left as data frame) ####

biom_kracken <- as.matrix(biom_kracken)
taxo_kracken <- as.matrix(taxo_kracken)

class(taxo_kracken)
class(biom_kracken)

# convert to phyloseq objects ####
OTU_kracken = otu_table(biom_kracken, taxa_are_rows = TRUE)
TAX_kracken = phyloseq::tax_table(taxo_kracken)


phylo_kracken <- phyloseq(OTU_kracken, TAX_kracken)
phylo_kracken

##############################
######## Metontiime ##########
##############################

phylo_metontiime <- qza_to_phyloseq(
  features="C:/Users/vanderheydenh/OneDrive - AGR-AGR/Projets/2023/Biovigilance/APOGEE/Apogee_mock/To_phyloseq/metontiime_table.qza",
  taxonomy="C:/Users/vanderheydenh/OneDrive - AGR-AGR/Projets/2023/Biovigilance/APOGEE/Apogee_mock/To_phyloseq/metontiime_taxonomy.qza"
)

sample_names(phylo_metontiime)

#######################################
########## Spaghetti ##################
#######################################


## load the biom_table, taxonomy and metadata ####
biom_spaghetti<- read.csv("C:/Users/vanderheydenh/OneDrive - AGR-AGR/Projets/2023/Biovigilance/APOGEE/Apogee_mock/To_phyloseq/spaghetti_otu_table.csv", header=TRUE, sep=";")
head(biom_spaghetti)

taxo_spaghetti<- read.csv("C:/Users/vanderheydenh/OneDrive - AGR-AGR/Projets/2023/Biovigilance/APOGEE/Apogee_mock/To_phyloseq/spaghetti_taxonomy.csv", header=TRUE, sep=";")
head(taxo_spaghetti)

# define the row names from the otu column ####

biom_spaghetti <- biom_spaghetti %>%
  tibble::column_to_rownames("otu") 

taxo_spaghetti <- taxo_spaghetti %>%
  tibble::column_to_rownames("otu") 

# Transform into matrixes otu and tax tables (sample table can be left as data frame) ####

biom_spaghetti <- as.matrix(biom_spaghetti)
taxo_spaghetti <- as.matrix(taxo_spaghetti)

class(taxo_kracken)
class(biom_spaghetti)

# convert to phyloseq objects ####
OTU_spaghetti = otu_table(biom_spaghetti, taxa_are_rows = TRUE)
TAX_spaghetti = phyloseq::tax_table(taxo_spaghetti)


phylo_spaghetti <- phyloseq(OTU_spaghetti, TAX_spaghetti)
phylo_spaghetti

###################################
# finally merge the 5 PS objects  #
###################################

# extract the tax tables from PS objects 1 to 5 
tax1  = phyloseq::tax_table(phylo_bracken)
tax2 = phyloseq::tax_table(phylo_kracken)
tax3 = phyloseq::tax_table(phylo_metontiime)
tax4  = phyloseq::tax_table(phylo_spaghetti)

# merge the tax tables 
tax <- merge_phyloseq(tax1, tax2, tax3, tax4)

# extract the otu tables from PS objects 1 to 5 
otu1  = otu_table(phylo_bracken)
otu2 = otu_table(phylo_kracken)
otu3 = otu_table(phylo_metontiime)
otu4  = otu_table(phylo_spaghetti)

# merge the otu tables 
otu <- merge_phyloseq(otu1, otu2, otu3, otu4)


sample_names(phylo_bracken)
sample_names(phylo_kracken)
sample_names(phylo_metontiime)
rank_names(phylo_metontiime)
sample_names(phylo_spaghetti)

# concatenate the meta data into a single file and import 
meta <- read.csv("C:/Users/vanderheydenh/OneDrive - AGR-AGR/Projets/2023/Biovigilance/APOGEE/Apogee_mock/To_phyloseq/meta.csv", header=TRUE, sep=",")
head(meta)

meta <- meta %>% 
  tibble::column_to_rownames("Sample_ID")

samples = sample_data(meta)

# merge the otu, tax and meta objects into a final PS object: 
mock <- merge_phyloseq(otu, tax, samples)

random_tree = rtree(ntaxa(mock), rooted=TRUE, tip.label=taxa_names(mock))

mock <- merge_phyloseq(mock, random_tree)
mock

# inspect your  phyloseq object: 
sample_names(mock)
rank_names(mock)
sample_variables(mock)
tax_table(mock)

saveRDS(mock, file= "C:/Users/vanderheydenh/OneDrive - AGR-AGR/Projets/2023/Biovigilance/APOGEE/Apogee_mock/To_phyloseq/mock.rds")

mock <- readRDS(file = "C:/Users/vanderheydenh/OneDrive - AGR-AGR/Projets/2023/Biovigilance/APOGEE/Apogee_mock/To_phyloseq/mock.rds")


