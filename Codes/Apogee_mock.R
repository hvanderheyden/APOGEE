



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

getwd()



#######################################
########## Bracken ####################
#######################################

## load the biom_table, taxonomy and metadata ####
biom_bracken<- read.csv("Data/mocks/bracken_table.biom", 
                        header=TRUE, sep="\t")
head(biom_bracken)

taxo_bracken<- read.csv("Data/mocks/bracken_taxo.tsv", 
                        header=TRUE, sep="\t")
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
biom_kracken<- read.csv("Data/mocks/kraken2_table.biom", 
                        header=TRUE, sep="\t")
head(biom_kracken)

taxo_kracken<- read.csv("Data/mocks/kraken2_taxo.tsv", 
                        header=TRUE, sep="\t")
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

phylo_qiime2 <- qza_to_phyloseq(
  features="Data/mocks/Qiime2_blast_table.qza",
  taxonomy="Data/mocks/Qiime2_blast_taxonomy.qza"
)

sample_names(phylo_qiime2)

#######################################
########## minimap2 ##################
#######################################


## load the biom_table, taxonomy and metadata ####
biom_minimap2<- read.csv("Data/mocks/minimap2_otu_table.csv", header=TRUE, sep=";")
head(biom_minimap2)

taxo_minimap2<- read.csv("Data/mocks/minimap2_taxo.csv", header=TRUE, sep=";")
head(taxo_minimap2)

# define the row names from the otu column ####

biom_minimap2 <- biom_minimap2 %>%
  tibble::column_to_rownames("otu") 

taxo_minimap2 <- taxo_minimap2 %>%
  tibble::column_to_rownames("otu") 

# Transform into matrixes otu and tax tables (sample table can be left as data frame) ####

biom_minimap2 <- as.matrix(biom_minimap2)
taxo_minimap2 <- as.matrix(taxo_minimap2)

class(taxo_minimap2)
class(biom_minimap2)

# convert to phyloseq objects ####
OTU_minimap2 = otu_table(biom_spaghetti, taxa_are_rows = TRUE)
TAX_minimap2 = phyloseq::tax_table(taxo_minimap2)


phylo_minimap2 <- phyloseq(OTU_minimap2, TAX_minimap2)
phylo_minimap2

###################################
###########   Theoretical #########
###################################

## load the biom_table, taxonomy and metadata ####
biom_theo<- read.csv("Data/mocks/theor_otu_table.csv", 
                     header=TRUE, sep=";")
head(biom_theo)

taxo_theo<- read.csv("Data/mocks/theor_taxo.csv", 
                     header=TRUE, sep=";")
head(taxo_theo)

# define the row names from the otu column ####

biom_theo <- biom_theo %>%
  tibble::column_to_rownames("otu") 

taxo_theo <- taxo_theo %>%
  tibble::column_to_rownames("otu") 

# Transform into matrixes otu and tax tables (sample table can be left as data frame) ####

biom_theo <- as.matrix(biom_theo)
taxo_theo <- as.matrix(taxo_theo)

class(taxo_theo)
class(biom_theo)

# convert to phyloseq objects ####
OTU_theo = otu_table(biom_theo, taxa_are_rows = TRUE)
TAX_theo = phyloseq::tax_table(taxo_theo)


phylo_theo <- phyloseq(OTU_theo, TAX_theo)
phylo_theo


###################################
# finally merge the 5 PS objects  #
###################################

# extract the tax tables from PS objects 1 to 5 
tax1  = phyloseq::tax_table(phylo_bracken)
tax2 = phyloseq::tax_table(phylo_kracken)
tax3 = phyloseq::tax_table(phylo_qiime2)
tax4  = phyloseq::tax_table(phylo_minimap2)
tax5  = phyloseq::tax_table(phylo_theo)

# merge the tax tables 
tax <- merge_phyloseq(tax1, tax2, tax3, tax4, tax5)
str(tax)

head(tax)

tax <- as.matrix(tax)

tax = phyloseq::tax_table(tax)

# extract the otu tables from PS objects 1 to 5 
otu1 = otu_table(phylo_bracken)
otu2 = otu_table(phylo_kracken)
otu3 = otu_table(phylo_qiime2)
otu4  = otu_table(phylo_minimap2)
otu5  = otu_table(phylo_theo)


# merge the otu tables 
otu <- merge_phyloseq(otu1, otu2, otu3, otu4, otu5)

sample_names(phylo_bracken)
sample_names(phylo_kracken)
sample_names(phylo_qiime2)
rank_names(phylo_qiime2)
sample_names(phylo_minimap2)

# concatenate the meta data into a single file and import 
meta <- read.csv("Data/mocks/meta.csv", 
                 header=TRUE, sep=";")
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

saveRDS(mock, file= "R_objects/mock.rds")



