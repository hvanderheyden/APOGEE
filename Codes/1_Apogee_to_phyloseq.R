

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

#######################################################################
############## build the phyloseq object for the mocks ################
#######################################################################

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
######## Qiime2 ##########
##############################

## load the biom_table, taxonomy and metadata ####
biom_qiime2<- read.csv("Data/mocks/Qiime2_vsearch_table.txt", 
                     header=TRUE, sep="\t")
head(biom_qiime2)

taxo_qiime2<- read.csv("Data/mocks/Qiime2_vsearch_taxonomy.tsv", 
                     header=TRUE, sep="\t")
head(taxo_qiime2)

# define the row names from the otu column ####

biom_qiime2 <- biom_qiime2 %>%
  tibble::column_to_rownames("otu") 

taxo_qiime2 <- taxo_qiime2 %>%
  tibble::column_to_rownames("otu") 

# Transform into matrixes otu and tax tables (sample table can be left as data frame) ####

biom_qiime2 <- as.matrix(biom_qiime2)
taxo_qiime2 <- as.matrix(taxo_qiime2)

class(taxo_qiime2)
class(biom_qiime2)

# convert to phyloseq objects ####
OTU_qiime2 = otu_table(biom_qiime2, taxa_are_rows = TRUE)
TAX_qiime2 = phyloseq::tax_table(taxo_qiime2)


phylo_qiime2 <- phyloseq(OTU_qiime2, TAX_qiime2)
phylo_qiime2

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
OTU_minimap2 = otu_table(biom_minimap2, taxa_are_rows = TRUE)
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
mock_raw <- merge_phyloseq(otu, tax, samples)

mock_raw

#create and merge 
library("ape")
random_tree_mock = rtree(ntaxa(mock_raw), 
                    rooted=TRUE, 
                    tip.label=taxa_names(mock_raw))

mock_raw <- phyloseq(otu, 
                           tax, 
                           samples, 
                           random_tree_mock)
mock_raw

# inspect your  phyloseq object: 
sample_names(mock_raw)
rank_names(mock_raw)
sample_variables(mock_raw)

saveRDS(mock_raw, file= "R_objects/mock_raw.rds")

#######################################################################
############## build the phyloseq object for the runs #################
#######################################################################


## load the biom_table, taxonomy and metadata ####
Apogee_biom<- read.csv("Data/runs/otu_table.csv", 
                       header=TRUE, 
                       sep=",")
head(Apogee_biom)

Apogee_taxo<- read.csv("Data/runs/phyloseq_taxonomy.csv", 
                       header=TRUE, 
                       sep=",")
head(Apogee_taxo)

Apogee_meta <- read.csv("Data/runs/Meta_full.tsv", 
                        header=TRUE, 
                        sep="\t")
head(Apogee_meta)

library("dplyr")
Apogee_meta$DOY<-as.factor(Apogee_meta$DOY)
Apogee_meta$Year<-as.factor(Apogee_meta$Year)

Apogee_meta <- Apogee_meta %>% 
  tibble::column_to_rownames("Sample_ID")

library("phyloseq")
samplesS = sample_data(Apogee_meta)

# define the row names from the otu column ####

Apogee_biom <- Apogee_biom %>%
  tibble::column_to_rownames("OTU") 

Apogee_taxo <- Apogee_taxo %>%
  tibble::column_to_rownames("OTU") 

# Transform into matrixes otu and tax tables (sample table can be left as data frame) ####

Apogee_biom <- as.matrix(Apogee_biom)
Apogee_taxo <- as.matrix(Apogee_taxo)

class(Apogee_biom)
class(Apogee_taxo)

# convert to phyloseq objects ####

Apogee_OTU = otu_table(Apogee_biom, taxa_are_rows = TRUE)
Apogee_TAX = phyloseq::tax_table(Apogee_taxo)


Apogee_run_raw <- phyloseq(Apogee_OTU, 
                      Apogee_TAX, 
                      samplesS)
Apogee_run_raw

sample_variables(Apogee_run_raw)

#create and merge 
library("ape")
random_tree = rtree(ntaxa(Apogee_run_raw), 
                    rooted=TRUE, 
                    tip.label=taxa_names(Apogee_run_raw))

Apogee_run_raw <- phyloseq(Apogee_OTU, 
                      Apogee_TAX, 
                      samplesS, 
                      random_tree)
Apogee_run_raw

saveRDS(Apogee_run_raw, file= "R_objects/Apogee_run_raw.rds")

