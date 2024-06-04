#trying to get the the heart of the data viz, not going through all of the decom

library(dplyr)
library(ggplot2)
install.packages("phyloseq")
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library(phyloseq)
install.packages("ranacapa")
library(ranacapa)
library(devtools)
library(vegan)

source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
       local = TRUE)
devtools::install_github("joey711/phyloseq")

tax_table18S <- read.csv("/home/shared/8TB_HDD_02/cnmntgna/GitHub/YellowIsland2023/data/Anacapa_Results/18S_taxonomy_tables/Summary_by_percent_confidence/90/18S_ASV_sum_by_taxonomy_90.txt", header= TRUE, sep = "\t")

head(tax_table18S)
summary(tax_table18S)

# Convert the matrix into a phyloseq otu_table object, with taxa as the rows
ana_taxon_table_physeq <- phyloseq::otu_table(tax_table18S, taxa_are_rows = TRUE)

# Extract the rownames of the matrix above- this has the full taxonomic path.
# Split the taxonomic path on semicolons, and turn the resulting matrix into
# a phyloseq tax_table object
taxon_names <- reshape2::colsplit(rownames(tax_table18S), ";",
                                  names = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix
rownames(taxon_names) <- rownames(tax_table18S)

tax_physeq <- phyloseq::tax_table(taxon_names)
colnames(tax_physeq) <- c("Domain","Phylum","Class","Order","Family","Genus","Species")




# Example: Summarize counts for each taxon
summarized_18S <- tax_table18S %>%
  group_by(Taxon) %>%
  summarize(total_counts = sum(Counts))
