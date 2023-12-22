---
title: "CIB_Bycatch"
author: "Rachel Meyer"
date: "2023-09-18"
output: html_document
---
#get decontam installed if you haven't 
library(devtools)


#devtools::install_github("benjjneb/decontam")

library(decontam)

#set working directory
setwd("/Users/rachelmeyer/Dropbox/Tronko_Outputs/Bycatch/assign")

#read in files
CIB_FITS_30<-read.csv("Bycatch_FITS_Max30_chosen.txt", header=TRUE, sep="\t")
CIB_FITS_5<-read.csv("Bycatch_FITS_Max5_chosen.txt", header=TRUE, sep="\t")
CIB_CO1_30<-read.csv("Bycatch_CO1_Max30_chosen.txt", header=TRUE, sep="\t")
CIB_CO1_5<-read.csv("Bycatch_CO1_Max5_chosen.txt", header=TRUE, sep="\t")


biom<-read.csv("ByCatch_MetadataOutput_Metabarcoding.txt", header=TRUE, sep="\t")

library(phyloseq)
library(ranacapa)
library(devtools)
library(ggplot2)
library(vegan)


#Import some ranacapa functions

convert_anacapa_to_phyloseq <- function(taxon_table, metadata_file) {
  
  # Validate the files
  validate_input_files(taxon_table, metadata_file)
  
  # Group the anacapa ouptut by taxonomy, if it has not yet happened, and turn it into a matrix
  taxon_table2 <- group_anacapa_by_taxonomy(taxon_table) %>%
    tibble::column_to_rownames("sum.taxonomy") %>%
    as.matrix
  # Reorder the columns (sites) for ease of displaying later
  taxon_table2 <- taxon_table2[ , order(colnames(taxon_table2))]
  
  # Convert the matrix into a phyloseq otu_table object, with taxa as the rows
  ana_taxon_table_physeq <- phyloseq::otu_table(taxon_table2, taxa_are_rows = TRUE)
  
  # Extract the rownames of the matrix above- this has the full taxonomic path.
  # Split the taxonomic path on semicolons, and turn the resulting matrix into
  # a phyloseq tax_table object
  taxon_names <- reshape2::colsplit(rownames(taxon_table2), ";",
                                    names = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
    as.matrix
  rownames(taxon_names) <- rownames(taxon_table2)
  
  tax_physeq <- phyloseq::tax_table(taxon_names)
  colnames(tax_physeq) <- c("Domain","Phylum","Class","Order","Family","Genus","Species")
  
  # Make a phyloseq object out of the otu_table and the tax_table objects
  physeq_object <- phyloseq::phyloseq(ana_taxon_table_physeq, tax_physeq)
  
  # Make sure the mapping file (ie the site metadata) is ordered according to site name
  rownames(metadata_file) <- metadata_file[, 1]
  metadata_file <- metadata_file[order(metadata_file[, 1]), ]
  
  # Convert the mapping file into a phyloseq sample_data object, and merge it with the
  # phyloseq object created above to make a phyloseq object with otu table, tax table, and sample data.
  sampledata <- phyloseq::sample_data(metadata_file)
  phyloseq::merge_phyloseq(physeq_object, sampledata)
}

vegan_otu <- function(physeq_object) {
  OTU <- phyloseq::otu_table(physeq_object)
  if (phyloseq::taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(methods::as(OTU, "matrix"))
}

#now make a phyloseq object
physeqCIB_FITS_30<-convert_anacapa_to_phyloseq(CIB_FITS_30,biom)
physeqCIB_FITS_5<-convert_anacapa_to_phyloseq(CIB_FITS_5,biom)
physeqCIB_CO1_30<-convert_anacapa_to_phyloseq(CIB_CO1_30,biom)
physeqCIB_CO1_5<-convert_anacapa_to_phyloseq(CIB_CO1_5,biom)


#summary_taxa_decontam

ntaxa(physeqCIB_FITS_30)
ntaxa(physeqCIB_FITS_5)
ntaxa(physeqCIB_CO1_30)
ntaxa(physeqCIB_CO1_5)


#decontaminate CO1_5
table(sample_data(physeqCIB_CO1_5)[,"Sample_or_Control"])


df <- as.data.frame(sample_data(physeqCIB_CO1_5))

library(ggplot2)

library(ggpubr)

library(decontam)

df$LibrarySize <- sample_sums(physeqCIB_CO1_5)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
p<-ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
ggsave("physeqCIB_CO1_5_libsizes_predecontam.pdf",p, width=6, height=5, units="cm", scale=3, limitsize = FALSE)
 
 #try 0.1 decontam
sample_data(physeqCIB_CO1_5)$Sample_or_Control == "Control"
sample_data(physeqCIB_CO1_5)$is.neg <- sample_data(physeqCIB_CO1_5)$Sample_or_Control == "Control"
contamdf.prev01 <- isContaminant(physeqCIB_CO1_5, method="prevalence", neg="is.neg", threshold=0.1)
head(which(contamdf.prev01$contaminant))
table(contamdf.prev01$contaminant)
physeqCIB_CO1_5_dc<-prune_taxa(!contamdf.prev01$contaminant, physeqCIB_CO1_5)
df <- as.data.frame(sample_data(physeqCIB_CO1_5_dc))
df$LibrarySize <- sample_sums(physeqCIB_CO1_5_dc)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
p<-ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
ggsave("physeqCIB_CO1_5_libsizes_postdecontam0.1.pdf",p, width=6, height=5, units="cm", scale=3, limitsize = FALSE)

physeqCIB_CO1_5_dct = filter_taxa(physeqCIB_CO1_5_dc,function(x) sum(x > 0) > 1, TRUE)

write.table(data.frame("sum.taxonomy"=rownames(otu_table(physeqCIB_CO1_5_dct)),otu_table(physeqCIB_CO1_5_dct)),"physeqCIB_CO1_5_dct_0.1.txt", row.names=FALSE, sep="\t",quote = FALSE)

#decontaminate FITS_5
table(sample_data(physeqCIB_FITS_5)[,"Sample_or_Control"])
df <- as.data.frame(sample_data(physeqCIB_FITS_5))
library(ggplot2)
library(ggpubr)
library(decontam)
df$LibrarySize <- sample_sums(physeqCIB_FITS_5)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
p<-ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
ggsave("physeqCIB_FITS_5_libsizes_predecontam.pdf",p, width=6, height=5, units="cm", scale=3, limitsize = FALSE)
 #try 0.1 decontam
sample_data(physeqCIB_FITS_5)$Sample_or_Control == "Control"
sample_data(physeqCIB_FITS_5)$is.neg <- sample_data(physeqCIB_FITS_5)$Sample_or_Control == "Control"
contamdf.prev01 <- isContaminant(physeqCIB_FITS_5, method="prevalence", neg="is.neg", threshold=0.1)
head(which(contamdf.prev01$contaminant))
table(contamdf.prev01$contaminant)
physeqCIB_FITS_5_dc<-prune_taxa(!contamdf.prev01$contaminant, physeqCIB_FITS_5)
df <- as.data.frame(sample_data(physeqCIB_FITS_5_dc))
df$LibrarySize <- sample_sums(physeqCIB_FITS_5_dc)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
p<-ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
ggsave("physeqCIB_FITS_5_libsizes_postdecontam0.1.pdf",p, width=6, height=5, units="cm", scale=3, limitsize = FALSE)
physeqCIB_FITS_5_dct = filter_taxa(physeqCIB_FITS_5_dc,function(x) sum(x > 0) > 1, TRUE)
write.table(data.frame("sum.taxonomy"=rownames(otu_table(physeqCIB_FITS_5_dct)),otu_table(physeqCIB_FITS_5_dct)),"physeqCIB_FITS_5_dct_0.1.txt", row.names=FALSE, sep="\t",quote = FALSE)

#decontaminate CO1_30
table(sample_data(physeqCIB_CO1_30)[,"Sample_or_Control"])
df <- as.data.frame(sample_data(physeqCIB_CO1_30))
library(ggplot2)
library(ggpubr)
library(decontam)
df$LibrarySize <- sample_sums(physeqCIB_CO1_30)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
p<-ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
ggsave("physeqCIB_CO1_30_libsizes_predecontam.pdf",p, width=6, height=5, units="cm", scale=3, limitsize = FALSE)
 #try 0.1 decontam
sample_data(physeqCIB_CO1_30)$Sample_or_Control == "Control"
sample_data(physeqCIB_CO1_30)$is.neg <- sample_data(physeqCIB_CO1_30)$Sample_or_Control == "Control"
contamdf.prev01 <- isContaminant(physeqCIB_CO1_30, method="prevalence", neg="is.neg", threshold=0.1)
head(which(contamdf.prev01$contaminant))
table(contamdf.prev01$contaminant)
physeqCIB_CO1_30_dc<-prune_taxa(!contamdf.prev01$contaminant, physeqCIB_CO1_30)
df <- as.data.frame(sample_data(physeqCIB_CO1_30_dc))
df$LibrarySize <- sample_sums(physeqCIB_CO1_30_dc)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
p<-ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
ggsave("physeqCIB_CO1_30_libsizes_postdecontam0.1.pdf",p, width=6, height=5, units="cm", scale=3, limitsize = FALSE)
physeqCIB_CO1_30_dct = filter_taxa(physeqCIB_CO1_30_dc,function(x) sum(x > 0) > 1, TRUE)
write.table(data.frame("sum.taxonomy"=rownames(otu_table(physeqCIB_CO1_30_dct)),otu_table(physeqCIB_CO1_30_dct)),"physeqCIB_CO1_30_dct_0.1.txt", row.names=FALSE, sep="\t",quote = FALSE)

#decontaminate FITS_30
table(sample_data(physeqCIB_FITS_30)[,"Sample_or_Control"])


df <- as.data.frame(sample_data(physeqCIB_FITS_30))

library(ggplot2)

library(ggpubr)

library(decontam)

df$LibrarySize <- sample_sums(physeqCIB_FITS_30)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
p<-ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
ggsave("physeqCIB_FITS_30_libsizes_predecontam.pdf",p, width=6, height=5, units="cm", scale=3, limitsize = FALSE)
 
 #try 0.1 decontam
sample_data(physeqCIB_FITS_30)$Sample_or_Control == "Control"
sample_data(physeqCIB_FITS_30)$is.neg <- sample_data(physeqCIB_FITS_30)$Sample_or_Control == "Control"
contamdf.prev01 <- isContaminant(physeqCIB_FITS_30, method="prevalence", neg="is.neg", threshold=0.1)
head(which(contamdf.prev01$contaminant))
table(contamdf.prev01$contaminant)
physeqCIB_FITS_30_dc<-prune_taxa(!contamdf.prev01$contaminant, physeqCIB_FITS_30)
df <- as.data.frame(sample_data(physeqCIB_FITS_30_dc))
df$LibrarySize <- sample_sums(physeqCIB_FITS_30_dc)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
p<-ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
ggsave("physeqCIB_FITS_30_libsizes_postdecontam0.1.pdf",p, width=6, height=5, units="cm", scale=3, limitsize = FALSE)

physeqCIB_FITS_30_dct = filter_taxa(physeqCIB_FITS_30_dc,function(x) sum(x > 0) > 1, TRUE)

write.table(data.frame("sum.taxonomy"=rownames(otu_table(physeqCIB_FITS_30_dct)),otu_table(physeqCIB_FITS_30_dct)),"physeqCIB_FITS_30_dct_0.1.txt", row.names=FALSE, sep="\t",quote = FALSE)





#remove blanks from decontaminated phyloseq objects
#remove blanks, make min 5

physeqCIB_CO1_30_dct_noblanks = subset_samples(physeqCIB_CO1_30_dct, Sample_or_Control == "Sample")
physeqCIB_CO1_5_dct_noblanks = subset_samples(physeqCIB_CO1_5_dct, Sample_or_Control == "Sample")
physeqCIB_FITS_30_dct_noblanks = subset_samples(physeqCIB_FITS_30_dct, Sample_or_Control == "Sample")
physeqCIB_FITS_5_dct_noblanks = subset_samples(physeqCIB_FITS_5_dct, Sample_or_Control == "Sample")

physeqCIB_CO1_30_dct_noblanks_min5 = prune_taxa(taxa_sums(physeqCIB_CO1_30_dct_noblanks) > 4, physeqCIB_CO1_30_dct_noblanks) 
write.table(data.frame("sum.taxonomy"=rownames(otu_table(physeqCIB_CO1_30_dct_noblanks_min5)),otu_table(physeqCIB_CO1_30_dct_noblanks_min5)),"physeqCIB_CO1_30_dct_noblanks_min5.txt", row.names=FALSE, sep="\t",quote = FALSE)
saveRDS(physeqCIB_CO1_30_dct_noblanks_min5, "physeqCIB_CO1_30_dct_noblanks_min5.RDS")

physeqCIB_CO1_5_dct_noblanks_min5 = prune_taxa(taxa_sums(physeqCIB_CO1_5_dct_noblanks) > 4, physeqCIB_CO1_5_dct_noblanks) 
write.table(data.frame("sum.taxonomy"=rownames(otu_table(physeqCIB_CO1_5_dct_noblanks_min5)),otu_table(physeqCIB_CO1_5_dct_noblanks_min5)),"physeqCIB_CO1_5_dct_noblanks_min5.txt", row.names=FALSE, sep="\t",quote = FALSE)
saveRDS(physeqCIB_CO1_5_dct_noblanks_min5, "physeqCIB_CO1_5_dct_noblanks_min5.RDS")

physeqCIB_FITS_5_dct_noblanks_min5 = prune_taxa(taxa_sums(physeqCIB_FITS_5_dct_noblanks) > 4, physeqCIB_FITS_5_dct_noblanks) 
write.table(data.frame("sum.taxonomy"=rownames(otu_table(physeqCIB_FITS_5_dct_noblanks_min5)),otu_table(physeqCIB_FITS_5_dct_noblanks_min5)),"physeqCIB_FITS_5_dct_noblanks_min5.txt", row.names=FALSE, sep="\t",quote = FALSE)
saveRDS(physeqCIB_FITS_5_dct_noblanks_min5, "physeqCIB_FITS_5_dct_noblanks_min5.RDS")

physeqCIB_FITS_30_dct_noblanks_min5 = prune_taxa(taxa_sums(physeqCIB_FITS_30_dct_noblanks) > 4, physeqCIB_FITS_30_dct_noblanks) 
write.table(data.frame("sum.taxonomy"=rownames(otu_table(physeqCIB_FITS_30_dct_noblanks_min5)),otu_table(physeqCIB_FITS_30_dct_noblanks_min5)),"physeqCIB_FITS_30_dct_noblanks_min5.txt", row.names=FALSE, sep="\t",quote = FALSE)
saveRDS(physeqCIB_FITS_30_dct_noblanks_min5, "physeqCIB_FITS_30_dct_noblanks_min5.RDS")

physeqCIB_CO1_5_dct_noblanks_min5_Arthropoda<-subset_taxa(physeqCIB_CO1_5_dct_noblanks_min5, Domain=="Arthropoda")

physeqCIB_CO1_30_dct_noblanks_min5_Arthropoda<-subset_taxa(physeqCIB_CO1_30_dct_noblanks_min5, Domain=="Arthropoda")


### rarefaction

##how to subset by a metadata column, in this case location
physeqCIB_CO1_5_dct_noblanks_min5_Arthropoda_CDFA = subset_samples(physeqCIB_CO1_5_dct_noblanks_min5_Arthropoda, Site == "CDFA")
physeqCIB_CO1_5_dct_noblanks_min5_Arthropoda_SDNHM = subset_samples(physeqCIB_CO1_5_dct_noblanks_min5_Arthropoda, Site == "SD Natural History Museum")

physeqCIB_CO1_5_dct_noblanks_min5_Arthropoda_SSI = subset_samples(physeqCIB_CO1_5_dct_noblanks_min5_Arthropoda, Site == "Sierra Streams Institute ")

write.table(data.frame("sum.taxonomy"=rownames(otu_table(physeqCIB_CO1_5_dct_noblanks_min5_Arthropoda_CDFA)),otu_table(physeqCIB_CO1_5_dct_noblanks_min5_Arthropoda_CDFA)),"physeqCIB_CO1_5_dct_noblanks_min5_Arthropoda_CDFA.txt", row.names=FALSE, sep="\t",quote = FALSE)

write.table(data.frame("sum.taxonomy"=rownames(otu_table(physeqCIB_CO1_5_dct_noblanks_min5_Arthropoda_SDNHM)),otu_table(physeqCIB_CO1_5_dct_noblanks_min5_Arthropoda_SDNHM)),"physeqCIB_CO1_5_dct_noblanks_min5_Arthropoda_SDNHM.txt", row.names=FALSE, sep="\t",quote = FALSE)

write.table(data.frame("sum.taxonomy"=rownames(otu_table(physeqCIB_CO1_5_dct_noblanks_min5_Arthropoda_SSI)),otu_table(physeqCIB_CO1_5_dct_noblanks_min5_Arthropoda_SSI)),"physeqCIB_CO1_5_dct_noblanks_min5_Arthropoda_SSI.txt", row.names=FALSE, sep="\t",quote = FALSE)

