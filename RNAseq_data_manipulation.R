# RNAseq_IMB211_R500.R
# R version 3.2.2 (2015-08-14)
# August 19, 2016. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Brassica RNA-seq data manipulation to import to Neo4j

# This "data manipulation" is meant to separate DP (Dense Planting) and 
# NP (Not Dense Planting) and organize the data into csv files that
# will aid its importation into Neo4j. Furthermore, csv's to be used for baySeq
# analysis will also be created. 

# Sources: collaborator Julin Maloof 
# http://jnmaloof.github.io/BIS180L_web/2016/05/19/RNAseq-edgeR/

library(dplyr)
library(data.table)

#------------------------------------------------------------------------------

# Read in greenhouse file "GH_merged_v1.5_mapping.tsv".
RNAseqcounts <- read.table(file = "GH_merged_v1.5_mapping.tsv", header = T)

# Remove first row "*" since these are reads that didn't map to a gene.
RNAseqcounts <- RNAseqcounts[-1, ]

# Replace NA's with 0. 
RNAseqcounts[is.na(RNAseqcounts)] <- 0

# Save gene row as vector and set genes to rownames in RNAseqcounts data frame.
genes <- RNAseqcounts[, 1]
RNAseqcounts <- RNAseqcounts[, -1]
rownames(RNAseqcounts) <- genes

# Drop endings .1_matched.merged.fq.bam from column names.
colnames(RNAseqcounts) <- sub(".1_matched.merged.fq.bam", "", 
                              colnames(RNAseqcounts))
genes <- rownames(RNAseqcounts)

# Drop FIELD data and leave only DP and NDP.
RNAseqcounts <- select(RNAseqcounts, -contains("FIELD"))

# Separate data into IMB211 and R500.
IMB211 <- select(RNAseqcounts, contains("IMB211"))
R500 <- select(RNAseqcounts, contains("R500"))
# Drop the IMB211 and R500 from the column names
colnames(IMB211) <- sub("IMB211_", "", colnames(IMB211))
colnames(R500) <- sub("R500_", "", colnames(R500))

################# BAYSEQ ANALYSIS#######################################
# For RNA seq data analysis, we'll only analyse genes with more than ten 
# reads in at least three samples. 
IMB211 <- IMB211[rowSums(IMB211 > 10) >= 3, ]
R500 <- R500[rowSums(R500 > 10) >= 3, ]
# Will reduce to only 27,023 genes for IMB211 and 26,536 for R500.

# Place data.frames into csv's for baySeq analysis.
write.csv(IMB211, file = "IMB211_DP_NDP_RNAseq.csv")
write.csv(R500, file = "R500_DP_NDP_RNAseq.csv")

################# NEO4J IMPORTATION######################################
# Turn data.frame into a data.table to melt. 
IMB211 <- data.table(IMB211, keep.rownames = T)

# Melt data.table to separate into treatment and tissue type.
IMB211 <- melt.data.table(IMB211, id = 1:25, measure.vars = 
                            patterns("1", "2", "3"), value.factor = T)

# Delete the columns that aren't the genes, variable, and value1-value3.
IMB211 <- IMB211[, grep("DP", colnames(IMB211)) := NULL]

# Rename the rn column to be "Gene," the value columns to be "Count1", etc.
colnames(IMB211) <- sub("rn", "Gene", colnames(IMB211))
colnames(IMB211) <- sub("value1", "Count1", colnames(IMB211), fixed = T)
colnames(IMB211) <- sub("value2", "Count2", colnames(IMB211), fixed = T)
colnames(IMB211) <- sub("value3", "Count3", colnames(IMB211), fixed = T)
# NOTE: This code is easier to read than the functions that would need to be
# written to do these all consecutively. 

# Variables 1-8 represent each treatment and tissue type. Add two columns that
# Separate these into a treatment column named "Treatment" and a tissue column
# named "Tissue"  by adding the appropriate columns for easy integration into 
# Neo4j. 
IMB211 <- IMB211[, Treatment := rep(c("DP", "NDP"), c(108092,108092))]
IMB211 <- IMB211[, Tissue := rep(c("INTERNODE", "LEAF", "PETIOLE", "SILIQUE"),
                                    c(27023, 27023, 27023, 27023))]
# TO DO: (M. Lai) Re-do in set syntax.

# Delete variable column.
IMB211 <- IMB211[, variable := NULL]

# Melt again with counts. 
IMB211 <- melt.data.table(IMB211, id = c("Gene", "Treatment", "Tissue"), 
                       measure.vars = patterns("Count"), value.factor = T)

# Rename columns value 1 and variable to be Count and Replicate.
colnames(IMB211) <- sub("value1", "Count", colnames(IMB211))
colnames(IMB211) <- sub("variable", "Replicate", colnames(IMB211))

# Change Count1 to 1, etc.
IMB211 <- IMB211[Replicate=="Count1", Replicate := "1"]
IMB211 <- IMB211[Replicate=="Count2", Replicate := "2"]
IMB211 <- IMB211[Replicate=="Count3", Replicate := "3"]

# Write csv for IMB211.
write.csv(IMB211, file = "IMB211_DP_NDP_RNAseq_Neo4j.csv")

# Now do the same with the R500 genotype:
# Turn data.frame into a data.table to melt. 
R500 <- data.table(R500, keep.rownames = T)

# Melt data.table to separate into treatment and tissue type.
R500 <- melt.data.table(R500, id = 1:25, measure.vars = 
                            patterns("1", "2", "3"), value.factor = T)

# Delete the columns that aren't the genes, variable, and value1-value3.
R500 <- R500[, grep("DP", colnames(R500)) := NULL]

# Rename the rn column to be "Gene," the value columns to be "Count1", etc.
colnames(R500) <- sub("rn", "Gene", colnames(R500))
colnames(R500) <- sub("value1", "Count1", colnames(R500), fixed = T)
colnames(R500) <- sub("value2", "Count2", colnames(R500), fixed = T)
colnames(R500) <- sub("value3", "Count3", colnames(R500), fixed = T)
# NOTE: This code is easier to read than the functions that would need to be
# written to do these all consecutively. 

# Variables 1-8 represent each treatment and tissue type. Add two columns that
# Separate these into a treatment column named "Treatment" and a tissue column
# named "Tissue"  by adding the appropriate columns for easy integration into 
# Neo4j. 
R500 <- R500[, Treatment := rep(c("DP", "NDP"), c(106144,106144))]
R500 <- R500[, Tissue := rep(c("INTERNODE", "LEAF", "PETIOLE", "SILIQUE"),
                                 c(26536, 26536, 26536, 26536))]
# TO DO: (M. Lai) Re-do in set syntax.

# Delete variable column.
R500 <- R500[, variable := NULL]

# Melt again with counts. 
R500 <- melt.data.table(R500, id = c("Gene", "Treatment", "Tissue"), 
                          measure.vars = patterns("Count"), value.factor = T)

# Rename columns value 1 and variable to be Count and Replicate.
colnames(R500) <- sub("value1", "Count", colnames(R500))
colnames(R500) <- sub("variable", "Replicate", colnames(R500))

# Change Count1 to 1, etc.
R500 <- R500[Replicate == "Count1", Replicate := "1"]
R500 <- R500[Replicate == "Count2", Replicate := "2"]
R500 <- R500[Replicate == "Count3", Replicate := "3"]

# Write csv for R500.
write.csv(R500, file = "R500_DP_NDP_RNAseq_Neo4j.csv")










