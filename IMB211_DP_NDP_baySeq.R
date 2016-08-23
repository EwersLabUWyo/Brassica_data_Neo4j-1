# IMB211_DP_NDP_baySeq.R
# R version 3.2.2 (2015-08-14)
# August 20, 2016. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Brassica IMB211 DP vs NDP RNA-seq data analysis with baySeq

# This RNA-seq data analysis compares DP (Dense Planting) and 
# NP (Not Dense Planting) in IMB211 using the baySeq package. 

# Sources: collaborator Julin Maloof 
# http://jnmaloof.github.io/BIS180L_web/2016/05/19/RNAseq-edgeR/
# baySeq author: Thomas J. Hardcastle
# https://www.bioconductor.org/packages/devel/bioc/vignettes/baySeq/inst/doc/baySeq.pdf

library(baySeq)
library(parallel)

#------------------------------------------------------------------------------

if(require("parallel")) cl <- makeCluster(8) else cl <- NULL

# Create replicate structure, used to estimate the prior distributions.
replicates <- c(rep("DP", 12), rep("NDP", 12))

# Define groups of no differential expression and differential expression.
# Since we expect there to be some genes expressed in all conditions at a   
# comparable level, the NDE group will cover all 24 columns. The DE group will
# be split into DP and NDP columns.
groups <- list(NDE = rep(1, 24), DE = c(rep(1,12), rep(2, 12)))

# Read in IMB211_DP_NDP_RNAseq csv.
IMB211_DP_NDP_RNAseq <- read.csv(file = "IMB211_DP_NDP_RNAseq", row.names = 1)
IMB211_DP_NDP_RNAseq <- as.matrix(IMB211_DP_NDP_RNAseq)

# Combine count data and groups into a countData object.
CD <- new("countData", data = IMB211_DP_NDP_RNAseq, 
          replicates = replicates, groups = groups)

# Allow library sizes to be inferred from data. 
libsizes(CD) <- getLibsizes(CD)

# Add annotation details into countData object.
CD@annotation <- data.frame(rownames(IMB211_DP_NDP_RNAseq))

# Estimate priors by bootstrapping from data. 
CD <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = cl)

# Get posterior likelihoods.
CD <- getLikelihoods(CD, cl = cl, bootStraps = 3, verbose = F)

CD@posteriors[1:10, ]





