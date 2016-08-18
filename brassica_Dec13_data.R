# brassica_Dec13_data.R
# R version 3.2.2 (2015-08-14)
# August 16, 2016. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Brassica data clean-up to import to Neo4j

# This "data clean-up" is meant to get rid of blank cells and data that's 
# not going to be imported to Neo4j.

library(data.table)

#------------------------------------------------------------------------------

# Read in file
brassica_Dec13_data <-read.csv(file.choose()) 

# Grab date/time from file
date <- as.character(brassica_Dec13_data[1, 1])

# Remove first 8 lines
output <- brassica_Dec13_data[9:dim(brassica_Dec13_data)[1], ] 

# Remove rows with "Remark="
output <- output[output$OPEN.6.1.4!="Remark=", ] 

# Turn first row into column names
colnames(output) <- as.character(unlist(output[1, ]))

# Delete row of column names
output <- output[-c(1,2), ]

# Keep only columns of data to import to Neo4j
output <- output[, c("Obs", "HHMMSS", "Photo", "Cond", "Fo", "Fm", "Fo'", 
                     "Fm'", "RedAbs", "BlueAbs", "VpdL", "Tair", "Tleaf", 
                     "TBlk", "CO2R", "CO2S", "H2OR", "H2OS", "RH_R", "RH_S",
                     "Flow", "PARabs")]

