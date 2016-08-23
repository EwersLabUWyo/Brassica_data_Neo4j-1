# 4Mall_data_edits.R
# R version 3.2.2 (2015-08-14)
# August 23, 2016. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# 4Mall* data consistency fixes.

# This "data clean-up" is meant to add consistency to our data files. 
# eg. '4Mall_Dried_Biomass_dec_2013.csv' uses "W" and "D" for treatment
# conditions whereas other files use "WW" (for "well-watered") and "Dry."
# Change "W's" and "D's" were to be consistent with the other files. The
# "Line" column should also be changed to R500 for consistency. Similar 
# edits will be made to other "4Mall" files. 

library(data.table)

#------------------------------------------------------------------------------

# Read in 4Mall_Dried_Biomass_dec_2013 file
DriedBiomass <- read.csv(file = "4Mall_Dried_Biomass_dec_2013.csv") 

# Change column names to get rid of periods. 
colnames(DriedBiomass) <- sub("Date.after.Oven", "DateAfterOven", colnames(DriedBiomass))
colnames(DriedBiomass) <- sub("Dry.Weight.Shoots", "DryWeightShoots", colnames(DriedBiomass))

# Change to data.table for following edits. 
DriedBiomass <- data.table(DriedBiomass)

# Change NA's in Line column to R500.
DriedBiomass <- DriedBiomass[, Line := as.character(Line)]
DriedBiomass <- DriedBiomass[, Line := "R500"]

# Change "W" in Treatment to "WW" and "D" to "Dry".
DriedBiomass <- DriedBiomass[Treatment == "W", Treatment := "WW"]
DriedBiomass <- DriedBiomass[Treatment == "D", Treatment := "Dry"]

write.csv(DriedBiomass, file = "4Mall_Dried_Biomass_dec_2013_edit.csv")

###############################################################################
# Read in 4Mall_NSC_Starchdec2013 file
Starch <- read.csv(file = "4Mall_NSC_Starchdec2013.csv")

# Turn into data.table.
Starch <- data.table(Starch)

# Remove rows AVG, SD, and SE
Starch <- Starch[ Tissue != 'AVG']
Starch <- Starch[ Tissue != 'SD']
Starch <- Starch[ Tissue != 'SE']

# Remove column X and X.1.
Starch <- Starch[, grep("X", colnames(Starch)) := NULL]

# Rename Time.Point Timepoint
colnames(Starch) <- sub("Time.Point", "Timepoint", colnames(Starch))





















