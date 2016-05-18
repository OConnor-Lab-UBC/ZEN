###################################################################################
#                                                                                ##
# ZEN 2014 data analysis (JED) V1                                                ##
# Data are current as of 17 May 2016                                             ##
# Emmett Duffy (duffye@si.edu)                                                   ##  
# Last updated 2016-05-17                                                        ##
#                                                                                ##
###################################################################################

# TO DO:
#   on main data spreadsheet: harmonize variable names (column heads) in data tabs
#   with associated variable names on metadata tab. The  remove the fossil 'tab' 
#   headers in metadata tab. 

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# METADATA                                                                        #
# LOAD PACKAGES                                                                   #
# READ IN AND PREPARE DATA                                                        #
# SUMMARIZE MESOGRAZER RICHNESS AND OUTPUT SPECIES X PLOT MATRIX                  #
# ASSEMBLE DATA FRAMES INTO A MASTER DATA SET                                     #
#                                                                                 #
###################################################################################

###################################################################################
# METADATA                                                                        #
###################################################################################

# Source data: ZEN_2014_Site&PlotData_2016_05_17_Released.xlsx                    
                                                                                
# This script analyzes data from the ZEN 2014 field sampling program in 2014. 
# The script started at ZEN workshop, Dvais, CA, May 2016.

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

# Load packages:
library(plyr)         
library(reshape)
library(reshape2)
library(plotrix)
library(ggplot2)
library(psych)
library(stringr)
library(nlme)
library(piecewiseSEM)
library(randomForest)
library(car)

###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

# Read in plot data: 
ZEN_2014_plot_data <- read.csv(file = "ZEN_2014_Site&PlotData_2016_05_17_Released.csv", header = TRUE)
names(ZEN_2014_plot_data)

# Remove redundant variables
ZEN_2014_plot_data <- ZEN_2014_plot_data[,!colnames(ZEN_2014_plot_data) %in% c("X")]
str(ZEN_2014_plot_data)

# Rename misspelled or confusing variables
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Mean.Sheath.Width.cm."] <- "Zostera.sheath.width"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Mean.Shealth.Length.cm."] <- "Zostera.sheath.length"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Mean.Longest.Leaft.Length.cm."] <- "Zostera.longest.leaf.length"

names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Pdn.Amphipod"] <- "amphipod.survival.24hr"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Pdn.Caprellid"] <- "caprellid.survival.24hr"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Pdn.Gastropod"] <- "gastropod.survival.24hr"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Pdn.Isopod"] <- "isopod.survival.24hr"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Pdn.Lettuce"] <- "lettuce.survival.24hr"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Pdn.Squid"] <- "squid.survival.24hr"

names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Mean.Above.Zmarina.g"] <- "Zostera.aboveground.mean.mass"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Mean.Above.Total.macrophytes.core.g"] <- "total.macrophytes.aboveground.mean.mass"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Mean.Below.Zmarina.g"] <- "Zostera.belowground.mean.mass"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Shoots.Zmarina.per.m2"] <- "Zostera.shoots.per.m2.core"

names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Epibiota.Periphyton"] <- "periphyton.mass.per.g.zostera"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Total.Abund.Mesograzers"] <- "mesograzer.abund.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Total.Abund.Molluscs"] <- "mollusc.abund.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Total.Abund.Crustaceans"] <- "crustacean.abund.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Total.Abund.Amphipods"] <- "amphipod.abund.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Abund.Gammarids"] <- "gammarid.abund.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Abund.Caprellids"] <- "caprellid.abund.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Abund.Gastropods"] <- "gastropod.abund.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Abund.Isopods"] <- "isopod.abund.per.g.plant"

names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Total.Biomass.Mesograzers"] <- "mesograzer.mass.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Total.Biomass.Molluscs"] <- "mollusc.mass.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Total.Biomass.Crustaceans"] <- "crustacean.mass.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Total.Biomass.Amphipods"] <- "amphipod.mass.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Biomass.Gammarids"] <- "gammarid.mass.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Biomass.Caprellids"] <- "caprellid.mass.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Biomass.Gastropods"] <- "gastropod.mass.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Biomass.Isopods"] <- "isopod.mass.per.g.plant"

names(ZEN_2014_plot_data)

# Read in plot-level (mesograzer) data to calculate site-level grazer richness.
epifauna <- read.csv(file = "ZENnverts.v8.csv", header = TRUE)

names(epifauna)

epifauna_clean <- epifauna[,!colnames(epifauna) %in% c("Data.Entered.by", "Measured.size",
  "Measured.Size","Mesured.Size","Mesured.Size.1","Mesured.Size.2","Mesured.Size.3",
  "Mesured.Size.4","Mesured.Size.5","Mesured.Size.6","Mesured.Size.7","Mesured.Size.8",
  "individuals.SizeCount","DW..g.","Vegetative.or.Flower","Shoot.ID", "Total.Abundance.to.1.mm",
  "Notes")]

names(epifauna_clean)

# Remove Sampling.Time = 2
epifauna_clean <- droplevels(subset(epifauna_clean, Sampling.Time == "1"))


###################################################################################
# SUMMARIZE MESOGRAZER RICHNESS AND OUTPUT SPECIES X PLOT MATRIX                  #
###################################################################################

# Create a dataframe with biomass of MESOGRAZER SPECIES as columns
# Melt the dataframe so that each SPECIES is a row:

names(epifauna_clean)

# Remove any non-grazers
grazers <- droplevels(subset(epifauna_clean, Type == "Mesograzer"))
levels(grazers$Type) # good

mesograzers <- melt(grazers, id.var = c("Site", "Plot.ID", "unique.ID","Species"), measure.var = c("Total.Abundance"))
# Cast the molten data frame so that each SPECIES is a column across the top:
mesograzers <- dcast(mesograzers, Site + Plot.ID + unique.ID ~ Species, sum)
str(mesograzers)
names(mesograzers)

# Export data frame as csv file: grazer species x plot ID matrix:
write.csv(mesograzers, "ZEN_2014_mesograzers_2016-05-17.csv", row.names = F)


# Create variable: PLOT-level grazer richness
mesograzers$richness.plot = 
  rowSums(mesograzers[, !colnames(mesograzers) %in% c("Site", "Plot.ID", "unique.ID")] > 0, na.rm = T)

# Create variable: SITE-level grazer richness 
richness.site <- ddply(mesograzers, c("Site"), function(x) {
  # Sum abundances 
  sum.x = colSums(x[, !colnames(x) %in% c("Site", "Plot.ID", "unique.ID")], na.rm = T)
  # Get richness
  sum(sum.x > 0, na.rm = T)
})

# change name of newly created variable 
colnames(richness.site)[2] <- c("richness.site")

# Create a dataframe with the richness values by site
temp <- mesograzers[c("Site", "Plot.ID", "unique.ID", "richness.plot")]
grazer.richness <- merge(temp, richness.site, by = c("Site"))
names(grazer.richness)


###################################################################################
# ASSEMBLE DATA FRAMES INTO A MASTER DATA SET                                     #
###################################################################################

# Ensure all different variables have unique, consistent, and unambiguous names

names(grazer.richness)
colnames(grazer.richness)[4:5] <- c("grazer.richness.plot", "grazer.richness.site") 

# Bind the data frames into a master data frame
ZEN_2014_master_data = ZEN_2014_plot_data
ZEN_2014_master_data$grazer.richness.plot <- grazer.richness$grazer.richness.plot[match(ZEN_2014_master_data$Site, grazer.richness$Site)]
ZEN_2014_master_data$grazer.richness.site <- grazer.richness$grazer.richness.site[match(ZEN_2014_master_data$Site, grazer.richness$Site)]

names(ZEN_2014_master_data)





