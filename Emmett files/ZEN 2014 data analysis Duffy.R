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

# Read in summary SUBSITE data: 
SubsiteSummary <- read.csv(file = "ZEN_2014_Site&PlotData_2016_05_17_Released_subsite_summary.csv", header = TRUE)
names(SubsiteSummary)
str(SubsiteSummary)
# NOTE: Variables of interest from this file are primarily site-level metadata, plus 
# genetic data (which are site-level metrics). The many site-level means are recalucklated below.  

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
# CREATE MESOGRAZER RICHNESS SPECIES X PLOT MATRIX                                #
###################################################################################

# PLOT level:

# Create a dataframe with biomass of MESOGRAZER SPECIES as columns
# Melt the dataframe so that each PLOT is a row:

names(epifauna_clean)

# Remove any non-grazers
grazers <- droplevels(subset(epifauna_clean, Type == "Mesograzer"))
levels(grazers$Type) # good

mesograzers_plot <- melt(grazers, id.var = c("Site", "Plot.ID", "unique.ID","Species"), measure.var = c("Total.Abundance"))
# Cast the molten data frame so that each SPECIES is a column across the top:
mesograzers_plot <- dcast(mesograzers_plot, Site + Plot.ID + unique.ID ~ Species, sum)
str(mesograzers_plot)
names(mesograzers_plot)

# Export data frame as csv file: grazer species x plot ID matrix:
write.csv(mesograzers_plot, "ZEN_2014_mesograzers_PLOT-LEVEL_2016-05-17.csv", row.names = F)


# SITE level:

# Create a dataframe with biomass of MESOGRAZER SPECIES as columns
# Melt the dataframe so that each SITE is a row:

mesograzers_site <- melt(grazers, id.var = c("Site","Species"), measure.var = c("Total.Abundance"))
# Cast the molten data frame so that each SPECIES is a column across the top:
mesograzers_site <- dcast(mesograzers_site, Site ~ Species, sum)
str(mesograzers_site)
names(mesograzers_site)

# Export data frame as csv file: grazer species x plot ID matrix:
write.csv(mesograzers_site, "ZEN_2014_mesograzers_SITE-LEVEL_2016-05-17.csv", row.names = F)


###################################################################################
# SUMMARIZE MESOGRAZER RICHNESS BY PLOT AND SITE                                  #
###################################################################################

# Create variable: PLOT-level grazer richness
mesograzers_plot$richness.plot = 
  rowSums(mesograzers_plot[, !colnames(mesograzers_plot) %in% c("Site", "Plot.ID", "unique.ID")] > 0, na.rm = T)

# Create variable: SITE-level grazer richness 
richness.site <- ddply(mesograzers_plot, c("Site"), function(x) {
  # Sum abundances 
  sum.x = colSums(x[, !colnames(x) %in% c("Site", "Plot.ID", "unique.ID")], na.rm = T)
  # Get richness
  sum(sum.x > 0, na.rm = T)
})

# change name of newly created variable 
colnames(richness.site)[2] <- c("richness.site")

# Create a dataframe with the richness values by site
temp <- mesograzers_plot[c("Site", "Plot.ID", "unique.ID", "richness.plot")]
grazer.richness <- merge(temp, richness.site, by = c("Site"))
names(grazer.richness)


###################################################################################
# ASSEMBLE DATA FRAMES INTO A MASTER DATA SET                                     #
###################################################################################

names(grazer.richness)
colnames(grazer.richness)[4:5] <- c("grazer.richness.plot", "grazer.richness.site") 

# Bind the data frames into a master data frame
ZEN_2014_master_data = ZEN_2014_plot_data
ZEN_2014_master_data$grazer.richness.plot <- grazer.richness$grazer.richness.plot[match(ZEN_2014_master_data$Site, grazer.richness$Site)]
ZEN_2014_master_data$grazer.richness.site <- grazer.richness$grazer.richness.site[match(ZEN_2014_master_data$Site, grazer.richness$Site)]

names(ZEN_2014_master_data)

# Select those environmental variables we want to include:
names(SubsiteSummary)

envdata <- SubsiteSummary[c("Site", "Subsite", "Subsite.Name", "Ocean", "Coast", "Basin", "Latitude", 
                            "Longitude", "Sampling.Time", "Month", "Date.Collected", "Temperature.C", "Salinity.ppt",
                            "Day.length.hours", "Depth.Categorical", "Depth.m", "Perc.Silt", "Perc.Sand", "Perc.Gravel", 
                            "GenotypicRichness", "AllelicRichness", "Inbreeding")]

ZEN_2014_master_data <- merge(ZEN_2014_master_data, envdata, by = c("Site"))
names(ZEN_2014_master_data)
nrow(ZEN_2014_master_data) # 1000


###################################################################################
# EXPLORE DISTRIBUTIONS OF VARIABLES (PLOT LEVEL)                                 #
###################################################################################

# Examine frequency distribution of sites by environmental factor
# par(mfrow = c(1,1))
par(mfrow = c(2,4))
hist(ZEN_2014_master_data$Latitude, col = "cyan", main = "Surveys by latitude")    
hist(ZEN_2014_master_data$Longitude, col = "cyan", main = "Surveys by longitude")    
hist(ZEN_2014_master_data$Temperature.C, col = "cyan", main = "Surveys by temperature")    
hist(ZEN_2014_master_data$Salinity.ppt, col = "cyan", main = "Surveys by salinity")    

hist(ZEN_2014_master_data$Zostera.aboveground.mean.mass, col = "cyan")    
hist(ZEN_2014_master_data$total.macrophytes.aboveground.mean.mass, col = "cyan")    
hist(ZEN_2014_master_data$periphyton.mass.per.g.zostera, col = "cyan")    
hist(ZEN_2014_master_data$crustacean.abund.per.g.plant, col = "cyan")    

hist(ZEN_2014_master_data$gastropod.abund.per.g.plant, col = "cyan")    
hist(ZEN_2014_master_data$grazer.richness.plot, col = "cyan")    
hist(ZEN_2014_master_data$amphipod.survival, col = "cyan")    


###################################################################################
# LOG TRANSFORMS                                                                  #
###################################################################################

names(ZEN_2014_master_data)

# NOTE: For many variables I add a constant roughly equal to the smallest value recorded 
ZEN_2014_master_data$log10.Zostera.AG.mass <- log10(ZEN_2014_master_data$Zostera.aboveground.mean.mass + 1) 
ZEN_2014_master_data$log10.macrophytes.total.AG.mass.core <- log10(ZEN_2014_master_data$total.macrophytes.aboveground.mean.mass + 1) 
ZEN_2014_master_data$log10.Zostera.BG.mass <- log10(ZEN_2014_master_data$Zostera.belowground.mean.mass + 1) 
ZEN_2014_master_data$log10.Zostera.shoots.core <- log10(ZEN_2014_master_data$Zostera.shoots.per.m2.core + 1) 
ZEN_2014_master_data$log10.Zostera.sheath.width <- log10(ZEN_2014_master_data$Zostera.sheath.width) 
ZEN_2014_master_data$log10.Zostera.sheath.length <- log10(ZEN_2014_master_data$Zostera.sheath.length) 
ZEN_2014_master_data$log10.Zostera.longest.leaf.length <- log10(ZEN_2014_master_data$Zostera.longest.leaf.length) 
ZEN_2014_master_data$log10.periphyton.mass.per.g.zostera <- log10(ZEN_2014_master_data$periphyton.mass.per.g.zostera + 0.0001) 

ZEN_2014_master_data$log10.mesograzer.abund.per.g.plant <- log10(ZEN_2014_master_data$mesograzer.abund.per.g.plant + 0.1) 
ZEN_2014_master_data$log10.crustacean.abund.per.g.plant <- log10(ZEN_2014_master_data$crustacean.abund.per.g.plant + 0.1) 
ZEN_2014_master_data$log10.amphipod.abund.per.g.plant <- log10(ZEN_2014_master_data$amphipod.abund.per.g.plant + 0.1) 
ZEN_2014_master_data$log10.gammarid.abund.per.g.plant <- log10(ZEN_2014_master_data$gammarid.abund.per.g.plant + 0.1) 
ZEN_2014_master_data$log10.gastropod.abund.per.g.plant <- log10(ZEN_2014_master_data$gastropod.abund.per.g.plant + 0.1) 

ZEN_2014_master_data$log10.mesograzer.mass.per.g.plant <- log10(ZEN_2014_master_data$mesograzer.mass.per.g.plant + 0.1) 
ZEN_2014_master_data$log10.crustacean.mass.per.g.plant <- log10(ZEN_2014_master_data$crustacean.mass.per.g.plant + 0.1) 
ZEN_2014_master_data$log10.amphipod.mass.per.g.plant <- log10(ZEN_2014_master_data$amphipod.mass.per.g.plant + 0.1) 
ZEN_2014_master_data$log10.gammarid.mass.per.g.plant <- log10(ZEN_2014_master_data$gammarid.mass.per.g.plant + 0.1) 
ZEN_2014_master_data$log10.gastropod.mass.per.g.plant <- log10(ZEN_2014_master_data$gastropod.mass.per.g.plant + 0.1) 

ZEN_2014_master_data$log10.grazer.richness.plot <- log10(ZEN_2014_master_data$grazer.richness.plot + 1) 
ZEN_2014_master_data$log10.grazer.richness.site <- log10(ZEN_2014_master_data$grazer.richness.site + 1) 

ZEN_2014_master_data$log10.Leaf.PercN <- log10(ZEN_2014_master_data$Leaf.PercN) 
ZEN_2014_master_data$log10.GenotypicRichness <- log10(ZEN_2014_master_data$GenotypicRichness) 
ZEN_2014_master_data$log10.AllelicRichness <- log10(ZEN_2014_master_data$AllelicRichness) 

# Obtain values for mean body mass per individual
ZEN_2014_master_data$crustacean.mean.body.mass <- ZEN_2014_master_data$crustacean.mass.per.g.plant / ZEN_2014_master_data$crustacean.abund.per.g.plant 
ZEN_2014_master_data$log10.crustacean.mean.body.mass <- log10(ZEN_2014_master_data$crustacean.mean.body.mass)

ZEN_2014_master_data$gastropod.mean.body.mass <- ZEN_2014_master_data$gastropod.mass.per.g.plant / ZEN_2014_master_data$gastropod.abund.per.g.plant 
ZEN_2014_master_data$log10.gastropod.mean.body.mass <- log10(ZEN_2014_master_data$gastropod.mean.body.mass)

ZEN_2014_master_data$mesograzer.mean.body.mass <- ZEN_2014_master_data$mesograzer.mass.per.g.plant / ZEN_2014_master_data$mesograzer.abund.per.g.plant 
ZEN_2014_master_data$log10.mesograzer.mean.body.mass <- log10(ZEN_2014_master_data$mesograzer.mean.body.mass)

ZEN_2014_master_data$basin <- ZEN_2014_master_data$Coast
levels(ZEN_2014_master_data$basin)[levels(ZEN_2014_master_data$basin)=='East.Atlantic'] <- 'East'
levels(ZEN_2014_master_data$basin)[levels(ZEN_2014_master_data$basin)=='East.Pacific'] <- 'East'
levels(ZEN_2014_master_data$basin)[levels(ZEN_2014_master_data$basin)=='West.Atlantic'] <- 'West'
levels(ZEN_2014_master_data$basin)[levels(ZEN_2014_master_data$basin)=='West.Pacific'] <- 'West'


# Change values of NaN to NA:
ZEN_2014_master_data[ZEN_2014_master_data == "NaN"] = NA 
ZEN_2014_master_data[ZEN_2014_master_data == "Inf"] = NA 
names(ZEN_2014_master_data)

# Export the data frame
write.csv(ZEN_2014_master_data, "ZEN_2014_master_data_2016-05-17.csv", row.names = F)


###################################################################################
# OBTAIN SITE MEANS                                                               #
###################################################################################






