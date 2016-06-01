# Load required libraries
# library(devtools)
# install_github("jslefche/divpart")
library(divpart)
library(ggplot2)
library(gridExtra)
library(iNEXT)
library(moments)
library(plyr)
library(reshape2)
library(vegan)

setwd("C:/Users/jslef/OneDrive/Documents/GitHub/ZEN/Diversity Partitioning")

# Create function to clean-up matrix
cleanup = function(mat, lvls) {
  
  # Convert NA to 0
  mat[is.na(mat)] = 0
  
  # Remove meta-data
  mat = mat[, -c(0:length(lvls))]
  
  # Convert to numeric matrix
  mat = data.matrix(mat)
  
  # Remove columns (species) where there are no data
  remove.cols = which(colSums(mat) == 0)
  
  if(length(remove.cols) > 0) mat = mat[, -remove.cols]
  
  # Remove rows (samples) where there are no data
  remove.rows <<- which(rowSums(mat) == 0)
  
  if(length(remove.rows) > 0) mat = mat[-remove.rows, ]
  
  return(mat) 
  
}

##########

# Load epifauna data
epifauna = read.table("ZENnverts.v8.txt", header = T)

# Subset epifauna to only include mesograzers
epifauna = subset(epifauna, Type == "Mesograzer")

# Load metadata 
metad = read.csv("zen_2014_metadata.csv")

# Bring metadata and epifaunal data together
zen = merge(
  metad[, c(1:6, 9:10, 14:23)],
  epifauna[, c(1:4, 7:9, 12, 15, 17, 20, 27, 31, 36)]# c(1:6, 11:12, 13, 17, 32)]
)

# Remove Sampling.Time > 1
zen = subset(zen, Sampling.Time == 1)

# Merge Subsite and subsite code
zen$Site = apply(zen, 1, function(x) paste0(x["Site.Code"], ".", x["Subsite"]))

# Drop Croatia from the analysis
# zen = subset(zen, Site.Code != "CR")

# Remove Portugal B from the analysis
# zen = subset(zen, Site != "PO.B")

# Creating nesting level
lvls = c("unique.ID", "Site", "Coast", "Ocean") 

# Check to see if unique entries in each level are listed in decreasing order
sapply(lvls, function(x) length(unique(zen[, x])))

##########

# Create sample-by-species matrix
D.mat = dcast(zen, as.formula(paste0(paste(lvls, collapse = " + "), " ~ Original.Species")), 
              fun.aggregate = sum, value.var = "Total.Abundance")

D.mat = cleanup(D.mat, lvls)

# Remove NAs
D.mat = D.mat[, -which(grepl("NA", colnames(D.mat)))]

##########

# Create sample-by-hierarchy matrix
G.mat = unique(zen[, lvls])

# Remove samples as above
G.mat = G.mat[-remove.rows, ]

# Convert to matrix
G.mat = as.matrix(G.mat)

# Check to see if nrow are the same
nrow(D.mat) == nrow(G.mat)

##########

# Partition diversity and store output in a list
zen.part.list = lapply(0:2, function(i) {
  
  part = divpart(D.mat, G.mat[, -1], q = i)
  
  write.csv(part, paste("Output/Diversity partitioning q = ", i, ".csv"))
  
} )

names(zen.part.list) = c("Richness", "Shannon", "Simpson")

# How do alpha, beta, and gamma diversities scale with increasing levels of the hierarchy?
part.plot.list = lapply(names(zen.part.list), function(i) {
  
  # Melt data.frame
  zen.part.melt = melt(zen.part.list[[i]], id.vars = c("level"), measure.vars = c("alpha", "gamma", "beta.add"))
  
  zen.part.melt$level = factor(zen.part.melt$level, levels = c("ID", "Site", "Coast", "Ocean"))
  
  # Re-level for better plotting
  levels(zen.part.melt$level) = c("Plot", "Site", "Coast", "Ocean")
  
  levels(zen.part.melt$variable) = c("Alpha", "Beta", "Gamma")
  
  # Plot results
  ggplot(zen.part.melt, aes(x = level, y = value, group = variable, shape = variable, col = variable)) +
    geom_point(col = "grey80", size = 1) +
    stat_summary(fun.data = "mean_cl_boot", size = 1.1, position = position_dodge(width = 1)) +
    geom_vline(xintercept = c(1.5, 2.5, 3.5), lwd = 1, col = "grey50") +
    scale_shape_manual(values = c(15:17), name = "") +
    scale_color_manual(values = c("grey50", "brown1", "black"), name = "") +
    labs(x = "", y = "", title = i) +
    theme_bw(base_size = 18) +
    theme(axis.ticks.x = element_blank())
  
} )

pdf(paste0("Figures/Diversity partitioning.PDF"),
    width = 6, height = 5)
bquiet = print(part.plot.list)
dev.off() 

# How does beta diversity scale with alpha and gamma diversities?
zen.part.plot.list = lapply(names(zen.part.list), function(i) {
  
  ggplot(zen.part.list[[i]], aes(x = alpha, y = beta.add)) +
    geom_point(size = 3) +
    geom_abline(slope = 1, lty = 2, lwd = 1) + 
    labs(x = "Alpha diversity", y = "Additive beta diversity", title = i) +
    theme_bw(base_size = 18)
  
} )

pdf(paste0("Figures/Diversity scaling.PDF"),
    width = 5, height = 5)
bquiet = print(zen.part.plot.list)
dev.off() 

##########

# Repeat for taxonomic dissimilarities

# Subset unique taxa
taxa = unique(zen[, c("Original.Species", "Genus", "Family", "Order", "Class")])

taxa = taxa[order(taxa$Original.Species), ]

taxa = taxa[!duplicated(taxa$Original.Species), ]

# Remove NAs
taxa = taxa[!is.na(taxa$Original.Species), ]

rownames(taxa) = taxa$Original.Species

# Subset to only include species in the abundance matrix
taxa = taxa[taxa$Original.Species %in% colnames(D.mat), ]

# Create taxonomic dissimilarity matrix
taxa.dist = taxa2dist(taxa)

# Partion diversity based on taxonomic dissimilarity
zen.taxo.part.list = lapply(0:2, function(i) {
  
  part = divpart(D.mat, G.mat[, -1], dissim = taxa.dist, q = i)
  
  write.csv(part, paste("Output/Diversity partitioning q = ", i, ".csv"))
  
} )

names(zen.taxo.part.list) = c("Richness", "Shannon", "Simpson")

# How do alpha, beta, and gamma diversities scale with increasing levels of the hierarchy?
taxo.part.plot.list = lapply(names(zen.taxo.part.list), function(i) {
  
  # Melt data.frame
  zen.part.melt = melt(zen.taxo.part.list[[i]], id.vars = c("level"), measure.vars = c("alpha", "gamma", "beta.add"))
  
  zen.part.melt$level = factor(zen.part.melt$level, levels = c("ID", "Site", "Coast", "Ocean"))
  
  # Re-level for better plotting
  levels(zen.part.melt$level) = c("Plot", "Site", "Coast", "Ocean")
  
  levels(zen.part.melt$variable) = c("Alpha", "Beta", "Gamma")
  
  # Plot results
  ggplot(zen.part.melt, aes(x = level, y = value, group = variable, shape = variable, col = variable)) +
    geom_point(col = "grey80", size = 1) +
    stat_summary(fun.data = "mean_cl_boot", size = 1.1, position = position_dodge(width = 1)) +
    geom_vline(xintercept = c(1.5, 2.5, 3.5), lwd = 1, col = "grey50") +
    scale_shape_manual(values = c(15:17), name = "") +
    scale_color_manual(values = c("grey50", "brown1", "black"), name = "") +
    labs(x = "", y = "", title = i) +
    theme_bw(base_size = 18) +
    theme(axis.ticks.x = element_blank())
  
} )

pdf(paste0("Figures/Taxonomic diversity partitioning.PDF"),
    width = 6, height = 5)
bquiet = print(taxo.part.plot.list)
dev.off() 

##########

# Non-metric Multidimensional Scaling

# Generate species, genus, family, class -by- site abundance matrix
site.mat.list = lapply(c("Species", "Genus", "Family", "Class"), function(i) {
  
  mat = dcast(subset(zen, Site.Code != "CR"), as.formula(paste0(paste(lvls[-1], collapse = " + "), " ~ ", i)), 
              fun.aggregate = sum, value.var = "Total.Abundance") 
  
  mat = cleanup(mat, lvls[-1])
  
  # Remove NA column
  mat = mat[, -which(grepl("NA", colnames(mat)))]
  
  return(mat)
  
} )

# Remove NA column
lapply(site.mat.list, dim)

names(site.mat.list) = c("Species", "Genus", "Family", "Class")

# Run NMDS and extract axes
nmds.list = lapply(names(site.mat.list), function(i) {
  
  # Run NMDS
  nmds = metaMDS(site.mat.list[[i]], trymax = 100)
  
  # Retrieve axes
  nmds. = nmds$points
  
  # Bind with metadata
  nmds. = data.frame(Site = unique(zen$Site)[-(5:6)], Stress = nmds$stress, nmds.)
  
  nmds. = cbind(nmds., zen[match(nmds.$Site, zen$Site), c("Site.Code", "Coast", "Ocean")])
  
  return(nmds.)
  
} )

names(nmds.list) = names(site.mat.list)

# Plot results
nmds.plot.list = llply(names(nmds.list), function(i) {
  
  # Use function chull() to derive hull vertices for different groupings
  hulls.df = ldply(unique(nmds.list[[i]]$Coast), function(j) {
    
    data = nmds.list[[i]][which(nmds.list[[i]]$Coast == j), c("MDS1", "MDS2")]
    
    data.frame(
      Coast = j,
      data[chull(data), ],
      row.names = NULL
    )
    
  } ) 
  
  # Re-level for better plotting
  levels(hulls.df$Coast) = gsub("\\.", " ", levels(hulls.df$Coast))
  
  # Plot results
  ggplot() +
    geom_polygon(data = hulls.df, aes(x = MDS1, y = MDS2, group = Coast, fill = Coast),
                 position = "identity",
                 alpha = 0.2, lwd = 0.8) +
    geom_text(data = nmds.list[[i]], aes(x = MDS1, y = MDS2, label = Site.Code)) +
    geom_text(data = nmds.list[[i]][1, ], aes(x = Inf, y = -Inf, label = paste0("Stress = ", round(Stress[1], 3))), 
              vjust = -1, hjust = 1.25) +
    labs(title = i) + 
    theme_bw(base_size = 18)
  
} )

pdf(paste0("Figures/NMDS by taxonomic level.PDF"),
    width = 7.5, height = 6)
bquiet = print(nmds.plot.list)
dev.off() 

# Output loadings
nmds.loadings = ldply(names(nmds.list), function(i) {
  
  # Get loadings for each species
  fit = envfit(nmds.list[[i]][, c("MDS1", "MDS2")], site.mat.list[[i]], perm = 1000)
  
  data.frame(
    level = i,
    name = rownames(fit$vectors$arrows),
    fit$vectors$arrows,
    p.value = fit$vectors$pvals
  )
  
} )

write.csv(nmds.loadings, "Output/NMDS loadings by taxonomic level.csv")

##########

# Rank abundance curves

# Create summary abundances by site
abund.site = ddply(zen, c("Site", "Original.Species"), summarize, abund = sum(Total.Abundance))

# Order species by their abundance within each site
abund.site.order = dlply(abund.site, "Site", function(x) {
  
  # Remove zero individuals
  x = subset(x, abund > 0)
  
  # Order species by abundance
  x = x[rev(order(x$abund)), ]

  x$Original.Species = factor(x$Original.Species, levels = x$Original.Species)
  
  return(x)
  
} )

# Plot results
abund.plot.list = llply(abund.site.order, function(x) 
  
  ggplot(x, aes(x = as.factor(as.numeric(Original.Species)), y = abund)) +
    geom_bar(stat = "identity", col = "black", fill = "white") + 
    labs(x = "Species", y = "Abundance", title = unique(x$Site)) 
  
)
  

# Output results
pdf(paste0("Figures/Rank abundance curves.PDF"),
    width = 30, height = 30)
bquiet = print(do.call(grid.arrange, abund.plot.list))
dev.off()

# Get kurtosis value for each site
kurtosis.df = ldply(abund.site.order, function(i) data.frame(kurtosis = kurtosis(i$abund)))
  
# write.csv(kurtosis.df, "Rank abundance kurtosis.csv")

##########

# Rarefaction curves

# Summarize species-by-site abundance matrix
D.mat.site = dcast(zen, Site ~ Original.Species, fun.aggregate = sum, value.var = "Total.Abundance")

rownames(D.mat.site) = unique(zen$Site)

# Clean up matrix
D.mat.site = cleanup(D.mat.site, 0)

# # Put each row in a list
# D.list.site = alply(D.mat.site, 1)
# 
# names(D.list.site) = rownames(D.mat.site)

# Run iNEXT
# zen.inext = llply(D.list.site, .progress = "text", iNEXT)

zen.inext = iNEXT(t(D.mat.site))

# Plot results and store in a list
# zen.inext.plot = ggiNEXT(zen.inext, facet.var = "site") + facet_grid(scales = "free")
# 
# # Output results 
# pdf(paste0("Figures/Rarefaction extrapolation curves.PDF"),
#     width = 25, height = 30)
# bquiet = print(zen.inext.plot)
# dev.off() 

# Get estimated richness at lowest level of coverage
zen.estD = estimateD(t(D.mat.site), base = "coverage")

# Unite with raw richness information
zen.estD = data.frame(zen.estD, raw = rowSums(D.mat.site > 0))

names(zen.estD)[5] = "rarefied"

# Bivariate correlation
cor(zen.estD$rarefied, zen.estD$raw)

# Plot results
ggplot(zen.estD, aes(x = raw, y = rarefied)) +
  geom_point() +
  stat_smooth(method = "lm") +
  labs(x = "Raw S", y = "Rarefied S")

# write.csv(zen.estD[, c("Site", "m", "SC", "rarefied", "raw")], "ZEN 2014 Rarefaction.csv")
