###----------------------------------------------------------------------###
### Paper: "The prevalence weapon damage: a proportional meta-analysis   ###
### Sarah M. Lane and Erin L. McCullough                                 ###                                                                    ###
### Code Author: Erin McCullough, Clark University                       ###
### Date: September 2024                                                 ###
###----------------------------------------------------------------------###

library(metafor)
library(ggplot2)
library(ape)
library(rotl)
library(dplyr)

# To install the orchaRd package:
install.packages("pacman")
pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, emmeans)
devtools::install_github("daniel1noble/orchaRd", force = TRUE)
library(orchaRd)

#------------------#
# 1. Phylogeny     #
#------------------#

data <- read.csv("meta_data.csv", h=T, stringsAsFactors = TRUE)
data$rotl_spp <- gsub(" ", "_", data$rotl_spp) # Replace space between names with underscore
data$rotl_spp <- as.factor(data$rotl_spp)
str(data)

# Get list of species from the OTL database
taxa <- levels(data$rotl_spp)
resolved_names <- tnrs_match_names(taxa,context_name = "Animals")
resolved_names # Check names in first two columns match

tree <- tol_induced_subtree(ott_ids = resolved_names$ott_id)
tree$tip.label <- strip_ott_ids(tree$tip.label) # Remove ott IDs for presentation
tree$node.label <- NULL # Remove node labels (these might be a problem later)
summary.phylo(tree) # Check number of tips match number of species
plot(tree, no.margin=TRUE, cex=0.5)

#--------------#
# 2. Setup     #
#--------------#

data <- read.csv("meta_data.csv", h=T, stringsAsFactors = TRUE)
data$Observation <- factor(data$Observation) # Unique identifier for each observation
data$Study <- factor(data$Study)
data$rotl_spp <- gsub(" ", "_", data$rotl_spp) # Replace space between names with underscore (species name in OTL database)
data$rotl_spp <- factor(data$rotl_spp)
data$Species <- gsub(" ", "_", data$Species) # Replace space between names with underscore (species name from study)
data$Species <- factor(data$Species)
data$Size <- factor(data$Size) # Two size categories: large and small
data$Regenerate <- factor(data$Regenerate) # Y/N can regenerate

# compute individual effect sizes using double arcsine transformation (because many values close to 0)
# add = 0 to prevent default adjustment because double arcsine transformation doesn't need to add 0.5 when have proportions equal to 0
data_es <- escalc(xi = Injured, ni = Total, data = data, measure = "PFT", add = 0)

str(data_es)
nlevels(data_es$Species) # Check number of species
nlevels(data_es$Observation) # Check number of studies

tree_grafen = compute.brlen(tree, method="Grafen", power=1)
plot(ladderize(tree_grafen), cex=0.6)
phylo_matrix <- vcv(tree_grafen, cor=TRUE, model="Brownian") # Make phylogenetic matrix

# Checking species names in tree and data file match
levels(data_es$rotl_spp) %in% row.names(phylo_matrix) # Should all say TRUE

#-----------------------#
# 3. Overall models     #
#-----------------------#

# Simple model (no random effects)
meta1 <- rma.uni(yi, vi, data= data_es, method= "REML")
summary(meta1) 

# Adding four random factors
meta2 <- rma.mv(yi, vi, random= list(~ 1|rotl_spp, ~ 1|Species, ~ 1|Study, ~1|Observation), 
                R= list(rotl_spp = phylo_matrix), data= data_es, method= "REML")
summary(meta2)
i2_ml(meta2, method=c("ratio")) # Heterogeneity at each random factor level

# Converting transformed proportion back to proportion
#Equation 13 from Wang 2023 (recommended)
proportion <- predict(meta2, transf = transf.ipft.hm, targ = list(ni = 1/(meta2$se)^2))

#-----------------------#
# 4. Meta-regressions   #
#-----------------------#

# Single categorical factor added as a fixed effect (Regeneration)
meta_regen <- rma.mv(yi, vi, mods = ~Regenerate, random= list(~ 1|rotl_spp, ~ 1|Species, ~ 1|Study, ~1|Observation), 
                     R= list(rotl_spp = phylo_matrix), data= data_es, method= "REML")
summary(meta_regen) 
# Look at "Test of Moderators"
# QM= test statistic, p-val= whether categories differ in the average correlation
r2_ml(meta_regen) # Marginal r squared (variance explained by fixed effect)

# To get means for each category- run the model without the intercept
meta_regen2 <- rma.mv(yi, vi, mod = ~Regenerate-1, random= list(~ 1|rotl_spp, ~ 1|Species, ~ 1|Study, ~1|Observation), 
                      R= list(rotl_spp = phylo_matrix), data= data_es, method= "REML")

pred_regen2 <- predict(meta_regen2, transf = transf.ipft.hm, targ = list(ni = 1/(meta_regen2$se)^2))

# Single categorical factor added as a fixed effect (Size)
meta_size <- rma.mv(yi, vi, mod = ~Size, random= list(~ 1|rotl_spp, ~ 1|Species, ~ 1|Study, ~1|Observation), 
                      R= list(rotl_spp = phylo_matrix), data= data_es, method= "REML")
summary(meta_size) 
r2_ml(meta_size) 

meta_size2 <- rma.mv(yi, vi, mod = ~Size-1, random= list(~ 1|rotl_spp, ~ 1|Species, ~ 1|Study, ~1|Observation), 
                      R= list(rotl_spp = phylo_matrix), data= data_es, method= "REML")

pred_size2 <- predict(meta_size2, transf = transf.ipft.hm, targ = list(ni = 1/(meta_size2$se)^2))

# Subsetting the data for the ramming behavior analysis
ram_data <- subset(data_es, Ram!="")
ram_data$Ram <- factor(ram_data$Ram)
ram_data$Observation <- factor(ram_data$Observation)
ram_data$Study <- factor(ram_data$Study)
ram_data$rotl_spp <- factor(ram_data$rotl_spp)
ram_data$Species <- factor(ram_data$Species)
ram_data_tips <- setdiff(levels(data_es$rotl_spp), levels(ram_data$rotl_spp)) # Species not found in full list
ram_tree <- drop.tip(tree_grafen, ram_data_tips) # Prune tree by removing species not found in the full list
ram_matrix <- vcv(ram_tree, cor=TRUE, model="Brownian")

meta_ram <- rma.mv(yi, vi, mod = ~Ram, random= list(~ 1|rotl_spp, ~ 1|Species, ~ 1|Study, ~1|Observation), 
                    R= list(rotl_spp = ram_matrix), data= ram_data, method= "REML",
                   control=list(rel.tol=1e-8))
#Needed to adjust threshold in order to achieve convergence
#Default is rel.tol = 1e-10, which was too strict

summary(meta_ram) 
r2_ml(meta_ram) 

meta_ram2 <- rma.mv(yi, vi, mod = ~Ram-1, random= list(~ 1|rotl_spp, ~ 1|Species, ~ 1|Study, ~1|Observation), 
                      R= list(rotl_spp = ram_matrix), data= ram_data, method= "REML")

pred_ram2 <- predict(meta_ram2, transf = transf.ipft.hm, targ = list(ni = 1/(meta_ram2$se)^2))

# Subsetting the data for the allometry analysis
allo_data <- subset(data_es, Allometry!="")
allo_data$Observation <- factor(allo_data$Observation)
allo_data$Study <- factor(allo_data$Study)
allo_data$rotl_spp <- factor(allo_data$rotl_spp)
allo_data$Species <- factor(allo_data$Species)
allo_data_tips <- setdiff(levels(data_es$rotl_spp), levels(allo_data$rotl_spp)) # Species not found in full list
allo_tree <- drop.tip(tree_grafen, allo_data_tips) # Prune tree by removing species not found in the full list
allo_matrix <- vcv(allo_tree, cor=TRUE, model="Brownian")

meta_allo <- rma.mv(yi, vi, mod = ~Allometry, random= list(~ 1|rotl_spp, ~ 1|Species, ~ 1|Study, ~1|Observation), 
                   R= list(rotl_spp = allo_matrix), data= allo_data, method= "REML")

summary(meta_allo) 
r2_ml(meta_allo) 

#--------------#
# 5. Figures   #
#--------------#

# Plot meta-regression results

dev.off()

#regeneration
plot(y=c(pred_regen2$pred[1],pred_regen2$pred[82]), x = 0:1, 
     xlim=c(-0.5, 1.5), ylim = c(-0.1,1), las=1, cex=2.5, 
     pch=21, bg="black", xlab = "Regeneration", xaxt='n', ylab = "Prevalence (%)", col=F)
axis(side = 1, at=0:1, labels = c("No", "Yes"), las=1)
segments(x0=0:1, y0 = c(pred_regen2$ci.lb[1],pred_regen2$ci.lb[82]), y1 = c(pred_regen2$ci.ub[1],pred_regen2$ci.ub[82]), lwd = 2, col = c("black"))
#[82] and [1] are rows for regen and non-regen observations, respectively
#saved as 4x5 (portrait)

#size
plot(y=c(pred_size2$pred[1],pred_size2$pred[2]), x = 0:1, xlim=c(-0.5, 1.5), ylim = c(-0.1,1), las=1, cex=2.5, pch=21, bg="black", xlab = "Relative weapon size", xaxt='n', ylab = "Prevalence (%)", col=F)
axis(side = 1, at=0:1, labels = c("Small", "Large"), las=1)
segments(x0=0:1, y0 = c(pred_size2$ci.lb[1],pred_size2$ci.lb[2]), y1 = c(pred_size2$ci.ub[1],pred_size2$ci.ub[2]), lwd = 2, col = c("black"))
#[1] and [2] are rows for small and large observations, respectively
#saved as 4x5 (portrait)

#ramming
plot(y=c(pred_ram2$pred[23],pred_ram2$pred[1]), x = 0:1, xlim=c(-0.5, 1.5), ylim = c(-0.1,1), las=1, cex=2.5, pch=21, bg="black", xlab = "Ramming", xaxt='n', ylab = "Prevalence (%)", col=F)
axis(side = 1, at=0:1, labels = c("No", "Yes"), las=1)
segments(x0=0:1, y0 = c(pred_ram2$ci.lb[23],pred_ram2$ci.lb[1]), y1 = c(pred_ram2$ci.ub[23],pred_ram2$ci.ub[1]), lwd = 2, col = c("black"))
#[23] and [1] are rows for non-ramming and ramming observations, respectively

# allometry
# point size scaled by the precision (inverse of standard error) of each estimate
regplot(meta_allo,mod = "Allometry",
        transf = transf.ipft.hm, 
        targ = list(ni = 1/(meta_allo$se)^2),
        psize = "seinv",
        shade = "white",
        bg = "gray90")

# Plot forest plot
data <- read.csv("meta_data_allo.csv", h=T, stringsAsFactors = TRUE)
meta.data <- data %>% select(Class, Order, Species, Weapon, Injured, Total, Regenerate, Study)
ordem = paste("expression(italic('", meta.data$Species, "'))", sep = "")
ital = sapply(ordem, function(x){eval(parse(text = x))})
ilab.text<-data.frame(meta.data$Class, meta.data$Order, meta.data$Weapon, meta.data$Study, meta.data$Regenerate, meta.data$Injured, meta.data$Total)

forest(meta2, header = TRUE, transf = transf.ipft.hm, targ = list(ni = 1/(meta2$se)^2),
       ilab=ilab.text,
       slab = ital, bg = "gray50",pch = 21,
       xlab=expression("Proportion damaged"))
#saved as 7x9 (landscape)

