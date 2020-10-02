## R script for analysing phylogenetic signal as done in:
## Engler et al. (revision). Niche evolution reveals disparate signatures of speciation in the 'great speciator' (white-eyes, Aves: Zosterops)
##
## see Supplement S2 for further details
##
## coded for R 3.6.3
##
## code by Jan O. Engler 
## https://github.com/joengler
## 
## further code can be accessed via:
## https://github.com/JOEngler/ZostiNicheEvol
##
## all necessary data available via FigShare:
## https://figshare.com/articles/dataset/Data_set_to_study_niche_evolution_in_White-eyes/13042031


## Contents
## 1) load and prepare data
##
## 2) testing and visualize phylogenetic signal estimates
##
## 3) assessing general and local phylogenetic correlation

## libs needed
library(phytools)
library(phylosignal)
library(phylobase)
setwd("D:/PATH/TO/FOLDER") #edit to your needs 


## 1) load and prepare data

# load full tree in newark format
mytree <-  read.tree(file = paste0(getwd(),"/FigshareData/MoreData/4.66%FinalTree_newark.txt"), text = NULL, tree.names = NULL, skip = 0, comment.char = "#", keep.multi = FALSE)

# read climate data (median values for five niche axes)
bclim <- read.csv(paste0(getwd(),"/FigshareData/MoreData/bclim.csv"))

# manipulating the tree
tipnames <- mytree$tip.label #extract species names
newtree <- drop.tip(mytree, tipnames[c(1:4,6,8:11,13:16,18:24,26:34,36,38:42,45)]) # drop species not needed


## 2) testing and visualize phylogenetic signal estimates

# create a phylo4d object from the data available
p4d <- phylo4d(newtree, tip.data = bclim[,-1])#, node.data = fitEB$ace) #this is only for one bioclimatic variable for testing ideally, there should be all predictors computed already! 

# visualize the trait(s)
barplot.phylo4d(p4d, tree.type = "phylo", tree.ladderize = FALSE)

# testing phyl signal across methods ((Blomberg’s K and K*, Abouheif’s Cmean, Moran’s I, and Pagel’s Lambda)
phySignal <- phyloSignal(p4d = p4d, method = "all")
phySignal

# assessing behavior of these methods using a simulation
phylosim <- phyloSim(tree = newtree, method = "all", nsim = 1000, reps = 99)
plot(phylosim, stacked.methods = FALSE, quantiles = c(0.05, 0.95))
plot.phylosim(phylosim, what = "pval", stacked.methods = TRUE)

# Phylogenetic signal estimation as a fraction of a Brownian Motion process.
phylosimsig <- phyloSimSignal(p4d, phylosim, quantiles = c(0.05, 0.95))
plot(phylosimsig)



## 3) assessing general and local phylogenetic correlation

# assessing phylogenetic correlation for each trait (Moran's I)
x.crlg <- phyloCorrelogram(p4d, trait = "bio15") #for each trait separately exchange niche axis by hand
plot(x.crlg)

# assessing phylogenetic correlation for multiple traita (Mantel's R)
t.crlg <- phyloCorrelogram(p4d, trait = c("bio04", "bio10")) #...for all thermal axes
plot(t.crlg)

p.crlg <- phyloCorrelogram(p4d, trait = c("bio12","bio15")) #...for all rainfall axes
plot(p.crlg)

# locating signal with Local Indicator of Phylogenetic Association (LIPA) analyses
carni.lipa <- lipaMoran(p4d, alternative = "two-sided")
carni.lipa.p4d <- lipaMoran(p4d, alternative = "two-sided", as.p4d = TRUE)

barplot.phylo4d(p4d, bar.col=(carni.lipa$p.value < 0.05) + 1, center = FALSE , scale = FALSE)

barplot.phylo4d(carni.lipa.p4d, bar.col = (carni.lipa$p.value < 0.05) + 1, center = FALSE, scale = FALSE)
