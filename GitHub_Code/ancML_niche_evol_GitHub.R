## R script for running phytools for ancestral niche reconstruction in e-space as done in:
## Engler et al. (revision). Niche evolution reveals disparate signatures of speciation in the 'great speciator' (white-eyes, Aves: Zosterops)
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
## 1) load and clean phylogenetic tree, distribution and climate data
##
## 2) compare multiple evolutionary models
##
## 3) visualize variation in the phenogram using subsets of occurrence information to estimate ancestral states

## libs needed
library(phytools)
library(raster)
library(scales)
setwd("D:/PATH/TO/FOLDER") #edit to your needs 

## 1) load and clean phylogenetic tree, distribution and climate data
# load full tree in newark format
mytree <-  read.tree(file = paste0(getwd(),"/FigshareData/MoreData/4.66%FinalTree_newark.txt"), text = NULL, tree.names = NULL, skip = 0, comment.char = "#", keep.multi = FALSE)

# manipulating the tree
tipnames <- mytree$tip.label #extract species names
newtree <- drop.tip(mytree, tipnames[c(1:4,6,8:11,13:16,18:24,26:34,36,38:42,45)]) # drop species not needed
tipnames <- newtree$tip.label #extract species names again for downstream analyses

# load occurrence and climate data
point_sp <- list.files(paste0(getwd(),"/FigshareData/Occ"))
env_files <- list.files(paste0(getwd(),"/FigshareData/Climate"), full.names = TRUE)


## 2) compare multiple evolutionary models
#
# These analyses are done for each climate niche axis separately
# the script requires manual edit of which environmental layer to load and it's name
# Legend:
# env_files[1] --> "bio01", env_files[2] --> "bio04", env_files[3] --> "bio10" for thermal niche axes
# env_files[4] --> "bio12", env_files[5] --> "bio15" for precipitation related niche axes

iras <- raster(env_files[1])
names(iras) <- "bio01"

conversion <- c(1, 10, 7, 9, 8, 5, 6, 2, 3, 4) #conversion between csv point localities and tip names

# extract climatic information at presence locations for each species
Species <- data.frame()
for(i in 1:length(tipnames)){
  S1occ <- read.csv(paste0(getwd(),"/FigshareData/Occ/",point_sp[conversion[i]]))
  coordinates(S1occ) <- ~lon+lat
  S1_vals <- extract(iras, S1occ) #climate at presence locations
  S1_vals <- S1_vals[complete.cases(S1_vals)] #delete NAs
  S1_vals <- data.frame(S1_vals)
  S1_vals$ID <- tipnames[i]
  
  Species <- rbind(Species, S1_vals)
}

# calculates the median value for each species
x <- tapply(Species[,1],Species$ID, median) 

# fits an evolutionary model 
fitEB <- anc.ML(newtree, x, model = "EB") #Early Burst
fitBM <- anc.ML(newtree, x, model = "BM") #Brownian Motion
fitOU <- anc.ML(newtree, x, model = "OU") #single optimum OU

# plots a phenogram of all three models compared
phenogram(newtree,c(x,fitEB$ace),col=alpha("black", 0.9),ftype="off") 
phenogram(newtree,c(x,fitBM$ace),col=alpha("red", 0.9),ftype="off", add = TRUE) 
phenogram(newtree,c(x,fitOU$ace),col=alpha("blue", 0.9),ftype="off", add = TRUE)  

# print all log-likelihoods in one object for better comparability
logLik <- data.frame(fitEB$logLik, fitBM$logLik, fitOU$logLik)
names(logLik)<-c(fitEB$model, fitBM$model, fitOU$model)
logLik #EB with the maximum logLik


## 3) visualize variation in the phenogram using subsets of occurrence information to estimate ancestral states

# continue with the EB model
phenogram(newtree,c(x,fitEB$ace),col=alpha("black", 0.5),ftype="reg", fsize = 0.8) # base tree with full data

for(k in 1:100){
  x <- c()
  for(o in 1:length(tipnames)){
    S1 <- Species[Species$ID == tipnames[o],]
    x <- c(x, median(sample(S1[,1], nrow(S1)*0.3))) #randomly pick 30% of data
  }
  names(x) <- tipnames
  fitEB <- anc.ML(newtree, x, model = "EB")
  
  phenogram(newtree,c(x,fitEB$ace),ftype="off",col=alpha("red", 0.01), add = TRUE)
}
phenogram(newtree,c(x,fitEB$ace),col=alpha("black", 0.9), add = TRUE) # base tree with full data

