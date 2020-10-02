## R script for performing niche comparisons in e-space:
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
## 1) prepare input data and objects
##
## 2) calculate niche overlaps and conduct test for niche similarity and equivalency

## libs needed
library(ecospat)
library(raster)
library(dismo)
setwd("D:/PATH/TO/FOLDER") #edit to your needs 

## 1) prepare input data and objects
# create background records from across the genus' range
ZGR <- raster(paste0(getwd(),"/FigshareData/MoreData/Zosterops_Genus_Range.asc")) #genus' range raster to extract background values from
env_files <- list.files(paste0(getwd(),"/FigshareData/Climate"), full.names = TRUE)
iras <- stack(env_files)

bg <- randomPoints(ZGR, n = 10000) #10000 background records
bg <- as.data.frame(bg)
coordinates(bg) <- ~x+y

bg_vals <- extract(iras, bg) #climate at presence locations
bg_vals <- data.frame(bg_vals)
bg_vals$ID <- 0

# load occurrence information
point_sp <- list.files(paste0(getwd(),"/FigshareData/Occ"))

# create labelled pairwise matrices
D_mat <- matrix(nrow = length(point_sp), ncol = length(point_sp), dimnames = list(strsplit(point_sp, split = ".csv"), strsplit(point_sp, split = ".csv")))
I_mat <- matrix(nrow = length(point_sp), ncol = length(point_sp), dimnames = list(strsplit(point_sp, split = ".csv"), strsplit(point_sp, split = ".csv")))
eq_p <- matrix(nrow = length(point_sp), ncol = length(point_sp), dimnames = list(strsplit(point_sp, split = ".csv"), strsplit(point_sp, split = ".csv")))
sim_p <- matrix(nrow = length(point_sp), ncol = length(point_sp), dimnames = list(strsplit(point_sp, split = ".csv"), strsplit(point_sp, split = ".csv")))


## 2) calculate niche overlaps and conduct test for niche similarity and equivalency
for(i in 1:length(point_sp)){
  S1occ <- read.csv(paste0(getwd(),"/FigshareData/Occ/",point_sp[i]))
  coordinates(S1occ) <- ~lon+lat
  S1_vals <- extract(iras, S1occ) #climate at presence locations
  S1_vals <- data.frame(S1_vals)
  S1_vals$ID <- 1
  S1 <- rbind(bg_vals, S1_vals)
  S1 <- S1[complete.cases(S1),]
  names(S1) <- c("bio01", "bio04", "bio10", "bio12", "bio15", "ID")
  
  for(k in 1:length(point_sp)){
    S2occ <- read.csv(paste0(getwd(),"/FigshareData/Occ/",point_lrsp[k]))
    coordinates(S2occ) <- ~lon+lat
    S2_vals <- extract(iras, S2occ) 
    S2_vals <- data.frame(S2_vals)
    S2_vals$ID <- 1
    S2 <- rbind(bg_vals, S2_vals)
    S2 <- S2[complete.cases(S2),]
    names(S2) <- c("bio01", "bio04", "bio10", "bio12", "bio15", "ID")
    pca.env <- dudi.pca(rbind(S1,S2)[,1:5],scannf=F,nf=2) 
    
    # PCA scores
    scores.globclim <- pca.env$li #for the whole study area
    
    scores.sp.S1 <- suprow(pca.env,S1[which(S1$ID==1),1:5])$li    # for species 1
    scores.sp.S2 <- suprow(pca.env,S2[which(S2$ID==1),1:5])$li  # for species 2
    scores.clim.S1 <- suprow(pca.env,S1[,1:5])$li #for the entire study area S1
    scores.clim.S2 <- suprow(pca.env,S2[,1:5])$li #for the entire study area S2
    
    # gridding the niches
    grid.clim.S1 <- ecospat.grid.clim.dyn(glob=scores.globclim, 
                                          glob1=scores.clim.S1,
                                          sp=scores.sp.S1, R=100,
                                          th.sp=0) 
    
    grid.clim.S2 <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                          glob1=scores.clim.S2,
                                          sp=scores.sp.S2, R=100,
                                          th.sp=0) 
    
    # Compute Schoener's D, index of niche overlap
    D_mat[i,k] <- ecospat.niche.overlap (grid.clim.S1, grid.clim.S2, cor=T)$D 
    I_mat[i,k] <- ecospat.niche.overlap (grid.clim.S1, grid.clim.S2, cor=T)$I 
    
    
    eq.test <- ecospat.niche.equivalency.test(grid.clim.S1, grid.clim.S2,
                                              rep=1000, alternative = "lower")
    
    eq_p[i,k] <- eq.test$p.D    
    
    
    sim.test <- ecospat.niche.similarity.test(grid.clim.S1, grid.clim.S2,
                                              rep=1000, alternative = "lower",
                                              rand.type=2) 
    sim_p[i,k] <- sim.test$p.D
    
    # save summary plots
    pdf(paste0(getwd(),"/FigshareData/NicheTests/",strsplit(point_lrsp[i], split = ".csv")," vs. ",strsplit(point_lrsp[k], split = ".csv"),".pdf"))
    par(mfrow=c(2,2))
    ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
    ecospat.plot.overlap.test(sim.test, "D", "Similarity")
    ecospat.plot.niche.dyn(grid.clim.S1, grid.clim.S2, quant=0.25, interest=2,
                           title= "Niche Overlap", name.axis1="PC1",
                           name.axis2="PC2")
    ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)
    dev.off()
    
  }
}

# write output files
write.csv(D_mat, paste0(getwd(),"/FigshareData/NicheTests/D_mat.csv"))
write.csv(I_mat, paste0(getwd(),"/FigshareData/NicheTests/I_mat.csv"))
write.csv(eq_p, paste0(getwd(),"/FigshareData/NicheTests/Equivalency_P.csv"))
write.csv(sim_p, paste0(getwd(),"/FigshareData/NicheTests/Similarity_P.csv"))



