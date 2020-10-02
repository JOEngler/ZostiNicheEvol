## R script for running basic ENMs for ancestral niche reconstruction in g-space as done in:
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
## 1) Do the Bioclim ENMs
##
## 2) Calculate AUC, COR, P/A thresholds for each ENM
##
## 3) ENM niche overlaps
##
## 4) Ancestral Niche Reconstruction with the phyloclim package


## libs needed
library(raster)
library(dismo)
setwd("D:/PATH/TO/FOLDER") #edit to your needs 


## 1) Do the Bioclim ENMs

# load in the data
point_sp <- list.files(paste0(getwd(),"/FigshareData/Occ"))
pred2.5 <- list.files(paste0(getwd(),"/FigshareData/Climate"), full.names = TRUE)
pred <- stack(pred2.5)


# Do the ENMs in a loop with projection to the background (genus's range)
for(i in 1:length(point_sp)){
  rec_i <- read.csv(paste0(getwd(),"/FigshareData/Occ/",point_sp[i]))
  
  coordinates(rec_i) <- ~lon+lat
  bc <- bioclim(pred, rec_i) 
  pd_i <- dismo::predict(pred, bc)
  writeRaster(pd_i, filename = paste0(getwd(),"/FigshareData/ENMs/",strsplit(point_sp[i], split = ".csv"),".asc"), overwrite = TRUE)
}


## 2) Calculate AUC, COR, P/A thresholds for each ENM

# list all prediction rasters
asc_lrsp <- paste0(getwd(),"/FigshareData/ENMs/",strsplit(point_sp, split = ".csv"),".asc")

# extract vals from presence and background
ZGR <- raster(paste0(getwd(),"/FigshareData/MoreData/Zosterops_Genus_Range.asc")) #genus' range raster to extract background values from

bg <- randomPoints(ZGR, n = 10000) #10000 background records; ignore warning message as assumption is correct
bg <- as.data.frame(bg)
coordinates(bg) <- ~x+y

zAUC <- c()
zCOR <- c()
zTR <- c()
for(i in 1:length(point_sp)){
  iras <- raster(paste0(getwd(),"/FigshareData/ENMs/",asc_lrsp[i]))
  iocc <- read.csv(paste0(getwd(),"/FigshareData/Occ/",point_sp[i]))
  coordinates(iocc) <- ~lon+lat
  p_vals <- extract(iras, iocc) #climate at presence locations
  bg_vals <- extract(iras, bg) #climate at background locations
  
  
  eval <- evaluate(p_vals, bg_vals)
  otr <- threshold(eval, 'spec_sens')
  
  zTR <- c(zTR, otr)
  zAUC <- c(zAUC, eval@auc)
  zCOR <- c(zCOR, eval@cor)
}
fit_res <- data.frame(unlist(strsplit(point_sp, split = ".csv")), zAUC, zCOR, zTR)
names(fit_res)[1]<- "species"
write.csv(fit_res, file = paste0(getwd(),"/FigshareData/MoreData/AUC_COR_SStresholds.csv"))



## 3) ENM niche overlaps

# calculate niche overlaps after applying the 'spec_sens' threshold (maxSS)
noD_matrix <- matrix(nrow = nrow(fit_res), ncol = nrow(fit_res), dimnames = list(fit_res$species,fit_res$species))
for(o in 1:(nrow(fit_res)-1)){
  for(u in 2:nrow(fit_res)){
    oras <- raster(paste0(getwd(),"/FigshareData/ENMs/",paste0(fit_res$species[o],".asc")))
    uras <- raster(paste0(getwd(),"/FigshareData/ENMs/",paste0(fit_res$species[u],".asc")))
    
    #sets vals < MaxSS to 0 to not inflate overlap values (see Rödder & Engler 2011)
    oras <- calc(oras, fun=function(x){ x[x < fit_res$zTR[o]] <- 0; return(x)} )
    uras <- calc(uras, fun=function(x){ x[x < fit_res$zTR[u]] <- 0; return(x)} )
    
    #calculate Schoener' D niche overlap and write in matrix
    noD <- nicheOverlap(oras, uras, stat = "D", mask = FALSE, checkNegatives = FALSE)
    noD_matrix[o,u]<- noD
  }
}
write.csv(noD_matrix, file = paste0(getwd(),"/FigshareData/MoreData/Niche_Overlap_D.csv"))


## 4) Ancestral Niche Reconstruction with the phyloclim package
library(phyloclim)

# load data (done for each niche axis manually)
mymodels <- paste0(getwd(),"/FigshareData/ENMs")     # folder with ascii grids from ENM
mytree <-  read.tree(file = paste0(getwd(),"/FigshareData/MoreData/4.66%FinalTree_newark.txt"), text = NULL, tree.names = NULL, skip = 0, comment.char = "#", keep.multi = FALSE)
tipnames <- mytree$tip.label #extract species names
newtree <- drop.tip(mytree, tipnames[c(1:4,6,8:11,13:16,18:24,26:34,36,38:42,45)]) # drop species not needed

# renaming the tips so they match with the grid names
newTipNames <- newtree$tip.label
x <- gsub("_", " ", newTipNames)
asc_names <- list.files(mymodels)
asc_names <- gsub(".asc", "", asc_names)

x_replace <- c(1,10,7,9,8,5,6,2,3,4) #to keep the order
x_new <- c()
for(i in 1:length(x)){
  x_new <- c(x_new,asc_names[x_replace[i]])
}
newtree$tip.label <- x_new

# create PNO profiles for each climate niche axis (these take a while to run)
myclimate <- paste0(getwd(),"/FigshareData/Climate/WC2_2.5m_wc2.0_bio_2.5m_01.asc")      
pno01 <- pno(path_bioclim=myclimate, path_model=mymodels, subset = NULL, bin_width = NULL, bin_number = 1000)
myclimate <- paste0(getwd(),"/FigshareData/Climate/WC2_2.5m_wc2.0_bio_2.5m_04.asc")     
pno04 <- pno(path_bioclim=myclimate, path_model=mymodels, subset = NULL, bin_width = NULL, bin_number = 1000)
myclimate <- paste0(getwd(),"/FigshareData/Climate/WC2_2.5m_wc2.0_bio_2.5m_10.asc")                                         
pno10 <- pno(path_bioclim=myclimate, path_model=mymodels, subset = NULL, bin_width = NULL, bin_number = 1000)
myclimate <- paste0(getwd(),"/FigshareData/Climate/WC2_2.5m_wc2.0_bio_2.5m_12.asc")       
pno12 <- pno(path_bioclim=myclimate, path_model=mymodels, subset = NULL, bin_width = NULL, bin_number = 1000)
myclimate <- paste0(getwd(),"/FigshareData/Climate/WC2_2.5m_wc2.0_bio_2.5m_15.asc")                    
pno15 <- pno(path_bioclim=myclimate, path_model=mymodels, subset = NULL, bin_width = NULL, bin_number = 1000)


# phylogenetic reconstruction of ancestral niches. Change 'pno' to the respective PNO profile
myresults <- anc.clim(target=newtree, posterior = NULL, pno=pno01, n = 1000, method = "GLS")

# plot results and save as PDF for later editing
plotAncClim(myresults, clades = NULL, density = TRUE, tipmode = 1, lwd = 1, ylab = "bio15")
savePlot(filename=paste0(getwd(),"/FigshareData/PNO/","AncClim_bio15_lrsp.pdf"), type = "pdf", device = dev.cur(), restoreConsole = T) #generates a PDF for later editing in vector graphic software


