#######################################
#
# identifying gap formation and closure areas
#
######################################

# --- libraries

library(dplyr)
library(tidyr)
library(terra)



# --- load gap layers ----

gap_stack <- rast("data/processed/gaps_final/gaps_masked.tif")

gaps2009 <- gap_stack[[1]]
gaps2017 <- gap_stack[[2]]
gaps2021 <- gap_stack[[3]]

#--- define functions ---

# function to identify gap formation and gap closure

gap_change_class <- function(gap_layer1, gap_layer2){
  exp_clo <- rast() #create empty raster to classify
  ext(exp_clo) <- ext(gap_layer1)
  res(exp_clo) <- res(gap_layer1)
  crs(exp_clo) <- crs(gap_layer1)
  # classify change group
  exp_clo[gap_layer1 >0  & is.na(gap_layer2) ] <- 1 #gap closure
  exp_clo[is.na(gap_layer1) & gap_layer2 >0] <- 2 #gap expansion
  return(exp_clo)
} 



exp_clo_917 <- gap_change_class(gaps2009, gaps2017) # between 2009 and 2017
exp_clo_1721 <- gap_change_class(gaps2017, gaps2021) # between 2017 and 2021


writeRaster(exp_clo_917, "data/processed/gap_change/formation_closure_917.tif")
writeRaster(exp_clo_1721, "data/processed/gap_change/formation_closure_1721.tif")


# --- extract vegetation growth in gap closure areas per time step ---

# CHM cannot be provided, but resulting growth values within gap area
# is provided in zenodo

chm9 <- rast("data/processed/CHM_data/chm9_artifacts_masked.tif")
chm17 <- rast("data/processed/CHM_data/chm17_artifacts_masked.tif")
chm21 <- rast("data/processed/CHM_data/chm21_artifacts_masked.tif")

chm9 <- crop(chm9, chm21)
chm17 <- crop(chm17, chm21)

# get vegetation changes
diff917 <- chm17-chm9
diff1721 <- chm21 - chm17

# extract only closure areas
clo_917 <- classify(exp_clo_917, cbind(2, NA)) #replace 2=expansion with NA to get only closure areas
clo_1721 <- classify(exp_clo_1721, cbind(2, NA)) #replace 2=expansion with NA to get only closure areas

# extract vegetation growth in closure areas for the two respective time steps

diff917 <- crop(diff917,clo_917 )
clo_growth_917 <-mask(diff917, clo_917) 

diff1721 <- crop(diff1721, clo_1721)
clo_growth_1721 <-mask(diff1721, clo_1721) 


writeRaster(clo_growth_917 , "data/processed/closure/closure_area_growth_917.tif") # data provided!
writeRaster(clo_growth_1721 , "data/processed/closure/closure_area_growth_1721.tif")
