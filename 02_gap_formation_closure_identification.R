#######################################
#
# identifying expansion and closure areas
#
######################################


library(dplyr)
library(tidyr)
library(terra)



# --- load gap layers ----

gaps2009 <- rast("data/processed/gaps_final/berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps2017 <- rast("data/processed/gaps_final/berchtesgaden_2017_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps2021 <- rast("data/processed/gaps_final/berchtesgaden_2021_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")


#--- define functions ---

# simplified function to identify gap expansion and gap closure

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


common_extent <- intersect(intersect(ext(gaps2009), ext(gaps2017)), ext(gaps2021))

gaps2009_aligned <- crop(gaps2009, common_extent)
gaps2017_aligned <- crop(gaps2017, common_extent)
gaps2021_aligned <- crop(gaps2021, common_extent)

# gaps2009 <- crop(gaps2009, gaps2021, snap="near",mask=TRUE)
# gaps2017 <-crop(gaps2017, gaps2021, snap="near",mask=TRUE)

exp_clo_917 <- gap_change_class(gaps2009_aligned, gaps2017_aligned)
exp_clo_1721 <- gap_change_class(gaps2017_aligned, gaps2021_aligned)
exp_clo_921 <- gap_change_class(gaps2009_aligned, gaps2021_aligned)


terra::writeRaster(exp_clo_917, "data/processed/gap_change/formation_closure_917_cn2cr2_mmu400n8_filtered.tif")
terra::writeRaster(exp_clo_1721, "data/processed/gap_change/formation_closure_1721_cn2cr2_mmu400n8_filtered.tif")
terra::writeRaster(exp_clo_921, "data/processed/gap_change/formation_closure_921_cn2cr2_mmu400n8_filtered.tif")

# --- extract vegetation growth in gap closure areas per time step ---

exp_clo_917 <- rast("data/processed/gap_change/formation_closure_917_cn2cr2_mmu400n8_filtered.tif")
exp_clo_1721 <- rast("data/processed/gap_change/formation_closure_1721_cn2cr2_mmu400n8_filtered.tif")

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

# extract vegetation growth in closure areas

diff917 <- crop(diff917,clo_917 )
clo_growth_917 <-mask(diff917, clo_917) 

diff1721 <- crop(diff1721, clo_1721)
clo_growth_1721 <-mask(diff1721, clo_1721) 


writeRaster(clo_growth_917 , "data/processed/closure/closure_area_growth_917.tif")
writeRaster(clo_growth_1721 , "data/processed/closure/closure_area_growth_1721.tif")
