###############################################################################
#
# Classifying gap closure area to ground regeneration and lateral canopy growth
#
#
##############################################################################

# Definition criteria

# Lateral canopy growth (due to crown plasticity) is classified when growth pixel
# is adjacent to the gap edge 
# and growth rate is > 0.5 m/yr 

# --- load libaries

library(terra)
library(sf)


# --- load layers ----

gaps2009 <- rast("data/processed/gaps_final/berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps2017 <- rast("data/processed/gaps_final/berchtesgaden_2017_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")


# create 1km grid across study area for gap boundary  delineation
aoi <- vect("data/raw/npb_zonierung_22_epsg25832.shp")

grid <- st_make_grid(aoi,
                     cellsize= 1000, # 1km
                     square=T)
grid <- vect(grid)
grid$id <- 1:nrow(grid)


# --- create forest/gap edge mask for lateral growth classification ---

# need to loop through a subset of gaps as R crashes when whole gap raster is processed

# --- 2009 ---

# loop where gap layer is cropped to grid cell and processed

grid_id <- grid$id

for (i in grid_id) {
  print(i)
  g.sub <- grid[grid$id==i]
  gaps <- gaps2009
  gaps.sub <- crop(gaps2009,g.sub)
  print("finished crop")
  bo <- boundaries(gaps.sub, directions=8, inner=TRUE)
  print("boundaries created")
  writeRaster(bo, paste0("data/processed/closure/gap_boundaries/2009/gap_boundaries_9.", i, ".tif"))
}


# --- create gap boundary layer for the whole National Park for 2009

bo2009.list <- list.files("data/processed/closure/gap_boundaries/2009/", full.names = TRUE)
r <- lapply(bo2009.list, raster)

r$overwrite <- TRUE
boundaries.2009 <- do.call(merge, r) # merge all tiles into one rasterlayer

writeRaster(boundaries.2009, "data/processed/closure/gap_boundaries_2009.tif")
#reprojecting the boundary layer in QGis to crs 25832


# --- 2017 ---

for (i in grid_id) {
  print(i)
  g.sub <- grid[grid$id==i]
  gaps <- gaps2017
  gaps.sub <- crop(gaps2017,g.sub)
  print("finished crop")
  bo <- boundaries(gaps.sub, directions=8, inner=TRUE)
  print("boundaries created")
  writeRaster(bo, paste0("data/processed/closure/gap_boundaries/2017/gap_boundaries_17.", i, ".tif"), overwrite=T)
}

# --- create gap boundary layer for the whole National Park for 2017

bo2017.list <- list.files("data/processed/closure/gap_boundaries/2017/gap_boundaries17", full.names = TRUE)
r <- lapply(bo2017.list, raster)

r$overwrite <- TRUE
boundaries.2017 <- do.call(merge, r) # merge all tiles into one rasterlayer

writeRaster(boundaries.2017, "data/processed/closure/gap_boundaries_2017.tif")
#reprojecting the boundary layer in QGis to crs 25832



# --- classify vertical and horizontal closure ---

#load closure areas with growth information

clo_growth_917 <- rast("data/processed/closure/closure_area_growth_917.tif")
clo_growth_1721 <- rast("data/processed/closure/closure_area_growth_1721.tif")

#load gap boundaries

boundaries.2009 <- rast("data/processed/closure/gap_boundaries_2009_25832.tif")
boundaries.2017 <- rast("data/processed/closure/gap_boundaries_2017_25832.tif")

#adjust extents
clo_growth_917 <- crop(clo_growth_917, boundaries.2009)
clo_growth_1721 <- crop(clo_growth_1721, boundaries.2017)


#  -- need to differ between processing both time steps due to different growing periods

gap_closure_mechanism917 <- function(diff_closure_layer, boundary_layer){       # 0.5 m * 8 = 4m height gain

  closure_mechanism <- diff_closure_layer

  # classify change group
  closure_mechanism[diff_closure_layer > 4  & boundary_layer ==1 ] <- 1 #horizontal closure (crown plasticity)
  closure_mechanism[diff_closure_layer >= 4  & boundary_layer ==0 ] <- 2 # above average vertical closure
  closure_mechanism[diff_closure_layer <= 4 & boundary_layer ==0 ] <- 2 #vertical closure (ground regeneration)
  closure_mechanism[diff_closure_layer <= 4 & boundary_layer ==1 ] <- 2 #vertical closure (ground regeneration)
  closure_mechanism[diff_closure_layer <=0 ] <- 0 # no closure
  return(closure_mechanism)
}

gap_closure_mechanism1721 <- function(diff_closure_layer, boundary_layer){      # 0.5 m * 4 = 2m height gain

  closure_mechanism <- diff_closure_layer

  # classify change group
  closure_mechanism[diff_closure_layer > 2  & boundary_layer ==1 ] <- 1 #horizontal closure (crown plasticity)
  closure_mechanism[diff_closure_layer >= 2  & boundary_layer ==0 ] <- 2 # above average vertical closure
  closure_mechanism[diff_closure_layer <= 2 & boundary_layer ==0 ] <- 2 #vertical closure (ground regeneration)
  closure_mechanism[diff_closure_layer <= 2 & boundary_layer ==1 ] <- 2 #vertical closure (ground regeneration)
  closure_mechanism[diff_closure_layer <=0 ] <- 0 # no closure
  return(closure_mechanism)
}

gap_closure_mechanism917 <- gap_closure_mechanism917(clo_growth_917, boundaries.2009)
gap_closure_mechanism1721 <- gap_closure_mechanism1721(clo_growth_1721, boundaries.2017)


terra::writeRaster(gap_closure_mechanism917, "data/processed/closure/gap_closure_mechanism917.tif")
terra::writeRaster(gap_closure_mechanism1721, "data/processed/closure/gap_closure_mechanism1721.tif")


