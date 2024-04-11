###########################################################################

# Identifying number of and area covered by gaps in reserach area per year

##########################################################################


library(dplyr)
library(tidyr)
library(terra)
library(ggplot2)


# --- load layers ----

gaps2009 <- rast("data/processed/gaps_final/berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps2017 <- rast("data/processed/gaps_final/berchtesgaden_2017_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps2021 <- rast("data/processed/gaps_final/berchtesgaden_2021_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")

# crop gaps to CHM 
gaps2009 <- crop(gaps2009, gaps2021, snap="near",mask=TRUE) 
gaps2017 <-crop(gaps2017, gaps2021, snap="near",mask=TRUE) 

gap.stack <- c(gaps2009, gaps2017, gaps2021)

#prepare layers to limit analysis to core zone, below 1800m and with forest type information

# --- load NP information 

foresttype <- rast("data/processed/environment_features/forest_type2020_reclass_1m.tif")
management <- vect("data/raw/npb_zonierung_22_epsg25832.shp")
elevation.below1800 <- rast("data/processed/environment_features/elevation_below1800_200steps.tif")
aspect<-  rast("data/processed/environment_features/aspect_2021_classified_1m.tif")


# exclude management zone
core.zone <- subset(management, management$zone_id == 4, c(1:2))

#mask down to reserach area
gap.stack <- mask(gap.stack, core.zone)
gap.stack <- mask(gap.stack, elevation.below1800)
gap.stack <- mask(gap.stack, foresttype)


writeRaster(gap.stack, "data/processed/gaps_final/gaps_masked_reserach_area_paper.tif")

# convert masked gap stack to data frame & write to disk
gap.stack.df <- as.data.frame(gap.stack)
saveRDS(gap.stack.df,"data/processed/gaps_final/gaps_masked_reserach_area_paper.df.rds" )

gap.stack.df <- readRDS("data/processed/gaps_final/gaps_masked_reserach_area_paper.df.rds")
names(gap.stack.df) <- c("gaps9", "gaps17", "gaps21")

# --- summarize number & area of gaps:

# Gather the data into long format
gathered_df <- gap.stack.df %>%
  pivot_longer(cols = everything(), names_to = "year", values_to = "gap_id") %>%
  filter(!if_all(.cols = everything(), is.na))%>% # Drop rows with all NAs
  count(year, gap_id) %>%
  filter(n >= 400) %>% # filter out gaps < 400m2, which emerged due to cropping of gaps layers to reserach area
  filter(!is.nan(gap_id))# Remove rows where gap_id is NaN

# overall gap area

gathered_df %>% group_by(year)%>%
  summarise(total_gap_area_ha = round(sum(n)/10000,2))


# Count the occurrences of each value in each column
gap_counts <- gathered_df %>%
  group_by(year) %>%
  summarize(number_gaps = length(unique(gap_id)))

# overall number of gaps
sum(gap_counts$number_gaps) # 11331 



# ---- assign Gap environmental features ----


# load & crop chm cropped to research area

chm9 <- rast("data/processed/CHM_data/chm9_artifacts_masked.tif")
chm17 <- rast("data/processed/CHM_data/chm17_artifacts_masked.tif")
chm21 <- rast("data/processed/CHM_data/chm21_artifacts_masked.tif")

names(gap.stack) <- c("gaps9", "gaps17", "gaps21")

crop_chm <- function(chm){
  chm <- crop(chm, gap.stack$gaps9)
  chm <- mask(chm, core.zone)
  chm <- mask(chm, elevation.below1800)
  chm <- mask(chm, foresttype)
}

chm9 <- crop_chm(chm9)
chm17 <- crop_chm(chm17)
chm21 <- crop_chm(chm21)


# --- define functions -----

# ---- forest type

#function to assign forest type
getForestType <- function(gap_chm) {
  x <-gap_chm %>% group_by(ID, ftype) %>% #count pixels per ftype per gap
    summarize(count = n())
  #identify dominating forest type in gap area
  xx <- data_frame()
  for (i in unique(x$ID)) {
    a <- x[x$ID == i,]        #subset to one ID
    a <- a[order(-a$count),]  #order descending according to pixel counts   
    if(nrow(a) == 1) {        #if only one entry = one forest type assign that one to ID
      xx <- rbind(xx,a[1,])
    }
    if(nrow(a) > 1) {         #if more than one entry = several ftypes or NAs
      if(is.na(a[1,2])) {xx <- rbind(xx,a[2,])}   #if is.na == TRUE = no forest type info assign NA
      if(is.na(a[1,2]) == FALSE) {xx <- rbind(xx,a[1,])}  #if there is other ftype info assign that one to ID
    }
  }
  return(xx)
}


Gap_Stats_ftype <- function (gap_layer, chm_layer, foresttype, year) 
{  t <- Sys.time()
gap_list <- data.frame(terra::freq(gap_layer))
gap_list$count <- gap_list$count * terra::res(chm_layer)[1]^2
gap_list <- gap_list[!is.na(gap_list[, 1]), ]
print("convert gaps to df + area: "); print(Sys.time()-t); t <- Sys.time()
# extract raster values per gap
gap_chm <- as.data.frame(c(gap_layer, chm_layer, foresttype), na.rm = FALSE)
names(gap_chm) <- c("ID", "chm_values", "ftype")
gap_chm <- gap_chm[!is.na(gap_chm$ID),]# delete pixels without any gap
print("extract values for gaps: "); print(Sys.time()-t); t <- Sys.time()
# calculate gap statistics
gap_list$chm_max <- as.data.frame(gap_chm %>% group_by(ID) %>% #create df and take second column
                                    summarize(chm_max = max(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_min <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                    summarize(chm_min = round(min(chm_values, na.rm=TRUE))))[,2]

gap_list$chm_mean <- as.data.frame(gap_chm %>%group_by(ID) %>% 
                                     summarize(chm_mean = mean(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_sd <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                   summarize(chm_mean = stats::sd(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_range <- round(gap_list$chm_max - gap_list$chm_min, 2)
print("calculate basic stats: "); print(Sys.time()-t); t <- Sys.time()

# gap_polygon <- as.polygons(gap_layer)
# gap_list$perimeter <- perim(gap_polygon)#add perimeter  
# 
# print("get perimeter: "); print(Sys.time()-t); t <- Sys.time()

gap_list$ftype <- as.data.frame(getForestType(gap_chm))[,2]
print("get forest type: "); print(Sys.time()-t); t <- Sys.time()

gap_list$year <- as.factor(year)

gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]

colnames(gap_list) <- c("gap_id", "gap_area", 
                        "chm_max", "chm_min", "chm_mean", "chm_sd"
                        ,"chm_range", "forest_type", "year") #"perimeter", 
print("finish df wrangling and labeling: "); print(Sys.time()-t); t <- Sys.time()
return(gap_list)
}


#---elevation 

#function to assign elevation class
getElevation <- function(gap_chm) {
  x <-gap_chm %>% group_by(ID, elevation) %>% #count pixels per elevation class per gap
    summarize(count = n())
  #identify dominating elevation in gap area
  xx <- data_frame()
  for (i in unique(x$ID)) {
    a <- x[x$ID == i,]        #subset to one ID
    a <- a[order(-a$count),]  #order descending according to pixel counts   
    if(nrow(a) == 1) {        #if only one entry = one forest type assign that one to ID
      xx <- rbind(xx,a[1,])
    }
    if(nrow(a) > 1) {         #if more than one entry = several ftypes or NAs
      if(is.na(a[1,2])) {xx <- rbind(xx,a[2,])}   #if is.na == TRUE = no forest type info assign NA
      if(is.na(a[1,2]) == FALSE) {xx <- rbind(xx,a[1,])}  #if there is other ftype info assign that one to ID
    }
  }
  xx<- xx %>% mutate(elevation = as.factor(recode(elevation,
                                                  `1`="600-800",
                                                  `2`="800-1000",
                                                  `3`="1000-1200",
                                                  `4`="1200-1400",
                                                  `5`="1400-1600",
                                                  `6`="1600-1800",
                                                  `7`="1800-2000",
                                                  `8`="2000-2800")))
  return(xx)
}


Gap_Stats_elevation <- function (gap_layer, chm_layer, dtm_class, year) 
{  t <- Sys.time()
gap_list <- data.frame(terra::freq(gap_layer))
gap_list$count <- gap_list$count * terra::res(chm_layer)[1]^2
gap_list <- gap_list[!is.na(gap_list[, 1]), ]
print("convert gaps to df + area: "); print(Sys.time()-t); t <- Sys.time()
# extract raster values per gap
gap_chm <- as.data.frame(c(gap_layer, chm_layer, dtm_class), na.rm = FALSE)
names(gap_chm) <- c("ID", "chm_values", "elevation")
gap_chm <- gap_chm[!is.na(gap_chm$ID),]# delete pixels without any gap
print("extract values for gaps: "); print(Sys.time()-t); t <- Sys.time()
# calculate gap statistics
gap_list$chm_max <- as.data.frame(gap_chm %>% group_by(ID) %>% #create df and take second column
                                    summarize(chm_max = max(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_min <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                    summarize(chm_min = round(min(chm_values, na.rm=TRUE))))[,2]

gap_list$chm_mean <- as.data.frame(gap_chm %>%group_by(ID) %>% 
                                     summarize(chm_mean = mean(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_sd <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                   summarize(chm_mean = stats::sd(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_range <- round(gap_list$chm_max - gap_list$chm_min, 2)
print("calculate basic stats: "); print(Sys.time()-t); t <- Sys.time()

# gap_polygon <- as.polygons(gap_layer)
# gap_list$perimeter <- perim(gap_polygon)#add perimeter  
# 
# print("get perimeter: "); print(Sys.time()-t); t <- Sys.time()

gap_list$elevation <- as.data.frame(getElevation(gap_chm))[,2]
print("get elevation: "); print(Sys.time()-t); t <- Sys.time()

gap_list$year <- as.factor(year)

gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]

colnames(gap_list) <- c("gap_id", "gap_area", 
                        "chm_max", "chm_min", "chm_mean", "chm_sd"
                        ,"chm_range",  "elevation", "year") #"perimeter",
print("finish df wrangling and labeling: "); print(Sys.time()-t); t <- Sys.time()
return(gap_list)
}


# ---- aspect 

#function to assign elevation class
getAspect <- function(gap_chm) {
  x <-gap_chm %>% group_by(ID, aspect) %>% #count pixels per aspect class per gap
    summarize(count = n())
  #identify dominating elevation in gap area
  xx <- data_frame()
  for (i in unique(x$ID)) {
    a <- x[x$ID == i,]        #subset to one ID
    a <- a[order(-a$count),]  #order descending according to pixel counts   
    if(nrow(a) == 1) {        #if only one entry = one aspect type assign that one to ID
      xx <- rbind(xx,a[1,])
    }
    if(nrow(a) > 1) {         #if more than one entry = several aspect classes or NAs
      if(is.na(a[1,2])) {xx <- rbind(xx,a[2,])}   #if is.na == TRUE = no aspect info assign NA
      if(is.na(a[1,2]) == FALSE) {xx <- rbind(xx,a[1,])}  #if there is other aspect info assign that one to ID
    }
  }
  xx<- xx %>% mutate(aspect = as.factor(recode(aspect,
                                               `1`="North",
                                               `2`="East",
                                               `3`="South",
                                               `4`="West")))
  return(xx)
}


Gap_Stats_aspect <- function (gap_layer, chm_layer, ascpect_class, year) 
{  t <- Sys.time()
gap_list <- data.frame(terra::freq(gap_layer))
gap_list$count <- gap_list$count * terra::res(chm_layer)[1]^2
gap_list <- gap_list[!is.na(gap_list[, 1]), ]
print("convert gaps to df + area: "); print(Sys.time()-t); t <- Sys.time()
# extract raster values per gap
gap_chm <- as.data.frame(c(gap_layer, chm_layer, ascpect_class), na.rm = FALSE)
names(gap_chm) <- c("ID", "chm_values", "aspect")
gap_chm <- gap_chm[!is.na(gap_chm$ID),]# delete pixels without any gap
print("extract values for gaps: "); print(Sys.time()-t); t <- Sys.time()
# calculate gap statistics
gap_list$chm_max <- as.data.frame(gap_chm %>% group_by(ID) %>% #create df and take second column
                                    summarize(chm_max = max(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_min <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                    summarize(chm_min = round(min(chm_values, na.rm=TRUE))))[,2]

gap_list$chm_mean <- as.data.frame(gap_chm %>%group_by(ID) %>% 
                                     summarize(chm_mean = mean(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_sd <- as.data.frame(gap_chm %>% group_by(ID) %>% 
                                   summarize(chm_mean = stats::sd(chm_values, na.rm=TRUE)))[,2]

gap_list$chm_range <- round(gap_list$chm_max - gap_list$chm_min, 2)
print("calculate basic stats: "); print(Sys.time()-t); t <- Sys.time()

# gap_polygon <- as.polygons(gap_layer)
# gap_list$perimeter <- perim(gap_polygon)#add perimeter  
# 
# print("get perimeter: "); print(Sys.time()-t); t <- Sys.time()

gap_list$aspect <- as.data.frame(getAspect(gap_chm))[,2]
print("get aspect: "); print(Sys.time()-t); t <- Sys.time()

gap_list$year <- as.factor(year)

gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]

colnames(gap_list) <- c("gap_id", "gap_area", 
                        "chm_max", "chm_min", "chm_mean", "chm_sd"
                        ,"chm_range",  "aspect", "year") #"perimeter",
print("finish df wrangling and labeling: "); print(Sys.time()-t); t <- Sys.time()
return(gap_list)
}


#------- calculate gap stats per forest type

stats_ftype_2009<- Gap_Stats_ftype(gap.stack$gaps9, chm9, foresttype, "2009")
stats_ftype_2017 <- Gap_Stats_ftype(gap.stack$gaps17, chm17, foresttype, "2017")
stats_ftype_2021 <- Gap_Stats_ftype(gap.stack$gaps21, chm21, foresttype, "2021")


stats_all_ftype <- rbind(stats_ftype_2009, stats_ftype_2017, stats_ftype_2021)
# stats_all_ftype$pa_ratio <- stats_all_ftype$perimeter/stats_all_ftype$gap_area #calculate perimeter/area ratio
# stats_all_ftype$shape.index <- stats_all_ftype$perimeter/ (2*(pi*stats_all_ftype$gap_area)^0.5)
stats_all_ftype$gap_area_ha <- stats_all_ftype$gap_area/10000 #convert area to ha
#drop gaps < 400 m2 (emerge through masking of research area, as large gaps transcending management area or elevation lines get cut)
stats_all_ftype <- stats_all_ftype[stats_all_ftype$gap_area >= 400,]

saveRDS(stats_all_ftype, "data/processed/gap_features/stats_all_ftype.rds")


#------- calculate gap stats per elevation

stats_elevation_2009<- Gap_Stats_elevation(gap.stack$gaps9, chm9, elevation.below1800, "2009")
stats_elevatione_2017 <- Gap_Stats_elevation(gap.stack$gaps17, chm17, elevation.below1800, "2017")
stats_elevation_2021 <- Gap_Stats_elevation(gap.stack$gaps21, chm21, elevation.below1800, "2021")

stats_all_elevation <- rbind(stats_elevation_2009, stats_elevatione_2017, stats_elevation_2021)
# stats_all_elevation$pa_ratio <- stats_all_elevation$perimeter/stats_all_elevation$gap_area #calculate perimeter/area ratio
# stats_all_elevation$shape.index <- stats_all_elevation$perimeter/ (2*(pi*stats_all_elevation$gap_area)^0.5)
stats_all_elevation$gap_area_ha <- stats_all_elevation$gap_area/10000 #convert area to ha
#drop gaps < 400 m2 (emerge through masking of research area, as large gaps transcending management area or elevation lines get cut)
stats_all_elevation <- stats_all_elevation[stats_all_elevation$gap_area >= 400,]

stats_all_elevation$elevation <- factor(stats_all_elevation$elevation , levels=c("600-800", "800-1000", "1000-1200",
                                                                                 "1200-1400", "1400-1600", "1600-1800"))

saveRDS(stats_all_elevation, "data/processed/gap_features/stats_all_elevation.rds")


#------- calculate gap stats per aspect

stats_aspect_2009<- Gap_Stats_aspect(gap.stack$gaps9, chm9, aspect, "2009")
stats_aspect_2017 <- Gap_Stats_aspect(gap.stack$gaps17, chm17, aspect, "2017")
stats_aspect_2021 <- Gap_Stats_aspect(gap.stack$gaps21, chm21, aspect, "2021")

stats_all_aspect <- rbind(stats_aspect_2009, stats_aspect_2017, stats_aspect_2021)
# stats_all_aspect$pa_ratio <- stats_all_aspect$perimeter/stats_all_aspect$gap_area #calculate perimeter/area ratio
# stats_all_aspect$shape.index <- stats_all_aspect$perimeter/ (2*(pi*stats_all_aspect$gap_area)^0.5)
stats_all_aspect$gap_area_ha <- stats_all_aspect$gap_area/10000 #convert area to ha
#drop gaps < 400 m2 (emerge through masking of research area, as large gaps transcending management area or elevation lines get cut)
stats_all_aspect <- stats_all_aspect[stats_all_aspect$gap_area >= 400,]

saveRDS(stats_all_aspect, "data/processed/gap_features/stats_all_aspect.rds")



# --- summary 

  
gap_stats_summary <- stats_all_ftype %>% group_by(year)%>%
  summarise(n_gaps = n_distinct(gap_id),
            gap_area_total = round(sum(gap_area)/10000,2),
            max_gap_size = max(gap_area)/10000,
            mean_gap_area = mean(gap_area)/10000)

