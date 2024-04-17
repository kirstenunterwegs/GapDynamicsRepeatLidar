###################################################

# Identifying number of and area covered by gaps 
# per environmental feature

##################################################

# --- libraries

library(dplyr)
library(tidyr)
library(terra)


# --- load layers ----

gap.stack <- rast("data/processed/gaps_final/gaps_masked.tif")

gaps2009 <- gap.stack[[1]]
gaps2017 <- gap.stack[[2]]
gaps2021 <- gap.stack[[3]]


# convert gap stack to data frame 

gap.stack.df <- as.data.frame(gap.stack)
names(gap.stack.df) <- c("gaps9", "gaps17", "gaps21")

# environmental features

environment <- rast("data/processed/environment_features/stack_environment_studyarea.tif")

foresttype <- environment[[1]]
elevation <- environment[[2]]
aspect <- environment[[3]]

# --- summarize number & area of gaps:

# Gather the data into long format
gathered_df <- gap.stack.df %>%
  pivot_longer(cols = everything(), names_to = "year", values_to = "gap_id") %>%
  filter(!if_all(.cols = everything(), is.na))%>% # Drop rows with all NAs
  count(year, gap_id) %>%
  filter(n >= 400) %>% # filter out gaps < 400m2, which emerged due to cropping of gaps layers to research area
  filter(!is.nan(gap_id))# Remove rows where gap_id is NaN

# overall gap area

gathered_df %>% group_by(year)%>%
  summarise(total_gap_area_ha = round(sum(n)/10000,2))


# Count the occurrences of each value in each column
gap_counts <- gathered_df %>%
  group_by(year) %>%
  summarize(number_gaps = length(unique(gap_id)))

# overall number of gaps
sum(gap_counts$number_gaps) 



# ---- assign Gap environmental features ----



# --- define functions -----

# ---- forest type

#function to assign forest type
getForestType <- function(gap_df) {
  x <-gap_df %>% group_by(ID, ftype) %>% #count pixels per ftype per gap
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


Gap_Stats_ftype <- function (gap_layer, foresttype, year) 
{  t <- Sys.time()

gap_list <- data.frame(terra::freq(gap_layer))
gap_list$count <- gap_list$count * terra::res(gap_layer)[1]^2
gap_list <- gap_list[!is.na(gap_list[, 1]), ]
print("convert gaps to df + area: "); print(Sys.time()-t); t <- Sys.time()

# extract raster values per gap
gap_df <- as.data.frame(c(gap_layer, foresttype), na.rm = FALSE) 
names(gap_df) <- c("ID", "ftype")

df_forest <- as.data.frame(getForestType(gap_df))
gap_list$ftype <- na.omit(df_forest$ftype)
print("get forest type: "); print(Sys.time()-t); t <- Sys.time()

gap_list$year <- as.factor(year)
gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]

colnames(gap_list) <- c("gap_id", "gap_area", "forest_type", "year")
print("finish df wrangling and labeling: "); print(Sys.time()-t); t <- Sys.time()
return(gap_list)
}


#---elevation 

#function to assign elevation class
getElevation <- function(gap_df) {
  x <-gap_df %>% group_by(ID, elevation) %>% #count pixels per elevation class per gap
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
                                                  `7`="1800-2000")))
  return(xx)
}


Gap_Stats_elevation <- function (gap_layer, dtm_class, year) 
{  t <- Sys.time()
gap_list <- data.frame(terra::freq(gap_layer))
gap_list$count <- gap_list$count * terra::res(gap_layer)[1]^2
gap_list <- gap_list[!is.na(gap_list[, 1]), ]
print("convert gaps to df + area: "); print(Sys.time()-t); t <- Sys.time()

# extract raster values per gap
gap_df <- as.data.frame(c(gap_layer, dtm_class), na.rm = FALSE)
names(gap_df) <- c("ID", "elevation")
gap_df <- gap_df[!is.na(gap_df$ID),]# delete pixels without any gap
print("extract values for gaps: "); print(Sys.time()-t); t <- Sys.time()

df_elevation <- as.data.frame(getElevation(gap_df))
gap_list$elevation <- na.omit(df_elevation$elevation)
print("get elevation: "); print(Sys.time()-t); t <- Sys.time()

gap_list$year <- as.factor(year)
gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]

colnames(gap_list) <- c("gap_id", "gap_area","elevation", "year") 
print("finish df wrangling and labeling: "); print(Sys.time()-t); t <- Sys.time()
return(gap_list)
}


# ---- aspect 

#function to assign elevation class
getAspect <- function(gap_df) {
  x <-gap_df %>% group_by(ID, aspect) %>% #count pixels per aspect class per gap
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


Gap_Stats_aspect <- function (gap_layer, ascpect_class, year) 
{  t <- Sys.time()
gap_list <- data.frame(terra::freq(gap_layer))
gap_list$count <- gap_list$count * terra::res(gap_layer)[1]^2
gap_list <- gap_list[!is.na(gap_list[, 1]), ]
print("convert gaps to df + area: "); print(Sys.time()-t); t <- Sys.time()
# extract raster values per gap
gap_df <- as.data.frame(c(gap_layer, ascpect_class), na.rm = FALSE)
names(gap_df) <- c("ID",  "aspect")
gap_df <- gap_df[!is.na(gap_df$ID),]# delete pixels without any gap
print("extract values for gaps: "); print(Sys.time()-t); t <- Sys.time()

df_aspect<- as.data.frame(getAspect(gap_df))
gap_list$aspect <- na.omit(df_aspect$aspect)
print("get aspect: "); print(Sys.time()-t); t <- Sys.time()

gap_list$year <- as.factor(year)
gap_list <- gap_list[ , !names(gap_list) %in% c("layer")]

colnames(gap_list) <- c("gap_id", "gap_area", "aspect", "year") 
print("finish df wrangling and labeling: "); print(Sys.time()-t); t <- Sys.time()
return(gap_list)
}


#------- calculate gap stats per forest type

stats_ftype_2009<- Gap_Stats_ftype(gaps2009, foresttype, "2009")
stats_ftype_2017 <- Gap_Stats_ftype(gaps2017, foresttype, "2017")
stats_ftype_2021 <- Gap_Stats_ftype(gaps2021, foresttype, "2021")


stats_all_ftype <- rbind(stats_ftype_2009, stats_ftype_2017, stats_ftype_2021)
stats_all_ftype$gap_area_ha <- stats_all_ftype$gap_area/10000 #convert area to ha
stats_all_ftype <- stats_all_ftype[stats_all_ftype$gap_area >= 400,]

saveRDS(stats_all_ftype, "data/processed/gap_features/stats_all_ftype.rds")


#------- calculate gap stats per elevation

stats_elevation_2009<- Gap_Stats_elevation(gaps2009, elevation, "2009")
stats_elevatione_2017 <- Gap_Stats_elevation(gaps2017, elevation, "2017")
stats_elevation_2021 <- Gap_Stats_elevation(gaps2021, elevation, "2021")

stats_all_elevation <- rbind(stats_elevation_2009, stats_elevatione_2017, stats_elevation_2021)
stats_all_elevation$gap_area_ha <- stats_all_elevation$gap_area/10000 #convert area to ha
stats_all_elevation <- stats_all_elevation[stats_all_elevation$gap_area >= 400,]

stats_all_elevation$elevation <- factor(stats_all_elevation$elevation , levels=c("600-800", "800-1000", "1000-1200",
                                                                                 "1200-1400", "1400-1600", "1600-1800", "1800-2000"))

saveRDS(stats_all_elevation, "data/processed/gap_features/stats_all_elevation.rds")


#------- calculate gap stats per aspect

stats_aspect_2009<- Gap_Stats_aspect(gaps2009, aspect, "2009")
stats_aspect_2017 <- Gap_Stats_aspect(gaps2017, aspect, "2017")
stats_aspect_2021 <- Gap_Stats_aspect(gaps2021, aspect, "2021")

stats_all_aspect <- rbind(stats_aspect_2009, stats_aspect_2017, stats_aspect_2021)
stats_all_aspect$gap_area_ha <- stats_all_aspect$gap_area/10000 #convert area to ha
stats_all_aspect <- stats_all_aspect[stats_all_aspect$gap_area >= 400,]

saveRDS(stats_all_aspect, "data/processed/gap_features/stats_all_aspect.rds")



# --- summary 


gap_stats_summary <- stats_all_ftype %>% group_by(year)%>%
  summarise(n_gaps = n_distinct(gap_id),
            gap_area_total = round(sum(gap_area)/10000,2),
            max_gap_size = max(gap_area)/10000,
            mean_gap_area = mean(gap_area)/10000)

