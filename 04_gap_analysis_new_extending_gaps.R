######################################
#
# Analyzing new and expanding Gaps
#
#####################################


library(plyr)
library(dplyr)
library(terra)
library(ggplot2)
require(scales)


# --- load classified Gap layers of new and expanding gaps ----

gaps2017.id <- rast("data/processed/gaps_final/berchtesgaden_2017_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")
gaps2021.id <- rast("data/processed/gaps_final/berchtesgaden_2021_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")

gaps2017 <- rast("data/processed/creation/gaps2017_new_extended_stable.tif")
gaps2021 <- rast("data/processed/creation/gaps2021_new_extended_stable.tif")

# crop gaps ID

gaps2017.id <-crop(gaps2017.id, gaps2017, snap="near",mask=TRUE) 
gaps2021.id <-crop(gaps2021.id, gaps2021, snap="near",mask=TRUE) 

#load expansion and closure layer and extract only expansion areas

exp_clo917 <- rast("data/processed/gap_change/formation_closure_917_cn2cr2_mmu400n8_filtered.tif") 
exp917 <- classify(exp_clo917, cbind(1, NA)) #replace 1=closure with NA to get only expansion areas

exp_clo1721 <- rast("data/processed/gap_change/formation_closure_1721_cn2cr2_mmu400n8_filtered.tif") # doesn't load - why?
exp1721 <- classify(exp_clo1721, cbind(1, NA)) #replace 1=closure with NA to get only expansion areas

# --- load NP information 

foresttype <- rast("data/processed/environment_features/forest_type2020_reclass_1m.tif")
management <- vect("data/raw/npb_zonierung_22_epsg25832.shp")
aspect<-  rast("data/processed/environment_features/aspect_2021_classified_1m.tif")
elevation.below1800 <- rast("data/processed/environment_features/elevation_below1800_200steps.tif")
closed.forest <- vect("data/raw/closed_forest_epsg25832.shp")

# extract core zone to exclude management zone
core.zone <- subset(management, management$zone_id == 4, c(1:2))


# crop all layers to same extent

gaps2017.id <- crop(gaps2017.id, elevation.below1800)
gaps2017<- crop(gaps2017, elevation.below1800)

gaps2021.id <- crop(gaps2021.id, elevation.below1800)
gaps2021<- crop(gaps2021, elevation.below1800)

# --- stack gap information and crop it

#2017

stack2017 <- c(gaps2017.id, gaps2017, exp917, foresttype , elevation.below1800, aspect)
names(stack2017) <- c("gap.id", "new_extended", "expansion", "forest_type", "elevation", "aspect")

stack2017 <- mask(stack2017, foresttype)
stack2017 <- mask(stack2017, core.zone)
stack2017 <- mask(stack2017, closed.forest)
stack2017 <- mask(stack2017, elevation.below1800)

writeRaster(stack2017, "data/processed/creation/stack.2017.all.gap.information.expansion.tif")
gap_stack_2017 <- rast("data/processed/creation/stack.2017.all.gap.information.expansion.tif")

df <- as.data.frame(gap_stack_2017, na.rm = FALSE) 

df1 <- df[!is.na(df$gap.id),] #expansion could only take place where there is a gap now, hence reduction of df to gap.id includes all expansion
# df1 <- df1[!is.na(df1$forest_type),] # exclude all areas with no forest type information
# df1 <- df1[!is.na(df1$elevation),]# exclude all areas > 1800 m (NA in this case, as it was re-coded above) 

write_rds(df1, "data/processed/creation/stack_2017_new_exp_df.rds")
df1 <- readRDS( "data/processed/creation/stack_2017_new_exp_df.rds")


#2021

stack21 <- c(gaps2021.id, gaps2021, exp1721, foresttype, elevation.below1800, aspect)
names(stack21) <- c("gap.id", "new_extended", "expansion", "forest_type", "elevation", "aspect")

stack21 <- mask(stack21, foresttype)
stack21 <- mask(stack21, core.zone)
stack21 <- mask(stack21, closed.forest)
stack21 <- mask(stack21, elevation.below1800)

writeRaster(stack21, "data/processed/creation/stack.2021.all.gap.information.expansion.tif")
gap_stack_2021 <- rast("data/processed/creation/stack.2021.all.gap.information.expansion.tif")

df <- as.data.frame(gap_stack_2021, na.rm = FALSE) 

df2 <- df[!is.na(df$gap.id),] #expansion could only take place where there is a gap now, hence reduction of df to gap.id includes all expansion
# df2 <- df2[!is.na(df2$forest_type),] # exclude all areas with no forest type information
# df2 <- df2[!is.na(df2$elevation),]# exclude all areas > 1800 m (na in this case, as it was recoded above) 

write_rds(df2, "data/processed/creation/stack_2021_new_exp_df.rds")
df2<- readRDS("data/processed/creation/stack_2021_new_exp_df.rds")


# --- calculate features per gap(.id)

df1 <- readRDS( "data/processed/creation/stack_2017_new_exp_df.rds")
df2<- readRDS("data/processed/creation/stack_2021_new_exp_df.rds")


gap_features_917 <- df1 %>% group_by(gap.id) %>%
  summarise(area.ha = n()/10000,
            exp.area.ha = (sum(expansion, na.rm = T)/2)/10000,
            exp.share = round(exp.area.ha/area.ha,3),
            new.exp = unique(new_extended)) %>% #add aspect, elevation and forest type
  mutate(new.exp = recode(new.exp, 
                          '0' = "new", 
                          '1' = "expanding",
                          '2' = "stable"))

gap_features_1721 <- df2 %>% group_by(gap.id) %>%
  summarise(area.ha = n()/10000,
            exp.area.ha = (sum(expansion, na.rm = T)/2)/10000,
            exp.share = round(exp.area.ha/area.ha,3),
            new.exp = unique(new_extended)) %>% #add aspect, elevation and forest type
  mutate(new.exp = recode(new.exp, 
                          '0' = "new", 
                          '1' = "expanding",
                          '2' = "stable"))


# --- functions to identify major forest type, elevation and aspect per gap and hence expansion ---

getForestType <- function(gap_chm) {
  x <-gap_chm %>% group_by(gap.id, forest_type) %>% #count pixels per ftype per gap
    summarize(count = n())
  #identify dominating forest type in gap area
  xx <- data_frame()
  for (i in unique(gap_chm$gap.id)) {
    a <- x[x$gap.id == i,]        #subset to one ID
    a <- a[order(-a$count),]  #order descending according to pixel counts   
    if(nrow(a) == 1) {        #if only one entry = one forest type assign that one to ID
      xx <- rbind(xx,a[1,])
    }
    if(nrow(a) > 1) {         #if more than one entry = several ftypes or NAs
      if(is.na(a[1,2])) {xx <- rbind(xx,a[2,])}   #if is.na == TRUE = no forest type info assign NA
      if(is.na(a[1,2]) == FALSE) {xx <- rbind(xx,a[1,])}  #if there is other ftype info assign that one to ID
    }
  }
  xx<- xx %>% mutate(forest_type = as.factor(recode(forest_type,
                                                  `1`="Beech",
                                                  `2`="Spruce-fir-beech",
                                                  `4`="Spruce",
                                                  `5`="Larch-pine")))
  return(xx)
}

#function to assign elevation class
getElevation <- function(gap_chm) {
  x <-gap_chm %>% group_by(gap.id, elevation) %>% #count pixels per elevation class per gap
    summarize(count = n())
  #identify dominating elevation in gap area
  xx <- data_frame()
  for (i in unique(x$gap.id)) {
    a <- x[x$gap.id == i,]        #subset to one ID
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


#function to assign elevation class
getAspect <- function(gap_chm) {
  x <-gap_chm %>% group_by(gap.id, aspect) %>% #count pixels per aspect class per gap
    summarize(count = n())
  #identify dominating elevation in gap area
  xx <- data_frame()
  for (i in unique(x$gap.id)) {
    a <- x[x$gap.id == i,]        #subset to one ID
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


#2017

ftype <- getForestType(df1)
elevation <- getElevation(df1)
aspect <- getAspect(df1)

gap_features_917 <- merge(gap_features_917, ftype[,c("gap.id","forest_type")], by = "gap.id", all.x = TRUE)
gap_features_917 <- merge(gap_features_917, elevation[,c("gap.id","elevation")], by = "gap.id", all.x = TRUE)
gap_features_917 <- merge(gap_features_917, aspect[,c("gap.id","aspect")], by = "gap.id", all.x = TRUE)

saveRDS(gap_features_917,"data/processed/creation/gap_features_new_expanding_917.rds")


#2021

ftype <- getForestType(df2)
elevation <- getElevation(df2)
aspect <- getAspect(df2)

gap_features_1721 <- merge(gap_features_1721, ftype[,c("gap.id","forest_type")], by = "gap.id", all.x = TRUE)
gap_features_1721 <- merge(gap_features_1721, elevation[,c("gap.id","elevation")], by = "gap.id", all.x = TRUE)
gap_features_1721 <- merge(gap_features_1721, aspect[,c("gap.id","aspect")], by = "gap.id", all.x = TRUE)

saveRDS(gap_features_1721,"data/processed/creation/gap_features_new_expanding_1721.rds")


# analyse new and expanding gaps ----------------------------------------------------------

gap_features_917 <- readRDS("data/processed/creation/gap_features_new_expanding_917.rds")
gap_features_1721 <- readRDS("data/processed/creation/gap_features_new_expanding_1721.rds")

gap_features_1721$year <- as.factor("17-21")
gap_features_917$year <- as.factor("9-17")

gap_features921 <- rbind(gap_features_917, gap_features_1721)
gap_features921 <- gap_features921[gap_features921$elevation != "1800-2000",]
gap_features921 <- gap_features921[gap_features921$area.ha >= 0.04,] #delete gaps smaller than 400m2, as they emerged out of the croping of the reserach area


# exclude stable/shrinking gaps for the analysis
gap_features921 <- subset(gap_features921, new.exp %in% c("new", "expanding")) 

gap_features921$new.exp <- as.factor(gap_features921$new.exp)
# renaming New and expansion labels
levels(gap_features921$new.exp)[levels(gap_features921$new.exp) == "new"] <- "New"
levels(gap_features921$new.exp)[levels(gap_features921$new.exp) == "expanding"] <- "Expanding"


# forest label ordering for plotting
gap_features921$forest_type <- ordered(gap_features921$forest_type, levels = c("Beech", "Spruce-fir-beech","Spruce","Larch-pine"))
gap_features921$forest_type <- factor(gap_features921$forest_type,levels=rev(levels(gap_features921$forest_type)))


# --- calculate area and share of new and expanding gaps

summary_exp_new <- gap_features921 %>%
  group_by(new.exp) %>%
  summarise(total_exp_area = sum(exp.area.ha),
            avg_exp_area = mean(exp.area.ha),
            num_observations = n())%>%
  mutate(share_expanding_gaps = (total_exp_area / sum(total_exp_area)) * 100)

summary_exp_new_time <- gap_features921 %>%
  group_by(new.exp,year) %>%
  summarise(total_exp_area = sum(exp.area.ha),
            avg_exp_area = mean(exp.area.ha),
            num_observations = n())

# ---- distance calculation between new and existing gaps ----

# convert raster gaps into polygons

# --- 2009 - 2017 ---

gaps_poly.917 <- as.polygons(gap_stack_2017$gap.id)

gaps_poly <- terra::extract(gap_stack_2017$new_extended, gaps_poly.917, fun=max, bind=TRUE)
gaps_poly_new <- gaps_poly[gaps_poly$new_extended == 0, ]
gaps_poly_rest <- gaps_poly[gaps_poly$new_extended != 0, ]

gap_distance <- distance(gaps_poly_new, gaps_poly_rest)


distance_gaps_m <- as.matrix(distance_gaps)
distance.df <- as.data.frame(distance_gaps_m)

# extract nearest gap distance
gap_dist <- as.data.frame(ncol(2), nrow(0))

for(i in 1:nrow(gap_distance)) {       # for-loop over rows
  dist_info <- c(i,sort(gap_distance[,i])[1])
  gap_dist <- rbind(gap_dist, dist_info)
  names(gap_dist) <- c("gap_id", "dist_near_gap")
}


quantile(gap_dist$dist_near_gap)
# 0%         25%         50%         75%        100% 
# 6.40   495.101      700.86      954.97      1967.067 

mean(gap_dist$dist_near_gap) # 738.7354

saveRDS(gap_dist, "data/processed/creation/dist_new_gap.917.rds")

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 30),
  axis.text.x = element_text(size = 24,angle = 45, hjust=1),
  axis.text.y = element_text(size = 24),
  axis.title.y = element_text(size = 30),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=30),
  legend.text = element_text(size=30),
  strip.text.x = element_text(size = 20),
  panel.spacing = unit(2, "lines"),
  legend.position="top")

tiff("data/results/gap_creation/new_gap_dist_917.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_dist, aes(x=dist_near_gap)) + geom_histogram(bins=100) +
  theme_classic() +
  geom_histogram(color = "#000000", fill = "lightgreen") +
  geom_vline(aes(xintercept = mean(dist_near_gap)), color = "#000000", size = 1.25) +
  geom_vline(aes(xintercept = mean(dist_near_gap) + sd(dist_near_gap)), color = "#000000", size = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = mean(dist_near_gap) - sd(dist_near_gap)), color = "#000000", size = 1, linetype = "dashed")+
  labs(x= "Distance from new gap to nearest existing gap (m)", y= "Density") + My_Theme
dev.off()

# --- 2017 - 2021 ---

gaps_poly.1721 <- as.polygons(gap_stack_2021$gap.id)

gaps_poly <- terra::extract(gap_stack_2021$new_extended, gaps_poly.1721, fun=max, bind=TRUE)
gaps_poly_new <- gaps_poly[gaps_poly$new_extended == 0, ]
gaps_poly_rest <- gaps_poly[gaps_poly$new_extended != 0, ]

gap_distance <- distance(gaps_poly_new, gaps_poly_rest)


# extract nearest gap distance
gap_dist <- as.data.frame(ncol(2), nrow(0))

for(i in 1:nrow(gap_distance)) {       # for-loop over rows
  dist_info <- c(i,sort(gap_distance[,i])[1])
  gap_dist <- rbind(gap_dist, dist_info)
  names(gap_dist) <- c("gap_id", "dist_near_gap")
}


quantile(gap_dist$dist_near_gap)
# 0%       25%       50%       75%      100% 
# 1.00000  64.05803 157.00225 257.05472 830.67021 

mean(gap_dist$dist_near_gap) # 179.6479

saveRDS(gap_dist, "data/processed/creation/dist_new_gap.1721.rds")

tiff("data/results/gap_creation/new_gap_dist_1721.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_dist, aes(x=dist_near_gap)) + geom_histogram(bins=100) +
  theme_classic() +
  geom_histogram(color = "#000000", fill = "lightgreen") +
  geom_vline(aes(xintercept = mean(dist_near_gap)), color = "#000000", size = 1.25) +
  geom_vline(aes(xintercept = mean(dist_near_gap) + sd(dist_near_gap)), color = "#000000", size = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = mean(dist_near_gap) - sd(dist_near_gap)), color = "#000000", size = 1, linetype = "dashed")+
  labs(x= "Distance from new gap to nearest existing gap (m)", y= "Density") + My_Theme
dev.off()

# --- combine both time steps

gap_dist1 <- readRDS("data/processed/creation/updated/dist_new_gap.917.rds")
gap_dist2 <- readRDS("data/processed/creation/updated/dist_new_gap.1721.rds")

gap_dist <- rbind(gap_dist1, gap_dist2)

quantile(gap_dist$dist_near_gap)
# 0%        25%        50%        75%       100% 
# 1.00000   88.13058  203.49447  417.66669 1967.06507 

mean(gap_dist$dist_near_gap) # 312.1125

tiff("data/results/gap_creation/new_gap_dist_921.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_dist, aes(x=dist_near_gap)) + geom_histogram(bins=100) +
  theme_classic() +
  geom_histogram(color = "#000000", fill = "lightgreen") +
  geom_vline(aes(xintercept = mean(dist_near_gap)), color = "#000000", size = 1.25) +
  geom_vline(aes(xintercept = mean(dist_near_gap) + sd(dist_near_gap)), color = "#000000", size = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = mean(dist_near_gap) - sd(dist_near_gap)), color = "#000000", size = 1, linetype = "dashed")+
  labs(x= "Distance from new gap to nearest existing gap (m)", y= "Density") + My_Theme
dev.off()


#------- check for randomness of new gap distance to existing gaps

# mask reserach area with existing gaps (research area = forest type, already adjusted to core zone)

gap_masked_research <- mask(gap_stack_2017.closed$forest_type, gap_stack_2017.closed$new_extended, maskvalue = 1)

# draw random points representing new gaps within the research area:

random_gaps <- spatSample(gap_masked_research, size=300,method = "random", replace=FALSE, as.points=TRUE, na.rm=TRUE )

writeRaster(gap_masked_research, "data/processed/creation/distance_newgap/sample_area_gapmasked_2017.tif", overwrite=T)
writeVector(random_gaps, "data/processed/creation/distance_newgap/random_new_gap.gpkg")

# calculate distances:

gaps_poly.917 <- as.polygons(gap_stack_2017$gap.id) 

gaps_poly <- terra::extract(gap_stack_2017$new_extended, gaps_poly.917, fun=max, bind=TRUE)

gaps_poly_rest <- gaps_poly[gaps_poly$new_extended != 0, ]

gap_distance <- distance(random_gaps, gaps_poly_rest)


# extract nearest gap distance
gap_dist <- as.data.frame(ncol(2), nrow(0))

for(i in 1:nrow(gap_distance)) {       # for-loop over rows
  dist_info <- c(i,sort(gap_distance[,i])[1])
  gap_dist <- rbind(gap_dist, dist_info)
  names(gap_dist) <- c("gap_id", "dist_near_gap")
}


quantile(gap_dist$dist_near_gap)
# 0%       25%       50%       75%      100% 
# 30.37269 313.80819 485.01804 635.44262 946.03620 

mean(gap_dist$dist_near_gap) # 477.7631

saveRDS(gap_dist, "data/processed/creation/distance_newgap/dist_new_gap.917_randomtest.rds")

tiff("data/results/gap_creation/new_gap_dist_917_randomtest.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_dist, aes(x=dist_near_gap)) + geom_histogram(bins=100) +
  theme_classic() +
  geom_histogram(color = "#000000", fill = "lightgreen") +
  geom_vline(aes(xintercept = mean(dist_near_gap)), color = "#000000", size = 1.25) +
  geom_vline(aes(xintercept = mean(dist_near_gap) + sd(dist_near_gap)), color = "#000000", size = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = mean(dist_near_gap) - sd(dist_near_gap)), color = "#000000", size = 1, linetype = "dashed")+
  labs(x= "Distance from random new gap to nearest existing gap (m)", y= "Density") + My_Theme
dev.off()

# --- end distance calculation ---


# ---- calculate annual gap creation rate ----

# -- all gaps

gap.creation <- gap_features921 %>% group_by(new.exp, year) %>%
  summarize(gap.creation.ha = sum(exp.area.ha),
            n.gaps = length(unique(gap.id)),
            ha.per.gap = round(gap.creation.ha / n.gaps, 3)) %>%
  mutate(time = recode(year,
                       `9-17`=8,
                       `17-21`=4),
         gap.creation.annual = round(gap.creation.ha/time,2)) %>%
  group_by(new.exp)%>%
  mutate(avg.gap.creation.annual = round(mean(gap.creation.annual),2),
         sd = sd(gap.creation.annual,2),
        median = round(median(gap.creation.annual),2),
        q2.5 = quantile(gap.creation.annual, 0.025),
        q5 = quantile(gap.creation.annual, 0.05),
        q95 = quantile(gap.creation.annual, 0.95),
        q97.5 = quantile(gap.creation.annual, 0.975),)

saveRDS(gap.creation, "data/processed/creation/gap_creation_final.rds" )

#--------- area scaling

# ---load  area shares for scaling
area_share_class <- readRDS("data/processed/environment_features/area_share_per_class_studyarea.rds")

gap.creation$area.scaling.factor <- 100/4000 #scaling to gap creation per 100 ha (percent!), total reserach area is 4000 ha (3999 ha)
gap.creation$avg.gap.creation.annual.scaled <- gap.creation$avg.gap.creation.annual*gap.creation$area.scaling.factor
gap.creation$sd_ascaled <- round(gap.creation$sd * gap.creation$area.scaling.factor,2)

# scale median and quantiles
gap.creation$median.scaled <- gap.creation$median*gap.creation$area.scaling.factor
gap.creation$q5_ascaled <- round(gap.creation$q5 * gap.creation$area.scaling.factor,2)
gap.creation$q95_ascaled <- round(gap.creation$q95 * gap.creation$area.scaling.factor,2)
gap.creation$q2.5_ascaled <- round(gap.creation$q2.5 * gap.creation$area.scaling.factor,2)
gap.creation$q97.5_ascaled <- round(gap.creation$q97.5 * gap.creation$area.scaling.factor,2)


# --- new and expanding gaps per forest type

gap.creation.ftype <- gap_features921 %>% group_by(new.exp, year, forest_type) %>%
  summarize(gap.creation.ha = sum(exp.area.ha),
            n.gaps = length(unique(gap.id)),
            ha.per.gap = round(gap.creation.ha / n.gaps, 3)) %>%
  mutate(time = recode(year,
                       `9-17`=8,
                       `17-21`=4),
         gap.creation.annual = round(gap.creation.ha/time,2)) %>%
  group_by(new.exp, forest_type)%>%
  mutate(avg.gap.creation.annual = round(mean(gap.creation.annual),2),
         sd = sd(gap.creation.annual),
         median = round(median(gap.creation.annual),2),
         q5 = quantile(gap.creation.annual, 0.05),
         q95 = quantile(gap.creation.annual, 0.95),
         q2.5 = quantile(gap.creation.annual, 0.025),
         q97.5 = quantile(gap.creation.annual, 0.975))


#----  area scaling median+quantiles

gap.creation.ftype.scaled <- gap.creation.ftype[,c("new.exp", "forest_type", "median", "q5", "q95", "q2.5", "q97.5")]
gap.creation.ftype.scaled <- gap.creation.ftype.scaled[!duplicated(gap.creation.ftype.scaled), ]

#merge with area share information
gap.creation.ftype.scaled <- merge(gap.creation.ftype.scaled, area_share_class[,c("class.name", "class_area_perc", "total_area", "total_area_category" )],
                                   by.x = "forest_type", by.y = "class.name", all.x = TRUE)

#scale annual gap creation by area share of subcategory to 100 ha
gap.creation.ftype.scaled$area.scaling.factor <- 100/gap.creation.ftype.scaled$total_area

gap.creation.ftype.scaled$median_ascaled <- round(gap.creation.ftype.scaled$median * gap.creation.ftype.scaled$area.scaling.factor,2)
gap.creation.ftype.scaled$q5_ascaled <- round(gap.creation.ftype.scaled$q5 * gap.creation.ftype.scaled$area.scaling.factor,2)
gap.creation.ftype.scaled$q95_ascaled <- round(gap.creation.ftype.scaled$q95 * gap.creation.ftype.scaled$area.scaling.factor,2)
gap.creation.ftype.scaled$q2.5_ascaled <- round(gap.creation.ftype.scaled$q2.5 * gap.creation.ftype.scaled$area.scaling.factor,2)
gap.creation.ftype.scaled$q97.5_ascaled <- round(gap.creation.ftype.scaled$q97.5 * gap.creation.ftype.scaled$area.scaling.factor,2)

forest_gap <- as.data.frame(gap.creation.ftype.scaled)
forest_gap <- forest_gap[,c("forest_type", "new.exp", "median_ascaled","q5_ascaled", "q95_ascaled", "q2.5_ascaled", "q97.5_ascaled") ]

ftype <- as.character(unique(forest_gap$forest_type))

# -- add up new and expanding gap area for plotting --

# change new.exp from factor level to character for the loop 
forest_gap$new.exp <- as.character(forest_gap$new.exp)

for(i in ftype) {
  sub <- subset(forest_gap, forest_type %in% i)
  k <- c(i, "Total", sum(sub$median_ascaled), sum(sub$q5_ascaled), sum(sub$q95_ascaled), 
                                                        sum(sub$q2.5_ascaled), sum(sub$q97.5_ascaled))
  forest_gap <- rbind(forest_gap, k)
  forest_gap$median_ascaled <- as.numeric(forest_gap$median_ascaled)
  forest_gap$q5_ascaled <- as.numeric(forest_gap$q5_ascaled)
  forest_gap$q95_ascaled <- as.numeric(forest_gap$q95_ascaled)
  forest_gap$q2.5_ascaled <- as.numeric(forest_gap$q2.5_ascaled)
  forest_gap$q97.5_ascaled <- as.numeric(forest_gap$q97.5_ascaled)
}


forest_gap$new.exp <- as.factor(forest_gap$new.exp)

# label ordering for forest type
forest_gap$forest_type <- ordered(forest_gap$forest_type, levels = c("Beech", "Spruce-fir-beech","Larch-pine","Spruce"))


# --- gap formation by elevation (for later contrasting formation and closure by elevation belt)

gap.creation.elevation<- gap_features921 %>% group_by(new.exp, year, elevation) %>%
  summarize(gap.creation.ha = sum(exp.area.ha),
            n.gaps = length(unique(gap.id)),
            ha.per.gap = round(gap.creation.ha / n.gaps, 3)) %>%
  mutate(time = recode(year,
                       `9-17`=8,
                       `17-21`=4),
         gap.creation.annual = round(gap.creation.ha/time,2)) %>%
  group_by(new.exp, elevation)%>%
  mutate(avg.gap.creation.annual = round(mean(gap.creation.annual),2),
         sd = sd(gap.creation.annual),
         median = round(median(gap.creation.annual),2),
         q5 = quantile(gap.creation.annual, 0.05),
         q95 = quantile(gap.creation.annual, 0.95),
         q2.5 = quantile(gap.creation.annual, 0.025),
         q97.5 = quantile(gap.creation.annual, 0.975))

gap.creation.elevation$elevation <- ordered(gap.creation.elevation$elevation, levels = c("1600-1800", "1400-1600","1200-1400","1000-1200", "800-1000",  "600-800" ))
gap.creation.elevation$elevation <- factor(gap.creation.elevation$elevation,levels=rev(levels(gap.creation.elevation$elevation)))

saveRDS(gap.creation.elevation, "data/processed/creation/gap_creation_elevation.rds")

# --------------------------------------------------- graphs


# plot distribution of new gaps and instances of gap expansion

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 26),
  axis.text.x = element_text(size = 24,angle = 45, hjust=1),
  axis.text.y = element_text(size = 24),
  axis.title.y = element_text(size = 26),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=26),
  legend.text = element_text(size=26),
  strip.text.x = element_text(size = 20),
  panel.spacing = unit(2, "lines"),
  legend.position = c(0.25, 0.8))

# only gap formation area

tiff("data/results/gap_creation/new_exp_density_creation.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_features921, aes(x=exp.area.ha, fill=factor(new.exp))) + geom_density(alpha=.5) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  guides(fill=guide_legend(title="Formation mechanism")) +
  theme_few() +
  My_Theme +  
  scale_fill_colorblind()+
  labs( x="Size of gap formation area in log10 (ha)", y ="Density")
dev.off()


# per time step 

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 26),
  axis.text.x = element_text(size = 24,angle = 45, hjust=1),
  axis.text.y = element_text(size = 24),
  axis.title.y = element_text(size = 26),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=26),
  legend.text = element_text(size=26),
  strip.text.x = element_text(size = 20),
  panel.spacing = unit(2, "lines"),
  legend.position = "top")

tiff("data/results/gap_creation/new_exp_density_formation_pertime.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_features921, aes(x=exp.area.ha, fill=factor(new.exp))) + geom_density(alpha=.5) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  guides(fill=guide_legend(title="Formation mechanism")) +
  theme_few() +
  My_Theme +  
  scale_fill_colorblind()+
  labs( x="Size of gap formation area in log10 (ha)", y ="Density") +facet_wrap(~year)
dev.off()




# range of new and expanding gaps

q <- c(0.01, 0.25, 0.5, 0.75, 0.99)

gap_features921 %>% group_by(new.exp) %>% 
  summarise(range = max(area.ha)-min(area.ha),
            mean.ha = round(mean(area.ha),2),
            quant1 = round(quantile(area.ha, probs = q[1]),3),
            quant25 = round(quantile(area.ha, probs = q[2]),3), 
            quant50 = round(quantile(area.ha, probs = q[3]),3),
            quant75 = round(quantile(area.ha, probs = q[4]),3),
            quant99 = round(quantile(area.ha, probs = q[5]),3))

gap_features921 %>% group_by(new.exp) %>% 
  summarise(range = max(exp.area.ha)-min(exp.area.ha),
            sum = sum(exp.area.ha),
            mean.ha = round(mean(exp.area.ha),2),
            quant1 = round(quantile(exp.area.ha, probs = q[1]),4),
            quant25 = round(quantile(exp.area.ha, probs = q[2]),3), 
            quant50 = round(quantile(exp.area.ha, probs = q[3]),3),
            quant75 = round(quantile(exp.area.ha, probs = q[4]),3),
            quant99 = round(quantile(exp.area.ha, probs = q[5]),3))


# plot new and expansion formation rate

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 60, margin = margin(t = 60), hjust = 0.5),
  axis.text.x = element_text(size = 60,angle = 45, hjust=1),
  axis.text.y = element_text(size = 0),
  axis.title.y = element_text(size = 40),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=30),
  legend.text = element_text(size=30),
  strip.text.x = element_text(size = 20),
  panel.spacing = unit(2, "lines"),
  legend.position = c(0.8, 0.8))

# 90th quantile data range

tiff("data/results/gap_creation/area_new_exp_90quantile.tiff", units="in", width=12, height=9, res=300)
ggplot(gap.creation, aes(x=new.exp , y=median.scaled, colour= new.exp, group= new.exp, fill=new.exp)) + 
  geom_point(shape = 21, size = 16) +
  theme_classic()+ coord_flip() +
  scale_color_manual(values = c("gray40", "#E69F00"), name = "Formation mechanism", guide = "none") +
  scale_fill_manual(values = c("gray40", "#E69F00"), name = "Formation mechanism", guide = "none") +
  My_Theme +
  labs(x = "", y=expression(atop("Rate of gap formation", "(" * ha * " " * 100 * ha^-1 * " " * year^-1 * ")"))) + 
  geom_pointrange(aes(ymin=q5_ascaled, ymax=q95_ascaled), linewidth = 4)
dev.off()



# 95th quantile data range

tiff("data/results/gap_creation/area_new_exp_95quantile.tiff", units="in", width=12, height=9, res=300)
ggplot(gap.creation, aes(x=new.exp , y=median.scaled, colour= new.exp, group= new.exp, fill=new.exp)) + 
  geom_point(shape = 21, size = 16) +
  theme_classic()+ coord_flip() +
  scale_color_manual(values = c("gray40", "#E69F00"), name = "Formation mechanism", guide = "none") +
  scale_fill_manual(values = c("gray40", "#E69F00"), name = "Formation mechanism", guide = "none") +
  My_Theme +
  #labs(x = "", y= expression( "Annual rate of gap formation \n (ha per 100 ha)"))+ 
  labs(x = "", y=expression(atop("Rate of gap formation", "(" * ha * " " * 100 * ha^-1 * " " * year^-1 * ")"))) + 
  geom_pointrange(aes(ymin=q2.5_ascaled, ymax=q97.5_ascaled), linewidth = 4)
dev.off()



# --- new and expanding gap area per forest type ----


My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 26),
  axis.text.x = element_text(size = 24,angle = 45, hjust=1),
  axis.text.y = element_text(size = 24),
  axis.title.y = element_text(size = 26),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=26),
  legend.text = element_text(size=26),
  strip.text.x = element_text(size = 20),
  panel.spacing = unit(2, "lines"),
  legend.position = c(0.75, 0.2))



forest_gap$new.exp <- ordered(forest_gap$new.exp, levels = c("New", "Expanding","Total"))


# 90th quantile data range

tiff("data/results/gap_creation/area_new_exp_ftype_90quantile.tiff", units="in", width=12, height=8, res=300)
ggplot(forest_gap, aes(x=forest_type , y=median_ascaled, fill=new.exp)) + 
  geom_point(aes(colour= new.exp), shape = 21, size = 8, position=position_dodge(width=0.7)) +
  scale_color_manual(values = c("#E69F00", "grey40", "#56B4E9"), name = "Formation mechanism")+
  scale_fill_manual(values = c("#E69F00", "grey40", "#56B4E9"), name = "Formation mechanism")+
  theme_minimal()+ coord_flip() +
  My_Theme +
  labs(x = "Forest type", y= expression("Rate of gap formation (" * ha * " " * 100 * ha^-1 * " " * year^-1 * ")"), fill= "Formation mechanism", shape = "Formation mechanism") + 
  geom_pointrange(aes(ymin=q5_ascaled, ymax=q95_ascaled, colour = new.exp), 
                  position = position_dodge(width=0.7), linewidth = 1.5)
dev.off()

# 95th quantile data range

tiff("data/results/gap_creation/area_new_exp_ftype_95quantile.tiff", units="in", width=12, height=8, res=300)
ggplot(forest_gap, aes(x=forest_type , y=median_ascaled, fill=new.exp)) + 
  geom_point(aes(colour= new.exp), shape = 21, size = 8, position=position_dodge(width=0.7)) +
  scale_color_manual(values = c("#E69F00", "grey40", "#56B4E9"), name = "Formation mechanism")+
  scale_fill_manual(values = c("#E69F00", "grey40", "#56B4E9"), name = "Formation mechanism")+
  theme_minimal()+ coord_flip() +
  My_Theme +
  labs(x = "Forest type", y= expression("Rate of gap formation (" * ha * " " * 100 * ha^-1 * " " * year^-1 * ")"), fill= "Formation mechanism", shape = "Formation mechanism") + 
  geom_pointrange(aes(ymin=q2.5_ascaled, ymax=q97.5_ascaled, colour = new.exp), 
                  position = position_dodge(width=0.7), linewidth = 1.5)
dev.off()


# --- Supporting information ---

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 28),
  axis.text.x = element_text(size = 24,angle = 45, hjust=1),
  axis.text.y = element_text(size = 24),
  axis.title.y = element_text(size = 28),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=24),
  legend.text = element_text(size=23),
  strip.text.x = element_text(size = 20),
  panel.spacing = unit(2, "lines"),
  legend.position = "top")

# --- aspect


gap.creation.aspect<- gap_features921 %>% group_by(new.exp, year, aspect) %>%
  summarize(gap.creation.ha = sum(exp.area.ha),
            n.gaps = length(unique(gap.id)),
            ha.per.gap = round(gap.creation.ha / n.gaps, 3)) %>%
  mutate(time = recode(year,
                       `9-17`=8,
                       `17-21`=4),
         gap.creation.annual = round(gap.creation.ha/time,2)) %>%
  group_by(new.exp, aspect)%>%
  mutate(avg.gap.creation.annual = round(mean(gap.creation.annual),2),
         sd = sd(gap.creation.annual),
         median = round(median(gap.creation.annual),2),
         q2.5 = quantile(gap.creation.annual, 0.025),
         q97.5 = quantile(gap.creation.annual, 0.975))

# area scaling 

gap.creation.aspect.scaled <- gap.creation.aspect[,c("new.exp", "aspect", "median", "q2.5","q97.5" )]
gap.creation.aspect.scaled <- gap.creation.aspect.scaled[!duplicated(gap.creation.aspect.scaled), ]

#merge with area share information
gap.creation.aspect.scaled <- merge(gap.creation.aspect.scaled, area_share_class[,c("class.name", "class_area_perc", "total_area", "total_area_category" )],
                                    by.x = "aspect", by.y = "class.name", all.x = TRUE)

#scale annual gap creation by area share of subcategory to 100 ha
gap.creation.aspect.scaled$area.scaling.factor <- 100/gap.creation.aspect.scaled$total_area


gap.creation.aspect.scaled$median_ascaled <- round(gap.creation.aspect.scaled$median * gap.creation.aspect.scaled$area.scaling.factor,2)
gap.creation.aspect.scaled$q2.5_ascaled <- round(gap.creation.aspect.scaled$q2.5 * gap.creation.aspect.scaled$area.scaling.factor,2)
gap.creation.aspect.scaled$q97.5_ascaled <- round(gap.creation.aspect.scaled$q97.5 * gap.creation.aspect.scaled$area.scaling.factor,2)

aspect_gap <- as.data.frame(gap.creation.aspect.scaled)
aspect_gap <- aspect_gap[,c("aspect", "new.exp", "median_ascaled","q2.5_ascaled", "q97.5_ascaled") ]


aspect <- as.character(unique(aspect_gap$aspect))
aspect_gap$new.exp <- as.character(aspect_gap$new.exp)

for(i in aspect) {
  sub <- subset(aspect_gap, aspect %in% i)
  k <- c(i, "Total", sum(sub$median_ascaled), sum(sub$q2.5_ascaled), sum(sub$q97.5_ascaled))
  aspect_gap <- rbind(aspect_gap, k)
  aspect_gap$median_ascaled <- as.numeric(aspect_gap$median_ascaled)
  aspect_gap$q2.5_ascaled <- as.numeric(aspect_gap$q2.5_ascaled)
  aspect_gap$q97.5_ascaled <- as.numeric(aspect_gap$q97.5_ascaled)
}

aspect_gap$new.exp <- as.factor(aspect_gap$new.exp)



aspect_gap$new.exp <- ordered(aspect_gap$new.exp, levels = c("New", "Expanding","Total"))


tiff("data/results/gap_creation/area_new_exp_aspect_95th.tiff", units="in", width=12, height=8, res=300)
ggplot(aspect_gap, aes(x=aspect , y=median_ascaled, colour= new.exp, group= new.exp, fill=new.exp)) + 
  geom_point(aes(colour= new.exp), shape = 21, size = 8, position=position_dodge(width=0.7)) +
  theme_minimal()+ coord_flip()  +
  scale_color_manual(values = c("#E69F00", "grey40", "#56B4E9"), name = "Formation mechanism")+
  scale_fill_manual(values = c("#E69F00", "grey40", "#56B4E9"), name = "Formation mechanism")+
  My_Theme +
  labs(x = "Aspect", y= expression("Rate of gap formation (" * ha * " " * 100 * ha^-1 * " " * year^-1 * ")"), fill= "Formation mechanism")+ 
  geom_pointrange(aes(ymin=q2.5_ascaled, ymax=q97.5_ascaled, colour = new.exp), 
                  position = position_dodge(width=0.7), linewidth = 1.5)
dev.off()



# --- elevation


# area -scaling

gap.creation.elevation.scaled <- gap.creation.elevation[,c("new.exp", "elevation", "median", "q2.5", "q97.5")]
gap.creation.elevation.scaled <- gap.creation.elevation.scaled[!duplicated(gap.creation.elevation.scaled), ]

#merge with area share information
gap.creation.elevation.scaled <- merge(gap.creation.elevation.scaled, area_share_class[,c("class.name", "class_area_perc", "total_area", "total_area_category" )],
                                       by.x = "elevation", by.y = "class.name", all.x = TRUE)

#scale annual gap creation by area share of subcategory to 100 ha
gap.creation.elevation.scaled$area.scaling.factor <- 100/gap.creation.elevation.scaled$total_area

gap.creation.elevation.scaled$median_ascaled <- round(gap.creation.elevation.scaled$median * gap.creation.elevation.scaled$area.scaling.factor,2)
gap.creation.elevation.scaled$q2.5_ascaled <- round(gap.creation.elevation.scaled$q2.5 * gap.creation.elevation.scaled$area.scaling.factor,2)
gap.creation.elevation.scaled$q97.5_ascaled <- round(gap.creation.elevation.scaled$q97.5 * gap.creation.elevation.scaled$area.scaling.factor,2)

elevation_gap <- as.data.frame(gap.creation.elevation.scaled)
elevation_gap <- elevation_gap[,c("elevation", "new.exp",  "median_ascaled","q2.5_ascaled", "q97.5_ascaled") ]

elev <- as.character(unique(elevation_gap$elevation))
elevation_gap$new.exp <- as.character(elevation_gap$new.exp)

for(i in elev) {
  sub <- subset(elevation_gap, elevation %in% i)
  k <- c(i, "Total", sum(sub$median_ascaled), sum(sub$q2.5_ascaled), sum(sub$q97.5_ascaled))
  elevation_gap <- rbind(elevation_gap, k)
  elevation_gap$median_ascaled <- as.numeric(elevation_gap$median_ascaled)
  elevation_gap$q2.5_ascaled <- as.numeric(elevation_gap$q2.5_ascaled)
  elevation_gap$q97.5_ascaled <- as.numeric(elevation_gap$q97.5_ascaled)
}

elevation_gap$new.exp <- as.factor(elevation_gap$new.exp)

elevation_gap$new.exp <- ordered(elevation_gap$new.exp, levels = c("New", "Expanding","Total"))


tiff("data/results/gap_creation/area_new_exp_elevation_95th.tiff", units="in", width=12, height=8, res=300)
ggplot(elevation_gap, aes(x=elevation , y=median_ascaled, colour= new.exp, group= new.exp, fill=new.exp)) + 
  geom_point(aes(colour= new.exp), shape = 21, size = 8, position=position_dodge(width=0.7)) +
  theme_minimal()+ coord_flip()  +
  scale_color_manual(values = c("#E69F00", "grey40", "#56B4E9"), name = "Formation mechanism")+
  scale_fill_manual(values = c("#E69F00", "grey40", "#56B4E9"), name = "Formation mechanism")+
  My_Theme +
  labs(x = "Elevation", y= expression("Rate of gap formation (" * ha * " " * 100 * ha^-1 * " " * year^-1 * ")"), fill= "Formation mechanism")+ 
  geom_pointrange(aes(ymin=q2.5_ascaled, ymax=q97.5_ascaled, colour = new.exp), 
                  position = position_dodge(width=0.7), linewidth = 1.5)
dev.off()

