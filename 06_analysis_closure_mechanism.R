###############################################################################
#
# Analyze gap closure (rate & mechanism)
#
#
##############################################################################


library(dplyr)
library(tidyr)
library(terra)
library(ggplot2)
library(RColorBrewer)


# --- prepare closure mechanism dataframes per time step ---


# prepare layers to limit analysis to core zone, below 1800m and with forest type information

# --- load NP information 

foresttype <- rast("data/processed/environment_features/forest_type2020_reclass_1m.tif")
management <- vect("data/raw/npb_zonierung_22_epsg25832.shp")
aspect<-  rast("data/processed/environment_features/aspect_2021_classified_1m.tif")
elevation.below1800 <- rast("data/processed/environment_features/elevation_below1800_200steps.tif")

# exclude management zone
core.zone <- subset(management, management$zone_id == 4, c(1:2))


# --- 2009 - 2017 ---

# --- merge closure mechanism with gaps

gap_closure_mechanism917 <- rast( "data/processed/closure/gap_closure_mechanism917.tif")
gaps2009 <- rast("data/processed/gaps_final/berchtesgaden_2009_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")

# crop and stack

gaps2009 <- crop(gaps2009, gap_closure_mechanism917)
gap_closure_mechanism_stack <- c(gap_closure_mechanism917, gaps2009)

# crop environmental features to same extent

elevation.below1800 <- crop(elevation.below1800, gap_closure_mechanism917 )
foresttype <- crop(foresttype, gap_closure_mechanism917 )

# mask stack to research area

gap_closure_mechanism_stack <- mask(gap_closure_mechanism_stack, elevation.below1800)
gap_closure_mechanism_stack <- mask(gap_closure_mechanism_stack, foresttype)
gap_closure_mechanism_stack <- mask(gap_closure_mechanism_stack, core.zone)

# - convert stack to a data frame for further processing and analysis

gap_closure_mechanism_stack.df <- as.data.frame(gap_closure_mechanism_stack, na.rm=FALSE)

# exclude pixels without gap (and hence closure):

gap_closure_mechanism_stack.df <- gap_closure_mechanism_stack.df[rowSums(is.na(gap_closure_mechanism_stack.df)) != ncol(gap_closure_mechanism_stack.df), ]
names(gap_closure_mechanism_stack.df) <- c("closure_mechanism", "gap_id")

# # save data frame in case RAM runs full
# saveRDS(gap_closure_mechanism_stack.df,"data/processed/closure/gap_closure_mechanism_pergap_917.rds" )
 gap_closure_mechanism_stack.df <- readRDS("data/processed/closure/gap_closure_mechanism_pergap_917.rds")

# aggregate closure and gap information 

gap_clo_per_id <-  gap_closure_mechanism_stack.df %>% group_by(gap_id) %>% 
  count(closure_mechanism) %>% # get amount of closure pixels per mechanism
  mutate(gap_area = sum(n))

#drop gaps < 400 m2 (emerge through masking of research area, as large gaps transcending management area or elevation lines get cut)
gap_clo_per_id <- gap_clo_per_id[gap_clo_per_id$gap_area >= 400,]


#calculate share of closure mechanism on whole gap area
gap_clo_per_id$closure_share = round(gap_clo_per_id$n/gap_clo_per_id$gap_area,2) 


# identify number of gaps not closing at all
gap_clo_per_id$contraction <- ifelse(is.na(gap_clo_per_id$closure_mechanism) & gap_clo_per_id$closure_share >= 0.99, 1,0 )
sum(gap_clo_per_id$contraction) # 5 gaps do not experience any closure from 2009-2017


# drop NAs, which is gap area which did not close - "no closure" closure area is here area, --!
# which falls within a context of closing gaps 
# - so they do not have a connection to the remaining gap area OR gather connected enough area for the gap definition

gap_clo_per_id_nona <- gap_clo_per_id %>% drop_na(closure_mechanism)  
gap_clo_per_id_nona$closure_mechanism <- as.factor(gap_clo_per_id_nona$closure_mechanism) # convert closure mechanism to factor

# recode closure mechanism
gap_clo_per_id_nona <- gap_clo_per_id_nona %>%
  mutate(gap_area_ha = round(gap_area/10000,2),
         closure_mechanism = as.factor(recode(closure_mechanism,
                                              `0`="no closure",
                                              `1`="lateral closure",
                                              `2`="vertical closure")))

# --- if I want to include the no closure pixels, you have to disable following line: !!!
# exclude no closure shares
gap_clo_per_id_nona <- subset(gap_clo_per_id_nona, closure_mechanism %in% c("vertical closure", "lateral closure"))


#aggregate closure share
gap_clo_per_id_nona <- gap_clo_per_id_nona %>% group_by(gap_id) %>%
  mutate(closure_share_sum = sum(closure_share),
         closure_area = n,
         closure_area_sum = sum(closure_area))

#gap_clo_per_id_nona$gap_area <- as.numeric(gap_clo_per_id_nona$gap_area)

# bin gap areas and get closure share per gap size bin
gap_clo_per_id_nona_917<-gap_clo_per_id_nona %>% 
  mutate(gap_area.ha = gap_area/10000,
         gap_area_bins = (cut(gap_area.ha, breaks = c(0.039,0.1,0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9,1,45))))%>% 
  mutate(gap.size = as.factor(recode(gap_area_bins,
                                     `(0.039,0.1]`="0.04-0.1",
                                     `(0.1,0.2]`="0.1-0.2",
                                     `(0.2,0.3]`="0.2-0.3",
                                     `(0.3,0.4]`="0.3-0.4",
                                     `(0.4,0.5]`="0.4-0.5",
                                     `(0.5,0.6]`="0.5-0.6",
                                     `(0.6,0.7]`="0.6-0.7",
                                     `(0.7,0.8]`="0.7-0.8",
                                     `(0.8,0.9]`="0.8-0.9",
                                     `(0.9,1]`="0.9-1",
                                     `(1,45]`=">1")))


#  --- 2017- 2021 ---


#merge closure mechanism with gaps
gap_closure_mechanism1721 <- rast("data/processed/closure/gap_closure_mechanism1721.tif")
gaps2017 <- rast("data/processed/gaps_final/berchtesgaden_2017_chm_1m_patchid_cn2cr2_mmu400n8_filtered_woheight.tif")

# crop and stack

gaps2017 <- crop(gaps2017, gap_closure_mechanism1721)
gap_closure_mechanism_stack_1721 <- c(gap_closure_mechanism1721, gaps2017)

elevation.below1800 <- crop(elevation.below1800, gap_closure_mechanism1721)
foresttype <- crop(foresttype, gap_closure_mechanism1721)
# mask stack to research area

gap_closure_mechanism_stack_1721 <- mask(gap_closure_mechanism_stack_1721, elevation.below1800)
gap_closure_mechanism_stack_1721 <- mask(gap_closure_mechanism_stack_1721, foresttype)
gap_closure_mechanism_stack_1721 <- mask(gap_closure_mechanism_stack_1721, core.zone)

gap_closure_mechanism_stack.df_1721 <- as.data.frame(gap_closure_mechanism_stack_1721, na.rm=FALSE)

#exclude pixels without gap (and hence closure):
gap_closure_mechanism_stack.df_1721 <- gap_closure_mechanism_stack.df_1721[rowSums(is.na(gap_closure_mechanism_stack.df_1721)) != ncol(gap_closure_mechanism_stack.df_1721), ]
names(gap_closure_mechanism_stack.df_1721) <- c("closure_mechanism", "gap_id")

# # save data frame in case RAM runs full
# saveRDS(gap_closure_mechanism_stack.df_1721,"processed/closure/gap_closure_mechanism_pergap_1721.rds" )
 gap_closure_mechanism_stack.df_1721 <- readRDS("data/processed/closure/gap_closure_mechanism_pergap_1721.rds" )


# aggregate closure and gap information 
gap_clo_per_id <-  gap_closure_mechanism_stack.df_1721 %>% group_by(gap_id) %>%
  count(closure_mechanism) %>% 
  mutate(gap_area = sum(n))

#drop gaps < 400 m2 (emerge through masking of research area, as large gaps transcending management area or elevation lines get cut)
gap_clo_per_id <- gap_clo_per_id[gap_clo_per_id$gap_area >= 400,]

#calculate share of closure mechanism on whole gap area
gap_clo_per_id$closure_share = round(gap_clo_per_id$n/gap_clo_per_id$gap_area,2) 

# identify number of gaps not closing
gap_clo_per_id$contraction <- ifelse(is.na(gap_clo_per_id$closure_mechanism) & gap_clo_per_id$closure_share >= 0.99, 1,0 )
sum(gap_clo_per_id$contraction) # 106 gaps do not experience any closure from 2009-2017

gaps_no_closing <- subset(gap_clo_per_id, contraction %in% 1)
range(gaps_no_closing$gap_area) # range 401-5678

# drop NAs, which is gap area which did not close - "no closure" closure area is here area, --!
# which falls within a context of closing gaps 
# - so they do not have a connection to the remaining gap area OR gather connected enough area for the gap definition
gap_clo_per_id_nona <- gap_clo_per_id %>% drop_na(closure_mechanism) #drop pixels not closing
gap_clo_per_id_nona$closure_mechanism <- as.factor(gap_clo_per_id_nona$closure_mechanism) #make closure mechanism as factor


# recode closure mechanism
gap_clo_per_id_nona <- gap_clo_per_id_nona %>%
  mutate(gap_area_ha = round(gap_area/10000,2),
         closure_mechanism = as.factor(recode(closure_mechanism,
                                              `0`="no closure", 
                                              `1`="lateral closure",
                                              `2`="vertical closure")))


# --- if I want to include the no closure pixels, you have to disable following line: !!!
# exclude no closure shares
gap_clo_per_id_nona.sub <- subset(gap_clo_per_id_nona, closure_mechanism %in% c("vertical closure", "lateral closure"))


#aggregate closure share
gap_clo_per_id_nona <- gap_clo_per_id_nona.sub %>% group_by(gap_id) %>%
  mutate(closure_share_sum = sum(closure_share),
         closure_area = n,
         closure_area_sum = sum(closure_area))


# bin gap areas and get closure share per gap size bin
gap_clo_per_id_nona_1721<-gap_clo_per_id_nona %>% 
  mutate(gap_area.ha = gap_area/10000,
         gap_area_bins = (cut(gap_area.ha, breaks = c(0.039,0.1,0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9,1,45))))%>% 
  mutate(gap.size = as.factor(recode(gap_area_bins,
                                     `(0.039,0.1]`="0.04-0.1",
                                     `(0.1,0.2]`="0.1-0.2",
                                     `(0.2,0.3]`="0.2-0.3",
                                     `(0.3,0.4]`="0.3-0.4",
                                     `(0.4,0.5]`="0.4-0.5",
                                     `(0.5,0.6]`="0.5-0.6",
                                     `(0.6,0.7]`="0.6-0.7",
                                     `(0.7,0.8]`="0.7-0.8",
                                     `(0.8,0.9]`="0.8-0.9",
                                     `(0.9,1]`="0.9-1",
                                     `(1,45]`=">1")))



# ------------merge 9-17 and 17-21 dfs for comparison -----------




gap_clo_per_id_nona_1721$timestep <- "17-21"
gap_clo_per_id_nona_917$timestep <- "9-17"


gap_clo_NP_91721 <- rbind(gap_clo_per_id_nona_917, gap_clo_per_id_nona_1721)
gap_clo_NP_91721$timestep <- factor(gap_clo_NP_91721$timestep , levels=c("9-17", "17-21"))# arrange time step labels


# -- add environmental feature information

stats_aspect <- readRDS("data/processed/gap_features/stats_all_aspect.rds")
# recode years to time steps for merges
stats_aspect <- stats_aspect %>% mutate( timestep = as.factor(recode(year,
                                                                     `2009`="9-17", 
                                                                     `2017`="17-21")))

#join with forest information
stats_ftype <- readRDS("data/processed/gap_features/stats_all_ftype.rds")
# recode years to time steps for merges
stats_ftype <- stats_ftype %>% mutate( timestep = as.factor(recode(year,
                                                                   `2009`="9-17", 
                                                                   `2017`="17-21")))

#join with management information
stats_elevation <- readRDS("data/processed/gap_features/stats_all_elevation.rds")
# recode years to time steps for merges
stats_elevation <- stats_elevation %>% mutate( timestep = as.factor(recode(year,
                                                                           `2009`="9-17", 
                                                                           `2017`="17-21")))

# merge environmental feature information with gap closure area
gap_clo_NP_91721 <- merge(x = gap_clo_NP_91721, y = stats_aspect[ , c("gap_id","aspect", "timestep")], by = c("gap_id", "timestep"), all.x=TRUE)
gap_clo_NP_91721 <- merge(x = gap_clo_NP_91721, y = stats_ftype[ , c("gap_id","forest_type", "timestep")], by = c("gap_id", "timestep"), all.x=TRUE)
gap_clo_NP_91721 <- merge(x = gap_clo_NP_91721, y = stats_elevation[ , c("gap_id","elevation", "timestep")], by = c("gap_id", "timestep"), all.x=TRUE)


# calculate annual closure rates

gap_clo_NP_91721<- gap_clo_NP_91721 %>% mutate(time = as.numeric(recode(timestep,
                                                                        `9-17`= 8, 
                                                                        `17-21`= 4)),
                                               clo_share_sum_annual = round(closure_share_sum/time,4)*100,
                                               clo_share_annual = round(closure_share/time,4)*100,
                                               clo_area_sum_annual = round(closure_area_sum/time,4),
                                               clo_area_annual = round(closure_area/time,4))

saveRDS(gap_clo_NP_91721, "data/processed/closure/gap_closure_elevation.rds") # store df to contrast gap formation and closure


# ---- relabel Larch-Swiss stone pine and Dwarf mountain pine to one forest type class ---

gap_clo_NP_91721 <- gap_clo_NP_91721 %>% mutate(forest_type = as.factor(recode(forest_type,
                                                                               `1`= "Beech",
                                                                               `2`= "Spruce-fir-beech",
                                                                               `4`= "Spruce",
                                                                               `5`= "Larch-pine")))

#label ordering
gap_clo_NP_91721$forest_type <- ordered(gap_clo_NP_91721$forest_type, levels = c("Beech", "Spruce-fir-beech","Spruce","Larch-pine"))
gap_clo_NP_91721$forest_type <- factor(gap_clo_NP_91721$forest_type,levels=rev(levels(gap_clo_NP_91721$forest_type)))

# Gap closure area summary
gap_clo_NP_91721 %>%
  summarise(sum_closure_area = sum(closure_area/10000),
            avg_clo_area_annual = mean(clo_area_sum_annual),
            sum_lateral_closure_area = sum(closure_area[closure_mechanism == "lateral closure"]/10000),
            sum_vertical_closure_area = sum(closure_area[closure_mechanism == "vertical closure"])/10000,
            share_lateral_closure = round(sum_lateral_closure_area / (sum_lateral_closure_area + sum_vertical_closure_area),2))


# Calculate the share of lateral closure on the sum of lateral and vertical closure per forest_type
lateral_share_per_forest_type <- gap_clo_NP_91721 %>%
  group_by(forest_type) %>%
  summarise(sum_lateral_closure_area = sum(closure_area[closure_mechanism == "lateral closure"]/10000),
            sum_vertical_closure_area = sum(closure_area[closure_mechanism == "vertical closure"])/10000,
            share_lateral_closure = round(sum_lateral_closure_area / (sum_lateral_closure_area + sum_vertical_closure_area),2))

# forest_type      sum_lateral_closure_area sum_vertical_closure_area share_lateral_closure

# 1 Larch-pine                          36.0                      139.                   0.21
# 2 Spruce                              24.7                       87.9                  0.22
# 3 Spruce-fir-beech                    11.2                       45.3                  0.2 
# 4 Beech                                2.78                      10.9                  0.2 



# --- append lateral + vertical closure info (overall gap closure) to main data frame  ----

gap_clo_NP_91721$id <- as.numeric(paste0(gap_clo_NP_91721$gap_id,gap_clo_NP_91721$time))

gap_clo <- as.data.frame(gap_clo_NP_91721[,c("id", "closure_mechanism",  "gap.size","aspect", "forest_type", "elevation", "clo_share_annual") ])
gap_clo$closure_mechanism <- as.character(gap_clo$closure_mechanism)
gap_clo$gap.size <- as.character(gap_clo$gap.size)
gap_clo$aspect <- as.character(gap_clo$aspect)
gap_clo$forest_type <- as.character(gap_clo$forest_type)
gap_clo$elevation <- as.character(gap_clo$elevation)

id <- as.character(unique(gap_clo$id))
#i=74
for(i in id) {
  sub <- subset(gap_clo, id %in% i)
  size <- unique(sub$gap.size)
  aspect <- unique(sub$aspect)
  ftype <- unique(sub$forest_type)
  elev <- unique(sub$elevation)
  k <- c(i,"Total", size, aspect, ftype, elev, sum(sub$clo_share_annual))
  gap_clo <- rbind(gap_clo, k)
  gap_clo$clo_share_annual <- as.numeric(gap_clo$clo_share_annual)
  gap_clo$id <- as.numeric(gap_clo$id)
}

gap_clo$closure_mechanism <- as.factor(gap_clo$closure_mechanism)
gap_clo$closure_mechanism <-  ordered(gap_clo$closure_mechanism, levels = c("lateral closure" , "vertical closure", "Total"))  
gap_clo$gap.size <- as.factor(gap_clo$gap.size)
gap_clo$aspect <- as.factor(gap_clo$aspect)
gap_clo$forest_type <- as.factor(gap_clo$forest_type)
gap_clo$elevation <- as.factor(gap_clo$elevation)

# order labels for plotting
gap_clo$forest_type <- ordered(gap_clo$forest_type, levels = c("Beech", "Spruce-fir-beech","Spruce","Larch-pine"))
gap_clo$forest_type <- factor(gap_clo$forest_type,levels=rev(levels(gap_clo$forest_type)))

gap_clo$gap.size <- ordered(gap_clo$gap.size, levels = c("0.04-0.1", "0.1-0.2",  "0.2-0.3",  "0.3-0.4",  "0.4-0.5",  "0.5-0.6",  "0.6-0.7",  "0.7-0.8",  "0.8-0.9",  "0.9-1", ">1" ))

saveRDS(gap_clo, "data/processed/closure/clo_analysis_ready.rds")



# --- Analyse gap closure ---

gap_clo <- readRDS("data/processed/closure/clo_analysis_ready.rds")


My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 28),
  axis.text.x = element_text(size = 24),
  axis.text.y = element_text(size = 24),
  axis.title.y = element_text(size = 28),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=18),
  legend.text = element_text(size=18),
  strip.text.x = element_text(size = 16),
  strip.text.y = element_text(size = 16),
  legend.position="top") #bottom


#  ---Calculate the average share of clo_share_annual for lateral closure and vertical closure in Beech forest
average_closures_share <- gap_clo %>%
  group_by(forest_type, closure_mechanism) %>%
  summarise(avg_share_closures = mean(clo_share_annual)) %>%
  pivot_wider(names_from = closure_mechanism, values_from = avg_share_closures) %>%
  mutate(share_lateral_on_total = `lateral closure`/`Total`)


# for Bechtesgaden core zone mean closure shares

# forest_type      `lateral closure` `vertical closure` Total share_lateral_on_total  --- new/updated

# 1 Larch-Pine                   0.690               2.70  3.39                  0.204
# 2 Spruce                       0.822               2.97  3.78                  0.218
# 3 Spruce-fir-beech             1.00                3.76  4.74                  0.211
# 4 Beech                        1.23                4.50  5.74                  0.215

# --- distribution of closure shares across forest types

closure_ftype_summary <- gap_clo %>%
  group_by(forest_type, closure_mechanism) %>%
  summarise(mean = mean(clo_share_annual),
            median = median(clo_share_annual),
            q2.5 = quantile(clo_share_annual, 0.025),
            q97.5 = quantile(clo_share_annual, 0.975))

gap_clo %>%
  group_by(closure_mechanism) %>%
  summarise(mean = mean(clo_share_annual),
            median = median(clo_share_annual),
            q2.5 = quantile(clo_share_annual, 0.025),
            q97.5 = quantile(clo_share_annual, 0.975))


#--- create panel view plot with closure per environmental feature ------

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 35),
  axis.text.x = element_text(size = 26,angle = 45, hjust=1),
  axis.text.y = element_text(size = 26),
  axis.title.y = element_text(size = 26),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=22),
  legend.text = element_text(size=22),
  strip.text.x = element_text(size = 25),
  strip.text.y = element_text(size = 25),
  legend.position="top") #bottom

forest_data <- select(gap_clo, -aspect) # Drop the "aspect" column
forest_data$closure_mechanism <- gsub("^(\\w)", "\\U\\1", forest_data$closure_mechanism, perl = TRUE) # change labels (capital letter first)
forest_data$closure_mechanism <- ordered(forest_data$closure_mechanism, levels = c("Lateral closure", "Vertical closure","Total")) # put labels into right order

new_df <- forest_data %>% 
  gather(category, feature, gap.size:elevation)


new_df$feature <-  ordered(new_df$feature, levels = c("0.04-0.1", "0.1-0.2",  "0.2-0.3",  "0.3-0.4",  "0.4-0.5",  "0.5-0.6",  "0.6-0.7",  "0.7-0.8",  "0.8-0.9",  "0.9-1", ">1",
                                                      "Beech", "Spruce-fir-beech","Spruce","Larch-pine",
                                                      "600-800", "800-1000","1000-1200","1200-1400", "1400-1600",  "1600-1800" ))

# Create a data frame for facet labels
facet_labels <- data.frame(category = unique(new_df$category), 
                           facet_label = rev(c("(A)", "(B)", "(C)")))

# Merge facet labels with the main data frame
new_df <- merge(new_df, facet_labels, by = "category")

tiff("data/results/gap_closure/gap_closure_mechanism_panel_box.tiff", units="in", width=16, height=8, res=300)
ggplot(new_df, aes(x = feature, y = clo_share_annual, fill = closure_mechanism)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  theme_minimal() +
  scale_fill_brewer(palette = "Dark2", name = "Closure mechanism") + 
  My_Theme +
  labs(x = "", y = "% of gap area closing annually", colour = "Closure mechanism") +
  facet_grid(~ category, scales = "free_x", 
             labeller = labeller(category = c("gap.size" = "Gap size [ha]", 
                                              "forest_type" = "Forest type", 
                                              "elevation" = "Elevation [m]"))) +
  theme(strip.text.y = element_blank(), strip.background = element_blank()) +
  geom_text(data = facet_labels, aes(x = -Inf, y = Inf, label = facet_label),
            hjust = 0, vjust = 1, size = 8, color = "black", inherit.aes = FALSE,
            nudge_x = -0.2)  
dev.off()


# --- bin elevation and gap size for better display

# Relabel the values in the feature column
new_df_agg <- new_df %>%
  mutate(feature = recode(feature,
                          "600-800" = "600-800",
                          "800-1000" = "800-1400",
                          "1000-1200" = "800-1400",
                          "1200-1400" = "800-1400",
                          "1400-1600" = "1400-1800",
                          "1600-1800" = "1400-1800",
                          "0.1-0.2" = "0.1-0.5",
                          "0.2-0.3" = "0.1-0.5",
                          "0.3-0.4" = "0.1-0.5",
                          "0.4-0.5" = "0.1-0.5",
                          "0.5-0.6" = "0.5-0.1",
                          "0.6-0.7" = "0.5-0.1",
                          "0.7-0.8" = "0.5-0.1",
                          "0.8-0.9" = "0.5-0.1",
                          "0.9-1" = "0.5-0.1"
  ))


tiff("data/results/gap_closure/gap_closure_mechanism_panel_box_aggregated.tiff", units="in", width=12, height=8, res=300)
ggplot(new_df_agg, aes(x = feature, y = clo_share_annual, fill = closure_mechanism)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  theme_minimal() +
  scale_fill_brewer(palette = "Dark2", name = "Closure mechanism") + 
  My_Theme +
  labs(x = "", y = "% of gap area closing annually", colour = "Closure mechanism") +
  facet_grid(~ category, scales = "free_x", 
             labeller = labeller(category = c("gap.size" = "Gap size [ha]", 
                                              "forest_type" = "Forest type", 
                                              "elevation" = "Elevation [m]"))) +
  theme(strip.text.y = element_blank(), strip.background = element_blank()) +
  geom_text(data = facet_labels, aes(x = -Inf, y = Inf, label = facet_label),
            hjust = 0, vjust = 1, size = 8, color = "black", inherit.aes = FALSE,
            nudge_x = -0.2)  
dev.off()


#--- panel plot with gap size and forest type without reduced bins

forest_data <- select(gap_clo, -aspect,-elevation) # Drop the "aspect" and "elevation" column


new_df <- forest_data %>% 
  gather(category, feature, gap.size:forest_type)


new_df$feature <-  ordered(new_df$feature, levels = c("0.04-0.1", "0.1-0.2",  "0.2-0.3",  "0.3-0.4",  "0.4-0.5",  "0.5-0.6",  "0.6-0.7",  "0.7-0.8",  "0.8-0.9",  "0.9-1", ">1",
                                                      "Beech", "Spruce-fir-beech","Spruce","Larch-pine"
))

tiff("data/results/gap_closure/gap_closure_mechanism_panel_box_v2.tiff", units="in", width=12, height=8, res=300)
ggplot(new_df, aes(x=feature , y=clo_share_annual, fill=closure_mechanism)) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  theme_minimal()  +  scale_fill_brewer(palette="Dark2", name = "Closure mechanism") + My_Theme +
  labs(x = "", y= "% of gap area closing annually", colour= "closure mechanism") +
  facet_grid(~category, scales="free_x", labeller = labeller(category = c("gap.size" = "Gap Size [ha]", "forest_type" = "Forest Type")))
dev.off()



#  --- single boxplots ---

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 32),
  axis.text.x = element_text(size = 26, hjust=1),
  axis.text.y = element_text(size = 26),
  axis.title.y = element_text(size = 32),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=23),
  legend.text = element_text(size=23),
  strip.text.x = element_text(size = 25),
  strip.text.y = element_text(size = 25),
  legend.position="top") 

#  --- closure share per gap size 

tiff("data/results/gap_closure/gap_closure_gap.size_box.tiff", units="in", width=12, height=8, res=300)
ggplot(subset(gap_clo, closure_mechanism %in% "Total"), aes(x=gap.size , y=clo_share_annual, fill="green")) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [ha]", y= "% of gap area closing annually", colour= "closure mechanism")+ guides(fill = FALSE)  
dev.off()


tiff("data/results/gap_closure/gap_closure_mechanism_gap.size_box.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo, aes(x=gap.size , y=clo_share_annual, fill=closure_mechanism)) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  theme_minimal()+ coord_flip()  +  scale_fill_brewer(palette="Dark2", name = "closure mechanism") + My_Theme +
  labs(x = "gap size [ha]", y= "% of gap area closing annually", colour= "closure mechanism")
dev.off()


#  --- closure share per aspect 

tiff("data/results/gap_closure/gap_closure_aspect.tiff", units="in", width=12, height=8, res=300)
ggplot(gap_clo, aes(x=clo_share_annual , y=aspect, fill= closure_mechanism)) + geom_boxplot() +
  theme_minimal()+ coord_flip()  + My_Theme +
  scale_fill_brewer(palette="Dark2", name = "closure mechanism")+
  labs(x = "% of gap area closing annually", y= "Aspect")
dev.off()


