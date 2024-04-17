###############################################
#
# Classify gaps as new or expanding gaps
#
##############################################

# --- libraries

library(terra)
library(raster)
library(doParallel)
library(dplyr)


# --- load gap layers ----

gap_stack <- rast("data/processed/gaps_final/gaps_masked.tif")

gaps2009 <- gap_stack[[1]]
gaps2017 <- gap_stack[[2]]
gaps2021 <- gap_stack[[3]]

# stack gaps and process for classification

gaps.df <- as.data.frame(gap_stack, na.rm=FALSE)
names(gaps.df) <- c("gaps2009", "gaps2017", "gaps2021")
gaps.df <- gaps.df[rowSums(is.na(gaps.df)) != ncol(gaps.df), ] # delete pixels without any gap at any moment in time
gaps.df[gaps.df == "NaN"] <- 0 # replace NaN with 0 to indicate vegetation pixel
gaps.df[is.na(gaps.df)] <- 0


###  ---- identify new and extending gaps ----

# concept: In the gap ID dataframe each row correspond to same pixel in space, column indicate different years.
# The foreach loop iterates over each unique value in the respective gap year column of the gaps.df data frame.
# gapID_t1 extracts unique gap IDs from the previous year column of this subset.

# The subsequent if conditions analyze the relationship between the gap IDs between the two years.
# If there is exactly one unique gap ID in the second time step which is 0 in the previous one, it indicates a new gap 
# If there is ONE unique gap ID from the previous year that is not 0, it is considered the same or a subset of an existing gap
# If there are exactly two or more unique IDs AND one of them is 0, it is considered as an extending gap 
#                                               (as previous vegetation pixels (==0), turned into a gap)

# cluster approach for whole National Park

# --- 2009-2017 ---

cl <- makeCluster(detectCores(-1)) # use all but one core for calculations (at least one free core necessary for operating system)
registerDoParallel(cl)

class_df <- foreach(i = unique(gaps.df$gaps2017), .combine = rbind) %dopar% { #"failed_polygons", "e" #, .packages = pkgs
  
  df <- subset(gaps.df, gaps2017 == i) 
  gapID_t1 <- unique(df$gaps2009)
  
  if (length(gapID_t1) == 1 & 0 %in% gapID_t1) {gap_info <- c(i,0)} # new gap
  if (length(gapID_t1) == 1 & !(0 %in% gapID_t1)) {gap_info <- c(i,1)} # gap same or sub of existing gap
  if (length(gapID_t1) == 2 & 0 %in% gapID_t1) {gap_info <- c(i,2)} # expanding gap
  # if (length(gapID_t1) == 2 & !(0 %in% gapID_t1)) {gap_info <- c(i,3)} # gap sub of several previously existing gaps (connecting several gaps without creating new gap area) > impossible
  # if (length(gapID_t1) >2 & !(0 %in% gapID_t1) ) {gap_info <- c(i,4)} # expansion of more than 1 gap (connecting several gaps without creating new gap area) > impossible
  if (length(gapID_t1) >2 &  0 %in% gapID_t1 ) {gap_info <- c(i,5)} # expansion of more than 1 gap (connecting several gaps)
  gap_info
  
}
stopCluster(cl)
class_df <- as.data.frame(class_df)

saveRDS(class_df, "data/processed/creation/new_exp_gap_class_917.rds")


class_df_917 <- class_df
names(class_df) <- c("gap_id", "class") #rename columns
gap_class_summary <- class_df %>% group_by(class) %>%  #get number of new , stable and expanding gaps
  summarise(n = n()) %>%
  mutate(perc = round(n/sum(n),2))

# 0 = new gap
# 1 = stable or shrinking gap
# 2 = expanding gap
# 5 = connection of several gaps


#subset raster layer to only new gaps
ID_vector_stablegap <- class_df$gap_id[class_df$class == 1] # get IDs of only stable/shrinking gaps
ID_vector_extendgap <- class_df$gap_id[class_df$class > 1] # get IDs of only extending gaps
ID_vector_newgap <- class_df$gap_id[class_df$class == 0] # get IDs of only new gaps

ID_vector_replace_stable <- rep(2, length(ID_vector_stablegap)) # create replace vector for stable gap       # 2 = stable gap
ID_vector_replace_extended <- rep(1, length(ID_vector_extendgap)) # create replace vector for extending gap  # 1 = expanding gap
ID_vector_replace_new <- rep(0, length(ID_vector_newgap)) # create replace vector for new gap                # 0 = new gap

rclmat1 <- cbind(ID_vector_extendgap,ID_vector_replace_extended) # create reclassification matrix  
rclamat2 <- cbind(ID_vector_newgap,ID_vector_replace_new )                                        
rclamat3 <- cbind(ID_vector_stablegap, ID_vector_replace_stable)                                
rclmat <- rbind(rclmat1, rclamat2, rclamat3) # create full reclassification matrix based on gap ID


gaps2017_class<- classify(gaps2017, rclmat, include.lowest=TRUE)
writeRaster(gaps2017_class, "data/processed/creation/gaps2017_new_extended_stable.tif")


# --- 2017-2021 ---

cl <- makeCluster(detectCores(-1)) # use all but one core for calculations (at least one free core necessary for operating system)
registerDoParallel(cl)

class_df <- foreach(i = unique(gaps.df$gaps2021), .combine = rbind) %dopar% { #"failed_polygons", "e" #, .packages = pkgs
  
  df <- subset(gaps.df, gaps2021 == i) 
  gapID_t1 <- unique(df$gaps2017)
  
  if (length(gapID_t1) == 1 & 0 %in% gapID_t1) {gap_info <- c(i,0)} # new gap
  if (length(gapID_t1) == 1 & !(0 %in% gapID_t1)) {gap_info <- c(i,1)} # gap same or sub of existing gap
  if (length(gapID_t1) == 2 & 0 %in% gapID_t1) {gap_info <- c(i,2)} # expanding gap
  #if (length(gapID_t1) == 2 & !(0 %in% gapID_t1)) {gap_info <- c(i,3)} # gap sub of several previously existing gaps
  # if (length(gapID_t1) >2 & !(0 %in% gapID_t1) ) {gap_info <- c(i,4)} # expansion of more than 1 gap (connecting several gaps without creating new gap area) > impossible
  if (length(gapID_t1) >2 &  0 %in% gapID_t1 ) {gap_info <- c(i,5)} # expansion of more than 1 gap (connecting several gaps)
  gap_info
  
}
stopCluster(cl)
class_df <- as.data.frame(class_df)

saveRDS(class_df, "data/processed/creation/new_exp_gap_class_1721.rds")

class_df_1721 <- class_df
names(class_df) <- c("gap_id", "class") #rename columns
gap_class_summary <- class_df %>% group_by(class) %>%  #get number of new , stable and expanding gaps
  summarise(n = n()) %>%
  mutate(perc = round(n/sum(n),2))

# 0 = new gap
# 1 = stable or shrinking gap
# 2 = expanding gap
# 5 = connection of several gaps


#subset raster layer to only new gaps
ID_vector_stablegap <- class_df$gap_id[class_df$class == 1] # get IDs of only stable/shrinking gaps
ID_vector_extendgap <- class_df$gap_id[class_df$class > 1] # get IDs of only extending gaps
ID_vector_newgap <- class_df$gap_id[class_df$class == 0] # get IDs of only new gaps

ID_vector_replace_stable <- rep(2, length(ID_vector_stablegap)) # create replace vector for stable gap       # 2 = stable gap
ID_vector_replace_extended <- rep(1, length(ID_vector_extendgap)) # create replace vector for extending gap  # 1 = expanding gap
ID_vector_replace_new <- rep(0, length(ID_vector_newgap)) # create replace vector for new gap                # 0 = new gap

rclmat1 <- cbind(ID_vector_extendgap,ID_vector_replace_extended) # create reclassification matrix 
rclamat2 <- cbind(ID_vector_newgap,ID_vector_replace_new )
rclamat3 <- cbind(ID_vector_stablegap, ID_vector_replace_stable)
rclmat <- rbind(rclmat1, rclamat2, rclamat3)  # create full reclassification matrix based on gap ID

gaps2021_class<- classify(gaps2021, rclmat, include.lowest=TRUE)
writeRaster(gaps2021_class, "data/processed/creation/gaps2021_new_extended_stable.tif")



