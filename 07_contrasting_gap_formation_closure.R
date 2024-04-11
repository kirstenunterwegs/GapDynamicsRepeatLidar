#######################################
#
# Comparing gap formation and closure rates
# in total and by elevation zone
#
######################################


# load libaries

library(dplyr)
library(tidyr)


# load data

area_share_class <- readRDS("data/processed/environment_features/area_share_per_class_studyarea.rds")
creation <- readRDS("data/processed/creation/gap_creation_elevation.rds")
closure <- readRDS("data/processed/closure/gap_closure_elevation.rds")


#--- process data for comparison

# --- gap creation

creation <- as.data.frame(creation[,c("elevation", "new.exp", "gap.creation.ha", "time") ])
creation$gap.creation.ha <- as.numeric(creation$gap.creation.ha)
creation$new.exp <- factor(creation$new.exp, levels = c(levels(creation$new.exp), "Total")) #pre-add new factor level

elev <- as.character(unique(creation$elevation))  

for(i in elev) {
  sub <- subset(creation, elevation %in% i)
  k <- c(i, "Total", sum(sub$gap.creation.ha), "12")
  creation <- rbind(creation, k)
  creation$gap.creation.ha <- as.numeric(creation$gap.creation.ha)
}

creation <- creation %>%
  filter(new.exp == "Total") %>% 
  mutate(annual.creation.area.ha = gap.creation.ha/as.numeric(time)) %>%
  select(elevation,gap.creation.ha, annual.creation.area.ha)



# --- gap closure


closure<- closure %>% 
  filter(closure_mechanism == "vertical closure") %>% # filter to one mechanism, as I use already Total closure area (avoid double counting)
  mutate(gap.closure.ha = round(closure_area_sum/10000,4))%>%
  group_by(elevation) %>%
  summarise(gap.closure.ha = sum(gap.closure.ha))


# merge gap formation and closure
creation.closure <- merge(creation, closure, by = "elevation")

# merge with area share information
creation.closure.scaled <- merge(creation.closure, area_share_class[,c("class.name", "class_area_perc", "total_area", "total_area_category" )],
                                 by.x = "elevation", by.y = "class.name", all.x = TRUE)


# --- calculate overall gap formation and closure for research area

sum_creation_closure <- creation.closure.scaled %>% 
  summarize(creation = sum(gap.creation.ha),
            closure = sum(gap.closure.ha))%>%
  mutate(area.scaling.factor = 100/3999, 
         obs.years = 12,
         creation.ha.yr = creation*area.scaling.factor/12,
         closure.ha.yr = closure*area.scaling.factor/12)


# creation  closure area.scaling.factor obs.years creation.ha.yr closure.ha.yr
# 278.7181 357.3844          0.02500625        12      0.5808079      0.744737


# --- calculate gap formation and closure rates per ecological elevation zone

creation.closure.elev.class <- creation.closure.scaled %>% mutate(elev.class = as.factor(recode(elevation,
                                                                                                '600-800' = "submontan",
                                                                                                '800-1000' = "submontan",
                                                                                                '1000-1200' = "montan",
                                                                                                '1200-1400' = "montan",
                                                                                                '1400-1600' = "montan",
                                                                                                '1600-1800' = "sub-alpine"))) %>%
  group_by(elev.class) %>%
  summarise(gap.creation.ha = sum(gap.creation.ha),
            gap.closure.ha = sum(gap.closure.ha),
            area = sum(total_area))

# area scaling and calculation of annual rates

creation.closure.elev.class <- creation.closure.elev.class %>%
  mutate(scaling.factor = 100/area, # scale rates to 100 ha
         creation.ha.scaled = gap.creation.ha*scaling.factor,
         closure.ha.scaled = gap.closure.ha*scaling.factor,
         creation.ha.yr = creation.ha.scaled/12, # 12 years observation period
         closure.ha.yr = closure.ha.scaled/12)


