

#### LOAD LIBRARIES & IMPORT DATA----
# install.packages("Tidyverse", "ggpmisc") #install packages if needed

# Load libraries
library("tidyverse")
library(ggplot2)
library(ggpmisc)
library(polynom)

# Import data from BOLD API
dfElysiaFull <- read_tsv(file = "http://www.boldsystems.org/index.php/API_Public/specimen?taxon=Elysia&format=tsv")



#### DATA REVIEW----

# Check summary info
summary(dfElysiaFull)

# Check number of records by country
table(dfElysiaFull$country)
hist(x = table(dfElysiaFull$country), xlab = "Count of BOLD Records per Country", ylab = "Frequency (No. Countries)")
# Note that 166 records taken at japan

# Check number of records by latitude
table(dfElysiaFull$lat)
hist(x = table(dfElysiaFull$lat), xlab = "Count of BOLD Records per Lat", ylab = "Frequency (No. at Same latitude)")
# Note that 120 records were taken at latitude of 36.7674 (Central Japan)

#### NOTE: From this, we can see potential bias due to the high number of records gathered from a single location



#### CLEANING DATA----

# NOTE: original species-based examination, do not mark this block
# # Add a column that indicates how many spaces are in the species name (1 indicates standard Linnaean species name)
# dfElysiaFull$spaces <- str_count(string = dfElysiaFull$species_name, pattern = "[\\s\\.\\d]")
# # Remove all entries with non-Linnaean species names, missing species names, missing bin numbers, or missing latitudes
# dfElysiaDataAvailable <- subset(dfElysiaFull, spaces == 1 & is.na(spaces) == F & is.na(species_name) == F & is.na(bin_uri) == F & is.na(lat) == F)
# unique(dfElysiaDataAvailable$spaces) #should only be "1" left
# unique(dfElysiaDataAvailable$species_name) #should look like normal species names
# print(n=45, unique(dfElysiaFull[dfElysiaFull$spaces != 1, "species_name"])) # Double check what weird species names look like



# Remove all entries with missing bin numbers or missing latitudes
dfElysiaDataAvailable <- subset(dfElysiaFull, is.na(bin_uri) == F & is.na(lat) == F)

# Check that data looks reasonable
unique(dfElysiaDataAvailable$lat) #should look like normal latitudes
unique(dfElysiaDataAvailable$bin_uri) #should have no NAs and look like good BINs
table(abs(dfElysiaDataAvailable$lat)>=90) #latitudes should not exceed +-90 (expect FALSE)



#### DATA MANIPULATION AND PROCESSING----

# Label species based on 9 ordered of latitude
# Skip commented section below, case statements are more readable/efficient
  # if (dfElysiaDataAvailable$lat <= 10 ) {dfElysiaDataAvailable$latCat <- 10}
  # else if (dfElysiaDataAvailable$lat <= 20 ) {dfElysiaDataAvailable$latCat <- 20}...

# Assign latitude category (new column 'latCat') based on absolute distance from equator (18 categories, +5 degrees from equator each)
# Categories: 5=[0,5], 10=(5,10], 15=(10,15], 20=(15,20], ... etc
dfElysiaDataAvailable <- dfElysiaDataAvailable %>%
  mutate (latCat = case_when (
    abs(lat) <= 5 ~ 5,
    abs(lat) <= 10 ~ 10,
    abs(lat) <= 15 ~ 15,
    abs(lat) <= 20 ~ 20,
    abs(lat) <= 25 ~ 25,
    abs(lat) <= 30 ~ 30,
    abs(lat) <= 35 ~ 35,
    abs(lat) <= 40 ~ 40,
    abs(lat) <= 45 ~ 45,
    abs(lat) <= 50 ~ 50,
    abs(lat) <= 55 ~ 55,
    abs(lat) <= 60 ~ 60,
    abs(lat) <= 65 ~ 65,
    abs(lat) <= 70 ~ 70,
    abs(lat) <= 75 ~ 75,
    abs(lat) <= 80 ~ 80,
    abs(lat) <= 85 ~ 85,
    TRUE ~ 90
  ))

# Check distribution of the number of records (may have sampling bias)
table(dfElysiaDataAvailable$latCat)
hist(dfElysiaDataAvailable$latCat, xlab = "Latitudinal Distance from Equator (°)", ylab = "Num. BOLD Records with BINs", main="Elysia Record Distribution Across Latitude", breaks=14)

# Find BIN richness in each latitude by removing duplicates at each latitude (remove duplicate BINs in the same latitude)
dfElysiaBINPerLat <- dfElysiaDataAvailable %>% distinct(latCat, bin_uri)

# Check distribution (looks skewed closer to equator but high at 40)
hist(dfElysiaBINPerLat$latCat, xlab = "Latitude", ylab = "Num. Unique BINs", main="BIN Richness by Latitude", breaks=14)

### NOTE: Distribution looks like it conforms to marine distribution. (referred region/climate?)



#### PLOTTING 1: LINEAR REGRESSION BINS BY LATITUDE----

# Count the number of BINs in each latitude category
lrBIN <- dfElysiaBINPerLat %>% group_by(latCat) %>% count()
lrBIN #check the dataset

# There are no samples from 70 degrees of latitude forth, therefore we should add in these missing datapoints
lrAdd <- data.frame(latCat=c(70, 75, 80, 85, 90), n=c(0, 0, 0, 0, 0))
lrBIN <- rbind(lrBIN, lrAdd) #append new data to lrBIN
lrBIN #check the dataset (looks good)
rm(lrAdd) #clean up

# PLOTTING & MODELS
set.seed(1000)

# Set up data to be plotted in new dataframe (latitudes on the x-axis and num unique species on y-axis)
dfPlot <- data.frame("x"=lrBIN$latCat, "y"=lrBIN$n)

# Set up formula & model for polynomial fit (for comparision)
my.formula <- y ~ poly(x, 3, raw = TRUE) #fit to 2nd degree (quadratic)
m <- lm(my.formula, dfPlot) #model formula on data
my.eq <- as.character(signif(as.polynomial(coef(m)), 3)) #set equation
label.text <- paste(gsub("x", "~italic(x)", my.eq, fixed = TRUE), #format equation text
                    paste("italic(R)^2", 
                          format(summary(m)$r.squared, digits = 2), 
                          sep = "~`=`~"),
                          sep = "~~~~")
label.text #check the text label for polynomial fit

ggplot(dfPlot, aes(x, y))+
  stat_poly_line(se=FALSE, color="blue") +
  stat_poly_eq(use_label(c("eq", "R2")), color="blue", family = "serif") +
  geom_point() +
  theme_minimal() +
  geom_smooth(method = "lm", se = FALSE, formula = my.formula, colour = "red") +
  annotate(geom = "text", x = 0.2, y = 0, label = label.text, family = "serif", hjust = 0, parse = TRUE, size = 4, color="red") +
  labs(x="Latitude", y="Unique BINs", title="Species Richness Poorly Correlated with Latitude") +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold")) +
  scale_x_continuous(breaks = round(seq(min(0), max(dfPlot$x), by = 5),1))

### NOTE: we see here that both regressions exhibit an outlier at approx. 30-40 degrees aways from equator (=where we found more than 100 records for japan). This may indicate that more sampling is required in other regions to determine whether Elysia are simply poorly sampled in other areas, leading to the poor species richness.  

### NOTE: here we check the distribution of records per country in the 40 degree latitude. Japan accounts for 120 of the 140 records available. 
subset(dfElysiaDataAvailable, latCat == 40) %>% group_by(country) %>% count()

### NOTE: here we check for unique BINs solely in Japan, and we see what is likely an over-representation of species diversity due to high sampling (only 12 unique BINs of 120 total samples). We likely have this high BIN number due to abnormally high sampling in this region. Thus, we can remove Japan and try the model again.
subset(dfElysiaDataAvailable, latCat == 40 & country =="Japan") %>% group_by(bin_uri) %>% count()

### NOTE: boxplot also displays the same outlier
ggplot(dfPlot, aes(x, y)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2)+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)



#### PLOTTING 2: OUTLIER REMOVED ----

#Lets remove Japan and try again
lrJapanRemoved <- subset(dfElysiaDataAvailable, country != "Japan") %>% distinct(latCat, bin_uri) %>% group_by(latCat) %>% count()

# Add back missing data again
lrAdd <- data.frame(latCat=c(70, 75, 80, 85, 90), n=c(0, 0, 0, 0, 0))
lrJapanRemoved <- rbind(lrJapanRemoved, lrAdd) #append new data to lrBIN
lrJapanRemoved #check the dataset (looks good)
rm(lrAdd) #clean up


# Plot again (same code, see comments above)
dfPlot <- data.frame("x"=lrJapanRemoved$latCat, "y"=lrJapanRemoved$n)

# Set up formula & model for polynomial fit (for comparision)
my.formula <- y ~ poly(x, 3, raw = TRUE) #fit to 2nd degree (quadratic)
m <- lm(my.formula, dfPlot) #model formula on data
my.eq <- as.character(signif(as.polynomial(coef(m)), 3)) #set equation
label.text <- paste(gsub("x", "~italic(x)", my.eq, fixed = TRUE), #format equation text
                    paste("italic(R)^2", 
                          format(summary(m)$r.squared, digits = 2), 
                          sep = "~`=`~"),
                    sep = "~~~~")
label.text #check the text label for polynomial fit

ggplot(dfPlot, aes(x, y))+
  stat_poly_line(se=FALSE, color="blue") +
  stat_poly_eq(use_label(c("eq", "R2")), color="blue", family = "serif") +
  geom_point() +
  theme_minimal() +
  geom_smooth(method = "lm", se = FALSE, formula = my.formula, colour = "red") +
  annotate(geom = "text", x = 0.2, y = 0, label = label.text, family = "serif", hjust = 0, parse = TRUE, size = 4, color="red") +
  labs(x="Latitude", y="Unique BINs", title="Species Richness Polynomially Correlated with Latitude") +
  theme(plot.title=element_text(hjust=0.5, size=12, face="bold")) +
  scale_x_continuous(breaks = round(seq(min(0), max(dfPlot$x), by = 5),1))

# Now we see a much better polynomial fit (3rd degree)



#### PLOTTING 3: LINEAR REGRESSION BINS BY SAMPLES----

# Obtain data
lrSamples <- dfElysiaDataAvailable %>% group_by(latCat) %>% count() #count number of samples per latitude
lrAdd <- data.frame(latCat=c(70, 75, 80, 85, 90), n=c(0, 0, 0, 0, 0)) #adjust data by adding 0s for latitudes with no samples
lrSamples <- rbind(lrSamples, lrAdd) 
lrSamples$nSpecies <- lrBIN$n #add bin count to the data
lrSamples

# Plot #species by #samples as linear regression
ggplot(lrSamples, aes(n, nSpecies)) +
  stat_poly_line(se=FALSE) +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  geom_point() +
  theme_minimal() +
  labs(x="Number of Samples", y="Number of Species", title="Correlation between #Samples and #Species") +
  theme(plot.title=element_text(hjust=0.5, size=20, face="bold"))

# Remove outlier and try again
outlierRemoved <- subset(lrSamples, n < 20 ) #remove outlier
ggplot(outlierRemoved, aes(n, nSpecies)) + #graph again
  stat_poly_line(se=FALSE) +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  geom_point() +
  theme_minimal() +
  labs(x="Number of Samples", y="Number of Species", title="Correlation between #Samples and #Species") +
  theme(plot.title=element_text(hjust=0.5, size=20, face="bold"))

# Both graphs indicate that there is a good correlation between sumber of samples and the number of species. Therefore, it is reasonable that more sampling is needed to to help determine species richness across latitudes. 



#End
stop()






#### EXTRA STUFF: HEAT MAP TESTING (NOT MARKED)----
# dfBINs.by.country <- dfElysiaDataAvailable %>% group_by(country, bin_uri) %>% count(bin_uri)
# dfBINs.spread <- pivot_wider(data = dfBINs.by.country, names_from = bin_uri,values_from = n)
# ggplot()
# library(ggrare) 
# library(vegan)
# ggrare(dfBINs.by.country, step = 20, se = TRUE)
# 
# ggplot(data = dfBINs.spread, aes(x = n, y = value, group = variable)) +
#   geom_line(aes(color = variable)) 


#install.packages("maps") #if not already installed
library(ggplot2) 
library(maps) 
world_map <- map_data("world") 
ggplot(world_map, aes(long, lat, group = group)) +  
  geom_polygon(fill = "gray90", color = "gray50") +  
  coord_map("mercator") +  
  ggtitle("World Map") +  
  theme_void()


#lat and number of species
dfElysiaSpeciesPerLat %>% group_by(latCat) %>% count()

mapRegions=data.frame(lon=c(-90,90,90,-90, -90,90,90,-90, -90,90,90,-90),
                      lat=c(90,90,-90,-90, 80,80,-80,-80, 70,70,-70,-70),
                      group=c(1,1,1,1, 2,2,2,2, 3,3,3,3),
                      id=c("Section 90","Section 90","Section 90","Section 90", "Section 80","Section 80","Section 80","Section 80", "Section 70","Section 70","Section 70","Section 70"))

mapRegions2=data.frame(lon=c(-90,90,90,-90, -90,90,90,-90, -90,90,90,-90 , -90,90,90,-90),
              lat=c(10,10,-10,-10, 20,20,-20,-20, 30,30,-30,-30, 10,10,-10,-10),
              group=c(9,9,9,9, 8,8,8,8, 7,7,7,7, 6,6,6,6),
              id=c("Section 10","Section 10","Section 10","Section 10", "Section 20","Section 20","Section 20","Section 20", "Section 30","Section 30","Section 30","Section 30", "Section 40","Section 40","Section 40","Section 40"))
mapRegions=rbind(mapRegions, mapRegions2)


ggplot(mapRegions, aes(lon, lat)) + 
  geom_point(size = .25, show.legend = FALSE) +
  coord_quickmap()

ggplot(mapRegions, aes(lon, lat, group = group)) +
  geom_polygon(fill = "white", colour = "grey50") + 
  coord_quickmap()


library(maps)
library(plyr)
library(gridExtra)

h2 <- dfElysiaSpeciesPerLat %>% group_by(latCat) %>% count()


world_map <- map_data("world")
#install.packages("sf")
install.packages(c("cowplot", "googleway", "ggrepel", "ggspatial", "libwgeom", "rnaturalearth", "rnaturalearthdata"))

# Load libraries
library("ggplot2")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
# Set up data for world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Set theme
theme_set(theme_bw())

ggplot(data = world) + 
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Species Distribution Map", subtitle = paste0("By latitude from equator")) +
  geom_sf(data = mapRegions, fill = NA)


gg <- ggplot(h2)
gg <- gg + geom_sf(data = world_map, mapping = aes(fill = NA), show.legend = FALSE)
# gg <- gg + geom_map(dat=world_map, map = world_map, aes(map_id=region), 
                    # fill="white", color="#7f7f7f", linewidth=0.25)
gg <- gg +   geom_sf(data = mapRegions, fill = NA)
# gg <- gg + geom_map(map = mapRegions, aes(map_id = id, fill = "#7f7f7f"), linewidth=0.25)
gg <- gg + scale_fill_gradient(low="#fff7bc", high="#cc4c02", name="Total Cases")
gg <- gg + expand_limits(x = world_map$long, y = world_map$lat)
gg <- gg + labs(x="", y="", title="World Hotspots")
gg <- gg + theme(panel.grid=element_blank(), panel.border=element_blank())
gg <- gg + theme(axis.ticks=element_blank(), axis.text=element_blank())
gg <- gg + theme(legend.position="top")
gg


# compare species richness vs latitued (linear regression)





