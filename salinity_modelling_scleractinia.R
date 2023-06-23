### Scleractinia Salinity Modelling for the Marshall Islands region

## RGBIF Coral Occurrence Data
library(rgbif)

scleractinia.marshall <- occ_search(scientificName = "Scleractinia",
                                        country = "MH", limit = 200000)

## Create subsets for coordinates and countries
scleractinia.marshall.coordinates <- scleractinia.marshall[,3:4]

write.csv(scleractinia.marshall.coordinates, file = "scleractinia_marshall_coord.csv")

scleractinia.marshall.coord.genus <- scleractinia.marshall[,c(3:4,31)]
write.csv(scleractinia.marshall.coord.genus, file = "scleractinia_marshall_coord_genus.csv")

## Data cleaning
## Remove missing coordinates
scleractinia.marshall.clean <- subset(scleractinia.marshall.coord.genus, !is.na(decimalLongitude) & !is.na(decimalLatitude))

## Removing duplicates (95 duplicate rows across all columns)
# distinct() in dplyr to keep only unique/distinct rows from a data frame

# library(tibble) # better to manage large data sets

# scleractinia.marshall.clean.df <- as_tibble(scleractinia.marshall.clean)

library(dplyr)
scleractinia.marshall.unique <- distinct(scleractinia.marshall.clean.df)

## Removing missing genera
scleractinia.marshall.unique.new <- subset(scleractinia.marshall.unique, !is.na(genus))
dim(scleractinia.marshall.unique.new)
# [1] 1611    3

## Removing erroneous and fossil genera
unique.genera <- unique(scleractinia.marshall.unique.new.coordcheck.genera, incomparables = FALSE)

## Export data
write.csv(unique.genera, file = "unique_genera.csv")

# Erroneous and fossil genera were removed from this list and the final occurrence data was saved (scleractinia_marshall_unique_new.csv)

##############
## Data import

# Coral data
occ.data <- read.csv("scleractinia_marshall_unique_new.csv")

# Geographic extent of occurrences
max.lat <- ceiling(max(occ.data$decimalLatitude))
min.lat <- floor(min(occ.data$decimalLatitude))
max.lon <- ceiling(max(occ.data$decimalLongitude))
min.lon <- floor(min(occ.data$decimalLongitude))
geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))

# Environmental data
library(sp)
library(raster)

sal.126.diff.data = raster("211007T114448_210824T225103_210726_sos_CMIP6_ssp126_2071_2101_diff-0p5deg.txt")
sal.585.diff.data = raster("211007T114700_210824T225103_210726_sos_CMIP6_ssp585_2071_2101_diff-0p5deg.txt")

sal.126.mean.data = raster("210927T094733_210824T225103_210726_sos_CMIP6_ssp126_2071_2101_mean-0p5deg.txt")
sal.585.mean.data = raster("210920T124039_210824T225103_210726_sos_CMIP6_ssp585_2071_2101_mean-0p5deg.txt")

## S4 method for Raster cropping
sal.126.diff.data.cropped <- crop(sal.126.diff.data, geographic.extent)
sal.585.diff.data.cropped <- crop(sal.585.diff.data, geographic.extent)

sal.126.mean.data.cropped <- crop(sal.126.mean.data, geographic.extent)
sal.585.mean.data.cropped <- crop(sal.585.mean.data, geographic.extent)

## Prepare occurrence data
names(occ.data)
# [1] "decimalLatitude"  "decimalLongitude" "genus"

dim(occ.data)
# [1] 1587    3

plot(occ.data$decimalLongitude, occ.data$decimalLatitude, pch = 20,
     cex = 0.5, col = 'blue')

## Plot occurrences on environmental layer
pdf("occurrences_Marshall.pdf") 
plot(sal.585.diff.data.cropped, legend = FALSE, axes = FALSE, box=FALSE)
points(occ.data$decimalLongitude, occ.data$decimalLatitude, pch = 20,
       cex = 0.5, col = 'blue')
dev.off()


#######
# Model
library(dismo)

coral.presence.locations <- SpatialPoints(cbind(occ.data$decimalLongitude, occ.data$decimalLatitude))

# MaxEnt
library(rJava)
maxent() # Version

fold <- kfold(coral.presence.locations, k=5)
occ.test <- coral.presence.locations[fold == 1, ]
occ.train <- coral.presence.locations[fold != 1, ]

####################
# Fitting the models

#######################################
# Salinity change data SSP scenario 126
xm.126.diff <- maxent(sal.126.diff.data.cropped, occ.train)
plot(sal.126.diff.data.cropped)

# Plot showing importance of each variable
plot(xm.126.diff)

# Response curves
response(xm.126.diff)

# Predict to entire dataset
xm.126.diff.pred <- predict(xm.126.diff, sal.126.diff.data.cropped)
plot(xm.126.diff.pred)

# Testing the model
# Background data
bg <- randomPoints(sal.126.diff.data.cropped, 10000) #background "pseudoabsences"

# Simplest way to use 'evaluate'
e1 <- evaluate(xm.126.diff, p=occ.test, a=bg, x=sal.126.diff.data.cropped)
plot(e1, 'ROC')

#######################################
# Salinity change data SSP scenario 585
xm.585.diff <- maxent(sal.585.diff.data.cropped, occ.train)
plot(sal.585.diff.data.cropped)

# Plot showing importance of each variable
plot(xm.585.diff)

# Response curves
response(xm.585.diff)

# Predict to entire dataset
xm.585.diff.pred <- predict(xm.585.diff, sal.585.diff.data.cropped)
plot(xm.585.diff.pred)

# Testing the model
# Background data
bg <- randomPoints(sal.585.diff.data.cropped, 10000) #background "pseudoabsences"

# Simplest way to use 'evaluate'
e1 <- evaluate(xm.585.diff, p=occ.test, a=bg, x=sal.585.diff.data.cropped)
plot(e1, 'ROC')

#####################################
# Salinity mean data SSP scenario 126
xm.126.mean <- maxent(sal.126.mean.data.cropped, occ.train)
plot(sal.126.mean.data.cropped)

# Plot showing importance of each variable
plot(xm.126.mean)

# Response curves
response(xm.126.mean)

# Predict to entire dataset
xm.126.mean.pred <- predict(xm.126.mean, sal.126.mean.data.cropped)
plot(xm.126.mean.pred)

# Testing the model
# Background data
bg <- randomPoints(sal.126.mean.data.cropped, 10000) #background "pseudoabsences"

# Simplest way to use 'evaluate'
e1 <- evaluate(xm.126.mean.pred, p=occ.test, a=bg, x=sal.126.mean.data.cropped)
plot(e1, 'ROC')

#####################################
# Salinity mean data SSP scenario 585
xm.585.mean <- maxent(sal.585.mean.data.cropped, occ.train)
plot(sal.585.mean.data.cropped)

# Plot showing importance of each variable
plot(xm.585.mean)

# Response curves
response(xm.585.mean)

# Predict to entire dataset
xm.585.mean.pred <- predict(xm.585.mean, sal.585.mean.data.cropped)
plot(xm.585.mean.pred)

# Testing the model
# Background data
bg <- randomPoints(sal.585.mean.data.cropped, 10000) #background "pseudoabsences"

# Simplest way to use 'evaluate'
e1 <- evaluate(xm.585.mean, p=occ.test, a=bg, x=sal.585.mean.data.cropped)
plot(e1, 'ROC')


###########
## Plotting
library(ggplot2)

# Define colour palette
cols = colorRampPalette(c("#C13127","#EDF0C0","#5E85B8"))

####################################################################
# Plot predictions for salinity change data SSP scenario 126 (Fig 4)
png(file="sos_change_126_marshall_prediction.png", width = 1230, height = 1200, res = 300)

#postscript("sos_change_126_marshall_prediction.eps", width = 1230, height = 1200)

plot(xm.126.diff.pred, col = cols(100), xaxt = "n", yaxt = "n")

map('worldHires', fill=FALSE, add=TRUE)
map('worldHires', fill=FALSE, add=TRUE)
map('worldHires', fill=FALSE, add=TRUE)

ylab="Latitude (Degrees)"
axis(side=1,at=c(140,150,160,170,180,190),labels=c("140E","150E","160E","170E","180E","190E"))
axis(side=2,at=c(-10,0,10,20,30),labels=c("10S","0","10N","20N","30N"))
legend(cex=0.8)
cex.axis=0.8
mtext("SSP1-2.6", side=3, font =2, adj=1, cex=0.8)

points(occ.data$decimalLongitude, occ.data$decimalLatitude, pch = 21, cex = 0.8, bg = "darkgrey", col = "black")

dev.off()

####################################################################
# Plot predictions for salinity change data SSP scenario 585 (Fig 4)

png(file="sos_change_585_marshall_prediction.png", width = 1230, height = 1200, res = 300)

#postscript("sos_change_585_marshall_prediction.eps", width = 1230, height = 1200)

plot(xm.585.diff.pred, col = cols(100), xaxt = "n", yaxt = "n")

map('worldHires', fill=FALSE, add=TRUE)
map('worldHires', fill=FALSE, add=TRUE)
map('worldHires', fill=FALSE, add=TRUE)

xlab="Longitude (Degrees)"
axis(side=1,at=c(140,150,160,170,180,190),labels=c("140E","150E","160E","170E","180E","190E"))
# axis(side=2,at=c(-10,0,10,20,30),labels=c("10S","0","10N","20N","30N"))
legend(cex=0.8)
cex.axis=0.8
mtext("SSP5-8.5", side=3, font =2, adj=1, cex=0.8)

points(occ.data$decimalLongitude, occ.data$decimalLatitude, pch = 21, cex = 0.8, bg = "darkgrey", col = "black")

dev.off()

###################################################################
# Plot predictions for salinity mean data SSP scenario 126 (Fig S2)

png(file="sos_mean_126_marshall_prediction.png", width = 1230, height = 1200, res = 300)

#postscript("sos_mean_126_marshall_prediction.eps", width = 1230, height = 1200)

plot(xm.126.mean.pred, col = cols(100),xaxt = "n", yaxt = "n")

map('worldHires', fill=FALSE, add=TRUE)
map('worldHires', fill=FALSE, add=TRUE)
map('worldHires', fill=FALSE, add=TRUE)

ylab="Latitude (Degrees)"
axis(side=1,at=c(140,150,160,170,180,190),labels=c("140E","150E","160E","170E","180E","190E"))
axis(side=2,at=c(-10,0,10,20,30),labels=c("10S","0","10N","20N","30N"))
legend(cex=0.8)
cex.axis=0.8
mtext("SSP1-2.6", side=3, font =2, adj=1, cex=0.8)

points(occ.data$decimalLongitude, occ.data$decimalLatitude, pch = 21, cex = 0.8, bg = "darkgrey", col = "black")

dev.off()

###################################################################
# Plot predictions for salinity mean data SSP scenario 585 (Fig S2)

png(file="sos_mean_585_marshall_prediction.png", width = 1230, height = 1200, res = 300)

#postscript("sos_mean_585_marshall_prediction.eps", width = 1230, height = 1200)

plot(xm.585.mean.pred, col = cols(100),xaxt = "n", yaxt = "n")

map('worldHires', fill=FALSE, add=TRUE)
map('worldHires', fill=FALSE, add=TRUE)
map('worldHires', fill=FALSE, add=TRUE)

xlab="Longitude (Degrees)"
axis(side=1,at=c(140,150,160,170,180,190),labels=c("140E","150E","160E","170E","180E","190E"))
# axis(side=2,at=c(-10,0,10,20,30),labels=c("10S","0","10N","20N","30N"))
legend(cex=0.8)
cex.axis=0.8
mtext("SSP5-8.5", side=3, font =2, adj=1, cex=0.8)

points(occ.data$decimalLongitude, occ.data$decimalLatitude, pch = 21, cex = 0.8, bg = "darkgrey", col = "black")

dev.off()
