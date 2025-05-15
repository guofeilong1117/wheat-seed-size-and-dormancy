###Extracting climate variables
library(raster)
setwd("F:/wheat/climate/wc2.1_2.5m_bio")
temp1 <- raster("wc2.1_2.5m_bio_1.tif")
temp12 <- raster("wc2.1_2.5m_bio_12.tif")
location <- read.table("F:/wheat/climate/wheat-lot.txt", sep="\t",header=T, stringsAsFactors = F)
head (location)
lon <- location$lon
lat <- location$lat
samples <- data.frame(lon, lat)
temp.data <- samples 
temp.data$BIO1 <- extract(temp1, samples)
temp.data$BIO12 <- extract(temp12, samples)
head(temp.dat
write.table(temp.data, "F:/wheat/climate/bioclim.txt", sep="\t", row.names = F,quote=F)

##### Cluster according to environmental variables
library(factoextra)
library(dplyr)
library(pacman)
library(cluster)
bc<-read.csv("F:/wheat/climate/wheat.csv",sep=',',header=TRUE)
names(bc)
bc.scaled<-scale(bc[2:3])
d<-dist(bc.scaled)
km <- eclust(bc[2:3], "kmeans", nstart = 25)
fviz_gap_stat(km$gap_stat) 
fviz_silhouette(km) 
set.seed(666)
kmeans1<-kmeans(bc.scaled,centers=3,nstart = 25)
fviz_cluster(object=kmeans1,data=bc[2:3],
             ellipse.type = "euclid",star.plot=T,repel=T,
             geom = ("point"),palette='jco',main="",
             ggtheme=theme_minimal())+
  theme(axis.title = element_blank())
summary(kmeans1)
kmeans1$cluster
kmeans1$size
aaa <- data.frame(bc[1:20], kmeans1$cluster)
aaa<-arrange(aaa,kmeans1.cluster)
write.table(aaa, "F:/wheat/climate/pop.xls", sep="\t", row.names = F,quote=F)
