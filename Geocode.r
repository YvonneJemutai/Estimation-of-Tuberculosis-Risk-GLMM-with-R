library(readr)
library(rgdal)
library(ggplot2)
#import coordinates from khis
nrb_hospitals <- read_csv("nrb_hospitals.csv")#coordinates
tb_summary <- read_csv("tb_summary.csv")
#rename columns



#merge the coordinates with tb records
tb_coords<-merge(tb_summary,nrb_hospitals,by="Health Facility" ,all =T )

tb_coords<-na.omit(tb_coords)
View(tb_coords)
write.csv(tb_coords,"tb_coordinates.csv")


#spatial distribution of health facilities in Nairobi county
Nairobi_sublocations<- readOGR("E:/JKUAT/shapefiles/Nairobi_sublocations","Nairobi_sublocations")#folder dir,shp name
plot(Nairobi_sublocations)

coords <- tb_coords[c("Longitude", "Latitude")]

# Making sure we are working with rows that don't have any blanks
coords <- coords[complete.cases(coords),]
# Letting R know that these are specifically spatial coordinates
sp <- SpatialPoints(coords)
plot(Nairobi_sublocations)
plot(sp, col="red", pch= 21,add=TRUE)
#delete points not within the polygon
##spatial join to get counts 

