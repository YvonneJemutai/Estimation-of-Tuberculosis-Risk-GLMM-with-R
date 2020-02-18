setwd("E:/GLMM")#define workspace
#======================================================================================
# Block 1: variable definitions, data import, preparation
#======================================================================================
library(cartools)
if (!require("pacman")) install.packages("pacman"); library(pacman) #package manager( loads required packages/libraries in list as below if not installed they will be installed
p_load(rgdal,raster,ggplot2,raster,readr,dplyr,rgdal,leaflet,CARBayes,INLA,spdep,shapefiles,sp,reshape2,data.table,leaflet)
#======================================================================================
#data loading
#======================================================================================
dataset<-read.csv("tb_data.csv")#tuberculosis count dataset
shp<-readOGR("E:/GLMM/boundary.shp")#boundary
tb_dataset<- merge(shp,dataset)#merge
#export dataset
#writeOGR(tb_dataset, layer = 'tb_dataset', 'E:/GLMM', driver="ESRI Shapefile")
train<-readOGR("E:/GLMM/Train_dataset.shp")#boundary
test<-readOGR("E:/GLMM/Test_dataset.shp")#boundary
#=======================================================================================
#======================================================================================
#Conventional SIR
#======================================================================================
#calculate IR
train$IR<-(sum(train$Observd)/sum(train$Popultn))
train$E<-train$IR*train$Popultn
train$SIR<-train$Observd/train$E
#==========================================================================================
#Neighbourhood definition
#============================================================================================
xy <- coordinates(train)

wr <- poly2nb(train, row.names=train$SLNAME, queen=FALSE)#rooks
plot(train, col='gray', border='blue')
plot(wr, xy, col='red', lwd=2, add=TRUE)
title(main="Neighborhood definition using Polygon Contiguity (ROOK)")

#converting it to a format readable by inla
nb2INLA("map.adj", wr)
g <- inla.read.graph(filename = "map.adj")
train$re_u <- 1:nrow(train@data)#spatial random effect
train$re_v <- 1:nrow(train@data)#noise
#=========================================================================================
#GLMM,GLM,RANDOM EFFECTS MODELS
#=========================================================================================
#glmmm
formula <-train$Observd ~ train$Popultn+train$pvrty_p+train$Relapse+train$HIV_Pos+f(re_u, model = "besag", graph = g) + f(re_v, model = "iid")
glmm_model <- inla(formula, family = "poisson", data = train@data, E = train$E,
                   control.predictor = list(compute = TRUE),control.compute=list(cpo=TRUE,dic=TRUE))
summary(glmm_model)
train$GLMM_train <- glmm_model$summary.fitted.values[, "mean"]
test$GLMM_Test<-(-0.3632 + -0.0148*(test$Popultn)+0.0175*(test$pvrty_p)+0.0020*(test$Relapse)+0.0015*(test$HIV_Pos))
test$GLMM_Test<-(-0.3632 + -0.0148*(test$Popultn)+0.0175*(test$pvrty_p)+0.0020*(test$Relapse)+0.0015*(test$HIV_Pos))

test$trial_test <- glmm_model$summary.fitted.values[, "mean"]                 


                 
                 