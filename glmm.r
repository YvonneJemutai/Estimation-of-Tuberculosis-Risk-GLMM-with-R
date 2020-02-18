setwd("E:/GLMM")
#======================================================================================
# Block 1: variable definitions, data import, preparation
#======================================================================================
if (!require("pacman")) install.packages("pacman"); library(pacman) #package manager( loads required packages/libraries in list as below if not installed they will be installed
p_load(rgdal,raster,ggplot2,raster,readr,dplyr,rgdal,leaflet,caTools,CARBayes,INLA,spdep,shapefiles,sp,reshape2,data.table,leaflet)
#======================================================================================
#data loading
#======================================================================================
#load dataset
dataset<-read.csv("tb_data.csv")#tuberculosis count dataset
shp<-readOGR("E:/GLMM/boundary.shp")#boundary
#=======================================================================================
#=======================================================================================
#data inspection 
#=======================================================================================
plot(dataset$Observed)
d<-density(dataset$Observed)#density plot
plot(d,main="Observed Cases")
polygon(d, col="red", border="blue")
boxplot(dataset$Observed,main="Observed Cases")#boxplot
#========================================================================================
#Correlation Plot
#========================================================================================
data<-dataset %>% select_if(is.numeric)

cormat <- round(cor(data),2)
melted_cormat <- melt(cormat)
View(melted_cormat)
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
upper_tri
View(melted_cormat)
#melted_cormat <- melt(upper_tri, na.rm = TRUE)
View(melted_cormat)
#create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
#plot
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
#======================================================================================
#Conventional SIR
#======================================================================================
#calculate IR
dataset$IR<-(sum(dataset$Observed)/sum(dataset$Population))
dataset$E<-dataset$IR*dataset$Population
dataset$SIR<-dataset$Observed/dataset$E
tb_dataset<- merge(shp,dataset)#merge
#=======================================================================================
#MAPPING Conventional SIR
#========================================================================================
l <- leaflet(tb_dataset) %>% addTiles()
pal <- colorNumeric(palette = "YlOrRd", domain = tb_dataset$SIR)
l <- leaflet(tb_dataset) %>% addTiles()
pal <- colorNumeric(palette = "YlOrRd", domain = tb_dataset$SIR)
labels<-sprintf(
  "<strong>%s</strong><br/>Observed: %g<br/>Expected: %g<br/>SIR: %g<br/>Population: %g",
  tb_dataset$SLNAME, tb_dataset$Observed,tb_dataset$E,tb_dataset$SIR,tb_dataset$Population) %>% lapply(htmltools::HTML)

l %>% addPolygons(color = "grey", weight = 1, fillColor = ~pal(SIR), fillOpacity = 0.5,
                  highlightOptions = highlightOptions(weight = 4), label = labels,
                  labelOptions = labelOptions(style = list("font-weight" = "normal",
                                                           padding = "3px 8px"),
                                              textsize = "15px",
                                              direction = "auto")) %>%
  addLegend(pal = pal, values = ~SIR, opacity = 0.5, title = "SIR",
            position = "bottomright")


#==========================================================================================
#Neighbourhood definition
#============================================================================================
xy <- coordinates(tb_dataset)

wr <- poly2nb(tb_dataset, row.names=tb_dataset$SLNAME, queen=FALSE)#rooks
plot(tb_dataset, col='gray', border='blue')
plot(wr, xy, col='red', lwd=2, add=TRUE)
title(main="Neighborhood definition using Polygon Contiguity (ROOK)")



wd1 <- dnearneigh(xy, 0, 1000)#distance=1000
plot(tb_dataset, col='gray', border='blue')
plot(wd1, xy, col='red', lwd=2, add=TRUE)
title(main="Neighborhood definition using Distance,d=1km")



wd25 <- dnearneigh(xy, 0, 2500)#distance =2500
plot(tb_dataset, col='gray', border='blue')
plot(wd25, xy, col='red', lwd=2, add=TRUE)
title(main="Neighborhood definition using Distance,d=2.5km")

#converting it to a format readable by inla
nb2INLA("map.adj", wr)
g <- inla.read.graph(filename = "map.adj")
tb_dataset$re_u <- 1:nrow(tb_dataset@data)#spatial random effect
tb_dataset$re_v <- 1:nrow(tb_dataset@data)#noise

#=========================================================================================
#GLMM,GLM,RANDOM EFFECTS MODELS
#=========================================================================================

#glmmm
formula <- tb_dataset$Observed ~ tb_dataset$Population+tb_dataset$poverty_perc+tb_dataset$Relapse+tb_dataset$HIV_Pos + f(re_u, model = "besag", graph = g) + f(re_v, model = "iid")
glmm_model <- inla(formula, family = "poisson", data = tb_dataset@data, E = tb_dataset$E,
                   control.predictor = list(compute = TRUE),control.compute=list(cpo=TRUE,waic=TRUE,dic=TRUE,pit=TRUE))
summary(glmm_model)
#glm 
glm_formula <- tb_dataset$Observed ~ tb_dataset$Population+tb_dataset$poverty_perc+tb_dataset$Relapse+tb_dataset$HIV_Pos 
glm_model <- inla(glm_formula, family = "poisson", data = tb_dataset@data, E = E,
                  control.predictor = list(compute = TRUE),control.compute=list(cpo=TRUE,dic=TRUE))
summary(glm_model)
tb_dataset$GLM_RR<- glm_model$summary.fitted.values[, "mean"]
#glm  with spatial structured
glmspatial_formula <- tb_dataset$Observed ~ tb_dataset$Population+tb_dataset$poverty_perc+tb_dataset$Relapse+tb_dataset$HIV_Pos + f(re_u, model = "besag", graph = g)
glmspatial_model <- inla(glmspatial_formula, family = "poisson", data = tb_dataset@data, E = E,
                         control.predictor = list(compute = TRUE),control.compute=list(cpo=TRUE,dic=TRUE))
summary(glmspatial_model)
#====================================================================================================
#Plotting
#=====================================================================================================
l <- leaflet(tb_dataset) %>% addTiles()
pal <- colorNumeric(palette = "YlOrRd", domain = tb_dataset$GLMM_RR)
l <- leaflet(tb_dataset) %>% addTiles()
pal <- colorNumeric(palette = "YlOrRd", domain = tb_dataset$GLMM_RR)
labels<-sprintf(
  "<strong>%s</strong><br/>Observed: %g<br/>Expected: %g<br/>SIR: %g<br/>Population: %g<br/>GLMM: %g",
  tb_dataset$SLNAME, tb_dataset$Observed,tb_dataset$E,tb_dataset$SIR,tb_dataset$Population,tb_dataset$GLMM_RR) %>% lapply(htmltools::HTML)

l %>% addPolygons(color = "grey", weight = 1, fillColor = ~pal(tb_dataset$GLMM_RR), fillOpacity = 0.5,
                  highlightOptions = highlightOptions(weight = 4), label = labels,
                  labelOptions = labelOptions(style = list("font-weight" = "normal",
                                                           padding = "3px 8px"),
                                              textsize = "15px",
                                              direction = "auto")) %>%
  addLegend(pal = pal, values = ~tb_dataset$GLMM_RR, opacity = 0.5, title = "GLMM RR",
            position = "bottomright")


#====================================================================================================
#Plotting GLMM
#=====================================================================================================
l <- leaflet(tb_dataset) %>% addTiles()
pal <- colorNumeric(palette = "YlOrRd", domain = tb_dataset$GLM_RR)
l <- leaflet(tb_dataset) %>% addTiles()
pal <- colorNumeric(palette = "YlOrRd", domain = tb_dataset$GLM_RR)
labels<-sprintf(
  "<strong>%s</strong><br/>Observed: %g<br/>Expected: %g<br/>SIR: %g<br/>Population: %g<br/>GLM: %g",
  tb_dataset$SLNAME, tb_dataset$Observed,tb_dataset$E,tb_dataset$SIR,tb_dataset$Population,tb_dataset$GLMM_RR) %>% lapply(htmltools::HTML)

l %>% addPolygons(color = "grey", weight = 1, fillColor = ~pal(tb_dataset$GLMM_RR), fillOpacity = 0.5,
                  highlightOptions = highlightOptions(weight = 4), label = labels,
                  labelOptions = labelOptions(style = list("font-weight" = "normal",
                                                           padding = "3px 8px"),
                                              textsize = "15px",
                                              direction = "auto")) %>%
  addLegend(pal = pal, values = ~tb_dataset$GLMM_RR, opacity = 0.5, title = "GLM RR",
            position = "bottomright")
#====================================================================================
#Pdf plots
#====================================================================================
#Population
marginal <- inla.smarginal(glmm_model$marginals.fixed$`tb_dataset|S|Population`)
marginal <- data.frame(marginal)
ggplot(marginal, aes(x = x, y = y)) + geom_line() +
labs(x = expression(beta[1](Population)), y = "Density") +
geom_vline(xintercept = 0, col = "blue") + theme_bw()
#poverty
marginal <- inla.smarginal(glmm_model$marginals.fixed$`tb_dataset|S|poverty_perc`)
marginal <- data.frame(marginal)
ggplot(marginal, aes(x = x, y = y)) + geom_line() +
labs(x = expression(beta[2](Poverty)), y = "Density") +
geom_vline(xintercept = 0, col = "blue") + theme_bw()

#Relapse 
marginal <- inla.smarginal(glmm_model$marginals.fixed$`tb_dataset|S|Relapse`)
marginal <- data.frame(marginal)
ggplot(marginal, aes(x = x, y = y)) + geom_line() +
labs(x = expression(beta[3](Relapse)), y = "Density") +
geom_vline(xintercept = 0, col = "blue") + theme_bw()
#HIV  Positive
marginal <- inla.smarginal(glmm_model$marginals.fixed$`tb_dataset|S|HIV_Pos`)
marginal <- data.frame(marginal)
ggplot(marginal, aes(x = x, y = y)) + geom_line() +
labs(x = expression(beta[4](HIV_Positive)), y = "Density") +
geom_vline(xintercept = 0, col = "blue") + theme_bw()

#===============================================================================
#Validation
#===============================================================================
#GLMM PLots
glmm_model$waic#WAIC
glmm_model$dic#DIC
tb_dataset$CPO <-glmm_model$cpo$cpo#CPO FOR GLMM MODEL
tb_dataset$pit <-glmm_model$cpo$pit#PIT FOR GLMM MODEL
tb_dataset$failure <-glmm_model$cpo$failure#FAILURE FOR GLMM MODEL
#glm MODEL
tb_dataset$GLM_CPO <-glm_model$cpo$cpo#CPO FOR GLM MODEL
tb_dataset$GLM_pit <-glm_model$cpo$pit#PIT FOR GLM MODEL
#glm  spatial model
tb_dataset$GLMSPATIAL_CPO <-glmspatial_model$cpo$cpo#CPO FOR GLM MODEL
tb_dataset$GLMSPATIAL_pit <-glmspatial_model$cpo$pit#PIT FOR GLM MODEL
#PIT PLOTS

hist(tb_dataset$pit,breaks=50,main="",xlab="PIT value",ylab="number")#Histogram of CPO
d<-density(tb_dataset$pit)#density plot
plot(d,main="GLM CPO Values")

hist(tb_dataset$GLMSPATIAL_pit,breaks=50,main="",xlab="CPO value",ylab="number")#Histogram of CPO
d<-density(tb_dataset$GLMSPATIAL_pit)#density plot
plot(d,main="GLM SPATIAL PIT Values")


hist(tb_dataset$GLM_pit,breaks=50,main="",xlab="PIT value",ylab="number")#Histogram of CPO
d<-density(tb_dataset$GLM_pit)#density plot
plot(d,main="GLM  PIT Values")

mean(tb_dataset$pit)
median(tb_dataset$pit)


mean(tb_dataset$GLMSPATIAL_pit)
median(tb_dataset$GLMSPATIAL_pit)

mean(tb_dataset$GLM_pit)
median(tb_dataset$GLM_pit)

shapiro.test(tb_dataset$pit)
summary(tb_dataset$GLM_pit)
summary(tb_dataset$GLMSPATIAL_pit)


shapiro.test(tb_dataset$GLMSPATIAL_pit)
shapiro.test(tb_dataset$GLM_pit)
shapiro.test(tb_dataset$pit)