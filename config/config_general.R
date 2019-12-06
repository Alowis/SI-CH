################################################################################
## file of R environment configuration
################################################################################


###########
## Paths ##
###########

dir.study = getwd()  # exemple: setwd("~/Rstudy_RFA")


if (!exists("dir.functions"))  dir.functions = paste0(dir.study,'/functions')
if (!exists("dir.input"))      dir.input     = paste0(dir.study,"/in")
if (!exists("dir.output"))     dir.output    = paste0(dir.study,"/out")
if (!exists("dir.obs"))        dir.obs       = paste0(dir.input,"/obs")
if (!exists("dir.figures"))    dir.figures   = paste0(dir.output,"/figures")
if (!exists("dir.gust"))       dir.gust      = paste0(dir.output,"/gust")

dir.create(dir.output,showWarnings = F,recursive = T)
dir.create(dir.figures,showWarnings = F,recursive = T)
dir.create(dir.gust,showWarnings = F,recursive = T)

## downloading library
## """"""""""""""""""""""""

library(igraph)
library(chron)
library(spatgraphs)
library(stats)
library(mapdata)
library(chron)
library(sp)
library(geosphere)
library(cartography)
library(chron)
library(EnvStats)
library(EnvStats)
library(geomorph)
# library(in2extRemes)
library(fExtremes)
library(extRemes)
library(ismev)
library(POT)
library(ks)
library(scales)
library(ggplot2)
library(reshape2)
library(plyr)
library(ggmap)
library(rgeos)
library(rgdal)
library(grid)
library(gridExtra)
library(lattice)
library(ggpubr)
library(dplyr)
library(texmex)
# library("fitdistrplus")
library(ismev)
library(maptools)
library(astsa)
library(evd)
library(ggmap)
library(kdevine)
library(kdecopula)
library(copula)
# library(lcopula)


# source functions directory
# """"""""""""""""""""""""""" 
for (func in list.files(dir.functions)) {
  print(func)
  try(source(paste(dir.functions,func,sep="/")))
}

# local functions
# """""""""""""""

fill.gap <- function (Sta){
  
  
  timeserie<-seq(Sta$time[1],Sta$time[length(na.omit((Sta$time)))], by = "hour") 
  dfts<-data.frame(list(time=timeserie))
  merged.data <- dplyr::full_join(dfts, Sta)
  
  return (merged.data)
}
fill.gap.2 <- function (Sta){
  
  Sta$time= as.POSIXlt(Sta$time,tz="GMT")
  Sta$time= as.POSIXct(Sta$time,tz="GMT")
  timeserie<-seq(Sta$time[1],Sta$time[length(na.omit((Sta$time)))], by = "day") 
  timeserie= strptime(timeserie ,format = "%Y-%m-%d")
  timeserie= as.POSIXct(timeserie)
  dfts<-data.frame(list(time=timeserie))
  merged.data <- dplyr::full_join(dfts, Sta)
  
  return (merged.data)
}
createScaleBar <- function(lon,lat,distanceLon,distanceLat,distanceLegend, dist.units = "km"){
  # First rectangle
  bottomRight <- gcDestination(lon = lon, lat = lat, bearing = 90, dist = distanceLon, dist.units = dist.units, model = "WGS84")
  
  topLeft <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distanceLat, dist.units = dist.units, model = "WGS84")
  rectangle <- cbind(lon=c(lon, lon, bottomRight[1,"long"], bottomRight[1,"long"], lon),
                     lat = c(lat, topLeft[1,"lat"], topLeft[1,"lat"],lat, lat))
  rectangle <- data.frame(rectangle, stringsAsFactors = FALSE)
  
  # Second rectangle t right of the first rectangle
  bottomRight2 <- gcDestination(lon = lon, lat = lat, bearing = 90, dist = distanceLon*2, dist.units = dist.units, model = "WGS84")
  rectangle2 <- cbind(lon = c(bottomRight[1,"long"], bottomRight[1,"long"], bottomRight2[1,"long"], bottomRight2[1,"long"], bottomRight[1,"long"]),
                      lat=c(lat, topLeft[1,"lat"], topLeft[1,"lat"], lat, lat))
  rectangle2 <- data.frame(rectangle2, stringsAsFactors = FALSE)
  
  # Now let's deal with the text
  onTop <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distanceLegend, dist.units = dist.units, model = "WGS84")
  onTop2 <- onTop3 <- onTop
  onTop2[1,"long"] <- bottomRight[1,"long"]
  onTop3[1,"long"] <- bottomRight2[1,"long"]
  
  legend <- rbind(onTop, onTop2, onTop3)
  legend <- data.frame(cbind(legend, text = c(0, distanceLon, distanceLon*2)), stringsAsFactors = FALSE, row.names = NULL)
  return(list(rectangle = rectangle, rectangle2 = rectangle2, legend = legend))
}

createOrientationArrow <- function(scaleBar, length, distance = 1, dist.units = "km"){
  lon <- scaleBar$rectangle2[1,1]
  lat <- scaleBar$rectangle2[1,2]
  
  # Bottom point of the arrow
  begPoint <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distance, dist.units = dist.units, model = "WGS84")
  lon <- begPoint[1,"long"]
  lat <- begPoint[1,"lat"]
  
  # Let us create the endpoint
  onTop <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = length, dist.units = dist.units, model = "WGS84")
  
  leftArrow <- gcDestination(lon = onTop[1,"long"], lat = onTop[1,"lat"], bearing = 225, dist = length/5, dist.units = dist.units, model = "WGS84")
  
  rightArrow <- gcDestination(lon = onTop[1,"long"], lat = onTop[1,"lat"], bearing = 135, dist = length/5, dist.units = dist.units, model = "WGS84")
  
  res <- rbind(
    cbind(x = lon, y = lat, xend = onTop[1,"long"], yend = onTop[1,"lat"]),
    cbind(x = leftArrow[1,"long"], y = leftArrow[1,"lat"], xend = onTop[1,"long"], yend = onTop[1,"lat"]),
    cbind(x = rightArrow[1,"long"], y = rightArrow[1,"lat"], xend = onTop[1,"long"], yend = onTop[1,"lat"]))
  
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  
  # Coordinates from which "N" will be plotted
  coordsN <- cbind(x = lon, y = (lat + onTop[1,"lat"])/2)
  
  return(list(res = res, coordsN = coordsN))
}

scaleBar <- function(lon, lat, distanceLon, distanceLat, distanceLegend, dist.unit = "km", rec.fill = "white", rec.colour = "black", rec2.fill = "black", rec2.colour = "black", legend.colour = "black", legend.size = 3, orientation = TRUE, arrow.length = 500, arrow.distance = 300, arrow.North.size = 6){
  laScaleBar <- createScaleBar(lon = lon, lat = lat, distanceLon = distanceLon, distanceLat = distanceLat, distanceLegend = distanceLegend, dist.unit = dist.unit)
  # First rectangle
  rectangle1 <- geom_polygon(data = laScaleBar$rectangle, aes(x = lon, y = lat), fill = rec.fill, colour = rec.colour,cex=0.1)
  
  # Second rectangle
  rectangle2 <- geom_polygon(data = laScaleBar$rectangle2, aes(x = lon, y = lat), fill = rec2.fill, colour = rec2.colour,cex=0.1)
  
  # Legend
  scaleBarLegend <- annotate("text", label = paste(laScaleBar$legend[,"text"], dist.unit, sep=""), x = laScaleBar$legend[,"long"], y = laScaleBar$legend[,"lat"], size = legend.size, colour = legend.colour)
  
  res <- list(rectangle1, rectangle2, scaleBarLegend)
  
  if(orientation){# Add an arrow pointing North
    coordsArrow <- createOrientationArrow(scaleBar = laScaleBar, length = arrow.length, distance = arrow.distance, dist.unit = dist.unit)
    arrow <- list(geom_segment(data = coordsArrow$res, aes(x = x, y = y, xend = xend, yend = yend)), annotate("text", label = "N", x = coordsArrow$coordsN[1,"x"], y = coordsArrow$coordsN[1,"y"], size = arrow.North.size, colour = "black"))
    res <- c(res, arrow)
  }
  return(res)
}

plotlag<- function (series1, series2,min.lag=0 , max.lag = 24, corr = TRUE, smooth = TRUE, xlab = NULL, ylab= NULL, main= NULL) 
{
  name1 = paste(deparse(substitute(series1)), "(t-", sep = "")
  name2 = paste(deparse(substitute(series2)), "(t)", sep = "")
  max.lag = as.integer(max.lag)
  m1 = max.lag + 1
  prow = 4
  pcol = 4
  a = stats::ccf(series1, series2, max.lag, plot = FALSE,na.action=na.pass)$acf
  old.par <- par(no.readonly = TRUE)
  par(mfrow = c(prow, pcol), mar = c(2.5, 4, 2.5, 1), cex.main = 1.1, font.main = 1)
  for (h in min.lag:max.lag) {
    plot(lag(series1,h), series2, main = paste(main,"(", h, " h)", sep = ""), ylab = ylab, xlab = xlab)
    if (smooth == TRUE) 
      lines(stats::lowess(ts.intersect(lag(series1, -h), 
                                       series2)[, 1], ts.intersect(lag(series1, -h), 
                                                                   series2)[, 2]), col = "red")
    if (corr == TRUE) 
      legend("topright", legend = paste0("R =",round(a[m1 - h], digits = 2)),  text.col = "blue",bty="n", cex=.6)
    on.exit(par(old.par))
  }
}


SpacializedMap=function(database="world",regions=".",...){
  #Create objects of class SpatialPolygons from geographical maps from map() function 
  #BBORGY 08/20/2014
  require(maps)
  require(sp)
  require(mapdata)
  k=map(database=database,regions=regions,plot=F,fill=T,...)
  pol=vector("list",length(unique(k$names)))
  w.start=c(1,which(is.na(k$x))+1)
  w.end=c(which(is.na(k$x))-1,length(k$x));w=1
  for(j in 1:length(unique(k$names))){
    pol[[j]]=vector("list",sum(k$names==unique(k$names)[j]))
    for(i in 1:sum(k$names==unique(k$names)[j])){
      pol[[j]][[i]]=Polygon(cbind(k$x[c(w.start[w]:w.end[w],w.start[w])],k$y[c(w.start[w]:w.end[w],w.start[w])]))
      w=w+1}
    pol[[j]]=Polygons(pol[[j]],unique(k$names)[j])}
  pol=SpatialPolygons(pol,proj4string=CRS("+init=epsg:4326"))
  return(pol)
}

cart2pol <- function(x, y)
{
  r <- sqrt(x^2 + y^2)
  t <- atan(y/x)
  
  cbind(r,t)
}
pol2cart <- function(r, t)
{
  x <- r*cos(t)
  y <- r*sin(t)
  
  cbind(x,y)
}