####First script which pre-processes input data######
####1. Import ncdf files
####2. Compute values of extreme thresholds
####3. Estimate parameters for DBSCAN
####4. Cluster extreme wind an precipitation
####5. Save intermediary data and metadata of single hazard cluster
rm(list=ls())  
gc()
getwd()

################Choose working directory here####################################

# setwd("C:/Users/PhD Student/OneDrive - King's College London/DATA/04_Spatiotemporal")
# setwd("C:/Users/k1638615/King's College London/OneDrive - King's College London/DATA/04_Spatiotemporal")

setwd("C:/Users/Alois/OneDrive - King's College London/DATA/04_Spatiotemporal/spatiotemporal_clustering")

##load libraries
config_file=paste0(getwd(),"/config/config_general.R")
source(config_file)
library(ncdf4)
library(dbscan)

#===========================================================================================================

############Simulation study to show DBSCAN skills#####################

#Synthetic data example for dbscan with data of the same format as the input data
library(mvtnorm)
library(dplyr)
diag(3)
generateGaussianData <- function(n, center, sigma, label) {
  data = rmvnorm(n, mean = center, sigma = sigma)
  data = data.frame(data)
  names(data) = c("x", "y","z")
  data = data %>% mutate(class=factor(label))
  data
}

dataset1 <- {
  # cluster 1
  n = 50
  center = c(-4, 58, 10)
  sigma = matrix(c(1, 0.5, 0.9, 0.5, 1, 0.3, 0.9, 0.3, 10), nrow = 3)
  
  data1 = generateGaussianData(n, center, sigma, 1)
  data1$z=round(data1$z)
  data1$x=round((data1$x + 0.25)*2) / 2 - 0.25
  data1$y=round((data1$y + 0.25)*2) / 2 - 0.25
  plot(data1)
  
  # cluster 2
  n = 1000
  center = c(0, 55, 40)
  sigma = matrix(c(1, 0.8, 0.9, 0.8, .5, 0, 0.9, 0, 1), nrow = 3)

  data2 = generateGaussianData(n, center, sigma, 2)
  data2$z=round(data2$z)
  data2$x=round((data2$x + 0.25)*2) / 2 - 0.25
  data2$y=round((data2$y + 0.25)*2) / 2 - 0.25
  plot(data2)
  # cluster 3
  n = 500
  center = c(0, 50, 150)
  sigma = matrix(c(2, 1, -10, 1, 1, 0, -10, 0, 50), nrow = 3)
  
  data3 = generateGaussianData(n, center, sigma, 3)
  data3$z=round(data3$z)
  data3$x=round((data3$x + 0.25)*2) / 2 - 0.25
  data3$y=round((data3$y + 0.25)*2) / 2 - 0.25
  plot(data3,col=data3$z)
 
  data = bind_rows(data1, data2, data3) 
  
  data$dataset = "1 - Mixture of Gaussians"
  
  data
}


dataset1 %>% ggplot(aes(x=x, y=y, color=class)) +
  geom_point() +
  coord_fixed() +
  scale_color_manual(values=c(1, 2, 3))

epcl<-2.5
coef<-4
dbdat<-dataset1[,c(1,2,3)]
dbdat$x=dataset1$x*coef
dbdat$y=dataset1$y*coef
rpip<-dbscan(dataset1[,c(1,2,3)], eps=epcl, minPts = 10)
length(unique(rpip$cluster))

dataset1<-cbind(dataset1,as.character(rpip$cluster))
names(dataset1)[6]="Cluster"
dataset1 %>% ggplot(aes(x=x, y=y, color=Cluster)) +
  geom_point() +
  coord_fixed() +
  theme_classic()
library(fpc)
eben<-cluster.stats(NULL,clustering=as.numeric(dataset1$class), alt.clustering = as.numeric(dataset1$Cluster))
library(igraph)
igraph::compare(as.numeric(dataset1$class), as.numeric(dataset1$Cluster), method = c("vi")) 
igraph::compare(as.numeric(dataset1$class), as.numeric(dataset1$Cluster), method = c("adjusted.rand")) 
library(mcclust)
vi.dist(as.numeric(dataset1$class), as.numeric(dataset1$Cluster))
library(clValid)
dunn(Data=dbdat, cluster=rpip$cluster)
table(dataset1$class, dataset1$Cluster)

newrun=FALSE



#################Tests and work on projection and study area definition#########

longlims=c(-6, 2)
longdom=seq(longlims[1],longlims[2],by=.25)
latlims=c(48,59)
latdom=seq(latlims[1],latlims[2],by=.25)
lagrid<-expand.grid(longdom,latdom)
logl<-c(-6,-5.75,1.75,2)
lagl<-c(48,48.25,58.75,59)
garea<-lagrid[which(!is.na(match(lagrid$Var1,logl)) | !is.na(match(lagrid$Var2,lagl))),]
gareax=garea
length(garea[,1])/length(lagrid[,1])

#test correctif de l'effet longlat
library(geosphere)
distGeo(c(longlims[2], latlims[1]), c(longlims[2], latlims[2]))/(1000*45)
distGeo(c(longlims[1], latlims[1]), c(longlims[2], latlims[1]))/(1000*33)

  la<-longdom[1]
  lag<-lagrid[which(lagrid$Var1==la),]
  log<-lagrid[which(lagrid$Var1==la+0.25),]
  laglog<-cbind(lag,log)
  laglog$Var1
  abc=c()
  fuck<-distGeo(c(laglog[1,1], laglog[1,2]), c(laglog[1,3], laglog[1,4]))/(1000)
  for(lo in 1:45){
    abc<-c(abc,(distGeo(c(laglog[lo,1], laglog[lo,2]), c(laglog[lo,3], laglog[lo,4]))/(1000)/fuck))
  }
plot(latdom,abc)


library(maps)

# Bonne equal-area projection with state abbreviations


library(mapproj)

optn<-mapproject(lagrid$Var1,lagrid$Var2,projection ="mercator")
plot(optn$x,optn$y)



###############Basemap - Important for all plots###########

mapUK = SpacializedMap(database="world",regions = c("UK","France","Ireland"))
plot(mapUK)
muki<-mapUK@polygons
ukk<- gSimplify(mapUK, tol=0.001, topologyPreserve=TRUE)
uk_fort <- ggplot2::fortify(ukk)
ukx<-uk_fort

lagrid.coord<-lagrid
coordinates(lagrid.coord) <- c("Var1", "Var2")
coordinates(ukx) <- c("long", "lat")
coordinates(garea)<-c("Var1","Var2")
ukxy<-as.SpatialPolygons.PolygonsList(muki, proj4string=CRS("+proj=longlat +datum=WGS84"))
proj4string(lagrid.coord) <- CRS("+proj=longlat +datum=WGS84")
proj4string(garea) <- CRS("+proj=longlat +datum=WGS84")
proj4string(ukx) <- CRS("+proj=longlat +datum=WGS84")
plot(ukxy)
ukxutm<-spTransform(ukxy, CRS("+proj=utm +zone=30U, +datum=WGS84"))

# end of important map settings

lagrid.utm <- spTransform(lagrid.coord, CRS("+proj=utm +zone=30U, +datum=WGS84"))
garea.utm<-spTransform(garea, CRS("+proj=utm +zone=30U, +datum=WGS84"))
plot(ukxutm)
points(lagrid.utm)

mamamia<-fortify(ukxutm)


lagfor<-data.frame(lagrid.utm@coords)
garea.utmf<-data.frame(garea.utm)
garea.longlat=spTransform(garea.utm, CRS("+proj=longlat +datum=WGS84"))
garea.llf<-data.frame(garea.longlat)
aziz<-data.frame(optn$x,optn$y)
33/(optn$range[2]-optn$range[1])
45/(optn$range[4]-optn$range[3])
aziz=466*aziz
plot(aziz)
diff(aziz$optn.y)



mapUK = SpacializedMap(database="world",regions = c("UK","France","Spain","Portugal","Italy","Ireland","Belgium","Netherland"))
plot(mapUK)
loli<-c(1.6e+5,9e+5)
lali<-c(5300000,6700000)
ukk<- gSimplify(mapUK, tol=0.001, topologyPreserve=TRUE)
uk_fort <- ggplot2::fortify(ukk)
lagfor$g=1

ggplot(mamamia, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill = "beige", color = "gray10", size = 1,alpha=.5) +
 geom_point(data=lagfor,aes(x=Var1,y=Var2,group=Var1),colour="grey60",size=1,shape=3)+
  theme_bw(16)+
  coord_fixed(xlim = loli,  ylim = lali, ratio = 1)+
  geom_point(data=garea.utmf,aes(x=Var1,y=Var2,group=Var1),alpha=.8,colour="red") +
scale_y_continuous("Latitude")+
  scale_x_continuous("Longitude")

########################Adding elevation background#########################
fileele="C:/Users/Alois/OneDrive - King's College London/DATA/04_Spatiotemporal/in/elev_0.1deg.nc"
elev<-nc_open(fileele)
name.var <- names(elev$var)
tdims=elev$ndims
tsize=c(1,1)

lonb<-c(which(round(elev$dim$longitude$vals,1)==longlims[1]-1),which(round(elev$dim$longitude$vals,1)==longlims[2]+1))  
latb<-c(which(round(elev$dim$latitude$vals,1)==latlims[1]-1),which(round(elev$dim$latitude$vals,1)==latlims[2]+1))
countlon<-lonb[2]-lonb[1]+1
countlat<-latb[2]-latb[1]+1
start <- rep(1,tdims) # begin with start=(1,1,1,...,1)
start[1]<-lonb[1]
start[2]<-latb[1]
count <- tsize # begin w/count=(nx,ny,nz,...,nt), reads entire var
count[1]<-countlon
count[2]<-countlat
elevation<-list()
elevation$data   = ncvar_get(elev,name.var,start,count)   
elevation$lon    = ncvar_get(elev,"longitude",start[1],count[1]) 
elevation$lat    = ncvar_get(elev,"latitude",start[2],count[2]) 
elonlat <- as.matrix(expand.grid(elevation$lon,elevation$lat))
ele1<-as.vector(elevation$data) 
alelu<-data.frame(elonlat,ele1)
alelu$gr=1
garea.llf$gr=1

ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill = "lemonchiffon1", color = "gray10", size = 1,alpha=.5) +
  theme_bw(16)+
  coord_map(xlim = longlims,  ylim = latlims, proj="albers",lat0=(50),lat1=(55))+
  # geom_tile(data=alelu,aes(x=Var1,y=Var2,fill=ele1,group=gr),alpha=.8,interpolate = F) +
  scale_fill_gradientn(colours = terrain.colors(100),na.value = "aliceblue")+
  geom_polygon(fill = "transparent", color = "gray10", size = 1.2)+
  theme(axis.title=element_text(size=18),
        panel.background = element_rect(fill = "aliceblue", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        panel.grid = element_line(colour="grey60", size=.5),
        panel.grid.minor  = element_line(colour="black", size=10),
        legend.position = "none",
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  geom_tile(data=lagrid,aes(x=Var1,y=Var2,group=Var1),alpha=.01, col="tomato4",fill="transparent") +
  geom_tile(data=gareax,aes(x=Var1,y=Var2,group=Var1),alpha=.8,fill="red") +
  scale_y_continuous(
    breaks = c(48,50,52,54,56,58),labels= c("48°N","50°N","52°N","54°N","56°N","58°N"),limits = c(40,70),"Latitude")+
  scale_x_continuous(
    breaks =c(-6,-4,-2,0,2),labels= c("-6°E","-4°E","-2°E","0°E","2°E"),"Longitude") 

 
################################Import ncdf files#####################################

filer=c(paste0(getwd(),"/in/raindat_7986.nc"),paste0(getwd(),"/in/raindat_8797.nc"),paste0(getwd(),"/in/raindat_9808.nc"),paste0(getwd(),"/in/raindat_0919.nc"))
filew=c(paste0(getwd(),"/in/windat_7986.nc"),paste0(getwd(),"/in/windat_8797.nc"),paste0(getwd(),"/in/windat_9808.nc"),paste0(getwd(),"/in/windat_0919.nc"))

Startdate=as.POSIXct("1979-01-01 10:00:00")
Enddate=as.POSIXct("1986-12-31 23:00:00")
ncr = nc_open(filer)
ncw = nc_open(filew)

haz<-c(filer,filew)
hazmat<-c()
timeco<-c()
hx=1

##############loop to open the 4 files for each variable###############

for (s in 1:4){
  if (hx==1) nc = nc_open(filer[s])
  if (hx==2) nc = nc_open(filew[s])

  nc$var$i10fg
  nc$var$p0001
  nc$var$p0005
  nc$var$tp
  file_time    <- ncvar_get(nc,'time')
  file_time.info=nc$dim$time$units   # recherche la date initial de la variable time
  file_time.origin=unlist(strsplit(file_time.info, " "))[[3]] 
  file_time.unit=unlist(strsplit(file_time.info, " "))[[1]] 
  if(length(nc$var)>1){
    name.var <- names(nc$var)[1]}else{
      name.var <- names(nc$var)
    }
  cat("Import parametre ", name.var,"\n")
  
  
  name.lon="longitude"
  name.lat="latitude"
  nc$data   = ncvar_get(nc,name.var)   
  nc$lon    = ncvar_get(nc,name.lon)
  nc$lat    = ncvar_get(nc,name.lat)
  dlname <- ncatt_get(nc,name.var,"long_name")
  dunits <- ncatt_get(nc,name.var,"units")
  fillvalue <- ncatt_get(nc,name.var,"_FillValue")
  
  dim(nc$lon)
  dim(nc$lat)
  time <- ncvar_get(nc,"time")
  time
  tunits <- ncatt_get(nc,"time","units")
  nt <- dim(time)
  nt
  
  
  library(chron)
  library(lattice)
  library(RColorBrewer)
  library(lubridate)
  # nc$lat <- rev(nc$lat)
  timestamp <- as_datetime(c(time*60*60),origin="1900-01-01")
  timeb<-c(which(timestamp==Startdate),which(timestamp==Enddate))
  timeb<-c(1,length(timestamp))
  lonb<-c(which(nc$lon==longlims[1]),which(nc$lon==longlims[2]))
  latb<-c(which(nc$lat==latlims[1]),which(nc$lat==latlims[2]))
  countime<-timeb[2]-timeb[1]+1
  countlon<-lonb[2]-lonb[1]+1
  countlat<-latb[2]-latb[1]+1
 
    t<-nc$var$tp
    if(length(t)==0) t=nc$var$i10fg
    if(length(t)==0){t<-nc$var$p0001}

    tsize<-t$varsize
    tdims<-t$ndims
    nt1<-tsize[tdims]
    
 
  # Initialize start and count to read one timestep of the variable.
  start <- rep(1,tdims) # begin with start=(1,1,1,...,1)
  start[tdims] <- timeb[1] # change to start=(1,1,1,...,i) to read timestep i
  start[1]<-lonb[1]
  start[2]<-latb[1]
  count <- tsize # begin w/count=(nx,ny,nz,...,nt), reads entire var
  count[1]<-countlon
  count[2]<-countlat
  count[tdims] <- countime # change to count=(nx,ny,nz,...,1) to read 1 tstep
  # 
  if(count[2]<0)
  {start[2]<-start[2]+count[2]-1
  count[2]=-count[2]+2}
  
  newformat<-list()
  
  

  newformat$data   = ncvar_get(nc,name.var, start = start, count= count)   
  newformat$lon    = ncvar_get(nc,name.lon,start = start[1], count= count[1]) 
  newformat$lat    = nc$lat[seq(start[2],start[2]+count[2]-1)]
  newformat$lat
  newformat$lon
  
  if(dlname$value=="Total precipitation") newformat$data<-newformat$data*1000
  
  dlname <- ncatt_get(nc,name.var,"long_name")
  dunits <- ncatt_get(nc,name.var,"units")
  fillvalue <- ncatt_get(nc,name.var,"_FillValue")
  time <- ncvar_get(nc,"time")
  time<-time[seq(start[3],start[3]+count[3]-1)]
  tunits <- ncatt_get(nc,"time","units")
  nt <- dim(time)
  which(timestamp== "2009-11-18 UTC")
  temp11 <- newformat$data[ , , 7705] 
  for (merd in 1:36){
    temp11=temp11+newformat$data[ , , 7705+merd] #Level is the third dimension and time the fourth.
  }

  lon=newformat$lon
  lat=newformat$lat
  max(temp11)
  rgb.palette=colorRampPalette(c("white", 
                                 "royalblue","orange","red","purple"),interpolate="linear",bias=1)
  grid <- expand.grid(lon=lon, lat=lat)
  cutpts <- c(1,5,10,15,20,25,30,35,40)
  levelplot(temp11 ~ lon * lat, data=grid,cuts=50, pretty=T, 
            col.regions=(rgb.palette(200)))
  
  hazmat<-c(hazmat,list(newformat))
  if (hx==2){
  if(s==1)save(newformat,file="interdat/windP1.Rdata")
  if(s==2)save(newformat,file="interdat/windP2.Rdata")
  if(s==3)save(newformat,file="interdat/windP3.Rdata")
  if(s==4)save(newformat,file="interdat/windP4.Rdata")
  }
  if(hx==1){
  if(s==1)save(newformat,file="interdat/rainP1.Rdata")
  if(s==2)save(newformat,file="interdat/rainP2.Rdata")
  if(s==3)save(newformat,file="interdat/rainP3.Rdata")
  if(s==4)save(newformat,file="interdat/rainP4.Rdata")
  }
  timeco<-c(timeco,list(time))
  if(s==1)save(time,file="interdat/timeP1.Rdata")
  if(s==2)save(time,file="interdat/timeP2.Rdata")
  if(s==3)save(time,file="interdat/timeP3.Rdata")
  if(s==4)save(time,file="interdat/timeP4.Rdata")
}


rbPal <- rbPal <- colorRampPalette(c('royalblue',"skyblue","gold","darkorange",'red',"purple"))
thx=T

#If required

# load(file="interdat/rainP4.Rdata")
# load(file="interdat/timeP4.Rdata")

timix<-as_datetime(c(time*60*60),origin="1900-01-01")
l1<-which(newformat$lon==-4.75)
l2<-which(newformat$lat==50.25)
ole<-which(month(timix)==11 & year(timix)==2010)
mierda<-newformat$data[l1,l2,which(newformat$data[l1,l2,]>=-1)]
plot(timix[ole],mierda[ole],type="o")
macs<-mierda[which(mierda>0)]
qtest<-quantile(macs,.99,na.rm=T)
abline(h=qtest)


mac<-c()
idd<-c()
x=0
for (id in c(x+1:(length(mierda)-x))){
mac[id]<-sum(mierda[c((id-x):(id+x))]) 
idd[id]<-id
}

which(mac==max(mac,na.rm=T))
plot(timix[ole],mac[ole],type="o")
oula<-mierda[ole]
tit<-timix[ole]
macs<-mac[which(mac>0)]
qtest<-quantile(macs,.99,na.rm=T)
abline(h=qtest)

aa<-which(mac[ole]>qtest)
length(aa)
diff(aa)
length(which(diff(aa)>1))
ot<-oula[aa]
tito<-tit[aa]
points(tito,ot,col=2)

mac=mierda
qtest<-quantile(mac,.99,na.rm=T)
abline(h=qtest)
mi<-sum(mac,na.rm=T)
mierdax<-idd[which(mac>qtest)]
id3=mierdax
id24=mierdax
id0=mierdax
id6=mierdax
plot(mac)
length(which(!is.na(match(id24,id3))))
mierdax/mierda
plot(sort(newformat$data[29,1,which(newformat$data[29,1,]>0)]))
if (thx==T){
pqt<-.99

#Set up a quantile for each grid cell

thr<-hazmat[[1]]$data[,,1]
for (i in (1:33)){
  for(j in 1:45){
    crappy<-c(hazmat[[1]]$data[i,j,],hazmat[[2]]$data[i,j,],hazmat[[3]]$data[i,j,],hazmat[[4]]$data[i,j,])
    thr[i,j]<-quantile(crappy[which(crappy>0)],pqt,na.rm=T)
  }
}

}


# coord_map(xlim = longlims,  ylim = latlims, proj="albers",lat0=(50),lat1=(55))


# save(thr,file="Wnd_99_AllP.Rdata")

#Now I have the threshold vectors.

#===============================================================#

#Plot for figure 4

load(paste0(getwd(),"/data/interdat/Wnd_99_AllP.Rdata"))
thw<-thr
load("data/interdat/Rain_99_AllP.Rdata")
rbPal <- rbPal <- colorRampPalette(c('royalblue',"skyblue","gold","darkorange",'red',"purple"))
thww<-as.vector(thw) 
thrr<-as.vector(thr) 
min(thrr)

#load one file for long lat
load(file=paste0(getwd(),"/data/interdat/windP1.Rdata"))
elonlat <- as.matrix(expand.grid(newformat$lon,newformat$lat))
thwbg<-data.frame(elonlat,thww)
thrbg<-data.frame(elonlat,thrr)
thwbg$gr=1
thrbg$gr=1

longlims=c(-5.7,1.7)
latlims=c(48.4,58.5)

ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  theme_bw(16)+
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3)+
  geom_raster(data=thrbg,aes(x=Var1,y=Var2,fill=thww,group=gr),alpha=1,interpolate = F) +
  #scale_fill_gradientn(name = expression(paste("w threshold [m s"^"-1","]")),colours = rbPal(100),na.value = "aliceblue")+
  scale_fill_distiller(palette = "YlOrBr",name = expression(paste("w threshold [m s"^"-1","]")),breaks=breaks_extended(6),guide="coloursteps",direction=1)+
  geom_polygon(fill = "transparent", color = "gray10", size = 1.2)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_y_continuous(
    breaks = c(48,50,52,54,56,58),labels= c("48°N","50°N","52°N","54°N","56°N","58°N"),limits = c(40,70),"Latitude")+
  scale_x_continuous(
    breaks =c(-6,-4,-2,0,2),labels= c("-6°E","-4°E","-2°E","0°E","2°E"),limits=c(-10,10),"Longitude") 


ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  theme_bw(16)+
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3)+
  geom_raster(data=thrbg,aes(x=Var1,y=Var2,fill=thrr,group=gr),alpha=1,interpolate = F) +
  #scale_fill_gradientn(name = expression(paste("p threshold [mm]")),colours = rbPal(100),na.value = "aliceblue")+
  scale_fill_distiller(palette = "YlGnBu",name = expression(paste("p threshold [mm]")),breaks=breaks_extended(6),guide="coloursteps",direction=1)+
  geom_polygon(fill = "transparent", color = "gray10", size = 1.2)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_y_continuous(
    breaks = c(48,50,52,54,56,58),labels= c("48°N","50°N","52°N","54°N","56°N","58°N"),limits = c(40,70),"Latitude")+
  scale_x_continuous(
    breaks =c(-6,-4,-2,0,2),labels= c("-6°E","-4°E","-2°E","0°E","2°E"),limits=c(-10,10),"Longitude") 

############Preparation of the data for clustering###########

rm(nc,newformat,hazmat)
gc()
sptdfx<-list()
metaHaz<-list()
metavHour<-list()
metavDaz<-list()
hx=1
wala=F

#################loop for clustering##############

for(hazard in 1:4){
  if(hx==2){
    if(hazard==1){load("data/interdat/windP1.Rdata");load(file="interdat/timeP1.Rdata")}
    if(hazard==2){load("interdat/windP2.Rdata");load(file="interdat/timeP2.Rdata")}
    if(hazard==3){load("interdat/windP3.Rdata");load(file="interdat/timeP3.Rdata")}
    if(hazard==4){load("interdat/windP4.Rdata");load(file="interdat/timeP4.Rdata")}
    th=thw
  }
  if(hx==1){
    if(hazard==1){load("data/interdat/rainP1.Rdata");load(file="data/interdat/timeP1.Rdata")}
    if(hazard==2){load("data/interdat/rainP2.Rdata");load(file="data/interdat/timeP2.Rdata")}
    if(hazard==3){load("data/interdat/rainP3.Rdata");load(file="data/interdat/timeP3.Rdata")}
    if(hazard==4){load("data/interdat/rainP4.Rdata");load(file="data/interdat/timeP4.Rdata")}
    th=thr
  }

  
  if (hx == 3){
    if(hazard==1){load("interdat/rainP1.Rdata");load(file="interdat/timeP1.Rdata");hazdat1<-newformat$data;load("interdat/windP1.Rdata")}
    if(hazard==2){load("interdat/rainP2.Rdata");load(file="interdat/timeP2.Rdata")}
    if(hazard==3){load("interdat/rainP3.Rdata");load(file="interdat/timeP3.Rdata")}
    if(hazard==4){load("interdat/rainP4.Rdata");load(file="interdat/timeP4.Rdata")}
  }
  wala=F
  hazdat=newformat
  rm(newformat)
  gc()
  lon=hazdat$lon
  lat=hazdat$lat
  oups<-length(time)-length(hazdat$data[1,1,])
  if (length(oups>0)){
    rm=seq(1:oups)+length(hazdat$data[1,1,])
    time<-time[-rm]
  }

  
  lonlatime <- expand.grid(lon, lat,time)
  
 if(hx==3){
   hazdat2<-newformat$data
   vecthouf1<-as.vector(rep(thr,length(hazdat2[1,1,])))
   vecthouf2<-as.vector(rep(thw,length(hazdat2[1,1,])))
   bolilos<-hazdat1
   
   bolilos[which(hazdat1<vecthouf1 | hazdat2<vecthouf2)] <- NA
   bolilos[which(hazdat1>=vecthouf1 & hazdat2>=vecthouf2)] <- 1  
   hazdat1[which(hazdat1<vecthouf1 | hazdat2<vecthouf2)] <- NA
   hazdat2[which(is.na(hazdat1))] <- NA
   
   formeta<-na.omit(as.vector(hazdat1))
   formeta2<-na.omit(as.vector(hazdat2))
   vecmeta<-data.frame(formeta,formeta2)
   
   vectouf<- as.vector(bolilos)
   length(na.omit(vectouf))/length(vectouf)
   
 }  else{
    vectouf<- as.vector(hazdat$data)
    vecthouf<-as.vector(rep(th,length(hazdat$data[1,1,])))
    vectouf[which(vectouf<vecthouf)]<-NA
    vectouf[which(vectouf>=vecthouf)]<-1
    vecmeta<-as.vector(hazdat$data)
    # vecmeta2<-as.vector(formeta2)
    if(hx==1){vecmeta[which(vecmeta<vecthouf)] <- NA}
    if(hx==2){vecmeta[which(vecmeta<vecthouf)]<- NA}
  
  metav<-data.frame(cbind(lonlatime,vecmeta))
  
  metav<-metav[which(!is.na(metav[,4])),]

}
  lonlatemp <- data.frame(cbind(lonlatime,vectouf))
  lonlatemp2<-lonlatemp[which(lonlatemp[,4]==1),]
  
  rm(lonlatemp,vecthouf,vectouf)
  gc()
  spdata<-lonlatemp2[,-4]
  
  if(hx==3){metav<-data.frame(cbind(spdata,vecmeta))}
  print(length(spdata$Var1))
  ep<-2
  coef<-4
  sampspd<-spdata
  sampspd$Var1<-sampspd$Var1*coef
  sampspd$Var2<-sampspd$Var2*coef
  sampspd$Var3<-sampspd$Var3-sampspd$Var3[1]+1
  sampspd$Var3<-sampspd$Var3
  samptt<-as.matrix(sampspd)
  

  
  ##############kNN plot for selection of eps################

  if (wala==T){
    walabibou<-kNNdist(sampspd,k=10,all=F)
    min(walabibou)
    walabibou<-walabibou[which(walabibou<10)]
    walabibof=jitter(as.vector(walabibou),15)
    walabibord<-seq(1:length(walabibou))
    walabibix<-data.frame(as.vector(walabibof[order(walabibof)]))
    walabibix<-data.frame(as.vector(walabibou[order(walabibou)]))
    
    walabibix$order=walabibord
    
    plot(walabibix)
    fatos<-walabibix[which(walabibix$order>1270000),]
    light<-seq(1,1270000,by=1000)
    #
    walabilight<-rbind(walabibix[which(!is.na(match(walabibix$order,light))),],fatos)
    lx<-length(walabilight[,1])
    lacurve=approx(walabilight[c(1,lx),2],walabilight[c(1,lx),1],n=200)
    lacurve<-as.data.frame(lacurve)
    xyf<-as.data.frame(approx(walabilight[,2],walabilight[,1],n=200))
    y_values=xyf[,2]
    x_values=xyf[,1]
    xcurve=lacurve[,1]
    ycurve=lacurve[,2]
    
    elbow_finder <- function(x_values, y_values,xcurve,ycurve) {
      # Max values to create line
      max_x_x <- min(x_values)
      max_x_y <- y_values[which.min(x_values)]
      max_y_y <- max(y_values)
      max_y_x <- x_values[which.max(y_values)]
      max_df <- data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))
      
      # Creating straight line between the max values
      fit <- lm(max_df$y ~ max_df$x)
      
      # Distance from point to line
      x_min_dist<-c()
      y_min_dist<-c()
      distances <- c()
      for(i in 1:length(x_values)) {
        distancex<-c()
        for(j in 1:length(xcurve)){
          # distancex <- c(distancex, abs(coef(fit)[2]*x_values[j] - y_values[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2))
          distancex<-c(distancex,sqrt((x_values[i]-xcurve[j])^2+(y_values[i]-ycurve[j])^2))
          # plot(distancex)
        }
        distances<-c(distances,min(distancex, na.rm=T))
        x_min_dist <- c(x_min_dist,xcurve[which.min(distancex)])
        y_min_dist <- c(y_min_dist,ycurve[which.min(distancex)])
      }
      
      # Max distance point
      x_max_dist <- x_values[which.max(distances)]
      y_max_dist <- y_values[which.max(distances)]
      xcmax <- x_min_dist[which.max(distances)]
      ycmax <- y_min_dist[which.max(distances)]
      
      
      poncu<-c()
      plot(x_values,distances)
      
      return(c(x_max_dist, y_max_dist,xcmax,ycmax))
    }
    
    
    merd2<-elbow_finder(x_values,y_values,xcurve,ycurve)

    
    plot(lacurve)
    points(x_values,y_values,col=3)
    points(merd2[1],merd2[2],col=2,pch=16)
    points(merd2[3],merd2[4],col=2,pch=16)
    lines(c(merd2[1],merd2[3]),c(merd2[2],merd2[4]))
    
    xtick=seq(0,1.5e6,by=2e5)
    par(mar=c(5,5,1,1))
    plot(walabilight[,2],walabilight[,1],ylim=c(0,10),ylab="10-NN distance",xlab="Points sorted by distance",cex.axis=1.5,cex.lab=1.8,xaxt="n")
    axis(side=1, at=xtick, labels = T,cex.axis=1.5)
    eps=round(merd2[2],2)

    txt<-bquote(epsilon == .(round(eps,2)))
    text(2e3,2.5,cex=1.8, labels=txt,pos=3,col="red")
    abline(h=merd2[2], col=2,lwd=3)
    
    
    # walabibou[order(walabibou)]
    #
  }
  
  if (hx==1)eps=2.45
  if (hx==2)eps=2.24
  if (hx==3)eps=2.45

  epcl=eps

  
  ##Update on epcl for rain events
 cluster=T

  rpip<-dbscan::dbscan(sampspd, eps=epcl, minPts = 10)
  
  length(unique(rpip$cluster))

    
  spdata<-cbind(spdata,rpip$cluster)
  sptdf<-cbind(metav,rpip$cluster)
  if (cluster==T){
  
  if(length(which(spdata$`rpip$cluster`==0))>0){
    metav<-metav[-which(spdata[,4]==0),]
    spdata<-spdata[-which(spdata[,4]==0),]
  }
  metav$cluster<-spdata[,4]
  

  spdata$Var3=spdata$Var3-spdata$Var3[1]+1
  length(metav$Var1)

  library(lubridate)
  metav$time<-as_datetime(c(metav$Var3*60*60),origin="1900-01-01")
  metav$month=month(metav$time)
  event<-metav
  charloc<-paste(event[,1],event[,2])
  event$cloc=charloc
  print(length(event$Var1))
  
  if (hx==3){
    testev<-aggregate(list(rf= event[,5],wg=event[,6]),
                      by = list(ev = event[,4],loc=event[,9]),
                      FUN = function(x) c(sum = sum(na.omit(x)),max=max(na.omit(x))))
    testev <- do.call(data.frame, testev)
    
    metamax<-aggregate(list(vir= testev[,3],viw=testev[,6]) ,
                       by = list(ev = testev[,1]),
                       FUN = function(x) c(sum = sum(na.omit(x)),max=max(na.omit(x)),mean=mean(x),sd=sd(x),surf=length(x)))
    
    metamax<- do.call(data.frame, metamax)
  }
  
  testev<-aggregate(list(vi= event[,4]),
                    by = list(ev = event[,5],loc=event[,8]),
                    FUN = function(x) c(sum = sum(na.omit(x)),max=max(na.omit(x))))
  testev <- do.call(data.frame, testev)
  
  thl<-aggregate(list(vi= event[,4]),
                 by = list(ev = event[,5]),
                 FUN = function(x) c(l= length(x)))
  thl<- do.call(data.frame, thl)
  small<-thl$ev[which(thl$vi<10)]
  if(length(small)>0){
  bip<-which(!is.na(match(testev$ev,small)))
  testev<-testev[-bip,]
  bop<-which(!is.na(match(event$cluster,small)))
  event<-event[-bop,]
  }
  length(unique(event$cluster))
  length(unique(testev$ev))
  
  metamax<-aggregate(list(vir= testev[,3],viw=testev[,4]) ,
                     by = list(ev = testev[,1]),
                     FUN = function(x) c(sum = sum(na.omit(x)),max=max(na.omit(x)),mean=mean(x),sd=sd(x),surf=length(x)))
  
  metamax<- do.call(data.frame, metamax)
  
  metave<-aggregate(list(vi= event[,4]) ,
                    by = list(ev = event[,5]),
                    FUN = function(x) c(surf=length(x)))
  
  metave <- do.call(data.frame, metave)
  length(unique(metave$ev))

  tempcom<-aggregate(event[,6] ,
                     by = list(ev = event[,5]),
                     FUN = function(x) c(dur=length(unique(x)),month=month(unique(x))[1],year=year(unique(x))[1]))
  
  tempcom<- do.call(data.frame, tempcom)
  if (hx==2) metave<-cbind(metamax[,c(1,7:11)],metave[,c(2)],tempcom[,c(2,3,4)])
  if (hx==1) metave<-cbind(metamax[,c(1:6)],metave[,c(2)],tempcom[,c(2,3,4)])
  maxR<-c()
  maxW<-c()
  evbk<-c()
  startev<-vector(length=length(metave$ev))
  endev<-vector(length=length(metave$ev))
  #Need to check all this
  for (eve in 1:length(metave$ev))
  {
    eves<-metave$ev[eve]
    eventr<-event$time[which(event$cluster==eves)]
    evint<-testev[which(testev$ev==eves),]
    met<-metamax[which(metamax$ev==eves),]
    startev[eve]<-min(eventr)
    endev[eve]<-max(eventr)
    elW<-as.character(evint$loc[which(evint$vi.max==met$viw.max)])[1]
    elR<-as.character(evint$loc[which(evint$vi.sum==met$vir.max)])[1]
    maxR<-c(maxR,elR)
    maxW<-c(maxW,elW)
    evbk<-c(evbk,eves)

  }
  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

  metave$season=metave$x.month
  
  metave$season[which(metave$x.month<3 | metave$x.month>11)]=1
  metave$season[which(metave$x.month<6 & metave$x.month>2)]=2
  metave$season[which(metave$x.month<10 & metave$x.month>5)]=3
  metave$season[which(metave$x.month<12& metave$x.month>9)]=4
  

  metave$startev<-startev
  metave$endev<-endev
  metaveN<-metave

  length(unique(metaveN$ev))
  metaHaz<-c(metaHaz,list(metaveN))
  metavHour<-c(metavHour,list(event))
  metavDaz<-c(metavDaz,list(testev))
}
  
  sptdfx<-c(sptdfx,list(sptdf))

}


sptdft<-rbind(sptdfx[[1]],sptdfx[[2]],sptdfx[[3]],sptdfx[[4]])
if(hx==1){
  save(sptdft,file="out/extremEvents_Rain.Rdata")}
if(hx==2){
  save(sptdft,file="out/extremEvents_Wind.Rdata")}
  



library(sp)
logl<-c(-6,-5.75,1.75,2)
lagl<-c(48,48.25,58.75,59)

coords = matrix(c(-5.75, 48.25,
                  -5.75, 58.75,
                  1.75, 58.75,
                  1.75, 48.25,
                  -5.75, 48.25), 
                ncol = 2, byrow = TRUE)


P1 = Polygon(coords)
Ps1 = SpatialPolygons(list(Polygons(list(P1), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
plot(Ps1, axes = TRUE)

pif <- ggplot2::fortify(Ps1)

############creating files containing metadata about clusters############
metaHax<-list()
metavHax<-list()
metavDax<-list()
for (hx in 1:4){
  idout<-c()
  for (cl in metaHaz[[hx]]$ev){ 
    clev<-metavHour[[hx]][which(metavHour[[hx]]$cluster==cl),]
    overlappy<-point.in.polygon(clev[,c(1)],clev[,c(2)],pif$long,pif$lat)
    if(length(which(overlappy>0))==0)idout<-c(idout,cl)
  }
  idout
  metax<-na.omit(match(idout,metaHaz[[hx]]$ev))
  metat<-as.numeric(which(!is.na((match(metavHour[[hx]]$cluster,idout)))))
  metas<-as.numeric(which(!is.na((match(metavDaz[[hx]]$ev,idout)))))
  metar<-metaHaz[[hx]][-metax,]
  metavr<-metavHour[[hx]][-metat,]
  metasr<-metavDaz[[hx]][-metas,]
  metaHax<-c(metaHax,list(metar))
  metavHax<-c(metavHax,list(metavr))
  metavDax<-c(metavDax,list(metasr))
  
}

rm(metaHaz,metavHour,metavDaz,lonlatemp,vecmeta,vecmeta2,vecthouf,vecthouf1,vecthouf2,vectouf,vrac,hazdat,hazdat2,newformat,bolilos,formeta,formeta2)
gc()

for (hx in 1:4){
  metaHax[[hx]]$startev<-as_datetime(c(metaHax[[hx]]$startev),origin="1970-01-01")
  metaHax[[hx]]$endev<-as_datetime(c(metaHax[[hx]]$endev),origin="1970-01-01")
}

for (hx in 2:4){
  metaHax[[hx]]$ev=metaHax[[hx]]$ev+max(metaHax[[hx-1]]$ev)
  metavHax[[hx]]$cluster=metavHax[[hx]]$cluster+max(metavHax[[hx-1]]$cluster,na.rm=T)
  metavDax[[hx]]$ev=metavDax[[hx]]$ev+max(metavDax[[hx-1]]$ev)
}
for (hx in 1:3){
difevent<-difftime(min(metaHax[[hx+1]]$startev),max(metaHax[[hx]]$endev),units = "hours")
print(difevent)
if(difevent<2){
event2<-metaHax[[hx+1]][which(metaHax[[hx+1]]$startev==min(metaHax[[hx+1]]$startev)  ),]
event1<-metaHax[[hx]][which(metaHax[[hx]]$endev==max(metaHax[[hx]]$endev)  ),]
cou1<-metavHax[[hx]][which(metavHax[[hx]]$cluster==event1$ev),]
cou2<-metavHax[[hx+1]][which(metavHax[[hx+1]]$cluster==event2$ev),] 
if(length(match(cou1$cloc,cou2$cloc))>0){
  metaHax[[hx+1]]$ev[which(metaHax[[hx+1]]$startev==min(metaHax[[hx+1]]$startev))]=event1$ev
  metavHax[[hx+1]]$cluster[which(metavHax[[hx+1]]$cluster==event2$ev)]=event1$ev
  metavDax[[hx+1]]$ev[which(metavDax[[hx+1]]$ev==event2$ev)]=event1$ev
}

Rcart=strsplit(cou1$cloc," ")
Rcar<-as.data.frame(matrix(unlist(Rcart),ncol=2,byrow=T))
Rcar[,1]<-as.numeric.factor((Rcar[,1]))
Rcar[,2]<-as.numeric.factor((Rcar[,2]))

Rcart=strsplit(cou2$cloc," ")
Rcar2<-as.data.frame(matrix(unlist(Rcart),ncol=2,byrow=T))
Rcar2[,1]<-as.numeric.factor((Rcar2[,1]))
Rcar2[,2]<-as.numeric.factor((Rcar2[,2]))
plot(Rcar,xlim=c(-6,3),ylim=c(47,55))
points(Rcar2,col=2,pch=3)

}
}

#########ha is the hazard to be saved 1 for rain and 2 for wind###########
ha=2

if(ha==1){
  save(metaHax,file="interdat/metaclusRain.Rdata")
  save(metavHax,file="interdat/rawclustRain.Rdata")
  save(metavDax,file="interdat/interclustRain.Rdata")
  
}
if(ha==2){
  save(metaHax,file="interdat/metaclustWind.Rdata")
  save(metavHax,file="interdat/rawclustWind.Rdata")
  save(metavDax,file="interdat/interclustWind.Rdata")
}


load(file="interdat/metaclustWind.Rdata")
load(file="interdat/rawclustWind.Rdata")
load(file="interdat/interclustWind.Rdata")
metaHaf<-rbind(metaHax[[1]],metaHax[[2]],metaHax[[3]],metaHax[[4]])
metavHaf<-rbind(metavHax[[1]],metavHax[[2]],metavHax[[3]],metavHax[[4]])
metavDaf<-rbind(metavDax[[1]],metavDax[[2]],metavDax[[3]],metavDax[[4]])
if(ha==2){
save(metaHaf,file="out/WindEv_meta_1979-2019.Rdata")
save(metavHaf,file="out/WindEv_hdat_1979-2019.Rdata")
save(metavDaf,file="out/WindEv_ldat_1979-2019.Rdata")
}
if(ha==1){
save(metaHaf,file="out/RainEv_meta_1979-2019.Rdata")
save(metavHaf,file="out/RainEv_hdat_1979-2019.Rdata")
save(metavDaf,file="out/RainEv_ldat_1979-2019.Rdata")
}



rm(lonlatime,ole,nvo,newformat,arain)
gc()
########### here ha is the file id (1,2,3,4) and hx is hazard id (rain,wind)#########
for(ha in 1:4){
if(hx==1){
  if(ha==1){load("interdat/rainP1.Rdata");time=timeco[[1]]}
  if(ha==2){load("interdat/rainP2.Rdata");time=timeco[[2]]}
  if(ha==3){load("interdat/rainP3.Rdata");time=timeco[[3]]}
  if(ha==4){load("interdat/rainP4.Rdata");time=timeco[[4]]}
}
arain<-as.vector(newformat$data)
max(arain)

lon=newformat$lon
lat=newformat$lat

lonlatime <- expand.grid(lon, lat,time)
lonlatime$time<-as_datetime(c(lonlatime$Var3*60*60),origin="1900-01-01")
head(lonlatime$time)
rm(newformat)
gc()

oc<-c()
nvo<-c()
nwo<-c()
nvo<-vector(mode="list",length=length(metaHax[[1]]$ev))

ptm <- proc.time()
for (c in 1: length(metaHax[[ha]]$ev)){ 
  
  cl=metaHax[[ha]]$ev[c]
  clev<-metavDax[[ha]][which(metavDax[[ha]]$ev==cl),]
  # dd<-unique(metavHax[[1]]$time[which(metavHax[[1]]$cluster==cl)])
  dd<-seq(metaHax[[ha]]$startev[c],metaHax[[ha]]$endev[c],by="hour")
  
  clex<-do.call(rbind, replicate(length(dd),clev[,c(1,2)], simplify=FALSE))
  cley<-do.call(rbind, replicate(length(clev[,1]),as.data.frame(dd),simplify=F))
  clex<-clex[order(clex$loc),]
  ouch<-cbind(clex,cley)
  
  nvo[[c]]<-ouch 
  print(round((c/length(metaHax[[ha]]$ev)*100),1))
  # nvo<-rbind(nvo,ouch)
  
}
proc.time()-ptm
nvo2<-ldply(nvo, data.frame)

newrun=T
if(newrun==T){
  lonlatime$loc<-paste(lonlatime$Var1,lonlatime$Var2)
  lonlatime<-lonlatime[,-c(1,2)]
  
  
  yoboys<-which(!is.na(fmatch(lonlatime$time,nvo2$dd)))
  ole<-data.frame(lonlatime[yoboys,],arain[yoboys])
  
  
  names(ole)[3] = "loc"
  names(ole)[2] = "dd"
  olel<-inner_join(ole,nvo2, by=c("loc","dd"))
  if(ha==1)save(olel,file="interdat/allraininclusters1.Rdata")
  if(ha==2)save(olel,file="interdat/allraininclusters2.Rdata")
  if(ha==3)save(olel,file="interdat/allraininclusters3.Rdata")
  if(ha==4)save(olel,file="interdat/allraininclusters4.Rdata")
}else{
  # load(file="allraininclusters.Rdata")
  print("chepa")
}
}

load(file="interdat/metaclustRain.Rdata")
load(file="out/RainEv_hdat_1979-2019.Rdata")

for(ha in 1:4){
  if(ha==1)load(file="interdat/allraininclusters1.Rdata")
  if(ha==2)load(file="interdat/allraininclusters2.Rdata")
  if(ha==3)load(file="interdat/allraininclusters3.Rdata")
  if(ha==4)load(file="interdat/allraininclusters4.Rdata")
rainacc<-aggregate(list(len=olel[,4]),
                   by = list(ev=olel[,5],loc=olel[,3]),
                   FUN = function(x) c(length=length(x),sum=sum(x),max=max(x)))
rainacc <- do.call(data.frame, rainacc)

rainacev<-aggregate(list(len=rainacc[,4]),
                    by = list(ev=rainacc[,1]),
                    FUN = function(x) c(length=length(x),sum=sum(x),max=max(x)))
rainacev <- do.call(data.frame, rainacev)
rmaxloc=c()
for (idev in rainacev$ev){
  racc<-rainacc[which(rainacc$ev==idev),]
  racev<-rainacev[which(rainacev$ev==idev),]
  locmax=which(racc$len.sum==racev$len.max)
  if(length(locmax)>1)locmax=locmax[1]
  locas=racc$loc[locmax]
  rmaxloc=c(rmaxloc,locas)
}

#last join
metaHax[[ha]]<-data.frame(metaHax[[ha]],rainacev[,c(3,4)],rmaxloc)
plot(metaHax[[ha]]$vir.max,metaHax[[ha]]$len.max) 
abline(a=0,b=1)
}

metaHar<-rbind(metaHax[[1]],metaHax[[2]],metaHax[[3]],metaHax[[4]])

save(metaHar,file="out/RainEv_metap2_1979-2019.Rdata")
evp<-metaHaf$ev[which.is.max(metaHaf$len.max)]
evento<-metavHaf[which(metavHaf$cluster==3250),]

#########################Quick event visualizer##################
library(ggalt)
ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill = "transparent", color = "gray10", size = 1) +
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
  geom_encircle(aes(x=Var1, y=Var2,group=cluster),
                 data = evento, size = 2, color = "blue",s_shape=.8, expand=0,fill="blue",alpha=.1)+
  geom_bin2d(data=evento,aes(x=Var1,y=Var2,group=cluster,weight=vecmeta),alpha=.7,binwidth=c(.25,.25))+
  # stat_density_2d(data=evento,geom = "polygon", aes(x=Var1,y=Var2,group=cluster,fill=..level..),alpha=.3)+
  scale_fill_viridis_c()



  geom_encircle(aes(x=Var1, y=Var2,group=gr),
                data = maxeventwx, size = 2, color = "orange",s_shape=.8, expand=0,fill="orange",alpha=.5)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="italic"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_y_continuous(
    breaks = c(48,50,52,54,56,58,60),limits = c(40,70),"Latitude")+
  scale_x_continuous(
    breaks =c(-6,-4,-2,0,2),limits=c(-10,10),"Longitiude") 








hx=2
if(hx==2){
  load(file="interdat/metaclustWind.Rdata")
  load(file="interdat/rawclustWind.Rdata")
  load(file="interdat/interclustWind.Rdata")
  nwo3<-c()
}

if(hx==1){
  load(file="interdat/metaclustRain.Rdata")
  load(file="interdat/rawclustRain.Rdata")
  load(file="interdat/interclustRain.Rdata")
  nvo3<-c()
}
# create progress bar
pb <- txtProgressBar(min = 0, max = length(metaHax[[ha]]$ev), style = 3)

##########Creation of files for compound hazard identification###############
for (ha in 1:4){
nwo<-vector(mode="list",length=length(metaHax[[ha]]$ev))
for (c in 1: length(metaHax[[ha]]$ev)){ 
  cl=metaHax[[ha]]$ev[c]
  clev<-metavDax[[ha]][which(metavDax[[ha]]$ev==cl),]
  # dd<-unique(metavHax[[2]]$time[which(metavHax[[2]]$cluster==cl)])
  dd<-seq(metaHax[[ha]]$startev[c],metaHax[[ha]]$endev[c],by="hour")
  
  clex<-do.call(rbind, replicate(length(dd),clev[,c(1,2)], simplify=FALSE))
  cley<-do.call(rbind, replicate(length(clev[,1]),as.data.frame(dd),simplify=F))
  clex<-clex[order(clex$loc),]
  ouch<-cbind(clex,cley)
  # nwo<-rbind(nwo,ouch)
  nwo[[c]]<-ouch
  Sys.sleep(0.1)
  # update progress bar
  setTxtProgressBar(pb, c)
  # print(round((c/length(metaHax[[ha]]$ev)*100),1))
}
close(pb)
nwo2<-ldply(nwo, data.frame)
rm(nwo)
gc()

if(hx==2)nwo3<-rbind(nwo3,nwo2)
if(hx==1)nvo3<-rbind(nvo3,nwo2)
}
mrdx<-unique(nwo3$ev)
plot(mrdx)
save(nwo3,file="out/Wind_stfprint2.Rdata")
save(nvo3,file="out/Rain_stfprint.Rdata")
