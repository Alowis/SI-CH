
rm(list=ls())  
gc()
getwd()
setwd("C:/Users/PhD Student/OneDrive - King's College London/DATA/04_Spatiotemporal")
setwd("C:/Users/k1638615/King's College London/OneDrive - King's College London/DATA/04_Spatiotemporal")
config_file=paste0(getwd(),"/config/config_general.R")
source(config_file)
library(ncdf4)
library(dbscan)

#===========================================================================================================

#I have to create my spatiotemporal cube


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

#Synthetic data example for dbscan with data of the same format as the input data


##############################################
#What i could do
# Do a sample for OR and find wind dominated, rain dominate and compound events
#1- for each grid cell: how many time within an event
#2- for each event: how long is it within the event (if not contimuous just add)
#3- for each event: what is its size
#4- Use some quantiles to display results in map



newrun=FALSE

longlims=c(-6, 2)
longdom=seq(longlims[1],longlims[2],by=.25)
latlims=c(48,59)
latdom=seq(latlims[1],latlims[2],by=.25)
lagrid<-expand.grid(longdom,latdom)
logl<-c(-6,-5.75,1.75,2)
lagl<-c(48,48.25,58.75,59)
garea<-lagrid[which(!is.na(match(lagrid$Var1,logl)) | !is.na(match(lagrid$Var2,lagl))),]
length(garea[,1])/length(lagrid[,1])

#conflict somewhere here
mapUK = SpacializedMap(database="world",regions = c("UK","France","Spain","Portugal","Italy","Ireland","Belgium","Netherland"))


ukk<- gSimplify(mapUK, tol=0.001, topologyPreserve=TRUE)
uk_fort <- ggplot2::fortify(ukk)
ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill = "white", color = "gray10", size = 1) +
  theme_bw(16)+
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3)

# file=paste0(getwd(),"/uerra_d1.nc")
# file=paste0(getwd(),"/Rain_era5_2010-2018.nc")
# filer=paste0(getwd(),"/PressurePrecipitation.nc")
# filew=paste0(getwd(),"/Wind_gust_ERA5_2018.nc")
fileele=paste0(getwd(),"/elev_0.1deg.nc")
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
min

ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  theme_bw(16)+
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.2)+
  geom_raster(data=alelu,aes(x=Var1,y=Var2,fill=ele1,group=gr),alpha=.8,interpolate = T) +
  scale_fill_gradientn(colours = terrain.colors(100),na.value = "aliceblue")+
  geom_polygon(fill = "transparent", color = "gray10", size = 1.2)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="italic"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.position = "none",
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  geom_raster(data=garea,aes(x=Var1,y=Var2,group=Var1),alpha=.8,interpolate = F,fill="red") +
  scale_y_continuous(
    breaks = c(48,50,52,54,56,58,60),limits = c(40,70),"Latitude")+
  scale_x_continuous(
    breaks =c(-6,-4,-2,0,2),limits=c(-10,10),"Longitiude") 
#################################################################################

filer=c(paste0(getwd(),"/in/raindat_7986.nc"),paste0(getwd(),"/in/raindat_8797.nc"),paste0(getwd(),"/in/raindat_9808.nc"),paste0(getwd(),"/in/raindat_0919.nc"))
filew=c(paste0(getwd(),"/in/windat_7986.nc"),paste0(getwd(),"/in/windat_8797.nc"),paste0(getwd(),"/in/windat_9808.nc"),paste0(getwd(),"/in/windat_0919.nc"))

Startdate=as.POSIXct("1979-01-01 10:00:00")
Enddate=as.POSIXct("1986-12-31 23:00:00")
ncr = nc_open(filer)
ncw = nc_open(filew)
# ncatt_get(ncr,name.var,"long_name")
# ncatt_get(ncw,name.var,"long_name")
haz<-c(filer,filew)
hazmat<-c()
timeco<-c()
hx=1
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
  # 
  # tustr <- strsplit(tunits$value, " ")
  # tdstr <- strsplit(unlist(tustr)[3], "-")
  # ori<-as.Date(as.character(tdstr))
  # tmonth <- as.integer(unlist(tdstr)[2])
  # tday <- as.integer(unlist(tdstr)[3])
  # tyear <- as.integer(unlist(tdstr)[1])
  # tcrap<-strsplit(unlist(tustr)[4], ":")
  # thour <- as.integer(unlist(tcrap)[1])
  # time=as.vector(time)
  # time<-time/24
  # chron(time,origin=c(tmonth, tday, tyear,thour))
  
  # nc$data[which(nc$data=="-32767")] <- NA
  
  # nc$data[nc$data==fillvalue$value] <- NA
  # length(na.omit(as.vector(nc$data[,,1])))
  
  
  
  # if(h==1){
    t<-nc$var$tp
    if(length(t)==0) t=nc$var$i10fg
    if(length(t)==0){t<-nc$var$p0001}

    tsize<-t$varsize
    tdims<-t$ndims
    nt1<-tsize[tdims]
    
  # }
  
  
  # if(h==2){
  #   t=nc$var$i10fg
  #   if(length(t)==0){t<-nc$var$p0001}
  #   tsize<-t$varsize
  #   tdims<-t$ndims
  #   nt1<-tsize[tdims]
  #   
  #   
  # }
  # 
  #CHOOSE MY SPATIAL COVERAGE
  #
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
  
  
  # tmp_array <- ncvar_get(nc,t, start = start, count= count)
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
  # temp11<-temp11[,seq(length(temp11[1,]),1)]
  # 
  # image(newformat$lon,newformat$lat,temp11,col=rev(brewer.pal(10,"RdBu")))
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
  # if (hx==2){ 
  # if(s==1)save(newformat,file="interdat/windP1.Rdata")
  # if(s==2)save(newformat,file="interdat/windP2.Rdata")
  # if(s==3)save(newformat,file="interdat/windP3.Rdata")
  # if(s==4)save(newformat,file="interdat/windP4.Rdata")
  # }
  # if(hx==1){
  # if(s==1)save(newformat,file="interdat/rainP1.Rdata")
  # if(s==2)save(newformat,file="interdat/rainP2.Rdata")
  # if(s==3)save(newformat,file="interdat/rainP3.Rdata")
  # if(s==4)save(newformat,file="interdat/rainP4.Rdata")
  # }
  timeco<-c(timeco,list(time))
  if(s==1)save(time,file="interdat/timeP1.Rdata")
  if(s==2)save(time,file="interdat/timeP2.Rdata")
  if(s==3)save(time,file="interdat/timeP3.Rdata")
  if(s==4)save(time,file="interdat/timeP4.Rdata")
}


rbPal <- colorRampPalette(c('darkgreen',"gold","darkorange",'red'))
thx=T
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
#thrx
}

thrr<-as.vector(thr) 
min(thrr)
elonlat <- as.matrix(expand.grid(hazmat[[1]]$lon,hazmat[[1]]$lat))
thrbg<-data.frame(elonlat,thrr)
thrbg$gr=1
ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  theme_bw(16)+
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.2)+
  geom_raster(data=thrbg,aes(x=Var1,y=Var2,fill=thrr,group=gr),alpha=.8,interpolate = F) +
  scale_fill_gradientn(colours = rbPal(100),na.value = "aliceblue")+
  geom_polygon(fill = "transparent", color = "gray10", size = 1.2)+
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


# save(thr,file="Wnd_99_AllP.Rdata")

#Now I have the threshold vectors. Now the tricky part is the selection of the events.
#Maybe I can do it one by one and then group everything 
#I will do one by one and then group and merge clusters if they are less the one hour diff

load("Wnd_99_AllP.Rdata")
thw<-thr
load("Rain_99_AllP.Rdata")
rm(nc,newformat,hazmat)
gc()
metaHaz<-list()
metavHour<-list()
metavDaz<-list()
hx=1
wala=F
for(hazard in 1:4){
  if(hx==2){
    if(hazard==1){load("interdat/windP1.Rdata");time=timeco[[1]]}
    if(hazard==2){load("interdat/windP2.Rdata");time=timeco[[2]]}
    if(hazard==3){load("interdat/windP3.Rdata");time=timeco[[3]]}
    if(hazard==4){load("interdat/windP4.Rdata");time=timeco[[4]]}
    th=thw
  }
  if(hx==1){
    if(hazard==1){load("interdat/rainP1.Rdata");time=timeco[[1]]}
    if(hazard==2){load("interdat/rainP2.Rdata");time=timeco[[2]]}
    if(hazard==3){load("interdat/rainP3.Rdata");time=timeco[[3]]}
    if(hazard==4){load("interdat/rainP4.Rdata");time=timeco[[4]]}
    th=thr
  }

  wala=F
  hazdat=newformat
  rm(newformat)
  gc()
  lon=hazdat$lon
  lat=hazdat$lat
  
  lonlatime <- expand.grid(lon, lat,time)
  
    
    vectouf<- as.vector(hazdat$data)
    vecthouf<-as.vector(rep(th,length(hazdat$data[1,1,])))
    vectouf[which(vectouf<vecthouf)]<-NA
    vectouf[which(vectouf>=vecthouf)]<-1
    vecmeta<-as.vector(hazdat$data)
    # vecmeta2<-as.vector(formeta2)
    if(hx==1){vecmeta[which(vecmeta<vecthouf)] <- NA}
    if(hx==2){vecmeta[which(vecmeta<vecthouf)]<- NA}
  
  lonlatemp <- data.frame(cbind(lonlatime,vectouf))
  metav<-data.frame(cbind(lonlatime,vecmeta))
  
  metav<-metav[which(!is.na(metav[,4])),]

  print(length(metav$Var1)) 
  head(metav)
  lonlatemp2<-lonlatemp[which(lonlatemp[,4]==1),]
  rm(lonlatemp,vecthouf,vectouf)
  gc()
  spdata<-lonlatemp2[,-4]
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
  # Need to figure out how to completely smooth this shit
  if (wala==T){
    walabibou<-kNNdist(sampspd,k=10,all=F)
    min(walabibou)
    walabibou<-walabibou[which(walabibou<10)]
    walabibof=jitter(as.vector(walabibou),15)
    walabibord<-seq(1:length(walabibou))
    walabibix<-data.frame(as.vector(walabibof[order(walabibof)]))
    walabibix<-data.frame(as.vector(walabibou[order(walabibou)]))
    
    walabibix$order=walabibord
    
    # plot(walabibix)
    # walaboss<-unique(walabibix[,1])
    # 
    # ord<-c()
    # for (wa in (1:length(walaboss))){
    #   wa<-walaboss[wa]
    #   id<-which(!is.na(match(walabibix$as.vector.walabibou.order.walabibou...,wa)))
    #   val<-mean(walabibix[id,2])
    #   ord<-c(ord,val)
    # }
    # walabibun<-cbind(ord,walaboss)
    # wlid<-which(walabibix[,1]==walaboss)
    # walabibun<-walabibix[wlid,]
    plot(walabibun)
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
    # lacurve2<-lacurve
    # lacurve2[,1]<-lacurve[,1]/max(lacurve[,1])
    # plot(lacurve)
    # points(x_values,y_values,col=3)
    # points(merd2[1],merd2[2],col=2,pch=16)
    # points(merd2[3],merd2[4],col=2,pch=16)
    # lines(c(merd2[1],merd2[3]),c(merd2[2],merd2[4]))
    # 
    # walabilight[,2]<-walabilight[,2]
    # lacurve[,1]<-lacurve[,1]
    
    
    
    
    
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
    # lines(lo, col='red', lwd=2)
    # points(lo$data$x, lo$data$y,pch=16,col="red",cex=1.5)
    txt<-bquote(epsilon == .(round(eps,2)))
    text(2e5,2.5,cex=1.8, labels=txt,pos=3,col="red")
    abline(h=merd2[2], col=2,lwd=3)
    
    
    # walabibou[order(walabibou)]
    #
  }
  
  if (hx==1)eps=2.45
  if (hx==2)eps=2.24

  epcl=eps

  
  
  # 
  # reservo <- optics(sampspd,eps=3, minPts = 15)
  # res=extractDBSCAN(reservo, eps_cl = epcl)
  # plot(res)
  # unique(res$cluster)
  # bip<-extractXi(object=res, xi=.2)
  # reach=as.reachability(bip)
  # dend<-as.dendrogram(reach)
  # 
  
  ##Update on epcl for rain events
 

  rpip<-dbscan(sampspd, eps=epcl, minPts = 10)
  
  length(unique(rpip$cluster))

  spdata<-cbind(spdata,rpip$cluster)
  if(length(which(spdata$`rpip$cluster`==0))>0){
    metav<-metav[-which(spdata[,4]==0),]
    spdata<-spdata[-which(spdata[,4]==0),]
  }
  metav$cluster<-spdata[,4]
  

  spdata$Var3=spdata$Var3-spdata$Var3[1]+1
  length(metav$Var1)


  metav$time<-as_datetime(c(metav$Var3*60*60),origin="1900-01-01")
  metav$month=month(metav$time)
  event<-metav
  charloc<-paste(event[,1],event[,2])
  event$cloc=charloc
  print(length(event$Var1))
  
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
  # Rcart=strsplit(maxR," ")
  # Rcar<-as.data.frame(matrix(unlist(Rcart),ncol=2,byrow=T))
  # Rcar[,1]<-as.numeric.factor((Rcar[,1]))
  # Rcar[,2]<-as.numeric.factor((Rcar[,2]))
  # 
  # Wcart=strsplit(maxW," ")
  # Wcar<-as.data.frame(matrix(unlist(Wcart),ncol=2,byrow=T))
  # Wcar[,1]<-as.numeric.factor((Wcar[,1]))
  # Wcar[,2]<-as.numeric.factor((Wcar[,2]))
  # EcWR=sqrt((Rcar[,1]-Wcar[,1])^2+(Rcar[,2]-Wcar[,2])^2)
  # geodist<-distGeo(Rcar,Wcar)/1000
  # plot(EcWR)
  # plot(geodist)
  metave$season=metave$x.month
  
  metave$season[which(metave$x.month<3 | metave$x.month>11)]=1
  metave$season[which(metave$x.month<6 & metave$x.month>2)]=2
  metave$season[which(metave$x.month<10 & metave$x.month>5)]=3
  metave$season[which(metave$x.month<12& metave$x.month>9)]=4
  
  # metave$cv=metave$rf.sd/metave$rf.mean
  # metave$Dmax=geodist
  # metave$lonMR<-Rcar[,1]
  # metave$latMR<-Rcar[,2]
  # metave$lonMW<-Wcar[,1]
  # metave$latMW<-Wcar[,2]
  metave$startev<-startev
  metave$endev<-endev
  metaveN<-metave
  # plot(metaveN$vir.surf,metaveN$x.dur)
  length(unique(metaveN$ev))
  
  
  
  metaHaz<-c(metaHaz,list(metaveN))
  metavHour<-c(metavHour,list(event))
  metavDaz<-c(metavDaz,list(testev))
}




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

hx=2

if(hx==1){
  save(metaHax,file="interdat/metaclusRain.Rdata")
  save(metavHax,file="interdat/rawclustRain.Rdata")
  save(metavDax,file="interdat/interclustRain.Rdata")
}
if(hx==2){
  save(metaHax,file="interdat/metaclustWind.Rdata")
  save(metavHax,file="interdat/rawclustWind.Rdata")
  save(metavDax,file="interdat/interclustWind.Rdata")
}
metaHaf<-rbind(metaHax[[1]],metaHax[[2]],metaHax[[3]],metaHax[[4]])
metavHaf<-rbind(metavHax[[1]],metavHax[[2]],metavHax[[3]],metavHax[[4]])
metavDaf<-rbind(metaHax[[1]],metaHax[[2]],metaHax[[3]],metaHax[[4]])
if(hx==2){
save(metaHaf,file="out/WindEv_meta_1979-2019.Rdata")
save(metavHaf,file="out/WindEv_hdat_1979-2019.Rdata")
save(metavDaf,file="out/WindEv_ldat_1979-2019.Rdata")
}
if(hx==1){
save(metaHaf,file="out/RainEv_meta_1979-2019.Rdata")
save(metavHaf,file="out/RainEv_hdat_1979-2019.Rdata")
save(metavDaf,file="out/RainEv_ldat_1979-2019.Rdata")
}



rm(lonlatime,ole,nvo,newformat,arain)
gc()

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
load(file="interdat/metaclusRain.Rdata")
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
#last join
metaHax[[ha]]<-data.frame(metaHax[[ha]],rainacev[,c(3,4)])
plot(metaHax[[ha]]$vir.max,metaHax[[ha]]$len.max) 
abline(a=0,b=1)
}

metaHaf<-rbind(metaHax[[1]],metaHax[[2]],metaHax[[3]],metaHax[[4]])
save(metaHaf,file="out/RainEv_metap_1979-2019.Rdata")
evp<-metaHaf$ev[which.is.max(metaHaf$len.max)]
evento<-metavHaf[which(metavHaf$cluster==3250),]

#Quick event visualizer
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








hx=1
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

save(nwo3,file="out/Wind_stfprint.Rdata")
save(nvo3,file="out/Rain_stfprint.Rdata")
