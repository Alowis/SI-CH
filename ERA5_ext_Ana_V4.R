
rm(list=ls())  
getwd()
setwd("C:/Users/PhD Student/OneDrive - King's College London/DATA/04_Spatiotemporal")
setwd("C:/Users/k1638615/King's College London/OneDrive - King's College London/DATA/04_Spatiotemporal")
config_file=paste0(getwd(),"/config/config_general.R")
source(config_file)
library(ncdf4)

#===========================================================================================================

#I have to create my spatiotemporal cube


##############################################
#What i could do
# Do a sample for OR and find wind dominated, rain dominate and compound events
#1- for each grid cell: how many time within an event
#2- for each event: how long is it within the event (if not contimuous just add)
#3- for each event: what is its size
#4- Use some quantiles to display results in map


Startdate=as.POSIXct("2009-01-01 10:00:00")
Enddate=as.POSIXct("2018-12-31 23:00:00")
longlims=c(-6, 2)
latlims=c(48,59)
mapUK = SpacializedMap(regions = c("UK","France","Spain","Portugal","Italy","Ireland","Belgium","Netherland"))


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
  scale_y_continuous(
    breaks = c(48,50,52,54,56,58,60),limits = c(40,70),"Latitude")+
  scale_x_continuous(
    breaks =c(-6,-4,-2,0,2),limits=c(-10,10),"Longitiude") 
#################################################################################

filer=paste0(getwd(),"/Rain_2009-2019.nc")
filew=paste0(getwd(),"/windG_2009-2019.nc")
ncr = nc_open(filer)
ncw = nc_open(filew)
ncatt_get(ncr,name.var,"long_name")
ncatt_get(ncw,name.var,"long_name")
haz<-c(filer,filew)
hazmat<-c()

for (h in 1:2){
nc = nc_open(haz[h])
nc$var$i10fg
nc$var$p0001
nc$var$p0005
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



if(h==1){
  t<-nc$var$tp
  t<-nc$var$p0001
  tsize<-t$varsize
  tdims<-t$ndims
  nt1<-tsize[tdims]

}


if(h==2){
  t=nc$var$i10fg
  t<-nc$var$p0001
  tsize<-t$varsize
  tdims<-t$ndims
  nt1<-tsize[tdims]


}

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
}



dim(hazmat[[1]]$data)

library(chron)
library(lattice)
library(RColorBrewer)

tustr <- strsplit(tunits$value, " ")
# tdstr<-strsplit(unlist(tustr)[3], "  ")
# 
# ori<-as.Date(as.character(tdstr))
# timi <- strsplit(unlist(tdstr),"-")
# tmonth<-as.integer(unlist(timi)[2])
# tday <- as.character(unlist(timi)[3])
# tday2<-strsplit(tday, "T")
# tday<-as.integer(unlist(tday2)[1])
# tyear <- as.integer(unlist(timi)[1])
# tcrap<-strsplit(unlist(tustr)[4], ":")
# thour1<-tday2[[1]][2]
# thour=0
# tmin=0
# tsec=0
# class(time) = c('POSIXt','POSIXct')
# time2=time
# time3<-which(time2=="2018-05-26 07:00:00 GMT")
# newformat$data[newformat$data=="-32767"] <- NA

library(dbscan)
library(dismo)
library(maps)
library(progress)

rbPal <- colorRampPalette(c('darkgreen',"gold","darkorange",'red'))
names(hazmat)=c("Pr","Wg")
#Set up a quantile for each grid cell
thr<-hazmat$Pr$data[,,1]
for (i in (1:33)){
  for(j in 1:45){
crappy<-hazmat$Pr$data[i,j,]
thr[i,j]<-quantile(crappy[which(crappy>0)],.99,na.rm=T)
  }
}

thw<-hazmat$Wg$data[,,1]
for (i in (1:33)){
  for(j in 1:45){
    crappy<-hazmat$Wg$data[i,j,]
    thw[i,j]<-quantile(crappy[which(crappy>0)],.99,na.rm=T)
  }
}
thrr<-as.vector(thw) 
elonlat <- as.matrix(expand.grid(hazmat$Pr$lon,hazmat$Pr$lat))
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

th1<-quantile(hazmat$Pr$data[which(hazmat[[1]]$data>0)],0.99,na.rm=T)
th2<-quantile(hazmat$Wg$data[which(hazmat[[2]]$data>0)],0.99,na.rm=T)

#####Spatiotemporal clustering#####

metaHaz<-list()
metavHour<-list()
metavDaz<-list()
for(hazard in 1:2){
if(hazard==1){hazdat=hazmat$Pr;th=thr}
if(hazard==2){hazdat=hazmat$Wg;th=thw}
if(hazard==3){hazdat1=hazmat$Pr;hazdat2=hazmat$Wg}

lon=hazmat$Pr$lon
lat=hazmat$Pr$lat
lonlatime <- expand.grid(lon, lat,time)

formeta<-hazmat$Pr$data
formeta2<-hazmat$Wg$data
if (hazard == 3){

  bolilos[hazdat1$data<th1 | hazdat2$data<th2] <- NA
  bolilos[hazdat1$data>=th1 & hazdat2$data>=th2] <- 1  
  formeta[hazdat1$data<th1 | hazdat2$data<th2] <- NA
  formeta2[hazdat1$data<th1 | hazdat2$data<th2] <- NA
}else
  {
  bolilos<-hazdat$data 

  # boli<-bolilos[,,t]
  # bol<-melt(boli)
  # bol[which(bol[,3]<th),3]<-NA
  # fmt<-array(as.matrix(bol$value), dim=c(33,45))
  # formeta[,,t]<-fmt
  # bol[which(bol[,3]>=th),3]<-1
  # bolo<-array(as.matrix(bol$value), dim=c(33,45))
  # bolilos[,,t]<-bolo
  # bolilos[hazdat$data[,,t]<th] <- NA
  vectouf<- as.vector(hazdat$data)
  vecthouf<-as.vector(rep(th,length(bolilos[1,1,])))
  vectouf[which(vectouf<vecthouf)]<-NA
  vectouf[which(vectouf>=vecthouf)]<-1
  vecmeta<-as.vector(formeta)
  vecmeta2<-as.vector(formeta2)
  if(hazard==1){vecmeta[which(vecmeta<vecthouf)] <- NA}
  if(hazard==2){vecmeta2[which(vecmeta2<vecthouf)]<- NA}
}
# vectouf<- as.vector(hazdat$data)
# vecthouf<-as.vector(rep(thr,length(bolilos[1,1,])))
# vectouf[which(vectouf<vecthouf)]<-NA
# length(na.omit(vectouf))
# vrac<-as.vector(bolilos)
# length(na.omit(vrac))
# vecmeta<-as.vector(formeta)
# vecmeta2<-as.vector(formeta2)
lonlatemp <- data.frame(cbind(lonlatime,vectouf))
metav<-data.frame(cbind(lonlatime,vecmeta,vecmeta2))


x=abs(rnorm(100*100,50,25))
x=matrix(x,nrow=100)
x1=melt(th)


if(hazard==1)metav<-metav[which(!is.na(metav[,4])),]
if(hazard==2)metav<-metav[which(!is.na(metav[,5])),]
if(hazard==3)metav<-metav[which(!is.na(metav[,5])),]

head(metav)
lonlatemp2<-lonlatemp[which(lonlatemp[,4]==1),]
spdata<-lonlatemp2[,-4]
ep<-2
sampspd<-spdata
sampspd$Var1<-sampspd$Var1*3
sampspd$Var2<-sampspd$Var2*3
sampspd$Var3<-sampspd$Var3-sampspd$Var3[1]+1
sampspd$Var3<-sampspd$Var3
samptt<-as.matrix(sampspd)

# walabibou<-kNNdist(sampspd,k=20,all=F)
# plot(walabibou[order(walabibou)],ylim=c(0,4))

sqrt(0.5^2+0.5^2+1^2)
epcl=2
epl<-sqrt(.75^2+.75^2+1)
# abline(h=epcl, col=2) 

# 
# reservo <- optics(sampspd,eps=100, minPts = 10)
# res=extractDBSCAN(reservo, eps_cl = ep)
# bip<-extractXi(object=res, xi=.2)
# reach=as.reachability(bip)
# dend<-as.dendrogram(reach)
# 

if(hazard==1)weightc=metav$vecmeta
if(hazard==2)weightc=metav$vecmeta2
if(hazard==3)weightc=rep(1,length(sampspd[,1]))
rpip<-dbscan(sampspd, eps=epcl, minPts = 15,weights = weightc)
unique(rpip$cluster)
# 
# nn<-frNN(sampspd,eps=epcl,sort=TRUE)
# 
# i <- 100
# nn$id[[i]]
# nn$dist[[i]]
# plot(sampspd[,c(1,2)], col = ifelse(1:nrow(sampspd) %in% nn$id[[i]], "red", "black"))
# 
# hist(sapply(adjacencylist(frchier), length),
#      xlab = "k", main="Number of Neighbors",
#      sub = paste("Neighborhood size eps =", frchier$eps))
# 
# nn <- list(ids = list(c(2,3), c(1,3), c(1,2,3), c(3,5), c(4,5)), eps = 1)
# class(nn) <- c("NN", "frNN")
# nn
# head(adjacencylist(nn))
# dbscan(nn, minPts = 2)
# library(cluster)
# library(subspace)
# res1 <- extractDBSCAN(reservo, eps_cl = epcl)
# unique(res1$cluster)
# plot(res1)

# length(unique(res1$cluster))
# # hullplot(yip,res)
# plot(spdata[which(res1$cluster==87),c(1,2) ],col=res1$cluster)

spdata<-cbind(spdata,rpip$cluster)
if(length(which(spdata$`rpip$cluster`==0))>0){
metav<-metav[-which(spdata[,4]==0),]
spdata<-spdata[-which(spdata[,4]==0),]
}
metav$cluster<-spdata[,4]

length(spdata[,4])
length(metav[,4])
spdata$Var3=spdata$Var3-spdata$Var3[1]+1


# for(m in unique(spdata[,4])){
metav$time<-as_datetime(c(metav$Var3*60*60),origin="1900-01-01")
metav$month=month(metav$time)
  event<-metav
  charloc<-paste(event[,1],event[,2])
  event$cloc=charloc

  testev<-aggregate(list(rf= event[,4],wg=event[,5]),
                    by = list(ev = event[,6],loc=event[,9]),
                    FUN = function(x) c(sum = sum(na.omit(x)),max=max(na.omit(x))))
  testev <- do.call(data.frame, testev)
  

  
  metamax<-aggregate(list(rf= testev[,3],wg=testev[,6]) ,
                     by = list(ev = testev[,1]),
                     FUN = function(x) c(sum = sum(na.omit(x)),max=max(na.omit(x)),mean=mean(x),sd=sd(x),surf=length(x)))
  
  metamax<- do.call(data.frame, metamax)
  
  metave<-aggregate(list(rf= event[,4],wg=event[,5]) ,
                      by = list(ev = event[,6]),
                      FUN = function(x) c(surf=length(x)))
  metave <- do.call(data.frame, metave)

  
  tempcom<-aggregate(event[,7] ,
                     by = list(ev = event[,6]),
                     FUN = function(x) c(dur=length(unique(x)),month=month(unique(x))[1],year=year(unique(x))[1]))

  tempic<- do.call(data.frame, tempcom)
  
  metave<-cbind(metamax,metave[,c(2,3)],tempic[,c(2,3,4)])
maxR<-c()
maxW<-c()
  for (eve in unique(testev$ev))
  {
    evint<-testev[which(testev$ev==eve),]
    met<-metamax[which(metamax$ev==eve),]
    elR<-as.character(evint$loc[which(evint$rf.sum==met$rf.max)])[1]
    elW<-as.character(evint$loc[which(evint$wg.max==met$wg.max)])[1]
    maxR<-c(maxR,elR)
    maxW<-c(maxW,elW)
    
  }
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
  Rcart=strsplit(maxR," ")
  Rcar<-as.data.frame(matrix(unlist(Rcart),ncol=2,byrow=T))
  Rcar[,1]<-as.numeric.factor((Rcar[,1]))
  Rcar[,2]<-as.numeric.factor((Rcar[,2]))
  
  Wcart=strsplit(maxW," ")
  Wcar<-as.data.frame(matrix(unlist(Wcart),ncol=2,byrow=T))
  Wcar[,1]<-as.numeric.factor((Wcar[,1]))
  Wcar[,2]<-as.numeric.factor((Wcar[,2]))
  EcWR=sqrt((Rcar[,1]-Wcar[,1])^2+(Rcar[,2]-Wcar[,2])^2)
  geodist<-distGeo(Rcar,Wcar)/1000
  plot(EcWR)
  plot(geodist)
  metave$season=metave$x.month

  metave$season[which(metave$x.month<3 | metave$x.month>11)]=1
  metave$season[which(metave$x.month<6 & metave$x.month>2)]=2
  metave$season[which(metave$x.month<10 & metave$x.month>5)]=3
  metave$season[which(metave$x.month<12& metave$x.month>9)]=4

  metave$cv=metave$rf.sd/metave$rf.mean
  metave$Dmax=geodist
  metave$lonMR<-Rcar[,1]
  metave$latMR<-Rcar[,2]
  metave$lonMW<-Wcar[,1]
  metave$latMW<-Wcar[,2]
  metaveN<-metave[which(metave$x.dur>2),]
  hist(metaveN$Dmax)
  plot(metaveN$Dmax,metaveN$x.dur)
  


  
plot(metaveN$rf.max,metaveN$wg.max)
metaHaz<-c(metaHaz,list(metaveN))
metavHour<-c(metavHour,list(metav))
metavDaz<-c(metavDaz,list(testev))
}


#Retain only events that affected lands

metaHax<-list()
metavHax<-list()
metavDax<-list()
for (hx in 1:2){
idout<-c()
for (cl in metaHaz[[hx]]$ev){ 
  clev<-metavHour[[hx]][which(metavHour[[hx]]$cluster==cl),]
  overlappy<-point.in.polygon(clev[,c(1)],clev[,c(2)],uk_fort$long,uk_fort$lat)
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
############################################
#Then select of pairs of events that overlap spatially
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
oc<-c()
nvo<-c()
nwo<-c()
for (cl in metaHax[[2]]$ev){ 
  clev<-metavDax[[2]][which(metavDax[[2]]$ev==cl),]
  dd<-unique(metavHax[[2]]$time[which(metavHax[[2]]$cluster==cl)])
  
  clex<-do.call(rbind, replicate(length(dd),clev[,c(1,2)], simplify=FALSE))
  cley<-do.call(rbind, replicate(length(clev[,1]),as.data.frame(dd),simplify=F))
  ouch<-cbind(clex,cley)
  nwo<-rbind(nwo,ouch)
  print(cl)
}


save(nvo,file="rainallclusters.Rdata")
save(nwo,file="windallclusters.Rdata")

#####MADE IT!!!####
  cclev<-metavHax[[2]]
  clev$charloc<-paste(clev$Var1,clev$Var2)
  cclev$charloc<-paste(metavHax[[2]]$Var1,metavHax[[2]]$Var2)
  sp10<-match_df(nvo,nwo,on=c("loc","dd"))
  sp03<-inner_join(nwo,nvo,by=c("loc","dd"))
  
  tesp<-aggregate(list(shit=sp03[,3]),
                    by = list(e1=sp03[,1],e2 = sp03[,4],l=sp03[,2]),
                    FUN = function(x) c(length=length(unique(x))))
  tesp <- do.call(data.frame, tesp)
  
  
  endgame<-aggregate(list(evd=tesp[,3]),
                     by = list(e1=tesp[,1],e2=tesp[,2]),
                     FUN = function(x) c(length=length(x)))
  endgame <- do.call(data.frame, endgame)
  tespx<-tesp[,c(1,2)]
  
  names(tespx)=c("dd","loc")
  sp11<-join(tespx,sp03,by=c("loc","dd"),match="first")
  matchev<-unique(sp11[c(3,4)])
  
  endgame<-aggregate(list(evd=sp03[,2]),
                  by = list(ev=sp03[,1],sp03[,4],sp03[,3]),
                  FUN = function(x) c(length=length(x)))
  endgame <- do.call(data.frame, endgame)
  max(endgame$evd)
  max(metaHax[[1]]$wg.surf)
  shi<-unique(sp03[,c(2,3)])
  krak<-(which(sp03[,c(2,3)]==shi))
  sp11<-match_df(sp03,tespx,on=c("loc","dd"))
  walla<-sp03[,c(7,15)]
  assad<-which(walla[,1]==walla[,2])
  wallaz<-sp03[assad,]
  length(unique(sp10$ev))

  sp07<-data.frame(clev$time,sp01)
  sp08<-data.frame(cclev$time,sp02)
  names(sp07)=names(sp08)


  crap<-unique(sp03[,6])
  crap2<-unique(clev$cluster)
  sp09<-unique(as.numeric(which(!is.na((match(sp07,sp08))))))
  sp04<-metavHax[[2]][sp03,]  
  sp05<-as.numeric(which(!is.na((match(sp04$time,clev$time)))))
  graal<-sp04[sp05,]
  graal$clr<-cl
  oc<-rbind(oc,graal)
  
}

names(metaHaz)=c("Pr-Events","Wind_Events","Compound")
names(metavHour)=c("Pr-Events","Wind_Events","Compound")

#############Analysis###################
library(viridis)
rbPal <- colorRampPalette(c('darkgreen',"gold","darkorange",'red'))

mapUK = SpacializedMap(regions = c("UK","France","Spain","Portugal","Italy","Ireland"))


ukk<- gSimplify(mapUK, tol=0.001, topologyPreserve=TRUE)
plot(ukk)

plot(metaHaz[[1]]$lonMR,metaHaz[[1]]$latMR)

rbPal2<-colorRampPalette(c('white',"blue","purple","red","orange"))
#Hotspot locations
ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill = "white", color = "gray10", size = 1) +
  theme_bw(16)+
  scale_fill_gradientn(colours = rbPal2(100))+
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.2)+
  # geom_bin2d(data=metaHaz[[2]],aes(x=lonMW,y=latMW,group=season),alpha=.7,binwidth=c(.5,.5))+
# geom_density2d(data=metaHaz[[2]],aes(x=lonMW,y=latMW,group=season),alpha=.2)
stat_density_2d(data=metaHaz[[1]],aes(x=lonMR,y=latMR,group=season,fill = stat(density)),geom = "raster", contour = FALSE,alpha=.6,n=150)
allprec<-metavHour[[2]]
allprec$season=1
allprec$season[which(allprec$month<6 & allprec$month>2)]=2
allprec$season[which(allprec$month<10 & allprec$month>5)]=3
allprec$season[which(allprec$month<12& allprec$month>9)]=4
ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill = "white", color = "gray10", size = 1) +
  theme_bw(16)+
  scale_fill_gradientn(colours = rbPal2(100))+
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.2)+
  # geom_bin2d(data=metaHaz[[2]],aes(x=lonMW,y=latMW,group=season),alpha=.7,binwidth=c(.5,.5))+
  # geom_density2d(data=metaHaz[[2]],aes(x=lonMW,y=latMW,group=season),alpha=.2)
  stat_density_2d(data=allprec,aes(x=Var1,y=Var2,group=season,fill = stat(density)),geom = "raster", contour = FALSE,alpha=.6,n=150)


hist(metaHaz[[3]]$lonMR)

rbPal2 <- colorRampPalette(c('red',"darkorange","gold",'darkgreen'))
Col2 <- rbPal2(1000)[as.numeric(cut(log(metaHaz[[1]]$Dmax+0.00001),breaks = 1000))] 
plot(metaHaz$`Pr-Events`$wg.max,metaHaz$`Pr-Events`$rf.max,col=alpha(Col2,.4),pch=16,log="x")




####Polat plot

polarEvent<-aggregate(metar[[1]]$season,
                 by = list(Month = metar[[1]]$x.month,Year=metar[[1]]$x.year),
                 FUN = function(x) c(n =length(x)))
polarEventR <- do.call(data.frame, polarEvent)
polarEventR$date=as.character(paste0(polarEventR$Year,"-",sprintf("%02d", as.numeric(polarEventR$Month))))
polarEventR$date=as.Date(paste(polarEventR$date,"-01",sep=""))
polarEventR$monthstr=months(polarEventR$date)

polarEventR$monthstr = factor(polarEventR$monthstr, levels=month.name)

coord_radar <- function (theta = "x", start = 0, direction = 1) {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}

rbPal2<-colorRampPalette(c('skyblue',"orange","darkorange","red","purple"))
polarEventR$gr=1
library(ggplot2)
data(beavers)
plot(polarEventR[,c(1,3)])
p <- ggplot() + 
  geom_path(data=polarEventR, aes(x=monthstr,y=x,group=gr,colour=Year),size=1)+
  geom_line(data=polarEventR, aes(x=monthstr,y=x,colour = Year,group=Year),size=1) + 
  geom_point(data=polarEventR, aes(x=monthstr,y=x,colour = Year,group=Year),size=2)+
  scale_colour_gradientn(colours = rbPal2(30))
p+ theme_bw() +coord_radar()



#Seasonal analysis

bibis<-aggregate(metaHaz[[1]]$season,
                 by = list(Month = metaHaz[[1]]$x.month),
                 FUN = function(x) c(n =length(x),x=x[1]))
bibis <- do.call(data.frame, bibis)
bibis$x.x[which(bibis$Month==9)]=4
dodge <- position_dodge(width = 1)

ggplot(bibis,aes(x=Month, y=x.n,fill=factor(x.x)))+geom_bar(stat = "identity", position = dodge,alpha=.6) +scale_fill_manual(values = c( "blue", "red","orange","skyblue")) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c("Jan","Fev","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))+
  labs(y = "# of occurences",x = "Month",colour = "Parameter",fill= "Season",size=20)+
  theme_classic() +
  theme(axis.text.x = element_text(angle=0, vjust = 0.5,size=16))+
  theme(text = element_text(size=16))


hist(metaHaz$`Pr-Events`$wg.max)
doubleExRevent<-metaHaz$`Pr-Events`[which(metaHaz$`Pr-Events`$wg.max>th2),]
doubleExWevent<-metaHaz[[2]][which(metaHaz[[2]]$rf.mean>th1),]

stseason<-metaHaz[[1]]
var.int=stseason$x.dur
hist(var.int,breaks=c(0,6,12,24,96))
stseason$sizegr<-1
stseason$sizegr[which(stseason$x.dur>6 & stseason$x.dur<13)]=2
stseason$sizegr[which(stseason$x.dur>12 & stseason$x.dur<25)]=3
stseason$sizegr[which(stseason$x.dur>24)]=4



gr1<-stseason[which(stseason$x.dur<7)]
gr2<-stseason[which(stseason$x.dur>6 & stseason$x.dur<13),]
gr3<-stseason[which(stseason$x.dur>12 & stseason$x.dur<25),]
gr4<-stseason[which(stseason$x.dur>24),]

plot(gr1$wg.max,gr1$rf.max,xlim=c(0,50),ylim=c(0,130))
points(gr2$wg.max,gr2$rf.max,col=2)
points(gr3$wg.max,gr3$rf.max,col=3)
points(gr4$wg.max,gr4$rf.max,col=4)


plot(gr1$rf.surf,gr1$x.dur,log="xy",xlim=c(1,10000),ylim=c(1,130))
points(gr2$rf.surf,gr2$x.dur,col=2)
points(gr3$rf.surf,gr3$x.dur,col=3)
points(gr4$rf.surf,gr4$x.dur,col=4)

########################Self organized map#####################
# Load the kohonen package 
require(kohonen)

# Create a training data set (rows are samples, columns are variables
# Here I am selecting a subset of my variables available in "data"
stseason<-metaHaz[[1]]
stseason$dispr<-stseason$rf.max/stseason$rf.sum
plot(stseason$dispr)
data_train <- (stseason[c(17,20,21,24)])
data_train$rf.surf=log(data_train$rf.surf)
data_train$x.dur=log(data_train$x.dur)
data_train$w=ifelse(data_train$season==1,1,0)
data_train$p=ifelse(data_train$season==2,1,0)
data_train$s=ifelse(data_train$season==3,1,0)
data_train$a=ifelse(data_train$season==4,1,0)
data_train<-data_train[,-1]
# Change the data frame with training data to a matrix
# Also center and scale all variables to give them equal importance during
# the SOM training process. 
data_train_matrix <- as.matrix(scale(data_train))

# Create the SOM Grid - you generally have to specify the size of the 
# training grid prior to training the SOM. Hexagonal and Circular 
# topologies are possible
dims<-sqrt(5*sqrt(length(data_train[,1])))
dims=10
som_grid <- somgrid(xdim = round(dims), ydim=round(dims), topo="hexagonal")

# Finally, train the SOM, options for the number of iterations,
# the learning rates, and the neighbourhood are available
som_model <- som(data_train_matrix, 
                 grid=som_grid, 
                 rlen=200, 
                 keep.data = T)
length(na.omit(data_train[,var]))
var <- 6 #define the variable to plot 
var_unscaled <- aggregate(as.numeric(data_train[,var]), by=list(som_model$unit.classif), FUN=mean, simplify=TRUE)[,2] 
plot(som_model, type = "property", property=var_unscaled, main=names(data_train)[var],palette.name = viridis)
var_unscaled <- aggregate(as.numeric(data_train[,var]), by=list(som_model$unit.classif), FUN=mean, simplify=TRUE)[,2] 
plot(som_model, type="changes")
plot(som_model, type = "property", property=var_unscaled, main=names(data_train)[var],palette.name = viridis)
max(var_unscaled)
plot(som_model, type="dist.neighbours")
plot(som_model, type = "codes", 
     property = som_model$codes[[1]][,6], main=names(som_model$data[[1]])[6])

c.class <- kmeans(som_model$codes[[1]], 4)
plot(som_model, type = "codes", 
     bgcol = rbPal2(4)[c.class$cluster])
add.cluster.boundaries(som_model, c.class$cluster)

som_model$unit.classif

## use hierarchical clustering to cluster the codebook vectors
som_cluster <- cutree(hclust(dist(som_model$codes[[1]])), 2)
# plot these results:
plot(som_model, type="codes", bgcol = som_cluster) 
add.cluster.boundaries(som_model, som_cluster)

# get vector with cluster value for each original data sample
cluster_assignment <- som_cluster[som_model$unit.classif]
stseason$somclust<-cluster_assignment
# for each of analysis, add the assignment as a column in the original data:
rain<-stseason[which(stseason$somclust==2),]
wind<-stseason[which(stseason$somclust==4),]
plot(wind$wg.max,wind$rf.max)

colorz<-c("blue", "red","orange","skyblue")
rbk <- colorRampPalette(c("blue", "red","orange","skyblue"))
stseason$sizegr<-1
stseason$sizegr[which(stseason$x.dur>6 & stseason$x.dur<13)]=2
stseason$sizegr[which(stseason$x.dur>12 & stseason$x.dur<25)]=3
stseason$sizegr[which(stseason$x.dur>24)]=4
Colx <- rbk(4)[stseason$season]
Coly <- rbk(4)[stseason$somclust]
Colz <-rbk(4)[stseason$sizegr]
plot(stseason$wg.max,stseason$rf.max,col=Colx)
plot(stseason$wg.max,stseason$rf.max,col=Coly)
plot(stseason$wg.max,stseason$rf.max,col=Colz)
plot(stseason$season,stseason$somclust)
plot(stseason$rf.surf,stseason$x.dur,col=alpha(Coly,.5),log="xy",xlim=c(1,10000),ylim=c(1,130),pch=16,cex=2)
plot(stseason$rf.surf,stseason$x.dur,col=Colx,log="xy",xlim=c(1,10000),ylim=c(1,130))

ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill = "white", color = "gray10", size = 1) +
  theme_bw(16)+
  scale_fill_gradientn(colours = rbPal2(100))+
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.2)+
  # geom_bin2d(data=metaHaz[[2]],aes(x=lonMW,y=latMW,group=season),alpha=.7,binwidth=c(.5,.5))+
  # geom_density2d(data=metaHaz[[2]],aes(x=lonMW,y=latMW,group=season),alpha=.2)
  stat_density_2d(data=stseason[which(stseason$season==4),],aes(x=lonMW,y=latMW,group=season,fill = stat(density)),geom = "raster", contour = FALSE,alpha=.6,n=150)

######################################################################################




mlev<-order(-var.int)
plot(var.int[mlev])
idmt<-stseason$ev[mlev]
plot(metaHaz[[1]]$x.dur,metaHaz[[1]]$wg.max,ylim=c(0.,52))
points(metaHaz[[2]]$x.dur,metaHaz[[2]]$wg.max,col=2)
dates<-c()
loca<-c()
totemp<-c()
haz=1
maxevcomp=c()
for(idi in 1:length(stseason[,1])){
  idm=idmt[idi]
  maxevent<-metavHour[[haz]][which(metavHour[[haz]]$cluster==idm),]
  
  uniquetstep<-aggregate(list(rf=maxevent$vecmeta,wg=maxevent$vecmeta2),
                   by = list(time = maxevent$time),
                   FUN = function(x) c(size = length(x),mx= max(x)))
  ustep <- do.call(data.frame, uniquetstep)
  
  plot(ustep$time, ustep$wg.mx,type="l")
  plot(ustep$time, ustep$rf.mx,typ="l",col=2)
  
  maxcomb<-maxevent[which(maxevent$vecmeta2>th2),]
  maxevcomp<-rbind(maxevcomp,maxcomb)
  maxevent$idloc<-paste(maxevent$Var1,maxevent$Var2)
}
  idloco<-unique(maxevent$idloc)
  metamax<-stseason[which(stseason$ev==idm),]
  
  temp<-(unique(maxevent$Var3))
  tempev<-data.frame(temp, rep(idm, length(temp)))
  totemp<-rbind(totemp, tempev)
  maxevent$time <- as_datetime(c(maxevent$Var3*60*60),origin="1900-01-01")
  dates<-rbind(dates,c(maxevent$time[1],maxevent$time[length(maxevent$time)]))
  loca<-rbind(loca,metamax)



chr <- chull(metaHaz[[1]]$rf,metaHaz[[1]]$x.dur)
chr<-data.frame(metaHaz[[1]]$rf[c(chr,chr[1])],metaHaz[[1]]$x.dur[c(chr,chr[1])])
chw <- chull(metaHaz[[2]]$rf,metaHaz[[2]]$x.dur)
chw<-data.frame(metaHaz[[2]]$rf[c(chw,chw[1])],metaHaz[[2]]$x.dur[c(chw,chw[1])])
chc<-chull(metaHaz[[3]]$rf,metaHaz[[3]]$x.dur)
chc<-data.frame(metaHaz[[3]]$rf[c(chc,chc[1])],metaHaz[[3]]$x.dur[c(chc,chc[1])])

rbPal <- c("blue","red","orange","skyblue")

plot(metaHaz[[1]]$rf,metaHaz[[1]]$x.dur,col=alpha(1,0.4),pch=16,log="xy", xlim=c(1,40000),ylim=c(2.9,120),cex= round(metaHaz[[1]]$rf.max,4)/10)
points(doubleExRevent$rf,doubleExRevent$x.dur,col=alpha(2,.4),pch=18,cex= round(doubleExRevent$rf.max,4)/10)


plot(metaHaz[[1]]$rf,metaHaz[[1]]$x.dur,col=alpha(1,0.6),pch=16,log="xy", xlim=c(1,40000),ylim=c(2.9,120),cex= round(metaHaz[[1]]$wg.max/mean(metaHaz[[1]]$wg.max),4))
points(metaHaz[[2]]$rf,metaHaz[[2]]$x.dur,col=alpha(2,0.6),pch=16,log="xy", xlim=c(1,40000),ylim=c(2.9,120),cex= round(metaHaz[[2]]$wg.max/mean(metaHaz[[2]]$wg.max),4))
points(metaHaz[[3]]$rf,metaHaz[[3]]$x.dur,col=alpha(3,0.6),pch=16,log="xy", xlim=c(1,40000),ylim=c(2.9,120),cex= round(metaHaz[[3]]$wg.max/mean(metaHaz[[3]]$wg.max),4))

points(doubleExWevent$rf,doubleExWevent$x.dur,col=alpha(2,.5),pch=16, xlim=c(1,40000),ylim=c(2.9,120),cex= round(doubleExWevent$rf.max)/20)

plot(metaHaz[[1]]$wg.max,metaHaz[[1]]$rf.max,cex=1.2,pch=16,col=alpha(1,0.6))
points(metaHaz[[2]]$wg.max,metaHaz[[2]]$rf.max,col=alpha(2,0.6),cex=1.1,pch=16)
points(metaHaz[[3]]$wg.max,metaHaz[[3]]$rf.max,col=alpha(3,0.6),cex=1,pch=16)
points(doubleExRevent$wg.max,doubleExRevent$rf.max,col=alpha(4,0.6),pch=17)
points(doubleExWevent$wg.max,doubleExWevent$rf.max,col=alpha(5,0.6),pch=18)


chwr <- chull(doubleExWevent$rf,doubleExWevent$x.dur)
chwr<-data.frame(doubleExWevent$rf[c(chwr,chwr[1])],doubleExWevent$x.dur[c(chwr,chwr[1])])
chrw<-chull(doubleExRevent$rf,doubleExRevent$x.dur)
chrw<-data.frame(doubleExRevent$rf[c(chrw,chrw[1])],doubleExRevent$x.dur[c(chrw,chrw[1])])


plot(chr,col=2,lwd=3,type="l",xlim=c(1,40000),ylim=c(2.9,120), log="xy")
lines(chw,col=3,lwd=3)
lines(chc,col=4,lwd=3)
lines(chrw,col=5,lwd=3)
lines(chwr,col=6,lwd=3)


rbPal <- colorRampPalette(c('lightskyblue',"skyblue","gold","darkorange",'red',"purple"))
Colw <- rbPal(20)[as.numeric(cut((metaHaz[[2]]$wg.max),breaks = 20))]
plot(metaHaz[[2]]$rf,metaHaz[[2]]$x.dur,col=alpha(Colw,0.8),pch=16, xlim=c(1,40000),ylim=c(2.9,120),cex= (metaHaz[[2]]$rf.max)/10,log="xy") 
Colr <- rbPal(20)[as.numeric(cut((metaHaz[[1]]$wg.max),breaks = 20))]
plot(metaHaz[[1]]$rf,metaHaz[[1]]$x.dur,col=alpha(Colr,0.8),pch=16, xlim=c(1,40000),ylim=c(2.9,120),cex= (metaHaz[[1]]$rf.max)/10,log="xy") 

Colc <- rbPal(20)[as.numeric(cut((metaHaz[[3]]$wg.max),breaks = 20))]
legend_image <- as.raster(matrix(rbPal(20)[seq(20,1)], ncol=1))
unique(Colr)
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
par(mar=rep(.5, 4), oma=rep(3, 4), las=1)
plot(metaHaz[[3]]$rf,metaHaz[[3]]$x.dur,col=alpha(Colr,0.8),pch=16, xlim=c(1,40000),ylim=c(2.9,120),cex= (metaHaz[[3]]$rf.max)/10,log="xy") 
AddGradientLegend(breaks, title = "Title", loc = "topright")
legend("topright", inset = -.0,legend=round(seq(min(metaHaz[[3]]$rf.max),max(metaHaz[[3]]$rf.max),l=5)),pt.cex=seq(min(metaHaz[[3]]$rf.max)/10,max(metaHaz[[3]]$rf.max)/10,l=5),pch=16,title="Maximum \n rainfall [mm]",col=alpha("grey",.8))



# plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
# text(x=1.5, y = seq(0,1,l=5), labels = round(seq(min(metaHaz[[3]]$wg.max),max(metaHaz[[3]]$wg.max),l=5)))
# rasterImage(legend_image, 0, 0, 1,1)
# 
# 
metaHaz[[1]]$rf.max[which(metaHaz[[1]]$rf.max<0)]=0
max(metaHaz[[1]]$rf.max)

ggplot(data=metaHaz[[1]],aes(x=rf.surf,y=x.dur,colour=wg.max,size=rf.max))+
geom_point(alpha=.5)+
  theme(axis.text=element_text(size=16),
                axis.title=element_text(size=18,face="bold"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
scale_colour_gradientn(trans=scales::boxcox_trans(1.5),colors=rbPal(100),"Max wind gust [m/s]",breaks=c(10,20,30,40,50,Inf),limits=c(1,51))+
  scale_size(trans=scales::boxcox_trans(1.3),range=c(.4,17),"Max precipitation [mm]",breaks=c(20,60,120,Inf),limits=c(0,140),guide = guide_legend(override.aes = list(colour = alpha("grey",.8))))+
  scale_y_continuous(trans = log_trans(),
                     breaks = c(3,6,12,24,48,72,148),limits = c(2.9,150),"Temporal scale",
                     labels=c("3h","6h","12h","1day","2days","3days","1week"))+
scale_x_continuous(trans = log_trans(),
                   breaks =c(1000,10000,100000,1000000)/450,limits=c(1,3000),"Spatial scale [km2]",
                   labels=c("1000","10000","100000","1000000")) 

boxcox_trans(2)
range(metaHaz[[2]]$wg)
c(1000,10000,100000,1000000,1e7)/961
45*33*961
#add event speed?



legend(5,65, legend=c("DJF", "MAM","JJA","SON"),
       col=rbPal,pch=16 ,cex=0.8)


#######Event cas par cas##############################

#############VALIDATION OF THE CLUSTERING METHOD##################
totdates=list()
for (haz in 1:3){
stseason<-metaHaz[[haz]]
var.int=stseason$rf.max
mlev<-order(-var.int)
plot(var.int[mlev])
idmt<-stseason$ev[mlev]

dates<-c()
for(idi in 1:10){
idm=idmt[idi]
maxevent<-metavHour[[haz]][which(metavHour[[haz]]$cluster==idm),]
maxevent$time <- as_datetime(c(maxevent$Var3*60*60),origin="1900-01-01")
print(idm)
print(c(maxevent$time[1],maxevent$time[length(maxevent$time)]))
dates<-rbind(dates,c(maxevent$time[1],maxevent$time[length(maxevent$time)]))
}
dates<-as.data.frame(dates)
dates[,1]<-as_datetime(dates[,1],origin="1970-01-01")
dates[,2]<-as_datetime(dates[,2],origin="1970-01-01")
names(dates)=c("start","end")
totdates<-c(totdates,list(dates))
}
names(totdates)=c("Pr-Events","Wind_Events","Compound")
plot(totdates[[1]][,1])
points(totdates[[3]][,1],col=2,pch=3)
points(totdates[[2]][,1],col=3,pch=4)

write.csv(totdates[[3]], file = "windcomptop10valid.csv")




idi=10
idm=idmt[idi]
maxevent<-metavHour[[haz]][which(metavHour[[haz]]$cluster==idm),]
maxevent$time <- as_datetime(c(maxevent$Var3*60*60),origin="1900-01-01")
print(idm)
print(c(maxevent$time[1],maxevent$time[length(maxevent$time)]))
idkik=c()
inti=seq(1,length(unique(maxevent$time)))
colvect<-heat.colors(max(inti), alpha = 1)




for (d in 1:length(unique(maxevent$time))){
  kik<-maxevent[which(maxevent$time==unique(maxevent$time)[d]),]
# points(kik$Var1,kik$Var2,col=d)

    centi<-as.numeric(kik[which(kik$vecmeta==max(kik$vecmeta))[1],c(1,2)])
  idkik<-rbind(idkik,c(d,max(kik$vecmeta),length(kik$vecmeta),centi))
  points(centi[1],centi[2],col=colvect[d],pch=16)
}

idkik=data.frame(idkik)


uk_fort <- ggplot2::fortify(ukk)
# ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
#   geom_polygon(fill = "white", color = "gray10", size = 1) +
#   theme_bw(16)+
#   scale_color_viridis(option="C", alpha=.8)+
#   geom_point(data=idkik,aes(x=X4,y=X5,group=X1,color=X1,size=X2))+
#   coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3)


idkik$gr<-1
idkik$rtime<-1
# ggplot(idkik, aes(X4, X5))+
#   geom_path(aes(group=gr), arrow = arrow())+
#   coord_fixed()

f1_mod <- loess(X4 ~ X1, data = idkik,span=0.5)
f2_mod <- loess(X5 ~ X1, data = idkik,span=0.4)
pred <- data.frame(X1 = seq(1,length(idkik$X1), length = 100))
pred$F1 <- predict(f1_mod, newdata = pred)
pred$F2 <- predict(f2_mod, newdata = pred)
pred$G<-1
indices <- seq(1,99, by = 5)
ggplot(pred, aes(F1, F2))+
  geom_path()+
  geom_point(data = pred[indices,])+
  coord_fixed()

# ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
#   geom_polygon(fill = "white", color = "gray10", size = 1) +
#   theme_bw(16)+
#   scale_color_viridis(option="D", alpha=.8)+
#   xlim(-10, 10)+
#   ylim(35,62)+
#   geom_point(data=idkik,aes(x=X4,y=X5,group=X1,color=X2,size=X3),alpha=0.8)+
#   geom_path(data=pred, aes(F1, F2,group=G),linejoin="mitre",colour="red",size=1.5,alpha=0.6)+
#   geom_point(data=pred, aes(F1, F2,group=G,color=X1),alpha=0.1)+
#   coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3)

temp01<- hazdat$data[ , ,1]
tmp_01 <- as.vector(temp01)

lon=hazdat$lon
lat=hazdat$lat
lonlat <- as.matrix(expand.grid(lon,lat))

tmp_02<-paste(maxevent[,1],maxevent[,2])
mxs<-c()
for (mx in unique(tmp_02)){
maxevsum<-sum(maxevent$vecmeta[which(tmp_02==mx)])
mxs<-c(mxs,maxevsum)
}

mxv<-c()
for (mv in unique(tmp_02)){
  maxevacc<-max(maxevent$vecmeta2[which(tmp_02==mv)])
  mxv<-c(mxv,maxevacc)
}
tmp_05<-unique(tmp_02)
tmp_03<-paste(lonlat[,1],lonlat[,2])
tmp_04<-match(tmp_05,tmp_03)
tmp_01[tmp_04]<-mxs
tmp_01[-tmp_04]<-0
nlat<-dim(lat)
nlon=dim(lon)
tmp_array <- array(tmp_01, dim=c(nlon,nlat))

  grid <- expand.grid(lon=lon, lat=lat)
  mapmat=tmp_array
  #column 1634 corresponding to Dec 2015
  #This command compresses numbers larger than 6 to 6
  # plot.new()
  tmp_06<-cbind(mxs,tmp_05)
  tmp_07<-maxevent[match(tmp_05,tmp_02),]
  elmaxou<-tmp_07
  elmaxou$maxx=mxs
  elmaxou$mwg=mxv
  
  ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
    theme_bw(16)+
    geom_raster(data=elmaxou,aes(x=Var1,y=Var2,fill=maxx,group=cluster),alpha=.9,interpolate = F)+
    geom_polygon(fill = "transparent", color = "gray10", size = 1) +
    coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
  scale_fill_gradientn(colours = rgb.palette(100))
  


  par(mar=c(4,5,3,0))
  int=seq(0,max(mxs),length.out=100)
  rgb.palette=colorRampPalette(c("white", 
                                 "royalblue","orange","red","purple"),interpolate="linear",bias=1)
  int2<-c(0,1)
  rgb.palette2=colorRampPalette(c("white", 
                                  "yellow","pink","red","maroon"))
  
  mapmat= mapmat[,seq(length(mapmat[1,]),1)]
  
  lat=rev(lat)
  
  filled.contour(lon, lat, mapmat, col = rgb.palette(100), levels=int,
                 plot.axes={axis(1, cex.axis=1.5);
                   axis(2, cex.axis=1.5);map("world", add=TRUE);grid()},
                 key.title={par(cex.main=1);title(main="[mm]")},
                 key.axes={axis(4, cex.axis=1.5,lwd.ticks = 1)})
  
##############################Temporal lag between events###########

  totdates=list()
  totloc<-list()
  tempclus1<-list()
  for (haz in 1:2){
    stseason<-metaHaz[[haz]]
    var.int=stseason$rf.max
    mlev<-order(-var.int)
    plot(var.int[mlev])
    idmt<-stseason$ev[mlev]
    
    dates<-c()
    loca<-c()
    totemp<-c()
    for(idi in 1:length(stseason[,1])){
      idm=idmt[idi]
      maxevent<-metavHour[[haz]][which(metavHour[[haz]]$cluster==idm),]
      maxevent$idloc<-paste(maxevent$Var1,maxevent$Var2)
      idloco<-unique(maxevent$idloc)
      metamax<-stseason[which(stseason$ev==idm),]
      
      temp<-(unique(maxevent$Var3))
      tempev<-data.frame(temp, rep(idm, length(temp)))
      totemp<-rbind(totemp, tempev)
      maxevent$time <- as_datetime(c(maxevent$Var3*60*60),origin="1900-01-01")
      dates<-rbind(dates,c(maxevent$time[1],maxevent$time[length(maxevent$time)]))
      loca<-rbind(loca,metamax)
    }
    dates<-as.data.frame(dates)
    dates[,1]<-as_datetime(dates[,1],origin="1970-01-01")
    dates[,2]<-as_datetime(dates[,2],origin="1970-01-01")
    names(dates)=c("start","end")
    totdates<-c(totdates,list(dates))
    totloc<-c(totloc,list(loca))
    tempclus1<-c(tempclus1,list(totemp))
  }
  
  
  names(tempclus1)=c("Pr-Events","Wind_Events","Compound")
  tc<-na.omit(match(tempclus1[[1]]$temp,tempclus1[[2]]$temp))
  tcx<-na.omit(match(tempclus1[[2]]$temp,tempclus1[[1]]$temp))
  tc1<-tempclus1[[2]][tc,]
  tc2<-tempclus1[[1]][tcx,]
  tcu0<-data.frame(tc1,tc2)
  qqev<-c()
  for (z in unique(tc1$temp)){ 
    qev<-length(tc1$temp[which(tc1$temp==z)])
    qqev<-c(qqev,qev)
  }
plot(qqev)


qqev2<-c()
for (z in unique(tc2$temp)){ 
  qev<-length(tc2$temp[which(tc2$temp==z)])
  qqev2<-c(qqev2,qev)
}
plot(qqev2)

  tcx<-unique(tcx)
  tc<-unique(tc)
  tc1<-tempclus1[[2]][tc,]
  tc2<-tempclus1[[1]][tcx,]

  tc1<-tc1[order(tc1[,1]),]
  tc2<-tc2[order(tc2[,1]),]
  tcu<-data.frame(tc1,tc2,qqev,qqev2)
  names(tcu)=c("timev","rainev","time","windev","nbrain","nbwind")
  tcu$time <- as_datetime(c(tcu$time*60*60),origin="1900-01-01")
  
  plot(tcu$time, tcu$nbwind,type="o")
  points(tcu$time, tcu$nbrain,type="o",col=2)
  didi<-diff(tcu$time)
  ltot<-c()
  metatc=c()
  l=1
  for (t in 1:length(diff(tcu$time))){
    if(l==1){
      metatc<-rbind(metatc,tcu[t,])}
    if (didi[t]==1){l=l+1
    }
    if(didi[t]>1){ lcb=l
   ltot<-c(ltot,lcb)
   l=1 }
    if(t==length(diff(tcu$time))) ltot<-c(ltot,l)
    
  }
  plot(ltot[order(ltot)]) 
  hist(ltot,breaks=50)
  metatc<-data.frame(metatc,ltot)
  plot(totdates[[1]][,1])
  points(totdates[[3]][,1],col=2,pch=3)
  points(totdates[[2]][,1],col=3,pch=4) 
  teslag<-totdates[[2]]
  teslag2<-totloc[[2]]
  
  teslag2<-teslag2[order(teslag[,1]),]
  teslag<-teslag[order(teslag[,1]),]
metav2<-data.frame(teslag2,teslag)
  plot(metav2$start,metav2$rf.max)
plot(metav2$start,metav2$rf.sum, type="h")  
ggplot(metav2, aes(start,x.dur))+
  geom_point(aes(colour=wg.max,size=rf.max))+
  scale_colour_gradientn(colours = rgb.palette(100))

  geom_rect(aes(xmin = start, xmax = end,ymin = (ev - .5) , ymax = (ev +.5) ,fill=wg.max))


  lagtot<-c()
  lagsp<-c()
  sRt<-c()
  sWt<-c()
  ev1<-c()
  ev2<-c()
  for(lag in 1:(length(teslag[,1])-1)){
    laggy<-teslag[lag+1,1]-teslag[lag,1]
    lagtot<-c(lagtot,laggy)
    interes<-as.matrix(rbind(teslag2[lag+1,c(20,21)],teslag2[lag,c(20,21)]))
    latsp<-distGeo(interes[1,],interes[2,])/1000
    lagsp<-c(lagsp,latsp)
    evname1<-teslag2$ev[lag]
    ev1<-c(ev1,evname1)
    evname2<-teslag2$ev[lag+1]
    ev2<-c(ev2,evname2)
    sumR<-sum(teslag2$rf.sum[lag],teslag2$rf.sum[lag+1])
    sRt<-c(sRt,sumR)
    sumW<-max(teslag2$wg.max[lag],teslag2$wg.max[lag+1])
    sWt<-c(sWt,sumW)
  }


plot(lagtot+0.1,log="y") 
plot(lagsp+0.1,lagtot+0.1,log="xy")
min(lagtot)
lagh1<-data.frame(lagtot,lagsp,sRt,sWt,ev1,ev2)

names(lagh1)=c("templag","spacelag","r","w","ev1","ev2")

plot(lagh1$templag)

sameplace<-lagh1[which(lagh1$spacelag==0),]

boxplot(sameplace$templag)

ggplot(data=lagh1,aes(x=spacelag,y=templag,colour=w,size=r))+
  geom_jitter(alpha=.5,width=0, height=0)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_colour_gradientn(trans=scales::boxcox_trans(1.5),colors=rbPal(100),"Max wind gust [m/s]",breaks=c(10,20,30,40,50,Inf),limits=c(1,51))+
  scale_size(trans=scales::boxcox_trans(1.3),range=c(.4,10),"Max precipitation [mm]",breaks=c(100,1000,10000,20000,Inf),limits=c(0,25000),guide = guide_legend(override.aes = list(colour = alpha("grey",.8))))+
  scale_y_continuous(
                     breaks = c(1,3,6,12,24,48,72,148),limits = c(.49,125),"Temporal lag",
                     labels=c("0h","3h","6h","12h","1day","2days","3days","1week"))+
  scale_x_continuous(
                     breaks =c(1,10,100,1000,10000,100000),limits=c(.49,1000),"Spatial lag [km2]",labels=c("0","10","100","1,000","10,000","100,000"))





pb <- progress_bar$new(total = count[3])
el<-c()
stepb<-c()
yipi<-list()
metaevent<-list()
stevent<-list()
spoltot<-c()
spoltoti<-c()
fams<-c()
aeraf<-c()

for (i in 1:count[3]){ 
  pb$tick()
  Sys.sleep(1 / 100)
  temp11 <- newformat$data[ , ,i] #Level is the third dimension and time the fourth.
  temp12<-expand.grid(tp=temp11)
  elq<-quantile(temp12,0.98,na.rm=T)
  el<-c(el,elq)
  tempon<-temp11
  tempon[temp11<th1] <- 0
  temp11[temp11<th1] <- 0
  temp11[temp11>=th1] <- 1
  # newformat$lat <- rev(newformat$lat)
  # newformat$lon <- ifelse(newformat$lon > 180, -360 + newformat$lon, newformat$lon)
  
  # newformat$lon <- rev(newformat$lon)
  # image(newformat$lon,newformat$lat,temp11,col=rev(brewer.pal(11,"RdBu")))
  lon=newformat$lon
  lat=newformat$lat
  quantile(temp11,na.rm=T)
  grid <- expand.grid(lon=lon, lat=lat)

  
  #Install maps package if not done before
  #install.packages("maps")

  Lat= seq(-87.5, 87.5, length=36)
  Lon=seq(2.5, 357.5, length=72)
  mapmat=tempon
  #column 1634 corresponding to Dec 2015
  #This command compresses numbers larger than 6 to 6
  # plot.new()
  par(mar=c(4,5,3,0))
  int=seq(0,5.4,length.out=10)
  rgb.palette=colorRampPalette(c("white", 
                                 "yellow","pink","red","maroon"),interpolate="spline")
  int2<-c(0,1)
  rgb.palette2=colorRampPalette(c("white", 
                                 "yellow","pink","red","maroon"))
  mapmat= mapmat[,seq(length(mapmat[1,]),1)]
  # 

  # # lon=rev(lon)
  # lat=rev(lat)
  # filled.contour(lon, lat, mapmat, col = rgb.palette(10), levels=int,
  #                plot.axes={axis(1, cex.axis=1.5);
  #                  axis(2, cex.axis=1.5);map("world", add=TRUE);grid()},
  #                key.title=title(main="[mm]"),
  #                key.axes={axis(4, cex.axis=1.5)})
  # lat=rev(lat)

  lonlat <- as.matrix(expand.grid(lon,lat))
  tmp_vec <- as.vector(temp11)
  tmp_met <- as.vector(tempon)
  yihou <- data.frame(cbind(lonlat,tmp_vec))
  metb<- data.frame(cbind(lonlat,tmp_met))
  names(yihou) <- c("lon","lat",paste(name.var[1]))
  head(na.omit(yihou), 10)
  
  names(metb) <- c("lon","lat",paste(name.var[1]))
  head(na.omit(metb), 10)
  
  yihou<-yihou[which(yihou$tp==1),]
  metb<- metb[which(metb$tp>0),]
  yip<-yihou[,c(1,2)]
  minpt<-35
  if(length(yihou[,1])>minpt){

  ep<-2
  res <- optics(yip,eps=ep, minPts = 8)
  ### get order
  res$order
  ### plot produces a reachability plot

  ### plot the order of points in the reachability plot
  plot(yip, col = "grey")
  polygon(yip[res$order,])
  # 
  # res <- extractDBSCAN(res, eps_cl = .065)


  ### re-cut at a higher eps threshold
  epcl=sqrt(.50^2*2)
  res <- extractDBSCAN(res, eps_cl = epcl)


  plot(res)
  # hullplot(yip,res)
  plot(yip,col=res$cluster)
  # hullplot(yip[,c(1,2)], res)
  yip<-cbind(yip,res$cluster)
  metbg<-cbind(metb,yip$`res$cluster`)

  lc<-c()
  for(uk in unique(res$cluster))lc<-c(lc,length(which(res$cluster==uk)))
  names(lc)<-unique(res$cluster)
  op<-lc[which(lc<3 & as.numeric(names(lc)!=0))]
  if(length(op)>0){
    if (length(op)==1){
  yip<-yip[-which(yip[,3]==as.numeric(names(op))),]
  metbg<-metbg[-which(metbg[,4]==as.numeric(names(op))),] }else{
    for(opp in as.numeric(names(op))) 
    yip<-yip[-which(yip[,3]==opp),]
    metbg<-metbg[-which(metbg[,4]==opp),]
    
  }

  }
  yip<-yip[-which(yip[,3]==0),]
  metbg<-metbg[-which(metbg[,4]==0),]
  
  yip$`res$cluster`=yip$`res$cluster`+10*i
  metbg[,4]<-  yip$`res$cluster`
  if(length(metbg[,1])>0){
  metevent<-aggregate(metbg[,3] ,
                      by = list(ev = metbg[,4]),
                      FUN = function(x) c(sum = sum(na.omit(x)),max =max(x),surf=30.25*length(x),mean=mean(x),sd=sd(x)))
  metevent <- do.call(data.frame, metevent)

  }else{metevent<-NA}
  st<-unique(yip[,3])
cent<-c()
  for (u in st){ 
    ab<-convHull(yip[which(yip[,3]==u),c(1,2)])
    ab1<-ab@polygons
    CRS.new <- CRS("+init=epsg:2056")
    proj4string(ab1) <- CRS("+init=epsg:4326")
    dab1<- spTransform(ab1, CRS.new)
  
    spol<-c(u,dab1@polygons[[1]]@area/1000000)
    spoltot<-rbind(spoltot,spol)
    if(ab1@polygons[[1]]@area==0){
      cent<-rbind(cent,c(NA,NA))
    }else{
      cent<-rbind(cent,centroid(ab1))}
    # print(u)
    
    if(length(stepb[,1])>minpt){
    for (v in unique(stepb[,3])){ 
      ac<-convHull(stepb[which(stepb[,3]==v),c(1,2)])
      ac1<-ac@polygons
      spola<-c(v,ac1@polygons[[1]]@area)
    
    
      # aa<-over(ab1,ac1)
      aa=0
      try(aa<-predict(ac,yip[which(yip[,3]==u),c(1,2)]),silent=T)
      if(sum(aa)>0){
        yip[which(yip[,3]==u),3]=v
        metbg[which(metbg[,4]==u),3]=v
        fams<-rbind(fams,c(u,v))
        surfc<-convHull(rbind(yip[which(yip[,3]==u),c(1,2)],stepb[which(stepb[,3]==v),c(1,2)]))@polygons
        CRS.new <- CRS("+init=epsg:2056")
        proj4string(surfc) <- CRS("+init=epsg:4326")
        ssurfc<-spTransform(surfc, CRS.new)
        spoli<-c(u,v,ssurfc@polygons[[1]]@area)
        spoltoti<-rbind(spoltoti,spoli)
        # print(v)
      }
    }
    }
  
  } 
metevent<-cbind(metevent,cent)
if (length(metevent[,1])>1) names(metevent)[c(7,8)]=c("centlon","centlat")
  yipi<-c(yipi,list(yip))
  metaevent<-c(metaevent,list(metevent))
  stevent<-c(stevent,list(metbg))
  stepb<-yip
}else{stepb=data.frame(0,0)}



}

save(yipi,file="rainevents_2018_ERA5.Rdata")
save(metaevent,file="metaevents_rain_2018.Rdata")
save(fams,file="Temporal_clust_R2018.Rdata")

###===========================================================================================================================

load("rainevents_2018_ERA5.Rdata")
load("metaevents_rain_2018.Rdata")

unique(yipi[[1]][,3])

# eventX<-c()
# plot(stevent[[1]][,c(1,2)],col=stevent[[1]][,3])
# for(ld in 1:length(yipi)){
#   eventX<-rbind(eventX,stevent[[ld]][which(stevent[[ld]][,4]==11),c(1,2,3)])
# points(yipi[[ld]][which(yipi[[ld]][,3]==11),c(1,2)],col=alpha(ld,0.1))
#   print(length(which(yipi[[ld]][,3]==11)))
# }
# plot(eventX[,c(1,2)])


# nlon<- length(eventX$lon)
# nlat<-length(eventX$lat)
# eventXt <- matrix(eventX$tp, nrow=nlon, ncol=nlat)
# grid <- expand.grid(lon=eventX$lon, lat=eventX$lat)
# levelplot(eventXt ~ eventX$lon * eventX$lat, data=grid, pretty=T, 
#           col.regions=(rev(brewer.pal(10,"RdBu"))), main="MAT (C)") 


totev<-c()

upi<-c()
allmet<-c()
for (iti in 1:length(yipi)){try(totev<-c(totev,(unique(yipi[[iti]][,3]))),silent=T)
  if(length(metaevent[[iti]])>1){
  names(metaevent[[iti]])[c(7,8)]=c("centlon","centlat")
  allmet<-rbind(allmet,metaevent[[iti]])
  }
  }

length(totev)
plot(allmet$x.sd,allmet$x.sum)
tempclus<-match(fams[,1],allmet[,1])

allmet[,9]<-allmet[,1]

allmet[tempclus,9]<-fams[,2]

nmet<-na.omit(allmet)

###STOPPED HERE


for (tc in 2:length(nmet[,1])){
  if((nmet[tc-1,5]-nmet[tc,5])==0){
    nmet[tc,c(2,4)]<-sum(allmet[tc,c(2,4)],allmet[tc-1,c(2,4)])
    nmet[tc,3]<-max(allmet[tc-1,3],allmet[tc-1,3])
    nmet[tc-1,c(2,3,4)]<-NA
  }
}
nmet[,1]<-nmet[,5]
fmet<-na.omit(nmet)

fmet
unique(fmet[,1])-unique(totev)
outf<-c()
for(kx in unique(fmet[,1]))
{
  temet<-fmet[which(fmet[,1]==kx),]
  id=kx
  vol<-sum(temet[,2])
  peak<-max(temet[,3])
  surf<-sum(temet[,4])
  dur<-length(temet[,1])
outp<-c(id,vol,peak,surf,dur)
outf<-rbind(outf,outp)
  
}
outf<-data.frame(outf)
plot(outf[,4],outf[,5])
match(totev,fmet[,1])
names(outf)=c("id","vol","peak","surf","dur")

plot(outf[,c(2,3,4,5)])
plot(outf$peak,outf$vol/outf$surf,log="xy",pch=19)
plot(outf$dur,outf$surf,cex=5*outf$vol/max(outf$vol),pch=19)
#   for (un in unique(yipi[[iti]][,3])){
#   area<-cbind(un,30.25*length(which(yipi[[iti]][,3]==un)))
#   upi<-rbind(upi,area)}
# }
length(yipi)
length(totev)
fams
plot(allmet$x.sum,allmet$x.surf)
totev
for (kk in 1:length(fams[,1])){
spoltot[which(spoltot[,1]==fams[kk,1]),1]=fams[kk,2]}

plot(spoltot[,1],log="xy")
points(upi[,1],col=2)
length(upi[,1])
length(spoltot[,1])
spoltot[,1]
length(spoltoti[,1])
ak47<-match(spoltot[,1],spoltoti[,1])
aa<-unique(spoltoti[,2])
wok<-spoltot[which(!is.na(ak47)),1]
length(wok)
wok



emerd<-match(upi[,1],spoltot[,1])
twist<-spoltot[emerd,1]
iya<-cbind(twist,upi[,1])
match(wok,twist)
spoltoti
totev<-unique(totev)
length(fams[,1])
aeraf
fams
spoltot
spoltoti
totfr<-as.data.frame(table(upi[,1]))
bibis<-aggregate(upi[,2] ,
                 by = list(ev = upi[,1]),
                 FUN = function(x) c(sum = sum(na.omit(x)),n =length(x)))

bibis <- do.call(data.frame, bibis)

bobos<-aggregate(spoltot[,2] ,
                 by = list(ev = spoltot[,1]),
                 FUN = function(x) c(sum = sum(na.omit(x)),n =length(x)))

bobos <- do.call(data.frame, bobos)

plot(bibis$x.n,bibis$x.sum)
hihi<-jitter(bibis$x.n, factor = 1, amount = NULL)
plot(bobos$x.sum,bobos$x.n)
totfr[,1]<-as.numeric(totfr[,1])
totfr<-totfr[order(-totfr$Freq),]
plot(totfr, type="h")
hist(totfr$Freq)
totun<-unique(totev)
library(animation)
saveGIF(for (i in 1:365){
temp11 <- newformat$data[ , ,i] #Level is the third dimension and time the fourth.
temp12<-expand.grid(tp=temp11)
elq<-quantile(temp12,0.98,na.rm=T)
el<-c(el,elq)
temp11[temp11<20] <- NA
# newformat$lat <- rev(newformat$lat)
newformat$lon <- ifelse(newformat$lon > 180, -360 + newformat$lon, newformat$lon)

# newformat$lon <- rev(newformat$lon)
# image(newformat$lon,newformat$lat,temp11,col=rev(brewer.pal(11,"RdBu")))
lon=newformat$lon
lat=newformat$lat
# quantile(temp11,na.rm=T)
# grid <- expand.grid(lon=lon, lat=lat)
# 
# cutpts <- c(0,10,20,30,40,50,60,70,80,90,100,110,120)
# levelplot(newformat$data ~ lon * lat, data=grid, at=cutpts, cuts=13, pretty=T, 
#           col.regions=(rev(brewer.pal(10,"RdBu"))))
# 

#Install maps package if not done before
#install.packages("maps")
library(maps)
Lat= seq(-87.5, 87.5, length=36)
Lon=seq(2.5, 357.5, length=72)
mapmat=temp11
#column 1634 corresponding to Dec 2015
#This command compresses numbers larger than 6 to 6
# plot.new()
par(mar=c(4,5,3,0))
int=seq(0,120,length.out=10)
rgb.palette=colorRampPalette(c("green", 
                                "yellow","pink","red","maroon"),interpolate="spline")
mapmat= mapmat[,seq(length(mapmat[1,]),1)]

filled.contour(lon, lat, mapmat, color.palette=rgb.palette, levels=int,
               plot.axes={axis(1, cex.axis=1.5);
                 axis(2, cex.axis=1.5);map("world", add=TRUE);grid()}, main=time2[i],
               key.title=title(main="[mm]"), 
               key.axes={axis(4, cex.axis=1.5)}) 

}
, movie.name = paste0(getwd(),"/champion.gif"), img.name = "cc", convert = "magick",cmd.fun, clean = TRUE, extra.opts = "") 
lonlat <- as.matrix(expand.grid(lon,lat))
dim(lonlat)  
# vector of `tmp` values
tmp_vec <- as.vector(temp11)
length(tmp_vec)



# create dataframe and add names
tmp_df01 <- data.frame(cbind(grid,tmp_vec))
names(tmp_df01) <- c("lon","lat",paste(name.var[1]))
head(na.omit(tmp_df01), 10)


coordinates(tmp_df01)= ~lon+lat
EXTEND<-bbox(tmp_df01)
register_google(key = "AIzaSyCaHpsh-Y3zXBwtc1tVfDt0zPwc-UUq4Hk")
mamap <- ggmap(get_map(location = EXTEND, zoom=3 ,maptype = "satellite"))

mamap
tmp_df01 <- data.frame(cbind(grid,tmp_vec))
names(tmp_df01) <- c("lon","lat",paste(name.var[1]))
mamap+ geom_sf(data = tmp_df01, aes(x=lon, y = lat, fill=tp),alpha=.7,interpolate=F)+ coord_cartesian()+
  scale_fill_distiller(palette = "Spectral")


install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", 
                   "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata")) 

world1 <- sf:::st_as_sf(map('europe', plot = FALSE, fill = TRUE)) 

ggplot() + geom_sf(data = world1)


ggplot(data = world1) +
  geom_sf() +
  coord_sf(crs = "+proj=laea +lat_0=30 +lon_0=-95")+
geom_tile(data = tmp_df01, aes(x=lon, y = lat, fill=tp),alpha=1)

  scale_fill_distiller(palette = "Spectral")


ggplot(data = world1) +
  geom_sf() +
  coord_sf(crs = st_crs(3035))




#YAAAAY
#NOW all time only on heathrow

long    = ncvar_get(nc,name.lon)
latt    = ncvar_get(nc,name.lat)
# Initialize start and count to read one timestep of the variable.
start <- rep(1,tdims) # begin with start=(1,1,1,...,1)
start[1] <- which(round(long,1)==-0.5) 
start[2] <- which(round(latt,1)==51.4) 
start[3] <-which(time2=="01/01/70") # change to start=(1,1,1,...,i) to read timestep i
count <- tsize # begin w/count=(nx,ny,nz,...,nt), reads entire var
count[1] <- 1 # change to count=(nx,ny,nz,...,1) to read 1 tstep
count[2] <-1
count[3] <- count[3]-start[3]+1
data3 <- ncvar_get( nc, t$name, start=start, count=count )


newformat<-list()


# tmp_array <- ncvar_get(nc,t, start = start, count= count)

newformat$data   = data3
newformat$lon    = ncvar_get(nc,name.lon,start=start[1],count=count[1])
newformat$lat    = ncvar_get(nc,name.lat,start=start[2],count=count[2])
dlname <- ncatt_get(nc,name.var,"long_name")
dunits <- ncatt_get(nc,name.var,"units")
fillvalue <- ncatt_get(nc,name.var,"_FillValue")
dim(tmp_array)
time <- ncvar_get(nc,"time")
time
tunits <- ncatt_get(nc,"time","units")
nt <- dim(time)

tustr <- strsplit(tunits$value, " ")
tdstr<-strsplit(unlist(tustr)[3], "  ")

ori<-as.Date(as.character(tdstr))
timi <- strsplit(unlist(tdstr),"-")
tmonth<-as.integer(unlist(timi)[2])
tday <- as.integer(unlist(timi)[3])
tyear <- as.integer(unlist(timi)[1])
tcrap<-strsplit(unlist(tustr)[4], ":")
thour <- as.integer(unlist(tcrap)[1])
time=as.vector(time)
time<-time
time2<-chron(time,origin=c(tmonth, tday, tyear))
time3<-time2[start[3]:length(time2)]

#this is a time serie for the London Heathrow pixel
timeserie<-data.frame(time3,newformat$data)

#Now compute the average temperature of the past week everyday

avgtemp<-c()
for (k in 8:length(timeserie$newformat.data)){
 avgweek<-mean(timeserie$newformat.data[(k-7):k]) 

 avgtemp<-c(avgtemp,avgweek)
}
timek<-timeserie$time3[8:length(timeserie$newformat.data)]
avgw<-data.frame(timek,avgtemp)

plot(timeserie,type="p",col="blue")  
lines(avgw,type="l",col="red")
plot(timeserie$newformat.data[-(1:7)],avgw$avgtemp)

#Same for mean last month temperature

avgmonth<-c()
for (k in 31:length(timeserie$newformat.data)){
  avgm<-mean(timeserie$newformat.data[(k-30):k]) 
  
  avgmonth<-c(avgmonth,avgm)
}
timek<-timeserie$time3[31:length(timeserie$newformat.data)]
avgmo<-data.frame(timek,avgmonth)

plot(timeserie,type="p",col="blue")  
lines(avgmo,type="l",col="red")
plot(timeserie$newformat.data[-(1:30)],avgmo$avgmonth)

names(timeserie)=c("date","temperature")
names(avgw)=c("date","mean_past_week_temp")
mm<-match(timeserie$date,avgw$date)
mm<-which(!is.na(mm)) 

Temp_Heathrouw<-timeserie[mm,]

Temp_Heathrouw<- data.frame(Temp_Heathrouw,avgw$mean_past_week_temp)

names(avgmo)=c("date","mean_past_month_temp")
mm<-match(Temp_Heathrouw$date,avgmo$date)
mm<-which(!is.na(mm)) 

Temp_Heathrow2<-Temp_Heathrouw[mm,]
Temp_Heathrow2<- data.frame(Temp_Heathrow2,avgmo$mean_past_month_temp)
save(Temp_Heathrow2, file= "Heathrowtemperature_2000-2018_N.Rdata")



