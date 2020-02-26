
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

newrun=FALSE
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


tustr <- strsplit(tunits$value, " ")

library(chron)
library(lattice)
library(RColorBrewer)
library(dbscan)
library(dismo)
library(maps)
library(progress)

rbPal <- colorRampPalette(c('darkgreen',"gold","darkorange",'red'))
names(hazmat)=c("Pr","Wg")
pqt<-.99
#Set up a quantile for each grid cell
thr<-hazmat$Pr$data[,,1]
for (i in (1:33)){
  for(j in 1:45){
crappy<-hazmat$Pr$data[i,j,]
thr[i,j]<-quantile(crappy[which(crappy>0)],pqt,na.rm=T)
  }
}

thw<-hazmat$Wg$data[,,1]
for (i in (1:33)){
  for(j in 1:45){
    crappy<-hazmat$Wg$data[i,j,]
    thw[i,j]<-quantile(crappy[which(crappy>0)],pqt,na.rm=T)
  }
}
thrr<-as.vector(thr) 
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
sum(hazmat[[1]]$data)/sum(hazmat[[1]]$data[which(hazmat[[1]]$data>th1)])
quantile(hazmat$Pr[[1]][16,4,],.99)
metaHaz<-list()
metavHour<-list()
metavDaz<-list()
for(hazard in 1:3){
if(hazard==1){hazdat=hazmat$Pr;th=thr}
if(hazard==2){hazdat=hazmat$Wg;th=thw}
if(hazard==3){hazdat1=hazmat$Pr;hazdat2=hazmat$Wg}

lon=hazmat$Pr$lon
lat=hazmat$Pr$lat

lonlatime <- expand.grid(lon, lat,time)

formeta<-hazmat$Pr$data
formeta2<-hazmat$Wg$data
if (hazard == 3){
  bolilos<-hazdat1$data
  vecthouf1<-as.vector(rep(thr,length(bolilos[1,1,])))
  vecthouf2<-as.vector(rep(thw,length(bolilos[1,1,])))

  bolilos[which(hazdat1$data<vecthouf1 | hazdat2$data<vecthouf2)] <- NA
  bolilos[which(hazdat1$data>=vecthouf1 & hazdat2$data>=vecthouf2)] <- 1  
  formeta[which(hazdat1$data<vecthouf1 | hazdat2$data<vecthouf2)] <- NA
  formeta2[which(hazdat1$data<vecthouf1 | hazdat2$data<vecthouf2)] <- NA
  
  vectouf<- as.vector(bolilos)
  length(na.omit(vectouf))
  vrac<-as.vector(bolilos)
  length(na.omit(vrac))
  vecmeta<-as.vector(formeta)
  vecmeta2<-as.vector(formeta2)
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

lonlatemp <- data.frame(cbind(lonlatime,vectouf))
metav<-data.frame(cbind(lonlatime,vecmeta,vecmeta2))

# 
# x=abs(rnorm(100*100,50,25))
# x=matrix(x,nrow=100)
# x1=melt(th)


if(hazard==1)metav<-metav[which(!is.na(metav[,4])),]
if(hazard==2)metav<-metav[which(!is.na(metav[,5])),]
if(hazard==3)metav<-metav[which(!is.na(metav[,5])),]
print(length(metav$Var1)) 
head(metav)
lonlatemp2<-lonlatemp[which(lonlatemp[,4]==1),]
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

# walabibou<-kNNdist(sampspd,k=10,all=F)
# min(walabibou)
# walabibou<-walabibou[which(walabibou<10)]
# walabibord<-seq(1:length(walabibou))
# walabibix<-data.frame(as.vector(walabibou[order(walabibou)]))
# walabibix$order=walabibord

# fatos<-walabibix[which(walabibix$order>1270000),]
# light<-seq(1,1270000,by=1000)
# 
# walabilight<-rbind(walabibix[which(!is.na(match(walabibix$order,light))),],fatos)
# xtick=seq(0,1.5e6,by=2e5)
# par(mar=c(5,5,1,1))
# plot(walabilight$order,walabilight$as.vector.walabibou.,ylim=c(0,10),ylab="10-NN distance",xlab="Points sorted by distance",cex.axis=1.5,cex.lab=1.8,xaxt="n")
# axis(side=1, at=xtick, labels = T,cex.axis=1.5)
# text(2e5,2.5,cex=1.8, labels=expression(paste( epsilon, " = 2.5")),pos=3,col="red")
# abline(h=2.5, col=2,lwd=3) 

# walabibou[order(walabibou)]
#
sqrt(0.5^2+0.5^2+1^2)
epcl=2.5
epl<-sqrt(1.4^2+1.4^2+2)


# 
# reservo <- optics(sampspd,eps=3, minPts = 15)
# res=extractDBSCAN(reservo, eps_cl = epcl)
# plot(res)
# unique(res$cluster)
# bip<-extractXi(object=res, xi=.2)
# reach=as.reachability(bip)
# dend<-as.dendrogram(reach)
# 

if(hazard==1)weightc=metav$vecmeta
if(hazard==2)weightc=metav$vecmeta2
if(hazard==3)weightc=rep(1,length(sampspd[,1]))
rpip<-dbscan(sampspd, eps=epcl, minPts = 10,weights = weightc)

length(unique(rpip$cluster))

# 
# nn<-frNN(sampspd,eps=epcl,sort=TRUE)

# Number of neighbors 
# hist(sapply(adjacencylist(nn), length), xlab = "k", main="Number of Neighbors", sub = paste("Neighborhood size eps =", nn$eps))
# 
# # Explore neighbors of point i = 10 
# i <- 10
# nn$id[[i]]
# nn$dist[[i]] 
# plot(sampspd, col = ifelse(1:nrow(iris) %in% nn$id[[i]], "red", "black"))
# # get an adjacency list
# 
# head(adjacencylist(nn))
# 
# #plot the fixed radius neighbors (and then reduced to a radius of .3) 
# plot(nn, sampspd) 
# plot(frNN(nn, .3), x)

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
length(metav$Var1)
length(spdata$Var1)
# for(m in unique(spdata[,4])){
metav$time<-as_datetime(c(metav$Var3*60*60),origin="1900-01-01")
metav$month=month(metav$time)
  event<-metav
  charloc<-paste(event[,1],event[,2])
  event$cloc=charloc
  print(length(event$Var1))

  testev<-aggregate(list(rf= event[,4],wg=event[,5]),
                    by = list(ev = event[,6],loc=event[,9]),
                    FUN = function(x) c(sum = sum(na.omit(x)),max=max(na.omit(x))))
  testev <- do.call(data.frame, testev)
  
  thl<-aggregate(list(rf= event[,4],wg=event[,5]),
                    by = list(ev = event[,6]),
                    FUN = function(x) c(l= length(x)))
  thl<- do.call(data.frame, thl)
  small<-thl$ev[which(thl$rf<15)]
  bip<-which(!is.na(match(testev$ev,small)))
  testev<-testev[-bip,]
  bop<-which(!is.na(match(event$cluster,small)))
  event<-event[-bop,]
  length(unique(event$cluster))
  length(unique(testev$ev))
  
  metamax<-aggregate(list(rf= testev[,3],wg=testev[,6]) ,
                     by = list(ev = testev[,1]),
                     FUN = function(x) c(sum = sum(na.omit(x)),max=max(na.omit(x)),mean=mean(x),sd=sd(x),surf=length(x)))
  
  metamax<- do.call(data.frame, metamax)
  
  metave<-aggregate(list(rf= event[,4],wg=event[,5]) ,
                      by = list(ev = event[,6]),
                      FUN = function(x) c(surf=length(x)))
  
  metave <- do.call(data.frame, metave)
  length(unique(metave$ev))
  
  tempcom<-aggregate(event[,7] ,
                     by = list(ev = event[,6]),
                     FUN = function(x) c(dur=length(unique(x)),month=month(unique(x))[1],year=year(unique(x))[1]))

  tempic<- do.call(data.frame, tempcom)
  
  metave<-cbind(metamax,metave[,c(2,3)],tempic[,c(2,3,4)])
maxR<-c()
maxW<-c()
evbk<-c()
startev<-vector(length=length(metave$ev))
endev<-vector(length=length(metave$ev))
  for (eve in 1:length(metave$ev))
  { 
    eves<-metave$ev[eve]
    eventr<-event$time[which(event$cluster==eves)]
    evint<-testev[which(testev$ev==eves),]
    met<-metamax[which(metamax$ev==eves),]
    startev[eve]<-min(eventr)
    endev[eve]<-max(eventr)
    elW<-as.character(evint$loc[which(evint$wg.max==met$wg.max)])[1]
    elR<-as.character(evint$loc[which(evint$rf.sum==met$rf.max)])[1]
    elW<-as.character(evint$loc[which(evint$wg.max==met$wg.max)])[1]
    maxR<-c(maxR,elR)
    maxW<-c(maxW,elW)
    evbk<-c(evbk,eves)
    
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
  metave$startev<-startev
  metave$endev<-endev
  metaveN<-metave
  hist(metaveN$Dmax)
  min(metaveN$rf)
  plot(metaveN$rf.surf,metaveN$x.dur)
  length(unique(metaveN$ev))


  
plot(metaveN$rf.max,metaveN$wg.max)
metaHaz<-c(metaHaz,list(metaveN))
metavHour<-c(metavHour,list(event))
metavDaz<-c(metavDaz,list(testev))
}


#Retain only events that affected lands

metaHax<-list()
metavHax<-list()
metavDax<-list()
for (hx in 1:3){
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
rm(metaHaz,metavHour,metavDaz,lonlatemp,vecmeta,vecmeta2,vecthouf,vecthouf1,vecthouf2,vectouf,vrac,hazdat,hazdat2,newformat,bolilos,formeta,formeta2)
gc()
############################################
#Then select of pairs of events that overlap spatially

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

metaHax[[1]]$startev<-as_datetime(c(metaHax[[1]]$startev),origin="1970-01-01")
metaHax[[1]]$endev<-as_datetime(c(metaHax[[1]]$endev),origin="1970-01-01")
metaHax[[2]]$startev<-as_datetime(c(metaHax[[2]]$startev),origin="1970-01-01")
metaHax[[2]]$endev<-as_datetime(c(metaHax[[2]]$endev),origin="1970-01-01")
lonlatime$time<-as_datetime(c(lonlatime$Var3*60*60),origin="1900-01-01")
arain<-as.vector(hazmat[[1]]$data)
head(lonlatime$time)
oc<-c()
nvo<-c()
nwo<-c()
nvo<-vector(mode="list",length=length(metaHax[[1]]$ev))
library(fastmatch)
ptm <- proc.time()
for (c in 1: length(metaHax[[1]]$ev)){ 

  cl=metaHax[[1]]$ev[c]
  clev<-metavDax[[1]][which(metavDax[[1]]$ev==cl),]
  # dd<-unique(metavHax[[1]]$time[which(metavHax[[1]]$cluster==cl)])
  dd<-seq(metaHax[[1]]$startev[c],metaHax[[1]]$endev[c],by="hour")
  
  clex<-do.call(rbind, replicate(length(dd),clev[,c(1,2)], simplify=FALSE))
  cley<-do.call(rbind, replicate(length(clev[,1]),as.data.frame(dd),simplify=F))
  clex<-clex[order(clex$loc),]
  ouch<-cbind(clex,cley)

  nvo[[c]]<-ouch 
  print(round((c/length(metaHax[[1]]$ev)*100),1))
  # nvo<-rbind(nvo,ouch)

}
proc.time()-ptm
nvo2<-ldply(nvo, data.frame)
lonlatime$loc<-paste(lonlatime$Var1,lonlatime$Var2)
lonlatime<-lonlatime[,-c(1,2)]
yoboys<-which(!is.na(fmatch(lonlatime$time,nvo2$dd)))
ole<-data.frame(lonlatime[yoboys,],arain[yoboys])

ole$loc=paste(ole$Var1,ole$Var2)
names(ole)[3] = "dd"
olel<-inner_join(ole,nvo2, by=c("loc","dd"))

rainacc<-aggregate(list(len=olel[,4]),
                by = list(ev=olel[,5],loc=olel[,2]),
                FUN = function(x) c(length=length(x),sum=sum(x),max=max(x)))
rainacc <- do.call(data.frame, rainacc)

rainacev<-aggregate(list(len=rainacc[,4]),
                   by = list(ev=rainacc[,1]),
                   FUN = function(x) c(length=length(x),mean=mean(x),max=max(x)))
rainacev <- do.call(data.frame, rainacev)
#last join
metaHax[[1]]<-data.frame(metaHax[[1]],rainacev[,c(3,4)])
plot(metaHax[[1]]$rf.max,metaHax[[1]]$len.max) 
abline(a=0,b=1)

nwo<-vector(mode="list",length=length(metaHax[[2]]$ev))
for (c in 1: length(metaHax[[2]]$ev)){ 
  cl=metaHax[[2]]$ev[c]
  clev<-metavDax[[2]][which(metavDax[[2]]$ev==cl),]
  # dd<-unique(metavHax[[2]]$time[which(metavHax[[2]]$cluster==cl)])
  dd<-seq(metaHax[[2]]$startev[c],metaHax[[2]]$endev[c],by="hour")
  
  clex<-do.call(rbind, replicate(length(dd),clev[,c(1,2)], simplify=FALSE))
  cley<-do.call(rbind, replicate(length(clev[,1]),as.data.frame(dd),simplify=F))
  clex<-clex[order(clex$loc),]
  ouch<-cbind(clex,cley)
  # nwo<-rbind(nwo,ouch)
  nwo[[c]]<-ouch
  print(round((c/length(metaHax[[2]]$ev)*100),1))
}
nwo2<-ldply(nwo, data.frame)
rm(nvo,nwo,ole1,lonlatime)



 #these data frame contains spatial id and temporal id over the duration of an event of each cell involved in that event. (e.g., if the event last 10 hours, each cell involved will be repeated 10 times even if the the cell is only impacted during 1 hour)

# save(nvo,file="rainallclustersF.Rdata")
# save(nwo,file="windallclustersF.Rdata")
# 
# load(file="rainallclustersF.Rdata")
# load(file="windallclustersF.Rdata")

  # clev<-metavHax[[1]]
  # cclev<-metavHax[[2]]
  # clev$charloc<-paste(clev$Var1,clev$Var2)
  # cclev$charloc<-paste(metavHax[[2]]$Var1,metavHax[[2]]$Var2)
  
  sp10<-match_df(nvo2,nwo2,on=c("loc","dd"))
  putin<-nwo[which(nwo$ev==498),]
  demerd<-nvo[which(nvo$ev==905),]
  #Inner join in R:  Return only the rows in which the left table have matching keys in the right table
  #this one set a pre-filter on events which have temporal overlap
  sp03<-inner_join(nwo2,nvo2,by=c("loc","dd"))
   #aggregation by events and space, each row correpond to one cell involved in both rain and wind event
   
 length(rainacc$ev)
  tesp<-aggregate(list(len=sp03[,3]),
                    by = list(ev2=sp03[,1],ev1 = sp03[,4],loc=sp03[,2]),
                    FUN = function(x) c(length=length(unique(x))))
  tesp <- do.call(data.frame, tesp)
  
  sesp<-aggregate(list(len=sp03[,2]),
                  by = list(ev2=sp03[,1],ev1 = sp03[,4],dd=sp03[,3]),
                  FUN = function(x) c(length=length(x)))
  sesp <- do.call(data.frame, sesp)

    
  
  #extract more metadata about the events Basically I will do everything from here: Spatial value (for combined  duration of the event), spatiotemporal value (only during the compound duration and in compound space), temporal (outside the compound space inside the compound time))

#1) the spatial compound (all time inside the compound space)
#2) the temporal compound (all space inside the compound time)
#3) the spatiotemporal compound ( inside the compound space and time)


###############################RAIN###########################################
  
  metatest1<-metavDax[[1]]
  names(metatest1)[c(1)]=c("ev1")
  
  metatest<-rainacc
  names(metatest)[c(1)]=c("ev1")
  
#1) Spatial compound
  metatestt<-join(metatest,metatest1,by=c("ev1","loc"))
  metatestt<-metatestt[,-c(6,7,8)]
  SCompR<-inner_join(metatestt,tesp,by=c("ev1","loc"))
  

  

  spSCR<-aggregate(list(rf=SCompR[,4],wg=SCompR[,6]),
                  by = list(ev2=SCompR[,7],ev1 = SCompR[,1]),
                  FUN = function(x) c(length=length((x)),max=max(x)))
  spSCR <- do.call(data.frame, spSCR)
  
  # spSCR2<-aggregate(list(rf=spSCR[,6],wg=spSCR[,8]), 
  #                    by = list(ev2=spSCR[,1],ev1=spSCR[,2]),  
  #                    FUN = function(x) c(length=length(x),max=max(x))) 
  # spSCR2<- do.call(data.frame, spSCR2)
  SCompR$cev=paste(SCompR$ev1,SCompR$ev2)
  for(rg in 1:length(spSCR$ev2)){
    lesev<-paste(spSCR[rg,2],spSCR[rg,1])
    mr<-spSCR$rf.max[rg]
    spSCR$maxlocsR[rg]<-as.character(SCompR$loc[which(SCompR$cev==lesev & SCompR$rf.sum==mr)])
  }

#3) Spatiotemporal compound
  names(sp03)<-c("ev2","loc","dd","ev1")
  spkik<-sp03 %>% group_by(ev1,dd,loc) %>% mutate(id = row_number())
  metatest<-metavHax[[1]]
  names(metatest)[c(6,7,9)]=c("ev1","dd","loc")
  names(olel)[5]="ev1"
  metatestx<-full_join(metatest,olel,by=c("ev1","dd","loc"))
  sp20<-inner_join(metatestx,spkik,by=c("ev1","loc","dd"))
  spSTCR<-sp20[which(sp20$id==1),]
  
  spSTCR2<-aggregate(list(rf=spSTCR[,11],wg=spSTCR[,5]),
                  by = list(ev2=spSTCR[,12],ev1 = spSTCR[,6],loc=spSTCR[,9]),
                  FUN = function(x) c(length=length(x),max=max(x,na.rm=T),sum=sum(x,na.rm=T)))
  spSTCR2 <- do.call(data.frame, spSTCR2)
  
  
  spSTFR<-aggregate(list(rf=spSTCR2[,6],wg=spSTCR2[,8]), 
                   by = list(ev2=spSTCR2[,1],ev1=spSTCR2[,2]),  
                   FUN = function(x) c(length=length(x),max=max(x))) 
  spSTFR<- do.call(data.frame, spSTFR)

  spSTCR2$cev=paste(spSTCR2$ev1,spSTCR2$ev2)
  spSTFR$maxlocR="00 00"
  for(rg in 1:length(spSTFR$ev2)){
    lesev<-paste(spSTFR[rg,2],spSTFR[rg,1])
    mr<-spSTFR$rf.max[rg]
    rtime=as.character(spSTCR2$loc[which(spSTCR2$cev==lesev & spSTCR2$rf.sum==mr)])
    if(length(rtime)>0) spSTFR$maxlocstR[rg]<-rtime
  } 
#2) Temporal compound
  # metateste<-aggregate(list(rf=metatest[,4],wg=metatest[,5]),
  #                      by = list(ev1 = metatest[,6],dd=metatest[,7]),
  #                      FUN = function(x) c(length=length((x)),max=max(x),sum=sum(x)))
  # metateste <- do.call(data.frame, metateste)

  TCompR<-inner_join(metatest,sesp,by=c("ev1","dd"))
  
  spTCR<-aggregate(list(rf=TCompR[,4],wg=TCompR[,5]),
                  by = list(ev2=TCompR[,10],ev1 = TCompR[,6],loc=TCompR[,9]),
                  FUN = function(x) c(length=length((x)),max=max(x),sum=sum(x)))
  spTCR <- do.call(data.frame, spTCR)
  
  spTCR2<-aggregate(list(rf=spTCR[,6],wg=spTCR[,8]),
                  by = list(ev2=spTCR[,1],ev1 = spTCR[,2]),
                  FUN = function(x) c(length=length((x)),max=max(x)))
  spTCR2 <- do.call(data.frame, spTCR2)

  spTCR$cev=paste(spTCR$ev1,spTCR$ev2)
  for(rg in 1:length(spTCR2$ev2)){
    lesev<-paste(spTCR2[rg,2],spTCR2[rg,1])
    mr<-spTCR2$rf.max[rg]
    spTCR2$maxloctR[rg]<-as.character(spTCR$loc[which(spTCR$cev==lesev & spTCR$rf.sum==mr)])
  }

 ########################################WIND######################################## 
metatest<-metavHax[[2]]
names(metatest)[c(6,7,9)]=c("ev2","dd","loc")

#1) Spatial compound

SCompW<-inner_join(metatest,tesp,by=c("ev2","loc"))



spSCW<-aggregate(list(rf=SCompW[,4],wg=SCompW[,5]),
                 by = list(ev2=SCompW[,10],ev1 = SCompW[,6],loc=SCompW[,9]),
                 FUN = function(x) c(length=length((x)),max=max(x),sum=sum(x)))
spSCW <- do.call(data.frame, spSCW)

spSCW2<-aggregate(list(rf=spSCW[,6],wg=spSCW[,8]), 
                  by = list(ev2=spSCW[,1],ev1=spSCW[,2]),  
                  FUN = function(x) c(length=length(x),max=max(x))) 
spSCW2<- do.call(data.frame, spSCW2)

spSCW$cev=paste(spSCW$ev1,spSCW$ev2)
for(rg in 1:length(spSCW2$ev2)){
  lesev<-paste(spSCW2[rg,2],spSCW2[rg,1])
  mr<-spSCW2$wg.max[rg]
  spSCW2$maxlocsW[rg]<-as.character(spSCW$loc[which(spSCW$cev==lesev & spSCW$wg.max==mr)])
}
#3) Spatiotemporal compound

spkik<-sp03 %>% group_by(ev2,dd,loc) %>% mutate(id = row_number())
sp20<-inner_join(metatest,spkik,by=c("ev2","loc","dd"))
spSTCW<-sp20[which(sp20$id==1),]

spSTCW2<-aggregate(list(rf=spSTCW[,4],wg=spSTCW[,5]),
                   by = list(ev2=spSTCW[,6],ev1 = spSTCW[,10],loc=spSTCW[,9]),
                   FUN = function(x) c(length=length((x)),max=max(x),sum=sum(x)))
spSTCW2 <- do.call(data.frame, spSTCW2)


spSTFW<-aggregate(list(rf=spSTCW2[,6],wg=spSTCW2[,8]), 
                  by = list(ev2=spSTCW2[,1],ev1=spSTCW2[,2]),  
                  FUN = function(x) c(length=length(x),max=max(x))) 
spSTFW<- do.call(data.frame, spSTFW)

spSTCW2$cev=paste(spSTCW2$ev1,spSTCW2$ev2)
for(rg in 1:length(spSTFW$ev2)){
  lesev<-paste(spSTFW[rg,2],spSTFW[rg,1])
  mr<-spSTFW$wg.max[rg]
  spSTFW$maxlocstW[rg]<-as.character(spSTCW2$loc[which(spSTCW2$cev==lesev & spSTCW2$wg.max==mr)])
}
#2) Temporal compound

TCompW<-inner_join(metatest,sesp,by=c("ev2","dd"))

spTCW<-aggregate(list(rf=TCompW[,4],wg=TCompW[,5]),
                 by = list(ev2=TCompW[,10],ev1 = TCompW[,6],loc=TCompW[,9]),
                 FUN = function(x) c(length=length((x)),max=max(x),sum=sum(x)))
spTCW <- do.call(data.frame, spTCW)

spTCW2<-aggregate(list(rf=spTCW[,6],wg=spTCW[,8]),
                  by = list(ev2=spTCW[,1],ev1 = spTCW[,2]),
                  FUN = function(x) c(length=length((x)),max=max(x)))
spTCW2 <- do.call(data.frame, spTCW2)

spTCW$cev=paste(spTCW$ev1,spTCW$ev2)
for(rg in 1:length(spTCW2$ev2)){
  lesev<-paste(spTCW2[rg,2],spTCW2[rg,1])
  mr<-spTCW2$wg.max[rg]
  spTCW2$maxlocTW[rg]<-as.character(spTCW$loc[which(spTCW$cev==lesev & spTCW$wg.max==mr)])
}



####################################################################
  spSTFR$gr="r"
  spSTFW$gr="w"
  spSTF<-full_join(spSTFR,spSTFW,by=c("ev1","ev2"))

  plot(spSTF$wg.max.y,spSTF$rf.max.x)
  for (mm in 1:length(spSTF$ev2)){
  spSTF$rf.max.c[mm]=max(spSTF$rf.max.x[mm],na.rm=T)
  spSTF$wg.max.c[mm]=max(spSTF$wg.max.y[mm],na.rm=T)
  }
  
  plot(spSTF$wg.max.c,spSTF$rf.max.c)
  spSTF<-spSTF[,c(2,1,3,8,9,10,14,15,16,17)]
  spSCW2<-spSCW2[,-3]
  spTCW2<-spTCW2[,-3]
  spSCR2<-spSCR2[,-3]
  spTCR2<-spTCR2[,-3]
  
  # names(metavDax[[2]])[1]="ev2"
  # sp11<-inner_join(tesp,metavDax[[2]],by=c("ev2","loc"))
  # sp11<-sp11[,c(1,2,3,4,8)]
  # 
  # names(metavDax[[1]])[1]="ev1"
  # sp12<-inner_join(tesp,metavDax[[1]],by=c("ev1","loc"))
  # sp12<-sp12[,c(1:6)]
  # tespp<-data.frame(sp12,sp11[,5])
  # 
  # 
  # endgame<-aggregate(list(rf=tespp[,5],wg=tespp[,7]), 
  #                    by = list(ev2=tespp[,1],ev1=tespp[,2]),  
  #                    FUN = function(x) c(length=length(x),max=max(x))) 
  # endgame <- do.call(data.frame, endgame)
  

  #this part looks useless, just add duration of the overlap
  
  # dude<-c()
  # for (clb in unique(spTCR2$ev1)){
  #   clev<-metavDax[[1]][which(metavDax[[1]]$ev1==clb),]
  #   dd<-unique(metavHax[[1]]$time[which(metavHax[[1]]$cluster==clb)])
  #   evr<-rep(clb,length(dd))
  #   durax<-data.frame(evr,dd)
  #   dude<-rbind(dude,durax)
  # }
  # dude2<-c()
  # for (clb in unique(spTCW2$ev1)){
  #   clev<-metavDax[[2]][which(metavDax[[2]]$ev2==clb),]
  #   dd<-unique(metavHax[[2]]$time[which(metavHax[[2]]$cluster==clb)])
  #   evr<-rep(clb,length(dd))
  #   durax<-data.frame(evr,dd)
  #   dude2<-rbind(dude2,durax)
  # }
  # spdude<-inner_join(dude,dude2,by=c("dd"))
  # spdude$combin=paste(spdude$evr.x,spdude$evr.y)
  # 
  # length(unique(spdude$combin))
  # 
  # temsp<-aggregate(list(shit=spdude[,2]),
  #                 by = list(ev1=spdude[,1],ev2 = spdude[,3]),
  #                 FUN = function(x) c(length=length(unique(x))))
  # temsp <- do.call(data.frame, temsp)
  # 
  # endgame$combin=paste(endgame$ev1,endgame$ev2)
  # temsp$combin=paste(temsp$ev1,temsp$ev2)
  # endgame2<-na.omit(match(endgame$combin,temsp$combin))
  # endgame3<-temsp[endgame2,]
  # 
  
  tempscale<-aggregate(list(len=sesp[,4]),
                       by = list(ev2=sesp[,1],ev1 = sesp[,2]),
                       FUN = function(x) c(length=length(x)))
  tempscale <- do.call(data.frame, tempscale)
  
  tempscale$combin=paste(tempscale$ev1,tempscale$ev2)
  
  spacescale<-aggregate(list(len=tesp[,4]),
                       by = list(ev2=tesp[,1],ev1 = tesp[,2]),
                       FUN = function(x) c(length=length(x)))
 spacescale <- do.call(data.frame, spacescale)
  
  tempscale$combin=paste(tempscale$ev1,tempscale$ev2)
  names(tempscale)[3]="timescale"
  names(spacescale)[3]="spacescale"
  tscale<-inner_join(tempscale,spacescale,by=c("ev1","ev2"))
  
  spTiming<-data.frame(tempscale$combin,tempscale$ev1,tempscale$ev2)
  names(spTiming)=c("combin","ev1","ev2") 
  sesp$combin=paste(sesp$ev1,sesp$ev2)
  
  for (kv in 1:length(spTiming$combin)){
    kev<-spTiming$combin[kv]
    idc<-which(!is.na(match(sesp$combin,kev)))
    inev<-sesp$dd[idc[1]]
    endev<-sesp$dd[idc[length(idc)]]
    spTiming$startime[kv]=inev
    spTiming$endtime[kv]=endev
  }
 
  for (kv in 1:length(spTiming$ev1)){
    kev<-spTiming$ev1[kv]
    idc<-which(!is.na(match(metavHax[[1]]$cluster,kev))) 
    inev<-metavHax[[1]]$time[idc[1]]
    endev<-metavHax[[1]]$time[idc[length(idc)]]
    spTiming$startimeRain[kv]=inev
    spTiming$endtimeRain[kv]=endev
  } 
  
  for (kv in 1:length(spTiming$ev2)){
    kev<-spTiming$ev2[kv]
    idc<-which(!is.na(match(metavHax[[2]]$cluster,kev))) 
    inev<-metavHax[[2]]$time[idc[1]]
    endev<-metavHax[[2]]$time[idc[length(idc)]]
    spTiming$startimeWind[kv]=inev
    spTiming$endtimeWind[kv]=endev
  }
  
  # for (kv in 1:length(endgame4$ev2)){
  #   kev<-endgame4$ev1[kv]
  #   idc<-which(!is.na(match(sp25$ev.y,kev))) 
  #   inev<-endgame4$startime[kv]
  #   endev<-endgame4$endtime[kv]
  #   rev<-sp25[idc,]
  #   # [which(sp25$dd[idc]>=inev & sp25$dd[idc]<=endev),] 
  #   if(length(rev$ev.y)>0){
  # 
  #     reve<-which(rev$vecmeta==max(rev$vecmeta))
  #     ttime<-rev$dd[reve]
  #     rx<-rev$vecmeta[reve]
  #     locx<-rev$loc[reve] 
  #     endgame4$rf.max.m[kv]=rx
  #     endgame4$MaxtimeRain[kv]=ttime
  #     endgame4$locmaxR[kv]=locx
  #     raincompound<-rbind(raincompound,rev)
  #   }else{
  #     endgame4$MaxtimeRain[kv]=NA
  #     endgame4$rf.max.m[kv]=NA
  #     endgame4$locmaxR[kv]=NA
  #   }
  # }
  # 
  # for (kv in 1:length(endgame4$ev2)){
  #   kev<-endgame4$ev2[kv]
  #   idc<-which(!is.na(match(sp251$ev.x,kev))) 
  #   idf<-which(!is.na(match(metavHax[[2]]$cluster,kev))) 
  #   inev<-endgame4$startime[kv]
  #   endev<-endgame4$endtime[kv]
  #   rev<-sp251[idc,]
  #   revr<-metavHax[[2]][idf,][which(metavHax[[2]]$time[idf]>=inev & metavHax[[2]]$time[idf]<=endev),]
  #   if(length(rev$ev.x)>0){
  #   reve<-which(rev$vecmeta2==max(rev$vecmeta2))
  #   ttime<-rev$dd[reve]
  #   locx<-rev$loc[reve]
  #   mgx<-rev$vecmeta2[reve]
  #   endgame4$MaxtimeWind[kv]=ttime
  #   endgame4$wg.max.m[kv]=mgx
  #   endgame4$locmaxW[kv]=locx
  #   }else{
  #     endgame4$MaxtimeWind[kv]=NA
  #     endgame4$wg.max.m[kv]=NA
  #     endgame4$locmaxW[kv]=NA
  #   }
  # }

  # spTiming$MaxtimeWind<-as_datetime(c(spTiming$MaxtimeWind),origin="1970-01-01")
  # spTiming$MaxtimeRain<-as_datetime(c(spTiming$MaxtimeRain),origin="1970-01-01")
  # 
  spTiming$startime<-as_datetime(c(spTiming$startime),origin="1970-01-01")
  spTiming$endtime<-as_datetime(c(spTiming$endtime),origin="1970-01-01")
  spTiming$startimeWind<-as_datetime(c(spTiming$startimeWind),origin="1970-01-01")
  spTiming$endtimeWind<-as_datetime(c(spTiming$endtimeWind),origin="1970-01-01")
  spTiming$startimeRain<-as_datetime(c(spTiming$startimeRain),origin="1970-01-01")
  spTiming$endtimeRain<-as_datetime(c(spTiming$endtimeRain),origin="1970-01-01")
  
  # plot((endgame4$endtime-endgame4$MaxtimeRain)/3600)
  
  spSTP<-full_join(spTiming,spSTF,by=c("ev1","ev2"),all=T)
  
#   #Finito endgame4 is the compound events set
# metaHaz[[1]][which(metaHaz[[1]]$ev==3443),]
# endgame[which(endgame4$ev1==3443),]
# 
#   #Max and mins for the compounds
# nmax2<-c()
# for(omg in 1:length(endgame4$ev1)){
#   perevent<-metavHax[[2]][which(metavHax[[2]]$cluster==endgame4$ev2[omg] & metavHax[[2]]$time>=endgame4$startime[omg] & metavHax[[2]]$time<=endgame4$endtime[omg]),]
# 
# # perevent2<-perevent[which(perevent$time>=endgame4$startime[omg] & perevent$time<=endgame4$endtime[omg]),]
# 
# newmax<-aggregate(list(rain=perevent$vecmeta,wind=perevent$vecmeta2),
#                   by = list(ev1=perevent$cloc),
#                   FUN = function(x) c(sum=sum(x),max=max(x)))
# newmax <- do.call(data.frame, newmax)
# nmax<-c(max(newmax$rain.sum),max(newmax$wind.max))
# endgame5$wg.max.c1[omg]<-nmax[2]
# }
# 
# nmax1<-c()
# for(omg in 1:length(endgame4$ev1)){
#     perevent<-metavHax[[1]][which(metavHax[[1]]$cluster==endgame4$ev1[omg] & metavHax[[1]]$time>=endgame4$startime[omg] & metavHax[[1]]$time<=endgame4$endtime[omg]),]
# 
#     # perevent2<-perevent[which(perevent$time>=endgame4$startime[omg] & perevent$time<=endgame4$endtime[omg]),]
# 
#     newmax<-aggregate(list(rain=perevent$vecmeta,wind=perevent$vecmeta2),
#                       by = list(ev1=perevent$cloc),
#                       FUN = function(x) c(sum=sum(x),max=max(x)))
#     newmax <- do.call(data.frame, newmax)
#     nmax<-c(max(newmax$rain.sum),max(newmax$wind.max))
#     endgame5$rf.max.c1[omg]<-nmax[1]
# }

names(metaHax[[1]])
allez<-metaHax[[1]][,c(1,3)][which(!is.na(match(metaHax[[1]]$ev,spSTP$ev1))),]
lesbleu<-metaHax[[2]][,c(1,8)][which(!is.na(match(metaHax[[2]]$ev,spSTP$ev2))),]
names(lesbleu)<-c("ev2","wg.max.a")
names(allez)<-c("ev1","rf.max.a")
spSTP<-inner_join(spSTP,lesbleu,by=c("ev2"),all=T)
spSTP<-inner_join(spSTP,allez,by=c("ev1"),all=T)
compound<-spSTP

# compound=data.frame(compound,nmax1[,1],nmax2[,2])
# names(compound)=c("windevent","rainevent","space","rf.max.s","wg.max.s","cev","time","startime",
#                   "endtime","startimeRain","endtimeRain","startimeWind","endtimeWind","rf.max.h",
#                   "maxtimeRain","locmaxRain","maxtimeWind","wg.max.h","locmaxWind","space.cr","olocmaxR",
#                   "space.cw","olocmaxW","rf.max.st","wg.max.st","wg.max.t","rf.max.t","wg.max.a","rf.max.a")
plot(compound$rf.max.c,compound$rf.max.a)
plot(metaHaz[[1]]$wg.max,metaHaz[[1]]$rf.max,pch=16,col=alpha("blue",.5))
points(metaHaz[[2]]$wg.max,metaHaz[[2]]$rf.max,pch=16,col=alpha("red",.5))
points(compound$wg.max.c,compound$rf.max.c,pch=16,col=alpha("green",.5))
abline(a=0,b=1)

comprain<-na.omit(match(metaHax[[1]]$ev,compound$ev1))
comprain2<-na.omit(match(compound$ev1,metaHax[[1]]$ev))
cocor<-metaHax[[1]][comprain,]
cocor2<-metaHax[[1]][comprain2,]
cocor3<-data.frame(compound,cocor2[,c(13:21)])

compwind<-na.omit(match(metaHax[[2]]$ev,compound$ev2))
compwind2<-na.omit(match(compound$ev2,metaHax[[2]]$ev))
cocow<-metaHax[[2]][compwind,]
cocow2<-metaHax[[2]][compwind2,]
# compwind<-inner_join(metaHax[[2]],compound, by=c("ev"="ev2"))
names(cocow2)
cocow3<-cocow2[,c(13,14,22,23)]
colnames(cocow3)


compoundfinalST<-data.frame(cocor3,cocow3)
compoundfinalST<-inner_join(compoundfinalST,tscale,by=c("ev1","ev2","combin"))

names(spSCW2)[c(1,2)]<-c("ev1","ev2")
compoundfinalS<-inner_join(spSCR2,spSCW2,by=c("ev1","ev2"))
names(spTCW2)[c(1,2)]<-c("ev1","ev2")
compoundfinalT<-inner_join(spTCR2,spTCW2,by=c("ev1","ev2"))


names(compoundfinalST)[c(20,21,29,30)]<-c("rf.space","rf.time","wg.space","wg.time")
compoundfinalST$ORspace=compoundfinalST$rf.space+compoundfinalST$wg.space-compoundfinalST$spacescale
for (ti in 1:length(compoundfinalST$ev2)){
compoundfinalST$ORtime[ti]<-difftime(max(compoundfinalST$endtimeRain[ti],compoundfinalST$endtimeWind[ti]),min(compoundfinalST$startimeRain[ti],compoundfinalST$startimeWind[ti]),unit="hours")+1}

compoundfinalST$ratiospace=compoundfinalST$spacescale/compoundfinalST$ORspace

compoundfinalST$ratiotime=compoundfinalST$timescale/compoundfinalST$ORtime
plot(compoundfinalST$ratiospace,compoundfinalST$ratiotime)
compoundfinalST$s="rain"
compoundfinalST$s[which(compoundfinalST$wg.space>compoundfinalST$rf.space)]="wind"

compoundfinalST$t="rain"
compoundfinalST$t[which(compoundfinalST$wg.time>compoundfinalST$rf.time)]="wind"

compoundfinalST$hstar="rain"
compoundfinalST$hstar[which(compoundfinalST$startimeRain>compoundfinalST$startimeWind)]="wind"

compoundfinalST$hend="rain"
compoundfinalST$hend[which(compoundfinalST$endtimeRain<compoundfinalST$endtimeWind)]="wind"

compoundfinalST$stdom<-paste(compoundfinalST$t,compoundfinalST$s)
compoundfinalST$stseq<-paste(compoundfinalST$hstar,compoundfinalST$hend)


plot(compoundfinalT$wg.max.y,compoundfinalT$rf.max.x)
plot(compoundfinalS$wg.max.y,compoundfinalS$rf.max.x)
plot(compoundfinalST$wg.max.c,compoundfinalST$rf.max.c)

compoundfinalST$maxlocstW[which(is.na(compoundfinalST$maxlocstW))]="00 00"
compoundfinalST$maxlocstR[which(is.na(compoundfinalST$maxlocstR))]="00 00"

maxwindc<-matrix(as.numeric(unlist(strsplit(compoundfinalST$maxlocstW," "))),ncol=2,byrow = T)
maxrainc<-matrix(as.numeric(unlist(strsplit(compoundfinalST$maxlocstR," "))),ncol=2,byrow = T)
compoundfinalST<-data.frame(compoundfinalST,maxrainc,maxwindc)

maxwinds<-matrix(as.numeric(unlist(strsplit(compoundfinalS$maxlocsW," "))),ncol=2,byrow = T)
maxrains<-matrix(as.numeric(unlist(strsplit(compoundfinalS$maxlocsR," "))),ncol=2,byrow = T)
compoundfinalS<-data.frame(compoundfinalS,maxrains,maxwinds)

maxwindt<-matrix(as.numeric(unlist(strsplit(compoundfinalT$maxlocTW," "))),ncol=2,byrow = T)
maxraint<-matrix(as.numeric(unlist(strsplit(compoundfinalT$maxloctR," "))),ncol=2,byrow = T)
compoundfinalT<-data.frame(compoundfinalT,maxraint,maxwindt)

plot(compoundfinalST$ORspace,compoundfinalST$ORtime,log="xy")
points(compoundfinalST$rf.space,compoundfinalST$rf.time,col=2)
points(compoundfinalST$spacescale,compoundfinalST$timescale,col=3)


maxspr<-maxrainc[-which(maxrainc[,2]==0),]
geodistc<-distGeo(compoundfinalST[,c(27,28)],compoundfinalST[,c(31,32)])/1000
ORandistRain<-distGeo(compoundfinalST[,c(27,28)],maxrains)/1000
hist(ORandistRain[-which(ORandistRain==0)],breaks=50)

geodistct<-distGeo(maxraint,maxwindt)/1000
geodistct[which(geodistct>2000)]=NA
geodistcs<-distGeo(maxrains,maxwinds)/1000
geodistcs[which(geodistcs>2000)]=NA
geodistcst<-distGeo(maxrainc,maxwindc)/1000
geodistcst[which(geodistcst>2000)]=NA
plot(geodistcs,geodistcst)

length(na.omit(geodistc2))
abline(a=0,b=1)
timedist=difftime(compoundfinalST$startimeWind, compoundfinalST$startimeRain,unit="hour")
hist(as.numeric(timedist),breaks=60)
geoshit<-data.frame(compoundfinalST$combin, geodistc,geodistcs,geodistct,geodistcst)



stdist<-data.frame(geoshit,as.numeric(timedist),as.character(compoundfinalST$season))
names(stdist)[c(1,6,7)]<-c("cev","timedist","season")
stdist$season=as.character(stdist$season)

#This is a result: classification if maximum OR are the same as maximum AND: 17% have both maximum of rain and wind in comp area
#Interesting stuff with space logged 
#I have to choose which metric to use, I have: 1. which is larger(space), longer(time). 2. which is first, which is las(seqence). 3. Is the total max in the compound area for both hazards.

compoundfinalST$maxin="none"
compoundfinalST$maxin[which(compoundfinalS$wg.max.y==compoundfinalST$wg.max.a)]="wind"
compoundfinalST$maxin[which(compoundfinalS$rf.max.x==compoundfinalST$rf.max.a)]="rain"
compoundfinalST$maxin[which(compoundfinalS$rf.max.x==compoundfinalST$rf.max.a & compoundfinalS$wg.max.y==compoundfinalST$wg.max.a)]="rain+wind"

length(which(compoundfinalS$rf.max.x==compoundfinalST$rf.max.a & compoundfinalS$wg.max.y==compoundfinalST$wg.max.a))/length(compoundfinalST$rf.max.a)

compoundfinalST$season=as.character(compoundfinalST$season)

ggplot(compoundfinalST, aes(x=timescale,y=..density..,fill=maxin)) + 
  # geom_histogram(binwidth=.4, position="identity", alpha=0.4,color="darkgreen")+
  geom_density(alpha=.6,size=.4,position = "fill",adjust=2)
  # geom_histogram(aes(x=geodistcst,y=..density..),binwidth=20,fill="skyblue",color="blue",position="dodge", alpha=0.4)+ 
  #   geom_density(aes(x=geodistcst,y=..density..),alpha=.1,fill="darkblue",size=.4)+
  geom_histogram(aes(x=geodistcs,y=..density..),binwidth=20,color="red",fill="orange",position="dodge", alpha=0.4)+
  geom_density(aes(x=geodistcs,y=..density..),alpha=.1,fill="pink",size=.4)


coord_radar <- function (theta = "x", start = 0, direction = 1) {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}
ggplot(compoundfinalST, aes(x=log(spacescale),group=season,fill=season))+
geom_density(adjust=1,position = "fill",alpha=1,n=500)

windst<-metaHax[[2]][,c(6,14)]
windst<-metaHax[[2]][,c(6,14)]
# hist(geodistc,breaks=30)
# hist(as.numeric(timedist),breaks=30)


plot(compoundfinal$time,compoundfinal$space,cex=2,pch=16,col=alpha(3,.6),log="xy")
abline(a=0,b=1)
points(compoundfinal$rf.time,compoundfinal$time,cex=2,pch=16,col=alpha(2,.6))
plot(stdist$geodistc,abs(stdist$timedist))
plot(as.numeric(timedist),geodistc)
xtracomp<-compoundfinal[which(geodistc<=200 & as.numeric(timedist)<=6 & as.numeric(timedist)>=-6),]

cacacomp<-compoundfinal[which(geodistc>300 & (as.numeric(timedist)>=12 | as.numeric(timedist)<=-12)),]

plot(compoundfinal$rf.max,compoundfinal$wg.max,pch=16)
points(xtracomp$rf.max,xtracomp$wg.max,pch=16,col=2)
points(cacacomp$rf.max,cacacomp$wg.max,pch=16,col=3)


m <- ggplot(compoundfinalS, aes(x = ratiospace,y = ratiotime )) +
  geom_point(aes(color=wg.max.a,size=rf.max.a),alpha=.5) 
m + scale_color_gradientn(colours = rgb.palette(30))

m + geom_bin2d(bins=15)
m + stat_density_2d(aes(fill = stat(level)), geom = "polygon")
m+stat_density_2d(geom = "raster", aes(fill = stat(density)),alpha=1, contour = FALSE,n=400)+
  theme_bw()+
  scale_fill_gradientn(colours = rgb.palette(30))


stseason<-compoundfinal
var.int=stseason$rf.max
qr<-quantile(var.int,.95)
mlev<-stseason[which(var.int>=qr),]
mlev$group="r"
plot(var.int[mlev])
idmr<-stseason$rainevent[mlev]
idmw<-stseason$windevent[mlev]

var.int=stseason$wg.max
qw<-quantile(var.int,.95)
mlev1<-stseason[which(var.int>=qw),]
mlev1$group="w"

bleu<-rbind(mlev,mlev1)

plot(mlev,mlev1)

e <- ggplot(compoundfinal, aes(y = wg.max, x = t))
e + geom_boxplot()



hist(compoundfinal$ratiotime[which(compoundfinal$t=="wind")],breaks=25)
hist(compoundfinal$ratiospace,breaks=25)
plot(compoundfinal$windevent,compoundfinal$rainevent)
length(unique(compoundfinal$windevent))

rainpw<-aggregate(list(rr=compoundfinal[,3]),
                 by = list(windev=compoundfinal$windevent), 
                 FUN = function(x) c(length=length(unique(x))))
rainpw <- do.call(data.frame, rainpw)

windpr<-aggregate(list(rr=compoundfinal[,3]),
                  by = list(windev=compoundfinal$rainevent), 
                  FUN = function(x) c(length=length(unique(x))))
windpr <- do.call(data.frame, windpr)
plot(windpr$rr)
plot(rainpw$rr)

######################################################################
# names(metaHaz)=c("Pr-Events","Wind_Events","Compound")
# names(metavHour)=c("Pr-Events","Wind_Events","Compound")

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
  scale_fill_gradientn(colours = rgb.palette(100))+
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.2)+
  # geom_bin2d(data=metaHaz[[2]],aes(x=lonMW,y=latMW,group=season),alpha=.7,binwidth=c(.5,.5))+
# geom_density2d(data=metaHaz[[2]],aes(x=lonMW,y=latMW,group=season),alpha=.2)
stat_density_2d(data=xtracomp,aes(x=lonMR,y=latMR,group=1,fill = stat(density)),geom = "raster", contour = FALSE,alpha=.6,n=300)


ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill = "white", color = "gray10", size = 1) +
  theme_bw(16)+
  scale_fill_gradientn(colours = rgb.palette(100))+
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.2)+
  # geom_bin2d(data=metaHaz[[2]],aes(x=lonMW,y=latMW,group=season),alpha=.7,binwidth=c(.5,.5))+
  # geom_density2d(data=metaHaz[[2]],aes(x=lonMW,y=latMW,group=season),alpha=.2)
  stat_density_2d(data=compoundfinal,aes(x=lonMW.1,y=latMW.1,group=1,fill = stat(density)),geom = "raster", contour = FALSE,alpha=.6,n=300)
sp03$loc=as.character(sp03$loc)
loco<-matrix(as.numeric(unlist(strsplit(sp03$loc," "))),ncol=2,byrow = T)
sp03<-data.frame(sp03,loco)
allprec<-sp03

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
  stat_density_2d(data=allprec,aes(x=X1,y=X2,group=season,fill = stat(density)),geom = "raster", contour = FALSE,alpha=.6,n=150)


rbPal2 <- colorRampPalette(c('red',"darkorange","gold",'darkgreen'))
Col2 <- rbPal2(1000)[as.numeric(cut(log(metaHaz[[1]]$Dmax+0.00001),breaks = 1000))] 
plot(metaHaz$`Pr-Events`$wg.max,metaHaz$`Pr-Events`$rf.max,col=alpha(Col2,.4),pch=16,log="x")




####Polat plot
tabchoos<-compoundfinalST

polarEvent<-aggregate(tabchoos$season,
                 by = list(Month = tabchoos$x.month,Year=tabchoos$x.year),
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


camembert<-aggregate(tabchoos$ev2,
                      by = list(season =tabchoos$season),
                      FUN = function(x) c(n =length(x)))
brie <- do.call(data.frame, camembert)
brie$seasonc<-c("Winter","Spring","Summer","Autum")

p1 <- ggplot(brie, aes(x = "", y =x , fill = seasonc)) + geom_bar(width = 1, stat = 'identity') + coord_polar("y", start = 0)
p1+ scale_fill_brewer("Blues") + 
  theme_bw() +
  geom_text(aes(y = x/3 + c(0, cumsum(x)[-length(x)]),label = percent(x/100)), size = 5)
print(p1)

#Seasonal analysis
bibcase<-metaHax[[2]]
bibis<-aggregate(bibcase$season,
                 by = list(Month = bibcase$x.month),
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

stseason<-compoundfinalST
var.int=stseason$time
hist(var.int,breaks=c(0,6,12,24,96))
stseason$sizegr<-1
stseason$sizegr[which(stseason$x.dur>6 & stseason$x.dur<13)]=2
stseason$sizegr[which(stseason$x.dur>12 & stseason$x.dur<25)]=3
stseason$sizegr[which(stseason$x.dur>24)]=4



gr1<-stseason[which(stseason$time<7),]
gr2<-stseason[which(stseason$time>6 & stseason$time<13),]
gr3<-stseason[which(stseason$time>12 & stseason$time<25),]
gr4<-stseason[which(stseason$time>24),]

plot(gr1$wg.max,gr1$rf.max,xlim=c(0,50),ylim=c(0,130))
points(gr2$wg.max,gr2$rf.max,col=2)
points(gr3$wg.max,gr3$rf.max,col=3)
points(gr4$wg.max,gr4$rf.max,col=4)


plot(gr1$space,gr1$time,log="xy",xlim=c(1,10000),ylim=c(1,130))
points(gr2$space,gr2$time,col=2)
points(gr3$space,gr3$time,col=3)
points(gr4$space,gr4$time,col=4)





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
stseason=metaHax[[1]]
stseason2=metaHax[[2]]
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
  idm=1
  metamax1<-stseason[which(stseason$ev==idm),]
  metamax2<-stseason2[which(stseason2$ev==idm),]
  maxevent1<-metavHax[[1]][which(metavHax[[1]]$cluster==idm),]
  maxevent2<-metavHax[[2]][which(metavHax[[2]]$cluster==idm),]
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
chc<-chull(compoundfinalST$spacescale,compoundfinalST$timescale)
chc<-data.frame(compoundfinalST$spacescale[c(chc,chc[1])],compoundfinalST$timescale[c(chc,chc[1])])


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


names(chr)=c("s","t")
names(chw)=c("s","t")
names(chc)=c("s","t")
plot(chr,col=2,lwd=3,type="l",xlim=c(1,40000),ylim=c(1,120),log="xy")
lines(chw,col=3,lwd=3)
lines(chc,col=4,lwd=3)

ggplot(chr,aes(x=s,y=t))+
  geom_polygon(fill="blue",alpha=.4)+
  geom_polygon(data=chw, aes(x=s,y=t),fill="red",alpha=.4)+
  geom_polygon(data=chc, aes(x=s,y=t),fill="green",alpha=.4)+
  scale_y_continuous(trans="log")+
  scale_x_continuous(trans="log")
  




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

compoundfinal$tdist<-stdist$timedist
compoundfinal$sdist<-stdist$geodistc

###################Spatiotemporal plot
compoundfinalST$wg.max.s=compoundfinalS$wg.max.y
compoundfinalST$rf.max.s=compoundfinalS$rf.max.x

rbPal <- colorRampPalette(c('lightskyblue',"skyblue","gold","darkorange",'red',"purple"))
# ggplot(data=compoundfinalST,aes(x=spacescale,y=ORtime,colour=wg.max.s,size=rf.max.s))+
max(metaHax[[1]]$rf.surf)
ggplot(data=metaHax[[1]],aes(x=rf.surf,y=x.dur,colour=wg.max,size=rf.max))+
geom_point(alpha=.5)+
  theme(axis.text=element_text(size=16),
                axis.title=element_text(size=18,face="bold"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  # geom_hline(yintercept=54)+
  # geom_vline(xintercept=1306)+
scale_colour_gradientn(trans=scales::modulus_trans(1.5),colors=rbPal(100),"Max wind gust [m/s]",breaks=c(10,20,30,40,50,Inf),limits=c(1,51))+
  scale_size(trans=scales::modulus_trans(1.3),range=c(.4,17),"Max precipitation [mm]",breaks=c(20,60,120,Inf),limits=c(0,140),guide = guide_legend(override.aes = list(colour = alpha("grey",.8))))+
  scale_y_continuous(trans = log_trans(),
                     breaks = c(3,6,12,24,48,72,148),limits = c(2.9,150),"Temporal scale",
                     labels=c("3h","6h","12h","1day","2days","3days","1week"))+
scale_x_continuous(trans = log_trans(),
                   breaks =c(1,10,100,1000),limits=c(1,3000),"Spatial scale [km2]",
                   labels=c("1","10","100","1000")) 

compoundfinal$tdist[which(compoundfinal$tdist==0)]=0.1
compoundfinal$timemin<-60*compoundfinal$time

ggplot(data=compoundfinal,aes(x=space,y=timemin,colour=abs(1/tdist),size=1/sdist))+
  geom_point(alpha=.5)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_colour_gradientn(trans=scales::boxcox_trans(.2),colors=rbPal(100),"1/timediff",breaks=c(0.01,.5,1,10),limits=c(0.01,10))+
  scale_size(trans=scales::boxcox_trans(.9),range=c(.4,25),"1/spacediff",breaks=c(0.01,0.05,.100),limits=c(0,.1),guide = guide_legend(override.aes = list(colour = alpha("grey",.8))))+
  scale_y_continuous(trans = log_trans(),
                     breaks = c(3*60,6*60,12*60,24*60,48*60,72*60,148*60),limits = c(0.9*60,150*60),"Temporal scale",
                     labels=c("3h","6h","12h","1day","2days","3days","1week"))+
  scale_x_continuous(trans = log_trans(),
                     breaks =c(1000,10000,100000,1000000)/450,limits=c(1,3000),"Spatial scale [km2]",
                     labels=c("1000","10000","100000","1000000")) 

min(abs(1/compoundfinal$tdist))
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
# for (haz in 1:3){
# stseason<-metaHaz[[haz]]


stseason<-compoundfinal
var.int=stseason$rf.max
mlev<-order(-var.int)
plot(var.int[mlev])
idmr<-stseason$rainevent[mlev]
idmw<-stseason$windevent[mlev]
datesr<-c()
datesw<-c()
for(idi in 1:10){
print(idi)
idm1=idmr[idi]
idm2=idmw[idi]
maxeventr<-metavHour[[1]][which(metavHour[[1]]$cluster==idm1),]
maxeventw<-metavHour[[2]][which(metavHour[[2]]$cluster==idm2),]

maxeventr$time <- as_datetime(c(maxeventr$Var3*60*60),origin="1900-01-01")
maxeventw$time <- as_datetime(c(maxeventw$Var3*60*60),origin="1900-01-01")
print(idm1)
print(idm2)
print(c(maxeventr$time[1],maxeventr$time[length(maxeventr$time)]))
print(c(maxeventw$time[1],maxeventw$time[length(maxeventw$time)]))
datesr<-rbind(datesr,c(maxeventr$time[1],maxeventr$time[length(maxeventr$time)]))
datesw<-rbind(datesw,c(maxeventw$time[1],maxeventw$time[length(maxeventr$time)]))
}
datesr<-as.data.frame(datesr)
datesw<-as.data.frame(datesw)
dates[,1]<-as_datetime(dates[,1],origin="1970-01-01")
dates[,2]<-as_datetime(dates[,2],origin="1970-01-01")
names(dates)=c("start","end")
totdates<-c(totdates,list(dates))


# names(totdates)=c("Pr-Events","Wind_Events","Compound")
# plot(totdates[[1]][,1])
# points(totdates[[3]][,1],col=2,pch=3)
# points(totdates[[2]][,1],col=3,pch=4)
# 
# write.csv(totdates[[3]], file = "windcomptop10valid.csv")





stseason<-compoundfinal
idi=770
idmr<-stseason$ev1
idmw<-stseason$ev2
print(idi)
idm1=idmr[idi]
idm2=idmw[idi]
maxeventr<-metavHour[[1]][which(metavHax[[1]]$cluster==idm1 & metavHax[[1]]$time>=stseason$startime[idi] & metavHax[[1]]$time<=stseason$endtime[idi] ),]
maxeventw<-metavHour[[2]][which(metavHax[[2]]$cluster==idm2 & metavHax[[2]]$time>=stseason$startime[idi] & metavHax[[2]]$time<=stseason$endtime[idi] ),]

maxeventr<-metavHax[[1]][which(metavHax[[1]]$cluster==idm1),]
maxeventw<-metavHax[[2]][which(metavHax[[2]]$cluster==idm2),]
compou<-compoundfinal[which(compoundfinal$combin==paste(idm1,idm2)),]

  
maxeventr$time <- as_datetime(c(maxeventr$Var3*60*60),origin="1900-01-01")
maxeventw$time <- as_datetime(c(maxeventw$Var3*60*60),origin="1900-01-01")

maxeventri<- maxeventr[which(maxeventr$time==compou$startime),]
maxeventrf<- maxeventr[which(maxeventr$time==compou$endtime),]

maxeventwi<- maxeventw[which(maxeventw$time==compou$startime),]
maxeventwf<- maxeventw[which(maxeventw$time==compou$endtime),]

# maxeventri<- maxeventr
# maxeventrf<- maxeventr
# 
# maxeventwi<- maxeventw
# maxeventwf<- maxeventw

compou$sr=min(maxeventr$time)
compou$er=max(maxeventr$time)
compou$sw=min(maxeventw$time)
compou$ew=max(maxeventw$time)
print(maxeventr$time[1])

#######################TRAJECTORY######################################
# 
# for (d in 1:length(unique(maxevent$time))){
#   kik<-maxevent[which(maxevent$time==unique(maxevent$time)[d]),]
# # points(kik$Var1,kik$Var2,col=d)
# 
#     centi<-as.numeric(kik[which(kik$vecmeta==max(kik$vecmeta))[1],c(1,2)])
#   idkik<-rbind(idkik,c(d,max(kik$vecmeta),length(kik$vecmeta),centi))
#   points(centi[1],centi[2],col=colvect[d],pch=16)
# }
# 
# idkik=data.frame(idkik)
# 
# 
# uk_fort <- ggplot2::fortify(ukk)
# # ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
# #   geom_polygon(fill = "white", color = "gray10", size = 1) +
# #   theme_bw(16)+
# #   scale_color_viridis(option="C", alpha=.8)+
# #   geom_point(data=idkik,aes(x=X4,y=X5,group=X1,color=X1,size=X2))+
# #   coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3)
# 
# 
# idkik$gr<-1
# idkik$rtime<-1
# # ggplot(idkik, aes(X4, X5))+
# #   geom_path(aes(group=gr), arrow = arrow())+
# #   coord_fixed()
# 
# f1_mod <- loess(X4 ~ X1, data = idkik,span=0.5)
# f2_mod <- loess(X5 ~ X1, data = idkik,span=0.4)
# pred <- data.frame(X1 = seq(1,length(idkik$X1), length = 100))
# pred$F1 <- predict(f1_mod, newdata = pred)
# pred$F2 <- predict(f2_mod, newdata = pred)
# pred$G<-1
# indices <- seq(1,99, by = 5)
# ggplot(pred, aes(F1, F2))+
#   geom_path()+
#   geom_point(data = pred[indices,])+
#   coord_fixed()

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
#####################################


temp01<- hazdat1$data[ , ,1]
tmp_01a <- as.vector(temp01)

lon=hazdat1$lon
lat=hazdat1$lat
lonlat <- as.matrix(expand.grid(lon,lat))

tmp_02a<-paste(maxeventr[,1],maxeventr[,2])
tmp_02b<-paste(maxeventw[,1],maxeventw[,2])
mxs<-c()
for (mx in unique(tmp_02a)){
maxevsum<-sum(maxeventr$vecmeta[which(tmp_02a==mx)])
mxs<-c(mxs,maxevsum)
}

mxv<-c()
for (mv in unique(tmp_02b)){
  maxevacc<-max(maxeventw$vecmeta2[which(tmp_02b==mv)])
  mxv<-c(mxv,maxevacc)
}
tmp_05a<-unique(tmp_02a)
tmp_03<-paste(lonlat[,1],lonlat[,2])
tmp_04a<-match(tmp_05a,tmp_03)
tmp_01a[tmp_04a]<-mxs
tmp_01a[-tmp_04a]<-0
nlat<-dim(lat)
nlon=dim(lon)
tmp_array <- array(tmp_01a, dim=c(nlon,nlat))

  grid <- expand.grid(lon=lon, lat=lat)
  mapmat=tmp_array
  #column 1634 corresponding to Dec 2015
  #This command compresses numbers larger than 6 to 6
  # plot.new()
  tmp_06a<-cbind(mxs,tmp_05a)
  tmp_07a<-maxeventr[match(tmp_05a,tmp_02a),]
  ohcrap<-which(is.na(match(tmp_03,tmp_02a)))
  ouh<-data.frame(lonlat[ohcrap,],0,0,0,0,tmp_07a$time[1],0,0,0,0)

  
  elmaxoua<-tmp_07a
  elmaxoua$maxx=mxs
  elmaxoua$mwg=0
  names(ouh)=names(elmaxoua)
  elmaxoua1=rbind(elmaxoua,ouh)
  elmaxoua$gr<-"0"
  elmaxoua$tf="0"
  timeco<-which(!is.na(match(elmaxoua$cloc,maxeventri$cloc)))
  elmaxoua$gr[timeco]<-"1"
  timeco<-which(!is.na(match(elmaxoua$cloc,maxeventrf$cloc)))
  elmaxoua$gr[timeco]<-"2"
  elmaxoua$tf[which(elmaxoua$time>=compou$startime&elmaxoua$time<=compou$endtime)]<-"1"

  tmp_01b <- as.vector(temp01)
  tmp_05b<-unique(tmp_02b)
  tmp_04b<-match(tmp_05b,tmp_03)
  tmp_01b[tmp_04b]<-mxv
  tmp_01b[-tmp_04b]<-0
  
  nlat<-dim(lat)
  nlon=dim(lon)
  tmp_array <- array(tmp_01b, dim=c(nlon,nlat))
  
  grid <- expand.grid(lon=lon, lat=lat)
  mapmat=tmp_array
  #column 1634 corresponding to Dec 2015
  #This command compresses numbers larger than 6 to 6
  # plot.new()
  # tmp_06b<-cbind(mxv,tmp_01b)
  tmp_07b<-maxeventw[(match(tmp_05b,tmp_02b)),]
  ohcrap<-which(is.na(match(tmp_03,tmp_02b)))
  ouh<-data.frame(lonlat[ohcrap,],0,0,0,0,tmp_07b$time[1],0,0,0,0)
 
  # tmp_07b<-tmp_07b[-which(is.na(tmp_07b$Var1)),]

  elmaxoub<-tmp_07b
  elmaxoub$maxx=NA
  elmaxoub$mwg=mxv
  names(ouh)=names(elmaxoub)
  elmaxoub1=rbind(elmaxoub,ouh)
  elmaxoub$gr<-"0"
  elmaxoub$tf<-"0"
  timeco<-which(!is.na(match(elmaxoub$cloc,maxeventwi$cloc)))
  elmaxoub$gr[timeco]<-"3"
  timeco<-which(!is.na(match(elmaxoub$cloc,maxeventwf$cloc)))
  elmaxoub$gr[timeco]<-"4"
  elmaxoub$tf[which(elmaxoub$time>=compou$startime&elmaxoub$time<=compou$endtime)]<-"2"

  
  bs.palette=colorRampPalette(c("white","aliceblue","skyblue","dodgerblue","dodgerblue4"),interpolate="spline",bias=1)
  bs.palette2=colorRampPalette(c("lightgoldenrod1","goldenrod1" ,"orange","red"),interpolate="spline",bias=1.5)
  bs.palette3=colorRampPalette(c("orange","blue","purple"),interpolate="linear",bias=1)
    library(ggnewscale)
    library(metR)
  metavHax[[2]][which(metavHax[[2]]$cluster==3155),]
  sp251[which(sp251$ev.x==3155),]
  
  ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
    theme_bw(16)+
    geom_raster(data=elmaxoub,aes(x=Var1,y=Var2,fill=mwg,group=cluster,alpha=mwg),interpolate = F)+
    # geom_contour_fill(data=elmaxoub1,aes(x=Var1,y=Var2,z=mwg,group=cluster,alpha=.1),na.fill=T)+
    # geom_contour(data=elmaxoub,colour="black",aes(x=Var1,y=Var2,z=mwg,group=cluster))+
    geom_polygon(fill = "transparent", color = "gray10", size = 1) +
    coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
    scale_fill_gradientn(colours = bs.palette2(100))+
    new_scale_fill() +
      geom_raster(data=elmaxoua,aes(x=Var1,y=Var2,fill=maxx,group=cluster,alpha=maxx),interpolate = F)+
    # geom_contour_fill(data=elmaxoua1,aes(x=Var1,y=Var2,z=maxx,group=cluster),alpha=0.3,na.fill=T)+
    scale_alpha_continuous(range = c(0.2, 0.7),trans="sqrt")+
    scale_fill_gradientn(colours = bs.palette(100))+
    geom_point(data=geoshit[idi,],aes(x=X1,y=X2,group=geodistc),colour="blue")+
    geom_point(data=geoshit[idi,],aes(x=X1.1,y=X2.1,group=geodistc),colour="red")
    # scale_colour_gradientn(colours = bs.palette(100))
  
    ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
      theme_bw(16)+
      geom_polygon(fill = "transparent", color = "gray10", size = 1) +
      coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
      geom_raster(data=elmaxoub,aes(x=Var1,y=Var2,fill="2",group=cluster),alpha=.6,interpolate = F)+
      geom_raster(data=elmaxoua,aes(x=Var1,y=Var2,fill="1",group=cluster),alpha=0.7,interpolate = F)+
      scale_fill_manual(values = c( "1"="blue", "2" ="orange"))+
      geom_point(data=compoundfinalST[idi,],aes(x=lonMR,y=latMR,group=ev1),colour="skyblue",alpha=.6,size=7)+
      geom_point(data=compoundfinalST[idi,],aes(x=lonMW,y=latMW,group=ev1),colour="red",alpha=.6,size=7)+
      geom_point(data=compoundfinalT[idi,],aes(x=X1,y=X2,group=ev1),colour="skyblue",alpha=.6,size=2,shape=2)+
      geom_point(data=compoundfinalT[idi,],aes(x=X1.1,y=X2.1,group=ev1),colour="red",alpha=.6,size=2,shape=2)+
      geom_point(data=compoundfinalS[idi,],aes(x=X1,y=X2,group=ev1),colour="skyblue",alpha=.6,size=2,shape=3)+
      geom_point(data=compoundfinalS[idi,],aes(x=X1.1,y=X2.1,group=ev1),colour="red",alpha=.6,size=2,shape=3)+
      geom_point(data=compoundfinalST[idi,],aes(x=X1,y=X2,group=ev1),colour="skyblue",alpha=.6,size=2,shape=4)+
      geom_point(data=compoundfinalST[idi,],aes(x=X1.1,y=X2.1,group=ev1),colour="red",alpha=.6,size=2,shape=4)

      
    
      
  
    
      scale_alpha_continuous(range = c(0.4, 0.7),trans="exp")+
      scale_fill_manual(values = c("0" = "lightgrey" ,"1"="blue", "2" ="red","3" = "gold","4" ="darkorange"))
  
  ohshit+   
    geom_raster(data=elmaxoub,aes(x=Var1,y=Var2,fill=2,group=cluster),alpha=.5,interpolate = F)+
    # new_scale("fill") +
    scale_fill_gradientn(colours = bs.palette2(10))
    
  
    geom_raster(data=verif,aes(x=Var1,y=Var2,fill=15,group=cluster),alpha=.1,interpolate = F)
  
  ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
    theme_bw(16)+
    geom_raster(data=elmaxoua,aes(x=Var1,y=Var2,fill=maxx,group=cluster),alpha=.9,interpolate = F)+
    geom_polygon(fill = "transparent", color = "gray10", size = 1) +
    coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
    scale_fill_gradientn(colours = rgb.palette(100))
  
  ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
    theme_bw(16)+
    geom_raster(data=verif,aes(x=Var1,y=Var2,fill=1,group=cluster),alpha=.9,interpolate = F)+
    geom_polygon(fill = "transparent", color = "gray10", size = 1) +
    coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
    scale_fill_gradientn(colours = rgb.palette(100))
  

elmaxou$charloc=paste(elmaxou$Var1,elmaxou$Var2) 
oyoy$charloc=paste(oyoy$Var1,oyoy$Var2) 
verif<-match_df(oyoy,elmaxou,on=c("charloc"))

oyoy<-elmaxou

  par(mar=c(4,5,3,0))
  int=seq(0,max(mxv),length.out=100)
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

##############################Cell by cell= where is the most prone to hazard comb?##############

cellcount<-metavDax[[1]][which(!is.na(match(metavDax[[1]]$ev,compoundfinalST$ev1))),]
cellshit<-metavHax[[1]][which(!is.na(match(metavHax[[1]]$cluster,compoundfinalST$ev1))),]

length(unique(cellcount$ev))

cellcount2<-metavDax[[2]][which(!is.na(match(metavDax[[2]]$ev,compoundfinalST$ev2))),]
cellshit2<-metavHax[[2]][which(!is.na(match(metavHax[[2]]$cluster,compoundfinalST$ev2))),]

ccop<-na.omit(match(metavDax[[1]]$ev,compoundfinalST$ev1))
ccop2<-na.omit(match(metavDax[[2]]$ev,compoundfinalST$ev2))
cellcount$megashit<-compound$combin[ccop]
cellcount2$megashit<-compound$combin[ccop2]

# cellshit$vecmeta2=NA
# cellshit2$vecmeta=NA
ccop<-na.omit(match(metavHax[[1]]$cluster,compoundfinalST$ev1))
ccop2<-na.omit(match(metavHax[[2]]$cluster,compoundfinalST$ev2))
cellshit$megashit<-compound$cev[ccop]
cellshit2$megashit<-compound$cev[ccop2]

sp03$combin=paste(sp03$ev1,sp03$ev2)

cellshit3<-as.data.frame(rbind(cellshit,cellshit2))
names(cellcount2)=names(cellcount)

cellshit$loc=paste(cellshit$Var1,cellshit$Var2)
cellshit2<-unique(cellshit$loc)

# compand<-inner_join(cellshit,cellshit2,by=c("megashit","cloc"))
# bib<-c()
# for(loci in 1: length(unique(compand$cloc))){
#   local<-unique(compand$cloc)[loci]
# trial<-compand[which(compand$cloc==local),]
# trial2<-trial[order(trial$cloc),]
# # trial2<-trial2[which(trial2$cloc=="-4 57"),]
# dt<-difftime(trial2$time.x,trial2$time.y,unit="hours")
# length(which(dt==0))
# bib<-c(bib,median(na.omit(dt)))
# }
# 
# bib<-data.frame(bib,unique(compand$cloc))
# mbib<-match(bib$unique.compand.cloc.,compand$cloc)
# bib<-data.frame(bib,compand[mbib,])


cellcorr<-cellshit[na.omit(match(cellshit2,cellshit$loc)),]
head(cellcount)


cellcount3<-as.data.frame(rbind(cellcount,cellcount2))




length(unique(cellcount3$ev))
ucount<-aggregate(cellcount3[,1] ,
                  by = list(ev=cellcount3[,7],loc=cellcount3[,2]),
                  FUN = function(x) c(n =length(x)))

ucount2<-ucount[which(ucount$x==2),]

ucount3<-aggregate(sp03[,1] ,
                  by = list(ev=sp03[,5],loc=sp03[,2]),
                  FUN = function(x) c(n =length(x)))

length(unique(cellcount3$megashit))
tesp$cev=paste(tesp$ev1,tesp$ev2)
ucount2<-tesp

ucount2<-spSTCR2

location<-"-0.5 51.5"

ucouthH<-ucount2[which(ucount2$loc==location),]
mates<-which(!is.na(match(ucount2$cev,ucouthH$cev)))
ucountmH<-ucount2[mates,]

heathrowext<-aggregate(ucountmH[,4] ,
                       by = list(loc=ucountmH[,3]),
                       FUN = function(x) c(n =sum(x)))

heathrowext<-do.call(data.frame,heathrowext)
# heathrowext$x=heathrowext$x/max(heathrowext$x)
heathrowext<-inner_join(heathrowext,cellcorr,by=c("loc"))
heathrowext$combin=heathrowext$megashit

ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  theme_bw(16)+
  geom_raster(data=heathrowext,aes(x=Var1,y=Var2,fill=x,group=cluster),interpolate = F)+
  # geom_contour_fill(data=elmaxoub1,aes(x=Var1,y=Var2,z=mwg,group=cluster,alpha=.1),na.fill=T)+
  # geom_contour(data=elmaxoub,colour="black",aes(x=Var1,y=Var2,z=mwg,group=cluster))+
  geom_polygon(fill = "transparent", color = "gray10", size = 1) +
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
  scale_fill_gradientn(colours = rgb.palette(100))

sptdurloc<-inner_join(heathrowext,compoundfinalST,by=c("combin"))
mean(compoundfinalST$spacescale)
mean(sptdurloc$time.y)
bib<-c()
for(loci in 1: length(unique(ucount2$loc))){
  local<-unique(ucount2$loc)[loci]
  ucouthH<-ucount2[which(ucount2$loc==local),]
  ucouthH$cev=ucouthH$ev
  trial<-semi_join(compoundfinalST,ucouthH,by=c("cev"))
  ohmerd<-c(median(trial$space),median(trial$time))
  bib<-rbind(bib,ohmerd)
}

nrow(bib)
bib<-data.frame(bib,unique(ucount2$loc))
mbib<-match(bib$unique.ucount2.loc.,cellshit$cloc)
bib<-data.frame(bib,cellshit[mbib,])
plot(bib$X1) 
hist(compoundfinalST$space)
median(compoundfinalST$space)

ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  theme_bw(16)+
  geom_raster(data=bib,aes(x=Var1,y=Var2,fill=X1,group=cluster),interpolate = F)+
  # geom_contour_fill(data=elmaxoub1,aes(x=Var1,y=Var2,z=mwg,group=cluster,alpha=.1),na.fill=T)+
  # geom_contour(data=elmaxoub,colour="black",aes(x=Var1,y=Var2,z=mwg,group=cluster))+
  geom_polygon(fill = "transparent", color = "gray10", size = 1) +
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
  scale_fill_gradientn(colours = rgb.palette(100))



bcount<-aggregate(ucount[,3] ,
                 by = list(loc=ucount[,2]),
                 FUN = function(x) c(n =length(x)))
bcount<-do.call(data.frame,bcount)

bcount2<-aggregate(ucount2[,3] ,
                  by = list(loc=ucount2[,2]),
                  FUN = function(x) c(n =length(x)))
bcount2<-do.call(data.frame,bcount2)

bcount$loc=as.character(bcount$loc)
bcoo1<-inner_join(bcount,cellcorr,by=c("loc"))

bcount2$loc=as.character(bcount2$loc)
bcoo2<-inner_join(bcount2,cellcorr,by=c("loc"))


ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  theme_bw(16)+
  geom_raster(data=bcoo2,aes(x=Var1,y=Var2,fill=x,group=cluster),interpolate = F)+
  # geom_contour_fill(data=elmaxoub1,aes(x=Var1,y=Var2,z=mwg,group=cluster,alpha=.1),na.fill=T)+
  # geom_contour(data=elmaxoub,colour="black",aes(x=Var1,y=Var2,z=mwg,group=cluster))+
  geom_polygon(fill = "transparent", color = "gray10", size = 1) +
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
  scale_fill_gradientn(colours = rgb.palette(100))

cellcount<-metavDax[[2]][which(!is.na(match(metavDax[[2]]$ev,compoundfinalST$ev2))),]
cellshit<-metavHax[[2]][which(!is.na(match(metavHax[[2]]$cluster,compoundfinalST$ev2))),]
cellshit$loc=paste(cellshit$Var1,cellshit$Var2)
cellshit2<-unique(cellshit$loc)
cellcorr<-cellshit[na.omit(match(cellshit2,cellshit$loc)),]
head(cellcount)
length(unique(cellcount$ev))
ucount<-aggregate(cellcount[,2] ,
                  by = list(ev=cellcount[,1]),
                  FUN = function(x) c(n =length(x)))

bcount<-aggregate(cellcount[,2] ,
                  by = list(loc=cellcount[,2]),
                  FUN = function(x) c(n =length(x)))
bcount<-do.call(data.frame,bcount)
bcoo<-inner_join(bcount,cellcorr,by=c("loc"))

bcootot<-data.frame(bcoo,bcoo1)
bcootot$xt<-bcootot$x+bcootot$x.1
ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  theme_bw(16)+
  geom_raster(data=bcootot,aes(x=Var1,y=Var2,fill=xt,group=cluster),interpolate = F)+
  # geom_contour_fill(data=elmaxoub1,aes(x=Var1,y=Var2,z=mwg,group=cluster,alpha=.1),na.fill=T)+
  # geom_contour(data=elmaxoub,colour="black",aes(x=Var1,y=Var2,z=mwg,group=cluster))+
  geom_polygon(fill = "transparent", color = "gray10", size = 1) +
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
  scale_fill_gradientn(colours = rgb.palette(100))


cellcount<-metavDax[[1]]
cellshit<-metavHax[[1]]
cellshit$loc=paste(cellshit$Var1,cellshit$Var2)
cellshit2<-unique(cellshit$loc)
cellcorr<-cellshit[na.omit(match(cellshit2,cellshit$loc)),]
head(cellcount)
length(unique(cellcount$ev))

ucount<-aggregate(cellcount[,3] ,
                  by = list(ev=cellcount[,1],loc=cellcount[,2]),
                  FUN = function(x) c(n =length(x)))
ucount<-do.call(data.frame,ucount)

bcount<-aggregate(cellcount[,2] ,
                  by = list(loc=cellcount[,2]),
                  FUN = function(x) c(n =length(x)))
bcount<-do.call(data.frame,bcount)
bcoos<-inner_join(bcount,cellcorr,by=c("loc"))

location<-"-0.5 51.5"
ucouthH<-ucount[which(ucount$loc==location),]
mates<-which(!is.na(match(ucount$ev,ucouthH$ev)))
ucountmH<-ucount[mates,]

heathrowext<-aggregate(ucountmH[,3] ,
                       by = list(loc=ucountmH[,2]),
                       FUN = function(x) c(n =length(x)))

heathrowext<-do.call(data.frame,heathrowext)
heathrowext$x=heathrowext$x/max(heathrowext$x)
heathrowext<-inner_join(heathrowext,cellcorr,by=c("loc"))

ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  theme_bw(16)+
  geom_raster(data=heathrowext,aes(x=Var1,y=Var2,fill=x,group=cluster),interpolate = F)+
  # geom_contour_fill(data=elmaxoub1,aes(x=Var1,y=Var2,z=mwg,group=cluster,alpha=.1),na.fill=T)+
  # geom_contour(data=elmaxoub,colour="black",aes(x=Var1,y=Var2,z=mwg,group=cluster))+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="italic"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  geom_polygon(fill = "transparent", color = "gray10", size = 1) +
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
  scale_fill_gradientn(colours = rgb.palette(100))



ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  theme_bw(16)+
  geom_raster(data=bcoos,aes(x=Var1,y=Var2,fill=x,group=cluster),interpolate = F)+
  # geom_contour_fill(data=elmaxoub1,aes(x=Var1,y=Var2,z=mwg,group=cluster,alpha=.1),na.fill=T)+
  # geom_contour(data=elmaxoub,colour="black",aes(x=Var1,y=Var2,z=mwg,group=cluster))+
  geom_polygon(fill = "transparent", color = "gray10", size = 1) +
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
  scale_fill_gradientn(colours = rgb.palette(100))


#Comparison with spatial dependence method, the chi measure

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
thrr<-as.vector(thr) 
elonlat <- as.matrix(expand.grid(hazmat$Pr$lon,hazmat$Pr$lat))
thrbg<-data.frame(elonlat,thrr)
thrbg$gr=1
thrbg$loc=paste(thrbg$Var1,thrbg$Var2)
lolon<--0.5 
lolat<-51.5
lla<-which(hazmat$Pr$lat==lolat)
llo<-which(hazmat$Pr$lon==lolon)
localdat<-hazmat$Pr$data[llo,lla,]
localth<-thrbg[which(thrbg$loc==location),]
exid<-which(localdat>=localth$thrr)
localex<-localdat[which(localdat>=localth$thrr)]

empichi<-hazmat$Wg$data[,,1]
for (i in (1:33)){
  for(j in 1:45){
    crappy<-hazmat$Pr$data[i,j,]
    thw<-quantile(crappy[which(crappy>0)],.99,na.rm=T)
    abo<-crappy[exid]
    lex<-length(abo[which(abo>=thw)])
    empichi[i,j]<-lex/length(exid)
  }
}

chiemp<-as.vector(empichi) 
elonlat <- as.matrix(expand.grid(hazmat$Pr$lon,hazmat$Pr$lat))
chim<-data.frame(elonlat,chiemp)
chim$gr=1

ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  theme_bw(16)+
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.2)+
  geom_raster(data=chim,aes(x=Var1,y=Var2,fill=chiemp,group=gr),alpha=.8,interpolate = F) +
  scale_fill_gradientn(colours = rgb.palette(100),na.value = "aliceblue")+
  geom_polygon(fill = "transparent", color = "gray10", size = 1.2)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="italic"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_y_continuous(
    breaks = c(48,50,52,54,56,58,60),limits = c(40,70),"Latitude")+
  scale_x_continuous(
    breaks =c(-6,-4,-2,0,2),limits=c(-10,10),"Longitiude") 


######Now try to look at temporal overlap for 1 cell= lagtime between rainfall and wind in the cell (histogram)



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



