rm(list=ls()) 
gc()
getwd()
setwd("C:/Users/PhD Student/OneDrive - King's College London/DATA/04_Spatiotemporal")
setwd("C:/Users/k1638615/King's College London/OneDrive - King's College London/DATA/04_Spatiotemporal")
setwd("C:/Users/Alois/OneDrive - King's College London/DATA/04_Spatiotemporal")
config_file=paste0(getwd(),"/config/config_general.R")
source(config_file)
library(ncdf4)
library(chron)
library(lattice)
library(RColorBrewer)
library(lubridate)
library(dbscan)
library(dismo)
library(maps)
library(progress)
library(fastmatch)
library(viridis)
library(ggplot2)
library(ggnewscale)
library(metR)
#===========================================================================================================





#I have to create my spatiotemporal cube


##############################################
#What i could do
# Do a sample for OR and find wind dominated, rain dominate and compound events
#1- for each grid cell: how many time within an event
#2- for each event: how long is it within the event (if not contimuous just add)
#3- for each event: what is its size
#4- Use some quantiles to display results in map

newrun=TRUE
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
fileele=paste0(getwd(),"/in/elev_0.1deg.nc")
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
  # coord_map(xlim = longlims,  ylim = latlims, "azequalarea")+
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.2)+
  geom_raster(data=alelu,aes(x=Var1,y=Var2,fill=ele1,group=gr),alpha=.8,interpolate = F) +
  scale_fill_gradientn(colours = terrain.colors(100),na.value = "aliceblue")+
  geom_polygon(fill = "transparent", color = "gray10", size = 1.2)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="italic"),
        panel.background = element_rect(fill = "transparent", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        panel.grid = element_line(colour="black", size=.5),
        legend.position = "none",
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_y_continuous(
    breaks = c(48,50,52,54,56,58,60),limits = c(40,70),"Latitude",  minor_breaks = seq(48, 60, 0.25))+
  scale_x_continuous(
    breaks =c(-6,-4,-2,0,2),limits=c(-10,10),"Longitiude", minor_breaks = seq(-6, 2, 0.25)) 




#################################################################################

filer=paste0(getwd(),"/in/Rain_2009-2019.nc")
filew=paste0(getwd(),"/in/windG_2009-2019.nc")
ncr = nc_open(filer)
ncw = nc_open(filew)
# ncatt_get(ncr,name.var,"long_name")
# ncatt_get(ncw,name.var,"long_name")
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


timix<-as_datetime(c(time*60*60),origin="1900-01-01")
l1<-which(hazmat[[1]]$lon==-4.75)
l2<-which(hazmat[[1]]$lat==50.25)
ole<-which(month(timix)==11 & year(timix)==2010)
mierda<-hazmat[[1]]$data[l1,l2,which(hazmat[[1]]$data[l1,l2,]>=-1)]
plot(timix[ole],mierda[ole],type="o")
macs<-mierda[which(mierda>0)]
qtest<-quantile(macs,.99,na.rm=T)
abline(h=qtest)


mac<-c()
idd<-c()
x=1
for (id in c(x+1:(length(mierda)-x))){
  mac[id]<-sum(mierda[c((id-x):(id+x))]) 
  idd[id]<-id
}
f3 <- rep(1/x,x)
mac<- stats::filter(mierda, f3, sides=2) *x
which(mac==max(mac,na.rm=T))
plot(timix[ole],mac[ole],type="o")
oula<-mierda[ole]
tit<-timix[ole]
macs<-mac[which(mac>0)]
qtest<-quantile(macs,.99,na.rm=T)
abline(h=qtest)

idsam<-which(year(timix)==2016)
hazmatx<-hazmat
hazmatx[[1]]$data<-hazmat[[1]]$data[,,idsam]
hazmatx[[2]]$data<-hazmat[[2]]$data[,,idsam]
timixsam<-timix[idsam]
times<-time[idsam]
rbPal <- colorRampPalette(c('darkgreen',"gold","darkorange",'red'))
names(hazmat)=c("Pr","Wg")
names(hazmatx)=c("Pr","Wg")

outix<-c()
pth<-c(.95,.96,.97,.98,.99)
minptv<-c(5,10,20,30)
acoef<-c(2,4,8)
totalcases<-length(pth)*length(minptv)*length(acoef)
aco<-rep(acoef,20)
mco<-rep(minptv,15)
pco<-rep(pth,12)
parset<-cbind(aco,mco,pco)

rm(nc,newformat,timestamp,crappy,alelu,file_time,hazmat)
gc()
for(m in 1:totalcases) {

mint=parset[m,2]
pqt=parset[m,3]
acof=parset[m,1]
#Step 1: Threshold and sampling of extreme events
thr<-hazmatx$Pr$data[,,1]
hazmatx$Pr2<-hazmatx$Pr
for (i in (1:33)){
  for(j in 1:45){
    crappy<-hazmatx$Pr$data[i,j,]
    x=1
    f3 <- rep(1/x,x)
    mac <- stats::filter(crappy, f3, sides=2) *x
    mac<-as.vector(mac)
    
    thr[i,j]<-quantile(mac[which(mac>0)],pqt,na.rm=T) 
    # hazmat$Pr2$data[i,j,]=mac
  } 
}
thrx6<-thr

thw<-hazmatx$Wg$data[,,1]
for (i in (1:33)){
  for(j in 1:45){
    crappy<-hazmatx$Wg$data[i,j,]
    thw[i,j]<-quantile(crappy[which(crappy>0)],pqt,na.rm=T)
  }
}
thrr<-as.vector(thr) 
elonlat <- as.matrix(expand.grid(hazmatx$Pr$lon,hazmatx$Pr$lat))
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

th1<-quantile(hazmatx$Pr$data[which(hazmatx[[1]]$data>0)],0.99,na.rm=T)
th2<-quantile(hazmatx$Wg$data[which(hazmatx[[2]]$data>0)],0.99,na.rm=T)


#####Step 2: Spatiotemporal clustering of wind and rain extremes

metaHaz<-list()
metavHour<-list()
metavDaz<-list()

for(hazard in 1:2){
  if(hazard==1){hazdat=hazmatx$Pr;th=thr}
  if(hazard==2){hazdat=hazmatx$Wg;th=thw}

  wala=F
  lon=hazmatx$Pr$lon
  lat=hazmatx$Pr$lat
  
  lonlatime <- expand.grid(lon, lat,times)
  
  formeta<-hazmatx[[hazard]]$data

    vectouf<- as.vector(hazdat$data)
    vecthouf<-as.vector(rep(th,length(hazdat$data[1,1,]))) 
    vectouf[which(vectouf<vecthouf)]<-NA
    vectouf[which(vectouf>=vecthouf)]<-1
    vecmeta<-as.vector(formeta)
    # vecmeta2<-as.vector(formeta2)
    if(hazard==1){vecmeta[which(is.na(vectouf))] <- NA}
    if(hazard==2){vecmeta[which(vecmeta<vecthouf)]<- NA}
  
  rm(formeta,vecthouf)
  gc()
  lonlatime <- data.frame(cbind(lonlatime,vectouf))
  metav<-data.frame(cbind(lonlatime[,c(1:3)],vecmeta))
  
  
  
  if(hazard==1)metav<-metav[which(!is.na(metav[,4])),]
  if(hazard==2)metav<-metav[which(!is.na(metav[,4])),]
  print(length(metav$Var1)) 
  head(metav)
  
  ep<-2
  
  # vary this coefficient !
  coef<-acof
  sampspd<-metav[,c(1,2,3,4)]
  sampspd$Var1<-sampspd$Var1*coef
  sampspd$Var2<-sampspd$Var2*coef
  sampspd$Var3<-sampspd$Var3-sampspd$Var3[1]+1
  sampspd$Var3<-sampspd$Var3
  
  
  ##############kNN plot for selection of eps################
  # Need to figure out how to completely smooth this shit
  wala=T
  if (wala==T){
    walabibou<-kNNdist(sampspd,k=10,all=F)
    min(walabibou)
    walabibou<-walabibou[which(walabibou<10)]
    walabibof=jitter(as.vector(walabibou),15)
    walabibord<-seq(1:length(walabibou))
    walabibix<-data.frame(as.vector(walabibof[order(walabibof)]))
    walabibix<-data.frame(as.vector(walabibou[order(walabibou)]))
    
    walabibix$order=walabibord
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
      
      
      plot(x_values,distances,type="l",lwd=3)
      points(x_max_dist,max(distances),col=2, cex=2,pch=16)
      
      return(c(x_max_dist, y_max_dist,xcmax,ycmax))
    }
    
    
    merd2<-elbow_finder(x_values,y_values,xcurve,ycurve)


    # plot(lacurve,type="l",lwd=4,col=3)
    # points(x_values,y_values,col=1,pch=16)
    # points(merd2[1],merd2[2],col=2,pch=16,cex=2)
    # points(merd2[3],merd2[4],col=2,pch=16,cex=2)
    # lines(c(merd2[1],merd2[3]),c(merd2[2],merd2[4]),lwd=4,lty=2,col=2)
    
    xtick=seq(0,1e5,by=2e4)
    par(mar=c(5,5,1,1))
    plot(walabilight[,2],walabilight[,1],ylim=c(0,10),ylab="10-NN distance",xlab="Points sorted by distance",cex.axis=1.5,cex.lab=1.8,xaxt="n")
    axis(side=1, at=xtick, labels = T,cex.axis=1.5)
    eps=round(merd2[2],2)
    txt<-bquote(epsilon == .(round(eps,2)))
    points(merd2[1],merd2[2],col=2,pch=16,cex=2)
    text(2e4,2.5,cex=1.8, labels=txt,pos=3,col="red")
    abline(h=merd2[2], col=2,lwd=3,lty=3)
    
    epsx=round(merd2[2],2)
  }
  print(epsx)
  # if (hazard==1){
  #   # eps=2.45
  #   eps=epsx
  # }
  # if (hazard==2)eps=epsx
  # sqrt(0.5^2+0.5^2+1^2)
  # epcl=eps
  # 
  
  
  ##Update on epcl for rain events
  if(hazard==1)epcl=epsx
  if(hazard==2)epcl=epsx
  
  #No weight
  #also need to vary minPts !
  rpip<-dbscan(sampspd, eps=epcl, minPts = mint)
  
  length(unique(rpip$cluster))
  
  
  spdata<-cbind(metav,rpip$cluster)
  noisepts<-length(which(spdata$`rpip$cluster`==0))
  if(length(which(spdata$`rpip$cluster`==0))>0){
    metav<-metav[-which(spdata[,5]==0),]
    spdata<-spdata[-which(spdata[,5]==0),]
  }
  metav$cluster<-spdata[,5]
  
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
  length(unique(event$cluster))
  
  testev<-aggregate(list(vr= event[,4]),
                    by = list(ev = event[,5],loc=event[,8]),
                    FUN = function(x) c(sum = sum(na.omit(x)),max=max(na.omit(x))))
  testev <- do.call(data.frame, testev)
  
  thl<-aggregate(list(vr= event[,4]),
                 by = list(ev = event[,5]),
                 FUN = function(x) c(l= length(x)))
  thl<- do.call(data.frame, thl)
  small<-thl$ev[which(thl$rf<10)]
  if (length(small)>1){
  bip<-which(!is.na(match(testev$ev,small)))
  testev<-testev[-bip,]
  bop<-which(!is.na(match(event$cluster,small)))
  event<-event[-bop,]
  length(unique(event$cluster))
  length(unique(testev$ev))
  }
  
  metamax<-aggregate(list(vr= testev[,3]) ,
                     by = list(ev = testev[,1]),
                     FUN = function(x) c(sum = sum(na.omit(x)),max=max(na.omit(x)),mean=mean(x),sd=sd(x),surf=length(x)))
  
  metamax<- do.call(data.frame, metamax)
  
  metave<-aggregate(list(vr= event[,4]) ,
                    by = list(ev = event[,5]),
                    FUN = function(x) c(surf=length(x)))
  
  
  metave <- do.call(data.frame, metave)
  length(unique(metave$ev))
  
  tempcom<-aggregate(event[,6] ,
                     by = list(ev = event[,5]),
                     FUN = function(x) c(dur=length(unique(x)),month=month(unique(x))[1],year=year(unique(x))[1]))
  
  tempic<- do.call(data.frame, tempcom)
  
  metave<-cbind(metamax,metave,tempic[,c(2,3,4)])
  maxV<-c()
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
    if (hazard==1)elV<-as.character(evint$loc[which(evint$vr.sum==met$vr.max)])[1]
    if (hazard==2)elV<-as.character(evint$loc[which(evint$vr.max==met$vr.max)])[1]
    maxV<-c(maxV,elV)
    evbk<-c(evbk,eves)
    
  }
  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
  Rcart=strsplit(maxV," ")
  Rcar<-as.data.frame(matrix(unlist(Rcart),ncol=2,byrow=T))
  Rcar[,1]<-as.numeric((Rcar[,1]))
  Rcar[,2]<-as.numeric((Rcar[,2]))
  
  metave$season=metave$x.month

  
  metave$startev<-as.POSIXct(startev,origin = "1970-01-01")
  metave$endev<-as.POSIXct(endev,origin = "1970-01-01")
  metaveN<-metave
  print(length(metamax$ev))
  
  metaHaz<-c(metaHaz,list(metaveN))
  metavHour<-c(metavHour,list(event))
  metavDaz<-c(metavDaz,list(testev))
}

## Step 3: Matchind hazard events

nwox<-list()
for(ha in c(1,2)){
  nwo<-vector(mode="list",length=length(metaHaz[[ha]]$ev))
for (c in 1: length(metaHaz[[ha]]$ev)){ 
  cl=metaHaz[[ha]]$ev[c]
  clev<-metavDaz[[ha]][which(metavDaz[[ha]]$ev==cl),]
  # dd<-unique(metavHax[[2]]$time[which(metavHax[[2]]$cluster==cl)])
  dd<-seq(metaHaz[[ha]]$startev[c],metaHaz[[ha]]$endev[c],by="hour")
  
  clex<-do.call(rbind, replicate(length(dd),clev[,c(1,2)], simplify=FALSE))
  cley<-do.call(rbind, replicate(length(clev[,1]),as.data.frame(dd),simplify=F))
  clex<-clex[order(clex$loc),]
  ouch<-cbind(clex,cley)
  # nwo<-rbind(nwo,ouch)
  nwo[[c]]<-ouch
  Sys.sleep(0.1)
  # update progress bar
}

nwo2<-ldply(nwo, data.frame)
rm(nwo)
nwox=c(nwox, list(nwo2))
}
gc()

sp03<-inner_join(nwox[[1]],nwox[[2]],by=c("loc","dd"))

tesp<-aggregate(list(len=sp03[,3]),
                by = list(ev2=sp03[,1],ev1 = sp03[,4]),
                FUN = function(x) c(length=length(unique(x))))
tesp <- do.call(data.frame, tesp)

outputs<-c(length(metaHaz[[1]]$ev),length(metaHaz[[2]]$ev),length(tesp$ev2))

outix<-rbind(outix,outputs)
}
outish<-cbind(parset,outix)
outish<-as.data.frame(outish)
names(outish)[c(4,5,6)]=c("rainc","windc","compc")

save(outish,file="sensitivity_out.Rdata")





########################################################
library(GGally)
kik<-cor(outish)

library(sensitivity)
srcplot=list()
for (id in 4:6){
merd<-src(X=outish[,c(1:3)], y=outish[,id], rank = F, logistic = FALSE, nboot = 200, conf = 0.95)
ggplot(merd)

srcplot<-c(srcplot, list(merd$SRC))
print(merd)
}
src1<-cbind(srcplot[[1]]$original,srcplot[[1]]$`min. c.i.`,srcplot[[1]]$`max. c.i.`)
src1<-data.frame(abs(src1))

src2<-cbind(srcplot[[2]]$original,srcplot[[2]]$`min. c.i.`,srcplot[[2]]$`max. c.i.`)
src2<-data.frame(abs(src2))

srcm<-c(srcplot[[1]]$original,srcplot[[2]]$original,srcplot[[3]]$original)
srcm<-data.frame(abs(srcm))
srcm$var=c(rep("Nr",3),rep("Nw",3),rep("Nc",3))
srcm$idv=c(rep(1,3),rep(2,3),rep(3,3))
srcm$par=c("r","mu","u","r","mu","u","r","mu","u")
src1<-srcm[which(srcm$par=="r"),]
src2<-srcm[which(srcm$par=="mu"),]
src3<-srcm[which(srcm$par=="u"),]
require(plotrix)
plotCI(src1$X1, ui=src1$X2, li=src1$X3,col="red",ylim=c(0,1))
plotCI(src2$X1, ui=src2$X2, li=src2$X3,col="blue",add=T)


par(mar=c(4.1, 4.1, 2.1, 10.1), xpd=TRUE)
plot(src1$idv,src1$abs.srcm.,col=2,type="o",ylim=c(0,1),pch=0,cex.axis=1.5,cex.lab=1.5,
     lwd=2,cex=1.5,ylab="abs(SRC)",xaxt="n",xlab="Output")
points(src2$idv,src2$abs.srcm.,col=3,type="o",pch=1,lwd=2,lty=2,cex=1.5)
points(src3$idv,src3$abs.srcm.,col=4,type="o",pch=2,lwd=2,lty=3,cex=1.5)
Axis(side=1, labels=c("Nr","Nw","Nc"),at=c(1,2,3),cex.axis=1.5)
legend(3.2,0.8,legend=c("r","mu","u"),col=c(2,3,4),pch=c(0,1,2),lwd=2,lty=c(1,2,3),cex=1.5)


require(ggplot2)
ggplot(src1, aes(x = X1, y = )) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L))



model <- lm(compc ~ pco + mco + aco, data = outish)
summary(model)

library(MASS)

# Stepwise regression model
step.model <- stepAIC(model, direction = "both", 
                      trace = FALSE)
summary(step.model)
