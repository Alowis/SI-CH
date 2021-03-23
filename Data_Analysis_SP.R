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

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
load(file="out/RainEv_metap_1979-2019.Rdata")
metaHar=metaHaf
load(file="out/WindEv_meta_1979-2019.Rdata")
boxplot(log(metaHar$vir.surf),log(metaHaf$viw.surf))
load(file="out/CompoundRW_79-19.v2x.Rdata")
load(file="out/messycompound_79-19.Rdata")

plot(metaHar$x.dur,metaHar$metave...c.2../metaHar$x.dur)
plot(metaHaf$x.dur,metaHaf$`metave[, c(2)]`/ metaHaf$x.dur,log="xy")
hist(log(metaHaf$`metave[, c(2)]`/ metaHaf$x.dur))
hist(log(metaHar$metave...c.2../metaHar$x.dur))
load(file="out/extremEvents_Rain.Rdata")
sptdfr<-sptdft
load(file="out/extremEvents_Wind.Rdata")

metaHar$spm=log(metaHar$x.dur*metaHar$vir.surf)
metaHaf$spm=log(metaHaf$x.dur*metaHaf$viw.surf)
compound$spm=log(compound$spacescale*compound$ORtime)
hist(metaHar$spm)
hist(metaHaf$spm)
hist(compound$spm,breaks=20)

load(file="out/Wind_stfprint2.Rdata")
load(file="out/Rain_stfprint.Rdata")
sp03<-inner_join(nwo3,nvo3,by=c("loc","dd"))

#######TESTS#########################
compoundl<-compound[which(compound$footprint<.01),]
compoundR<-compound[which(compound$footprint>=.01 & compound$footprint<.5),]
compoundM<-compound[which(compound$footprint>=.5),]

ggplot(data=metaHaf, aes(x=x.dur, group=con, fill=con)) +
  geom_density(adjust=5, alpha=.4,position="fill") +
  scale_fill_viridis()

plot(compoundl$wg.max.cs,compoundl$rf.max.cs,xlim=c(15,50),ylim=c(0,140))
points(compoundM$wg.max.cs,compoundM$rf.max.cs,col=2)
points(compoundR$wg.max.cs,compoundR$rf.max.cs,col=3)
orspace<-compoundfinalST$vir.surf+compoundfinalST$wg.space-compoundfinalST$spacescale.x
compoundfinalST$ORspace=orspace

save(compoundfinalST,file ="out/CompoundRW_79-19.v3x.Rdata")
plot(orspace,compoundfinalST$spacescale.x)
max(orspace)
rati<-data.frame(compoundfinalST$timescale.x/compoundfinalST$ORtime,compoundfinalST$spacescale.x/orspace)

plot(orspace,compoundfinalST$spacescale.x)
abline(a=0,b=1)
names(rati)=c("ratiotime","ratiospace")
plot(compoundfinalST$spacescale.x/orspace,rati$ratiospace)
hist(rati$ratiospace,breaks=20)

compoundfinalST$ratitime=rati[,1]
compoundfinalST$ratispace=rati[,2]

ratix<-aggregate(list(sp=compoundfinalST[,55],ti=compoundfinalST[,54],ar=compoundfinalST[,51]),
                by = list(mo=compoundfinalST[,30],y=compoundfinalST[,31]),
                FUN = function(x) c(m=mean(unique(x)))) 

ratix<- do.call(data.frame, ratix)

plot(ratix$ti,type="h")
plot(ratix$sp,type="h")
plot(ratix$ar,type="h")

plot(ratix$sp,ratix$ti)

ggplot(ratix, aes(x=factor(mo), y=ar)) +
  scale_fill_manual(values = c(  "darkorange","dodgerblue","forestgreen"),labels = c("Wind", "Rain", "Compound"))+
  scale_color_manual(values = c(  "chocolate4","dodgerblue4","darkgreen"),labels = c("Wind", "Rain", "Compound"))+
  stat_boxplot(geom = "errorbar")+
  geom_boxplot(position="dodge2",size=1) +
  theme(axis.text.x = element_text(angle=0, vjust = 0.5,size=16))+
  theme(axis.text=element_text(size=16),
        panel.grid.major = element_line(colour = "lightgrey") ,
        panel.border = element_rect(size=1.2, fill = NA),
        
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  theme(text = element_text(size=16))+
  labs(y = "# of occurences",x = "Month",colour = "Parameter",fill= "Hazard events",size=20)+
  scale_x_discrete(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c("Jan","Fev","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))+
  scale_y_continuous(name = "monthly occurences")

wala<-aggregate(list(len=metaHar[,1]),
                by = list(mo=metaHar[,9],y=metaHar[,10]),
                FUN = function(x) c(length=length(unique(x)))) 
wala<- do.call(data.frame, wala)
hist(wala$len,breaks=50)

wala<-aggregate(list(len=metaHaf[,1]),
                by = list(mo=metaHaf[,9],y=metaHaf[,10]),
                FUN = function(x) c(length=length(unique(x)))) 
wala<- do.call(data.frame, wala)
hist(wala$len,breaks=50)
difrr<-c()
for(i in 1:(length(metaHaf$ev)-1)){
  
difrr[i]<-difftime(metaHaf$startev[i+1],metaHaf$endev[i],unit="hours")

}
difrr

plot(difrr)
mean(difrr)
hist(difrr,breaks=130)
compound$sgroup=1
compound$month=month(compound$startime)
compound$sgroup[which(compound$month<10 & compound$month>3)]=2

comp1<-compound[which(compound$sgroup==1),]
comp2<-compound[which(compound$sgroup==2),]

plot(comp2$wg.max.cs,comp2$rf.max.cs)
plot(comp1$wg.max.cs,comp1$rf.max.cs)

chit<-chi(data.frame(comp1$wg.max.cs,comp1$rf.max.cs))
plot(chit$quantile, chit$chi[,2])


###############NEW SEction###########################

hex<-aggregate(list(len=sptdfr$Var1),
               by = list(tstep<-sptdfr$Var3),
               FUN = function(x) c(length=length(x)))
hexr <- do.call(data.frame, hex)
library(lubridate)
hexr$time<-as_datetime(c(unique(sptdfr$Var3*60*60)),origin="1900-01-01")
days<-as.Date(as_datetime(c(sptdfr$Var3*60*60),origin="1900-01-01"),format="%Y-%m")
mo<-format(as.Date(days),"%Y-%m")
mon<-month(days)

sptdfr$loc=paste(sptdfr$Var1,sptdfr$Var2,sep=" ")
length(sptdfr$Var1[which(sptdfr$`rpip$cluster`==0)])/length(sptdfr$Var1)
hex2<-aggregate(list(len=sptdfr$loc),
                by = list(tstep<-mo),
                FUN = function(x) c(length=length(x)))
hex2 <- do.call(data.frame, hex2)

hex3<-aggregate(list(len=sampspd$Var3),
                by = list(spstep<-sampspd$loc),
                FUN = function(x) c(length=length(unique(x))))
hexw3 <- do.call(data.frame, hex3)

hexx<-fill.gap(hexr)
plot(hex2$len,type="o",ylim=c(0,50000))
plot(hexx$time, hexx$len,type="h")
plot(sort(hex2$len),log="y")

acf (hexx$len, lag = 48)
abline(h=quantile(hexr$len,.95))
hexx$len[which(is.na(hexx$len))]=0
exi<-extremalindex(hexx$len, threshold=100, method = "intervals")
exi

length(hexr$Group.1)/length(hexx$time)
rbPal <- colorRampPalette(c('darkgreen',"gold","darkorange",'red'))

mapUK = SpacializedMap(regions = c("UK","France","Spain","Portugal","Italy","Ireland"))


ukk<- gSimplify(mapUK, tol=0.001, topologyPreserve=TRUE)
plot(ukk)
uk_fort <- ggplot2::fortify(ukk)
longlims=c(-6, 2)
longdom=seq(longlims[1],longlims[2],by=.25)
latlims=c(48,59)
latdom=seq(latlims[1],latlims[2],by=.25)

lagrid<-expand.grid(longdom,latdom)
lagrid2<-SpatialEpi::latlong2grid(lagrid)
lagrid2$x[1]-lagrid2$x[2]

logl<-c(-6,-5.75,1.75,2)
lagl<-c(48,48.25,58.75,59)
garea<-lagrid[which(!is.na(match(lagrid$Var1,logl)) | !is.na(match(lagrid$Var2,lagl))),]
lown<-c(lagrid[1,1]-0.125,lagrid[1,1]-0.125,lagrid[1,1]+0.125,lagrid[1,1]+0.125)
lawn<-c(lagrid[1,2]-0.125,lagrid[1,2]+0.125,lagrid[1,2]+0.125,lagrid[1,2]-0.125)

cell1<-data.frame(lown,lawn)
areaPolygon(cell1)/1000000
distm(lagrid[1,], lagrid[2,], fun = distGeo)

rgb.palette=colorRampPalette(c("white", 
                               "royalblue","orange","red","purple"),interpolate="linear",bias=1)
rbPal2<-colorRampPalette(c('white',"blue","purple","red","orange"))


Rcart=strsplit(compound$maxlocsR," ")
Rcar<-as.data.frame(matrix(unlist(Rcart),ncol=2,byrow=T))
Rcar[,1]<-as.numeric((Rcar[,1]))
Rcar[,2]<-as.numeric((Rcar[,2]))

Rcort=strsplit(compound$maxlocsW[which(compound$sgroup==1)]," ")
Rcor<-as.data.frame(matrix(unlist(Rcort),ncol=2,byrow=T))
Rcor[,1]<-as.numeric((Rcor[,1]))
Rcor[,2]<-as.numeric((Rcor[,2]))

############################Event visualiser###########################
load(file="out/RainEv_hdat_1979-2019.Rdata")
load(file="out/RainEv_ldat_1979-2019.Rdata")
metavDar=metavDaf
metavHar=metavHaf
load(file="out/WindEv_hdat_1979-2019.Rdata")
load(file="out/WindEv_ldat_1979-2019.Rdata")
load(file="out/Rainev_hdatp_1979-2019.Rdata")

load(file="out/CompoundRW_79-19.v2x.Rdata")
load(file="out/messycompound_79-19.Rdata")
load(file="out/CompoundRW_79-19.RP.Rdata")

roups<-quantile(compound$rf.max,.95) 
compmer<-compound[which(compound$rf.max>roups),]
plot(compmer$wg.max.cs,compmer$rf.max.cs)
plot(compmer$footprint,compmer$ORtime)
plot(compound$rf.max[order(compound$rf.max)])
plot(metatest$len.sum[order(metatest$len.sum)])
length(unique(metatix$ev1))
#Quick event visualizer
mapUK = SpacializedMap(database="world",regions = c("UK","France","Spain","Portugal","Italy","Ireland","Belgium","Netherland"))

library(nnet)
#visualizing storm Angus
rev<-stplot$ev1[2]

wev<-stplot$ev2[2]

evp<-metaHar$ev
eventa<-metavHar[which(metavHar$cluster==rev),]

com8<-compound[which(compoundfinalST$ev1==rev),]
com9<-com8[c(16,21)]
comr=strsplit(as.character(com9$maxlocsR)," ")
comRri<-as.data.frame(matrix(unlist(comr),ncol=2,byrow=T))
com9$lonR<-as.numeric((comRri[,1]))
com9$latR<-as.numeric((comRri[,2]))

comw=strsplit(as.character(com9$maxlocsW)," ")
comWi<-as.data.frame(matrix(unlist(comw),ncol=2,byrow=T))
com9$lonW<-as.numeric((comWi[,1]))
com9$latW<-as.numeric((comWi[,2]))

ppr<-aggregate(list(ri=eventa[,4]),
               by = list(loc=eventa[,8]),
               FUN = function(x) c(length=length(x)))
ppr <- do.call(data.frame, ppr)

evento<-metatest[which(metatest$ev1==rev),]  
evento$len=ppr$ri
Rr=strsplit(as.character(evento$loc)," ")
Rri<-as.data.frame(matrix(unlist(Rr),ncol=2,byrow=T))
evento$lon<-as.numeric((Rri[,1]))
evento$lat<-as.numeric((Rri[,2]))
  
cps<-metaHaf$ev
eventi<-metavHaf[which(metavHaf$cluster==wev),]

ppw<-aggregate(list(ri=eventi[,4]),
               by = list(loc=eventi[,8]),
               FUN = function(x) c(length=length(x),max=max(x)))
ppw <- do.call(data.frame, ppw)
Rw=strsplit(as.character(ppw$loc)," ")
Rwi<-as.data.frame(matrix(unlist(Rw),ncol=2,byrow=T))
ppw$lon<-as.numeric((Rwi[,1]))
ppw$lat<-as.numeric((Rwi[,2]))



# hollyshit<-inner_join(evento, maxeventwx,by="cloc")
# 
# txx<-aggregate(list(len=hollyshit[,3]),
#                 by = list(loc=hollyshit[,8]),
#                 FUN = function(x) c(length=length((x))))
# txx <- do.call(data.frame, txx)
# plot(txx$len,rank(txx$len))
# sp03$loc=as.character(sp03$loc)
# loco<-matrix(as.numeric(unlist(strsplit(sp03$loc," "))),ncol=2,byrow = T)
# sp03<-data.frame(sp03,loco)
# 
oula<-sp03[which(sp03$ev.y== rev & sp03$ev.x==wev),]
oulai<-unique(oula$loc)


Ro=strsplit(as.character(oulai)," ")
Roi<-as.data.frame(matrix(unlist(Ro),ncol=2,byrow=T))
oulai<-as.data.frame(oulai)
oulai$lon<-as.numeric((Roi[,1]))
oulai$lat<-as.numeric((Roi[,2]))
oulai$gr=1
cof<-c()
com<-c()
for (lat in unique(oulai$lat)){
elpmin<-min(oulai$lon[which(oulai$lat==lat)])
elpmax<-max(oulai$lon[which(oulai$lat==lat)])
coo<-c(elpmin,lat)
coa<-c(elpmax,lat)
cof<-rbind(cof,coo)
com<-rbind(com,coa)
}

cof[,1]
cof<-data.frame(cof)
cof$gr=1
com<-data.frame(com)
com$gr=1
cop<-c(-2,48)
cok<-c(-6,48)
comic<-rbind(cop,com[order(com[,1]),],cof[order(cof[,1],decreasing = T),],cok)



##############################Plots################### 

ukk<- gSimplify(mapUK, tol=0.001, topologyPreserve=TRUE)
uk_fort <- ggplot2::fortify(ukk)
longlims=c(-6, 2)
longdom=seq(longlims[1],longlims[2],by=.25)
latlims=c(48,59)
latdom=seq(latlims[1],latlims[2],by=.25)
lagrid<-expand.grid(longdom,latdom)
logl<-c(-6,-5.75,1.75,2)
lagl<-c(48,48.25,58.75,59)
garea<-lagrid[which(!is.na(match(lagrid$Var1,logl)) | !is.na(match(lagrid$Var2,lagl))),]
length(garea[,1])/length(lagrid[,1])

solorain<-evento[-match(oulai$oulai,evento$loc),]
solowind<-ppw[-match(oulai$oulai,ppw$loc),]
library(ggalt)
library(ggnewscale)

bs.palette=colorRampPalette(c("skyblue","dodgerblue","blue","darkblue"),interpolate="spline",bias=1)
bs.palette2=colorRampPalette(c("lightgreen","yellowgreen","gold","orangered","deeppink"),interpolate="spline",bias=1.5)

longlims=c(-5.7,1.7)
latlims=c(48.4,58.5)


ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill = "lemonchiffon1", color = "gray10", size = 1,alpha=.9) +
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
  # geom_tile(data=oulai,aes(x=lon,y=lat,group=gr),alpha=.7,fill="grey")+
  # geom_encircle(aes(x=Var1, y=Var2,group=cluster),
                # data = evento, size = 2, color = "blue",s_shape=.8, expand=0,fill="blue",alpha=.1)+

  geom_point(data=evento,aes(x=lon,y=lat,size=len,col=len.sum, group=loc),alpha=.9)+
  scale_color_gradientn(trans=scales::modulus_trans(1),colors=bs.palette(100),"ra [mm]",limits=c(1,80))+
  new_scale_color() +

  geom_point(data=ppw,aes(x=lon,y=lat,size=ri.length,col=ri.max, group=loc),shape=18,alpha=.8)+
  # geom_polygon(data=comic,aes(x=X1,y=X2,group=NA),alpha=1,size=1.5,fill="transparent",color="darkgreen")+

  # geom_bin2d(data=maxeventwx,aes(x=Var1,y=Var2,group=cluster),alpha=.5,binwidth=c(.25,.25))+
  # stat_density_2d(data=hollyshit,geom = "polygon", aes(x=Var1.x,y=Var2.x,group=cluster.x,fill=..level..),alpha=.3,n=c(100,45))+
  scale_color_gradientn(trans=scales::modulus_trans(.5),colors=bs.palette2(100),"w a [m/s]",limits=c(17,50))+
  
  scale_size(trans=scales::modulus_trans(1),limits=c(1,30),range=c(.5,5),"Duration [h]",
             guide = guide_legend(override.aes = list(colour = alpha("grey",.8))))+
  # 
  # geom_encircle(aes(x=Var1, y=Var2,group=cluster),
  #               data = maxeventwx, size = 2, color = "orange",s_shape=.8, expand=0,fill="orange",alpha=.5)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="italic"),
        panel.background = element_rect(fill = "transparent", colour = "grey50"),
        legend.title = element_text(size=18),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_y_continuous(
    breaks = c(48,50,52,54,56,58,60),limits = c(40,70),"Latitude")+
  scale_x_continuous(
    breaks =c(-6,-4,-2,0,2),limits=c(-10,10),"Longitiude") 



ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill = "lemonchiffon1", color = "gray10", size = 1,alpha=.9) +
  coord_map(xlim = longlims,  ylim = latlims, proj="albers",lat0=(50),lat1=(55))+
  # coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
  # geom_tile(data=oulai,aes(x=lon,y=lat,group=gr),alpha=.7,fill="grey")+
  # geom_encircle(aes(x=Var1, y=Var2,group=cluster),
  # data = evento, size = 2, color = "blue",s_shape=.8, expand=0,fill="blue",alpha=.1)+
  
 
  scale_color_gradientn(trans=scales::modulus_trans(1),colors=bs.palette(100),"ra [mm]",limits=c(1,50))+
  new_scale_color() +
  
  geom_point(data=ppw,aes(x=lon,y=lat,size=ri.length,col=ri.max, group=loc),shape=18,alpha=.8)+
  # geom_polygon(data=comic,aes(x=X1,y=X2,group=NA),alpha=1,size=1.5,fill="transparent",color="darkgreen")+
  
  # geom_bin2d(data=maxeventwx,aes(x=Var1,y=Var2,group=cluster),alpha=.5,binwidth=c(.25,.25))+
  # stat_density_2d(data=hollyshit,geom = "polygon", aes(x=Var1.x,y=Var2.x,group=cluster.x,fill=..level..),alpha=.3,n=c(100,45))+
  scale_color_gradientn(trans=scales::modulus_trans(.5),colors=bs.palette2(100),expression(paste("wa [m s"^"-1","]")),limits=c(17,40))+
  
  scale_size(trans=scales::modulus_trans(1),limits=c(1,30),range=c(.5,5),"Duration [h]",
             guide = guide_legend(override.aes = list(colour = alpha("grey",.8))))+
  # 
  # geom_encircle(aes(x=Var1, y=Var2,group=cluster),
  #               data = maxeventwx, size = 2, color = "orange",s_shape=.8, expand=0,fill="orange",alpha=.5)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="italic"),
        panel.background = element_rect(fill = "transparent", colour = "grey50"),
        legend.title = element_text(size=18),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_y_continuous(
    breaks = c(48,50,52,54,56,58,60),limits = c(40,70),"Latitude")+
  scale_x_continuous(
    breaks =c(-6,-4,-2,0,2),limits=c(-10,10),"Longitiude")





#################################################################


ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill = "lemonchiffon1", color = "gray10", size = 1,alpha=.9) +
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
  geom_raster(data=oulai,aes(x=lon,y=lat,group=gr),alpha=.7,fill="chartreuse4")+
  geom_raster(data=solorain,aes(x=lon,y=lat,group=ev1),alpha=.7,fill="blue")+
  geom_raster(data=solowind,aes(x=lon,y=lat,group=loc),alpha=.7,fill="orange")+
  # geom_encircle(aes(x=Var1, y=Var2,group=cluster),
  # data = evento, size = 2, color = "blue",s_shape=.8, expand=0,fill="blue",alpha=.1)+
  
  # geom_point(data=evento,aes(x=lon,y=lat,size=len,col=len.sum, group=loc),alpha=.9)+
  geom_point(data=com9,aes(x=lonR,y=latR,size=1, group=maxlocsW),colour="darkblue",alpha=1, size=4)+
  geom_point(data=com9,aes(x=lonW,y=latW,size=1, group=maxlocsW),colour="darkred",fill="darkred",alpha=1, size=4,shape=23)+
  scale_color_gradientn(trans=scales::modulus_trans(1),colors=bs.palette(100),"ra [mm]",limits=c(1,50))+
  new_scale_color() +
  
  # geom_point(data=ppw,aes(x=lon,y=lat,size=ri.length,col=ri.max, group=loc),shape=18,alpha=.8)+
  # geom_polygon(data=comic,aes(x=X1,y=X2,group=NA),alpha=1,size=1.5,fill="transparent",color="darkgreen")+
  
  # geom_bin2d(data=maxeventwx,aes(x=Var1,y=Var2,group=cluster),alpha=.5,binwidth=c(.25,.25))+
  # stat_density_2d(data=hollyshit,geom = "polygon", aes(x=Var1.x,y=Var2.x,group=cluster.x,fill=..level..),alpha=.3,n=c(100,45))+
  scale_color_gradientn(trans=scales::modulus_trans(.5),colors=bs.palette2(100),expression(paste("wa [m s"^"-1","]")),limits=c(17,40))+
  
  scale_size(trans=scales::modulus_trans(1),limits=c(1,30),range=c(.5,5),"Duration [h]",
             guide = guide_legend(override.aes = list(colour = alpha("grey",.8))))+
  # 
  # geom_encircle(aes(x=Var1, y=Var2,group=cluster),
  #               data = maxeventwx, size = 2, color = "orange",s_shape=.8, expand=0,fill="orange",alpha=.5)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="italic"),
        panel.background = element_rect(fill = "transparent", colour = "grey50"),
        legend.title = element_text(size=18),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_y_continuous(
    breaks = c(48,50,52,54,56,58,60),limits = c(40,70),"Latitude")+
  scale_x_continuous(
    breaks =c(-6,-4,-2,0,2),limits=c(-10,10),"Longitiude") 




Rcar$V3=jitter(Rcar$V1,factor=3)
Rcar$V4=jitter(Rcar$V2,factor=3)
Rcor$V3=jitter(Rcor$V1,factor=5)
Rcor$V4=jitter(Rcor$V2,factor=5)
library(ggpointdensity)

#####################Hotspot locations###########################

ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill = "white", color = "gray10", size = 1) +
  theme_bw(16)+
  scale_fill_gradientn(colours = rgb.palette(100))+
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.2)+
  # geom_bin2d(data=metaHaz[[2]],aes(x=lonMW,y=latMW,group=season),alpha=.7,binwidth=c(.5,.5))+
  # geom_point(data=Rcar,aes(x=V3,y=V4,group=1),alpha=.3)+
  # geom_point(data=Rcor,aes(x=V3,y=V4,group=1,color="red"),alpha=.3)
geom_pointdensity(size=2,adjust = 1,data=Rcor,aes(x=V3,y=V4,group=1),shape=16,method="kde2d") + scale_color_viridis_c(option="magma")
  stat_density_2d(data=Rcor,aes(x=V1,y=V2,group=1,fill = stat(density)),geom = "raster", contour = FALSE,alpha=.6,n=300)

  
  distus<-dist(data.frame(Rcor[,c(1,2)]))
  hc<-hclust(distus,method="ward.D2")
  plot(hc)

  hcx<-hclust(distus^2,method="centroid")
  plot(hcx)
  memb <- cutree(hc, k = 4)
  ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
    geom_polygon(fill = "white", color = "gray10", size = 1) +
    theme_bw(16)+
    scale_fill_gradientn(colours = rgb.palette(100))+
    coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.2)+
    geom_point(data=Rcor,aes(x=V3,y=V4,group=1,color=factor(memb)),alpha=.4,size=3)
  
  hc2 <- agnes(distus, method = "complete")
  hc2$ac
  plot(hc2)
  library(factoextra)
  fviz_nbclust(as.matrix(distus), FUN = hcut, method = "wss")
  
  reservo <- optics(Rcar[,c(3,4)],eps=10, minPts = 100)
  res <- extractXi(reservo, xi = 0.)
  res=extractDBSCAN(reservo, eps_cl =.7)
  plot(res)
  reach <- as.reachability(res)
  dend <- as.dendrogram(reach)
  dend
  plot(dend)
  
  res <- extractXi(res, xi = .07)
  hullplot(Rcar[,c(1,2)], res)
  
Rcast=strsplit(compound$maxlocsW," ")
Rcas<-as.data.frame(matrix(unlist(Rcast),ncol=2,byrow=T))
Rcas[,1]<-as.numeric.factor((Rcas[,1]))
Rcas[,2]<-as.numeric.factor((Rcas[,2]))

ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill = "white", color = "gray10", size = 1) +
  theme_bw(16)+
  scale_fill_gradientn(colours = rgb.palette(100))+
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.2)+
  # geom_bin2d(data=metaHaz[[2]],aes(x=lonMW,y=latMW,group=season),alpha=.7,binwidth=c(.5,.5))+
  # geom_density2d(data=metaHaz[[2]],aes(x=lonMW,y=latMW,group=season),alpha=.2)
  stat_density_2d(data=Rcas,aes(x=V1,y=V2,group=1,fill = stat(density)),geom = "raster", contour = FALSE,alpha=.6,n=300)

# 


#####################traditional hotspot################

load("out/CompoundRW_79-19.ME2.Rdata")
sp03<-inner_join(nwo3,nvo3,by=c("loc","dd"))
sp03$combin=paste(sp03$ev.y,sp03$ev.x)
kaka<-which(!is.na(match(sp03$combin,compoundX$combin)))
length(kaka)/length(sp03$combin)
sp04<-sp03[kaka,]
focus="major"
if(focus=="major"){
  sp03=sp04
}
tesp<-aggregate(list(len=sp03[,3]),
                by = list(loc=sp03[,2]),
                FUN = function(x) c(length=length(unique(x))))
tesp <- do.call(data.frame, tesp)

aggregR<-aggregate(list(len=nvo3[,3]),
                   by = list(loc=nvo3[,2]),
                   FUN = function(x) c(length=length(unique(x))))
aggregR <- do.call(data.frame, aggregR) 

aggregW<-aggregate(list(len=nwo3[,3]),
                   by = list(loc=nwo3[,2]),
                   FUN = function(x) c(length=length(unique(x))))
aggregW <- do.call(data.frame, aggregW) 


sp03$loc=as.character(sp03$loc)
loco<-matrix(as.numeric(unlist(strsplit(sp03$loc," "))),ncol=2,byrow = T) 
sp03<-data.frame(sp03,loco) 
allprec<-sp03 

allprec$month=month(allprec$dd) 
allprec$year=year(allprec$dd) 
allprec$season=1
allprec$season[which(allprec$month<6 & allprec$month>2)]=2
allprec$season[which(allprec$month<10 & allprec$month>5)]=3
allprec$season[which(allprec$month<12& allprec$month>9)]=4

allprec$sgroup=1
allprec$sgroup[which(allprec$month<10 & allprec$month>3)]=2

nvo3$sgroup=1
nvo3$sgroup[which(month(nvo3$dd)<10 & month(nvo3$dd)>3)]=2

nwo3$sgroup=1
nwo3$sgroup[which(month(nwo3$dd)<10 & month(nwo3$dd)>3)]=2

allprec$evxy=paste(allprec$ev.x,allprec$ev.y,sep=" ")
a=1

tesp<-aggregate(list(len=allprec[which(allprec$sgroup==a),3]),
                by = list(loc=allprec[which(allprec$sgroup==a),2]),
                FUN = function(x) c(length=length(unique(x))))
tesp <- do.call(data.frame, tesp)

aggregR<-aggregate(list(len=nvo3[which(nvo3$sgroup==a),3]),
                   by = list(loc=nvo3[which(nvo3$sgroup==a),2]),
                   FUN = function(x) c(length=length(unique(x))))
aggregR <- do.call(data.frame, aggregR) 

aggregW<-aggregate(list(len=nwo3[which(nwo3$sgroup==a),3]),
                   by = list(loc=nwo3[which(nwo3$sgroup==a),2]),
                   FUN = function(x) c(length=length(unique(x))))
aggregW <- do.call(data.frame, aggregW) 




windval<-inner_join(allprec, meta)

length(allprec$sgroup[which(allprec$sgroup==a)])/length(allprec$sgroup)

omshit<-aggregate(allprec[which(allprec$sgroup==a),1] ,
                         by = list(loc=allprec[which(allprec$sgroup==a),2],y=allprec$year[which(allprec$sgroup==a)]),
                         FUN = function(x) c(n =length(unique(x)),l=length(x)))
omshit <- do.call(data.frame, omshit)

###Need to comment all these... wtf

omgNew<-aggregate(allprec$evxy ,
                  by = list(loc=allprec[,2],m=allprec$month),
                  FUN = function(x) c(n =length(unique(x)),l=length(x)))
omgNew <- do.call(data.frame, omgNew)

omgNew2<-aggregate(list(a=omgNew$x.l) ,
                  by = list(loc=omgNew$loc),
                  FUN = function(x) c(n =max(x)))
omgNew2 <- do.call(data.frame, omgNew2)

ccr<-match(omgNew2$a,omgNew$x.l)
moismax=c()
for (xx in unique(omgNew2$loc)){
  laloc<-omgNew[which(omgNew$loc==xx),]
  lemois<-unique(laloc$m[which(laloc$x.l==max(laloc$x.l))])
  print(lemois)
  moismax<-c(moismax,lemois[1])
}

maxm<-data.frame(omgNew2,moismax)

locus<-strsplit(maxm$loc," ")
locus<-as.data.frame(matrix(unlist(locus),ncol=2,byrow=T))
maxm$lon<-as.numeric((locus[,1]))
maxm$lat<-as.numeric((locus[,2]))
maxm$grp=1
rgb.palette=colorRampPalette(c("lightblue", 
                               "royalblue","orange","red","purple"),interpolate="linear",bias=1)

ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill = "white", color = "gray10", size = 1) +
  theme_bw(16)+
  # scale_fill_gradientn(colours = rgb.palette(100),"Total hours in a CE",breaks=breaks_extended(5),limits=c(min(locav$x)-0.1*min(locav$x),max(locav$x)+.1*min(locav$x)))+
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3)+
  geom_raster(data=locav,aes(x=lon,y=lat,fill=tcom,group=loc),interpolate = F,alpha=.7)+
  scale_fill_gradientn(trans=scales::modulus_trans(1),colours = rgb.palette(100)," ",breaks=breaks_extended(5),limits=c(0,16))+
  # scale_fill_gradientn(trans=scales::modulus_trans(1),colours = rgb.palette(100),"LMF",breaks=breaks_extended(5),limits=c(0,.015))+
  labs(y = "Latitude",x = "Longitude",fill= "Total hours in a CE",size=20)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="italic"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(angle=0,size=18),
        legend.text = element_text(angle=0,size=14),
        legend.position = "right",
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_y_continuous(
    breaks = c(48,50,52,54,56,58,60),limits = c(40,70),"Latitude")+
  scale_x_continuous(
    breaks =c(-6,-4,-2,0,2),limits=c(-10,10),"Longitiude") 
guides(fill = guide_legend())




worstY<-aggregate(list(n=omshit$x.l) ,
                  by = list(loc=omshit$loc), 
                  FUN = function(x) c(n =length(unique(x)),m=max(x)))
worstY <- do.call(data.frame,worstY)

omrd<-match(paste(worstY$n.m,worstY$loc) ,paste(omshit$x.l,omshit$loc))
worstY$year<-omshit$y[omrd]
locav=worstY


locav<-aggregate(list(me=omshit[,4]) ,
                 by = list(loc=omshit[,1]),
                 FUN = function(x) c(n =sum(x),m=mean(x)))
locav <- do.call(data.frame, locav)


locav<-aggregate(allprec[,1] ,
                       by = list(loc=allprec[,2]),
                       FUN = function(x) c(me =length(x)))

########Locav for hotspots#####################

locav$tcom1=tesp$len/ttstep
locav$cdens1=locav$tcom1/max(locav$tcom1)
locav$densdif=locav$tcom/locav$tcom1

####work on this######

locav=tesp
locav$loc<-as.character(locav$loc)
locus<-strsplit(locav$loc," ")
locus<-as.data.frame(matrix(unlist(locus),ncol=2,byrow=T))
locav$x=locav$me.n
locav$lon<-as.numeric((locus[,1]))
locav$lat<-as.numeric((locus[,2]))
ttstep<-41
locav$p<-locav$x/ttstep
locav$grp<-1
locav$train<-aggregR$len/ttstep
locav$twin<-aggregW$len/ttstep
locav$tcom<-tesp$len/ttstep
locav$pcw<-locav$tcom/locav$twin
locav$pcr<-locav$tcom/locav$train
locav$pr<-locav$train/(24*365.25)
locav$pw<-locav$twin/(24*365.25)
locav$prc=locav$pcr*locav$pr
locav$pwc=locav$pcw*locav$pw
locav$pnaif<-locav$pr*locav$pw
locav$mult<-locav$pwc/locav$pnaif
max(locav$pwc)
longlims=c(-5.7,1.7)
latlims=c(48.4,58.5)
max(locav$mult)
locav$cdens=locav$tcom/max(locav$tcom)


rgb.palette=colorRampPalette(c("lightblue", 
                               "royalblue","orange","red","purple"),interpolate="linear",bias=1)

ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill = "white", color = "gray10", size = 1) +
  theme_bw(16)+
  # scale_fill_gradientn(colours = rgb.palette(100),"Total hours in a CE",breaks=breaks_extended(5),limits=c(min(locav$x)-0.1*min(locav$x),max(locav$x)+.1*min(locav$x)))+
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3)+
  # geom_raster(data=locav,aes(x=lon,y=lat,fill=tcom,group=loc),interpolate = F,alpha=.7)+
  geom_contour_fill(data=locav,aes(x=lon,y=lat,z = densdif,group=grp),alpha=0.5) +
  scale_fill_gradientn(trans=scales::modulus_trans(1),colours = rgb.palette(100)," ",breaks=breaks_extended(5),limits=c(0,0.6))+
  # scale_fill_gradientn(trans=scales::modulus_trans(1),colours = rgb.palette(100),"LMF",breaks=breaks_extended(5),limits=c(0,.015))+
  labs(y = "Latitude",x = "Longitude",fill= "Total hours in a CE",size=20)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="italic"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(angle=0,size=18),
        legend.text = element_text(angle=0,size=14),
        legend.position = "right",
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_y_continuous(
    breaks = c(48,50,52,54,56,58,60),limits = c(40,70),"Latitude")+
  scale_x_continuous(
    breaks =c(-6,-4,-2,0,2),limits=c(-10,10),"Longitiude") 
  guides(fill = guide_legend())

  # geom_bin2d(data=metaHaz[[2]],aes(x=lonMW,y=latMW,group=season),alpha=.7,binwidth=c(.5,.5))+
  # geom_density2d(data=metaHaz[[2]],aes(x=lonMW,y=latMW,group=season),alpha=.2)
  # stat_density_2d(data=allprec,aes(x=X1,y=X2,group=season,fill = stat(density)),geom = "raster", contour = FALSE,alpha=.6,n=150)

#################################SOMETHING Happens here##########################
# load(file="out/RainEv_hdat_1979-2019.Rdata")
load(file="out/RainEv_ldat_1979-2019.Rdata")
load(file="out/RainEv_hdat_1979-2019.Rdata")
# load(file="out/RainEv_metap_1979-2019.Rdata")
metavDar=metavDaf
metavHar=metavHaf
# load(file="out/WindEv_hdat_1979-2019.Rdata")
load(file="out/WindEv_ldat_1979-2019.Rdata")
load(file="out/WindEv_hdat_1979-2019.Rdata")

load(file="out/Rainev_hdatp_1979-2019.Rdata")

qx<-quantile(metatest$len.sum,.9)
metates2<-metatest[which(metatest$len.sum>=qx),]
qw<-quantile(metavHaf$vecmeta,.9)
metates3<-metavHaf[which(metavHaf$vecmeta>=qw),]
proprain<-c()
propwind<-c()
rmm<-c()
wmm<-c()
for(bx in 1:1485)
{
location=locav$loc[bx]
compbx<-sp03[which(sp03$loc==location),]
comcr<-unique(compbx$ev.y)
comcw<-unique(compbx$ev.x)
rainc<-metatest[which(metatest$loc==location),]
cx<-compound[which(!is.na(match(compound$ev1,comcr))),]
cxx<-cx[which(!is.na(match(cx$ev2,comcw))),]
cy<-rainc$len.sum[na.omit(match(cxx$ev1,rainc$ev))]
rmax<-cy/cxx$rf.max.cs
# plot(rmax)
rmm[bx]<-mean(rmax)
winc<-metavDaf[which(metavDaf$loc==location),]
cw<-winc$vi.max[na.omit(match(cxx$ev2,winc$ev))]

wmax<-cw/cxx$wg.max.cs
wmm[bx]<-mean(wmax)
proprain[bx]<-length(comcr)/length(rainc$ev)
propwind[bx]<-length(comcw)/length(winc$ev)
}

locav$prain<-proprain
locav$pwind=propwind
locav$rmm<-rmm
locav$wmm<-wmm


locav$fwind=locav$pwind*100
locav$frain=locav$prain*100
locav$frmm<-locav$rmm*100
locav$fwmm<-locav$wmm*100
locav$rwmm<-rmm*wmm*100
min(locav$fwmm)

ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill = "white", color = "gray10", size = 1) +
  theme_bw(16)+

  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.2)+
  geom_contour_fill(data=locav,aes(x=lon,y=lat,z = fwind,group=grp),binwidth=5,alpha=0.5) +
  scale_fill_gradientn(colours = rgb.palette(100),"Frequency [%]",breaks=breaks_extended(4),limits=c(0,60))+
  labs(y = "Latitude",x = "Longitude",fill= "Frequency [%]",size=20)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="italic"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(angle=0,size=14),
        legend.position = "right",
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_y_continuous(
    breaks = c(48,50,52,54,56,58,60),limits = c(40,70),"Latitude")+
  scale_x_continuous(
    breaks =c(-6,-4,-2,0,2),limits=c(-10,10),"Longitiude") 


max(metavHaf$vecmeta)
min(metavHaf$vecmeta)
compourain<-which(!is.na(match(metavHaf$cluster,compound$ev1)))
submetavHaf<-metavHaf[compourain,]
metest<-metavHaf[which(metavHaf$cluster==10),]
qcon<-quantile(metest$vecmeta,.8)

les20<-sum(metest$vecmeta[which(metest$vecmeta>qcon)])
lesall<-sum(metest$vecmeta)
xratio<-les20/lesall
plot(sort(metest$vecmeta,decreasing=T),ylim=c(0,20),xlim=c(1,10000),type="l",log="x")

les20<-c()
lesall<-c()
xr<-c()
qcon<-c()
for (idx in 1:2000){
# lines(sort(submetavHaf$vecmeta[which(submetavHaf$cluster==evx)],decreasing = T))
  evx<-metaHaf$ev[idx]
  metest<-metavDaf[which(metavDaf$ev==evx),]
  qcon[idx]<-max(metest$vi.max)
  # les20[idx]<-sum(metest$vecmeta[which(metest$vecmeta>qcon)])
  # lesall[idx]<-sum(metest$vecmeta)
}
wix<-which(metaHaf$ev==evx)
dates<-metaHaf$season[1:2000]

rat<-les20/lesall
boxplot(qcon~dates)
plot(dates,qcon)

max(compoundfinalST$ORtime)/24
24*4
max(compoundfinalST$spacescale.x)/1485*100
compound$month<-month(compound$startime)
compound$sgroup<-1
compound$sgroup[which(compound$month<10 & compound$month>3)]=2

a=2
rbPal2 <- colorRampPalette(c('red',"darkorange","gold",'darkgreen'))
Col2 <- rbPal2(1000)[as.numeric(cut(log(metaHaz[[1]]$Dmax+0.00001),breaks = 1000))] 
plot(compound$wg.max.cs[which(compound$sgroup==a)],compound$rf.max.cs[which(compound$sgroup==a)],pch=16)


compoundfinalST$sgroup=compound$sgroup
mean(log(compoundfinalST$spacescale.y[which(compoundfinalST$sgroup==2)]))

vart1<-log(compoundfinalST$ORtime[which(compoundfinalST$sgroup==1)])
vart2<-log(compoundfinalST$ORtime[which(compoundfinalST$sgroup==2)])

t.test(vart1,vart2)

boxplot(vart1 , vart2)
#################Polat plot#####################################
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


GeomSegment2 <- ggproto("GeomSegment2",
                        GeomSegment,
                        draw_panel = function (data, panel_params, coord, arrow = NULL,
                                               arrow.fill = NULL, lineend = "butt", 
                                               linejoin = "round", na.rm = FALSE) {
                          data <- remove_missing(data, na.rm = na.rm, 
                                                 c("x", "y", "xend", "yend", "linetype", 
                                                   "size", "shape"), 
                                                 name = "geom_segment")                          
                          if (ggplot2:::empty(data)) 
                            return(zeroGrob())
                          # remove option for linear coordinate system
                          data$group <- 1:nrow(data)
                          starts <- subset(data, select = c(-xend, -yend))
                          ends <- plyr::rename(subset(data, select = c(-x, -y)), 
                                               c(xend = "x", yend = "y"), 
                                               warn_missing = FALSE)
                          pieces <- rbind(starts, ends)
                          pieces <- pieces[order(pieces$group), ]
                          GeomPath$draw_panel(pieces, panel_params, coord, arrow = arrow, 
                                              lineend = lineend)
                        })

geom_segment2 <- function (mapping = NULL, data = NULL, stat = "identity", 
                           position = "identity", ..., arrow = NULL, arrow.fill = NULL, 
                           lineend = "butt", linejoin = "round", na.rm = FALSE, 
                           show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, 
        geom = GeomSegment2, # instead of GeomSegment
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(arrow = arrow, arrow.fill = arrow.fill, 
                      lineend = lineend, linejoin = linejoin, na.rm = na.rm, 
                      ...))
}
rbPal2<-colorRampPalette(c('skyblue',"orange","darkorange","red","purple"))
polarEventR$gr=1

data(beavers)
plot(polarEventR[,c(1,3)])
p <- ggplot() + 
  geom_path(data=polarEventR, aes(x=monthstr,y=x,group=gr,colour=Year),size=.8)+
  geom_line(data=polarEventR, aes(x=monthstr,y=x,colour = Year,group=Year),size=.8) + 
  geom_point(data=polarEventR, aes(x=monthstr,y=x,colour = Year,group=Year),size=2)+
  # geom_segment2(data=polarEventR,aes(x=monthstr, 
  #                                    xend=monthstr, 
  #                                    y=0, 
  #                                    yend=x),color="blue")+
  scale_colour_gradientn(colours = rbPal2(30))+
  scale_y_continuous(name = "# of occurences")+
  scale_x_discrete(name = "Months")
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

##############################Seasonal analysis#######################################

matchingr<-which(!is.na(match(nvo3$ev,compoundfinalST$ev1)))
matchingw<-which(!is.na(match(nwo3$ev,compoundfinalST$ev2)))

ttm<-aggregate(list(len=allprec$dd), 
                by = list(m=allprec$month,y=allprec$year),
                FUN = function(x) c(length=length(unique(x))))
ttm <- do.call(data.frame, ttm)
ttm$prob=ttm$len/(24*30) 
ggplot(ttm,aes(x= m,y=y)) +
  theme_bw(16)+
  geom_raster(aes(fill=prob),interpolate = F)+
  scale_fill_gradientn(trans=scales::modulus_trans(.6),colors=rgb.palette(100))

nvo3$month=month(nvo3$dd) 
nvo3$year=year(nvo3$dd) 

nwo3$month=month(nwo3$dd) 
nwo3$year=year(nwo3$dd) 

nvo4<-nvo3[-matchingr,]
aggretR<-aggregate(list(len=nvo3[,3]),
                   by = list(m=nvo3$month,y=nvo3$year),
                   FUN = function(x) c(length=length(unique(x))))
aggretR <- do.call(data.frame, aggretR) 

aggretR$prob=aggretR$len/(24*30) 

rgb.palette=colorRampPalette(c("lightblue", 
                               "royalblue","orange","red","purple"),interpolate="linear",bias=1)

ggplot(aggretW,aes(x= m,y=y)) +
  theme_bw(16)+
  geom_raster(aes(fill=prob),interpolate = F)+
  scale_fill_gradient(trans=scales::modulus_trans(.6),low="white",high="red")

aggretW<-aggregate(list(len=nwo3[,3]),
                   by = list(m=nwo3$month,y=nwo3$year),
                   FUN = function(x) c(length=length(unique(x))))
aggretW <- do.call(data.frame, aggretW) 

aggretW$prob=aggretW$len/(24*30)*100

aggreTenp<-full_join(aggretR,aggretW,by=c("m","y"))
aggreTemp<-full_join(aggreTenp,ttm,by=c("m","y"))
aggreTemp$prob[which(is.na(aggreTemp$prob))]=0
aggreTemp$prob.y[which(is.na(aggreTemp$prob.y))]=0
aggreTemp$prob=aggreTemp$prob*100
aggreTemp$pinde<-aggreTemp$prob.x*aggreTemp$prob.y
aggreTemp$lmf<-aggreTemp$prob/aggreTemp$pinde
aggreTemp$len[which(is.na(aggreTemp$len))]=0
aggreTemp$pi<-aggreTemp$len/vivix$x.n
aggreTemp$nb<-vivix$x.n
aggreTemp$dur=aggreTemp$len/aggreTemp$nb

ggplot(aggreTemp,aes(x= factor(m),y=y)) +
  theme_bw(16)+
  geom_raster(aes(fill=prob),interpolate = F,alpha=.9)+
  scale_fill_gradientn(trans=scales::modulus_trans(1),colors=rgb.palette(100),na.value="white",limits=c(0,60))+
  theme(axis.text.x = element_text(angle=0, vjust = 0.5,size=16))+
  theme(axis.text=element_text(size=16),
        panel.grid.major = element_line(colour = NA) ,
        panel.border = element_rect(size=1.2, fill = NA),
        
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  theme(text = element_text(size=16))+
  labs(y = "# of occurences",x = "Month",colour = "Parameter",fill= "Frequency [%]",size=20)+
  scale_x_discrete(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c("Jan","Fev","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))+
  scale_y_continuous(name = "Years")


wamr<-aggregate(aggreTemp$prob.x,
                by=list(aggreTemp$m),
                FUN = function(x) c(m=median(x)))
wamr<- do.call(data.frame, wamr)
wamr$c=2
wamw<-aggregate(aggreTemp$prob.y,
                by=list(aggreTemp$m),
                FUN = function(x) c(m=median(x)))
wamw<- do.call(data.frame, wamw)
wamw$c=1
wamc<-aggregate(aggreTemp$prob,
                by=list(aggreTemp$m),
                FUN = function(x) c(m=median(x)))
wamc<- do.call(data.frame, wamc)
wamc$c=3

wam=rbind(wamr,wamw,wamc)


ggplot(wam,aes(x= factor(Group.1),y=c)) +
  theme_bw(16)+
  geom_raster(aes(fill=x),interpolate = F,alpha=.9)+
  scale_fill_gradientn(trans=scales::modulus_trans(1),colors=rgb.palette(100))+
  theme(axis.text.x = element_text(angle=0, vjust = 0.5,size=16))+
  theme(axis.text=element_text(size=16),
        panel.grid.major = element_line(colour = NA) ,
        panel.border = element_rect(size=1.2, fill = NA),
        
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  theme(text = element_text(size=16))+
  labs(y = "# of occurences",x = "Month",colour = "Parameter",fill= "Frequency [%]",size=20)+
  scale_x_discrete(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c("Jan","Fev","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))+
  scale_y_continuous(name = "Years")

# extract rain only and wind only events

matchingr<-unique(match(compoundfinalST$ev1,metaHar$ev))
matchingw<-unique(match(compoundfinalST$ev2,metaHaf$ev))

rainevi<-metaHar[-matchingr,c(1,9,10,11,8,6)]
rainevi$c="r"
windevi<-metaHaf[-matchingw,c(1,9,10,11,8,6)]
windevi$c="w"

rainevf<-metaHar[-matchingr,]

hist(rainevf$vir.max)
rainotro<-metaHar[matchingr,]
hist(rainotro$vir.max)
rainevf$c="r"
rainotro$c="c"



###############Quantile-space plot#############################

quantile(metatest$len.sum,.5)
qt=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,.95,.99)
qlist=c(0,quantile(metatest$len.sum,.1),quantile(metatest$len.sum,.2),
        quantile(metatest$len.sum,.3),quantile(metatest$len.sum,.4),
        quantile(metatest$len.sum,.5),quantile(metatest$len.sum,.6),
        quantile(metatest$len.sum,.7),quantile(metatest$len.sum,.8),
        quantile(metatest$len.sum,.9),quantile(metatest$len.sum,.95),
        quantile(metatest$len.sum,.99))

bb<-metaHar$ev[which(metaHar$vir.max>10)]
metatest2<-metatest[which(!is.na(match(metatest$ev1,bb))),]
plot(qlist)

ccl=c()
clo=c()
for(evi in bb){
  metev<-metatest[which(metatest$ev1==evi),]
  lan=unique(metev$len.length)
  cl=c()
  for(q in 1:12){
    l=length(metev$len.sum[which(metev$len.sum>qlist[q])])
    cl=c(cl,l)
  }
  clo<-cbind(qt,rep(evi,11),cl,rep(lan,11))
  ccl<-rbind(ccl,clo)
}

rev=12593 #midlands floods 2007
wev=5423  #storm xaver 2013

cev=4172 #storm angus 4286

compound$id=c(1:length(compound$combin))
cevc=compound$combin[which(compound$id==cev)]


rltim=ultim[which(ultim$Dominant.hazard=="Extreme rainfall"),]
rainid=as.vector(rltim$ID_mod)
rix<-  strsplit(rainid, ",")
irc=c()
listev<-c()
ism=c()
comidr=c()
for(idr in 1:length(rix)){
  rmu=c()
  bof=rix[[idr]]
  bofv<-as.data.frame(matrix(unlist(bof),ncol=length(bof),byrow=T))
  
  if(length(bof)>1)rmu=rep(1,length(bof))
  if(length(bof)==1)rmu=2

  iscom=(match(as.numeric(bof),compound$ev1))
  icomid=na.omit(compound$combin[iscom])
  if(length(iscom)>1)rcom=T
  if(length(iscom)<1)rcom=F
  comidr=c(comidr,icomid)
  irc=c(irc,rcom)
  ism=c(ism,rmu)
  listev=c(listev,bofv)
}
length(irc[which(irc==T)])

rltim$iscompound=irc

listev=as.numeric(listev)
listevp=cbind(listev,ism)
listevp=as.data.frame(listevp)
names(listevp)=c("V2","linetype")

rlen=length(listevp$V2)/length(rltim$Year)

coui=as.data.frame(ccl)
coui$cl=(coui$cl/1485)*100
max(coui$V4)
cofl=coui[which(coui$V2==rev),]
cool=inner_join(coui,listevp,by="V2")

rgb.palette=colorRampPalette(c("#ffffcc","#ffeda0", "#fed976", "#feb24c","#fd8d3c","#fc4e2a",
  "#e31a1c","#bd0026","#800026"),interpolate="linear",bias=1)


ggplot(coui, aes(x = qt, y = cl, group=as.character(V2))) +
  geom_line(alpha=.5,color="grey",size=1.2)+
  
  scale_x_continuous(limits = c(0,1), expand = c(0, 0),"Intensity") +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0),"Spatial footprint [%]") +
  scale_size(range = c(0.3, 1.2)) +
  scale_linetype_identity() +
  
  geom_line(data=cool, aes(x = qt, y = cl, group=as.character(V2),color=V4),alpha=.7,size=1.2,linemitre=10) +
  geom_line(data=cofl,aes(x = qt, y = cl, group=as.character(V2)),color="darkblue",size=1.2,linemitre=10) +
  # scale_alpha_continuous(c(0.7, 0.8),guide = FALSE)+
  theme(axis.text=element_text(size=16),
        panel.grid.major = element_line(colour = "lightgrey") ,
        panel.border = element_rect(size=1.2, fill = NA),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_color_gradientn(colors=rgb.palette(100),na.value="white",limits=c(1,100),"Duration [h]")+
  theme(text = element_text(size=16))

##########################################################


load(file="ident_events_v5.Rdata")
qltim=ultim[which(ultim$Dominant.hazard=="Extreme wind"),]
windid=as.vector(qltim$ID_mod)
wix<-  strsplit(windid, ",")
isc=c()
listev<-c()
ism=c()
comid=c()
for(idw in 1:length(wix)){
  wmu=c()
  wcom=c()
  bof=wix[[idw]]
  bofv<-as.data.frame(matrix(unlist(bof),ncol=length(bof),byrow=T))

  if(length(bof)>1)wmu=rep(1,length(bof))
  if(length(bof)==1)wmu=2
  
  iscom=(match(as.numeric(bof),compound$ev2))
  icomid=na.omit(compound$combin[iscom])
  if(length(iscom)>1)wcom=T
  if(length(iscom)<=1)wcom=F
  comid=c(icomid,comid)
  isc=c(isc,wcom)
  ism=c(ism,wmu)
  listev=c(listev,bofv)
  
}

listev=as.numeric(listev)
listevp=cbind(listev,ism)
listevp=as.data.frame(listevp)
names(listevp)=c("V2","linetype")

qmin=quantile(metaHaf$viw.max,.5)
qt=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,.95,.99)
qlist=c(0,quantile(metaHaf$viw.max,.1),quantile(metaHaf$viw.max,.2),
        quantile(metaHaf$viw.max,.3),quantile(metaHaf$viw.max,.4),
        quantile(metaHaf$viw.max,.5),quantile(metaHaf$viw.max,.6),
        quantile(metaHaf$viw.max,.7),quantile(metaHaf$viw.max,.8),
        quantile(metaHaf$viw.max,.9),quantile(metaHaf$viw.max,.95),
        quantile(metaHaf$viw.max,.99))

bb<-metaHaf$ev[which(metaHaf$viw.max>qmin)]
plot(qlist)

ccl=c()
clo=c()
for(evi in bb){
  metev<-metavDaf[which(metavDaf$ev==evi),]
  lan=metaHaf$x.dur[which(metaHaf$ev==evi)]
  cl=c()
  for(q in 1:12){
    l=length(metev$vi.max[which(metev$vi.max>qlist[q])])
    cl=c(cl,l)
  }
  clo<-cbind(qt,rep(evi,11),cl,rep(lan,11))
  ccl<-rbind(ccl,clo)
}

wlen=length(listevp$V2)/length(qltim$Year)
coui=as.data.frame(ccl)
coui$cl=(coui$cl/1485)*100
cost=coui[which(coui$V2==wev),]
cool=inner_join(coui,listevp,by="V2")
rgb.palette=colorRampPalette(c("#ffffcc","#ffeda0", "#fed976", "#feb24c","#fd8d3c","#fc4e2a",
                               "#e31a1c","#bd0026","#800026"),interpolate="linear",bias=1)


ggplot(coui, aes(x = qt, y = cl, group=as.character(V2))) +
  geom_line(alpha=.5,color="grey",size=1.2)+

  scale_x_continuous(limits = c(0,1), expand = c(0, 0),"Intensity") +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0),"Spatial footprint [%]") +
  scale_size(range = c(0.3, 1.2)) +
  scale_linetype_identity() +
  
  geom_line(data=cool, aes(x = qt, y = cl, group=as.character(V2),color=V4),alpha=.7,size=1.2,linemitre=10) +
  geom_line(data=cost,aes(x = qt, y = cl, group=as.character(V2)),color="darkblue",size=1.2,linemitre=10) +
  # scale_alpha_continuous(c(0.7, 0.8),guide = FALSE)+
  theme(axis.text=element_text(size=16),
        panel.grid.major = element_line(colour = "lightgrey") ,
        panel.border = element_rect(size=1.2, fill = NA),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_color_gradientn(colors=rgb.palette(100),na.value="white",limits=c(1,100),"Duration [h]")+
  theme(text = element_text(size=16))
  # geom_point(size = 4, shape = 21) 



#########################################Compound########################
comidf=c(comidr,comid)


listev=as.numeric(listev)
listevp=cbind(listev,ism)
listevp=as.data.frame(listevp)
names(listevp)=c("V2","linetype")


tesp<-aggregate(list(len=sp03[,3]),
                by = list(ev2=sp03[,1],ev1 = sp03[,4],loc=sp03[,2]),
                FUN = function(x) c(length=length(unique(x))))
tesp <- do.call(data.frame, tesp)

tesp$cev=paste(tesp$ev1,tesp$ev2)

caca<-inner_join(tesp,metatest,by=c("loc","ev1"))
names(metavDaf)[1]="ev2"
caca2=inner_join(caca,metavDaf,by=c("loc","ev2"))
bivar=caca2[,c(10,7)]

## Evaluate empirical copula


ep=as.data.frame(pobs(bivar))

tesp=cbind(tesp,ep)
u <- epdata # random points were to evaluate the empirical copula
ecb <- C.n(u, X = bivar)
cbprob=1-qt-qt+ec

plot(as.numeric(cbprob))

qt=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,.95,.99)
qlist=c(0,quantile(metaHaf$viw.max,.1),quantile(metaHaf$viw.max,.2),
        quantile(metaHaf$viw.max,.3),quantile(metaHaf$viw.max,.4),
        quantile(metaHaf$viw.max,.5),quantile(metaHaf$viw.max,.6),
        quantile(metaHaf$viw.max,.7),quantile(metaHaf$viw.max,.8),
        quantile(metaHaf$viw.max,.9),quantile(metaHaf$viw.max,.95),
        quantile(metaHaf$viw.max,.99))

bb<-metaHaf$ev[which(metaHaf$viw.max>qmin)]
plot(qlist)
bij=((match(tesp$cev, compound$combin)))
ortime=compound$ORtime[bij]
tesp$ortime=ortime
ccl=c()
clo=c()
for(evi in unique(tesp$cev)){
  metev<-tesp[which(tesp$cev==evi),]
  lan=unique(tesp$ortime[which(tesp$cev==evi)])
  cl=c()
  for(q in 1:12){
    l=length(metev$cev[which(metev$vi.max>qt[q] & metev$len.sum>qt[q])])
    cl=c(cl,l)
  }
  clo<-cbind(qt,rep(evi,12),cl,rep(lan,12))
  ccl<-rbind(ccl,clo)
}

clen=length(comidf)/length(ulticom$Year)

coui=as.data.frame(ccl)
coui$V4=as.numeric(coui$V4)
coui$cl=(as.numeric(coui$cl)/1485)*100
coui$qt=as.numeric(coui$qt)
cool=coui[which(!is.na(match(coui$V2,comidf))),]
coc=coui[which(coui$V2==cevc),]
rgb.palette=colorRampPalette(c("#ffffcc","#ffeda0", "#fed976", "#feb24c","#fd8d3c","#fc4e2a",
                               "#e31a1c","#bd0026","#800026"),interpolate="linear",bias=1)


ggplot(coui, aes(x = qt, y = cl, group=as.character(V2))) +
  geom_line(alpha=.5,color="grey",size=1.2)+
  
  scale_x_continuous(limits = c(0,1), expand = c(0, 0),"Intensity") +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0),"Spatial footprint [%]") +
  scale_size(range = c(0.3, 1.2)) +
  scale_linetype_identity() +
  
  geom_line(data=cool, aes(x = qt, y = cl, group=as.character(V2),color=V4),alpha=.7,size=1.2,linemitre=10) +
  geom_line(data=coc,aes(x = qt, y = cl, group=as.character(V2)),color="darkblue",size=1.2,linemitre=10) +
  # scale_alpha_continuous(c(0.7, 0.8),guide = FALSE)+
  theme(axis.text=element_text(size=16),
        panel.grid.major = element_line(colour = "lightgrey") ,
        panel.border = element_rect(size=1.2, fill = NA),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_color_gradientn(colors=rgb.palette(100),na.value="white",limits=c(1,100),"Duration [h]")+
  theme(text = element_text(size=16))
# geom_point(size = 4, shape = 21) 



###################timeline##################
coki=rbind(rainevf,rainotro)

ggplot( data= coki, aes(x=vir.surf, fill=c,y=..density..)) +
  # geom_histogram( color="#e9ecef", alpha=0.4, position = 'identity') +
  geom_density(aes(x = vir.surf),alpha=.2,adjust=2) +
  geom_rug(aes(x = vir.surf, y = 0), position = position_jitter(height = 0))+
  # scale_fill_manual(values=c("#69b3a2", "#404080")) +
  scale_x_continuous(trans = log_trans(),
                     breaks =c(1,10,100,1000),limits=c(1,1500),"Spatial footprint [%]",
                     labels=c("0.1","1","10","100")) 
  # scale_x_continuous(trans = log_trans(),
  #                      breaks = c(3,6,12,24,48,72,148),limits = c(1.3,150),"Duration",
  #                      labels=c("3h","6h","12h","1day","2days","3days","1week"))+
  labs(fill="")


  
  ggplot( data= compound, aes(x= startime)) +
    geom_histogram( color="black", alpha=0.4, position = 'identity',bins=2000) +
    # geom_density(aes(x = startev),alpha=.4) +
    geom_rug(data=validationdays, aes(x = startime, y = 0), position = position_jitter(height = 0))
  
windevf<-metaHaf[-matchingw,]
windotro<-metaHaf[matchingw,]
windevf$c="w"
windotro$c="c"

cowi=rbind(windevf,windotro)

ggplot( data= cowi, aes(x=x.dur, fill=c,y=..density..)) +
  geom_histogram( color="#e9ecef", alpha=0.4, position = 'identity') +
  # geom_density(aes(x = viw.surf),alpha=.2,adjust=1) +
  geom_rug(aes(x = x.dur, y = 0), position = position_jitter(height = 0))+
  # scale_fill_manual(values=c("#69b3a2", "#404080")) +
  # scale_x_continuous(trans = log_trans(),
  #                    breaks =c(1,10,100,1000),limits=c(1,1500),"Spatial footprint [%]",
  #                    labels=c("0.1","1","10","100")) 
scale_x_continuous(trans = log_trans(),
                     breaks = c(3,6,12,24,48,72,148),limits = c(1.3,150),"Duration",
                     labels=c("3h","6h","12h","1day","2days","3days","1week"))+
labs(fill="")



write.csv(rainevf,'rainonly_79-19.csv')
write.csv(windevf,'windonly_79-19.csv')
write.csv(compoundfinalST,'compoundevents_79-19.csv')


rainevi2<-metaHar
rainevi2$c="r"
windevi2<-metaHaf
windevi2$c="w"

compevi<-compoundfinalST[,c(2,30,31,32,53,51)]
compevi$c="c"
names(compevi)[1]="ev"
names(rainevi)=names(compevi)
names(windevi)=names(compevi)
compound$year=year(compound$startime)

compevi2<-compound[,c(2,33,34,30,28,23,24,32,25,26)]
totevi<-rbind(compevi,rainevi,windevi)

bibcase<-rainevi2

bibis1<-aggregate(bibcase$season,
                  by = list(Month = bibcase$x.month),
                  FUN = function(x) c(n =length(x),x=x[1]))
bibis1<- do.call(data.frame, bibis1)
bibis1$x.x[which(bibis1$Month==9)]=4
bibis1$c=1
bibis1$prop<-bibis1$x.n/(bibis1$x.n+bibis2$x.n+bibis3$x.n)*100

vivis1<-aggregate(bibcase$season,
                  by = list(Month = bibcase$x.month,Year=bibcase$x.year),
                  FUN = function(x) c(n =length(x),x=x[1]))
vivis1<- do.call(data.frame, vivis1)
vivis1$x.x[which(vivis1$Month==9)]=4
vivis1$c=2
vivis1$date<-as.yearmon(paste(vivis1$Year, vivis1$Month, sep = "-"))
vivis1$pm<-vivis1$x.n/sum(vivis1$x.n)




bibcase<-windevi2
bibis2<-aggregate(bibcase$season,
                  by = list(Month = bibcase$x.month),
                  FUN = function(x) c(n =length(x),x=x[1]))
bibis2<- do.call(data.frame, bibis2)
bibis2$x.x[which(bibis2$Month==9)]=4
bibis2$c=2
bibis2$prop<-bibis2$x.n/(bibis1$x.n+bibis2$x.n+bibis3$x.n)*100


vivis2<-aggregate(bibcase$season,
                  by = list(Month = bibcase$x.month,Year=bibcase$x.year),
                  FUN = function(x) c(n =length(x),x=x[1]))
vivis2<- do.call(data.frame, vivis2)
vivis2$x.x[which(vivis2$Month==9)]=4
vivis2$c=1
vivis2$date<-as.yearmon(paste(vivis2$Year, vivis2$Month, sep = "-"))


bibcase<-compevi
bibis3<-aggregate(bibcase$season,
                  by = list(Month = bibcase$x.month),
                  FUN = function(x) c(n =length(x),x=x[1]))
bibis3<- do.call(data.frame, bibis3)
bibis3$x.x[which(bibis3$Month==9)]=4
bibis3$c=3



vivis3<-aggregate(bibcase$season,
                  by = list(Month = bibcase$x.month,Year=bibcase$x.year),
                  FUN = function(x) c(n =length(x),x=x[1]))
vivis3<- do.call(data.frame, vivis3)
vivis3$c=3



vivis2$pm<-vivis2$x.n/sum(vivis2$x.n)
vivis3$pm<-vivis3$x.n/sum(vivis3$x.n)

# bibcase<-totevi
# bibis4<-aggregate(bibcase$season,
#                   by = list(Month = bibcase$x.month),
#                   FUN = function(x) c(n =length(x),x=x[1]))
# bibis4<- do.call(data.frame, bibis3)
# bibis4$x.x[which(bibis4$Month==9)]=4
# bibis4$c=4


library(zoo)

vivis3$date<-as.Date( paste( vivis3$Month , vivis3$Year , sep = "." )  , format = "%m.%Y" )
vivis3$date<-as.yearmon(paste(vivis3$Year, vivis3$Month, sep = "-"))


plot(vivis1$pm)

ggplot(vivis3,aes(x= Month,y=Year)) +
  theme_bw(16)+
  geom_raster(aes(fill=x.n,group=c),interpolate = F)+
  scale_fill_gradientn(trans=scales::modulus_trans(.8),colors=rgb.palette(100))



vivix<-full_join(vivis1,vivis2,by=c("date"))

vivix<-full_join(vivix,vivis3,by=c("date"))

vivix$x.n[which(is.na(vivix$x.n))]=0
vivix$x.n.y[which(is.na(vivix$x.n.y))]=0
vivix$pmm<-vivix$pm.x*vivix$pm.y
vivix$mult<-vivix$pm/vivix$pmm
vivix$rat<-vivix$x.n/(vivix$x.n.x+vivix$x.n.y+vivix$x.n)*100
mean(vivix$rat)
ggplot(vivix,aes(x= factor(Month.x),y=Year.x)) +
  theme_bw(16)+
  geom_raster(aes(fill=rat,group=factor(c)),interpolate = F)+
scale_fill_gradient2(trans=scales::modulus_trans(1),low="white",mid ="orangered", high = "black",na.value="white",midpoint=20)+
  theme(axis.text.x = element_text(angle=0, vjust = 0.5,size=16))+
  theme(axis.text=element_text(size=16),
        panel.grid.major = element_line(colour = NA) ,
        panel.border = element_rect(size=1.2, fill = NA),
        
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  theme(text = element_text(size=16))+
  labs(y = "# of occurences",x = "Month",colour = "Parameter",fill= "Frequency [%]",size=20)+
  scale_x_discrete(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))+
  scale_y_continuous(name = "Years")

vivix$propr<-vivix$x.n.x/(vivix$x.n.x+vivix$x.n.y+vivix$x.n)*100
vivix$propw<-vivix$x.n.y/(vivix$x.n.x+vivix$x.n.y+vivix$x.n)*100
vivix$propc<-vivix$x.n/(vivix$x.n.x+vivix$x.n.y+vivix$x.n)*100
vivix$c.y=1
vivix$c=3


vivis4<-aggregate(list(x=compevi2$rf.max.cs,y=compevi2$wg.max.cs),
                  by = list(Month = compevi2$month,Year=compevi2$year),
                  FUN = function(x,y) c(n =cor.test(x,y,method="kendall")))
vivis4<- do.call(data.frame, vivis4)

tg<-cor.test(compevi2$rf.max.cs,compevi2$wg.max.cs,method="kendall")$estimate

omerde<-ddply(compevi2,c("sgroup","year"),function(x) cor.test(x$rf.max.cs,x$wg.max.cs,method="kendall")$estimate)
plot(omerde$sgroup,omerde$tau)

vivis4$g="r"
vivis6<-aggregate(list(m=compevi2$wg.max.cs),
                  by = list(Month = compevi2$month,Year=compevi2$year),
                  FUN = function(x) c(n =mean(x)))
vivis6<- do.call(data.frame, vivis6)
vivis6$g="w"
vivis7<-rbind(vivis4,vivis6)

vivis5<-aggregate(list(t=compevi2$ORtime,s=compevi2$spacescale,w=compevi2$wg.max.cs,r=compevi2$rf.max.cs),
                  by = list(season=compevi2$sgroup,year=compevi2$year,month=compevi2$month), 
                  FUN = function(x) c(n =mean(x)))
vivis5<- do.call(data.frame, vivis5)

vivis3$x.x[which(vivis3$Month==9)]=4

bibis3$prop<-bibis3$x.n/(bibis1$x.n+bibis2$x.n+bibis3$x.n)*100

bibis<-rbind(bibis1,bibis2,bibis3)
vivis<-rbind(vivis1,vivis2,vivis3)

names(vivix)[c(3,5,17)]=names(vivix)[c(14,16,19)]
names(vivix)[c(9,11,18)]=names(vivix)[c(14,16,19)]
vivisp<-rbind(vivix[,c(1,2,14,16,19)],vivix[,c(1,2,3,5,17)],vivix[,c(1,2,9,11,18)])

mis1<-aggregate(vivisp$propc[which(vivisp$c==1)],
                by = list(Month =  vivisp$Month.x[which(vivisp$c==1)]),
                FUN = function(x) c(n =median(x)))
mis1<- do.call(data.frame, mis1)
mis1$g=2

mis2<-aggregate(vivisp$propc[which(vivisp$c==2)],
                by = list(Month =  vivisp$Month.x[which(vivisp$c==2)]),
                FUN = function(x) c(n =median(x)))
mis2<- do.call(data.frame, mis2)
mis2$g=1

mis3<-aggregate(vivisp$propc[which(vivisp$c==3)],
                by = list(Month = vivisp$Month.x[which(vivisp$c==3)]),
                FUN = function(x) c(n =median(x)))
mis3<- do.call(data.frame, mis3)
mis3$g=3

ggplot(vivisp)+
  geom_path(aes(x=factor(Month.x), y=propc,group=factor(paste(Year.x,c)),color=factor(c)),alpha=.3,size=1)+
  geom_line(data=mis1,aes(x=factor(Month), y=x,group=g,color=factor(1)),size=2)+
  
  geom_line(data=mis2,aes(x=factor(Month), y=x,color=factor(2),group=g),size=2)+
  geom_line(data=mis3,aes(x=factor(Month), y=x,color=factor(3),group=g),size=2)+
  scale_color_manual(values = c(  "chocolate4","dodgerblue4","darkgreen"),labels = c("Wind", "Rain", "Compound"))+
  theme(axis.text.x = element_text(angle=0, vjust = 0.5,size=16))+
  theme(axis.text=element_text(size=16),
        panel.grid.major = element_line(colour = "lightgrey") ,
        panel.border = element_rect(size=1.2, fill = NA),
        
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  theme(text = element_text(size=16))+
  labs(y = "# of occurences",x = "Month",colour = "Parameter",fill= "Hazard events",size=20)+
  scale_x_discrete(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c("Jan","Fev","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))+
  scale_y_continuous(name = "Frequency [%]")

dodge <- position_dodge(width = 1)

hiv1l<-c(1,1,3.5,3.5)
hiv1L<-c(100,1000,1000,100)

h1<-data.frame(hiv1l,hiv1L)
h1$g<-1


ggplot(bibis,aes(x=Month, y=x.n,group=c,order=c))+geom_col(aes(fill=factor(c)),size=1,position = "dodge",alpha=.9) +
  
  scale_fill_manual(values = c(  "dodgerblue","darkorange","forestgreen"),labels = c( "Rain", "Wind","Compound")) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c("Jan","Fev","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))+
  scale_y_continuous(name = "# of occurences",sec.axis = sec_axis(~./20, name = "% of compound events"))+
  labs(y = "# of occurences",x = "Month",colour = "Parameter",fill= "Hazard events",size=20)+
  geom_line(data=bibis, mapping = aes(x = Month, y = prop*20,color=factor(c)), size = 1.5, linetype=1)+
  theme_classic() +
  theme(axis.text.x = element_text(angle=0, vjust = 0.5,size=16))+
  theme(text = element_text(size=16))
  
ggplot(vivis, aes(x=factor(Month), y=x.n,fill=factor(c),color=factor(c))) +
  scale_fill_manual(values = c(  "darkorange","dodgerblue","forestgreen"),labels = c("Wind", "Rain", "Compound"))+
  scale_color_manual(values = c(  "chocolate4","dodgerblue4","darkgreen"),labels = c("Wind", "Rain", "Compound"))+
  stat_boxplot(geom = "errorbar")+
  geom_boxplot(position="dodge2",size=1) +
  theme(axis.text.x = element_text(angle=0, vjust = 0.5,size=16))+
  theme(axis.text=element_text(size=16),
        panel.grid.major = element_line(colour = "lightgrey") ,
        panel.border = element_rect(size=1.2, fill = NA),

        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  theme(text = element_text(size=16))+
  labs(y = "# of occurences",x = "Month",colour = "Parameter",fill= "Hazard events",size=20)+
  scale_x_discrete(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c("Jan","Fev","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))+
  scale_y_continuous(name = "monthly occurences")

metaHar$con=1
metaHar$con[matchingr]=2

metaHaf$con=1
metaHaf$con[matchingw]=2

metaHar$footprint=metaHar$vir.surf/1485
metaHaf$footprint=metaHaf$viw.surf/1485

maprob<-c()
laprob<-c()
diff<-c()
las<-seq(.01,100,by=0.1)
for(i in 1:length(las)){
 
  l=las[i]
  noco<-length(metaHar$ev[which(metaHar$footprint>l)])
  wco<-length(metaHar$ev[which(metaHar$con==2 & metaHar$footprint>l)])
  
  moco<-length(metaHaf$ev[which(metaHaf$footprint>l)])
  zco<-length(metaHaf$ev[which(metaHaf$con==2 & metaHaf$footprint>l)])

  maprob[i]=wco/noco
  laprob[i]<-zco/moco
}
plot(las,maprob,ylim=c(0,1),log="x",type="p",pch=1,lwd=2,col=alpha("dodgerblue",.9),ylab="Proportion of CE",
     xlab="Footprint [%]",xaxt="n",cex.lab=1.6,cex.axis=1.5)
points(las,laprob,col=alpha("darkorange",.9),pch=4,lwd=2,type="p")
axis(side=1,at=c(0.1, 1, 10, 25, 50, 100), labels= c(0.1, 1, 10, 25, 50, 100),cex.axis=1.5)
abline(h=50)
merdicus<-data.frame(maprob,las)
maprob<-c()
laprob<-c()
diff<-c()
max(metaHaf$x.dur)
las<-seq(1,100,by=1)
for(i in 1:length(las)){
  
  l=las[i]
  noco<-length(metaHar$ev[which(metaHar$x.dur>l)])
  wco<-length(metaHar$ev[which(metaHar$con==2 & metaHar$x.dur>l)])
  
  moco<-length(metaHaf$ev[which(metaHaf$x.dur>l)])
  zco<-length(metaHaf$ev[which(metaHaf$con==2 & metaHaf$x.dur>l)])
  
  maprob[i]=wco/noco
  laprob[i]<-zco/moco
}
plot(las,maprob,ylim=c(0,1),type="p",pch=1,lwd=2,col=alpha("dodgerblue",.9),ylab="Proportion of CE",
     xlab="Duration [h]",xaxt="n",cex.axis=1.5,cex.lab=1.6)
points(las,laprob,col=alpha("darkorange",.9),pch=4,lwd=2,type="p")
axis(side=1,at=c( 1, 12, 24, 48, 72, 96), labels= c( 1, 12, 24, 48, 72, 96),cex.axis=1.5)
abline(v=50)

length(matchingr)/length(metaHar$ev)
length(matchingw)/length(metaHaf$ev)

aaah<-metaHar[order(metaHar$vir.max,decreasing = T),]
aah<-aaah[1:100,]
length(aah$ev[which(aah$con==2)])/100 

aaah<-metaHaf[order(metaHaf$viw.max,decreasing = T),]
aah<-aaah[1:100,]
length(aah$ev[which(aah$con==2)])/100 

                 
ggplot(data=metaHar,aes(x=vir.surf/1485,y=x.dur,group=factor(con),color=factor(con))) +
  geom_point(size=3,alpha=.5) + 
  scale_color_manual(values = c(  "lightslateblue","lightgoldenrod"),labels = c("chaz", "shaz"))+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24),
        panel.background = element_rect(fill = "white", colour = "black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size=24),
        legend.text = element_text(size=18),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_y_continuous(trans = log_trans(),
                     breaks = c(3,6,12,24,48,72,148),limits = c(1.9,150),"Duration",
                     labels=c("3h","6h","12h","1day","2days","3days","1week"))+
  scale_x_continuous(trans = log_trans(),
                     breaks =c(0.001,0.01,0.1,1),limits=c(0.0006,1.2),"Spatial footprint [%]",
                     labels=c("0.1","1","10","100")) 
  
  

ggplot(compevi2, aes(x=factor(sgroup), y=ORtime, fill=factor(sgroup))) +
  scale_fill_manual(values = c(  "lightslateblue","lightgoldenrod"),labels = c("ONDJFM", "AMJJAS"))+
  scale_color_manual(values = c(  "chocolate4","dodgerblue4","darkgreen"),labels = c("Wind", "Rain", "Compound"))+
  stat_boxplot(geom = "errorbar",width = 0.5,size=1)+
  geom_boxplot(position="dodge",size=1) +
  theme(axis.text.x = element_text(angle=0, vjust = 0.5,size=16))+
  theme(axis.text=element_text(size=16),
        panel.grid.major = element_line(colour = "lightgrey") ,
        panel.border = element_rect(size=1.2, fill = NA),
        
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  # geom_hline(yintercept =mean(vivis5$t),size=2, col="red")+
  theme(text = element_text(size=16))+
  labs(y = "# of occurences",x = "Season",colour = "Parameter",fill= "Season",size=20)+
  scale_x_discrete(breaks=c(1,2),labels=c("ONDJFM", "AMJJAS"))+
  scale_y_continuous(name = "Duration",trans="log",breaks = c(1,6,12,24,36,48,72))

totevi$spacescale.y=totevi$spacescale.y/1485


ggplot(totevi, aes(x=factor(c), y=spacescale.y, fill=factor(c))) +
  scale_fill_manual(values = c( "forestgreen", "darkorange","dodgerblue"),labels = c("Compound","Wind", "Rain"))+
  scale_color_manual(values = c(  "chocolate4","dodgerblue4","darkgreen"),labels = c("Wind", "Rain", "Compound"))+
  stat_boxplot(geom = "errorbar",width = 0.5,size=1)+
  geom_boxplot(position="dodge",size=1) +
  theme(axis.text.x = element_text(angle=0, vjust = 0.5,size=16))+
  theme(axis.text=element_text(size=16),
        panel.grid.major = element_line(colour = "lightgrey") ,
        panel.border = element_rect(size=1.2, fill = NA),
        
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  # geom_hline(yintercept =mean(vivis5$t),size=2, col="red")+
  theme(text = element_text(size=16))+
  labs(y = "# of occurences",x = "Season",colour = "Parameter",fill= "Season",size=20)+
  scale_x_discrete(breaks=c(1,2),labels=c("ONDJFM", "AMJJAS"))+
  # scale_y_continuous(name = "Duration",trans="log",breaks = c(1,6,12,24,36,48,72))
  scale_y_continuous(name = "Footprint",trans="log",breaks = c(.001,.01,.1,.2,.5,.75,1),labels=c(0.1,1,10,20,50,75,100))





mis1<-aggregate(vivis1$x.n,
                by = list(Month = vivis1$Month),
                FUN = function(x) c(n =mean(x)))
mis1<- do.call(data.frame, mis1)
mis1$g=1

mis2<-aggregate(vivis2$x.n,
                by = list(Month = vivis2$Month),
                FUN = function(x) c(n =mean(x)))
mis2<- do.call(data.frame, mis2)
mis2$g=1

mis3<-aggregate(vivis3$x.n,
                  by = list(Month = vivis3$Month),
                  FUN = function(x) c(n =mean(x)))
mis3<- do.call(data.frame, mis3)
mis3$g=1

ggplot(vivis)+
  geom_path(aes(x=factor(Month), y=x.n,group=factor(paste(Year,c)),alpha=Year,color=factor(c)),size=1)+
  geom_point(data=mis3,aes(x=factor(Month), y=x),size=4, col="black")+
  geom_line(data=mis3,aes(x=factor(Month), y=x,group=g),size=2, col="black")
  
  


geom_line(data=bibis, mapping = aes(x = Month, y = prop,color=factor(c)), size = 2, linetype=1)
mean(metaHax[[1]]$x.dur)
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

min(compound$wg.mean.cs)


days<-seq(as.Date('1979-01-01'), as.Date('2019-09-01'), by = "7 days")

interarival<-diff(compound$startime)/3600

hist(compound$ORtime,breaks=50)
hist(windevf$x.dur ,breaks=50)

cor(compound$ORtime,compound$footprint)
cor(rainevf$x.dur,rainevf$footprint)

tt<-sum(metaHaf$viw.surf)+sum(metaHar$vir.surf)
vv<-sum(windevf$viw.surf)+sum(rainevf$vir.surf)
cc<-sum(compound$spacescale)
cc/tt
######################Spatiotemporal plot####################################
compoundfinalST$wg.max.s=compoundfinalS$wg.max.y
compoundfinalST$rf.max.s=compoundfinalS$rf.max.x
compoundfinalST$footprint=compoundfinalST$spacescale.x/1485
compound$footprint=compound$spacescale/1485
metaHar$footprint=metaHar$vir.surf/1485
bestof<-metaHaf[which(metaHaf$viw.max>45),]
rbPal <- colorRampPalette(c('lightskyblue',"skyblue","gold","darkorange",'red',"purple"))
# ggplot(data=compoundfinalST,aes(x=spacescale,y=ORtime,colour=wg.max.s,size=rf.max.s))+
max(metaHax[[1]]$rf.surf)
a=2

compound$shi=1
plot(rati$ratiospace,rati$ratiotime)
compound2<-compound[which(rati$ratiospace>0.5|rati$ratiotime>0.5),]
length(compound2[which(compound2$sgroup==1),1])/906
plot(compound2$rf.max.cs,compound2$wg.max.cs)
rspace<-1/rank(compound$spacescale)
plot(sort(rspace))
library(poweRlaw)
m_pl<- conpl$new(compound$ORtime)
est = estimate_xmin(m_pl)
m_pl$setXmin(est)
plot(m_pl)
lines(m_pl, col = 2, lwd = 2)
m_bl_ln = conlnorm$new(compound$spacescale)
est = estimate_xmin(m_bl_ln)
m_bl_ln$setXmin(est)
lines(m_bl_ln, col = 3, lwd = 2)


meerd<-compound$ORtime
meerd<-metaHaf$viw.surf
plot(sort(meerd))
descdist(meerd,boot=1000)
shit<-fitdist(meerd,"lnorm")
plot(shit)
plot(log(compound$spacescale),rspace)
plot(sort(compound$footprint))
plot(sort(metaHaf$x.dur))
max(compound$ORtime)
compound$spt<-log(compound$spacescale*compound$ORtime)
hist(compound$spt,breaks=20)

fepeter<-aggregate(compound$combin,
                  by = list(sp=compound$footprint,t=compound$ORtime),
                  FUN = function(x) c(n =length(x)))
fepeter <- do.call(data.frame, fepeter)

rfepeter<-aggregate(metaHar$ev,
                   by = list(sp=metaHar$footprint,t=metaHar$x.dur),
                   FUN = function(x) c(n =length(x)))
rfepeter <- do.call(data.frame, rfepeter)

wfepeter<-aggregate(metaHaf$ev,
                   by = list(sp=metaHaf$footprint,t=metaHaf$x.dur),
                   FUN = function(x) c(n =length(x)))
wfepeter <- do.call(data.frame, wfepeter)

mod <- lm(ORtime ~ footprint, data = compoundmax)

plot(metaHar$metave...c.2..)

compoundmax<-compound[which(compound$spm>1),]
rainemax<-rainevf[which(rainevf$spm>1),]
windemax<-windevf[which(windevf$spm>1),]
level=c(0.1,0.01)

com.dens<-kde2d(log(compoundmax[,c(33)]),log(compoundmax[,30]),n=1000)

contour_95 <- with(com.dens, contourLines(x=x, y=y,
                                    z=z, levels=c(0.05)))
contour_95 <- data.frame(contour_95)
contour_95=exp(contour_95)

contour(x=com.dens$x, y=com.dens$y,
             z=com.dens$z, levels=c(.1,.01,0.15,.4,.5,.6,.7,.99))


compoundmax$footprint=compoundmax$spacescale/1485



ggplot(data=compoundmax,aes(x=footprint,y=timescale))+
  # geom_jitter(shape=21,alpha=0.5,fill="transparent", width = 0.3, height= 0.2,size=2)+
  geom_point(col="green", alpha=.2) +
  # geom_path(aes(x, y), data=contour_95,color="blue", size=1) +
  # stat_density_2d(data=compoundmax,aes(x=footprint,y=ORtime,weight=rf.max,fill=..level..),alpha=.1,geom="polygon",contour_var="ndensity",color="darkgreen",size=1)+
  # stat_density_2d(data=rainemax,aes(x=footprint,y=x.dur,fill=..level..),alpha=.1,geom="polygon",contour_var = "ndensity",color="royalblue",size=1)+
  # stat_density_2d(data=windemax,aes(x=footprint,y=x.dur,fill=..level..),alpha=.1,geom="polygon",contour_var = "ndensity",color="orange",size=1)+
  # geom_jitter(data= rfepeter,aes(x=sp,y=t,size=x), shape=21,alpha=0.3,fill="blue") +
  # geom_jitter(data= wfepeter,aes(x=sp,y=t, size=x), shape=21,alpha=0.3,fill="orange") +
  # geom_jitter(data=rainemax,aes(x=footprint,y=x.dur),shape=21,alpha=0.1,fill="transparent",width = 0.3, height= 0.3,size=1,color="blue")+
  # geom_jitter(data=rainotro,aes(x=footprint,y=x.dur),shape=21,alpha=0.2,fill="grey",width = 0.3, height= 0.3,size=3)+
  # geom_jitter(data=windemax,aes(x=footprint,y=x.dur),shape=21,alpha=0.1,fill="orange",width = 0.3, height= 0.3,size=1)+
  # geom_point(shape=21)+
  # geom_abline(intercept = 4,size=2,slope= .6)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  # geom_vline(xintercept = 1,size=2,alpha=.5)+
  # geom_hline(yintercept=54)+
  # geom_vline(xintercept=1306)+
 
  scale_y_continuous(trans = log_trans(),
                   breaks = c(3,6,12,24,48,72,148),limits = c(1.3,150),"Duration",
                    labels=c("3h","6h","12h","1day","2days","3days","1week"))+
  scale_x_continuous(trans = log_trans(),
                     breaks =c(0.001,0.01,0.1,1),limits=c(0.0004,1.2),"Spatial footprint [%]",
                     labels=c("0.1","1","10","100")) 
  scale_fill_gradientn(trans=scales::modulus_trans(.2),colors=rgb.palette(100),"Duration [h]",limits=c(0,15))
 
   # scale_colour_manual(values = c( "blue", "red","orange","skyblue"))+
  # scale_size(trans=scales::modulus_trans(.6),range=c(1,15),"# of occurencces",limits=c(0,300))
# guide = guide_legend(override.aes = list(colour = alpha("grey",.8))))
  
  # scale_fill_gradientn(trans=scales::modulus_trans(2),colors=rbPal(100),"Mean wa [m/s]",breaks=c(20,30,40,Inf),limits=c(1,43))+
  #   # scale_colour_manual(values = c( "blue", "red","orange","skyblue"))+
  #   scale_size(trans=scales::modulus_trans(1.8),range=c(1,14),"Mean ra [mm]",breaks=c(25,50,75,Inf),limits=c(1,75),guide = guide_legend(override.aes = list(colour = alpha("grey",.8))))

compoundmax$sgr=3
rainemax$sgr=1
windemax$sgr=2

ilmonstro=rbind(as.matrix(compoundmax[,c(33,30,35)]),as.matrix(rainemax[c(17,8,19)]),as.matrix(windemax[c(15,8,17)]))
ilmonstro=data.frame(ilmonstro)
ilmonstro$ORtime=jitter(ilmonstro$ORtime,2)
ilmonstro$footprint=jitter(ilmonstro$footprint,2)

d=ggplot(data=ilmonstro,aes(x=footprint,y=ORtime))
d +     geom_point(aes(color=as.character(sgr)),alpha=0.1)+
  scale_color_manual(values=c("royalblue","orange","darkgreen"))+
  stat_density_2d(contour_var = "ndensity",aes(fill=..level..),alpha=0.3,geom="polygon") + 

    facet_wrap(vars(sgr))+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  
  scale_y_continuous(trans = log_trans(),
                     breaks = c(3,6,12,24,48,72,148),limits = c(1.3,150),"Duration",
                     labels=c("3h","6h","12h","1day","2days","3days","1week"))+
  scale_x_continuous(trans = log_trans(),
                     breaks =c(0.001,0.01,0.1,1),limits=c(0.0004,1.2),"Spatial footprint [%]",
                     labels=c("0.1","1","10","100")) +
  scale_fill_viridis_b()

tr<-(compound$timescale/ compound$ORtime)
min(tr)
max(tr)
plot(tr,compound$wg.mean.cs)

library(ggpointdensity)
library(ggalt)
compound$sgr=1
compoundX$sgr=1

ggplot(data=compound,aes(x=footprint,y=ORtime,group=sgr)) +
  geom_pointdensity(size=3,adjust = 1,shape=16,method="kde2d") + scale_color_viridis_c(option="magma")+
  # geom_point(size=4,alpha=.6 ,shape=16,color="darkgrey") + 
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24),
        panel.background = element_rect(fill = "white", colour = "black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_text(angle=90,size=24),
        legend.text = element_text(size=18),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_y_continuous(trans = log_trans(),
                     breaks = c(3,6,12,24,48,72,148),limits = c(1,150),"Duration",
                     labels=c("3h","6h","12h","1day","2days","3days","1week"))+
  scale_x_continuous(trans = log_trans(),
                     breaks =c(0.1,1,10,100),limits=c(0.01,1.2),"Spatial footprint [%]",
                     labels=c("0.1","1","10","100")) 
  # geom_encircle(data=compound,aes(x=footprint,y=ORtime,group=1), size = 5, color = "darkgrey",s_shape=1.3, expand=0.03,fill=NA,alpha=.6)+
  # geom_point(data=compoundX,(aes(x=footprint,y=ORtime,group=sgr,color=1/(rp*111))), size=3,stroke=4,shape=1,alpha=.8) +
  # scale_colour_gradientn(trans="log",colors=rev(rbPal2(100)),"Return period (Y)",breaks=c(1,10,100))+
  # geom_encircle(data=compoundX,aes(x=footprint,y=ORtime,group=1), size = 5, color = "red",s_shape=1.3, expand=0.03,fill=NA,alpha=.6)


compound$rpa<-compound$rp
compound$rpa[which(is.na(compound$rpa))]=.9

ggplot(data=compound,aes(x=footprint/100,y=ORtime,fill=rp))+
  geom_point(shape=21,stroke=1,size=3,alpha=.001+1/(111*compound$rpa)/250)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  # geom_vline(xintercept = 1,size=2,alpha=.5)+
  # geom_hline(yintercept=54)+
  # geom_vline(xintercept=1306)+
  scale_fill_gradientn(trans="log",colors=rbPal2(100),"Joint Return Period (Y)",breaks=c(1/111,.1/111,.01/111),labels = c(1,10,100))+
  # scale_colour_manual(values = c( "blue", "red","orange","skyblue"))+
  # scale_colour_gradientn(trans=scales::modulus_trans(1.2),colors=rbPal2(100),"Max precipitation [mm]",breaks=c(40,80,120,Inf),limits=c(1,141))+
  scale_y_continuous(trans = log_trans(),
                     breaks = c(3,6,12,24,48,72,148),limits = c(2.9,150),"Duration",
                     labels=c("3h","6h","12h","1day","2days","3days","1week"))+
  scale_x_continuous(trans = log_trans(),
                     breaks =c(0.001,0.01,0.1,1),limits=c(0.0006,1.2),"Spatial footprint [%]",
                     labels=c("0.1","1","10","100")) 

plot(compound$rpa/max(compound$rpa))

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


compound$sgr<-1
compound$sgr[which(month(compound$startime)>3 & month(compound$startime)<10)]=2
length(compound$combin[which(compound$sgr==1)])/length(compound$combin)

##############################Cell by cell= where is the most prone to hazard comb?#########################

load(file="out/RainEv_metap_1979-2019.Rdata")
metaHar=metaHaf
load(file="out/WindEv_meta_1979-2019.Rdata")

load(file="out/RainEv_ldat_1979-2019.Rdata")
metavDar=metavDaf
load(file="out/WindEv_ldat_1979-2019.Rdata")

cellcount<-metavDar[which(!is.na(match(metavDax[[1]]$ev,compoundfinalST$ev1))),]
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
sp03$sgr<-1
sp03$cev=paste(sp03$ev.x, sp03$ev.y,sep=" ")

sp03$sgr[which(month(sp03$dd)>3 & month(sp03$dd)<10)]=2
sp03w<-sp03[which(sp03$sgr==1),]
sp03s<-sp03[which(sp03$sgr==2),]

tesp<-aggregate(list(len=sp03[,3]),
                by = list(ev2=sp03[,1],ev1 = sp03[,4],loc=sp03[,2]),
                FUN = function(x) c(length=length(unique(x))))
tesp <- do.call(data.frame, tesp)

tesp$cev=paste(tesp$ev1,tesp$ev2)

caca<-inner_join(tesp,metatest,by=c("loc","ev1"))
names(metavDaf)[1]="ev2"
caca2=inner_join(caca,metavDaf,,by=c("loc","ev2"))
bivar=caca2[,c(10,7)]
plot(bivar)



ucount2<-sp03
ucount2$cev=paste(ucount2$ev.x, ucount2$ev.y,sep=" ")
sp03$loc<-as.character(sp03$loc)

plotmap<-list()
locationHr<-"-0.5 51.5"
locationSh<-"-1.5 53.5"
locationCu<-"-2.75 54.25"
locationGl<-"-4.25 56"
locas<-c(locationCu,locationGl,locationSh,locationHr)
namesloc<-c("Cumbria","Glasgow","Sheffield","London")
for(lo in 1:4){

location=locas[lo]

ucouthH<-sp03[which(sp03$loc==location),]
mates<-which(!is.na(match(sp03$cev,ucouthH$cev)))
ucountmH<-ucount2[mates,]

tip<-aggregate(list(len=ucountmH[,3]),
                by = list(ev2=ucountmH[,1],ev1 = ucountmH[,4],loc=ucountmH[,2]),
                FUN = function(x) c(length=length(unique(x))))
tip <- do.call(data.frame, tip)

heathrowext<-aggregate(tip[,4] ,
                       by = list(loc=tip[,3]),
                       FUN = function(x) c(n =sum(x), l=length(x)))

heathrowext<-do.call(data.frame,heathrowext)
heathrowext$loc=as.character(heathrowext$loc)


Rcast=strsplit(heathrowext$loc," ")

Rcas<-as.data.frame(matrix(unlist(Rcast),ncol=2,byrow=T))
heathrowext$lon<-as.numeric((Rcas[,1]))
heathrowext$lat<-as.numeric((Rcas[,2]))
heathrowext$gr=1
heathrowext$l.p<-heathrowext$x.l/max(heathrowext$x.l)
heathrowext$n.p<-heathrowext$x.n/max(heathrowext$x.n)

plotmap[[lo]]<-ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  theme_bw(16)+
  geom_raster(data=heathrowext,aes(x=lon,y=lat,fill=l.p,group=gr),interpolate = F)+
  # geom_contour_fill(data=elmaxoub1,aes(x=Var1,y=Var2,z=mwg,group=cluster,alpha=.1),na.fill=T)+
  # geom_contour(data=elmaxoub,colour="black",aes(x=Var1,y=Var2,z=mwg,group=cluster))+
  geom_polygon(fill = "transparent", color = "gray10", size = 1) +
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
  labs(y = "Latitude",x = "Longitude",fill= "Ps",size=20,title=namesloc[lo])+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="italic"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=14),
        legend.text = element_text(angle=0,size=13),
        legend.position = "right",
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_y_continuous(
    breaks = c(48,50,52,54,56,58,60),limits = c(40,70),"Latitude")+
  scale_x_continuous(
    breaks =c(-6,-4,-2,0,2),limits=c(-10,10),"Longitiude")+
  # scale_fill_gradientn(trans=scales::modulus_trans(1.5),colors=rbPal(100),"Max wind gust [m/s]",breaks=c(10,20,30,40,50,Inf),limits=c(5,51))+
  scale_fill_gradientn(trans=scales::modulus_trans(1.6), colours = rgb.palette(100),limits=c(0,1))
}

lay <- rbind(c(1,1,2,2),
             c(3,3,4,4))

grid.arrange(plotmap[[1]],plotmap[[2]],plotmap[[3]],plotmap[[4]], layout_matrix = lay, top="test") 

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

#########################WEATHER TYPES#################
compound$day<-as.Date(compound$startime)
WeatherTypes=read.csv(paste0(dir.study,"/in/ERA5_1979_2018_WT_UK.csv"),sep=",",dec=".",header=T)  

WeatherTypes$Date<- paste(WeatherTypes$year,WeatherTypes$month,WeatherTypes$ï..day,sep= " ")
WeatherTypes$Date = as.Date(strptime(WeatherTypes$Date ,format = "%Y%m%d"))
wahou<-na.omit(match(compound$day,WeatherTypes$Date))
owow<-na.omit(match(WeatherTypes$Date,compound$day))
int<-WeatherTypes[wahou,]
ziz<-compound[owow,]

int$day=int$Date
data2<-inner_join(int,compound,by=c("day"))
data3<-data2[na.omit(unique(match(compound$combin,data2$combin))),]


datowow<-aggregate(list(ev=data3$combin),
                   by = list(lwt=data3$LWT),
                   FUN = function(x) c(length=length(x)))
datowow<- do.call(data.frame, datowow)
datowow$freq=datowow$ev/length(data3$ï..day)
rowow<-aggregate(list(ev=WeatherTypes$Date),
                   by = list(lwt=WeatherTypes$LWT),
                   FUN = function(x) c(length=length(x)))
rowow<- do.call(data.frame, rowow)
rowow$freq<-rowow$ev/length(WeatherTypes$month)
laverif<-match(data3$combin,compound$combin)
compoundpls<-compound[match(data3$combin,compound$combin),]
compoundpls$lwt=as.character(data3$LWT)
hist(data3$LWT)
hist(WeatherTypes$LWT)


sp<-data.frame(data2$`Speed_m/s`)
length(unique(WeatherTypes$LWT))
summary(data2)
lwt<-data2$`int$LWT`
lwt<-sort(lwt)
lwt<-lwt[-which(lwt==-1)]
plot(lwt)
hist(data2$`int$LWT`, breaks = 27)
names(data2)[8]="LWT"

cols <- c("8" = "red", "4" = "blue", "6" = "darkgreen", "10" = "orange")

load(file="out/CompoundRW_79-19.RP.Rdata")
load(file="out/CompoundRW_79-19.ME.Rdata")
cols= c("-1"="white","0"="skyblue","1"="grey","2"="grey","3"="grey","4"="grey","5"="grey","6"="grey","7"="grey","8"="grey","11"="grey","12"="orange","13"="gold","14"="yellow","15"="yellowgreen","16"="green","17"="darkgreen","18"="black","20"="tomato3","21"="tomato","22"="tomato","23"="tomato", "24"="tomato","25"="tomato","26"="tomato","27"="tomato","28"="tomato") 
ggplot(data=xtremdat,aes(x=wg.max.cs,y=rf.max.cs,color=lwt))+
  geom_point(size=4,shape=19,alpha=0.7)+
  
  # scale_fill_gradientn(colours = rbPal2(30))+
  # scale_fill_gradientn(na.value = "gray95" ,trans="log",colors=rbPal2(100),"Joint Return Period (Y)",breaks=c(1/111,.1/111,.01/111),labels = c(1,10,100))+
  scale_colour_manual(values = cols)+
  # ggplot(data=metaHar,aes(x=footprint,y=x.dur,fill=vir.max,size=len.max))+
  # geom_point(shape=21)+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="italic"),
        panel.background = element_rect(fill = "white", colour = "black",size=1),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18),
        axis.line = element_line(colour = "black"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.key.size = unit(1, "cm"))+
  geom_line(data=inil, aes(x=X1,y=X2),size=2,color="lightblue4",alpha=.7)+
  scale_x_continuous(limits=c(17,50),expression(paste("w" ['a']," [m s "^-1,"]")))+
  scale_y_continuous(expression(paste("r" [a]," [mm]")),breaks=c(0,25,50,75,100,125,150))


xtremdat<-compoundpls[which(!is.na(compoundpls$rp)),]

powow<-aggregate(list(ev=xtremdat$combin),
                 by = list(lwt=xtremdat$lwt),
                 FUN = function(x) c(length=length(x)))
powow<- do.call(data.frame, powow)
powow$freq=powow$ev/222
ggplot(data3, aes(x=LWT)) + geom_histogram(binwidth=1)
plot(data3$wg.max.cs,data3$rf.max.cs,col=data3$LWT[which(data3$LWT>=0)], pch=16,col=rgb.palette())
compoundpls$lwtn<-as.numeric(compoundpls$lwt)
compoundpls$mainfam=0
compoundpls$mainfam[which(compoundpls$lwtn>=20)]=2
compoundpls$mainfam[which(compoundpls$lwtn<18 & compoundpls$lwtn>14)]=1
compoundpls$mainfam[which(compoundpls$lwtn<10)]=3

ggplot(data=xtremdat,aes(x=footprint,y=ORtime,color=lwt)) +
  # geom_pointdensity(size=3,adjust = 1,shape=16,method="kde2d") + scale_color_viridis_c(option="magma")+
  geom_point(size=4,alpha=.6 ,shape=16) + 
  scale_colour_manual(values = cols)+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24),
        panel.background = element_rect(fill = "white", colour = "black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_text(angle=90,size=24),
        legend.text = element_text(size=18),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_y_continuous(trans = log_trans(),
                     breaks = c(3,6,12,24,48,72,148),limits = c(1.9,150),"Duration",
                     labels=c("3h","6h","12h","1day","2days","3days","1week"))+
  scale_x_continuous(trans = log_trans(),
                     breaks =c(0.1,1,10,100),limits=c(0.06,120),"Spatial footprint [%]",
                     labels=c("0.1","1","10","100")) 
  # geom_encircle(data=compound,aes(x=footprint,y=ORtime,group=1), size = 5, color = "darkgrey",s_shape=1.3, expand=0.03,fill=NA,alpha=.6)+
  # geom_point(data=xtremdat,(aes(x=footprint,y=ORtime,color=lwt)), size=3,stroke=4,shape=1,alpha=.8) +
  # # scale_colour_gradientn(trans="log",colors=rev(rbPal2(100)),"Return period (Y)",breaks=c(1,10,100))+
  # geom_encircle(data=compoundX,aes(x=footprint,y=ORtime,group=1), size = 5, color = "red",s_shape=1.3, expand=0.03,fill=NA,alpha=.6)

plot(powow$lwt,powow$freq,type="h")

plot(rowow$lwt,rowow$freq,type="h")

plot(datowow$lwt,datowow$freq,type="h")


proprain<-c()
propwind<-c()
rmm<-c()
wmm<-c()
wti=c()
for(bx in 1:1485)
{
  location=locav$loc[bx]
  compbx<-sp03[which(sp03$loc==location),]
  comcr<-unique(compbx$ev.y)
  comcw<-unique(compbx$ev.x)
  wt<-compoundpls[which(!is.na(match(compoundpls$combin,compbx$combin))),]
  wtbg<-aggregate(list(ev=wt$combin),
                     by = list(lwt=wt$lwt),
                     FUN = function(x) c(length=length(x)))
  wtbg<- do.call(data.frame, wtbg)
  wtd=wtbg$lwt[which(wtbg$ev==max(wtbg$ev))]
  wti[bx]=wtd
  rainc<-metatest[which(metatest$loc==location),]
  cx<-compound[which(!is.na(match(compound$ev1,comcr))),]
  cxx<-cx[which(!is.na(match(cx$ev2,comcw))),]
  cy<-rainc$len.sum[na.omit(match(cxx$ev1,rainc$ev))]
  rmax<-cy/cxx$rf.max.cs
  # plot(rmax)
  rmm[bx]<-mean(rmax)
  winc<-metavDaf[which(metavDaf$loc==location),]
  cw<-winc$vi.max[na.omit(match(cxx$ev2,winc$ev))]
  
  wmax<-cw/cxx$wg.max.cs
  wmm[bx]<-mean(wmax)
  proprain[bx]<-length(comcr)/length(rainc$ev)
  propwind[bx]<-length(comcw)/length(winc$ev)
}

locav$prain<-proprain
locav$pwind=propwind
locav$rmm<-rmm
locav$wmm<-wmm
locav$lwt=wti

locav$fwind=locav$pwind*100
locav$frain=locav$prain*100
locav$frmm<-locav$rmm*100
locav$fwmm<-locav$wmm*100
locav$rwmm<-rmm*wmm*100
min(locav$fwmm)

ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill = "white", color = "gray10", size = 1) +
  theme_bw(16)+
  
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3)+
  geom_raster(data=locav,aes(x=lon,y=lat,fill = lwt,group=grp),binwidth=5,alpha=0.5) +
  scale_fill_manual(values=cols)+
  labs(y = "Latitude",x = "Longitude",fill= "LWT",size=20)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="italic"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(angle=0,size=14),
        legend.position = "right",
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  scale_y_continuous(
    breaks = c(48,50,52,54,56,58,60),limits = c(40,70),"Latitude")+
  scale_x_continuous(
    breaks =c(-6,-4,-2,0,2),limits=c(-10,10),"Longitiude") 

# Plot 2: Same plot with custom features
colors_border=c(rgb(0.8,0.2,0.5,0.9) ,  rgb(0.7,0.5,0.1,0.9), rgb(0.2,0.5,0.5,0.9) )
colors_in=c( rgb(0.8,0.2,0.5,0.4), rgb(0.7,0.5,0.1,0.4) , rgb(0.2,0.5,0.5,0.4))
datox<-datowow[,3]
pox<-powow[,3]
rox<-rowow[,3]
powow$lwt=as.numeric(powow$lwt)
prox<-full_join(powow,rowow,by="lwt")
prox<-full_join(prox,datowow,by="lwt")
datmort<-t(prox[,c(3,5,7)])
datamort<-as.data.frame(datmort)
colnames(datamort)=as.character(prox$lwt)
datamort<-rbind(rep(.4,27) , rep(0,27) , datamort)

radarchart(datamort, axistype=4 , 
           #custom polygon
           pcol=colors_border , pfcol=colors_in , plwd=1 , plty=1,pty=32,
           #custom the grid
           cglcol="darkgrey", cglty=1, axislabcol="darkgrey", caxislabels=seq(0,60,10),calcex=1.2, cglwd=1,seg=6,
           #custom labels
           vlcex=1.6 
)
legend(x=1.2, y=0.4, legend =c("MCE","all","CE"), bty = "n", pch=20 , col=colors_in , text.col = "black", cex=1.6, pt.cex=5, y.intersp = 0.5)

# Create data: note in High school for several students
set.seed(99)
data <- as.data.frame(matrix( sample( 0:20 , 15 , replace=F) , ncol=5))
colnames(data) <- c("math" , "english" , "biology" , "music" , "R-coding" )
rownames(data) <- paste("mister" , letters[1:3] , sep="-")

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
data <- rbind(rep(20,5) , rep(0,5) , data)



#############sensitivity analysis#####################
library(multisensi)
# NOT RUN {
## Test case : the Winter Wheat Dynamic Models (WWDM)
#  input factors design
data(biomasseX)
# input climate variable
data(Climat)
# output variables (precalculated to speed up the example)
data(biomasseY)

# to do dynsi process
# argument reduction=NULL
resD <- multisensi(design=biomasseX, model=biomasseY, reduction=NULL,
                   dimension=NULL, analysis=analysis.anoasg,
                   analysis.args=list(formula=2,keep.outputs = FALSE))
summary(resD)
