####Second script which creates compound hazard clusters######
####1. Import outputs from 01_PC
####2. Estimate spatial and temporal overlaps of wind and precipitation clusters
####3. Save results

rm(list=ls()) 
gc()
getwd()
###################Choose working directory#####################
# setwd("C:/Users/PhD Student/OneDrive - King's College London/DATA/04_Spatiotemporal")
# setwd("C:/Users/k1638615/King's College London/OneDrive - King's College London/DATA/04_Spatiotemporal")

setwd("C:/Users/Alois/OneDrive - King's College London/DATA/04_Spatiotemporal")

config_file=paste0(getwd(),"/config/config_general.R")
source(config_file)




############Load outputs from 01_preprocess_clustering##############

load(file="out/Wind_stfprint2.Rdata")
load(file="out/Rain_stfprint.Rdata")
load(file="out/RainEv_metap_1979-2019.Rdata")
metaHar=metaHaf
load(file="out/WindEv_meta_1979-2019.Rdata")


load(file="out/RainEv_hdat_1979-2019.Rdata")
load(file="out/RainEv_ldat_1979-2019.Rdata")
metavDar=metavDaf
metavHar=metavHaf
load(file="out/WindEv_hdat_1979-2019.Rdata")
load(file="out/WindEv_ldat_1979-2019.Rdata")

library(nnet)

###################Quick event visualizer#################


#choose rain event to display
evp<-metaHar$ev
evento<-metavHar[which(metavHar$cluster==3802),]

#choose wind event to display
cps<-metaHaf$ev
maxeventwx<-metavHaf[which(metavHaf$cluster==1323),]

mapUK = SpacializedMap(database="world",regions = c("UK","France","Spain","Portugal","Italy","Ireland","Belgium","Netherland"))


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
library(ggalt)

ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill = "transparent", color = "gray10", size = 1) +
  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
  # geom_encircle(aes(x=Var1, y=Var2,group=cluster),
  #               data = evento, size = 2, color = "blue",s_shape=.8, expand=0,fill="blue",alpha=.1)+
  geom_bin2d(data=evento,aes(x=Var1,y=Var2,group=cluster,weight=vecmeta),alpha=.7,binwidth=c(.25,.25))+
  # stat_density_2d(data=evento,geom = "polygon", aes(x=Var1,y=Var2,group=cluster,fill=..level..),alpha=.3)+
  scale_fill_viridis_c()+
  geom_bin2d(data=maxeventwx,aes(x=Var1,y=Var2,group=cluster,weight=vecmeta),alpha=.7,binwidth=c(.25,.25))+
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
    breaks =c(-6,-4,-2,0,2),limits=c(-10,10),"Longitude")







#these data frame contains spatial id and 
#temporal id over the duration of an event of each cell involved in that event. 
#(e.g., if the event last 10 hours, each cell involved will be repeated 10 times even if the the cell is only impacted during 1 hour)


#Inner join in R:  Return only the rows in which the left table have matching keys in the right table
#this one set a pre-filter on events which have temporal overlap

sp03<-inner_join(nwo3,nvo3,by=c("loc","dd"))

#aggregation by events and space, each row correpond to one cell involved in both rain and wind event

tesp<-aggregate(list(len=sp03[,3]),
                by = list(ev2=sp03[,1],ev1 = sp03[,4],loc=sp03[,2]),
                FUN = function(x) c(length=length(unique(x))))
tesp <- do.call(data.frame, tesp)

sesp<-aggregate(list(len=sp03[,2]),
                by = list(ev2=sp03[,1],ev1 = sp03[,4],dd=sp03[,3]),
                FUN = function(x) c(length=length(x)))
sesp <- do.call(data.frame, sesp)

###################
# Extract more metadata about the events 
# Spatial value (for combined  duration of the event), 
# Spatiotemporal value (only during the compound duration and in compound space), 
# Temporal (outside the compound space inside the compound time))
#Here are computed:
#1) the spatial compound (all time inside the compound space)
#2) the temporal compound (all space inside the compound time)
#3) the spatiotemporal compound ( inside the compound space and time)


#########compound clusters based on precipitation clusters##########

metatest1<-metavDar
names(metatest1)[c(1)]=c("ev1")

rrainacc<-c()
olei<-c()
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
  rrainacc<-rbind(rrainacc,rainacc)
  olei<-rbind(olei, olel)
  #last join

}


metatest<-rrainacc
names(metatest)[c(1)]=c("ev1")
sptotal=F

#1) Spatial compound
metatestt<-join(metatest,metatest1,by=c("ev1","loc"))
metatestt<-metatestt[,-c(6,7,8)]
SCompR<-inner_join(metatestt,tesp,by=c("ev1","loc"))

cevcx=unique(SCompR$cev[order(SCompR$len.sum,decreasing=T)])
cevc=cevcx[13]
oula<-SCompR[which(SCompR$cev==cevc),]
ziv<-SCompR$len.sum[which(SCompR$cev==cevc)]
plot(ziv[order(ziv)])


spSCR<-aggregate(list(rf=SCompR[,4],wg=SCompR[,6]),
                 by = list(ev2=SCompR[,6],ev1 = SCompR[,1]),
                 FUN = function(x) c(length=length((x)),max=max(x),mean=mean(x)))
spSCR <- do.call(data.frame, spSCR)


SCompR$cev=paste(SCompR$ev1,SCompR$ev2)
for(rg in 1:length(spSCR$ev2)){
  lesev<-paste(spSCR[rg,2],spSCR[rg,1])
  mr<-spSCR$rf.max[rg]
  spSCR$maxlocsR[rg]<-as.character(SCompR$loc[which(SCompR$cev==lesev & SCompR$len.sum==mr)])
}



#3) Spatiotemporal compound

names(sp03)<-c("ev2","loc","dd","ev1")

spkik<-sp03 %>% group_by(ev1,dd,loc) %>% mutate(id = row_number())
metatest<-metavHar
names(metatest)[c(5,6,8)]=c("ev1","dd","loc")
names(olei)[5]="ev1"
metatestx<-full_join(metatest,olei,by=c("ev1","dd","loc"))
sp20<-inner_join(metatestx,spkik,by=c("ev1","loc","dd"))
spSTCR<-sp20[which(sp20$id==1),]

spSTCR2<-aggregate(list(rf=spSTCR[,10]),
                   by = list(ev2=spSTCR[,11],ev1 = spSTCR[,5],loc=spSTCR[,8]),
                   FUN = function(x) c(length=length(x),max=max(x,na.rm=T),sum=sum(x,na.rm=T)))
spSTCR2 <- do.call(data.frame, spSTCR2)


spSTFR<-aggregate(list(rf=spSTCR2[,6]), 
                  by = list(ev2=spSTCR2[,1],ev1=spSTCR2[,2]),  
                  FUN = function(x) c(length=length(x),max=max(x),mean=mean(x))) 
spSTFR<- do.call(data.frame, spSTFR)

spSTCR2$cev=paste(spSTCR2$ev1,spSTCR2$ev2)
spSTFR$maxlocR="00 00"
for(rg in 1:length(spSTFR$ev2)){
  lesev<-paste(spSTFR[rg,2],spSTFR[rg,1])
  mr<-spSTFR$rf.max[rg]
  rtime=as.character(spSTCR2$loc[which(spSTCR2$cev==lesev & spSTCR2$rf.sum==mr)])
  if(length(rtime)>0) spSTFR$maxlocstR[rg]<-rtime
} 
if(sptotal==T){
#2) Temporal compound
  

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
}
rm(metatestx,metatest,olei,olel)
gc()

####################compound clusters based on wind gust clusters############ 
metatest<-metavHaf
names(metatest)[c(5,6,8)]=c("ev2","dd","loc")

#1) Spatial compound

SCompW<-inner_join(metatest,tesp,by=c("ev2","loc"))



spSCW<-aggregate(list(wg=SCompW[,4]),
                 by = list(ev2=SCompW[,9],ev1 = SCompW[,5],loc=SCompW[,8]),
                 FUN = function(x) c(length=length((x)),max=max(x),sum=sum(x)))
spSCW <- do.call(data.frame, spSCW)

spSCW2<-aggregate(list(wg=spSCW[,5]), 
                  by = list(ev2=spSCW[,1],ev1=spSCW[,2]),  
                  FUN = function(x) c(length=length(x),max=max(x),mean=mean(x))) 
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

spSTCW2<-aggregate(list(wg=spSTCW[,4]),
                   by = list(ev2=spSTCW[,5],ev1 = spSTCW[,9],loc=spSTCW[,8]),
                   FUN = function(x) c(length=length((x)),max=max(x),sum=sum(x)))
spSTCW2 <- do.call(data.frame, spSTCW2)


spSTFW<-aggregate(list(wg=spSTCW2[,5]), 
                  by = list(ev2=spSTCW2[,1],ev1=spSTCW2[,2]),  
                  FUN = function(x) c(length=length(x),max=max(x))) 
spSTFW<- do.call(data.frame, spSTFW)

spSTCW2$cev=paste(spSTCW2$ev1,spSTCW2$ev2)
for(rg in 1:length(spSTFW$ev2)){
  lesev<-paste(spSTFW[rg,2],spSTFW[rg,1])
  mr<-spSTFW$wg.max[rg]
  spSTFW$maxlocstW[rg]<-as.character(spSTCW2$loc[which(spSTCW2$cev==lesev & spSTCW2$wg.max==mr)])
}
if(sptotal==T){
  
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


}


#################Creation of compound clusters database#############

if(sptotal==F){
  spSTFR<-spSCR
  spSTFW<-spSCW2
  names(spSTFW)[c(1,2)]=c("ev1","ev2")
}
spSTFR$gr="r"
spSTFW$gr="w"
spSTF<-full_join(spSTFR,spSTFW,by=c("ev1","ev2"))
max(spSCR$rf.max)
max(metaHar$len.max)

plot(spSTF$wg.max.y,spSTF$rf.max)
for (mm in 1:length(spSTF$ev2)){
  spSTF$rf.max.cs[mm]=max(spSTF$rf.max[mm],na.rm=T)
  spSTF$wg.max.cs[mm]=max(spSTF$wg.max.y[mm],na.rm=T)
  spSTF$rf.mean.cs[mm]=max(spSTF$rf.mean[mm],na.rm=T)
  spSTF$wg.mean.cs[mm]=max(spSTF$wg.mean.y[mm],na.rm=T)
  
}

plot(spSTF$wg.mean.cs,spSTF$rf.mean.cs)

spSCW2<-spSCW2[,-3]
spTCW2<-spTCW2[,-3]

spTCR2<-spTCR2[,-3]


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

rm(metatest, spkik)
gc()
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
  idc<-which(!is.na(match(metavHar$cluster,kev))) 
  inev<-metavHar$time[idc[1]]
  endev<-metavHar$time[idc[length(idc)]]
  spTiming$startimeRain[kv]=inev
  spTiming$endtimeRain[kv]=endev
} 

for (kv in 1:length(spTiming$ev2)){
  kev<-spTiming$ev2[kv]
  idc<-which(!is.na(match(metavHaf$cluster,kev))) 
  inev<-metavHaf$time[idc[1]]
  endev<-metavHaf$time[idc[length(idc)]]
  spTiming$startimeWind[kv]=inev
  spTiming$endtimeWind[kv]=endev
}


spTiming$startime<-as_datetime(c(spTiming$startime),origin="1970-01-01")
spTiming$endtime<-as_datetime(c(spTiming$endtime),origin="1970-01-01")
spTiming$startimeWind<-as_datetime(c(spTiming$startimeWind),origin="1970-01-01")
spTiming$endtimeWind<-as_datetime(c(spTiming$endtimeWind),origin="1970-01-01")
spTiming$startimeRain<-as_datetime(c(spTiming$startimeRain),origin="1970-01-01")
spTiming$endtimeRain<-as_datetime(c(spTiming$endtimeRain),origin="1970-01-01")

spSTP<-full_join(spTiming,spSTF,by=c("ev1","ev2"),all=T)


names(metaHax[[1]])
allez<-metaHax[[1]][,c(1,27)][which(!is.na(match(metaHax[[1]]$ev,spSTP$ev1))),]
lesbleu<-metaHax[[2]][,c(1,8)][which(!is.na(match(metaHax[[2]]$ev,spSTP$ev2))),]
names(lesbleu)<-c("ev2","wg.max.a")
names(allez)<-c("ev1","rf.max.a")
spSTP<-inner_join(spSTP,lesbleu,by=c("ev2"),all=T)
spSTP<-inner_join(spSTP,allez,by=c("ev1"),all=T)



compound<-inner_join(spSTP,tscale,by=c("ev1","ev2","combin"))
for (ti in 1:length(compound$ev2)){
  compound$ORtime[ti]<-difftime(max(compound$endtimeRain[ti],compound$endtimeWind[ti]),min(compound$startimeRain[ti],compound$startimeWind[ti]),unit="hours")+1}

compound$rf.mean.cs= spSTF$rf.mean.cs


compound$rf.mean.h=compound$rf.mean.cs/compound$ORtime
compound$rf.max.h=compound$rf.max.cs/compound$ORtime
plot(compound$ORtime,compound$rf.max.h)

#Compound clusters database
save(compound,file="out/CompoundRW_79-19.v2x.Rdata")
write.csv(compoundfinalST,file="out/CompoundRW_79-19.v3.csv")

comprain<-na.omit(match(metaHar$ev,compound$ev1))
comprain2<-na.omit(match(compound$ev1,metaHar$ev))
cocor<-metaHar[comprain,]
cocor2<-metaHar[comprain2,]
cocor3<-data.frame(compound,cocor2)

compwind<-na.omit(match(metaHaf$ev,compound$ev2))
compwind2<-na.omit(match(compound$ev2,metaHaf$ev))
cocow<-metaHaf[compwind,]
cocow2<-metaHaf[compwind2,]

names(cocow2)
cocow3<-cocow2
colnames(cocow3)


compoundfinalST<-data.frame(cocor3,cocow3)
compoundfinalST<-inner_join(compoundfinalST,tscale,by=c("ev1","ev2","combin"))

compoundBG<-compoundfinalST[c(1:9,37,39,23,25,46,52,54,24,26,49,60,61)]
compoundBG$ORspace=compoundBG$vir.surf+compoundBG$viw.surf-compoundBG$spacescale.y

for (ti in 1:length(compoundBG$ev2)){
  compoundBG$ORtime[ti]<-difftime(max(compoundBG$endtimeRain[ti],compoundBG$endtimeWind[ti]),min(compoundBG$startimeRain[ti],compoundBG$startimeWind[ti]),unit="hours")+1}

write.csv(compoundBG,file="out/ch5.compoundclusters.csv")

names(spSCW2)[c(1,2)]<-c("ev1","ev2")
compoundfinalS<-inner_join(spSCR,spSCW2,by=c("ev1","ev2"))
names(spTCW2)[c(1,2)]<-c("ev1","ev2")
compoundfinalT<-inner_join(spTCR2,spTCW2,by=c("ev1","ev2"))


names(compoundfinalST)[c(28,29,42,44)]<-c("rf.space","rf.time","wg.space","wg.time")
compoundfinalST$ORspace=compoundfinalST$vir.surf+compoundfinalST$wg.space-compoundfinalST$spacescale.x
for (ti in 1:length(compoundfinalST$ev2)){
  compoundfinalST$ORtime[ti]<-difftime(max(compoundfinalST$endtimeRain[ti],compoundfinalST$endtimeWind[ti]),min(compoundfinalST$startimeRain[ti],compoundfinalST$startimeWind[ti]),unit="hours")+1}


rm(cocor,cocor2,cocor3,cocow,cocow2,cocow3,spSTCR,spSTCR2,spSTCW,spSTCW2,spSTFW,spSTP,spSCR,spSCW,spSCW2,spkik,sp20,spSTF,spSTFR,samptt)
gc()

save (compoundfinalST,file="messycompound_79-19.Rdata")


compoundfinalST$ratiospace=compoundfinalST$spacescale/compoundfinalST$ORspace

compoundfinalST$ratiotime=compoundfinalST$timescale/compoundfinalST$ORtime



###########################Supplementary analysis####################
moreana=F
if(moreana==T){
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
}

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


rm(lonlatime)
gc()
if(moreana==T){
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
  
  ggplot(compoundfinalST, aes(x=timescale,y=..density..,fill=season)) + 
    geom_density(alpha=.6,size=.4,position = "fill",adjust=2)
  geom_histogram(aes(x=geodistcs,y=..density..),binwidth=20,color="red",fill="orange",position="dodge", alpha=0.4)+
    geom_density(aes(x=geodistcs,y=..density..),alpha=.1,fill="pink",size=.4)
  
  
  coord_radar <- function (theta = "x", start = 0, direction = 1) {
    theta <- match.arg(theta, c("x", "y"))
    r <- if (theta == "x") "y" else "x"
    ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start, 
            direction = sign(direction),
            is_linear = function(coord) TRUE)
  }
  
  ggplot(metaHax[[1]], aes(x=log(x.dur),group=season,fill=as.character(season)))+
    geom_density(adjust=2,position = "fill",alpha=.6,n=500)+
    scale_fill_brewer(palette="Accent")
  
  windst<-metaHax[[2]][,c(6,14)]
  windst<-metaHax[[2]][,c(6,14)]

  
  
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
}
