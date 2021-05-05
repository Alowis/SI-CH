library(sf)
library(lubridate)
library(vistime)
rm(list=ls())  
getwd()


setwd("C:/Users/PhD Student/OneDrive - King's College London/DATA/04_Spatiotemporal")
setwd("C:/Users/k1638615/King's College London/OneDrive - King's College London/DATA/04_Spatiotemporal")
setwd("C:/Users/Alois/OneDrive - King's College London/DATA/04_Spatiotemporal")
config_file=paste0(getwd(),"/config/config_general.R")
source(config_file)


#load boundaries of NUTS1 regions in the UK
uk_reg <- st_read("in/NUTS_Level_1__January_2018__Boundaries.shp")
areaval<-data.frame(uk_reg$nuts118nm, uk_reg$st_areasha/1000000)
class(uk_reg)
uk_reg=uk_reg[-12,]
uk_reg$nuts118nm
sf<- readOGR( "in/NUTS_Level_1__January_2018__Boundaries.shp")

ukg<-st_transform(uk_reg, "+proj=utm +zone=30U, +datum=WGS84")
class(ukg)


mapUK = SpacializedMap(regions = c("UK","France","Spain","Portugal","Italy","Ireland"))


ukk<- gSimplify(mapUK, tol=0.001, topologyPreserve=TRUE)
plot(ukk)
uk_fort <- ggplot2::fortify(ukk)
longlims=c(-6, 2)
longdom=seq(longlims[1],longlims[2],by=.25)
latlims=c(48,59)
latdom=seq(latlims[1],latlims[2],by=.25)

#############VALIDATION OF THE CLUSTERING METHOD##################

#identifies temporal (days) and spatial matches 
#(NUTS1 regions) between catalog and cluster databases

load("out/CompoundRW_79-19.v2x.Rdata")
load(file="out/messycompound_79-19.Rdata")

#days with an event in the catalog happening
validationdays<-read.csv(file=paste0(getwd(),"/csvs/identifix.csv"),header = T)
stardate=as.Date(validationdays$Startdate,format='%d/%m/%Y')
endate=as.Date(validationdays$Enddate,format='%d/%m/%Y')
stardate-endate

totsek=data.frame(lasek=integer(), id=integer())
for(dd in 1: length(stardate)){
  sd=stardate[dd]
  ed=endate[dd]
  lasek<-seq(sd, ed, by="days")
  lasek=as.data.frame(lasek)
  lasek$id=dd
  totsek=rbind(totsek,lasek)
}

totsek$type=validationdays$Dominant.hazard[match(totsek$id,validationdays$ï..Dis.No)] 
validationdays$id=validationdays$ï..Dis.No

validationdays$startime=as.POSIXct(stardate)

#data visualisation
ggplot( data= compound, aes(x= startime)) +
  geom_histogram( color="grey", alpha=0.4, position = 'identity',bins=2000) +
  geom_rug(data=validationdays, aes(x = startime, y = 0), position = position_jitter(height = 0))


stseason<-data.frame(compound)

load(file="out/RainEv_metap2_1979-2019.Rdata")
load(file="out/Rainev_hdatp_1979-2019.Rdata")

load(file="out/WindEv_meta_1979-2019.Rdata")
stseason$startimeD<-as.Date(stseason$startime,format='%d/%m/%Y')
stseason$endtimeD<-as.Date(stseason$endtime,format='%d/%m/%Y')
stclust<-c()
tipi<-c()

# load(file="out/Wind_stfprint2.Rdata")
# load(file="out/Rain_stfprint.Rdata")
load(file= "out/WindEv_ldat_1979-2019.Rdata")
wmaxloc<-c()
for (idev in metaHaf$ev){
  wacc<-metavDaf[which(metavDaf$ev==idev),]
  wacev<-metaHaf[which(metaHaf$ev==idev),]
  locmax=which(wacc$vi.max==wacev$viw.max)
  if(length(locmax)>1)locmax=locmax[1]
  locas=as.character(wacc$loc[locmax])
  wmaxloc=c(wmaxloc,locas)
}
metaHaf$locmax=wmaxloc
names(metaHar)[16]="locmax"




longlims=c(-6, 2)
longdom=seq(longlims[1],longlims[2],by=.25)
latlims=c(48,59)
latdom=seq(latlims[1],latlims[2],by=.25)
lagrid<-expand.grid(longdom,latdom)
lagrid$loc=paste(lagrid$Var1,lagrid$Var2,sep=" ")
lagrid.coord<-lagrid
coordinates(lagrid.coord) <- c("Var1", "Var2")


proj4string(lagrid.coord) <- CRS("+proj=longlat +datum=WGS84")

lagrid.utm <- spTransform(lagrid.coord, CRS("+proj=utm +zone=30U, +datum=WGS84"))

sf.utm <- spTransform(sf, CRS("+proj=utm +zone=30U, +datum=WGS84"))
sf.utm@proj4string
mamamia<-over(lagrid.utm,sf.utm)
mamamia$nuts118cd[which(mamamia$nuts118cd=="UKN")]=NA
mam<-cbind(lagrid,mamamia)


########Spationtemporal match################
returnout<-c()
for (idiwant in 1:length(validationdays$Year) ){
  
  print(idiwant)
  # idiwant=100
  laval<-validationdays[which(validationdays$id==idiwant),]
  sekance=totsek[which(totsek$id==idiwant),]

  if (laval$Dominant.hazard=="Extreme rainfall") hazevent<-metaHar
  if (laval$Dominant.hazard=="Extreme wind") hazevent<-metaHaf
  
  hazevent$startimeD<-as.Date(hazevent$startev,format='%d/%m/%Y')
  hazevent$endtimeD<-as.Date(hazevent$endev,format='%d/%m/%Y')
  
  totsev=data.frame(lasev=integer(), id=integer())
  for(dd in 1: length(hazevent$ev)){
    sd=hazevent$startimeD[dd]
    ed=hazevent$endtimeD[dd]
    lasev<-seq(sd, ed, by="days")
    lasev=as.data.frame(lasev)
    lasev$id=hazevent$ev[dd]
    totsev=rbind(totsev,lasev)
  }
  
  chosenones<-which(!is.na(match(totsev$lasev,sekance[,1])))
  

  stconfx<-totsev[chosenones,]
  stconfr<-hazevent[which(!is.na(match(hazevent$ev,stconfx$id))),]
  # stconfr$daymatch=stconfx$lasev
  
  
  
    rev<-unique(stconfx$id)
    wev=rev

  if (laval$Dominant.hazard=="Extreme rainfall") {

    ppr<-metatest[which(!is.na(match(metatest$ev1,rev))),]
    # eventa<-metavHar[which(!is.na(match(metavHar$cluster,rev))),]
  
  # ppr<-aggregate(list(ri=eventa[,4]),
  #                by = list(loc=eventa[,8]),
  #                FUN = function(x) c(length=length(x)))
  # ppr <- do.call(data.frame, ppr)
    if(length(ppr$ev)>0){
      px=strsplit(as.character(ppr$loc)," ")
      px<-as.data.frame(matrix(unlist(px),ncol=2,byrow=T))
      ppr$lon<-as.numeric((px[,1]))
      ppr$lat<-as.numeric((px[,2]))
      ppr$clust=ppr$ev1
    }
  
    }
    
  if (laval$Dominant.hazard=="Extreme wind") {

  ppr<-metavDaf[which(!is.na(match(metavDaf$ev,wev))),]
  
  # ppr<-aggregate(list(ri=eventi[,4]),
  #                by = list(loc=eventi[,8]),
  #                FUN = function(x) c(length=length(x),max=max(x)))
  # ppr <- do.call(data.frame, ppr)
  if(length(ppr$ev)>0){
  Rw=strsplit(as.character(ppr$loc)," ")
  Rwi<-as.data.frame(matrix(unlist(Rw),ncol=2,byrow=T))
  ppr$lon<-as.numeric((Rwi[,1]))
  ppr$lat<-as.numeric((Rwi[,2]))
  ppr$clust=ppr$ev
  }
  }
  
  god<-mam[which(!is.na(match(mam$loc,ppr$loc))),]
  god1<-god[which(!is.na(god$objectid)),]
  
  god2<-inner_join(ppr,god1,by=c("loc"))
  
  
  unique(god2$ev1)
  reg=laval[,10:20]
  regf<-names(reg[which(reg==1)])
  clf<-unique(god2$nuts118cd)
  guys<-regf[which(!is.na(match(regf,clf)))]
  # print(guys)
  rat=length(guys)/length(regf)
  
  nclust=god2[which(!is.na(match(god2$nuts118cd,guys))),]
  ncl<-unique(nclust$clust)
  cwap<-god2[which(!is.na(match(god2$clust,ncl))),]
  ptn<-unique(cwap$nuts118cd)
  ptd<-unique(nclust$nuts118cd)
  
   # 
   # ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
   #  geom_polygon(fill = "transparent", color = "gray10", size = 1) +
   #  coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
   #  geom_raster(data=ppr,aes(x=lon,y=lat,group=clust,fill=len.sum),alpha=.7)+
   #  scale_fill_gradientn(trans=scales::modulus_trans(1),colours = rgb.palette(12)," ",breaks=breaks_extended(10),guide= "colorbar")+
   #  # scale_fill_manual(values=rgb.palette(12))+
   #  theme(axis.text=element_text(size=16),
   #        axis.title=element_text(size=18,face="italic"),
   #        panel.background = element_rect(fill = "white", colour = "grey50"),
   #        legend.title = element_text(size=18),
   #        legend.text = element_text(size=14),
   #        legend.key = element_rect(fill = "transparent", colour = "transparent"),
   #        legend.key.size = unit(1, "cm"))+
   #  scale_y_continuous(
   #    breaks = c(48,50,52,54,56,58,60),limits = c(40,70),"Latitude")+
   #  scale_x_continuous(
   #    breaks =c(-6,-4,-2,0,2),limits=c(-10,10),"Longitiude")
   # 
  # print(ncl)
  merd<-paste(ncl,collapse=",")
  outi<-c(idiwant, paste(ncl,collapse=","), length(ncl),paste(guys, collapse=","),rat,length(ptn),length(clf),length(regf))
  
  returnout<-rbind(returnout,outi)
}

retout<-as.data.frame(returnout)
names(retout)= c("ID_obs","ID_mod","N_ev","regions","spatial_match","N_r_mod","N_rt_mod","N_r_obs")
retout$N_r_mod[which(retout$N_r_mod==12)]=11
retout$N_rt_mod[which(retout$N_rt_mod==12)]=11

ze<-as.numeric(retout$N_r_obs)
zr<-as.numeric(retout$N_r_mod)
cor(ze,zr)
ze=data.frame(ze)
ze$scale="regional"
ze$scale[which(ze$ze<3)]="local"
ze$scale[which(ze$ze>6)]="multi-regional"
zr=data.frame(zr)
zr$scale="regional"
zr$scale[which(zr$zr<3)]="local"
zr$scale[which(zr$zr>6)]="multi-regional"
zr$scale[which(zr$scale==ze$scale)]
length(zr$zr[which(zr$zr>ze$ze)])




#save the data frame with the severe events caltalog
save(retout,file="ident_events_v5.Rdata")


#############Now plot the results#################
load(file="ident_events_v5.Rdata")

ultim<-cbind(validationdays,retout)
retout$N_ev=as.numeric(retout$N_ev)
mean(retout$N_ev[which(retout$N_ev>0)])
ptin<-ultim[which(ultim$Dominant.hazard=="Extreme rainfall"),]
ratio=1-length(ptin$Year[which(ptin$ID_mod=="")])/length(ptin$Year)
ccm<-ultim[which(ultim$N_ev==1),]

oups<-ultim$Dominant.hazard[which(ultim$ID_mod=="")]
oups
m1<-month(ultim$Startdate[which(ultim$ID_mod=="")])
m2<-month(ultim$Startdate[which(ultim$ID_mod!="")])
hist(m1,breaks=12)
hist(m2,breaks=12)

metad=c()
for (ul in 1:11){
ola<-ultim[which(ultim[,(9+ul)]==1),]
ouy<-length(ola$Startdate[which(ola$ID_mod=="")])
oc=length(ola$Startdate)
oc2=c(oc,ouy)
metad=rbind(metad,oc2)
}

metad=as.data.frame(metad)
metad$rati=1-metad$V2/metad$V1

#Overlap of grid cells and NUTS1 regions
rgb.palette=colorRampPalette(c("lightblue","royalblue","darkolivegreen3","gold","orange","red","darkred")
                               ,interpolate="linear",bias=1)
  
ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
    geom_polygon(fill = "transparent", color = "gray10", size = 1) +
    coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3) +
    geom_raster(data=mam,aes(x=Var1,y=Var2,group=objectid,fill=nuts118cd),alpha=.7)+
    scale_fill_manual(values=rgb.palette(12))+
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

metad$objectid=c(1,2,3,4,5,6,7,8,9,10,11)
ukg2<-left_join(ukg,metad, by="objectid")

hist(as.numeric(ultim$N_ev),breaks=20)

dur=as.Date(ultim$Enddate,format='%d/%m/%Y')-as.Date(ultim$Startdate,format='%d/%m/%Y')

plot(dur,as.numeric(ultim$N_ev))
length(which(ultim$N_ev==1))
class(ukg2)
class(sf.utm)

ukg2$V1=as.numeric(ukg2$V1)
rgb.palette=colorRampPalette(rev(c("#d53e4f","#fc8d59","#fee08b","#ffffbf","#e6f598",
  "#99d594","#3288bd")),interpolate="linear",bias=1)

#Plots showing number of events per region and hit rate per region

#Figure 5.11
ggplot(data = ukg2) + 
  geom_sf(aes(fill = V1)) +
  scale_fill_gradientn(colours = rgb.palette(20)," ",
                       guide= "colourbar",breaks = c(30,40,50,60,70),limits=c(20,70))+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="italic"),
        panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14)) +
  ggtitle("Events per NUTS1 Region") + 
  coord_sf()

ukg2$perc=ukg2$rati*100

#Figure 5.12
ggplot(data = ukg2) + 
  geom_sf(aes(fill = perc)) +
  scale_fill_gradient2(guide= "colourbar",low="khaki",mid ="#addd8e", high = "#31a354",na.value="white",midpoint=95,"%")+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="italic"),
        panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14)) +
  ggtitle("Percentage of events detected") + 
  coord_sf()

qltim=ultim[which(ultim$Dominant.hazard=="Extreme wind"),]
windid=as.vector(qltim$ID_mod)
wix<-  strsplit(windid, ",")
isc=c()
for(idw in 1:length(wix)){
  bof=wix[[idw]]
  wcom=F
  iscom=(match(as.numeric(bof),compound$ev2))
  if(length(iscom)>1)wcom=T
  if(length(iscom)<1)wcom=F
  isc=c(isc,wcom)
  
}
length(isc[which(isc==T)])
qltim$iscompound=isc

rltim=ultim[which(ultim$Dominant.hazard=="Extreme rainfall"),]
rainid=as.vector(rltim$ID_mod)
rix<-  strsplit(rainid, ",")
irc=c()
for(idw in 1:length(rix)){
  bof=rix[[idw]]
  wcom=F
  iscom=(match(as.numeric(bof),compound$ev1))
  if(length(iscom)>1)wcom=T
  if(length(iscom)<1)wcom=F
  irc=c(irc,wcom)
  
}
length(irc[which(irc==T)])

rltim$iscompound=irc
  


alultim=rbind(rltim,qltim)


alultim$start=as.Date(alultim$Startdate,format='%d/%m/%Y')
alultim$end=as.Date(alultim$Enddate,format='%d/%m/%Y')
alultim$event=alultim$ï..Dis.No
alultim$event_type="rain"
alultim$event_type[which(alultim$Dominant.hazard=="Extreme wind")]="wind"
ulticom=alultim[which(alultim$iscompound==T),]

positions <- rep(c(0.5, -0.5, 1.0, -1.0, 1.25, -1.25, 1.5, -1.5),20)
positions<-positions[-c(1:3)]
alultim$position=positions



#########Spatial and temporal locations of the 151 severe events from the catalog#############
ev.tot=c() 
for (ite in 1:length(alultim$Year)){
  mom=alultim[ite,]
  regl<-which(mom[,10:20]==1)
  ev.rep=rep(mom$Startdate,length(regl))
  ev.tip=rep(mom$event_type,length(regl))
  ev.com=rep(mom$iscompound,length(regl))
  ev.repi=cbind(ev.rep,regl,ev.tip,ev.com)
  ev.tot=rbind(ev.tot,ev.repi)
}

ev.tot=as.data.frame(ev.tot)
ev.tot$ev.rep=as.Date(ev.tot$ev.rep,format='%d/%m/%Y')
ev.tot$regl=as.numeric(ev.tot$regl)


ev.comp=ev.tot[which(ev.tot$ev.com==T),]
insert_minor <- function(major_labs, n_minor) {
  labs <-c( sapply( major_labs, function(x) c(x, rep("", 9) ) ) )
  labs[1:(length(labs)-n_minor)]}

start_date=as.Date("01/01/1980",format='%d/%m/%Y')
end_date=as.Date("01/01/2020",format='%d/%m/%Y')
date_br10 <- seq(from = start_date, to = end_date, by = "10 years")
date_br5 <- seq(from = start_date, to = end_date, by = "1 year")
bra<-c(insert_minor(format(date_br10, "%Y"), 9))

ev.tot$pos=ev.tot$regl
ev.tot$pos[which(ev.tot$ev.tip=="rain")]=ev.tot$pos[which(ev.tot$ev.tip=="rain")]+0.2
ev.tot$pos[which(ev.tot$ev.tip=="wind")]=ev.tot$pos[which(ev.tot$ev.tip=="wind")]-0.2
ev.comp$pos=ev.comp$regl

locs=c("North East","North West","Yorkshire","East Midlands","West Midlands","East of England",
       "London","South East","South West","Wales","Scotland")


#Figure 5.10
ggplot( data= ev.tot, aes(x= ev.rep,y=pos,shape=ev.tip,color=ev.tip)) +
  # geom_histogram( color="black", alpha=0.4, position = 'identity') +
  # geom_density(aes(x = startev),alpha=.4) +
  geom_point(size=3,alpha=1,stroke=2)+
  geom_point(data=ev.comp,size=5,alpha=.5,stroke=1.5,shape=13,color="darkgreen",fill="darkgreen")+
  scale_shape_manual(values=c("rain"=1,"wind"=4))+
  theme(axis.text=element_text(size=16),
        axis.text.y=element_text(size=16),
        axis.title=element_text(size=18,face="italic"),
        panel.grid.major.x = element_blank() ,
        panel.border = element_rect(size=1.2, fill = NA),
        # panel.grid.minor.x = element_line( size=.8, color="black" ) ,
        panel.grid.major.y = element_line( size=1, color="black" ), 
        panel.background = element_rect(fill = "transparent",color="black"),
        legend.title = element_blank(),
        legend.text = element_blank()) +
  scale_color_manual(values=c("rain"="royalblue","wind"="orange"),"dominant hazard")+
  scale_x_date(breaks= date_br5,labels = bra, "Date")+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11),labels=locs,"UK NUTS1 Regions")
