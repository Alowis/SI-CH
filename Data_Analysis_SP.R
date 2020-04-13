
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
  stat_density_2d(data=compoundfinalST,aes(x=lonMR,y=latMR,group=1,fill = stat(density)),geom = "raster", contour = FALSE,alpha=.6,n=300)


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

data(beavers)
plot(polarEventR[,c(1,3)])
p <- ggplot() + 
  geom_path(data=polarEventR, aes(x=monthstr,y=x,group=gr,colour=Year),size=1)+
  geom_line(data=polarEventR, aes(x=monthstr,y=x,colour = Year,group=Year),size=1) + 
  geom_point(data=polarEventR, aes(x=monthstr,y=x,colour = Year,group=Year),size=2)+
  geom_segment2(data=polarEventR,aes(x=monthstr, 
                                     xend=monthstr, 
                                     y=0, 
                                     yend=x),color="blue")+
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

matchingr<-unique(match(compoundfinalST$ev1,metaHax[[1]]$ev))
matchingw<-unique(match(compoundfinalST$ev2,metaHax[[2]]$ev))
rainevi<-metaHax[[1]][-matchingr,c(1,15,17)]
rainevi$c="r"
windevi<-metaHax[[2]][-matchingw,c(1,15,17)]
windevi$c="w"
compevi<-compoundfinalST[,c(2,22,24)]
compevi$c="c"
names(compevi)[1]="ev"
totevi<-rbind(compevi,rainevi,windevi)

bibcase<-rainevi

bibis1<-aggregate(bibcase$season,
                  by = list(Month = bibcase$x.month),
                  FUN = function(x) c(n =length(x),x=x[1]))
bibis1<- do.call(data.frame, bibis1)
bibis1$x.x[which(bibis1$Month==9)]=4
bibis1$c=1

bibcase<-windevi
bibis2<-aggregate(bibcase$season,
                  by = list(Month = bibcase$x.month),
                  FUN = function(x) c(n =length(x),x=x[1]))
bibis2<- do.call(data.frame, bibis2)
bibis2$x.x[which(bibis2$Month==9)]=4
bibis2$c=2

bibcase<-compevi
bibis3<-aggregate(bibcase$season,
                  by = list(Month = bibcase$x.month),
                  FUN = function(x) c(n =length(x),x=x[1]))
bibis3<- do.call(data.frame, bibis3)
bibis3$x.x[which(bibis3$Month==9)]=4
bibis3$c=3

bibis<-rbind(bibis1,bibis2,bibis3)
dodge <- position_dodge(width = 1)

ggplot(bibis,aes(x=Month, y=x.n,group=c))+geom_bar(stat = "identity", aes(fill=factor(c)),color="black",size=1,position = "dodge",alpha=.9) +
  scale_fill_manual(values = c( "blue", "red","orange","skyblue")) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c("Jan","Fev","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))+
  labs(y = "# of occurences",x = "Month",colour = "Parameter",fill= "Season",size=20)+
  theme_classic() +
  theme(axis.text.x = element_text(angle=0, vjust = 0.5,size=16))+
  theme(text = element_text(size=16))

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

###################Spatiotemporal plot
compoundfinalST$wg.max.s=compoundfinalS$wg.max.y
compoundfinalST$rf.max.s=compoundfinalS$rf.max.x

bestof<-metaHaf[which(metaHaf$viw.max>45),]
rbPal <- colorRampPalette(c('lightskyblue',"skyblue","gold","darkorange",'red',"purple"))
# ggplot(data=compoundfinalST,aes(x=spacescale,y=ORtime,colour=wg.max.s,size=rf.max.s))+
max(metaHax[[1]]$rf.surf)
ggplot(data=metaHaf[which(metaHaf$viw.max>10),],aes(x=viw.surf,y=x.dur,colour=viw.max,size=20))+
  geom_point(alpha=.5)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  geom_vline(xintercept = 1485,size=2,alpha=.5)+
  # geom_hline(yintercept=54)+
  # geom_vline(xintercept=1306)+
  scale_colour_gradientn(trans=scales::modulus_trans(2.5),colors=rbPal(100),"Max wind gust [m/s]",breaks=c(10,20,30,40,50,Inf),limits=c(15,51))+
  # scale_colour_manual(values = c( "blue", "red","orange","skyblue"))+
  # scale_size(trans=scales::modulus_trans(4),range=c(.1,10),"Max precipitation [mm]",breaks=c(10,20,30,40,50,Inf),limits=c(1,51),guide = guide_legend(override.aes = list(colour = alpha("grey",.8))))+
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
