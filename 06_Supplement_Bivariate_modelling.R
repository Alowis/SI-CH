
rm(list=ls())  
getwd()

setwd("C:/Users/PhD Student/OneDrive - King's College London/DATA/04_Spatiotemporal")
setwd("C:/Users/k1638615/King's College London/OneDrive - King's College London/DATA/04_Spatiotemporal")
setwd("C:/Users/Alois/OneDrive - King's College London/DATA/04_Spatiotemporal")
config_file=paste0(getwd(),"/config/config_general.R")
source(config_file)

functions=paste0(getwd(),"/spatiotemporal_clustering/function/Functions_Appli-Script.R")
source(functions)

JointExcureve_file=paste0(getwd(),"/spatiotemporal_clustering/function/Joint_excurves.R")
source (JointExcureve_file)

interh<-"comb" #"comb or "casc"
upobj=0.001


load(file="out/CompoundRW_79-19.v2x.Rdata")
load(file="out/messycompound_79-19.Rdata")

u2<-data.frame(compound$wg.max.cs,compound$rf.max.cs)

u=u2


##################Marginal distributions#################

tcplot(u[,2],tlim=c(6,40),n=100)
tr1=.92
tr2=.92

th1=quantile(u[,1],tr1)
th2=quantile(u[,2],tr2)

idr<-u[which(u[,2]>th2),2]

# exD<-extremalIndex(u[,1],threshold =  th1)
# 
# d<-declust(exD)
# 
# 
# 
# exD2<-extremalIndex(u[,2],threshold =  th2)
# 
# d2<-declust(exD2)
# plot(d2)

which(is.na(u[,1]))
gpdN1<-evm(u[,1],u,family=gpd,qu=tr1,start=c(1,1))
plot(gpdN1)
gpdN2<-evm(u[,2],u,family=gpd,qu=tr2,start=c(1,1))
plot(gpdN2)
gpdN2$par


id<-seq(1,length(u[,2]))
u3<-cbind(u,id)


parv1n<-gpdN1$par
parv1n[1]<-exp(parv1n[1])
parv2n<-gpdN2$par
parv2n[1]<-exp(parv2n[1])

pxt1<-texmex::pgpd(u[which(u[,1]>th1+1),1], parv1n[1], xi=parv1n[2], u = as.numeric(th1), lower.tail = T, log.p = F)
pxt1=pxt1*(1-tr1)+tr1


pyt1<-texmex::pgpd(u[which(u[,2]>th2+1),2], parv2n[1], xi=parv2n[2], u = as.numeric(th2), lower.tail = T, log.p = F)
pyt1=pyt1*(1-tr2)+tr2


epdata <- apply(u, 2, rank, na.last = "keep")
nasm <- apply(u, 2, function(x) sum(!is.na(x)))
epdata <-  epdata/rep(nasm + 1, each = nrow(epdata))
tooLow1 <- which(u[,1] <= min(u[,1]))
tooLow2 <- which(u[,2] <= min(u[,2]))
epdata[tooLow1,1] <- 0
epdata[tooLow2,2] <- 0
px<-epdata[,1]
tr1<-approx(u[,1], px, xout = th1, method = "linear",yleft = max(px),yright = min(px), rule = 1)$y
py<-epdata[,2]
tr2<-approx(u[,2], py, xout = th2, method = "linear",yleft = max(py),yright = min(py), rule = 1)$y



interw<-1-(1-tr1)*(1+parv1n[2]*((u[which(u[,1]>th1 & u[,1]<=th1+1),1]-th1)/parv1n[1]))^(-1/parv1n[2])
interr<-1-(1-tr2)*(1+parv2n[2]*((u[which(u[,2]>th2 & u[,2]<=th2+1),2]-th2)/parv2n[1]))^(-1/parv2n[2])
ww=seq(0,1,length=length(interw))
wr=seq(0,1,length=length(interr))

empix<-px[which(u[,1]>th1 & u[,1]<=th1+1)]
empiy<-py[which(u[,2]>th2 & u[,2]<=th2+1)]

interdw<-((1-ww)*empix+ww*interw)
interdr<-((1-wr)*empiy+wr*interr)

pxx<-px[-which(u[,1]>th1)]
pyy<-py[-which(u[,2]>th2)]

pxf<-as.numeric(c(pxx,interdw,pxt1))
pyf<-as.numeric(c(pyy,interdr,pyt1))
plot(pyt1)

pxp<-c(seq(0.999,0.999999,by=0.000001))

u1p<-SpatialExtremes::qgpd(pxp, loc=th1, scale=parv1n[1], shape=parv1n[2], lower.tail = T,lambda=tr1)

pyp<-c(seq(0.999,0.999999,by=0.000001))
u2p<-SpatialExtremes::qgpd(pyp, loc=th2, scale=parv2n[1], shape=parv2n[2], lower.tail = TRUE,lambda=tr2) 

u1b<-c(u[,1],u1p)
u2b<-c(u[,2],u2p)


pxo<-c(u3[which(u3[,1]<=th1),3],u3[which(u3[,1]>th1 & u3[,1]<=th1+1),3],u3[which(u3[,1]>th1+1),3]) 
pxf<-pxf[order(pxo)] 
px<-px[order(pxo)] 
plot(pxf)

pyo<-c(u3[which(u3[,2]<=th2),3],u3[which(u3[,2]>th2 & u3[,2]<=th2+1),3],u3[which(u3[,2]>th2+1),3])
pyf<-pyf[order(pyo)]
plot(pyf)

plot(px[order(px)],pxf[order(pxf)],xlim=c(0.9,1),ylim=c(0.9,1),type="p")
segments(x0=0,y0=0,x1=1,y1=1,col=2)
plot(pxf[order(pxf)])


pxfp<-as.numeric(c(pxf,pxp))
pyfp<-as.numeric(c(pyf,pyp))

plot(pxf,u[,1])
plot(pyf,u[,2])

plot(pxfp,u1b)
plot(pyfp,u2b)

kk<-data.frame(pxf,pyf)
kk1<-data.frame(px,py)
vtau<-cor.test(x=u[,1],y=u[,2],method="kendall")$estimate


u2=u
pbas=  1/111
4555/41
111*0.01

pobj=0.001
beta=100

vtau=vtau
mar1=u1b
mar2=u2b
px=pxfp
py=pyfp

e1<-seq(0,1.2*max(u2[,1]),length.out=200)
e2<-seq(0,1.2*max(u2[,2]),length.out=200)
evp<-as.data.frame(cbind(e1,e2))
bwaa=expand.grid(e1,e2)
aa<-kcde(u2, gridsize=300,tail.flag = "upper.tail",xmin=c(0,0),xmax=c(1.5*max(u2[,1]),1.5*max(u2[,2])))

lox<-aa$eval.points[[1]]
loy<-aa$eval.points[[2]]
ngx=100000

godx<-approx(mar1,px, n = ngx, method = "linear",
             yleft = min(px), yright = max(px))
gody<-approx(mar2,py, n = ngx, method = "linear",
             yleft = min(py), yright = max(py))


pxe<-approx(godx$x, godx$y, xout = lox, method = "linear",yleft = min(px),yright = max(px), rule = 1)$y
pye<-approx(gody$x, gody$y, xout = loy, method = "linear",yleft = min(py),yright = max(py), rule = 1)$y

if(devplot==T){
  plot(aa,cont = c(0.05,0.1),display="filled.contour" ,col=viridis(10))
  
}
az<-aa$estimate
wq<-contourLines(lox,loy, az, levels = c(0.7,0.5,0.2,0.1,pbas))
wq0ri<-cbind(wq[[5]]$x,wq[[5]]$y)

var1<-u2[,1]
var2<-u2[,2]
wu<-data.frame(u2)

ovar1<-var1[order(var1)]
rx1<-rank(ovar1)
ovar2<-var2[order(var2)]
rx2<-rank(ovar2)
p01=0.01
p02=0.01
q0=0.92





if(length(which(wq0ri[,1]<0))>0)wq0ri<-wq0ri[-which(wq0ri[,1]<0),]
if(length(which(wq0ri[,2]<0))>0)wq0ri<-wq0ri[-which(wq0ri[,2]<0),]


wqUnitx<- approx(mar1, px, xout = wq0ri[,1], method = "linear",yleft = min(px),yright = max(px), rule = 1)$y
wqUnity<- approx(mar2, py, xout = wq0ri[,2], method = "linear",yleft = min(py),yright = max(py), rule = 1)$y


#Transforming Pbase to Frechet Margins
wqxFrechet<--1/log(wqUnitx)
wqyFrechet<--1/log(wqUnity)


pb<-1-pbas

#Setting up Pobjective


s=pbas/pobj


xFrechet<--1/log(px)
yFrechet<--1/log(py)

#############################Use JT_KDE model to assess joint return period of every event###########
ModelHugo_file="C:/Users/k1638615/King's College London/OneDrive - King's College London/DATA/CODE ANALYSIS/Model comparison/LF/BveLTDep.R"
ModelHugo_file="C:/Users/Alois/OneDrive - King's College London/DATA/CODE ANALYSIS/Model comparison/LF/BveLTDep.R"
source (ModelHugo_file)
if(vtau<0)q0=0.9
qq<-c()
cm<-c()
cq<-c()
ccd<-c()
qc<-.95
rq0<-seq(0.75,0.95,by=0.01)
for(q0 in rq0){
  estims<-try(BveLTDep (data= kk,mod.thresh.u = q0,crit.lev.u = qc,sig.lev=0.05,ci.meth='se',marg.inf=T),silent = T)
  qd<-try((estims$par[2]),silent=T)
  cd<-try((estims$chiCIs),silent=T)
  cc<-try((estims$chi),silent=T)
  
  if(is.numeric(qd)){
    qq<-c(qq,qd)
    ccd<-c(ccd,cc)
  }
  if(is.numeric(cd)){
    cq<-c(cq,cd)
  }
}
plot(rq0,ccd)
plot(rq0,qq)
sumd<-0
for(i in 1:length(diff(qq))){
  sdiff<-diff(qq)[i]
  sumd<-c(sumd,sumd[i]+sdiff)}
sh<-which(sumd<=-0.05|sumd>=0.05) [1]
sh=sh[length(sh)]
q0<-rq0[sh-1]
plot(sumd)
q0= estims<-try(BveLTDep (data= kk,mod.thresh.u = q0,crit.lev.u = qc,sig.lev=0.05,ci.meth='se',marg.inf=T),silent = T)
chat=NA
etahat=NA
Chilow=NA
Chimed=NA
try(chat<-estims$par[1],silent=T)
try(etahat<-estims$par[2],silent=T)
try(Chilow<-estims$chiCIs[1],silent=T)
try(Chimed<-estims$chi,silent=T)

#Loop for asymptotic dependence

if (!is.na(Chilow) & (Chilow<0.05 & etahat<0.75| etahat<0.6)){
  print("AI")
  bet=beta
  m1= 1- (wqxFrechet/(wqxFrechet+wqyFrechet))^bet
  m2<-1- (wqyFrechet/(wqxFrechet+wqyFrechet))^bet
  eta1<-m1*etahat + (1-m1)
  eta2<-m2*etahat + (1-m2)
  projx<-s^(eta1)*wqxFrechet
  projy<-s^(eta2)*wqyFrechet 
}else{ 
  projx<-s*wqxFrechet
  projy<-s*wqyFrechet}



projbackx<-exp(-(1/projx))
projbacky<-exp(-(1/projy))




objx<-approx(px,mar1, xout = projbackx, method = "linear",yleft = min(mar1),yright = max(mar1), rule = 1)$y
objy<-approx(py,mar2, xout = projbacky, method = "linear",yleft = min(mar2),yright = max(mar2), rule = 1)$y 

wqobj<-cbind(objx, objy)


  pobj<-c(seq(0.00001,0.0002,by=0.000001),seq(0.0002,pbas,by=0.0001))
  
  s=pbas/pobj
  
  wqobjf<-c()
  for (sl in 1:length(s)){
    
    if (!is.na(Chilow) & (Chilow<0.05 & etahat<0.75| etahat<0.6)){
      # print("AI")
      bet=beta
      m1= 1- (wqxFrechet/(wqxFrechet+wqyFrechet))^bet
      m2<-1- (wqyFrechet/(wqxFrechet+wqyFrechet))^bet
      eta1<-m1*etahat + (1-m1)
      eta2<-m2*etahat + (1-m2)
      projx<-s[sl]^(eta1)*wqxFrechet
      projy<-s[sl]^(eta2)*wqyFrechet 
    }else{ 
      projx<-s[sl]*wqxFrechet
      projy<-s[sl]*wqyFrechet}
    
    projbackx<-exp(-(1/projx))
    projbacky<-exp(-(1/projy))
    
    
    
    objx<-approx(px,mar1, xout = projbackx, method = "linear",yleft = min(mar1),yright = max(mar1), rule = 1)$y
    objy<-approx(py,mar2, xout = projbacky, method = "linear",yleft = min(mar2),yright = max(mar2), rule = 1)$y 
    
    wqobj<-cbind(objx, objy)
    wqobj<-removeNA(wqobj)
    wqobj<-data.frame(wqobj)
    wqobj[,1]<-jitter(wqobj[,1])
    xlt=seq(min(wqobj[,1]),max(wqobj[,1])-0.1,length.out = 120)
    xlto=seq(min(wq0ri[,1]),max(wq0ri[,1])-0.1,length.out = 120)
    
    
    if(length(wqobj[,1])<102){
      repeat{
        mirror<-wqobj[c(length(wqobj[,1]):1),]
        wqobj<-rbind(wqobj,mirror)
        if(length(wqobj[,1])>=102) break
      }
      wqobj<-round(wqobj,8)
      wqobj<-wqobj[order(wqobj[,1],-wqobj[,2]),]
      wqobj[,1]<-jitter(wqobj[,1])
    }
    
    ltl<-digit.curves.p(start=wqobj[1,], curve=as.matrix(wqobj), nPoints=98, closed = FALSE)
    ltl1<-ltl
    ltl2<-ltl
    if (pobj[sl]>0.0000001){
      gridx<-(seq(min(ltl[,1],na.rm=T),max(mar1),length.out=100))
      gridy<-(seq(min(ltl[,2],na.rm=T),max(mar2),length.out=100))
      ltl1[,1]<-approx(ltl[,1], ltl[,2], xout = gridx, method = "linear", rule = 1)$x
      ltl1[,2]<-approx(ltl[,1], ltl[,2], xout = gridx, method = "linear", rule = 1)$y
      
      ltl2[,1]<-approx(ltl[,2], ltl[,1], xout = gridy, method = "linear", rule = 1)$y
      ltl2[,2]<-approx(ltl[,2], ltl[,1], xout = gridy, method = "linear", rule = 1)$x
      
    }
    ltl1<-data.frame(ltl1,rep(pobj[sl],100),rep(sl,100))
    ltl2<-data.frame(ltl2,rep(pobj[sl],100),rep(sl,100))
    wqobjf<-rbind(wqobjf,ltl1,ltl2)
    
  }
  
  
  plot(wqobjf[,1],wqobjf[,2],col=wqobjf[,4])
  points(u2, col=alpha("black",1),pch=16)
  
  tg=100
  

  gridx<-seq(min(wqobjf[,1],na.rm=T),max(mar1),length.out=tg)
  gridy<-seq(min(wqobjf[,2],na.rm=T),max(mar2),length.out=tg)
  
  
  pxg<-approx(mar1, px, xout = gridx, method = "linear",yleft = min(px),yright = max(px), rule = 1)$y
  pyg<-approx(mar2, py, xout = gridy, method = "linear",yleft = min(py),yright = max(py), rule = 1)$y
  
  
  matjt<-matrix(ncol=tg,nrow=tg)
  for (k in 1:(length(pxg)-1)){
    colx<-which(wqobjf[,1]>gridx[k] & wqobjf[,1]<=gridx[k+1])
    for (j in 1:(length(pyg)-1)){
      coly<-wqobjf[colx,3][which(wqobjf[colx,2]>gridy[j] & wqobjf[colx,2]<=gridy[j+1])]
      if(length(coly)==0){matjt[j,k]=NA}else{
        matjt[k,j]=mean(coly,na.rm=T)}
    }
  }
  
  for (k in 1:(length(pxg)-1)){
    colx<-which(wqobjf[,2]>gridy[k] & wqobjf[,2]<=gridy[k+1])
    for (j in 1:(length(pyg)-1)){
      coly<-wqobjf[colx,3][which(wqobjf[colx,1]>gridx[j] & wqobjf[colx,1]<=gridx[j+1])]
      if(length(coly)==0){matjt[j,k]=matjt[j,k]}else{
        matjt[j,k]=mean(coly,na.rm=T)}
    }
  }
  
  grid <- expand.grid(lon=gridx, lat=gridy)
  
  
  aac<-kcde(u2[,c(1,2)], gridsize=100,tail.flag = "upper.tail",xmin=c(min(gridx),min(gridy)),xmax=c(max(gridx),max(gridy)))
  plot(aac,"filled.contour")
  matjt
  
  
  library(zoo)
  
  for (nap in 1: length(pxg)){ matjt[,nap]<-na.approx(matjt[,nap],maxgap = 5,na.rm=F)}
  
  for (nap in 1: length(pxg)){ matjt[nap,]<-na.approx(matjt[nap,],maxgap = 5,na.rm=F)}
  
  levelplot(matjt ~ lon * lat, data=grid, cuts=20, pretty=T,contour=T) 
  
  
  contour(gridx,gridy,matjt,levels=seq(0.01,0.0001,by=-.0001))
  obj<- list( x= gridx, y=gridy, z= matjt)
  levelplot(obj$z ~ lon * lat, data=grid, cuts=20, pretty=T,contour=T) 
  points(jtres$levelcurve,pch=16)
  length(unique(u2$Var2))
  points(u2,pch=16,col="red")

  krap<-as.vector(matjt) 
  kdrap<-as.vector(aac$estimate)
  krap[which(is.na(krap))]=kdrap[which(is.na(krap))]
  locas<-expand.grid(gridx,gridy)
  locas$idl<-c(1:length(locas$Var1))
  names(u2)<-c("Var1","Var2")

  lokrap<-cbind(locas,krap)
  obl<- array(krap, dim=c(100,100))
  objl<- list( x= gridx, y=gridy, z= obl)
  levelplot(objl$z ~ lon * lat, data=grid, cuts=20, pretty=T,contour=T) 
  
  #important function computing the 2d density field
  omerd<-fields::interp.surface(obj, u2)
  plot(omerd)
  
  
  abline(h=1e-3)
  4555/41
  111*0.01
  1/111
  u2<-data.frame(u2[,c(1,2)],omerd)
  plot(u2[,c(1,2)],pch=16,col=factor(u2$omerd))
  rbPal2<-colorRampPalette(c("purple","red","darkorange","orange",'skyblue'))
  inil<-data.frame(wq0ri)

  ux<-u2[which(!is.na(u2$omerd)),]
  cfx<-compound[which(!is.na(u2$omerd)),]
  library(lubridate)
  plot(year(cfx$startime),month(cfx$startime))
  cfx$year=year(cfx$startime)
  pyear<-aggregate(cfx$combin ,
                  by = list(year= cfx$year),
                  FUN = function(x) c(len = length(na.omit(x))))
  pyear <- do.call(data.frame, pyear)
  plot(pyear,type="h")
  mean(pyear$x)
  
  #plot of compound hazard clusters with 1 year level curve (See appendix H)
  ggplot(data=u2,aes(x=Var1,y=Var2,fill=omerd))+
    geom_point(size=4,shape=21,alpha=0.7)+
    scale_fill_gradientn(na.value = "gray95" ,trans="log",colors=rbPal2(100),"Joint Return Period (Y)",breaks=c(1/111,.1/111,.01/111),labels = c(1,10,100))+
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=20,face="italic"),
          panel.background = element_rect(fill = "white", colour = "black",size=1),
          legend.title = element_text(size=20),
          legend.text = element_text(size=18),
          axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.key.size = unit(1, "cm"))+
  geom_line(data=inil, aes(x=X1,y=X2,fill=.01),size=2,color="lightblue4",alpha=.7)+
    scale_x_continuous(limits=c(17,50),expression(paste("w" ['a']," [m s "^-1,"]")))+
  scale_y_continuous(expression(paste("r" [a]," [mm]")),breaks=c(0,25,50,75,100,125,150))
  
  
  
  #save joint return period
  compound$rp<-u2$omerd
  compound$footprint=compound$spacescale/1485*100
  plot(1/(compound$rp*111),compound$ORtime,log="x")
  plot(compound$footprint[which(!is.na(compound$rp))],compound$ORtime[which(!is.na(compound$rp))],log="xy")
  compoundX<-compound[which(!is.na(compound$rp)),]
  compoundY<-compound[which(!is.na(compound$rp)),]
  mean(compoundX$footprint)/  mean(compound$footprint)
  mean(compoundX$ORtime)
  save(compound,file="out/CompoundRW_79-19.RP.Rdata")
  write.csv(compound,"out/CompoundRW_79-19.RP.csv")
  save(compoundX,file="out/CompoundRW_79-19.ME2.Rdata")
  tr<-aggregate(list(len=compoundX$combin),
                by = list(loc=compoundX$maxlocsR),
                FUN = function(x) c(length=length(unique(x))))
  tr <- do.call(data.frame, tr)
  tr$loc=as.character(tr$loc)
  Rcart=strsplit(tr$loc," ")
  Rcar<-as.data.frame(matrix(unlist(Rcart),ncol=2,byrow=T))
  Rcar[,1]<-as.numeric.factor((Rcar[,1]))
  Rcar[,2]<-as.numeric.factor((Rcar[,2]))
  tr<-data.frame(tr,Rcar)

  
  tw<-aggregate(list(len=compoundX$combin),
                by = list(loc=compoundX$maxlocsW),
                FUN = function(x) c(length=length(unique(x))))
  tw <- do.call(data.frame, tw)
  tw$loc=as.character(tw$loc)
  Rcort=strsplit(tw$loc," ")
  Rcor<-as.data.frame(matrix(unlist(Rcort),ncol=2,byrow=T))
  Rcor[,1]<-as.numeric.factor((Rcor[,1]))
  Rcor[,2]<-as.numeric.factor((Rcor[,2]))
  tw<-data.frame(tw,Rcor)
  
  
  #############################clustering of areas impacted by the same events######
  minifun<-function(sequence){
  
  out<-as.numeric(unlist((strsplit(sequence, " "))))
  return(out)
  }
  minifun(secul$sequence[1])
  seq <- lapply(secul$sequence, minifun) #I am assuming that get1grams produces a vector
  names(seq) <- paste0("seq_", secul$local)

  
  library(arules)
  seqTrans <- as(seq, "transactions") #create a transactions object
  seqMat <- as(seqTrans, "matrix") #turn the transactions object into an incidence matrix each row represents a sequence and each column a 1gram each cell presence/absence of the 1gram
  seqMat <- +(seqMat) #convert boolean to 0/1
  j.dist <- dist(seqMat, method = "binary") #make use of base R's distance function
  
  ##Matrix multiplication to calculate the jaccard distance
  tseqMat <- t(seqMat)
  a <- t(tseqMat) %*% tseqMat
  b <- t(matrix(rep(1, length(tseqMat)), nrow = nrow(tseqMat), ncol = ncol(tseqMat))) %*% tseqMat
  b <- b - a
  c <- t(b)
  j <- as.dist(1-a/(a+b+c))
  ja<-data.frame(as.matrix(j))
  
  # methods to assess
  m <- c( "average", "single", "complete", "ward")
  names(m) <- c( "average", "single", "complete", "ward")
  
  # function to compute coefficient
  ac <- function(x) {
    agnes(ja, method = x)$ac
  }
  
  map_dbl(m, ac)
  ##   average    single  complete      ward 
  ## 0.7379371 0.6276128 0.8531583 0.9346210
  
  clusters <- hclust(j,  method = "ward")
  km.res <- eclust(ja, "kmeans", k = 3, nstart = 25, graph = FALSE)
  km.res$silinfo
  tespx$clust<-km.res$cluster
    
  fviz_cluster(km.res, geom = "point", ellipse.type = "norm",
               palette = "jco", ggtheme = theme_minimal())
  
  clusterCut <- cutree(clusters, h=20)
  table(clusterCut)
  plot(clusters, cex = 0.6)
  rect.hclust(clusters, k = 8, border = 2:5)
  sub_grp <- cutree(clusters, k = 8)
  tespx$clust<-sub_grp
  
  library (factoextra)
  fviz_cluster(list(data = j, cluster = sub_grp))
  fviz_nbclust(ja, FUN = hcut, method = "silhouette")
  fviz_nbclust(ja, FUN = hcut, method = "wss",k.max= 20)
  gap_stat <- clusGap(ja, FUN = hcut, nstart = 5, K.max = 10, B = 2)
  fviz_gap_stat(gap_stat)
  
  length(!is.na(ja))/length(ja)
  library(NbClust)
  dunix<-c()
  for (csize in 2:15){
  print(csize)
  km.res <- eclust(ja, "kmeans", k = csize, nstart = 25, graph = FALSE)
  dunix<-c(dunix,dunn(ja, as.numeric(km.res$cluster)))
  }
  
  dunix
  tg<-cluster.stats(d = ja, clustering=sub_grp)
  tg
  tespx$loc<-as.character(tespx$loc)
  locus<-strsplit(tespx$loc," ")
  locus<-as.data.frame(matrix(unlist(locus),ncol=2,byrow=T))
  tespx$lon<-as.numeric((locus[,1]))
  tespx$lat<-as.numeric((locus[,2]))
  tespx$grp=1
  ttstep<-41
  tespx$tcomb<-tespx$len/ttstep
  tespx$xratio<-tespx$tcomb/locav$tcom
  tespx$mday=dayof
  tespx$rxx<-rxt
  
  rgb.palette=colorRampPalette(c("lightblue","royalblue","darkolivegreen3","gold","orange","red","darkred")
                               ,interpolate="linear",bias=1)
  
  rgb.palette2=colorRampPalette(c("darkred","orange","gold","darkolivegreen3","lightblue", 
                                  "royalblue","darkred"),interpolate="linear",bias=1)
  rbPal <- colorRampPalette(c('darkgreen',"gold","darkorange",'red'))
  
  mapUK = SpacializedMap(regions = c("UK","France","Spain","Portugal","Italy","Ireland"))
  
  library(ggalt)
  library(metR)
  
  ukk<- gSimplify(mapUK, tol=0.001, topologyPreserve=TRUE)
  plot(ukk)
  uk_fort <- ggplot2::fortify(ukk)
  longlims=c(-6, 2)
  longdom=seq(longlims[1],longlims[2],by=.25)
  latlims=c(48,59)
  latdom=seq(latlims[1],latlims[2],by=.25)
  
  tespx$mbins<-cut(tespx$mday, 12)
  levels(tespx$mbins)
  levels(tespx$mbins) <- list(A="alpha", B="beta", C="gamma")
  
  ggplot(uk_fort, aes(x=long,y=lat,group=group)) +
    geom_polygon(fill = "white", color = "gray10", size = 1) +
    theme_bw(16)+
    # scale_fill_gradientn(colours = rgb.palette(100),"Total hours in a CE",breaks=breaks_extended(5),limits=c(min(locav$x)-0.1*min(locav$x),max(locav$x)+.1*min(locav$x)))+
    coord_fixed(xlim = longlims,  ylim = latlims, ratio = 1.3)+
    geom_raster(data=tespx,aes(x=lon,y=lat,fill=clust,group=loc),interpolate = F,alpha=.7)+
    # geom_contour_fill(data=tespx,aes(x=lon,y=lat,z = clust,group=grp),alpha=0.6) +
    scale_fill_gradientn(trans=scales::modulus_trans(1),colours = rgb.palette(10)," ",breaks=breaks_extended(5),guide= "colorbar")+
    # scale_fill_manual(values=rgb.palette2(12),breaks = 12)+
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
  #   geom_point(data=tr,aes(x=V1,y=V2,size=len,group=loc),alpha=.5)+
  # geom_point(data=tw,aes(x=V1,y=V2,size=len,group=loc),alpha=.5,shape=18)
  # 
  guides(fill = guide_legend())

  