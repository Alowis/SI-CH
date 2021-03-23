
######################################################################################################
######################################################################################################
############################APPLICATION OF MULTIPLE MODELS TO REAL DATA###############################
######################################################################################################
######################################################################################################

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
hazard="TempFire" #"WindRain" or "TempFire"
hazard="WindRain"

load(file="out/CompoundRW_79-19.v2x.Rdata")
load(file="out/messycompound_79-19.Rdata")

u2<-data.frame(compound$wg.max.cs,compound$rf.max.cs)
##############################Rain + wind##############################
if(hazard=="WindRain"){
filei<-paste0(getwd(),"/data/Winf_Rain_anaH.Rdata")
filech<-paste0(getwd(),"/data/windRainprobHaethrow.Rdata")



load(filei)
load(filech)
u2$Rainmm[which(u2$Rainmm<=0)]=NA

u2<-u2[-which(is.na(u2$Speedms)|is.na(u2$Rainmm)),]
plot(u2)
}
###########################################Temperature and fire####################
if(hazard=="TempFire"){


filef<-paste0(getwd(),"/data/fire01meantemp.Rdata")
load(filef)
u2<-na.omit(fire01meantemp)


dir.study = getwd() 
if (!exists("dir.day"))     dir.day    = paste0(dir.study,"/out/",Sys.Date())
dir.create(dir.day,showWarnings = F,recursive = T)
dir.coef    = paste0(dir.study,"/out/",Sys.Date(),"/bootstraplife")
dir.create(dir.coef,showWarnings = F,recursive = T)

}

#################### Semi-Automatic analogous selection##############
chiest<-chi(u2,qlim=c(0.75,0.98),nq=100)
summary(chiest)
plot(chiest)
plot(chiest$quantile,chiest$chi[,2])
abline(h=mean(chiest$chi[,2]))
plot(chiest$chi[,1])
etaest=(chiest$chibar+1)/2

chint<-round(c(min(chiest$chi[,1]),mean(chiest$chi[,2]),max(chiest$chi[,3])),2)
etint<-round(c(min(etaest[,1]),mean(etaest[,2]),max(etaest[,3])),2)

etalist<-c(0.25,0.5,0.75,0.9)
chilist<-c(0.05,0.1,0.3,0.5,0.9)
et<-c()
ct=c()
et<-c(etalist[which(etalist>=etint[1]&etalist<=etint[3])])
ct<-c(chilist[which(chilist>=chint[1]&chilist<=chint[3])])
if(chint[1]>0.1)et<-c()
cat(paste(c("analogous datasets to be tested: \n eta = ",et," \n chi = ",ct),collapse= " "))

#################Nonparametric bootsrap (because nonparametric method)###################
R=100
library(boot)
inut<-function(ts){
  return(ts)
}
e2<-tsboot(as.matrix(u2), statistic=inut, R=R, l = 1, sim = "fixed")


oup<-t(e2$t)
ux<-c()
uy<-c()
for (i in (1:R)){
  ux<-cbind(ux,oup[1:length(u2[,1]),i])
  uy<-cbind(uy,oup[(length(u2[,1])+1):(2*length(u2[,1])),i])
}
plot(u2,col="grey")
ugig<-list()

for (i in (1:R)){
  ugig<-c(ugig,list(cbind(ux[,i],uy[,i])))}

######################Declare variables########################
totalchiT<-c()
totaltauT<-c()
totaletaT<-c()

totalchi<-as.data.frame(matrix(0, ncol = 6, nrow = R))
totaleta<-as.data.frame(matrix(0, ncol = 6, nrow = R))

HTlinex <- as.data.frame(matrix(0, ncol = R, nrow = 100))
LTlinex <- as.data.frame(matrix(0, ncol = R, nrow = 100))
Gumlinex <- as.data.frame(matrix(0, ncol = R, nrow = 100))
Normlinex <- as.data.frame(matrix(0, ncol = R, nrow = 100))
FGMlinex <- as.data.frame(matrix(0, ncol = R, nrow = 100))
Gallinex<- as.data.frame(matrix(0, ncol = R, nrow = 100))

HTliney <- as.data.frame(matrix(0, ncol = R, nrow = 100))
LTliney <- as.data.frame(matrix(0, ncol = R, nrow = 100))
Gumliney <- as.data.frame(matrix(0, ncol = R, nrow = 100))
Normliney <- as.data.frame(matrix(0, ncol = R, nrow = 100))
FGMliney <- as.data.frame(matrix(0, ncol = R, nrow = 100))
Galliney<- as.data.frame(matrix(0, ncol = R, nrow = 100))

HTlined <- as.data.frame(matrix(0, ncol = R, nrow = 100))
LTlined <- as.data.frame(matrix(0, ncol = R, nrow = 100))
Gumlined <- as.data.frame(matrix(0, ncol = R, nrow = 100))
Normlined <- as.data.frame(matrix(0, ncol = R, nrow = 100))
FGMlined <- as.data.frame(matrix(0, ncol = R, nrow = 100))
Gallined<- as.data.frame(matrix(0, ncol = R, nrow = 100))

#####################################################################################
################################Bootstrap loop#######################################
#####################################################################################

for (jo in 1:R){
  print(paste0("loop # ", jo))
  if(jo==0){
    u2[,2]<-as.numeric(u2[,2])
    u2<-data.frame(u2)
    u<-u2

  }else{
u<-data.frame(ugig[[jo]])
}
##########################################################
##################Marginal distributions#################
##########################################################
tr1=.95
tr2=.95

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

gpdN2<-evm(u[,2],u,family=gpd,qu=tr2,start=c(1,1))

# 
# gpdbg1<-fevd(u[,1],threshold =as.numeric(th1), location.fun = ~1,
#              scale.fun = ~1, shape.fun = ~1, use.phi = FALSE,
#              type = c("GP"),
#              method = c("MLE"), time.units = "days",na.action = na.omit)
# 

id<-seq(1,length(u[,2]))
u3<-cbind(u,id)

# gpdbg2<-fevd(u[,2],threshold =as.numeric(th2), location.fun = ~1,
#              scale.fun = ~1, shape.fun = ~1, use.phi = FALSE,
#              type = c("GP"),
#              method = c("MLE"), time.units = "days",na.action = na.omit)
# plot(gpdbg1)
# 
# parv1<-gpdbg1$results$par
# parv2<-gpdbg2$results$par

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

##########################################################
####################Joint Tail Model######################
##########################################################
if(length(u[,1])<5000){pb=0.02;t.sim=0.98}
if(length(u[,1])>=5000){pb=0.01;t.sim=0.99}

print("Joint tail model fitting")
jtres<-JT_KDEap(u2=u,pbas=0.01,pobj=upobj,beta=100,vtau=vtau,devplot=F,mar1=u1b,mar2=u2b,px=pxfp,py=pyfp,interh=interh)
  plot(jtres$levelcurve) 
  wq0ri=jtres$wq0ri
  wqobj<-removeNA(jtres$levelcurve)
  wqobj<-data.frame(wqobj)
  wq0ri<-jtres$wq0ri
  wqobj[,1]=jitter(wqobj[,1])
  plot(wqobj)
  jtres$etaJT
  jtres$chiJT
  
  wqobj<-removeNA(wqobj)
  wqobj<-data.frame(wqobj)
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
  plot(u)
  lines(ltl,col=2)

##########################################################
####################Conditional model#####################
##########################################################

  print("Conditional model fitting")
 
  
  condexres<-Cond_modap(u2=u,tr1,tr2,tsim=t.sim,num.sim=10000,pobj=0.001,interh=interh)

  jline<-data.frame(condexres$jline)
  xht=seq(ceiling(min(jline[,1])),round(max(jline[,1]),2)-0.1,length.out = 120)
  
  library(ks)
  cem.dens<-kde(condexres$onlysim ,gridsize = 200)
  jt.dens<-kde(u,gridsize = 200)
  
  #new method for curve drawing

  
  if(length(wqobj[,1])<102){
    idz<-which(diff(wqobj[,1])==0)
    zr<-as.matrix(wqobj[idz,])
    if(length(zr)>0)wqobj<-wqobj[-idz,]
    ltlq=approx(wqobj[,1], wqobj[,2], xout = xlt, method = "linear",yleft = max(wqobj[,2]),yright = min(wqobj[,2]), rule = 1)
    ltlq<-cbind(ltlq$x,ltlq$y)
    ltlq<-data.frame(rbind(ltlq,zr))
    wqobj<-ltlq[order(ltlq[,1]),]
    
  }
  if(length(wq0ri[,1])<102){
    idz<-which(diff(wq0ri[,1])==0)
    zr<-wq0ri[idz,]
    if(length(zr)>0)wq0ri<-wq0ri[-idz,]
    ltlqo=approx(wq0ri[,1], wq0ri[,2], xout = xlto, method = "linear",yleft = max(wq0ri[,2]),yright = min(wq0ri[,2]), rule = 1)
    ltlqo<-cbind(ltlqo$x,ltlqo$y)
    ltlqo<-data.frame(rbind(ltlqo,zr))
    wq0ri<-ltlqo[order(ltlqo[,1]),]
    
  }
  ltl<-digit.curves.p(start=wqobj[1,], curve=as.matrix(wqobj), nPoints=98, closed = FALSE)
  ltlo<-digit.curves.p(start=wq0ri[1,], as.matrix(wq0ri), nPoints=98, closed = FALSE)
  htl<-digit.curves.p(start=jline[1,], as.matrix(jline), nPoints=98, closed = FALSE)
  
  densicurvem<-function(kdetab,lines,tl,lines2){
    densim<-c()
    aba<-kdetab$estimate/sum(kdetab$estimate)
    for (stro in 1:100){
      pwin<-as.matrix(lines[stro,])
      onarray<-aba[which(abs(kdetab$eval.points[[1]]-pwin[1])==min(abs(kdetab$eval.points[[1]]-pwin[1]))),which(abs(kdetab$eval.points[[2]]-pwin[2])==min(abs(kdetab$eval.points[[2]]-pwin[2])))]
      if(length(onarray)==0)onarray<-aba[which(round(kdetab$eval.points[[1]],1)==round(pwin[1],1)),which(round(kdetab$eval.points[[2]],1)==round(pwin[2],1))]
      if(length(onarray)==0)print(stro)
      if(length(onarray)>1)onarray=mean(onarray)
      densim<-c(densim,onarray)
      
    }
    dm<-sum(densim)
    densim<-densim/dm
    posoncurve<-seq(1/length(densim),1,length=length(densim))
    plot(posoncurve,densim, type="h")
    if (tl=="l"){
      exit<-data.frame(lines2,densim)}
    if (tl=="h"){
      exit<-data.frame(lines,densim)}
    names(exit)=c("x","y","dens")
    exit
  }
  densicurvcop<-function(lines, copi){
    linep<-lines
    linep[,1]<-spline(u[,1],pxf, n = 300, method = "fmm",
                 xmin = min(lines[,1]), xmax = max(lines[,1]), ties = mean,xout=lines[,1])$y
    linep[,2]<-spline(u[,2],pyf, n = 300, method = "fmm",
                  xmin = min(lines[,2]), xmax = max(lines[,2]), ties = mean,xout=lines[,2])$y
    shit<-copula::dCopula(u=as.matrix(data.frame(linep)), copula =copi)
    shma<-sum(shit)
    shit<-shit/shma
    posonshit<-seq(1/length(shit),1,length=length(shit))
    plot(posonshit,shit, type="h")
    shuba<-data.frame(lines,shit)
    names(shuba)=c("x","y","dens")
    shuba
  }
  
  htl<-densicurvem(cem.dens,htl, tl="h")
  ltl<-densicurvem(jt.dens,ltlo, tl="l", ltl)

  denslim<-1e-9
  
  
  
##########################################################
####################Copula models#########################
##########################################################
  
  print("Copula fitting")
  
  library(VineCopula)
  library(Metrics)
  
  copulas<-list()
  copulas[[1]]<-fitCopula(gumbelCopula(),data=kk,method="itau")
  print("one passed")
  copulas[[2]]<-fitCopula(normalCopula(),data=kk,method="itau")
  print("one passed")
  copulas[[3]]<-fitCopula(fgmCopula(),data=kk, method="itau")
  print("one passed")
  copulas[[4]]<-NA
  galsh=F
  try(copulas[[4]]<-fitCopula(galambosCopula(),data=kk,method="itau"),silent=T)
  if(is.na(copulas[[4]])){copulas[[4]]<-indepCopula()
  galsh=T}
  print("one passed")
  
  library(acopula)
  library(copBasic)
  NORMcop <- function(u,v, para=NULL) {
    if(length(u) == 1) u <- rep(u, length(v)) # see asCOP() for reasoning of
    if(length(v) == 1) v <- rep(v, length(u)) # this "vectorization" hack
    return(copula::pCopula(matrix(c(u,v), ncol=2),normalCopula(para)))
  }
  copbask<-c(GHcop,NORMcop,FGMcop,GLcop)
  clines<-list()
  for (v in 1:4){

    if( v!=4 | galsh==F){
      copulas[[v]]@copula
      c1<-copulas[[v]]@estimate
    }
    if( v==1){o<-gumbelCopula(c1,dim=2)}
    if( v==2){o1<-copNormal(parameters=c1,dim=2)
    o<-normalCopula(c1,dim=2)}
    if( v==3)o<-fgmCopula(c1,dim=2)
    if( v==4 & galsh==F)o<-galambosCopula(c1)
    if( v==4 & galsh==T)o<-indepCopula(dim=2)
    
    # if( v==4 & vtau<=0)o<-indepCopula(dim=2)
    coco<-copbask[[v]]
    
    
    cl1<-curve.funct.b(pxfp,pyfp,u1b,u2b,u,pos="l",pobje=upobj,ng=100,inter=interh) 
    cl2<-curve.funct.b(pxfp,pyfp,u1b,u2b,u,pos="m",pobje=upobj,ng=100,inter=interh) 
    cl3<-curve.funct.b(pxfp,pyfp,u1b,u2b,u,pos="r",pobje=upobj,ng=100,inter=interh) 
    # contour(godx,gody,acp2)
    
    plot(u,xlim=c(0,60),ylim=c(0,300))
    points(cl1,col=2,type="o")
    points(cl2,col=3,type="o")
    points(cl3,col=4,type="o")
    

    
    cl<-rbind(cl1,cl2,cl3)
    cl<-round(cl,8)
    cl<-cl[order(cl[,1],-cl[,2]),]
    cl[,1]<-jitter(cl[,1])
    if(length(cl[,1])<102){
      repeat{
        mirror<-cl[c(length(cl[,1]):1),]
        cl<-rbind(cl,mirror)
        if(length(cl[,1])>=102) break
      }
      cl<-round(cl,8)
      cl<-cl[order(cl[,1],-cl[,2]),]
      cl[,1]<-jitter(cl[,1])
    }
    
    cli<-digit.curves.p(start=c(cl[1,1],cl[1,2]), as.matrix(cl), nPoints=98, closed = FALSE)
    cli<-densicurvcop(cli,o)
    clines<-c(clines,list(cli))
    
  }

  
  names(clines)=c("gumbel","normal","fgm","galambos")

  
  

######################################################################
####################Compute and store results and diagnostics#########
######################################################################
  if(jo==0){
  
  tauN<-tau(copulas[[2]]@copula)
  tauGu<-tau(copulas[[1]]@copula)
  taufgm<-tau(copulas[[3]]@copula)
  if(galsh==F)tauGal<-tau(copulas[[4]]@copula)
  if(galsh==T)tauGal=NA
  
  totaltau<-c(vtau,tauGu,tauN,taufgm,tauGal)
  names(totaltau)=c("empirical","gumbel","normal","fgm","galambos")
  
  #chi for gumbel 2-2^alpha(hugo)
  #fgm copula chi = 0
  #for galambos = 2-2*A(1/2)
  #For normal copula BicpoPaer2Tail =0
  #trouver des moyens d'estimer chi
  # 
  etaht=condexres$etaht
  etalt<-jtres$etaJT
  etaGum<-1
  etafgm=0.5
  etaNorm<-(1+copulas[[2]]@estimate)/2
  
  # An<-abvevd(dep=param[4],model="neglog")
  etaGal=1
  # etaemp<-(1+rho)/2
  etadep<-c(etaht,etalt,etaGum,etaNorm,etafgm,etaGal)
  names(etadep)= c("HT","LT","Gumcop","Normalcop","FGMcop","GalambosCop")
  
  chiht=condexres$chiHT
  chilt<-jtres$chiJT
  chiGum<-2-2^(1/copulas[[1]]@estimate)
  chifgm=0
  chiNorm<-BiCopPar2TailDep(1, copulas[[2]]@estimate)$upper
  if(vtau>0){
    An<-abvevd(dep=copulas[[4]]@estimate,model="neglog")
    chiGal=2-2*An}
  if(vtau<=0){
    chiGal=NA
  }
  
  # chireal<-2-2^(1/theta)
  chidep<-c(chiht,chilt,chiGum,chiNorm,chifgm,chiGal)
  names(chidep)= c("HT","LT","Gumcop","Normalcop","FGMcop","GalambosCop")
  
  linesInit<-list(htl,ltl,clines[[1]],clines[[2]],clines[[3]],clines[[4]])
  names(linesInit)<-c("HT","LT","Gumcop","Normalcop","FGMcop","GalambosCop")
  coefs<-list(totaltau,chidep,etadep)
  names(coefs)<-c("tau","chi","eta")
  if (hazard=="WindRain"){
    xl=expression(paste("Daily wind gust [m ",s^-1,"]"))
    yl="Daily rainfall [mm]" 
  }
  if (hazard=="TempFire"){
  xl="Mean daily temperature [°C]"
  yl="# of wildfire/day"
  }
  mods<-c("Cond.Ex","JT.KDE","Gumcop","Normalcop","FGMcop","GalambosCop")
  cols <- c("Cond-Ex" = "purple", "JT-KDE" = "green", "Gumcop" = "red", "Normalcop" = "skyblue", "FGMcop" = "orange", "GalambosCop" = "chocolate4") 
  #
  coland<-c( "purple","green", "red","skyblue","orange","chocolate4")
  
  plot(u, col="grey",xlim=c(0,1.3*max(ltl[,1])),ylim=c(0,1.3*max(ltl[,2])),main=paste0("Joint exceedance curves for p = 0.0001 \n RP = 240 years @Heathrow"),xlab="Hazard 1",ylab="Hazard 2")
  # legend(30,35,legend=c("JTKDE", "CONdEx","GumbelCop","NormalCop","FGMCop","GalambosCop"),
         # col=c("green", "purple", "red","skyblue","orange","chocolate4"), lty=c(1,1,1,1,1,1), cex=.9, lwd=3)
  lines(ltl[,1],ltl[,2],col="green", lwd=3)
  lines(htl[,1],htl[,2],col="purple", lwd=3,type="l")
  lines(clines[[1]],col="red", lwd=3,type="l")
  lines(clines[[2]],col="skyblue", lwd=3,type="l")
  lines(clines[[3]],col="orange", lwd=3,type="l")
  lines(clines[[4]],col="chocolate4",lwd=3,type="l")
  

  print(coefs)
  print("results from the raw data /n If you want to continie strat loop from 1")
  
  break
  
  
  # save(linesInit,file=paste0(dir.coef,"/lineinit_T_F_Porto_0.001.Rdata"))
  # save(coefs,file=paste0(dir.coef,"/coefsInit_T_F_Porto_0.001.Rdata"))
  # 
  # save(linesInit,file=paste0(dir.coef,"/lineinit_W_R_Heathrow_0.001.Rdata"))
  # save(coefs,file=paste0(dir.coef,"/coefsInit_W_R_Heathrow_0.001.Rdata"))
  }else{
  
  etaht=condexres$etaht
  etalt<-jtres$etaJT
  etaGum<-1
  etafgm=0.5
  etaNorm<-(1+copulas[[2]]@estimate)/2
  

  etaGal=1
  totaleta[jo,]<-c(etaht,etalt,etaGum,etaNorm,etafgm,etaGal)
  names(totaleta)= c("HT","LT","Gumcop","Normalcop","FGMcop","GalambosCop")
  
  chiht=condexres$chiHT
  chilt<-jtres$chiJT
  chiGum<-2-2^(1/copulas[[1]]@estimate)
  chifgm=0
  chiNorm<-BiCopPar2TailDep(1, copulas[[2]]@estimate)$upper
  if(vtau>0){
    An<-abvevd(dep=copulas[[4]]@estimate,model="neglog")
    chiGal=2-2*An}
  if(vtau<=0){
    chiGal=NA
  }
  
  totalchi[jo,]<-c(chiht,chilt,chiGum,chiNorm,chifgm,chiGal)
  names(totalchi)= c("HT","LT","Gumcop","Normalcop","FGMcop","GalambosCop")
  
  
  names(u)= c("V1","V2")
  plot(u, col="grey",xlim=c(0,(max(u$V1)+2*sd(u$V1))),ylim=c(0,(max(u$V2)+2*sd(u$V2))),main=paste0("Joint exceedance curves for p = 0.001"))
  legend((max(u$V1)),(max(u$V2)),legend=c("True","L&T", "H&T","GumbelCop","NormalCop","FGMCop","GalambosCop"),
         col=c("blue", "green", "purple","red","red","red","red"), lty=c(1,1,1,3,4,5,6), cex=.9, lwd=2)
  # contour(lox,loy,az,levels = c(0.1,0.05), add=T, col="darkgreen",delta=0,xlim=c(0,100),ylim=c(0,100),box=F, lwd=4)
  lines(wqobj,col="green", lwd=2)
  lines(jline,col="purple", lwd=2,type="l")
  lines(clines[[1]],col="red",lty=3, lwd=2,type="l")
  lines(clines[[2]],col="red",lty=4, lwd=2,type="l")
  lines(clines[[3]],col="red",lty=5, lwd=2,type="l")
  lines(clines[[4]],col="red",lty=6, lwd=2,type="l")
  
  
  
  
  HTlinex[,jo]<-htl[,1]
  LTlinex[,jo]<-ltl[,1]
  Gumlinex[,jo]<-clines[[1]][,1]
  Normlinex[jo]<-clines[[2]][,1]
  FGMlinex[,jo]<-clines[[3]][,1]
  Gallinex[,jo]<-clines[[4]][,1]
  
  HTliney[,jo]<-htl[,2]
  LTliney[,jo]<-ltl[,2]
  Gumliney[,jo]<-clines[[1]][,2]
  Normliney[jo]<-clines[[2]][,2]
  FGMliney[,jo]<-clines[[3]][,2]
  Galliney[,jo]<-clines[[4]][,2]
  
  HTlined[,jo]<-htl[,3]
  LTlined[,jo]<-ltl[,3]
  Gumlined[,jo]<-clines[[1]][,3]
  Normlined[jo]<-clines[[2]][,3]
  FGMlined[,jo]<-clines[[3]][,3]
  Gallined[,jo]<-clines[[4]][,3]
}
}

##########################################################################################
##################Computing and storing CIs for level curves and dep measure estimations##
##########################################################################################

linesfinalx<-list(HTlinex,LTlinex, Gumlinex, Normlinex, FGMlinex,Gallinex)
linesfinaly<-list(HTliney,LTliney, Gumliney, Normliney, FGMliney,Galliney)
linesfinald<-list(HTlined,LTlined, Gumlined, Normlined, FGMlined,Gallined)
# load(file="C:/Users/PhD Student/OneDrive - King's College London/DATA/CODE ANALYSIS/Model comparison/out/2019-11-12/bootstraplife/lines_T_F_PC2_0.001.Rdata")
rbPal <- colorRampPalette(c("green","orange","red"))
adpd<-ltl
plot(u, col="grey",xlim=c(0,1.5*max(adpd$x)),ylim=c(0,1.5*max(adpd$y)),main=paste0("Joint exceedance curves for p = 0.001"))
legend(1.1*max(adpd$x),1.2*max(adpd$y),legend=c("True","L&T", "H&T","GumbelCop","NormalCop","FGMCop","GalambosCop"),
       col=c("blue", "green", "purple","red","red","red","red"), lty=c(1,1,1,3,4,5,6), cex=.9, lwd=2)
for(tt in 1:jo){
  HTlinec <- rbPal(10)[as.numeric(cut(100*HTlined[,tt],breaks = 10))]
  points(HTlinex[,tt],HTliney[,tt],type="l", col="purple",lwd=1, pch=16)
  lines(LTlinex[,tt],LTliney[,tt],type="l", col="darkgreen",lwd=1)
  Glinec <- rbPal(10)[as.numeric(cut(100*Gumlined[,tt],breaks = 10))]
  points(Gumlinex[,tt],Gumliney[,tt],type="l", col="red")
  points(Normlinex[,tt],Normliney[,tt],type="l", col="darkred")
  points(FGMlinex[,tt],FGMliney[,tt],type="l", col="darkred")
  lines(Gallinex[,tt],Galliney[,tt],type="l", col="darkred")
}

formL1<-matrix(ncol=3,nrow=10000)
pfinal<-list(formL1, formL1,formL1,formL1, formL1,formL1)
nsim=R
for (lis in 1:6){
  htpf<-c()  
  for(tt in 1:nsim){
    htp<-cbind(linesfinalx[[lis]][,tt],linesfinaly[[lis]][,tt],linesfinald[[lis]][,tt])
    htpf<-rbind(htpf,htp)
  }
  HTlinec <- rbPal(6)[as.numeric(cut(htpf[,3],breaks = 6))]
  plot(htpf[,1],htpf[,2],col=HTlinec,lwd=3, type="p")
  pfinal[[lis]]<-htpf
}


transflines<-function (linesx=linesfinalx,linesy=linesfinaly){
  craps<-list()
  for (ll in 1:length(linesx)){
    linex<-linesx[[ll]]
    liney<-linesy[[ll]]
    linescore<-linex*liney
    craps<-c(craps,list(min(linex),max(linex),min(liney),max(liney)))
  }
  craps
}
thelines<-transflines(linesx=linesfinalx,linesy=linesfinaly)


formL<-matrix(ncol=3,nrow=100)
linesICx<-list(formL, formL,formL,formL, formL,formL)
linesICy<-list(formL, formL,formL,formL, formL,formL)
linesICd<-list(formL, formL,formL,formL, formL,formL)


EtaM<-matrix(ncol=3,nrow=6)
ChiM<-matrix(ncol=3,nrow=6)

for (los in 1:6){ 


  EtaM[los,1]<-median(as.numeric(totaleta[,los]))
  EtaM[los,2]<-quantile(as.numeric(totaleta[,los]),.025,na.rm=T)
  EtaM[los,3]<-quantile(as.numeric(totaleta[,los]),.975,na.rm=T)
  

  ChiM[los,1]<-median(as.numeric(totalchi[,los]))
  ChiM[los,2]<-quantile(as.numeric(totalchi[,los]),.025,na.rm=T)
  ChiM[los,3]<-quantile(as.numeric(totalchi[,los]),.975,na.rm=T)
  
  
  x1<-c()
  x2<-c()
  x3<-c()
  for (ios in 1:100)
  {
    x1<-c(x1,as.numeric(linesfinalx[[los]][,ios]))
    x2<-c(x2,as.numeric(linesfinaly[[los]][,ios]))
    x3<-c(x3,as.numeric(linesfinald[[los]][,ios]))
    

  }
if(los==6 | los==3 & length(which(x1>31))>0){
  x2<-x2[-which(x1>31)]
  x3<-x3[-which(x1>31)]
  x1<-x1[-which(x1>31)]}
  
  scs1<-scale(as.numeric(x1),center=F)
  ss1<-attr(scs1,'scaled:scale')
  scs2<-scale(as.numeric(x2),center=F)
  ss2<-attr(scs2,'scaled:scale')
  x1p<-as.vector(scale(as.numeric(x1),center=F))
  x2p<-as.vector(scale(as.numeric(x2),center=F))
  poco<-as.data.frame(cbind(x1p,x2p))
  
  if (interh=="comb"){
  polco<-cart2pol(x1p, x2p)
  lala<-seq(max(min(polco[,2],na.rm=T),0),max(polco[,2],na.rm=T),length.out = 100)
  quantile(polco[,2],0.95,na.rm=T)
  for(al in 1:100)
  {
    alo<-lala[al]
    
    pre<-matrix(ncol=3,nrow=100)
    for (ibs in 1:100){
      
      y1<-as.numeric(linesfinalx[[los]][,ibs])/ss1
      y2<-as.numeric(linesfinaly[[los]][,ibs])/ss2
      y3<-as.numeric(linesfinald[[los]][,ibs])
      if(min(y2,na.rm=T)<=10){
      fshit<-cbind(y1,y2,y3)
      polco<-cart2pol(y1, y2)
      lapol<-polco[,2]
      sel<-abs(lapol-alo)
      if(length(which(sel==min(sel,na.rm = T)))>1){
        pre[ibs,]<-c(pol2cart(polco[which(sel==min(sel,na.rm=T))[1],1],alo),fshit[which(sel==min(sel,na.rm=T))[1],3])
      }else{pre[ibs,]<-c(pol2cart(polco[which(sel==min(sel,na.rm=T)),1],alo),fshit[which(sel==min(sel,na.rm=T)),3])}
      pre[ibs,1]<-pre[ibs,1]*ss1
      pre[ibs,2]<-pre[ibs,2]*ss2}
    }
    alex<-as.data.frame(pre)
    linesICx[[los]][al,1]<-median(alex[,1],na.rm=T)
    linesICy[[los]][al,1]<-median(alex[,2],na.rm=T)
    linesICd[[los]][al,1]<-median(alex[,3],na.rm=T)
    linesICx[[los]][al,2]<-quantile(alex[,1],0.025,na.rm=T)
    linesICx[[los]][al,3]<-quantile(alex[,1],0.975,na.rm=T)
    linesICy[[los]][al,2]<-quantile(alex[,2],0.025,na.rm=T)
    linesICy[[los]][al,3]<-quantile(alex[,2],0.975,na.rm=T)
  }
  } 
  if(interh=="casc"){
    polco<-data.frame(x1p, x2p)
  lala<-seq(max(min(polco[,1],na.rm=T),0),max(polco[,1],na.rm=T),length.out = 100)
  quantile(polco[,1],0.95,na.rm=T)
  for(al in 1:100)
  {
    alo<-lala[al]
    
    pre<-matrix(ncol=3,nrow=100)
    for (ibs in 1:100){
      
      y1<-as.numeric(linesfinalx[[los]][,ibs])/ss1
      y2<-as.numeric(linesfinaly[[los]][,ibs])/ss2
      y3<-as.numeric(linesfinald[[los]][,ibs])
        fshit<-cbind(y1,y2,y3)
        polco<-cbind(y1, y2)
        lapol<-polco[,1]
        sel<-abs(lapol-alo)
        if(length(which(sel==min(sel,na.rm = T)))>1){
          pre[ibs,]<-c(alo,polco[which(sel==min(sel,na.rm=T))[1],2],fshit[which(sel==min(sel,na.rm=T))[1],3])
        }else{pre[ibs,]<-c(alo,polco[which(sel==min(sel,na.rm=T)),2],fshit[which(sel==min(sel,na.rm=T)),3])}
        pre[ibs,1]<-pre[ibs,1]*ss1
        pre[ibs,2]<-pre[ibs,2]*ss2}
    
    alex<-as.data.frame(pre)
    linesICx[[los]][al,1]<-median(alex[,1],na.rm=T)
    linesICy[[los]][al,1]<-median(alex[,2],na.rm=T)
    linesICd[[los]][al,1]<-median(alex[,3],na.rm=T)
    linesICx[[los]][al,2]<-quantile(alex[,1],0.025,na.rm=T)
    linesICx[[los]][al,3]<-quantile(alex[,1],0.975,na.rm=T)
    linesICy[[los]][al,2]<-quantile(alex[,2],0.025,na.rm=T)
    linesICy[[los]][al,3]<-quantile(alex[,2],0.975,na.rm=T)
  }
  }
  
} 



totaletaT<-rbind(totaletaT,EtaM)
totalchiT<-rbind(totalchiT,ChiM)
linesfinal<-list(linesfinalx,linesfinaly,linesfinald)
names(linesfinal)=c("x","y","d")
linesIC<-list(linesICx,linesICy, linesICd)
names(linesIC)=c("x","y","d")

# save(linesfinal,file=paste0(dir.coef,"/lines_W_R_HC_0.001.Rdata"))
# save(linesIC,file=paste0(dir.coef,"/linesIC95_W_R_HC_0.001.Rdata"))

save(linesfinal,file=paste0(dir.coef,"/lines_T_F_PC2_0.001.Rdata"))
save(linesIC,file=paste0(dir.coef,"/linesIC95_T_F_PC2_0.001.Rdata"))

totalchiT=as.data.frame(totalchiT)
names(totalchiT)=c("median","Q2.5","Q97.5")

totaletaT=as.data.frame(totaletaT)
names(totaletaT)=c("median","Q2.5","Q97.5")
graphout<-list (t(totalchiT),t(totaletaT))
names(graphout)=c("Chival","Etaval")

# save(graphout, file=paste0(dir.day,"/stat_output_W_R_HC_0.001.Rdata"))
save(graphout, file=paste0(dir.day,"/stat_output_T_F_PC2_0.001.Rdata"))


###################################FINAL PLOT #####################################
#loader les functions
# load(file="C:/Users/PhD Student/OneDrive - King's College London/DATA/CODE ANALYSIS/Model comparison/out/2019-11-12/stat_output_T_F_PC2_0.001.Rdata")
# load(file="C:/Users/PhD Student/OneDrive - King's College London/DATA/CODE ANALYSIS/Model comparison/out/2019-11-12/bootstraplife/linesIC95_T_F_PC2_0.001.Rdata")
linesICx=linesIC$x
linesICy=linesIC$y
linesICd=linesIC$d

idmods<-c(1,2,4,5)
bl<-linesInit$LT
epy<-length(fire01meantemp$temp2.temperature)/26
1/(0.0001*epy)

plot(u, col="grey",xlim=c(0,1.2*max(bl$x)),ylim=c(0,1.4*max(bl$y)),xlab=xl,ylab=yl,cex.axis=1.4,cex.lab=1.5)


# grid (NULL,NULL, lty = 1, col = "lightgray")
legend(.75*max(bl$x),1.5*max(bl$y),legend=mods[idmods], bg=NULL,
       col=cols[idmods], cex=1.3, lwd=3,bty="n",lty=idmods+9)
for (zg in idmods){
points(linesInit[[zg]]$x,linesInit[[zg]]$y,type="l" ,lwd=3,lty=zg+9, col=cols[zg])
}
cols2<-c("purple","green","red")
rbPal <- colorRampPalette(c("green","orange","red"))
# lines(linesIC$x[[zg]][,3],linesIC$y[[zg]][,3],type="l", lwd=1,lty=2,col=cols[zg])




for (zg in idmods){
  # coli <- rbPal(6)[as.numeric(cut(linesICd[[zg]][,1],breaks = 6))]
  # lines(linesICx[[zg]][,2],linesICy[[zg]][,2],type="l" ,lwd=1,lty=2,col=coland[zg])
  # lines(linesICx[[zg]][,3],linesICy[[zg]][,3],type="l" ,lwd=1,lty=2,col=coland[zg])
  polygon(c(linesICx[[zg]][,2],rev(linesICx[[zg]][,3])),c(linesICy[[zg]][,2],rev(linesICy[[zg]][,3])),col = alpha(coland[zg],.2), border = FALSE)
}




1/((3442/26)*0.001)
1/1000*365
