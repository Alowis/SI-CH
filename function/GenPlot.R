
gendist<-function(AIndep=F, coef=0.3, nbpt=5000, margin="lnorm",pm1=c(2,0.5),pm2=c(2,.5),ctype="comb"){

loc1=18
scale1=1.353
shape1=-0.1
loc2=13.9
scale2=1.865
shape2=0.0506
meanlog1=pm1[1]
sdlog1=pm1[2]
meanlog2=pm2[1]
sdlog2=pm2[2]
if (AIndep==T){
rhx=2*coef-1
  eta=coef
  rho=rhx
  vtau<-2*asin(rho)/pi
  gc<-normalCopula(rho, dim=2)
}
if (AIndep==F){
  thx=1/log2(2-coef) 
  chi=coef
  theta=thx
  vtau<-1 - 1/theta
  gc<-gumbelCopula(theta, dim=2)
}
  num.points=nbpt
if (margin=="lnorm") {
    omg<-mvdc(gc, margins = c("lnorm","lnorm"), paramMargins= list(list(meanlog = meanlog1, sdlog = sdlog1), list(meanlog = meanlog2, sdlog = sdlog2)), marginsIdentical = F, check = TRUE, fixupNames = TRUE)}
if (margin=="other"){  
    shp=1.5
    sca=4.2
    rt=1.5
    omg<-mvdc(gc, margins = c("weibull","exp"), paramMargins= list(list(shape = 1.5, scale = 4.2), list(rate=1.5)), marginsIdentical = F, check = TRUE, fixupNames = TRUE)
  } 
    u2<-rMvdc(num.points,omg) #synthetic data
    u2<-as.data.frame(u2)
    id<-seq(1,length(u2[,2]))
    u3<-cbind(u2,id)
    dev.new()
    plot(u2)
    tau(gc)
    data=u2
    
    if(margin=="lnorm"){
      px<-plnorm(u2[,1],meanlog = meanlog1, sdlog = sdlog1)
      py<-plnorm(u2[,2],meanlog = meanlog2, sdlog = sdlog2)
    }
    
    if(margin=="other"){
      px<-pweibull(u2[,1],shape= shp, scale = sca)
      py<-pexp(u2[,2],rate = rt)
    }

    kk<-data.frame(px,py)


    n.grid=100
    
    #Utiliser Mvdc pour ca car sinon les contours deviennes tout petit dans les extremes
    cont<-contour(omg,pMvdc,levels = c(0.3,0.5,0.8,0.9,0.995,0.999), add=T, delta=0.01,box=F,nlevels=50,xlim=c(0,2*max(u2$V1)), ylim=c(0,2*max(u2$V2)), lwd=.01,col="#0000ff00",n.grid=300)
    dev.off()
    cont$z[which(cont$z==0.3)]
    ad<-contourLines(cont$x,cont$y, cont$z, levels = c(0.3))
    
    cox=cont$x
    coy<-cont$y

    if(margin=="lnorm"){
      pxc<-plnorm(cont$x,meanlog = meanlog1, sdlog = sdlog1)
      pyc<-plnorm(cont$y,meanlog = meanlog2, sdlog = sdlog2)
    }
    
    if(margin=="other"){
      pxc<-pweibull(cont$x,shape= shp, scale = sca)
      pyc<-pexp(cont$y,rate = rt)
    }

    ax<-cont$z
    if(ctype=="comb"){ 
    for (k in 1:length(pyc)){
      for (j in 1:length(pxc)){
        ax[j,k]=1-(pxc[j]+pyc[k]-(cont$z[j,k]))}}
    }
    if(ctype=="casc"){
    for (k in 1:length(pyc)){
      for (j in 1:length(pxc)){
        ax[j,k]=(1-(pxc[j]+pyc[k]-(cont$z[j,k])))/(1-pxc[j])}}
    }
    
    #ax is the joint probability from the copula
  sortie<- list(u2,cont,ax)
  # wupwup<-qcbvnonpar(p = seq(0.991, 0.999, 0.001), data=u2, epmar = T, nsloc1 =
  #              NULL, nsloc2 = NULL, mint = 1, method = "cfg", convex = FALSE, plot=T,add=T)
  names(sortie)=c("u2","cont","ax")
  return(sortie)
} 
   