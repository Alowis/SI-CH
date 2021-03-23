
################## DIAGNOSTIC FUNCTION ######################
Diagnostix<-function(obs,model,devplot=F){
  obsxN<-obs[which(obs[,3]>denslim),1]/max(obs[,1],na.rm = T)
  obsyN<-obs[which(obs[,3]>denslim),2]/max(obs[,2],na.rm = T)
  obsyd<-obs[which(obs[,3]>denslim),3]
  obsN<-cbind(obsxN,obsyN)
  modxN<-model[,1]/max(obs[,1],na.rm = T)
  modyN<-model[,2]/max(obs[,2],na.rm = T)
  modN<-cbind(modxN,modyN)
  obsp<-cart2pol(obsN[,1],obsN[,2])
  modp<-cart2pol(modN[,1],modN[,2])
  lala<-seq(max(min(modp[,2]),min(obsp[,2])),min(max(modp[,2]),max(obsp[,2])),length.out = 80)
  denra=obs[,3]
  preo<-matrix(ncol=2,nrow=80)
  prem<-matrix(ncol=2,nrow=80)
  denw<-matrix(ncol=1,nrow=80)
  rob<-matrix(ncol=2,nrow=80)
  for(al in 1:80)
  {
    alo<-lala[al]
    x1<-modN[,1]
    x2<-modN[,2]
    mshit<-cbind(x1,x2)
    pmo<-cart2pol(x1,x2)
    lapom<-pmo[,2]
    prem[al,]<-c(approx(pmo[,2],pmo[,1],xout=alo, n=1)$y,approx(pmo[,2],pmo[,1],xout=alo, n=1)$x)
    # sem<-abs(lapom-alo)
    # if(length(which(sem==min(sem,na.rm = T)))>1){
    #   prem[al,]<-mshit[which(sem==min(sem,na.rm = T))[1],] 
    #   print("kikoo")
    # }else{prem[al,]<-mshit[which(sem==min(sem,na.rm = T)),]}
    y1<-obsN[,1]
    y2<-obsN[,2]
    did<-obsyd
    pob<-cart2pol(y1, y2)
    preo[al,]<-c(approx(pob[,2],pob[,1],xout=alo, n=1)$y,approx(pob[,2],pob[,1],xout=alo, n=1)$x)
    denw[al,]<-approx(pob[,2],did,xout=alo, n=1)$y
    # lapob<-pob[,2]
    # seo<-abs(lapob-alo)
    # rob[al]<-pob[which(seo==min(seo)),]
    # preo[al,]<-obsN[which(seo==min(seo)),]
  }
  obs1<-pol2cart(preo[,1],preo[,2])
  model1<-pol2cart(prem[,1],prem[,2])
  if (devplot==T){
  plot(obs1, pch=16, xlim=c(0,1), ylim=c(0,1),col="royalblue",type="o",xlab="X1",ylab="X2")
  points(model1,pch=17,col="red",type="o")
  points(obs1[40,1],obs1[40,2],col="gold",pch=16)
  points(model1[40,1],model1[40,2],col="gold",pch=17)
  for(al in 1:80){
    tix<-cbind(c(0 ,cos(lala[al])),c(0,sin(lala[al])))
    rlm<-lm(tix[,2] ~ tix[,1])
    abline(rlm,col="darkgrey")
  }
  }
  sd1<-(obs1[,1])
  sd2<-(obs1[,2])
  sd3<-(model1[,1]) 
  sd4<-(model1[,2])
  # qqq1=(sqrt((sd3^2+sd4^2))-sqrt((sd1^2+sd2^2)))
  # plot(qqq1,type="h")
  # qq2<-sqrt(mean(qqq1))
  bidis<-mean(sqrt((sd3^2+sd4^2))-sqrt((sd1^2+sd2^2)))
  ecdist=mean(sqrt((sd1-sd3)^2+(sd2-sd4)^2))
  ecdistW<-sum(denw*(sqrt((sd1-sd3)^2+(sd2-sd4)^2)))
  rmsb<-sqrt(mean((sqrt(sd3^2+sd4^2)-sqrt(sd1^2+sd2^2))^2))
  cmvt<-cvm.test(prem[,1],ecdf(preo[,1]))$statistic
  kst<-ks.test(preo[,1],prem[,1])$statistic
  outbig<-c(bidis,ecdist,ecdistW,rmsb,cmvt,kst)
  names(outbig)=c("bias_N","ecdist_N","Wecdist_N","RMS_bias_N","Cvm_N","KS_N")
  return(outbig)
}

################## Curve functions #########################

xcurvDraw<-function(mvd,u2,meanlog1,meanlog2,sdlog1,sdlog2,n.grid1,n.grid2,ms,cop,par,pos,pobj){
  
  cont<-contour(mvd,pMvdc,levels = c(1-pobj), add=T,box=F,col="orange",xlim=c(0,1.2*max(u2$V1)), ylim=c(0,1.2*max(u2$V2)), lwd=3,n.grid=c(n.grid1,n.grid1))
  # ad<-contourLines(cont$x,cont$y, cont$z, levels = c(0.999))
  # adol<-data.frame(ad[[1]]$x,ad[[1]]$y)
  cox=cont$x
  coy<-cont$y
  ax<-cont$z
  mar1=u2[order(u2[,1]),1]
  mar2=u2[order(u2[,2]),2]

  pxc<-plnorm(cont$x,meanlog = meanlog1, sdlog = sdlog1)
  pyc<-plnorm(cont$y,meanlog = meanlog2, sdlog = sdlog2)

  for (k in 1:length(pyc)){
    for (j in 1:length(pxc)){
      ax[j,k]=(1-(pxc[j]+pyc[k]-(cont$z[j,k])))}}
  
  
  contour(cox,coy,ax,level=pobj)
  r<-contourLines(cox,coy, ax, levels = c(pobj)) 
  rs<-data.frame(r[[1]]$x,r[[1]]$y)
  q1<-c(quantile(rs[,1],0.3),quantile(rs[,1],0.8))
  q2<-c(quantile(rs[,2],0.3),quantile(rs[,2],0.8))
  
  
  pxc<-plnorm(mar1,meanlog = meanlog1, sdlog = sdlog1)
  pyc<-plnorm(mar2,meanlog = meanlog2, sdlog = sdlog2)
  
  
  godx<-pxc
  
  gody<-pyc
  
  coxi<-cox
  
  coyi<-coy

  laperle<-surfuncCOP(godx, gody, cop=cop, para=par)
  idloose<-which(laperle<0.01)
  id3<-idloose[1]-1
  v1<-mar1[idloose][1]
  v2<-mar2[idloose][1]
  v3<-mar1[id3]
  v4<-mar2[id3]
  
  
  if(pos=="l"){
    xmin=0
    xmax=v1
    ymin=v2
    ymax=1.2*max(coyi)
  }
  if(pos=="m"){
    xmin=v1
    xmax=1.2*max(coxi)
    ymin=0
    ymax=v2
  }
  if(pos=="t"){
    xmin=v3
    xmax=q1[2]
    ymin=v4
    ymax=q2[2]
  }
  
  cont<-contour(mvd,pMvdc,levels = c(1-pobj), add=F,box=F,col="orange",xlim=c(xmin,xmax), ylim=c(ymin,ymax), lwd=3,n.grid=n.grid2)
  cox<-cont$x
  coy<-cont$y
  ax<-cont$z
  pxc<-plnorm(cont$x,meanlog = meanlog1, sdlog = sdlog1)
  pyc<-plnorm(cont$y,meanlog = meanlog2, sdlog = sdlog2)
  for (k in 1:length(pyc)){
    for (j in 1:length(pxc)){
      ax[j,k]=(1-(pxc[j]+pyc[k]-(cont$z[j,k])))}}
  
  
  contour(cox,coy,ax,level=pobj)
  reference<-contourLines(cox,coy, ax, levels = c(pobj)) 
  
  refc<-data.frame(reference[[1]]$x,reference[[1]]$y)
  plot(refc)
  x=seq(min(refc[,1]),max(refc[,1]),length.out = 120)
  xf<-x
# if(ms>3){
#   print("shit")
#   xmin<-
#   conti<-contour(mvd,pMvdc,levels = c(0.999), add=F,box=F,col="orange",xlim=c(0,median(refc[,1])-0.1), ylim=c(0,1.2*max(refc[,2])), lwd=3,n.grid=c(n.grid2,2*n.grid2))
#   cox=conti$x
#   coy<-conti$y
#   
#   if(margins=="lnorm"){
#     pxc<-plnorm(conti$x,meanlog = meanlog1, sdlog = sdlog1)
#     pyc<-plnorm(conti$y,meanlog = meanlog2, sdlog = sdlog2)
#   }
#   
#   ax<-conti$z
#   for (k in 1:length(pyc)){
#     for (j in 1:length(pxc)){
#       ax[j,k]=(1-(pxc[j]+pyc[k]-(conti$z[j,k])))}}
#   
#   
#   reference1<-contourLines(cox,coy, ax, levels = c(0.001)) 
#   refc1<-data.frame(reference1[[1]]$x,reference1[[1]]$y)
# 
#   
#   contj<-contour(mvd,pMvdc,levels = c(0.999), add=F,box=F,col="orange",xlim=c(max(refc1[,1])+0.1,1.2*max(refc[,1])), ylim=c(0,median(refc[,2])-0.1), lwd=3,n.grid=c(n.grid2,n.grid2))
#   cox=contj$x
#   coy<-contj$y
#   
#   if(margins=="lnorm"){
#     pxc<-plnorm(contj$x,meanlog = meanlog1, sdlog = sdlog1)
#     pyc<-plnorm(contj$y,meanlog = meanlog2, sdlog = sdlog2)
#   }
#   
#   ax<-contj$z
#   for (k in 1:length(pyc)){
#     for (j in 1:length(pxc)){
#       ax[j,k]=(1-(pxc[j]+pyc[k]-(contj$z[j,k])))}}
#   
#   
#   reference2<-contourLines(cox,coy, ax, levels = c(0.001)) 
#   refc2<-data.frame(reference2[[1]]$x,reference2[[1]]$y)
# plot(refc)
# rel<-c(first(refc1[,1]),last(refc1[,1]),first(refc2[,1]),last(refc2[,1]))
# 
#   refcf<-rbind(as.matrix(refc2),as.matrix(refc1))
#   refc<-round(refcf,6)
#   refc=refcf[order(refc[,1],-refc[,2]),]
#   plot(refc,type="l")
#   }
# 
#   rab<-digit.curves.p(start=c(refc[1,1],refc[1,2]), as.matrix(refc), nPoints=98, closed = FALSE)
#   rab<-densicurvcop1(rab,omg)
#   
  rab=refc
  return(rab)
}


curve.funct.a<-function(px,py,mar1,mar2,u,pos,pobje,ng=100,inter="comb",logm=F){
 if(logm==T) {
  mar1<-log(mar1)
  mar2<-log(mar2)
 }
  xmin=0
  xmax=1
  ymin=0
  ymax=1
  tieef<-mean
  ngx=100000


  godx<-approx(px,mar1, n = ngx, method = "linear",
               yleft = min(mar1), yright = max(mar1), ties = tieef)
  
  gody<-approx(py,mar2, n = ngx, method = "linear",
               yleft = min(mar1), yright = max(mar1), ties = tieef)


  # coxi<-approx(godx$x,godx$y, n = ngx, method = "linear",
  #              yleft = min(godx$y), yright= max(godx$y), ties = tieef)
  # 
  # coyi<-approx(gody$x,gody$y, n = ngx, method = "linear",
  #              yleft = min(gody$y), yright = max(gody$y), ties = tieef)


  coxi<-approx(godx$y,godx$x, n = ng, method = "linear",
               yleft = 0, yright = 1, ties = mean)
  
  coyi<-approx(gody$y,gody$x, n = ng, method = "linear",
               yleft = 0, yright = 1, ties = mean)

  godx<-coxi$y
  gody<-coyi$y
  coxi<-coxi$x
  coyi<-coyi$x
  repeat{
    idx<-which(diff(coxi)<=0)
    # print("cl")
    if(length(idx)<1){
      break
    }
    coxi[idx+1]=coxi[idx]+0.001
    
  }
  repeat{
    idy<-which(diff(coyi)<=0)
    if(length(idy)<1){
      break
    }
    coyi[idy+1]=coyi[idy]+0.001
    
  }
  
  acp3<-matrix(NA, nrow = ng, ncol = ng)
  for (k in 1:length(gody)){
    for (j in 1:length(godx)){
      if (inter=="comb"){
        acp3[j,k]=surfuncCOP(godx[j], gody[k], cop=coco, para=c1)}
      if (inter=="casc"){
        acp3[j,k]= surfuncCOP(godx[j], gody[k], cop=coco, para=c1)/(1-godx[j])} 
      
    }
  }
  
  clf<-contourLines(coxi,coyi, acp3, levels = pobje) 
  if(is.list(clf)|length(clf)>1){
    clf<-as.matrix(data.frame(clf[[1]]$x,clf[[1]]$y))} else{clf<-"niketamere"}
  if(logm==T) clf<-exp(clf)
  return(clf)
  
}


################ Multivariate extreme models ################

JT_KDE<-function(u2,pbas ,pobj,beta,margins,vtau,devplot=F,mar1,mar2,px,py){
  
  kk<-data.frame(px,py)
  aa<-kcde(u2,gridsize=400, tail.flag = "upper.tail")
  
  lox<-aa$eval.points[[1]]
  loy<-aa$eval.points[[2]]
  if(margins=="lnorm"){
    pxe<-plnorm(lox,meanlog = meanlog1, sdlog = sdlog1)
    pye<-plnorm(loy,meanlog = meanlog2, sdlog = sdlog2)
  }
  
  if(margins=="emp"){
    pxc<-approx(mar1, px, xout = lox, method = "linear",yleft = min(px),yright = max(px), rule = 1,ties=mean)$y
    pyc<-approx(mar2, py, xout = loy, method = "linear",yleft = min(py),yright = max(py), rule = 1)$y
  }
  
  if(devplot==T){
  plot(aa,cont = c(0.01,0.05,0.1,1),display="filled.contour" ,col=viridis(10))
  
  }
  az<-aa$estimate
  wq<-contourLines(lox,loy, az, levels = c(0.7,0.5,0.2,0.1,pbas))
  wq0ri<-cbind(wq[[5]]$x,wq[[5]]$y)
  
  var1<-u2[,1]
  var2<-u2[,2]
  wu<-data.frame(u2)
  
  p01=0.01
  p02=0.01
  q0=0.9
  
  
  
  
  
  if(length(which(wq0ri[,1]<0))>0)wq0ri<-wq0ri[-which(wq0ri[,1]<0),]
  if(length(which(wq0ri[,2]<0))>0)wq0ri<-wq0ri[-which(wq0ri[,2]<0),]
  if(margins=="lnorm"){
    wqUnitx<- plnorm(wq0ri[,1],meanlog = meanlog1, sdlog = sdlog1)
    wqUnity<- plnorm(wq0ri[,2],meanlog = meanlog2, sdlog = sdlog2)
  }
  
  if(margins=="emp"){
    wqUnitx<- approx(mar1, px, xout = wq0ri[,1], method = "linear",yleft = min(px),yright = max(px), rule = 1)$y
    wqUnity<- approx(mar2, py, xout = wq0ri[,2], method = "linear",yleft = min(py),yright = max(py), rule = 1)$y
  }
  
  #Transforming Pbase to Frechet Margins
  wqxFrechet<--1/log(wqUnitx)
  wqyFrechet<--1/log(wqUnity)
  
  
  pb<-1-pbas
  
  #Setting up Pobjective
  
  
  s=pbas/pobj
  
  
  xFrechet<--1/log(px)
  yFrechet<--1/log(py)
  ModelHugo_file=paste0(getwd(),"/LF/BveLTDep.R")
  source (ModelHugo_file)
  if(vtau<0)q0=0.9
  qq<-c()
  cm<-c()
  qc<-.95
  rq0<-seq(0.7,0.94,by=0.02)
  for(q0 in rq0){
  estims<-try(BveLTDep (data= kk,mod.thresh.u = q0,crit.lev.u = qc,sig.lev=0.05,ci.meth='se',marg.inf=F),silent = T)
  qd<-try((estims$par[2]),silent=T)
  if(is.numeric(qd)){
    qq<-c(qq,qd)
  }
  }
  sumd<-0
  for(i in 1:length(diff(qq))){
    sdiff<-diff(qq)[i]
    sumd<-c(sumd,sumd[i]+sdiff)}
  sh<-which(sumd<=-0.02|sumd>=0.02 )[1]
  if(!is.na(sh)) q0<-rq0[sh-1]
  
  estims<-try(BveLTDep (data= kk,mod.thresh.u = q0,crit.lev.u = qc,sig.lev=0.05,ci.meth='se',marg.inf=F),silent = T)
  chat=NA
  etahat=NA
  Chilow=NA
  Chimed=NA
  try(chat<-estims$par[1],silent=T)
  try(etahat<-estims$par[2],silent=T)
  try(Chilow<-estims$chiCIs[1],silent=T)
  try(Chimed<-estims$chi,silent=T)
  
  #Loop for asymptotic dependence
  
  if (!is.na(Chilow) & (etahat<0.9)){
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
  
  
  
  if(margins=="lnorm"){
    objx<-qlnorm(projbackx,meanlog = meanlog1, sdlog = sdlog1)
    objy<-qlnorm(projbacky,meanlog = meanlog2, sdlog = sdlog2)
  }
  
  if(margins=="emp"){
    objx<-approx(px,mar1, xout = projbackx, method = "linear",yleft = min(mar1),yright = max(mar1), rule = 1)$y
    objy<-approx(py,mar2, xout = projbacky, method = "linear",yleft = min(mar2),yright = max(mar2), rule = 1)$y 
  }
  
  wqobj<-cbind(objx, objy)
  wqobj<-wqobj[order(wqobj[,1]),]
  wqd<-diff(round(wqobj,10))
  swqd<-wqd[,1]+wqd[,2]
  weirdp<-which(abs(swqd)==0)
  if(length(weirdp)>0) wqobj<-wqobj[which(abs(swqd)>0)[-1]-1,]
  # wqobjr<-cbind(objx, objy)
  # frejc<-cbind(projx,projy)
  # wqfr<-cbind(wqxFrechet,wqyFrechet)
  # m1<-min(u2[,1])
  # m2<-min(u2[,2])
  # ad1<-NA
  # ad2<-NA
  # if (m1<min(wqobj[,1]))ad1<-cbind(m1,wqobj[1,2])
  # if (m2<min(wqobj[,2]))ad2<-cbind(wqobj[length(wqobj[,1]),1],m2)
  # wqobj<-rbind(ad1,wqobj,ad2)
  
  
  res<-list(levelcurve=wqobj,etaJT=etahat,chiJT= Chimed,wq0ri=wq0ri)
}

Cond_mod<-function(u2,tr1,tr2,tsim,num.sim,pobj=0.001){
  names(u2)= c("V1","V2")
  thresh1 <- tr1
  thresh2 <- tr2
  
  
  
  if(num.points>50000){
    ext.q=0.995}
  if(num.points<=50000 & num.points>5000){
    ext.q=0.95
  }
  if(num.points<=5000){
    ext.q=0.9
  }
  
  mex.fit <- mex(data = u2 , which = 1, mqu = thresh1, dqu = ext.q, margins = "Laplace", constrain = F) 
  mex.fit2 <- mex(data = u2, which = 2, mqu =c(thresh1), dqu = ext.q, constrain=F) 
  
  
  mex.pred <-predict(mex.fit, pqu = tsim, nsim = num.sim,smoothZdistribution=TRUE)
  mex.pred2<-predict(mex.fit2, pqu = tsim, nsim = num.sim,smoothZdistribution=T)
  
  #Estimation of the H&T chi
  qexp1<-quantile(u2[,2],tsim)
  l1<-length(mex.pred$data$simulated[which(mex.pred$data$simulated[,2]>qexp1),2])
  chilt1<-l1/num.sim
  tault1<-cor.test(x=mex.pred$data$simulated$V1,y=mex.pred$data$simulated$V2,method="kendall")$estimate
  qexp2<-quantile(u2[,1],tsim)
  l2<-length(mex.pred2$data$simulated[which(mex.pred2$data$simulated[,2]>qexp2),2])
  chilt2<-l2/num.sim
  chiht<-mean(chilt1,chilt2)
  tault2<-cor.test(x=mex.pred2$data$simulated$V1,y=mex.pred2$data$simulated$V2,method="kendall")$estimate
  tauht<-mean(tault1,tault2)
  chibarht=2*log(1-tsim)/log(chiht*(1-tsim))-1
  etaht=(chibarht+1)/2
  
  
  #Not used here
  q2<-quantile(u2[,2],tsim)
  q1<-quantile(u2[,1],tsim)
  totsim=num.sim*3/2
  #
  sim1<-mex.pred$data$simulated
  #
  sim2<-mex.pred2$data$simulated
  sim2<-data.frame(mex.pred2$data$simulated$V1,mex.pred2$data$simulated$V2)
  names(sim2)<-names(u2)
  
  sim1<-sim1[which(sim1$V2<(sim1$V1+q2-q1)),]
  sim2<-sim2[which(sim2$V2>(sim2$V1+q2-q1)),]
  # mex.pred2$data$simulated=sim2
  # mex.pred$data$simulated=sim1
  onlysim<-rbind(sim1,sim2)
  u2<-as.data.frame(u2)
  names(onlysim)=names(u2)
  simu<-rbind(u2,onlysim)
  
  
  
  
  JointExcureve_file=paste0(getwd(),"/Multimodel_Sim/Joint_excurves.R")
  source (JointExcureve_file)
  pobj1=pobj/(1-tsim)
  j1 <- JointExceedanceCurve(mex.pred,pobj,n=120) 
  j2 <- JointExceedanceCurve(mex.pred2,pobj,n=120)
  
  j1o<-data.frame(j1[[1]],j1[[2]])
  j2o<-data.frame(j2[[1]],j2[[2]])
  
  ##########Level curves drawing########################
  
  names(j1o)=names(j2o)
  jline<-rbind(j1o[which(j1o[,2]<j1o[,1]+q2-q1),],j2o[which(j2o[,2]>j2o[,1]+q2-q1),])
  jline<-jline[order(jline$j2..1..),]
  jj1<-as.matrix(j1o[which(j1o[,2]<j1o[,1]+q2-q1),])
  jj2<-as.matrix(j2o[which(j2o[,2]>j2o[,1]+q2-q1),])
  mj1<-min(jj1[,1])
  mj2<-max(jj2[,1])
  mj3<-max(jj1[,1])
  if(mj2>mj1){
    t1<- jj1[which(jj1[,1]>mj1&jj1[,1]<mj2),]
    t2<- jj2[which(jj2[,1]>mj1&jj2[,1]<mj2),]
    if(mj3<mj2)mj2=mj3
    xjj<-seq(mj1,mj2,length.out = 50)
    cj1=approx(jj1[,1], jj1[,2], xout = xjj, method = "linear",yleft = max(jj1[,2]),yright = min(jj1[,2]), rule = 1)$y
    cj2=approx(jj2[,1], jj2[,2], xout = xjj, method = "linear",yleft = max(jj2[,2]),yright = min(jj2[,2]), rule = 1)$y
    cjf<-(cj1+cj2)/2
    jj1<-jj1[order(jj1[,2]),]
    jjcrop<-cbind(xjj,cjf)
    jline<-rbind(jj1[1,],jj1[which(jj1[,1]>=(mj2)),],jjcrop,jj2[which(jj2[,1]<mj1),])
    
  } 
  jline<-jline[order(jline[,1]),]
  repeat{
    idm<-which(diff(jline[,1])<0.001 & diff(jline[,2])>0)
    if(length(idm)<1){
      break
    }
    jline<-jline[-(idm),]
    
  }
  difj<-diff(jline[,2])+diff(jline[,1])
  if(length(which(difj==0))>0)jline<-jline[-which(difj==0),]
  
  m1<-min(u2[,1])
  m2<-min(u2[,2])
  ad1<-NA
  ad2<-NA
  if (m1<min(jline[,1]))ad1<-cbind(m1,jline[1,2])
  if (m2<min(jline[,2]))ad2<-cbind(jline[length(jline[,1]),1],m2)
  jline<-as.matrix(jline)
  jline<-rbind(ad1,jline,ad2)
  jline<-removeNA(jline)
  jline<-data.frame(jline)
  xht=seq(ceiling(min(jline[,1])),round(max(jline[,1]),2)-0.1,length.out = 120)
  
  res=list(jline=jline,chiHT=chiht,etaht=etaht,onlysim=onlysim)
}

################# Others ##################################
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

densicurvem<-function(kdetab,lines,tl,lines2,devplot=F){
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
  if(devplot==T)plot(posoncurve,densim, type="h")
  if (tl=="l"){
    exit<-data.frame(lines2,densim)}
  if (tl=="h"){
    exit<-data.frame(lines,densim)}
  names(exit)=c("x","y","dens")
  exit
}

densicurvcop1<-function(lines, mvf,devplot=F){
  shit<-dMvdc(x=as.matrix(data.frame(lines)), mvdc=mvf)
  shma<-sum(shit)
  shit<-shit/shma
  posonshit<-seq(1/length(shit),1,length=length(shit))
  if(devplot==T)plot(posonshit,shit, type="h")
  shuba<-data.frame(lines,shit)
  names(shuba)=c("x","y","dens")
  shuba
}

densicurvcop2<-function(lines, copi){
  linep<-lines
  linep[,1]<-spline(u2[,1],pxf, n = 300, method = "fmm",
                    xmin = min(lines[,1]), xmax = max(lines[,1]), ties = mean,xout=lines[,1])$y
  linep[,2]<-spline(u2[,2],pyf, n = 300, method = "fmm",
                    xmin = min(lines[,2]), xmax = max(lines[,2]), ties = mean,xout=lines[,2])$y
  shit<-copula::dCopula(u=as.matrix(data.frame(linep)), copula =copi)
  shma<-sum(shit)
  shit<-shit/shma
  posonshit<-seq(1/length(shit),1,length=length(shit))
  shuba<-data.frame(lines,shit)
  names(shuba)=c("x","y","dens")
  shuba
}

digit.curves.p <- function(start, curve, nPoints, closed=TRUE){
  nPoints <- nPoints+2
  if(!is.matrix(curve)) stop("Input must be a p-x-k matrix of curve coordinates")
  nCurvePoints = NROW(curve)
  if(nCurvePoints < 2) stop("curve matrix does not have enough points to estimate any interior points")
  if(nPoints > (nCurvePoints - 1)) {
    if((nCurvePoints - 1) == 1) nPoints = 1
    if((nCurvePoints - 1) > 1) nPoints = nCurvePoints - 2
    cat("\nWarning: because the number of desired points exceeds the number of curve points,")
    cat("\nthe number of points will be truncated to", nPoints, "\n\n")
  }
  start <- as.numeric(start)
  if(!setequal(start, curve[1,])) curve <- rbind(start, curve)
  if(closed) curve <- rbind(curve, curve[1,])
  res <- evenPts.p(curve, nPoints)
  if(closed) res <- res[-NROW(res),]
  res
}

evenPts.p <- function(x, n){
  x <- as.matrix(na.omit(x))
  # x<-round(x,6)
  # x <-x[order(x[,1],-x[,2]),]
  at1<-scale(x[,1])
  at2<-scale(x[,2])
  x[,1]<-as.vector(scale(x[,1]))
  x[,2]<-as.vector(scale(x[,2]))
  N <- NROW(x); p <- NCOL(x)
  if(N == 1) stop("x must be a matrix")
  if(n < 3) {
    n <- 2
    nn <- 3 # so lapply function works
  } else nn <- n
  
  if(N == 2) {
    x <- rbind(x, x[2,])
    N <- 3 # third row cut off later
  }
  xx <- x[2:N, ] - x[1:(N - 1), ]
  ds <- as.numeric(sqrt(xx[,1]^2+xx[,2]^2))
  cds <- c(0, cumsum(ds))
  cuts <- cumsum(rep(cds[N]/(n-1), n-1))
  targets <- lapply(1:(nn-2), function(j){
    dtar <- cuts[j]
    ll <- which.max(cds[cds < dtar])
    ul <- ll + 1
    adj <- (dtar- cds[ll])/(cds[[ul]] - cds[ll])
    x[ll,] + adj * (x[ul,] - x[ll,])
  })
  
  out <- matrix(c(x[1,], unlist(targets), x[N,]), n, p, byrow = TRUE)
  out[,1]<- out[,1]*attr(at1,'scaled:scale')+attr(at1, 'scaled:center')
  out[,2]<-out[,2]*attr(at2,'scaled:scale')+attr(at2, 'scaled:center')
  out
}

curve.funct.aold<-function(px,py,mar1,mar2,u,pos,pobje,ng=100,inter="comb",logm=F){
  if(logm==T) {
    mar1<-log(mar1)
    mar2<-log(mar2)
  }
  xmin=0
  xmax=1
  ymin=0
  ymax=1
  # clx<-c()
  # for(dd in 1:(length(disc1)-1))
  # {
  #   xmin=disc1[dd]
  #   xmax=disc1[dd+1]
  #   ymin=0
  #   ymax=1
  #   ngx=1000000
  #   tieef<-mean
  #   godx<-spline(px,mar1, n = ngx, method = "natural",
  #                xmin = xmin, xmax = xmax, ties = tieef)
  #   
  #   gody<-spline(py,mar2, n = ngx, method = "natural",
  #                xmin = ymin, xmax = ymax, ties = tieef)
  #   
  #   # coxi<-spline(godx$y,godx$x, n = ng, method = "natural",
  #   #              xmin = min(godx$y), xmax = max(godx$y), ties = tieef)
  #   # 
  #   # coyi<-spline(gody$y,gody$x, n = ng, method = "natural",
  #   #              xmin = min(gody$y), xmax = max(gody$y), ties = tieef)
  # 
  #   coxi<-approx(godx$y,godx$x, n = ng, method = "linear",
  #                yleft = 0, yright = 1, ties = mean)
  #   
  #   coyi<-approx(gody$y,gody$x, n = ng, method = "linear",
  #                yleft = 0, yright = 1, ties = mean)
  #   # coxi<-godx$y
  #   # coyi<-gody$y
  #   # godx<-godx$x
  #   # gody<-gody$x
  # 
  #   godx<-coxi$y
  #   gody<-coyi$y
  #   coxi<-coxi$x
  #   coyi<-coyi$x
  #   plot(godx,coxi)
  #   acp3<-matrix(NA, nrow = ng, ncol = ng)
  #   for (k in 1:length(gody)){
  #     for (j in 1:length(godx)){
  #       if (inter=="comb"){
  #         acp3[j,k]=surfuncCOP(godx[j], gody[k], cop=coco, para=c1)}
  #       if (inter=="casc"){
  #         acp3[j,k]= surfuncCOP(godx[j], gody[k], cop=coco, para=c1)/(1-godx[j])}
  # 
  #     }
  #   }
  # 
  #   clf<-contourLines(coxi,coyi, acp3, levels = pobje)
  #   if(is.list(clf)|length(clf)>1){
  #     clf<-as.matrix(data.frame(clf[[1]]$x,clf[[1]]$y))} else{clf<-"niketamere"}
  #   plot(clf)
  #   min(clf[,1])
  #   clx<-rbind(clx,clf)
  # }
  # for(dd in 1:(length(disc1)-1))
  # {
  #   ymin=disc1[dd]
  #   ymax=disc1[dd+1]
  #   xmin=0
  #   xmax=1
  #   ngx=1000000
  #   tieef<-mean
  #   godx<-spline(px,mar1, n = ngx, method = "natural",
  #                xmin = xmin, xmax = xmax, ties = tieef)
  #   
  #   gody<-spline(py,mar2, n = ngx, method = "natural",
  #                xmin = ymin, xmax = ymax, ties = tieef)
  #   
  #   coxi<-spline(godx$y,godx$x, n = ng, method = "natural",
  #                xmin = min(godx$y), xmax = max(godx$y), ties = tieef)
  #   
  #   coyi<-spline(gody$y,gody$x, n = ng, method = "natural",
  #                xmin = min(gody$y), xmax = max(gody$y), ties = tieef)
  #   
  #   
  #   # coxi<-godx$y
  #   # coyi<-gody$y
  #   # godx<-godx$x
  #   # gody<-gody$x
  #   
  #   godx<-coxi$y
  #   gody<-coyi$y
  #   coxi<-coxi$x
  #   coyi<-coyi$x
  #   plot(godx,coxi)
  #   acp3<-matrix(NA, nrow = ng, ncol = ng)
  #   for (k in 1:length(gody)){
  #     for (j in 1:length(godx)){
  #       if (inter=="comb"){
  #         acp3[j,k]=surfuncCOP(godx[j], gody[k], cop=coco, para=c1)}
  #       if (inter=="casc"){
  #         acp3[j,k]= surfuncCOP(godx[j], gody[k], cop=coco, para=c1)/(1-godx[j])}
  #       
  #     }
  #   }
  #   
  #   clf<-contourLines(coxi,coyi, acp3, levels = pobje)
  #   if(is.list(clf)|length(clf)>1){
  #     clf<-as.matrix(data.frame(clf[[1]]$x,clf[[1]]$y))} else{clf<-"niketamere"}
  #   plot(clf)
  #   clx<-rbind(clx,clf)
  # }
  # plot(clx)
  ngx=10000
  tieef<-mean
  # godx<-spline(px,mar1, n = ngx, method = "natural",
  # #              xmin = xmin, xmax = xmax, ties = tieef)
  # # 
  # # gody<-spline(py,mar2, n = ngx, method = "natural",
  # #              xmin = ymin, xmax = ymax, ties = tieef)
  # godx<-approx(px,mar1, n = ngx, method = "linear",
  #              yleft = min(mar1), yright = max(mar1), ties = tieef)
  # 
  # gody<-approx(py,mar2, n = ngx, method = "linear",
  #              yleft = min(mar1), yright = max(mar1), ties = tieef)
  # 
  # for(god in 1:100){
  # 
  # godx<-approx(godx$x,godx$y, n = ngx, method = "linear",
  #              yleft = min(mar1), yright = max(mar1), ties = tieef)
  # 
  # gody<-approx(gody$x,gody$y, n = ngx, method = "linear",
  #              yleft = min(mar1), yright = max(mar1), ties = tieef)
  # }
  # godx<-approx(godx$x,godx$y, n = ngx, method = "linear",
  #              yleft = min(mar1), yright = max(mar1), ties = tieef)
  # 
  # gody<-approx(gody$x,gody$y, n = ngx, method = "linear",
  #              yleft = min(mar1), yright = max(mar1), ties = tieef)
  # coxi<-spline(godx$y,godx$x, n = ng, method = "natural",
  #              xmin = min(godx$y), xmax = max(godx$y), ties = tieef)
  # 
  # coyi<-spline(gody$y,gody$x, n = ng, method = "natural",
  #              xmin = min(gody$y), xmax = max(gody$y), ties = tieef)
  # 
  # godx<-coxi$y
  # gody<-coyi$y
  # coxi<-coxi$x
  # coyi<-coyi$x
  # laperle<-surfuncCOP(godx, gody, cop=coco, para=c1)
  # idloose<-which(laperle<0.005)
  # id3<-idloose[1]-1
  # v1<-godx[idloose][1]
  # v2<-gody[idloose][1]
  # v3<-godx[id3]
  # v4<-gody[id3]
  # 
  # 
  # if(pos=="l"){
  #   xmin=0
  #   xmax=v1
  #   ymin=v2
  #   ymax=max(gody)
  # }
  # if(pos=="m"){
  #   xmin=v1
  #   xmax=max(godx)
  #   ymin=0
  #   ymax=v2
  # }
  # if(pos=="t"){
  #   xmin=0
  #   xmax=1
  #   ymin=0
  #   ymax=1
  # }
  ngx=100000
  
  
  godx<-approx(px,mar1, n = ngx, method = "linear",
               yleft = min(mar1), yright = max(mar1), ties = tieef)
  
  gody<-approx(py,mar2, n = ngx, method = "linear",
               yleft = min(mar1), yright = max(mar1), ties = tieef)
  
  
  coxi<-approx(godx$x,godx$y, n = ngx, method = "linear",
               yleft = min(godx$y), yright= max(godx$y), ties = tieef)
  
  coyi<-approx(gody$x,gody$y, n = ngx, method = "linear",
               yleft = min(gody$y), yright = max(gody$y), ties = tieef)
  
  
  coxi<-approx(godx$y,godx$x, n = ng, method = "linear",
               yleft = 0, yright = 1, ties = mean)
  
  coyi<-approx(gody$y,gody$x, n = ng, method = "linear",
               yleft = 0, yright = 1, ties = mean)
  
  godx<-coxi$y
  gody<-coyi$y
  coxi<-coxi$x
  coyi<-coyi$x
  repeat{
    idx<-which(diff(coxi)<=0)
    # print("cl")
    if(length(idx)<1){
      break
    }
    coxi[idx+1]=coxi[idx]+0.001
    
  }
  repeat{
    idy<-which(diff(coyi)<=0)
    if(length(idy)<1){
      break
    }
    coyi[idy+1]=coyi[idy]+0.001
    
  }
  
  acp3<-matrix(NA, nrow = ng, ncol = ng)
  for (k in 1:length(gody)){
    for (j in 1:length(godx)){
      if (inter=="comb"){
        acp3[j,k]=surfuncCOP(godx[j], gody[k], cop=coco, para=c1)}
      if (inter=="casc"){
        acp3[j,k]= surfuncCOP(godx[j], gody[k], cop=coco, para=c1)/(1-godx[j])} 
      
    }
  }
  
  clf<-contourLines(coxi,coyi, acp3, levels = pobje) 
  # contour(coxi,coyi, acp3)
  if(is.list(clf)|length(clf)>1){
    clf<-as.matrix(data.frame(clf[[1]]$x,clf[[1]]$y))} else{clf<-"niketamere"}
  if(logm==T) clf<-exp(clf)
  return(clf)
  
}
