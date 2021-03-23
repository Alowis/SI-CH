curve.funct<-function(pxf,pyf,mar1,mar2,u,pos,pobje,ng=100,inter="comb"){
  if(pos=="l"){
    xmin=0
    xmax=0.9
    ymin=0.99
    ymax=1
  }
  if(pos=="m"){
    xmin=0.9
    xmax=0.995
    ymin=0.9
    ymax=.995
  }
  if(pos=="r"){
    xmin=0.99
    xmax=1
    ymin=0
    ymax=.9
  }
  godx<-spline(pxf,mar1, n = ng, method = "fmm",
               xmin = xmin, xmax = xmax, ties = mean)
  
  gody<-spline(pyf,mar2, n = ng, method = "fmm",
               xmin = ymin, xmax = ymax, ties = mean)
  
  coxi<-approx(godx$y,godx$x, n = ng, method = "linear",
               yleft = 0, yright = 1, ties = mean)
  
  coyi<-approx(gody$y,gody$x, n = ng, method = "linear",
               yleft = 0, yright = 1, ties = mean)
  plot(coxi$x,coxi$y)
  godx<-coxi$y
  gody<-coyi$y
  coxi<-coxi$x
  coyi<-coyi$x
  
  repeat{
    idx<-which(diff(coxi)<=0)
    print("cl")
    if(length(idx)<1){
      break
    }
    coxi[idx+1]=coxi[idx]+0.001
    
  }
  repeat{
    idy<-which(diff(coyi)<=0)
    print("cl")
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
  
  cl2<-contourLines(coxi,coyi, acp3, levels = pobje) 
  if(length(cl2)>0){
    cl2<-as.matrix(data.frame(cl2[[1]]$x,cl2[[1]]$y))} else{cl2<-NA}
  
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

curve.funct.b<-function(pxf,pyf,mar1,mar2,u,pos,pobje,ng=100,inter="comb"){
  if(pos=="l"){
    xmin=0
    xmax=0.99
    ymin=0.99
    ymax=1
  }
  if(pos=="m"){
    xmin=0.98
    xmax=1
    ymin=0.98
    ymax=1
  }
  if(pos=="r"){
    xmin=0.99
    xmax=1
    ymin=0
    ymax=0.99
  }
  ngx=10000
  godx<-spline(pxfp,mar1, n = ngx, method = "natural",
               xmin = xmin, xmax = xmax, ties = mean)
  
  gody<-spline(pyfp,mar2, n = ngx, method = "natural",
               xmin = ymin, xmax = ymax, ties = mean)
  
  coxi<-approx(godx$y,godx$x, n = ng, method = "linear",
               yleft = 0, yright = 1, ties = mean)
  
  coyi<-approx(gody$y,gody$x, n = ng, method = "linear",
               yleft = 0, yright = 1, ties = mean)
  plot(coxi$x,coxi$y)
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
  grid <- expand.grid(lon=godx, lat=gody)
  levelplot(acp3 ~ lon * lat, data=grid, cuts=20, pretty=T) 
  cl2<-contourLines(coxi,coyi, acp3, levels = pobje) 
  if(length(cl2)>0){
    cl2<-as.matrix(data.frame(cl2[[1]]$x,cl2[[1]]$y))} else{cl2<-NA}
  
}

JT_KDEap<-function(u2,pbas ,pobj,beta,margins,vtau,devplot=F,mar1,mar2,px,py,interh=NA){
  
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
  ModelHugo_file="C:/Users/k1638615/King's College London/OneDrive - King's College London/DATA/CODE ANALYSIS/Model comparison/LF/BveLTDep.R"
  source (ModelHugo_file)
  if(vtau<0)q0=0.9
  qq<-c()
  cm<-c()
  cq<-c()
  ccd<-c()
  qc<-.95
  rq0<-seq(0.75,0.95,by=0.01)
  for(q0 in rq0){
    estims<-try(BveLTDep (data= kk,mod.thresh.u = q0,crit.lev.u = qc,sig.lev=0.05,ci.meth='pl',marg.inf=T),silent = T)
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
  plot(ccd)
  plot(qq)
  sumd<-0
  for(i in 1:length(diff(qq))){
    sdiff<-diff(qq)[i]
    sumd<-c(sumd,sumd[i]+sdiff)}
  sh<-which(sumd<=-0.02|sumd>=0.02 )[1]
  q0<-rq0[sh-1]
  plot(sumd)
  estims<-try(BveLTDep (data= kk,mod.thresh.u = q0,crit.lev.u = qc,sig.lev=0.05,ci.meth='se',marg.inf=T),silent = T)
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
  
  if(interh=="casc"){
    pobj<-c(seq(0.00000005,0.00005,by=0.0000001),seq(0.00002,pbas/5,by=0.00005))
    pobj<-c(seq(0.0001,0.01,by=0.00001))
    s=pbas/pobj
    
    wqobjf<-c()
    for (sl in 1:length(s)){
      
      if (!is.na(Chilow) & (Chilow<0.05 & etahat<0.75| etahat<0.6)){
        print("AI")
        beta=200
        m1= 1- (wqxFrechet/(wqxFrechet+wqyFrechet))^beta
        m2<-1- (wqyFrechet/(wqxFrechet+wqyFrechet))^beta
        eta1<-m1*etahat + (1-m1)
        eta2<-m2*etahat + (1-m2)
        projx<-s[sl]^(eta1)*wqxFrechet
        projy<-s[sl]^(eta2)*wqyFrechet
      }else{ 
        # plot(xFrechet,yFrechet)
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
    
    
    plot(wqobjf[which(wqobjf$rep.pobj.sl...100.==0.001),1],wqobjf[which(wqobjf$rep.pobj.sl...100.==0.001),2],col=wqobjf[which(wqobjf$rep.pobj.sl...100.==0.001),4])
    points(u2, col=alpha("black",1),pch=16)
    
    tg=50
    
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
          matjt[k,j]=mean(coly,na.rm=T)/(1-pxg[k])}
      }
    }
    
    for (k in 1:(length(pxg)-1)){
      colx<-which(wqobjf[,2]>gridy[k] & wqobjf[,2]<=gridy[k+1])
      for (j in 1:(length(pyg)-1)){
        coly<-wqobjf[colx,3][which(wqobjf[colx,1]>gridx[j] & wqobjf[colx,1]<=gridx[j+1])]
        if(length(coly)==0){matjt[j,k]=matjt[j,k]}else{
          matjt[j,k]=mean(coly,na.rm=T)/(1-pxg[j])}
      }
    }
    
    
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
    
    library(zoo)
    
    for (nap in 1: length(pxg)){ matjt[,nap]<-na.approx(matjt[,nap],maxgap = 5,na.rm=F)}
    
    for (nap in 1: length(pxg)){ matjt[nap,]<-na.approx(matjt[nap,],maxgap = 5,na.rm=F)}
    
    levelplot(matjt ~ lon * lat, data=grid, cuts=20, pretty=T,contour=T) 
    
    contour(gridx,gridy,matjt,levels=0.001)
    points(jtres,pch=16)
    
    sh<-contourLines(gridx,gridy,matjt,levels=upobj)
    obx<-c()
    oby<-c()
    for (ssh in 1:length(sh)){
      obx<-c(obx,sh[[ssh]]$x)
      oby<-c(oby,sh[[ssh]]$y)
    }
    plot(obx,oby)
    
    
    
    wqobj<-data.frame(obx,oby)
  }
  
  res<-list(levelcurve=wqobj,etaJT=etahat,chiJT= Chimed,wq0ri=wq0ri)
}

Cond_modap<-function(u2,tr1,tr2,tsim,num.sim,pobj=0.001,interh="comb"){
  names(u2)= c("V1","V2")
  thresh1 <- tr1
  thresh2 <- tr2
  
ext.q=0.95
  
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
  names(sim1)<-names(u2)
  
  sim1<-sim1[which(sim1$V2<(sim1$V1+q2-q1)),]
  sim2<-sim2[which(sim2$V2>(sim2$V1+q2-q1)),]
  # mex.pred2$data$simulated=sim2
  # mex.pred$data$simulated=sim1
  onlysim<-rbind(sim1,sim2)
  u2<-as.data.frame(u2)
  names(onlysim)=names(u2)
  simu<-rbind(u2,onlysim)
  plot(onlysim)
  
  
  
  # JointExcureve_file=paste0(getwd(),"/Multimodel_Sim/Joint_excurves.R")
  # source (JointExcureve_file)
  pobj1=pobj/(1-tsim)
  
  if(interh=="casc"){
    simv=sim2#Can also condition on y being gretaer than a threshole
    ngx=10000
    
    
    godx<-approx(u1b,pxfp, n = ngx, method = "linear",
                 yleft = min(pxfp), yright = max(pxfp))
    cem.dens<-kcde(simv,gridsize=200, tail.flag = "upper.tail",xmin=c(min(simv[,1]),min(simv[,2])),xmax=c(max(godx$x),max(simv[,2])))

    
    jesus<-cem.dens$estimate
    lesx<-cem.dens$eval.points[[1]]
    lesy<-cem.dens$eval.points[[2]]


    lexp<-approx(godx$x, godx$y, xout = lesx, method = "linear",yleft = min(pxfp),yright = max(pxfp), rule = 1)$y
    ale<-jesus
    for (k in 1:length(lesy)){
      for (j in 1:length(lesx)){
        ale[j,k]=jesus[j,k]/(1-lexp[j])}}
    
    plot(u2,xlim=c(0,60),ylim=c(0,300))
    contour(lesx,lesy,jesus,levels = c(0.1),col=3,ylim=c(0,200),add=T)
    
    sh2<-contourLines(lesx,lesy,ale,levels=upobj/(1-tsim))
    obx<-c()
    oby<-c()
    for (ssh in 1:length(sh2)){
      obx<-c(obx,sh2[[ssh]]$x)
      oby<-c(oby,sh2[[ssh]]$y)
    }
    plot(obx,oby)
    
    
    jline<-data.frame(obx,oby)
  }
  
  
  if(interh=="comb"){
    j1 <- JointExceedanceCurve(mex.pred,pobj,n=100,meth="J") 
    j2 <- JointExceedanceCurve(mex.pred2,pobj,n=100,meth="J")
    summary(mex.pred2)
    #For confidence interval bootstrap
    # jb<-JointExceedanceCurve(mex.pred$replicates[[1]] ,pobj1,n=10,which=1) 
    # jbo<-data.frame(jb[[1]],jb[[2]])
    ##############################
    
    
    j1o<-data.frame(j1[[1]],j1[[2]])
    j2o<-data.frame(j2[[1]],j2[[2]])
    
    
    names(j1o)=names(j2o)
    # sim1<-sim1[which(sim1$V2<(sim1$V1+q2-q1)),]
    # sim2<-sim2[which(sim2$V2>(sim2$V1+q2-q1)),]
    jline<-rbind(j1o[which(j1o[,2]<j1o[,1]+q2-q1),],j2o[which(j2o[,2]>j2o[,1]+q2-q1),])
    jline<-jline[order(jline$j2..1..),]
      jj1<-as.matrix(j1o[which(j1o[,2]<j1o[,1]+q2-q1),])
      jj2<-as.matrix(j2o[which(j2o[,2]>j2o[,1]+q2-q1),])
      mj1<-min(jj1[,1],na.rm=T)
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
        jline<-jline[order(jline[,1]),]
      } 
      
      idm<-which(jline[,2]==min(jline[,2]))
      jline<-jline[-((idm+1):length(jline[,1])),]
    }
    plot(jline)
  
  res=list(jline=jline,chiHT=chiht,etaht=etaht,onlysim=onlysim)
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
