#THIS SCRIPT WILL BE FOR ANALYSING OUTPUTS



setwd("C:/Users/PhD Student/OneDrive - King's College London/DATA/CODE ANALYSIS/Model comparison")
setwd("C:/Users/k1638615/King's College London/OneDrive - King's College London/DATA/CODE ANALYSIS/Model comparison/Bivariate_Models/bivariate_models")
config_file=paste0(getwd(),"/config/config_general.R")
source(config_file)
source(paste0(getwd(),"/function/GenPlot.R"))
# source(config_file)
library(copula)
library(MASS)
library(animation)
chiv1<-c(0.05,0.1,0.5,0.7,0.9)
thx=1/log2(2-chiv1)
nbpt<-c(500,1000,2000,5000,10000)
theta=round(thx[1],2)
realchi=chiv1[1]
mars=c(0.5,0.5)
pc<-2
mama=4
print=F
td=3
nbpo=5000
if (pc==1){d="C:/Users/PhD Student/"}
if (pc==2){d="C:/Users/k1638615/King's College London/"}

generoutputs<-function(mars,mama,td,print=F,nbpo){
  
if(mama==2){locf<-"Nmar"}
if(mama==3){locf<-"Newrun"}
if(mama==4){locf<-"NewRunConf"}
if(mama==1){locf<-"Clean"}  
if(td==1)type="nb"
if(td==2)type="chi"
if(td==3)type="eta"
if(td==4)type="perso"
mods<-c("Cond-Ex","JT-KDE","Gumcop","Normalcop","FGMcop","GalambosCop")
cols <- c("Cond-Ex" = "purple", "JT-KDE" = "green", "Gumcop" = "red", "Normalcop" = "skyblue", "FGMcop" = "darkorange", "GalambosCop" = "darkred") 
#Analysis stats

nbpt<-c(500,1000,2000,5000,10000)  

 if (type=="nb"){
    directory=paste0(d,"OneDrive - King's College London/DATA/CODE ANALYSIS/Model comparison/out/",locf,"/nbpt")
  load(paste0(directory,"/stat_output_F1_chi-0.3_nbv_0.5_0.5.Rdata"))
  sep=c(1,7,13,19,25)
  labs=c("500","1000","2000","5000","10000")
  realchi=0.3
  br=c(1:5)
  euclidist<-t(graphout$eucl_dist_Norm_curve)
  RChi<-t(graphout$RMSE_Chi)
  Chiv<-t(graphout$Chival)
  BiasF<-t(graphout$Bias_Norm)
  
  model<-c(rep(NA,30))
  wup<-c(rep(NA,30))
  pointsN<-c(rep(NA,30))
  
  euclidist<-data.frame(cbind(euclidist,model,wup))
  for (rs in 0:5)
  {euclidist[sep+rs,4]<-mods[rs+1]}
  
  for (rs in 1:5)
  { rl<-sep[rs]
  euclidist[seq(rl,rl+5),5]<-rs
  euclidist[seq(rl,rl+5),6]<-nbpt[rs]
  names(euclidist)[6]="PointNames"
 
  RChi<-data.frame(cbind(RChi,model,wup))
  for (rs in 0:5)
  {RChi[sep+rs,4]<-mods[rs+1]}
  
  for (rs in 1:5)
  { rl<-sep[rs]
  RChi[seq(rl,rl+5),5]<-rs
  RChi[seq(rl,rl+5),6]<-nbpt[rs]
  }
  names(RChi)[6]="PointNames"
  
  Chiv<-data.frame(cbind(Chiv,model,wup))
  for (rs in 0:5)
  {Chiv[sep+rs,4]<-mods[rs+1]}
  
  for (rs in 1:5)
  { rl<-sep[rs]
  Chiv[seq(rl,rl+5),5]<-rs
  Chiv[seq(rl,rl+5),6]<-nbpt[rs]
  }
  names(Chiv)[6]="PointNames"
  

  
  BiasF<-data.frame(cbind(BiasF,model,wup))
  
  for (rs in 0:5)
  {BiasF[sep+rs,4]<-mods[rs+1]}
  
  for (rs in 1:5)
  { rl<-sep[rs]
  BiasF[seq(rl,rl+5),5]<-rs
  BiasF[seq(rl,rl+5),6]<-nbpt[rs]
  }
  names(BiasF)[6]="PointNames"
  
  }
  
 }else if(type=="chi")
   {
    directory=paste0(d,"OneDrive - King's College London/DATA/CODE ANALYSIS/Model comparison/out/",locf,"/chi/",nbpo)
    # load(paste0(directory,"/Stat_final-Chivar_nbtotal -2019-03-07.Rdata"))
    # load(paste0(directory,"/vals4chi0.3n5000.Rdata"))
    
    load(paste0(directory,"/stat_output_F1_chi-0.9_chiv_",mars[1],"_",mars[2],".Rdata"))
    chiv1<-c(0.05,0.1,0.3,0.5,0.7,0.9)
    sep=c(1,7,13,19,25,31)
    realchi=chiv1
    br=c(1:6)
  labs=c("0.05","0.1","0.3","0.5","0.7","0.9") 
    euclidist<-t(graphout$eucl_dist_Norm_curve)
    RChi<-t(graphout$RMSE_Chi)
    Chiv<-t(graphout$Chival)
    BiasF<-t(graphout$Bias_Norm)
    #change 30 to 36
    model<-c(rep(NA,36))
    wup<-c(rep(NA,36))
    chivalN<-c(rep(NA,36))
    
    euclidist<-data.frame(cbind(euclidist,model,wup,chivalN))
    for (rs in 0:5)
    {euclidist[sep+rs,4]<-mods[rs+1]}
    # names(addchi03[[1]])[5:6]=c("wup","chivalN")
    # euclidist<-data.frame(rbind(euclidist[1:12,],addchi03$eulci,euclidist[13:30,]))  
    for (rs in 1:6)
    { rl<-sep[rs]
    euclidist[seq(rl,rl+5),5]<-rs
    euclidist[seq(rl,rl+5),6]<-chiv1[rs]
    }
    names(euclidist)[6]="PointNames"
    
    RChi<-data.frame(cbind(RChi,model,wup,chivalN))
    for (rs in 0:5)
    {RChi[sep+rs,4]<-mods[rs+1]}
    # names(addchi03[[2]])[5:6]=c("wup","chivalN")
    # RChi<-data.frame(rbind(RChi[1:12,],addchi03$rmsechi,RChi[13:30,]))    
    for (rs in 1:6)
    { rl<-sep[rs]
    RChi[seq(rl,rl+5),5]<-rs
    RChi[seq(rl,rl+5),6]<-chiv1[rs]
    }
    names(RChi)[6]="PointNames"
    
    Chiv<-data.frame(cbind(Chiv,model,wup,chivalN))
    for (rs in 0:5)
    {Chiv[sep+rs,4]<-mods[rs+1]}
    # names(addchi03[[3]])[5:6]=c("wup","chivalN")
    # Chiv<-data.frame(rbind(Chiv[1:12,],addchi03$chi,Chiv[13:30,]))     
    for (rs in 1:6)
    { rl<-sep[rs]
    Chiv[seq(rl,rl+5),5]<-rs
    Chiv[seq(rl,rl+5),6]<-chiv1[rs]
    }
    names(Chiv)[6]="PointNames" 
    

    
    BiasF<-data.frame(cbind(BiasF,model,wup, chivalN))
    
    for (rs in 0:5)
    {BiasF[sep+rs,4]<-mods[rs+1]}
    # names(addchi03[[4]])[5:6]=c("wup","chivalN")
    # BiasF<-data.frame(rbind(BiasF[1:12,],addchi03$bias,BiasF[13:30,]))     
    for (rs in 1:6)
    { rl<-sep[rs]
    BiasF[seq(rl,rl+5),5]<-rs
    BiasF[seq(rl,rl+5),6]<-chiv1[rs]
    }
    names(BiasF)[6]="PointNames"
    pd <- position_dodge(.2)
  }else if(type=="eta")
    {
    directory=paste0(d,"OneDrive - King's College London/DATA/CODE ANALYSIS/Model comparison/out/",locf,"/eta/",nbpo)
    # load(paste0(directory,"/Stat_final-Chivar_nbtotal -2019-03-07.Rdata"))
    # load(paste0(directory,"/vals4chi0.3n5000.Rdata"))
    
    load(paste0(directory,"/stat_output_F1_chi-0.9_chiv_",mars[1],"_",mars[2],".Rdata"))
    Etav1<-c(0.25,0.5,0.75,0.9)
    realeta=Etav1
    labs=c("0.25","0.5","0.75","0.9")
    br=c(1:4)

    sep=c(1,7,13,19)

    euclidist<-t(graphout[[2]])
    Etav<-t(graphout[[3]])
    BiasF<-t(graphout[[6]])
    
    model<-c(rep(NA,24))
    wup<-c(rep(NA,24))
    EtavalN<-c(rep(NA,24))
    euclidist<-data.frame(cbind(euclidist,model,wup,EtavalN))
    for (rs in 0:5)
    {euclidist[sep+rs,4]<-mods[rs+1]}
    
    for (rs in 1:4)
    { rl<-sep[rs]
    euclidist[seq(rl,rl+5),5]<-rs
    euclidist[seq(rl,rl+5),6]<-Etav1[rs]
    
    }
    names(euclidist)[6]="PointNames"
    Etav<-data.frame(cbind(Etav,model,wup,EtavalN))
    
    for (rs in 0:5)
    {Etav[sep+rs,4]<-mods[rs+1]}
    
    for (rs in 1:4)
    { rl<-sep[rs]
    Etav[seq(rl,rl+5),5]<-rs
    Etav[seq(rl,rl+5),6]<-Etav1[rs]
    }
    names(Etav)[6]="PointNames"
    
    BiasF<-data.frame(cbind(BiasF,model,wup,EtavalN))
    
    for (rs in 0:5)
    {BiasF[sep+rs,4]<-mods[rs+1]}
    
    for (rs in 1:4)
    { rl<-sep[rs]
    BiasF[seq(rl,rl+5),5]<-rs
    BiasF[seq(rl,rl+5),6]<-Etav1[rs]
    }
    names(BiasF)[6]="PointNames"
    
  }

  
    
  


  #   #Sacred file
  #   load(paste0(getwd(),"/out/Clean/Stat_final-Chivar_nbtotal -2019-03-07.Rdata"))
  # }else if(type=="eta"){
  # load(paste0(directory,"Stat_final-qtyPoints_nbtotal -2019-03-06.Rdata")) 
  # }


# for (rs in 0:5)
# {ecdistF<-c(ecdistF,list(graphout[[1]][,sep+rs]))}
# 
# brx<-unique(euclidist$PointNames)
# sh<-(unique(euclidist$model))
# fp<-euclidist[which(euclidist$wup==5),]
# newvar<-c()
# newvarm<-c()
# newvarM<-c()
# modix<-c()
# for (st in 1:6){
# newvar<-c(newvar,euclidist[which(euclidist$model==fp$model[st]),(1)] -  fp[which(fp$model==fp$model[st]),(1)])
# newvarm<-c(newvarm,euclidist[which(euclidist$model==fp$model[st]),(2)] -  fp[which(fp$model==fp$model[st]),(2)])
# newvarM<-c(newvarM,euclidist[which(euclidist$model==fp$model[st]),(3)] -  fp[which(fp$model==fp$model[st]),(3)])
# modix<-c(modix,rep(fp$model[st],6))
# }
# merde<-data.frame(modix,newvar,newvarm,newvarM)
# merde<-merde[order(merde$modix),]
# euclidist<-euclidist[order(euclidist$model),]
# euclidist$difmed<-merde$newvar
# euclidist$difm<-merde$newvarm
# euclidist$difM<-merde$newvarM
# plot(euclidist$difmed,col=cols)
# 
# euclidistprout<-euclidist[which(euclidist$model=="Gumcop"),]
# euclidist$ci<-euclidist$Q97.5-euclidist$Q2.5
# pd <- position_dodge(.100) # move them .05 to the left and right
# 
# ggplot(euclidist, aes(x=model, y=median, colour=PointNames)) +
#   geom_point() +
#   xlab("sd of margins")+
#   ylab("Eucliddean distance to the reference curve")   +    
#   ggtitle(paste0("Euclidean distance level curves - p = 0.001 - n = ",nbpo)) +
#   theme_bw() 


p1<-ggplot(euclidist, aes(y=median, x=PointNames, colour=model)) +
   # geom_errorbar(aes(ymin=Q2.5, ymax=Q97.5), width=1000,size=1, alpha=0.3,position=pd) +
  geom_line(position=pd,size=1) +
  geom_point(position=pd, size=3,shape=21, aes(fill=model))+
xlab("# of realisations") +
  ylab("Eucliddean distance to the reference curve") +
  scale_colour_manual(values = cols,breaks=mods,labels=mods)   +    
  # scale_x_continuous(breaks=brx,labels=labs)+
  ggtitle(paste0("Euclidean distance level curves - p = 0.001 - n = ",nbpo)) +
  # expand_limits(y=1) +                        # Expand y range
  # scale_y_continuous(limits=c(0,1.2*max(euclidist$Q97.5))) +         # Set tick every 4
  theme_bw() +scale_y_continuous(breaks=seq(0.001,0.3,by=0.01))
p1



ggplot(euclidist, aes(x=PointNames, y=Q97.5, fill=model)) + 
  geom_area()+
  geom_area(alpha=0.6 , size=1, colour="black")

ggplot(euclidist, aes(x=PointNames, y=Q2.5, fill=model)) + 
  geom_area()+
  geom_area(alpha=0.6 , size=1, colour="black")


# 
# nRes<-nbpt[-1]
# specgum<-euclidist[which(RChi$model=="Gumcop"),]
# spcchig<-euclidist[which(RChi$model=="Gumcop"),]
# # seng<--diff(specgum$median)/diff(specgum$PointNames)*1000
# seng=specgum
# par(mar = c(5,5,3,1))
# plot(seng$PointNames,seng$median,type="n",lwd=2,ylim=c(0,0.3), xaxt="n",xlab= " n (# of realisations)"
#      ,ylab=expression(paste("-",delta,"(RMSE"[chi],")/",delta,"n x 1000")),cex.lab=1.5,cex.axis=1.4)
# axis(1, at=nRes,labels=nRes, las=1,cex.axis=1.4)
# # legend(3.3,0.035,legend=c("Cond_Ex", "JT-KDE","GumbelCop","NormalCop","FGMCop","GalambosCop"), bg=NULL,
#        # col=cols, cex=.9, lwd=2)
# spec<-list()
# for (ik in (1:6)){
# spec[[ik]]<-euclidist[which(euclidist$model==mods[ik]),]
# seng<-spec[[ik]]
# # seng<--diff(spec[[ik]]$median)/diff(spec[[ik]]$PointNames)*1000
# # sengm<--diff(spec[[ik]]$Q2.5)/diff(spec[[ik]]$PointNames)*1000
# # sengM<--diff(spec[[ik]]$Q97.5)/diff(spec[[ik]]$PointNames)*1000
# points(seng$PointNames,seng$median,type="b",lwd=2,col=cols[ik]) 
# # points (sengm,type="b",lwd=2,lty=2,col=cols[ik])
# # points (sengM,type="b",lwd=2,lty=2,col=cols[ik])
# }
# #RMSE Chi
if (type=="nb"|type=="chi"|type=="perso"){
# 
pd <- position_dodge(.2) # move them .05 to the left and right
RChi2<-RChi[,c(1,2,3,5,6)]
plot(RChi2[,5],RChi2[,1])

p2<-ggplot(RChi, aes(x=PointNames, y=median, colour=model)) +
  geom_errorbar(aes(ymin=Q2.5, ymax=Q97.5), width=.3,size=1, position=pd) +
  geom_line(position=pd,size=1) +
  geom_point(position=pd, size=3,shape=21, fill="white")+
  xlab("number of observations") +
  ylab("RMSE Chi value") +
 scale_colour_manual(values = cols,breaks=mods,labels=mods)  +
  ggtitle("RMSE Chi value estimated by the models -  Chi = 0.3") +
  expand_limits(y=0.05) +                        # Expand y range        # Set tick every 4
  theme_bw() 

p2
# scale_colour_hue(name="Models",    # Legend label, use darker colors
# breaks=mods,
# l=50) 

}
#Chi Value
if (type=="nb"|type=="chi"|type=="perso"){

pd <- position_dodge(.2) # move them .05 to the left and right

p3<-ggplot(Chiv, aes(x=wup, y=median, colour=model)) +
  geom_errorbar(aes(ymin=Q2.5, ymax=Q97.5), width=.7,size=0.7, position=pd) +
  geom_line(position=pd,size=.7) +
  geom_point(position=pd, size=3,shape=21, fill="white")+
  xlab("number of observations") +
  ylab("Chi value") +
  scale_colour_manual(values = cols,breaks=mods,labels=mods)  +    
  scale_x_continuous(breaks=br,labels=labs)+
  ggtitle("Chi values estimated by the models") +
  expand_limits(y=0.05) +                        # Expand y range        # Set tick every 4
  theme_bw() 
p3
}

if (type=="eta"){
  
  pd <- position_dodge(.2) # move them .05 to the left and right
  
  p2<-ggplot(Etav, aes(x=wup, y=median, colour=model)) +
    geom_errorbar(aes(ymin=Q2.5, ymax=Q97.5), width=.7,size=0.7, position=pd) +
    geom_line(position=pd,size=.7) +
    geom_point(position=pd, size=3,shape=21, fill="white")+
    xlab("number of observations") +
    ylab("Eta value") +
    scale_colour_manual(values = cols,breaks=mods,labels=mods)  +    
    scale_x_continuous(breaks=br,labels=labs)+
    ggtitle("Eta values estimated by the models") +
    expand_limits(y=0.05) +                        # Expand y range        # Set tick every 4
    theme_bw() 
}
#Bias
p2

pd <- position_dodge(.02)
# move them .05 to the left and right

p4<-ggplot(BiasF, aes(x=PointNames, y=median, colour=model)) +
  geom_errorbar(aes(ymin=Q2.5, ymax=Q97.5), width=.1,size=0.7, position=pd) +
  geom_line(position=pd,size=1) +
  geom_point(position=pd, size=3,shape=21, fill="white")+
  geom_hline(yintercept=0, col="darkgrey", size=1,linetype="dashed")+
  xlab("value of Chi measure") +
  ylab("Bias") +
  scale_colour_manual(values = cols,breaks=mods,labels=mods)  +   
  scale_x_continuous(breaks=br,labels=labs)+
  ggtitle("Biases level curves - p = 0.001 - Chi = 0.3") +
  expand_limits(y=0.05) +                        # Expand y range        # Set tick every 4
  theme_bw()
p4

if(print==T)pdf(paste0(d,"OneDrive - King's College London/DATA/CODE ANALYSIS/Model comparison/out/figures/Clean/results_comp_study_",type,"_",mars[1],"_",mars[2],"_",nbpo,".pdf"),29 / 3.5,25 / 3.5)
if (type=="nb"|type=="chi"){
  print(p1)
  print(p2)
  print(p3)
  print(p4)
}
if (type=="eta"){
  print(p1)
  print(p2)
  print(p4)
  # ggsave(filename=paste0("results_comp_study_",type,"_01.pdf"), plot = p1, device = pdf, path = "C:/Users/PhD Student/OneDrive - King's College London/DATA/CODE ANALYSIS/Model comparison/out/figures/Clean/", scale = 1, width = 29, height = 25, units = "cm", dpi = 300, limitsize = TRUE)
}

graphics.off()

# addchi03<-list()
# addchi03<-c(addchi03,list(euclidist[which(euclidist$PointNames==5000),]))
# addchi03<-c(addchi03,list(RChi[which(RChi$PointNames==5000),]))
# addchi03<-c(addchi03,list(Chiv[which(Chiv$PointNames==5000),]))
# addchi03<-c(addchi03,list(BiasF[which(BiasF$PointNames==5000),]))
# names(addchi03)=c("eulci","rmsechi","chi","bias")


if(print==T)pdf(paste0(d,"OneDrive - King's College London/DATA/CODE ANALYSIS/Model comparison/out/figures/Clean/Levelcurves_IC_comp_study_",type,"_",mars[1],"_",mars[2],"_",nbpo,".pdf"),29 / 3.5,25 / 3.5)

if (type=="nb"|type=="perso"){
  ch<-c(500,1000,2000,5000,10000)}
if(type=="chi"){
  ch<-c(0.05,0.1,0.3,0.5,0.7,0.9)
  thx=1/log2(2-ch) }
if(type=="eta"){
  ch<-c(0.25,0.5,0.75,0.9)
  }

for (is in 1:length(ch)){ 

if(type=="nb"){
  chireal=0.3
  theta=1/log2(2-chireal)
  real=ch[is]
  num.points=real
  load(paste0(directory,"/","nb_",real,"_sd_",mars[1],"_",mars[2],"/lines_chi",chireal,"_nb",real,".Rdata"))
  load(paste0(directory,"/","nb_",real,"_sd_",mars[1],"_",mars[2],"/linesIC95_chi",chireal,"_nb",real,".Rdata"))
  dat<-gendist(AIndep = F, coef=chireal,nbpt=real, margin="lnorm",pm1 =c(2,mars[1]) ,pm2=c(2,mars[2]))
}else if(type=="perso"){
  chireal=0.3
  theta=1/log2(2-chireal)
  real=ch[is]
  load(paste0(directory,"/",real,"/lines_chi",round(chireal,2),"_nb",real,".Rdata"))
  load(paste0(directory,"/",real,"/lines_IC95",round(chireal,2),"_nb",real,".Rdata"))
  dat<-gendist(AIndep = F, coef=chireal,nbpt=real, margin="lnorm",pm1 =c(2,mars[1]) ,pm2=c(2,mars[1]))
}else if (type=="chi"){
    num.points=nbpo
    theta=thx[is]
    real=ch[is]
    load(paste0(directory,"/chi_",real,"_sd_",mars[1],"_",mars[2],"/lines_dep",real,"_nb",num.points,".Rdata"))
    load(paste0(directory,"/chi_",real,"_sd_",mars[1],"_",mars[2],"/linesIC95_dep",real,"_nb",num.points,".Rdata"))
    dat<-gendist(AIndep = F, coef=real,nbpt=5000, margin="lnorm",pm1 =c(2,mars[1]) ,pm2=c(2,mars[2]),ctype="comb")
  }else if(type=="eta"){
    real=ch[is]
    num.points=nbpo
load(paste0(directory,"/eta ",real,"_sd_",mars[1],"_",mars[2],"/lines_dep",real,"_nb",num.points,".Rdata"))
load(paste0(directory,"/eta ",real,"_sd_",mars[1],"_",mars[2],"/linesIC95_dep",real,"_nb",num.points,".Rdata"))
dat<-gendist(AIndep = T, coef=real,nbpt=5000, margin="lnorm",pm1 =c(2,mars[1]) ,pm2=c(2,mars[2]))
  }

# if (type=="perso"){ 
# plot(dat$u2, col="grey",xlim=c(0,max(dat$cont$x)),ylim=c(0,max(dat$cont$y)),main=paste0("Joint exceedance curves for p = 0.001 \n ",type, " = ",round(real,3), " n = ",num.points," with 95% CI"), xlab="Var 1", ylab="Var 2") 
# }else{
  adpd<-contourLines(dat$cont$x,dat$cont$y, dat$ax, levels = c(0.001))[[1]]
   plot(dat$u2, col="grey",xlim=c(0,1.5*max(adpd$x)),ylim=c(0,1.5*max(adpd$y)),main=paste0("Joint exceedance curves for p = 0.001 \n ",type, " = ",round(real,3), " n = ",num.points," with 95% CI"), xlab="Var 1", ylab="Var 2")
   
   # plot(dat$cont$x,dat$cont$y, dat$ax, levels = c(0.1,0.001,0.00001,0.000001),add =F,lwd=8 ,drawlabels=F,col=alpha(c("grey","gold","darkorange","darkred"),0.6),xaxt="n",yaxt="n")
   # grid<-expand.grid(lon=dat$cont$x, lat=dat$cont$y)
   # library(viridis)
   # u<-dat$u2
   # cols <- brewer.pal(4, "YlOrBr")
   # 
   # afu<-kcde(u,gridsize=1000, xmin=c(min(u[,1]),min(u[,2])), xmax=c(15,15),tail.flag = "upper.tail")
   # plot(afu,display="image" ,cont=c(0.001,0.05,0.2,0.5),col=cols,xlim=c(2,10),ylim=c(2,10))
   # 

grid (NULL,NULL, lty = 1, col = "lightgray")
legend(1.1*max(adpd$x),1.1*max(adpd$y),legend=c("Cond_Ex", "JT-KDE","GumbelCop","NormalCop","FGMCop","GalambosCop"), bg=NULL,
       col=cols, cex=.9, lwd=2)
contour(dat$cont$x,dat$cont$y, dat$ax, levels = c(0.1,0.05,0.01,0.001), add=T, col="black",delta=0.1,xlim=c(0,100),ylim=c(0,100),box=F, lwd=4,n.grid=200) 

cols2<-c("purple","green","red")
rbPal <- colorRampPalette(c("green","orange","red"))
# for(tt in 1:100){
#   lines(linesfinal$x[[1]][,tt],linesfinal$y[[1]][,tt],type="p" ,lwd=1,lty=1, col=1)
#   lines(linesfinal$x[[2]][,tt],linesfinal$y[[2]][,tt],type="p" ,lwd=1,lty=1, col=2)
#   lines(linesfinal$x[[3]][,tt],linesfinal$y[[3]][,tt],type="p" ,lwd=1,lty=1, col=3)
#   lines(linesfinal$x[[4]][,tt],linesfinal$y[[4]][,tt],type="p" ,lwd=1,lty=1, col=4)
#   lines(linesfinal$x[[5]][,tt],linesfinal$y[[5]][,tt],type="p" ,lwd=1,lty=1, col=5)
# }
# for (zg in 1:6){
# 
#   colo <- rbPal(10)[as.numeric(cut(linesIC$d[[zg]][,1],breaks = 10))]
#   points(linesIC$x[[zg]][,1],linesIC$y[[zg]][,1],type="p" ,lwd=3,lty=1, col=colo)
# 
# }
# 
# lic<-data.frame(linesIC$x[[4]][,1],linesIC$y[[4]][,1],linesIC$d[[4]][,1])
# names(lic)=c("x","y","d")
# 
# lic2<-data.frame(linesIC$x[[1]][,1],linesIC$y[[1]][,1],linesIC$d[[1]][,1])
# names(lic2)=c("x","y","d")
# cumd=lic$d
# for(n in 2:100){
#   cumd[n]= sum(cumd[n-1],lic$d[n])
# }
# 
# lim1<-lic[max(which(cumd<0.025)),]
# lim2<-lic[min(which(cumd>0.975)),]
# 
# idlim<-rbind(lim1,lim2)
# 
# pq<-ggplot() +
#   geom_point(data=dat$u2,aes(x=V1,y=V2,alpha=0.5),shape=21)
# 
# pq+geom_line(data=lic,aes(x=x,y=y,colour=d,size=2))+
# scale_color_gradient2(midpoint=mean(lic$d), low="white", mid="green",
#                         high="red", space ="Lab" )+
#   geom_point(data=lic[which(lic$d==max(lic$d)),],aes(x=x,y=y),shape=19,colour="black",size=4)+
#   geom_point(data=idlim,aes(x=x,y=y),shape=18,size=4)+
#   theme_bw()+
# theme(legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=14),axis.title=element_text(size=15))+
# 
#   labs(x =  expression('X'[1]), y =expression('X'[2]),size=10)
# 
# plot(lic$d)
# 
# 
# points(idlim$d,col=2)

for (zg in 1:6){
   polygon(c(linesIC$x[[zg]][,2],rev(linesIC$x[[zg]][,3])),c(linesIC$y[[zg]][,2],rev(linesIC$y[[zg]][,3])),col = alpha(cols[zg], 0.2), border = FALSE)
  colo <- rbPal(10)[as.numeric(cut(linesIC$d[[zg]][,1],breaks = 10))]
  points(linesIC$x[[zg]][,1],linesIC$y[[zg]][,1],type="l" ,lwd=3,lty=1, col=cols[zg])
  lines(linesIC$x[[zg]][,3],linesIC$y[[zg]][,3],type="l", lwd=3,lty=2,col=cols[zg])
  lines(linesIC$x[[zg]][,2],linesIC$y[[zg]][,2],type="l" ,lwd=3,lty=2,col=cols[zg])
}
}
dev.off()

# saveGIF(for (is in 1:length(ch)){ 
#   
#   if(type=="nb"){
#     chireal=0.3
#     theta=1/log2(2-chireal)
#     real=ch[is]
#     load(paste0(directory,"/","nb_",real,"_sd_",mars[1],"_",mars[2],"/lines_chi",chireal,"_nb",real,".Rdata"))
#     load(paste0(directory,"/","nb_",real,"_sd_",mars[1],"_",mars[2],"/linesIC95_chi",chireal,"_nb",real,".Rdata"))
#     dat<-gendist(AIndep = F, coef=chireal,nbpt=real, margin="lnorm",pm1 =c(2,mars[1]) ,pm2=c(2,mars[2]))
#   }else if (type=="chi"){
#     num.points=nbpo
#     theta=thx[is]
#     real=ch[is]
#     load(paste0(directory,"/chi_",real,"_sd_",mars[1],"_",mars[2],"/lines_chi",real,"_nb",num.points,".Rdata"))
#     load(paste0(directory,"/chi_",real,"_sd_",mars[1],"_",mars[2],"/linesIC95_chi",real,"_nb",num.points,".Rdata"))
#     dat<-gendist(AIndep = F, coef=real,nbpt=5000, margin="lnorm",pm1 =c(2,mars[1]) ,pm2=c(2,mars[2]))
#   }else if(type=="eta"){
#     num.points=nbpo
#     real=ch[is]
#     load(paste0(directory,"/eta ",real,"_sd_",mars[1],"_",mars[2],"/lines_eta",real,"_nb",num.points,".Rdata"))
#     load(paste0(directory,"/eta ",real,"_sd_",mars[1],"_",mars[2],"/linesIC95_eta",real,"_nb",num.points,".Rdata"))
#     dat<-gendist(AIndep = T, coef=real,nbpt=5000, margin="lnorm",pm1 =c(2,mars[1]) ,pm2=c(2,mars[2]))
#   }
#   adpd<-contourLines(dat$cont$x,dat$cont$y, dat$ax, levels = c(0.001))[[1]]
#   plot(dat$u2, col="grey",xlim=c(0,1.5*max(adpd$x)),ylim=c(0,1.5*max(adpd$y)),main=paste0("Joint exceedance curves for p = 0.001 \n ",type, " = ",round(real,3), " n = ",num.points," with 95% CI"), xlab="Var 1", ylab="Var 2")
# 
#   
#   grid (NULL,NULL, lty = 1, col = "lightgray")
#   legend(1.1*max(adpd$x),1.1*max(adpd$y),legend=c("Cond_Ex", "JT-KDE","GumbelCop","NormalCop","FGMCop","GalambosCop"), bg=NULL,
#          col=cols, cex=.9, lwd=2)
#   contour(dat$cont$x,dat$cont$y, dat$ax, levels = c(0.1,0.05,0.01,0.001), add=T, col="black",delta=0.1,xlim=c(0,100),ylim=c(0,100),box=F, lwd=4) 
#   
#   cols2<-c("purple","green","red")
#   rbPal <- colorRampPalette(c("green","orange","red"))
#   
#   for (zg in 1:6){
#     polygon(c(linesIC$x[[zg]][,2],rev(linesIC$x[[zg]][,3])),c(linesIC$y[[zg]][,2],rev(linesIC$y[[zg]][,3])),col = alpha(cols[zg], 0.4), border = FALSE)
#     lines(linesIC$x[[zg]][,1],linesIC$y[[zg]][,1],type="l" ,lwd=3,lty=1, col=cols[zg])
#     lines(linesIC$x[[zg]][,3],linesIC$y[[zg]][,3],type="l", lwd=3,lty=2,col=cols[zg])
#     lines(linesIC$x[[zg]][,2],linesIC$y[[zg]][,2],type="l" ,lwd=3,lty=2,col=cols[zg])
#   }
# }, movie.name = paste0(d,"OneDrive - King's College London/DATA/CODE ANALYSIS/Model comparison/out/figures/Clean/animation",type,"_",mars[1],"_",mars[2],"_",nbpo,".gif"), img.name = "Rplot", convert = "magick",cmd.fun, clean = TRUE, extra.opts = "")

}


sdl1<-c(0.25,0.5,0.5,0.25,0.5,1.5)
sdl2<-c(0.25,0.25,0.5,1.5,1.5,1.5)


for (m in 1:6){
mars=c(sdl1[m],sdl2[m])
generoutputs(mama=4,mars=mars,td=2,print=F,nbpo=5000)
}


graphoutT<-list()
if(mama==2){locf<-"Nmar"}
if(mama==1){locf<-"Clean"} 
if(mama==3){locf<-"Newrun"}
if(mama==4){locf<-"NewRunConf"}
if(td==1)type="nb"
if(td==2)type="chi"
if(td==3)type="eta"
if(td==4)type="perso"
mods<-c("Cond-Ex","JT-KDE","Gumcop","Normalcop","FGMcop","GalambosCop")
cols <- c("Cond-Ex" = "purple", "JT-KDE" = "green", "Gumcop" = "red", "Normalcop" = "orange", "FGMcop" = "darkorange", "GalambosCop" = "darkred") 







#Analysis stats
for (m in 1:6){
  mars=c(sdl1[m],sdl2[m])
  if(type=="chi"){
    directory=paste0(d,"OneDrive - King's College London/DATA/CODE ANALYSIS/Model comparison/out/",locf,"/chi/",nbpo)
  # load(paste0(directory,"/Stat_final-Chivar_nbtotal -2019-03-07.Rdata"))
  # load(paste0(directory,"/vals4chi0.3n5000.Rdata"))
  
  load(paste0(directory,"/stat_output_F1_chi-0.9_chiv_",mars[1],"_",mars[2],".Rdata"))
  }
  if(type=="eta"){
  directory=paste0(d,"OneDrive - King's College London/DATA/CODE ANALYSIS/Model comparison/out/",locf,"/eta/",nbpo)
  # load(paste0(directory,"/Stat_final-Chivar_nbtotal -2019-03-07.Rdata"))
  # load(paste0(directory,"/vals4chi0.3n5000.Rdata"))
  
  load(paste0(directory,"/stat_output_F1_chi-0.9_chiv_",mars[1],"_",mars[2],".Rdata"))
  }
graphoutT[[m]]<-graphout

}

iter<-c("0.25_0.25","0.25_0.5","0.5_0.5","0.25_1.5","0.5_1.5","1.5_1.5")
euclidisti<-as.data.frame(matrix(0, ncol = 6, nrow = 0))
RMSEi<-as.data.frame(matrix(0, ncol = 6, nrow = 0))
weighecli<-as.data.frame(matrix(0, ncol = 6, nrow = 0))

mods<-c("Cond-Ex","JT-KDE","Gumcop","Normalcop","FGMcop","GalambosCop")
cols <- c("Cond-Ex" = "purple", "JT-KDE" = "green", "Gumcop" = "red", "Normalcop" = "orange", "FGMcop" = "darkorange", "GalambosCop" = "darkred") 
#Analysis stats

nbpt<-c(500,1000,2000,5000,10000)  
if (type=="chi"){
for (m in 1:6){
  
chiv1<-c(0.05,0.1,0.3,0.5,0.7,0.9)
sep=c(1,7,13,19,25,31)
iti<-iter[m]
realchi=chiv1
br=c(1:6)
graphout<-graphoutT[[m]]
labs=c("0.05","0.1","0.3","0.5","0.7","0.9") 
euclidist<-t(graphout[[2]])
Weucli<-t(graphout$eucl_dist_curve_WeighNorm)
RMSE<-t(graphout$RMSE_Chi)

# BiasF<-t(graphout[[6]])
#change 30 to 36
model<-c(rep(NA,36))
wup<-c(rep(NA,36))
chivalN<-c(rep(NA,36))
iteri<-c(rep(NA,36))
euclidist<-cbind(euclidist,model,wup,chivalN,iteri)
euclidist<-as.data.frame(euclidist)
for (rs in 0:5)
{euclidist[sep+rs,4]<-mods[rs+1]}
# names(addchi03[[1]])[5:6]=c("wup","chivalN")
# euclidist<-data.frame(rbind(euclidist[1:12,],addchi03$eulci,euclidist[13:30,]))  
for (rs in 1:6)
{ rl<-sep[rs]
euclidist[seq(rl,rl+5),5]<-rs
euclidist[seq(rl,rl+5),6]<-chiv1[rs]
}
euclidist[,7]<-iti
names(euclidist)[6]="Chinames"

euclidisti=rbind(euclidisti,euclidist)


Weucli<-data.frame(cbind(Weucli,model,wup,chivalN,iteri))
for (rs in 0:5)
{Weucli[sep+rs,4]<-mods[rs+1]}
# names(addchi03[[2]])[5:6]=c("wup","chivalN")
# Weucli<-data.frame(rbind(Weucli[1:12,],addchi03$rmsechi,Weucli[13:30,]))    
for (rs in 1:6)
{ rl<-sep[rs]
Weucli[seq(rl,rl+5),5]<-rs
Weucli[seq(rl,rl+5),6]<-chiv1[rs]
}
Weucli[,7]<-iti
names(Weucli)[6]="PointNames"
weighecli=rbind(weighecli,Weucli)

RMSE<-data.frame(cbind(RMSE,model,wup,chivalN,iteri))
for (rs in 0:5)
{RMSE[sep+rs,4]<-mods[rs+1]}
# names(addchi03[[3]])[5:6]=c("wup","RMSEalN")
# RMSE<-data.frame(rbind(RMSE[1:12,],addchi03$chi,RMSE[13:30,]))     
for (rs in 1:6)
{ rl<-sep[rs]
RMSE[seq(rl,rl+5),5]<-rs
RMSE[seq(rl,rl+5),6]<-chiv1[rs]
}
RMSE[,7]<-iti
names(RMSE)[6]="chivals" 
RMSEi=rbind(RMSEi,RMSE)
# 
# 
# BiasF<-data.frame(cbind(BiasF,model,wup, chivalN))
# 
# for (rs in 0:5)
# {BiasF[sep+rs,4]<-mods[rs+1]}
# # names(addchi03[[4]])[5:6]=c("wup","chivalN")
# # BiasF<-data.frame(rbind(BiasF[1:12,],addchi03$bias,BiasF[13:30,]))     
# for (rs in 1:6)
# { rl<-sep[rs]
# BiasF[seq(rl,rl+5),5]<-rs
# BiasF[seq(rl,rl+5),6]<-chiv1[rs]
# }
# names(BiasF)[6]="chivals"
# pd <- position_dodge(.2)
 }
}
euclidistiCh<-euclidisti
RMSEich<-RMSEi

if(type=="eta"){
for (m in 1:6){
  
  etav1<-c(0.25,0.5,0.75,0.9)
  sep=c(1,7,13,19)
  iti<-iter[m]
  realchi=etav1
  br=c(1:4)
  graphout<-graphoutT[[m]]
  labs=c("0.25","0.5","0.75","0.9") 
  euclidist<-t(graphout[[2]])
  Weucli<-t(graphout$eucl_dist_curve_WeighNorm)
  RMSE<-t(graphout$RMSE_Eta) 
  Etav<-t(graphout$Etaval)
  BiasF<-t(graphout$Bias_Norm)
  #change 30 to 36
  model<-c(rep(NA,24))
  wup<-c(rep(NA,24))
  etavalN<-c(rep(NA,24))
  iteri<-c(rep(NA,24))
  euclidist<-cbind(euclidist,model,wup,etavalN,iteri)
  euclidist<-as.data.frame(euclidist)
  for (rs in 0:5)
  {euclidist[sep+rs,4]<-mods[rs+1]}
  
  for (rs in 1:4)
  { rl<-sep[rs]
  euclidist[seq(rl,rl+5),5]<-rs
  euclidist[seq(rl,rl+5),6]<-etav1[rs]
  
  }
  euclidist[,7]<-iti
  names(euclidist)[6]="PointNames"
  
  euclidisti=rbind(euclidisti,euclidist)

  
  
  Weucli<-data.frame(cbind(Weucli,model,wup,etavalN,iteri))
  for (rs in 0:5)
  {Weucli[sep+rs,4]<-mods[rs+1]}
  # names(addchi03[[2]])[5:6]=c("wup","chivalN")
  # Weucli<-data.frame(rbind(Weucli[1:12,],addchi03$rmsechi,Weucli[13:30,]))    
  for (rs in 1:4)
  { rl<-sep[rs]
  Weucli[seq(rl,rl+5),5]<-rs
  Weucli[seq(rl,rl+5),6]<-etav1[rs]
  }
  Weucli[,7]<-iti
  names(Weucli)[6]="PointNames"
  weighecli=rbind(weighecli,Weucli)
  
  RMSE<-data.frame(cbind(RMSE,model,wup,etavalN,iteri))
  for (rs in 0:5)
  {RMSE[sep+rs,4]<-mods[rs+1]}
  # names(addchi03[[3]])[5:6]=c("wup","RMSEalN")
  # RMSE<-data.frame(rbind(RMSE[1:12,],addchi03$chi,RMSE[13:30,]))     
  for (rs in 1:4)
  { rl<-sep[rs]
  RMSE[seq(rl,rl+5),5]<-rs
  RMSE[seq(rl,rl+5),6]<-etav1[rs]
  }
  RMSE[,7]<-iti
  names(RMSE)[6]="PointNames" 
  RMSEi=rbind(RMSEi,RMSE)
  
  Etav<-data.frame(cbind(Etav,model,wup,etavalN,iteri))
  
  for (rs in 0:5)
  {Etav[sep+rs,4]<-mods[rs+1]}
  
  for (rs in 1:4)
  { rl<-sep[rs]
  Etav[seq(rl,rl+5),5]<-rs
  Etav[seq(rl,rl+5),6]<-etav1[rs]
  }
  Etav[,7]<-iti
  names(Etav)[6]="PointNames"
  
  BiasF<-data.frame(cbind(BiasF,model,wup,etavalN,iteri))
  
  for (rs in 0:5)
  {BiasF[sep+rs,4]<-mods[rs+1]}
  
  for (rs in 1:4)
  { rl<-sep[rs]
  BiasF[seq(rl,rl+5),5]<-rs
  BiasF[seq(rl,rl+5),6]<-nbpt[rs]
  }
  names(BiasF)[6]="PointNames"
}
}


cols3<-c("red","green","blue","darkorange","purple","pink")
pd <- position_dodge(.2) # move them .05 to the left and right

wgalt<-weighecli[which(weighecli$model=="GalambosCop"),]

ggplot(weighecli[which(weighecli$model=="GalambosCop"),], aes(x=iteri, y=median, colour=iteri)) +
  geom_boxplot(aes(ymin=Q2.5, ymax=Q97.5), width=.7,size=0.7, position=pd) +
  xlab("sd of margins")+
ylab("Eucliddean distance to the reference curve")   +    
  ggtitle(paste0("Euclidean distance level curves - p = 0.001 - n = ",nbpo)) +
  theme_bw() 

# euclidisti$difference<-euclidisti$median-weighecli$median

weighecli$sizeCI<-weighecli$Q97.5-weighecli$Q2.5
euclidisti$sizeCI<-euclidisti$Q97.5-euclidisti$Q2.5
RMSEi$sizeCI<-RMSEi$Q97.5-RMSEi$Q2.5
# 
# write.csv(euclidisti, file="d_AI_Val.csv")
# write.csv(RMSEi, file="R_AI_Val.csv")

# RMSEi[,c(1,8)][which(RMSEi$model=="Gumcop"| RMSEi$model=="GalambosCop"),]=NA
# RMSEi$median<-RMSEi$median/RMSEi$chivals
# write.csv(RMSEit, file="rmsecoefs_all.csv")

wupishit<-c(rep(5,6),rep(6,6),rep(7,6),rep(8,6),rep(9,6),rep(10,6))
wupis<-rep(wupishit,6)
wupishiti<-c(rep(1,6),rep(2,6),rep(3,6),rep(4,6))
wupisi<-rep(wupishiti,6)
RMSEit$wup[c(1:216)]=wupis
RMSEit$wup[c(217:360)]=wupisi
# 
RMSEit<-read.csv(file="rmsecoefs_all.csv")[,-1]
# RMSEit$median<-RMSEit$median/RMSEit$chivals
# ######################################################################
# 
# 
#Verifying consistancy
#

weighecliAD_Cal<-read.csv(file="data/Model_comp_metadata/wd_AD_Cal.csv")[,-1]
weighecliAI_Cal<-read.csv(file="data/Model_comp_metadata/wd_AI_Cal.csv")[,-1]
weighecliAD_Val<-read.csv(file="data/Model_comp_metadata/wd_AD_Val.csv")[,-1]
weighecliAI_Val<-read.csv(file="data/Model_comp_metadata/wd_AI_Val.csv")[,-1]

euclidistAD_Cal<-read.csv(file="data/Model_comp_metadata/d_AD_Cal.csv")[,-1]
euclidistAI_Cal<-read.csv(file="data/Model_comp_metadata/d_AI_Cal.csv")[,-1]
euclidistAD_Val<-read.csv(file="data/Model_comp_metadata/d_AD_Val.csv")[,-1]
euclidistAI_Val<-read.csv(file="data/Model_comp_metadata/d_AI_Val.csv")[,-1]

RMSEiAD_Cal<-read.csv(file="data/Model_comp_metadata/R_AD_Cal.csv")[,-1]
RMSEiAI_Cal<-read.csv(file="data/Model_comp_metadata/R_AI_Cal.csv")[,-1]
RMSEiAD_Val<-read.csv(file="data/Model_comp_metadata/R_AD_Val.csv")[,-1]
RMSEiAI_Val<-read.csv(file="data/Model_comp_metadata/R_AI_Val.csv")[,-1]

weighecliAI_Cal[which(weighecliAI_Cal$model=="Gumcop"& weighecliAI_Cal$PointNames<0.75),c(1,8)]=NA 
weighecliAI_Cal[which(weighecliAI_Cal$model=="GalambosCop"& weighecliAI_Cal$PointNames<0.75),c(1,8)]=NA 
weighecliAI_Val[which(weighecliAI_Val$model=="Gumcop"& weighecliAI_Val$PointNames<0.75),c(1,8)]=NA 
weighecliAI_Val[which(weighecliAI_Val$model=="GalambosCop"& weighecliAI_Val$PointNames<0.75),c(1,8)]=NA 

euclidistAI_Cal[which(euclidistAI_Cal$model=="Gumcop"& euclidistAI_Cal$PointNames<0.75),c(1,8)]=NA 
euclidistAI_Cal[which(euclidistAI_Cal$model=="GalambosCop"& euclidistAI_Cal$PointNames<0.75),c(1,8)]=NA 
euclidistAI_Val[which(euclidistAI_Val$model=="Gumcop"& euclidistAI_Val$PointNames<0.75),c(1,8)]=NA 
euclidistAI_Val[which(euclidistAI_Val$model=="GalambosCop"& euclidistAI_Val$PointNames<0.75),c(1,8)]=NA 

weighecliAD_Cal$wup=weighecliAD_Cal$wup+4
weighecliAD_Val$wup=weighecliAD_Val$wup+4

euclidistAD_Cal$wup=euclidistAD_Cal$wup+4
euclidistAD_Val$wup=euclidistAD_Val$wup+4

RMSEiAD_Cal$wup=RMSEiAD_Cal$wup+4
RMSEiAD_Val$wup=RMSEiAD_Val$wup+4

weighecli_Cal<-rbind(weighecliAD_Cal,weighecliAI_Cal)
weighecli_Val<-rbind(weighecliAD_Val,weighecliAI_Val)

names(euclidistAD_Cal)=names(euclidistAI_Cal)
names(euclidistAD_Val)=names(euclidistAI_Val)

euclidist_Cal<-rbind(euclidistAD_Cal,euclidistAI_Cal)
euclidist_Val<-rbind(euclidistAD_Val,euclidistAI_Val)

names(RMSEiAD_Cal)=names(RMSEiAI_Cal)
names(RMSEiAD_Val)=names(RMSEiAI_Val)

RMSEi_Cal<-rbind(RMSEiAD_Cal,RMSEiAI_Cal)
RMSEi_Val<-rbind(RMSEiAD_Val,RMSEiAI_Val)

plot(abs((weighecli_Cal[,1]-weighecli_Val[,1])/weighecli_Cal[,1]*100))
mean(abs((weighecli_Cal[,1]-weighecli_Val[,1])/weighecli_Cal[,1]*100))

veri<-c()
for (id in 1:360){
veri<-c(veri,inside.range(weighecli_Val[id,1],c(weighecli_Cal[id,2],weighecli_Cal[id,3])))
}
veri



weighecli<-list(weighecli_Cal,weighecli_Val)
names(weighecli)=c("Cal","Val")

plot(weighecli$sizeCI,col=weighecli$wup)
plot(euclidisti$sizeCI,col=weighecli$wup)

mods<-c("Cond-Ex","JT-KDE","Gumcop","Normalcop","FGMcop","GalambosCop")
cols <- c("Cond-Ex" = "purple", "JT-KDE" = "green", "Gumcop" = "red", "Normalcop" = "skyblue", "FGMcop" = "darkorange", "GalambosCop" = "darkred") 
type="eta"

myweight=weighecliAI_Cal
euclidisti=euclidistAI_Cal
RMSEit=RMSEiAI_Cal

weighecliAI_emp<-read.csv(file="out/empirix/weucli_AI.csv")
names(weighecliAI_emp)=c("id","median","Q2.5","Q97.5")
weighecliAI_emp$wup=rep(c(1,2,3,4),6)
weighecliAI_emp$iteri<-as.factor(c(rep("0.25_0.25",4),rep("0.25_0.5",4),rep("0.5_0.5",4),
rep("0.25_1.5",4),rep("0.5_1.5",4),rep("1.5_1.5",4)))
weighecliAI_emp$sizeCI<-weighecliAI_emp$Q97.5-weighecliAI_emp$Q2.5

weighecliAD_emp<-read.csv(file="out/empirix/weucli_AD.csv")
names(weighecliAD_emp)=c("id","median","Q2.5","Q97.5")
weighecliAD_emp$wup=rep(c(5,6,7,8,9,10),6)
weighecliAD_emp$iteri<-as.factor(c(rep("0.25_0.25",6),rep("0.25_0.5",6),rep("0.5_0.5",6),
                                   rep("0.25_1.5",6),rep("0.5_1.5",6),rep("1.5_1.5",6)))
weighecliAD_emp$sizeCI<-weighecliAD_emp$Q97.5-weighecliAD_emp$Q2.5

weighecliAD_Cal$id<-sort(c(rep(c(1:36),6)))
weighecliAD_Val$id<-sort(c(rep(c(1:36),6)))

weighecliAI_Cal$id<-sort(c(rep(c(1:24),6)))
weighecliAI_Val$id<-sort(c(rep(c(1:24),6)))

bmex<-c()
bmex2<-c()
bmcx<-c()
bmcx2<-c()
mode="AI"
if(mode=="AI"){
  for (id in 1:24){
    samp<-weighecliAI_Cal[which(weighecliAI_Cal$id==id),]
    saml<-weighecliAI_Val[which(weighecliAI_Val$id==id),]
    
    bme<-1-samp$median/weighecliAI_emp$median[id]
    bme2<-1-saml$median/weighecliAI_emp$median[id]
    
    bmc<-1-samp$sizeCI/weighecliAI_emp$sizeCI[id]
    bmc2<-1-saml$median/weighecliAI_emp$sizeCI[id]
  
    bmex<-c(bmex,bme)
    bmex2<-c(bmex2,bme2)
    bmcx<-c(bmcx,bmc)
    bmcx2<-c(bmcx2,bmc2)
  }
}
if(mode=="AD"){
  for (id in 1:36){
    samp<-weighecliAD_Cal[which(weighecliAD_Cal$id==id),]
    saml<-weighecliAD_Val[which(weighecliAD_Val$id==id),]
    
    bme<-1-samp$median/weighecliAD_emp$median[id]
    bme2<-1-saml$median/weighecliAD_emp$median[id]
    
    bmc<-1-samp$sizeCI/weighecliAD_emp$sizeCI[id]
    bmc2<-1-saml$median/weighecliAD_emp$sizeCI[id]
    
    bmex<-c(bmex,bme)
    bmex2<-c(bmex2,bme2)
    bmcx<-c(bmcx,bmc)
    bmcx2<-c(bmcx2,bmc2)
  }
}
length(bmcx[which(bmcx>=0)])
length(bmcx2[which(bmcx2>=0)])
bmex[which(bmex<0)]=NA
bmex2[which(bmex2<0)]=NA
plot(bmex, bmex2)
abline(a=0,b=1)
bmcx[which(bmcx<0)]=NA
bmcx2[which(bmcx2<0)]=NA

if(mode=="AI"){
weighecliAI_Cal$bme<-round(bmex,2)
weighecliAI_Cal$bms<-0
weighecliAI_Cal$bms[which(!is.na(bmex))]<-1

weighecliAI_Val$bme<-round(bmex2,2)
weighecliAI_Val$bms<-0
weighecliAI_Val$bms[which(!is.na(bmex2))]<-1
weighecliAI_Cal$bmc<-round(bmcx,2)
weighecliAI_Val$bmc<-round(bmcx2,2)
}

if(mode=="AD"){
  weighecliAD_Cal$bme<-round(bmex,2)
  weighecliAD_Cal$bms<-0
  weighecliAD_Cal$bms[which(!is.na(bmex))]<-1
  weighecliAD_Val$bme<-round(bmex2,2)
  weighecliAD_Val$bms<-0
  weighecliAD_Val$bms[which(!is.na(bmex2))]<-1
  weighecliAD_Cal$bmc<-round(bmcx,2)
  weighecliAD_Val$bmc<-round(bmcx2,2)
}


myweight=rbind(weighecliAI_Val,weighecliAD_Val)
euclidisti=euclidistAD_Cal
RMSEit=RMSEiAD_Cal

length(weighecliAD_Cal$median[which(weighecliAD_Cal$bms==1 & as.character(weighecliAD_Cal$model)=="Cond-Ex")])
length(weighecliAD_Val$median[which(weighecliAD_Val$bms==1)])
length(weighecliAI_Cal$median[which(weighecliAI_Cal$bms==1)])
length(weighecliAI_Val$median[which(weighecliAI_Val$bms==1)])

bymodW<-list()
bymodE<-list()
bymodR<-list()

for (imo in 1:6){
 bymodW<-c(bymodW,list(myweight[which(myweight$model==mods[imo]),] ))
 bymodE<-c(bymodE,list(euclidisti[which(euclidisti$model==mods[imo]),] ))
 bymodR<-c(bymodR,list(RMSEit[which(RMSEit$model==mods[imo]),] ))


}
names(bymodW)=mods
names(bymodE)=mods
names(bymodR)=mods

maxiP<-list()
maxiP1<-list()
maxiP2<-list()
for (kak in 1:6){
if(type=="chi"){
  br<-c(5:10)
  labs=c("0.05","0.1","0.3","0.5","0.7","0.9")
}
if(type=="eta"){
  br<-c(1:4)
  labs=c("0.25","0.5","0.75","0.9") 
}
if (type=="all"){
  br<-c(1:10) 
  labs=c("0.25","0.5","0.75","0.9","0.05","0.1","0.3","0.5","0.7","0.9") 
}

ylabs<-c("A-A","A-B","A-C","B-B","B-C","C-C")
colo <- c("1" = "blue", "0" = "black")
maxiP[[kak]]<-ggplot(data = bymodE[[kak]], aes(x=iteri, y=wup, fill=median)) + 
  geom_tile(aes(fill= median,width=0.8, height=0.8),color ="black",size=10*(bymodE[[kak]]$sizeCI))+
  scale_fill_gradient2(low = "darkgreen", high = "red", mid = "yellow", 
                       midpoint = log(mean(euclidisti$median,na.rm=T)), limit = c((min(euclidisti$median)),(quantile(euclidisti$median,.98,na.rm=T))), space = "Lab", 
                       name="Euclidean\nDistance",trans = "log") +
  # scale_color_ma(low = "skyblue", high = "red", 
  #                        limit = c(0,quantile(euclidisti$sizeCI,.975,na.rm=T)), space = "Lab", 
  #                       name="Euclidean\nDistance")+

  theme_minimal()+ 
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=br,labels=labs)+
  scale_x_discrete(labels=ylabs)+
  labs(title=mods[kak],
       x =  "Marginal distributions", y =expression(paste(eta," value")))+
  geom_text(aes(iteri, wup, label = round(median,2)), color = "black", size = 3)+
  coord_fixed()


maxiP1[[kak]]<-ggplot(data = bymodW[[kak]], aes(x=iteri, y=wup, fill=median)) + 
  geom_tile(aes(fill= median,width=0.8, height=0.8),color="black",size=log(20*(bymodW[[kak]]$sizeCI)+1))+
  scale_fill_gradient2(low = "darkgreen", high = "red", mid = "yellow", 
                       midpoint = log(mean(myweight$median,na.rm=T)), 
                       name="Euclidean\nDistance",trans="log") +
  # scale_color_gradient(low = "skyblue", high = "red", 
  #                      limit = c(0,quantile(myweight$sizeCI,.975,na.rm=T)), space = "Lab", 
  #                      name="Euclidean\nDistance")+
  scale_colour_manual(values = colo)+
  theme_minimal()+ 
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=br,labels=labs)+
  scale_x_discrete(labels=ylabs)+
  labs(title=mods[kak],
       x =  "Marginal distributions", y =expression(paste(eta," value")))+
  geom_text(aes(iteri, wup, label = round(median,2),color = factor(bms),fontface=bms+1), size = 4)+
  coord_fixed()

maxiP2[[kak]]<-ggplot(data = bymodR[[kak]], aes(x=iteri, y=wup, fill=median)) + 
  geom_tile(aes(fill= median,width=0.8, height=0.8),color ="black",size=10*(bymodR[[kak]]$sizeCI))+
  scale_fill_gradient2(low = "darkgreen", high = "red", mid = "yellow", 
                       midpoint =-1,space = "Lab", 
                       name="Euclidean\nDistance",trans="log10") +
  # scale_color_gradient2(low = "skyblue", high = "red", mid = "mediumpurple2", 
  #                       midpoint = 0.05, limit = c(0,0.14) space = "Lab", 
  #                       name="Euclidean\nDistance")+
  theme_minimal()+ 
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=br,labels=labs)+
  scale_x_discrete(labels=ylabs)+
  labs(title=mods[kak],
       x =  "Marginal distributions", y =expression(paste(eta," value")))+
  geom_text(aes(iteri, wup, label = round(median,2)), color = "black", size = 3)+
  coord_fixed()

}

lay <- rbind(c(1,1,2,2,3,3),
             c(4,4,5,5,6,6))


pdf(paste0(d,"OneDrive - King's College London/DATA/CODE ANALYSIS/Model comparison/Bivariate_Models/bivariate_models/out/MatrixCIsAIf_CoefsVa_",type,"_",mars[1],"_",mars[2],"_",nbpo,".pdf"),width = 11.6, height = 8.2)

grid.arrange(maxiP[[1]],maxiP[[2]],maxiP[[3]],
             maxiP[[4]],maxiP[[5]],maxiP[[6]], layout_matrix = lay, top="Normalized Euclidean Distance")
grid.arrange(maxiP1[[1]],maxiP1[[2]],maxiP1[[3]],
             maxiP1[[4]],maxiP1[[5]],maxiP1[[6]], layout_matrix = lay, top="Weighted Normalized Euclidean Distance")
grid.arrange(maxiP2[[1]],maxiP2[[2]],maxiP2[[3]],
             maxiP2[[4]],maxiP2[[5]],maxiP2[[6]], layout_matrix = lay, top="Root Mean Squared Normalized Bias")
dev.off()

grid <- expand.grid(lon=unique(bymodE[[1]]$Chinames), lat=unique(bymodE[[1]]$iteri))
x<-euclidistit[order(euclidistit$model,euclidistit$Chinames),]

Sales <- array(NA, c(36, 10))
for (i in 1:36) 
  for(j in (1:10))
  Sales[i,j] <- x[(10*(i-1)+(j)), "median"]
rownames(Sales)=rep(iter,6)
  
colnames(Sales)= c("0.05","0.1","0.25","0.3","0.5","0.5","0.7","0.75","0.9","0.9")
  
write.csv(euclidistit, file="euclidist_Totolit.csv")


# }



  
matB<-bymodW
mean(matB$`Cond-Ex`$median[which(matB$`Cond-Ex`$wup<=4)])
mean(matB$`JT-KDE`$median[which(matB$`Cond-Ex`$wup>4)])
pmM<-c()
pmQM<-c()
pmQm<-c()
pm1<-c()
for (ziz in 1:6){
# {idi<-seq(0.01,0.1,by=0.01)
# pm<-idi
# for(zg in 1:10){
  # pmQM[zg]<-(length(matB[[ziz]]$Q97.5[which(matB[[ziz]]$Q97.5<idi[zg])])/length(matB[[ziz]]$Q97.5))*100
  pmM[ziz]<-(length(matB[[ziz]]$median[which(matB[[ziz]]$bms>0)])/length(matB[[ziz]]$median))*100
  # pmQm[zg]<-(length(matB[[ziz]]$Q2.5[which(matB[[ziz]]$Q2.5<idi[zg])])/length(matB[[ziz]]$Q2.5))*100
# }
# pmb<-data.frame(idi,pmQm,pmM,pmQM,mods[ziz])
# pm1<-as.data.frame(rbind(pm1,pmb))
}

# colnames(pm1)=c("id",mods)
# plot(pm1$`Cond-Ex`,type="o",col=1)
# points(pm1[,2],type="o",col=2)
# points(pm1[,3],type="o",col=3)
# points(pm1[,4],type="o",col=4)
# points(pm1[,5],type="o",col=5)
# points(pm1[,6],type="o",col=6)

ggplot(pm1, aes(x=idi, y=pmM, colour=mods.ziz.))  +
  geom_line(size=.7) +
  geom_point(size=3,shape=21, fill="white") +
  xlab("wd threshold") +
  ylab("% of occurence below threshold") +
  scale_colour_manual(values = cols,breaks=mods,labels=mods)   +    
  scale_x_continuous(breaks=seq(0.01,0.1,by=0.01),labels=idi)+
  ggtitle(paste0("Euclidean distance level curves - p = 0.001 - n = ",nbpo)) +
  # expand_limits(y=1) +                        # Expand y range
  scale_y_continuous(limits=c(0,100)) +         # Set tick every 4
  theme_bw()

pm1[which(pm1$idi>0.099),]

pm_Cal[which(pm_Cal$idi>0.099),]
pm_Val[which(pm_Val$idi>0.099),]

pmu<-c()
pmu1<-c()
for (ziz in 1:6)
{idi<-seq(0.01,0.1,by=0.01)
pmu<-idi
for(zg in 1:10){
  pmu[zg]<-(length(bymodE[[ziz]]$median[which(bymodE[[ziz]]$median<=idi[zg])])/length(bymodE[[ziz]]$median))*100
}
pmb1<-data.frame(idi,pmu,mods[ziz])
pmu1<-as.data.frame(rbind(pmu1,pmb1))
}

ggplot(pmu1, aes(x=idi, y=pmu, colour=mods.ziz.))  +
  geom_line(size=.7) +
  geom_point(size=3,shape=21, fill="white")+
  xlab("wd threshold") +
  ylab("% of occurence below threshold") +
  scale_colour_manual(values = cols,breaks=mods,labels=mods)   +    
  scale_x_continuous(breaks=seq(0.01,0.1,by=0.01),labels=idi)+
  ggtitle(paste0("Euclidean distance level curves - p = 0.001 - n = ",nbpo)) +
  # expand_limits(y=1) +                        # Expand y range
  scale_y_continuous(limits=c(0,100)) +         # Set tick every 4
  theme_bw() 

