#clear variables and console
rm(list=ls())
cat("\014") 

# bias analysis for incomplete genotyping & LOH
# NLM project prototype
# Started: April 2022
# Edits: Switched order so LOH adjusted b4 incompleteness
#        reamed file from Lash et al 2011 v4.r
#        9/20/22: redoing imputation for missing snp phenotype
#        2/28/23: checking, using new data, finally wrapping this up
#        4/13/23: Altering phenotype imputation. Now, impute the phenotype
#                  using study phenotype, if available.
#       12/29/23: Final analysis and meta analysis code
#        2/15/24: Final check - minor changes in phenotype algorithms
#       8/11/24: updated cpic GS algorithm for *9 and *41
    
#Read in needed libraries
library(MASS)
library(mvtnorm)
library(mc2d)
library(parallel)
library(matrixStats)
library(meta)
library(readxl)
#library(rstan)
#library(rjags)
#library(coda)
library(ggplot2)
library(haven)
  
#Read in additional functions
#bias functions.r contains most of the functions for the main bias analysis
#data.r contains functions to read in and process data
#xtra fxns.r contains sundry simple functions
#genetic fxns.r contans function used for processing the genetic data (phenotyping)
  
    
source("~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/cyp2d6 bias functions.r")
source("~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/cyp2d6 data.r")
source("~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/xtra fxns.r")
source("~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/genetic fxns.r")
    
#######################################
# Parameters for what data to use and # of iterations
niter <- 10^5 #number of iterations
  
#which concordance data to use
  # 1= all
  # 2= non goetz
  # 3= goetz
concord=1

# Set up data matrices
Data.meta.pm.vs.nm=as.data.frame(matrix(data=0,nrow=40,ncol=16))
colnames(Data.meta.pm.vs.nm)<-c("Effect","Low","Up","Eff.c","Low.c","Up.c","Year","Region","Samp.Type","vv.1","vw.1","ww.1","vv.0","vw.0","ww.0","niter")
Data.meta.pm.vs.nm$Name=c("Temp")
Data.meta.pm.vs.nm$eff.loh=0
Data.meta.pm.vs.nm$low.loh=0
Data.meta.pm.vs.nm$up.loh=0
Data.meta.pmim.vs.nm<-Data.meta.im.vs.nm<-Data.meta.pm.vs.nm
#same data frames for asian studies
Data.meta.im.vs.nm.as<-Data.meta.pm.vs.nm
  
  
#run adjustment for EUROPEAN data with summary level data
#Studies with cell counts 
eur.data=c(1:6,8:10,25,26,29)
#studies with only aggregate data
eur.data.eff=c(7,11,12,15,16,18:23)
as.data=c(27,28,30,33)
as.data.eff=c(24,31,32,34,36,38)
#conduct bias analysis for studies with individual level european data
for (i in eur.data){
  Study.id=i
  #read in study data
  sdl<-bias.setup(Study.id)
  res<-cyp.loh.inc(sdl$vv.1,sdl$vw.1,sdl$ww.1,sdl$vv.0,sdl$vw.0,sdl$ww.0,sdl$cohort,sdl$samp.type,concord,sdl$sesp.prior.mat,niter)
  #vv.1<-sdl$vv.1;vw.1<-sdl$vw.1;ww.1<-sdl$ww.1;vv.0<-sdl$vv.0;vw.0<-sdl$vw.0;ww.0<-sdl$ww.0;cohort<-sdl$cohort;samp.type<-sdl$samp.type;sesp.prior.mat<-sdl$sesp.prior.mat;niter<-100
  if (sdl$samp.type=="n"){
    Data.meta.pm.vs.nm[i,]<-c(res$Hom[3,],res$Hom[1,],sdl$year,"Europe",sdl$samp.type,sdl$vv.1,sdl$vw.1,sdl$ww.1,sdl$vv.0,sdl$vw.0,sdl$ww.0,res$Drop[2],sdl$author,res$Hom[1,])
    Data.meta.im.vs.nm[i,]<-c(res$Het[3,],res$Het[1,],sdl$year,"Europe",sdl$samp.type,sdl$vv.1,sdl$vw.1,sdl$ww.1,sdl$vv.0,sdl$vw.0,sdl$ww.0,res$Drop[2],sdl$author,res$Het[1,])
    Data.meta.pmim.vs.nm[i,]<-c(res$Any[3,],res$Any[1,],sdl$year,"Europe",sdl$samp.type,sdl$vv.1,sdl$vw.1,sdl$ww.1,sdl$vv.0,sdl$vw.0,sdl$ww.0,res$Drop[2],sdl$author,res$Any[1,])
  }
  if (sdl$samp.type=="t"){
    Data.meta.pm.vs.nm[i,]<-c(res$Hom[3,],res$Hom[1,],sdl$year,"Europe",sdl$samp.type,sdl$vv.1,sdl$vw.1,sdl$ww.1,sdl$vv.0,sdl$vw.0,sdl$ww.0,res$Drop[2],sdl$author,res$Hom[2,])
    Data.meta.im.vs.nm[i,]<-c(res$Het[3,],res$Het[1,],sdl$year,"Europe",sdl$samp.type,sdl$vv.1,sdl$vw.1,sdl$ww.1,sdl$vv.0,sdl$vw.0,sdl$ww.0,res$Drop[2],sdl$author,res$Het[2,])
    Data.meta.pmim.vs.nm[i,]<-c(res$Any[3,],res$Any[1,],sdl$year,"Europe",sdl$samp.type,sdl$vv.1,sdl$vw.1,sdl$ww.1,sdl$vv.0,sdl$vw.0,sdl$ww.0,res$Drop[2],sdl$author,res$Any[2,])
  }
  cat(paste(i, sdl$author, round(res$Drop[2]/res$Drop[1],2), sep = " "), "\n")
  rm(sdl,res)
}
    
#Abstract EUROPEAN studies without summary level data:
for (i in eur.data.eff){
  #Read in study data
  study.data.2<-extract.cyp.data.2(i)
  #put study data in pm vs nm if ecode=1
  if (study.data.2[8]=="1"){
    if (study.data.2[11]=="n"){
      Data.meta.pm.vs.nm[i,]<-c(rep("NA",3),study.data.2[1:3],study.data.2[14],"Europe",study.data.2[11],rep("NA",7),study.data.2[13],study.data.2[1:3])
    }
    if (study.data.2[11]=="t"){
      Data.meta.pm.vs.nm[i,]<-c(rep("NA",3),study.data.2[1:3],study.data.2[14],"Europe",study.data.2[11],rep("NA",7),study.data.2[13],rep(NA,3))
    }
  }
  # put study data in pm.im vs nm if ecode=2
  if (study.data.2[8]=="2"){
    if (study.data.2[11]=="n"){
      Data.meta.pmim.vs.nm[i,]<-c(rep("NA",3),study.data.2[1:3],study.data.2[14],"Europe",study.data.2[11],rep("NA",7),study.data.2[13],study.data.2[1:3])
    }
    if (study.data.2[11]=="t"){
      Data.meta.pmim.vs.nm[i,]<-c(rep("NA",3),study.data.2[1:3],study.data.2[14],"Europe",study.data.2[11],rep("NA",7),study.data.2[13],rep(NA,3))
    }
  }
  if (study.data.2[9]=="2"){
    if (study.data.2[11]=="n"){
      Data.meta.pmim.vs.nm[i,]<-c(rep("NA",3),study.data.2[4:6],study.data.2[14],"Europe",study.data.2[11],rep("NA",7),study.data.2[13],study.data.2[4:6])
    }
    if (study.data.2[11]=="t"){
      Data.meta.pmim.vs.nm[i,]<-c(rep("NA",3),study.data.2[4:6],study.data.2[14],"Europe",study.data.2[11],rep("NA",7),study.data.2[13],rep(NA,3))
    }
  }
  
  # put study data in im vs nm if ecode==4
  #only happens for second effect estimates (4:6 in dataframe)
  if (study.data.2[9]=="4"){
    if (study.data.2[11]=="n"){
      Data.meta.im.vs.nm[i,]<-c(rep("NA",3),study.data.2[4:6],study.data.2[14],"Europe",study.data.2[11],rep("NA",7),study.data.2[13],study.data.2[4:6])
    }
    if (study.data.2[11]=="t"){
      Data.meta.im.vs.nm[i,]<-c(rep("NA",3),study.data.2[4:6],study.data.2[14],"Europe",study.data.2[11],rep("NA",7),study.data.2[13],rep(NA,3))
    }
  }
}

#Run adjustment for ASIAN studies
for (i in as.data){
  Study.id=i
  sdl<-bias.setup(Study.id)
  res<-cyp.loh.inc.asian(sdl$vv.1,sdl$vw.1,sdl$ww.1,sdl$vv.0,sdl$vw.0,sdl$ww.0,sdl$cohort,sdl$samp.type,concord,sdl$sesp.prior.mat,niter)
  #vv.1<-sdl$vv.1;vw.1<-sdl$vw.1;ww.1<-sdl$ww.1;vv.0<-sdl$vv.0;vw.0<-sdl$vw.0;ww.0<-sdl$ww.0;cohort<-sdl$cohort;samp.type<-sdl$samp.type;sesp.prior.mat<-sdl$sesp.prior.mat;
  if (sdl$samp.type=="n"){
    Data.meta.im.vs.nm.as[i,]<-c(res$Het[3,],res$Het[1,],sdl$year,"Asia",sdl$samp.type,sdl$vv.1,sdl$vw.1,sdl$ww.1,sdl$vv.0,sdl$vw.0,sdl$ww.0,res$Drop[2],sdl$author,res$Het[1,])
  }
  if (sdl$samp.type=="t"){
    Data.meta.im.vs.nm.as[i,]<-c(res$Het[3,],res$Het[1,],sdl$year,"Asia",sdl$samp.type,sdl$vv.1,sdl$vw.1,sdl$ww.1,sdl$vv.0,sdl$vw.0,sdl$ww.0,res$Drop[2],sdl$author,rep(NA,3))
  }
  cat(paste(i, sdl$author, round(res$Drop[2]/res$Drop[1],2), sep = " "), "\n")
  rm(res,sdl)
}
  
#Abstract ASIAN studies without summary level data:
for (i in as.data.eff){
  study.data.2<-extract.cyp.data.2(i)
  if (study.data.2[11]=="n"&study.data.2[8]!=5){
    Data.meta.im.vs.nm.as[i,]<-c(rep("NA",3),study.data.2[1:3],study.data.2[14],"Asia",study.data.2[11],rep("NA",7),study.data.2[13],study.data.2[1:3])
  }
  if (study.data.2[11]=="t"&study.data.2[8]!=5){
    Data.meta.im.vs.nm.as[i,]<-c(rep("NA",3),study.data.2[1:3],study.data.2[14],"Asia",study.data.2[11],rep("NA",7),study.data.2[13],rep(NA,3))
  }
}
    
    
#recompute the crude effect for Asian studies that have gtype/ptype data
#  In these studies, the IM vs NM effect can't be calculated directly from ptype/gtype
#  e.g., in studies with *10, the 'vv' should be IM and the 'vw' should be NM according to Caudle
#        studies with ptype have similar issue.
#        reclassification scheme by study: 
#           Chamn: fine as is. do nothing
#           Kiyatoni: PM-->IM and IM-->NM (most of their IM category are 10/1; most of PM is 10/10)
#           Park: PM-->IM and IM-->NM (same as above)
#           Sirach: has *10 data.  10/10->IM  10/1-->NM
#           Lan: same as Sirach

for (index in c(30,33,27,28)){
  abcd<-c(as.numeric(Data.meta.im.vs.nm.as$vv.1[index]),
          as.numeric(Data.meta.im.vs.nm.as$vw.1[index])+
          as.numeric(Data.meta.im.vs.nm.as$ww.1[index]),
          as.numeric(Data.meta.im.vs.nm.as$vv.0[index]),
          as.numeric(Data.meta.im.vs.nm.as$vw.0[index])+
          as.numeric(Data.meta.im.vs.nm.as$ww.0[index]))
  if (index==30){
    Data.meta.im.vs.nm.as[index,4:6]<-Data.meta.im.vs.nm.as[index,18:20]<-
      do.call(crude.or.ci, as.list(abcd))[c(1,3,4)]
  }
  if (index !=30){
    Data.meta.im.vs.nm.as[index,4:6]<-Data.meta.im.vs.nm.as[index,18:20]<-
      do.call(crude.rr.ci, as.list(abcd))[c(1,3,4)]
  }
  Data.meta.im.vs.nm.as$ww.1[index]<-as.numeric(Data.meta.im.vs.nm.as$ww.1[index])+as.numeric(Data.meta.im.vs.nm.as$vw.1[index])
  Data.meta.im.vs.nm.as$ww.0[index]<-as.numeric(Data.meta.im.vs.nm.as$ww.0[index])+as.numeric(Data.meta.im.vs.nm.as$vw.0[index])
  Data.meta.im.vs.nm.as$vw.1[index]<-Data.meta.im.vs.nm.as$vv.1[index]
  Data.meta.im.vs.nm.as$vw.0[index]<-Data.meta.im.vs.nm.as$vv.0[index]
}

#set Data up for meta analysis and introduce bias adjustment for studies with only point estimates
Data.meta.pm.vs.nm.2<-pt.est.bias.adj(meta.data.setup(Data.meta.pm.vs.nm))
Data.meta.pmim.vs.nm.2<-pt.est.bias.adj(meta.data.setup(Data.meta.pmim.vs.nm))
Data.meta.im.vs.nm.2<-pt.est.bias.adj(meta.data.setup(Data.meta.im.vs.nm))
Data.meta.im.vs.nm.as.2<-pt.est.bias.adj(meta.data.setup(Data.meta.im.vs.nm.as))
Data.meta.im.vs.nm.all.2<-rbind(Data.meta.im.vs.nm.2,Data.meta.im.vs.nm.as.2)
    
#export data for tables
Data.meta.im.vs.nm.all.2$type<-"IM vs NM"
Data.meta.pm.vs.nm.2$type<-"PM vs NM"
Data.meta.pmim.vs.nm.2$type<-"PM/IM vs NM"
merged<-rbind(Data.meta.pm.vs.nm.2,Data.meta.im.vs.nm.all.2,Data.meta.pmim.vs.nm.2)
merged<-merged[,which(colnames(merged) %in% c("Name","Year","Region","Samp.Type","Eff.c","Low.c","Up.c","vv.1","vw.1","ww.1","vv.0","vw.0","ww.0","type"))]
merged[,c(1:3)]<-round(merged[,c(1:3)],digits=1)
merged$rrc <- paste(merged$Eff.c, " (", merged$Low.c, ", ", merged$Up.c, ")", sep = "")
merged.out<-merged[,c("Name","Year","Region","Samp.Type","type","rrc","vv.1","vw.1","ww.1","vv.0","vw.0","ww.0")]
merged.out$Samp.Type<-ifelse(merged.out$Samp.Type=="n","Non","Tumor")
table1.out<-merged.out[merged.out$type=="PM vs NM",]
table2.out<-merged.out[merged.out$type=="IM vs NM",]
table3.out<-merged.out[merged.out$type=="PM/IM vs NM",]
table1.out<-table1.out[,c("Name","Year","Region","Samp.Type","type","rrc","vv.1","ww.1","vv.0","ww.0")]
table2.out<-table2.out[,c("Name","Year","Region","Samp.Type","type","rrc","vw.1","ww.1","vw.0","ww.0")]
table3.out$vv.vw.1<-as.numeric(table3.out$vv.1)+as.numeric(table3.out$vw.1)
table3.out$vv.vw.0<-as.numeric(table3.out$vv.0)+as.numeric(table3.out$vw.0)
    
table3.out<-table3.out[,c("Name","Year","Region","Samp.Type","type","rrc","vv.vw.1","ww.1","vv.vw.0","ww.0")]

#Fix crude data in Ismail...For that study +1 was added to cell counts to aid in model fit
# remove the +1 here. Still reporting the RR based on +1...little difference between that and exact
table1.out[17,7:10]<-c(4,45,0,33)
table2.out[15,7:10]<-c(8,45,5,33)
table3.out[17,7:10]<-c(12,45,5,33)    
write.csv(table1.out, "~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/results/table1.csv")
write.csv(table2.out, "~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/results/table2.csv")
write.csv(table3.out, "~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/results/table3.csv")
    
# Meta Analysis: PM vs NM
#Note: these meta analyses are for double checking. Stata analyses are main ones
# Stata uses a different formula for I2 than R
# The stata formula is more common (See Higgins and Thompson 2002)
meta.out.f<-metagen(log(as.numeric(Effect)),lower=log(as.numeric(Low)),upper=log(as.numeric(Up)),data=Data.meta.pm.vs.nm.2,sm="RR",studlab=Name)
meta.out.l<-metagen(log(as.numeric(eff.loh)),lower=log(as.numeric(low.loh)),upper=log(as.numeric(up.loh)),data=Data.meta.pm.vs.nm.2,sm="RR",studlab=Name)
meta.out.c<-metagen(log(as.numeric(Eff.c)),lower=log(as.numeric(Low.c)),upper=log(as.numeric(Up.c)),data=Data.meta.pm.vs.nm.2,sm="RR",studlab=Name)
table.pm<-matrix(0,nrow=3,ncol=5)
colnames(table.pm)<-c("RR","LCI","UCI","p-val","I2")
rownames(table.pm)<-c("Crude","LOH","Full")
table.pm[3,]<-c(exp(meta.out.f$TE.random),exp(meta.out.f$lower.random),exp(meta.out.f$upper.random),meta.out.f$pval.Q,meta.out.f$I2)
table.pm[2,]<-c(exp(meta.out.l$TE.random),exp(meta.out.l$lower.random),exp(meta.out.l$upper.random),meta.out.l$pval.Q,meta.out.l$I2)
table.pm[1,]<-c(exp(meta.out.c$TE.random),exp(meta.out.c$lower.random),exp(meta.out.c$upper.random),meta.out.c$pval.Q,meta.out.c$I2)

    
#sens analyses: drop studies with larger number of missing 
to.drop=which(Data.meta.pm.vs.nm.2$Name %in% c("Thompson","Wegman","Ismail Al khalil","Chamnanphon","Teh"))
Data.meta.pm.vs.nm.2.sens<-Data.meta.pm.vs.nm.2[-c(to.drop),]
meta.out.f<-metagen(log(as.numeric(Effect)),lower=log(as.numeric(Low)),upper=log(as.numeric(Up)),data=Data.meta.pm.vs.nm.2.sens,sm="RR",studlab=Name)
meta.out.l<-metagen(log(as.numeric(eff.loh)),lower=log(as.numeric(low.loh)),upper=log(as.numeric(up.loh)),data=Data.meta.pm.vs.nm.2.sens,sm="RR",studlab=Name)
meta.out.c<-metagen(log(as.numeric(Eff.c)),lower=log(as.numeric(Low.c)),upper=log(as.numeric(Up.c)),data=Data.meta.pm.vs.nm.2.sens,sm="RR",studlab=Name)
table.pm.sens1<-matrix(0,nrow=3,ncol=5)
colnames(table.pm.sens1)<-c("RR","LCI","UCI","p-val","I2")
rownames(table.pm.sens1)<-c("Crude","LOH","Full")
table.pm.sens1[3,]<-c(exp(meta.out.f$TE.random),exp(meta.out.f$lower.random),exp(meta.out.f$upper.random),meta.out.f$pval.Q,meta.out.f$I2)
table.pm.sens1[2,]<-c(exp(meta.out.l$TE.random),exp(meta.out.l$lower.random),exp(meta.out.l$upper.random),meta.out.l$pval.Q,meta.out.l$I2)
table.pm.sens1[1,]<-c(exp(meta.out.c$TE.random),exp(meta.out.c$lower.random),exp(meta.out.c$upper.random),meta.out.c$pval.Q,meta.out.c$I2)
    

#sens analyses: drop studies with unclear algorithm
to.drop=which(Data.meta.pm.vs.nm.2$Name %in% c("Thompson","De Ameida Melo","Argalacsova"))
Data.meta.pm.vs.nm.2.sens<-Data.meta.pm.vs.nm.2[-c(to.drop),]
meta.out.f<-metagen(log(as.numeric(Effect)),lower=log(as.numeric(Low)),upper=log(as.numeric(Up)),data=Data.meta.pm.vs.nm.2.sens,sm="RR",studlab=Name)
meta.out.l<-metagen(log(as.numeric(eff.loh)),lower=log(as.numeric(low.loh)),upper=log(as.numeric(up.loh)),data=Data.meta.pm.vs.nm.2.sens,sm="RR",studlab=Name)
meta.out.c<-metagen(log(as.numeric(Eff.c)),lower=log(as.numeric(Low.c)),upper=log(as.numeric(Up.c)),data=Data.meta.pm.vs.nm.2.sens,sm="RR",studlab=Name)
table.pm.sens2<-matrix(0,nrow=3,ncol=5)
colnames(table.pm.sens2)<-c("RR","LCI","UCI","p-val","I2")
rownames(table.pm.sens2)<-c("Crude","LOH","Full")
table.pm.sens2[3,]<-c(exp(meta.out.f$TE.random),exp(meta.out.f$lower.random),exp(meta.out.f$upper.random),meta.out.f$pval.Q,meta.out.f$I2)
table.pm.sens2[2,]<-c(exp(meta.out.l$TE.random),exp(meta.out.l$lower.random),exp(meta.out.l$upper.random),meta.out.l$pval.Q,meta.out.l$I2)
table.pm.sens2[1,]<-c(exp(meta.out.c$TE.random),exp(meta.out.c$lower.random),exp(meta.out.c$upper.random),meta.out.c$pval.Q,meta.out.c$I2)


#sens analyses: drop studies without genotype phenotype data 
to.drop=which(Data.meta.pm.vs.nm.2$Name %in% c("Thompson","Wegman","Goetzb","Gor","Rae"))
Data.meta.pm.vs.nm.2.sens<-Data.meta.pm.vs.nm.2[-c(to.drop),]
meta.out.f<-metagen(log(as.numeric(Effect)),lower=log(as.numeric(Low)),upper=log(as.numeric(Up)),data=Data.meta.pm.vs.nm.2.sens,sm="RR",studlab=Name)
meta.out.l<-metagen(log(as.numeric(eff.loh)),lower=log(as.numeric(low.loh)),upper=log(as.numeric(up.loh)),data=Data.meta.pm.vs.nm.2.sens,sm="RR",studlab=Name)
meta.out.c<-metagen(log(as.numeric(Eff.c)),lower=log(as.numeric(Low.c)),upper=log(as.numeric(Up.c)),data=Data.meta.pm.vs.nm.2.sens,sm="RR",studlab=Name)
table.pm.sens3<-matrix(0,nrow=3,ncol=5)
colnames(table.pm.sens3)<-c("RR","LCI","UCI","p-val","I2")
rownames(table.pm.sens3)<-c("Crude","LOH","Full")
table.pm.sens3[3,]<-c(exp(meta.out.f$TE.random),exp(meta.out.f$lower.random),exp(meta.out.f$upper.random),meta.out.f$pval.Q,meta.out.f$I2)
table.pm.sens3[2,]<-c(exp(meta.out.l$TE.random),exp(meta.out.l$lower.random),exp(meta.out.l$upper.random),meta.out.l$pval.Q,meta.out.l$I2)
table.pm.sens3[1,]<-c(exp(meta.out.c$TE.random),exp(meta.out.c$lower.random),exp(meta.out.c$upper.random),meta.out.c$pval.Q,meta.out.c$I2)

    
# Meta Analysis: IM vs NM
#All (Eur and Asian)
meta.out.f<-metagen(log(as.numeric(Effect)),lower=log(as.numeric(Low)),upper=log(as.numeric(Up)),data=Data.meta.im.vs.nm.all.2,sm="RR",studlab=Name)
meta.out.l<-metagen(log(as.numeric(eff.loh)),lower=log(as.numeric(low.loh)),upper=log(as.numeric(up.loh)),data=Data.meta.im.vs.nm.all.2,sm="RR",studlab=Name)
meta.out.c<-metagen(log(as.numeric(Eff.c)),lower=log(as.numeric(Low.c)),upper=log(as.numeric(Up.c)),data=Data.meta.im.vs.nm.all.2,sm="RR",studlab=Name)
table.im.all<-matrix(0,nrow=3,ncol=5)
colnames(table.im.all)<-c("RR","LCI","UCI","p-val","I2")
table.im.all[3,]<-c(exp(meta.out.f$TE.random),exp(meta.out.f$lower.random),exp(meta.out.f$upper.random),meta.out.f$pval.Q,meta.out.f$I2)
table.im.all[2,]<-c(exp(meta.out.l$TE.random),exp(meta.out.l$lower.random),exp(meta.out.l$upper.random),meta.out.l$pval.Q,meta.out.l$I2)
table.im.all[1,]<-c(exp(meta.out.c$TE.random),exp(meta.out.c$lower.random),exp(meta.out.c$upper.random),meta.out.c$pval.Q,meta.out.c$I2)

#Asian only    
meta.out.f<-metagen(log(as.numeric(Effect)),lower=log(as.numeric(Low)),upper=log(as.numeric(Up)),data=Data.meta.im.vs.nm.as.2,sm="RR",studlab=Name)
meta.out.l<-metagen(log(as.numeric(eff.loh)),lower=log(as.numeric(low.loh)),upper=log(as.numeric(up.loh)),data=Data.meta.im.vs.nm.as.2,sm="RR",studlab=Name)
meta.out.c<-metagen(log(as.numeric(Eff.c)),lower=log(as.numeric(Low.c)),upper=log(as.numeric(Up.c)),data=Data.meta.im.vs.nm.as.2,sm="RR",studlab=Name)
table.im.as<-matrix(0,nrow=3,ncol=5)
colnames(table.im.as)<-c("RR","LCI","UCI","p-val","I2")
table.im.as[3,]<-c(exp(meta.out.f$TE.random),exp(meta.out.f$lower.random),exp(meta.out.f$upper.random),meta.out.f$pval.Q,meta.out.f$I2)
table.im.as[2,]<-c(exp(meta.out.l$TE.random),exp(meta.out.l$lower.random),exp(meta.out.l$upper.random),meta.out.l$pval.Q,meta.out.l$I2)
table.im.as[1,]<-c(exp(meta.out.c$TE.random),exp(meta.out.c$lower.random),exp(meta.out.c$upper.random),meta.out.c$pval.Q,meta.out.c$I2)
    
#European only    
meta.out.f<-metagen(log(as.numeric(Effect)),lower=log(as.numeric(Low)),upper=log(as.numeric(Up)),data=Data.meta.im.vs.nm.2,sm="RR",studlab=Name)
meta.out.l<-metagen(log(as.numeric(eff.loh)),lower=log(as.numeric(low.loh)),upper=log(as.numeric(up.loh)),data=Data.meta.im.vs.nm.2,sm="RR",studlab=Name)
meta.out.c<-metagen(log(as.numeric(Eff.c)),lower=log(as.numeric(Low.c)),upper=log(as.numeric(Up.c)),data=Data.meta.im.vs.nm.2,sm="RR",studlab=Name)
table.im<-matrix(0,nrow=3,ncol=5)
colnames(table.im)<-c("RR","LCI","UCI","p-val","I2")
table.im[3,]<-c(exp(meta.out.f$TE.random),exp(meta.out.f$lower.random),exp(meta.out.f$upper.random),meta.out.f$pval.Q,meta.out.f$I2)
table.im[2,]<-c(exp(meta.out.l$TE.random),exp(meta.out.l$lower.random),exp(meta.out.l$upper.random),meta.out.l$pval.Q,meta.out.l$I2)
table.im[1,]<-c(exp(meta.out.c$TE.random),exp(meta.out.c$lower.random),exp(meta.out.c$upper.random),meta.out.c$pval.Q,meta.out.c$I2)

# Meta Analysis: IM/PM vs NM
meta.out.f<-metagen(log(as.numeric(Effect)),lower=log(as.numeric(Low)),upper=log(as.numeric(Up)),data=Data.meta.pmim.vs.nm.2,sm="RR",studlab=Name)
meta.out.l<-metagen(log(as.numeric(eff.loh)),lower=log(as.numeric(low.loh)),upper=log(as.numeric(up.loh)),data=Data.meta.pmim.vs.nm.2,sm="RR",studlab=Name)
meta.out.c<-metagen(log(as.numeric(Eff.c)),lower=log(as.numeric(Low.c)),upper=log(as.numeric(Up.c)),data=Data.meta.pmim.vs.nm.2,sm="RR",studlab=Name)
table.pmim<-matrix(0,nrow=3,ncol=5)
colnames(table.pmim)<-c("RR","LCI","UCI","p-val","I2")
table.pmim[3,]<-c(exp(meta.out.f$TE.random),exp(meta.out.f$lower.random),exp(meta.out.f$upper.random),meta.out.f$pval.Q,meta.out.f$I2)
table.pmim[2,]<-c(exp(meta.out.l$TE.random),exp(meta.out.l$lower.random),exp(meta.out.l$upper.random),meta.out.l$pval.Q,meta.out.l$I2)
table.pmim[1,]<-c(exp(meta.out.c$TE.random),exp(meta.out.c$lower.random),exp(meta.out.c$upper.random),meta.out.c$pval.Q,meta.out.c$I2)

    
round(table.pm,digits=2)
round(table.im,digits=2)
round(table.im.as,digits=2)
round(table.pmim,digits=2)
  
round(table.pm.sens1,digits=2)
round(table.pm.sens2,digits=2)
round(table.pm.sens3,digits=2)
  


  
#output effects for meta regression analysis in Stata
study.data.excel<-read_excel("~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/meta_data final.xlsx")
study.data.excel$Name<- substr(study.data.excel$author,1,nchar(study.data.excel$author)-6)
  
merged_pm <- merge(study.data.excel, Data.meta.pm.vs.nm.2[,c("Effect","Low","Up","Eff.c","Low.c","Up.c","eff.loh","low.loh","up.loh","Name")], by = "Name", all.x = TRUE)
merged_im <- merge(study.data.excel, Data.meta.im.vs.nm.2[,c("Effect","Low","Up","Eff.c","Low.c","Up.c","eff.loh","low.loh","up.loh","Name")], by = "Name", all.x = TRUE)
merged_im.as <- merge(study.data.excel, Data.meta.im.vs.nm.as.2[,c("Effect","Low","Up","Eff.c","Low.c","Up.c","eff.loh","low.loh","up.loh","Name")], by = "Name", all.x = TRUE)
merged_pmim <- merge(study.data.excel, Data.meta.pmim.vs.nm.2[,c("Effect","Low","Up","Eff.c","Low.c","Up.c","eff.loh","low.loh","up.loh","Name")], by = "Name", all.x = TRUE)
  
  
write.csv(merged_pm, "~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/meta_data_pm.csv")
write.csv(merged_im, "~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/meta_data_im.csv")
write.csv(merged_im.as, "~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/meta_data_imas.csv")
write.csv(merged_pmim, "~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/meta_data_pmim.csv")
  
  
  
  
  
  
  
  
  
  
  
    
#for display
Disp<-Data.meta.pm.vs.nm.2
for (i in c(1:6)){
  Disp[,i]<-round(Disp[,i],2)
}
Disp$Name[Disp$Name=="Goetzb"] <- "Goetz"
#Inv rank jpgs

# Log-transform the 'Low' and 'Up' columns in your dataframe
Disp$log_Effect <- log(Disp$Effect)
Disp$log_Low <- log(Disp$Low)
Disp$log_Up <- log(Disp$Up)

# Create a ggplot object: FIGURE 2
Disp$Inverse_Normal_Rank <- qnorm((rank(Disp$Effect) - 0.5) / length(Disp$Effect))
Disp$weights <- 1 / ((Disp$log_Up - Disp$log_Low) / (2 * qnorm(0.975)))^2

weighted_lm <- lm(log_Effect ~ Inverse_Normal_Rank, data = Disp, weights = weights)
Disp$predicted <- predict(weighted_lm)
# Create a ggplot object
p <- ggplot(Disp, aes(x = Inverse_Normal_Rank)) +
  geom_point(aes(y = log_Effect), color = "blue") +  # Add points for the effect
  geom_errorbar(aes(ymin = log_Low, ymax = log_Up), width = 0.1) +  # Add error bars
  geom_line(aes(y = predicted), linetype = "dashed", color = "red") +  # Add dashed line for predicted values
  labs(x = "Inverse Normal Rank Percentile", y = "Log(Relative Risk)") +  # Label axes
  theme_minimal() +  # Apply a minimal theme
  theme(axis.text.x = element_text(size = 14),  # Increase font size of x-axis labels
        axis.text.y = element_text(size = 14),axis.title = element_text(size = 14))  # Increase font size of y-axis labels

jpeg("/Users/maclehose/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/results/invnormplotpmvsnm.jpg",
     quality=100,width=600,height=500,units="px")
print(p)
dev.off()
  



jpeg("/Users/maclehose/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/results/change effects pm v nm .jpg",
     quality=100,width=600,height=500,units="px")
ggplot() +
  geom_point(data = Disp, aes(x = Name, y = Effect), position = position_nudge(x = 0.1), color = "red") +
  geom_errorbar(data = Disp, aes(x = Name, ymin = Low, ymax = Up), width = 0.2, position = position_nudge(x = 0.1), color = "red") +
  geom_point(data = Disp, aes(x = Name, y = Eff.c), position = position_nudge(x = -0.1)) +
  geom_errorbar(data = Disp, aes(x = Name, ymin = Low.c, ymax = Up.c), width = 0.2, position = position_nudge(x = -0.1)) +
  theme_light() +  # Set background plane to white
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increase font size on x-axis labels
        axis.text.y = element_text(size = 12),  # Set font size on y-axis labels
        axis.title = element_text(size = 14)) +  # Increase font size for axis titles
  scale_y_continuous(trans = 'log', breaks = c(0.1, 0.5, 1, 2, 5, 10)) +
  labs(y = "RR", x = NULL)  # Set y-axis label to "Log(RR)" and remove x-axis label

  dev.off()
  
  
#for display
Disp<-Data.meta.im.vs.nm.2
for (i in c(1:6)){
  Disp[,i]<-round(Disp[,i],2)
}
Disp$log_Effect <- log(Disp$Effect)
Disp$log_Low <- log(Disp$Low)
Disp$log_Up <- log(Disp$Up)
  
# Create a ggplot object
Disp$Inverse_Normal_Rank <- qnorm((rank(Disp$Effect) - 0.5) / length(Disp$Effect))
Disp$weights <- 1 / ((Disp$log_Up - Disp$log_Low) / (2 * qnorm(0.975))^2)
  
weighted_lm <- lm(log_Effect ~ Inverse_Normal_Rank, data = Disp, weights = weights)
Disp$predicted <- predict(weighted_lm)
# Create a ggplot object
p <- ggplot(Disp, aes(x = Inverse_Normal_Rank)) +
  geom_point(aes(y = log_Effect), color = "blue") +  # Add points for the effect
  geom_errorbar(aes(ymin = log_Low, ymax = log_Up), width = 0.1) +  # Add error bars
  geom_line(aes(y = predicted), linetype = "dashed", color = "red") +  # Add dashed line for predicted values
  labs(x = "Inverse Normal Rank Percentile", y = "Log(Relative Risk)") +  # Label axes
  theme_minimal() +  # Apply a minimal theme
  theme(axis.text.x = element_text(size = 12),  # Increase font size of x-axis labels
        axis.text.y = element_text(size = 12),axis.title = element_text(size = 14))  # Increase font size of y-axis labels

jpeg("/Users/maclehose/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/results/invnormplotimvsnm.jpg",
     quality=100,width=600,height=500,units="px")
print(p)
dev.off()
  
  
  
jpeg("/Users/maclehose/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/results/change effects im v nm .jpg",
     quality=100,width=800,height=600,units="px")
ggplot() +
  geom_point(data = Disp, aes(x = Name, y = Effect), position = position_nudge(x = 0.1), color = "red") +
  geom_errorbar(data = Disp, aes(x = Name, ymin = Low, ymax = Up), width = 0.2, position = position_nudge(x = 0.1), color = "red") +
  geom_point(data = Disp, aes(x = Name, y = Eff.c), position = position_nudge(x = -0.1)) +
  geom_errorbar(data = Disp, aes(x = Name, ymin = Low.c, ymax = Up.c), width = 0.2, position = position_nudge(x = -0.1)) +
  theme_light() +  # Set background plane to white
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increase font size on x-axis labels
        axis.text.y = element_text(size = 12),  # Set font size on y-axis labels
        axis.title = element_text(size = 14)) +  # Increase font size for axis titles
  scale_y_continuous(trans = 'log', breaks = c(0.1, 0.5, 1, 2, 5, 10)) +
  labs(y = "RR", x = NULL)  # Set y-axis label to "Log(RR)" and remove x-axis label

dev.off()
  
  
  #for display
  Disp<-Data.meta.im.vs.nm.as.2
  for (i in c(1:6)){
    Disp[,i]<-round(Disp[,i],2)
  }
  
  jpeg("/Users/maclehose/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/results/change effects im v nm as .jpg",
       quality=100,width=800,height=600,units="px")
  ggplot() +
    geom_point(data = Disp, aes(x = Name, y = Effect), position = position_nudge(x = 0.1), color = "red") +
    geom_errorbar(data = Disp, aes(x = Name, ymin = Low, ymax = Up), width = 0.2, position = position_nudge(x = 0.1), color = "red") +
    geom_point(data = Disp, aes(x = Name, y = Eff.c), position = position_nudge(x = -0.1)) +
    geom_errorbar(data = Disp, aes(x = Name, ymin = Low.c, ymax = Up.c), width = 0.2, position = position_nudge(x = -0.1)) +
    theme_light() +  # Set background plane to white
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increase font size on x-axis labels
          axis.text.y = element_text(size = 12),  # Set font size on y-axis labels
          axis.title = element_text(size = 14)) +  # Increase font size for axis titles
    scale_y_continuous(trans = 'log', breaks = c(0.1, 0.5, 1, 2, 5, 10)) +
    labs(y = "RR", x = NULL)  # Set y-axis label to "Log(RR)" and remove x-axis label
  
  dev.off()
  
  
  
  
  
  
