# CYP2D6 functions to load data
# 1. study data
# 2. validation data
  
  
##################
#Load Study Data for studies WITH summary leve data
extract.cyp.data<-function(idnum){
  study.data<-read_excel("~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/meta_data final.xlsx")
  study.data<-study.data[idnum,]
  if (study.data$genotype==1){
    vv.1<-study.data$VV_1
    vw.1<-study.data$VW_1
    ww.1<-study.data$WW_1
    vv.0<-study.data$VV_0
    vw.0<-study.data$VW_0
    ww.0<-study.data$WW_0
  }
  if (study.data$genotype==0){
    vv.1<-study.data$PM_1
    vw.1<-study.data$IM_1
    ww.1<-study.data$NM_1
    vv.0<-study.data$PM_0
    vw.0<-study.data$IM_0
    ww.0<-study.data$NM_0
  }
      
  if (study.data$european==1){
    region <-"Eur"
  }
  if (study.data$asian==1){
    region<-"As"
  }
  if (study.data$dnasource=="Tumor"){
    samp.type="t"
  }
  if (study.data$dnasource=="Non-neoplastic"){
    samp.type="n"
  }
  if(study.data$`Study type`=="cohort"){
    cohort=1
  }
  if(study.data$`Study type`=="cact"){
    cohort=0
  }
  author<-substr(study.data$author,1,nchar(study.data$author)-6)
  year<-study.data$`year published`
  c(vv.1,vw.1,ww.1,vv.0,vw.0,ww.0,region,samp.type,cohort,author,year)
}
    
##################
#Extract data from studies without summary level data
extract.cyp.data.2<-function(idnum){
  study.data<-read_excel("~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/meta_data final.xlsx")
  study.data<-study.data[idnum,]
  #abstract up to 4 effects from each: 2 contrasts each has crude and adjusted
  eff1<-study.data$`Point Est crude 1`
  up1<-study.data$`Upper CI crude 1`
  low1<-study.data$`Lower CI crude 1`
  eff1.type<-study.data$`ecode1`
  eff2<-study.data$`Point Est crude 2`
  up2<-study.data$`Upper CI crude 2`
  low2<-study.data$`Lower CI crude 2`
  eff2.type<-study.data$`ecode2`
  eff.adj.type<-"crude"
      
  if (is.na(eff1)==1){
    eff1<-study.data$`Point Est adj 1`
    up1<-study.data$`Upper CI adj 1`
    low1<-study.data$`Lower CI adj 1`
    eff.adj.type="adj"
  }
      
  if (is.na(eff2)==1){
    eff2<-study.data$`Point Est adj 2`
    up2<-study.data$`Upper CI adj 2`
    low2<-study.data$`Lower CI adj 2`
  }
      
  if (study.data$european==1){
    region <-"Eur"
  }
  if (study.data$asian==1){
    region<-"As"
  }
  if (study.data$dnasource=="Tumor"){
    samp.type="t"
  }
  if (study.data$dnasource=="Non-neoplastic"){
    samp.type="n"
  }
  if(study.data$`Study type`=="cohort"){
    cohort=1
  }
  if(study.data$`Study type`=="cact"){
    cohort=0
  }
  author<-substr(study.data$author,1,nchar(study.data$author)-6)
  year<-study.data$`year published`
  c(eff1,low1,up1,eff2,low2,up2,eff.adj.type,eff1.type,eff2.type,region,samp.type,cohort,author,year)
}
    
###################
#Load SNP data
#Read in data from Human Genome provided by Katherine Li, 11/2022
#European Data
load.gen.data<-function(region){
  eur.snp.data<-as.data.frame(read.csv("~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/3.1EA_long.csv"))
  eur.snp.data<-eur.snp.data[,c(1,c(1,4,5,6,7,9,10,12,13,14,16,20)+1)]
  colnames(eur.snp.data)<-c("ID","*2a","*41","*2b","*9","*3","*4","*8","*6","*29a","*29b","*17","*10")
  aa.snp.data<-as.data.frame(read.csv("~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/3.1Asian_long.csv"))
  aa.snp.data<-aa.snp.data[,c(1,c(1,4,5,6,7,9,10,12,13,14,16,20)+1)]
  colnames(aa.snp.data)<-c("ID","*2a","*41","*2b","*9","*3","*4","*8","*6","*29a","*29b","*17","*10")
  #number of obs
  n.people.eur<-dim(eur.snp.data)[1]
  n.people.aa<-dim(aa.snp.data)[1]  
  #Augment data so each person has two obs (1 for each copy of chrom 22)
  #  Data comes as "1|0" characters. want make it long dataset with 
  #  first row =1 and second =0
  snps.eur<-rbind(eur.snp.data,eur.snp.data)
  #assign snp 1/0 values for each rsid
  for (i in 1:n.people.eur){
    snps.eur[i,2:13]=substr(snps.eur[i,2:13],1,1)
    snps.eur[(i+n.people.eur),2:13]=substr(snps.eur[(i+n.people.eur),2:13],3,3)
  }
  snps.eur[,2:13]<-mcmapply(as.numeric,snps.eur[,2:13])
  snps.eur<-snps.eur[order(snps.eur$ID),]
  snps.aa<-rbind(aa.snp.data,aa.snp.data)
  #assign snp 1/0 values for each rsid
  for (i in 1:n.people.aa){
    snps.aa[i,2:13]=substr(snps.aa[i,2:13],1,1)
    snps.aa[(i+n.people.aa),2:13]=substr(snps.aa[(i+n.people.aa),2:13],3,3)
  }
  snps.aa[,2:13]<-mcmapply(as.numeric,snps.aa[,2:13])
  snps.aa<-snps.aa[order(snps.aa$ID),]
  if (region=="Eur"){
    snp.out<-snps.eur
  }
  if (region=="As"){
    snp.out<-snps.aa
  }
  snp.out
}
    
    
#Set up data for meta analysis
#compute variances
meta.data.setup<-function(dataname){
  #remove empty rows
  dataname<-dataname[dataname[,1]!=0,]
  #convert to numeric
  i<-c(1:7,18:20)
  dataname[ , i] <- apply(dataname[ , i], 2,function(x) as.numeric(as.character(x)))
  #calc
  dataname$v<-((log(dataname$Up)-log(dataname$Low))/(2*1.96))^2
  dataname$v.c<-((log(dataname$Up.c)-log(dataname$Low.c))/(2*1.96))^2
  dataname$v.loh<-((log(dataname$up.loh)-log(dataname$low.loh))/(2*1.96))^2
  dataname$bias.eff<-log(dataname$Effect)-log(dataname$Eff.c)
  #approximate increase in variance from crude to adjusted 
  dataname$bias.v<-dataname$v-dataname$v.c
  dataname$bias.v[dataname$bias.v<0]<-0
  dataname$bias.eff.loh<-log(dataname$eff.loh)-log(dataname$Eff.c)
  dataname$bias.v.loh<-dataname$v.loh-dataname$v.c
  dataname$bias.v.loh[dataname$bias.v.loh<0]<-0
  dataname
}
    
    
#introduce bias adjustment for studies that have only effect estimates, not summary level tables
pt.est.bias.adj<-function(dataname){
  #var weighted avg difference between crude and bias adjusted estimate
  avg.miss.bias<-sum(dataname$bias.eff/dataname$v,na.rm=TRUE)/sum(1/dataname$v,na.rm=TRUE)
  #avg increase in variance between crude and bias adjsuted estiamte
  avg.miss.v<-sum(dataname$bias.v/dataname$v,na.rm=TRUE)/sum(1/dataname$v,na.rm=TRUE)
  #row with no adjusted effect
  forwhom<-is.na(dataname$Effect)
  #add in bias
  dataname$Effect[forwhom]<-exp(log(dataname$Eff.c[forwhom])+avg.miss.bias)
  #add in variance
  dataname$v[forwhom]<-dataname$v.c[forwhom]+avg.miss.v
  #compute CIs
  dataname$Low[forwhom]<-exp(log(dataname$Effect[forwhom])-1.96*sqrt(dataname$v[forwhom]))
  dataname$Up[forwhom]<-exp(log(dataname$Effect[forwhom])+1.96*sqrt(dataname$v[forwhom]))
      
  #repeat for loh
  #avg difference between crude and bias adjusted estimate
  avg.miss.bias.loh<-sum(dataname$bias.eff.loh[dataname$Samp.Type=="t"]/dataname$v.loh[dataname$Samp.Type=="t"],na.rm=TRUE)/sum(1/dataname$v.loh[dataname$Samp.Type=="t"],na.rm=TRUE)
  #avg increase in variance between crude and bias adjsuted estiamte
  avg.miss.v.loh<-sum(dataname$bias.v.loh[dataname$Samp.Type=="t"]/dataname$v.loh[dataname$Samp.Type=="t"],na.rm=TRUE)/sum(1/dataname$v.loh[dataname$Samp.Type=="t"],na.rm=TRUE)
  #row with no adjusted effect
  forwhom<-is.na(dataname$eff.loh)
  #add in bias
  dataname$eff.loh[forwhom]<-exp(log(dataname$Eff.c[forwhom])+avg.miss.bias.loh)
  #add in variance
  dataname$v.loh[forwhom]<-dataname$v.c[forwhom]+avg.miss.v.loh
  #compute CIs
  dataname$low.loh[forwhom]<-exp(log(dataname$eff.loh[forwhom])-1.96*sqrt(dataname$v.loh[forwhom]))
  dataname$up.loh[forwhom]<-exp(log(dataname$eff.loh[forwhom])+1.96*sqrt(dataname$v.loh[forwhom]))
  dataname
}
