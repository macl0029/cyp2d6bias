    
#function to assign activity scores to SNP data (long form)
assign.act.score<-function(X){
  t1<-X
  t1.n.snp<-dim(t1)[2]-1
  #create a column for wild type (lack of any other *'s)
  t1$wt<-rowSums(t1[,2:(t1.n.snp+1)])==0
  #Give activity score to each snp
  t1$`*2a`<-ifelse(t1$`*2a`==1,1,NA)
  t1$`*2b`<-ifelse(t1$`*2b`==1,1,NA)
  t1$`*3`<-ifelse(t1$`*3`==1,0,NA)
  t1$`*4`<-ifelse(t1$`*4`==1,0,NA)
  t1$`*6`<-ifelse(t1$`*6`==1,0,NA)
  t1$`*8`<-ifelse(t1$`*8`==1,0,NA)
  t1$`*9`<-ifelse(t1$`*9`==1,0.25,NA)
  t1$`*10`<-ifelse(t1$`*10`==1,0.25,NA)
  t1$`*17`<-ifelse(t1$`*17`==1,0.5,NA)
  t1$`*29a`<-ifelse(t1$`*29a`==1,0.5,NA)
  t1$`*29b`<-ifelse(t1$`*29b`==1,0.5,NA)
  t1$`*41`<-ifelse(t1$`*41`==1,0.25,NA)
  #Activity score is min of snp-specific scores
  t1$as<-rowMins(as.matrix(t1[,2:(t1.n.snp+1)]),na.rm = TRUE)
  #reassign wild type's to AS=1
  t1$as[t1$wt=="TRUE"]=1
  X$AS<-t1$as
  X
}
assign.act.score.regan<-function(X){
  t1<-X
  col.get<-which(names(t1)%in%c("*3","*4","*6","*41"))
  #create a column for wild type (lack of any other *'s)
  t1$wt<-rowSums(t1[,col.get])==0
  #Give activity score to each snp
  t1$`*3`<-ifelse(t1$`*3`==1,0,NA)
  t1$`*4`<-ifelse(t1$`*4`==1,0,NA)
  t1$`*6`<-ifelse(t1$`*6`==1,0,NA)
  t1$`*41`<-ifelse(t1$`*41`==1,0.5,NA)
  #Activity score is min of snp-specific scores
  t1$as<-rowMins(as.matrix(t1[,col.get]),na.rm = TRUE)
  #reassign wild type's to AS=1
  t1$as[t1$wt=="TRUE"]=1
  X$AS.alt<-t1$as
  X
}
assign.act.score.caudle<-function(X){
  t1<-X
  col.get<-which(names(t1)%in%c("*4","*10","*41"))
  #create a column for wild type (lack of any other *'s)
  t1$wt<-rowSums(t1[,col.get])==0
  #Give activity score to each snp
  t1$`*4`<-ifelse(t1$`*4`==1,0,NA)
  t1$`*10`<-ifelse(t1$`*10`==1,.25,NA)
  t1$`*41`<-ifelse(t1$`*41`==1,0.5,NA)
  #Activity score is min of snp-specific scores
  t1$as<-rowMins(as.matrix(t1[,col.get]),na.rm = TRUE)
  #reassign wild type's to AS=1
  t1$as[t1$wt=="TRUE"]=1
  X$AS.alt<-t1$as
  X
}
assign.act.score.goetz<-function(X){
  t1<-X[,-14]
  t1.n.snp<-dim(t1)[2]-1
  #create a column for wild type (lack of any other *'s)
  t1$wt<-rowSums(t1[,2:(t1.n.snp+1)])==0
  #Give activity score to each snp
  t1$`*2a`<-ifelse(t1$`*2a`==1,1,NA)
  t1$`*2b`<-ifelse(t1$`*2b`==1,1,NA)
  t1$`*3`<-ifelse(t1$`*3`==1,0,NA)
  t1$`*4`<-ifelse(t1$`*4`==1,0,NA)
  t1$`*6`<-ifelse(t1$`*6`==1,0,NA)
  t1$`*8`<-ifelse(t1$`*8`==1,0,NA)
  t1$`*9`<-ifelse(t1$`*9`==1,0.5,NA)
  t1$`*10`<-ifelse(t1$`*10`==1,0.5,NA)
  t1$`*17`<-ifelse(t1$`*17`==1,0.5,NA)
  t1$`*29a`<-ifelse(t1$`*29a`==1,0.5,NA)
  t1$`*29b`<-ifelse(t1$`*29b`==1,0.5,NA)
  t1$`*41`<-ifelse(t1$`*41`==1,0.5,NA)
  #Activity score is min of snp-specific scores
  t1$as<-rowMins(as.matrix(t1[,2:(t1.n.snp+1)]),na.rm = TRUE)
  #reassign wild type's to AS=1
  t1$as[t1$wt=="TRUE"]=1
  X$AS.alt<-t1$as
  X
}
assign.act.score.markkula<-function(X){
  t1<-X
  col.get<-which(names(t1)%in%c("*3","*4","*6","*10","*41"))
  #create a column for wild type (lack of any other *'s)
  t1$wt<-rowSums(t1[,col.get])==0
  #Give activity score to each snp
  t1$`*3`<-ifelse(t1$`*3`==1,0,NA)
  t1$`*4`<-ifelse(t1$`*4`==1,0,NA)
  t1$`*6`<-ifelse(t1$`*6`==1,0,NA)
  t1$`*10`<-ifelse(t1$`*10`==1,0.5,NA)
  t1$`*41`<-ifelse(t1$`*41`==1,0.5,NA)
  #Activity score is min of snp-specific scores
  t1$as<-rowMins(as.matrix(t1[,col.get]),na.rm = TRUE)
  #reassign wild type's to AS=1
  t1$as[t1$wt=="TRUE"]=1
  X$AS.alt<-t1$as
  X
}
assign.act.score.schroth<-function(X){
  t1<-X
  col.get<-which(names(t1)%in%c("*4","*10","*41"))
  #create a column for wild type (lack of any other *'s)
  t1$wt<-rowSums(t1[,c(col.get)])==0
  #Give activity score to each snp
  t1$`*4`<-ifelse(t1$`*4`==1,0,NA)
  t1$`*10`<-ifelse(t1$`*10`==1,0.5,NA)
  t1$`*41`<-ifelse(t1$`*41`==1,0.5,NA)
  #Activity score is min of snp-specific scores
  t1$as<-rowMins(as.matrix(t1[,col.get]),na.rm = TRUE)
  #reassign wild type's to AS=1
  t1$as[t1$wt=="TRUE"]=1
  X$AS.alt<-t1$as
  X
}
assign.act.score.cham<-function(X){
  t1<-X
  col.get<-which(names(t1)%in%c("*3","*4","*6","*8","*9","*10","*17","*29a","*41"))
  #create a column for wild type (lack of any other *'s)
  t1$wt<-rowSums(t1[,c(col.get)])==0
  #Give activity score to each snp
  t1$`*3`<-ifelse(t1$`*3`==1,0,NA)
  t1$`*4`<-ifelse(t1$`*4`==1,0,NA)
  t1$`*6`<-ifelse(t1$`*6`==1,0,NA)
  t1$`*8`<-ifelse(t1$`*8`==1,0,NA)
  t1$`*9`<-ifelse(t1$`*9`==1,.5,NA)
  t1$`*10`<-ifelse(t1$`*10`==1,0.5,NA)
  t1$`*17`<-ifelse(t1$`*17`==1,0.5,NA)
  t1$`*29a`<-ifelse(t1$`*29a`==1,.5,NA)
  t1$`*41`<-ifelse(t1$`*41`==1,0.5,NA)
  #Activity score is min of snp-specific scores
  t1$as<-rowMins(as.matrix(t1[,col.get]),na.rm = TRUE)
  #reassign wild type's to AS=1
  t1$as[t1$wt=="TRUE"]=1
  X$AS.alt<-t1$as
  X
}
assign.act.score.mwinyi<-function(X){
  t1<-X
  col.get<-which(names(t1)%in%c("*3","*4","*6","*10","*41"))
  #create a column for wild type (lack of any other *'s)
  t1$wt<-rowSums(t1[,col.get])==0
  #Give activity score to each snp
  t1$`*3`<-ifelse(t1$`*3`==1,0,NA)
  t1$`*4`<-ifelse(t1$`*4`==1,0,NA)
  t1$`*6`<-ifelse(t1$`*6`==1,0,NA)
  t1$`*10`<-ifelse(t1$`*10`==1,0.5,NA)
  t1$`*41`<-ifelse(t1$`*41`==1,0.5,NA)
  #Activity score is min of snp-specific scores
  t1$as<-rowMins(as.matrix(t1[,col.get]),na.rm = TRUE)
  #reassign wild type's to AS=1
  t1$as[t1$wt=="TRUE"]=1
  X$AS.alt<-t1$as
  X
}
assign.act.score.he<-function(X){
  t1<-X
  col.get<-which(names(t1)%in%c("*3","*4","*6","*9","*10","*17","*29a","*41"))
  #create a column for wild type (lack of any other *'s)
  t1$wt<-rowSums(t1[,col.get])==0
  #Give activity score to each snp
  t1$`*3`<-ifelse(t1$`*3`==1,0,NA)
  t1$`*4`<-ifelse(t1$`*4`==1,0,NA)
  t1$`*6`<-ifelse(t1$`*6`==1,0,NA)
  t1$`*9`<-ifelse(t1$`*9`==1,.5,NA)
  t1$`*10`<-ifelse(t1$`*10`==1,.5,NA)
  t1$`*17`<-ifelse(t1$`*17`==1,.5,NA)
  t1$`*29a`<-ifelse(t1$`*29a`==1,.5,NA)
  t1$`*41`<-ifelse(t1$`*41`==1,0.5,NA)
  #Activity score is min of snp-specific scores
  t1$as<-rowMins(as.matrix(t1[,col.get]),na.rm = TRUE)
  #reassign wild type's to AS=1
  t1$as[t1$wt=="TRUE"]=1
  X$AS.alt<-t1$as
  X
}
assign.act.score.park<-function(X){
  t1<-X
  col.get<-which(names(t1)%in%c("*3","*4","*10","*41"))
  #create a column for wild type (lack of any other *'s)
  t1$wt<-rowSums(t1[,col.get])==0
  #Give activity score to each snp
  t1$`*3`<-ifelse(t1$`*3`==1,0,NA)
  t1$`*4`<-ifelse(t1$`*4`==1,0,NA)
  t1$`*10`<-ifelse(t1$`*10`==1,0.5,NA)
  t1$`*41`<-ifelse(t1$`*41`==1,0.5,NA)
  #Activity score is min of snp-specific scores
  t1$as<-rowMins(as.matrix(t1[,col.get]),na.rm = TRUE)
  #reassign wild type's to AS=1
  t1$as[t1$wt=="TRUE"]=1
  X$AS.alt<-t1$as
  X
}
assign.act.score.melo<-function(X){
  t1<-X
  col.get<-which(names(t1)%in%c("*4","*10","*17"))
  #create a column for wild type (lack of any other *'s)
  t1$wt<-rowSums(t1[,col.get])==0
  #Give activity score to each snp
  t1$`*4`<-ifelse(t1$`*4`==1,0,NA)
  t1$`*10`<-ifelse(t1$`*10`==1,0.5,NA)
  t1$`*17`<-ifelse(t1$`*17`==1,0.5,NA)
  #Activity score is min of snp-specific scores
  t1$as<-rowMins(as.matrix(t1[,col.get]),na.rm = TRUE)
  #reassign wild type's to AS=1
  t1$as[t1$wt=="TRUE"]=1
  X$AS.alt<-t1$as
  X
}
assign.act.score.kiyotani<-function(X){
  t1<-X
  col.get<-which(names(t1)%in%c("*4","*10","*41"))
  #create a column for wild type (lack of any other *'s)
  t1$wt<-rowSums(t1[,col.get])==0
  #Give activity score to each snp
  t1$`*4`<-ifelse(t1$`*4`==1,0,NA)
  t1$`*10`<-ifelse(t1$`*10`==1,0,NA)
  t1$`*41`<-ifelse(t1$`*41`==1,0,NA)
  #Activity score is min of snp-specific scores
  t1$as<-rowMins(as.matrix(t1[,col.get]),na.rm = TRUE)
  #reassign wild type's to AS=1
  t1$as[t1$wt=="TRUE"]=1
  X$AS.alt<-t1$as
  X
}
assign.act.score.wrang<-function(X){
  t1<-X
  col.get<-which(names(t1)%in%c("*3","*4","*10","*17","*29a"))
  #create a column for wild type (lack of any other *'s)
  t1$wt<-rowSums(t1[,c(col.get)])==0
  #Give activity score to each snp
  t1$`*3`<-ifelse(t1$`*3`==1,0,NA)
  t1$`*4`<-ifelse(t1$`*4`==1,0,NA)
  t1$`*10`<-ifelse(t1$`*10`==1,0.25,NA)
  t1$`*17`<-ifelse(t1$`*17`==1,0.5,NA)
  t1$`*29a`<-ifelse(t1$`*29a`==1,.5,NA)
  #Activity score is min of snp-specific scores
  t1$as<-rowMins(as.matrix(t1[,col.get]),na.rm = TRUE)
  #reassign wild type's to AS=1
  t1$as[t1$wt=="TRUE"]=1
  X$AS.alt<-t1$as
  X
}
assign.act.score.teh<-function(X){
  t1<-X
  col.get<-which(names(t1)%in%c("*4","*10"))
  #create a column for wild type (lack of any other *'s)
  t1$wt<-rowSums(t1[,c(col.get)])==0
  #Give activity score to each snp
  t1$`*4`<-ifelse(t1$`*4`==1,0,NA)
  t1$`*10`<-ifelse(t1$`*10`==1,0.5,NA)
  #Activity score is min of snp-specific scores
  t1$as<-rowMins(as.matrix(t1[,col.get]),na.rm = TRUE)
  #reassign wild type's to AS=1
  t1$as[t1$wt=="TRUE"]=1
  X$AS.alt<-t1$as
  X
}

#Functions to reshape to wide format, create total AS, and genotype =w/w, w/v, v/v for *4
val.phenotype.star4<-function(X){
  long.ptype<-X
  len<-dim(long.ptype)[1]
  #create variable for first vs second chromosome
  long.ptype$chr<-rep(c(1,2),len/2)
  #reshape long to wide
  wide.ptype<-reshape(long.ptype,idvar="ID",timevar="chr",direction="wide")
  #Create 'exposure': 2=W/W; 1= *4/W (w/v); 0= *4/*4 (v/v)
  wide.ptype$Ptype.exp<-NA
  wide.ptype$Ptype.exp[wide.ptype$`*4.1`==1&wide.ptype$`*4.2`==1]=0
  wide.ptype$Ptype.exp[wide.ptype$`*4.1`==0&wide.ptype$`*4.2`==1]=1
  wide.ptype$Ptype.exp[wide.ptype$`*4.1`==1&wide.ptype$`*4.2`==0]=1
  wide.ptype$Ptype.exp[wide.ptype$`*4.1`==0&wide.ptype$`*4.2`==0]=2
  #Create Activity score
  wide.ptype$AS<-wide.ptype$AS.1+wide.ptype$AS.2
  #Create Ptype: 0=PM; 1=IM; 2=NM
  #from Caudle 2020
  wide.ptype$Ptype<-NA
  wide.ptype$Ptype[wide.ptype$AS==0]=0
  wide.ptype$Ptype[wide.ptype$AS>0&wide.ptype$AS<1.25]=1
  wide.ptype$Ptype[wide.ptype$AS>=1.25&wide.ptype$AS<=2.25]=2
  Out<-wide.ptype[,c(28:30)]
  Out
}
val.phenotype.star10<-function(X){
  long.ptype<-X
  len<-dim(long.ptype)[1]
  #create variable for first vs second chromosome
  long.ptype$chr<-rep(c(1,2),len/2)
  #reshape long to wide
  wide.ptype<-reshape(long.ptype,idvar="ID",timevar="chr",direction="wide")
  #Create 'exposure': 2=W/W; 1= *10/W (w/v); 0= *10/*10 (v/v)
  wide.ptype$Ptype.exp<-NA
  wide.ptype$Ptype.exp[wide.ptype$`*10.1`==1&wide.ptype$`*10.2`==1]=0
  wide.ptype$Ptype.exp[wide.ptype$`*10.1`==0&wide.ptype$`*10.2`==1]=1
  wide.ptype$Ptype.exp[wide.ptype$`*10.1`==1&wide.ptype$`*10.2`==0]=1
  wide.ptype$Ptype.exp[wide.ptype$`*10.1`==0&wide.ptype$`*10.2`==0]=2
  #Create Activity score
  wide.ptype$AS<-wide.ptype$AS.1+wide.ptype$AS.2
  #Create Ptype: 0=PM; 1=IM; 2=NM
  #from Caudle 2020
  wide.ptype$Ptype<-NA
  wide.ptype$Ptype[wide.ptype$AS==0]=0
  wide.ptype$Ptype[wide.ptype$AS>0&wide.ptype$AS<1.25]=1
  wide.ptype$Ptype[wide.ptype$AS>=1.25&wide.ptype$AS<=2.25]=2
  Out<-wide.ptype[,c(28:30)]
  Out
}
val.phenotype.goetz<-function(X){
  long.ptype<-X
  len<-dim(long.ptype)[1]
  #create variable for first vs second chromosome
  long.ptype$chr<-rep(c(1,2),len/2)
  #reshape long to wide
  wide.ptype<-reshape(long.ptype,idvar="ID",timevar="chr",direction="wide")
  
  #Create Activity score
  wide.ptype$AS<-wide.ptype$AS.1+wide.ptype$AS.2
  #Create Ptype: 0=PM; 1=IM; 2=NM
  wide.ptype$Ptype<-NA
  wide.ptype$Ptype[wide.ptype$AS==0]=0
  wide.ptype$Ptype[wide.ptype$AS>0&wide.ptype$AS<1.25]=1
  wide.ptype$Ptype[wide.ptype$AS>=1.25&wide.ptype$AS<=2.25]=2
  
  #Create Activity score - ALT
  wide.ptype$AS.alt<-wide.ptype$AS.alt.1+wide.ptype$AS.alt.2
  #Create Ptype: 0=PM; 1=IM; 2=NM
  wide.ptype$Ptype.exp<-NA
  wide.ptype$Ptype.exp[wide.ptype$AS.alt==0]=0
  wide.ptype$Ptype.exp[wide.ptype$AS.alt>0&wide.ptype$AS.alt<1.5]=1
  wide.ptype$Ptype.exp[wide.ptype$AS.alt>=1.5]=2
    
  Out<-wide.ptype[,c(31,33)]
    
  Out
}
val.phenotype.alt<-function(X){
  long.ptype<-X
  len<-dim(long.ptype)[1]
  #create variable for first vs second chromosome
  long.ptype$chr<-rep(c(1,2),len/2)
  #reshape long to wide
  wide.ptype<-reshape(long.ptype,idvar="ID",timevar="chr",direction="wide")
  
  #Create Activity score
  wide.ptype$AS<-wide.ptype$AS.1+wide.ptype$AS.2
  #Create Ptype: 0=PM; 1=IM; 2=NM
  wide.ptype$Ptype<-NA
  wide.ptype$Ptype[wide.ptype$AS==0]=0
  wide.ptype$Ptype[wide.ptype$AS>0&wide.ptype$AS<1.25]=1
  wide.ptype$Ptype[wide.ptype$AS>=1.25&wide.ptype$AS<=2.25]=2
    
  #Create Activity score - ALT
  wide.ptype$AS.alt<-wide.ptype$AS.alt.1+wide.ptype$AS.alt.2
  #Create Ptype: 0=PM; 1=IM; 2=NM
  wide.ptype$Ptype.exp<-NA
  wide.ptype$Ptype.exp[wide.ptype$AS.alt==0]=0
  wide.ptype$Ptype.exp[wide.ptype$AS.alt>0&wide.ptype$AS.alt<2]=1
  wide.ptype$Ptype.exp[wide.ptype$AS.alt==2]=2
      
  Out<-wide.ptype[,c(31,33)]
  Out
}
val.phenotype.ism<-function(X){
  long.ptype<-X
  len<-dim(long.ptype)[1]
  #create variable for first vs second chromosome
  long.ptype$chr<-rep(c(1,2),len/2)
  #reshape long to wide
  wide.ptype<-reshape(long.ptype,idvar="ID",timevar="chr",direction="wide")
  
  #Create Activity score
  wide.ptype$AS<-wide.ptype$AS.1+wide.ptype$AS.2
  #Create Ptype: 0=PM; 1=IM; 2=NM
  wide.ptype$Ptype<-NA
  wide.ptype$Ptype[wide.ptype$AS==0]=0
  wide.ptype$Ptype[wide.ptype$AS>0&wide.ptype$AS<1.25]=1
  wide.ptype$Ptype[wide.ptype$AS>=1.25&wide.ptype$AS<=2.25]=2
  
  #Create Activity score - ALT
  wide.ptype$AS.alt<-wide.ptype$AS.alt.1+wide.ptype$AS.alt.2
  #Create Ptype: 0=PM; 1=IM; 2=NM
  wide.ptype$Ptype.exp<-NA
  wide.ptype$Ptype.exp[wide.ptype$AS.alt==0]=0
  wide.ptype$Ptype.exp[wide.ptype$AS.alt>0&wide.ptype$AS.alt<=1]=1
  wide.ptype$Ptype.exp[wide.ptype$AS.alt>1]=2
    
  Out<-wide.ptype[,c(31,33)]
  Out
}
val.phenotype.he<-function(X){
  long.ptype<-X
  len<-dim(long.ptype)[1]
  #create variable for first vs second chromosome
  long.ptype$chr<-rep(c(1,2),len/2)
  #reshape long to wide
  wide.ptype<-reshape(long.ptype,idvar="ID",timevar="chr",direction="wide")
  
  #Create Activity score
  wide.ptype$AS<-wide.ptype$AS.1+wide.ptype$AS.2
  #Create Ptype: 0=PM; 1=IM; 2=NM
  wide.ptype$Ptype<-NA
  wide.ptype$Ptype[wide.ptype$AS==0]=0
  wide.ptype$Ptype[wide.ptype$AS>0&wide.ptype$AS<1.25]=1
  wide.ptype$Ptype[wide.ptype$AS>=1.25&wide.ptype$AS<=2.25]=2
    
  #Create Activity score - ALT
  wide.ptype$AS.alt<-wide.ptype$AS.alt.1+wide.ptype$AS.alt.2
  #Create Ptype: 0=PM; 1=IM; 2=NM
  wide.ptype$Ptype.exp<-NA
  wide.ptype$Ptype.exp[wide.ptype$AS.alt==0]=0
  wide.ptype$Ptype.exp[wide.ptype$AS.alt>0&wide.ptype$AS.alt<1.5]=1
  wide.ptype$Ptype.exp[wide.ptype$AS.alt>=1.5]=2
    
  Out<-wide.ptype[,c(31,33)]
  Out
}
val.phenotype.park<-function(X){
  long.ptype<-X
  len<-dim(long.ptype)[1]
  #create variable for first vs second chromosome
  long.ptype$chr<-rep(c(1,2),len/2)
  #reshape long to wide
  wide.ptype<-reshape(long.ptype,idvar="ID",timevar="chr",direction="wide")
  
  #Create Activity score
  wide.ptype$AS<-wide.ptype$AS.1+wide.ptype$AS.2
  #Create Ptype: 0=PM; 1=IM; 2=NM
  wide.ptype$Ptype<-NA
  wide.ptype$Ptype[wide.ptype$AS==0]=0
  wide.ptype$Ptype[wide.ptype$AS>0&wide.ptype$AS<1.25]=1
  wide.ptype$Ptype[wide.ptype$AS>=1.25&wide.ptype$AS<=2.25]=2
  
  #Create Activity score - ALT
  wide.ptype$AS.alt<-wide.ptype$AS.alt.1+wide.ptype$AS.alt.2
  #Create Ptype: 0=PM; 1=IM; 2=NM
  wide.ptype$Ptype.exp<-NA
  wide.ptype$Ptype.exp[wide.ptype$AS.alt<1]=0
  wide.ptype$Ptype.exp[wide.ptype$AS.alt>=1&wide.ptype$AS.alt<2]=1
  wide.ptype$Ptype.exp[wide.ptype$AS.alt>=2]=2
  
  Out<-wide.ptype[,c(31,33)]
  Out
}
val.phenotype.schroth<-function(X){
  long.ptype<-X
  len<-dim(long.ptype)[1]
  #create variable for first vs second chromosome
  long.ptype$chr<-rep(c(1,2),len/2)
  #reshape long to wide
  wide.ptype<-reshape(long.ptype,idvar="ID",timevar="chr",direction="wide")
  
  #Create Activity score
  wide.ptype$AS<-wide.ptype$AS.1+wide.ptype$AS.2
  #Create Ptype: 0=PM; 1=IM; 2=NM
  wide.ptype$Ptype<-NA
  wide.ptype$Ptype[wide.ptype$AS==0]=0
  wide.ptype$Ptype[wide.ptype$AS>0&wide.ptype$AS<1.25]=1
  wide.ptype$Ptype[wide.ptype$AS>=1.25&wide.ptype$AS<=2.25]=2
  
  #Create Activity score - ALT
  wide.ptype$AS.alt<-wide.ptype$AS.alt.1+wide.ptype$AS.alt.2
  #Create Ptype: 0=PM; 1=IM; 2=NM
  wide.ptype$Ptype.exp<-NA
  wide.ptype$Ptype.exp[wide.ptype$AS.alt==0]=0
  wide.ptype$Ptype.exp[wide.ptype$AS.alt>0&wide.ptype$AS.alt<1.5]=1
  wide.ptype$Ptype.exp[wide.ptype$AS.alt>=1.5]=2
    
  Out<-wide.ptype[,c(31,33)]
  Out
} 
val.phenotype.wrang<-function(X){
  long.ptype<-X
  len<-dim(long.ptype)[1]
  #create variable for first vs second chromosome
  long.ptype$chr<-rep(c(1,2),len/2)
  #reshape long to wide
  wide.ptype<-reshape(long.ptype,idvar="ID",timevar="chr",direction="wide")
  
  #Create Activity score
  wide.ptype$AS<-wide.ptype$AS.1+wide.ptype$AS.2
  #Create Ptype: 0=PM; 1=IM; 2=NM
  wide.ptype$Ptype<-NA
  wide.ptype$Ptype[wide.ptype$AS==0]=0
  wide.ptype$Ptype[wide.ptype$AS>0&wide.ptype$AS<1.25]=1
  wide.ptype$Ptype[wide.ptype$AS>=1.25&wide.ptype$AS<=2.25]=2
  
  #Create Activity score - ALT
  wide.ptype$AS.alt<-wide.ptype$AS.alt.1+wide.ptype$AS.alt.2
  #Create Ptype: 0=PM; 1=IM; 2=NM
  wide.ptype$Ptype.exp<-NA
  wide.ptype$Ptype.exp[wide.ptype$AS.alt==0]=0
  wide.ptype$Ptype.exp[wide.ptype$AS.alt>0&wide.ptype$AS.alt<1.25]=1
  wide.ptype$Ptype.exp[wide.ptype$AS.alt>=1.25]=2
  
  Out<-wide.ptype[,c(31,33)]
  Out
}
val.phenotype.teh<-function(X){
  long.ptype<-X
  len<-dim(long.ptype)[1]
  #create variable for first vs second chromosome
  long.ptype$chr<-rep(c(1,2),len/2)
  #reshape long to wide
  wide.ptype<-reshape(long.ptype,idvar="ID",timevar="chr",direction="wide")
  
  #Create Activity score
  wide.ptype$AS<-wide.ptype$AS.1+wide.ptype$AS.2
  #Create Ptype: 0=PM; 1=IM; 2=NM
  wide.ptype$Ptype<-NA
  wide.ptype$Ptype[wide.ptype$AS==0]=0
  wide.ptype$Ptype[wide.ptype$AS>0&wide.ptype$AS<1.25]=1
  wide.ptype$Ptype[wide.ptype$AS>=1.25&wide.ptype$AS<=2.25]=2
  
  #Create Activity score - ALT
  wide.ptype$AS.alt<-wide.ptype$AS.alt.1+wide.ptype$AS.alt.2
  #Create Ptype: 0=PM; 1=IM; 2=NM
  wide.ptype$Ptype.exp<-NA
  wide.ptype$Ptype.exp[wide.ptype$AS.alt>=1.5]=2
  wide.ptype$Ptype.exp[wide.ptype$AS.alt<1.5]=1

  Out<-wide.ptype[,c(31,33)]
  Out
} 
