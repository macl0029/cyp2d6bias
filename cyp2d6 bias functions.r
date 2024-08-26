###################
#This function reads in study and validation data and formats it up for analysis
bias.setup<-function(Study.id){
  #Load Cyp2D6 data from chosen study  
  cyp.data<-extract.cyp.data(Study.id)
  vv.1 <- as.numeric(cyp.data[1])  #cases, no functional allele
  vw.1 <- as.numeric(cyp.data[2]) #cases, one functional allele
  ww.1 <- as.numeric(cyp.data[3]) #cases, two functional alleles
  vv.0 <- as.numeric(cyp.data[4])  #controls, no functional allele
  vw.0 <- as.numeric(cyp.data[5]) #controls, one functional allele
  ww.0 <- as.numeric(cyp.data[6]) #controls, no functional allele 
  region <- cyp.data[7]
  samp.type<-cyp.data[8]
  cohort=cyp.data[9]
  author=cyp.data[10]
  year=cyp.data[11]
      
    
  #load genotype data for imputation data for given region
  val.data<-load.gen.data(region)
  #Add activity score to validation data
  val.data<-assign.act.score(val.data)
    
  #functions to assign an activity score for the algorithm used in article:
  #Star 4 articles
  if (author=="Lash"|author=="Bijl"|author=="Dezentje"|author=="Abraham"|author=="Wegman"){
    star.phenotype.data<-val.phenotype.star4(val.data)
    #Create bias parameter matrix for SE/SP for validation data
    # 3 x 3 matrix. i=row; j=column
    # i,j cell is Sum(misclass genotype=j , true phenotype=i)
      sesp.prior.mat<-table(star.phenotype.data$Ptype,star.phenotype.data$Ptype.exp)
  }
  #star 10 articles
  if (author=="Sirachainan"){
    star.phenotype.data<-val.phenotype.star10(val.data)
    #Create bias parameter matrix for SE/SP for validation data
    # 3 x 3 matrix. i=row; j=column
    # i,j cell is Sum(misclass genotype=j , true phenotype=i)
    sesp.prior.mat<-table(star.phenotype.data$Ptype,star.phenotype.data$Ptype.exp)
    #no 0 genotype so add 0 row
    sesp.prior.mat<-rbind(c(0,0,0),sesp.prior.mat)
  }
  #phenotype algorithms    
  if (author=="Argalacsova"){
    val.data<-assign.act.score.goetz(val.data)
    star.phenotype.data<-val.phenotype.goetz(val.data)
    #Create bias parameter matrix for SE/SP for validation data
    # 3 x 3 matrix. i=row; j=column
    # i,j cell is Sum(genotype=j | phenotype=i)
    sesp.prior.mat<-table(star.phenotype.data$Ptype,star.phenotype.data$Ptype.exp)
  }
  if (author=="Regan"){
    val.data<-assign.act.score.regan(val.data)
    star.phenotype.data<-val.phenotype.alt(val.data)
    #Create bias parameter matrix for SE/SP for validation data
    # 3 x 3 matrix. i=row; j=column
    # i,j cell is Sum(genotype=j | phenotype=i)
    sesp.prior.mat<-table(star.phenotype.data$Ptype,star.phenotype.data$Ptype.exp)
  }
  if (author=="Markkula"){
    val.data<-assign.act.score.markkula(val.data)
    star.phenotype.data<-val.phenotype.alt(val.data)
    #Create bias parameter matrix for SE/SP for validation data
    # 3 x 3 matrix. i=row; j=column
    # i,j cell is Sum(genotype=j | phenotype=i)
    sesp.prior.mat<-table(star.phenotype.data$Ptype,star.phenotype.data$Ptype.exp)
  }
  if (author=="Chamnanphon"){
    val.data<-assign.act.score.cham(val.data)
    star.phenotype.data<-val.phenotype.alt(val.data)
    #Create bias parameter matrix for SE/SP for validation data
    # 3 x 3 matrix. i=row; j=column
    # i,j cell is Sum(genotype=j | phenotype=i)
    sesp.prior.mat<-table(star.phenotype.data$Ptype,star.phenotype.data$Ptype.exp)
  }
  if (author=="Schroth"){
    val.data<-assign.act.score.schroth(val.data)
    star.phenotype.data<-val.phenotype.schroth(val.data)
    #Create bias parameter matrix for SE/SP for validation data
    # 3 x 3 matrix. i=row; j=column
    # i,j cell is Sum(genotype=j | phenotype=i)
    sesp.prior.mat<-table(star.phenotype.data$Ptype,star.phenotype.data$Ptype.exp)
  }
  if (author=="Mwinyi"){
    val.data<-assign.act.score.mwinyi(val.data)
    star.phenotype.data<-val.phenotype.alt(val.data)
      #Create bias parameter matrix for SE/SP for validation data
      # 3 x 3 matrix. i=row; j=column
      # i,j cell is Sum(genotype=j | phenotype=i)
    sesp.prior.mat<-table(star.phenotype.data$Ptype,star.phenotype.data$Ptype.exp)
  }
  if (author=="Park"){
    val.data<-assign.act.score.park(val.data)
      star.phenotype.data<-val.phenotype.park(val.data)
      #Create bias parameter matrix for SE/SP for validation data
      # 3 x 3 matrix. i=row; j=column
      # i,j cell is Sum(genotype=j | phenotype=i)
      sesp.prior.mat<-table(star.phenotype.data$Ptype,star.phenotype.data$Ptype.exp)
    }
  if (author=="He"){
    val.data<-assign.act.score.he(val.data)
    star.phenotype.data<-val.phenotype.he(val.data)
      #Create bias parameter matrix for SE/SP for validation data
      # 3 x 3 matrix. i=row; j=column
      # i,j cell is Sum(genotype=j | phenotype=i)
    sesp.prior.mat<-table(star.phenotype.data$Ptype,star.phenotype.data$Ptype.exp)
  }
  if (author=="Ismail Al khalil"){
      val.data<-assign.act.score.caudle(val.data)
      star.phenotype.data<-val.phenotype.ism(val.data)
      #Create bias parameter matrix for SE/SP for validation data
      # 3 x 3 matrix. i=row; j=column
      # i,j cell is Sum(genotype=j | phenotype=i)
      sesp.prior.mat<-table(star.phenotype.data$Ptype,star.phenotype.data$Ptype.exp)
  }
  if (author=="De Ameida Melo"){
    val.data<-assign.act.score.melo(val.data)
    star.phenotype.data<-val.phenotype.alt(val.data)
    #Create bias parameter matrix for SE/SP for validation data
    # 3 x 3 matrix. i=row; j=column
    # i,j cell is Sum(genotype=j | phenotype=i)
    sesp.prior.mat<-table(star.phenotype.data$Ptype,star.phenotype.data$Ptype.exp)
  }
  if (author=="Kiyotani"){
    val.data<-assign.act.score.kiyotani(val.data)
    star.phenotype.data<-val.phenotype.alt(val.data)
    #Create bias parameter matrix for SE/SP for validation data
    # 3 x 3 matrix. i=row; j=column
    # i,j cell is Sum(misclass phenotype=j & true phenotype=i)
    sesp.prior.mat<-table(star.phenotype.data$Ptype,star.phenotype.data$Ptype.exp)
    #no 0 genotype so add 0 row
    sesp.prior.mat<-rbind(c(0,0,0),sesp.prior.mat)
  }
  if (author=="Rangel-Mendez"){
    val.data<-assign.act.score.wrang(val.data)
    star.phenotype.data<-val.phenotype.wrang(val.data)
    #Create bias parameter matrix for SE/SP for validation data
    # 3 x 3 matrix. i=row; j=column
    # i,j cell is Sum(genotype=j | phenotype=i)
    sesp.prior.mat<-table(star.phenotype.data$Ptype,star.phenotype.data$Ptype.exp)
  }
      
  #The next few are for studies with impossible values including for completeness but not using
  if (author=="Thompson"){
      val.data<-assign.act.score.schroth(val.data)
      star.phenotype.data<-val.phenotype.schroth(val.data)
      sesp.prior.mat<-table(star.phenotype.data$Ptype,star.phenotype.data$Ptype.exp)
  }
  if (author=="Teh"){
    val.data<-assign.act.score.teh(val.data)
    star.phenotype.data<-val.phenotype.teh(val.data)
    sesp.prior.mat<-table(star.phenotype.data$Ptype,star.phenotype.data$Ptype.exp)
    sesp.prior.mat<-cbind(c(0,0),sesp.prior.mat)
  }
  if (author=="Lan"){
    star.phenotype.data<-val.phenotype.star10(val.data)
    #Create bias parameter matrix for SE/SP for validation data
    # 3 x 3 matrix. i=row; j=column
    # i,j cell is Sum(misclass genotype=j , true phenotype=i)
    sesp.prior.mat<-table(star.phenotype.data$Ptype,star.phenotype.data$Ptype.exp)
    #no 0 genotype so add 0 row
    sesp.prior.mat<-rbind(c(0,0,0),sesp.prior.mat)
  }
  
  
    #Fix sesp.prior for Chamnanphon/Park. no PM and no v/v  
    #including zeros for PM row (none of those in data)
    #including 1 for V/V ---makes code run easier and likely has no impact. 
  if (author=="Chamnanphon"){
      if (!("0" %in% rownames(sesp.prior.mat))){
        temp<-sesp.prior.mat
        temp<-cbind(c(0,0),sesp.prior.mat)
        temp<-rbind(c(0,0,0),temp)
        colnames(temp)<-c("0","1","2")
        rownames(temp)<-c("0","1","2")
        sesp.prior.mat<-temp
        }
    }
  if (author=="Park"){
      if (!("0" %in% rownames(sesp.prior.mat))){
        temp<-sesp.prior.mat
        temp<-rbind(c(0,0,0),temp)
        colnames(temp)<-c("0","1","2")
        rownames(temp)<-c("0","1","2")
        sesp.prior.mat<-temp
      }
    }
  if (author=="Teh"){
    if (!("0" %in% rownames(sesp.prior.mat))){
      temp<-sesp.prior.mat
      temp<-rbind(c(0,0,0),temp)
      colnames(temp)<-c("0","1","2")
      rownames(temp)<-c("0","1","2")
      sesp.prior.mat<-temp
    }
  }

  #Fix 0 cell in data for Ismail by adding 1 to each cell
  if (author=="Ismail Al khalil"){
    vv.1=vv.1+1
    vw.1=vw.1+1
    ww.1=ww.1+1
    vv.0=vv.0+1
    vw.0=vw.0+1
    ww.0=ww.0+1
  }
  
    
  #note that v/v, v/w, w/w is short hand for misclassified phenotype (could be genotype or phenotype depending on study)
  colnames(sesp.prior.mat)<-c("V/V","V/W","W/W")
  rownames(sesp.prior.mat)<-c("PM","IM","NM")
  #Add 1 to validation data to stabilize estimates
  if (region=="As"){
    sesp.prior.mat[2:3,]<-sesp.prior.mat[2:3,]+0
  }
  if (region=="Eur"){
    sesp.prior.mat=sesp.prior.mat+0
  }

  #display Probability of VV/VW/WW given PM,IM,NM
  sesp.prior.mat/rowSums(sesp.prior.mat)
  datalist=list("vv.1"=vv.1,"vw.1"=vw.1,"ww.1"=ww.1,"vv.0"=vv.0,"vw.0"=vw.0,"ww.0"=ww.0,"cohort"=cohort,"samp.type"=samp.type,"sesp.prior.mat"=sesp.prior.mat,"author"=author,"year"=year)
  return(datalist)
}
  
    
###################
#This function adjusts for LOH and incomplete genotyping
cyp.loh.inc<-function(vv.1,vw.1,ww.1,vv.0,vw.0,ww.0,cohort,samp.type,concord,sesp.prior.mat,niter){
  niter.init=niter
  if (cohort==0){
    # crude association for one functional allele versus two
    res_1vs2<-crude.or.ci(vw.1,ww.1,vw.0,ww.0)
    # crude association for no functional allele versus two
    res_0vs2<-crude.or.ci(vv.1,ww.1,vv.0,ww.0)
    # crude association for no or one functional allele versus two
    res_01vs2<-crude.or.ci((vv.1+vw.1),ww.1,(vv.0+vw.0),ww.0)
  }
  if (cohort==1){
    # crude association for one functional allele versus two
    res_1vs2<-crude.rr.ci(vw.1,ww.1,vw.0,ww.0)
    # crude association for no functional allele versus two
    res_0vs2<-crude.rr.ci(vv.1,ww.1,vv.0,ww.0)
    # crude association for no or one functional allele versus two
    res_01vs2<-crude.rr.ci((vv.1+vw.1),ww.1,(vv.0+vw.0),ww.0)
  }
    
  ############################
  #Begin LOH adjustment
  #Load data on concordance
  if (samp.type=="t"){
    adj=1
    if (concord==1){
      vv_1_alpha=11+adj
      vv_2_alpha=0+adj
      vv_3_alpha=1+adj
      vw_1_alpha=3+adj
      vw_2_alpha=72+adj
      vw_3_alpha=15+adj
      ww_1_alpha=2+adj
      ww_2_alpha=2+adj
      ww_3_alpha=195+adj
    }
    if (concord==2){
      vv_1_alpha=3+adj
      vv_2_alpha=0+adj
      vv_3_alpha=0+adj
      vw_1_alpha=0+adj
      vw_2_alpha=38+adj
      vw_3_alpha=0+adj
      ww_1_alpha=0+adj
      ww_2_alpha=0+adj
      ww_3_alpha=74+adj
    }
    if (concord==3){
      vv_1_alpha=8+adj
      vv_2_alpha=0+adj
      vv_3_alpha=1+adj
      vw_1_alpha=3+adj
      vw_2_alpha=34+adj
      vw_3_alpha=15+adj
      ww_1_alpha=2+adj
      ww_2_alpha=2+adj
      ww_3_alpha=121+adj
    }
      
      
    #draw misclassification probabilities given true: Sens/spec analogs
    #These probabilities are prob of tumor gtype [vv,vw,ww] given benign gtype
    # e.g.: p_vv=[pr(vv_tumor|vv_benign),pr(vw_tumor|vv_benign),pr(ww_tumor|vv_benign)]
    # sample from Dirichlet by first sampling ind Gamma(alpha,1)
    p_vv_1<-rgamma(niter,vv_1_alpha,1)
    p_vv_2<-rgamma(niter,vv_2_alpha,1)
    p_vv_3<-rgamma(niter,vv_3_alpha,1)
    all_vv<-cbind(p_vv_1,p_vv_2,p_vv_3)
    p_vv=all_vv/(rowSums(all_vv))
        
    p_vw_1<-rgamma(niter,vw_1_alpha,1)
    p_vw_2<-rgamma(niter,vw_2_alpha,1)
    p_vw_3<-rgamma(niter,vw_3_alpha,1)
    all_vw<-cbind(p_vw_1,p_vw_2,p_vw_3)
    p_vw=all_vw/(rowSums(all_vw))
        
    p_ww_1<-rgamma(niter,ww_1_alpha,1)
    p_ww_2<-rgamma(niter,ww_2_alpha,1)
    p_ww_3<-rgamma(niter,ww_3_alpha,1)
    all_ww<-cbind(p_ww_1,p_ww_2,p_ww_3)
    p_ww=all_ww/(rowSums(all_ww))
        
    #Create 3x3xniter array of misclassification probabilities (se/sp)
    #for [,,i] matrix in array, cells are sens/spec: eg. pr(vv*|vv), pr(vw*|vv), etc 
    #     columns are 1=vv; 2=vw; 3=ww
    #     rows are 1=vv*; 2=vw*; 3=ww*
    P.array=array(t(cbind(p_vv,p_vw,p_ww)),c(3,3,niter))
    #Generate adjusted cell counts for each 3x3 misclass matrix in the array
    # returns matrix: 3x niter. each column and row 1=vv; row2=vw; row3=ww
    Adj.1_LOH_1=apply(P.array,3,mis.adj.counts,counts=c(vv.1,vw.1,ww.1))
    Adj.0_LOH_1=apply(P.array,3,mis.adj.counts,counts=c(vv.0,vw.0,ww.0))
    #find impossible values & eliminate those iterations
    Trap1<-((Adj.1_LOH_1[1,]<0)|(Adj.1_LOH_1[2,]<0)|(Adj.1_LOH_1[3,]<0)|(Adj.0_LOH_1[1,]<0)|(Adj.0_LOH_1[2,]<0)|(Adj.0_LOH_1[3,]<0))
    Adj.1_LOH_1<-Adj.1_LOH_1[,Trap1=="FALSE"]
    Adj.0_LOH_1<-Adj.0_LOH_1[,Trap1=="FALSE"]
    P.array<-P.array[,,Trap1=="FALSE"]
    niter=dim(Adj.0_LOH_1)[2]
        
    # calculate prevalence of true (benign) exposure (vv,vw,ww) in cases and controls, accounting for sampling error
    # Draw from Dirichlet Dist'n by first sampling from gammas
    vv.1.gamma<-rgamma(niter,Adj.1_LOH_1[1,],1)
    vw.1.gamma<-rgamma(niter,Adj.1_LOH_1[2,],1)
    ww.1.gamma<-rgamma(niter,Adj.1_LOH_1[3,],1)
    all.1.gamma<-cbind(vv.1.gamma,vw.1.gamma,ww.1.gamma)
    #standardized gammas are a Dirichlet draw
    PrevE.1<-all.1.gamma/(rowSums(all.1.gamma))
        
    vv.0.gamma<-rgamma(niter,Adj.0_LOH_1[1,],1)
    vw.0.gamma<-rgamma(niter,Adj.0_LOH_1[2,],1)
    ww.0.gamma<-rgamma(niter,Adj.0_LOH_1[3,],1)
    all.0.gamma<-cbind(vv.0.gamma,vw.0.gamma,ww.0.gamma)
    PrevE.0<-all.0.gamma/(rowSums(all.0.gamma))
        
    # calculate Predictive Values of exposure classification in cases and controls
    # naming: PV_XX_xx_z: PV=predictive value;of XX= true genotype;given xx=tumor genotype;and given z=case/control
    PPV_1_loh<-array(unlist(lapply(1:niter,function(i)pv.fun(P.array[,,i],PrevE.1[i,]))),dim=c(3,3,niter))
    PPV_0_loh<-array(unlist(lapply(1:niter,function(i)pv.fun(P.array[,,i],PrevE.0[i,]))),dim=c(3,3,niter))
        
    #calculate the expected number of exposed cases and controls conditional on 
    #   Observed allele status
    count_vv_1=apply(t(PPV_1_loh[1,,]),1,rmultinom,n=1,size=vv.1)
    count_vw_1=apply(t(PPV_1_loh[2,,]),1,rmultinom,n=1,size=vw.1)
    count_ww_1=apply(t(PPV_1_loh[3,,]),1,rmultinom,n=1,size=ww.1)
        
    count_vv_0=apply(t(PPV_0_loh[1,,]),1,rmultinom,n=1,size=vv.0)
    count_vw_0=apply(t(PPV_0_loh[2,,]),1,rmultinom,n=1,size=vw.0)
    count_ww_0=apply(t(PPV_0_loh[3,,]),1,rmultinom,n=1,size=ww.0)
        
    #aggregate across observe allele, to get total expected number of true allele
    vv.1_loh=count_vv_1[1,]+count_vw_1[1,]+count_ww_1[1,]
    vw.1_loh=count_vv_1[2,]+count_vw_1[2,]+count_ww_1[2,]
    ww.1_loh=count_vv_1[3,]+count_vw_1[3,]+count_ww_1[3,]
    vv.0_loh=count_vv_0[1,]+count_vw_0[1,]+count_ww_0[1,]
    vw.0_loh=count_vv_0[2,]+count_vw_0[2,]+count_ww_0[2,]
    ww.0_loh=count_vv_0[3,]+count_vw_0[3,]+count_ww_0[3,]
        
    #dropping zero cells.
    Trap0=vv.1_loh==0|vw.1_loh==0|ww.1_loh==0|vv.0_loh==0|vw.0_loh==0|ww.0_loh==0
    vv.1_loh=vv.1_loh[Trap0=="FALSE"]
    vw.1_loh=vw.1_loh[Trap0=="FALSE"]
    ww.1_loh=ww.1_loh[Trap0=="FALSE"]
    vv.0_loh=vv.0_loh[Trap0=="FALSE"]
    vw.0_loh=vw.0_loh[Trap0=="FALSE"]
    ww.0_loh=ww.0_loh[Trap0=="FALSE"]
    niter=table(Trap0)[1]
        
    if (cohort==0){
      # calculate bias-adjusted OR (1 functional allele vs 2) with LOH uncertainty
      or_1vs2_loh <- (vw.1_loh/ww.1_loh)/(vw.0_loh/ww.0_loh)
      se_or_1vs2_loh <- sqrt(1/vw.1_loh +1/ww.1_loh+ 1/vw.0_loh+ 1/ww.0_loh)
      or_1vs2_loh_r <- exp(rnorm(length(or_1vs2_loh),log(or_1vs2_loh),se_or_1vs2_loh))
          
      # calculate bias-adjusted OR (0 functional allele vs 2) with LOH uncertainty
      or_0vs2_loh <- (vv.1_loh/ww.1_loh)/(vv.0_loh/ww.0_loh)
      se_or_0vs2_loh <- sqrt(1/vv.1_loh +1/ww.1_loh+ 1/vv.0_loh+ 1/ww.0_loh)
      or_0vs2_loh_r <- exp(rnorm(length(or_0vs2_loh),log(or_0vs2_loh),se_or_0vs2_loh))
          
      # calculate bias-adjusted OR (0/1 functional allele vs 2) with LOH uncertainty
      or_01vs2_loh <- ((vw.1_loh+vv.1_loh)/ww.1_loh)/((vw.0_loh+vv.0_loh)/ww.0_loh)
      se_or_01vs2_loh <- sqrt(1/(vw.1_loh+vv.1_loh) +1/ww.1_loh+ 1/(vw.0_loh+vv.0_loh)+ 1/ww.0_loh)
      or_01vs2_loh_r <- exp(rnorm(length(or_01vs2_loh),log(or_01vs2_loh),se_or_01vs2_loh))
    }
        
    if (cohort==1){
      # calculate bias-adjusted OR (1 functional allele vs 2) with LOH uncertainty
      or_1vs2_loh <- (vw.1_loh/(vw.1_loh+vw.0_loh))/(ww.1_loh/(ww.1_loh+ww.0_loh))
      se_or_1vs2_loh <- sqrt(1/vw.1_loh +1/ww.1_loh - 1/(vw.1_loh+vw.0_loh) - 1/(ww.1_loh+ww.0_loh))
      or_1vs2_loh_r <- exp(rnorm(length(or_1vs2_loh),log(or_1vs2_loh),se_or_1vs2_loh))
          
      # calculate bias-adjusted OR (0 functional allele vs 2) with LOH uncertainty
      or_0vs2_loh <- (vv.1_loh/(vv.1_loh+vv.0_loh))/(ww.1_loh/(ww.1_loh+ww.0_loh))
      se_or_0vs2_loh <- sqrt(1/vv.1_loh +1/ww.1_loh - 1/(vv.1_loh+vv.0_loh)+ 1/(ww.1_loh+ww.0_loh))
      or_0vs2_loh_r <- exp(rnorm(length(or_0vs2_loh),log(or_0vs2_loh),se_or_0vs2_loh))
        
      # calculate bias-adjusted OR (0/1 functional allele vs 2) with LOH uncertainty
      or_01vs2_loh <- ((vw.1_loh+vv.1_loh)/(vw.1_loh+vv.1_loh+vw.0_loh+vv.0_loh))/(ww.1_loh/(ww.1_loh+ww.0_loh))
      se_or_01vs2_loh <- sqrt(1/(vw.1_loh+vv.1_loh) + 1/ww.1_loh - 1/(vw.1_loh+vv.1_loh+vw.0_loh+vv.0_loh)- 1/(ww.1_loh+ww.0_loh))
      or_01vs2_loh_r <- exp(log(or_01vs2_loh) - rnorm(length(or_01vs2_loh),0,se_or_01vs2_loh))
    }
        
        
  }
    #End LOH adjustment
      
      
      
      
      
  ##############################################################################
  #Begin misclassification adjustment for incomplete genotyping
    
  #draw misclassification probabilities given true: Sens/spec analogs
  #These probabilities are prob of *4 genotype (or ptype) [vv,vw,ww] given true 
  #   phenotype (PM,IM,NM)
  # e.g.: p_pm=[pr(vv|PM),pr(vw|PM),pr(ww|PM)]
  # sample from Dirichlet by first sampling ind Gamma(alpha,1)
  p_pm_1<-rgamma(niter,sesp.prior.mat[1,1],1)
  p_pm_2<-rgamma(niter,sesp.prior.mat[1,2],1)
  p_pm_3<-rgamma(niter,sesp.prior.mat[1,3],1)
  all_pm<-cbind(p_pm_1,p_pm_2,p_pm_3)
  p_pm=all_pm/(rowSums(all_pm))
      
  p_im_1<-rgamma(niter,sesp.prior.mat[2,1],1)
  p_im_2<-rgamma(niter,sesp.prior.mat[2,2],1)
  p_im_3<-rgamma(niter,sesp.prior.mat[2,3],1)
  all_im<-cbind(p_im_1,p_im_2,p_im_3)
  p_im=all_im/(rowSums(all_im))
      
  p_nm_1<-rgamma(niter,sesp.prior.mat[3,1],1)
  p_nm_2<-rgamma(niter,sesp.prior.mat[3,2],1)
  p_nm_3<-rgamma(niter,sesp.prior.mat[3,3],1)
  all_nm<-cbind(p_nm_1,p_nm_2,p_nm_3)
  p_nm=all_nm/(rowSums(all_nm))
      
  #Create 3x3xniter array of misclassification probabilities (se/sp)
  #for [,,i] matrix in array, cells are sens/spec: eg. pr(vv|pm), pr(vw|pm), etc 
  #     columns are 1=pm; 2=im; 3=nm
  #     rows are 1=vv; 2=vw; 3=ww
  P.inc.array=array(t(cbind(p_pm,p_im,p_nm)),c(3,3,niter))
  #Generate adjusted cell counts for each 3x3 misclass matrix in the array
  # returns matrix: 3x niter. each column and row 1=vv; row2=vw; row3=ww
  if (samp.type=="t"){
    gtype.1_loh=cbind(vv.1_loh,vw.1_loh,ww.1_loh)
    gtype.0_loh=cbind(vv.0_loh,vw.0_loh,ww.0_loh)
    Adj.1_inc_1<-matrix(unlist(lapply(1:niter,function(i)mis.adj.counts(P.inc.array[,,i],gtype.1_loh[i,]))),ncol=niter)
    Adj.0_inc_1<-matrix(unlist(lapply(1:niter,function(i)mis.adj.counts(P.inc.array[,,i],gtype.0_loh[i,]))),ncol=niter)
  }
  if (samp.type=="n"){
    Adj.1_inc_1=apply(P.inc.array,3,mis.adj.counts,counts=c(vv.1,vw.1,ww.1))
    Adj.0_inc_1=apply(P.inc.array,3,mis.adj.counts,counts=c(vv.0,vw.0,ww.0))
  }
      
  #find impossible values & eliminate those iterations
  Trap1.inc<-((Adj.1_inc_1[1,]<=0)|(Adj.1_inc_1[2,]<=0)|(Adj.1_inc_1[3,]<=0)|(Adj.0_inc_1[1,]<=0)|(Adj.0_inc_1[2,]<=0)|(Adj.0_inc_1[3,]<=0))
  Adj.1_inc_1<-Adj.1_inc_1[,Trap1.inc=="FALSE"]
  Adj.0_inc_1<-Adj.0_inc_1[,Trap1.inc=="FALSE"]
  P.inc.array<-P.inc.array[,,Trap1.inc=="FALSE"]
  niter=dim(Adj.1_inc_1)[2]
      
  # calculate prevalence of exposure (pm,im,nm) in cases and controls, accounting for sampling error
  # Draw from Dirichlet Dist'n by first sampling from gammas
  vv.1.gamma<-rgamma(niter,Adj.1_inc_1[1,],1)
  vw.1.gamma<-rgamma(niter,Adj.1_inc_1[2,],1)
  ww.1.gamma<-rgamma(niter,Adj.1_inc_1[3,],1)
  all.1.gamma<-cbind(vv.1.gamma,vw.1.gamma,ww.1.gamma)
  #standardized gammas are a Dirichlet draw
  PrevE.1<-all.1.gamma/(rowSums(all.1.gamma))
      
  vv.0.gamma<-rgamma(niter,Adj.0_inc_1[1,],1)
  vw.0.gamma<-rgamma(niter,Adj.0_inc_1[2,],1)
  ww.0.gamma<-rgamma(niter,Adj.0_inc_1[3,],1)
  all.0.gamma<-cbind(vv.0.gamma,vw.0.gamma,ww.0.gamma)
  PrevE.0<-all.0.gamma/(rowSums(all.0.gamma))
    
  P.array=P.inc.array
  # calculate Predictive Values of exposure classification in cases and controls
  # naming: PV_XX_xx_z: PV=predictive value;of XX= true Phenotype; 
  #                     given xx=genotype/obs phenotype;and given z=case/control
  # arrays are indexed by i,j,k.  K is the number of iterations
  #                               i (row) is observed genotype/phenothype xx
  #                               j (col) is true phenotype XX

  PPV_1_inc<-array(unlist(lapply(1:niter,function(i)pv.fun(P.array[,,i],PrevE.1[i,]))),dim=c(3,3,niter))
  PPV_0_inc<-array(unlist(lapply(1:niter,function(i)pv.fun(P.array[,,i],PrevE.0[i,]))),dim=c(3,3,niter))
      
  #calculate the expected number of true exposed cases and controls conditional on 
  #   Observed allele status
  if (samp.type=="t"){
    count_vv_1=t(rmultinomial(n=niter,size=vv.1_loh,prob=t(PPV_1_inc[1,,])))
    count_vw_1=t(rmultinomial(n=niter,size=vw.1_loh,prob=t(PPV_1_inc[2,,])))
    count_ww_1=t(rmultinomial(n=niter,size=ww.1_loh,prob=t(PPV_1_inc[3,,])))
      
    count_vv_0=t(rmultinomial(n=niter,size=vv.0_loh,t(PPV_0_inc[1,,])))
    count_vw_0=t(rmultinomial(n=niter,size=vw.0_loh,t(PPV_0_inc[2,,])))
    count_ww_0=t(rmultinomial(n=niter,size=ww.0_loh,t(PPV_0_inc[3,,])))
  }
  if (samp.type=="n"){
    count_vv_1=apply(t(PPV_1_inc[1,,]),1,rmultinom,n=1,size=vv.1)
    count_vw_1=apply(t(PPV_1_inc[2,,]),1,rmultinom,n=1,size=vw.1)
    count_ww_1=apply(t(PPV_1_inc[3,,]),1,rmultinom,n=1,size=ww.1)
        
    count_vv_0=apply(t(PPV_0_inc[1,,]),1,rmultinom,n=1,size=vv.0)
    count_vw_0=apply(t(PPV_0_inc[2,,]),1,rmultinom,n=1,size=vw.0)
    count_ww_0=apply(t(PPV_0_inc[3,,]),1,rmultinom,n=1,size=ww.0)
  }
    
  #aggregate across observe allele, to get total expected number of true allele
  pm.1_inc=count_vv_1[1,]+count_vw_1[1,]+count_ww_1[1,]
  im.1_inc=count_vv_1[2,]+count_vw_1[2,]+count_ww_1[2,]
  nm.1_inc=count_vv_1[3,]+count_vw_1[3,]+count_ww_1[3,]
  pm.0_inc=count_vv_0[1,]+count_vw_0[1,]+count_ww_0[1,]
  im.0_inc=count_vv_0[2,]+count_vw_0[2,]+count_ww_0[2,]
  nm.0_inc=count_vv_0[3,]+count_vw_0[3,]+count_ww_0[3,]
      
  #remove zeros
  trap2=(pm.1_inc<=0)|(im.1_inc<=0)|(nm.1_inc<=0)|(pm.0_inc<=0)|(im.0_inc<=0)|(nm.0_inc<=0)
  pm.1_inc=pm.1_inc[trap2=="FALSE"]
  im.1_inc=im.1_inc[trap2=="FALSE"]
  nm.1_inc=nm.1_inc[trap2=="FALSE"]
  pm.0_inc=pm.0_inc[trap2=="FALSE"]
  im.0_inc=im.0_inc[trap2=="FALSE"]
  nm.0_inc=nm.0_inc[trap2=="FALSE"]
  niter=length(pm.1_inc)
      
  if (cohort==0){
    # calculate bias-adjusted OR (IM functional allele vs nm) with inc uncertainty
    res_im_vs_nm <- t(matrix(unlist(lapply(1:niter,function(i)crude.or.ci(im.1_inc[i],nm.1_inc[i],im.0_inc[i],nm.0_inc[i]))),nrow=4))
    or_im_vs_nm_r <- exp(rnorm(niter,log(res_im_vs_nm[,1]),res_im_vs_nm[,2]))
        
    # calculate bias-adjusted OR (pm functional allele vs nm) with inc uncertainty
    res_pm_vs_nm <- t(matrix(unlist(lapply(1:niter,function(i)crude.or.ci(pm.1_inc[i],nm.1_inc[i],pm.0_inc[i],nm.0_inc[i]))),nrow=4))
    or_pm_vs_nm_r <- exp(rnorm(niter,log(res_pm_vs_nm[,1]),res_pm_vs_nm[,2]))
        
    # calculate bias-adjusted OR (pm/im functional allele vs nm) with inc uncertainty
    res_pmim_vs_nm <- t(matrix(unlist(lapply(1:niter,function(i)crude.or.ci((im.1_inc[i]+pm.1_inc[i]),nm.1_inc[i],(im.0_inc[i]+pm.0_inc[i]),nm.0_inc[i]))),nrow=4))
    or_pmim_vs_nm_r <- exp(rnorm(niter,log(res_pmim_vs_nm[,1]),res_pmim_vs_nm[,2]))
  }
      
  if (cohort==1){
    # calculate bias-adjusted RR (IM functional allele vs nm) with inc uncertainty
    res_im_vs_nm <- t(matrix(unlist(lapply(1:niter,function(i)crude.rr.ci(im.1_inc[i],nm.1_inc[i],im.0_inc[i],nm.0_inc[i]))),nrow=4))
    or_im_vs_nm_r <- exp(rnorm(niter,log(res_im_vs_nm[,1]),res_im_vs_nm[,2]))
        
    # calculate bias-adjusted RR (pm functional allele vs nm) with inc uncertainty
    res_pm_vs_nm <- t(matrix(unlist(lapply(1:niter,function(i)crude.rr.ci(pm.1_inc[i],nm.1_inc[i],pm.0_inc[i],nm.0_inc[i]))),nrow=4))
    or_pm_vs_nm_r <- exp(rnorm(niter,log(res_pm_vs_nm[,1]),res_pm_vs_nm[,2]))
        
    # calculate bias-adjusted RR (pm/im functional allele vs nm) with inc uncertainty
    res_pmim_vs_nm <- t(matrix(unlist(lapply(1:niter,function(i)crude.rr.ci((im.1_inc[i]+pm.1_inc[i]),nm.1_inc[i],(im.0_inc[i]+pm.0_inc[i]),nm.0_inc[i]))),nrow=4))
    or_pmim_vs_nm_r <- exp(rnorm(niter,log(res_pmim_vs_nm[,1]),res_pmim_vs_nm[,2]))
  }
      
  if (samp.type=="t"){
    LOH<-quantile(or_1vs2_loh_r,c(.5,.025,.975))
    LOH_Inc<-quantile(or_im_vs_nm_r,c(.5,.025,.975))
    or1vs2<-rbind(res_1vs2[c(1,3,4)],LOH,LOH_Inc)
    rownames(or1vs2)<-c("Crude","LOH only","LOH & INC")
      
    LOH<-quantile(or_0vs2_loh_r,c(.5,.025,.975))
    LOH_Inc<-quantile(or_pm_vs_nm_r,c(.5,.025,.975))
    or0vs2<-rbind(res_0vs2[c(1,3,4)],LOH,LOH_Inc)
    rownames(or0vs2)<-c("Crude","LOH only","LOH & INC")
        
    LOH<-quantile(or_01vs2_loh_r,c(.5,.025,.975))
    LOH_Inc<-quantile(or_pmim_vs_nm_r,c(.5,.025,.975))
    or01vs2<-rbind(res_01vs2[c(1,3,4)],LOH,LOH_Inc)
    rownames(or01vs2)<-c("Crude","LOH only","LOH & INC")
  }
  if (samp.type=="n"){
    LOH<-c("NA","NA","NA")
    LOH_Inc<-quantile(or_im_vs_nm_r,c(.5,.025,.975))
    or1vs2<-rbind(res_1vs2[c(1,3,4)],LOH,LOH_Inc)
    rownames(or1vs2)<-c("Crude","LOH only","INC")
        
    LOH<-c("NA","NA","NA")
    LOH_Inc<-quantile(or_pm_vs_nm_r,c(.5,.025,.975))
    or0vs2<-rbind(res_0vs2[c(1,3,4)],LOH,LOH_Inc)
    rownames(or0vs2)<-c("Crude","LOH only","INC")
        
    LOH<-c("NA","NA","NA")
    LOH_Inc<-quantile(or_pmim_vs_nm_r,c(.5,.025,.975))
    or01vs2<-rbind(res_01vs2[c(1,3,4)],LOH,LOH_Inc)
    rownames(or01vs2)<-c("Crude","LOH only","INC")
  }
      
      
    dropped<-c(niter.init,niter)
    dropped[is.na(dropped)]=0
    names(dropped)=c("Initial","Final")
    return.list=list("Het"=or1vs2,"Hom"=or0vs2,"Any"=or01vs2,"Drop"=dropped)
  }
    

###################
#This function adjusts for LOH and incomplete genotyping in ASIAN Populations
cyp.loh.inc.asian<-function(vv.1,vw.1,ww.1,vv.0,vw.0,ww.0,cohort,samp.type,concord,sesp.prior.mat,niter){
    niter.init=niter
    if (cohort==0){
      # crude association for one functional allele versus two
      res_1vs2<-crude.or.ci(vw.1,ww.1,vw.0,ww.0)
      # crude association for no functional allele versus two
      res_0vs2<-crude.or.ci(vv.1,ww.1,vv.0,ww.0)
      # crude association for no or one functional allele versus two
      res_01vs2<-crude.or.ci((vv.1+vw.1),ww.1,(vv.0+vw.0),ww.0)
    }
    if (cohort==1){
      # crude association for one functional allele versus two
      res_1vs2<-crude.rr.ci(vw.1,ww.1,vw.0,ww.0)
      # crude association for no functional allele versus two
      res_0vs2<-crude.rr.ci(vv.1,ww.1,vv.0,ww.0)
      # crude association for no or one functional allele versus two
      res_01vs2<-crude.rr.ci((vv.1+vw.1),ww.1,(vv.0+vw.0),ww.0)
    }
    
    ############################
    #Begin LOH adjustment
    #Load data on concordance
    if (samp.type=="t"){
      if (concord==1){
        vv_1_alpha=11.5
        vv_2_alpha=0.5
        vv_3_alpha=1.5
        vw_1_alpha=3.5
        vw_2_alpha=72.5
        vw_3_alpha=15.5
        ww_1_alpha=2.5
        ww_2_alpha=2.5
        ww_3_alpha=195.5
      }
      if (concord==2){
        vv_1_alpha=3.5
        vv_2_alpha=0.5
        vv_3_alpha=0.5
        vw_1_alpha=0.5
        vw_2_alpha=38.5
        vw_3_alpha=0.5
        ww_1_alpha=0.5
        ww_2_alpha=0.5
        ww_3_alpha=74.5
      }
      if (concord==3){
        vv_1_alpha=8.5
        vv_2_alpha=0.5
        vv_3_alpha=1.5
        vw_1_alpha=3.5
        vw_2_alpha=34.5
        vw_3_alpha=15.5
        ww_1_alpha=2.5
        ww_2_alpha=2.5
        ww_3_alpha=121.5
      }
      #draw misclassification probabilities given true: Sens/spec analogs
      #These probabilities are prob of tumor [vv,vw,ww] given benign type
      # e.g.: p_vv=[pr(vv_tumor|vv_benign),pr(vw_tumor|vv_benign),pr(ww_tumor|vv_benign)]
      # sample from Dirichlet by first sampling ind Gamma(alpha,1)
      p_vv_1<-rgamma(niter,vv_1_alpha,1)
      p_vv_2<-rgamma(niter,vv_2_alpha,1)
      p_vv_3<-rgamma(niter,vv_3_alpha,1)
      all_vv<-cbind(p_vv_1,p_vv_2,p_vv_3)
      p_vv=all_vv/(rowSums(all_vv))
      
      p_vw_1<-rgamma(niter,vw_1_alpha,1)
      p_vw_2<-rgamma(niter,vw_2_alpha,1)
      p_vw_3<-rgamma(niter,vw_3_alpha,1)
      all_vw<-cbind(p_vw_1,p_vw_2,p_vw_3)
      p_vw=all_vw/(rowSums(all_vw))
      
      p_ww_1<-rgamma(niter,ww_1_alpha,1)
      p_ww_2<-rgamma(niter,ww_2_alpha,1)
      p_ww_3<-rgamma(niter,ww_3_alpha,1)
      all_ww<-cbind(p_ww_1,p_ww_2,p_ww_3)
      p_ww=all_ww/(rowSums(all_ww))
      
      #Create 3x3xniter array of misclassification probabilities (se/sp)
      #for [,,i] matrix in array, cells are sens/spec: eg. pr(vv*|vv), pr(vw*|vv), etc 
      #     columns are 1=vv; 2=vw; 3=ww
      #     rows are 1=vv*; 2=vw*; 3=ww*
      P.array=array(t(cbind(p_vv,p_vw,p_ww)),c(3,3,niter))
      #Generate adjusted cell counts for each 3x3 misclass matrix in the array
      # returns matrix: 3x niter. each column and row 1=vv; row2=vw; row3=ww
      Adj.1_LOH_1=apply(P.array,3,mis.adj.counts,counts=c(vv.1,vw.1,ww.1))
      Adj.0_LOH_1=apply(P.array,3,mis.adj.counts,counts=c(vv.0,vw.0,ww.0))
      #find impossible values & eliminate those iterations
      Trap1<-((Adj.1_LOH_1[1,]<0)|(Adj.1_LOH_1[2,]<0)|(Adj.1_LOH_1[3,]<0)|(Adj.0_LOH_1[1,]<0)|(Adj.0_LOH_1[2,]<0)|(Adj.0_LOH_1[3,]<0))
      Adj.1_LOH_1<-Adj.1_LOH_1[,Trap1=="FALSE"]
      Adj.0_LOH_1<-Adj.0_LOH_1[,Trap1=="FALSE"]
      P.array<-P.array[,,Trap1=="FALSE"]
      niter=dim(Adj.0_LOH_1)[2]
      
      # calculate prevalence of exposure (vv,vw,ww) in cases and controls, accounting for sampling error
      # Draw from Dirichlet Dist'n by first sampling from gammas
      vv.1.gamma<-rgamma(niter,Adj.1_LOH_1[1,],1)
      vw.1.gamma<-rgamma(niter,Adj.1_LOH_1[2,],1)
      ww.1.gamma<-rgamma(niter,Adj.1_LOH_1[3,],1)
      all.1.gamma<-cbind(vv.1.gamma,vw.1.gamma,ww.1.gamma)
      #standardized gammas are a Dirichlet draw
      PrevE.1<-all.1.gamma/(rowSums(all.1.gamma))
      
      vv.0.gamma<-rgamma(niter,Adj.0_LOH_1[1,],1)
      vw.0.gamma<-rgamma(niter,Adj.0_LOH_1[2,],1)
      ww.0.gamma<-rgamma(niter,Adj.0_LOH_1[3,],1)
      all.0.gamma<-cbind(vv.0.gamma,vw.0.gamma,ww.0.gamma)
      PrevE.0<-all.0.gamma/(rowSums(all.0.gamma))
      
      # calculate Predictive Values of exposure classification in cases and controls
      # namin: PV_XX_xx_z: PV=predictive value;of XX= true genotype;given xx=tumor genotype;and given z=case/control
      
      PPV_1_loh<-array(unlist(lapply(1:niter,function(i)pv.fun(P.array[,,i],PrevE.1[i,]))),dim=c(3,3,niter))
      PPV_0_loh<-array(unlist(lapply(1:niter,function(i)pv.fun(P.array[,,i],PrevE.0[i,]))),dim=c(3,3,niter))
      
      
      #calculate the expected number of exposed cases and controls conditional on 
      #   Observed allele status
      count_vv_1=apply(t(PPV_1_loh[1,,]),1,rmultinom,n=1,size=vv.1)
      count_vw_1=apply(t(PPV_1_loh[2,,]),1,rmultinom,n=1,size=vw.1)
      count_ww_1=apply(t(PPV_1_loh[3,,]),1,rmultinom,n=1,size=ww.1)
      
      count_vv_0=apply(t(PPV_0_loh[1,,]),1,rmultinom,n=1,size=vv.0)
      count_vw_0=apply(t(PPV_0_loh[2,,]),1,rmultinom,n=1,size=vw.0)
      count_ww_0=apply(t(PPV_0_loh[3,,]),1,rmultinom,n=1,size=ww.0)
      
      #aggregate across observe allele, to get total expected number of true allele
      vv.1_loh=count_vv_1[1,]+count_vw_1[1,]+count_ww_1[1,]
      vw.1_loh=count_vv_1[2,]+count_vw_1[2,]+count_ww_1[2,]
      ww.1_loh=count_vv_1[3,]+count_vw_1[3,]+count_ww_1[3,]
      vv.0_loh=count_vv_0[1,]+count_vw_0[1,]+count_ww_0[1,]
      vw.0_loh=count_vv_0[2,]+count_vw_0[2,]+count_ww_0[2,]
      ww.0_loh=count_vv_0[3,]+count_vw_0[3,]+count_ww_0[3,]
      
      #dropping zero cells. maybe add 1?
      Trap0=vv.1_loh==0|vw.1_loh==0|ww.1_loh==0|vv.0_loh==0|vw.0_loh==0|ww.0_loh==0
      vv.1_loh=vv.1_loh[Trap0=="FALSE"]
      vw.1_loh=vw.1_loh[Trap0=="FALSE"]
      ww.1_loh=ww.1_loh[Trap0=="FALSE"]
      vv.0_loh=vv.0_loh[Trap0=="FALSE"]
      vw.0_loh=vw.0_loh[Trap0=="FALSE"]
      ww.0_loh=ww.0_loh[Trap0=="FALSE"]
      niter=table(Trap0)[1]
      
      if (cohort==0){
        # calculate bias-adjusted OR (1 functional allele vs 2) with LOH uncertainty
        or_1vs2_loh <- (vw.1_loh/ww.1_loh)/(vw.0_loh/ww.0_loh)
        se_or_1vs2_loh <- sqrt(1/vw.1_loh +1/ww.1_loh+ 1/vw.0_loh+ 1/ww.0_loh)
        or_1vs2_loh_r <- exp(rnorm(length(or_1vs2_loh),log(or_1vs2_loh),se_or_1vs2_loh))
        
        # calculate bias-adjusted OR (0 functional allele vs 2) with LOH uncertainty
        or_0vs2_loh <- (vv.1_loh/ww.1_loh)/(vv.0_loh/ww.0_loh)
        se_or_0vs2_loh <- sqrt(1/vv.1_loh +1/ww.1_loh+ 1/vv.0_loh+ 1/ww.0_loh)
        or_0vs2_loh_r <- exp(rnorm(length(or_0vs2_loh),log(or_0vs2_loh),se_or_0vs2_loh))
        
        # calculate bias-adjusted OR (0/1 functional allele vs 2) with LOH uncertainty
        or_01vs2_loh <- ((vw.1_loh+vv.1_loh)/ww.1_loh)/((vw.0_loh+vv.0_loh)/ww.0_loh)
        se_or_01vs2_loh <- sqrt(1/(vw.1_loh+vv.1_loh) +1/ww.1_loh+ 1/(vw.0_loh+vv.0_loh)+ 1/ww.0_loh)
        or_01vs2_loh_r <- exp(rnorm(length(or_01vs2_loh),log(or_01vs2_loh),se_or_01vs2_loh))
      }
        
      if (cohort==1){
        # calculate bias-adjusted OR (1 functional allele vs 2) with LOH uncertainty
        or_1vs2_loh <- (vw.1_loh/(vw.1_loh+vw.0_loh))/(ww.1_loh/(ww.1_loh+ww.0_loh))
        se_or_1vs2_loh <- sqrt(1/vw.1_loh +1/ww.1_loh - 1/(ww.1_loh+vw.0_loh) - 1/(ww.1_loh+ww.0_loh))
        or_1vs2_loh_r <- exp(rnorm(length(or_1vs2_loh),log(or_1vs2_loh),se_or_1vs2_loh))
        
        # calculate bias-adjusted OR (0 functional allele vs 2) with LOH uncertainty
        or_0vs2_loh <- (vv.1_loh/(vv.1_loh+vv.0_loh))/(ww.1_loh/(ww.1_loh+ww.0_loh))
        se_or_0vs2_loh <- sqrt(1/vv.1_loh +1/ww.1_loh - 1/(vv.1_loh+vv.0_loh)+ 1/(ww.1_loh+ww.0_loh))
        or_0vs2_loh_r <- exp(rnorm(length(or_0vs2_loh),log(or_0vs2_loh),se_or_0vs2_loh))
        
        # calculate bias-adjusted OR (0/1 functional allele vs 2) with LOH uncertainty
        or_01vs2_loh <- ((vw.1_loh+vv.1_loh)/(vw.1_loh+vv.1_loh+vw.0_loh+vv.0_loh))/(ww.1_loh/(ww.1_loh+ww.0_loh))
        se_or_01vs2_loh <- sqrt(1/(vw.1_loh+vv.1_loh) + 1/ww.1_loh - 1/(vw.1_loh+vv.1_loh+vw.0_loh+vv.0_loh)- 1/(ww.1_loh+ww.0_loh))
        or_01vs2_loh_r <- exp(log(or_01vs2_loh) - rnorm(length(or_01vs2_loh),0,se_or_01vs2_loh))
      }
        
        
    }
    #End LOH adjustment
      
      
      
      
      
    ##############################################################################
    #Begin misclassification adjustment for incomplete genotyping
      
    #draw misclassification probabilities given true: Sens/spec analogs
    #These probabilities are prob of *4 genotype [vv,vw,ww] given true phenotype (PM,IM,NM)
    # e.g.: p_im=[pr(vv|IM),pr(vw|IM),pr(ww|IM)]
    # sample from Dirichlet by first sampling ind Gamma(alpha,1)
    # of note: there are no PM in the validation data so we are not sampling them
    
    p_im_1<-rgamma(niter,sesp.prior.mat[2,1],1)
    p_im_2<-rgamma(niter,sesp.prior.mat[2,2],1)
    p_im_3<-rgamma(niter,sesp.prior.mat[2,3],1)
    all_im<-cbind(p_im_1,p_im_2,p_im_3)
    p_im=all_im/(rowSums(all_im))
      
    p_nm_1<-rgamma(niter,sesp.prior.mat[3,1],1)
    p_nm_2<-rgamma(niter,sesp.prior.mat[3,2],1)
    p_nm_3<-rgamma(niter,sesp.prior.mat[3,3],1)
    all_nm<-cbind(p_nm_1,p_nm_2,p_nm_3)
    p_nm=all_nm/(rowSums(all_nm))
      
    #Create 3x2xniter array of misclassification probabilities (se/sp)
    #for [,,i] matrix in array, cells are sens/spec: eg. pr(vv|im), pr(vw|im), etc 
    #     columns are  1=im; 2=nm
    #     rows are 1=vv; 2=vw; 3=ww
    P.inc.array=array(t(cbind(p_im,p_nm)),c(3,2,niter))
    #Generate adjusted cell counts for each 3x3 misclass matrix in the array
    # returns matrix: 3x niter. each column and row 1=vv; row2=vw; row3=ww
    if (samp.type=="t"){
      gtype.1_loh=cbind(vv.1_loh,vw.1_loh,ww.1_loh)
      gtype.0_loh=cbind(vv.0_loh,vw.0_loh,ww.0_loh)
      Adj.1_inc_1<-matrix(unlist(lapply(1:niter,function(i)mis.adj.counts.ls(P.inc.array[,,i],gtype.1_loh[i,]))),ncol=niter)
      Adj.0_inc_1<-matrix(unlist(lapply(1:niter,function(i)mis.adj.counts.ls(P.inc.array[,,i],gtype.0_loh[i,]))),ncol=niter)
    }
    if (samp.type=="n"){
      Adj.1_inc_1=apply(P.inc.array,3,mis.adj.counts.ls,counts=c(vv.1,vw.1,ww.1))
      Adj.0_inc_1=apply(P.inc.array,3,mis.adj.counts.ls,counts=c(vv.0,vw.0,ww.0))
    }
      
    #find impossible values & eliminate those iterations
    # or iterations with probabilities so small the matrix won't invert later in program (rare)
    #low.val.flag<-(P.inc.array[1,1,]<1e-15)|(P.inc.array[2,1,]<1e-15)|(P.inc.array[3,1,]<1e-15)|(P.inc.array[1,2,]<1e-15)|(P.inc.array[2,2,]<1e-15)|(P.inc.array[3,2,]<1e-15)
    #Trap1.inc<-((Adj.1_inc_1[1,]<=0)|(Adj.1_inc_1[2,]<=0)|(Adj.0_inc_1[1,]<=0)|(Adj.0_inc_1[2,]<=0)|low.val.flag==1)
    Trap1.inc<-((Adj.1_inc_1[1,]<=0)|(Adj.1_inc_1[2,]<=0)|(Adj.0_inc_1[1,]<=0)|(Adj.0_inc_1[2,]<=0))
    Adj.1_inc_1<-Adj.1_inc_1[,Trap1.inc=="FALSE"]
    Adj.0_inc_1<-Adj.0_inc_1[,Trap1.inc=="FALSE"]
    P.inc.array<-P.inc.array[,,Trap1.inc=="FALSE"]
      
    niter=dim(Adj.1_inc_1)[2]
      
    # calculate prevalence of exposure (im NOTE THIS IS DIFFERENT FOR ASIAN POPS WHERE WE ASSUME NO PM)
    #   in cases and controls, accounting for sampling error

    
    PrevE.1<-rbeta(niter,Adj.1_inc_1[1,],Adj.1_inc_1[2,])
    PrevE.0<-rbeta(niter,Adj.0_inc_1[1,],Adj.0_inc_1[2,])
    #Occasionally prevalences are too close to 1 or 0 and it causes a problem 
    #  in the next few lines where. This replaces values too close to 1 or 0 
    #  with values that are a little further from 1 or 0
    PrevE.1[PrevE.1>.99999]=.99999
    PrevE.1[PrevE.1<.00001]=.00001
    PrevE.0[PrevE.0>.99999]=.99999
    PrevE.0[PrevE.0<.00001]=.00001
    
    PrevE.0<-cbind(PrevE.0,1-PrevE.0)
    PrevE.1<-cbind(PrevE.1,1-PrevE.1)
    
    
      
      
    P.array=P.inc.array
    # calculate Predictive Values of exposure classification in cases and controls
    # namin: PV_XX_xx_z: PV=predictive value;of XX= Phenotype; given xx=genotype;and given z=case/control
      
    PPV_1_inc<-array(unlist(lapply(1:niter,function(i)pv.fun(P.array[,,i],PrevE.1[i,]))),dim=c(3,2,niter))
    PPV_0_inc<-array(unlist(lapply(1:niter,function(i)pv.fun(P.array[,,i],PrevE.0[i,]))),dim=c(3,2,niter))
      
    #calculate the expected number of exposed cases and controls conditional on 
    #   Observed allele status
    if (samp.type=="t"){
      count_vv_1=t(rmultinomial(n=niter,size=vv.1_loh,prob=t(PPV_1_inc[1,,])))
      count_vw_1=t(rmultinomial(n=niter,size=vw.1_loh,prob=t(PPV_1_inc[2,,])))
      count_ww_1=t(rmultinomial(n=niter,size=ww.1_loh,prob=t(PPV_1_inc[3,,])))
        
      count_vv_0=t(rmultinomial(n=niter,size=vv.0_loh,t(PPV_0_inc[1,,])))
      count_vw_0=t(rmultinomial(n=niter,size=vw.0_loh,t(PPV_0_inc[2,,])))
      count_ww_0=t(rmultinomial(n=niter,size=ww.0_loh,t(PPV_0_inc[3,,])))
    }
    if (samp.type=="n"){
      count_vv_1=apply(t(PPV_1_inc[1,,]),1,rmultinom,n=1,size=vv.1)
      count_vw_1=apply(t(PPV_1_inc[2,,]),1,rmultinom,n=1,size=vw.1)
      count_ww_1=apply(t(PPV_1_inc[3,,]),1,rmultinom,n=1,size=ww.1)
        
      count_vv_0=apply(t(PPV_0_inc[1,,]),1,rmultinom,n=1,size=vv.0)
      count_vw_0=apply(t(PPV_0_inc[2,,]),1,rmultinom,n=1,size=vw.0)
      count_ww_0=apply(t(PPV_0_inc[3,,]),1,rmultinom,n=1,size=ww.0)
    }
      
    #aggregate across observe allele, to get total expected number of true allele
    im.1_inc=count_vv_1[1,]+count_vw_1[1,]+count_ww_1[1,]
    nm.1_inc=count_vv_1[2,]+count_vw_1[2,]+count_ww_1[2,]
    im.0_inc=count_vv_0[1,]+count_vw_0[1,]+count_ww_0[1,]
    nm.0_inc=count_vv_0[2,]+count_vw_0[2,]+count_ww_0[2,]
    
    #remove zeros
    trap2=(im.1_inc<=0)|(nm.1_inc<=0)|(im.0_inc<=0)|(nm.0_inc<=0)
    im.1_inc=im.1_inc[trap2=="FALSE"]
    nm.1_inc=nm.1_inc[trap2=="FALSE"]
    im.0_inc=im.0_inc[trap2=="FALSE"]
    nm.0_inc=nm.0_inc[trap2=="FALSE"]
    niter=length(im.1_inc)
      
    if (cohort==0){
      # calculate bias-adjusted OR (IM functional allele vs nm) with inc uncertainty
      res_im_vs_nm <- t(matrix(unlist(lapply(1:niter,function(i)crude.or.ci(im.1_inc[i],nm.1_inc[i],im.0_inc[i],nm.0_inc[i]))),nrow=4))
      or_im_vs_nm_r <- exp(rnorm(niter,log(res_im_vs_nm[,1]),res_im_vs_nm[,2]))
    }
      
    if (cohort==1){
      # calculate bias-adjusted OR (IM functional allele vs nm) with inc uncertainty
      res_im_vs_nm <- t(matrix(unlist(lapply(1:niter,function(i)crude.rr.ci(im.1_inc[i],nm.1_inc[i],im.0_inc[i],nm.0_inc[i]))),nrow=4))
      or_im_vs_nm_r <- exp(rnorm(niter,log(res_im_vs_nm[,1]),res_im_vs_nm[,2]))
    }
      
    if (samp.type=="t"){
      LOH<-quantile(or_1vs2_loh_r,c(.5,.025,.975))
      LOH_Inc<-quantile(or_im_vs_nm_r,c(.5,.025,.975))
      or1vs2<-rbind(res_1vs2[c(1,3,4)],LOH,LOH_Inc)
      rownames(or1vs2)<-c("Crude","LOH only","LOH & INC")
    }
    
    if (samp.type=="n"){
      LOH<-c("NA","NA","NA")
      LOH_Inc<-quantile(or_im_vs_nm_r,c(.5,.025,.975))
      or1vs2<-rbind(res_1vs2[c(1,3,4)],LOH,LOH_Inc)
      rownames(or1vs2)<-c("Crude","LOH only","INC")
    }
      
      
    dropped<-c(niter.init,niter)
    dropped[is.na(dropped)]=0
    names(dropped)=c("Initial","Final")
    return.list=list("Het"=or1vs2,"Drop"=dropped)
  }