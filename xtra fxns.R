#Assorted functions for cyp2d6 analysis
  
  
#function to compute adjusted cell counts using matrix method
mis.adj.counts<-function(sesp.matrix,counts){
  solve(sesp.matrix)%*%counts
}
    
    
#function to compute adjusted cell counts using matrix method
#using least squares for over-identified solution in Asian pop
mis.adj.counts.ls<-function(sesp.matrix,counts){
  adjusted.counts<-solve(t(sesp.matrix)%*%sesp.matrix)%*%t(sesp.matrix)%*%counts
  #LS solution can (will) result in more adjusted people than in "counts." next line rescales the adjusted counts
  adjusted.counts*sum(counts)/sum(adjusted.counts)
}
    
    
#function to compute pred values using matrix method
pv.fun<-function(Prob,Prev){
  solve(diag(as.vector(Prob%*%Prev)))%*%Prob%*%diag(Prev)
}
    
    
    
    
    
crude.or.ci<-function(a,b,c,d){
  fun.or<-a*d/(b*c)
  fun.or.se<-sqrt(1/a+1/b+1/c+1/d)
  fun.ci.low<-exp(log(fun.or)-1.96*fun.or.se)
  fun.ci.hi<-exp(log(fun.or)+1.96*fun.or.se)
  c(fun.or,fun.or.se,fun.ci.low,fun.ci.hi)
}
    
    
crude.rr.ci<-function(a,b,c,d){
  fun.rr<-(a/(a+c))/(b/(b+d))
  fun.rr.se<-sqrt(1/a+1/b-1/(a+c)-1/(b+d))
  fun.ci.low<-exp(log(fun.rr)-1.96*fun.rr.se)
  fun.ci.hi<-exp(log(fun.rr)+1.96*fun.rr.se)
  c(fun.rr,fun.rr.se,fun.ci.low,fun.ci.hi)
}
