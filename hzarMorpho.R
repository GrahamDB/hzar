source("hzarClasses.R");

logLikeCLTEstPopMean=function(muEst,muObs,varObs,nObs){
  #kappa=nObs/(2*varObs);
t.model<-sqrt(nObs/varObs)*(muObs-muEst);

  return( dt(t.model,df=nObs-1,log=TRUE)-dt(0,df=nObs-1,log=TRUE));
}


#log(kappa/pi)/2
## logLikeCLTEstPopMeanVar=function(muEst,varEst,muObs,varObs,nObs){
##   kappa=nObs/(2*varEst);
##   return(log(kappa/pi)/2-log(nObs/(pi*2*varObs))/2 -kappa*(muObs-muEst)^2);
## }

#Perhaps I should estimate the variance as proportianal to p(1-p).
#Given that I could consider hybrid zone as a mix of 2 distribution,
#with the mixing factor following a Bernoulli distribution, I would
#expect a p(1-p) to be involved.


#Make cline data frame
doCLTData1DPops<-function(distance,muObs,varObs,nSamples){
  if((length(distance) != length(muObs))     ||
     (length(distance) != length(nSamples)) ||
     (length(distance) != length(varObs))  ){
    stop("Distance, muObs, varObs and nSamples are not all of the same length.");
  }
  if(sum(nSamples<2)>0)
    stop("There must be at least two samples per population!");
  if(sum(nSamples<10)>0)
    warning("Some populations have less than 10 samples.  More samples would be nice.");
  ## Check for 0 variance singularity
   if(sum(varObs==0)>0){
    warning("Some of the population samples have a variance of 0.  Adding in estimated variance due to measurement error.");
    
    m.err<-5/3*10^-quantile(1+as.numeric(lapply(muObs,
             function(nV) {min((-1:12)[round(nV,digits=-1:12)==nV])}
                                                    )),probs=0.75)[[1]];
    ##For consistancy, error estimate should applied accross the board.
    varObs<-varObs+m.err*m.err;
  
  
  }
  ## if(sum(varObs==0)>0)
  ##   stop("varObs is 0 for one of the samples! varObs should at least be greater than the variance of the measurement error.");
  obj<-list(frame=data.frame(dist=distance,obsMean=muObs,obsVariance=varObs,n=nSamples));
  obj$model.LL <- function(model.func){
    muEst=model.func(obj$frame$dist);
##res<-numeric(length(pEst));
  ##  for(iter in 1:(length(pEst))){
      res<-logLikeCLTEstPopMeanVar(muEst=as.numeric(muEst),
                                      muObs=as.numeric(obj$frame$obsMean),
                                      varObs=as.numeric(obj$frame$obsVariance),
                                      N=as.numeric(obj$frame$n));
    ##}
                                                
    result<-sum(res);
    if(is.na(result))
      return(-1e8);
    return(result);
  }
class(obj)<-"clineSampleData1D";
  return(obj);
}

doCLTData1DRaw<-function(distance,traitValue){
  if(length(distance)!=length(traitValue))
    stop("Distance and traitValue vectors not of equal length.");
  ## determine population groups
  dist.group<-unique(distance);
  ## determine sample counts
  group.nSamp<-as.numeric(lapply(X=dist.group,
                                 function(a) sum(distance==a)));
  ## Make sure there are enough samples
  if(sum(group.nSamp<3)>0)
    stop("There are not enough samples in discrete populations.  If the data submitted is correct, you should either drop low sample populations or use one of the population sample interpolation methods.");
  if(sum(group.nSamp<10)>0)
    warning("There are very few samples in discrete populations.  You should consider dropping low sample populations or using one of the population sample interpolation methods.");
  
  ## Calculate sample means
  group.mean<-as.numeric(lapply(X=dist.group,
                                function(a)
                                mean(traitValue[distance==a])));
  ## Calculate sample variance
  group.var<-as.numeric(lapply(X=dist.group,
                                function(a)
                                var(traitValue[distance==a])));
  ## Check for 0 variance singularity
 if(sum(group.var==0)>0){
    warning("Some of the population samples have a variance of 0.  Adding in estimated variance due to measurement error.");
    m.err<-5/3*10^-quantile(1+as.numeric(lapply(traitValue,
             function(nV) {min((-1:12)[round(nV,digits=-1:12)==nV])}
                                                    )),probs=0.75)[[1]];
    ##For consistancy, error estimate should applied accross the board.
    group.var<-group.var+m.err*m.err;
  
  
  }
  return(data.frame(dist=dist.group,
                    mu=group.mean,
                    sigma2=group.var,
                    nSamp=group.nSamp));
}
