#source("hzarClasses.R");

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
      res<-logLikeCLTEstPopMean(muEst=as.numeric(muEst),
                                      muObs=as.numeric(obj$frame$obsMean),
                                      varObs=as.numeric(obj$frame$obsVariance),
                                      nObs=as.numeric(obj$frame$n));
    ##}
                                                
    result<-sum(res);
    if(is.na(result))
      return(-1e8);
    return(result);
  }
class(obj)<-c("clineSampleData1DCLT","hzar.obsData");
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

hzar.meta.CLTnA =
  list(
       prior=function(center,width,xMin,xMax){
  return(0); },
      func=function(center,width,xMin,xMax){
         pCline<- function(x) {
           u <- (x - center) * 4/width
           return(xMin+(xMax-xMin)* (1/(1+ exp(-u)))) }
         
         return(pCline)
       },
       req=function(center,width,xMin,xMax){
         return(xMin <xMax & width>0)},
       parameterTypes=CLINEPARAMETERS[c("center","width","xMin","xMax")]
       );
hzar.meta.CLTnD =
  list(
       prior=function(center,width,xMin,xMax){
  return(0); },
      func=function(center,width,xMin,xMax){
         pCline<- function(x) {
           u <- (x - center) * -4/width
           return(xMin+(xMax-xMin)* (1/(1+ exp(-u)))) }
         
         return(pCline)
       },
       req=function(center,width,xMin,xMax){
         return(xMin <xMax &width>0)},
       parameterTypes=CLINEPARAMETERS[c("center","width","xMin","xMax")]
       );
class(hzar.meta.CLTnA)<-"clineMetaModel";
class(hzar.meta.CLTnD)<-"clineMetaModel";
hzar.meta.CLTrA =
  list(req= function(center,width,xMin,xMax,deltaR,tauR)
       {
         return(width>0 & deltaR>=0 &
                xMin <xMax &
                tauR>=0 & tauR<=1 )
       },
       prior=function(center,width,xMin,xMax,deltaR,tauR){
  return(0); },
      func=function(center,width,xMin,xMax,deltaR,tauR)
       {
         gamma=4/width;
        # tail.LO=meta.tail.lower(gamma=gamma,d1=deltaL,tau1=tauL);
         tail.HI=meta.tail.upper(gamma=gamma,d2=deltaR,tau2=tauR);
         clineComposite=
           meta.cline.func.upStep(center=center,
                                    direction=1,
                                    gamma=gamma,
                              #      lowerTail=tail.LO,
                                    upperTail=tail.HI);
         return(meta.cline.func.pScale(xMin,xMax,clineComposite));
       },
       parameterTypes=CLINEPARAMETERS[c("center","width","xMin","xMax","deltaR","tauR")]
       );
hzar.meta.CLTrD =
  list(req= function(center,width,xMin,xMax,deltaR,tauR)
       {
         return(width>0 & deltaR>=0 &
                xMin <xMax &
                tauR>=0 & tauR<=1 )
       },
       prior=function(center,width,xMin,xMax,deltaR,tauR){
  return(0); },
      func=function(center,width,xMin,xMax,deltaR,tauR)
       {
         gamma=4/width;
         tail.LO=meta.tail.lower(gamma=gamma,d1=deltaR,tau1=tauR);
        # tail.HI=meta.tail.upper(gamma=gamma,d2=deltaL,tau2=tauL);
         clineComposite=
           meta.cline.func.lowStep(center=center,
                                    direction=-1,
                                    gamma=gamma,
                                    lowerTail=tail.LO);
         return(meta.cline.func.pScale(xMin,xMax,clineComposite));
       },
       parameterTypes=CLINEPARAMETERS[c("center","width","xMin","xMax","deltaR","tauR")]
       );
class(hzar.meta.rtail.scaled.ascending)<-"clineMetaModel";
class(hzar.meta.rtail.scaled.descending)<-"clineMetaModel";


setupCLTCenterClineParameters<-function(myModel,scaling,x=NULL,y=NULL) {
   if(scaling=="fixed") {
    attr(myModel$parameterTypes$xMin,"fixed")<-TRUE;
    attr(myModel$parameterTypes$xMax,"fixed")<-TRUE;
    if(!is.null(y)){
      myModel$parameterTypes$xMin$val<-min(y);
      myModel$parameterTypes$xMax$val<-max(y);
    }
  } else if(scaling=="free") {
    attr(myModel$parameterTypes$xMin,"fixed")<-FALSE;
    attr(myModel$parameterTypes$xMax,"fixed")<-FALSE;
    if(!is.null(y)){
      myModel$parameterTypes$xMin$val<-min(y);
      myModel$parameterTypes$xMax$val<-max(y);
      hzar.suggestionFunc1D$xMin(x,y)->junk;
    attr(myModel$parameterTypes$xMin,"limit.lower")<-junk[[1]];
    attr(myModel$parameterTypes$xMin,"limit.upper")<-junk[[2]];
      hzar.suggestionFunc1D$xMax(x,y)->junk;
    attr(myModel$parameterTypes$xMax,"limit.lower")<-junk[[1]];
    attr(myModel$parameterTypes$xMax,"limit.upper")<-junk[[2]];
    }
  } else {
    stop(paste("Scaling type",scaling,"unrecignized. Please use none, fixed, or free."));
  }
   pTnames<-names(myModel$parameterTypes);
   
  if(!is.null(x)){
    qX<-quantile(x,probs=c(0.25,0.5,0.75));
    myModel$parameterTypes$center$val<-qX[[2]]; 
    myModel$parameterTypes$width$val<-qX[[3]]-qX[[1]];
    hzar.suggestionFunc1D$center(x,y)->junk;
    attr(myModel$parameterTypes$center,"limit.lower")<-junk[[1]];
    attr(myModel$parameterTypes$center,"limit.upper")<-junk[[2]];
    
    hzar.suggestionFunc1D$width(x,y)->junk;
    attr(myModel$parameterTypes$width,"limit.lower")<-junk[[1]];
    attr(myModel$parameterTypes$width,"limit.upper")<-junk[[2]];

    index<-"deltaR";
    if(index %in% pTnames){
      hzar.suggestionFunc1D[[index]](x,y)->junk;
      attr(myModel$parameterTypes[[index]],"limit.lower")<-junk[[1]];
      attr(myModel$parameterTypes[[index]],"limit.upper")<-junk[[2]];
    }
    index<-"deltaM";
    if(index %in% pTnames){
      hzar.suggestionFunc1D[[index]](x,y)->junk;
      attr(myModel$parameterTypes[[index]],"limit.lower")<-junk[[1]];
      attr(myModel$parameterTypes[[index]],"limit.upper")<-junk[[2]];
    }
    index<-"deltaL";
    if(index %in% pTnames){
      hzar.suggestionFunc1D[[index]](x,y)->junk;
      attr(myModel$parameterTypes[[index]],"limit.lower")<-junk[[1]];
      attr(myModel$parameterTypes[[index]],"limit.upper")<-junk[[2]];
    }
    
    
  }
   index<-"tauR";
   if(index %in% pTnames){
     hzar.suggestionFunc1D[[index]](x,y)->junk;
     attr(myModel$parameterTypes[[index]],"limit.lower")<-junk[[1]];
     attr(myModel$parameterTypes[[index]],"limit.upper")<-junk[[2]];
   }
   index<-"tauM";
   if(index %in% pTnames){
     hzar.suggestionFunc1D[[index]](x,y)->junk;
     attr(myModel$parameterTypes[[index]],"limit.lower")<-junk[[1]];
     attr(myModel$parameterTypes[[index]],"limit.upper")<-junk[[2]];
   }
   index<-"tauL";
   if(index %in% pTnames){
     hzar.suggestionFunc1D[[index]](x,y)->junk;
     attr(myModel$parameterTypes[[index]],"limit.lower")<-junk[[1]];
     attr(myModel$parameterTypes[[index]],"limit.upper")<-junk[[2]];
   }
    
  return(myModel);
}

makeCline1DCLT<- function(data=NULL,scaling="free",tails="none",direction=NULL){
  if(identical(tolower(tails),"none")){
    return(buildCline1D(data,scaling,direction,
                        hzar.meta.CLTnA,hzar.meta.CLTnD));
 ##  }else if(identical(tolower(tails),"both")) {
##     return(makeTailedCline1D(data,scaling,direction));
   }else if(identical(tolower(tails),"right")) {
     myRightCline<-buildCline1D(data,scaling,direction,
                                hzar.meta.CLTrA, hzar.meta.CLTrD);
     
     attr(myRightCline,"tails")<-"right";
     return(myRightCline);
##   }else if(identical(tolower(tails),"left")) {
## myLeftCline<-buildCline1D(data,scaling,direction,
##                               hzar.meta.ltail.scaled.ascending,
##                               hzar.meta.ltail.scaled.descending);
##   attr(myLeftCline,"tails")<-"left";
##     return(myLeftCline);
##   }else if(identical(tolower(tails),"mirror")) {
## myMirrorCline<-buildCline1D(data,scaling,direction,
##                               hzar.meta.mtail.scaled.ascending,
##                               hzar.meta.mtail.scaled.descending);
##   attr(myMirrorCline,"tails")<-"mirror";
##     return(myMirrorCline);
  }
  stop(paste("Cline with",tails,"tail(s) not available.")); 
}


