## So I need this file to handle the setup of the data and the model.
source("modelEquationPrototypeGrouping.R");
require(MCMCpack);

## First, i will define the sample level likelihood Function.

sampleLikelihoodMolecularPop=function(pEst,pObs,N) {
  p1=N*log(pEst);
  p0=N*log(1-pEst);
  lP=((1-pObs)*log((1-pEst)/(1-pObs))+pObs*log(pEst/pObs))*N;
  tmp=ifelse(pObs==1,p1,lP);
  res=ifelse(pObs==0,p0,tmp);
  return(res);
}

## Second, i will setup the data.
doMolecularData1DPops<-function(distance,pObs,nSamples){
  if((length(distance) != length(pObs))     ||
     (length(distance) != length(nSamples)) ||
     (length(pObs) != length(nSamples))  ){
    stop("Distance, pObs and nSamples are not of the same length.");
  }
  obj<-list(frame=data.frame(dist=distance,obsFreq=pObs,n=nSamples));
  obj$model.LL <- function(model.func){
    pEst=model.func(obj$frame$dist);
res<-numeric(length(pEst));
    for(iter in 1:(length(pEst))){
      res[[iter]]<-sampleLikelihoodMolecularPop(pEst=pEst[[iter]],
                                                pObs=obj$frame$obsFreq[[iter]],
                                                N=obj$frame$n[[iter]]);
    }
                                                
    result<-sum(res);
    if(is.na(result))
      return(-1e8);
    return(result);
  }
class(obj)<-"clineSampleData1D";
  return(obj);
}

##formate of parameter data frame?
mkParam<-function(name,value,weight,lower,upper,isFixed=FALSE){
  objCParamC<-list(val=value,w=weight);
  class(objCParamC)<-"clineParameter";
  attr(objCParamC,"param")<-name
  attr(objCParamC,"fixed")<-isFixed;
  attr(objCParamC,"limit.lower")<- lower;
  attr(objCParamC,"limit.upper")<- upper;
  return(objCParamC);
}
CLINEPARAMETERS<-list(center=mkParam("center",10,1.5,-1e8,1e8),
                      width =mkParam("width", 10,1.5,1e-8,1e8),
                      pMin  =mkParam("pMin",  0, 1.1,   0,  1),
                      pMax  =mkParam("pMax",  1, 1.1,   0,  1),
                      deltaL=mkParam("deltaL",1,1.5,  1e-8, 1e8),
                      deltaR=mkParam("deltaR",1,1.5,  1e-8, 1e8),
                      deltaM=mkParam("deltaM",1,1.5,  1e-8, 1e8),
                      tauL  =mkParam("tauL", 0.5,1.1,  0,   1),
                      tauR  =mkParam("tauR", 0.5,1.1,  0,   1),
                      tauM  =mkParam("tauM", 0.5,1.1,  0,   1));

splitParameters<-function(paramList){
  param.fixed=list();
  param.free.init =list();
  param.free.weight =list();
  param.free.lower =list();
  param.free.upper =list();
  for(iter in seq(along=paramList)){
    if(class(paramList[[iter]])!="clineParameter")
      stop(paste("Item",paramList[[iter]],"not a cline parameter!"));
    parName<-attr(paramList[[iter]],"param");
    parValue<-paramList[[iter]]$val;
    
    if(attr(paramList[[iter]],"fixed")){
      param.fixed[[parName]]<-parValue;
    }else{
      param.free.init[[parName]]<-parValue;
      param.free.weight[[parName]]<-paramList[[iter]]$w;
      param.free.lower[[parName]]<-attr(paramList[[iter]],"limit.lower");
      param.free.upper[[parName]]<-attr(paramList[[iter]],"limit.upper");
    }
  }
  param<-list(init=param.free.init,tune=param.free.weight,
              lower=param.free.lower,upper=param.free.upper,
              fixed=param.fixed);
  return(param);
}

## ## format of the cline meta model
## objCMeta<-list();
## objCMeta$parameterTypes<-list("center","width");
## objCMeta$req<-function(center,width);
## objCMeta$func<-function(center,width);
## class(objCMeta)<-"clineMetaModel";

cline.unscale <- function(cline.data, pMin,pMax){

  cline.data$obsFreq <- (cline.data$obsFreq - pMin)/(pMax-pMin);
  return(cline.data);
}

cline.logit <- function(cline.data){
  attach(cline.data);
  cline.data$obsFreq <- log( obsFreq / ( 1- obsFreq));
  detach();
  return(cline.data);
}

cline.samp.dist <- function(cline.data) {
  dist.pool <- unique(sort(as.numeric(cline.data$dist)));
  dist.init= quantile(dist.pool, probs=c(0.25,0.75));
  potLeftOut=length(dist.pool[dist.pool <= dist.init[[1]]]);
  potRightOut=length(dist.pool[dist.pool >= dist.init[[2]]])-1;
  potInner=length(dist.pool[(dist.pool < dist.init[[2]])&(dist.pool > dist.init[[1]])]);
##   print(dist.init);
##   print(potLeftOut);
##   print(potRightOut);
##   print(potInner);
  dist.len<-length(dist.pool);
  distLeftSeries<-numeric(1);
  distLeftSeries[1]<-dist.init[[1]];
  distRightSeries<-numeric(1);
  
  distRightSeries[1]<-dist.init[[2]];
  if(potLeftOut>2)  distLeftSeries  <- c(dist.pool[2:potLeftOut],
                                         distLeftSeries);
  if(potRightOut>1) distRightSeries <- c(distRightSeries        ,
                                         dist.pool[(dist.len-potRightOut):(dist.len-1)]);
  if(potInner > 1) {
    distLeftSeries  <- c(distLeftSeries,
                         dist.pool[potLeftOut+1]);
    distRightSeries <- c(dist.pool[dist.len-potRightOut-1],
                         distRightSeries);
  }
  return(list(unique(distLeftSeries),unique(distRightSeries)));
}

cssp<-function(cline.data,do.tail="none",pMinS=NULL,pMaxS=NULL){
  do.left<-FALSE;
  do.right<-FALSE;
  do.mirror<-FALSE;
 
  keyValues=c("center","width","pMin","pMax")
    
  attach(cline.data);
  pOQ<-quantile(unique(sort(obsFreq)),probs=c(0,0.1,0.9,1));
  
  detach(cline.data);
  if(is.null(pMinS)){
    pMinS<-as.numeric(unique(c(0,pOQ[[1]]/2,pOQ[[1]],pOQ[[2]])));
    if(length(pMaxS)<3) {
      pMinS<-as.numeric(unique(c(pMinS,(pOQ[[1]]+pOQ[[2]])/2,0.01+pOQ[[2]])));
    }
  } else {
    pMinS<-as.numeric(pMinS);
    if(length(pMinS)==0) {
      warning("Malformed pMinS, using 0 instead");
      pMinS<-as.numeric(c(0));
    }
  }
  if(is.null(pMaxS)){
    pMaxS<-as.numeric(unique(c(1,(1+pOQ[[4]])/2,pOQ[[4]],pOQ[[3]])));
    if(length(pMaxS)<3) {
      pMaxS<-as.numeric(unique(c(pMaxS,(pOQ[[3]]+pOQ[[4]])/2,pOQ[[3]]-0.01)));
    }
  } else {
    pMaxS<-as.numeric(pMaxS);
    if(length(pMaxS)==0) {
      warning("Malformed pMaxS, using 1 instead");
      pMaxS<-as.numeric(c(1));
    }
  }
  
  bothDist<-cline.samp.dist(cline.data);
  
  if(identical(tolower(do.tail),"mirror")){
    do.left<-TRUE;
    do.right<-TRUE;
    do.mirror<-TRUE;
    result<-expand.grid(pMin=pMinS,pMax=pMaxS,dist.left=bothDist[[1]],dist.right=bothDist[[2]]);
    keyValues=c(keyValues,"deltaM","tauM");
  } else if (identical(tolower(do.tail),"both")) {
    do.left<-TRUE;
    do.right<-TRUE;
    result<-expand.grid(pMin=pMinS,pMax=pMaxS,dist.left=bothDist[[1]],dist.right=bothDist[[2]]);
    keyValues=c(keyValues,"deltaL","tauL","deltaR","tauR");
  } else if (identical(tolower(do.tail),"left")) {
    do.left<-TRUE;
    result<-expand.grid(pMin=pMinS,pMax=pMaxS,dist.left=bothDist[[1]]);
    keyValues=c(keyValues,"deltaL","tauL");
  } else if (identical(tolower(do.tail),"right")) {
    do.right<-TRUE;
    result<-expand.grid(pMin=pMinS,pMax=pMaxS,dist.right=bothDist[[2]]);
    keyValues=c(keyValues,"deltaR","tauR");
  } else {
    result<-expand.grid(pMin=pMinS,pMax=pMaxS)
  }
  
  nThetas<-dim(result)[[1]];
  result$center<-numeric(nThetas);
  result$width <-numeric(nThetas);
  result$deltaL <-numeric(nThetas);
  result$tauL <-numeric(nThetas);
  result$deltaR <-numeric(nThetas);
  result$tauR <-numeric(nThetas);
  result$deltaM <-numeric(nThetas);
  result$tauM <-numeric(nThetas);
  for(iter in 1:nThetas){
    tPMin<- result$pMin[[iter]];
    tPMax<- result$pMax[[iter]];
##     myClineData=subset(cline.data,
##       (cline.data$obsFreq >tPMin)&
##       (cline.data$obsFreq <tPMax));
    myClineData<-cline.data;
    myClineData$obsFreq[ myClineData$obsFreq< (1e-4+tPMin) ] <-
      rep( 1e-4+tPMin , length(which( myClineData$obsFreq< (1e-4+tPMin))));
    myClineData$obsFreq[ myClineData$obsFreq> (-1e-4+tPMax )]<-
      rep( tPMax-1e-4 , length(which( myClineData$obsFreq > (tPMax-1e-4)))); #=  0.9999*tPMax ;
## print(myClineData);
    myLogitData<-cline.logit(cline.unscale(myClineData,tPMin,tPMax));
    myFullLogitData<-myLogitData;
    if(do.left){
      myLeftData=subset(myLogitData,myLogitData$dist <= result$dist.left[[iter]]);
      myLogitData=subset(myLogitData,myLogitData$dist >= result$dist.left[[iter]]);
    }
    if(do.right){
      myRightData=subset(myLogitData,myLogitData$dist >= result$dist.right[[iter]]);
      myLogitData=subset(myLogitData,myLogitData$dist <= result$dist.right[[iter]]);
    }
    lmEst<-NULL;
    try(lmEst<-lm(obsFreq~dist,data=myLogitData));
    if(is.null(lmEst)) {
      try(lmEst<-lm(obsFreq~dist,data=myFullLogitData));
      if(is.null(lmEst)) {
        result$width[[iter]]<-NA;
        next;
      }
    }
    lambda=lmEst$coefficients[[2]]
    result$width[[iter]]<- abs(0.25/lambda);
    result$center[[iter]]<- -lmEst$coefficients[[1]]/lmEst$coefficients[[2]];
   
    if(do.mirror){
      result$center[[iter]]<- ( result$dist.left[[iter]]+ result$dist.right[[iter]])/2;
    }
    if(do.left){
       if(result$center[[iter]]< result$dist.left[[iter]])
      result$center[[iter]]=result$dist.left[[iter]];
      lmEst<-NULL;
      try(lmEst<-lm(obsFreq~dist,data=myLeftData));
      if(is.null(lmEst)) {
        result$tauL[[iter]]<-0; 
      }else{
        propose.tau<- lmEst$coefficients[[2]]/lambda;
        propose.tau=ifelse(propose.tau<0,1e-4+propose.tau*1e-5,propose.tau);
        propose.tau=ifelse(propose.tau>=1,1+(propose.tau-1)*1e-5-1e-4,propose.tau);
        propose.tau=ifelse(propose.tau<0,0,propose.tau);
        propose.tau=ifelse(propose.tau>1,1,propose.tau);
        result$tauL[[iter]]<-propose.tau;
        
      }
      result$deltaL[[iter]]<- result$center[[iter]]- result$dist.left[[iter]];
    }
    if(do.right){
      
    if(result$center[[iter]]> result$dist.right[[iter]])
      result$center[[iter]]=result$dist.right[[iter]];
    
      lmEst<-NULL;
      try(lmEst<-lm(obsFreq~dist,data=myRightData));
      if(is.null(lmEst)) {
        result$tauR[[iter]]<-0;
      }else{
        propose.tau<- lmEst$coefficients[[2]]/lambda;
        propose.tau=ifelse(propose.tau<0,1e-4+propose.tau*1e-5,propose.tau);
        propose.tau=ifelse(propose.tau>=1,1+(propose.tau-1)*1e-5-1e-4,propose.tau);
        propose.tau=ifelse(propose.tau<0,0,propose.tau);
        propose.tau=ifelse(propose.tau>1,1,propose.tau);
        result$tauR[[iter]]<- propose.tau;
      }
      result$deltaR[[iter]]<-result$dist.right[[iter]]- result$center[[iter]]
    }
    if(do.mirror){
      tau.mirror=( result$tauL[[iter]]+ result$tauR[[iter]])/2;
      result$tauR[[iter]]<-tau.mirror;
      result$tauL[[iter]]<-tau.mirror;
      result$tauM[[iter]]<-tau.mirror;
      result$deltaM[[iter]]<- result$deltaL[[iter]]
    }
    
  }
  
  result=result[,keyValues]
  result=na.omit(result);
  return(result);
} 

hzar.meta.simple.scaled.ascending =
  list(
       func=function(center,width,pMin,pMax){
         pCline<- function(x) {
           u <- (x - center) * 4/width
           return(pMin+(pMax-pMin)* (1/(1+ exp(-u)))) }
         
         return(pCline)
       },
       req=function(center,width,pMin,pMax){
         return(pMin>=0 & pMax<=1 & pMin<pMax & width>0)},
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax")]
       );
hzar.meta.simple.scaled.descending =
  list(
       func=function(center,width,pMin,pMax){
         pCline<- function(x) {
           u <- (x - center) * -4/width
           return(pMin+(pMax-pMin)* (1/(1+ exp(-u)))) }
         
         return(pCline)
       },
       req=function(center,width,pMin,pMax){
         return(pMin>=0 & pMax<=1 & pMin<pMax & width>0)},
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax")]
       );
class(hzar.meta.simple.scaled.ascending)<-"clineMetaModel";
class(hzar.meta.simple.scaled.descending)<-"clineMetaModel";
hzar.meta.tailed.scaled.ascending =
  list(req= function(center,width,pMin,pMax,deltaL,tauL,deltaR,tauR)
       {
         return(width>0 & deltaL>=0 & deltaR>=0 &
                pMin>=0 & pMax<=1 & pMin<pMax & 
                tauL>=0 & tauL<=1 &
                tauR>=0 & tauR<=1 )
       },
       func=function(center,width,pMin,pMax,deltaL,tauL,deltaR,tauR)
       {
         gamma=4/width;
         tail.LO=meta.tail.lower(gamma=gamma,d1=deltaL,tau1=tauL);
         tail.HI=meta.tail.upper(gamma=gamma,d2=deltaR,tau2=tauR);
         clineComposite=
           meta.cline.func.stepBoth(center=center,
                                    direction=1,
                                    gamma=gamma,
                                    lowerTail=tail.LO,
                                    upperTail=tail.HI);
         return(meta.cline.func.pScale(pMin,pMax,clineComposite));
       },
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax","deltaL","tauL","deltaR","tauR")]
       );
hzar.meta.tailed.scaled.descending =
  list(req= function(center,width,pMin,pMax,deltaL,tauL,deltaR,tauR)
       {
         return(width>0 & deltaL>=0 & deltaR>=0 &
                pMin>=0 & pMax<=1 & pMin<pMax & 
                tauL>=0 & tauL<=1 &
                tauR>=0 & tauR<=1 )
       },
       func=function(center,width,pMin,pMax,deltaL,tauL,deltaR,tauR)
       {
         gamma=4/width;
         tail.LO=meta.tail.lower(gamma=gamma,d1=deltaR,tau1=tauR);
         tail.HI=meta.tail.upper(gamma=gamma,d2=deltaL,tau2=tauL);
         clineComposite=
           meta.cline.func.stepBoth(center=center,
                                    direction=-1,
                                    gamma=gamma,
                                    lowerTail=tail.LO,
                                    upperTail=tail.HI);
         return(meta.cline.func.pScale(pMin,pMax,clineComposite));
       },
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax","deltaL","tauL","deltaR","tauR")]
       );
class(hzar.meta.tailed.scaled.ascending)<-"clineMetaModel";
class(hzar.meta.tailed.scaled.descending)<-"clineMetaModel";

hzar.meta.mtail.scaled.descending =
  list(req= function(center,width,pMin,pMax,deltaM,tauM)
       {
         return(width>0 & deltaM>=0 &
                pMin>=0 & pMax<=1 & pMin<pMax & 
                tauM>=0 & tauM<=1 )
       },
       func=function(center,width,pMin,pMax,deltaM,tauM)
       {
         gamma=4/width;
         tail.LO=meta.tail.lower(gamma=gamma,d1=deltaM,tau1=tauM);
         tail.HI=meta.tail.upper(gamma=gamma,d2=deltaM,tau2=tauM);
         clineComposite=
           meta.cline.func.stepBoth(center=center,
                                    direction=-1,
                                    gamma=gamma,
                                   lowerTail=tail.LO,
                                    upperTail=tail.HI);
         return(meta.cline.func.pScale(pMin,pMax,clineComposite));
       },
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax","deltaM","tauM")]
       );
hzar.meta.mtail.scaled.ascending =
  list(req= function(center,width,pMin,pMax,deltaM,tauM)
       {
         return(width>0 & deltaM>=0 &
                pMin>=0 & pMax<=1 & pMin<pMax & 
                tauM>=0 & tauM<=1 )
       },
       func=function(center,width,pMin,pMax,deltaM,tauM)
       {
         gamma=4/width;
         tail.LO=meta.tail.lower(gamma=gamma,d1=deltaM,tau1=tauM);
        tail.HI=meta.tail.upper(gamma=gamma,d2=deltaM,tau2=tauM);
         clineComposite=
           meta.cline.func.stepBoth(center=center,
                                    direction=1,
                                    gamma=gamma,
                                    lowerTail=tail.LO,
                                    upperTail=tail.HI);
         return(meta.cline.func.pScale(pMin,pMax,clineComposite));
       },
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax","deltaM","tauM")]
       );
class(hzar.meta.mtail.scaled.ascending)<-"clineMetaModel";
class(hzar.meta.mtail.scaled.descending)<-"clineMetaModel";
hzar.meta.ltail.scaled.descending =
  list(req= function(center,width,pMin,pMax,deltaL,tauL)
       {
         return(width>0 & deltaL>=0 &
                pMin>=0 & pMax<=1 & pMin<pMax & 
                tauL>=0 & tauL<=1 )
       },
       func=function(center,width,pMin,pMax,deltaL,tauL)
       {
         gamma=4/width;
        # tail.LO=meta.tail.lower(gamma=gamma,d1=deltaL,tau1=tauL);
         tail.HI=meta.tail.upper(gamma=gamma,d2=deltaL,tau2=tauL);
         clineComposite=
           meta.cline.func.upStep(center=center,
                                    direction=-1,
                                    gamma=gamma,
                              #      lowerTail=tail.LO,
                                    upperTail=tail.HI);
         return(meta.cline.func.pScale(pMin,pMax,clineComposite));
       },
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax","deltaL","tauL")]
       );
hzar.meta.ltail.scaled.ascending =
  list(req= function(center,width,pMin,pMax,deltaL,tauL)
       {
         return(width>0 & deltaL>=0 &
                pMin>=0 & pMax<=1 & pMin<pMax & 
                tauL>=0 & tauL<=1 )
       },
       func=function(center,width,pMin,pMax,deltaL,tauL)
       {
         gamma=4/width;
         tail.LO=meta.tail.lower(gamma=gamma,d1=deltaL,tau1=tauL);
        # tail.HI=meta.tail.upper(gamma=gamma,d2=deltaL,tau2=tauL);
         clineComposite=
           meta.cline.func.lowStep(center=center,
                                    direction=1,
                                    gamma=gamma,
                                    lowerTail=tail.LO);
         return(meta.cline.func.pScale(pMin,pMax,clineComposite));
       },
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax","deltaL","tauL")]
       );
class(hzar.meta.ltail.scaled.ascending)<-"clineMetaModel";
class(hzar.meta.ltail.scaled.descending)<-"clineMetaModel";
hzar.meta.rtail.scaled.ascending =
  list(req= function(center,width,pMin,pMax,deltaR,tauR)
       {
         return(width>0 & deltaR>=0 &
                pMin>=0 & pMax<=1 & pMin<pMax & 
                tauR>=0 & tauR<=1 )
       },
       func=function(center,width,pMin,pMax,deltaR,tauR)
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
         return(meta.cline.func.pScale(pMin,pMax,clineComposite));
       },
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax","deltaR","tauR")]
       );
hzar.meta.rtail.scaled.descending =
  list(req= function(center,width,pMin,pMax,deltaR,tauR)
       {
         return(width>0 & deltaR>=0 &
                pMin>=0 & pMax<=1 & pMin<pMax & 
                tauR>=0 & tauR<=1 )
       },
       func=function(center,width,pMin,pMax,deltaR,tauR)
       {
         gamma=4/width;
         tail.LO=meta.tail.lower(gamma=gamma,d1=deltaR,tau1=tauR);
        # tail.HI=meta.tail.upper(gamma=gamma,d2=deltaL,tau2=tauL);
         clineComposite=
           meta.cline.func.lowStep(center=center,
                                    direction=-1,
                                    gamma=gamma,
                                    lowerTail=tail.LO);
         return(meta.cline.func.pScale(pMin,pMax,clineComposite));
       },
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax","deltaR","tauR")]
       );
class(hzar.meta.rtail.scaled.ascending)<-"clineMetaModel";
class(hzar.meta.rtail.scaled.descending)<-"clineMetaModel";
setupMoleCenterClineParameters<-function(myModel,scaling,x=NULL,y=NULL) {
  if(scaling=="none"){
    attr(myModel$parameterTypes$pMin,"fixed")<-TRUE;
    attr(myModel$parameterTypes$pMax,"fixed")<-TRUE;
    myModel$parameterTypes$pMin$val<-0;
    myModel$parameterTypes$pMax$val<-1;
    
  } else if(scaling=="fixed") {
    attr(myModel$parameterTypes$pMin,"fixed")<-TRUE;
    attr(myModel$parameterTypes$pMax,"fixed")<-TRUE;
    if(!is.null(y)){
      myModel$parameterTypes$pMin$val<-min(y);
      myModel$parameterTypes$pMax$val<-max(y);
    }
  } else if(scaling=="free") {
    attr(myModel$parameterTypes$pMin,"fixed")<-FALSE;
    attr(myModel$parameterTypes$pMax,"fixed")<-FALSE;
    if(!is.null(data)){
      myModel$parameterTypes$pMin$val<-min(y);
      myModel$parameterTypes$pMax$val<-max(y);
    }
  } else {
    stop(paste("Scaling type",scaling,"unrecignized. Please use none, fixed, or free."));
  }
  if(!is.null(x)){
    qX<-quantile(x,probs=c(0.25,0.5,0.75));
    myModel$parameterTypes$center$val<-qX[[2]]; 
    myModel$parameterTypes$width$val<-qX[[3]]-qX[[1]];
  }
  return(myModel);
}

buildCline1D<-function(data,scaling,direction,ascending.cline,descending.cline){
  if(is.null(direction) ){
    if(is.null(data))
      stop("Either provide observation data, or set direction.");
    if(class(data)=="clineSampleData1D"){
      if(mean(data$frame$dist) <
         (sum(data$frame$dist*data$frame$obsFreq)/sum(data$frame$obsFreq))){
        direction="ascending";
      } else {
        direction="descending";
      }
    } else {
      stop("Class of data unknown.");
    }
  }
  if(identical(tolower(direction),"ascending")){
    myModel=ascending.cline;
  }else if(identical(tolower(direction),"descending")){
    myModel=descending.cline;
  }else {
    stop(paste("Direction",direction,"not recognized. Specify either ascending or descending."));
  }
  myModel=setupMoleCenterClineParameters(myModel,scaling,data$frame$dist,data$frame$obsFreq);
  attr(myModel,"tails")<-"none";
  return(myModel);
}

makeSimpleCline1D<-function(data=NULL,scaling="none",direction=NULL){
  return(buildCline1D(data,scaling,direction,
                      hzar.meta.simple.scaled.ascending,
                      hzar.meta.simple.scaled.descending));
}

makeTailedCline1D<-function(data=NULL,scaling="none",direction=NULL){
  myTailedModel<-buildCline1D(data,scaling,direction,
                              hzar.meta.tailed.scaled.ascending,
                              hzar.meta.tailed.scaled.descending);
  Width<-myTailedModel$parameterTypes$width$val;
  attr(myTailedModel,"tails")<-"both";
##   myTailedModel$parameterTypes$deltaL$val<-Width;
##   myTailedModel$parameterTypes$deltaR$val<-Width;
  return(myTailedModel);
}

makeCline1D<- function(data=NULL,scaling="none",tails="none",direction=NULL){
  if(identical(tolower(tails),"none")){
    return(makeSimpleCline1D(data,scaling,direction));
  }else if(identical(tolower(tails),"both")) {
    return(makeTailedCline1D(data,scaling,direction));
  }else if(identical(tolower(tails),"right")) {
myRightCline<-buildCline1D(data,scaling,direction,
                              hzar.meta.rtail.scaled.ascending,
                              hzar.meta.rtail.scaled.descending);
  attr(myRightCline,"tails")<-"right";
    return(myRightCline);
  }else if(identical(tolower(tails),"left")) {
myLeftCline<-buildCline1D(data,scaling,direction,
                              hzar.meta.ltail.scaled.ascending,
                              hzar.meta.ltail.scaled.descending);
  attr(myLeftCline,"tails")<-"left";
    return(myLeftCline);
  }else if(identical(tolower(tails),"mirror")) {
myMirrorCline<-buildCline1D(data,scaling,direction,
                              hzar.meta.mtail.scaled.ascending,
                              hzar.meta.mtail.scaled.descending);
  attr(myMirrorCline,"tails")<-"mirror";
    return(myMirrorCline);
  }
  stop(paste("Cline with",tails,"tail(s) not available.")); 
}

## ## format of the cline model Frame
## objCFrame<-list();
## objCFrame$meta.model <- objCMeta;
## objCFrame$data <-objSampleData;  ##Keep Seperate?
## objCFrame$parameters <-objParameterDataFrame;
## objCFrame$mcmc <- mcmcObj;
## class(objCframe)<-"clineModelFrame";

fitClineModel <- function(model,sampleData, verbose=10000,
                          burnin=1000,mcmc=1e6,thin=100,
                          seed=list(NA,1),rejectLL=-1e8){
  myClineFrame=list(data=sampleData);
##  myModel<-model;
  attr(myClineFrame,"mcmc.count")<-0;
  if(identical(is.null(seed),TRUE)){
    seed=list(NA, attr(clineFrame,"mcmc.count")+1)
  }
  if(is.list(seed)&&length(seed)>1)
    attr(myClineFrame,"mcmc.count")<-seed[[2]];
  
  myParams<-splitParameters(model$parameterTypes);
##   param<-list(init=param.free.init,tune=param.free.weight,
##               lower=param.free.lower,upper=param.free.upper,
##               fixed=param.fixed);
  oldFormals<-formals(model$func);
  tttFormals<-oldFormals[names(myParams$init)];
  names(tttFormals)<-names(myParams$init);
  formals(model$func)<-c(tttFormals,myParams$fixed);
  formals(model$req)<-c(tttFormals,myParams$fixed);
## print(c(tttFormals,myParams$fixed)); 
  myRejectionLL<-rejectLL;
  clineLLFunc <- function(theta,meta.model,obsData){
##      print("D");
##      print(list(theta,meta.model,obsData));
    if(! do.call(meta.model$req,as.list(theta))) return(myRejectionLL);
    model=do.call(meta.model$func,as.list(theta));
    return(obsData$model.LL(model));
  }
  VMATRIX<-NULL;
  pMinS.cssp<-NULL;
  pMaxS.cssp<-NULL;
  if(sum(names(myParams$fixed) %in% "pMin")>0) pMinS.cssp<-myParams$fixed$pMin;
  if(sum(names(myParams$fixed) %in% "pMax")>0) pMaxS.cssp<-myParams$fixed$pMax;
  sampleModels<-cssp(sampleData$frame,
                     attr(model,"tails"),
                     pMinS=NULL,
                     pMaxS=NULL
                     )[,names(myParams$init)];
  if(dim(sampleModels)[[1]]>2){
    sampleModels<-subset(sampleModels,do.call(model$req,as.list(sampleModels[names(myParams$init)])))
    try(VMATRIX<-cov(sampleModels));
   
      
  }
  print("C");
   print(VMATRIX);
## print(  clineLLFunc (myParams$init,meta.model=model,obsData=sampleData));
##  print(list(fun=clineLLFunc, logfun="TRUE",
##                  burnin=burnin, mcmc=mcmc, thin=thin,
##                  theta.init=myParams$init,
##                  tune=myParams$tune,
##                  meta.model=model,
##                  obsData=sampleData,
##                  seed=seed
             ##     optim.control=list(fnscale = -1, trace = 0,
##                    REPORT = 10,maxit = 5000),
##                  optim.method= "L-BFGS-B",
##                  optim.lower=myParams$lower,
##                  optim.upper=myParams$upper
##              ))
 
##   print(optim(myParams$init, function(ttt) clineLLFunc(ttt,model,sampleData),
##               control = list(fnscale = -1, trace = 0, REPORT = 10,maxit = 5000),
##             ##lower = myParams$lower, upper =myParams$upper ,
##               method =  "SANN", 
##             hessian = FALSE));
## stop("readyToRun");
  myClineFrame$mcmc <-
    MCMCmetrop1R(fun=clineLLFunc, logfun="TRUE",
                 burnin=burnin, mcmc=mcmc, thin=thin,
                 theta.init=myParams$init,
                 tune=myParams$tune,
                 meta.model=model,
                 force.samp=TRUE,
                 obsData=sampleData,
                 seed=seed,
                 V=VMATRIX,
                 verbose=verbose,
                 optim.control=list(fnscale = -1, trace = 0,
                   REPORT = 10,maxit = 5000),
                 optim.method= "L-BFGS-B",
                 optim.lower=myParams$lower,
                 optim.upper=myParams$upper);

  oldParams<-model$parameterTypes;
  colnames(myClineFrame$mcmc)<-names(myParams$init)
##   print(summary(myClineFrame$mcmc));
## plot(myClineFrame$mcmc);
  model.info=as.data.frame(myClineFrame$mcmc);
  nSamples=length(model.info[[names(myParams$init)[[1]]]]);
  if(length(myParams$fixed>0)){
    for(iter in names(myParams$fixed)){
      model.info[[iter]]=rep(myParams$fixed[[iter]],nSamples);
    }
  }
  for(iter in 1:nSamples){
    model.info[iter,"model.LL"] <-
      clineLLFunc(model.info[iter,names(myParams$init)],
                  model,
                  sampleData);
  }
  
  myClineFrame$model=model;
  myClineFrame$allClines=model.info;
  myClineFrame$maxLL=max(model.info$model.LL);
  myClineFrame$maxLL.theta=subset(model.info,model.info$model.LL==max(model.info$model.LL))[1,];
   myClineFrame$maxLL.theta<-myClineFrame$maxLL.theta[names( formals(model$func))];
  myClineFrame$maxLL.cline=do.call(model$func,as.list( myClineFrame$maxLL.theta))
  class(myClineFrame)<-"clineModelFrame";
  return(myClineFrame);
}


reFitClineFunc<-function(clineFrame,mcmc=1e6, verbose=10000,thin=NULL,
                         seed=NULL,rejectLL=-1e8){
  myClineFrame=list(data=clineFrame$data);
  model<-clineFrame$model;
  sampleData<-clineFrame$data;
  attr(myClineFrame,"mcmc.count")<-0; 
  if(identical(is.null(seed),TRUE)){
    if(is.numeric(attr(clineFrame,"mcmc.count"))&&
       length(attr(clineFrame,"mcmc.count")==1)){
      seed=list(NA, attr(clineFrame,"mcmc.count")+1)
    } else {
      seed=list(NA, attr(myClineFrame,"mcmc.count")+1)
    }
  }
  if(is.list(seed)&&length(seed)>1)
    attr(myClineFrame,"mcmc.count")<-seed[[2]];
  if(is.null(thin)) thin<-thin(clineFrame$mcmc);
  for(typeName in names(model$parameterTypes[lapply(model$parameterTypes,attr,"fixed")==FALSE])) #names(clineFrame$maxLL.theta))
    model$parameterTypes[[typeName]]$val<-
      weighted.mean(clineFrame$allClines[clineFrame$allClines$model.LL>(clineFrame$maxLL-8),typeName],exp(clineFrame$allClines$model.LL[clineFrame$allClines$model.LL>(clineFrame$maxLL-8)]));

                                        
## clineFrame$maxLL.theta[[typeName]][[1]];
  
  myParams<-splitParameters(model$parameterTypes);

  oldFormals<-formals(model$func);
  tttFormals<-oldFormals[names(myParams$init)];
  names(tttFormals)<-names(myParams$init);
  formals(model$func)<-c(tttFormals,myParams$fixed);
  formals(model$req)<-c(tttFormals,myParams$fixed);

  myRejectionLL<-rejectLL;
  clineLLFunc <- function(theta,meta.model,obsData){
    if(! do.call(meta.model$req,as.list(theta))) return(myRejectionLL);
    model=do.call(meta.model$func,as.list(theta));
    return(obsData$model.LL(model));
  }
  attach(clineFrame$allClines);
  credibleLLspace<-data.frame(LL=sort(model.LL), percentile=cumsum(exp(sort(model.LL)))/sum(exp(sort(model.LL))));
  detach(clineFrame$allClines);
  credible.LLcut<-min(subset(credibleLLspace,credibleLLspace$percentile>0.05)$LL);
  VMATRIX<-cov.wt(subset(clineFrame$allClines,clineFrame$allClines$model.LL>=clineFrame$maxLL-8)[,c(names(myParams$init),"model.LL")],wt=exp(clineFrame$allClines$model.LL[clineFrame$allClines$model.LL>=(clineFrame$maxLL-8)]))$cov;
  counter.inv<-diag(1/VMATRIX["model.LL",])
  dimnames(counter.inv)<-list(rownames(VMATRIX),colnames(VMATRIX));
  VMATRIX<-(counter.inv%*%VMATRIX%*%counter.inv)
  VMATRIX<-VMATRIX[names(myParams$init),names(myParams$init)]/VMATRIX["model.LL","model.LL"];
  print("C");
   print(VMATRIX);
 ##  print(list(fun=clineLLFunc, logfun="TRUE",
##                  burnin=0, mcmc=mcmc, thin=thin,
##                  theta.init=myParams$init,
##                  tune=myParams$tune,
##                  meta.model=model,
##                  force.samp=TRUE,
##                  obsData=clineFrame$data,
##                  seed=seed,
##                  V=VMATRIX,
##                  verbose=verbose,
##                  optim.control=list(fnscale = -1, trace = 0,
##                    REPORT = 10,maxit = 5000),
##                  optim.method= "L-BFGS-B",
##                  optim.lower=myParams$lower,
##                  optim.upper=myParams$upper));
  myClineFrame$mcmc <-
    MCMCmetrop1R(fun=clineLLFunc, logfun="TRUE",
                 burnin=0, mcmc=mcmc, thin=thin,
                 theta.init=myParams$init,
                 tune=myParams$tune,
                 meta.model=model,
                 force.samp=TRUE,
                 obsData=sampleData,
                 seed=seed,
                 V=VMATRIX,
                 verbose=verbose,
                 optim.control=list(fnscale = -1, trace = 0,
                   REPORT = 10,maxit = 5000),
                 optim.method= "L-BFGS-B",
                 optim.lower=myParams$lower,
                 optim.upper=myParams$upper);

  oldParams<-model$parameterTypes;
  colnames(myClineFrame$mcmc)<-names(myParams$init)
    model.info=as.data.frame(myClineFrame$mcmc);
  nSamples=length(model.info[[names(myParams$init)[[1]]]]);
  if(length(myParams$fixed>0)){
    for(iter in names(myParams$fixed)){
      model.info[[iter]]=rep(myParams$fixed[[iter]],nSamples);
    }
  }
  for(iter in 1:nSamples){
    model.info[iter,"model.LL"] <-
      clineLLFunc(model.info[iter,names(myParams$init)],
                  model,
                  myClineFrame$data);
  }
  myClineFrame$mcmc.old=clineFrame$mcmc;
  myClineFrame$mcmc.new=myClineFrame$mcmc;
  
  if((thin( myClineFrame$mcmc)==thin(clineFrame$mcmc))&
     (dim( myClineFrame$mcmc)[[2]]==dim(clineFrame$mcmc)[[2]])){
    myClineFrame$mcmc<-mcmc(rbind(clineFrame$mcmc, myClineFrame$mcmc),
                            start(clineFrame$mcmc),
                            end(clineFrame$mcmc)+end(myClineFrame$mcmc),
                            thin(myClineFrame$mcmc));
  }else {
    myClineFrame$mcmc<-mcmc(myClineFrame$mcmc,
                            start(myClineFrame$mcmc)+ end(clineFrame$mcmc),
                            end(clineFrame$mcmc)+end(myClineFrame$mcmc),
                            thin(myClineFrame$mcmc));
  }
  myClineFrame$model=model;
  myClineFrame$allClines.old=clineFrame$allClines;
  myClineFrame$allClines.new=model.info;
  myClineFrame$allClines=rbind(clineFrame$allClines,model.info);
  myClineFrame$maxLL=max( myClineFrame$allClines$model.LL);
  myClineFrame$maxLL.theta=subset( myClineFrame$allClines, myClineFrame$allClines$model.LL==max( myClineFrame$allClines$model.LL))[1,];
   myClineFrame$maxLL.theta<-myClineFrame$maxLL.theta[names( formals(model$func))];
  myClineFrame$maxLL.cline=do.call(model$func,as.list( myClineFrame$maxLL.theta))
  class(myClineFrame)<-"clineModelFrame";
  return(myClineFrame);
}

plotClineFunc <- function(clineFrame){

  attach(clineFrame$data$frame);
  x<-seq(min(dist),max(dist),length.out=101);
  detach();
  plot(obsFreq ~ dist, data=clineFrame$data$frame);
  lines(clineFrame$maxLL.cline(x)~x);

}





## ## format of the cline Multi Model Frame
## objCMMFrame<-list();
## objCMMFrame$modelFrames <- list(objCFrame1,objCFrame2);
## objCMMFrame$data <-objSampleData;  ##Keep Seperate?
## objCMMFrame$parameters <-objParameterDataFrame;  ##Attach to frames?
## class(objCMMFrame)<-"clineMultiModelFrame";
