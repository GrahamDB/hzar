

require(foreach);
source("hzarClasses.R")

getCredibleCut<-function(clineFrame,rejectPercent=0.05){

model.LL=clineFrame$allClines$model.LL
  credibleLLspace<-data.frame(LL=sort(model.LL), percentile=cumsum(exp(sort(model.LL)))/sum(exp(sort(model.LL))));

  credible.LLcut<-min(subset(credibleLLspace,credibleLLspace$percentile>rejectionPercent)$LL);
return(credible.LLcut)
}
getCredibleLLspace<-function(clineFrame){

model.LL=clineFrame$allClines$model.LL
  credibleLLspace<-data.frame(LL=sort(model.LL), percentile=cumsum(exp(sort(model.LL)))/sum(exp(sort(model.LL))));

return(credibleLLspace)
}

makeMultiCline1D<-function(data=NULL,scaling=c("fixed","free","free"),tails=c("none","none","both"),direction=NULL){
  if((!identical(is.null(direction),FALSE))&&
     (!identical(is.null(data),FALSE))) {
    stop("Either direction or data must be specified.");
  }
  if((!identical(is.null(scaling),FALSE))||
     (!identical(is.null(tails),FALSE))) {
    stop("The parameters 'scaling' _and_ 'tails' must be specified.");
  }
  if(length(scaling) != length(tails)) {
    if(length(scaling)==1) scaling=rep(scaling,length(tails));
    if(length(tails)==1)   tails=  rep(tails,  length(scaling));
    if(length(scaling) != length(tails))
      stop(paste("  Scaling:",scaling,"\nand tails:",tails,"\nnot of equal length.",sep="\t"));
  }
  clineRequests<-data.frame(c.s=I(scaling),c.t=I(tails));
  foreach(c.rqst=iter(clineRequests)) %do% {
    makeCline1D(data=data,scaling=c.rqst$c.s,tails=c.rqst$c.t,direction=direction);
  } -> clineMultiModel;
  class(clineMuliModel)<-"clineMultiModel";
  return(clineMultiModel);
}


## fitting multiple clines
## Need to assess package of all the clines

 class(myClineFrame)<-"clineMultiModelFrame"

## External Functions to summarize / compare internals?

