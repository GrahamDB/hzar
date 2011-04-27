

require(foreach);
source("hzarClasses.R")

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
 # clineRequests<-data.frame(c.s=scaling,c.t=tails);
  foreach(c.rqst.c.s=scaling,c.rqst.c.t=tails) %do% {
    makeCline1D(data=data,scaling=c.rqst.c.s,tails=c.rqst.c.t,direction=direction);
  } -> clineMultiModel;
  class(clineMultiModel)<-"clineMultiModel";
  return(clineMultiModel);
}


fitClineFunc.staged<-function(model,data,seedChannel=1,verbose=0,mcmc=1e6,rejectLL=-1e8,preBurnin=1e3,mcmcB1=5e5,mcmcB2=1e6,thin=100){
  clineFrame.B<-NULL;
  clineFrame.C<-NULL;
  for(iter in names(model$parameterTypes)){
    model$parameterTypes[[iter]]$w<-1.5;
  }
  
  clineFrame.A<-fitClineModel(model,
                              sampleData=data,
                              verbose=verbose,
                              burnin=preBurnin,
                              mcmc=mcmcB1,
                              thin=thin,
                              seed=list(NA,seedChannel),
                              rejectLL=rejectLL);
  ##try(  print(clineFrame.A$allClines.new$model.LL[1:100]));
##  lapply(clineFrame.A$model$parameterTypes,function(x) {x$w<<-1.2;})
  for(iter in names(model$parameterTypes)){
    clineFrame.A$model$parameterTypes[[iter]]$w<-0.5;
  }
  try( clineFrame.B<-reFitClineFunc( clineFrame.A,
                                    verbose=verbose,
                                    mcmc=mcmcB2,
                                    thin=thin,
                                    rejectLL=rejectLL));
  if(identical(is.null(clineFrame.B),TRUE))
    return(clineFrame.A);
  ## try(  plot(clineFrame.B$mcmc,density=FALSE,ask=FALSE));
##try(  print(clineFrame.B$allClines.new$model.LL[1:100]));
  ##  lapply(clineFrame.B$model$parameterTypes,function(x) {x$w<<-1.1;})
  for(iter in names(model$parameterTypes)){
    clineFrame.B$model$parameterTypes[[iter]]$w<-1.5;
  }
  try( clineFrame.C<-reFitClineFunc( clineFrame.B,
                                    verbose=verbose,
                                    mcmc=mcmc,
                                    thin=thin,
                                    rejectLL=rejectLL));
  if(identical(is.null(clineFrame.C),TRUE))
    return(clineFrame.B);
  ## try(  plot(clineFrame.C$mcmc,density=FALSE,ask=FALSE));
##try(  print(clineFrame.C$allClines.new$model.LL[1:100]));

  if(getCredibleCut(clineFrame.C,0.01)>clineFrame.C$maxLL-4)
    return(clineFrame.C)
for(iter in names(model$parameterTypes)){
    clineFrame.C$model$parameterTypes[[iter]]$w<-1.1;
  }
  ##lapply(clineFrame.C$model$parameterTypes,function(x) {x$w<<-0.5;})
   try( clineFrame.C<-reFitClineFunc( clineFrame.C,
                              verbose=verbose,
                              mcmc=mcmc*2,
                              thin=thin,
                              rejectLL=rejectLL));
  
  if(getCredibleCut(clineFrame.C,0.01)>clineFrame.C$maxLL-4)
    return(clineFrame.C)
for(iter in names(model$parameterTypes)){
    clineFrame.C$model$parameterTypes[[iter]]$w<-0.9;
  }
  ##lapply(clineFrame.C$model$parameterTypes,function(x) {x$w<<-0.5;})
   try( clineFrame.C<-reFitClineFunc( clineFrame.C,
                              verbose=verbose,
                              mcmc=mcmc*2,
                              thin=thin,
                              rejectLL=rejectLL));
   return(clineFrame.C)
}


## fitting multiple clines
## Need to assess package of all the clines
fitMultipleModels<-function(multiModel,data,mcmc=5e6,seedChannelSpacing=1,rejectLL=-1e8,verbose=0){
  myClineFrame<-list();
  class(myClineFrame)<-"clineMultiModelFrame"
  myClineFrame$modelList<-multiModel
  sChannels=seq(from=1,by=seedChannelSpacing,length.out=length(multiModel));
  foreach(clineModel=multiModel,seedChannel=sChannels) %dopar% {
    fitClineFunc.staged(model=clineModel,
                        data=data,
                        seedChannel=seedChannel,
                        mcmc=mcmc,
                        rejectLL=rejectLL,
                        verbose=verbose);
  } -> myClineFrame$clineFrames
   ## myClineFrame$mcmc.comp
  return(myClineFrame);
}
## External Functions to summarize / compare internals?

