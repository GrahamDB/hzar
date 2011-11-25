## This is the new way of fitting the cline model, in an attempt to
## streamline the process. I want to break the code into smaller
## chunks by taking advantage of internal methods.  That way I can
## provide complex behavior to the end user with minimal code.
## Although technically the code will be more complex in structure,
## each method should have a clear meaning.

## In writing this, I plan on starting with biggest atomic slice of
## the fitting process, the call to MCMCmetrop1R.

hzar.doFit <- function(fitRequest){
  result<-NULL;                        
  useParam=fitRequest$mcmcParam;
  mdlParam=fitRequest$modelParam;
  try({
    result<-
      MCMCmetrop1R(fun=fitRequest$llFunc, logfun="true", force.samp=TRUE,
                   mcmc=useParam$chainLength, burnin=useParam$burnin,
                   verbose=useParam$verbosity, thin=useParam$thin,
                   theta.init=mdlParam$init, tune=mdlParam$tune,
                   V=fitRequest$cM, seed=useParam$seed,
                   optim.control=list(fnscale=-1,trace=0,REPORT=10,maxit=5000),
                   optim.method="L-BFGS-B",
                   optim.lower=useParam$lower,
                   optim.upper=useParam$upper );
    
    colnames(result)<-names(mdlParam$init);
  });
  fitRequest$mcmcRaw<-result;
  attr(fitRequest,"fit.run")<-TRUE;
  attr(fitRequest,"fit.success")<-FALSE;
  if(identical(is.null(result),TRUE)){
    warning("Fitting failed.");
  } else {
    attr(fitRequest,"fit.success")<-TRUE;
  }
  return(fitRequest);
}

## fitRequest has the following base structure:
hzar.make.fitRequest <-
  function(modelParameters,             #output of splitParameters
           covMatrix,                   #covariance matrix to use
           clineLLfunc,                 #method to get LL given theta
           mcmcParameters,              #mcmc attributes & controls
           mcmcRaw=NULL,                #what MCMCmetrop1R returns
           fit.run=FALSE,               #has this request been run?
           fit.success=FALSE            #has the run succeeded?
           ){
    fitRequest<-list(modelParam=modelParameters,cM=covMatrix,
                     llFunc=clineLLfunc,mcmcParam=mcmcParameters,
                     mcmcRaw=mcmcRaw);
    class(fitRequest)<-"hzar.fitRequest";
    attr(fitRequest,"fit.run")<-fit.run;
    attr(fitRequest,"fit.success")<-fit.success;
    return(fitRequest);
  }

## I already have splitParameters to make modelParameters.

## mcmcParameters should have a simple structure:
##
## chainLength, burnin, verbosity, thin, seed(seedStreamChannel,
## useSeedStream, mersenneSeed, lecuyerSeed)
hzar.make.mcmcParam <-
  function(chainLength, burnin, verbosity, thin,
           seedStreamChannel=1, useSeedStream=TRUE,
           mersenneSeed=12345,
           lecuyerSeed=rep(12345,6)){
    mcmcParam<-list(chainLength=chainLength,burnin=burnin,
                    verbosity=verbosity,thin=thin);
    mcmcParam$seed<-mersenneSeed;
    if(useSeedStream){
      mcmcParam$seed<-list(lecuyerSeed,seedStreamChannel)
    }
    return(mcmcParam);
  }



## I need a method to make clineLLfunc.  I will start with the
## original cline log likelihood function, wrapped with a basic method
## call.

## I need two llFuncs, one for Maximum Likelihood and one for Bayesian
## analysis. Also, I need to recognize that under the old process, the
## fixed parameters are magically inserted during engine startup to
## simplify the method calls.
hzar.make.clineLLfunc.old.ML <-
  function(param.free.names, param.fixed, #lists made by splitParam
           param.check.func,              #func to check parameters
           meta.cline.func,               #func to generate cline func
           model.LL,                      #returns LL of cline func
           LLrejectedModel=-1E8           #when rejecting, return this
           ){
    ## First, assign references locally.
    model.req<-param.check.func;
    model.gen<-meta.cline.func;
    eval.clineLL<- model.LL;
    myRejectionLL<-LLrejectedModel;
    ## Second, modify function signitures to account for fixed
    ## parameters.
    old.formals<-formals(model.gen);
    if(length(old.formals)!=(length(param.free.names)+
                             length(param.fixed))){
      warning("The length of the method formals does not match the length of the parameters supplied.");
    }
    ttt.formals<-old.formals[param.free.names];
    names(ttt.formals)<-param.free.names;
    new.formals<-c(ttt.formals,param.fixed);
    formals(model.req)<-new.formals;
    formals(model.gen)<-new.formals;
    
    ## references: theta, meta.model, obsData, myRejectionLL
    ##
    ## in-depth:  theta, meta.model$req, myRejectionLL, meta.model$func, 
    ## obsData$model.LL
    
    llFunc<-function(theta){
      if(! do.call(model.req,as.list(theta))) return(myRejectionLL);
      model=do.call(model.gen,as.list(theta));
      result<-eval.clineLL(model);
      if(identical(is.finite(result),TRUE))
        return(result);
      return(myRejectionLL);
    }
    
    return(llFunc);
  }

hzar.make.clineLLfunc.old.bayes <-
  function(param.free.names, param.fixed, #lists made by splitParam
           param.check.func,              #func to check parameters
           meta.cline.func,               #func to generate cline func
           model.LL,                      #returns LL of cline func
           prior.LL,                      #returns LL of priors
           LLrejectedModel=-1E8           #when rejecting, return this
           ){

    ## First, assign references locally.
    model.req<-param.check.func;
    model.gen<-meta.cline.func;
    eval.priorLL<-prior.LL;
    eval.clineLL<- model.LL;
    myRejectionLL<-LLrejectedModel;
    ## Second, modify function signitures to account for fixed
    ## parameters.
    old.formals<-formals(model.gen);
    if(length(old.formals)!=(length(param.free.names)+
                             length(param.fixed))){
      warning("The length of the method formals does not match the length of the parameters supplied.");
    }
    
    ttt.formals<-old.formals[param.free.names];
    names(ttt.formals)<-param.free.names;
    new.formals<-c(ttt.formals,param.fixed);
    formals(model.req)<-new.formals;
    formals(model.gen)<-new.formals;
    formals(eval.priorLL)<-new.formals;
    
    ## references: theta, meta.model, obsData, myRejectionLL
    ##
    ## in-depth: theta, meta.model$req, myRejectionLL, meta.model$func,
    ## meta.model$prior, obsData$model.LL
    llFunc<-function(theta){
      if(! do.call(model.req,as.list(theta))) return(myRejectionLL);
      model=do.call(model.gen,as.list(theta));
      thetaLL=do.call(eval.priorLL,as.list(theta));
      
      result<-eval.clineLL(model)+thetaLL;
      if(identical(is.finite(result),TRUE))
        return(result);
      return(myRejectionLL);
    }
    
    return(llFunc);
  }


## I need method(s) to generate covMatrix.

## I will need helper functions, specifically one to generate a vector
## of likelihood values given a data frame of parameter values. I want
## this to be capable of executing in parallel.

## require(foreach);

fitter.wedgeSlice <- function (count,slice.size=1000){
  res=list();
  
  if(count>slice.size)
    res<-lapply(1:((count-1)%/%slice.size)-1,
                function(x,y=slice.size) { 1:y + x*y });
  return(c(res,list(1:(1+(count-1)%%slice.size)
                    +slice.size*((count-1)%/%slice.size))));
}

hzar.eval.clineLL <- function(data, llFunc,doPar=FALSE){
  ## print("A");
  slices<-fitter.wedgeSlice(dim(data)[[1]]);
  ## cat("Eval wedge size:",object.size(slices),"\n");
  
  useFunc=llFunc;
  useData=data;
  if(doPar){
    
    result<-foreach(tttIndex=slices,
                    .combine=rbind) %dopar% {
                      data.frame(model.LL=as.numeric(lapply(tttIndex,function(x)
                                   {useFunc(useData[x,]); })));}
  }else{
    
    result<-foreach(tttIndex=slices,
                    .combine=rbind) %do% {
                      data.frame(model.LL=as.numeric(lapply(tttIndex,function(x)
                                   {useFunc(useData[x,]); })));}
  }
  ## cat("Eval result size:",object.size(result),"\n");
  return(result);
}
           
fitter.gen.rParam.uniform<-function(param.lower,param.upper,count=1000){
  raw<-foreach(low=param.lower,
               high=param.upper,
               .combine=cbind) %do% {runif(count,low,high)};
  result<-as.data.frame(raw);
  colnames(result)<-names(param.lower);
  names(result)<-names(param.lower);
  return(result);
}

## So, what do I need to know to generate the matrix?
## Work backwards?

## cov.wt requires a data frame of sampled parameters, weights
## (Likelihood, intregration differentials, scaling), and possibly a
## pre-calculated center.

## data    := sampled parameters

## data.wt := weights needs: (data, clineLLfunc, data.dP

## data.LL := Likelihood; needs (clineLLfunc)
## data.LL =  hzar.eval.clineLL(data, clineLLfunc)
## data.dP := area given to each sampled point
## scaling := 1/(sum(exp(data.LL)*data.dP))

fitter.getCovWeights <- function(data,clineLLfunc, data.dP){
  data.LL <- hzar.eval.clineLL(data,clineLLfunc);
  return(exp(data.LL)*data.dP/sum(exp(data.LL)*data.dP));
}

## I can assume I have clineLLfunc, so I just need data and data.dP. I
## could just sample a rectangular area.

fitter.gen.samples.rect <- function(param.lower, param.upper, pDiv=11){
  param.names<-names(param.lower);
  nParam<-length(param.names);
  deltas=(as.numeric(param.upper[param.names])-
          as.numeric(param.lower[param.names]))/(pDiv-1);
 ##  print(deltas);print(param.names);
  grid.formals<-matrix(nrow=pDiv,ncol=nParam,
                       data=rep(deltas,each=pDiv))* rep(0:(pDiv-1),nParam)+
                         rep(as.numeric(param.lower[param.names]),each=pDiv);
  colnames(grid.formals)<-param.names;
  ## names(grid.formals)<-param.names;
  ## print(as.list(as.data.frame(grid.formals)));
  result<-list(dTheta=prod(abs(deltas)));
  result$data<-do.call(expand.grid,as.list(as.data.frame(grid.formals)));
  ## print(class(result$data))
  ## grid.formals<-names(param.lower
  return(result);
}

hzar.cov.rect<-function(clineLLfunc,param.lower,param.upper,pDiv=11,random=0,passCenter=FALSE){
  ## print("A");

  ## Check to make sure we aren't about to do something stupid.
  ## Yes, that means that models with over ten million parameters
  ## will fail spectacularly. We hope that you would have
  ## considered writing your own software for problems of such
  ## scale.
  if(random>1e9){
    stop("Covariance matrix calculation with random sampling requested with far too many samples.  Stopping.")
  }
  if(random>0){
    data.mat<-list(dTheta=prod(abs(as.numeric(param.upper)
                     -as.numeric(param.lower)))
                   / random,
                   data=fitter.gen.rParam.uniform(param.lower,
                     param.upper,
                     random));
  }else{
    ## Stupidity check.  See above
    if((pDiv^length(param.lower))>1e6){
      ## This is recoverable by switching to random sampling.
      warning("Covariance matrix calculation requested for too complex of a lattice structure. Switching to random sampling using ten thousand samples.");
      return(hzar.cov.rect(clineLLfunc,
                           param.lower,
                           param.upper,
                           random=1e4,
                           passCenter=passCenter));
    }
    data.mat<-fitter.gen.samples.rect(param.lower,param.upper,pDiv);
  }
  param.names<-names(data.mat$data);##print(names(data.mat));
  ## print("A");
  ##data.wt<-fitter.getCovWeights(data.mat$data,clineLLfunc,data.mat$dTheta);
  data.wt<-hzar.eval.clineLL(data.mat$data,clineLLfunc);
  data.mat$data<-data.mat$data[data.wt>-1e6,];
  data.wt<-data.wt[data.wt>-1e6];
  MIN.DATA<-(1+length(param.upper))
  if(length(data.wt)<MIN.DATA){
    ## need more samples... recurse.
    if(random>0){ 
      ## Double the amount of sampling
      return(hzar.cov.rect(clineLLfunc,
                           param.lower,
                           param.upper,
                           random=2*random,
                           passCenter=passCenter));
    }else{
      ##Increase the resolution slightly.
      return(hzar.cov.rect(clineLLfunc,
                           param.lower,
                           param.upper,
                           pDiv=pDiv+1,
                           passCenter=passCenter));
    }
  }
  while(sum(data.wt>-723)<MIN.DATA){
    ## insufficient data in finite range
    if(sum(data.wt>609)>0){
      ## scaling won't fix the problem, so more samples needed.
      ## Recurse.
      if(random>0){ 
        ## Double the amount of sampling
        return(hzar.cov.rect(clineLLfunc,
                             param.lower,
                             param.upper,
                             random=2*random,
                             passCenter=passCenter));
      }else{
        ##Increase the resolution slightly.
        return(hzar.cov.rect(clineLLfunc,
                             param.lower,
                             param.upper,
                             pDiv=pDiv+1,
                             passCenter=passCenter));
      }
    }
    ## iteratively shift likelihood space to bring samples into finite range.
    data.wt<-data.wt+100;
  }
  ## print("A");
  VDATA<-cov.wt(x=cbind(data.mat$data,model.LL=data.wt),wt=exp(data.wt))
  ##   *data.mat$dTheta)
  VMATRIX<-VDATA$cov;
  ## diag(1/(VMATRIX["model.LL",]))->counter.inv;
  ## diag(sign(VMATRIX["model.LL",]))->counter.inv;
  ## dimnames(counter.inv)<-list(rownames(VMATRIX),colnames(VMATRIX));
  ## counter.inv2<-counter.inv/sqrt(counter.inv["model.LL","model.LL"]);
  ## mat.scaled<-counter.inv%*%VMATRIX%*%counter.inv;
  mat.scaled<-VMATRIX;
  if(passCenter)
    return(list(cov=mat.scaled[param.names,param.names],center=VDATA$center[param.names]));
  return(mat.scaled[param.names,param.names]);
}

cfg.hzar.default.mcmc <- hzar.make.mcmcParam(chainLength=1e6,
                                         burnin=1e4,
                                         verbosity=5e5,
                                         thin=100);

cfg.hzar.quiet.mcmc   <- hzar.make.mcmcParam(chainLength=1e6,
                                         burnin=1e4,
                                         verbosity=0,
                                         thin=100);

##refitting
hzar.cov.mcmc<-function(clineLLfunc,mcmcRaw,pDiv=15,random=1e4,passCenter=FALSE){
  mcmc.nm<-colnames(mcmcRaw);
  pL<-lapply(mcmc.nm,function(x){min(mcmcRaw[,x])});
  names(pL)<-mcmc.nm;
  pU<-lapply(mcmc.nm,function(x){max(mcmcRaw[,x])});
  names(pU)<-mcmc.nm;
  return(hzar.cov.rect(clineLLfunc,pL,pU,pDiv,random,passCenter));
}
hzar.next.fitRequest <- function(oldFitRequest){
  seedChannel<-1;
  if(is.list(oldFitRequest$mcmcParam$seed)){
    seedChannel=oldFitRequest$mcmcParam$seed[[2]];
  }
  if(identical( attr(oldFitRequest,"fit.run") , TRUE)){
    seedChannel<-seedChannel+1;
  } else {
    seedChannel<-seedChannel+10;
  }
  mcmcParam=hzar.make.mcmcParam(oldFitRequest$mcmcParam$chainLength,
    oldFitRequest$mcmcParam$burnin,
    oldFitRequest$mcmcParam$verbosity,
    oldFitRequest$mcmcParam$thin,
    seedChannel);
  mdlParam<-oldFitRequest$modelParam;
  covMatrix<-oldFitRequest$cM
  if(identical( attr(oldFitRequest,"fit.success") , TRUE)){
    mcmcSubset<-oldFitRequest$mcmcRaw[sample(dim(oldFitRequest$mcmcRaw)[[1]]),];
    subLL<-hzar.eval.clineLL(mcmcSubset,oldFitRequest$llFunc);
    covData<-hzar.cov.mcmc(oldFitRequest$llFunc,mcmcSubset[subLL>max(subLL-4),],passCenter=TRUE);
    covMatrix<-covData$cov;
    new.center<-covData$center[names(mdlParam$init)];
    if(oldFitRequest$llFunc(new.center)>1e-6)
      mdlParam$init <- new.center;
  }
  return(hzar.make.fitRequest(mdlParam,
                              covMatrix,
                              oldFitRequest$llFunc,
                              mcmcParam));
}
  
hzar.first.fitRequest.old.ML <-function(model,obsData,verbose=TRUE){
  
  if(verbose){
    mcmcParam<-cfg.hzar.default.mcmc;
  }else {
    mcmcParam<-cfg.hzar.quiet.mcmc;
  }
   ## print("A");
  
  modelParam<-splitParameters(model$parameterTypes);
   ## print("A");
  clineLLfunc<-hzar.make.clineLLfunc.old.ML(names(modelParam$init),
                                            modelParam$fixed,
                                            model$req,
                                            model$func,
                                            obsData$model.LL);
   ## print("A");
  covMatrix<-NULL;
  try(  covMatrix<-hzar.cov.rect(clineLLfunc,modelParam$lower,modelParam$upper,random=1e4));
   ## print("A");
  return(hzar.make.fitRequest(modelParam,covMatrix,clineLLfunc,mcmcParam));
}


## hzar.multiFit.doNext

hzar.chain.doSeq <- function(hzar.request, count=3,collapse=FALSE,announce.complete="Chain Complete"){
  if(collapse){
    mcmcParam=hzar.request$mcmcParam;
    mcmcParam$chainLength= mcmcParam$chainLength*count;
    mcmcParam$burnin= mcmcParam$burnin*count;
    mdlParam=hzar.request$modelParam;
    
  }
  hzar.results<-list();
  for(iter in 1:count){
    print(iter);
    if(iter>1){
      hzar.request<-NULL;
      print(summary(try(hzar.request<-
                        hzar.next.fitRequest(hzar.results[[iter-1]] ))));
    }
    if(!inherits(hzar.request,c("hzar.fitRequest"))){
      warning("Failed to generate next fit request. Returning successful runs.");
      
      return(hzar.results);
    }
    hzar.results[[iter]] <- hzar.request;
    ##format(hzar.request$cM,width=8);
    print(summary(try(hzar.results[[iter]] <- hzar.doFit(hzar.request))));
    ##print(attributes(hzar.results[[iter]]))
  }
  print(announce.complete);
  if(collapse){
    raw.data<-do.call(rbind,lapply(hzar.results,function(x)x$mcmcRaw));
    
    rawMCMC<-mcmc(data=raw.data,thin=thin(hzar.results[[count]]$mcmcRaw),start=mcmcParam$burnin+1);
    names(rawMCMC)<-names(hzar.results[[count]]$mcmcRaw);
    return(list(hzar.make.fitRequest(mdlParam,
                                     hzar.results[[count]]$cM,
                                     hzar.results[[count]]$llFunc,
                                     mcmcParam,
                                     mcmcRaw=rawMCMC,
                                     TRUE,
                                     prod(as.logical(lapply(hzar.results,attr,"fit.success"))))));
  }
  return(hzar.results);
}
                                                                                                              
