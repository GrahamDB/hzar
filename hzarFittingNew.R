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
  try(result<-
    MCMCmetrop1R(fun=fitRequest$llFunc, logfun="true", force.samp=TRUE,
                 mcmc=useParam$chainLength, burnin=useParam$burnin,
                 verbose=useParam$verbosity, thin=useParam$thin,
                 theta.init=mdlParam$init, tune=mdlParam$tune,
                 V=fitRequest$cM, seed=useParam$seed,
                 optim.control=list(fnscale=-1,trace=0,REPORT=10,maxit=5000),
                 optim.method="L-BFGS-B",
                 optim.lower=useParam$lower,
                 optim.upper=useParam$upper )
      )
                 
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

require(foreach);

hzar.eval.clineLL <- function(data, llFunc){
  result<-foreach(ttt=iter(data,by='row'),.combine=c) %dopar% { llFunc(ttt); };
  return(result);
}
           
hzar.gen.rParam.uniform<-function(param.lower,param.upper,count=1000){
  raw<-foreach(low=param.lower,
               high=param.upper,
               .combine=cbind) %dopar% {runif(count,low,high)};
  result<-as.data.frame(raw);
  colnames(result)<-names(param.lower);
  names(result)<-names(param.lower);
  return(result);
}

## So, what do I need to know to generate the matrix?
## Work backwards?
