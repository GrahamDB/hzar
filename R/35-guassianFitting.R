

g.LLfuncA <- function(obsData, reqExp, muExp, varExp, tFunc, tArgs, tFixed,
           LLrejectedModel = -1e+08){
    baseFunc <- function(theta) 0;
    
    gLL <-
      guassianThetaLLExpF(distance=obsData$frame$dist,
                          sampleMean=obsData$frame$mu,
                          sampleVar=obsData$frame$var,
                          nEff=obsData$frame$nEff,
                          muExp=muExp,
                          varExp=varExp)
    gLL <-  step1VectorExpF(reqExp,gLL,LLrejectedModel)
    gLL <- eval(substitute(substitute(LLfunc,tF),list(LLfunc=gLL,tF=tFixed)))
    lapply(1:length(tArgs),function(x) bquote(theta[.(x)]))->tMap
    names(tMap) <- tArgs;
    body(baseFunc) <- eval(substitute(substitute(LLfunc,eL),
                                      list(LLfunc=gLL,eL=tMap)))
    environment(baseFunc) <- .GlobalEnv
    baseFunc
}
g.LLfunc <- function(obsData, model,
           tFixed=model$init[as.logical(model$fixed)],
           tArgs=model$init[!as.logical(model$fixed)],
           LLrejectedModel = -1e+08){
    baseFunc <- function(theta) 0;
    model.req=model$req;
    model.gen=model$mFunc
    model.obsData=obsData
    param.fixed=tFixed
    
    new.formals=tArgs
    
    gLL <-
      guassianThetaLLExpF(distance=obsData$frame$dist,
                          sampleMean=obsData$frame$mu,
                          sampleVar=obsData$frame$var,
                          nEff=obsData$frame$nEff,
                          muExp=model$mu,
                          varExp=model$var)
    gLL <-  step1VectorExpF(body(model$req)[[2]],gLL,LLrejectedModel)
    gLL <- eval(substitute(substitute(LLfunc,tF),list(LLfunc=gLL,tF=tFixed)))
    
    lapply(1:length(tArgs),function(x) bquote(theta[[.(x)]]))->tMap
    names(tMap) <- names(tArgs);
    gLL <- eval(substitute(substitute(LLfunc,tF),
                           list(LLfunc=gLL,tF=c(tFixed,tMap))))
    body(baseFunc) <-as.call(c(as.name("{"),
                               bquote(res <- .(gLL)),
                               expression(if(any(is.na(res))) print(theta)),
                               bquote(ifelse(is.na(res),.(LLrejectedModel),res))));
    llFunc=baseFunc
    ## eval(substitute(substitute(LLfunc,eL),
    ##                                 list(LLfunc=gLL,eL=tMap)))
    ## body(baseFunc) <- substitute(evalq(LLfunc,envir=theta),list(LLfunc=gLL))
    ## environment(baseFunc) <- .GlobalEnv 
    baseFunc
}


naiveHessian <- function(theta,zF,k=0.05){
  zFunc <- function(a,b,c,d)
    zF(a)+zF(d)-zF(b)-zF(c);
  n=length(theta)
  res <- matrix(0,nrow=n,ncol=n);
  for(i in 1:n){
    for(j in i:n){
      
      if(i==j){
        t0=t2=theta
        t0[[i]]=theta[[i]]*(1-k)
        t2[[i]]=theta[[i]]*(1+k)
        res[i,i]=zFunc(t0,theta,theta,t2)/(k^2*theta[[i]]^2);
    
      }else{
        t00=t02=t20=t22=theta
        t00[[i]]=theta[[i]]*(1-k)
        t02[[i]]=theta[[i]]*(1-k)
        t20[[i]]=theta[[i]]*(1+k)
        t22[[i]]=theta[[i]]*(1+k)
        t00[[i]]=theta[[j]]*(1-k)
        t20[[i]]=theta[[j]]*(1-k)
        t02[[i]]=theta[[j]]*(1+k)
        t22[[i]]=theta[[j]]*(1+k)
        res[j,i]=res[i,j]=
          zFunc(t00,t20,t02,t22)/(k^2*theta[[i]]*theta[[j]]);
      }
    }
  }
  rownames(res) <- names(theta)
  colnames(res) <- names(theta)
  
  
  return(res);

}

## gC.V <- function(model,LLfunc){
##   solve(-naiveHessian(model$init,LLfunc));
## }

hzar.first.fitRequest.gC <- function(gModel,obsData,verbose=TRUE){
  if (verbose) {
    mcmcParam <- hzar:::cfg.hzar.default.mcmc
  } else {
    mcmcParam <- hzar:::cfg.hzar.quiet.mcmc
  }
  modelParam <- list();
  modelParam$init <- gModel$init[!as.logical(gModel$fixed)];
  modelParam$tune <- gModel$tune[!as.logical(gModel$fixed)];
  
  
  LLfunc <- g.LLfunc(obsData,gModel)
  cV <-  solve(-naiveHessian(modelParam$init,LLfunc));
  return(hzar.make.fitRequest(modelParam, cV, LLfunc, 
        mcmcParam))
}
