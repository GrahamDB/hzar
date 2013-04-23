
thEvalFExp <- function( target, exp, val)
  bquote(.(target) <- .(eval(substitute(substitute(a,b),
                                        list(a=exp,b=val)))))
optMuPExpF <- function( sampleMean, nEff,distance,
                       fExp, varExp,
                       tMuL=quote(muL),
                       tMuR=quote(muR)){
  fTa <- quote(fXA)
  fTd <- quote(fXD)
  varT <- quote(vX)
  cLT <- quote(cL)
  cRT <- quote(cR)
  aLT <- quote(aL)
  aRT <- quote(aR)
  aMT <- quote(aM)
  tDen <- quote(aLaRaMM)
  
  
  c(thEvalFExp(fTa,fExp,list(x=distance)),
    bquote(.(fTd) <- 1 - .(fTa)),
    thEvalFExp(varT,varExp,list(x=distance)),
    bquote(.(cLT) <- sum( .(sampleMean)*.(nEff)*.(fTd)/ .(varT))),
    bquote(.(cRT) <- sum( .(sampleMean)*.(nEff)*.(fTa)/ .(varT))),
    bquote(.(aLT) <- sum( .(nEff)*.(fTd)^2/ .(varT))),
    bquote(.(aRT) <- sum( .(nEff)*.(fTa)^2/ .(varT))),
    bquote(.(aMT) <- sum( .(nEff)*.(fTd)*.(fTa)/ .(varT))),
    bquote(.(tDen) <- .(aRT)*.(aLT)-.(aMT)^2),
    bquote(.(tMuL) <- (.(cLT)*.(aRT)-.(cRT)*.(aMT))/.(tDen)),
    bquote(.(tMuR) <- (.(cRT)*.(aLT)-.(cLT)*.(aMT))/.(tDen))
    )
}

mu2fExp <- function(muExp, tMuL=quote(muL), tMuR=quote(muR))
  simplify.exp(eval(substitute(substitute(a,b),list(a=muExp,b=list(muL=0,muR=1)))))

optMuPExpF.fit <- function(fitRequest,
                           model=cline.extract.modelPrep(fitRequest)$model,...){
  if(!all(c("mu","var") %in% names(model)))
    return(list());
  fixed=cline.extract.modelPrep(fitRequest)$modelParam$fixed
  
  optMuPExpF(distance=quote(frame$dist),
             sampleMean=quote(frame$mu),
             nEff=quote(frame$nEff),
             fExp=mu2fExp(eval(substitute(substitute(a,fixed),
               list(a=model$mu)))),
             varExp=simplify.exp(eval(substitute(substitute(a,fixed),
               list(a=model$var)))),
             ...)
}

optimizeTraceF <- function(fitRequest){
  resF <- function(mcmcRaw) {as.data.frame(mcmcRaw)}
  params <- names(cline.extract.modelPrep(fitRequest)$modelParam$init)
  params <- params[!(params %in% c("muL","muR"))]
  repTheta <- lapply(params,function(x) bquote(data.param[[.(x)]][iter]))
  names(repTheta) <- params
  ## str(repTheta)
  optMuPExpF.fit(fitRequest,
                 tMuL=quote(data.param$muL[iter]),
                 tMuR=quote(data.param$muR[iter]))->optExp
  optExp <- lapply(optExp,
                   function(x) eval(substitute(substitute(a,b),
                            list(a=x,b=repTheta))))
  as.call(c(as.name("{"),
            c(quote(data.param<-as.data.frame(mcmcRaw)),
              bquote(for(iter in 1:nrow(data.param))
                     .(as.call(c(as.name("{"),optExp)))),
              quote(data.param)) ) ) -> body(resF)
  environment(resF) <- list2env(list(frame=hzar.extract.obsData(fitRequest)$frame),
                                parent=.GlobalEnv)
  resF
}

optimizeTrace <- function(fitRequest,...){
  ##print(
    traceF <- optimizeTraceF(fitRequest)
  ##  )
  
  if(inherits(fitRequest,"hzar.fitRequest")){
    if(!all(c("mcmcRaw") %in% names(fitRequest)))
      stop("No trace to optimize")
    
    tempPar <- attr(fitRequest$mcmcRaw, "mcpar")
    fitRequest$mcmcRaw <- as.mcmc(traceF(fitRequest$mcmcRaw))
    tempPar -> attr(fitRequest$mcmcRaw, "mcpar")
    return(fitRequest)
  } else if (inherits(fitRequest,"hzar.dataGroup")){
    tempPar <- attr(fitRequest$data.mcmc, "mcpar")
    
    fitRequest$data.param <-traceF(fitRequest$data.mcmc)
    fitRequest$data.mcmc <- as.mcmc(fitRequest$data.param )
    tempPar -> attr(fitRequest$data.mcmc, "mcpar")
    return(hzar.make.dataGroup(data.mcmc=fitRequest$data.mcmc,
                               llFunc=fitRequest$llFunc,
                               data.param=fitRequest$data.param,
                               obsData=fitRequest$obsData,
                               ...))
  }
}
