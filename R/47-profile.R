
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
    bquote(.(cLT) <- sum( .(sampleMean)*.(nEff)*.(fTa)/ .(varT))),
    bquote(.(cRT) <- sum( .(sampleMean)*.(nEff)*.(fTd)/ .(varT))),
    bquote(.(aLT) <- sum( .(nEff)*.(fTa)^2/ .(varT))),
    bquote(.(aRT) <- sum( .(nEff)*.(fTd)^2/ .(varT))),
    bquote(.(aMT) <- sum( .(nEff)*.(fTd)*.(fTa)/ .(varT))),
    bquote(.(tDen) <- .(aRT)*.(aLT)-.(aMT)^2),
    bquote(.(tMuL) <- (.(cLT)*.(aRT)-.(cRT)*.(aMT))/.(tDen)),
    bquote(.(tMuR) <- (.(cRT)*.(aLT)-.(cLT)*.(aMT))/.(tDen))
    )
}

mu2fExp <- function(muExp, tMuL=quote(muL), tMuR=quote(muR))
  simplify.exp(eval(substitute(substitute(a,b),list(a=muExp,b=list(muL=0,muR=1)))))

