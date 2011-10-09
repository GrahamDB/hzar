\name{hzar.fz.getCheckFunc}
\alias{hzar.fz.getCheckFunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.fz.getCheckFunc(paramNames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{paramNames}{
%%     ~~Describe \code{paramNames} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (paramNames) 
{
    accLeft <- quote(acc$center <= can$center)
    accRight <- quote(acc$center >= can$center)
    accWIn <- quote(acc$width >= can$width)
    accWOut <- quote(acc$width <= can$width)
    accMMOut <- quote(acc$pMin <= can$pMin & acc$pMax >= can$pMax)
    accMMIn <- quote(acc$pMin >= can$pMin & acc$pMax <= can$pMax)
    accLTIn <- quote(acc$deltaL <= can$deltaL & acc$tauL <= can$tauL)
    accLTOut <- quote(acc$deltaL >= can$deltaL & acc$tauL >= 
        can$tauL)
    accRTIn <- quote(acc$deltaR <= can$deltaR & acc$tauR <= can$tauR)
    accRTOut <- quote(acc$deltaR >= can$deltaR & acc$tauR >= 
        can$tauR)
    accMTIn <- quote(acc$deltaM <= can$deltaM & acc$tauM <= can$tauM)
    accMTOut <- quote(acc$deltaM >= can$deltaM & acc$tauM >= 
        can$tauM)
    check <- quote(if (all(!((accScale) & (accLoc)))) return(TRUE))
    checkFunc <- function(acc, can) return(TRUE)
    if ((!("center" \%in\% paramNames)) || (!("width" \%in\% paramNames))) 
        return(checkFunc)
    accOut <- accWOut
    accIn <- accWIn
    if (("pMin" \%in\% paramNames) && ("pMax" \%in\% paramNames)) {
        accOut <- substitute(accOld & accEx, list(accOld = accOut, 
            accEx = accMMOut))
        accIn <- substitute(accOld & accEx, list(accOld = accIn, 
            accEx = accMMIn))
    }
    else if (("pMin" \%in\% paramNames) && (!("pMax" \%in\% paramNames)) || 
        ("pMax" \%in\% paramNames) && (!("pMin" \%in\% paramNames))) {
        return(checkFunc)
    }
    if (("deltaM" \%in\% paramNames) && ("tauM" \%in\% paramNames)) {
        accOut <- substitute(accOld & accEx, list(accOld = accOut, 
            accEx = accMTOut))
        accIn <- substitute(accOld & accEx, list(accOld = accIn, 
            accEx = accMTIn))
    }
    else if (("deltaM" \%in\% paramNames) && (!("tauM" \%in\% paramNames)) || 
        ("tauM" \%in\% paramNames) && (!("deltaM" \%in\% paramNames))) {
        return(checkFunc)
    }
    if (("deltaL" \%in\% paramNames) && ("tauL" \%in\% paramNames)) {
        accLOut <- substitute(accOld & accEx, list(accOld = accOut, 
            accEx = accLTOut))
        accLIn <- substitute(accOld & accEx, list(accOld = accIn, 
            accEx = accLTIn))
    }
    else if (("deltaL" \%in\% paramNames) && (!("tauL" \%in\% paramNames)) || 
        ("tauL" \%in\% paramNames) && (!("deltaL" \%in\% paramNames))) {
        return(checkFunc)
    }
    else {
        accLOut <- accOut
        accLIn <- accIn
    }
    cLO <- substitute(substitute(c, list(accLoc = accRight, accScale = accLOut)), 
        list(c = check))
    cLI <- substitute(substitute(c, list(accLoc = accLeft, accScale = accLIn)), 
        list(c = check))
    if (("deltaR" \%in\% paramNames) && ("tauR" \%in\% paramNames)) {
        accROut <- substitute(accOld & accEx, list(accOld = accOut, 
            accEx = accRTOut))
        accRIn <- substitute(accOld & accEx, list(accOld = accIn, 
            accEx = accRTIn))
    }
    else if (("deltaR" \%in\% paramNames) && (!("tauR" \%in\% paramNames)) || 
        ("tauR" \%in\% paramNames) && (!("deltaR" \%in\% paramNames))) {
        return(checkFunc)
    }
    else {
        accROut <- accOut
        accRIn <- accIn
    }
    cRO <- substitute(substitute(c, list(accLoc = accLeft, accScale = accROut)), 
        list(c = check))
    cRI <- substitute(substitute(c, list(accLoc = accRight, accScale = accRIn)), 
        list(c = check))
    checkE <- do.call(expression, list(eval(cLI), eval(cLO), 
        eval(cRI), eval(cRO), quote(return(FALSE))))
    body(checkFunc) <- as.call(c(as.name("{"), checkE))
    return(checkFunc)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
