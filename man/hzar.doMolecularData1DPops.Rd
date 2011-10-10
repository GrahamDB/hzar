\name{doMolecularData1DPops}
\alias{doMolecularData1DPops}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
doMolecularData1DPops(distance, pObs, nSamples)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distance}{
%%     ~~Describe \code{distance} here~~
}
  \item{pObs}{
%%     ~~Describe \code{pObs} here~~
}
  \item{nSamples}{
%%     ~~Describe \code{nSamples} here~~
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
function (distance, pObs, nSamples) 
{
    if ((length(distance) != length(pObs)) || (length(distance) != 
        length(nSamples)) || (length(pObs) != length(nSamples))) {
        stop("Distance, pObs and nSamples are not of the same length.")
    }
    obj <- list(frame = data.frame(dist = distance, obsFreq = pObs, 
        n = nSamples))
    obj$model.LL <- function(model.func) {
        pEst = model.func(obj$frame$dist)
        res <- sampleLikelihoodMolecularPop(pEst = as.numeric(pEst), 
            pObs = as.numeric(obj$frame$obsFreq), N = as.numeric(obj$frame$n))
        result <- sum(res)
        if (is.na(result)) 
            return(-1e+08)
        return(result)
    }
    class(obj) <- c("clineSampleData1D", "hzar.obsData")
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
