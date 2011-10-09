\name{doCLTData1DPops}
\alias{doCLTData1DPops}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
doCLTData1DPops(distance, muObs, varObs, nSamples)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distance}{
%%     ~~Describe \code{distance} here~~
}
  \item{muObs}{
%%     ~~Describe \code{muObs} here~~
}
  \item{varObs}{
%%     ~~Describe \code{varObs} here~~
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
function (distance, muObs, varObs, nSamples) 
{
    if ((length(distance) != length(muObs)) || (length(distance) != 
        length(nSamples)) || (length(distance) != length(varObs))) {
        stop("Distance, muObs, varObs and nSamples are not all of the same length.")
    }
    if (sum(nSamples < 2) > 0) 
        stop("There must be at least two samples per population!")
    if (sum(nSamples < 10) > 0) 
        warning("Some populations have less than 10 samples.  More samples would be nice.")
    if (sum(varObs == 0) > 0) {
        warning("Some of the population samples have a variance of 0.  Adding in estimated variance due to measurement error.")
        m.err <- 5/3 * 10^-quantile(1 + as.numeric(lapply(muObs, 
            function(nV) {
                min((-1:12)[round(nV, digits = -1:12) == nV])
            })), probs = 0.75)[[1]]
        varObs <- varObs + m.err * m.err
    }
    obj <- list(frame = data.frame(dist = distance, obsMean = muObs, 
        obsVariance = varObs, n = nSamples))
    obj$model.LL <- function(model.func) {
        muEst = model.func(obj$frame$dist)
        res <- logLikeCLTEstPopMean(muEst = as.numeric(muEst), 
            muObs = as.numeric(obj$frame$obsMean), varObs = as.numeric(obj$frame$obsVariance), 
            nObs = as.numeric(obj$frame$n))
        result <- sum(res)
        if (is.na(result)) 
            return(-1e+08)
        return(result)
    }
    class(obj) <- c("clineSampleData1DCLT", "hzar.obsData")
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
