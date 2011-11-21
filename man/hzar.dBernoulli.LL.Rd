\name{hzar.dBernoulli.LL}
\alias{hzar.dBernoulli.LL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Transformations of Scalar Data into Bernoulli Trials 
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.dBernoulli.LL(values, locations, getMax = FALSE, getProbs = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{values}{
%%     ~~Describe \code{values} here~~
}
  \item{locations}{
%%     ~~Describe \code{locations} here~~
}
  \item{getMax}{
%%     ~~Describe \code{getMax} here~~
}
  \item{getProbs}{
%%     ~~Describe \code{getProbs} here~~
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
function (values, locations, getMax = FALSE, getProbs = FALSE) 
{
    midpoints <- function(x) {
        y <- as.numeric(sort(unique(x)))
        return((y[y > min(y)] + y[y < max(y)])/2)
    }
    mVal <- midpoints(values)
    bernoulli.LL <- function(nSuccess, nTotal) {
        nFailures <- nTotal - nSuccess
        hSuccess <- ifelse(nSuccess > 0, nSuccess * log(nTotal/nSuccess), 
            0)
        hFailure <- ifelse(nFailures > 0, nFailures * log(nTotal/nFailures), 
            0)
        return(hSuccess + hFailure)
    }
    pL <- lapply(unique(locations), function(x) as.numeric(values[locations == 
        x]))
    bLLSum <- function(junk, procList) {
        nSamples = as.numeric(lapply(procList, length))
        return(sum(bernoulli.LL(as.numeric(lapply(procList, function(x, 
            cutValue) sum(x < cutValue), cutValue = junk)), nSamples)))
    }
    if (getProbs) {
        junk <- as.numeric(lapply(mVal, bLLSum, pL))
        cutValue <- mVal[junk == max(junk)][1]
        nSamples <- as.numeric(lapply(pL, length))
        nLow <- as.numeric(lapply(pL, function(x, cutValue) sum(x < 
            cutValue), cutValue))
        return(data.frame(locID = unique(locations), pObs = nLow/nSamples, 
            nLow = nLow, nSamples = nSamples, cutValue = cutValue))
    }
    if (getMax) {
        junk <- as.numeric(lapply(mVal, bLLSum, pL))
        return(mVal[junk == max(junk)][1])
    }
    return(as.numeric(lapply(mVal, bLLSum, pL)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
