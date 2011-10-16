\name{hzar.next.fitRequest}
\alias{hzar.next.fitRequest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a new fitRequest using data from another fitRequest.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.next.fitRequest(oldFitRequest)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oldFitRequest}{
%%     ~~Describe \code{oldFitRequest} here~~
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
function (oldFitRequest) 
{
    seedChannel <- 1
    if (is.list(oldFitRequest$mcmcParam$seed)) {
        seedChannel = oldFitRequest$mcmcParam$seed[[2]]
    }
    if (identical(attr(oldFitRequest, "fit.run"), TRUE)) {
        seedChannel <- seedChannel + 1
    }
    else {
        seedChannel <- seedChannel + 10
    }
    mcmcParam = hzar.make.mcmcParam(oldFitRequest$mcmcParam$chainLength, 
        oldFitRequest$mcmcParam$burnin, oldFitRequest$mcmcParam$verbosity, 
        oldFitRequest$mcmcParam$thin, seedChannel)
    mdlParam <- oldFitRequest$modelParam
    covMatrix <- oldFitRequest$cM
    if (identical(attr(oldFitRequest, "fit.success"), TRUE)) {
        mcmcSubset <- oldFitRequest$mcmcRaw[sample(dim(oldFitRequest$mcmcRaw)[[1]]), 
            ]
        subLL <- hzar.eval.clineLL(mcmcSubset, oldFitRequest$llFunc)
        covData <- hzar.cov.mcmc(oldFitRequest$llFunc, mcmcSubset[subLL > 
            max(subLL - 4), ], passCenter = TRUE)
        covMatrix <- covData$cov
        new.center <- covData$center[names(mdlParam$init)]
        if (oldFitRequest$llFunc(new.center) > 1e-06) 
            mdlParam$init <- new.center
    }
    return(hzar.make.fitRequest(mdlParam, covMatrix, oldFitRequest$llFunc, 
        mcmcParam))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
