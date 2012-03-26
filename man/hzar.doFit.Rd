\name{hzar.doFit}
\alias{hzar.doFit}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Run the optimizer.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.doFit(fitRequest)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitRequest}{
%%     ~~Describe \code{fitRequest} here~~
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
  Graham Derryberry \email{asterion@alum.mit.edu}

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
function (fitRequest) 
{
    result <- NULL
    useParam = fitRequest$mcmcParam
    mdlParam = fitRequest$modelParam
    try({
        result <- MCMCmetrop1R(fun = fitRequest$llFunc, logfun = "true", 
            force.samp = TRUE, mcmc = useParam$chainLength, burnin = useParam$burnin, 
            verbose = useParam$verbosity, thin = useParam$thin, 
            theta.init = mdlParam$init, tune = mdlParam$tune, 
            V = fitRequest$cM, seed = useParam$seed, optim.control = list(fnscale = -1, 
                trace = 0, REPORT = 10, maxit = 5000), optim.method = "L-BFGS-B", 
            optim.lower = useParam$lower, optim.upper = useParam$upper)
        colnames(result) <- names(mdlParam$init)
    })
    fitRequest$mcmcRaw <- result
    attr(fitRequest, "fit.run") <- TRUE
    attr(fitRequest, "fit.success") <- FALSE
    if (identical(is.null(result), TRUE)) {
        warning("Fitting failed.")
    }
    else {
        attr(fitRequest, "fit.success") <- TRUE
    }
    return(fitRequest)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
