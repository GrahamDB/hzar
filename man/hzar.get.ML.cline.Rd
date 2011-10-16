\name{hzar.get.ML.cline}
\alias{hzar.get.ML.cline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extract the maximum likelihood cline.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.get.ML.cline(fitRequest)
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
    if (inherits(fitRequest, "hzar.fitRequest")) {
        if (!identical(attr(fitRequest, "fit.success"), TRUE)) {
            stop("fitRequest needs a successful MCMC sampling chain.")
        }
        data.LL <- hzar.eval.clineLL(llFunc = fitRequest$llFunc, 
            data = fitRequest$mcmcRaw)
        data.param <- as.data.frame(fitRequest$mcmcRaw)
    }
    else if (inherits(fitRequest, "hzar.dataGroup")) {
        if (!identical(is.null(fitRequest$ML.cline), TRUE)) 
            return(fitRequest$ML.cline)
        data.LL <- fitRequest$data.LL
        data.param <- fitRequest$data.param
    }
    else {
        stop(paste("hzar does not understand objects of type ", 
            class(fitRequest)))
    }
    param.free = data.param[data.LL == max(data.LL), ][1, ]
    return(hzar.gen.cline(param.free, fitRequest))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
