\name{hzar.make.fitRequest}
\alias{hzar.make.fitRequest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a hzar.fitRequest with the specified attributes.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.make.fitRequest(modelParameters, covMatrix, clineLLfunc, mcmcParameters, mcmcRaw = NULL, fit.run = FALSE, fit.success = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelParameters}{
%%     ~~Describe \code{modelParameters} here~~
}
  \item{covMatrix}{
%%     ~~Describe \code{covMatrix} here~~
}
  \item{clineLLfunc}{
%%     ~~Describe \code{clineLLfunc} here~~
}
  \item{mcmcParameters}{
%%     ~~Describe \code{mcmcParameters} here~~
}
  \item{mcmcRaw}{
%%     ~~Describe \code{mcmcRaw} here~~
}
  \item{fit.run}{
%%     ~~Describe \code{fit.run} here~~
}
  \item{fit.success}{
%%     ~~Describe \code{fit.success} here~~
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
function (modelParameters, covMatrix, clineLLfunc, mcmcParameters, 
    mcmcRaw = NULL, fit.run = FALSE, fit.success = FALSE) 
{
    fitRequest <- list(modelParam = modelParameters, cM = covMatrix, 
        llFunc = clineLLfunc, mcmcParam = mcmcParameters, mcmcRaw = mcmcRaw)
    class(fitRequest) <- "hzar.fitRequest"
    attr(fitRequest, "fit.run") <- fit.run
    attr(fitRequest, "fit.success") <- fit.success
    return(fitRequest)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
