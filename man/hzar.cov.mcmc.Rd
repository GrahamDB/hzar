\name{hzar.cov.mcmc}
\alias{hzar.cov.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.cov.mcmc(clineLLfunc, mcmcRaw, pDiv = 15, random = 10000, passCenter = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clineLLfunc}{
%%     ~~Describe \code{clineLLfunc} here~~
}
  \item{mcmcRaw}{
%%     ~~Describe \code{mcmcRaw} here~~
}
  \item{pDiv}{
%%     ~~Describe \code{pDiv} here~~
}
  \item{random}{
%%     ~~Describe \code{random} here~~
}
  \item{passCenter}{
%%     ~~Describe \code{passCenter} here~~
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
function (clineLLfunc, mcmcRaw, pDiv = 15, random = 10000, passCenter = FALSE) 
{
    mcmc.nm <- colnames(mcmcRaw)
    pL <- lapply(mcmc.nm, function(x) {
        min(mcmcRaw[, x])
    })
    names(pL) <- mcmc.nm
    pU <- lapply(mcmc.nm, function(x) {
        max(mcmcRaw[, x])
    })
    names(pU) <- mcmc.nm
    return(hzar.cov.rect(clineLLfunc, pL, pU, pDiv, random, passCenter))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
