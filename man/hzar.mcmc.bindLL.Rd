\name{hzar.mcmc.bindLL}
\alias{hzar.mcmc.bindLL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a data frame with sampled parameters and log likelihoods.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.mcmc.bindLL(fitRequest, dataGroup = hzar.fit2DataGroup(fitRequest), mcmcData = fitRequest$mcmcRaw, llFunc = dataGroup$llFunc, llData = dataGroup$data.LL, t0 = start(mcmcData), tF = thin(mcmcData))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitRequest}{
%%     ~~Describe \code{fitRequest} here~~
}
  \item{dataGroup}{
%%     ~~Describe \code{dataGroup} here~~
}
  \item{mcmcData}{
%%     ~~Describe \code{mcmcData} here~~
}
  \item{llFunc}{
%%     ~~Describe \code{llFunc} here~~
}
  \item{llData}{
%%     ~~Describe \code{llData} here~~
}
  \item{t0}{
%%     ~~Describe \code{t0} here~~
}
  \item{tF}{
%%     ~~Describe \code{tF} here~~
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
function (fitRequest, dataGroup = hzar.fit2DataGroup(fitRequest), 
    mcmcData = fitRequest$mcmcRaw, llFunc = dataGroup$llFunc, 
    llData = dataGroup$data.LL, t0 = start(mcmcData), tF = thin(mcmcData)) 
{
    data <- cbind(mcmcData, llData)
    result <- mcmc(data = as.matrix(data), start = t0, thin = tF)
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
