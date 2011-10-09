\name{hzar.make.dataGroup}
\alias{hzar.make.dataGroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.make.dataGroup(data.mcmc, llFunc, ML.cline = NULL, doPar = FALSE, data.LL = hzar.eval.clineLL(llFunc = llFunc, data = data.mcmc, doPar = doPar), data.param = as.data.frame(data.mcmc), obsData = hzar.extract.obsData(llFunc))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.mcmc}{
%%     ~~Describe \code{data.mcmc} here~~
}
  \item{llFunc}{
%%     ~~Describe \code{llFunc} here~~
}
  \item{ML.cline}{
%%     ~~Describe \code{ML.cline} here~~
}
  \item{doPar}{
%%     ~~Describe \code{doPar} here~~
}
  \item{data.LL}{
%%     ~~Describe \code{data.LL} here~~
}
  \item{data.param}{
%%     ~~Describe \code{data.param} here~~
}
  \item{obsData}{
%%     ~~Describe \code{obsData} here~~
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
function (data.mcmc, llFunc, ML.cline = NULL, doPar = FALSE, 
    data.LL = hzar.eval.clineLL(llFunc = llFunc, data = data.mcmc, 
        doPar = doPar), data.param = as.data.frame(data.mcmc), 
    obsData = hzar.extract.obsData(llFunc)) 
{
    result <- list(llFunc = llFunc, data.mcmc = data.mcmc, data.param = data.param, 
        data.LL = data.LL, ML.cline = ML.cline, obsData = obsData)
    class(result) <- "hzar.dataGroup"
    result$ML.cline <- hzar.get.ML.cline(result)
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
