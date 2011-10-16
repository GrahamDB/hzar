\name{hzar.fit2DataGroup}
\alias{hzar.fit2DataGroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Prepare output from optimizer for analysis.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.fit2DataGroup(fitRequest, doPar = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitRequest}{
%%     ~~Describe \code{fitRequest} here~~
}
  \item{doPar}{
%%     ~~Describe \code{doPar} here~~
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
function (fitRequest, doPar = FALSE) 
{
    if (inherits(fitRequest, "hzar.fitRequest")) {
        if (!identical(attr(fitRequest, "fit.success"), TRUE)) {
            mcmcPsuedo <- mcmc(as.data.frame(fitRequest$modelParam$init))
            return(hzar.make.dataGroup(data.mcmc = mcmcPsuedo, 
                llFunc <- fitRequest$llFunc))
        }
        return(hzar.make.dataGroup(data.mcmc = fitRequest$mcmcRaw, 
            llFunc <- fitRequest$llFunc, doPar = doPar))
    }
    else if (inherits(fitRequest, "hzar.dataGroup")) {
        return(fitRequest)
    }
    else {
        stop(paste("hzar does not understand objects of type ", 
            class(fitRequest)))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
