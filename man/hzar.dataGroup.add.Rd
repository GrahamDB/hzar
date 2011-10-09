\name{hzar.dataGroup.add}
\alias{hzar.dataGroup.add}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.dataGroup.add(dataGroup, fitRequest = list(), doPar = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataGroup}{
%%     ~~Describe \code{dataGroup} here~~
}
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
function (dataGroup, fitRequest = list(), doPar = FALSE) 
{
    if (inherits(dataGroup, c("hzar.fitRequest", "hzar.dataGroup"))) {
        dataGroup <- hzar.fit2DataGroup(dataGroup, doPar = doPar)
    }
    else if (is.list(dataGroup) && (length(fitRequest) == 0)) {
        if (length(dataGroup) == 0) 
            return(NULL)
        oldFitRequest <- dataGroup
        dataGroup <- hzar.fit2DataGroup(dataGroup[[1]], doPar = doPar)
        if (length(oldFitRequest) > 1) 
            fitRequest <- lapply(2:length(oldFitRequest), function(x) oldFitRequest[[x]])
    }
    else {
        stop("dataGroup not of apropriate type.")
    }
    if (inherits(fitRequest, c("hzar.fitRequest", "hzar.dataGroup"))) {
        fitRequest <- hzar.fit2DataGroup(fitRequest, doPar = doPar)
    }
    else if (is.list(fitRequest)) {
        lapply(fitRequest, function(x) dataGroup <<- hzar.dataGroup.add(dataGroup, 
            x, doPar = doPar))
        return(dataGroup)
    }
    else {
        stop("fitRequest not of appropriate type.")
    }
    if (!hzar.sameModel(dataGroup, fitRequest)) 
        stop("dataGroup and fitRequest must use the same model.")
    Left.Cline <- hzar.get.ML.cline(dataGroup)
    Right.Cline <- hzar.get.ML.cline(fitRequest)
    if (Right.Cline$logLike > Left.Cline$logLike) 
        Left.Cline <- Right.Cline
    return(hzar.make.dataGroup(data.mcmc = rbind(dataGroup$data.mcmc, 
        fitRequest$data.mcmc), llFunc = dataGroup$llFunc, ML.cline = Left.Cline, 
        data.LL = rbind(dataGroup$data.LL, fitRequest$data.LL), 
        data.param = rbind(dataGroup$data.param, fitRequest$data.param), 
        obsData = dataGroup$obsData))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
