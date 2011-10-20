\name{hzar.make.fzCline}
\alias{hzar.make.fzCline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compile a list of hzar.cline objects into a hzar.fzCline object
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.make.fzCline(clineList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clineList}{
%%     ~~Describe \code{clineList} here~~
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
function (clineList) 
{
    res <- list(clines = clineList, listFuncInt = function(xVal, 
        funcList = res$clines) {
        yList <- as.numeric(lapply(funcList, function(x, u) x$clineFunc(u), 
            u = xVal))
        return(data.frame(x = xVal, yMin = min(yList), yMax = max(yList)))
    }, fzCline = function(xVals, listFunc = res$listFuncInt) {
        xVals <- as.numeric(xVals)
        if (length(xVals > 1)) return(do.call(rbind, lapply(xVals, 
            listFunc)))
        if (length(xVals < 1)) return(numeric(0))
        return(listFunc(xVals))
    })
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
