\name{splitParameters}
\alias{splitParameters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
splitParameters(paramList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{paramList}{
%%     ~~Describe \code{paramList} here~~
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
function (paramList) 
{
    param.fixed = list()
    param.free.init = list()
    param.free.weight = list()
    param.free.lower = list()
    param.free.upper = list()
    for (iter in seq(along = paramList)) {
        if (class(paramList[[iter]]) != "clineParameter") 
            stop(paste("Item", paramList[[iter]], "not a cline parameter!"))
        parName <- attr(paramList[[iter]], "param")
        parValue <- paramList[[iter]]$val
        if (attr(paramList[[iter]], "fixed")) {
            param.fixed[[parName]] <- parValue
        }
        else {
            param.free.init[[parName]] <- parValue
            param.free.weight[[parName]] <- paramList[[iter]]$w
            param.free.lower[[parName]] <- attr(paramList[[iter]], 
                "limit.lower")
            param.free.upper[[parName]] <- attr(paramList[[iter]], 
                "limit.upper")
        }
    }
    param <- list(init = param.free.init, tune = param.free.weight, 
        lower = param.free.lower, upper = param.free.upper, fixed = param.fixed)
    return(param)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
