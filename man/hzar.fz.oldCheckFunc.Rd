\name{hzar.fz.oldCheckFunc}
\alias{hzar.fz.oldCheckFunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.fz.oldCheckFunc(paramNames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{paramNames}{
%%     ~~Describe \code{paramNames} here~~
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
function (paramNames) 
{
    checkCW <- function(accepted, candidate) {
        if (sum(accepted$center <= candidate$center & accepted$width >= 
            candidate$width) == 0) 
            return(TRUE)
        if (sum(accepted$center >= candidate$center & accepted$width <= 
            candidate$width) == 0) 
            return(TRUE)
        if (sum(accepted$center >= candidate$center & accepted$width >= 
            candidate$width) == 0) 
            return(TRUE)
        if (sum(accepted$center <= candidate$center & accepted$width <= 
            candidate$width) == 0) 
            return(TRUE)
        return(FALSE)
    }
    checkCWpMM <- function(accepted, candidate) {
        if (sum(accepted$center <= candidate$center & accepted$width >= 
            candidate$width & accepted$pMin >= candidate$pMin & 
            accepted$pMax <= candidate$pMax) == 0) 
            return(TRUE)
        if (sum(accepted$center >= candidate$center & accepted$width <= 
            candidate$width & accepted$pMin <= candidate$pMin & 
            accepted$pMax >= candidate$pMax) == 0) 
            return(TRUE)
        if (sum(accepted$center >= candidate$center & accepted$width >= 
            candidate$width & accepted$pMin >= candidate$pMin & 
            accepted$pMax <= candidate$pMax) == 0) 
            return(TRUE)
        if (sum(accepted$center <= candidate$center & accepted$width <= 
            candidate$width & accepted$pMin <= candidate$pMin & 
            accepted$pMax >= candidate$pMax) == 0) 
            return(TRUE)
        return(FALSE)
    }
    if (identical(sort(paramNames), c("center", "width"))) {
        return(checkCW)
    }
    if (identical(sort(paramNames), c("center", "pMax", "pMin", 
        "width"))) {
        return(checkCWpMM)
    }
    return(function(accepted, candidate) return(TRUE))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
