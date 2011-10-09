\name{meta.tail.upper}
\alias{meta.tail.upper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
meta.tail.upper(gamma, d2, tau2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{d2}{
%%     ~~Describe \code{d2} here~~
}
  \item{tau2}{
%%     ~~Describe \code{tau2} here~~
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
function (gamma, d2, tau2) 
{
    upperStep <- gamma * d2
    upperA <- eval.intercept(upperStep)
    upperK <- eval.kappa(tau2, upperStep)
    upper.tail <- function(u) {
        upperU <- u - upperStep
        return((upperU > 0) * (1 - upperA * exp(-upperK * upperU)))
    }
    return(list(upperStep = upperStep, upperA = upperA, upperK = upperK, 
        upper.tail = upper.tail))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
