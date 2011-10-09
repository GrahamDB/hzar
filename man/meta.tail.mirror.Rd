\name{meta.tail.mirror}
\alias{meta.tail.mirror}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
meta.tail.mirror(gamma, d, tau)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
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
function (gamma, d, tau) 
{
    step <- gamma * d
    A <- eval.intercept(step)
    K <- eval.kappa(tau, step)
    lower.tail <- function(u) {
        lowerU <- u + step
        return((lowerU < 0) * (A * exp(K * lowerU)))
    }
    upper.tail <- function(u) {
        upperU <- u - step
        return((upperU > 0) * (1 - A * exp(-K * upperU)))
    }
    return(list(lowerStep = step, lowerA = A, lowerK = K, lower.tail = lower.tail, 
        upperStep = step, upperA = A, upperK = K, upper.tail = upper.tail))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
