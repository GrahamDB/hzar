\name{meta.tail.lower}
\alias{meta.tail.lower}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
meta.tail.lower(gamma, d1, tau1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{d1}{
%%     ~~Describe \code{d1} here~~
}
  \item{tau1}{
%%     ~~Describe \code{tau1} here~~
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
function (gamma, d1, tau1) 
{
    lowerStep <- gamma * d1
    lowerA <- eval.intercept(lowerStep)
    lowerK <- eval.kappa(tau1, lowerStep)
    lower.tail <- function(u) {
        lowerU <- u + lowerStep
        return((lowerU < 0) * (lowerA * exp(lowerK * lowerU)))
    }
    return(list(lowerStep = lowerStep, lowerA = lowerA, lowerK = lowerK, 
        lower.tail = lower.tail))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
