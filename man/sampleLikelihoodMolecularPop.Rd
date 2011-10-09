\name{sampleLikelihoodMolecularPop}
\alias{sampleLikelihoodMolecularPop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampleLikelihoodMolecularPop(pEst, pObs, N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pEst}{
%%     ~~Describe \code{pEst} here~~
}
  \item{pObs}{
%%     ~~Describe \code{pObs} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
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
function (pEst, pObs, N) 
{
    p1 = N * log(pEst)
    p0 = N * log(1 - pEst)
    lP = ((1 - pObs) * log((1 - pEst)/(1 - pObs)) + pObs * log(pEst/pObs)) * 
        N
    tmp = ifelse(pObs == 1, p1, lP)
    res = ifelse(pObs == 0, p0, tmp)
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
