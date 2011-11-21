\name{hzar.doMorphoSets}
\alias{hzar.doMorphoSets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Make hzar.obsData objects from scalar observations using reference tables 
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.doMorphoSets(traitNames, tDist, tDLocCol, tDDistCol, tValues, tVLocCol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{traitNames}{
%%     ~~Describe \code{traitNames} here~~
}
  \item{tDist}{
%%     ~~Describe \code{tDist} here~~
}
  \item{tDLocCol}{
%%     ~~Describe \code{tDLocCol} here~~
}
  \item{tDDistCol}{
%%     ~~Describe \code{tDDistCol} here~~
}
  \item{tValues}{
%%     ~~Describe \code{tValues} here~~
}
  \item{tVLocCol}{
%%     ~~Describe \code{tVLocCol} here~~
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
function (traitNames, tDist, tDLocCol, tDDistCol, tValues, tVLocCol) 
{
    distOfLocation <- tDist[[tDDistCol]]
    names(distOfLocation) <- tDist[[tDLocCol]]
    res <- lapply(traitNames, function(id) {
        tVal <- tValues[!is.na(tValues[[id]]), c(tVLocCol, id)]
        return(hzar.makeTraitObsData(distOfLocation, locationOfValue = tVal[[tVLocCol]], 
            values = tVal[[id]]))
    })
    names(res) <- traitNames
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
