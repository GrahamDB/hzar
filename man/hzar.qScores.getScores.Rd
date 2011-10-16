\name{hzar.qScores.getScores}
\alias{hzar.qScores.getScores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate credibility scores.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.qScores.getScores(xSeries, dSeries, raw.x, raw.wt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xSeries}{
%%     ~~Describe \code{xSeries} here~~
}
  \item{dSeries}{
%%     ~~Describe \code{dSeries} here~~
}
  \item{raw.x}{
%%     ~~Describe \code{raw.x} here~~
}
  \item{raw.wt}{
%%     ~~Describe \code{raw.wt} here~~
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
function (xSeries, dSeries, raw.x, raw.wt) 
{
    if (length(xSeries) > 1000) {
        cutValue <- as.integer(length(xSeries)/2) + 1
        return(c(hzar.qScores.getScores(xSeries[1:cutValue], 
            dSeries[1:cutValue], raw.x[raw.x <= xSeries[[cutValue]]], 
            raw.wt[raw.x <= xSeries[[cutValue]]]), hzar.qScores.getScores(xSeries[cutValue:length(xSeries)], 
            dSeries[cutValue:length(xSeries)], raw.x[raw.x >= 
                xSeries[[cutValue]]], raw.wt[raw.x >= xSeries[[cutValue]]])))
    }
    getScore <- function(x, xS = xSeries, dS = dSeries, rx = raw.x, 
        rwt = raw.wt) {
        data = c(rwt[xS[[x]] == rx], rwt[xS[[x + 1]] == rx])
        return(mean(exp(data)) * dS[[x]])
    }
    return(sapply(1:(length(xSeries) - 1), getScore))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
