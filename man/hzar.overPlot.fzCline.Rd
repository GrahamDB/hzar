\name{hzar.overPlot.fzCline}
\alias{hzar.overPlot.fzCline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot a cline region for multiple models and / or loci. 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.overPlot.fzCline(dataGroupSet,
fzClineSet = sapply(
dataGroupSet, hzar.getCredParamRed, simplify = FALSE),
type = "p",
fzDens = 8,
fzShadeAngle = ((1:length(dataGroupSet)) * 180)
\%/\% (1 + length(dataGroupSet)),
...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataGroupSet}{
%%     ~~Describe \code{dataGroupSet} here~~
}
  \item{fzClineSet}{
%%     ~~Describe \code{fzClineSet} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{fzDens}{
%%     ~~Describe \code{fzDens} here~~
}
  \item{fzShadeAngle}{
%%     ~~Describe \code{fzShadeAngle} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (dataGroupSet, fzClineSet = sapply(dataGroupSet, hzar.getCredParamRed, 
    simplify = FALSE), type = "p", fzDens = 8, fzShadeAngle = ((1:length(dataGroupSet)) * 
    180)\%/\%(1 + length(dataGroupSet)), ...) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
