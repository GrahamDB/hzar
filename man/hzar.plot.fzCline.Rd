\name{hzar.plot.fzCline}
\alias{hzar.plot.fzCline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot the 95\% credible cline region for the given locus model.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.plot.fzCline(dataGroup, fzCline = hzar.getCredParamRed(dataGroup), type = "p", pch = "+", col = "black", fzCol = "gray", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataGroup}{
%%     ~~Describe \code{dataGroup} here~~
}
  \item{fzCline}{
%%     ~~Describe \code{fzCline} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{pch}{
%%     ~~Describe \code{pch} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{fzCol}{
%%     ~~Describe \code{fzCol} here~~
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
function (dataGroup, fzCline = hzar.getCredParamRed(dataGroup), 
    type = "p", pch = "+", col = "black", fzCol = "gray", ...) 
{
    hzar.plot.obsData(dataGroup, col = "transparent", ...)
    xSeries <- seq(from = par("usr")[1], to = par("usr")[2], 
        length.out = 109)
    if (par("xaxs") == "r") 
        xSeries <- xSeries[2:108]
    fzCoor <- fzCline$fzCline(xSeries)
    polygon(x = c(fzCoor$x, rev(fzCoor$x)), y = c(fzCoor$yMin, 
        rev(fzCoor$yMax)), border = fzCol, col = fzCol)
    lines(x = xSeries, y = dataGroup$ML.cline$clineFunc(xSeries), 
        col = col)
    hzar.plot.obsData(dataGroup, col = col, type = type, pch = pch, 
        add = TRUE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
