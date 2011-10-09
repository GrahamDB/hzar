\name{hzar.plot.cline}
\alias{hzar.plot.cline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.plot.cline(cline, add = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cline}{
%%     ~~Describe \code{cline} here~~
}
  \item{add}{
%%     ~~Describe \code{add} here~~
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
function (cline, add = FALSE, ...) 
{
    if (inherits(cline, "hzar.cline")) 
        curve(cline$clineFunc(x), add = add, ...)
    if (inherits(cline, c("hzar.dataGroup", "hzar.fitRequest"))) {
        dataGroup <- hzar.fit2DataGroup(cline)
        hzar.plot.obsData(dataGroup, add = add, ...)
        hzar.plot.cline(hzar.get.ML.cline(dataGroup), add = TRUE, 
            ...)
    }
    if (inherits(cline, c("hzar.obsDataGroup"))) {
        hzar.plot.obsData(cline, add = add, ...)
        lapply(cline$data.groups, function(dataGroup) hzar.plot.cline(hzar.get.ML.cline(dataGroup), 
            add = TRUE, ...))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
