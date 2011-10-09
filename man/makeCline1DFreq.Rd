\name{makeCline1DFreq}
\alias{makeCline1DFreq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeCline1DFreq(data = NULL, scaling = "none", tails = "none", direction = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{scaling}{
%%     ~~Describe \code{scaling} here~~
}
  \item{tails}{
%%     ~~Describe \code{tails} here~~
}
  \item{direction}{
%%     ~~Describe \code{direction} here~~
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
function (data = NULL, scaling = "none", tails = "none", direction = NULL) 
{
    if (identical(tolower(tails), "none")) {
        return(makeSimpleCline1D(data, scaling, direction))
    }
    else if (identical(tolower(tails), "both")) {
        return(makeTailedCline1D(data, scaling, direction))
    }
    else if (identical(tolower(tails), "right")) {
        myRightCline <- buildCline1D(data, scaling, direction, 
            hzar.meta.rtail.scaled.ascending, hzar.meta.rtail.scaled.descending)
        attr(myRightCline, "tails") <- "right"
        return(myRightCline)
    }
    else if (identical(tolower(tails), "left")) {
        myLeftCline <- buildCline1D(data, scaling, direction, 
            hzar.meta.ltail.scaled.ascending, hzar.meta.ltail.scaled.descending)
        attr(myLeftCline, "tails") <- "left"
        return(myLeftCline)
    }
    else if (identical(tolower(tails), "mirror")) {
        myMirrorCline <- buildCline1D(data, scaling, direction, 
            hzar.meta.mtail.scaled.ascending, hzar.meta.mtail.scaled.descending)
        attr(myMirrorCline, "tails") <- "mirror"
        return(myMirrorCline)
    }
    stop(paste("Cline with", tails, "tail(s) not available."))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
