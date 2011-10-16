\name{hzar.makeCline1DCLT}
\alias{hzar.makeCline1DCLT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Make a cline model for the prediction of population means.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.makeCline1DCLT(data = NULL, scaling = "free", tails = "none", direction = NULL)
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
function (data = NULL, scaling = "free", tails = "none", direction = NULL) 
{
    if (identical(tolower(tails), "none")) {
        return(buildCline1D(data, scaling, direction, hzar.meta.CLTnA, 
            hzar.meta.CLTnD))
    }
    else if (identical(tolower(tails), "right")) {
        myRightCline <- buildCline1D(data, scaling, direction, 
            hzar.meta.CLTrA, hzar.meta.CLTrD)
        attr(myRightCline, "tails") <- "right"
        return(myRightCline)
    }
    stop(paste("Cline with", tails, "tail(s) not available."))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
