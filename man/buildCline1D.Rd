\name{buildCline1D}
\alias{buildCline1D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
buildCline1D(data, scaling, direction, ascending.cline, descending.cline)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{scaling}{
%%     ~~Describe \code{scaling} here~~
}
  \item{direction}{
%%     ~~Describe \code{direction} here~~
}
  \item{ascending.cline}{
%%     ~~Describe \code{ascending.cline} here~~
}
  \item{descending.cline}{
%%     ~~Describe \code{descending.cline} here~~
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
function (data, scaling, direction, ascending.cline, descending.cline) 
{
    dist = NULL
    obs = NULL
    if (is.null(direction)) {
        if (is.null(data)) 
            stop("Either provide observation data, or set direction.")
        if (inherits(data, "clineSampleData1D")) {
            if (mean(data$frame$dist) < (sum(data$frame$dist * 
                data$frame$obsFreq)/sum(data$frame$obsFreq))) {
                direction = "ascending"
            }
            else {
                direction = "descending"
            }
            dist = data$frame$dist
            obs = data$frame$obsFreq
        }
        else if (inherits(data, "clineSampleData1DCLT")) {
            if (mean(data$frame$dist) < (sum(data$frame$dist * 
                data$frame$obsMean)/sum(data$frame$obsMean))) {
                direction = "ascending"
            }
            else {
                direction = "descending"
            }
            dist = data$frame$dist
            obs = data$frame$obsMean
        }
        else {
            stop("Class of data unknown.")
        }
    }
    if (identical(tolower(direction), "ascending")) {
        myModel = ascending.cline
    }
    else if (identical(tolower(direction), "descending")) {
        myModel = descending.cline
    }
    else {
        stop(paste("Direction", direction, "not recognized. Specify either ascending or descending."))
    }
    myModel = setupMoleCenterClineParameters(myModel, scaling, 
        dist, obs)
    attr(myModel, "tails") <- "none"
    return(myModel)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
