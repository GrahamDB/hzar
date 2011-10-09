\name{cline.samp.dist}
\alias{cline.samp.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cline.samp.dist(cline.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cline.data}{
%%     ~~Describe \code{cline.data} here~~
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
function (cline.data) 
{
    dist.pool <- unique(sort(as.numeric(cline.data$dist)))
    dist.init = quantile(dist.pool, probs = c(0.25, 0.75))
    potLeftOut = length(dist.pool[dist.pool <= dist.init[[1]]])
    potRightOut = length(dist.pool[dist.pool >= dist.init[[2]]]) - 
        1
    potInner = length(dist.pool[(dist.pool < dist.init[[2]]) & 
        (dist.pool > dist.init[[1]])])
    dist.len <- length(dist.pool)
    distLeftSeries <- numeric(1)
    distLeftSeries[1] <- dist.init[[1]]
    distRightSeries <- numeric(1)
    distRightSeries[1] <- dist.init[[2]]
    if (potLeftOut > 2) 
        distLeftSeries <- c(dist.pool[2:potLeftOut], distLeftSeries)
    if (potRightOut > 1) 
        distRightSeries <- c(distRightSeries, dist.pool[(dist.len - 
            potRightOut):(dist.len - 1)])
    if (potInner > 1) {
        distLeftSeries <- c(distLeftSeries, dist.pool[potLeftOut + 
            1])
        distRightSeries <- c(dist.pool[dist.len - potRightOut - 
            1], distRightSeries)
    }
    return(list(unique(distLeftSeries), unique(distRightSeries)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
