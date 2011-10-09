\name{hzar.copyModelLabels}
\alias{hzar.copyModelLabels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.copyModelLabels(group1, group2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{group1}{
%%     ~~Describe \code{group1} here~~
}
  \item{group2}{
%%     ~~Describe \code{group2} here~~
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
function (group1, group2) 
{
    if (inherits(group1, c("hzar.obsDataGroup"))) 
        group1 <- group1$data.groups
    if (inherits(group2, c("hzar.obsDataGroup"))) {
        group2$data.groups <- hzar.copyModelLabels(group1, group2$data.groups)
        return(group2)
    }
    if (!is.list(group1) || !is.list(group2)) 
        stop("group1 or group2 of malformed type")
    names(group2) <- as.character(lapply(group2, function(x) names(group1)[which(as.logical(lapply(group1, 
        hzar.sameModel, x)))[1]]))
    return(group2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
