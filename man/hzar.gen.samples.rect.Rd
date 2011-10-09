\name{hzar.gen.samples.rect}
\alias{hzar.gen.samples.rect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.gen.samples.rect(param.lower, param.upper, pDiv = 11)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param.lower}{
%%     ~~Describe \code{param.lower} here~~
}
  \item{param.upper}{
%%     ~~Describe \code{param.upper} here~~
}
  \item{pDiv}{
%%     ~~Describe \code{pDiv} here~~
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
function (param.lower, param.upper, pDiv = 11) 
{
    param.names <- names(param.lower)
    nParam <- length(param.names)
    deltas = (as.numeric(param.upper[param.names]) - as.numeric(param.lower[param.names]))/(pDiv - 
        1)
    grid.formals <- matrix(nrow = pDiv, ncol = nParam, data = rep(deltas, 
        each = pDiv)) * rep(0:(pDiv - 1), nParam) + rep(as.numeric(param.lower[param.names]), 
        each = pDiv)
    colnames(grid.formals) <- param.names
    result <- list(dTheta = prod(abs(deltas)))
    result$data <- do.call(expand.grid, as.list(as.data.frame(grid.formals)))
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
