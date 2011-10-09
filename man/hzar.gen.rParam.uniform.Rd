\name{hzar.gen.rParam.uniform}
\alias{hzar.gen.rParam.uniform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.gen.rParam.uniform(param.lower, param.upper, count = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param.lower}{
%%     ~~Describe \code{param.lower} here~~
}
  \item{param.upper}{
%%     ~~Describe \code{param.upper} here~~
}
  \item{count}{
%%     ~~Describe \code{count} here~~
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
function (param.lower, param.upper, count = 1000) 
{
    raw <- foreach(low = param.lower, high = param.upper, .combine = cbind) \%do\% 
        {
            runif(count, low, high)
        }
    result <- as.data.frame(raw)
    colnames(result) <- names(param.lower)
    names(result) <- names(param.lower)
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
