\name{hzar.make.cline}
\alias{hzar.make.cline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.make.cline(free.parameters, parameters, func, LL, isValid = is.function(func))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{free.parameters}{
%%     ~~Describe \code{free.parameters} here~~
}
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
}
  \item{func}{
%%     ~~Describe \code{func} here~~
}
  \item{LL}{
%%     ~~Describe \code{LL} here~~
}
  \item{isValid}{
%%     ~~Describe \code{isValid} here~~
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
function (free.parameters, parameters, func, LL, isValid = is.function(func)) 
{
    obj <- list(param.free = free.parameters, param.all = parameters, 
        clineFunc = func, logLike = LL, isValid = isValid)
    class(obj) <- "hzar.cline"
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
