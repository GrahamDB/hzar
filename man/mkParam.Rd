\name{mkParam}
\alias{mkParam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mkParam(name, value, weight, lower, upper, isFixed = FALSE, isRB01 = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
%%     ~~Describe \code{name} here~~
}
  \item{value}{
%%     ~~Describe \code{value} here~~
}
  \item{weight}{
%%     ~~Describe \code{weight} here~~
}
  \item{lower}{
%%     ~~Describe \code{lower} here~~
}
  \item{upper}{
%%     ~~Describe \code{upper} here~~
}
  \item{isFixed}{
%%     ~~Describe \code{isFixed} here~~
}
  \item{isRB01}{
%%     ~~Describe \code{isRB01} here~~
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
function (name, value, weight, lower, upper, isFixed = FALSE, 
    isRB01 = FALSE) 
{
    objCParamC <- list(val = value, w = weight)
    class(objCParamC) <- "clineParameter"
    attr(objCParamC, "param") <- name
    attr(objCParamC, "fixed") <- isFixed
    attr(objCParamC, "limit.lower") <- lower
    attr(objCParamC, "limit.upper") <- upper
    attr(objCParamC, "realBTWN01") <- isRB01
    return(objCParamC)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
