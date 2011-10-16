\name{hzar.AIC.hzar.dataGroup}
\alias{hzar.AIC.hzar.dataGroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the AIC score for the given hzar.dataGroup object.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.AIC.hzar.dataGroup(dataGroup)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataGroup}{
%%     ~~Describe \code{dataGroup} here~~
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
function (dataGroup) 
{
    return(hzar.AIC.hzar.cline(hzar.get.ML.cline(dataGroup)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
