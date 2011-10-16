\name{hzar.sameModel}
\alias{hzar.sameModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Test hzar objects for identical associated cline models.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.sameModel(fitA, fitB)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitA}{
%%     ~~Describe \code{fitA} here~~
}
  \item{fitB}{
%%     ~~Describe \code{fitB} here~~
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
function (fitA, fitB) 
{
    if (!inherits(fitA, c("hzar.fitRequest", "hzar.dataGroup")) || 
        !inherits(fitB, c("hzar.fitRequest", "hzar.dataGroup"))) {
        stop("Can only compare hzar.fitRequest or hzar.dataGroup objects.")
    }
    return(isTRUE(all.equal(hzar.extract.modelPrep(fitA), hzar.extract.modelPrep(fitB))))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
