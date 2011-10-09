\name{hzar.fzCline.addOne}
\alias{hzar.fzCline.addOne}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.fzCline.addOne(accepted, candidate, checkFirst = 1, checkFunc = hzar.fz.getCheckFunc(names(accepted)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{accepted}{
%%     ~~Describe \code{accepted} here~~
}
  \item{candidate}{
%%     ~~Describe \code{candidate} here~~
}
  \item{checkFirst}{
%%     ~~Describe \code{checkFirst} here~~
}
  \item{checkFunc}{
%%     ~~Describe \code{checkFunc} here~~
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
function (accepted, candidate, checkFirst = 1, checkFunc = hzar.fz.getCheckFunc(names(accepted))) 
{
    accepted <- as.data.frame(accepted)
    candidate <- as.data.frame(candidate)
    numAccepted <- dim(accepted)[[1]]
    if (dim(candidate)[[1]] > 1) 
        return(hzar.fzCline.addOne(hzar.fzCline.addOne(accepted, 
            candidate[1, ], checkFunc = checkFunc, checkFirst = checkFirst), 
            candidate[2:(dim(candidate)[[1]]), ], checkFunc = checkFunc, 
            checkFirst = (checkFirst + 1)\%\%4))
    if (checkFunc(accepted, candidate)) 
        accepted <- rbind(accepted, candidate)
    numAccepted <- dim(accepted)[[1]]
    return(accepted)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
