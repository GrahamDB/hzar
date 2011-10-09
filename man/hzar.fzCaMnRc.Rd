\name{hzar.fzCaMnRc}
\alias{hzar.fzCaMnRc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.fzCaMnRc(accepted, candidates, checkFunc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{accepted}{
%%     ~~Describe \code{accepted} here~~
}
  \item{candidates}{
%%     ~~Describe \code{candidates} here~~
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
function (accepted, candidates, checkFunc) 
{
    numCan <- dim(candidates)[[1]]
    if (numCan < 4) 
        return(hzar.fzCline.addOne(accepted, candidate = candidates, 
            checkFunc = checkFunc))
    for (canI in 4 * (1:(numCan\%/\%4)) - 3) {
        if (checkFunc(accepted, candidates[canI, ])) 
            accepted <- rbind(accepted, candidates[canI, ])
        if (checkFunc(accepted, candidates[canI + 1, ])) 
            accepted <- rbind(accepted, candidates[canI + 1, 
                ])
        if (checkFunc(accepted, candidates[canI + 2, ])) 
            accepted <- rbind(accepted, candidates[canI + 2, 
                ])
        if (checkFunc(accepted, candidates[canI + 3, ])) 
            accepted <- rbind(accepted, candidates[canI + 3, 
                ])
    }
    canI = 4 * (numCan\%/\%4) - 3
    if ((!numCan < (canI + 4)) && checkFunc(accepted, candidates[canI + 
        4, ])) 
        accepted <- rbind(accepted, candidates[canI + 4, ])
    if ((!numCan < (canI + 5)) && checkFunc(accepted, candidates[canI + 
        5, ])) 
        accepted <- rbind(accepted, candidates[canI + 5, ])
    if ((!numCan < (canI + 6)) && checkFunc(accepted, candidates[canI + 
        6, ])) 
        accepted <- rbind(accepted, candidates[canI + 6, ])
    return(accepted)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
