\name{hzar.fzCline.reduce}
\alias{hzar.fzCline.reduce}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.fzCline.reduce(candidates, checkFunc = hzar.fz.getCheckFunc(names(candidates)), noReCheck = FALSE, glitz = character(0), wedgeCut = 8000 * length(names(candidates)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{candidates}{
%%     ~~Describe \code{candidates} here~~
}
  \item{checkFunc}{
%%     ~~Describe \code{checkFunc} here~~
}
  \item{noReCheck}{
%%     ~~Describe \code{noReCheck} here~~
}
  \item{glitz}{
%%     ~~Describe \code{glitz} here~~
}
  \item{wedgeCut}{
%%     ~~Describe \code{wedgeCut} here~~
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
function (candidates, checkFunc = hzar.fz.getCheckFunc(names(candidates)), 
    noReCheck = FALSE, glitz = character(0), wedgeCut = 8000 * 
        length(names(candidates))) 
{
    candidates <- as.data.frame(candidates)
    numCandidates <- dim(candidates)[[1]]
    if (numCandidates > wedgeCut) {
        a <- hzar.fzCline.reduce(candidates[1:as.integer(numCandidates/2), 
            ], checkFunc = checkFunc, noReCheck = noReCheck, 
            glitz = c(glitz, " "))
        b <- hzar.fzCline.reduce(candidates[as.integer(1 + (numCandidates/2)):numCandidates, 
            ], checkFunc = checkFunc, noReCheck = noReCheck, 
            glitz = c(glitz, "|"))
        numA <- dim(a)[[1]]
        numB <- dim(b)[[1]]
        if ((length(glitz) == 0) || ((rev(glitz)[1]) == " ")) {
            rN <- 1
            rD <- 2
            if (numA < numB) {
                cutV <- (rN * numA)\%/\%(rD)
                res <- (hzar.fzCaMnRc(rbind(b[numB, ], a[cutV:numA, 
                  ]), rbind(a[rev(1:(cutV - 1)), ], b[rev(1:(numB - 
                  1)), ]), checkFunc = checkFunc))
            }
            else {
                cutV <- (rN * numB)\%/\%(rD)
                res <- (hzar.fzCaMnRc(rbind(a[numA, ], b[cutV:numB, 
                  ]), rbind(b[rev(1:(cutV - 1)), ], a[rev(1:(numA - 
                  1)), ]), checkFunc = checkFunc))
            }
            cat("+")
            if (numCandidates > 16 * wedgeCut) 
                cat("\n")
        }
        else {
            if (numA < numB) {
                res <- (hzar.fzCaMnRc(a, b[rev(1:(numB)), ], 
                  checkFunc = checkFunc))
            }
            else {
                res <- (hzar.fzCaMnRc(b, a[rev(1:(numA)), ], 
                  checkFunc = checkFunc))
            }
            cat("-")
        }
        return(res)
    }
    if (numCandidates < 8) 
        return(candidates)
    return(hzar.fzCaMnRc(candidates[1:4, ], candidates[5:numCandidates, 
        ], checkFunc = checkFunc))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
