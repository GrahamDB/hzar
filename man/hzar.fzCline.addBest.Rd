\name{hzar.fzCline.addBest}
\alias{hzar.fzCline.addBest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.fzCline.addBest(accepted, candidates, checkFunc)
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
    addFz1 <- hzar.fzCline.addOne
    addFzN <- hzar.fzCline.addMany
    pluckSet <- function(pluckTF) {
        pluckI <- which(pluckTF)
        pluckO <- which(!(pluckTF))
        if (length(pluckI) == 0) 
            return()
        if (length(pluckO) == 0) {
            if (length(pluckI) == 1) {
                accepted <<- rbind(accepted, candidates)
                candidates <<- NULL
                return()
            }
            accepted <<- rbind(accepted, addFz1(candidates[1, 
                ], candidates[2:(dim(candidates)[[1]]), ], checkFunc = checkFunc))
            candidates <<- NULL
            return()
        }
        if (length(pluckI) == 1) {
            accepted <<- rbind(accepted, candidates[pluckI, ])
            candidates <<- candidates[pluckO, ]
            return()
        }
        pkFirst <- pluckI[1]
        pluckI <- pluckI[pluckI != pkFirst]
        accepted <<- rbind(accepted, addFz1(candidates[pkFirst, 
            ], candidates[pluckI, ], checkFunc = checkFunc))
        candidates <<- candidates[pluckO, ]
        return()
    }
    doPluck <- function(param) {
        if (!is.data.frame(candidates)) 
            return()
        if (!all(param \%in\% names(accepted), param \%in\% names(candidates))) 
            return()
        if (max(accepted[, param]) < (pVal <- max(candidates[, 
            param]))) {
            pluckSet(candidates[, param] == pVal)
            if (!is.data.frame(candidates)) 
                return()
        }
        if (min(accepted[, param]) > (pVal <- min(candidates[, 
            param]))) {
            pluckSet(candidates[, param] == pVal)
            if (!is.data.frame(candidates)) 
                return()
        }
        return()
    }
    for (paramI in names(accepted)) {
        doPluck(paramI)
    }
    if (!is.data.frame(candidates)) 
        return(accepted)
    return(addFzN(accepted, candidates, checkFunc = checkFunc))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
