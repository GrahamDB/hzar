\name{hzar.eval.clineLL}
\alias{hzar.eval.clineLL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.eval.clineLL(data, llFunc, doPar = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{llFunc}{
%%     ~~Describe \code{llFunc} here~~
}
  \item{doPar}{
%%     ~~Describe \code{doPar} here~~
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
function (data, llFunc, doPar = FALSE) 
{
    slices <- hzar.wedgeSlice(dim(data)[[1]])
    useFunc = llFunc
    useData = data
    if (doPar) {
        result <- foreach(tttIndex = slices, .combine = rbind) \%dopar\% 
            {
                data.frame(model.LL = as.numeric(lapply(tttIndex, 
                  function(x) {
                    useFunc(useData[x, ])
                  })))
            }
    }
    else {
        result <- foreach(tttIndex = slices, .combine = rbind) \%do\% 
            {
                data.frame(model.LL = as.numeric(lapply(tttIndex, 
                  function(x) {
                    useFunc(useData[x, ])
                  })))
            }
    }
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
