\name{hzar.cov.rect}
\alias{hzar.cov.rect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a covariance matrix for the cline optimizer.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.cov.rect(clineLLfunc, param.lower, param.upper, pDiv = 11, random = 0, passCenter = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clineLLfunc}{
%%     ~~Describe \code{clineLLfunc} here~~
}
  \item{param.lower}{
%%     ~~Describe \code{param.lower} here~~
}
  \item{param.upper}{
%%     ~~Describe \code{param.upper} here~~
}
  \item{pDiv}{
%%     ~~Describe \code{pDiv} here~~
}
  \item{random}{
%%     ~~Describe \code{random} here~~
}
  \item{passCenter}{
%%     ~~Describe \code{passCenter} here~~
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
function (clineLLfunc, param.lower, param.upper, pDiv = 11, random = 0, 
    passCenter = FALSE) 
{
    if (random > 1e+09) {
        stop("Covariance matrix calculation with random sampling requested with far too many samples.  Stopping.")
    }
    if (random > 0) {
        data.mat <- list(dTheta = prod(abs(as.numeric(param.upper) - 
            as.numeric(param.lower)))/random, data = hzar.gen.rParam.uniform(param.lower, 
            param.upper, random))
    }
    else {
        if ((pDiv^length(param.lower)) > 1e+06) {
            warning("Covariance matrix calculation requested for too complex of a lattice structure. Switching to random sampling using ten thousand samples.")
            return(hzar.cov.rect(clineLLfunc, param.lower, param.upper, 
                random = 10000, passCenter = passCenter))
        }
        data.mat <- hzar.gen.samples.rect(param.lower, param.upper, 
            pDiv)
    }
    param.names <- names(data.mat$data)
    data.wt <- hzar.eval.clineLL(data.mat$data, clineLLfunc)
    data.mat$data <- data.mat$data[data.wt > -1e+06, ]
    data.wt <- data.wt[data.wt > -1e+06]
    MIN.DATA <- (1 + length(param.upper))
    if (length(data.wt) < MIN.DATA) {
        if (random > 0) {
            return(hzar.cov.rect(clineLLfunc, param.lower, param.upper, 
                random = 2 * random, passCenter = passCenter))
        }
        else {
            return(hzar.cov.rect(clineLLfunc, param.lower, param.upper, 
                pDiv = pDiv + 1, passCenter = passCenter))
        }
    }
    while (sum(data.wt > -723) < MIN.DATA) {
        if (sum(data.wt > 609) > 0) {
            if (random > 0) {
                return(hzar.cov.rect(clineLLfunc, param.lower, 
                  param.upper, random = 2 * random, passCenter = passCenter))
            }
            else {
                return(hzar.cov.rect(clineLLfunc, param.lower, 
                  param.upper, pDiv = pDiv + 1, passCenter = passCenter))
            }
        }
        data.wt <- data.wt + 100
    }
    VDATA <- cov.wt(x = cbind(data.mat$data, model.LL = data.wt), 
        wt = exp(data.wt))
    VMATRIX <- VDATA$cov
    mat.scaled <- VMATRIX
    if (passCenter) 
        return(list(cov = mat.scaled[param.names, param.names], 
            center = VDATA$center[param.names]))
    return(mat.scaled[param.names, param.names])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
