\name{hzar.doCLTData1DRaw}
\alias{hzar.doCLTData1DRaw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create an hzar.obsData object using a table of individual traits.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.doCLTData1DRaw(distance, traitValue)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distance}{
%%     ~~Describe \code{distance} here~~
}
  \item{traitValue}{
%%     ~~Describe \code{traitValue} here~~
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
function (distance, traitValue) 
{
    if (length(distance) != length(traitValue)) 
        stop("Distance and traitValue vectors not of equal length.")
    dist.group <- unique(distance)
    group.nSamp <- as.numeric(lapply(X = dist.group, function(a) sum(distance == 
        a)))
    if (sum(group.nSamp < 3) > 0) 
        stop("There are not enough samples in discrete populations.  If the data submitted is correct, you should either drop low sample populations or use one of the population sample interpolation methods.")
    if (sum(group.nSamp < 10) > 0) 
        warning("There are very few samples in discrete populations.  You should consider dropping low sample populations or using one of the population sample interpolation methods.")
    group.mean <- as.numeric(lapply(X = dist.group, function(a) mean(traitValue[distance == 
        a])))
    group.var <- as.numeric(lapply(X = dist.group, function(a) var(traitValue[distance == 
        a])))
    if (sum(group.var == 0) > 0) {
        warning("Some of the population samples have a variance of 0.  Adding in estimated variance due to measurement error.")
        m.err <- 5/3 * 10^-quantile(1 + as.numeric(lapply(traitValue, 
            function(nV) {
                min((-1:12)[round(nV, digits = -1:12) == nV])
            })), probs = 0.75)[[1]]
        group.var <- group.var + m.err * m.err
    }
    return(data.frame(dist = dist.group, mu = group.mean, sigma2 = group.var, 
        nSamp = group.nSamp))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
