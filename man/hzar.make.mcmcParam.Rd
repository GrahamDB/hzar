\name{hzar.make.mcmcParam}
\alias{hzar.make.mcmcParam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.make.mcmcParam(chainLength, burnin, verbosity, thin, seedStreamChannel = 1, useSeedStream = TRUE, mersenneSeed = 12345, lecuyerSeed = rep(12345, 6))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chainLength}{
%%     ~~Describe \code{chainLength} here~~
}
  \item{burnin}{
%%     ~~Describe \code{burnin} here~~
}
  \item{verbosity}{
%%     ~~Describe \code{verbosity} here~~
}
  \item{thin}{
%%     ~~Describe \code{thin} here~~
}
  \item{seedStreamChannel}{
%%     ~~Describe \code{seedStreamChannel} here~~
}
  \item{useSeedStream}{
%%     ~~Describe \code{useSeedStream} here~~
}
  \item{mersenneSeed}{
%%     ~~Describe \code{mersenneSeed} here~~
}
  \item{lecuyerSeed}{
%%     ~~Describe \code{lecuyerSeed} here~~
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
function (chainLength, burnin, verbosity, thin, seedStreamChannel = 1, 
    useSeedStream = TRUE, mersenneSeed = 12345, lecuyerSeed = rep(12345, 
        6)) 
{
    mcmcParam <- list(chainLength = chainLength, burnin = burnin, 
        verbosity = verbosity, thin = thin)
    mcmcParam$seed <- mersenneSeed
    if (useSeedStream) {
        mcmcParam$seed <- list(lecuyerSeed, seedStreamChannel)
    }
    return(mcmcParam)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
