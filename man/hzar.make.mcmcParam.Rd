\name{hzar.make.mcmcParam}
\alias{hzar.make.mcmcParam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Set the parameters controling the optimizer environment.

}

\description{
  Defines the general parameters controlling an MCMC process. This
  includes the burnin, chainLength, diangostic output, and the random
  number generator.
  
}
\usage{
hzar.make.mcmcParam(chainLength, burnin, verbosity, thin, seedStreamChannel = 1, useSeedStream = TRUE, mersenneSeed = 12345, lecuyerSeed = rep(12345, 6))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chainLength}{
    Defines how many generations to run the MCMC process after completing burnin. 
  }
  \item{burnin}{
    Defines how many states to generate and discard at the begining of
    the chain. 
  }
  \item{verbosity}{
    If greater than zero, output a diagnostic message every verbosity'th
    generation after completing burnin.
  }
  \item{thin}{ 
    Instruct the MCMC process to return every state that is n * thin
    generations after burnin, where n is a natural number.
  }
  \item{seedStreamChannel}{
    Set the channel number used in the lecuyer random number generator.
  }
  \item{useSeedStream}{
    If TRUE, use the lecuyer random number generator in the MCMC process.
    If FALSE, use the mersenne twister random generator in the MCMC process.
  }
  \item{mersenneSeed}{
    Sets the seed value for mersenne twister.  Expects a numeric of length 1.
  }
  \item{lecuyerSeed}{
    Sets the seed value for lecuyer random number generator.  Expects a
    numeric of length 6. 
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
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
