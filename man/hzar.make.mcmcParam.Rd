\name{hzar.make.mcmcParam}
\alias{hzar.make.mcmcParam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Set the parameters controling the optimizer environment.
%%  ~~function to do ... ~~
}
\description{
  This method assists in assigning appropriate values for controller \code{\link{MCMCmetrop1R}} optimizer.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.make.mcmcParam(chainLength, burnin, verbosity, thin, seedStreamChannel = 1, useSeedStream = TRUE, mersenneSeed = 12345, lecuyerSeed = rep(12345, 6))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chainLength}{
    How many generations to run the mcmc process after burn-in.
    %%     ~~Describe \code{chainLength} here~~
  }
  \item{burnin}{
    How many initial generations are to be considered burn-in and discarded.
    %%     ~~Describe \code{burnin} here~~
  }
  \item{verbosity}{
    If 0, output nothing while mcmc process is running.
    If positive, print the model state every verbosity generations after burn-in.
  }
  \item{thin}{
    Keep only states whose number of generations after burn-in is evenly
    divisible by thin.
  }
  \item{seedStreamChannel}{
    In the lecuyer random number generator, use channel \code{seedStreamChannel}.
  }
  \item{useSeedStream}{
    If TRUE, use the lecuyer random number generator.
    If FALSE, use the mersenne random number generator.
  }
  \item{mersenneSeed}{
    Seed for the mersenne random number generator.
  }
  \item{lecuyerSeed}{
    Seed for the lecuyer random number generator.
  }
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  A list of 5 values:
  \item{chainLength}{The value of the mcmc argument for \code{\link{MCMCmetrop1R}}.}
  \item{burnin}{The value of the burnin argument for \code{\link{MCMCmetrop1R}}.}
  \item{verbosity}{The value of the verbose argument for \code{\link{MCMCmetrop1R}}.}
  \item{thin}{The value of the thin argument for \code{\link{MCMCmetrop1R}}.}
  \item{seed}{The value of the seed argument for \code{\link{MCMCmetrop1R}}.}
  
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
  Graham Derryberry \email{asterion@alum.mit.edu}
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{MCMCmetrop1R}}
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
