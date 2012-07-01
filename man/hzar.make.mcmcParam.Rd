\name{hzar.make.mcmcParam}
\alias{hzar.make.mcmcParam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Set the parameters controling the optimizer environment.
}
\description{
  Defines the general parameters controlling an MCMC process. This
  includes the burnin, chainLength, diangostic output, and the random
  number generator seed, as used by the \code{\link{MCMCmetrop1R}} 
  optimizer. 
  
}
\usage{
hzar.make.mcmcParam(chainLength, burnin, verbosity, thin,
  seedStreamChannel = 1, useSeedStream = TRUE, mersenneSeed = 12345,
  lecuyerSeed = rep(12345, 6))
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
    If 0, output nothing while mcmc process is running.
    If positive, print the model state every \code{verbosity} generations after burn-in.
  }
  \item{thin}{
    Keep only states whose number of generations after burn-in is evenly
    divisible by thin.
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
}
\examples{
##TODO
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
