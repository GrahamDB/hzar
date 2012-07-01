\name{hzar.make.clineLLfunc.old.ML}
\alias{hzar.make.clineLLfunc.old.ML}
\alias{hzar.make.clineLLfunc.old.bayes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a Log Likelihood method for the cline model.
%%  ~~function to do ... ~~
}
\description{
  A low-level method for assembling a cline likelihood method.  See
  \code{\link{hzar.first.fitRequest.old.ML}} for an easier to use method
  for working with clines.
  
  Using the arguments given, assemble either a maximum likelihood (ML)
  or bayesian (bayes) function with a single argument theta.  The result
  is suitable for \code{\link{MCMCmetrop1R}}.
}
\usage{
hzar.make.clineLLfunc.old.ML(param.free.names, param.fixed,
  param.check.func, meta.cline.func, model.LL, LLrejectedModel = -1e+08) 
hzar.make.clineLLfunc.old.bayes(param.free.names, param.fixed,
  param.check.func, meta.cline.func, model.LL, prior.LL,
  LLrejectedModel = -1e+08) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param.free.names}{
    The names of the free parameters.
  }
  \item{param.fixed}{
    A named list of the fixed parameters.  
  }
  \item{param.check.func}{
    A boolean function of the parameters, which returns true if the
  model described is valid.
  }
  \item{meta.cline.func}{
    A function of the parameters, which returns a function of distance
  estimating frequency.
  }
  \item{model.LL}{
    A function which takes one argument and returns a log
  likelihood. The argument is a function of distance estimating
  frequency.    
  }
  \item{LLrejectedModel}{
    A finite log likelihood of an invalid model.
  }
  
  \item{prior.LL}{
    A function of the parameters, which returns the log likelihood of
    those parameters.
  }

}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  A function with a single argument theta that returns a finite log
  likelihood.

  Theta is a named list of the free parameters for a cline.

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
  \code{\link{hzar.first.fitRequest.old.ML}}
}
\examples{
##TODO
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
