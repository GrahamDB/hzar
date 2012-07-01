\name{hzar.dBernoulli.LL}
\alias{hzar.dBernoulli.LL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Transformations of Scalar Data into Bernoulli Trials 
}
\description{
  This method operates on a collection of sampled scalar values and the
  sample location site factor.  It calculates a score for each
  "cut" value that can split the samples into unique groups, and
  either returns those scores, the cut value with the best score, or a
  table of the frequencies of the sample values that are less than the
  cut value with the best score.
}
\usage{
hzar.dBernoulli.LL(values, locations, getMax = FALSE, getProbs = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{values}{
    The sample values to use.
  }
  \item{locations}{
    The factor grouping sample values by sample site.
  }
  \item{getMax}{
    Should this method return the best cut value?
  }
  \item{getProbs}{
    Should this method return the table of frequencies of sample values
    which are less than the best cut value?
  }
}
\details{
  The score for a cut value is the information of learning the location
  of a sample conditioned on the knowledge of the whether the sample
  value is greater or less than the cut value.  This score indicates how
  little a cut value distinguishes between populations.

  The best scoring cut value is the one that determines the greatest
  amount of information about the sample location and therefore it is
  the one with lowest score.
}
\value{
  Either a vector of scores, a single cut value, or a table of
  frequencies of the sample values that are less than the cut value with
  the best score. 
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
  \code{\link{hzar.doMorphoSets}}
  \code{\link{hzar.obsData}}
}
% \examples{
% ##TODO
% }
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
