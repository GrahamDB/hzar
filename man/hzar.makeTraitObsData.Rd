\name{hzar.makeTraitObsData}
\alias{hzar.makeTraitObsData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compile and transform raw scalar data
}
\description{
  To be written.
  Compile and transform raw scalar data into a frequency based
  \code{\link{hzar.obsData}} object. 
}
\usage{
hzar.makeTraitObsData(distOfLocation, locationOfValue, values)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distOfLocation}{
    A vector mapping population ids to distances of each population. 
  }
  \item{locationOfValue}{
    A vector of population ids of each sample.
  }
  \item{values}{
    A vector of trait values for each sample. 
  }
}
% \details{
%   TODO
% }
\value{
  A \code{\link{hzar.obsData}} object.
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
  \code{\link{hzar.obsData}}
  \code{\link{hzar.doMorphoSets}}
}
\examples{
##TODO
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
