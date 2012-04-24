\name{hzar.qScores}
\alias{hzar.qScores}
\alias{hzar.qScores.dataGroup}
\alias{hzar.qScores.obsDataGroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate credibility intervals.
}
\description{
  To be written.
}
\usage{
hzar.qScores(x, wt, probs = c(0, 0.25, 0.5, 0.75, 1))
hzar.qScores.dataGroup(dataGroup, probs = c(0.025, 0.5, 0.975))
hzar.qScores.obsDataGroup(oDG, probs = c(0.025, 0.5, 0.975))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    The series of values to analyze.
  }
  \item{wt}{
    The log of the weight applied to each value.
  }
  \item{probs}{
    The frequencies to calculate.
  }
  \item{dataGroup}{
    The \code{\link{hzar.dataGroup}} to analyze.
  }
  \item{oDG}{
    The \code{\link{hzar.obsDataGroup}} to analyze.
  }
}
\details{
  TODO
  
  hzar.qScores.dataGroup generates quantiles for all of the free
  parameters.
  
  hzar.qScores.obsDataGroup generates quantiles for the cline model with
  the best AICc score. 

}
\value{
  TODO

 \item{comp1 }{Description of 'comp1'}
 \item{comp2 }{Description of 'comp2'}
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
  \code{\link{hzar.dataGroup}}
  \code{\link{hzar.obsDataGroup}}
  \code{\link{hzar.AICc.default}}
}
\examples{
##TODO
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
