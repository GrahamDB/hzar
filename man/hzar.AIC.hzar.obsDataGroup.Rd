\name{hzar.AIC.hzar.obsDataGroup}
\alias{hzar.AIC.hzar.obsDataGroup}
\alias{hzar.AICc.hzar.obsDataGroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate an AIC score table.
}
\description{
  Calculate the AIC or corrected AIC score table for the given
  hzar.obsDataGroup object. There will be one score generated for each
  model associated with this object.
}
\usage{
  hzar.AIC.hzar.obsDataGroup(obsDataGroup, label = "AIC",
    show.count = FALSE, show.param = FALSE)
  hzar.AICc.hzar.obsDataGroup(obsDataGroup, label = "AICc",
    show.count = FALSE, show.param = FALSE) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obsDataGroup}{
    The hzar.obsDataGroup object for which to generate the score table.
  }
  \item{label}{
    The name to use for the score column. 
  }
  \item{show.count}{
    Include an addition column with a count of the number of free parameters.
  }
  \item{show.param}{
    Currently does nothing.
    
    Include additional columns with the parameter values of the maximum
    likelhood model for each score.
  }
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  A data frame with at least one column, with the label specified above.

  If the models in the obsDataGroup have names, then the rownames of the
  data frame are the models' names.

  If show.count is TRUE, the result will have one additional column with
  the label "count".  For each score this column will have the number of
  free parameters.

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
  \code{\link{AIC}}
  \code{\link{hzar.AIC.default}}
}
\examples{
##TODO
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
