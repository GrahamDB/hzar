\name{hzar.doMorphoSets}
\alias{hzar.doMorphoSets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Make hzar.obsData objects from scalar observations using reference tables 
}
\description{
  Perform a Bernoulli transform on a table of scalar traits of sampled
  individuals while using a separate table of localities.
}
\usage{
hzar.doMorphoSets(traitNames, tDist, tDLocCol, tDDistCol, tValues, tVLocCol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{traitNames}{
    The columns of \code{tValues} to transform.
  }
  \item{tDist}{
    A data.frame with a column of locality identifiers and a
    column of locality distances. See
    \code{\link{manakinLocations}} for an example. 
  }
  \item{tDLocCol}{
    The name of the column of \code{tDist} with locality identifiers.
  }
  \item{tDDistCol}{
    The name of the column of \code{tDist} with locality distances.
  }
  \item{tValues}{
    A data.frame of observed scalar traits of individuals.  See
    \code{\link{manakinMorphological}} for an example.
  }
  \item{tVLocCol}{
    The name of the column of \code{tValues} with locality
    identifiers where the sample was taken.
  }
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  A list of \code{\link{hzar.obsData}} objects. The values of
  \code{traitNames} are used as names.  The Bernoulli likelihood
  function is for each \code{\link{hzar.obsData}} object.
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
  \code{\link{manakinLocations}}
  \code{\link{manakinMorphological}}
  \code{\link{hzar.plot.obsData}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(manakinMorphological);
data(manakinLocations);
mkn <-
  hzar.doMorphoSets("beard.length",
                    tDist=manakinLocations,
                    tDLocCol="LocalityID",
                    tDDistCol="distance",
                    tValues=manakinMorphological,
                    tVLocCol="Locality")
print(str(mkn));
hzar.plot.obsData(mkn$beard.length);
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
