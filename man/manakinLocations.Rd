\name{manakinLocations}
\alias{manakinLocations}
\docType{data}
\title{
  Distance from locality A for each locality.
  %%   ~~ data name/kind ... ~~
}
\description{
  Distance from locality A for each locality sampled accross the
  Manakin Cline.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(manakinLocations)}
\format{
  A data frame with 12 observations on the following 2 variables.
  \describe{
    \item{\code{LocalityID}}{a factor with levels \code{A} \code{B} \code{C} \code{D} \code{E} \code{F} \code{G} \code{H} \code{I} \code{J} \code{K} \code{L}}
    \item{\code{distance}}{a numeric vector}
  }
}
% \details{
% %%  ~~ If necessary, more details than the __description__ above ~~
% }
\source{
  Brumfield, R. T., R. W. Jernigan, D. B. McDonald, and M. J. Braun. 2001. Evolutionary implications of divergent clines in an avian (Manacus: Aves) hybrid zone. Evolution 55:2070-2087.
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
% \references{
% %%  ~~ possibly secondary sources and usages ~~
% }
\examples{
data(manakinLocations)
print(manakinLocations)
## maybe str(manakinLocations) ; plot(manakinLocations) ...
}
\keyword{datasets}
