\name{manakinMorphological}
\alias{manakinMorphological}
\docType{data}
\title{
  Manakin observed Morphological Traits
%%   ~~ data name/kind ... ~~
}
\description{
  Manakin observed Morphological Traits
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(manakinMorphological)}
\format{
  A data frame with 213 observations on the following 14 variables.
  \describe{
   
    \item{\code{PopulationID}}{a factor}
    \item{\code{Location}}{a factor with levels \code{Chiriqui Grande} \code{Costa Rica} \code{Quebrada Pastores} \code{Rio Changuinola} \code{Rio Oeste} \code{Rio Robalo} \code{Rio Sixaola} \code{Rio Teribe} \code{Rio Uyama} \code{Soberania} \code{Tierra Oscura} \code{Valiente Peninsula}}
    \item{\code{TissueID}}{a factor}
    \item{\code{Mass}}{a numeric vector}
    \item{\code{Wing.length}}{a numeric vector}
    \item{\code{Culmen.length}}{a numeric vector}
    \item{\code{Bill.depth}}{a numeric vector}
    \item{\code{Bill.width}}{a numeric vector}
    \item{\code{Tail.length}}{a numeric vector}
    \item{\code{Tarsus.length}}{a numeric vector}
    \item{\code{Epaulet.width}}{a numeric vector}
    \item{\code{Beard.length}}{a numeric vector}
    \item{\code{Collar.color}}{a numeric vector}
    \item{\code{Belly.color}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
  Brumfield, R. T., R. W. Jernigan, D. B. McDonald, and M. J. Braun. 2001. Evolutionary implications of divergent clines in an avian (Manacus: Aves) hybrid zone. Evolution 55:2070-2087.
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(manakinMorphological)
str(manakinMorphological) ;

}
\keyword{datasets}
