\name{manakinMorphological}
\alias{manakinMorphological}
\docType{data}
\title{
   Manakin observed Morphological Traits
}
\description{
   Morphological Traits observations of individuals sampled from the
   Manakin Cline.
}
\usage{data(manakinMorphological)}
\format{
  A data frame with 165 observations on the following 7 variables.
  \describe{
    \item{\code{Population}}{The id code of the population sample site,
      a factor with levels \code{A} \code{B} \code{C} \code{D} \code{E}
      \code{F} \code{G} \code{H} \code{I} \code{J} \code{K} \code{L}} 
    \item{\code{ID}}{The ID code of the individual sampled}
    \item{\code{Name}}{The population site name, a factor with levels
      \code{Chiriqui_Grande} \code{Costa_Rica} \code{Quebrada_Pastores}
      \code{Rio_Changuinola} \code{Rio_Oeste} \code{Rio_Robalo}
      \code{Rio_Sixaola} \code{Rio_Teribe} \code{Rio_Uyama}
      \code{Soberania} \code{Tierra_Oscura} \code{Valiente_Peninsula}} 
    \item{\code{collar.color}}{The collar color of the individual sampled}
    \item{\code{belly.color}}{The belly color of the individual sampled}
    \item{\code{epaulet.width}}{The width of the epaulet of the individual sampled}
    \item{\code{beard.length}}{The length of the beard of the individual sampled}
  }
}
% \details{
% %%  ~~ If necessary, more details than the __description__ above ~~
% }
\source{
  Brumfield, R. T., R. W. Jernigan, D. B. McDonald, and
  M. J. Braun. 2001. Evolutionary implications of divergent clines in an
  avian (Manacus: Aves) hybrid zone. Evolution 55:2070-2087. 

}
% \references{
% %%  ~~ possibly secondary sources and usages ~~
% }
\examples{
data(manakinMorphological)
str(manakinMorphological)
}
\keyword{datasets}
