\name{manakinMolecular}
\alias{manakinMolecular}
\docType{data}
\title{
  Molecular data samples for multiple loci accross the Manakin Cline.
%%   ~~ data name/kind ... ~~
}
\description{
  Molecular data samples for multiple loci accross the Manakin Cline.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(manakinMolecular)}
\format{
  A data frame with 11 observations on the following 30 variables.
  \describe{
    \item{\code{locationID}}{a factor with levels \code{B} \code{C} \code{D} \code{E} \code{F} \code{G} \code{H} \code{I} \code{J} \code{K} \code{L}}
    \item{\code{distance}}{a numeric vector}
    \item{\code{ada.A}}{a numeric vector}
    \item{\code{ada.B}}{a numeric vector}
    \item{\code{ada.nSamples}}{a numeric vector}
    \item{\code{ak2.A}}{a numeric vector}
    \item{\code{ak2.B}}{a numeric vector}
    \item{\code{ak2.nSamples}}{a numeric vector}
    \item{\code{gsr.A}}{a numeric vector}
    \item{\code{gsr.B}}{a numeric vector}
    \item{\code{gsr.C}}{a numeric vector}
    \item{\code{gsr.D}}{a numeric vector}
    \item{\code{gsr.E}}{a numeric vector}
    \item{\code{gsr.nSamples}}{a numeric vector}
    \item{\code{pgm2.A}}{a numeric vector}
    \item{\code{pgm2.B}}{a numeric vector}
    \item{\code{pgm2.C}}{a numeric vector}
    \item{\code{pgm2.D}}{a numeric vector}
    \item{\code{pgm2.nSamples}}{a numeric vector}
    \item{\code{l5.A}}{a numeric vector}
    \item{\code{l5.B}}{a numeric vector}
    \item{\code{l5.nSamples}}{a numeric vector}
    \item{\code{pscn3.A}}{a numeric vector}
    \item{\code{pscn3.B}}{a numeric vector}
    \item{\code{pscn3.nSamples}}{a numeric vector}
    \item{\code{mtDNA.A}}{a numeric vector}
    \item{\code{mtDNA.B}}{a numeric vector}
    \item{\code{mtDNA.nSamples}}{a numeric vector}
    \item{\code{geneticHybridIndex.mu}}{a numeric vector}
    \item{\code{geneticHybridIndex.sigma}}{a numeric vector}
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
data(manakinMolecular)
str(manakinMolecular) ;
%% ##plot(manakinMolecular) ...
}
\keyword{datasets}
