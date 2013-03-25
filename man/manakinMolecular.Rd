\name{manakinMolecular}
\alias{manakinMolecular}
\docType{data}
\title{
  Molecular data samples for multiple loci accross the Manakin Cline.
}
\description{
  Allele frequencies for multiple loci from localitys sampled in the Manakin Cline.
}
\usage{data(manakinMolecular)}
\format{
  A data frame with 11 observations on the following 30 variables.
  \describe{
    \item{\code{locationID}}{ID code for the locality, a factor
      with levels \code{B} \code{C} \code{D} \code{E} \code{F} \code{G}
      \code{H} \code{I} \code{J} \code{K} \code{L}} 
    \item{\code{distance}}{The distance from locality \code{A}}
    \item{\code{ada.A}}{The frequency of allelle A of locus Ada}
    \item{\code{ada.B}}{The frequency of allelle B of locus Ada}
    \item{\code{ada.nSamples}}{The number of allelles sampled from locus
      Ada}
    \item{\code{ak2.A}}{The frequency of allelle A of locus Ak2}
    \item{\code{ak2.B}}{The frequency of allelle B of locus Ak2}
    \item{\code{ak2.nSamples}}{The number of allelles sampled from locus
      Ak2}
    \item{\code{gsr.A}}{The frequency of allelle A of locus Gsr}
    \item{\code{gsr.B}}{The frequency of allelle B of locus Gsr}
    \item{\code{gsr.C}}{The frequency of allelle C of locus Gsr}
    \item{\code{gsr.D}}{The frequency of allelle D of locus Gsr}
    \item{\code{gsr.E}}{The frequency of allelle E of locus Gsr}
    \item{\code{gsr.nSamples}}{The number of allelles sampled from locus
      Gsr}
    \item{\code{pgm2.A}}{The frequency of allelle A of locus Pgm2}
    \item{\code{pgm2.B}}{The frequency of allelle B of locus Pgm2}
    \item{\code{pgm2.C}}{The frequency of allelle C of locus Pgm2}
    \item{\code{pgm2.D}}{The frequency of allelle D of locus Pgm2}
    \item{\code{pgm2.nSamples}}{The number of allelles sampled from locus
      Pgm2}
    \item{\code{l5.A}}{The frequency of allelle A of locus L5}
    \item{\code{l5.B}}{The frequency of allelle B of locus L5}
    \item{\code{l5.nSamples}}{The number of allelles sampled from locus
      L5}
    \item{\code{pscn3.A}}{The frequency of allelle A of locus Pscn3}
    \item{\code{pscn3.B}}{The frequency of allelle B of locus Pscn3}
    \item{\code{pscn3.nSamples}}{The number of allelles sampled from locus
      Pscn3}
    \item{\code{mtDNA.A}}{The frequency of allelle A of mitochondrial DNA}
    \item{\code{mtDNA.B}}{The frequency of allelle B of mitochondrial DNA}
    \item{\code{mtDNA.nSamples}}{The number of allelles sampled from
      mitochondrial DNA} 
    \item{\code{geneticHybridIndex.mu}}{The mean genetic hybrid index
      value} 
    \item{\code{geneticHybridIndex.sigma}}{The standard deviation of the
      genetic hybrid index value}
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
data(manakinMolecular)
str(manakinMolecular) ;
%% ##plot(manakinMolecular) ...
}
\keyword{datasets}
