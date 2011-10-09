\name{hzar.meta.CLTnA}
\alias{hzar.meta.CLTnA}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(hzar.meta.CLTnA)}
\format{
  The format is:
List of 4
 $ prior         :function (center, width, xMin, xMax)  
 $ func          :function (center, width, xMin, xMax)  
 $ req           :function (center, width, xMin, xMax)  
 $ parameterTypes:List of 4
  ..$ center:List of 2
  .. ..$ val: num 10
  .. ..$ w  : num 1.5
  .. ..- attr(*, "class")= chr "clineParameter"
  .. ..- attr(*, "param")= chr "center"
  .. ..- attr(*, "fixed")= logi FALSE
  .. ..- attr(*, "limit.lower")= num -1e+08
  .. ..- attr(*, "limit.upper")= num 1e+08
  .. ..- attr(*, "realBTWN01")= logi FALSE
  ..$ width :List of 2
  .. ..$ val: num 10
  .. ..$ w  : num 1.5
  .. ..- attr(*, "class")= chr "clineParameter"
  .. ..- attr(*, "param")= chr "width"
  .. ..- attr(*, "fixed")= logi FALSE
  .. ..- attr(*, "limit.lower")= num 1e-08
  .. ..- attr(*, "limit.upper")= num 1e+08
  .. ..- attr(*, "realBTWN01")= logi FALSE
  ..$ xMin  :List of 2
  .. ..$ val: num 0
  .. ..$ w  : num 1.5
  .. ..- attr(*, "class")= chr "clineParameter"
  .. ..- attr(*, "param")= chr "xMin"
  .. ..- attr(*, "fixed")= logi FALSE
  .. ..- attr(*, "limit.lower")= num -1e+08
  .. ..- attr(*, "limit.upper")= num 1e+08
  .. ..- attr(*, "realBTWN01")= logi FALSE
  ..$ xMax  :List of 2
  .. ..$ val: num 10
  .. ..$ w  : num 1.5
  .. ..- attr(*, "class")= chr "clineParameter"
  .. ..- attr(*, "param")= chr "xMax"
  .. ..- attr(*, "fixed")= logi FALSE
  .. ..- attr(*, "limit.lower")= num -1e+08
  .. ..- attr(*, "limit.upper")= num 1e+08
  .. ..- attr(*, "realBTWN01")= logi FALSE
 - attr(*, "class")= chr "clineMetaModel"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(hzar.meta.CLTnA)
## maybe str(hzar.meta.CLTnA) ; plot(hzar.meta.CLTnA) ...
}
\keyword{datasets}
