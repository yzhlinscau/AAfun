\name{MET}
\alias{MET}
\docType{data}
\title{
multi-environment trial dataset
}
\description{
A dataset of one plant in multi-environment trials(MET)
}
\usage{data(MET)}
\format{
  A data frame with 648 observations on the following 9 variables.
  \describe{
    \item{\code{Genotype}}{a factor with levels \code{G1} \code{G10} \code{G11} \code{G12} \code{G13} \code{G14} \code{G15} \code{G16} \code{G17} \code{G18} \code{G19} \code{G2} \code{G20} \code{G21} \code{G22} \code{G23} \code{G24} \code{G25} \code{G26} \code{G27} \code{G28} \code{G29} \code{G3} \code{G30} \code{G31} \code{G32} \code{G33} \code{G34} \code{G35} \code{G36} \code{G4} \code{G5} \code{G6} \code{G7} \code{G8} \code{G9}}
    \item{\code{Loc}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6}}
    \item{\code{Loc2}}{a factor with levels \code{Balcarce} \code{CMuerto} \code{CSuarez} \code{Miramar} \code{Orense} \code{TArroyos}}
    \item{\code{Row}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6} \code{7} \code{8} \code{9}}
    \item{\code{Col}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6} \code{7} \code{8} \code{9} \code{10} \code{11} \code{12} \code{13} \code{14} \code{15} \code{16} \code{17} \code{18}}
    \item{\code{Rep}}{a factor with levels \code{1} \code{2} \code{3}}
    \item{\code{Block}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6}}
    \item{\code{Plot}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6}}
    \item{\code{yield}}{a numeric vector}
  }
}

\references{
Yuanzhen Lin. R & ASReml-R Statistics. China Forestry Publishing House. 2016
}
\seealso{
Website for instant update: yzhlin-asreml.ys168.com
}
\examples{
data(MET)
## maybe str(MET) ; plot(MET) ...
}
\keyword{datasets}
