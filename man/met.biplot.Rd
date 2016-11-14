\name{met.biplot}
\alias{met.biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Biplot asreml-MET results
}
\description{
This function biplots MET factor analytic results from asreml to find the relation of trial sites and the best variety suitable to trial sites. 
}
\usage{
met.biplot(object, siteN, VarietyN, faN,dSco.u = NULL,dLam.u = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Asreml factor analytic results for MET
}
  \item{siteN}{
Total trial site numbers
}
  \item{VarietyN}{
Total variety numbers
}
  \item{faN}{
Total factor numbers
}
  \item{dSco.u}{
Least score of Variety breeding value
}
  \item{dLam.u}{
Least distance from center
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\author{
Yuanzhen Lin <yzhlinscau@163.com>
}
\references{
Yuanzhen Lin. R & ASReml-R Statistics. China Forestry Publishing House. 2016
}
\seealso{
 \code{\link{met.plot}} 
 \code{\link{met.corr}}
}

\examples{
library(asreml)
library(AAfun)
data(MET)

MET$yield<-0.01*MET$yield
met.asr<-asreml(yield~Loc, random=~ Genotype:fa(Loc,2), 
              rcov=~ at(Loc):ar1(Col):ar1(Row), 
              data=MET, maxiter=50)
              
met.biplot(met.asr,6,36,2)
met.biplot(met.asr,siteN=6,VarietyN=36,faN=2,dSco.u=1.8,dLam.u=1.0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ asreml.MET }
\keyword{ met.biplot }% __ONLY ONE__ keyword per line
