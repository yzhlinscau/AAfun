\name{met.plot}
\alias{met.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot asreml-MET data
}
\description{
This function plots MET data for further factor analytic by asreml to find the relation of trial sites, etc.
}
\usage{
met.plot(object, plot.title = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
MET data
}
  \item{plot.title}{
MET plot title
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
 \code{\link{met.biplot}} 
 \code{\link{met.corr}}
}

\examples{
library(AAfun)
data(MET)

# example 1
# variable order: genotype,yield,site,row,col
MET2<-MET[,c(1,9,2,4:5)] 
met.plot(MET2)

# example 2
MET3<-MET[,c(1,9,2,4:7)] # add variable order on MET2: Rep, Block
met.plot(MET3,"My met trials")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ asreml.MET.PLOT }
\keyword{ met.plot }% __ONLY ONE__ keyword per line
