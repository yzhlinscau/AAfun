\name{met.corr}
\alias{met.corr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
asreml-MET corr matrix
}
\description{
This function calculate var/cov/corr from asreml MET factor analytic results to further research the relation of trial sites.
}
\usage{
met.corr(object, site, kn = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Asreml factor analytic results for MET
}
  \item{site}{
variable site of data
}
  \item{kn}{
Site cluster group numbers
}

}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{object }{Asreml factor analytic results for MET, met.asr}
  \item{site }{variable site of data, MET$Site}
  \item{kn }{cluster group numbers,3}
%% ...
}
\author{
Yuanzhen Lin <yzhlinscau@163.com>
}
\references{
Yuanzhen Lin, Xiaoyang Chen. R & ASReml-R Statistical Analysis Tutorial. China Forestry Publishing House. 2014
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{met.plot}} 
 \code{\link{met.biplot}}
}
\examples{
library(asreml)
data(MET)
MET$yield<-0.01*MET$yield
met.asr<-asreml(yield~Loc, random=~ Genotype:fa(Loc,2), 
              rcov=~ at(Loc):ar1(Col):ar1(Row), 
              data=MET, maxiter=50)
met.corr(met.asr,MET$Loc,2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ asreml.MET.CORR }
\keyword{ met.corr }% __ONLY ONE__ keyword per line
