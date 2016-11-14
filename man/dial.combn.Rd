\name{dial.combn}
\alias{dial.combn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Recipro and family code of diallel mating.
}
\description{
This function generates the recipro and family code from a diallel mating.
}
\usage{
dial.combn(lbls, N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lbls}{
An english character for genotype
}
  \item{N}{
diallel mating number
}

}
\details{
The input data is a character with a total number for diallel mate.
}
\value{
%%  If it is a LIST, use
\item{lbls }{An english character for genotype,such as "A"}
\item{N }{positive integer}
%% ...
}
\author{
Yuanzhen Lin <yzhlinscau@163.com>
}
\references{
Yuanzhen Lin. R & ASReml-R Statistics. China Forestry Publishing House. 2016
}
\seealso{
Website for instant update: yzhlin-asreml.ys168.com
}

\examples{
library(AAfun)
  
dial.combn("G",8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dial.combn }
\keyword{ AAfun.diallel }% __ONLY ONE__ keyword per line
