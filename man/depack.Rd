\name{depack}
\alias{depack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Unload one package with its dependencies
}
\description{
This function will unload one package with all of its dependencies at one time.
}
\usage{
depack(library(AAfun))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{
Expression to unload package.
}

}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{expr }{Expression to unload package, format: library(AAfun)}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Yuanzhen Lin <yzhlinscau@163.com>
}
\references{
Yuanzhen Lin, Xiaoyang Chen. R & ASReml-R Statistical Analysis Tutorial. China Forestry Publishing House. 2014
}
\seealso{
Website for instant update: yzhlin-asreml.ys168.com
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

depack(library(AAfun))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ depack }% __ONLY ONE__ keyword per line
