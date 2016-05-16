\name{snp.batch}
\alias{snp.batch}
\alias{sig.level2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SNP batch analysis
}
\description{
This function carries out batch analysis for SNP effects with trait in ASReml-R package.
}
\usage{
snp.batch(data,SNP.N,traitN,FMod=NULL,maxit=NULL,SNP.detail=FALSE,SNP.figure=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
asreml dataset
}
  \item{SNP.N}{
A vector with sites of all SNPs.
}
  \item{traitN}{
A vector with sites of all traits
}
  \item{FMod}{
Fixed mode, y~x(default).
}
  \item{maxit}{
Maximum number of iterations, 20(default).
}
  \item{SNP.detail}{
Index to output information of every SNP type, F(default) for non-SNP-detail output.
}
  \item{SNP.figure}{
Index to output boxplot of every SNP type, F(default) for non-figure output.
}
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

\examples{

library(asreml)
library(AAfun)
data(SNP)
df<-SNP #str(df)

snp.batch(df,SNP.N=2:4,traitN=c(5:6))
snp.batch(df,SNP.N=2:4,traitN=c(5:6), SNP.detail=T)
snp.batch(df,SNP.N=2:4,traitN=c(5:6), SNP.figure=T)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ snp.batch }% __ONLY ONE__ keyword per line
