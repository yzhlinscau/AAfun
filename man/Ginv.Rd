\name{Ginv}
\alias{Ginv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
G inverse matrices output
}
\description{
This function generates G matrices using different methods
}
\usage{
Ginv(marker.file, ped.file, aped.rowNames, path = NULL, Goptions = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{marker.file}{
The marker file, first column is id then SNPs, SNPs must be 0,1,2.
}
  \item{ped.file}{
The pedigree file, contains id, sire(or female) and dam (or male).
}
  \item{aped.rowNames}{
The rowNames of pedigree inverse matrix.
}
  \item{path}{
The path of all datasets for GBLUP.
}
  \item{Goptions}{
G matrices for different methods,1(default) using frequencies,2 using regression.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{marker.file }{marker.file="Genotype.csv"}
  \item{ped.file }{ped.file="pedigree.csv"}
  \item{path}{path="G:\\Users\\yzhlin\\GBLUP"}
  \item{Goptions}{Goptions=1, or 2}
}
\references{
1. VanRaden PM. Efficient methods to compute genomic predictions. J Dairy Sci,2008, 91: 4414-4423.
2. Yuanzhen Lin. R & ASReml-R Statistics. China Forestry Publishing House. 2016
}
\author{
Yuanzhen Lin <yzhlinscau@163.com>
}

\seealso{
Website for instant update: yzhlin-asreml.ys168.com
}

\examples{

library(asreml)
library(AAfun)

#### for demo examples
data(G.data)
data(G.pedigree)
data(G.marker)

gpedinv=asreml.Ainverse(G.pedigree)$ginv
gpnames=attributes(gpedinv)$rowNames

Ginv1=Ginv(marker.file=G.marker,ped.file=G.pedigree,aped.rowNames=gpnames,Goptions=1)

gblup<-asreml(t1~1+Site,random=~ giv(ID), ## ped(ID)
               ginverse=list(ID=Ginv1),data=G.data)

summary(gblup)$varcomp
c<-coef(gblup)$random

#### for user

df<-asreml.read.table("data.csv",T,",")
aped<-asreml.read.table("pedigree.csv",T,",")
apedinv=asreml.Ainverse(aped)$ginv
apnames=attributes(apedinv)$rowNames 

path="G:/Users/yzhlin/Desktop/GBLUP"
Ginv1=Ginv(marker.file="Genotype.csv",ped.file="pedigree.csv",
                       aped.rowNames=apnames,path=path,Goptions=1)

## GBLUP
df.asr<-asreml(t1~1+Site,random=~ giv(ID),
               ginverse=list(ID=Ginv1),data=df)

## ABLUP
df2.asr<-asreml(t1~1+Site,random=~ ped(ID),
               ginverse=list(ID=apedinv),data=df)               

summary(df.asr)$varcomp
summary(df2.asr)$varcomp

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ asreml.Ginv }
\keyword{ GBLUP }% __ONLY ONE__ keyword per line
