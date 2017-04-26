\name{makeA1}
\alias{makeA1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates the additive genetic relationship matrix
}
\description{
This function Creates the additive genetic relationship matrix from MCMCglmm's pedigree format, when the pedigree file without parents details, the pedigree should be transformed with mc.ped() at first.
}
\usage{
makeA1(pedigree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{
Pedigree file from MCMCglmm's format data.
}

}

\value{
  \item{ped }{Pedigree file with parents details}
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
data(PrSpa)
df<-PrSpa

### create an asreml-format pedigree
ped<-df[,c(1,4)]
names(ped)[2]<-"Mum"
ped$Dad<-NA
for(i in 1:3) ped[,i]=as.factor(ped[,i])

### transform to MCMCglmm-format pedigree
ped1<-mc.ped(ped)


### Create the additive genetic relationship matrix
makeA1(ped1) # Am=makeA1(ped1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ makeA1 }
\keyword{ pedigree }% __ONLY ONE__ keyword per line
