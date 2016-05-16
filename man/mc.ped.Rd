\name{mc.ped}
\alias{mc.ped}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pedigree format transformation
}
\description{
This function changes Pedigree format from asreml to MCMCglmm, when the pedigree file without parents details.
}
\usage{
mc.ped(ped)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{
Pedigree file from asreml format data.
}

}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{ped }{Pedigree file without parents details}
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
data(PrSpa)
df<-PrSpa

### create an asreml-format pedigree
ped<-df[,c(1,4)]
names(ped)[2]<-"Mum"
ped$Dad<-NA

### transform to MCMCglmm-format pedigree
ped2<-mc.ped(ped)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mc.ped }
\keyword{ pedigree }% __ONLY ONE__ keyword per line
