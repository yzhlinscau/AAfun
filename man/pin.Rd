\name{pin}
\alias{pin}
\alias{sig.level}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Count error for h2 and corr
}
\description{
This function counts standard error(se) for heritability(h2) and corr value and also outputs significent level for corr value in ASReml-R package.
}
\usage{
pin(object, formula=NULL, signif=NULL, corN = NULL, Rdf=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
asreml results
}
  \item{formula}{
formula for h2 or corr
}
  \item{signif}{
Index to output signif levels, F(default) for non-signif.
}
  \item{corN}{
Number of corr, 1(default).
}
  \item{Rdf}{
Index to output results to vector, F(default) for non-vector output.
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
data(PrSpa)
df<-PrSpa

# exmaple 1 for sigle trait model
fm<-asreml(h5~1+Rep, random=~Fam, subset=Spacing=='3',data=df)
summary(fm)$varcomp[,1:3]
pin(fm, h2 ~4*V1/(V1+V2))
pin(fm, h2 ~4*V1/(V1+V2),Rdf=TRUE)

# exmaple 2 for us model
fm2<-asreml(cbind(dj,h5)~ trait+trait:Rep, 
              random=~ us(trait):Fam, rcov=~units:us(trait),
              subset=Spacing=='3',data=df,maxit=40)

summary(fm2)$varcomp[,1:3]

pin(fm2, h2_A ~ 4 * V1/(V1+V5)) # heritability for trait A
pin(fm2, h2_B ~ 4 * V3/(V3+V7)) # heritability for trait B

pin(fm2, gCORR ~ V2/sqrt(V1*V3),signif=T) # genetic corr
pin(fm2, pCORR ~ (V2+V6)/sqrt((V1+V5)*(V3+V7)),signif=T) # phenotype corr

# exmaple 3 for corr model
fm3<-asreml(cbind(dj,h3,h5)~ trait+trait:Rep, 
              random=~ corgh(trait):Fam, rcov=~units:us(trait),
              subset=Spacing=='3',data=df,maxit=40)

summary(fm3)$varcomp[,1:3]
pin(fm3,corN=3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ asreml.pin }
\keyword{ corr.signif }% __ONLY ONE__ keyword per line
