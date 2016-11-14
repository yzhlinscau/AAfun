\name{mc.se}
\alias{mc.se}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Count error for h2 and corr
}
\description{
This function counts standard error(se) for heritability(h2) and corr value for MCMCglmm package.
}
\usage{
mc.se(object = NULL, Nmc = NULL, confinterval = NULL, lv = NULL, 
        uv = NULL, n = NULL, conf.level = NULL, sigf = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
MCMCglmm model results or h2/corr formula results
}
  \item{Nmc}{
Use MCMCglmm result directly (T) or not (F)
}
  \item{confinterval}{
Confidence interval for heritability or corr
}
  \item{lv}{
Lower value from confidence interval
}
  \item{uv}{
Upper value from confidence interval
}
  \item{n}{
Total number of aim trait observation
}
  \item{conf.level}{
Confidence level
}
  \item{sigf}{
Output significent level(T) or not(F,default)
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Nmc }{Default is T to use MCMCglmm results directly, F for not.}
  \item{n }{Total number of aim trait observation,default value is 1000.}
  \item{conf.level }{Confidence level, default value is 0.95.}
  \item{sigf}{Output significant level (T) or not (default F).}
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
library(AAfun)
data(PrSpa)

df<-subset(PrSpa,Spacing=='3')

######## single trait model ######## 
m1.glmm <- MCMCglmm(h5 ~ 1 + Rep, 
                    random = ~ Fam, pr = TRUE,
                    family = 'gaussian',
                    data = df)

h2.glmm<-4*m1.glmm$VCV[,'Fam']/(m1.glmm$VCV[,'Fam']+m1.glmm$VCV[,'units'])
posterior.mode(h2.glmm)

mc.se(h2.glmm)

confinterval<-HPDinterval(h2.glmm)
mc.se(confinterval=confinterval,Nmc=F)
mc.se(confinterval=confinterval,Nmc=F,n=559,conf.level=0.95)

### Second method has the same result ###
lv<-HPDinterval(h2.glmm)[1]
uv<-HPDinterval(h2.glmm)[2]
mc.se(lv=lv,uv=uv,Nmc=F)

######## bi-trait model ########

df$dj<-1000*df$dj
phen.var<-matrix(c(var(df$dj,na.rm=TRUE),0,0,
                   var(df$h5,na.rm=TRUE)),2,2)

prior<-list(G=list(G1=list(V=phen.var/2,n=2)),
            R=list(V=phen.var/2,n=2))

set.seed(1234)
m2.glmm <- MCMCglmm(cbind(dj,h5)~ trait-1+trait:Rep, random=~ us(trait):Fam, 
                    rcov=~ us(trait):units, data=df, family=c("gaussian", "gaussian"), 
                    nitt=130000,thin=100,burnin=30000,
                    prior=prior,verbose=FALSE,pr=TRUE)

posterior.mode(m2.glmm$VCV)
HPDinterval(m2.glmm$VCV)
#### count se for variance component
mc.se(m2.glmm$VCV)

#### count se for fixed and randomed effects
# mc.se(m2.glmm$Sol)
posterior.mode(m2.glmm$Sol)[c(1:5,40:45,80:85)]
mc.se(m2.glmm$Sol)[c(1:5,40:45,80:85),]

#### count se for heritability
A.h2.glmm<-4*m2.glmm$VCV[,'dj:dj.Fam']/(m2.glmm$VCV[,'dj:dj.Fam']+m2.glmm$VCV[,'dj:dj.units'])
posterior.mode(A.h2.glmm)

mc.se(A.h2.glmm)

confinterval<-HPDinterval(A.h2.glmm)
mc.se(confinterval=confinterval,Nmc=F)
mc.se(confinterval=confinterval,n=559,conf.level=0.95,Nmc=F)

### Second method has the same result for h2 ###
lv<-HPDinterval(A.h2.glmm)[,1]
uv<-HPDinterval(A.h2.glmm)[,2]
mc.se(lv=lv,uv=uv,Nmc=F)

#### count se for corr
gCorr.glmm<-m2.glmm$VCV[,'h5:dj.Fam']/sqrt(m2.glmm$VCV[,'dj:dj.Fam']*m2.glmm$VCV[,'h5:h5.Fam'])
mc.se(gCorr.glmm,sigf=T)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MCMCglmm.se }
\keyword{ mc.se }% __ONLY ONE__ keyword per line
