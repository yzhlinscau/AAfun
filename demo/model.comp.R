##
## demo file for model comparison. 
##

library(AAfun)
df<-PrSpa
library(asreml)

fm<-asreml(h5~ 1+Rep,random=~ Fam, 
           subset=Spacing=='3',data=df,maxit=40,trace=F)

fm1a<-asreml(cbind(dj,h5)~ trait+trait:Rep, 
             random=~ us(trait):Fam, rcov=~units:us(trait),
             subset=Spacing=='3',data=df,maxit=40,trace=F)
fm1b<-asreml(cbind(dj,h5)~ trait+trait:Rep, 
             random=~ diag(trait):Fam, rcov=~units:us(trait),
             subset=Spacing=='3',data=df,maxit=40,trace=F)

fm3a<-asreml(cbind(dj,h3,h5)~ trait+trait:Rep, 
             random=~ diag(trait):Fam, rcov=~units:us(trait),
             subset=Spacing=='3',data=df,maxit=40,trace=F)

fm3b<-asreml(cbind(dj,h3,h5)~ trait+trait:Rep, 
             random=~ diag(trait):Fam, rcov=~units:us(trait),
             subset=Spacing=='3',data=df,maxit=40,trace=F)

#####   model comparison    #####
model.comp(m1=fm1a,m2=fm1b)
model.comp(m1=fm1a,m2=fm1b,LRT=TRUE)
model.comp(m1=fm1a,m2=fm1b,LRT=TRUE,rdDF=TRUE)

model.comp(Nml=c(fm3a,fm3b,fm1a,fm1b,fm),mulM=TRUE)
model.comp(Nml=c(fm3a,fm3b,fm1a,fm1b,fm),mulM=TRUE,LRT=TRUE)


