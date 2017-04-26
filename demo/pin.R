##
## demo file for pin() function. 
##

library(asreml)
library(AAfun)
data(PrSpa)
df<-PrSpa

# exmaple 1 for sigle trait model
fm<-asreml(h5~1+Rep, random=~Fam, subset=Spacing=='3',data=df,trace=F)
summary(fm)$varcomp[,1:3]
pin(fm, h2 ~4*V1/(V1+V2))
pin(fm, h2 ~4*V1/(V1+V2),Rdf=TRUE)

# exmaple 2 for us model
fm2<-asreml(cbind(dj,h5)~ trait+trait:Rep, 
            random=~ us(trait):Fam, rcov=~units:us(trait),
            subset=Spacing=='3',data=df,maxit=40,trace=F)

summary(fm2)$varcomp[,1:3]

pin(fm2, h2_A ~ 4 * V1/(V1+V5)) # heritability for trait A
pin(fm2, h2_B ~ 4 * V3/(V3+V7)) # heritability for trait B

pin(fm2, gCORR ~ V2/sqrt(V1*V3),signif=T) # genetic corr
pin(fm2, pCORR ~ (V2+V6)/sqrt((V1+V5)*(V3+V7)),signif=T) # phenotype corr

# exmaple 3 for corr model
fm3<-asreml(cbind(dj,h3,h5)~ trait+trait:Rep, 
            random=~ corgh(trait):Fam, rcov=~units:us(trait),
            subset=Spacing=='3',data=df,maxit=40,trace=F)

summary(fm3)$varcomp[,1:3]
pin(fm3,corN=3)