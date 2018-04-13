# AAfun
ASReml-R Added functions

# Attension: AAfun is now available and works for ASReml-R V3.0 or ASReml-R V4.0. The online version AAfun is only for ASReml-R V3.0. The new version AAfun4 works for ASReml-R V4.0, but not supplied online in present. If readers feel interesting in AAfun4, you can get the package AAfun4 by sending email to me (yzhlinscau@163.com).

## INSTALL
source("http://bioconductor.org/biocLite.R") 
biocLite(c("GeneticsPed","genetics"))

install.packages(c('amap',"agricolae","agridat","corrgram","grid",'ggplot2',"gplots","MCMCglmm",'nadiv',"plyr","reshape2","sqldf"))

devtools::install_github('yzhlinscau/AAfun')

##DEMO
library(asreml)

library(AAfun)

data(PrSpa)

df<-PrSpa

### function 1 asreml.pin 
##### exmaple 1 for sigle trait model

fm<-asreml(h5~1+Rep, random=~Fam, subset=Spacing=='3',data=df,trace=F)

summary(fm)$varcomp[,1:3]

pin(fm, h2 ~4*V1/(V1+V2),signif=T) # asreml result, h2 formlua

##### exmaple 2 for us model for bi-trait
fm2<-asreml(cbind(dj,h5)~ trait+trait:Rep,
               random=~ us(trait):Fam, rcov=~units:us(trait),
               subset=Spacing=='3',data=df,maxit=40,trace=F)

summary(fm2)$varcomp[,1:3]

pin(fm2, h2_A ~ 4 * V1/(V1+V5)) # heritability for trait A

pin(fm2, h2_B ~ 4 * V3/(V3+V7)) # heritability for trait B

pin(fm2, gCORR ~ V2/sqrt(V1*V3),signif=TRUE) # genetic corr: asreml result, corr formlua, T for signif

pin(fm2, pCORR ~ (V2+V6)/sqrt((V1+V5)*(V3+V7)),signif=TRUE) # phenotype corr

##### exmaple 3 for corr model for multi-trait
fm3<-asreml(cbind(dj,h3,h5)~ trait+trait:Rep, 
            random=~ corgh(trait):Fam, rcov=~units:us(trait),
            subset=Spacing=='3',data=df,maxit=40,trace=F)
            
summary(fm3)$varcomp[,1:3]

pin(fm3,corN=3) #  asreml result, corr coef number

pin(fm3) # corr coef default N is 1, for the first corr
