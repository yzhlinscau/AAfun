# AAfun
ASReml-R Added functions

# Attension: AAfun is now available and works for ASReml-R V3.0 or ASReml-R V4.0. The online version AAfun is only for ASReml-R V3.0. The new version AAfun4 works for ASReml-R V4.0, but not supplied online in present. If readers feel interesting in AAfun4, you can get the package AAfun4 by sending email to me (yzhlinscau@163.com).

## INSTALL
``` r
install.packages(c('amap',"agricolae","agridat","grid",'ggplot2',"gplots",
"devtools","MCMCglmm",'nadiv',"plyr","reshape2","sqldf"))

devtools::install_github('yzhlinscau/AAfun')
``` 
## DEMO functions
``` r
library(asreml)
library(AAfun)
```
## DEMO data
#### data1 with no pedigree
``` r
data(PrSpa)
df<-PrSpa
```
#### data2 with pedigree
``` r
data(dfm2)
df2<-dfm2

``` 
## function 1 pin():calculate se for h2 or corr
##### exmaple 1.1 for sigle trait model
``` r
fm<-asreml(h5~1+Rep, random=~Fam, subset=Spacing=='3',data=df,trace=F)

summary(fm)$varcomp[,1:3]
``` 
calculate heritability:
``` r
pin(fm, h2 ~4*V1/(V1+V2),signif=T) 
``` 
##### exmaple 1.2 for us model for bi-trait
``` r
fm2<-asreml(cbind(dj,h5)~ trait+trait:Rep,
               random=~ us(trait):Fam, rcov=~units:us(trait),
               subset=Spacing=='3',data=df,maxit=40,trace=F)

summary(fm2)$varcomp[,1:3]
``` 
calculate heritability for both traits:
``` r
pin(fm2, h2_A ~ 4 * V1/(V1+V5)) 

pin(fm2, h2_B ~ 4 * V3/(V3+V7)) 
``` 
calculate genetic and phenotypic corr between traits:
``` r
pin(fm2, gCORR ~ V2/sqrt(V1*V3),signif=TRUE) 

pin(fm2, pCORR ~ (V2+V6)/sqrt((V1+V5)*(V3+V7)),signif=TRUE) 
``` 
##### exmaple 1.3 for corr model for multi-trait
``` r
fm3<-asreml(cbind(dj,h3,h5)~ trait+trait:Rep, 
            random=~ corgh(trait):Fam, rcov=~units:us(trait),
            subset=Spacing=='3',data=df,maxit=40,trace=F)
            
summary(fm3)$varcomp[,1:3]
``` 
return corr results:
``` r
pin(fm3,corN=3) 
```
just return the first corr:
``` r
pin(fm3) 
``` 
pin() also works for data with pedigree files.

## function 2 asreml.batch():batch analysis
#### exmaple 2.1 for sigle trait model
``` r
df1=subset(df,Spacing==3)
asreml.batch(data=df1,factorN=1:5,traitN=c(9:13),
             FMod=y~1+Rep+Plot,RMod=~Fam,
             pformula=h2 ~ 4 * V1/(V1+V2))
```
#### exmaple 2.2 for us model
``` r
asreml.batch(data=df1,factorN=1:5,traitN=c(10:13),
             FMod=cbind(y1,y2)~trait+trait:Rep,
             RMod=~us(trait):Fam,
             EMod=~units:us(trait),
             mulT=TRUE,mulN=2,mulR=TRUE,corMout=T,
             pformula=r.g ~ V2/sqrt(V1*V3),
             pformula1=h2.A ~ 4*V1/(V1+V5),
             pformula2=h2.B ~ 4*V3/sqrt(V3+V7))
```

#### exmaple 2.3 for corr model
``` r
asreml.batch(data=df1,factorN=1:5,traitN=c(10:13),
             FMod=cbind(y1,y2,y3)~trait+trait:Rep,
             RMod=~corgh(trait):Fam,
             EMod=~units:us(trait), maxit=30,
             mulT=TRUE,mulN=3,mulR=TRUE,corM=TRUE)
```
### example II for dataset with pedigree
``` r
ped<-df2[,1:3]  
pedinv<-asreml.Ainverse(ped)$ginv
df2a=subset(df2,Spacing==3)
```
#### example 2.4 sigle trait model
``` r
asreml.batch(data=df2a,factorN=1:6,traitN=c(7:14),
             FMod=y~1+Rep,RMod=~ped(TreeID),
             ped=T,pedinv=pedinv,                                   
             ginverse=list(TreeID=pedinv),
             pformula=h2 ~ V1/(V1+V2))
```

#### exmaple 2.5 us model 
``` r
asreml.batch(data=df2a,factorN=1:6,traitN=c(10:14),
                    FMod=cbind(y1,y2)~trait+trait:Rep,
                    RMod=~us(trait):ped(TreeID),
                    EMod=~units:us(trait),maxit=40,
                    mulT=TRUE,mulN=2,mulR=TRUE,corMout=F,
                    ped=T,pedinv=pedinv,
                    ginverse=list(TreeID=pedinv),
                    pformula=r.g ~ V2/sqrt(V1*V3),
                    pformula1=h2.A ~ V1/(V1+V5),
                    pformula2=h2.B ~ V3/(V3+V7))
```

#### exmaple 2.6  corr model
``` r
asreml.batch(data=df2a,factorN=1:6,traitN=c(10:14),
             FMod=cbind(y1,y2)~trait+trait:Rep,
             RMod=~corgh(trait):ped(TreeID),
             EMod=~units:us(trait),maxit=40,
             mulT=TRUE,mulN=2,corM=TRUE,
             ped=T,pedinv=pedinv,
             ginverse=list(TreeID=pedinv)) 
```
## function 3 model.comp(): model comparison
``` r
fm1<-asreml(h5~ 1+Rep,random=~ Fam, 
            subset=Spacing=='3',data=df,maxit=40)
            
fm1a<-asreml(cbind(dj,h5)~ trait+trait:Rep, 
            random=~ us(trait):Fam, rcov=~units:us(trait),
            subset=Spacing=='3',data=df,maxit=40)
fm1b<-asreml(cbind(dj,h5)~ trait+trait:Rep, 
            random=~ diag(trait):Fam, rcov=~units:us(trait),
            subset=Spacing=='3',data=df,maxit=40)

fm3a<-asreml(cbind(dj,h3,h5)~ trait+trait:Rep, 
            random=~ diag(trait):Fam, rcov=~units:us(trait),
            subset=Spacing=='3',data=df,maxit=40)
            
fm3b<-asreml(cbind(dj,h3,h5)~ trait+trait:Rep, 
            random=~ diag(trait):Fam, rcov=~units:us(trait),
            subset=Spacing=='3',data=df,maxit=40)
```
comparison between two models:
``` r
model.comp(m1=fm1a,m2=fm1b)
model.comp(m1=fm1a,m2=fm1b,LRT=TRUE)
model.comp(m1=fm1a,m2=fm1b,LRT=TRUE,rdDF=TRUE)
```
comparison among more than two models:
``` r
model.comp(Nml=c(fm3a,fm3b,fm1a,fm1b,fm1),mulM=TRUE)
model.comp(Nml=c(fm3a,fm3b,fm1a,fm1b,fm1),mulM=TRUE,LRT=TRUE)
```
### More details will be updated in the future.
