##
## demo file for asreml.batch(). 
##

library(asreml)
library(AAfun)

##### example I for dataset without pedigree
data(PrSpa)
df<-PrSpa

# exmaple 1.1 for sigle trait model
df1=subset(df,Spacing==3)
asreml.batch(data=df1,factorN=1:5,traitN=c(9:13),
             FMod=y~1+Rep+Plot,RMod=~Fam,
             pformula=h2 ~ 4 * V1/(V1+V2))

# exmaple 1.2 for us model 
asreml.batch(data=df1,factorN=1:5,traitN=c(10:13),
             FMod=cbind(y1,y2)~trait+trait:Rep,
             RMod=~us(trait):Fam,
             EMod=~units:us(trait),
             mulT=TRUE,mulN=2,mulR=TRUE,corMout=T,
             pformula=r.g ~ V2/sqrt(V1*V3),
             pformula1=h2.A ~ 4*V1/(V1+V5),
             pformula2=h2.B ~ 4*V3/sqrt(V3+V7))


# exmaple 1.3 for corr model
asreml.batch(data=df1,factorN=1:5,traitN=c(10:13),
             FMod=cbind(y1,y2,y3)~trait+trait:Rep,
             RMod=~corgh(trait):Fam,
             EMod=~units:us(trait), maxit=30,
             mulT=TRUE,mulN=3,mulR=TRUE,corM=TRUE)

##### example II for dataset with pedigree
data(dfm2);df=dfm2

ped<-df[,1:3]  
pedinv<-asreml.Ainverse(ped)$ginv
df1=subset(df,Spacing==3)

# example 2.1 sigle trait model
asreml.batch(data=df1,factorN=1:6,traitN=c(7:14),
             FMod=y~1+Rep,RMod=~ped(TreeID),
             ped=T,pedinv=pedinv,                                   
             ginverse=list(TreeID=pedinv),
             pformula=h2 ~ V1/(V1+V2))


# exmaple 2.2 us model 
asreml.batch(data=df1,factorN=1:6,traitN=c(10:14),
             FMod=cbind(y1,y2)~trait+trait:Rep,
             RMod=~us(trait):ped(TreeID),
             EMod=~units:us(trait),maxit=40,
             mulT=TRUE,mulN=2,mulR=TRUE,corMout=F,
             ped=T,pedinv=pedinv,
             ginverse=list(TreeID=pedinv),
             pformula=r.g ~ V2/sqrt(V1*V3),
             pformula1=h2.A ~ V1/(V1+V5),
             pformula2=h2.B ~ V3/(V3+V7))


# exmaple 2.3  corr model
asreml.batch(data=df1,factorN=1:6,traitN=c(10:14),
             FMod=cbind(y1,y2)~trait+trait:Rep,
             RMod=~corgh(trait):ped(TreeID),
             EMod=~units:us(trait),maxit=40,
             mulT=TRUE,mulN=2,corM=TRUE,
             ped=T,pedinv=pedinv,
             ginverse=list(TreeID=pedinv))

