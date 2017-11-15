snp.batch <- function (data,SNP.N,traitN,FMod=NULL,maxit=NULL,SNP.detail=FALSE,
                       SNP.figure=FALSE,SNP.signif=TRUE,alpha=0.05){
  options(digits=3)
  
  ## theme for ggplot
  theme.1<-theme(
    axis.title=element_text(face="bold",colour="red",size=20),
    axis.text=element_text(face="bold",colour="black",size=10),
    axis.line=element_line(colour="black",size=1.5),
    panel.grid.major=element_line(colour="grey"),
    panel.grid.minor=element_line(colour="grey",size=.5), #linetype="dashed",
    panel.background=element_rect(fill="grey",colour="blue",size=1.5)    
  )
  
  if(is.null(FMod)) FMod=y~x
  if(is.null(maxit)) maxit=20
  
  aa=SNP.N;cc=traitN
  aaN=length(aa);ccN=length(cc)
  
  mm<-data.frame();tt=SNP.infor=list()
  NTrait=N.Snp=SNP.p=vector() #=NULL
  
  for(i in 1:ccN){
    df1=data[,c(aa,cc[i])]
    nn=length(df1)
    
    NTrait[i]=names(df1)[nn]
    for(j in 1:aaN){
      df2=df1[,c(j,nn)]
      N.Snp[j]=names(df2)[1]
      names(df2)[1:2]=c("x","y")
      
      fm<-asreml(fixed=FMod, maxit=maxit, data=df2,trace=FALSE)
      
      SNP.p[j]=round(wald(fm)[2,4],3)
      
      mm[j,1]=N.Snp[j]
      mm[j,2]=SNP.p[j]
      mm[j,3]=sig.level2(mm[j,2])
      snpinfor=summarise(group_by(df2, x), N=length(y), mean = mean(y), 
                         sd = sd(y), min=min(y),max=max(y))
      SNP.infor[[j]]=list(SNP=N.Snp[j],SNP.infor=snpinfor)
      
      ggfigure=ggplot(df2, aes(x=x, y=y)) + 
        labs(x = N.Snp[j], y = NTrait[i])+
        geom_boxplot(fill="forestgreen")+theme.1
      
      if(SNP.figure==TRUE) print(ggfigure)
    }
    names(mm)=c("SNP","SNP.p","Sig.level")
    
    if(SNP.signif==TRUE) {
      mm=subset(mm,SNP.p<=alpha)
    } else {mm=mm}
    
    if(SNP.detail==FALSE) {tt[[i]]=list(Trait=NTrait[i],SNP.p=mm)} else{
      tt[[i]]=list(Trait=NTrait[i],SNP.p=mm,SNP.infor=SNP.infor)
    }
  }  
  print(tt)   
  cat("=================\n")
  cat("Sig.level: 0'***' 0.001 '**' 0.01 '*' 0.05 '.' 0.10 'Not signif' 1\n\n")
}
