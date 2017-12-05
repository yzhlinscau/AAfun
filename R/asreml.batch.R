asreml.batch <- function (data,factorN,traitN,FMod=NULL,RMod=NULL, EMod=NULL,
                           mulT=NULL, mulN=NULL,mulR=NULL,corM=NULL,corMout=FALSE,
                           pformula=NULL,pformula1=NULL,pformula2=NULL,
                           pformula3=NULL,pformula4=NULL,sigFormula=NULL,
                           maxit=NULL,#SPmodel=FALSE,
                           ped=NULL,pedinv=NULL,ginverse=NULL) {

  #options(digits=3)
  if(is.null(mulT)) mulT=FALSE
  if(is.null(mulR)) mulR=FALSE
  if(is.null(corM)) corM=FALSE
  if(is.null(maxit)) maxit=20
  if(is.null(ped)) ped=FALSE
  tr=list()
  
  pin2<-function (object, transform)
  {
    pframe <- as.list(object$gammas)
    names(pframe) <- paste("V", seq(1, length(pframe)), sep = "")
    tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)),
                   pframe)
    X <- as.vector(attr(tvalue, "gradient"))
    X[object$gammas.type == 1] <- 0
    tname <- if (length(transform) == 3) transform[[2]] else ""
    n <- length(pframe)
    i <- rep(1:n, 1:n)
    j <- sequence(1:n)
    k <- 1 + (i > j)
    Vmat <- object$ai
    se <- sqrt(sum(Vmat * X[i] * X[j] * k))
    #data.frame(row.names = tname, Estimate = tvalue, SE = se)
    vv=vector() #<-NULL
    vv[1]=tvalue;vv[2]=se
    names(vv)=c(tname,paste(tname,"se",sep="."))
    return(vv)
  }
  
  aa=factorN;cc=traitN
  aaN=length(aa);ccN=length(cc)
  NTrait=pedinv=ginverse2=vector() #=NULL
  
  fm <-list()
  mm1<-data.frame();mm2<-data.frame()
  H2=RS=RS2=H2.se=vector() #=NULL
  Nvar=Nvar1=Nvar2=Nvar3=vector() #=NULL
  H3=H4=H5=H6=vector() #=NULL
  H3.se=H4.se=H5.se=H6.se=vector() #=NULL
  vv2=vv3=vv4=vv5=vv6=vector() #=NULL
  ARV=RV=RV.se=vector() #=NULL
  
  if(mulT==FALSE){
    for(i in 1:ccN){
      df1=data[,c(aa,cc[i])]
      nn=length(df1)
      
      NTrait[i]=names(df1)[nn]
      names(df1)[nn]="y"
      if(is.null(EMod)) EMod=~units
      
      if(ped==FALSE) fm<-asreml(fixed=FMod, random=RMod,rcov=EMod,maxit=maxit, data=df1,trace=FALSE)
      if(ped==TRUE) {
        #fm<-asreml(fixed=FMod, random=RMod,maxit=maxit, ginverse=ginverse2, data=df1,trace=F)
        fm<- do.call(asreml,list(fixed=FMod, random=RMod,rcov=EMod,maxit=maxit,ginverse=ginverse, data=quote(df1),trace=FALSE))
      }
      
      Var=summary(fm)$varcomp
      Nvar=row.names(Var)
      Nvar1=strsplit(Nvar,"!")
      nn=length(Nvar1)
      #if(SPmodel==TRUE) {
      #  for(x in 1:(nn-2)) Nvar2[x]=Nvar1[[x]][1]
      #  for(x in (nn-1):nn) Nvar2[x]=Nvar1[[x]][2]
      #}else {
        for(x in 1:nn) Nvar2[x]=Nvar1[[x]][1]
        #}  ## random factors
      
      ff=row.names(wald(fm)) 
      nf=length(ff)
      ffa=ff[c(-1,-nf)] ## fixed factors
      
      vvN=nrow(Var)
      for(jj in 1:vvN) {
        mm1[i,jj]=round(Var[jj,2],4)
        mm2[i,jj]=round(Var[jj,3],4)
      }

      if(!is.null(pformula)){
        vv2=pin2(fm, pformula) #h2 ~ 4 * V1/(V1+V2)
        H2[i]=round(vv2[1],3);H2.se[i]=round(vv2[2],3)
      }
      if(!is.null(pformula1)){
        vv3=pin2(fm, pformula1) 
        H3[i]=round(vv3[1],3);H3.se[i]=round(vv3[2],3)
      }
      if(!is.null(pformula2)){
        vv4=pin2(fm, pformula2) 
        H4[i]=round(vv4[1],3);H4.se[i]=round(vv4[2],3)
      }
      if(!is.null(pformula3)){
        vv5=pin2(fm, pformula3) 
        H5[i]=round(vv5[1],3);H5.se[i]=round(vv5[2],3)
      }
      if(!is.null(pformula4)){
        vv6=pin2(fm, pformula4) 
        H6[i]=round(vv6[1],3);H6.se[i]=round(vv6[2],3)
      }
      RS[i]=fm$converge
      RS2[i]=ncol(fm$monitor)-2
    } 
  }
   
  if(mulT==TRUE){
    if(is.null(mulN)) mulN=2
    
    if((ccN/mulN)>1) {bb=combn(cc,mulN);bbn=ncol(bb)}
    if((ccN/mulN)==1){bb=cc;bbn=1}
    if((ccN/mulN<1)){cat("\nThe trait No is less than in the model!\n");break}
    
    vvN=vector() #=NULL
    
    for(n in 1:bbn){
      if((ccN/mulN)>1) df1=data[,c(aa,bb[,n])] else df1=data[,c(aa,bb)]
      nn=length(df1)
      
      if((ccN/mulN)>1) {NTrait[n]=paste(names(data)[bb[,n]],collapse="-") 
                        names(df1)[(nn-mulN+1):nn]=paste("y",1:mulN,sep="")
      }else {NTrait[n]=paste(names(data)[bb],collapse="-")
                       names(df1)[(nn-mulN+1):nn]=paste("y",1:mulN,sep="")
      }
      
      if(ped==FALSE) fm<-asreml(fixed=FMod, random=RMod,rcov=EMod,maxit=maxit, data=df1,trace=FALSE)
      if(ped==TRUE) {
        fm<- do.call(asreml,list(fixed=FMod, random=RMod,rcov=EMod,maxit=maxit,ginverse=ginverse, data=quote(df1),trace=FALSE))
        }
      
      Var=summary(fm)$varcomp
      
      Nvar=row.names(Var)
      vvN=nrow(Var)
      dd=0.5*mulN*(mulN+1)
      
      Nvar1=strsplit(Nvar,"!")
      for(x in 1:length(Nvar1)) Nvar2[x]=Nvar1[[x]][1]  ## random factors
      
      if(corM==FALSE){
        Nvar3=sub("R!trait","R",Nvar)
        Nvar3=sub("!trait","",Nvar3)
        Nvar3=sub("trait:","",Nvar3)
        Nvar3=sub("R!variance",NA,Nvar3)
        Nvar3=na.omit(Nvar3)
      }else{
        Nvar3=sub("R!trait","R",Nvar)
        Nvar3=sub("!trait","",Nvar3)
        Nvar3=sub(":!trait","",Nvar3)
        Nvar3=sub("trait:","",Nvar3)
        Nvar3=sub("R!variance",NA,Nvar3)
        Nvar3=na.omit(Nvar3)
      }
      
      ff=row.names(wald(fm)) 
      nf=length(ff)
      ffa=ff[c(-1,-nf)] ## fixed factors
      
      for(jj in 1:vvN) {
        mm1[n,jj]=round(Var[jj,2],4) # var
        mm2[n,jj]=round(Var[jj,3],4) # se
      }
      
      ## just only gcorr or add parameter -- sigFormula
      #if(is.null(sigFormula)) ARV=pin2(fm, rg~V2/sqrt(V1 * V3)) else ARV=pin2(fm, sigFormula)
      ARV=pin2(fm, rg~V2/sqrt(V1 * V3))
      RV[n]=round(ARV[1],3);RV.se[n]=round(ARV[2],3)
      
      if(!is.null(pformula)){
        vv2=pin2(fm, pformula) #h2 ~ 4 * V1/(V1+V2)
        H2[n]=round(vv2[1],3);H2.se[n]=round(vv2[2],3)
      }
      if(!is.null(pformula1)){
        vv3=pin2(fm, pformula1) 
        H3[n]=round(vv3[1],3);H3.se[n]=round(vv3[2],3)
      }
      if(!is.null(pformula2)){
        vv4=pin2(fm, pformula2) 
        H4[n]=round(vv4[1],3);H4.se[n]=round(vv4[2],3)
      }
      if(!is.null(pformula3)){
        vv5=pin2(fm, pformula3) 
        H5[n]=round(vv5[1],3);H5.se[n]=round(vv5[2],3)
      }
      if(!is.null(pformula4)){
        vv6=pin2(fm, pformula4) 
        H6[n]=round(vv6[1],3);H6.se[n]=round(vv6[2],3)
      }
      RS[n]=fm$converge
      RS2[n]=ncol(fm$monitor)-2
    }
      
    ffa=sub("trait:","",ffa)
    Nvar2=sub("trait:","",Nvar2)
  }
    
  if(mulT==TRUE) {
    nmm1=ncol(mm1)
    aaa=0.5*(nmm1+1)
    mm1[,aaa]<-NULL 
    mm2[,aaa]<-NULL 
  }
  
  nmm2=ncol(mm1)
  names(mm1)<-paste("V",1:nmm2,sep="")
  names(mm2)<-paste("V",1:nmm2,".se",sep="")
  mm1[,(1+nmm2):(2*nmm2)]=mm2

  nmm=ncol(mm1)  
  if(!is.null(pformula)){
      mm1[,(nmm+1)]=H2;mm1[,(nmm+2)]=H2.se
      names(mm1)[(nmm+1)]=names(vv2[1]);names(mm1)[(nmm+2)]=names(vv2[2])
      if(!is.null(pformula1)){
        mm1[,(nmm+3)]=H3;mm1[,(nmm+4)]=H3.se
        names(mm1)[(nmm+3)]=names(vv3[1]);names(mm1)[(nmm+4)]=names(vv3[2])
        if(!is.null(pformula2)){
          mm1[,(nmm+5)]=H4;mm1[,(nmm+6)]=H4.se
          names(mm1)[(nmm+5)]=names(vv4[1]);names(mm1)[(nmm+6)]=names(vv4[2])
          if(!is.null(pformula3)){
            mm1[,(nmm+7)]=H5;mm1[,(nmm+8)]=H5.se
            names(mm1)[(nmm+7)]=names(vv5[1]);names(mm1)[(nmm+8)]=names(vv5[2])
            if(!is.null(pformula4)){
              mm1[,(nmm+9)]=H6;mm1[,(nmm+10)]=H6.se
              names(mm1)[(nmm+9)]=names(vv6[1]);names(mm1)[(nmm+10)]=names(vv6[2])
            }else{ mm1[,(nmm+9)]=NTrait
            names(mm1)[(nmm+9)]="Trait"
            mm1<-mm1[,c((nmm+9),1:(nmm+8))]}
          }else{ mm1[,(nmm+7)]=NTrait
          names(mm1)[(nmm+7)]="Trait"
          mm1<-mm1[,c((nmm+7),1:(nmm+6))]}
        }else{ mm1[,(nmm+5)]=NTrait
        names(mm1)[(nmm+5)]="Trait"
        mm1<-mm1[,c((nmm+5),1:(nmm+4))]}
      }else{ mm1[,(nmm+3)]=NTrait
       names(mm1)[(nmm+3)]="Trait"
       mm1<-mm1[,c((nmm+3),1:(nmm+2))]}
  }else{
        mm1[,(nmm+1)]=NTrait 
        names(mm1)[(nmm+1)]="Trait"
        mm1<-mm1[,c((nmm+1),1:nmm)]
  }
  
  if(mulR==TRUE&&mulN==2){
    nrr=length(traitN)
    rr<-diag(nrr)
    nn=ncol(mm1) #?
    #rv=mm1[,(nn-1)] #?
    #re=mm1[,nn] #?
    
    rr[lower.tri(rr)]=RV
    rr=t(rr)
    rr[lower.tri(rr)]=RV.se
    row.names(rr)=colnames(rr)=names(data)[traitN]
    
    siglevel<-sig.level(RV,RV.se)
    siglevel=sub("Not signif","",siglevel)
    rr2=rr
    rr2[lower.tri(rr2)]=siglevel
    rr2=as.data.frame(rr2)
  }
  
  cat("\n\nASReml-R batch analysis results:\n")  
  cat("\nFixed Factors --", ffa,"\n")
  cat("Randomed Factors --", unique(Nvar2),"\n") 
  if(!is.null(pformula)){cat("Index formula -- ");print(pformula)}
  if(!is.null(pformula1)){cat("Index formula1 -- ");print(pformula1)}
  if(!is.null(pformula2)){cat("Index formula2 -- ");print(pformula2)}
  if(!is.null(pformula3)){cat("Index formula3 -- ");print(pformula3)}
  if(!is.null(pformula4)){cat("Index formula4 -- ");print(pformula4)}
  if(mulT==TRUE) cat("\nVariance order:",paste(Nvar3,collapse=", ")) else cat("\nVariance order:",paste(unique(Nvar2),collapse=", "))
  cat("\n\n")
  
  nn=ncol(mm1) # var,v.se--nmm
  mm1[,nn+1]=RS
  names(mm1)[nn+1]="Converge"
  mm1[,nn+2]=RS2
  names(mm1)[nn+2]="Maxit"
  
  nn2=ncol(mm1)
  if(is.null(sigFormula)) sigFormula=FALSE
  if(sigFormula==TRUE){
    dsig=mm1[,c(1,(nmm+2):(nn2-2))]
    ndsig=ncol(dsig)
    nsig=(ndsig-1)/2
    nmsig=c("pformula",paste("pformula",nsig-1,sep=""))
    for(i in 1:nsig){
      RV2=dsig[,2*i]
      RV2.se=dsig[,(2*i+1)]
      dsiglevel<-sig.level(RV2,RV2.se)
      dsiglevel=sub("Not signif","",dsiglevel)
      dsig[,(ndsig+i)]=dsiglevel
      #colnames(dsig[(ndsig+i)])=nmsig[i]#names(dsig[2*i])
    }
    ndsig2=ncol(dsig)
    names(dsig[,(ndsig2-nsig+1):ndsig2])=nmsig
    print(dsig)
  }
  
  if(mulR==TRUE&&mulN==2){
    tr=list(Varcomp=mm1, Corr.erro.matrix=rr,Corr.sig.matrix=rr2)
    #cat("\n\nCorr and error matrix:\n")
    return(tr)
    cat("=================\n")
    cat("upper is corr and lower is error (or sig.level) for corr matrix.\n")
    cat("Sig.level: 0'***' 0.001 '**' 0.01 '*' 0.05 'Not signif' 1\n\n")
  }else {return(mm1)}
  
  #if(corMout==TRUE) return(rr)
}
