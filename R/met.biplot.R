met.biplot <-
  function(object,siteN,VarietyN,faN,dSco.u=NULL,dLam.u=NULL,...){
    # faRS=1
    par(mar=c(4,2,4.5,2))  
    n<-siteN*(1+faN)
    #arr<-summary(object)$varcomp$gamma[1:n] # fa loading, site n., 1(Psi) + fa n.
    gamma<-summary(object)$varcomp[1] 
    rownm<-row.names(gamma)
    aimnm<-rownm[grep('fa\\(.*,.*\\)',rownm)]
    arr<-gamma[aimnm,]
    Xfam<-matrix(arr,siteN,(1+faN))
    
    fa.name<-paste("FA",1:faN,sep="")
    dimnames(Xfam)<-list(paste("S",1:siteN,sep=""),c("Psi",fa.name))
    #windows(8,8)
    pairs(Xfam,main="Fig 1 pairs of Psi with FAs")
    
    ss<-svd(Xfam[,-1])
    Lam<-Xfam[,-1] %*% ss$v
    colnames(Lam)<-c(paste("FA",1:faN,sep="")) # c("Fa1","Fa2")
    Gvar<-Lam %*% t(Lam)+diag(Xfam[,1])
    cLam<-diag(1/sqrt(diag(Gvar))) %*% Lam
    varp<-round(mean(diag(Lam %*% t(Lam))/diag(Gvar))*100,2) # %variance explained
    cat("\nFA number is: ",faN,",\t\t%Var.explained is:\t",varp,"\n")
    
    if(faN==1){cat("\nAttension: biplot worked when more than 2 FAs!\n\n")}
    
    if(faN>1){
      
      bv<-coef(object)$random # here would be complexed!!
      alln<-row.names(bv)
      aimn<-alln[grep('fa\\(.*,.*\\)',alln)]
      Xfasln<-bv[aimn,]
      
      Xfa2<-matrix(Xfasln,nrow=VarietyN) # effects, row=Variety n., col=site n. + fa n.
      colN<-length(Xfasln)/VarietyN
      scores<-Xfa2[,-1:-(colN-faN)]
      dimnames(scores)<-list(paste("V",1:VarietyN,sep=""),paste("Fa",1:faN,sep=""))
      
      acb<-combn(1:faN,2)
      bl<-faN*(faN-1)/2
      
      mLam<-rep(1/siteN,siteN) %*% Lam # get loading means
      sLam<-Lam-rep(mLam,rep(siteN,faN)) # center loadings 
      dLam<-sqrt((sLam*sLam) %*% rep(1,faN)) # distance from center 
      dSco<-sqrt((scores*scores) %*% rep(1,faN)) 
      
      dSco.a<-0.65*max(dSco,rm=T)
      dLam.a<-0.7*max(dLam,rm=T)
      if(is.null(dSco.u)) dSco.u=round(dSco.a,1) # 2
      if(is.null(dLam.u)) dLam.u=round(dLam.a,1) #0.1
      
      if(faN>2){
        for(i in 1:bl){
          #windows(18,8)
          par(mfrow=c(1,2))
          biplot(scores[,acb[,i]],Lam[,acb[,i]],cex=0.75,
                 main=paste("Fig 2-",i, " biplot with all genotypes",sep=""))
          
          biplot(scores[dSco>dSco.u,acb[,i]],Lam[dLam>dLam.u,acb[,i]],cex=0.75,           
                 main=paste("Fig 3-",i, " biplot when dSco>",dSco.u,sep="")) # dSco>2
          abline(h=0,lty=3)
          abline(v=0,lty=3)        
        }
      }else {
        #windows(18,8)
        par(mfrow=c(1,2))
        biplot(scores[,1:2],Lam[,1:2],cex=0.75,
               main="Fig 2 biplot with all genotypes")
        
        biplot(scores[dSco>dSco.u,1:2],Lam[dLam>dLam.u,1:2],cex=0.75,           
               main=paste("Fig 3 biplot when dSco>",dSco.u,sep="")) # dSco>2
        abline(h=0,lty=3)
        abline(v=0,lty=3)
      }
      
      dscores<-data.frame(scores[dSco>dSco.u,],Scores=dSco[dSco>dSco.u]) #2
      ddLam<-data.frame(Lam[dLam>dLam.u,],distFC=dLam[dLam>dLam.u]) # 0.1
      
      cat("\nScores.u is:",dSco.u,"\n")
      print(round(dscores,3)) 
      cat("\ndistFC.u is:",dLam.u,"\n")
      print(round(ddLam,3))
      cat("\n")
    }
}
