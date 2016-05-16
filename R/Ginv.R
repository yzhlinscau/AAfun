Ginv <-
function (marker.file, ped.file, aped.rowNames, path=NULL, Goptions=1) {
  #library(genetics)   
  #library(GeneticsPed)  
  #library(asreml)
  #library(MASS)
   
  if(!is.null(path)) setwd(path)
  # marker.file="Genotype.csv";ped.file="pedigree.csv"
  # Reading pedigree file
  if(is.character(ped.file)) pedigree<-asreml.read.table(ped.file,header=TRUE, sep=",") else pedigree=ped.file
  #ainv<-asreml.Ainverse(pedigree)$ginv
  
  # Reading marker file
  if(is.character(marker.file)) genos <-read.table(marker.file, header=TRUE, sep=",") else genos=marker.file
  
  # this function generates G matrices using different methods
  # option ( 1=genomic rel using frequencies,2=genomic rel using regression),  
  # sort (not implemented yet), data (first column id then SNPs, SNPs must be 0,1,2) no header, 
  # ped (id, sire, dam) no header
  #source("RelFxn.R")
  GenomicRel = function(option,data,ped,sort=FALSE){

    if(option==1){
      M1=data
      
      M= M1[,2:ncol(M1)]-1
      p1=round((apply(M,2,sum)+nrow(M))/(nrow(M)*2),3)
      p=2*(p1-.5)
      P = matrix(p,byrow=T,nrow=nrow(M),ncol=ncol(M))
      Z = as.matrix(M-P)
      
      b=1-p1
      c=p1*b
      d=2*(sum(c))
      
      ZZt = Z %*% t(Z)
      G = (ZZt/d)
      invG=solve(G)
    }
    if(option==2){
      names(ped)=c("id","father","mother")
      ped=as.Pedigree(ped)
      fullA=relationshipAdditive(ped)
      rowName=as.numeric(rownames(fullA))
      A=fullA[match(data[,1],rowName),match(data[,1],rowName)]
      n=(nrow(A))^2
      sumA=sum(A)
      SqA=A*A
      sumSqA=sum(SqA)
      
      M=as.matrix((data[,2:ncol(data)])-1)
      MtM=M%*%t(M)
      
      rhs1=sum(MtM)
      rhs2=sum(MtM*A)
      
      lhs=cbind(rbind(n,sumA),rbind(sumA,sumSqA))
      rhs=rbind(rhs1,rhs2)
      g=solve(lhs,rhs)
      #g=ginv(lhs)%*%rhs
      
      one=matrix(1,nrow=nrow(data))
      
      G=(MtM-(g[1,1]*(one%*%t(one))))/g[2,1]
      invG=solve(G)
    }
    #giv = matrix(NA,nrow=((nrow(G)*nrow(G))/2),ncol=3)
    col1=NA
    col2=NA
    col3=NA
    
    print("Still working, this takes awhile")
    
    for (i in 1:nrow(invG)){
      for (j in 1:ncol(invG)){
        if (i >= j){
          col1=cbind(col1,i)
          col2=cbind(col2,j)
          col3=cbind(col3,invG[i,j])
        }	
      }	
      
    }
    Ginv=cbind(t(col1),t(col2),t(col3))
    row.names(Ginv)=c(0:(nrow(Ginv)-1))
    #write.table(Ginv[-1,],"GinvTEST.giv",sep=" ", row.names=F, col.names=F,quote = F)
    return(Ginv[-1,])
  }
  
  # Generating G matrix based on frequencies
  if(Goptions==1) Ginv <-GenomicRel(option=1,data=genos, sort=FALSE)
  # Generating G matrix based on the regression method (option=2)
  if(Goptions==2) Ginv <-GenomicRel(option=2,data=genos, ped=pedigree, sort=FALSE)
  
  colnames(Ginv)<-c("row","column","value")
  attributes(Ginv)$rowNames<- aped.rowNames

  return(Ginv)
}
