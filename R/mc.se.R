mc.se <-
  function(object=NULL,Nmc=NULL,confinterval=NULL,lv=NULL,
           uv=NULL,n=NULL,conf.level=NULL,sigf=NULL){
    #object<-mc.model$VCV or h2/corr formula 
    if(is.null(sigf)) sigf=FALSE
    if(is.null(Nmc)) Nmc=TRUE
    if(Nmc==TRUE){         
      if(is.null(confinterval)) confinterval=HPDinterval(object) 
    }
    
    if(is.null(lv)) lv=confinterval[,1]
    if(is.null(uv)) uv=confinterval[,2]  
    if(is.null(n)) n=1000
    if(is.null(conf.level)) conf.level=0.95
    a=1-conf.level
    ta=qt(1-a/2,n-1)
    se=vector() #NULL
    for(i in 1:length(lv)){
      se[i]=(uv[i]-lv[i])/(2*ta)
    }
    se<-matrix(se,nrow=length(lv))
    if(length(lv)>1){rownames(se)<-rownames(confinterval)
    colnames(se)<-"se"} else {colnames(se)<-"se"}
    
    sig.level<-function(tvalue,se,...){
      n<-length(se)
      siglevel<-1:n
      for(i in 1:n){    
        sig.se<-c(se[i]*1.450,se[i]*2.326,se[i]*3.090)  
        if(abs(tvalue[i])>sig.se[1]){siglevel[i]<-"*"}else{siglevel[i]<-"Not signif"}
        if(abs(tvalue[i])>sig.se[2]){siglevel[i]<-"**"}
        if(abs(tvalue[i])>sig.se[3]){siglevel[i]<-"***"}
      }
      siglevel
    }  
    
    if(Nmc==TRUE){  
      var<-round(posterior.mode(object),3)
      se<-round(se[,1],3) 
      z.ratio<-round((var/se),3)
      df<-data.frame(var=var,se=se,z.ratio=z.ratio)
      
      tvalue<-as.vector(df[,1])
      se<-as.vector(df[,2])
      tname<-rownames(df)    
      siglevel<-sig.level(tvalue,se)
      
      if(length(df[,1])==1) rownames(df)<-deparse(substitute(object))
      if(sigf==FALSE){        
        return(df)
      }
      
      if(sigf==TRUE){
        options(digits=3) 
        print(data.frame(df, sig.level=siglevel))
        cat("---------------")
        cat("\nSig.level: 0'***' 0.001 '**' 0.01 '*' 0.05 'Not signif' 1\n\n")            
      }
      
    }
    if(Nmc==FALSE){
      return(round(se,3))
    }
  }
