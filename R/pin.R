pin <-
function(object, formula=NULL,signif=NULL, corN=NULL,Rdf=NULL){
  if(is.null(signif)) signif=FALSE
  if(is.null(Rdf)) Rdf=FALSE
  

  
  if(!is.null(formula)){
    transform<-formula
    #if(is.null(N)) N<-0
    pframe <- as.list(object$gammas)
    names(pframe) <- paste("V", seq(1, length(pframe)), sep = "")
    tvalue<-eval(deriv(transform[[length(transform)]], names(pframe)),pframe)
    X <- as.vector(attr(tvalue, "gradient"))
    X[object$gammas.type == 1] <- 0
    tname <- if(length(transform)==3){transform[[2]]}else ""
    n <- length(pframe)
    i <- rep(1:n, 1:n)
    j <- sequence(1:n)
    k <- 1 + (i > j)
    Vmat <- object$ai
    se <- sqrt(sum(Vmat * X[i] * X[j] * k))
    
    vv=vector() #NULL
    vv[1]=tvalue;vv[2]=se
    
    result<-data.frame(row.names=tname, Estimate=tvalue, SE=se)
    result1<-result
    result1$sig.level<-sig.level(tvalue,se)
  
    cat("\n")
    options(digits=3)
    if(signif==TRUE){ 
      print(result1)
      cat("---------------")
      cat("\nSig.level: 0'***' 0.001 '**' 0.01 '*' 0.05 'Not signif' 1\n")    
    }else{
      if(Rdf==TRUE) print(vv) else print(result)
      }
    cat("\n")
  }
  
  if(is.null(formula)){
    
    if(is.null(corN)){cat("\nAttension: since no N value, here just show fisrt one corr!!\n\n")
                    corN<-1} 
    n<-corN
    df<-summary(object)$varcomp
    tvalue<-as.vector(df[1:n,2])
    se<-as.vector(df[1:n,3])
    tname<-rownames(summary(object)$varcomp)[1:n]    
    siglevel<-sig.level(tvalue,se)
    
    options(digits=3) 
    print(data.frame(row.names=tname,Estimate=tvalue, SE=se, sig.level=siglevel))
    cat("---------------")
    cat("\nSig.level: 0'***' 0.001 '**' 0.01 '*' 0.05 'Not signif' 1\n\n")    
  }
}

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# sig.level functions

sig.level<-function(tvalue,se,...){
  n<-length(se)
  siglevel<-1:n
  for(i in 1:n){    
    sig.se<-c(se[i]*1.450,se[i]*2.326,se[i]*3.090)  
    
    if(abs(tvalue[i])>sig.se[3]) {siglevel[i]<-"***"}
     else if(abs(tvalue[i])>sig.se[2]) {siglevel[i]<-"**"}
     else if(abs(tvalue[i])>sig.se[1]) {siglevel[i]<-"*"}
     else {siglevel[i]<-"Not signif"}
  }
  siglevel
}

sig.level2=function(x){
  tt=vector() #NULL
  n=length(x)
  for(i in 1:n){
    if(abs(x[i])<0.001) tt[i]='***'
    else if(abs(x[i])<0.01) tt[i]='**'
    else if(abs(x[i])<0.05) tt[i]='*'
    else if(abs(x[i])<0.10) tt[i]='.' 
    else tt[i]=''
  }
  return(tt)
}
