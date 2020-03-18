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

