dial.combn <-
function(lbls,N){
  n<-N  # dial mate number
  t=vector() #<-NULL
  if(n<10){ t<-paste(lbls,1:n,sep="")}
  if(n<100){ t[1:9]<-paste(lbls,1:9,sep="00")
             t[10:99]<-paste(lbls,10:99,sep="0")}  
  t1<-paste(lbls,1:n,sep="")
  
  mm<-function(t,n){
    tb<-diag(n)
    for(i in 1:n){
      for(j in 1:n){
        tb[i,j]<-paste(t[i],t[j],sep="")
        j<-j+1}
      i<-i+1
    }
    tb
  }
  tb1.v<-as.vector(mm(t1,n))
  f1.v<-rep(t1,rep(n,n))
  m1.v<-rep(t1,n)
  
  df<-mm(t,n)
  for(i in 1:(n-1)){
    for(j in 2:n){
      if(i<j){df[j,i]<-df[i,j];j<-j+1}
    }
    i<-i+1
  }
  df.v<-as.vector(df)
  
  dial<-data.frame(Male=m1.v,Female=f1.v,Recipro=tb1.v,Fam=df.v)
  
  b<-nlevels(dial[,4])
  if(b<10){levels(dial[,4])[1:b]<-paste("F",1:b,sep="00")}
  if(b>=10){
    levels(dial[,4])[1:9]<-paste("F",1:9,sep="00")
    levels(dial[,4])[10:b]<-paste("F",10:b,sep="0")
  }
  if(b>=100){
    levels(dial[,4])[1:9]<-paste("F",1:9,sep="00")
    levels(dial[,4])[10:99]<-paste("F",10:99,sep="0")
    levels(dial[,4])[100:b]<-paste("F",100:b,sep="")
  }  
  dial
  
}
