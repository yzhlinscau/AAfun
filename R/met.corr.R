met.corr <-
function(object,site,kn=NULL){
  #require(amap)
  #require(corrgram)
  
  if(is.null(kn)) kn=3
  n<-nlevels(site)
  vect1<-summary(object)$varcomp$component[1:n]
  w.var<-diag(vect1)
  vect2<-summary(object)$varcomp$component[(n+1):(3*n)]
  t.var<-matrix(vect2,nrow=n)
  
  wt.var<-t.var%*%t(t.var)+w.var
  df<-wt.var
  
  for(i in 1:(n-1)){
    for(j in 2:n){
      if(i<j){df[i,j]<-df[j,i]/(sqrt(df[i,i]*df[j,j]))
              j<-j+1}
    }
    i<-i+1
  }
  
  df.2<-df
  
  for(i in 1:(n-1)){
    for(j in 2:n){
      if(i<j){df[j,i]<-df[i,j]
              j<-j+1}
    }
    i<-i+1
  }
  diag(df)<-1
  rownames(df)<-c(paste("S",levels(site),sep=''))
  colnames(df)<-c(paste("S",levels(site),sep=''))  
  
  chcluster <- hclusterpar(na.omit(df), method="manhattan")
  windows(10,8)
  plot(chcluster, main="Fig.1 Cluster of different sites",hang=-1)  #  labels=F
  rect.hclust(chcluster, k=kn)
  cat("Site cluster results:\n")
  print(tree.id<-cutree(chcluster,k=kn))
  
  if(n<16){
    windows(10,8)
    corrgram(df, type="cor",order=T, lower.panel=panel.pie,
               upper.panel=panel.conf, text.panel=panel.txt,
               main="Fig.2 Correlogram of different sites")
  }
  
  cat("\nCov\\Var\\Corr matrix\n\n")
  rownames(df.2)<-c(paste("S",levels(site),sep=''))
  colnames(df.2)<-c(paste("S",levels(site),sep=''))
  print(round(df.2,3))
  cat("\n")
}
