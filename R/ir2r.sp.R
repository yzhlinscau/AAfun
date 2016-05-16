ir2r.sp<-
function(object,row.max=NULL,col.max=NULL){
  # Row=NULL,Col=NULL,
  #require(sqldf)
  
  df<-object
  #names(df)[1:2]<-c("Row","Col")
  #attach(object)
  #if(is.null(Row)) df$Row=object$Row
  #if(is.null(Col)) df$Col=object$Col
  
  df$Row<-as.numeric(df$Row) # df$Row
  df$Col<-as.numeric(df$Col) # df$Col
  rmax<-max(df$Row,na.rm=TRUE)
  cmax<-max(df$Col,na.rm=TRUE)
  if(is.null(row.max)) row.max=rmax
  if(is.null(col.max)) col.max=cmax
  
  nr<-nrow(df)
  k=0
  for(i in 1:row.max){
    df2=subset(df,Row==i)
    a<-as.data.frame(1:col.max) # for all col
    b<-as.data.frame(df2$Col) # for data col
    nc<-length(b[,1])
    if(nc<col.max){
      aNotInb<-sqldf('SELECT * FROM a EXCEPT SELECT * FROM b')
      dc<-nrow(aNotInb)
      df[(nr+1+k):(nr+dc+k),]<-NA
      df[(nr+1+k):(nr+dc+k),1]<-i # for Row--1
      df[(nr+1+k):(nr+dc+k),2]<-aNotInb[,1] # for Col--2      
    } else dc=0
    k=k+dc
  }
  if(row.max<rmax) df<-subset(df,Row<(row.max+1))
  if(col.max<cmax) df<-subset(df,Col<(col.max+1))
  df<-arrange(df,Row,Col)
  for(i in 1:2) df[,i]<-as.factor(df[,i])
  
  cat("Row Number is:",row.max,"Col Number is:",col.max)
  #detach(df)
  
  df
}
