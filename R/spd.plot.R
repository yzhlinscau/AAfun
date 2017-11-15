spd.plot<-
function(object,type="data",p.lbls=NULL,key.unit=NULL,   
             x.unit=NULL,y.unit=NULL,na=NULL,
             color.p=NULL,mtitle=NULL,...){
  #require(plyr)  
  par(mar=c(4,4,2,2), cex.main=1)
  if(is.null(color.p)) color.p=terrain.colors
  
  if(type=="data"){
    for(i in 1:2){object[,i]<-as.numeric(object[,i])}
    #object<-arrange(object,object[2],object[1])
    object<-object[order(object[2],object[1]),]
  }
  if(type=="variogram"){object=object[,-4];for(i in 1:2) object[,i]=object[,i]+1} 
  
  ncol<-max(object[2])
  lbls<-names(object[3])
  lbls2<-paste(lbls,key.unit,sep="\n(")
  lbls2<-paste(lbls2,"",sep=")")
  object.1<-(object[,3])
  df <-matrix(object.1,nrow=ncol,byrow=TRUE)
  if(is.null(na)) na=1
  if(na==0)  df[is.na(df)]<-0.0001  # reduce effects of NA value in data 
  
  x = 1 : nrow(df) 
  y = 1 : ncol(df)
  
  #if(max(x)>20){
  #  if(max(x)%%2!=0){if(is.null(x.unit)) seq(1,max(x),by=2)
  #                    else{x.axis<-seq(1,max(x),by=x.unit)} 
  #                          }else{
  #                              if(is.null(x.unit)) x.axis<-seq(2,max(y),by=2)
  #                               else{x.axis<-seq(2,max(x),by=x.unit)}}
  #}else{if(is.null(x.unit)) x.axis<-seq(1,max(x),by=1)
  #       else{x.axis<-seq(1,max(x),by=x.unit)}}
  
  #if(max(y)>15){
  #  if(max(y)%%2!=0){if(is.null(y.unit))seq(1,max(y),by=2)
  #                    else{y.axis<-seq(1,max(y),by=y.unit)}
  #                           }else{
  #                                  if(is.null(y.unit)) y.axis<-seq(2,max(y),by=2)
  #                                   else{y.axis<-seq(2,max(y),by=y.unit)}}
  #}else{if(is.null(y.unit)) y.axis<-seq(1,max(y),by=1)
  #      else{y.axis<-seq(1,max(y),by=y.unit)}}
  
  if(is.null(x.unit)) x.unit=1
  if(is.null(y.unit)) y.unit=1
  x.axis<-seq(x.unit,max(x),by=x.unit)
  y.axis<-seq(y.unit,max(y),by=y.unit)
  
  if(is.null(p.lbls)) p.lbls=" " else p.lbls=paste(": ",p.lbls)
  if(is.null(key.unit)) lbls2=lbls else lbls2=lbls2
  if(is.null(mtitle)) mtitle=paste("The Topography of ",lbls, p.lbls) else mtitle=mtitle
  #windows(10,8)
  filled.contour(x, y, df, color.palette=color.p, 
                 plot.title = title( main=mtitle, 
                                     xlab="Col", ylab="Row"), 
                 plot.axes = { 
                   axis(1,x.axis) 
                   axis(2,y.axis) 
                 }, 
                 key.title = title(main=lbls2, cex.main=1.0) 
  )
  #abline(v=0, h=seq(1, max(y), by=1),lty=2,col="grey75") 
  #abline(h=0, v=seq(1, max(x), by=1),lty=2,col="gray75")
  
}
