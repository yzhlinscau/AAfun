group.plot <-
function(object,x.lbls,y.lbls=NULL,y.zero=NULL,alpha=0.05, ...){
   #par(mar=c(5,6.5,4,2))
   
   #require(agricolae)
   #require(gplots)
   #require(plyr)
   
   if(is.null(y.zero)) y.zero=0

   object$groups2<-arrange(object$groups,object$groups$trt)
   lbls<-object$groups2[,3]
   #lbls<-toupper(lbls) # tolower()
   if(is.null(y.lbls)) y.lbls<-names(object$means[1])  
   df=sum(object$means[,3])-1
   
   mu.i <- object$means[,1]
   se.i <- qt(1-0.5*alpha, df) * object$means[,2] 
   
   #if(y.zero==TRUE) y.min<- ifelse(min(mu.i-se.i)>10,(min(mu.i-se.i)-5),0)
   #if(y.zero==FALSE) y.min=0
   
   y.min=y.zero
   y.max<-max(mu.i+se.i)+10
   y.max2<- mu.i + se.i +2
   
   #windows(10,8)
   bp <- barplot2(mu.i, names.arg=rownames(object$means),  
                  ylab=list(y.lbls, cex=1.5,font=2),xpd=FALSE,                  
                  ylim=c(min(y.min),y.max),density=10,font=2, col="red",
                  plot.ci=TRUE, ci.l=mu.i-se.i, ci.u=mu.i+se.i)
   text(bp, y.max2, lbls, cex=1.5,font=2, pos=3,col="blue")
   title(cex.main=1.5,font=2,main="Comparison between\ntreatment means",
         xlab=list(x.lbls, cex=1.5,font=2))
   box()
}
