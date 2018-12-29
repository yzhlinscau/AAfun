met.plot <-
function(object,plot.title=NULL,...){
  #require(desplot) # V1.4
  #require(grid)
  #require(reshape2)
  
  if(is.null(plot.title)) plot.title<-"MET data plot"
  dat<-object
  levels(dat[,3])<-paste("S",1:nlevels(dat[,3]),sep="")
  names(dat)[1:5]<-c("genotype","yield","site","row","col")
  for(i in 4:5){dat[,i]<-as.numeric(dat[,i])}
  #windows(10,8)
  # desplot(yield~ col*row|site, dat, main=plot.title)
  if(length(dat)==5){  
    desplot::desplot(yield~ col*row|site, dat, main=plot.title)
    }else{    
      names(dat)[6:7]<-c("Rep","Blk")  
      desplot::desplot(yield ~ col*row|site, dat, main=plot.title,
          out1=Rep, out2=Blk,strip.cex=1.5,
          out1.gpar=list(col="blue", lwd=4),
          out2.gpar=list(col="red", lwd=1, lty=1),
          par.settings = list(layout.heights=list(strip=2)))
  }
  #windows(10,8)
  #desplot(genotype~ col*row|site, dat,main="Genotype plot placement")  
}
