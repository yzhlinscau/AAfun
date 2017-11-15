mc.ped <-
function(ped){
  ped[ped==0]=NA
  ped[ped=='.']=NA
  
  for(i in 1:3){ped[,i]<-as.factor(ped[,i])}
  p1a<-levels(ped[,2])
  p1b<-levels(ped[,3])
  p1<-c(p1a,p1b)
  ped1<-data.frame(p1,p2=0,p3=0)
  names(ped1)<-names(ped)
  ped2<-rbind(ped1,ped)
  ped2[is.na(ped2)]=0
  ped2   
}
