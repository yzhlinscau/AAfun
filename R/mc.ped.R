mc.ped <-
function(ped){
  for(i in 2:3){ped[,i]<-as.factor(ped[,i])}
  p1a<-levels(ped[,2])
  p1b<-levels(ped[,3])
  p1<-c(p1a,p1b)
  ped1<-data.frame(p1,p2=0,p3=0)
  names(ped1)<-names(ped)
  ped2<-rbind(ped1,ped)
  ped2   
}
