##
## demo file for spd.plot. 
##

library(AAfun)
######## example 1 plot regular spatial data
data(barley)


aim.trait<-subset(barley,select=c(Row,Column,yield))
spd.plot(aim.trait)
spd.plot(aim.trait,color.p=topo.colors)
spd.plot(aim.trait,key.unit="Kg")
spd.plot(aim.trait,p.lbls="barley",x.unit=2,y.unit=1)

#AR1*AR1
barley1.asr<-asreml(yield~Variety, rcov =~ ar1(Row):ar1(Column), data=barley)

summary(barley1.asr)$varcomp
plot(variogram(barley1.asr),main="M1")

aa=variogram(barley1.asr)
spd.plot(aa,type="variogram",color.p=topo.colors)

######## example 2 plot spatial data with NA's
data(ir.sp)

ir.sp2<-ir.sp[,5:16] # order: Row,Col,h05,cw05,...
#ir.sp2<-subset(ir.sp,select=c(Row,Col,h05,cw05))

sp1<-ir2r.sp(ir.sp2,row.max=10,col.max=20)

aim.trait=subset(sp1,select=c(Row,Col,d10))
spd.plot(aim.trait,key.unit="cm")
spd.plot(aim.trait,color.p=topo.colors,na=0)
spd.plot(aim.trait,na=0,x.unit=3)


