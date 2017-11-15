\name{spd.plot}
\alias{spd.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot asreml-spatial data or variogram
}
\description{
This function plots spatial data for further spatial model or variogram result from spatial model by asreml, etc.
}
\usage{
spd.plot(object, type="data",p.lbls = NULL, key.unit = NULL,
          x.unit = NULL, y.unit = NULL, 
          na = NULL, color.p = NULL, mtitle=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Aim dataset
}
  \item{type}{
Type of dataset, "data" or "variogram"
}
  \item{p.lbls}{
Extra labels in figure title
}
  \item{key.unit}{
The unit of key
}
  \item{x.unit}{
Axis x least unit
}
  \item{y.unit}{
Axis y least unit
}
  \item{na}{
Transform NA to 0(na=0) or keep NA (default)
}
  \item{color.p}{
Parameters of the colors for figures
}
  \item{\mtitle}{
Parameters of the function filled.contour()
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{type }{Dataset type, default value is "data", when "variogram" for variogram.plot in spatial analysis in ASReml-R.}
  \item{key.unit }{Key unit, default value is empty.}
  \item{x.unit }{Axis x least unit, default value is 1.}
  \item{y.unit }{Axis y least unit, default value is 1.}
  \item{na }{Control NA value, default value is 1 to keep NA.}
  \item{color.p }{Parameters of the colors for figures, default value is terrain.colors, it could be rainbow, heat.colors, cm.colors and topo.colors, ...}
%% ...
}
\author{
Yuanzhen Lin <yzhlinscau@163.com>
}
\references{
Yuanzhen Lin. R & ASReml-R Statistics. China Forestry Publishing House. 2016
}
\seealso{
Website for instant update: yzhlin-asreml.ys168.com
}

\examples{
library(AAfun)
######## example 1 plot regular spatial data
data(barley)

aim.trait<-subset(barley,select=c(Row,Column,yield))
spd.plot(aim.trait)
spd.plot(aim.trait,color.p=topo.colors)
spd.plot(aim.trait,key.unit="Kg")
spd.plot(aim.trait,p.lbls="barley",x.unit=2,y.unit=1)

#AR1×AR1 model to plot variogram
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ asreml.spatial }
\keyword{ spd.plot }% __ONLY ONE__ keyword per line
