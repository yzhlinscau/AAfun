\name{ir2r.sp}
\alias{ir2r.sp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Irregular to regular spatial data
}
\description{
This function changes irregular spatial data to regular data for further plot with spd.plot().
}
\usage{
ir2r.sp(object, row.max = NULL, col.max = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Irregular spatial data with Row and Col
}
  \item{row.max}{
Total row numbers for regular data
}
  \item{col.max}{
Total column numbers for regular data
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{object }{Irregular spatial data, Must have 'Row' and 'Col' as the first two column.}
%%  \item{comp2 }{Description of 'comp2'}
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

data(ir.sp)
ir.sp2<-ir.sp[,5:16] # order: Row,Col,h05,cw05,...
#ir.sp2<-subset(ir.sp,select=c(Row,Col,h05,cw05))

sp1<-ir2r.sp(ir.sp2,row.max=10)
sp2<-ir2r.sp(ir.sp2,col.max=20)
sp3<-ir2r.sp(ir.sp2,row.max=10,col.max=20)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ irregular.spatial }
\keyword{ ir2r.sp }% __ONLY ONE__ keyword per line
