\name{corrgram2}
\alias{corrgram2}
\alias{panel.bar}
\alias{panel.conf}
\alias{panel.density}
\alias{panel.ellipse}
\alias{panel.minmax}
\alias{panel.pie}
\alias{panel.pts}
\alias{panel.shade}
\alias{panel.txt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw a correlogram
}
\description{
The corrgram2 function produces a graphical display of a correlation matrix, called a correlogram. This function is rechanged from corrgram.
}
\usage{
corrgram2(x, type=NULL, order = FALSE, labels, panel=panel.shade,
  lower.panel=panel, upper.panel=panel, diag.panel=NULL,
  text.panel=textPanel, label.pos=c(0.5, 0.5), label.srt=0, cex.labels=NULL,
  font.labels=1, row1attop=TRUE, dir="", gap=0, abs=FALSE,
  col.regions=colorRampPalette(c("red","salmon","white","royalblue","navy")),
  cor.method="pearson",
  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A tall data frame with one observation per row, or a correlation matrix.
}
\item{type}{Use 'data' or 'cor'/'corr' to explicitly specify that 'x' is data or a correlation matrix.
}
\item{order}{Should variables be re-ordered? Use TRUE/"PCA" for PCA-based re-ordering. Options from the 'seriate' package include "OLO" for optimal leaf ordering, "GW", and "HC".
}
\item{...}{ Additional arguments samed to corrgram.
}
}

\value{
No value is returned. A plot is created.
}

\seealso{
Website for instant update: yzhlin-asreml.ys168.com
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

library(AAfun)
df=PrSpa
#names(df)

corrgram2(df[,6:13], type="data",order=T, lower.panel=panel.shade,
          upper.panel=panel.conf, text.panel=panel.txt,
          main="Correlogram Figures")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ corrgram2 }% __ONLY ONE__ keyword per line
