\name{ir.sp}
\alias{ir.sp}
\docType{data}
\title{
Irregular spatial dataset
}
\description{
An irregular spatial dataset example to be transformed to regular spatial data.
}
\usage{data(ir.sp)}
\format{
  A data frame with 324 observations on the following 16 variables.
  \describe{
    \item{\code{Provence}}{a factor with levels \code{21} \code{22} \code{23}}
    \item{\code{Plustree}}{a factor with levels \code{21001} \code{21002} \code{21003} \code{21004} \code{21007} \code{21007-5} \code{21009} \code{21010} \code{21012} \code{21013} \code{21014} \code{22001} \code{22002} \code{22003} \code{22004} \code{22005} \code{22006} \code{22007} \code{22008} \code{22009} \code{22010} \code{22011} \code{22013} \code{22014} \code{22015} \code{22016} \code{22017} \code{22018} \code{22019} \code{22020} \code{22021} \code{L1} \code{L11} \code{L2} \code{L21} \code{L22} \code{L24} \code{L26} \code{L27} \code{L3} \code{L32} \code{L33} \code{L34} \code{L4} \code{L5} \code{L52} \code{L54} \code{L55} \code{L56} \code{L57} \code{L6} \code{L60} \code{L61} \code{L62} \code{L8} \code{L9}}
    \item{\code{Family}}{a factor with levels \code{201} \code{202} \code{203} \code{204} \code{205} \code{206} \code{207} \code{208} \code{209} \code{210} \code{211} \code{212} \code{213} \code{214} \code{215} \code{216} \code{217} \code{218} \code{219} \code{220} \code{221} \code{222} \code{224} \code{225} \code{226} \code{227} \code{228} \code{229} \code{230} \code{231} \code{232} \code{233} \code{234} \code{235} \code{236} \code{237} \code{238} \code{240} \code{241} \code{243} \code{244} \code{245} \code{247} \code{249} \code{250} \code{252} \code{253} \code{254} \code{256} \code{258} \code{259} \code{260} \code{261} \code{263} \code{264} \code{265}}
    \item{\code{Block}}{a factor with levels \code{2}}
    \item{\code{Row}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6} \code{7} \code{8} \code{9} \code{10} \code{11} \code{12} \code{13} \code{14} \code{15}}
    \item{\code{Col}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6} \code{7} \code{8} \code{9} \code{10} \code{11} \code{12} \code{13} \code{14} \code{15} \code{16} \code{17} \code{18} \code{19} \code{20} \code{21} \code{22} \code{23} \code{24} \code{25} \code{26} \code{27} \code{28} \code{29} \code{30}}
    \item{\code{h05}}{a numeric vector}
    \item{\code{cw05}}{a numeric vector}
    \item{\code{gd05}}{a numeric vector}
    \item{\code{d05}}{a numeric vector}
    \item{\code{h07}}{a numeric vector}
    \item{\code{cw07}}{a numeric vector}
    \item{\code{d07}}{a numeric vector}
    \item{\code{h10}}{a numeric vector}
    \item{\code{d10}}{a numeric vector}
    \item{\code{v10}}{a numeric vector}
  }
}

\references{
Yuanzhen Lin, Xiaoyang Chen. R & ASReml-R Statistical Analysis Tutorial. China Forestry Publishing House. 2014
}
\examples{
data(ir.sp)
## maybe str(ir.sp) ; plot(ir.sp) ...
}
\keyword{datasets}
