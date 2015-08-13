\name{print.ChIPComp}
\alias{print.ChIPComp}
\title{
	Print top ranked differential binding sites
}
\description{
	Print top differential binding sites ranked by posterior probability in a decreasing order.
}

\usage{	
	\method{print}{ChIPComp}(x,topK=10,...)
}

\arguments{
  \item{x}{ A \code{ChIPComp} object. }
  \item{topK}{top K differential binding sites. Default is 10.}
  \item{\dots}{Other parameters to \code{print}}
}

\value{
	Print differential binding sites ranked by posterior probability
}


\author{Hao Wu<hao.wu@emory.edu>, Li Chen <li.chen@emory.edu>}

\examples{
	data(seqData)
	seqData=ChIPComp(seqData)
	print(seqData)
}







