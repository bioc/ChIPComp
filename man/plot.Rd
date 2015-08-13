\name{plot.ChIPComp}
\alias{plot.ChIPComp}
\title{plot correlation between log ChIP read counts and smoothing control counts in common binding sites.}

\description{
	plot correlation between log ChIP counts and smoothing control counts in common binding sites.
}

\usage{
	\method{plot}{ChIPComp}(x,...)
}

\arguments{
  \item{x}{ A \code{ChIPComp} object.}
  \item{\dots}{Other graphical parameters to \code{plot}}
}

\value{
	Plot the correlation between ChIP sample and control sample
}

\author{Hao Wu<hao.wu@emory.edu>, Li Chen <li.chen@emory.edu>}

\examples{
	data(seqData)
	plot(seqData)
}

