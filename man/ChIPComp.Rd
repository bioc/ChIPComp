\name{ChIPComp}
\alias{ChIPComp}
\title{ 
	Perform hypothesis testing to detect differential binding sites 
}
\description{
	Perform hypothesis testing to detect differential binding sites 
}
\usage{
	ChIPComp(countSet,A,threshold=1)
}
\arguments{
  \item{countSet}{ A \code{ChIPComp} object.}
  \item{A}{
  	User-specified regions to fit the model. 
  	It is a bed file with three columns, named ("chr","start","end"), 
  	could be separated by space or tab.
  }
  \item{threshold}{
  	User specified posterior probability threshold. Default is 1.
  }
}

\value{
	A object \code{ChIPComp} contains
	Column \code{chr},\code{start},\code{end} are the binding site genomic coordinate;
	Column \code{ip_c(\#condition)_r(\#replicate)} indicates ChIP counts in \#replicate in \#condition;
	Column \code{ct_c(\#condition)_r(\#replicate)} indicates smoothing control counts in \#replicate in \#condition; 
	Column \code{commonPeak} 1s indicate common binding sites;
	Column \code{prob.post} is the posterior probability for each binding site.
}

\author{Hao Wu<hao.wu@emory.edu>, Li Chen <li.chen@emory.edu>}

\examples{

	data(seqData)
	seqData=ChIPComp(seqData)

}








