\name{makeConf}
\alias{makeConf}
\title{make configurations for experimental design written in csv sheet}
\description{
	 Make a list with two elements. 
	 The first element is a data frame containing two group comparison study information.
	 The second element is the design matrix.
}
\usage{
	makeConf(sampleSheet)
}
\arguments{
	
  \item{sampleSheet}{
  	A csv sheet represents ChIP experiments design. 
  	It contains 6 columns,\code{sampleID},\code{condition},\code{factor},\code{ipReads},\code{ctReads},\code{peaks}.
  	\code{condition} refers to treatment condition or cell line; 
  	\code{factor} refers to transcription factor or histone modification; 
  	\code{ipReads} is the ChIP sequence data in bam or bed format; 
  	\code{ctReads} is the control sequence data in bam or bed format; 
  	\code{peaks} is the called peaks from existing peak-calling software.
  }
}

\value{
 	A list with two elements. 
	The first element is a data frame containing two group comparison study information.
	The second element is the design matrix.
}

\author{Hao Wu<hao.wu@emory.edu>, Li Chen <li.chen@emory.edu>}

\examples{

	confs=makeConf(system.file("extdata", "conf.csv", package="ChIPComp"))
	conf=confs$conf
	design=confs$design

}








