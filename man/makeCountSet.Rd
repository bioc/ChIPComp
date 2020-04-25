\name{makeCountSet}
\alias{makeCountSet}
\title{make differential binding sites data frame}
\description{
	This is an utility function to create a data frame. 
	The data frame contains binding sites merged by peaks from two conditions, 
	count ChIP read counts, smoothing control counts for each candidate region, 
	and indicate the common peaks from two conditions.
}

\usage{
	makeCountSet(conf,design,filetype,species,peak.center=FALSE,peak.ext=0,binsize=50,mva.span=c(1000,5000,10000))
}

\arguments{
	\item{conf}{
		A data frame that represents the ChIP experiments information. 
		It contains 6 columns,\code{sampleID},\code{condition},\code{factor},\code{ipReads},\code{ctReads},\code{peaks}.
  		\code{condition} refers to treatment condition or cell line; 
  		\code{factor} refers to transcription factor or histone modification; 
  		\code{ipReads} is the ChIP sequence data in bam or bed format; 
  		\code{ctReads} is the control sequence data in bam or bed format; 
  		\code{peaks} is the called peaks from existing peak-calling software.
	}
	\item{design}{
		Two column design matrix. The number of rows equals number of ChIP samples from two conditions.
		The first column are all 1s, which indicates intercept in regression model.
		The second column are 1s for one condition and 0s for another condition.
	}
	\item{filetype}{
		Two sequence file types are supported (bed or bam).
	}
	\item{species}{
		Two species are supported (hg19 or mm9). Other species are supported by specifying other.
	}
	\item{peak.center}{
		This argument is coupled with \code{peak.ext}. Default is FALSE. The argument 
		is used when centered regions of peaks are more of interest.
	}
	\item{peak.ext}{
		This argument is coupled with \code{peak.center}. Default is 0.
	}
	\item{binsize}{
		binsize in bp to calculate the smooth local lambda in poisson distribution. The default is 50bp.
	}
	\item{mva.span}{
		1 kb, 5 kb or 10 kb window centered at the peak location in the control sample.
	}

}


\value{
	A object \code{ChIPComp}.
	Column \code{chr},\code{start},\code{end} are the binding site genomic coordinate;
	Column \code{ip_c(\#condition)_r(\#replicate)} indicates the ChIP counts in \#replicate in \#condition;
	Column \code{ct_c(\#condition)_r(\#replicate)} indicates the smoothing control counts in \#replicate in \#condition; 
	Column \code{commonPeak} indicates the common binding sites.
}




\examples{
	conf=data.frame(
		SampleID=1:4,
		condition=c("Helas3","Helas3","K562","K562"),
		factor=c("H3k27ac","H3k27ac","H3k27ac","H3k27ac"),
		ipReads=system.file("extdata",c("Helas3.ip1.bed","Helas3.ip2.bed","K562.ip1.bed","K562.ip2.bed"),package="ChIPComp"),
		ctReads=system.file("extdata",c("Helas3.ct.bed","Helas3.ct.bed","K562.ct.bed","K562.ct.bed"),package="ChIPComp"),
		peaks=system.file("extdata",c("Helas3.peak.bed","Helas3.peak.bed","K562.peak.bed","K562.peak.bed"),package="ChIPComp")
	)
	conf$condition=factor(conf$condition)
  conf$factor=factor(conf$factor)
	design=as.data.frame(lapply(conf[,c("condition","factor")],as.numeric))-1
	design=as.data.frame(model.matrix(~condition,design))
	countSet=makeCountSet(conf,design,filetype="bed", species="hg19",binsize=1000)
}



