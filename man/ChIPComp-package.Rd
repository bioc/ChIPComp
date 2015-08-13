\name{ChIPComp-package}
\docType{package}
\alias{ChIPComp-package}
\title{
   Detect differential binding sites for ChIP sequencing data
}
\description{
  ChIPComp is an R library performing the differential binding analysis for
  ChIP-seq count data. Compared with other similar packages (DBChIP,
  DIME), ChIPComp considers the control samples in the process of detecting
  the differential binding sites. Extensive simulation results
  showed that ChIPComp performs favorably compared to DBChIP and DIME when
  the control samples are ignored.
  ChIPComp only works for two group comparison at this time, that is, to detect
  the differential binding sites for one transcription factor(histone)
  between two conditions (cell lines). We plan to extend the functionalities
  and make it work for more general experimental designs in the near future. 
}

\author{
  Hao Wu <hao.wu@emory.edu>, Li Chen <li.chen@emory.edu>
}

\keyword{package}
