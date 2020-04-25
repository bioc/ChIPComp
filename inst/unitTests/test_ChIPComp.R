
test.makeConf<-function(){
   	x=makeConf(system.file("extdata", "conf.csv", package="ChIPComp"))
    checkIdentical("list",as.character(class(x)))
    checkIdentical("data.frame",as.character(class(x$conf)))
    checkIdentical("data.frame",as.character(class(x$design)))
    checkEquals(length(x),2)
	checkEquals(ncol(x$conf),6)
	checkEquals(ncol(x$design),2)
}



test.makeCountSet<-function(){
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
 	x=makeCountSet(conf,design,filetype="bed", species="hg19",binsize=1000)
	checkIdentical("ChIPComp",as.character(class(x)))
  checkEquals(length(x),2)
	checkEquals(ncol(x$db),12)
	checkEquals(ncol(x$design),2)
  checkEquals(names(x), c("db","design"))
  checkEquals(names(x$db),c("chr","start","end","ip_c0_r1","ip_c0_r2","ip_c1_r1","ip_c1_r2","ct_c0_r1","ct_c0_r2","ct_c1_r1","ct_c1_r2","commonPeak"))
	checkException(makeCountSet(conf,design,filetype="bed", species="hg18"))
}
 
 
 
 test.ChIPComp<-function(){
   data(seqData)
 	x=ChIPComp(seqData)
    checkIdentical("ChIPComp",as.character(class(x)))
    checkEquals(length(x),2)
	checkEquals(ncol(x$db),14)
	checkEquals(ncol(x$design),2)
	checkEquals(names(x), c("db","design"))
    checkEquals(names(x$db),c("chr","start","end","ip_c0_r1","ip_c0_r2","ip_c1_r1","ip_c1_r2","ct_c0_r1","ct_c0_r2","ct_c1_r1","ct_c1_r2","commonPeak","prob.post","pvalue.wald"))
 }


 
 
 
