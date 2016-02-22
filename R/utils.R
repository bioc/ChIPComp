read.BAM<-function(fn){
	param=ScanBamParam(what=c("rname","strand","pos","qwidth"))
	bam=scanBam(fn,param=param)[[1]]
	ix=!is.na(bam$rname) & !is.na(bam$pos)
	qwidth=bam$qwidth[ix]
	IRange.reads=GRanges(seqnames=Rle(bam$rname[ix]),ranges=IRanges(bam$pos[ix], width=bam$qwidth[ix]))
	IRange.reads
}

getWinCounts<-function(files,wins,filetype=c("bed","bam")){
	if(class(wins)!="data.frame" & class(wins)!="GRanges")
		stop("Input genomic intervals must be a GRanges or data frame!")
	if(class(wins)=="data.frame") wins=import(wins)
	counts=matrix(0,length(wins),length(files))
	for(i in 1:length(files)) {
  		if(filetype=="bam") reads=read.BAM(files[i])
    	else if(filetype=="bed") reads=import(files[i])
    	counts[,i]=countOverlaps(wins,reads)
  	}
  	colnames(counts)=files
  	counts
}


getCTCounts<-function(files,peak.gr,filetype=c("bed","bam"),species=c("hg19","mm9"),binsize,mva.span){

	n=length(files)
	p=length(peak.gr)
	if(species=="hg19"){
		wins=tileGenome(seqinfo(BSgenome.Hsapiens.UCSC.hg19), tilewidth=binsize,cut.last.tile.in.chrom=TRUE)
	}else if(species=="mm9"){ 
		wins=tileGenome(seqinfo(BSgenome.Mmusculus.UCSC.mm9), tilewidth=binsize,cut.last.tile.in.chrom=TRUE)
	}
	tmp=findOverlaps(peak.gr,wins)
  	counts.peak=matrix(0,p,n)
  	chrs=seqlengths(wins)
  	
  	for(i in seq_len(n)){
  		if(filetype=="bam") reads=read.BAM(files[i])
    	else if(filetype=="bed") reads=import(files[i])
  		binCounts=countOverlaps(wins,reads)
  		binCounts.rm=rmdup(binCounts,sum(binCounts),sum(chrs/1000),binsize/1000)
  		binCounts.sm=MACS.mva(binCounts.rm, mva.span=mva.span,binsize)
  		counts=tapply(binCounts.sm[subjectHits(tmp)],queryHits(tmp), sum)
		id1=as.numeric(intersect(names(counts),seq_len(p)))
  		id2=as.numeric(setdiff(names(counts),seq_len(p)))
  		counts.peak[id1,i]=counts
  		counts.peak[id2,i]=1
  		counts.peak[,i][counts.peak[,i]>quantile(counts.peak[,i],0.95)]=quantile(counts.peak[,i],0.95)
  	}
  	counts.peak
}



regress<-function(ip, ct, ix.commonPeak,method=c("spline","lm")){
  method=match.arg(method)
  xx=log(ct+1)
  yy=log(ip+1)
  if(method=="lm") {
    fit=lm(yy~xx, subset=ix.commonPeak)
    cc=coef(fit)
    std=sd(resid(fit), na.rm=TRUE)
    yfit=cc[1] + cc[2]*xx
    res=(yy-yfit)/ std
  }else if(method=="spline"){
    fit=smooth.spline(xx[ix.commonPeak], yy[ix.commonPeak],df=3)
    std=sd(resid(fit), na.rm=TRUE)
    yfit=predict(fit, xx)$y
    res=yy-yfit
  }
  list(yfit=yfit, resid=res)
}


findCommonPeak<-function(ipmat,ctmat){
  ix.peak=rep(TRUE,nrow(ipmat))
  for(i in 1:ncol(ipmat)){
    k=sum(ipmat[,i])/sum(ctmat[,i])
    pp=1- ppois(ipmat[,i], ctmat[,i]*k)
    ix.peak=ix.peak & pp<0.01
  }
  ix.peak
}


MACS.mva<-function(x,mva.span=c(1000,5000,10000),binsize=50){
  nmvg=mva.span/binsize
  loc.lambda=mean(x)
  for(i in 1:length(nmvg)){
  	#dyn.load("~/Dropbox/dbind/bioconductor/ChIPComp/src/mva.so")
    out=.Call("mva",x=as.double(x),z=as.integer(nmvg[i]))
    loc.lambda=pmax(loc.lambda,out)
  }
  loc.lambda
}


rmdup<-function(counts,countsN,Lg,Lw,pval=1e-5){
  if(is.vector(counts)){
    cutoff=qbinom(1-pval,countsN,Lw/Lg)
    counts[counts>cutoff]=cutoff
  }else{
    for(i in seq(ncol(counts))){
      cutoff=qbinom(1-pval,countsN[i],Lw/Lg)
      counts[,i][counts[,i]>cutoff]=cutoff
    }
  }
  counts
}


findReplicate<-function(design){
	reps=double(nrow(design))
	ix=duplicated(design)
	for(i in seq(nrow(design))){
		if(!ix[i]){
			reps[i]=1
			k=1
		}else{
			reps[i]=k+1
			k=k+1
		}
	}
	reps
}


