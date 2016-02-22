
makeConf=function(sampleSheet){
	conf=read.csv(sampleSheet,stringsAsFactors=TRUE)
	design=as.data.frame(lapply(conf[,c("condition","factor")],as.numeric))-1
	design=as.data.frame(model.matrix(~condition,design))
	list(conf=conf,design=design)
}


makePeakSet=function(peaks,peak.center,peak.ext){

	message("Making peak list......\n")
	n=length(peaks)
	peak.list=GRangesList()
	pmat=import(peaks[1])
	if(peak.center){
		start(pmat)=(start(pmat)+end(pmat))/2-peak.ext
		end(pmat)=(start(pmat)+end(pmat))/2+peak.ext
		start(pmat)[start(pmat)<=0]=1
	}	
	peak.list[[1]]=pmat
	for(i in 2:n){
		mat=import(peaks[i])
		peak.list[[i]]=mat
		if(peak.center){
			start(pmat)=(start(pmat)+end(pmat))/2-peak.ext
			end(pmat)=(start(pmat)+end(pmat))/2+peak.ext
			start(pmat)[start(pmat)<=0]=1
		}
		pmat=union(pmat,mat)
	}
	tmp=findOverlaps(peak.list[[1]],pmat)
	oidx=unique(subjectHits(tmp))	
	for(i in 2:n){
		tmp=findOverlaps(peak.list[[i]],pmat)
		oidx=intersect(oidx, unique(subjectHits(tmp)))	
	}
	peakSet=list(peak.list=peak.list,pmat=pmat,oidx=oidx)
	peakSet
}


makeIPSet=function(ips,peakSet,filetype){

	message("Making ip counts......")
	ip.allpeak=getWinCounts(ips,peakSet$pmat,filetype)
	ip.opeak=ip.allpeak[peakSet$oidx,]
	ipSet=list(ip.allpeak=ip.allpeak,ip.opeak=ip.opeak)
	ipSet
}


makeCTSet=function(cts,peakSet,filetype,species,binsize,mva.span){

	message("Making control counts......\n")
	ct.allpeak=getCTCounts(cts,peakSet$pmat,filetype,species,binsize,mva.span)	
	ct.opeak=ct.allpeak[peakSet$oidx,]
	controlSet=list(ct.allpeak=ct.allpeak,ct.opeak=ct.opeak)
	controlSet
}



makeCountSet=function(conf,design,filetype=c("bed","bam"),species=c("hg19","mm9"),peak.center=FALSE,peak.ext=0,binsize=50,mva.span=c(1000,5000,10000)){
				
		if(missing(conf))
			stop("The configuration data frame should be provided!")
		if(missing(design))
			stop("The design matrix should be provided!")
		filetype=match.arg(filetype)
		species=match.arg(species)
		
		peaks=as.vector(conf$peaks)
		ips=as.vector(conf$ipReads)
		cts=as.vector(conf$ctReads)
		peakSet=makePeakSet(peaks,peak.center,peak.ext)
		ipSet=makeIPSet(ips,peakSet,filetype)
		controlSet=makeCTSet(cts,peakSet,filetype,species,binsize,mva.span)
		reps=findReplicate(design)
		
		if(is.null(design$factor)){
			ipNames=paste("ip_c",design[,2],"_r",reps,sep="")
			ctNames=paste("ct_c",design[,2],"_r",reps,sep="")
		}else{
			ipNames=paste("ip_c",design[,2],"_f",design[,3],"_r",reps,sep="")
			ctNames=paste("ct_c",design[,2],"_f",design[,3],"_r",reps,sep="")
		}
		
		commonPeak=numeric(length(peakSet$pmat))
		commonPeak[peakSet$oidx]=1
		db=cbind( as.data.frame(peakSet$pmat)[,seq(3)],ipSet$ip.allpeak,controlSet$ct.allpeak,commonPeak)	
		colnames(db)= c("chr","start","end", ipNames, ctNames,"commonPeak")
		countSet=list(db=db,design=design)
		structure(countSet,class="ChIPComp")
}














