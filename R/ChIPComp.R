ChIPComp<-function(countSet,A,threshold=1){
	
	ip=countSet$db[,grep("ip",colnames(countSet$db))]
	ct=countSet$db[,grep("ct",colnames(countSet$db))]
	if(missing(A)) ix.commonPeak=countSet$db$commonPeak == 1
	else{
		m1=makeGRangesFromDataFrame(countSet$db)
		m2=import(A)
		a=findOverlaps(m1,m2)
		ix.commonPeak=unique(queryHits(a))	
	}
	
	signals=glambda=ip
	X=as.matrix(countSet$design)
	for(i in seq_len(ncol(ip))) {
    	tmp=regress(ip[,i], ct[,i], ix.commonPeak)
    	signals[,i]=tmp$resid
    	glambda[,i]=tmp$yfit
 	}
 	limmaFit=lmFit(signals, design=X)
  	limma.res=eBayes(limmaFit)
    beta.est=limma.res$coefficients
  	XXX=solve(t(X)%*%X)%*%t(X)
  	mu0=-X%*%t(beta.est)- t(glambda)
  	emu0=sweep(mu0,2,limma.res$s2.post/2,FUN='+')
  	emu0=exp(emu0)
  	vy=sweep(emu0,2,limma.res$s2.post,FUN='+')
  	vbeta.cond=apply(vy, 2, function(x){
   	 					cc=XXX%*%diag(x)%*%t(XXX)
   	 					cc[2,2]
   	 			})			
  	wstat=limma.res$coefficients[,2]/sqrt(vbeta.cond)
	countSet$db$pvalue.wald=2*(1-pnorm(abs(wstat)))

	prob.post=1-pnorm(threshold,beta.est[,2],sqrt(vbeta.cond))+pnorm(-threshold,beta.est[,2],sqrt(vbeta.cond))
  	countSet$db$prob.post=prob.post
	  	
	structure(countSet,class="ChIPComp")
}




print.ChIPComp<-function(x,topK=10,...){
	print(x$db[order(x$db$prob.post,decreasing=TRUE),][seq(topK),])
	
}


plot.ChIPComp<-function(x,...){
	mypar=function (a = 2, b = 2, brewer.n = 8, brewer.name = "Dark2", ...) {
  		par(mar = c(2.5, 2.5, 1.6, 1.1), mgp = c(1.5, 0.5, 0))
  		par(mfrow = c(a, b), ...)
	}
	len=x$db$end-x$db$start
	ix.len=len>500 & len<10000
	ip=x$db[,grep("ip",colnames(x$db))]
	ct=x$db[,grep("ct",colnames(x$db))]
	mains=gsub("ip_","",colnames(x$db)[grep("ip",colnames(x$db))])
	ix.commonPeak=x$db$commonPeak == 1
	ix.commonPeak=ix.commonPeak & ix.len
	n=ncol(ip[,,drop=FALSE])
	mypar(sqrt(n),sqrt(n))
	for(i in seq_len(n)){
		x=log(ct[,i]+1); y=log(ip[,i]+1)
  		fit=lm(y~x)
 		beta=round(coef(fit)[2], 2)
  		if(length(x)>5000) {
    		ix=sample(1:length(x), 5000)
    		x=x[ix]; y=y[ix]
  		}
  		plot(x, y, log="", col="#00000015",cex=0.4, xlab="log(background)", ylab="log(IP)",main=mains[i])
  		fit=smooth.spline(x, y, df=5)
    	lines(fit, lwd=2, lty=2, col=2)
	}
	on.exit()
}







