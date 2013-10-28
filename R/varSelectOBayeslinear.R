varSelectOBayeslinear<-function(y,X,X0,
	type="IP",prior.par,EFF=TRUE,
	model.prior.type="Beta-Binomial",model.prior.par=c(1,1),
	RWflag=TRUE,N.draws=10^4,shuf=.2,start,RWpropflag=FALSE,
	inference="both",alpha=0.05,m.considered=100){

	if(missing(y)){stop(" y was not specified")}	
	if(class(y)!="numeric"){stop(" y must be a numerical vector")}
	n<-length(as.vector(y))
	if(missing(X0)){X0<-matrix(1,n,1)}
	if(!is.numeric(X0)){stop(" X0 must be numerical")}
	if(class(X0)=="numeric"){X0<-matrix(X0,nrow=n)}
	if(class(X0)=="matrix"){if(dim(X0)[1]!=n){stop(" number of rows of X0 must be the same as the length of y")}}
	p0<-dim(X0)[2]
	if(min(svd(crossprod(X0[1:p0,]))$d)<1e-12){
		stop(" X0 has collinear columns")
	}
	if(min(svd(cbind(c(1),X0)[1:(p0+1),])$d)>1e-12){
		X0<-cbind(c(1),X0);
		p0<-p0+1
		warning(" Intercept is not collinear with X0, it has been added")
	}
	colnames(X0)<-paste("Base",1:p0,sep="_")
	if(missing(X)){stop(" X was not specified")}
	if(!is.numeric(X)){stop(" X must be numerical")}
	if(class(X)=="numeric"){X<-matrix(X,nrow=n)}
	if(class(X)=="matrix"){
		if(dim(X)[1]!=n){
			stop(" number of rows of X must be the same as the length of y")
		}
		pmax<-dim(X)[2]
		if(min(svd(cbind(X0,X)[1:(p0+pmax),])$d)<1e-12){
			stop(" cbind(X0,X) has collinear columns")
		}
	}
	colnames(X)<-paste("Test",c(1:pmax),sep="_")
	namesvec<-c(colnames(X0),colnames(X))
	if(!(type%in%c("NP","IP","ZS","HG"))){stop(" type is not in accepted values of 'NP','IP','ZS','HG'")}
	if(missing(prior.par)){
		if(type=="NP"){
			warning(" prior.par was missing and has been set to its default value of 1")
			prior.par<-1
		}else if(type=="HG"){
			warning(" prior.par was missing and has been set to its default value of c(2,1)")
			prior.par<-c(2,1)
		}else{prior.par<-1}
	}
	if(class(prior.par)!="numeric"){stop(" prior.par must be a positive numerical vector")}
	if(min(prior.par)<=0||min(prior.par)==Inf){stop(" prior.par must be a positive numerical vector")}
	if(type=="HG"&&length(prior.par)<2){
		warning(" type='HG' and length(prior.par)<2, values have been repeated")
		prior.par<-rep(prior.par,2)
	}
	if(type=="NP"&&length(prior.par)>1){
		warning("type='NP' and length(prior.par)>1, only first value will been used")
		prior.par<-prior.par[1]
	}
	if(!(EFF%in%c(FALSE,TRUE))){
		warning(" EFF was not TRUE or FALSE, it has been set to the default of TRUE")
		EFF<-TRUE
	}
	if(!(model.prior.type)%in%c("Beta-Binomial","Binomial","Uniform")){
		warning(" model.prior.type misspecified, set to default value of Beta-Binomial")
		model.prior.type<-"Beta-Binomial"
	}
	if(class(model.prior.par)!="numeric"){
		warning(" model.prior.par misspecified, set to default value")
		if(model.prior.type=="Beta-Binomial"){
			model.prior.par<-c(1,1)
		}else{model.prior.par<-1}
	}else if(model.prior.type=="Beta-Binomial" && length(model.prior.par)<2){
		warning(" model.prior.par misspecified, set to default value")
		model.prior.par<-c(2,1)
	}else if(model.prior.type=="Beta-Binomial"){
		model.prior.par<-abs(model.prior.par[1:2])
		if(min(model.prior.par)==0){
			warning(" model.prior.par misspecified, set to default value")
			model.prior.par<-c(2,1)
		}
	}else if(model.prior.type=="Binomial" && length(model.prior.par)<1){
		warning(" model.prior.par misspecified, set to default value")
		model.prior.par<-1
	}else if(model.prior.type=="Binomial"){
		model.prior.par<-abs(model.prior.par[1])
		if(model.prior.par==0 ||model.prior.par>=1){
			warning(" model.prior.par misspecified, set to default value")
			model.prior.par<-.5
		}
	}
	if(min(model.prior.par)==0){ind<-which(model.prior.par==0)}
	if(!(RWflag%in%c(FALSE,TRUE))){
		warning(" RWflag was not TRUE or FALSE, it has been set to the default of TRUE")
		RWflag<-TRUE
	}
	if(pmax>20&&RWflag==FALSE){
		message(" Random Walk is not selected and the number of covariates is large.  This may result in very slow enumeration of the model space.  Do you want to continue?")
		iwanttoleave<-2
		while(!(iwanttoleave%in%c("y","n"))){
			message("Yes (y) or No (n):")
			iwanttoleave<-scan(n=1,what=character(),quiet=TRUE)
			if(length(iwanttoleave)==0){iwanttoleave<-2}
		}
		if(iwanttoleave=="n"){stop(" try RWflag=TRUE")}
	} 
	if(!(RWpropflag%in%c(FALSE,TRUE))){
		warning(" RWpropflag was not TRUE or FALSE, it has been set to the default of FALSE")
		RWpropflag<-FALSE
	}
	if(!(inference%in%c("none","both","selected","averaged"))){
		return("Errpr: infernece is not in accepted values of 'none','both','selected','averaged'")
	}
	if(class(alpha)!="numeric" || alpha<=0 || alpha>=1){
		warning(" alpha was misspecified.  It has been set to its default value of 0.05")
		alpha<-0.05
	}
	if(class(shuf)!="numeric" || shuf<0 || shuf>1){
		warning(" shuf was misspecified.  It has been set to its default value of 0.2")
		shuf<-0.2
	}
	if(missing(start)){
		start<-paste(rep(0,pmax),collapse="")
	}else if(!(class(start)%in%c("numeric","character"))){
		start<-paste(rep(0,pmax),collapse="")
		warning(" start was misspecified, it has been set to the default value of ",start)
	}else if(class(start)=="numeric"){
		if(min(start)<0){
			start<-paste(rep(0,pmax),collapse="")
			warning("start was misspecified, it has been set to the default value of ",start)
		}
		if(max(start)>pmax){
			start<-paste(rep(0,pmax),collapse="")
			warning("start was misspecified, it has been set to the default value of ",start)
		}
		if(length(start)<pmax){
			start2<-rep(0,pmax)
			start2[unique(start[start>0])]<-1
			start<-paste(start2,collapse=TRUE)
		}
	}else{
		start.split<-as.integer(strsplit(start,split="")[[1]])
		if(length(start.split)!=pmax && min(start.split)==0 && max(start.split)==1){
			start<-paste(rep(0,pmax),collapse="")
			warning("start was misspecified, it has been set to the default value of ",start)
		}
	}
	if(class(N.draws)!="numeric"){
			N.draws<-10^4
			warning(" N.draws was misspecified, it has been set to the default value of 10^4")
	}else if((N.draws[1])<0){
			N.draws<-10^4
			warning(" N.draws was misspecified, it has been set to the default value of 10^4")
	}else{N.draws<-ceiling(N.draws[1])}
	if(class(m.considered)!="numeric"){
			m.considered<-100
			warning(" m.considered was misspecified, it has been set to the default value of 100")
	}else if(m.considered[1]<=0){
			m.considered<-100
			warning(" m.considered was misspecified, it has been set to the default value of 100")
	}else if(m.considered[1]>1){
		m.considered<-ceiling(m.considered[1])
	}

		
	
	
#############################################################################
## These are the parameters that are needed   																						   ##
#############################################################################

#############################################################################
#############################################################################


#############################################################################
## Setting function for log of model space prior																						   ##
## model.prior.type can be "Beta-Binomial" with hyperparameters (alpha,beta)=model.prior.par    ##
## default value is model.prior.par=c(1,1)                                                                                                     ##
## if model.prior.type is not "Beta-Binomial" then the model space prior is assumed to be uniform ##
#############################################################################
if(model.prior.type=="Beta-Binomial"){
	log.prior<-function(gamma){
		gamma<-as.integer(strsplit(gamma,split="")[[1]])
		s<-sum(gamma)
		l<-length(gamma)
		lbeta(s+model.prior.par[1],l-s+model.prior.par[2])-lbeta(model.prior.par[1],model.prior.par[2])
	}
}else if(model.prior.type=="Binomial"){
	log.prior<-function(gamma){
		gamma<-as.integer(strsplit(gamma,split="")[[1]])
		s<-sum(gamma)
		l<-length(gamma)
		s*log(model.prior.par[1])+(l-s)*log(1-model.prior.par[1])
	}	
}else{
	log.prior<-function(gamma){
		gamma<-as.integer(strsplit(gamma,split="")[[1]])
		-length(gamma)*log(2)	
	}
}
#############################################################################
#############################################################################

#############################################################################
## Setting function for log Bayes' factor																										   ##
## Penalization can come from one of four families					      															   ##
## type="IP" for intrinsic priors, "ZS" for Zellner-Siow prior, "HG" for hyper-g prior, "NP" for normal  ##
## type="HG" has an additional hyper-parameter (a,b)=prior.par        								                       ##
## default value is hg.prior.par=c(2,1)																											   ##
## type="NP" has an additional hyper-parameter w=prior.par															       ##
## default value is np.prior.par=1     																											   ##
## Additionally, each prior is scaled by n, the number of observations, providing consistency         ##
## The prior can also be chosen to scale by the minimal sample size by setting EFF=TRUE           ##
## EFF=FALSE is the usual g/n framework																							   	  ## 
## default value is EFF=TRUE          																											  ##
#############################################################################
if(EFF){
	qfun<-function(pA){(pA+1)*(pA>p0)}
}else{
	qfun<-function(pA){1}
}
if(type=="IP"){
	w.trans<-function(theta){
		(sin(pi/2*theta))^2
	}
	jacobian.factor<-function(theta){
		0
	}
}else if(type=="ZS"){
	w.trans<-function(theta){
		(sin(theta*pi/2)^(-2)-1)
	}
	jacobian.factor<-function(theta){
		w<-w.trans(theta)
		.5*log(pi/2)+log((1+w)^((sign(theta))))-w/2
	}
}else if(type=="HG"){
	w.trans<-function(theta){
		k<-ceiling(1/prior.par[1])
		prior.par[2]*(sin(theta*pi/2)^(-2*k)-1)
	}
	jacobian.factor<-function(theta){
		k<-ceiling(1/prior.par[1])
		v<-apply(matrix(sin(pi/2*theta)^2,length(0:(k-1)),length(theta),byrow=TRUE)^(0:(k-1)),2,sum)
		log(sin(theta*pi/2)^(k*prior.par[1]-1))-.5*log(v)+lgamma((prior.par[1]+1)/2)-lgamma(prior.par[1]/2)-lgamma(1/2)+log(k*pi)
	}
}else if(type=="NP"){
	w.trans<-function(theta){
		rep(prior.par[1],length(theta))
	}
	jacobian.factor<-function(theta){
		rep(0,length(theta))
	}
}

f1.BF<-function(theta,pA,Rsq){	
		multn1<-qfun(pA)/n
		multn2<-(n-pA)/2
		multRsq0<-1-Rsq
		multRsq1<-1/multRsq0
		multRsq2<-Rsq*multRsq1
		multp<-(pA-p0)/2
		w<-w.trans(theta)
		fw1<-1/(w*multn1)
		fw2<-((fw1)+1)
		fw3<-((fw1)+multRsq1)
		(-multn2*log(1+multRsq2/fw2)-log(fw3^(multp>0))*multp+jacobian.factor(theta))
}

f1.shrink<-function(theta,pA,Rsq){	
		w<-w.trans(theta)
		multn1<-qfun(pA)/n
		f1.BF(theta=theta,pA=pA,Rsq=Rsq) - log((1+w*multn1))
}
	
f2<-function(theta,pA,Rsq,val1,f1){
		exp(f1(theta,pA,Rsq)-val1)
}

iii<-seq(0,1,length.out=5)#max(c(ceiling(sqrt(n)),20)))
ind1<-iii[1:(length(iii)-1)]
ind2<-iii[2:length(iii)]		
	
log.BF.and.shrink<-function(pA,Rsq){
	add<- -.5*(n-p0)*log(1-Rsq)
	val1<-optimize(f1.BF,c(0,1),maximum=TRUE,tol=1E-12,pA=pA,Rsq=Rsq)$objective
	val2<-0
	for(i in 1:length(ind1)){
		val2<-val2+integrate(f2,ind1[i],ind2[i],abs.tol=1E-12,f1=f1.BF,pA=pA,Rsq=Rsq,val1=val1)$value
	}
	log.NC<-as.vector(log(val2)+val1)
	log.BF<-as.vector(log(val2)+val1+add)
	val1<-optimize(f1.shrink,c(0,1),maximum=TRUE,tol=1E-12,pA=pA,Rsq=Rsq)$objective
	val2<-0
	for(i in 1:length(ind1)){
		val2<-val2+integrate(f2,ind1[i],ind2[i],abs.tol=1E-12,f1=f1.shrink,pA=pA,Rsq=Rsq,val1=val1)$value
	}
	shrink<-as.vector(exp(log(val2)+val1-log.NC))
	return(list(logBF=log.BF,shrink=shrink,logNC=log.NC))
}

#############################################################################
#############################################################################

#############################################################################
## These are the functions that are needed to get model info  																   ##
#############################################################################
X0X0inv<-solve(crossprod(X0))
X0X<-X0X0inv%*%crossprod(X0,X)
Xorth<-X-(X0%*%X0X0inv)%*%crossprod(X0,X)

baseinfo<-function(){
	beta0hat<-X0X0inv%*%crossprod(X0,y)
	BaseRSS<-as.vector(crossprod(y)-crossprod(y,X0)%*%(beta0hat))
	sigmasq0<-as.vector(BaseRSS/(n-p0))
	logPriorA<-log.prior(paste(rep("0",pmax),collapse=""))
	nums<-cbind(c(diag(X0X0inv),rep(0,pmax)),c(0),c(beta0hat,rep(0,pmax)),c(0))
	rownames(nums)<-namesvec
	inf.ind<-1:p0
	return(list(gammaA=paste(rep("0",pmax),collapse=""),pA=p0,RsqA=0,logPriorA=logPriorA,logBFA=0,logNCA=0,shrinkA=1,nums=nums,inf.ind=inf.ind,sigmasq0=sigmasq0))
}

Base.Model<-baseinfo()
sigmasq0<-Base.Model$sigmasq0
gammaA.out<-Base.Model$gammaA
pA.out<-Base.Model$pA
RsqA.out<-Base.Model$RsqA
logPriorA.out<-Base.Model$logPriorA
logBFA.out<-Base.Model$logBFA
logNCA.out<-Base.Model$logNCA
shrinkA.out<-Base.Model$shrinkA
inf.ind.out<-list(Base.Model$inf.ind)
nums.out<-list(Base.Model$nums)
rm(Base.Model)

modelinfo<-function(gamma){
	A.ind<-which((strsplit(gamma,split="")[[1]]=="1"))
	inf.ind<-c(1:p0,A.ind+p0)
	size.A<-length(A.ind)
	XAorth<-matrix(Xorth[,A.ind],ncol=size.A)
	X0XA<-matrix(X0X[,A.ind],ncol=size.A)
	XA<-matrix(X[,A.ind],ncol=size.A)
	XAorthXAorthinv<-solve(crossprod(XAorth))
	yproj<-crossprod(XAorth,y)
	betaAhat<-XAorthXAorthinv%*%yproj
	ExtraRSS<-crossprod(yproj,betaAhat)
	RsqA<-as.vector(ExtraRSS/sigmasq0/(n-p0))
	logPriorA<-log.prior(gamma)
	beta0add<- -X0X0inv%*%crossprod(X0,XA)%*%betaAhat
	cc<-cccc<-rep(0,pmax)
	cc[A.ind]<-betaAhat
	logBFinfo<-log.BF.and.shrink(size.A+p0,RsqA)
	bb<-X0X0inv%*%X0XA
	Sigma0r<-bb%*%XAorthXAorthinv
	Sigma00add<-diag(tcrossprod(Sigma0r,bb))
	cccc[A.ind]<-diag(XAorthXAorthinv)
	nums<-cbind(c(nums.out[[1]][1:p0,1],rep(0,pmax)),c(Sigma00add,cccc),c(nums.out[[1]][1:p0,3],rep(0,pmax)),c(beta0add,cc))
	rownames(nums)<-namesvec
	return(list(gammaA=gamma,pA=size.A+p0,RsqA=RsqA,logPriorA=logPriorA,logBFA=logBFinfo[[1]],logNCA=logBFinfo[[3]],shrinkA=logBFinfo[[2]],nums=nums,inf.ind=inf.ind))
}


#############################################################################
#############################################################################

#############################################################################
## Function to enumerate models                                                   																   ##
#############################################################################
if(!RWflag){
	bit.increase<-function(x){
		if(prod(x)==TRUE){return(c(rep(FALSE,length(x)),TRUE))}
		index<-1
		flag<-1
		while(flag){
			if(x[index]==FALSE){x[index]<-TRUE;flag<-0}else{x[index]<-FALSE;index<-index+1}
		}
		return(x)
	}
	gamma<-matrix(0,2^pmax,pmax)
	for(i in 2:2^pmax){
		gamma[i,]<-bit.increase(gamma[i-1,])
	}
	gamma<-apply(gamma,1,paste,collapse="")
	gammaA.out<-c(gammaA.out,gamma[2:(2^pmax)])
	rm(gamma)
	for(i in 2:(2^pmax)){
		Model<-modelinfo(gammaA.out[i])
		pA.out[i]<-Model$pA
		RsqA.out[i]<-Model$RsqA
		logPriorA.out[i]<-Model$logPriorA
		logBFA.out[i]<-Model$logBFA
		logNCA.out[i]<-Model$logNCA
		shrinkA.out[i]<-Model$shrinkA
		inf.ind.out[[i]]<-Model$inf.ind
		nums.out[[i]]<-Model$nums
	}
	rm(Model)
	log.probs<-logBFA.out+logPriorA.out
	log.probs<-log.probs-max(log.probs)
	log.probs<-log.probs-log(sum(exp(log.probs)))
	freq.visit<-(exp(log.probs))
	visit.order<-gammaA.out
	names(freq.visit)<-names(log.probs)<-names(RsqA.out)<-names(logPriorA.out)<-names(pA.out)<-names(logBFA.out)<-names(nums.out)<-names(inf.ind.out)<-gammaA.out
}else{
#############################################################################
#############################################################################

#############################################################################
## Function to random walk models                                                   															   	##
#############################################################################
	log.probs<-logBFA.out+logPriorA.out
	if(missing(start)){
		mod.num.current<-1
	}else if(start!=gammaA.out[1]){
		i<-2
		gammaA.out[i]<-start
		Model<-modelinfo(gammaA.out[i])
		pA.out[i]<-Model$pA
		RsqA.out[i]<-Model$RsqA
		logPriorA.out[i]<-Model$logPriorA
		logBFA.out[i]<-Model$logBFA
		log.probs[i]<-logBFA.out[i]+logPriorA.out[i]
		logNCA.out[i]<-Model$logNCA
		shrinkA.out[i]<-Model$shrinkA
		inf.ind.out[[i]]<-Model$inf.ind
		nums.out[[i]]<-Model$nums
		rm(Model)
		mod.num.current<-2
	}else{
		mod.num.current<-1
	}
	gamma.current<-gammaA.out[mod.num.current]
	gamma.current.vec<-as.integer(strsplit(gamma.current,split="")[[1]])
	M<-length(gammaA.out)
	visit.order<-rep("",N.draws)
	freq.visit<-rep(0,M)
	mod.num.current<-M
	
for(draws in 1:N.draws){
	gamma.current<-gammaA.out[mod.num.current]
	gamma.current.vec<-as.integer(strsplit(gamma.current,split="")[[1]])
	to.shuf<-rbinom(1,1,shuf)
	if(to.shuf){
		if(model.prior.type=="Beta-Binomial"){
			p.draw<-rbeta(2,model.prior.par[1],model.prior.par[2])
		}else if(model.prior.type=="Beta-Binomial"){
			p.draw<-rep(model.prior.par[1],2)
		}else{
			p.draw<-rep(0.5,2)
		}
		vec.proposal<-rbinom(pmax,1,p.draw[1])
		gamma.proposals<-paste(vec.proposal,collapse="")
		gamma.comp<-function(gg){
			max(c(0,which(gammaA.out==gg)))
		}
		mod.num<-gamma.comp(gg=gamma.proposals)
		MM<-0
		if(mod.num==0){
			j<-1
			MM<-1
			i<-M+j
			gammaA.out[i]<-gamma.proposals
			Model<-modelinfo(gammaA.out[i])
			pA.out[i]<-Model$pA
			RsqA.out[i]<-Model$RsqA
			logPriorA.out[i]<-Model$logPriorA
			logBFA.out[i]<-Model$logBFA
			log.probs[i]<-logBFA.out[i]+logPriorA.out[i]
			logNCA.out[i]<-Model$logNCA
			shrinkA.out[i]<-Model$shrinkA
			inf.ind.out[[i]]<-Model$inf.ind
			nums.out[[i]]<-Model$nums
			freq.visit[i]<-0
			rm(Model)
			mod.num<-i

		}
		M<-M+MM
		log.prob.forwards<-(sum(vec.proposal*log(p.draw[1]))+sum((1-vec.proposal)*(1-p.draw[1])))
		log.prob.backwards<-(sum(gamma.current.vec*log(p.draw[2]))+sum((1-gamma.current.vec)*(1-p.draw[2])))
		log.prob.targ.forwards<-log.probs[mod.num]
		log.prob.targ.backwards<-log.probs[mod.num.current]
		if(log(runif(1))<(log.prob.targ.forwards-log.prob.forwards+log.prob.backwards-log.prob.targ.backwards)){
			gamma.current<-gamma.proposals
			mod.num.current<-mod.num
		}else{
			gamma.current<-gamma.current
		}
	}else{
		gamma.proposals<-as.integer(strsplit(gamma.current,split="")[[1]])
		gamma.proposals<-matrix(gamma.proposals,pmax,pmax,byrow=TRUE)
		diag(gamma.proposals)<-(diag(gamma.proposals)+1)%%2
		gamma.proposals<-apply(gamma.proposals,1,paste,collapse="")
		gamma.comp<-function(gg){
			max(c(0,which(gammaA.out==gg)))
		}
		Vgc<-Vectorize(gamma.comp,vectorize.args="gg")
		mod.num<-Vgc(gg=gamma.proposals)
		ind.new<-which(mod.num==0)
		MM<-length(ind.new)
		if(MM>0){
			for(j in 1:MM){
				i<-M+j
				gammaA.out[i]<-gamma.proposals[ind.new[j]]
				Model<-modelinfo(gammaA.out[i])
				pA.out[i]<-Model$pA
				RsqA.out[i]<-Model$RsqA
				logPriorA.out[i]<-Model$logPriorA
				logBFA.out[i]<-Model$logBFA
				log.probs[i]<-logBFA.out[i]+logPriorA.out[i]
				logNCA.out[i]<-Model$logNCA
				shrinkA.out[i]<-Model$shrinkA
				inf.ind.out[[i]]<-Model$inf.ind
				nums.out[[i]]<-Model$nums
				freq.visit[i]<-0
				rm(Model)
				mod.num[ind.new[j]]<-i
			}
		}
		M<-M+MM	
		log.prob.draw.forwards<-log.probs[mod.num]
		ml<-max(log.prob.draw.forwards)
		log.prob.draw.forwards<-log.prob.draw.forwards-ml
		log.prob.draw.forwards<-log.prob.draw.forwards-log(sum(exp(log.prob.draw.forwards)))
		log.prob.draw.forwards<-log(exp(log.prob.draw.forwards)/2+1/2/pmax)
		j.draw.forwards<-sample(1:pmax,1,prob=exp(log.prob.draw.forwards))
		mod.num.draw.forwards<-mod.num[j.draw.forwards]
		log.prob.draw.forwards<-log.prob.draw.forwards[j.draw.forwards]
		gamma.draw.forwards<-gammaA.out[mod.num.draw.forwards]
		log.prob.targ.forwards<-log.probs[mod.num.draw.forwards]
		gamma.proposals<-as.integer(strsplit(gamma.draw.forwards,split="")[[1]])
		gamma.proposals<-matrix(gamma.proposals,pmax,pmax,byrow=TRUE)
		diag(gamma.proposals)<-(diag(gamma.proposals)+1)%%2
		gamma.proposals<-apply(gamma.proposals,1,paste,collapse="")
		gamma.comp<-function(gg){
			max(c(0,which(gammaA.out==gg)))
		}
		Vgc<-Vectorize(gamma.comp,vectorize.args="gg")
		mod.num<-Vgc(gg=gamma.proposals)
		ind.new<-which(mod.num==0)
		MM<-length(ind.new)
		if(MM>0){
			for(j in 1:MM){
				i<-M+j
				gammaA.out[i]<-gamma.proposals[ind.new[j]]
				Model<-modelinfo(gammaA.out[i])
				pA.out[i]<-Model$pA
				RsqA.out[i]<-Model$RsqA
				logPriorA.out[i]<-Model$logPriorA
				logBFA.out[i]<-Model$logBFA
				log.probs[i]<-logBFA.out[i]+logPriorA.out[i]
				logNCA.out[i]<-Model$logNCA
				shrinkA.out[i]<-Model$shrinkA
				inf.ind.out[[i]]<-Model$inf.ind
				nums.out[[i]]<-Model$nums
				freq.visit[i]<-0
				rm(Model)
				mod.num[ind.new[j]]<-i
			}
		}
		M<-M+MM	
		log.prob.draw.backwards<-log.probs[mod.num]
		ml<-max(log.prob.draw.backwards)
		log.prob.draw.backwards<-log.prob.draw.backwards-ml
		log.prob.draw.backwards<-log.prob.draw.backwards-log(sum(exp(log.prob.draw.backwards)))
		log.prob.draw.backwards<-log(exp(log.prob.draw.backwards)/2+1/2/pmax)
		j.draw.backwards<-j.draw.forwards
		mod.num.draw.backwards<-mod.num[j.draw.backwards]
		log.prob.draw.backwards<-log.prob.draw.backwards[j.draw.backwards]
		gamma.draw.backwards<-gammaA.out[mod.num.draw.backwards]
		log.prob.targ.backwards<-log.probs[mod.num.draw.backwards]
		if(log(runif(1))<(log.prob.targ.forwards-log.prob.draw.forwards+log.prob.draw.backwards-log.prob.targ.backwards)){
			gamma.current<-gamma.draw.forwards
			mod.num.current<-mod.num.draw.forwards
		}else{
			gamma.current<-gamma.current
		}
	}
	freq.visit[mod.num.current]<-freq.visit[mod.num.current]+1
	visit.order[draws]<-gamma.current
}
log.probs<-log.probs-max(log.probs)
log.probs<-log.probs-log(sum(exp(log.probs)))
}
names(freq.visit)<-names(log.probs)<-names(RsqA.out)<-names(logPriorA.out)<-names(pA.out)<-names(logBFA.out)<-names(logNCA.out)<-names(shrinkA.out)<-names(nums.out)<-names(inf.ind.out)<-gammaA.out

rm(X0X0inv)
rm(X0X)
rm(Xorth)

if(RWpropflag){probs<-freq.visit/sum(freq.visit)}else{probs=exp(log.probs)}


#############################################################################
#############################################################################

#############################################################################
## Extractng Model space information ##
#############################################################################
mean.fun<-function(gg){
	nums.out[[gg]][,3]+shrinkA.out[gg]*nums.out[[gg]][,4]
}
Vmf<-Vectorize(mean.fun)
means<-Vmf(gammaA.out)
colnames(means)<-gammaA.out
scales.fun<-function(gg){
	sigmasq0*(1-RsqA.out[gg]*shrinkA.out[gg])*(nums.out[[gg]][,1]+shrinkA.out[gg]*nums.out[[gg]][,2])
}
Vsf<-Vectorize(scales.fun)
scales<-Vmf(gammaA.out)
colnames(scales)<-gammaA.out

MA.mean<-apply(t(means)*probs,2,sum)
SM.index<-which(probs==max(probs))
SM.mean<-means[,SM.index]
SM.gamma<-gammaA.out[SM.index]
SM.prob<-probs[SM.index]
gamma.mat<-matrix(as.integer(unlist(strsplit(gammaA.out,split=""))),nrow=length(gammaA.out),byrow=TRUE)
prob0<-c(rep(0,p0),apply((1-gamma.mat)*probs,2,sum))
names(prob0)<-namesvec

probs<-sort(probs,decreasing=TRUE,index.return=TRUE)
gammaA.out<-gammaA.out[probs$ix]
logBFA.out<-logBFA.out[probs$ix]
logNCA.out<-logNCA.out[probs$ix]
logPriorA.out<-logPriorA.out[probs$ix]
pA.out<-pA.out[probs$ix]
RsqA.out<-RsqA.out[probs$ix]
shrinkA.out<-shrinkA.out[probs$ix]
nums.out<-nums.out[probs$ix]
inf.ind.out<-inf.ind.out[probs$ix]
freq.visit<-freq.visit[probs$ix]
means<-means[,probs$ix]
scales<-scales[,probs$ix]
probs<-probs$x
out.1<-list(gamma=gammaA.out,p=pA.out,prob=probs,logBF=logBFA.out,logPrior=logPriorA.out,Rsquared=RsqA.out,prob0=prob0,freq.visit=freq.visit)

means.out<-cbind(MA.mean,SM.mean)
rownames(means)<-c(colnames(X0),colnames(X))
if(m.considered>1){m.considered<-min(c(m.considered,length(probs)))}else{m.considered<-min(which(cumsum(probs)>=m.considered))}

gammaA.out<-gammaA.out[1:m.considered]
logBFA.out<-logBFA.out[1:m.considered]
logPriorA.out<-logPriorA.out[1:m.considered]
logNCA.out<-logNCA.out[1:m.considered]
pA.out<-pA.out[1:m.considered]
RsqA.out<-RsqA.out[1:m.considered]
shrinkA.out<-shrinkA.out[1:m.considered]
nums.out<-nums.out[1:m.considered]
inf.ind.out<-inf.ind.out[1:m.considered]
freq.visit<-freq.visit[1:m.considered]
means<-means[,1:m.considered]
scales<-scales[,1:m.considered]
probs<-probs[1:m.considered]
gamma.mat<-matrix(as.integer(unlist(strsplit(gammaA.out,split=""))),nrow=length(gammaA.out),byrow=TRUE)
prob0<-c(rep(0,p0),apply((1-gamma.mat)*probs,2,sum))
names(prob0)<-namesvec
#############################################################################
#############################################################################

#############################################################################
## Functions for computing credible sets ##
#############################################################################

prob.fun.model<-function(gg,u,index){
if(abs(u)==Inf){return(log(.5*(1+sign(u))))}
	if(!(index%in%inf.ind.out[[gg]])){
		return(-Inf)
	}else{
		Rsq<-RsqA.out[gg]
		pA<-pA.out[gg]
		logNCA<-logNCA.out[gg]
		nums<-nums.out[[gg]][index,]
		f1.prob<-function(theta,pA,Rsq){
			w<-w.trans(theta)
			weight<-1/(1+qfun(pA)*w/n)
			sigmasqA<-sigmasq0*(1-weight*Rsq)
			pt((u-nums[3]-weight*nums[4])/sqrt(sigmasqA*(nums[1]+weight*nums[2])),df=n-p0,log.p=TRUE)+f1.BF(theta,pA=pA,Rsq=Rsq)-logNCA
		}
	}
	val1<-optimize(f1.prob,c(0,1),maximum=TRUE,tol=1E-12,pA=pA,Rsq=Rsq)$objective
	val2<-0
	for(i in 1:length(ind1)){
		val2<-val2+integrate(f2,ind1[i],ind2[i],abs.tol=1E-12,f1=f1.prob,pA=pA,Rsq=Rsq,val1=val1)$value
	}
	return(log(val2)+val1)
}

Vpfm<-Vectorize(prob.fun.model,vectorize.args="gg")
prob.fun.models<-function(u,index){	
	log(sum(exp(Vpfm(gammaA.out,u=u,index=index))*probs))
}
V.prob.fun.model<-Vectorize(prob.fun.model,vectorize.args="u")
V.prob.fun.models<-Vectorize(prob.fun.models,vectorize.args="u")


credsingleQuantile<-function(gg,alpha){
shrink<-shrinkA.out[gg]
Rsq<-RsqA.out[gg]
pA<-pA.out[gg]
sigmasqA<-sigmasq0*(1-shrink*Rsq)
numsmat<-nums.out[[gg]]
compprob<-function(index,q){
	nums<-numsmat[index,]
	ggfun<-function(v){
		u<-qt(v,n-p0)*sqrt(sigmasqA*(nums[1]+shrink*nums[2]))+nums[3]+shrink*nums[4]
		log(q)-V.prob.fun.model(u=u,gg=gg,index=index)
	}
	obj<-(uniroot(ggfun,c(0,1),tol=q*1e-20,maxiter=10000)$root)
	qt(obj,n-p0)*sqrt(sigmasqA*(nums[1]+shrink*nums[2]))+nums[3]+shrink*nums[4]
}
Vcp<-Vectorize(compprob,vectorize.args="index")
	lowers<-rep(0,p0+pmax)
	uppers<-rep(0,p0+pmax)
	lowers[inf.ind.out[[gg]]]<-Vcp(inf.ind.out[[gg]],q=alpha/2)
	uppers[inf.ind.out[[gg]]]<-Vcp(inf.ind.out[[gg]],q=1-alpha/2)
	out<-cbind(lowers,uppers)
return(out)
}

if(inference%in%c("both","selected")){
credsets.selected<-credsingleQuantile(SM.gamma,alpha)
rownames(credsets.selected)<-namesvec
}else{credsets.selected<-NULL}

credmixQuantile<-function(prob0,alpha){
get.cred.quantile<-function(index,alpha){
	prob0loc<-prob0[index]
	if(prob0loc>=(1-alpha)){return(list(intervals=matrix(0,2,1),probabilities=prob0loc))}
	quantile.approx<-function(v){
		qt((v)/(1-prob0loc),df=(n-p0))*sqrt(sum(probs*scales)+sum(probs*(means-sum(means*probs))^2))+sum(probs*means)
	}
	probl0<-exp(prob.fun.models(u=0,index=index))
	probg0<-1-prob0loc-probl0
	

	intervals<-NULL
	probabilities<-NULL
	
	if(prob0loc>=alpha){
		alphause<-alpha
		intervals<-matrix(rep(0,2),2,1)
		probabilities<-c(prob0loc)
	}else if((probl0>(alpha-prob0loc)/2) && (probg0>(alpha-prob0loc)/2)){
		alphause<-alpha
		intervals<-matrix(rep(0,2),2,1)
		probabilities<-c(prob0loc)
	}else{
		alphause<-alpha-prob0loc
	}

	f<-function(v,val){
		x<-quantile.approx(v)
		exp(V.prob.fun.models(u=x,index=index))-val
	}
	lb<-quantile.approx(uniroot(f,c(0,1-prob0loc),val=alphause/2,tol=1e-10)$root)
	ub<-quantile.approx(uniroot(f,c(0,1-prob0loc),val=1-prob0loc-alphause/2,tol=1e-10)$root)
	problb<-exp(prob.fun.models(lb,index))
	probub<-exp(prob.fun.models(ub,index))
	probincont<-probub-problb
	if(lb<0 && ub>0){
		intervals<-matrix(c(lb,ub),2,1)
		probabilities<-prob0loc+probincont
	}else{
		intervals<-cbind(intervals,matrix(c(lb,ub),2,1))
		probabilities<-c(probabilities,probincont)
	}
	intervals<-t(intervals)
	colnames(intervals)<-c("lowers","uppers")
	names(probabilities)<-NULL
	return(list(intervals=intervals,probabilities=probabilities))
}

out<-list(NA)
for(i in 1:(p0+pmax)){
	out[[i]]<-get.cred.quantile(i,alpha)
}
return(out)
}

if(inference%in%c("both","averaged")){
credsets.mixture<-credmixQuantile(prob0,alpha)
names(credsets.mixture)<-c(colnames(X0),colnames(X))
}else{credsets.mixture<-NULL}

MA.resids<-y-cbind(X0,X)%*%means.out[,1]
SM.resids<-y-cbind(X0,X)%*%means.out[,2]
resids<-cbind(MA.resids,SM.resids)
return(list(gamma=out.1$gamma,p=out.1$p,logPrior=out.1$logPrior,Rsquared=out.1$Rsquared,logBF=out.1$logBF,prob=out.1$prob,selected.model=SM.gamma,selected.model.prob=out.1$prob[1],prob0=out.1$prob0,means=means.out,credsets.selected=credsets.selected,credsets.averaged=credsets.mixture,resids=resids))
}