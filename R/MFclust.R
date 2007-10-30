MFclust.compute<-function(data, minG, maxG, nchain, thres, iter.max, my.alpha , ...){
	if(is.null(nchain)) nchain = 5
	if(is.null(thres)) thres = 0.5
	if(is.null(iter.max)) iter.max=10
	data=as.matrix(data)
	out=list();count=1
	for(nclust in minG:maxG){
		if(is.null(my.alpha)) alpha=rep(1,nclust) else alpha=my.alpha[nclust,]
		likeli=NULL;likeli[1]=10^20     
		for(chain in 1:nchain){
			clust=kmeans(data,nclust)$cluster
			if(check.start(clust)==TRUE){
				initial=compute.M.all(data, weight=NULL, clust)
				EM=em.clust(data,clust,mu=initial$mu,zeta=initial$zeta,varht=initial$varht,thres=thres,iter.max=iter.max,alpha=alpha)
				current.likeli=EM$loglikelihood
				if(current.likeli<min(likeli)) my.EM=EM
				likeli[chain]=current.likeli
				rm(EM)
				rm(initial)
			}else{
				print("warning: there is an empty cluster,reduce the number of cluster") 
			}			
		}
	out[[count]]=my.EM
	count=count+1
	}
	return(out=out)	
}

check.start<-function(clust){
	my.level<-unique(clust)
	my.check=TRUE
	for(k in 1:length(my.level)){
		if(sum(clust==my.level[k])<3) my.check=FALSE
	}
	if(length(my.level)!=max(clust)) my.check=FALSE	
	return(my.check)
}

MFclust<-function(data,minG,maxG,nchain=NULL,thres=NULL,iter.max=NULL,my.alpha=NULL,...){
	mf<-match.call(expand.dots=FALSE)
	m <- match(c("data", "minG", "maxG"), names(mf), 0)
	n=NROW(data)
	m=NCOL(data)
	my.obj=MFclust.compute(data,minG=minG,maxG=maxG,nchain=nchain,thres=thres,iter.max=iter.max,my.alpha=my.alpha)
	my.bic=NULL;count=maxG-minG+1
	for(k in 1:count){my.bic[k]=my.obj[[k]]$BIC}
	out.k<-which(my.bic==min(my.bic))
	OUTPUT=my.obj[[out.k]]
	row.names(OUTPUT$mu)<-paste("clust",1:(minG+out.k-1),sep="",collapse=NULL)
	colnames(OUTPUT$mu)<-paste("t",1:m,sep="",collapse=NULL)
	my.out=list(BIC=my.bic, nclust=minG+out.k-1, clust=OUTPUT$clust, clust.center=OUTPUT$mu)
	my.out
}


Estep.tik=function(xx.i,old){
   temp=ind.p=NULL
   m=length(xx.i)
   nclust=NROW(old$mu)	
   for(k in 1:nclust){
           clust.sigma=matrix(nrow=m,ncol=m)
           tempvar=10^old$zeta[k]*old$varht[k]
           rho=tempvar/(tempvar+old$varht[k])
           clust.sigma=(matrix(rho,m,m)+diag(rep((1-rho),m)))*(tempvar+old$varht[k]) 
           temp[k]=old$pk[k]*dmvnorm(xx.i,old$mu[k,],clust.sigma)
   }  		 
   t.ik=temp*1000/sum(temp*1000)
   loglike=log(sum(temp))	
   ind.p=which(t.ik==max(t.ik))[1]
   return(list(clust=ind.p,t.ik=t.ik,loglike=loglike))
}

Estep.tk<-function(x,old,nclust){
	if(is.matrix(x)){
		n=NROW(x)
	}else{
	stop("x is not a matrix in computing likelihood")
	}
	tk=matrix(nrow=n,ncol=nclust)
	clust=loglike=NULL
	for(i in 1:n){
         	temp=NULL
         	temp=Estep.tik(xx.i=x[i,],old)
         	tk[i,]=temp$t.ik
		clust[i]=temp$clust
		loglike[i]=temp$loglike
         	rm(temp)
         }	
	rm(old)
	loglikelihood=-sum(loglike)
	return(list(clust=clust,tk=tk,loglikelihood=loglikelihood))
}


compute.reject<-function(tk,k,thres){
   	n=NROW(tk)
	em.select=NULL  
      em.select=which(tk[,k]>=thres)
      label1=NULL
      label1=which(tk[,k]<thres)
      my.remain=NULL
      my.remain=tk[label1,k]
      sem.select=NULL
      sem.select=as.vector(tapply(my.remain,1:length(my.remain),
           		function(x,p=thres){sample(c(1,0),1,prob=c(x/p,(1-x/p)))}))
      label2=NULL
      label2=label1[which(sem.select==1)]
	return(list(em.select=em.select,sem.select=label2))
}


compute.weight<-function(tk,thres){
	nclust=NCOL(tk)
	n=NROW(tk)
	weight=matrix(0,nrow=n,ncol=nclust)
	my.label=list()
 	for(k in 1:nclust){
		 label=compute.reject(tk,k,thres)
		 weight[label$em.select,k]=tk[label$em.select,k]
		 weight[label$sem.select,k]=rep(thres,length(label$sem.select))
		 my.label[[k]]=sort(c(label$em.select,label$sem.select))
	}
	return(list(weight=weight,label=my.label))
}



compute.M.pk<-function(tk,nclust,alpha){
        n=NROW(tk)
        pk=NULL
        temp=NULL
        temp=apply(tk,2,sum)
        pk=(temp+alpha)/(sum(temp)+sum(alpha))
        return(pk)
}



compute.center<-function(x,weight){
	ni=NROW(x)
	m=NCOL(x)
	tt=1:m
	if(is.null(weight)){		
		temp.data=data.frame(gi=as.numeric(as.vector(t(x))),tm=rep(tt,ni),my.geneid=as.factor(sort(rep(1:ni,m))))	
	      my.ran=mkran(~1|my.geneid,data=temp.data)
		temp=ssanova(gi~tm, random=my.ran,id.basis=tt,alpha=1,data=temp.data)
	}else{
		my.weight=weight%x%rep(1,m)		 	
      	temp.data=data.frame(gi=as.numeric(as.vector(t(x))),tm=rep(tt,ni), my.weight=as.numeric(my.weight), my.geneid=as.factor(sort(rep(1:ni,m)))) 
  		my.ran=mkran(~1|my.geneid,data=temp.data)
		temp=ssanova(gi~tm, random=my.ran,weights=my.weight,id.basis=1:m,alpha=1,data=temp.data)
	}			
      temp.pred=predict(temp,data.frame(tm=tt))
      mu=temp.pred
	zeta=temp$zeta
	varht=temp$varht
    	trc=temp$trc
	return(list(mu=mu,zeta=zeta,varht=varht,trc=trc))
}



compute.M.all<-function(x,weight,my.label){
	if(is.null(weight)){
		 nclust=length(unique(my.label))
	}else{
		 nclust=NCOL(weight)
	}
	m=NCOL(x)
	mu=matrix(nrow=nclust,ncol=m)
	zeta=varht=trc=NULL
	for(k in 1:nclust){
		if(!is.null(weight))label=my.label[[k]] else label=which(my.label==k)
		gi=x[label,]
		if(is.null(weight)){
			temp=compute.center(x=gi,weight=NULL)
		}else{
			my.weight=weight[label,k]	
                  temp=compute.center(x=gi,weight=my.weight)
		}
		mu[k,]=temp$mu
		zeta[k]=temp$zeta
		varht[k]=temp$varht
		trc[k]=temp$trc
		rm(temp)		
	}
	return(list(mu=mu,zeta=zeta,varht=varht,trc=trc))
}

em.bic<-function(likelihood,trc,n,m){
	bic=2*likelihood+sum(trc)*log(n*m)
}



em.clust=function(x,clust,mu,zeta,varht,thres,iter.max,alpha){
  loglikelihood=bic=NULL;loglikelihood[1]=10^20
  n=NROW(x);m=NCOL(x);nclust=max(clust)	
  pk=NULL
  for(k in 1:nclust){		  
  	pk[k]=sum(clust==k)/n
  }	
  old=list(pk=pk,mu=mu,zeta=zeta,varht=varht) 		  
  for(iter in 1:iter.max){ 
	my.tk=Estep.tk(x=x,old=old,nclust=nclust)
	tk=my.tk$tk
      like.current=my.tk$loglikelihood	
	rc=compute.weight(tk,thres)			
	my.weight=rc$weight					
	label=rc$label
	rm(rc)
	pk=compute.M.pk(my.weight,nclust,alpha)
	obj=compute.M.all(x,weight=my.weight,my.label=label)
	if(like.current<min(loglikelihood)){
		 out.clust=my.tk$clust
		 out.mu=obj$mu
		 out.trc=obj$trc
	}
	mu=obj$mu;zeta=obj$zeta;varht=obj$varht
	rm(obj)
	loglikelihood[iter]=like.current
	rm(old)
	old=list(pk=pk,mu=mu,zeta=zeta,varht=varht) 
  }
  rm(old) 	
  BIC=em.bic(min(loglikelihood),out.trc,n,m)
  return(list(clust=out.clust,mu=out.mu,BIC=BIC,loglikelihood=min(loglikelihood)))
}
