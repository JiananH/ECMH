#new strategy with bivariate normal case
pvec=c(0,0.05,0.125,0.25,0.375,0.5,0.7,0.9)
cvec=c(0.05,0.1,0.5,1,1.5,2,2.5,3,3.5)
itr=1000
N=1000

arate=function(x,c,mu,sigma,p){
    (1-p*(pnorm(x+c,mu,sigma)-pnorm(x-c,mu,sigma)))/(1-pnorm(x+c,mu,sigma)+pnorm(x-c,mu,sigma))
}

for (i in 1:length(pvec)){
	esjd.mcmhc=matrix(0,itr,length(cvec))
	mse.mcmhc=matrix(0,itr,length(cvec))
	draw=matrix(0,itr,length(cvec))
	acceptrate=matrix(0,itr,length(cvec))
	p=pvec[i]
	for (l in 1:itr){
  	  	for (k in 1:length(cvec)){
   	     	x1=x2=double(N)
   	     	x1[1]=x2[1]=0   
    	    c=cvec[k]
    	    draws=double(itr)
        	acceptance=update=matrix(0,itr,2)
        	for (i in 2:N){
           	 	judge=runif(1)
            	#update x1 first
            	if (judge>=0.5){
                	judge1=runif(1)
                #draw samples from the original proposal density
                if (judge1<p){
                	draws[i]=1
                	temp=rnorm(1,x2[i-1],1)
                	alpha1=min(1,arate(temp,c,x2[i-1],1,p)/arate(x1[i-1],c,x2[i-1],1,p))
                    judge2=runif(1)
                    x1[i]=ifelse(judge2<=alpha1,temp,x1[i-1])
                    accept=ifelse(judge2<=alpha1,1,0)
                #draw samples from the truncated proposal density
                } else{
                	repeat{
                   	draws[i]=draws[i]+1
                    temp=rnorm(1,x2[i-1],1)                  
                    if (abs(temp-x1[i-1])>c){break}
                }
					alpha1=min(1,arate(temp,c,x2[i-1],1,p)/arate(x1[i-1],c,x2[i-1],1,p))
                    judge2=runif(1)
                    x1[i]=ifelse(judge2<=alpha1,temp,x1[i-1])
                    accept=ifelse(judge2<=alpha1,1,0)
                }                	
                update[i,]=c(1,0)
                acceptance[i,]=c(accept,0)
                x2[i]=x2[i-1]    
                #update x2 first    
            	} else {
                	judge1=runif(1)
                	#draw samples from the original proposal density
               		if (judge1<p){
                		draws[i]=1
                		temp=rnorm(1,x1[i-1]/2,sqrt(1/2))
                		alpha2=min(1,arate(temp,c*sqrt(1/2),x1[i-1]/2,sqrt(1/2),p)/arate(x2[i-1],c*sqrt(1/2),x1[i-1]/2,sqrt(1/2),p))
                    	judge2=runif(1)
                    	x2[i]=ifelse(judge2<=alpha2,temp,x2[i-1])
                    	accept=ifelse(judge2<=alpha2,1,0)
                    #draw samples from the truncated proposal density
                	}else{
                		repeat{
		 	        	draws[i]=draws[i]+1
                    	temp=rnorm(1,x1[i-1]/2,sqrt(1/2))                 
                    	if (abs(temp-x2[i-1])>c*sqrt(1/2)){break}
                		}
                		alpha2=min(1,arate(temp,c*sqrt(1/2),x1[i-1]/2,sqrt(1/2),p)/arate(x2[i-1],c*sqrt(1/2),x1[i-1]/2,sqrt(1/2),p))
                    	judge2=runif(1)
                    	x2[i]=ifelse(judge2<=alpha2,temp,x2[i-1])
                    	accept=ifelse(judge2<=alpha2,1,0)
                	}
                	x1[i]=x1[i-1] 
                	update[i,]=c(0,1)
                	acceptance[i,]=c(0,accept)                 
            	}   
        	}
        	draw[l,k]=mean(draws)
        	updates = c(sum(update[,1]), sum(update[,2]))/N
        	acceptrate[l,k]=mean(c(sum(acceptance[,1]), sum(acceptance[,2]))/(updates*N))
			#evaluate the jump size and draws per iteration
			esjd.mcmhc[l,k]=(sum(diff(x1)^2+diff(x2)^2))/(N-1)
			mse.mcmhc[l,k]=mean(x1)^2
    	}
	}
	write.table(acceptrate,paste(p,"acceptrate.txt",sep=''))
	write.table(draw,paste(p,"draw.txt",sep=''))
	write.table(esjd.mcmhc,paste(p,"esjd.txt",sep=''))
	write.table(mse.mcmhc,paste(p,"mse.txt",sep=''))
}





