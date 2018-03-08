## Loading auxillary custom functions ####
trapz_simple<-function(x,y){
  nx=length(x)
  
  xmesh=x[(2:nx)]-x[(1:(nx-1))]
  ymesh=0.5*(y[(2:nx)]+y[(1:(nx-1))])
  
  sum(xmesh*ymesh)
}


# Michealis-Menten function
MM=function(x,K){
  x/K/(1+x/K)
}

# Binomial distribution functions with parameters
PYx=function(Y,x,N,g,H){
  dbinom(Y,N,g(x,H))
}

# generating output probability with exponential noise
PYxTexp=function(Y,xT,N,g,H,sd,signal_lambda,Nseq){
  
  lambda_rate=1/sd
  x=c(0,exp(seq(from=log(10^(-6)),to=-(log( 10^(-6)  ))/( lambda_rate),length.out=Nseq)))
  ppr=dexp(x,rate =  lambda_rate) 
  
  pl=lapply(Y,PYx,(xT+signal_lambda*x),N,g,H)
  
  Z=sapply(pl,function(x1) trapz_simple(x,x1*ppr) )
  #Z=Z/sum(Z)
  #Z#P=as.numeric(lapply(Z,"sum"))
}

# generating output probability with gamma noise
PYxTgamma=function(Y,xT,N,g,H,sd,signal_lambda,Nseq){
  cvar=sd^2
  k_rate=(2*cvar+1+sqrt(4*cvar+1))/(2*cvar)
  theta_rate=1/(k_rate-1)
  prob_tol=10^(-6)
  
  x=c(0,exp(seq(from=log(prob_tol),to=qgamma((1-prob_tol),shape=k_rate,scale=theta_rate),length.out=Nseq)))
  ppr=dgamma(x,shape=k_rate,scale=theta_rate) 
  
  pl=lapply(Y,PYx,(xT+signal_lambda*x),N,g,H)
  
  Z=sapply(pl,function(x1) trapz_simple(x,x1*ppr) )
  #Z=Z/sum(Z)
  #Z#P=as.numeric(lapply(Z,"sum"))
}


# generating output probability with lognormal noise
PYxTLN2=function(Y,xT,N,g,H,signal_lambda,sigma=1,mu=0,Nseq){
  
  x=seq(from=-10*sigma,to=10*sigma,length.out=Nseq)
  # x=xT+c(0,(seq(from=(qlnorm(0.0001, meanlog = mu+log(lambda), sdlog = sigma)),
  #                 to=(qlnorm(0.9999, meanlog = mu+log(lambda), sdlog = sigma)),
  #                 length.out=Nseq)))
  ppr=dnorm(x,mean=mu,sd=sigma) 
  
  pl=lapply(Y,PYx,(xT+signal_lambda*exp(x)),N,g,H)
  
  Z=sapply(pl,function(x1) trapz_simple(x,x1*ppr) )
  #Z=Z/sum(Z)
  #Z#P=as.numeric(lapply(Z,"sum"))
}



# 
# C=function(N){
#   0.5*log2(N/(2*exp(1)/pi))
# }
# C_Num=function(lambda,N,g,H,RL=-2,RR=4,sigma,mu){
#   gridl=0.01
#   x_T=10^seq(RL,RR,gridl) 
#   #print(c(sigma,mu))
#   SS=melt((lapply(x_T,FI_num, N,g,H,lambda,sigma,mu)))
#   S=t(matrix(SS[,1],3,length(x_T)))
#   Q=sqrt(S[,1])*x_T*log(10)
#   dmu=diff(S[,2])/diff(x_T)
#   QSN=sqrt(dmu^2/S[-1,3])*x_T[-1]*log(10)
#   R={}
#   list(x=x_T,C=log2(sum(Q)/sqrt(2*exp(1)*pi)*gridl),CSN=log2(sum(QSN)/sqrt(2*exp(1)*pi)*gridl),FIs=S[,1],FISNs=dmu^2/S[-1,3],mu=S[,2],var=S[,3])
# }
# 
# FI_num=function(xT,N,g,H,lambda,sigma,mu){
#   Y=0:(N)
#   E=as.numeric(lapply(Y,ExYxT,N,g,H,xT,lambda,sigma,mu))
#   E[E=='NaN']=10^(-6)
#   P=PYxT(Y,xT,N,g,H,lambda,sigma,mu)
#   P[P=='NaN']=0
#   m=sum(P*Y)
#   list(FI=sum(P*(E^2)),mu=m, var= sum(P*(Y-m)^2))
# }
# 
# 
# PxxT=function(x,xT,lambda,sigma,mu){
#   #dlnorm(x,xT+mu*lambda,sigma*lambda) 
#   dlnorm(x-xT,mu+log(lambda),sigma)
# }
# 
# PxYxT_UN=function(Y,x,N,g,H,xT,lambda,sigma,mu){
#   dbinom(Y,N,g(x,H))*dlnorm(x-xT,mu+log(lambda),sigma)
# }
# 
# 
# 
# Pdlnorm=function(Y,xT,N,g,H,lambda,sigma,mu){
#   x=xT+exp(log(lambda)+mu+seq(-sigma*10,sigma*10 ,length.out = Nseq))
#   x=x[x>xT]
#   ppr=dlnorm(x-xT,mu+log(lambda),sigma) 
#   list(xs=x,P=ppr)
# }
# 
# 
# PYxTnoNoise=function(Y,xT,N,g,H){
#   
#   Z=PYx(Y,xT,N,g,H)
#   #Z=Z/sum(Z)
#   #Z#P=as.numeric(lapply(Z,"sum"))
# }
# 
# 
# 
# PYxTLN=function(Y,xT,N,g,H,signal_lambda,sigma,mu=0,Nseq){
#   
#   sigma2=sigma
#   sigma=sqrt(log( (1+sqrt(1+4*sigma2^2))/2 ))
#   x=seq(from=-10*sigma,to=10*sigma,length.out=Nseq)
#   # x=xT+c(0,(seq(from=(qlnorm(0.0001, meanlog = mu+log(lambda), sdlog = sigma)),
#   #                 to=(qlnorm(0.9999, meanlog = mu+log(lambda), sdlog = sigma)),
#   #                 length.out=Nseq)))
#   ppr=dnorm(x,mean=mu,sd=sigma) 
#   
#   pl=lapply(Y,PYx,(xT+signal_lambda*exp(x)),N,g,H)
#   
#   Z=sapply(pl,function(x1) trapz_simple(x,x1*ppr) )
#   #Z=Z/sum(Z)
#   #Z#P=as.numeric(lapply(Z,"sum"))
# }
# 
# 
# 
# PYxT=function(Y,xT,N,g,H,lambda,sigma,mu){
#   x=xT+c(0,exp(seq(from=log(qlnorm(0.0001, meanlog = mu+log(lambda), sdlog = sigma)),
#                    to=log(qlnorm(0.9999, meanlog = mu+log(lambda), sdlog = sigma)),
#                    length.out=Nseq)))
#   # x=xT+c(0,(seq(from=(qlnorm(0.0001, meanlog = mu+log(lambda), sdlog = sigma)),
#   #                 to=(qlnorm(0.9999, meanlog = mu+log(lambda), sdlog = sigma)),
#   #                 length.out=Nseq)))
#   ppr=dlnorm(x-xT,meanlog=mu+log(lambda),sdlog=sigma) 
#   
#   
#   
#   pl=lapply(Y,PYx,x,N,g,H)
#   
#   Z=sapply(pl,function(x1) trapz_simple(x,x1*ppr) )
#   #Z=Z/sum(Z)
#   #Z#P=as.numeric(lapply(Z,"sum"))
# }
# 
# 
# PYxT2=function(Y,xT,N,g,H,lambda,sigma,mu){
#   x=xT*c(exp(seq(from=log(qlnorm(0.0001, meanlog = mu+log(lambda), sdlog = sigma)),
#                  to=log(qlnorm(0.9999, meanlog = mu+log(lambda), sdlog = sigma)),
#                  length.out=Nseq)))
#   # x=xT+c(0,(seq(from=(qlnorm(0.0001, meanlog = mu+log(lambda), sdlog = sigma)),
#   #                 to=(qlnorm(0.9999, meanlog = mu+log(lambda), sdlog = sigma)),
#   #                 length.out=Nseq)))
#   ppr=dlnorm(x/xT,meanlog=mu+log(lambda),sdlog=sigma) 
#   
#   
#   
#   pl=lapply(Y,PYx,x,N,g,H)
#   
#   Z=sapply(pl,function(x1) trapz_simple(x,x1*ppr) )
#   #Z=Z/sum(Z)
#   #Z#P=as.numeric(lapply(Z,"sum"))
# }
# 
# 
# PYxT3=function(Y,xT,N,g,H,lambda,sigma,mu){
#   x=xT+c(0,exp(seq(from=log(0.0001),
#                    to=log(0.0001+5*sigma),
#                    length.out=Nseq)))
#   # x=xT+c(0,(seq(from=(qlnorm(0.0001, meanlog = mu+log(lambda), sdlog = sigma)),
#   #                 to=(qlnorm(0.9999, meanlog = mu+log(lambda), sdlog = sigma)),
#   #                 length.out=Nseq)))
#   ppr=2*dnorm(x-xT,mean=0,sd=sigma) 
#   
#   
#   
#   pl=lapply(Y,PYx,x,N,g,H)
#   
#   Z=sapply(pl,function(x1) trapz_simple(x,x1*ppr) )
#   #Z=Z/sum(Z)
#   #Z#P=as.numeric(lapply(Z,"sum"))
# }
# 
# 
# ExYxT=function(Y,N,g,H,xT,lambda,sigma,mu){
#   L=qlnorm(0.001,meanlog=mu+log(lambda),sdlog=sigma)
#   R=qlnorm(1-0.001,meanlog=mu+log(lambda),sdlog=sigma)
#   x=xT+(seq((L),(R),length.out = Nseq))
#   
#   x=xT+exp(seq(log(L),log(R),length.out = Nseq))
#   dd=diff(x)
#   d=c(dd[1],dd)
#   p=PxYxT_UN(Y,x,N,g,H,xT,lambda,sigma,mu)*d
#   #  p
#   # print(Y)
#   # print("* ")
#   #  print(c(dbinom(Y,N,g(x,H)),Y,  g(x,H))) 
#   #  print(" ** ")
#   sum(  p* ( -(-sigma^2 + mu - log(x - xT) + log(lambda))/(sigma^2*(x - xT)) ))/sum(p)
# }
# 
# ExYxT=function(Y,N,g,H,xT,lambda,sigma,mu){
#   x=xT+exp(log(lambda)+mu+seq(-sigma*10,sigma*10 ,length.out = Nseq))
#   x=x[x>xT]
#   dd=diff(x)
#   d=c(dd[1],dd)
#   p=PxYxT_UN(Y,x,N,g,H,xT,lambda,sigma,mu)*d
#   #  p
#   # print(Y)
#   # print("* ")
#   #  print(c(dbinom(Y,N,g(x,H)),Y,  g(x,H))) 
#   #  print(" ** ")
#   sum(  p* ( -(-sigma^2 + mu - log(x - xT) + log(lambda))/(sigma^2*(x - xT)) ))/sum(p)
# }
# # 1   2   3   4   5  10  50 100 200 300
# 
# ##
# 
# binombetaVar<-function(n,a,b){
#   VarB=(n*a*b*(a+b+n))/((a+b)^2*(a+b+1))
#   VarB
# }
# 
# binombetaVarDerA<-function(n,a,b){
#   VarBderA=(n*b)*( (2*a+b+n)*(a+b)*(a+b+1)  - a*(a+b+n)*(3*a+3*b+2)  )/( ((a+b)^3)*((a+b+1)^2) ) 
#   VarBderA
# }