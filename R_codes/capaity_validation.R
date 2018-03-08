# Custom code accompanying manuscript NCOMMS-18-04749:
# "An information-theoretic framework for deciphering pleiotropic and noisy biochemical signaling"

# Code I: Validation of capacity approximation by using Jeffrey's prior - R scripts
# R-script for calculating small noise approximation and C_N and C_A from Main Paper Results section
# for the model of ligand-receptor binding with false signal

## Set your workings directory ####
#setwd("~/Documents/Science/Paper_codes/R_codes/")
setwd("F:\\TJ\\Simulations\\Paper_codes\\R_codes")
dir.create("output/pyx_matr/",recursive = TRUE)

## Libraries needed for calculations ####
rm(list=ls())
library('latex2exp')
library(reshape2)
library(ggplot2)
library(ggthemes)
library(gridExtra)
source("aux_functions_capacity.R")


## Setting up variables and costants ####
nys=c(10,50,100,500) # Setting up number of molecules
isds=c(0.01,0.1,0.5,1,5,10,25,50) # Setting up standard variations
ilambdas=c(1,0.5,0.1) # Setting up lambdas
driv_step=0.0001 # shift parameter for derivative calculation
is=200  # Setting up density of signal
xTs=10^seq(-8,8,length.out = is)  # Setting up range of signal
Nseq=2000 # Setting up number sample size used to approximate single marginal distribution

outputList=list()    # Large object to store capacity values
outputListOpt=list() # Large object to store optimal distributions


## Calculations - Exponential distribution ####
outputList[["exp"]]=list()
outputListOpt[["exp"]]=list()
for (ny in nys){
  outputList[["exp"]][[as.character(ny)]]=list()
  outputListOpt[["exp"]][[as.character(ny)]]=list()
  for (isd in isds){
    outputList[["exp"]][[as.character(ny)]][[as.character(isd)]]=list()
    outputListOpt[["exp"]][[as.character(ny)]][[as.character(isd)]]=list()
    for (ilambda in ilambdas){
      
      FI=c()
      FIlist=list()
      FISN=c() #small noise approximation
      Plist=list()
      k=0
      for (ix in xTs){ #loop over all signals
        k=k+1
        if (k%%20 == 0){print(k)}
        FIix=c()
        Pix=c()
        # Probabilities at subsequent output states
        pf=c() # 'forward' signal shift (for derivative calculation)
        pc=c() # no signal shift (for derivative calculation)
        pb=c() # 'backward' signal shift (for derivative calculation)
        for (iy in 0:ny){ #loop over all output states
          pf[iy+1]=PYxTexp(iy,(ix+ix*driv_step),N=ny,g=MM,H=1,sd=isd,signal_lambda=ilambda,Nseq = Nseq)
          pc[iy+1]=PYxTexp(iy,ix,N=ny,g=MM,H=1,sd=isd,signal_lambda=ilambda,Nseq = Nseq)
          pb[iy+1]=PYxTexp(iy,(ix-ix*driv_step),N=ny,g=MM,H=1,sd=isd,signal_lambda=ilambda,Nseq = Nseq)
        }
        
        #normalising probabilities to 1 (for corrections of any numerical discrepancies)
        pf=pf/sum(pf) 
        pc=pc/sum(pc)
        pb=pb/sum(pb)
        
        #Fisher information at output states
        for (iy in 0:ny){
          Pix[iy+1]=pc[iy+1]
          FIix[iy+1]=(1/(ix*driv_step)^2)*(
            log(pf[iy+1]^pc[iy+1])-
              2*log(pc[iy+1]^pc[iy+1])+
              log(pb[iy+1]^pc[iy+1])
          )
        }
        
        dermu_x=(sum((0:ny)*pf)-sum((0:ny)*pb))/(2*ix*driv_step)
        var_x=sum(((0:ny)^2)*pc)-(sum((0:ny)*pc))^2
        
        Plist[[k]]=Pix
        FIlist[[k]]=FIix
        FISN[k]= (dermu_x^2)/var_x
        FI[k]=max(-sum(FIix[is.finite(FIix)]),0)
      }
      
      #Probability matrix (rows - signals, columns - output states)
      Pmatr=do.call(rbind,Plist)
      
      FI[!is.finite(FI)]=0
      CCC=trapz_simple(xTs,sqrt(FI)) 
      CCCSN=trapz_simple(xTs,sqrt(FISN)) # Channel capacity via small noise approximation
      JP=sqrt(FI)/CCC # Optimal probability via Jeffrey's prior
      JPSN=sqrt(FISN)/CCCSN # Optimal probability via small noise approximation
      vecFI=1:is
      
      # Channel capacity via Jeffrey's prior
      CCC=trapz_simple(xTs[vecFI],sqrt(FI)[vecFI])
      CCCC=log( (1/(2*pi*exp(1))^0.5) * CCC)
      CCCCb=log2(exp(CCCC)) 
      
      
      ProbX =  apply(Pmatr[vecFI,],2,function(x) trapz_simple(xTs[vecFI],x*((JP[vecFI])) ) )
      HX = -sum(log(ProbX^ProbX));
      HXS = trapz_simple(xTs[vecFI],((JP[vecFI]))*apply(Pmatr[vecFI,],1,function(x) (-sum(log(x^x))) )  )
      
      output_m=(apply(Pmatr,1,function(x) sum(0:ny*x)))
      output_sd=sqrt(apply(Pmatr,1,function(x) sum((0:ny)^2*x)-(sum(0:ny*x))^2 ))
      
      # Channel capacity via small noise approximation (two approaches for numerical error check)
      Cbialek=trapz_simple(output_m,1/output_sd)
      Cbialek2=trapz_simple(xTs,sqrt(FISN))
      
      outputList[["exp"]][[as.character(ny)]][[as.character(isd)]][[as.character(ilambda)]]=data.frame(jpexact=log2(exp(HX-HXS)),
                                                                                                       jpapp=CCCCb,
                                                                                                       sn1=log2( Cbialek/sqrt(2*pi*exp(1)) ),
                                                                                                       sn2=log2( Cbialek2/sqrt(2*pi*exp(1)) ))
      
      outputListOpt[["exp"]][[as.character(ny)]][[as.character(isd)]][[as.character(ilambda)]]$JP=data.frame(s=xTs,popt=JP)
      outputListOpt[["exp"]][[as.character(ny)]][[as.character(isd)]][[as.character(ilambda)]]$SN=data.frame(s=xTs,popt=JPSN)
      
      # c(log2(exp(HX-HXS)),CCCCb,log2( Cbialek/sqrt(2*pi*exp(1)) ),log2( Cbialek2/sqrt(2*pi*exp(1)) ))
      
      
      # Saving probability matrix of the single case in txt file
      write.table(Pmatr,paste("output/pyx_matr/ver7_exp_Ps",isd,"_sig",is,"_n",ny,"_lambda",ilambda,".csv",sep=""),
                  col.names = FALSE,
                  row.names = FALSE,
                  sep=",")
      
    }
  }
}



## Calculations - Gamma distribution ####
outputList[["gamma"]]=list()
outputListOpt[["gamma"]]=list()
for (ny in nys){
  outputList[["gamma"]][[as.character(ny)]]=list()
  outputListOpt[["gamma"]][[as.character(ny)]]=list()
  for (isd in isds){
    outputList[["gamma"]][[as.character(ny)]][[as.character(isd)]]=list()
    outputListOpt[["gamma"]][[as.character(ny)]][[as.character(isd)]]=list()
    for (ilambda in ilambdas){
      
      # ny=500
      # isd=10
      # 
      #xTs=seq(10^-8,10^8,length.out = is)
      #xTs=c(10^seq(-4,0,length.out = is/4),seq(from=1.01,to=10^2,length.out=is/4),
      #      seq(from=10^2+1,to=10^4,length.out=is/4),seq(from=10^4+1,to=10^6,length.out=is/4))
      
      FI=c()
      FIlist=list()
      FISN=c()
      Plist=list()
      k=0
      for (ix in xTs){
        k=k+1
        if (k%%20 == 0){print(k)}
        FIix=c()
        Pix=c()
        pf=c()
        pc=c()
        pb=c()
        for (iy in 0:ny){
          pf[iy+1]=PYxTgamma(iy,(ix+ix*driv_step),N=ny,g=MM,H=1,sd=isd,signal_lambda=ilambda,Nseq = Nseq)
          pc[iy+1]=PYxTgamma(iy,ix,N=ny,g=MM,H=1,sd=isd,signal_lambda=ilambda,Nseq = Nseq)
          pb[iy+1]=PYxTgamma(iy,(ix-ix*driv_step),N=ny,g=MM,H=1,sd=isd,signal_lambda=ilambda,Nseq = Nseq)
        }
        
        pf=pf/sum(pf)
        pc=pc/sum(pc)
        pb=pb/sum(pb)
        
        for (iy in 0:ny){
          Pix[iy+1]=pc[iy+1]
          FIix[iy+1]=(1/(ix*driv_step)^2)*(
            log(pf[iy+1]^pc[iy+1])-
              2*log(pc[iy+1]^pc[iy+1])+
              log(pb[iy+1]^pc[iy+1])
          )
        }
        
        dermu_x=(sum((0:ny)*pf)-sum((0:ny)*pb))/(2*ix*driv_step)
        var_x=sum(((0:ny)^2)*pc)-(sum((0:ny)*pc))^2
        
        Plist[[k]]=Pix
        FIlist[[k]]=FIix
        FISN[k]= (dermu_x^2)/var_x
        FI[k]=max(-sum(FIix[is.finite(FIix)]),0)
      }
      
      
      Pmatr=do.call(rbind,Plist)
      
      FI[!is.finite(FI)]=0
      CCC=trapz_simple(xTs,sqrt(FI))
      CCCSN=trapz_simple(xTs,sqrt(FISN))
      JP=sqrt(FI)/CCC
      JPSN=sqrt(FISN)/CCCSN
      vecFI=1:is
      
      CCC=trapz_simple(xTs[vecFI],sqrt(FI)[vecFI])
      CCCC=log( (1/(2*pi*exp(1))^0.5) * CCC)
      CCCCb=log2(exp(CCCC))
      
      
      ProbX =  apply(Pmatr[vecFI,],2,function(x) trapz_simple(xTs[vecFI],x*((JP[vecFI])) ) )
      HX = -sum(log(ProbX^ProbX));
      HXS = trapz_simple(xTs[vecFI],((JP[vecFI]))*apply(Pmatr[vecFI,],1,function(x) (-sum(log(x^x))) )  )
      
      output_m=(apply(Pmatr,1,function(x) sum(0:ny*x)))
      output_sd=sqrt(apply(Pmatr,1,function(x) sum((0:ny)^2*x)-(sum(0:ny*x))^2 ))
      Cbialek=trapz_simple(output_m,1/output_sd)
      Cbialek2=trapz_simple(xTs,sqrt(FISN))
      
      outputList[["gamma"]][[as.character(ny)]][[as.character(isd)]][[as.character(ilambda)]]=data.frame(jpexact=log2(exp(HX-HXS)),
                                                                                                         jpapp=CCCCb,
                                                                                                         sn1=log2( Cbialek/sqrt(2*pi*exp(1)) ),
                                                                                                         sn2=log2( Cbialek2/sqrt(2*pi*exp(1)) ))
      
      outputListOpt[["gamma"]][[as.character(ny)]][[as.character(isd)]][[as.character(ilambda)]]$JP=data.frame(s=xTs,popt=JP)
      outputListOpt[["gamma"]][[as.character(ny)]][[as.character(isd)]][[as.character(ilambda)]]$SN=data.frame(s=xTs,popt=JPSN)
      
      
      # c(log2(exp(HX-HXS)),CCCCb,log2( Cbialek/sqrt(2*pi*exp(1)) ),log2( Cbialek2/sqrt(2*pi*exp(1)) ))
      
      
      
      write.table(Pmatr,paste("output/pyx_matr/ver7_gamma_Ps",isd,"_sig",is,"_n",ny,"_lambda",ilambda,".csv",sep=""),
                  col.names = FALSE,
                  row.names = FALSE,
                  sep=",")
      
    }
  }
}


# # Lognormal distribution
# outputList[["lognormal"]]=list()
# outputListOpt[["lognormal"]]=list()
# for (ny in nys){
#   outputList[["lognormal"]][[as.character(ny)]]=list()
#   outputListOpt[["lognormal"]][[as.character(ny)]]=list()
#   for (isd in isds){
#     outputList[["lognormal"]][[as.character(ny)]][[as.character(isd)]]=list()
#     outputListOpt[["lognormal"]][[as.character(ny)]][[as.character(isd)]]=list()
#     for (ilambda in ilambdas){
#       
#       # ny=500
#       # isd=10
#       # 
#       #xTs=seq(10^-8,10^8,length.out = is)
#       #xTs=c(10^seq(-4,0,length.out = is/4),seq(from=1.01,to=10^2,length.out=is/4),
#       #      seq(from=10^2+1,to=10^4,length.out=is/4),seq(from=10^4+1,to=10^6,length.out=is/4))
#       
#       FI=c()
#       FIlist=list()
#       FISN=c()
#       Plist=list()
#       k=0
#       for (ix in xTs){
#         k=k+1
#         if (k%%20 == 0){print(k)}
#         FIix=c()
#         Pix=c()
#         pf=c()
#         pc=c()
#         pb=c()
#         for (iy in 0:ny){
#           pf[iy+1]=PYxTLN(iy,(ix+ix*driv_step),N=ny,g=MM,H=1,sigma=isd,signal_lambda=ilambda,Nseq = Nseq,mu=0)
#           pc[iy+1]=PYxTLN(iy,ix,N=ny,g=MM,H=1,sigma=isd,signal_lambda=ilambda,Nseq = Nseq,mu=0)
#           pb[iy+1]=PYxTLN(iy,(ix-ix*driv_step),N=ny,g=MM,H=1,sigma=isd,signal_lambda=ilambda,Nseq = Nseq,mu=0)
#         }
#         
#         pf=pf/sum(pf)
#         pc=pc/sum(pc)
#         pb=pb/sum(pb)
#         
#         for (iy in 0:ny){
#           Pix[iy+1]=pc[iy+1]
#           FIix[iy+1]=(1/(ix*driv_step)^2)*(
#             log(pf[iy+1]^pc[iy+1])-
#               2*log(pc[iy+1]^pc[iy+1])+
#               log(pb[iy+1]^pc[iy+1])
#           )
#         }
#         
#         dermu_x=(sum((0:ny)*pf)-sum((0:ny)*pb))/(2*ix*driv_step)
#         var_x=sum(((0:ny)^2)*pc)-(sum((0:ny)*pc))^2
#         
#         Plist[[k]]=Pix
#         FIlist[[k]]=FIix
#         FISN[k]= (dermu_x^2)/var_x
#         FI[k]=max(-sum(FIix[is.finite(FIix)]),0)
#       }
#       
#       
#       Pmatr=do.call(rbind,Plist)
#       
#       FI[!is.finite(FI)]=0
#       CCC=trapz_simple(xTs,sqrt(FI))
#       CCCSN=trapz_simple(xTs,sqrt(FISN))
#       JP=sqrt(FI)/CCC
#       JPSN=sqrt(FISN)/CCCSN
#       vecFI=1:is
#       
#       CCC=trapz_simple(xTs[vecFI],sqrt(FI)[vecFI])
#       CCCC=log( (1/(2*pi*exp(1))^0.5) * CCC)
#       CCCCb=log2(exp(CCCC))
#       
#       
#       ProbX =  apply(Pmatr[vecFI,],2,function(x) trapz_simple(xTs[vecFI],x*((JP[vecFI])) ) )
#       HX = -sum(log(ProbX^ProbX));
#       HXS = trapz_simple(xTs[vecFI],((JP[vecFI]))*apply(Pmatr[vecFI,],1,function(x) (-sum(log(x^x))) )  )
#       
#       output_m=(apply(Pmatr,1,function(x) sum(0:ny*x)))
#       output_sd=sqrt(apply(Pmatr,1,function(x) sum((0:ny)^2*x)-(sum(0:ny*x))^2 ))
#       Cbialek=trapz_simple(output_m,1/output_sd)
#       Cbialek2=trapz_simple(xTs,sqrt(FISN))
#       
#       outputList[["lognormal"]][[as.character(ny)]][[as.character(isd)]][[as.character(ilambda)]]=data.frame(jpexact=log2(exp(HX-HXS)),
#                                                                                                              jpapp=CCCCb,
#                                                                                                              sn1=log2( Cbialek/sqrt(2*pi*exp(1)) ),
#                                                                                                              sn2=log2( Cbialek2/sqrt(2*pi*exp(1)) ))
#       
#       outputListOpt[["lognormal"]][[as.character(ny)]][[as.character(isd)]][[as.character(ilambda)]]$JP=data.frame(s=xTs,popt=JP)
#       outputListOpt[["lognormal"]][[as.character(ny)]][[as.character(isd)]][[as.character(ilambda)]]$SN=data.frame(s=xTs,popt=JPSN)
#       
#       
#       # c(log2(exp(HX-HXS)),CCCCb,log2( Cbialek/sqrt(2*pi*exp(1)) ),log2( Cbialek2/sqrt(2*pi*exp(1)) ))
#       
#       
#       
#       write.table(Pmatr,paste("ver7_lnorm_Ps",isd,"_sig",is,"_n",ny,"_lambda",ilambda,".csv",sep=""),
#                   col.names = FALSE,
#                   row.names = FALSE,
#                   sep=",")
#       
#     }
#   }
# }


#saveRDS(outputList,file="output_list_ver7.rds")
#saveRDS(outputListOpt,file="outputPopt_list_ver7.rds")
# outputList=readRDS(file="output_list_ver7.rds")
# outputListOpt=readRDS(file="outputPopt_list_ver7.rds")


## Calculations - Lognormal distribution ####
outputList[["lognormal2"]]=list()
outputListOpt[["lognormal2"]]=list()
for (ny in nys){
  outputList[["lognormal2"]][[as.character(ny)]]=list()
  outputListOpt[["lognormal2"]][[as.character(ny)]]=list()
  for (isd in isds){
    outputList[["lognormal2"]][[as.character(ny)]][[as.character(isd)]]=list()
    outputListOpt[["lognormal2"]][[as.character(ny)]][[as.character(isd)]]=list()
    for (ilambda in ilambdas){
      
      # ny=500
      # isd=10
      # 
      #xTs=seq(10^-8,10^8,length.out = is)
      #xTs=c(10^seq(-4,0,length.out = is/4),seq(from=1.01,to=10^2,length.out=is/4),
      #      seq(from=10^2+1,to=10^4,length.out=is/4),seq(from=10^4+1,to=10^6,length.out=is/4))
      
      tempv=isd^2
      tempy=polyroot(c(-tempv,0,0,-1,1))
      tempy=Re(tempy[abs(Im(tempy))<1e-6&Re(tempy)>0])
      tempsd=sqrt(log(tempy))
      tempmu=tempsd^2
      # testsamp=rlnorm(100000,meanlog=tempmu,sdlog=tempsd)
      # tempdensity=density(testsamp,from=0,to=exp(tempmu+5*tempsd),n=10000)
      # 
      # 
      # c(mean(testsamp),exp(1.5*tempsd^2))
      # c(tempdensity$x[which.max(tempdensity$y)],1)
      # c(sd(testsamp),isd)
      
      FI=c()
      FIlist=list()
      FISN=c()
      Plist=list()
      k=0
      for (ix in xTs){
        k=k+1
        if (k%%20 == 0){print(k)}
        FIix=c()
        Pix=c()
        pf=c()
        pc=c()
        pb=c()
        for (iy in 0:ny){
          pf[iy+1]=PYxTLN2(iy,(ix+ix*driv_step),N=ny,g=MM,H=1,sigma=tempsd,signal_lambda=ilambda,Nseq = Nseq,mu=tempmu)
          pc[iy+1]=PYxTLN2(iy,ix,N=ny,g=MM,H=1,sigma=tempsd,signal_lambda=ilambda,Nseq = Nseq,mu=tempmu)
          pb[iy+1]=PYxTLN2(iy,(ix-ix*driv_step),N=ny,g=MM,H=1,sigma=tempsd,signal_lambda=ilambda,Nseq = Nseq,mu=tempmu)
        }
        
        pf=pf/sum(pf)
        pc=pc/sum(pc)
        pb=pb/sum(pb)
        
        for (iy in 0:ny){
          Pix[iy+1]=pc[iy+1]
          FIix[iy+1]=(1/(ix*driv_step)^2)*(
            log(pf[iy+1]^pc[iy+1])-
              2*log(pc[iy+1]^pc[iy+1])+
              log(pb[iy+1]^pc[iy+1])
          )
        }
        
        dermu_x=(sum((0:ny)*pf)-sum((0:ny)*pb))/(2*ix*driv_step)
        var_x=sum(((0:ny)^2)*pc)-(sum((0:ny)*pc))^2
        
        Plist[[k]]=Pix
        FIlist[[k]]=FIix
        FISN[k]= (dermu_x^2)/var_x
        FI[k]=max(-sum(FIix[is.finite(FIix)]),0)
      }
      
      
      Pmatr=do.call(rbind,Plist)
      
      FI[!is.finite(FI)]=0
      CCC=trapz_simple(xTs,sqrt(FI))
      CCCSN=trapz_simple(xTs,sqrt(FISN))
      JP=sqrt(FI)/CCC
      JPSN=sqrt(FISN)/CCCSN
      vecFI=1:is
      
      CCC=trapz_simple(xTs[vecFI],sqrt(FI)[vecFI])
      CCCC=log( (1/(2*pi*exp(1))^0.5) * CCC)
      CCCCb=log2(exp(CCCC))
      
      
      ProbX =  apply(Pmatr[vecFI,],2,function(x) trapz_simple(xTs[vecFI],x*((JP[vecFI])) ) )
      HX = -sum(log(ProbX^ProbX));
      HXS = trapz_simple(xTs[vecFI],((JP[vecFI]))*apply(Pmatr[vecFI,],1,function(x) (-sum(log(x^x))) )  )
      
      output_m=(apply(Pmatr,1,function(x) sum(0:ny*x)))
      output_sd=sqrt(apply(Pmatr,1,function(x) sum((0:ny)^2*x)-(sum(0:ny*x))^2 ))
      Cbialek=trapz_simple(output_m,1/output_sd)
      Cbialek2=trapz_simple(xTs,sqrt(FISN))
      
      outputList[["lognormal2"]][[as.character(ny)]][[as.character(isd)]][[as.character(ilambda)]]=data.frame(jpexact=log2(exp(HX-HXS)),
                                                                                                              jpapp=CCCCb,
                                                                                                              sn1=log2( Cbialek/sqrt(2*pi*exp(1)) ),
                                                                                                              sn2=log2( Cbialek2/sqrt(2*pi*exp(1)) ))
      
      outputListOpt[["lognormal2"]][[as.character(ny)]][[as.character(isd)]][[as.character(ilambda)]]$JP=data.frame(s=xTs,popt=JP)
      outputListOpt[["lognormal2"]][[as.character(ny)]][[as.character(isd)]][[as.character(ilambda)]]$SN=data.frame(s=xTs,popt=JPSN)
      
      
      # c(log2(exp(HX-HXS)),CCCCb,log2( Cbialek/sqrt(2*pi*exp(1)) ),log2( Cbialek2/sqrt(2*pi*exp(1)) ))
      
      
      
      write.table(Pmatr,paste("output/pyx_matr/ver7_lnorm2_Ps",isd,"_sig",is,"_n",ny,"_lambda",ilambda,".csv",sep=""),
                  col.names = FALSE,
                  row.names = FALSE,
                  sep=",")
      
    }
  }
}


# Writing large objects storing results to RDS files
saveRDS(outputList,file="output/output_list_ver7.rds")
saveRDS(outputListOpt,file="output/outputPopt_list_ver7.rds")







