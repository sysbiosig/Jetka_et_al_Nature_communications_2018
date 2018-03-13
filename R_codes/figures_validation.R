# Custom code accompanying manuscript NCOMMS-18-04749:
# "An information-theoretic framework for deciphering pleiotropic and noisy biochemical signaling"

# Code I: Validation of capacity approximation by using Jeffrey's prior
# R-script for generating figures in Main Paper Results section

# Generating figures
# Use after runing R_codes/capacity_validation.R and Matlab_codes_BA
# Assumes that Blahut-Arimoto calculations are in directory:
# "../Matlab_codes_BA/falsesignal/"

## Set your workings directory ####
#setwd("~/Documents/Science/Paper_codes/R_codes/")
setwd("F:\\TJ\\Simulations\\Paper_codes\\R_codes")

## Libraries needed for calculations ####
rm(list=ls())
library('latex2exp')
library(reshape2)
library(ggplot2)
library(ggthemes)
library(gridExtra)
source("aux_functions_figures.R")


## Setting up variables and costants ####
nys=c(10,50,100,500) # Setting up number of molecules
isds=c(0.01,0.1,0.5,1,5,10,25,50) # Setting up standard variations
ilambdas=c(1,0.5,0.1) # Setting up lambdas
driv_step=0.0001 # shift parameter for derivative calculation
is=200  # Setting up density of signal
xTs=10^seq(-8,8,length.out = is)  # Setting up range of signal
Nseq=2000 # Setting up number sample size used to approximate single marginal distribution

dir.create("output")
dir.create("paper_output")

## Generating exploratory plots ####
outputList=readRDS(file="output/output_list_ver7.rds")
outputListOpt=readRDS(file="output/outputPopt_list_ver7.rds")


outputDF=melt(outputList,id.vars=c("jpexact",
                                   "jpapp",
                                   "sn1",
                                   "sn2"))
colnames(outputDF)[5:8]<-c("lambda","sd","N","dist")
outputDF2=melt(outputDF,id.vars=c("lambda","sd","N","dist"))
colnames(outputDF2)<-c("lambda","sd","N","dist","method","capacityBit")

#saveRDS(outputDF,file="output_df_ver7.rds")
#saveRDS(outputDF2,file="output_df2_ver7.rds")

outputDF2$lambda<-as.numeric(outputDF2$lambda)
outputDF2$sd<-as.numeric(outputDF2$sd)
outputDF2$N<-as.numeric(outputDF2$N)

outputDF2$distlambda=paste(outputDF2$dist,outputDF2$lambda,sep="_")
plot=ggplot(data=outputDF2[!outputDF2$dist=="nonoise",],aes(x=sd,y=capacityBit,colour=method))+geom_line(size=2)+
  facet_grid(distlambda~N,scales="free_y")+theme_publ(version=2)

ggsave(plot,file="output/plot_withoutBA_ver7.pdf",width=10,height=20)


outputDF2$distlambda=paste(outputDF2$dist,outputDF2$lambda,sep="_")
plot=ggplot(data=outputDF2[outputDF2$dist=="lognormal2",],aes(x=sd,y=capacityBit,colour=method))+geom_line(size=2)+
  facet_grid(distlambda~N,scales="free_y")+theme_publ(version=2)

ggsave(plot,file="output/plot_withoutBA_ver7_lognormal2.pdf",width=10,height=20)

dists=c("exp","lnorm2","gamma")
outputDF3=outputDF2
outputDF4=list()
for (idist in dists){
  tempdf=read.table(file=paste("../Matlab_codes_BA/output/falsesignal/",idist,"_capacity.csv",sep=""),
                    header=FALSE,sep=",")
  colnames(tempdf)<-c("sd","N","lambda","capacityBit")
  if (idist=="lnorm"){
    idist="lognormal"
  }
  if (idist=="lnorm2"){
    idist="lognormal2"
  }
  
  tempdf=tempdf[!(tempdf$lambda==0&tempdf$N==0),]
  tempdf=data.frame(tempdf,dist=idist,method="blahut")
  tempdf$distlambda=paste(tempdf$dist,tempdf$lambda,sep="_")
  outputDF3=rbind(outputDF3,tempdf[,c(3,1,2,5,6,4,7)])
  outputDF4[[idist]]=tempdf
}

outputDF4=do.call(rbind,outputDF4)
outputDF4=merge(outputDF4,outputDF2,by=c( "sd", "N","lambda","dist","distlambda" ))
outputDF4$error=abs(outputDF4$capacityBit.y-outputDF4$capacityBit.x)/outputDF4$capacityBit.x

plot=ggplot(data=outputDF3[!outputDF3$dist=="nonoise",],aes(x=sd,y=capacityBit,colour=method))+geom_line(size=2)+
  facet_grid(distlambda~N,scales="free_y")+theme_publ(version=2)

ggsave(plot,file="output/plot_withBA_ver7.pdf",width=10,height=30)


plot=ggplot(data=outputDF4[!outputDF4$dist=="nonoise",],aes(x=sd,y=error,colour=method.y))+geom_line(size=2)+
  facet_grid(distlambda~N,scales="free_y")+theme_publ(version=2)+scale_y_log10()+scale_x_log10()

ggsave(plot,file="output/plot_Errors_ver7.pdf",width=10,height=30)


dataToPlot=aggregate(data=outputDF4[outputDF4$dist%in%c("exp","gamma"),],error~method.y+sd,mean)
plot=ggplot(data=dataToPlot,aes(x=sd,y=error,colour=method.y))+geom_line(size=2)+
  theme_publ(version=2)
plot

ggsave(plot,file="output/plot_MeanErrors_ver7.pdf",width=10,height=30)

dataToPlot=aggregate(data=outputDF4[outputDF4$dist%in%c("exp","gamma","lognormal2"),],error~method.y+sd,mean)
plot=ggplot(data=dataToPlot,aes(x=sd,y=error,colour=method.y))+geom_line(size=2)+
  theme_publ(version=2)
plot

ggsave(plot,file="output/plot_MeanErrors_ver7_withLN2.pdf",width=10,height=30)

dists=c("exp","lnorm2","gamma")
for (idist in dists){
  print(idist)
  if (idist=="lnorm"){
    idist2="lognormal"
  } else {
    idist2=idist
  }
  if (idist=="lnorm2"){
    idist2="lognormal2"
  } else {
    idist2=idist
  }
  for (ny in nys){
    for (isd in isds){
      for (ilambda in ilambdas){
        
        temppath=paste('../Matlab_codes_BA/output/falsesignal/ver7_',idist,'_sd_',isd,'inputnum_',is,'_n',ny,'_lambda',ilambda,'/',sep="")
        
        tempdf=dename(read.table(file=paste(temppath,"Q.csv",sep=""),
                                 header=FALSE,sep=","))
        tempBA=tempdf[(length(tempdf)-length(xTs)+1):(length(tempdf)) ]
        outputListOpt[[idist2]][[as.character(ny)]][[as.character(isd)]][[as.character(ilambda)]]$BA=data.frame(s=xTs,popt=tempBA)
        
        divn=4
        for (idiv in 1:(length(xTs)/divn)){
          tempBA[((idiv-1)*divn+1):((idiv)*divn)]=sum(tempBA[((idiv-1)*divn+1):((idiv)*divn)])/divn
        }
        outputListOpt[[idist2]][[as.character(ny)]][[as.character(isd)]][[as.character(ilambda)]]$BA2=data.frame(s=xTs,popt=tempBA)
        
      }
    }
  }
}


# Saving obejct after renaming
saveRDS(outputListOpt,file="output/outputPopt_listFull_ver7.rds")


outputListOpt=readRDS(file="output/outputPopt_listFull_ver7.rds")

temp_outputpath="output/ver7_optPlots/"
dir.create(temp_outputpath)

dists=c("exp","lognormal2","gamma")
for (idist in dists){
  for (ny in nys){
    for (isd in isds){
      for (ilambda in ilambdas){
        
        tempdata0=outputListOpt[[idist]][[as.character(ny)]][[as.character(isd)]][[as.character(ilambda)]]
        tempdata0$JP$popt=tempdata0$JP$popt*tempdata0$JP$s
        tempdata0$JP$popt=tempdata0$JP$popt/sum(tempdata0$JP$popt)
        tempdata0$SN$popt=tempdata0$SN$popt*tempdata0$SN$s
        tempdata0$SN$popt=tempdata0$SN$popt/sum(tempdata0$SN$popt)
        tempdata=melt(tempdata0,id.vars=c("s","popt"))
        
        temppath=paste(temp_outputpath,'1_ver7_',idist,'_sd_',isd,'inputnum_',is,'_n',ny,'_lambda',ilambda,'.pdf',sep="")
        plot=ggplot(data=tempdata,aes(x=s,y=popt,colour=L1))+geom_line()+scale_x_log10()+theme_publ(version=2)
        ggsave(plot,file=temppath)
        
        temppath=paste(temp_outputpath,'2_ver7_',idist,'_sd_',isd,'inputnum_',is,'_n',ny,'_lambda',ilambda,'.pdf',sep="")
        plot=ggplot(data=tempdata[!tempdata$L1=="BA2",],aes(x=s,y=popt,colour=L1))+geom_smooth(method="loess",n=100,span=0.08,alpha=0.5,se=FALSE,size=2)+
          scale_x_log10()+scale_y_continuous(limits=c(0,NA))+theme_publ(version=2)
        ggsave(plot,file=temppath)
        
        temppath=paste(temp_outputpath,'3_ver7_',idist,'_sd_',isd,'inputnum_',is,'_n',ny,'_lambda',ilambda,'.pdf',sep="")
        plot=ggplot(data=tempdata,aes(x=s,y=popt,colour=L1))+geom_smooth(method="loess",n=100,span=0.08,alpha=0.5,se=FALSE,size=1)+
          geom_point(size=2)+
          scale_x_log10()+scale_y_continuous(limits=c(0,NA))+theme_publ(version=2)
        ggsave(plot,file=temppath)
        
      }
    }
  }
}


## Generating figures for paper####
outputList=readRDS(file="output/output_list_ver7.rds")
outputListOpt=readRDS(file="output/outputPopt_list_ver7.rds")

outputDF=melt(outputList,id.vars=c("jpexact",
                                   "jpapp",
                                   "sn1",
                                   "sn2"))
colnames(outputDF)[5:8]<-c("lambda","sd","N","dist")
outputDF2=melt(outputDF,id.vars=c("lambda","sd","N","dist"))
colnames(outputDF2)<-c("lambda","sd","N","dist","method","capacityBit")


outputDF2$lambda<-as.numeric(outputDF2$lambda)
outputDF2$sd<-as.numeric(outputDF2$sd)
outputDF2$N<-as.numeric(outputDF2$N)
outputDF2$distlambda=paste(outputDF2$dist,outputDF2$lambda,sep="_")

dists=c("exp","lnorm2","gamma")
outputDF3=outputDF2
outputDF4=list()
#Loading Blahut Arimoto
for (idist in dists){
  tempdf=read.table(file=paste("../Matlab_codes_BA/output/falsesignal/",idist,"_capacity.csv",sep=""),
                    header=FALSE,sep=",")
  colnames(tempdf)<-c("sd","N","lambda","capacityBit")
  if (idist=="lnorm"){
    idist="lognormal"
  }
  if (idist=="lnorm2"){
    idist="lognormal2"
  }
  
  tempdf=tempdf[!(tempdf$lambda==0&tempdf$N==0),]
  tempdf=data.frame(tempdf,dist=idist,method="blahut")
  tempdf$distlambda=paste(tempdf$dist,tempdf$lambda,sep="_")
  outputDF3=rbind(outputDF3,tempdf[,c(3,1,2,5,6,4,7)])
  outputDF4[[idist]]=tempdf
}

outputDF4=do.call(rbind,outputDF4)
outputDF4=merge(outputDF4,outputDF2,by=c( "sd", "N","lambda","dist","distlambda" ))
outputDF4$error=abs(outputDF4$capacityBit.y-outputDF4$capacityBit.x)/outputDF4$capacityBit.x


outputListOpt=readRDS(file="output/outputPopt_listFull_ver7.rds")



plot=ggplot(data=outputDF3[(outputDF3$dist%in%c("exp","gamma")&
                              outputDF3$lambda%in%c(0.5)&
                              outputDF3$N%in%c(50,500)&
                              outputDF3$method%in%c("jpexact","jpapp","sn2","blahut")),],
            aes(x=sd,y=capacityBit,colour=method))+geom_line(size=2)+
  facet_grid(distlambda~N,scales="free_y")+theme_publ(version=2)#+scale_y_continuous(limits=c(0,NA))

ggsave(plot,file="paper_output/plot_withBA_ver7.pdf",width=8,height=6)


plot=ggplot(data=outputDF3[(outputDF3$dist%in%c("gamma","lognormal2")&
                              outputDF3$lambda%in%c(0.5)&
                              outputDF3$N%in%c(50,500)&
                              outputDF3$method%in%c("jpexact","jpapp","sn2","blahut")),],
            aes(x=sd,y=capacityBit,colour=method))+geom_line(size=2)+
  facet_grid(distlambda~N,scales="free_y")+theme_publ(version=2)+scale_y_continuous(limits=c(0,NA))

ggsave(plot,file="paper_output/plot_withBAwithLN_ver7.pdf",width=8,height=6)

outputDF4$capacityBit.y0=outputDF4$capacityBit.y
outputDF4$capacityBit.y0[outputDF4$capacityBit.y0<0]=0
outputDF4$error0=abs(outputDF4$capacityBit.y0-outputDF4$capacityBit.x)/outputDF4$capacityBit.x

plot=ggplot(data=outputDF4[!outputDF4$dist=="nonoise",],aes(x=sd,y=error0,colour=method.y))+geom_line(size=2)+
  facet_grid(distlambda~N,scales="free_y")+theme_publ(version=2)+scale_y_log10()+scale_x_log10()

ggsave(plot,file="paper_output/plot_Errors1_ver7.pdf",width=10,height=20)


plot=ggplot(data=outputDF4[!outputDF4$dist=="nonoise",],aes(x=sd,y=error0,colour=method.y))+geom_line(size=2)+
  facet_grid(distlambda~N,scales="free_y")+theme_publ(version=2)#+scale_y_log10()+scale_x_log10()

ggsave(plot,file="paper_output/plot_Errors2_ver7.pdf",width=10,height=20)


plot=ggplot(data=outputDF3[!outputDF3$dist=="nonoise",],aes(x=sd,y=capacityBit,colour=method))+geom_line(size=2)+
  facet_grid(distlambda~N,scales="free_y")+theme_publ(version=2)

ggsave(plot,file="paper_output/plot_withBAAll_ver7.pdf",width=10,height=20)


tempcond=(outputDF4$dist%in%c("exp","gamma","lognormal2"))&(outputDF4$lambda%in%c(0.1,0.5))&(outputDF4$N%in%c(50,500))
dataToPlot=aggregate(data=outputDF4[tempcond,],error0~method.y+sd,mean)
plot=ggplot(data=dataToPlot,aes(x=sd,y=error0,colour=method.y))+geom_line(size=2)+
  theme_publ(version=2)
plot

ggsave(plot,file="paper_output/plot_MeanErrors1_ver7.pdf",width=8,height=6)

tempcond=(outputDF4$dist%in%c("exp","gamma","lognormal2"))
dataToPlot=aggregate(data=outputDF4[tempcond,],error0~method.y+sd,mean)
plot=ggplot(data=dataToPlot,aes(x=sd,y=error0,colour=method.y))+geom_line(size=2)+
  theme_publ(version=2)
plot

ggsave(plot,file="paper_output/plot_MeanErrors2_ver7.pdf",width=8,height=6)


tempcond=(outputDF4$dist%in%c("exp","gamma","lognormal2"))&(outputDF4$lambda%in%c(0.1,0.5,1))&(outputDF4$N%in%c(50,100,500))
dataToPlot=aggregate(data=outputDF4[tempcond,],error0~method.y+sd,mean)
plot=ggplot(data=dataToPlot,aes(x=sd,y=error0,colour=method.y))+geom_line(size=2)+
  theme_publ(version=2)
plot

ggsave(plot,file="paper_output/plot_MeanErrors3_ver7.pdf",width=8,height=8)



tempcond=(outputDF4$dist%in%c("exp"))
dataToPlot=aggregate(data=outputDF4[tempcond,],error0~method.y+sd,mean)
plot=ggplot(data=dataToPlot,aes(x=sd,y=error0,colour=method.y))+geom_line(size=2)+
  theme_publ(version=2)
plot

ggsave(plot,file="paper_output/plot_MeanErrors_exp_ver7.pdf",width=8,height=6)

tempcond=(outputDF4$dist%in%c("gamma"))
dataToPlot=aggregate(data=outputDF4[tempcond,],error0~method.y+sd,mean)
plot=ggplot(data=dataToPlot,aes(x=sd,y=error0,colour=method.y))+geom_line(size=2)+
  theme_publ(version=2)
plot

ggsave(plot,file="paper_output/plot_MeanErrors_gamma_ver7.pdf",width=8,height=6)

tempcond=(outputDF4$dist%in%c("lognormal2"))
dataToPlot=aggregate(data=outputDF4[tempcond,],error0~method.y+sd,mean)
plot=ggplot(data=dataToPlot,aes(x=sd,y=error0,colour=method.y))+geom_line(size=2)+
  theme_publ(version=2)
plot

ggsave(plot,file="paper_output/plot_MeanErrors_lognormal_ver7.pdf",width=8,height=6)


tempcond=(outputDF4$dist%in%c("exp"))&(outputDF4$lambda%in%c(0.1,0.5))&(outputDF4$N%in%c(50,100,500))
dataToPlot=aggregate(data=outputDF4[tempcond,],error0~method.y+sd,mean)
plot=ggplot(data=dataToPlot,aes(x=sd,y=error0,colour=method.y))+geom_line(size=2)+
  theme_publ(version=2)
plot

ggsave(plot,file="paper_output/plot_MeanErrors_exp_chosen_ver7.pdf",width=8,height=6)

tempcond=(outputDF4$dist%in%c("gamma"))&(outputDF4$lambda%in%c(0.1,0.5))&(outputDF4$N%in%c(50,100,500))
dataToPlot=aggregate(data=outputDF4[tempcond,],error0~method.y+sd,mean)
plot=ggplot(data=dataToPlot,aes(x=sd,y=error0,colour=method.y))+geom_line(size=2)+
  theme_publ(version=2)
plot

ggsave(plot,file="paper_output/plot_MeanErrors_gamma_chosen_ver7.pdf",width=8,height=6)

tempcond=(outputDF4$dist%in%c("lognormal2"))&(outputDF4$lambda%in%c(0.1,0.5))&(outputDF4$N%in%c(50,100,500))
dataToPlot=aggregate(data=outputDF4[tempcond,],error0~method.y+sd,mean)
plot=ggplot(data=dataToPlot,aes(x=sd,y=error0,colour=method.y))+geom_line(size=2)+
  theme_publ(version=2)
plot

ggsave(plot,file="paper_output/plot_MeanErrors_lognormal_chosen_ver7.pdf",width=8,height=6)

