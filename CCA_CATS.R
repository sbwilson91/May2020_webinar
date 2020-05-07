##################################################################################
# Dummy R script for HOPC demo                                                  #                                                           #
#  Rushani Wijesuriya 03 May  2020
#  
##################################################################################
library(lme4)
library(DataCombine)
library(xlsx)

rm(list = ls())

#import data from directory to a list
temp = list.files(pattern="*.csv")


data=list()

for (i in 1:length(temp)){
  dataL=read.csv(temp[i])
  data[[i]]= dataL
}


comp_results.est=matrix(NA,nrow=7,ncol=length(temp))
comp_results.sd=matrix(NA,nrow=7,ncol=length(temp))
comp_results.RE=matrix(NA,nrow=3,ncol=length(temp))
comp_results.ICC=matrix(NA,nrow=2,ncol=length(temp))
comp_results.CI=c()


for (i in 1:length(temp)){
  print(i)
  simdataL= data[[i]]
  
  simdataL$c_id=as.numeric(simdataL$c_id)
  
  #create previous wave depression
  simdataL_lag<- slide(simdataL, Var = "c_dep", GroupVar = "c_id",
                       slideBy = -1)
  
  colnames(simdataL_lag)[colnames(simdataL_lag)=="c_dep-1"] <- "prev_dep"
  
  #remove unwanted waves
  simdataL_lag=subset(simdataL_lag, wave!= 2 & wave!=4 & wave!=6)
  
  simdataL_lag$p_sdq=NULL
  simdataL_lag$c_dep=NULL
  
  
  ##centre the exposure
  
  
  fit<- lmer(napscore_z~prev_dep+wave+prev_dep*wave+c_age+
               c_gender+c_nap1_z+c_ses
             +(1|school/child), data=simdataL_lag)
  
  x=as.data.frame(VarCorr(fit),comp=c("Variance"))
  
  
  comp_results.est[,i]=coef(summary(fit))[2:8,1]
  comp_results.sd[,i]<-coef(summary(fit))[2:8,2]
  comp_results.RE[1,i]=x[2,5]
  comp_results.RE[2,i]=x[1,5]
  comp_results.RE[3,i]=x[3,5]
  comp_results.ICC[1,i]=x[2,4]/(x[2,4]+x[1,4]+x[3,4])
  comp_results.ICC[2,i]=(x[2,4]+x[1,4])/(x[2,4]+x[1,4]+x[3,4])
 
  
}

rownames(comp_results.est)=c("prev_dep","time","c_age","c_gender","c_nap1_z","c_ses","inter")
colnames(comp_results.est)=c(1:length(temp))

rownames(comp_results.sd)=c("prev_dep","time","c_age","c_gender","c_nap1_z","c_ses","inter")
colnames(comp_results.sd)=c(1:length(temp))

rownames(comp_results.RE)=c("level 3","level 2","level 1")
colnames(comp_results.RE)=c(1:length(temp))

rownames(comp_results.ICC)=c("level 3","level 2")
colnames(comp_results.ICC)=c(1:length(temp))



write.xlsx(comp_results.est,"CCA.compResults_est.xlsx")
write.xlsx(comp_results.sd,"CCA.compResults_sd.xlsx")
write.xlsx(comp_results.RE,"CCA.compResults_RE.xlsx")
write.xlsx(comp_results.ICC,"comp_results.ICC.xlsx")


