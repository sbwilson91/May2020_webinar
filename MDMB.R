##################################################################################
#  mdmb -three-level MI                                                          #
#  Simulation Study 2                                                           #
#  Rushani Wijesuriya 3 Dec 2019                                                 #
##################################################################################
rm(list=ls())

library(miceadds)
library(mdmb)
library(DataCombine)
library(caret)

library(mitml)
library(lme4)
require(mice)
require(xlsx)




 args = commandArgs(trailingOnly=TRUE)  #these arguments are passed on from a command line that sends script to a HPC server
# # Alternatively, can comment out and set parameters as below
#
 datnum<-as.numeric(args[1])

T1 = list.files(pattern="*.csv")

 s2=seq(2,6,by=2)
 s1=s2-1
 
#s=seq(20,1000,by=20)

temp=T1[s1[datnum]:s2[datnum]]


#temp = list.files(pattern="*.csv")

data=list()

for (i in 1:length(temp)){
  dataL=read.csv(temp[i])
  data[[i]]= dataL
}

#---------------------------------------------------------------------------------

MDMB_results.est=matrix(NA,nrow=7,ncol=length(temp))
MDMB_results.sd=matrix(NA,nrow=7,ncol=length(temp))
MDMB_results.RE=matrix(NA,nrow=3,ncol=length(temp))
MDMB_results.ICC=matrix(NA,nrow=2,ncol=length(temp))

MDMB_results.CI=c()

for (i in 1:length(temp)){
  ST=Sys.time()
  print(i)
  simdataL= data[[i]]
  simdataL <- simdataL[order(simdataL$school,simdataL$child),]
  simdataL$c_id=as.numeric(simdataL$c_id)
  
  #---mdmb for three-level data
  
  simdataL<- slide(simdataL, Var = "c_dep", GroupVar = "c_id",
                    slideBy = -1)
  
  colnames(simdataL)[colnames(simdataL)=="c_dep-1"] <- "prev_dep"
  
  simdataL$c_dep=NULL
  
  
  #create previous wave SDQ variable
  
  simdataL<- slide(simdataL, Var = "p_sdq", GroupVar = "c_id",
                    slideBy = -1)
  
  colnames(simdataL)[colnames(simdataL)=="p_sdq-1"] <- "prev_sdq"
  
  simdataL$p_sdq=NULL
  
  #remove unwanted waves
  simdataL=subset(simdataL, wave!= 2 & wave!=4 & wave!=6)
  
  simdataL$c_gender=ifelse(simdataL$c_gender=="male",1,0)
  
  
  simdataL=simdataL[,!names(simdataL)%in%c("child")]
  

  ##create DI for schools
  school_DI=data.frame(model.matrix(simdataL$c_id~as.factor(simdataL$school)-1,
                                    simdataL))
  
  
  names(school_DI)[1:ncol(school_DI)] <- unlist(mapply(function(x,y) paste(x, seq(1,y), sep="_"), 
                                                       "schoo_Ind",40))
  school_DI=school_DI[,1:39]
  
  
  ##combine all data
  simdataL1=cbind(simdataL[,!names(simdataL)%in%c("school")],school_DI)
  
  
  #*** model specification
  #mcmc_iter <- 4 # number of MCMC iterations for model parameter sampling
  
  iter <-20 ; burnin <- 1   ##burn in for 1000 iteration and one imputed dataset is saved every 100th iteration=20imps? 
  Nimp <- 2
  
  
  # miceadds::cwc() and miceadds::gm() are functions for computing group-mean
  # centered variables and the group mean, respectively.
  # miceadds::cwc() and miceadds::gm() are functions for computing group-mean
  # centered variables and the group mean, respectively.
  Y_formula <- napscore_z ~prev_dep+wave+prev_dep*wave+c_age+c_ses+c_gender+c_nap1_z+prev_sdq+schoo_Ind_1+
    schoo_Ind_2+schoo_Ind_3+schoo_Ind_4+schoo_Ind_5+schoo_Ind_6+schoo_Ind_7+schoo_Ind_8+
    schoo_Ind_9+schoo_Ind_10+schoo_Ind_11+schoo_Ind_12+schoo_Ind_13+schoo_Ind_14+schoo_Ind_15+
    schoo_Ind_16+schoo_Ind_17+schoo_Ind_18+schoo_Ind_19+schoo_Ind_20+schoo_Ind_21+schoo_Ind_22+
    schoo_Ind_23+schoo_Ind_24+schoo_Ind_25+schoo_Ind_26+schoo_Ind_27+schoo_Ind_28+schoo_Ind_29+
    schoo_Ind_30+schoo_Ind_31+schoo_Ind_32+schoo_Ind_33+schoo_Ind_34+schoo_Ind_35+schoo_Ind_36+schoo_Ind_37+
    schoo_Ind_38+schoo_Ind_39+(1|c_id)
  
  # model for dependent variable
  dep <- list("model"="mlreg", "formula"=Y_formula )
  # dep <- list("model"="mlreg", "formula"=model_formula,
  #             R_args=list(iter=mcmc_iter, outcome="normal") )
  
  
  
  ##covariate models
  X1_formula=prev_dep~wave+c_age+c_ses+c_gender+c_nap1_z+prev_sdq+schoo_Ind_1+
    schoo_Ind_2+schoo_Ind_3+schoo_Ind_4+schoo_Ind_5+schoo_Ind_6+schoo_Ind_7+schoo_Ind_8+
    schoo_Ind_9+schoo_Ind_10+schoo_Ind_11+schoo_Ind_12+schoo_Ind_13+schoo_Ind_14+schoo_Ind_15+
    schoo_Ind_16+schoo_Ind_17+schoo_Ind_18+schoo_Ind_19+schoo_Ind_20+schoo_Ind_21+schoo_Ind_22+
    schoo_Ind_23+schoo_Ind_24+schoo_Ind_25+schoo_Ind_26+schoo_Ind_27+schoo_Ind_28+schoo_Ind_29+
    schoo_Ind_30+schoo_Ind_31+schoo_Ind_32+schoo_Ind_33+schoo_Ind_34+schoo_Ind_35+schoo_Ind_36+schoo_Ind_37+
    schoo_Ind_38+schoo_Ind_39+(1|c_id)
  
  ind_x <- list( "model"="mlreg", "formula"=X1_formula,
                 sampling_level="c_id" )
  
  
  
  X2_formula=c_ses~c_age+c_gender+c_nap1_z
  
  ind_x2<-list( "model"="linreg", "formula"=X2_formula,
                variable_level="c_id")
  
  #ind <- list(prev_dep=ind_x)
  ind <- list(c_ses=ind_x2,prev_dep=ind_x)
  
  # --- estimate model
  mod1 <- mdmb::frm_fb(simdataL1, dep, ind,iter=iter,burnin=burnin, Nimp =20,aggregation=TRUE)
  
  #** convert output into list of imputed datasets
  datlist <- mdmb::frm2datlist(mod1)
  
  
  
  mylist=list()
  
  for(j in 1:Nimp)
  {
    
    #extract the mth imputed dataset
     datL<-datlist[[j]]
    
    
    datL$school=simdataL$school
    
    datL=datL[,names(datL)%in% c("c_id","wave","napscore_z","c_age","c_gender","c_nap1_z","prev_dep",
                                 "prev_sdq","c_ses","school")]
    
    

    
    #save the dataset in a list
    mylist[[j]]= datL
    
    
  }
  
  #fit the analysis of interest on the imputed datasets 
  mods <- lapply(mylist,function(d) {lmer( napscore_z~prev_dep+wave+prev_dep*wave+c_age+
                                             as.factor(c_gender)+c_nap1_z+c_ses
                                           +(1|school/c_id), data = d)} )
  
  
  MI_est=testEstimates(mods, var.comp=TRUE,df.com=NULL)
  
  
  
  
  #store the estimates
  MDMB_results.est[,i]=MI_est$estimates[2:8,1]
  MDMB_results.sd[,i]=MI_est$estimates[2:8,2]
  MDMB_results.RE[1,i]=sqrt(MI_est$var.comp[2,1])
  MDMB_results.RE[2,i]=sqrt(MI_est$var.comp[1,1])
  MDMB_results.RE[3,i]=sqrt(MI_est$var.comp[3,1])
  
  
  MDMB_results.ICC[1,i]=MI_est$var.comp[2,1]/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  MDMB_results.ICC[2,i]=(MI_est$var.comp[2,1]+MI_est$var.comp[1,1])/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  ET=Sys.time()
  
}

write.xlsx(MDMB_results.est,paste0("MDMB_results.est",datnum,".xlsx"))
write.xlsx(MDMB_results.sd,paste0("MDMB_results.sd",datnum,".xlsx"))
write.xlsx(MDMB_results.RE,paste0("MDMB_results.RE",datnum,".xlsx"))
write.xlsx(MDMB_results.CI,paste0("MDMB_results.CI",datnum,".xlsx"))
write.xlsx(MDMB_results.ICC,paste0("MDMB_results.ICC",datnum,".xlsx"))


# write.xlsx(MDMB_results.est,"MDMB_results.est.xlsx")
# write.xlsx(MDMB_results.sd,"MDMB_results.sd.xlsx")
# write.xlsx(MDMB_results.RE,"MDMB_results.RE.xlsx")
# write.xlsx(MDMB_results.CI,"MDMB_results.CI.xlsx")
# write.xlsx(MDMB_results.ICC,"MDMB_results.ICC.xlsx")
