
set.seed(1)


source("Functions_matching_common.R")


#load("matching_data_20112017.Rdata")

data<-data_model

data$G_pred<-factor(data$G_pred,levels = c("Always_taker","Compliers","Never_taker"),labels=c("A","C","N"))


Datalist=list(Y=data$Y, D=data$D, Z=data$Z,X=DX.mat,G=data$G_pred)


prior_val=list(A=2,proposal_rw_width=0.5,proposal_rt_df=2)

ini_val=list(alpha_A_ini=0.1,alpha_C_ini=0.1,alpha_N_ini=0.1, alpha_ini=0.1,
             psi_A_ini=c(0.1,0.2,0.1), psi_C_ini=c(0.1,0.2,0.1),  
             mu_alpha_ini=0.01,sigma_alpha_ini=1,Tau_ini=0.1)

model_Matching_common<-Cross_binary_mathching(Datalist=Datalist,Potential.outcome=T,Iteration=20000,
                                            temporal.assumption="common",
                                            prior_values=prior_val,
                                            ini_values=ini_val,
                                            burn_in=5000,thin=2,printYes=T)


#save(model_Matching_common,file="Matching_common_2t_20112017_res.Rdata")
