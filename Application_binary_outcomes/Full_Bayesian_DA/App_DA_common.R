source("Functions_Full_Bayesian.R")

#load("data_for_run_2017_latest.Rdata")

##############################################3

burn_in=5000
Iteration=30000
thin=5

prior_val=list(prior_A=10,prior_C=10,prior_N=10,
               A=2,proposal_rw_width=5,proposal_rt_df=5)

ini_val=list(w_A_ini=1/3,w_C_ini=1/3,w_N_ini=1/3,Sigma_beta=diag(0.05,ncol(Datalist$X)),
             #beta_A_ini=coef(fit0)[-2],beta_N_ini=coef(fit0)[-2],beta_C_ini=coef(fit0)[-2],
             alpha_A_ini=0.01,alpha_C_ini=0.01,alpha_N_ini=0.1, alpha_ini=0.1,
             #psi_A_ini=c(0.1,0.2,0.1), psi_C_ini=c(0.1,0.2,0.1),  
             mu_alpha_ini=0.01,sigma_alpha_ini=10,Tau_ini=0.5)

model_DA_common.2t.probit.c1<-App_Cross_temporal_DA(Iteration, Data=Datalist,
                                                  prior_values=prior_val,
                                                  ini_values=ini_val,
                                                  Strata.model=T, burn_in,thin,printYes=T,Model.name="DA_common_2t_probit_c1")


#save(model_DA_common.2t.probit.c1,file="App_DA_common_tx_2ts_20112017_probit_c1.Rdata")



