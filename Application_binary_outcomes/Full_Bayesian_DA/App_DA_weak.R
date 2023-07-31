


source("Functions_DA_param_expand.R")

#load("data_for_run_2yrs_20112017_v3.Rdata")

burn_in=5000
Iteration=30000
thin=10


prior_val=list(prior_A=2,prior_C=2,prior_N=2,prior.sd=1,
               A=2,proposal_rw_width=4,proposal_rt_df=5,prior_step_size=21^(-1/3),
               prior.gamma.shape=1,
               prior.gamma.scale=0.5,
               prior.gamma.mu=0,
               prior.gamma.sd=2,
               prior.epsilon.shape=1,
               prior.epsilon.scale=1,
               proposal_gamma_width=0.5)

ini_val=list(w_A_ini=1/3,w_C_ini=1/3,w_N_ini=1/3,Sigma_beta=diag(0.05,ncol(Datalist$X)),
             alpha_A_ini=0.1,alpha_C_ini=0.1,alpha_N_ini=0.1, alpha_ini=1, gamma_ini=0.5,
             mu_alpha_ini=0.1,sigma_alpha_ini=5,Tau_ini=0.1)

model_DA_weak<-MCMC_Cross_temporal_binary(Iteration=Iteration,Data=Datalist, random.ini = T,
                                          prior_values=prior_val,
                                          ini_values=ini_val,
                                          burn_in=burn_in,thin=thin,printYes=T)

#save(model_DA_weak,file="App_DA_weak_20112017_res.Rdata")
