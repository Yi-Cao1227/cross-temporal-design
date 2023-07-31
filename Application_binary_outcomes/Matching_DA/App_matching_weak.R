
source("Functions_matching_weak.R")



#load("matching_data_20112017.Rdata")

data$G_pred<-factor(data$G_pred,levels = c("Always_taker","Compliers","Never_taker"),labels=c("A","C","N"))

Datalist=list(Y=data$Y, D=data$D, Z=data$Z,X=DX.mat,G=data$G_pred)

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
             #beta_A_ini=coef(fit0)[-2],beta_N_ini=coef(fit0)[-2],beta_C_ini=coef(fit0)[-2],
             alpha_A_ini=0.1,alpha_C_ini=0.1,alpha_N_ini=0.1, alpha_ini=0, gamma_ini=0.5,
             #psi_A_ini=c(0.1,0.2,0.1), psi_C_ini=c(0.1,0.2,0.1),  
             mu_alpha_ini=1,sigma_alpha_ini=4,Tau_ini=0)



model_Matching_weak<-Cross_binary_mathching(Datalist=Datalist,Potential.outcome=T,Iteration=30000,random.ini=T,
                                            temporal.assumption="weak",
                                            prior_values=prior_val,
                                            ini_values=ini_val,
                                            burn_in=5000,thin=5,printYes=T)


save(model_Matching_weak,file="Matching_weak_res.Rdata")

