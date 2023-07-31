
###########################################################################
###     Matching Algorithms for Cross-temporal Design                  ###
###########################################################################

library(brms)
library(extraDistr)
library(invgamma)
library(MCMCglmm) 
library(MCMCpack)
library(BayesLogit)
library(mvtnorm)
#library(truncnorm)


### Only update parameters, G known ######
MCMC_run_Hier_matching<-function(Iteration, Data=list(Y,D,Z,X,G), 
                                 prior_values=list(prior_A=10,prior_C=10,prior_N=10, prior_a=1, prior_b=1,hyperparam=list(a_l=1,b_l=1,a_k=1,b_k=1)),
                                 ini_values=list(beta_A0_ini, beta_N0_ini, beta_C0_ini, 
                                                 beta_A_ini, beta_N_ini, beta_C_ini, 
                                                 alpha_A_ini, alpha_C_ini, alpha_N_ini,
                                                 sigma_A_ini, sigma_C_ini, sigma_N_ini,
                                                 mu_alpha_ini,sigma_alpha_ini,Tau_ini),
                                 burn_in,thin,printYes=F){
  
  t=2
  
  if(T){  
    
    # Data 
    Y=Data$Y
    D=Data$D
    Z=Data$Z
    X=Data$X
    G=Data$G
    N=length(Y)
    
    
    # imputations
    Y_imp_04<-matrix(0,nrow=Iteration,ncol=N)
    Y_imp_09<-matrix(0,nrow=Iteration,ncol=N)
    
    
    beta_A0<-rep(0,Iteration)
    beta_C0<-rep(0,Iteration)
    beta_N0<-rep(0,Iteration)
    
    beta_A<-array(0,dim=c(Iteration,ncol(X)))
    beta_C<-array(0,dim=c(Iteration,ncol(X)))
    beta_N<-array(0,dim=c(Iteration,ncol(X)))
    
    beta_A0[1]=ini_values$beta_A0_ini
    beta_C0[1]=ini_values$beta_C0_ini
    beta_N0[1]=ini_values$beta_N0_ini
    
    beta_A[1,]=as.matrix(ini_values$beta_A_ini,ncol=1)
    beta_C[1,]=as.matrix(ini_values$beta_C_ini,ncol=1)
    beta_N[1,]=as.matrix(ini_values$beta_N_ini,ncol=1)
    
    
    
    sigma2<-rep(10,Iteration)
    sigma_A<-rep(10,Iteration)
    sigma_C<-rep(10,Iteration)
    sigma_N<-rep(10,Iteration)
    
    alpha_A<-rep(0,Iteration)
    alpha_C<-rep(0,Iteration)
    alpha_N<-rep(0,Iteration)
    
    alpha_A[1]=ini_values$alpha_A_ini
    alpha_C[1]=ini_values$alpha_C_ini
    alpha_N[1]=ini_values$alpha_N_ini
    
    mu_alpha<-rep(0,Iteration)
    sigma_alpha<-rep(1,Iteration)
    Tau<-rep(0,Iteration)
    
    mu_alpha[1]<-ini_values$mu_alpha_ini
    sigma_alpha[1]<-ini_values$sigma_alpha_ini
    Tau[1]<-ini_values$Tau_ini
    
    # place holders
    DID<-rep(0,Iteration)
    DID_C04<-rep(0,Iteration)
    DID_C09<-rep(0,Iteration)
    accept_sigma<-0
    
    
    # Priors 
    a_l=prior_values$hyperparam$a_l #prior for mu.alpha.var
    b_l=prior_values$hyperparam$b_l
    a_k=prior_values$hyperparam$a_k #prior for sigma_alpha
    b_k=prior_values$hyperparam$b_k
    
    prior_beta_mu=prior_values$prior_beta_mu
    prior_beta_sd2=prior_values$prior_beta_sd2
    
    prior_a=prior_values$prior_a
    prior_b=prior_values$prior_b
    prior_A=prior_values$prior_A
    prior_C=prior_values$prior_C
    prior_N=prior_values$prior_N
    prior_tau_mu=prior_values$prior_tau_mu
    prior_tau_sigma=prior_values$prior_tau_sigma
    prior_A_scale=prior_values$A
    proposal_width=prior_values$proposal_width
    
    id_A=which(G=="A")
    id_C=which(G=="C")
    id_N=which(G=="N")
    
    Y_A=Y[id_A]
    Y_C=Y[id_C]
    Y_N=Y[id_N]
    Yj=list(Y_A=Y_A,Y_N=Y_N)
    
    Z_A=Z[id_A]
    Z_C=Z[id_C]
    Z_N=Z[id_N]
    Zj=list(Z_A=Z_A,Z_N=Z_N)
    
    N_A=length(Y_A)
    N_C=length(Y_C)
    N_N=length(Y_N)
    
    X_A<-X[id_A,]
    X_C<-X[id_C,]
    X_N<-X[id_N,]
    
    mu_A_ini=mean(Y_A)
    mu_C_ini=mean(Y_C)
    mu_N_ini=mean(Y_N)
    
  }
  
  while(t <= Iteration){
    
    if(printYes==T){
      if(t%%5000==0){
        print(paste("Iteration=",t,sep=""))
      }
    }
    
    X.mat<-cbind(rep(1),X)
    q<-ncol(X.mat)
    Bigmat<-as.matrix(cbind(1*(G=="A")*X.mat,1*(G=="N")*X.mat,1*(G=="C")*X.mat,1*(G=="C")*D))
    colnames(Bigmat)<-c(paste("beta",c(1:q),"A",sep="_"),paste("beta",c(1:q),"N",sep="_"),paste("beta",c(1:q),"C",sep="_"),"tau")
    p<-ncol(Bigmat)
    beta_ini<-matrix(rep(0,p),ncol=1) #prior mu
    T0_beta<-diag(10000,p) #prior sigma
    
    ## calculate Ystar=Y-Z*temporal coef
    Y_star<-Y
    Y_star[id_A]<-Y[id_A]-Z[id_A]*alpha_A[t-1]
    Y_star[id_N]<-Y[id_N]-Z[id_N]*alpha_N[t-1]
    Y_star[id_C]<-Y[id_C]-Z[id_C]*alpha_C[t-1]
    
    ## update betas ##
    inv.Sigma0<-solve(T0_beta)
    M<-solve(inv.Sigma0+1/sigma2[t-1]*t(Bigmat)%*%Bigmat)%*%(inv.Sigma0%*%beta_ini+1/sigma2[t-1]*t(Bigmat)%*%Y_star)
    V<-solve(inv.Sigma0+1/sigma2[t-1]*t(Bigmat)%*%Bigmat)
    beta<-c(rmvnorm(1,M,V))
    beta_A0[t]<-beta[1]
    beta_A[t,]<-beta[2:3]
    beta_N0[t]<-beta[4]
    beta_N[t,]<-beta[5:6]
    beta_C0[t]<-beta[7]
    beta_C[t,]<-beta[8:9]
    Tau[t]<-beta[10]
    
    mu_alpha_var_prior=100
    
    ### Update temporal coefficients ####
    ## Update  mu_alpha ## (M=4, four groups)
    inv.sigma_alpha2<-solve(sigma_alpha[t-1])
    alpha.sum<-alpha_A[t-1]+alpha_N[t-1]+alpha_C[t-1]
    
    mu.alpha<-solve(3*inv.sigma_alpha2+1/mu_alpha_var_prior)*inv.sigma_alpha2*alpha.sum
    sigma.alpha<-solve(3*inv.sigma_alpha2+1/mu_alpha_var_prior)
    mu_alpha[t]<-rmvnorm(1,mean=mu.alpha,sigma = sigma.alpha)
    
    Y_star<-Y-Bigmat%*%matrix(beta,ncol=1)
    
    temp.A<-solve(t(Z_A)%*%Z_A/sigma2[t-1]+inv.sigma_alpha2)%*%(t(Z_A)%*%Y_star[id_A]/sigma2[t-1]+inv.sigma_alpha2*mu_alpha[t])
    temp.N<-solve(t(Z_N)%*%Z_N/sigma2[t-1]+inv.sigma_alpha2)%*%(t(Z_N)%*%Y_star[id_N]/sigma2[t-1]+inv.sigma_alpha2*mu_alpha[t])
    temp.C<-solve(t(Z_C)%*%Z_C/sigma2[t-1]+inv.sigma_alpha2)%*%(t(Z_C)%*%Y_star[id_C]/sigma2[t-1]+inv.sigma_alpha2*mu_alpha[t])
    
    temp.sigma.A<-1/(t(Z_A)%*%Z_A/sigma2[t-1]+inv.sigma_alpha2) 
    temp.sigma.N<-1/(t(Z_N)%*%Z_N/sigma2[t-1]+inv.sigma_alpha2)
    temp.sigma.C<-1/(t(Z_C)%*%Z_C/sigma2[t-1]+inv.sigma_alpha2)
    
    alpha_A[t]<-rnorm(1,mean=temp.A,sd=sqrt(temp.sigma.A))
    alpha_N[t]<-rnorm(1,mean=temp.N,sd=sqrt(temp.sigma.N))
    alpha_C[t]<-rnorm(1,mean=temp.C,sd=sqrt(temp.sigma.C))
    
    
    ### update Sigma_alpha
    sigma_alpha[t]<-invgamma::rinvgamma(1,shape=a_k+3/2,rate=b_k+1/2*
                                          ((alpha_A[t]-mu_alpha[t])^2+(alpha_N[t]-mu_alpha[t])^2+(alpha_C[t]-mu_alpha[t])^2))
    
    Y.res<-Y
    Y.res[id_A]<-Y[id_A]-Bigmat[id_A,]%*%matrix(beta,ncol=1)-Z_A*alpha_A[t]
    Y.res[id_N]<-Y[id_N]-Bigmat[id_N,]%*%matrix(beta,ncol=1)-Z_N*alpha_N[t]
    Y.res[id_C]<-Y[id_C]-Bigmat[id_C,]%*%matrix(beta,ncol=1)-Z_C*alpha_C[t]
    
    PP<-crossprod(Y.res)
    post_a<-length(Y)/2+prior_a
    post_b<-PP/2+prior_b
    sigma2[t]<-invgamma::rinvgamma(1,shape=post_a,rate=post_b)
    
    sigma_A[t]<-sigma_N[t]<-sigma_C[t]<-sigma2[t]
    
    ### imputation ##
    mu_C0_04<-beta_C0[t]+X_C%*%as.matrix(beta_C[t,],ncol=1)
    mu_C1_04<-beta_C0[t]+X_C%*%as.matrix(beta_C[t,],ncol=1)+Tau[t]
    
    mu_C0_09<-beta_C0[t]+X_C%*%as.matrix(beta_C[t,],ncol=1)+alpha_C[t]
    mu_C1_09<-beta_C0[t]+X_C%*%as.matrix(beta_C[t,],ncol=1)+alpha_C[t]+Tau[t]
    
    mu_A.p<-beta_A0[t]+X%*%as.matrix(beta_A[t,],ncol=1)
    mu_N.p<-beta_N0[t]+X%*%as.matrix(beta_N[t,],ncol=1)
    mu_C.p<-beta_C0[t]+X%*%as.matrix(beta_C[t,],ncol=1)
    

      Y_imp_C0_04<-rnorm(N_C,mu_C0_04,sd=sqrt(sigma_C[t]))
      Y_imp_C1_04<-rnorm(N_C,mu_C1_04,sd=sqrt(sigma_C[t]))
      
      Y_imp_C0_09<-rnorm(N_C,mu_C0_09,sd=sqrt(sigma_C[t]))
      Y_imp_C1_09<-rnorm(N_C,mu_C1_09,sd=sqrt(sigma_C[t]))
      
      
      mu_pred_04=mu_A.p*(G=="A")+ mu_C.p*(G=="C")+ mu_N.p*(G=="N")  
      mu_pred_09=(mu_A.p)*(G=="A")+ (mu_C.p+Tau[t])*(G=="C")+ (mu_N.p)*(G=="N") 
      sigma_pred=sigma_A[t]*(G=="A")+ sigma_C[t]*(G=="C")+ sigma_N[t]*(G=="N") 
      
      Y_imp_04[t,]<-rnorm(N,mu_pred_04,sd=sqrt(sigma_pred))
      Y_imp_09[t,]<-rnorm(N,mu_pred_09,sd=sqrt(sigma_pred))
      
      DID[t]=mean(Y_imp_09[t,which(G=="C")])-mean(Y_imp_04[t,which(G=="C")])#-alpha_C[t]
      DID_C04[t]<-mean(Y_imp_C1_04)-mean(Y_imp_C0_04)
      DID_C09[t]<-mean(Y_imp_C1_09)-mean(Y_imp_C0_09)
      

    
    t=t+1
    
    
  }
  
  
  MCMC_list<-data.frame(alpha_A=alpha_A, alpha_N=alpha_N, alpha_C=alpha_C,
                        mu_alpha=mu_alpha, sigma_alpha=sigma_alpha,
                        Tau=Tau,beta0_A=beta_A0,beta0_C=beta_C0,beta0_N=beta_N0, 
                        beta_A=beta_A,beta_C=beta_C,beta_N=beta_N, 
                        sigma_A=sigma_A,sigma_N=sigma_N,sigma_C=sigma_C)
  
  return(list(param=MCMC_list[seq(burn_in,Iteration,thin),],
              DID=DID[seq(burn_in,Iteration,thin)],
              DID_C04=DID_C04[seq(burn_in,Iteration,thin)],
              DID_C09=DID_C09[seq(burn_in,Iteration,thin)]
  ))
  
}


MCMC_run_NH_matching<-function(Iteration, Data=list(Y,D,Z,X,G),
                               prior_values=list(prior_mu,prior_sd2,prior_beta_mu,prior_beta_sd2,prior_a,prior_b,prior_A,prior_C,prior_N),
                               ini_values=list(mu_A_ini, mu_C_ini, mu_N_ini, 
                                               beta0_A_ini, beta0_C_ini, beta0_N_ini, 
                                               beta_A_ini, beta_C_ini, beta_N_ini, 
                                               psi_A_ini,psi_C_ini,
                                               alpha_ini, alpha_C_ini,
                                               sigma_A_ini, sigma_C_ini, sigma_N_ini),
                               burn_in,thin,printYes=F){
  
  t=1
  # Data 
  {
    Y=Data$Y
    D=Data$D
    Z=Data$Z
    X=Data$X
    G<-Data$G
    N=length(Y)
    
    
    # imputations
    Y_imp_04<-matrix(0,nrow=Iteration,ncol=N)
    Y_imp_09<-matrix(0,nrow=Iteration,ncol=N)
    DID<-rep(0,Iteration)
    DID_C04<-rep(0,Iteration)
    DID_C09<-rep(0,Iteration)
    
    # initial values
    mu_A_ini=ini_values$mu_A_ini
    mu_C_ini=ini_values$mu_C_ini
    mu_N_ini=ini_values$mu_N_ini
    
    alpha_ini=ini_values$alpha_ini
    alpha_C_ini=ini_values$alpha_C_ini
    
    sigma_A_ini=ini_values$sigma_A_ini
    sigma_C_ini=ini_values$sigma_C_ini
    sigma_N_ini=ini_values$sigma_N_ini
    
    
    
    beta_A0_ini=ini_values$beta0_A_ini
    beta_C0_ini=ini_values$beta0_C_ini
    beta_N0_ini=ini_values$beta0_N_ini
    
    beta_A_ini=as.matrix(ini_values$beta_A_ini,ncol=1)
    beta_C_ini=as.matrix(ini_values$beta_C_ini,ncol=1)
    beta_N_ini=as.matrix(ini_values$beta_N_ini,ncol=1)
    
    prior_beta_mu=prior_values$prior_beta_mu
    prior_beta_sd2=prior_values$prior_beta_sd2
    

    # place holders
    beta_A0<-rep(0,Iteration)
    beta_C0<-rep(0,Iteration)
    beta_N0<-rep(0,Iteration)
    
    beta_A<-array(0,dim=c(Iteration,ncol(X)))
    beta_C<-array(0,dim=c(Iteration,ncol(X)))
    beta_N<-array(0,dim=c(Iteration,ncol(X)))
    
    sigma_A<-rep(1,Iteration)
    sigma_C<-rep(1,Iteration)
    sigma_N<-rep(1,Iteration)
    
    alpha<-rep(0,Iteration)
    alpha_C<-rep(0,Iteration)
    
    # Priors 
    prior_mu=prior_values$prior_mu
    prior_sd2=prior_values$prior_sd2
    prior_a=prior_values$prior_a
    prior_b=prior_values$prior_b
    prior_A=prior_values$prior_A
    prior_C=prior_values$prior_C
    prior_N=prior_values$prior_N
    
    id_A=which(G=="A")
    id_C=which(G=="C")
    id_N=which(G=="N")
    
    Y_A=Y[id_A]
    Y_C=Y[id_C]
    Y_N=Y[id_N]
    Yj=list(Y_A=Y_A,Y_N=Y_N,Y_C=Y_C)
    
    Z_A=Z[id_A]
    Z_C=Z[id_C]
    Z_N=Z[id_N]
    Zj=list(Z_A=Z_A,Z_N=Z_N,Z_C=Z_C)
    
    N_A=length(Y_A)
    N_C=length(Y_C)
    N_N=length(Y_N)
    
    X_A<-X[id_A,]
    X_C<-X[id_C,]
    X_N<-X[id_N,]
    X.mat<-cbind(rep(1),X)
    
    
  }
  
  while(t <= Iteration){
    
    if(printYes==T){
      if(t%%1000==0){
        print(paste("Iteration=",t,sep=""))
      }
    }
    
    
    if(t==1){
      
 
          ## Update coefficients 
          ## Big design matrix ##
          Bigmat<-as.matrix(cbind(1*(G=="A")*X.mat,1*(G=="N")*X.mat,1*(G=="C")*cbind(X.mat,Z),Z))
          colnames(Bigmat)<-c("beta0A","beta1A","beta2A","beta0N","beta1N","beta2N","beta0C","beta1C","beta2C","alpha_C","alpha")
          p<-ncol(Bigmat)
          beta_ini<-rep(0,p)
          mu_beta<-rep(0,p)
          T0_beta<-diag(0.001,p)
          
          ## update betas ##
          v<-solve(crossprod(Bigmat)*1/sigma_A_ini+T0_beta) 
          m<-v%*%(T0_beta%*%mu_beta+1/sigma_A_ini*t(Bigmat)%*%Y) 
          beta<-c(rmvnorm(1,m,v))
          
          beta_A0[t]<-beta[1]
          beta_A[t,]<-beta[2:3]
          beta_N0[t]<-beta[4]
          beta_N[t,]<-beta[5:6]
          beta_C0[t]<-beta[7]
          beta_C[t,]<-beta[8:9]
          alpha_C[t]<-beta[10]
          alpha[t]<-beta[11]
          
          beta<-matrix(beta,ncol=1)
          vec.y<-Y-Bigmat%*%beta
          
          P<-t(vec.y)%*%vec.y
          post_a<-N/2+prior_a
          post_b<-P/2+prior_b
          sigma_A[t]<-extraDistr::rinvgamma(1,alpha=post_a,beta=post_b)
          sigma_N[t]<-sigma_C[t]<-sigma_A[t]
          
        
          }else{
      
      
    
          Bigmat<-as.matrix(cbind(1*(G=="A")*X.mat,1*(G=="N")*X.mat,1*(G=="C")*cbind(X.mat,Z),Z))
          colnames(Bigmat)<-c("beta0A","beta1A","beta2A","beta0N","beta1N","beta2N","beta0C","beta1C","beta2C","alpha_C","alpha")
          p<-ncol(Bigmat)
          mu_beta<-rep(0,p)
          T0_beta<-diag(0.001,p)
          
          v<-solve(crossprod(Bigmat)*1/sigma_A[t-1]+T0_beta) 
          m<-v%*%(T0_beta%*%mu_beta+1/sigma_A[t-1]*t(Bigmat)%*%Y) 
          beta<-c(rmvnorm(1,m,v))
          
          
          beta_A0[t]<-beta[1]
          beta_A[t,]<-beta[2:3]
          beta_N0[t]<-beta[4]
          beta_N[t,]<-beta[5:6]
          beta_C0[t]<-beta[7]
          beta_C[t,]<-beta[8:9]
          alpha_C[t]<-beta[10]
          alpha[t]<-beta[11]
          
          beta<-matrix(beta,ncol=1)
          vec.y<-Y-Bigmat%*%beta
          
          P<-t(vec.y)%*%vec.y
          post_a<-N/2+prior_a
          post_b<-P/2+prior_b
          sigma_A[t]<-extraDistr::rinvgamma(1,alpha=post_a,beta=post_b)
          sigma_N[t]<-sigma_C[t]<-sigma_A[t]

      
    }
    
    ### imputation ##
    mu_C0_04<-beta_C0[t]+X_C%*%as.matrix(beta_C[t,],ncol=1)
    mu_C1_04<-beta_C0[t]+X_C%*%as.matrix(beta_C[t,],ncol=1)+alpha_C[t]
    
    mu_C0_09<-beta_C0[t]+X_C%*%as.matrix(beta_C[t,],ncol=1)+alpha[t]
    mu_C1_09<-beta_C0[t]+X_C%*%as.matrix(beta_C[t,],ncol=1)+alpha[t]+alpha_C[t]
    
    mu_A.p<-beta_A0[t]+X%*%as.matrix(beta_A[t,],ncol=1)
    mu_N.p<-beta_N0[t]+X%*%as.matrix(beta_N[t,],ncol=1)
    mu_C.p<-beta_C0[t]+X%*%as.matrix(beta_C[t,],ncol=1)
    
  
      Y_imp_C0_04<-rnorm(N_C,mu_C0_04,sd=sqrt(sigma_C[t]))
      Y_imp_C1_04<-rnorm(N_C,mu_C1_04,sd=sqrt(sigma_C[t]))
      
      Y_imp_C0_09<-rnorm(N_C,mu_C0_09,sd=sqrt(sigma_C[t]))
      Y_imp_C1_09<-rnorm(N_C,mu_C1_09,sd=sqrt(sigma_C[t]))
      
      
      mu_pred_04=mu_A.p*(G=="A")+ mu_C.p*(G=="C")+ mu_N.p*(G=="N")  
      mu_pred_09=(mu_A.p)*(G=="A")+ (mu_C.p+alpha_C[t])*(G=="C")+ (mu_N.p)*(G=="N") 
      sigma_pred=sigma_A[t]*(G=="A")+ sigma_C[t]*(G=="C")+ sigma_N[t]*(G=="N") 
      
      Y_imp_04[t,]<-rnorm(N,mu_pred_04,sd=sqrt(sigma_pred))
      Y_imp_09[t,]<-rnorm(N,mu_pred_09,sd=sqrt(sigma_pred))
      

      DID[t]=mean(Y_imp_09[t,which(G=="C")])-mean(Y_imp_04[t,which(G=="C")])#-alpha_C[t]
      DID_C04[t]<-mean(Y_imp_C1_04)-mean(Y_imp_C0_04)
      DID_C09[t]<-mean(Y_imp_C1_09)-mean(Y_imp_C0_09)

  
    
    t=t+1
    
  }
  

    MCMC_list<-data.frame(alpha=alpha, alpha_C=alpha_C,
                          beta_A0=beta_A0,beta_C0=beta_C0,beta_N0=beta_N0,
                          beta_A=beta_A,beta_C=beta_C,beta_N=beta_N,
                          sigma_A=sigma_A,sigma_N=sigma_N,sigma_C=sigma_C)

  
  return(list(param=MCMC_list[seq(burn_in,Iteration,thin),],
              DID=DID[seq(burn_in,Iteration,thin)],
              DID_C04=DID_C04[seq(burn_in,Iteration,thin)],
              DID_C09=DID_C09[seq(burn_in,Iteration,thin)]))

}


################ Function for the matching algorithm ############33333
#Temporal.assumption: Parallel, Common, Weaker
Cross_mathching<-function(data,Iteration,Temporal.assumption="Parallel",
                          prior_values=list(prior_beta_mu,prior_beta_sd2,
                                            prior_a,prior_b,
                                            prior_A,prior_C,prior_N,
                                            prior_tau_mu,prior_tau_sigma,
                                            A, proposal_width),
                          burn_in,thin,printYes=F){
  
  #### Cross-temporal matching ####
  {
  ### Step 1: Match people from 2004 in hospice to people in 2009 use hospice ##
  ##########  S1_2004: Z=0, D=1  ---> S1_2009: Z=1, D=1 #######################
  data_2004<-data[which(data$Z==0),]
  data_2009<-data[which(data$Z==1),]
  
  ## propensity score model using cohort in 2004 ##
  PS2004_model<-glm(D~X1+X2,data=data_2004,family = binomial())
  
  ## Apply propensity score model on cohort 2009 ##
  Pred_2009<-predict(PS2004_model,newdata = data_2009,type="response")
  Pred_2004<-predict(PS2004_model,type="response")
  
  ## one-to-one match group 2004 with D=1 and 2009 with D=1 with replacement ##
  N01<-length(which(data_2004$D==1))
  data_2004_1<-data_2004[which(data_2004$D==1),]
  data_2009_1<-data_2009[which(data_2009$D==1),]
  
  Pred_2004_1<-Pred_2004[which(data_2004$D==1)]
  Pred_2009_1<-Pred_2009[which(data_2009$D==1)]
  
  data_2004_1$Matched_ID<-rep(0)
  
  for(i in 1:N01){
    
    dist_i<-abs(Pred_2004_1[i]-Pred_2009_1)  
    
    if(min(dist_i)<=0.01){
      matched_id<-which.min(dist_i) 
    }else{
      matched_id<-NA 
    }
    
    data_2004_1$Matched_ID[i]<-matched_id
    
  }
  
  
  ##remove non-matched ##
  data_2004_1<-na.omit(data_2004_1)
  data_2009_1_matched<-data_2009_1[data_2004_1$Matched_ID,]
  
  data_always_takers<-rbind(data_2004_1[,-7],data_2009_1_matched)
  data_always_takers$G_pred<-"Always_taker"  
  
  ##### Step 2: use 2009 propensity score model to match non-hospice decedent in 2009 ###
  ## propensity score model using cohort in 2009 ##
  PS2009_model<-glm(D~X1+X2,data=data_2009,family = binomial())
  
  ## Apply propensity score model on cohort 2009 ##
  Pred_2004<-predict(PS2009_model,newdata = data_2004,type="response")
  Pred_2009<-predict(PS2009_model,type="response")
  
  ## one-to-one match group 2004 with D=1 and 2009 with D=1 with replacement ##
  N10<-length(which(data_2009$D==0))
  data_2004_0<-data_2004[which(data_2004$D==0),]
  data_2009_0<-data_2009[which(data_2009$D==0),]
  
  Pred_2009_0<-Pred_2009[which(data_2009$D==0)]
  Pred_2004_0<-Pred_2004[which(data_2004$D==0)]
  
  data_2009_0$Matched_ID<-rep(0)
  
  for(i in 1:N10){
    
    dist_i<-abs(Pred_2009_0[i]-Pred_2004_0)  
    
    if(min(dist_i)<=0.01){
      matched_id<-which.min(dist_i) 
    }else{
      matched_id<-NA 
    }
    
    data_2009_0$Matched_ID[i]<-matched_id
    
  }
  
  ##remove non-matched ##
  data_2009_0<-na.omit(data_2009_0)
  data_2004_0_matched<-data_2004_0[data_2009_0$Matched_ID,]
  
  data_never_takers<-rbind(data_2009_0[,-7],data_2004_0_matched)
  data_never_takers$G_pred<-"Never_taker"  
  
  
  ######Step 2: Left unmatched group 2009 who use hospice are set to compliers ###
  data_2009_unmatched<-data_2009_1[-data_2004_1$Matched_ID,]

  Pred_2009_unmatched<-predict(PS2009_model,newdata=data_2009_unmatched,type="response")
  
  N<-nrow(data_2009_unmatched)
  
  for(i in 1:N){
    
    dist_i<-abs(Pred_2009_unmatched[i]-Pred_2004_0)  
    
    if(min(dist_i)<=0.01){
      matched_id<-which.min(dist_i) 
    }else{
      matched_id<-NA 
    }
    
    data_2009_unmatched$Matched_ID[i]<-matched_id
    
    
  }
  
  
  
  ##remove non-matched ##
  data_2009_unmatched<-na.omit(data_2009_unmatched)
  data_2004_0_matched_c<-data_2004_0[data_2009_unmatched$Matched_ID,]
  
  data_compliers<-rbind(data_2009_unmatched[,-7],data_2004_0_matched_c)
  data_compliers$G_pred<-"Compliers" 
  
  ########## Form dataset for estimating causal estimand ######
  data_model<-rbind(data_always_takers,data_never_takers,data_compliers)
  
  }
  
  ### Matching with DID under parallel trends assumption ######
  if(Temporal.assumption=="Parallel"){
    
    new.data<-rbind(data_compliers,data_never_takers)
    new.data$G2<-ifelse(new.data$G_pred=="Compliers",1,0)
    new.data$Post<-new.data$Z
    
    fit.DID<-lm(Y~G2+Post+X1+X2+G2*Post,data=new.data)
    DID<-coef(fit.DID)[6]
    CI<-confint(fit.DID)[6,]
    
    Result<-list(CACE=DID,CACE_CI=CI,cat.number=cat.number)
  }
  
  ### Matching with DA under common trends assumption ######
  if(Temporal.assumption=="Common"){
    
    data_PO<-data_model[,-6]
    data_PO$G_pred<-factor(data_PO$G_pred,levels = c("Always_taker","Compliers","Never_taker"),labels=c("A","C","N"))
    Datalist=list(Y=data_PO$Y, D=data_PO$D, Z=data_PO$Z,X=cbind(data_PO$X1,data_PO$X2),G=data_PO$G_pred)
    
      
      prior_val=list(prior_mu=0,prior_sd2=10000,prior_beta_mu=0,prior_beta_sd2=10000,prior_a=1,prior_b=1,prior_A=10,prior_C=10,prior_N=10)
      
      fit.A<-lm(Y~X1+X2+Z,data=data_PO[which(data_PO$G_pred=="A"),])
      fit.N<-lm(Y~X1+X2+Z,data=data_PO[which(data_PO$G_pred=="N"),])
      fit.C<-lm(Y~X1+X2+Z,data=data_PO[which(data_PO$G_pred=="C"),])
      
      ini_val=list(w_A_ini=1/3,w_C_ini=1/3,w_N_ini=1/3,
                   mu_A_ini=mean(data_PO$Y[which(data_PO$G_pred=="A")]), 
                   mu_C_ini=mean(data_PO$Y[which(data_PO$G_pred=="C")]), 
                   mu_N_ini=mean(data_PO$Y[which(data_PO$G_pred=="N")]),
                   beta0_A_ini=coef(fit.A)[1], 
                   beta0_C_ini=coef(fit.C)[1], 
                   beta0_N_ini=coef(fit.N)[1], 
                   beta_A_ini=coef(fit.A)[2:3],beta_C_ini=coef(fit.C)[2:3],beta_N_ini=coef(fit.N)[2:3],
                   psi_A_ini=c(1,2), psi_C_ini=c(1,2), 
                   alpha_ini=(coef(fit.A)[4]+coef(fit.N)[4])/2,
                   alpha_C_ini=coef(fit.C)[4],
                   sigma_A_ini=10, sigma_C_ini=10, sigma_N_ini=10) 
      
      MCMC_model<-MCMC_run_NH_matching(Iteration=Iteration, Data=Datalist,Intercept.model = F,Outcome.type=Outcome.type,
                                       prior_values=prior_val,
                                       ini_values=ini_val,
                                       burn_in=burn_in,printYes=printYes)
      tau<-MCMC_model$param$alpha_C
      DID<-MCMC_model$DID
      DID_C04<-MCMC_model$DID_C04
      DID_C09<-MCMC_model$DID_C09
      
      CACE.all<-DID
      CACE<-mean(CACE.all,na.rm=T)
      CI<-quantile(CACE.all,c(0.025,0.975),na.rm=T)
      
      CACE.all.04<-DID_C04
      CACE.04<-mean(CACE.all.04,na.rm=T)
      CI.04<-quantile(CACE.all.04,c(0.025,0.975),na.rm=T)
      
      CACE.all.09<-DID_C09
      CACE.09<-mean(CACE.all.09,na.rm=T)
      CI.09<-quantile(CACE.all.09,c(0.025,0.975),na.rm=T)
      
      Result<-list(Tau=tau,CACE=CACE,CACE_CI=CI,
                   CACE.04=CACE.04,CACE_CI.04=CI.04,
                   CACE.09=CACE.09,CACE_CI.09=CI.09,
                   cat.number=cat.number)
    
  }
  
  
  ### Matching with DA under weaker trends assumption ######
  if(Temporal.assumption=="Weaker"){
    
    data_PO<-data_model[,-6]
    data_PO$G_pred<-factor(data_PO$G_pred,levels = c("Always_taker","Compliers","Never_taker"),labels=c("A","C","N"))
    Datalist=list(Y=data_PO$Y, D=data_PO$D, Z=data_PO$Z,X=cbind(data_PO$X1,data_PO$X2),X1=data_PO$X1,X2=data_PO$X2,G=data_PO$G_pred)
    
    fit.A<-lm(Y~X1+X2+Z,data=data_PO[which(data_PO$G_pred=="A"),])
    fit.N<-lm(Y~X1+X2+Z,data=data_PO[which(data_PO$G_pred=="N"),])
    fit.C<-lm(Y~X1+X2+Z,data=data_PO[which(data_PO$G_pred=="C"),])
    
    ini_val=list(mu_A_ini=mean(data_PO$Y[which(data_PO$G_pred=="A")]), 
                 mu_C_ini=mean(data_PO$Y[which(data_PO$G_pred=="C")]), 
                 mu_N_ini=mean(data_PO$Y[which(data_PO$G_pred=="N")]),
                 beta0_A_ini=coef(fit.A)[1], 
                 beta0_C_ini=coef(fit.C)[1], 
                 beta0_N_ini=coef(fit.N)[1], 
                 beta_A_ini=coef(fit.A)[2:3],beta_C_ini=coef(fit.C)[2:3],beta_N_ini=coef(fit.N)[2:3],
                 psi_A_ini=c(1,2), psi_C_ini=c(1,2), 
                 alpha_A_ini=coef(fit.A)[4],
                 alpha_N_ini=coef(fit.N)[4],
                 alpha_C_ini=coef(fit.C)[4],
                 sigma_A_ini=10, sigma_C_ini=10, sigma_N_ini=10,
                 mu_alpha_ini=4,sigma_alpha_ini=10,Tau_ini=1)
    
    
    MCMC_model<-MCMC_run_Hier_matching(Iteration=Iteration, Data=Datalist, Outcome.type=Outcome.type,Intercept.model = F,
                                       prior_values=prior_values,
                                       ini_values=ini_val,
                                       burn_in=burn_in,printYes=printYes)
    tau<-MCMC_model$param$Tau
    DID<-MCMC_model$DID
    DID_C04<-MCMC_model$DID_C04
    DID_C09<-MCMC_model$DID_C09

    CACE.all<-DID
    CACE<-mean(CACE.all,na.rm=T)
    CI<-quantile(CACE.all,c(0.025,0.975),na.rm=T)
    
    CACE.all.04<-DID_C04
    CACE.04<-mean(CACE.all.04,na.rm=T)
    CI.04<-quantile(CACE.all.04,c(0.025,0.975),na.rm=T)
    
    CACE.all.09<-DID_C09
    CACE.09<-mean(CACE.all.09,na.rm=T)
    CI.09<-quantile(CACE.all.09,c(0.025,0.975),na.rm=T)

    Result<-list(Tau=tau,CACE=CACE,CACE_CI=CI,
                 CACE.04=CACE.04,CACE_CI.04=CI.04,
                 CACE.09=CACE.09,CACE_CI.09=CI.09,
                 cat.number=cat.number)
  }
  

  return(Result)
  
}

### Example: Matching with DID under parallel trends assumption ######
model_matching_did<-Cross_mathching(data=data,Iteration=Iteration,Temporal.assumption = "Parallel",
                                   prior_values=prior_val,
                                   burn_in,thin=thin,printYes=T)

### Example: Matching with DA under common trends assumption ######
prior_val=list(prior_mu=0,prior_sd2=10000,prior_beta_mu=0,prior_beta_sd2=10000,prior_a=1,prior_b=1,prior_A=10,prior_C=10,prior_N=10)

model_matching_common<-Cross_mathching(data=data,Temporal.assumption = "Common",
                                 Iteration=Iteration,
                                 prior_values=prior_val,
                                 burn_in,thin=thin,printYes=F)

### Example: Matching with DA under weaker trends assumption ######
prior_val= prior_values=list(prior_A=10,prior_C=10,prior_N=10, prior_a=1, prior_b=1,hyperparam=list(a_l=1,b_l=1,a_k=1,b_k=1))

matching_weak<-Cross_mathching(data=data,Iteration=Iteration,Temporal.assumption = "Weaker",
                               prior_values=prior_val,
                               burn_in=burn_in,thin=thin.in,printYes=T)

