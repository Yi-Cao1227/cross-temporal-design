
################### Topic2: Simulation Final for Binary outcome #####################
##################################################################

library(brms)
library(extraDistr)
library(invgamma)
library(MCMCglmm) 
library(MCMCpack)
library(mvtnorm)
library(maxLik)
library(BayesLogit)
library(rootSolve)


#library(truncnorm)


######### Functions file ################
# sd_alpha^2/sd_A^2=var.Ratio
invlogit<-function(x){
  
  exp(x)/(1+exp(x))
}

logit<-function(p){
  
  log(p/(1-p))
}



######### Matching Method ###############
PS_matching<-function(data,replace.sample=F){
  
  
  ### Step 1: Match people from 2004 in hospice to people in 2009 use hospice ##
  ##########  S1_2004: Z=0, D=1  ---> S1_2009: Z=1, D=1 #######################
  data_2004<-data[which(data$Z==0),]
  data_2009<-data[which(data$Z==1),]
  
  ## propensity score model using cohort in 2004 ##
  PS2004_model<-glm(D~X1+X2,data=data_2004,family = binomial())
  PS2009_model<-glm(D~X1+X2,data=data_2009,family = binomial())
  
  ## Apply propensity score model on cohort 2009 ##
  Pred_2009<-predict(PS2004_model,newdata = data_2009,type="response")
  Pred_2004<-predict(PS2004_model,type="response")
  
  
  ## one-to-one match group 2004 with D=1 and 2009 with D=1 with replacement ##
  N01<-length(which(data_2004$D==1))
  data_2004_1<-data_2004[which(data_2004$D==1),]
  data_2009_1<-data_2009[which(data_2009$D==1),]
  
  Pred_2004_1<-Pred_2004[which(data_2004$D==1)]
  Pred_2009_1<-Pred_2009[which(data_2009$D==1)]
  
  
  if(replace.sample==T){
    
    data_2004_1$Matched_ID<-rep(0)
    
    for(i in 1:N01){
      
      dist_i<-abs(Pred_2004_1[i]-Pred_2009_1)  
      
      #if(min(dist_i)<=0.01){
      matched_id<-which.min(dist_i) 
      #}else{
      # matched_id<-NA 
      #}
      
      data_2004_1$Matched_ID[i]<-matched_id
      
      
      
    }
    
    ##remove non-matched ##
    data_2004_1<-na.omit(data_2004_1)
    data_2009_1_matched<-data_2009_1[data_2004_1$Matched_ID,]
    
    data_always_takers<-rbind(data_2004_1[,-7],data_2009_1_matched)
    data_always_takers$G_pred<-"Always_taker" 
    
    data_2009_unmatched<-data_2009_1[-data_2004_1$Matched_ID,] #unique matched id <300
    Pred_2009_unmatched<-predict(PS2009_model,newdata=data_2009_unmatched,type="response")
    
    
    
  }else{
    
    Pred_2009_1.new<-Pred_2009_1
    data_2009_1_matched<-data_2004_1
    
    for(i in 1:N01){
      
      dist_i<-abs(Pred_2004_1[i]-Pred_2009_1.new)  #the length of dist_i decreases at every iteration
      
      #if(min(dist_i)<=0.01){
      
      matched_id<-which.min(dist_i) 
      
      # }else{
      #   
      #   matched_id<-NA 
      #   id.na<-c(id.na,i)
      # }
      
      data_2009_1_matched[i,]<-data_2009_1[matched_id,]
      
      Pred_2009_1.new<-Pred_2009_1.new[-matched_id]
      
      data_2009_1<-data_2009_1[-matched_id,]
      
      
    }
    
    data_always_takers<-rbind(data_2004_1,data_2009_1_matched)
    data_always_takers$G_pred<-"Always_taker" 
    
    data_2009_unmatched<-data_2009_1
    Pred_2009_unmatched<-predict(PS2009_model,newdata=data_2009_unmatched,type="response")
    
    
  }
  
  
  
  ##### Step 2: use 2009 propensity score model to match non-hospice decedent in 2009 ###
  ## propensity score model using cohort in 2009 ##
  
  ## Apply propensity score model on cohort 2009 ##
  Pred_2004<-predict(PS2009_model,newdata = data_2004,type="response")
  Pred_2009<-predict(PS2009_model,type="response")
  
  ## one-to-one match group 2009 with D=0 and 2004 with D=0 with replacement ##
  N10<-length(which(data_2009$D==0))
  data_2004_0<-data_2004[which(data_2004$D==0),]
  data_2009_0<-data_2009[which(data_2009$D==0),]
  
  Pred_2009_0<-Pred_2009[which(data_2009$D==0)]
  Pred_2004_0<-Pred_2004[which(data_2004$D==0)]
  
  if(replace.sample==T){
    
    data_2009_0$Matched_ID<-rep(0)
    
    for(i in 1:N10){
      
      dist_i<-abs(Pred_2009_0[i]-Pred_2004_0)  
      
      # if(min(dist_i)<=0.01){
      matched_id<-which.min(dist_i) 
      #}else{
      #  matched_id<-NA 
      #}
      
      data_2009_0$Matched_ID[i]<-matched_id
      
    }
    
    ##remove non-matched ##
    data_2009_0<-na.omit(data_2009_0)
    data_2004_0_matched<-data_2004_0[data_2009_0$Matched_ID,]
    
    data_never_takers<-rbind(data_2009_0[,-7],data_2004_0_matched)
    data_never_takers$G_pred<-"Never_taker"  
    
    
  }else{
    
    Pred_2004_0.new<-Pred_2004_0
    data_2004_0_matched<-data_2009_0
    
    for(i in 1:N10){
      
      dist_i<-abs(Pred_2009_0[i]-Pred_2004_0.new)  
      
      #if(min(dist_i)<=0.01){
      matched_id<-which.min(dist_i) 
      #}else{
      #  matched_id<-NA 
      #}
      
      #data_2009_0$Matched_ID[i]<-matched_id
      
      data_2004_0_matched[i,]<-data_2004_0[matched_id,]
      
      Pred_2004_0.new<-Pred_2004_0.new[-matched_id]
      
      data_2004_0<-data_2004_0[-matched_id,]
      
    }
    
    data_never_takers<-rbind(data_2009_0,data_2004_0_matched)
    data_never_takers$G_pred<-"Never_taker" 
    
    
    
  }
  
  
  
  ######Step 2: Left unmatched group 2009 who use hospice are set to compliers ###
  N<-nrow(data_2009_unmatched)
  
  if(replace.sample==T){
    
    for(i in 1:N){
      
      dist_i<-abs(Pred_2009_unmatched[i]-Pred_2004_0)  
      
      #if(min(dist_i)<=0.01){
      matched_id<-which.min(dist_i) 
      #}else{
      #  matched_id<-NA 
      #}
      
      data_2009_unmatched$Matched_ID[i]<-matched_id
      
      
    }
    
    ##remove non-matched ##
    data_2009_unmatched<-na.omit(data_2009_unmatched)
    data_2004_0_matched_c<-data_2004_0[data_2009_unmatched$Matched_ID,]
    
    data_compliers<-rbind(data_2009_unmatched[,-7],data_2004_0_matched_c)
    data_compliers$G_pred<-"Compliers" 
    
    
  }else{
    
    Pred_2004_0.new<-Pred_2004_0
    data_2004_0_matched_c<-data_2009_0
    
    for(i in 1:N){
      
      dist_i<-abs(Pred_2009_unmatched[i]-Pred_2004_0.new)  
      
      #if(min(dist_i)<=0.01){
      matched_id<-which.min(dist_i) 
      #}else{
      #  matched_id<-NA 
      #}
      data_2004_0_matched_c[i,]<-data_2004_0[matched_id,]
      
      Pred_2004_0.new<-Pred_2004_0.new[-matched_id]
      
      data_2004_0<-data_2004_0[-matched_id,]
      
      
    }
    
    
    data_compliers<-rbind(data_2009_unmatched,data_2004_0_matched_c)
    data_compliers$G_pred<-"Compliers" 
    
    
  }
  
  
  ########## Form dataset for estimating causal estimand ######
  data_model<-rbind(data_always_takers,data_never_takers,data_compliers)
  data_model$G_pred<-factor(data_model$G_pred,levels = c("Always_taker","Compliers","Never_taker"),labels=c("A","C","N"))
  
  return(data_model)
  
}


MH_mcmc_v2<-function(curr_X,args, accept,proposal_rt_df,proposal_rw_width,proposal_Hessian=NA){
  
  if(args$sample.param=="beta"){
    
    proposed_X=rmvnorm(1,mean=curr_X,sigma=proposal_Hessian)
    proposed_X<-matrix(proposed_X,ncol=1)
    
    a=exp(target_param_v2(proposed_X,args)-target_param_v2(curr_X,args))
    
    if(is.na(a)){
      
      x = curr_X        # otherwise "reject" move, and stay where we are
      accept=accept+0
      
      
      
    }else{
      
      if(runif(1)<a){
        x = proposed_X      # accept move with probabily min(1,A)
        accept=accept+1
        
        
      } else {
        x = curr_X      # otherwise "reject" move, and stay where we are
        accept=accept+0
        
      }   
      
    }
    
  }
  
  
  if(args$sample.param%in% c("alpha","tau","mu_alpha")){
    
    #proposed_X=rt.scaled(1,df=proposal_rt_df,mean=curr_X,sd=1)
    proposed_X=rnorm(1,mean=curr_X,sd=proposal_rw_width)
    
    a=exp(target_param_v2(proposed_X,args)-target_param_v2(curr_X,args))
    
    if(is.na(a)){
      
      x = curr_X        # otherwise "reject" move, and stay where we are
      accept=accept+0
      
      
    }else{
      
      if(runif(1)<a){
        x = proposed_X      # accept move with probabily min(1,A)
        accept=accept+1
        
        
      } else {
        x = curr_X      # otherwise "reject" move, and stay where we are
        accept=accept+0
        
      }   
      
    }
    
  }
  
  if(args$sample.param=="sigma_alpha"){
    
    curr_X<-sqrt(curr_X)
    proposed_X=rnorm(1,mean=curr_X,sd=proposal_rw_width)
    
    a=exp(target_param_v2(abs(proposed_X),args)-target_param_v2(curr_X,args))
    
    if(runif(1)<a){
      x = proposed_X^2      # accept move with probabily min(1,A)
      accept=accept+1
      
      
    } else{
      x = curr_X^2        # otherwise "reject" move, and stay where we are
      accept=accept+0
      
    }  
  }
  
  
  
  
  return(list(new.param=x,accept=accept))
  
}

# target function for logistic regression, beta, alpha, tau #  
target_param_v2<-function(param,args=list(sample.param="beta",Grp="A", temporal.assumption=temporal.assumption,prior.sd,
                                          A,Z,Y,X,G,beta_A,beta_C,beta_N,
                                          alpha_A,alpha_C,alpha_N,alpha,mu_alpha,sigma_alpha,Tau)){
  
  sample.param<-args$sample.param
  Grp=args$Grp
  temporal.assumption<-args$temporal.assumption
  prior.sd=args$prior.sd
  
  if(temporal.assumption=="weak"){
    
    A=args$A
    alpha_A=args$alpha_A
    alpha_C=args$alpha_C
    alpha_N=args$alpha_N
    mu_alpha=args$mu_alpha
    sigma_alpha=args$sigma_alpha
    
    
  }else if(temporal.assumption=="common"){
    
    alpha<-args$alpha  
    
  }
  
  beta_A=matrix(args$beta_A,ncol=1)
  beta_C=matrix(args$beta_C,ncol=1)
  beta_N=matrix(args$beta_N,ncol=1)
  
  Tau=args$Tau
  
  Y=args$Y
  G=args$G
  Z=args$Z
  X=args$X
  
  Y_A=Y[G=="A"]
  Y_C=Y[G=="C"]
  Y_N=Y[G=="N"]
  
  X_A.mat=cbind(1,X[G=="A",])
  X_N.mat=cbind(1,X[G=="N",])
  X_C.mat=cbind(1,X[G=="C",])
  
  Z_A<-Z[G=="A"]
  Z_N<-Z[G=="N"]
  Z_C<-Z[G=="C"]
  
  
  if(temporal.assumption=="weak"){
    
    
    if(sample.param=="beta"){
      
      beta_x<-param
      
      if(Grp=="A"){
        
        lp.A <- X_A.mat%*% beta_x+alpha_A*Z_A
        lp.N <- X_N.mat %*% beta_N+alpha_N*Z_N
        lp.C <- X_C.mat %*% beta_C+alpha_C*Z_C+Tau*Z_C
        
      }else if(Grp=="N"){
        
        lp.N <- X_N.mat%*% beta_x+alpha_N*Z_N
        lp.A <- X_A.mat %*% beta_A+alpha_A*Z_A
        lp.C <- X_C.mat %*% beta_C+alpha_C*Z_C+Tau*Z_C
        
        
      }else if(Grp=="C"){
        
        lp.C <- X_C.mat%*% beta_x+alpha_C*Z_C+Tau*Z_C
        lp.N <- X_N.mat %*% beta_N+alpha_N*Z_N
        lp.A <- X_A.mat %*% beta_A+alpha_A*Z_A
        
        
      }
      
      p.A<-1/(1+exp(-lp.A))
      p.N<-1/(1+exp(-lp.N))
      p.C<-1/(1+exp(-lp.C))
      
      f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A))+dmvnorm(c(beta_x),mean=rep(0,length(beta_x)),sigma=diag(prior.sd,length(beta_x)),log=T)+
        sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N))+
        sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C))
      
      
    }
    
    
    if(sample.param=="tau"){
      
      tau_x<-param
      
      lp.A <- X_A.mat %*% beta_A+alpha_A*Z_A
      lp.N <- X_N.mat %*% beta_N+alpha_N*Z_N
      lp.C <- X_C.mat %*% beta_C+alpha_C*Z_C+tau_x*Z_C
      
      p.A<-1/(1+exp(-lp.A))
      p.N<-1/(1+exp(-lp.N))
      p.C<-1/(1+exp(-lp.C))
      
      f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A))+
        sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N))+
        sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C))+dnorm(tau_x,mean=0,sd=prior.sd,log=T)
      
      
    }
    
    if(sample.param=="alpha"){
      
      alpha_x<-param
      
      if(Grp=="A"){
        lp.A <- X_A.mat %*% beta_A+alpha_x*Z_A
        lp.N <- X_N.mat %*% beta_N+alpha_N*Z_N
        lp.C <- X_C.mat %*% beta_C+alpha_C*Z_C+Tau*Z_C
        
        
      }else if(Grp=="N"){
        lp.A <- X_A.mat %*% beta_A+alpha_A*Z_A
        lp.N <- X_N.mat %*% beta_N+alpha_x*Z_N
        lp.C <- X_C.mat %*% beta_C+alpha_C*Z_C+Tau*Z_C   
      }else{
        lp.A <- X_A.mat %*% beta_A+alpha_A*Z_A
        lp.N <- X_N.mat %*% beta_N+alpha_N*Z_N
        lp.C <- X_C.mat %*% beta_C+alpha_x*Z_C+Tau*Z_C
      }
      
      p.A<-1/(1+exp(-lp.A))
      p.N<-1/(1+exp(-lp.N))
      p.C<-1/(1+exp(-lp.C))
      
      if(Grp=="A"){
        f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A))+dnorm(alpha_x,mean=mu_alpha,sd=sqrt(sigma_alpha),log=T)+
          sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N))+dnorm(alpha_N,mean=mu_alpha,sd=sqrt(sigma_alpha),log=T)+
          sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C))+dnorm(alpha_C,mean=mu_alpha,sd=sqrt(sigma_alpha),log=T)
        
      } 
      if(Grp=="N"){
        f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A))+dnorm(alpha_A,mean=mu_alpha,sd=sqrt(sigma_alpha),log=T)+
          sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N))+dnorm(alpha_x,mean=mu_alpha,sd=sqrt(sigma_alpha),log=T)+
          sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C))+dnorm(alpha_C,mean=mu_alpha,sd=sqrt(sigma_alpha),log=T)
        
        
        
      }
      if(Grp=="C"){
        
        f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A))+dnorm(alpha_A,mean=mu_alpha,sd=sqrt(sigma_alpha),log=T)+
          sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N))+dnorm(alpha_N,mean=mu_alpha,sd=sqrt(sigma_alpha),log=T)+
          sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C))+dnorm(alpha_x,mean=mu_alpha,sd=sqrt(sigma_alpha),log=T)
        
        
      }
    }
    
    
    if(sample.param=="mu_alpha"){
      
      mu_alpha_x<-param
      
      f=dnorm(alpha_A,mean=mu_alpha_x,sd=sqrt(sigma_alpha),log=T)+
        dnorm(alpha_N,mean=mu_alpha_x,sd=sqrt(sigma_alpha),log=T)+
        dnorm(alpha_C,mean=mu_alpha_x,sd=sqrt(sigma_alpha),log=T)+
        dnorm(mu_alpha_x,mean=0,sd=prior.sd,log=T)
      
    }
    
    
    if(sample.param=="sigma_alpha"){
      
      sigma_alpha_x<-param
      
      f=dhcauchy(sigma_alpha_x,sigma=A,log=T)+dnorm(alpha_A,mean=mu_alpha,sd=sigma_alpha_x,log=T)+
        dnorm(alpha_C,mean=mu_alpha,sd=sigma_alpha_x,log=T)+dnorm(alpha_N,mean=mu_alpha,sd=sigma_alpha_x,log=T)
    }
    
    
    
  }else if(temporal.assumption=="common"){
    
    
    
    if(sample.param=="beta"){
      
      beta_x<-param
      
      if(Grp=="A"){
        
        lp.A <- X_A.mat %*% beta_x+alpha*Z_A
        lp.N <- X_N.mat %*% beta_N+alpha*Z_N
        lp.C <- X_C.mat %*% beta_C+alpha*Z_C+Tau*Z_C
        
      }else if(Grp=="N"){
        
        lp.N <- X_N.mat%*%beta_x+alpha*Z_N
        lp.A <- X_A.mat %*% beta_A+alpha*Z_A
        lp.C <- X_C.mat %*% beta_C+alpha*Z_C+Tau*Z_C
        
        
      }else if(Grp=="C"){
        
        lp.C <-X_C.mat%*%beta_x+alpha*Z_C+Tau*Z_C
        lp.N <- X_N.mat %*% beta_N+alpha*Z_N
        lp.A <- X_A.mat %*% beta_A+alpha*Z_A
        
        
      }
      
      p.A<-1/(1+exp(-lp.A))
      p.N<-1/(1+exp(-lp.N))
      p.C<-1/(1+exp(-lp.C))
      
      
      f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A))+dmvnorm(c(beta_x),mean=rep(0,length(beta_x)),sigma=diag(1,length(beta_x)),log=T)+
        sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N))+
        sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C))
      
      
      
    }
    
    
    if(sample.param=="tau"){
      
      tau_x<-param
      
      lp.A <- X_A.mat %*% beta_A+alpha*Z_A
      lp.N <- X_N.mat %*% beta_N+alpha*Z_N
      lp.C <- X_C.mat %*% beta_C+alpha*Z_C+tau_x*Z_C
      
      p.A<-1/(1+exp(-lp.A))
      p.N<-1/(1+exp(-lp.N))
      p.C<-1/(1+exp(-lp.C))
      
      f=sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C))+dnorm(tau_x,mean=0,sd=prior.sd,log=T)
      sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A))+
        sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N))
      
    }
    
    if(sample.param=="alpha"){
      
      alpha_x<-param
      
      lp.A <- X_A.mat %*% beta_A+alpha_x*Z_A
      lp.N <- X_N.mat %*% beta_N+alpha_x*Z_N
      lp.C <- X_C.mat %*% beta_C+alpha_x*Z_C+Tau*Z_C
      
      p.A<-1/(1+exp(-lp.A))
      p.N<-1/(1+exp(-lp.N))
      p.C<-1/(1+exp(-lp.C))
      
      
      f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A))+dnorm(alpha_x,mean=0,sd=prior.sd,log=T)+
        sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N))+
        sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C))
      
      
      
    }
    
  }
  
  
  return(f)
  
}

target_param_expand_v2<-function(param,args=list(sample.param="beta",Grp="A",
                                                 prior.sd,
                                                 prior.gamma.shape,
                                                 prior.gamma.scale,
                                                 prior.epsilon.shape,
                                                 prior.epsilon.scale,
                                                 A,Z,Y,X,G,beta_A,beta_C,beta_N,
                                                 epsilon_A,epsilon_C,epsilon_N,mu_epsilon,sigma_epsilon,
                                                 gamma,Tau)){
  prior.sd=args$prior.sd
  prior.gamma.shape=args$prior.gamma.shape
  prior.gamma.scale=args$prior.gamma.scale
  prior.gamma.mu=args$prior.gamma.mu
  prior.gamma.sd=args$prior.gamma.sd
  prior.epsilon.shape=args$prior.epsilon.shape
  prior.epsilon.scale=args$prior.epsilon.scale
  
  sample.param<-args$sample.param
  Grp=args$Grp
  temporal.assumption<-args$temporal.assumption
  
  A=args$A
  epsilon_A=args$epsilon_A
  epsilon_C=args$epsilon_C
  epsilon_N=args$epsilon_N
  mu_epsilon=args$mu_epsilon
  sigma_epsilon=args$sigma_epsilon
  gamma<-args$gamma
  
  beta_A=matrix(args$beta_A,ncol=1)
  beta_C=matrix(args$beta_C,ncol=1)
  beta_N=matrix(args$beta_N,ncol=1)
  
  Tau=args$Tau
  
  Y=args$Y
  G=args$G
  Z=args$Z
  X=args$X
  
  Y_A=Y[G=="A"]
  Y_C=Y[G=="C"]
  Y_N=Y[G=="N"]
  
  X_A.mat=cbind(1,X[G=="A",])
  X_N.mat=cbind(1,X[G=="N",])
  X_C.mat=cbind(1,X[G=="C",])
  
  N_A=length(Y_A)
  N_C=length(Y_C)
  N_N=length(Y_N)
  
  Z_A<-Z[G=="A"]
  Z_N<-Z[G=="N"]
  Z_C<-Z[G=="C"]
  
  if(sample.param=="beta"){
    
    beta_x<-param
    
    if(Grp=="A"){
      
      lp.A <- X_A.mat%*% beta_x+epsilon_A*Z_A*gamma
      lp.N <- X_N.mat %*% beta_N+epsilon_N*Z_N*gamma
      lp.C <- X_C.mat %*% beta_C+epsilon_C*Z_C*gamma+Tau*Z_C
      
    }else if(Grp=="N"){
      
      lp.N <- X_N.mat%*% beta_x+epsilon_N*Z_N*gamma
      lp.A <- X_A.mat %*% beta_A+epsilon_A*Z_A*gamma
      lp.C <- X_C.mat %*% beta_C+epsilon_C*Z_C*gamma+Tau*Z_C
      
      
    }else if(Grp=="C"){
      
      lp.C <- X_C.mat%*% beta_x+epsilon_C*Z_C*gamma+Tau*Z_C
      lp.N <- X_N.mat %*% beta_N+epsilon_N*gamma*Z_N
      lp.A <- X_A.mat %*% beta_A+epsilon_A*Z_A*gamma
      
      
    }
    
    p.A<-1/(1+exp(-lp.A))
    p.N<-1/(1+exp(-lp.N))
    p.C<-1/(1+exp(-lp.C))
    
    
    
    f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A))+dmvnorm(c(beta_x),mean=rep(0,length(beta_x)),sigma=diag(prior.sd,length(beta_x)),log=T)+
      sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N))+
      sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C))
    
    
  }
  
  
  
  if(sample.param=="tau"){
    
    tau_x<-param
    
    lp.A <- X_A.mat %*% beta_A+epsilon_A*Z_A*gamma
    lp.N <- X_N.mat %*% beta_N+epsilon_N*Z_N*gamma
    lp.C <- X_C.mat %*% beta_C+epsilon_C*Z_C*gamma+tau_x*Z_C
    
    p.A<-1/(1+exp(-lp.A))
    p.N<-1/(1+exp(-lp.N))
    p.C<-1/(1+exp(-lp.C))
    
    
    f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A))+
      sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N))+
      sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C))+dnorm(tau_x,mean=0,sd=prior.sd,log=T)
    
  }
  
  if(sample.param=="epsilon"){
    
    epsilon_x<-param
    sigma_epsilon<-sqrt(sigma_epsilon)
    
    if(Grp=="A"){
      lp.A <- X_A.mat %*% beta_A+epsilon_x*Z_A*gamma
      lp.N <- X_N.mat %*% beta_N+epsilon_N*Z_N*gamma
      lp.C <- X_C.mat %*% beta_C+epsilon_C*Z_C*gamma+Tau*Z_C
      
      
    }else if(Grp=="N"){
      lp.A <- X_A.mat %*% beta_A+epsilon_A*Z_A*gamma
      lp.N <- X_N.mat %*% beta_N+epsilon_x*Z_N*gamma
      lp.C <- X_C.mat %*% beta_C+epsilon_C*Z_C*gamma+Tau*Z_C   
      
    }else{
      lp.A <- X_A.mat %*% beta_A+epsilon_A*Z_A*gamma
      lp.N <- X_N.mat %*% beta_N+epsilon_N*Z_N*gamma
      lp.C <- X_C.mat %*% beta_C+epsilon_x*Z_C*gamma+Tau*Z_C
    }
    
    p.A<-1/(1+exp(-lp.A))
    p.N<-1/(1+exp(-lp.N))
    p.C<-1/(1+exp(-lp.C))
    
    if(Grp=="A"){
      
      
      f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A))+dnorm(epsilon_x,mean=mu_epsilon,sd=sigma_epsilon,log=T)+
        sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N))+dnorm(epsilon_N,mean=mu_epsilon,sd=sigma_epsilon,log=T)+
        sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C))+dnorm(epsilon_C,mean=mu_epsilon,sd=sigma_epsilon,log=T)
      
      
    }
    
    if(Grp=="N"){
      
      f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A))+dnorm(epsilon_A,mean=mu_epsilon,sd=sigma_epsilon,log=T)+
        sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N))+dnorm(epsilon_x,mean=mu_epsilon,sd=sigma_epsilon,log=T)+
        sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C))+dnorm(epsilon_C,mean=mu_epsilon,sd=sigma_epsilon,log=T)
      
    }
    
    if(Grp=="C"){
      
      f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A))+dnorm(epsilon_A,mean=mu_epsilon,sd=sigma_epsilon,log=T)+
        sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N))+dnorm(epsilon_N,mean=mu_epsilon,sd=sigma_epsilon,log=T)+
        sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C))+dnorm(epsilon_x,mean=mu_epsilon,sd=sigma_epsilon,log=T)
      
    }
  }
  
  
  if(sample.param=="mu_epsilon"){
    
    mu_epsilon_x<-param
    sigma_epsilon<-sqrt(sigma_epsilon)
    
    f=N_A*dnorm(epsilon_A,mean=mu_epsilon_x,sd=sigma_epsilon,log=T)+
      N_N*dnorm(epsilon_N,mean=mu_epsilon_x,sd=sigma_epsilon,log=T)+
      N_C*dnorm(epsilon_C,mean=mu_epsilon_x,sd=sigma_epsilon,log=T)+
      dnorm(mu_epsilon_x,mean=0,sd=prior.sd,log=T)
    
  }
  
  
  if(sample.param=="sigma_epsilon"){
    
    sigma_epsilon_x<-param
    
    f=invgamma::dinvgamma(sigma_epsilon_x,shape=prior.epsilon.shape,scale=prior.epsilon.scale,log=T)+dnorm(epsilon_A,mean=mu_epsilon,sd=sqrt(sigma_epsilon_x),log=T)+
      dnorm(epsilon_C,mean=mu_epsilon,sd=sqrt(sigma_epsilon_x),log=T)+dnorm(epsilon_N,mean=mu_epsilon,sd=sqrt(sigma_epsilon_x),log=T)
    
    #f=dhcauchy(sigma_epsilon_x,sigma=A,log=T)+dnorm(epsilon_A,mean=mu_epsilon,sd=sigma_epsilon_x,log=T)+
    #  dnorm(epsilon_C,mean=mu_epsilon,sd=sigma_epsilon_x,log=T)+dnorm(epsilon_N,mean=mu_epsilon,sd=sigma_epsilon_x,log=T)
    #dhcauchy(sigma_epsilon_x,sigma=A,log=T)
  }
  
  #alpha~N()
  if(sample.param=="gamma"){
    
    gamma_x<-param
    
    lp.A <- X_A.mat %*% beta_A+epsilon_A*Z_A*gamma_x
    lp.N <- X_N.mat %*% beta_N+epsilon_N*Z_N*gamma_x
    lp.C <- X_C.mat %*% beta_C+epsilon_C*Z_C*gamma_x+Tau*Z_C 
    
    p.A<-1/(1+exp(-lp.A))
    p.N<-1/(1+exp(-lp.N))
    p.C<-1/(1+exp(-lp.C))
    
    f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A))+invgamma::dinvgamma(gamma_x,shape=prior.gamma.shape,scale=prior.gamma.scale,log=T)+
      sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N))+
      sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C))
    
    
  }
  
  
  
  
  
  
  
  return(f)
  
}


MH_mcmc_expand_v2<-function(curr_X,args, accept,proposal_rt_df,proposal_rw_width,proposal_gamma_width,proposal_Hessian=NA){
  
  
  
  if(args$sample.param=="beta"){
    
    proposed_X=rmvnorm(1,mean=curr_X,sigma=proposal_Hessian)
    
    #lp_curr<-dmvnorm(c(curr_X),mean=curr_X,sigma=proposal_Hessian,log=T)
    #lp_prop<-dmvnorm(c(proposed_X),mean=curr_X,sigma=proposal_Hessian,log=T)
    
    proposed_X<-matrix(proposed_X,ncol=1)
    
    a=min(c(1,exp(target_param_expand_v2(proposed_X,args)-target_param_expand_v2(curr_X,args))))#+lp_curr-lp_prop))) #random walk 
    
    if(is.na(a)){
      
      x = curr_X        # otherwise "reject" move, and stay where we are
      accept=accept+0
      
      
      
    }else{
      
      if(runif(1)<a){
        x = proposed_X      # accept move with probabily min(1,A)
        accept=accept+1
        
        
      } else {
        x = curr_X      # otherwise "reject" move, and stay where we are
        accept=accept+0
        
      }   
      
    }
    
  }
  
  if(args$sample.param%in%c("epsilon","gamma","tau","mu_epsilon")){
    
    # proposed_X=rt.scaled(1,df=proposal_rt_df,mean=curr_X,sd=1) #student t 
    
    if(args$sample.param=="gamma"){
      proposed_X=rnorm(1,mean=curr_X,sd=proposal_gamma_width)
      
      
    }else{
      proposed_X=rnorm(1,mean=curr_X,sd=proposal_gamma_width)
    }
    lp_curr<-dnorm(curr_X,mean=curr_X,sd=proposal_gamma_width,log=T)
    lp_prop<-dnorm(proposed_X,mean=curr_X,sd=proposal_gamma_width,log=T)  
    
    a=min(c(1,exp(target_param_expand_v2(proposed_X,args)-target_param_expand_v2(curr_X,args)+lp_curr-lp_prop)))
    
    
    if(is.na(a)){
      
      x = curr_X        # otherwise "reject" move, and stay where we are
      accept=accept+0
      
      
      
      
    }else{
      
      if(runif(1)<a){
        x = proposed_X      # accept move with probabily min(1,A)
        accept=accept+1
        
        
      } else {
        x = curr_X      # otherwise "reject" move, and stay where we are
        accept=accept+0
        
      }   
      
    }
    
  }
  
  if(args$sample.param=="sigma_epsilon"){
    
    
    curr_X<-curr_X
    
    proposed_X=abs(rnorm(1,mean=curr_X,sd=proposal_rw_width))
    
    #a=exp(target_param_expand_v2(proposed_X,args)-target_param_expand_v2(curr_X,args))
    #   
    #   rchisq(1,df=1,ncp=curr_X)
    
    #lp_curr<-dnorm(curr_X,mean=curr_X,sd=proposal_rw_width,log=T)
    #lp_prop<-dnorm(proposed_X,mean=curr_X,sd=proposal_rw_width,log=T)
    
    
    #lp_curr<-dunif(curr_X,min=min(c(0,curr_X-proposal_rw_width)),max=curr_X+proposal_rw_width,log=T)
    #lp_prop<-dunif(proposed_X,min=min(c(0,curr_X-proposal_rw_width)),max=curr_X+proposal_rw_width,log=T)
    
    #lp_curr<-dchisq(curr_X,df=1,ncp=curr_X,log=T)#dnorm(curr_X,mean=curr_X,sd=proposal_rw_width,log=T)
    #lp_prop<-dchisq(proposed_X,df=1,ncp=curr_X,log=T)#dnorm(proposed_X,mean=curr_X,sd=proposal_rw_width,log=T)
    
    a=min(c(1,exp(target_param_expand_v2(proposed_X,args)-target_param_expand_v2(curr_X,args))))#+lp_curr-lp_prop)))#+dinvgamma(curr_X,shape=3,scale=1)-dinvgamma(proposed_X,shape=3,scale=1)))
    
    if(runif(1)<a){
      x = proposed_X      # accept move with probabily min(1,A)
      accept=accept+1
      
      
    } else {
      x = curr_X     # otherwise "reject" move, and stay where we are
      accept=accept+0
      
    }  
  }
  
  # if(args$sample.param=="gamma"){
  #   
  # 
  #   proposed_X= rchisq(1,df=1,ncp=curr_X)
  #   
  #   #lp_curr<-dunif(curr_X,min=min(c(0,curr_X-proposal_rw_width)),max=curr_X+proposal_rw_width,log=T)
  #   #lp_prop<-dunif(proposed_X,min=min(c(0,curr_X-proposal_rw_width)),max=curr_X+proposal_rw_width,log=T)
  #   
  #   lp_curr<-dchisq(curr_X,df=1,ncp=curr_X,log=T)#dnorm(curr_X,mean=curr_X,sd=proposal_rw_width,log=T)
  #   lp_prop<-dchisq(proposed_X,df=1,ncp=curr_X,log=T)#dnorm(proposed_X,mean=curr_X,sd=proposal_rw_width,log=T)
  #   
  #   a=min(na.omit(c(1,exp(target_param_expand(proposed_X,args)-target_param_expand(curr_X,args)+lp_curr-lp_prop))))#+dinvgamma(curr_X,shape=3,scale=1)-dinvgamma(proposed_X,shape=3,scale=1)))
  #   
  #   if(runif(1)<a){
  #     x = proposed_X      # accept move with probabily min(1,A)
  #     accept=accept+1
  #     
  #     
  #   } else {
  #     x = curr_X     # otherwise "reject" move, and stay where we are
  #     accept=accept+0
  #     
  #   }  
  # }
  # 
  #  }
  
  return(list(new.param=x,accept=accept))
  
}

### Only update parameters, G known ######
MCMC_binary_Cross_matching<-function(Iteration, Data=list(Y,D,Z,X,G),random.ini=F,
                                     prior_values=list(prior.sd,
                                                       A, proposal_rt_df,proposal_rw_width,
                                                       prior.gamma.shape,
                                                       prior.gamma.scale,
                                                       prior.epsilon.shape,
                                                       prior.epsilon.scale,
                                                       proposal_gamma_width),
                                     ini_values=list(beta_A_ini, beta_N_ini, beta_C_ini, 
                                                     gamma_ini,
                                                     alpha_A_ini, alpha_C_ini, alpha_N_ini,alpha_ini,
                                                     mu_alpha_ini,sigma_alpha_ini,Tau_ini),
                                     burn_in,thin,printYes=T){
  
  t=1
  
  if(T){  
    # Data 
    Y=Data$Y
    D=Data$D
    Z=Data$Z
    X=Data$X
    G=Data$G
    N=length(Y)
    
    
    Tau_ini<-ini_values$Tau_ini
    
    
    # place holders
    DID<-rep(0,Iteration)
    DID_C04<-rep(0,Iteration)
    DID_C09<-rep(0,Iteration)
    
    accept_sigma<-0
    accept_beta<-c(0,0,0)
    accept_alpha.A<-0
    accept_alpha.N<-0
    accept_alpha.C<-0      
    accept_alpha<-0
    accept_tau<-0
    accept_mu.alpha<-0
    accept_sigma.alpha<-0
    accept_gamma<-0
    
    beta_A<-array(0,dim=c(Iteration,ncol(X)+1))
    beta_C<-array(0,dim=c(Iteration,ncol(X)+1))
    beta_N<-array(0,dim=c(Iteration,ncol(X)+1))
    
    sigma_A<-rep(1,Iteration)
    sigma_C<-rep(1,Iteration)
    sigma_N<-rep(1,Iteration)
    
      alpha_A_ini=ini_values$alpha_A_ini
      alpha_C_ini=ini_values$alpha_C_ini
      alpha_N_ini=ini_values$alpha_N_ini
      
      gamma_ini=ini_values$gamma_ini
      
      epsilon_A_ini=alpha_A_ini/gamma_ini
      epsilon_N_ini=alpha_N_ini/gamma_ini
      epsilon_C_ini=alpha_C_ini/gamma_ini
      
      mu_alpha_ini<-ini_values$mu_alpha_ini
      sigma_alpha_ini<-ini_values$sigma_alpha_ini
      
      mu_epsilon_ini<-mu_alpha_ini/gamma_ini
      sigma_epsilon_ini<-sigma_alpha_ini/(gamma_ini^2)
      
      
      alpha_A<-rep(0,Iteration)
      alpha_C<-rep(0,Iteration)
      alpha_N<-rep(0,Iteration)
      mu_alpha<-rep(0,Iteration)
      sigma_alpha<-rep(1,Iteration)
      
      epsilon_A<-rep(0,Iteration)
      epsilon_C<-rep(0,Iteration)
      epsilon_N<-rep(0,Iteration)
      mu_epsilon<-rep(0,Iteration)
      sigma_epsilon<-rep(1,Iteration)
      
      gamma<-rep(1,Iteration)
      
    Tau<-rep(0,Iteration)
    
    # Priors 
    prior_A_scale=prior_values$A
    prior.sd=prior_values$prior.sd
    proposal_rw_width=prior_values$proposal_rw_width
    proposal_rt_df=prior_values$proposal_rt_df
    
    proposal_gamma_width=prior_values$proposal_gamma_width
    
    prior.gamma.shape=prior_values$prior.gamma.shape
    prior.gamma.scale=prior_values$prior.gamma.scale
    prior.epsilon.shape=prior_values$prior.epsilon.shape
    prior.epsilon.scale=prior_values$prior.epsilon.scale
    
    
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
    
    X_A.mat<-cbind(rep(1),X_A)
    X_C.mat<-cbind(rep(1),X_C)
    X_N.mat<-cbind(rep(1),X_N)
    
    X.mat<-cbind(1,X) 
    
    
    ## Get initial values for coefficients 
    {
      ## Get initial values for betas
      ini_dataA<-data.frame(cbind(Y_A,X_A))
      fit0.A<-glm(Y_A~.,data=ini_dataA,family = binomial)
      
      ini_dataN<-data.frame(cbind(Y_N,X_N))
      fit0.N<-glm(Y_N~.,data=ini_dataN,family = binomial)
      
      ini_dataC<-data.frame(cbind(Y_C,X_C))
      fit0.C<-glm(Y_C~.,data=ini_dataC,family = binomial)
      
      beta_A_ini<-matrix(coef(fit0.A),ncol=1)
      beta_N_ini<-matrix(coef(fit0.N),ncol=1)
      beta_C_ini<-matrix(coef(fit0.C),ncol=1)
      
      
    }
  }
  
  while(t <= Iteration){
    
    if(printYes==T){
      
      if(t%%1000==0){
        print(paste("Iteration=",t,sep=""))
      }
    }
    
    if(t==1){
      
      

      llikfun<-function(param,Grp,beta_A,beta_C,beta_N,epsilon_A,epsilon_C,epsilon_N,gamma.val,Tau){
        
        
        if(Grp=="A"){  
          
          beta_A<-matrix(param,ncol=1)
          beta_C<-matrix(c(beta_C),ncol=1)
          beta_N<-matrix(c(beta_N),ncol=1)
          
        }else if(Grp=="N"){
          
          beta_N<-matrix(param,ncol=1)
          beta_C<-matrix(c(beta_C),ncol=1)
          beta_A<-matrix(c(beta_A),ncol=1)
          
        }else{
          
          beta_C<-matrix(param,ncol=1)
          beta_A<-matrix(c(beta_A),ncol=1)
          beta_N<-matrix(c(beta_N),ncol=1)  
          
        }
        
        lp.A <- X_A.mat %*% beta_A+epsilon_A*Z_A*gamma.val
        lp.N <- X_N.mat %*% beta_N+epsilon_N*Z_N*gamma.val
        lp.C <- X_C.mat %*% beta_C+epsilon_C*Z_C*gamma.val+Tau*Z_C
        
        p.A<-1/(1+exp(-lp.A))
        p.N<-1/(1+exp(-lp.N))
        p.C<-1/(1+exp(-lp.C))
        
        
        f=(sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A))+dnorm(epsilon_A,mean= mu_epsilon_ini,sd=sqrt(sigma_epsilon_ini),log=T)+
             sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N))+dnorm(epsilon_N,mean= mu_epsilon_ini,sd=sqrt(sigma_epsilon_ini),log=T)+
             sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C))+dnorm(epsilon_C,mean= mu_epsilon_ini,sd=sqrt(sigma_epsilon_ini),log=T)
           + invgamma::dinvgamma(gamma_ini,shape=prior.gamma.shape,scale=prior.gamma.scale,log=T)
           +dnorm(mu_epsilon_ini,mean=0,sd=prior.sd,log=T)+invgamma::dinvgamma(sigma_epsilon,shape=prior.epsilon.shape,scale=prior.epsilon.scale,log=T))
        
        
        
        
        return(f)
      }
      

        if(random.ini==F){
          ml.A <- maxLik(llikfun,start=beta_A_ini,
                         Grp="A",beta_A=beta_A_ini,
                         beta_C=beta_C_ini,
                         beta_N=beta_N_ini,
                         epsilon_A=epsilon_A_ini,
                         epsilon_C=epsilon_C_ini,epsilon_N=epsilon_N_ini,gamma.val=gamma_ini,Tau=Tau_ini)
          
          Sigma_hes.A<-solve(maxLik::hessian(ml.A))
          Sigma_hes.A<-make.positive.definite(Sigma_hes.A)
          
        }else{
          
          Sigma_hes.A<-diag(0.00001,ncol(X.mat))
        }
        
        tmp.beta.A<-MH_mcmc_expand_v2(curr_X=beta_A_ini,
                                      args=list(sample.param="beta",Grp="A",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                prior.gamma.shape=prior.gamma.shape,
                                                prior.gamma.scale=prior.gamma.scale,
                                                prior.epsilon.shape=prior.epsilon.shape,
                                                prior.epsilon.scale=prior.epsilon.scale,
                                                A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                                beta_A=beta_A_ini,beta_C=beta_C_ini,beta_N=beta_N_ini, 
                                                epsilon_A=epsilon_A_ini,epsilon_C=epsilon_C_ini,epsilon_N=epsilon_N_ini,gamma=gamma_ini,
                                                mu_epsilon=mu_epsilon_ini,sigma_epsilon=sigma_epsilon_ini,Tau=Tau_ini),
                                      accept = accept_beta[1],proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,
                                      proposal_gamma_width = proposal_gamma_width,proposal_Hessian=Sigma_hes.A)
        
        beta_A[t,]<-tmp.beta.A$new.param
        accept_beta[1]<-tmp.beta.A$accept
        
        if(random.ini==F){
          
          ml.N <- maxLik(llikfun, start=beta_N_ini,
                         Grp="N",beta_A=beta_A[t,],
                         beta_C=beta_C_ini,
                         beta_N=beta_N_ini,
                         epsilon_A=epsilon_A_ini,
                         epsilon_C=epsilon_C_ini,epsilon_N=epsilon_N_ini,gamma.val=gamma_ini,Tau=Tau_ini)
          
          Sigma_hes.N<--solve(maxLik::hessian(ml.N))
          Sigma_hes.N<-make.positive.definite(Sigma_hes.N)
          
        }else{
          
          Sigma_hes.N<-diag(0.00001,ncol(X.mat))
        }
        
        tmp.beta.N<-MH_mcmc_expand_v2(curr_X=beta_N_ini,
                                      args=list(sample.param="beta",Grp="N",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                prior.gamma.shape=prior.gamma.shape,
                                                prior.gamma.scale=prior.gamma.scale,
                                                prior.epsilon.shape=prior.epsilon.shape,
                                                prior.epsilon.scale=prior.epsilon.scale,
                                                A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                                beta_A=beta_A[t,],beta_C=beta_C_ini,beta_N=beta_N_ini, 
                                                epsilon_A=epsilon_A_ini,epsilon_C=epsilon_C_ini,epsilon_N=epsilon_N_ini,gamma=gamma_ini,
                                                mu_epsilon=mu_epsilon_ini,sigma_epsilon=sigma_epsilon_ini,Tau=Tau_ini),
                                      accept = accept_beta[2],proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,
                                      proposal_gamma_width = proposal_gamma_width,proposal_Hessian=Sigma_hes.N)
        beta_N[t,]<-tmp.beta.N$new.param
        accept_beta[2]<-tmp.beta.N$accept
        
        if(random.ini==F){
          ml.C <- maxLik(llikfun, start=beta_C_ini,
                         Grp="C",beta_A=beta_A[t,],
                         beta_C=beta_C_ini,
                         beta_N=beta_N[t,],
                         epsilon_A=epsilon_A_ini,
                         epsilon_C=epsilon_C_ini,epsilon_N=epsilon_N_ini,gamma.val=gamma_ini,Tau=Tau_ini)
          
          Sigma_hes.C<--solve(maxLik::hessian(ml.C))
          Sigma_hes.C<-make.positive.definite(Sigma_hes.C)
          
        }else{
          
          Sigma_hes.C<-diag(0.00001,ncol(X.mat))
        }
        
        tmp.beta.C<-MH_mcmc_expand_v2(curr_X=beta_C_ini,
                                      args=list(sample.param="beta",Grp="C",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                prior.gamma.shape=prior.gamma.shape,
                                                prior.gamma.scale=prior.gamma.scale,
                                                prior.epsilon.shape=prior.epsilon.shape,
                                                prior.epsilon.scale=prior.epsilon.scale,
                                                A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                                beta_A=beta_A[t,],beta_C=beta_C_ini,beta_N=beta_N[t,], 
                                                epsilon_A=epsilon_A_ini,epsilon_C=epsilon_C_ini,epsilon_N=epsilon_N_ini,gamma=gamma_ini,
                                                mu_epsilon=mu_epsilon_ini,sigma_epsilon=sigma_epsilon_ini,Tau=Tau_ini),
                                      accept = accept_beta[3],proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,
                                      proposal_gamma_width = proposal_gamma_width,proposal_Hessian=Sigma_hes.C)
        
        beta_C[t,]<-tmp.beta.C$new.param
        accept_beta[3]<-tmp.beta.C$accept
        
        
        
        ## Update tau ##
        tmp.tau<-MH_mcmc_expand_v2(curr_X=Tau_ini,
                                   args=list(sample.param="tau",Grp=NA,temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                             prior.gamma.shape=prior.gamma.shape,
                                             prior.gamma.scale=prior.gamma.scale,
                                             prior.epsilon.shape=prior.epsilon.shape,
                                             prior.epsilon.scale=prior.epsilon.scale,
                                             A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                             beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                             epsilon_A=epsilon_A_ini,epsilon_C=epsilon_C_ini,epsilon_N=epsilon_N_ini,
                                             mu_epsilon=mu_epsilon_ini,sigma_epsilon=sigma_epsilon_ini,Tau=Tau_ini,gamma=gamma_ini  ),
                                   accept = accept_tau,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
        
        Tau[t]<-tmp.tau$new.param
        accept_tau<-tmp.tau$accept
        
        
        ### Update alpha using parameter expansion ###
        ## alpha_g=gamma*epsilon_g
        ## update epsilon~metropolis
        ## update sigma_epsilon~p(sigma_alpha|epsilon, y)
        ## upadte gamma~
        ## compute alpha_g=gamma*epsilon_g
        
        
        ## Update epsilon_g ##
        tmp.epsilon.A<-MH_mcmc_expand_v2(curr_X=epsilon_A_ini,
                                         args=list(sample.param="epsilon",Grp="A",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                   prior.gamma.shape=prior.gamma.shape,
                                                   prior.gamma.scale=prior.gamma.scale,
                                                   prior.epsilon.shape=prior.epsilon.shape,
                                                   prior.epsilon.scale=prior.epsilon.scale,
                                                   A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                                   beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                   epsilon_A=epsilon_A_ini,epsilon_C=epsilon_C_ini,epsilon_N=epsilon_N_ini,
                                                   mu_epsilon=mu_epsilon_ini,sigma_epsilon=sigma_epsilon_ini,Tau=Tau[t],gamma=gamma_ini),
                                         accept = accept_alpha.A,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
        
        epsilon_A[t]<-tmp.epsilon.A$new.param
        alpha_A[t]<-tmp.epsilon.A$new.param*gamma_ini
        accept_alpha.A<-tmp.epsilon.A$accept
        
        tmp.epsilon.N<-MH_mcmc_expand_v2(curr_X=epsilon_N_ini,
                                         args=list(sample.param="epsilon",Grp="N",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                   prior.gamma.shape=prior.gamma.shape,
                                                   prior.gamma.scale=prior.gamma.scale,
                                                   prior.epsilon.shape=prior.epsilon.shape,
                                                   prior.epsilon.scale=prior.epsilon.scale,
                                                   A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                                   beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                   epsilon_A=epsilon_A[t],epsilon_C=epsilon_C_ini,epsilon_N=epsilon_N_ini,gamma=gamma_ini,
                                                   mu_epsilon=mu_epsilon_ini,sigma_epsilon=sigma_epsilon_ini,Tau=Tau[t]),
                                         accept = accept_alpha.N,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
        
        epsilon_N[t]<-tmp.epsilon.N$new.param
        alpha_N[t]<-tmp.epsilon.N$new.param*gamma_ini
        accept_alpha.N<-tmp.epsilon.N$accept
        
        tmp.epsilon.C<-MH_mcmc_expand_v2(curr_X=epsilon_C_ini,
                                         args=list(sample.param="epsilon",Grp="C",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                   prior.gamma.shape=prior.gamma.shape,
                                                   prior.gamma.scale=prior.gamma.scale,
                                                   prior.epsilon.shape=prior.epsilon.shape,
                                                   prior.epsilon.scale=prior.epsilon.scale,
                                                   A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                                   beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                   epsilon_A=epsilon_A[t],epsilon_C=epsilon_C_ini,epsilon_N=epsilon_N[t],gamma=gamma_ini,
                                                   mu_epsilon=mu_epsilon_ini,sigma_epsilon=sigma_epsilon_ini,Tau=Tau[t]),
                                         accept = accept_alpha.C,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
        
        epsilon_C[t]<-tmp.epsilon.C$new.param
        alpha_C[t]<-tmp.epsilon.C$new.param*gamma_ini
        accept_alpha.C<-tmp.epsilon.C$accept
        
        ##update gamma ##
        tmp.gamma<-MH_mcmc_expand_v2(curr_X=gamma_ini,
                                     args=list(sample.param="gamma",Grp=NA,temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                               prior.gamma.shape=prior.gamma.shape,
                                               prior.gamma.scale=prior.gamma.scale,
                                               prior.epsilon.shape=prior.epsilon.shape,
                                               prior.epsilon.scale=prior.epsilon.scale,
                                               A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                               beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                               epsilon_A=epsilon_A[t],epsilon_C=epsilon_C[t],epsilon_N=epsilon_N[t],
                                               mu_epsilon=mu_epsilon_ini,sigma_epsilon=sigma_epsilon_ini,Tau=Tau[t],gamma=gamma_ini),
                                     accept = accept_gamma,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,
                                     proposal_gamma_width=proposal_gamma_width)
        
        gamma[t]<-tmp.gamma$new.param
        accept_gamma<-tmp.gamma$accept
        
        
        ## update mu_epsilon ##
        tmp.mu.epsilon<-MH_mcmc_expand_v2(curr_X=mu_epsilon_ini,
                                          args=list(sample.param="mu_epsilon",Grp=NA,temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                    prior.gamma.shape=prior.gamma.shape,
                                                    prior.gamma.scale=prior.gamma.scale,
                                                    prior.epsilon.shape=prior.epsilon.shape,
                                                    prior.epsilon.scale=prior.epsilon.scale,
                                                    A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                                    beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                    epsilon_A=epsilon_A[t],epsilon_C=epsilon_C[t],epsilon_N=epsilon_N[t],gamma=gamma[t],
                                                    mu_epsilon=mu_epsilon_ini,sigma_epsilon=sigma_epsilon_ini,Tau=Tau[t]),
                                          accept = accept_mu.alpha,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
        
        mu_epsilon[t]<-tmp.mu.epsilon$new.param
        mu_alpha[t]<-mu_epsilon[t]*gamma[t]
        accept_mu.alpha<-tmp.mu.epsilon$accept
        
        ## update sigma_epsilon ##
        tmp.sigma.epsilon<-MH_mcmc_expand_v2(curr_X=sigma_epsilon_ini,
                                             args=list(sample.param="sigma_epsilon",Grp=NA,temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                       prior.gamma.shape=prior.gamma.shape,
                                                       prior.gamma.scale=prior.gamma.scale,
                                                       prior.epsilon.shape=prior.epsilon.shape,
                                                       prior.epsilon.scale=prior.epsilon.scale,
                                                       A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                                       beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                       epsilon_A=epsilon_A[t],epsilon_C=epsilon_C[t],epsilon_N=epsilon_N[t],
                                                       gamma=gamma[t],
                                                       mu_epsilon=mu_epsilon[t],sigma_epsilon=sigma_epsilon_ini,Tau=Tau[t]),
                                             accept = accept_sigma.alpha,proposal_rw_width = proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
        
        sigma_epsilon[t]<-tmp.sigma.epsilon$new.param
        sigma_alpha[t]<-sigma_epsilon[t]*(gamma[t]^2)
        accept_sigma.alpha<-tmp.sigma.epsilon$accept
        

    }else{
     

        X_A.mat<-cbind(rep(1),X_A)
        X_N.mat<-cbind(rep(1),X_N)
        X_C.mat<-cbind(rep(1),X_C)
        

          tmp.beta.A<-MH_mcmc_expand_v2(curr_X=beta_A[t-1,],
                                        args=list(sample.param="beta",Grp="A",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                  prior.gamma.shape=prior.gamma.shape,
                                                  prior.gamma.scale=prior.gamma.scale,
                                                  prior.epsilon.shape=prior.epsilon.shape,
                                                  prior.epsilon.scale=prior.epsilon.scale,
                                                  A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                                  beta_A=beta_A[t-1,],beta_C=beta_C[t-1,],beta_N=beta_N[t-1,], 
                                                  epsilon_A=epsilon_A[t-1],epsilon_C=epsilon_C[t-1],epsilon_N=epsilon_N[t-1],gamma=gamma[t-1],
                                                  mu_epsilon=mu_epsilon[t-1],sigma_epsilon=sigma_epsilon[t-1],Tau=Tau[t-1]),
                                        accept = accept_beta[1],proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,
                                        proposal_gamma_width = proposal_gamma_width,proposal_Hessian=Sigma_hes.A)
          beta_A[t,]<-tmp.beta.A$new.param
          accept_beta[1]<-tmp.beta.A$accept
          
          
          tmp.beta.N<-MH_mcmc_expand_v2(curr_X=beta_N[t-1,],
                                        args=list(sample.param="beta",Grp="N",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                  prior.gamma.shape=prior.gamma.shape,
                                                  prior.gamma.scale=prior.gamma.scale,
                                                  prior.epsilon.shape=prior.epsilon.shape,
                                                  prior.epsilon.scale=prior.epsilon.scale,
                                                  A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                                  beta_A=beta_A[t,],beta_C=beta_C[t-1,],beta_N=beta_N[t-1,], 
                                                  epsilon_A=epsilon_A[t-1],epsilon_C=epsilon_C[t-1],epsilon_N=epsilon_N[t-1],gamma=gamma[t-1],
                                                  mu_epsilon=mu_epsilon[t-1],sigma_epsilon=sigma_epsilon[t-1],Tau=Tau[t-1]),
                                        accept = accept_beta[2],proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,
                                        proposal_gamma_width = proposal_gamma_width,proposal_Hessian=Sigma_hes.N)
          beta_N[t,]<-tmp.beta.N$new.param
          accept_beta[2]<-tmp.beta.N$accept
          
          tmp.beta.C<-MH_mcmc_expand_v2(curr_X=beta_C[t-1,],
                                        args=list(sample.param="beta",Grp="C",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                  prior.gamma.shape=prior.gamma.shape,
                                                  prior.gamma.scale=prior.gamma.scale,
                                                  prior.epsilon.shape=prior.epsilon.shape,
                                                  prior.epsilon.scale=prior.epsilon.scale,
                                                  A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                                  beta_A=beta_A[t,],beta_C=beta_C[t-1,],beta_N=beta_N[t,], 
                                                  epsilon_A=epsilon_A[t-1],epsilon_C=epsilon_C[t-1],epsilon_N=epsilon_N[t-1],gamma=gamma[t-1],
                                                  mu_epsilon=mu_epsilon[t-1],sigma_epsilon=sigma_epsilon[t-1],Tau=Tau[t-1]),
                                        accept = accept_beta[3],proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,
                                        proposal_gamma_width = proposal_gamma_width,proposal_Hessian=Sigma_hes.C)
          
          
          beta_C[t,]<-tmp.beta.C$new.param
          accept_beta[3]<-tmp.beta.C$accept  
          
          
          ## Update tau ##
          tmp.tau<-MH_mcmc_expand_v2(curr_X=Tau[t-1],
                                     args=list(sample.param="tau",Grp=NA,temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                               prior.gamma.shape=prior.gamma.shape,
                                               prior.gamma.scale=prior.gamma.scale,
                                               prior.epsilon.shape=prior.epsilon.shape,
                                               prior.epsilon.scale=prior.epsilon.scale,
                                               A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                               beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                               epsilon_A=epsilon_A[t-1],epsilon_C=epsilon_C[t-1],epsilon_N=epsilon_N[t-1],gamma=gamma[t-1],
                                               mu_epsilon=mu_epsilon[t-1],sigma_epsilon=sigma_epsilon[t-1],Tau=Tau[t-1]),
                                     accept = accept_tau,proposal_rt_df = proposal_rt_df, proposal_gamma_width = proposal_gamma_width,proposal_rw_width=proposal_rw_width)
          
          Tau[t]<-tmp.tau$new.param
          accept_tau<-tmp.tau$accept
          
          ## Update epsilon_g ##
          tmp.epsilon.A<-MH_mcmc_expand_v2(curr_X=epsilon_A[t-1],
                                           args=list(sample.param="epsilon",Grp="A",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                     prior.gamma.shape=prior.gamma.shape,
                                                     prior.gamma.scale=prior.gamma.scale,
                                                     prior.epsilon.shape=prior.epsilon.shape,
                                                     prior.epsilon.scale=prior.epsilon.scale,
                                                     A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                                     beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                     epsilon_A=epsilon_A[t-1],epsilon_C=epsilon_C[t-1],epsilon_N=epsilon_N[t-1],gamma=gamma[t-1],
                                                     mu_epsilon=mu_epsilon[t-1],sigma_epsilon=sigma_epsilon[t-1],Tau=Tau[t]),
                                           accept = accept_alpha.A,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
          
          epsilon_A[t]<-tmp.epsilon.A$new.param
          alpha_A[t]<-epsilon_A[t]*gamma[t-1]
          accept_alpha.A<-tmp.epsilon.A$accept
          
          tmp.epsilon.N<-MH_mcmc_expand_v2(curr_X=epsilon_N[t-1],
                                           args=list(sample.param="epsilon",Grp="N",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                     prior.gamma.shape=prior.gamma.shape,
                                                     prior.gamma.scale=prior.gamma.scale,
                                                     prior.epsilon.shape=prior.epsilon.shape,
                                                     prior.epsilon.scale=prior.epsilon.scale,
                                                     A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                                     beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                     epsilon_A=epsilon_A[t],epsilon_C=epsilon_C[t-1],epsilon_N=epsilon_N[t-1],gamma=gamma[t-1],
                                                     mu_epsilon=mu_epsilon[t-1],sigma_epsilon=sigma_epsilon[t-1],Tau=Tau[t]),
                                           accept = accept_alpha.N,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
          
          
          epsilon_N[t]<-tmp.epsilon.N$new.param
          alpha_N[t]<-epsilon_N[t]*gamma[t-1]
          accept_alpha.N<-tmp.epsilon.N$accept
          
          tmp.epsilon.C<-MH_mcmc_expand_v2(curr_X=epsilon_C[t-1],
                                           args=list(sample.param="epsilon",Grp="C",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                     prior.gamma.shape=prior.gamma.shape,
                                                     prior.gamma.scale=prior.gamma.scale,
                                                     prior.epsilon.shape=prior.epsilon.shape,
                                                     prior.epsilon.scale=prior.epsilon.scale,
                                                     A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                                     beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                     epsilon_A=epsilon_A[t],epsilon_C=epsilon_C[t-1],epsilon_N=epsilon_N[t],gamma=gamma[t-1],
                                                     mu_epsilon=mu_epsilon[t-1],sigma_epsilon=sigma_epsilon[t-1],Tau=Tau[t]),
                                           accept = accept_alpha.C,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
          
          epsilon_C[t]<-tmp.epsilon.C$new.param
          alpha_C[t]<-epsilon_C[t]*gamma[t-1]
          accept_alpha.C<-tmp.epsilon.C$accept
          
          ##update gamma ##
          tmp.gamma<-MH_mcmc_expand_v2(curr_X=gamma[t-1],
                                       args=list(sample.param="gamma",Grp=NA,temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                 prior.gamma.shape=prior.gamma.shape,
                                                 prior.gamma.scale=prior.gamma.scale,
                                                 prior.epsilon.shape=prior.epsilon.shape,
                                                 prior.epsilon.scale=prior.epsilon.scale,
                                                 A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                                 beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                 epsilon_A=epsilon_A[t],epsilon_C=epsilon_C[t],epsilon_N=epsilon_N[t],
                                                 mu_epsilon=mu_epsilon_ini,sigma_epsilon=sigma_epsilon_ini,Tau=Tau[t],gamma=gamma[t-1]),
                                       accept = accept_gamma,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
          
          gamma[t]<-tmp.gamma$new.param
          accept_gamma<-tmp.gamma$accept
          
          
          ## update mu_epsilon ##
          tmp.mu.epsilon<-MH_mcmc_expand_v2(curr_X=mu_epsilon[t-1],
                                            args=list(sample.param="mu_epsilon",Grp=NA,temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                      prior.gamma.shape=prior.gamma.shape,
                                                      prior.gamma.scale=prior.gamma.scale,
                                                      prior.epsilon.shape=prior.epsilon.shape,
                                                      prior.epsilon.scale=prior.epsilon.scale,
                                                      A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                                      beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                      epsilon_A=epsilon_A[t],epsilon_C=epsilon_C[t],epsilon_N=epsilon_N[t],gamma=gamma[t],
                                                      mu_epsilon=mu_epsilon[t-1],sigma_epsilon=sigma_epsilon[t-1],Tau=Tau[t]),
                                            accept = accept_mu.alpha,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
          
          
          mu_epsilon[t]<-tmp.mu.epsilon$new.param
          mu_alpha[t]<-mu_epsilon[t]*gamma[t]
          accept_mu.alpha<-tmp.mu.epsilon$accept
          
          ## update sigma_epsilon ##
          tmp.sigma.epsilon<-MH_mcmc_expand_v2(curr_X=sigma_epsilon[t-1],
                                               args=list(sample.param="sigma_epsilon",Grp=NA,temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                         prior.gamma.shape=prior.gamma.shape,
                                                         prior.gamma.scale=prior.gamma.scale,
                                                         prior.epsilon.shape=prior.epsilon.shape,
                                                         prior.epsilon.scale=prior.epsilon.scale,
                                                         A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,
                                                         beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                         epsilon_A=epsilon_A[t],epsilon_C=epsilon_C[t],epsilon_N=epsilon_N[t],gamma=gamma[t],
                                                         mu_epsilon=mu_epsilon[t],sigma_epsilon=sigma_epsilon[t-1],Tau=Tau[t]),
                                               accept = accept_sigma.alpha,proposal_rw_width = proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
          
          
          sigma_epsilon[t]<-tmp.sigma.epsilon$new.param
          sigma_alpha[t]<-sigma_epsilon[t]*gamma[t]^2
          accept_sigma.alpha<-tmp.sigma.epsilon$accept
          
      
    }
    
    
    ### imputation ###
    X_A.mat<-cbind(1,X_A)
    X_N.mat<-cbind(1,X_N)
    X_C.mat<-cbind(1,X_C)
    
    mu_A.pred<-X.mat%*%as.matrix(beta_A[t,],ncol=1)
    mu_N.pred<-X.mat%*%as.matrix(beta_N[t,],ncol=1)
    mu_C.pred<-X.mat%*%as.matrix(beta_C[t,],ncol=1)
    

      mu_C0_04<-X_C.mat%*%as.matrix(beta_C[t,],ncol=1)
      mu_C1_04<-X_C.mat%*%as.matrix(beta_C[t,],ncol=1)+Tau[t]
      
      mu_C0_09<-X_C.mat%*%as.matrix(beta_C[t,],ncol=1)+alpha_C[t]
      mu_C1_09<-X_C.mat%*%as.matrix(beta_C[t,],ncol=1)+alpha_C[t]+Tau[t]
      
      mu_pred_04=mu_A.pred*(G=="A")+ mu_C.pred*(G=="C")+ mu_N.pred*(G=="N")  
      mu_pred_09=mu_A.pred*(G=="A")+ (mu_C.pred+Tau[t])*(G=="C")+ mu_N.pred*(G=="N") 
      
    Y_imp_C0_04<-rbinom(N_C,1,invlogit(mu_C0_04))
    Y_imp_C1_04<-rbinom(N_C,1,invlogit(mu_C1_04))
    #
    Y_imp_C0_09<-rbinom(N_C,1,invlogit(mu_C0_09))
    Y_imp_C1_09<-rbinom(N_C,1,invlogit(mu_C1_09))
  
    
    Y_imp_04<-rbinom(N,1,invlogit(mu_pred_04))
    Y_imp_09<-rbinom(N,1,invlogit(mu_pred_09))
    
    #predicted odds ratio
      DID[t]=(mean(Y_imp_09[which(G=="C")])/(1-mean(Y_imp_09[which(G=="C")])))/(mean(Y_imp_04[which(G=="C")])/(1-mean(Y_imp_04[which(G=="C")])))
      DID_C04[t]<-(mean(Y_imp_C1_04)/(1-mean(Y_imp_C1_04)))/(mean(Y_imp_C0_04)/(1-mean(Y_imp_C0_04)))
      DID_C09[t]<-(mean(Y_imp_C1_09)/(1-mean(Y_imp_C1_09)))/(mean(Y_imp_C0_09)/(1-mean(Y_imp_C0_09)))
      
  
    
    t=t+1
    
    
  }
  
  
  ## results
    MCMC_list<-data.frame(alpha_A=alpha_A, alpha_N=alpha_N, alpha_C=alpha_C,
                          mu_alpha=mu_alpha, sigma_alpha=sigma_alpha,
                          Tau=Tau,
                          beta_A=beta_A,beta_C=beta_C,beta_N=beta_N)
    
    accept_rates=list(sigma_alpha=accept_sigma.alpha/Iteration,
                      mu_alpha=accept_mu.alpha/Iteration,
                      beta=accept_beta/Iteration,
                      gamma=accept_gamma/Iteration,
                      alpha_A=accept_alpha.A/Iteration,
                      alpha_N=accept_alpha.N/Iteration,
                      alpha_C=accept_alpha.C/Iteration,
                      tau=accept_tau/Iteration)
  
  return(list(param=MCMC_list,
              DID=DID,
              DID.04=DID_C04,
              DID.09=DID_C09,
              accept_rates=accept_rates))
  
  
  
  
}




Cross_binary_mathching<-function(Datalist=Datalist,Potential.outcome=F,Iteration,temporal.assumption="common", binary.model="logit",random.ini=T,
                                 prior_values=list(prior_A,prior_C,prior_N,prior.sd,
                                                   A, proposal_width),
                                 ini_values=list(beta_A_ini, beta_N_ini, beta_C_ini, 
                                                 alpha_A_ini, alpha_C_ini, alpha_N_ini,alpha_ini,
                                                 mu_alpha_ini,sigma_alpha_ini,Tau_ini),
                                 burn_in,thin,printYes=T){
  
  if(F){  
    ### Step 1: Match people from 2011 to people in 2017 use MA ##
    data_2011<-data[which(data$Z==0),]
    data_2017<-data[which(data$Z==1),]
    
    ## propensity score model using cohort in 2011 ##
    Model.str.2011<-formula(paste("D ~ ", paste0(var.choose, collapse = " + ")))
    PS2011_model<-glm(Model.str.2011,data=data_2011,family = binomial())
    
    ## Apply propensity score model on cohort 2017 ##
    Pred_2017<-predict(PS2011_model,newdata = data_2017,type="response")
    Pred_2011<-predict(PS2011_model,type="response")
    
    ## one-to-one match group 2011 with D=1 and 2017 with D=1 with replacement ##
    N01<-length(which(data_2011$D==1))
    data_2011_1<-data_2011[which(data_2011$D==1),]
    data_2017_1<-data_2017[which(data_2017$D==1),]
    
    Pred_2011_1<-Pred_2011[which(data_2011$D==1)]
    Pred_2017_1<-Pred_2017[which(data_2017$D==1)]
    
    data_2011_1$Matched_ID<-rep(0)
    
    for(i in 1:N01){
      
      dist_i<-abs(Pred_2011_1[i]-Pred_2017_1)  
      
      if(min(dist_i)<=0.01){
        matched_id<-which.min(dist_i) 
      }else{
        matched_id<-NA 
      }
      
      data_2011_1$Matched_ID[i]<-matched_id
      
    }
    
    
    ##remove non-matched ##
    data_2011_1<-na.omit(data_2011_1)
    data_2017_1_matched<-data_2017_1[data_2011_1$Matched_ID,]
    
    data_always_takers<-rbind(data_2011_1[,-ncol(data_2011_1)],data_2017_1_matched)
    data_always_takers$G_pred<-"Always_taker"  
    
    ##### Step 2: use 2017 propensity score model to match non-hospice decedent in 2017 ###
    ## propensity score model using cohort in 2017 ##
    Model.str.2017<-formula(paste("D ~ ", paste0(var.choose, collapse = " + ")))
    PS2017_model<-glm(Model.str.2017,data=data_2017,family = binomial())
    
    ## Apply propensity score model on cohort 2017 ##
    Pred_2011<-predict(PS2017_model,newdata = data_2011,type="response")
    Pred_2017<-predict(PS2017_model,type="response")
    
    ## one-to-one match group 2011 with D=1 and 2017 with D=1 with replacement ##
    N10<-length(which(data_2017$D==0))
    data_2011_0<-data_2011[which(data_2011$D==0),]
    data_2017_0<-data_2017[which(data_2017$D==0),]
    
    Pred_2017_0<-Pred_2017[which(data_2017$D==0)]
    Pred_2011_0<-Pred_2011[which(data_2011$D==0)]
    
    data_2017_0$Matched_ID<-rep(0)
    
    for(i in 1:N10){
      
      dist_i<-abs(Pred_2017_0[i]-Pred_2011_0)  
      
      if(min(dist_i)<=0.01){
        matched_id<-which.min(dist_i) 
      }else{
        matched_id<-NA 
      }
      
      data_2017_0$Matched_ID[i]<-matched_id
      
    }
    
    ##remove non-matched ##
    data_2017_0<-na.omit(data_2017_0)
    data_2011_0_matched<-data_2011_0[data_2017_0$Matched_ID,]
    
    data_never_takers<-rbind(data_2017_0[,-ncol(data_2017_0)],data_2011_0_matched)
    data_never_takers$G_pred<-"Never_taker"  
    
    
    ######Step 2: Left unmatched group 2017 who use hospice are set to compliers ###
    data_2017_unmatched<-data_2017_1[-data_2011_1$Matched_ID,]
    
    Pred_2017_unmatched<-predict(PS2017_model,newdata=data_2017_unmatched,type="response")
    
    N<-nrow(data_2017_unmatched)
    
    for(i in 1:N){
      
      dist_i<-abs(Pred_2017_unmatched[i]-Pred_2011_0)  
      
      if(min(dist_i)<=0.01){
        matched_id<-which.min(dist_i) 
      }else{
        matched_id<-NA 
      }
      
      data_2017_unmatched$Matched_ID[i]<-matched_id
      
      
    }
    
    ##remove non-matched ##
    data_2017_unmatched<-na.omit(data_2017_unmatched)
    data_2011_0_matched_c<-data_2011_0[data_2017_unmatched$Matched_ID,]
    
    data_compliers<-rbind(data_2017_unmatched[,-ncol(data_2011_1)],data_2011_0_matched_c)
    data_compliers$G_pred<-"Compliers" 
    
    data_model<-rbind(data_always_takers,data_never_takers,data_compliers)
  }
  
  
  if(Potential.outcome==T){
    
    # data_model$G_pred<-factor(data_model$G_pred,levels = c("Always_taker","Compliers","Never_taker"),labels=c("A","C","N"))
    # 
    # data_PO<-data_model
    # 
    # DX.mat<-model.matrix(~.-1,data=data_PO[,which(colnames(data_PO)%in%var.choose)])[,-3]
    # 
    # Datalist=list(Y=data_PO$Y, D=data_PO$D, Z=data_PO$Z,X=DX.mat,G=data_PO$G_pred)
    
    if(temporal.assumption=="common"){
      

      #### Method3.1: Potential outcome framework + common temporal assumption #####
      MCMC_model<-MCMC_binary_Cross_matching(Iteration=Iteration, Data=Datalist,random.ini=random.ini,
                                             temporal.assumption="common",
                                             prior_values=prior_values,
                                             ini_values=ini_values,
                                             burn_in=burn_in,thin=thin,printYes=T)
      
    }else if(temporal.assumption=="weak"){
      
      MCMC_model<-MCMC_binary_Cross_matching(Iteration=Iteration, Data=Datalist,random.ini=random.ini,
                                             temporal.assumption="weak",
                                             prior_values=prior_values,
                                             ini_values=ini_values,
                                             burn_in=burn_in,thin=thin,printYes=T)
      
    }
    
    
    Result<-MCMC_model
    

  }else{
    
    new.data<-rbind(data_compliers,data_never_takers)
    new.data$G2<-ifelse(new.data$G_pred=="Compliers",1,0)
    new.data$Post<-new.data$Z
    
    fit.DID<-glm(Y~G2+Post+X1+X2+G2*Post,data=new.data,family = binomial("logit"))
    DID<-exp(coef(fit.DID)[6]) #odds ratio
    CI<-exp(confint(fit.DID)[6,])
    
    Result<-list(CACE=DID,CACE_CI=CI,err=err)
    
  }
  
  return(Result)
  
}














