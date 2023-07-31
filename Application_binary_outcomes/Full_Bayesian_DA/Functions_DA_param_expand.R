

####### Algorithm for DA under weaker temporal assumption (parameter-expansion Gibbs) ###########

library(brms)
library(extraDistr)
library(invgamma)
library(MCMCglmm) 
library(MCMCpack)
library(mvtnorm)
library(maxLik)
library(BayesLogit)
library(rootSolve)
library(corpcor)



######### Functions file ################
# sd_alpha^2/sd_A^2=var.Ratio
invlogit<-function(x){
  
  exp(x)/(1+exp(x))
}

logit<-function(p){
  
  log(p/(1-p))
}

{
  Update.G<-function(Y,D,Z,G, w_A, w_N, w_C, mu_A, mu_C, mu_N, 
                     alpha_A, alpha_N, alpha_C, Tau){
    
    
    N<-length(Y)
    
    ## Draw compliers from (1,1) (mixture of always-takers and compliers)
    id_a_or_c<-which(Z==1&D==1)
    Y_a_or_c<-Y[id_a_or_c]
    
    Prob.muA1=invlogit(mu_A[id_a_or_c]+alpha_A)
    Prob.muC1=invlogit(mu_C[id_a_or_c]+alpha_C+Tau)
    
    Prob_YA1<-dbinom(Y_a_or_c,size=1,prob=Prob.muA1)
    Prob_YC1<-dbinom(Y_a_or_c,size=1,prob=Prob.muC1)
    
    Prob_C1<-w_C[id_a_or_c]*Prob_YC1/(w_A[id_a_or_c]*Prob_YA1+w_C[id_a_or_c]*Prob_YC1)
      
    G_a_or_c<-ifelse(rbinom(length(id_a_or_c),1,Prob_C1)==1,"C","A")
    
    ## Draw compliers from (0,0)(mixture of neve-takers and compliers)
    id_n_or_c<-which(Z==0&D==0)
    Y_n_or_c<-Y[id_n_or_c]
    Prob.muN0=invlogit(mu_N[id_n_or_c])
    Prob.muC0=invlogit(mu_C[id_n_or_c])
    Prob_YN0<-dbinom(Y_n_or_c,size=1,prob=Prob.muN0)
    Prob_YC0<-dbinom(Y_n_or_c,size=1,prob=Prob.muC0)  
    
    Prob_C0<-w_C[id_n_or_c]*Prob_YC0/(w_N[id_n_or_c]*Prob_YN0+w_C[id_n_or_c]*Prob_YC0)

    G_n_or_c<-ifelse(rbinom(length(id_n_or_c),1,Prob_C0),"C","N")
    
    G[id_a_or_c]<-G_a_or_c
    G[id_n_or_c]<-G_n_or_c
    
    return(G)
    
  }
  

  ## Metropolis-hastings algorithm ##
  #proposal density: scaled t_5 
  MH_mcmc<-function(curr_X,args, accept,proposal_rt_df,proposal_rw_width,gradient.mu,proposal_Hessian=NA){
    
    if(args$sample.param=="beta"){
      
      proposal_step<-proposal_Hessian
      proposed_X<-rnorm(1,mean=curr_X+proposal_step/2*gradient.mu,sd=sqrt(proposal_step))

      lp_curr<-dnorm(c(curr_X),mean=curr_X+proposal_step/2*gradient.mu,sd= sqrt(proposal_step),log=T)
      lp_prop<-dnorm(c(proposed_X),mean=curr_X+proposal_step/2*gradient.mu,sd= sqrt(proposal_step),log=T)
      
      a=min(c(1,exp(target_param(proposed_X,args)-target_param(curr_X,args)-lp_prop+lp_curr))) #random walk 
      
      if(is.na(a)){
        
        x = curr_X        # otherwise "reject" move, and stay where we are
        accept=accept+0
        
        
        
      }else{
        
        if(runif(1)<a){
          x = proposed_X      # accept move with probability min(1,A)
          accept=accept+1
          
          
        } else {
          x = curr_X      # otherwise "reject" move, and stay where we are
          accept=accept+0
          
        }   
        
      }
      
    }
    
    if(args$sample.param%in%c("alpha","tau","mu_alpha")){
      

      
      if(args$sample.param=="mu_alpha"){
        
        proposed_X=rnorm(1,mean=curr_X,sd=proposal_rw_width)
        
        lp_curr<-dnorm(c(curr_X),mean=curr_X,sd= proposal_step,log=T)
        lp_prop<-dnorm(c(proposed_X),mean=curr_X,sd= proposal_step,log=T)
        
      }else{
        
        proposed_X=rnorm(1,mean=curr_X+proposal_rw_width/2*gradient.mu,sd=sqrt(proposal_rw_width))
        
        lp_curr<-dnorm(c(curr_X),mean=curr_X+proposal_step/2*gradient.mu,sd= sqrt(proposal_step),log=T)
        lp_prop<-dnorm(c(proposed_X),mean=curr_X+proposal_step/2*gradient.mu,sd= sqrt(proposal_step),log=T)
        
      }
      
      
      a=min(c(1,exp(target_param(proposed_X,args)-target_param(curr_X,args)-lp_prop+lp_curr))) #random walk 
      
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
      
      
      a=min(c(1,exp(target_param(abs(proposed_X),args)-target_param(curr_X,args))))#+dinvgamma(curr_X,shape=3,scale=1)-dinvgamma(proposed_X,shape=3,scale=1)))
      
      if(runif(1)<a){
        x = proposed_X^2      # accept move with probabily min(1,A)
        accept=accept+1
        
        
      } else {
        x = curr_X^2        # otherwise "reject" move, and stay where we are
        accept=accept+0
        
      }  
    }
    
    #  }
    
    return(list(new.param=x,accept=accept))
    
  }
  
  # target function for logistic regression, beta, alpha, tau #  
  target_param<-function(param,args=list(sample.param="beta",Grp="A",k, temporal.assumption=temporal.assumption,prior.sd,
                                         A,Z,Y,X,G,w_A, w_C, w_N,beta_A,beta_C,beta_N,
                                         alpha_A,alpha_C,alpha_N,mu_alpha,sigma_alpha,alpha,Tau,Strata.model)){
    
    sample.param<-args$sample.param
    Grp=args$Grp
    temporal.assumption<-args$temporal.assumption
    Strata.model=args$Strata.model
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
    
    beta_A=args$beta_A
    beta_C=args$beta_C
    beta_N=args$beta_N
    
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
    
    w_A=args$w_A
    w_N=args$w_N
    w_C=args$w_C
    
 
      if(sample.param=="beta"){
        
        beta_x<-c(beta_A)
        beta_x[k]<-param
        beta_x<-matrix(beta_x,ncol=1)
        
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
        

          
          f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A)+log(w_A[G=="A"]))+dnorm(beta_x[k],mean=0,sd=prior.sd,log=T)
          sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N)+log(w_N[G=="N"]))+ sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C)+log(w_C[G=="C"]))
       
        
      }
      
      
      
      
      if(sample.param=="tau"){
        
        tau_x<-param
        
        lp.A <- X_A.mat %*% beta_A+alpha_A*Z_A
        lp.N <- X_N.mat %*% beta_N+alpha_N*Z_N
        lp.C <- X_C.mat %*% beta_C+alpha_C*Z_C+tau_x*Z_C
        
        p.A<-1/(1+exp(-lp.A))
        p.N<-1/(1+exp(-lp.N))
        p.C<-1/(1+exp(-lp.C))
        
        if(Strata.model==T){
          f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A)+log(w_A[G=="A"]))+
            sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N)+log(w_N[G=="N"]))+
            sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C)+log(w_C[G=="C"]))+dnorm(tau_x,mean=0,sd=prior.sd,log=T)
        }else{
          f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A)+log(w_A))+
            sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N)+log(w_N))+
            sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C)+log(w_C))+dnorm(tau_x,mean=0,sd=prior.sd,log=T)
        }
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
          
  
            
            f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A)+log(w_A[G=="A"]))+dnorm(alpha_x,mean=mu_alpha,sd=sqrt(sigma_alpha),log=T)+
              sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N)+log(w_N[G=="N"]))+dnorm(alpha_N,mean=mu_alpha,sd=sqrt(sigma_alpha),log=T)+
              sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C)+log(w_C[G=="C"]))+dnorm(alpha_C,mean=mu_alpha,sd=sqrt(sigma_alpha),log=T)
            
          
        }
        
        if(Grp=="N"){
          
      
            f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A)+log(w_A[G=="A"]))+dnorm(alpha_A,mean=mu_alpha,sd=sqrt(sigma_alpha),log=T)+
              sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N)+log(w_N[G=="N"]))+dnorm(alpha_x,mean=mu_alpha,sd=sqrt(sigma_alpha),log=T)+
              sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C)+log(w_C[G=="C"]))+dnorm(alpha_C,mean=mu_alpha,sd=sqrt(sigma_alpha),log=T)
            
          
        }
        
        if(Grp=="C"){
          
        
            f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A)+log(w_A[G=="A"]))+dnorm(alpha_A,mean=mu_alpha,sd=sqrt(sigma_alpha),log=T)+
              sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N)+log(w_N[G=="N"]))+dnorm(alpha_N,mean=mu_alpha,sd=sqrt(sigma_alpha),log=T)+
              sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C)+log(w_C[G=="C"]))+dnorm(alpha_x,mean=mu_alpha,sd=sqrt(sigma_alpha),log=T)
            
         
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
        
        sigma_alpha_x<-sqrt(param)
        
        f=dhcauchy(sigma_alpha_x,sigma=A,log=T)+dnorm(alpha_A,mean=mu_alpha,sd=sigma_alpha_x,log=T)+
          dnorm(alpha_C,mean=mu_alpha,sd=sigma_alpha_x,log=T)+dnorm(alpha_N,mean=mu_alpha,sd=sigma_alpha_x,log=T)
      }
      

    return(f)
    
  }
  

  target_param_expand<-function(param,args=list(sample.param="beta",Grp="A",
                                                prior.sd,
                                                prior.gamma.shape,
                                                prior.gamma.scale,
                                                prior.epsilon.shape,
                                                prior.epsilon.scale,
                                                A,Z,Y,X,G,w_A, w_C, w_N,beta_A,beta_C,beta_N, psi_A,psi_C,
                                                epsilon_A,epsilon_C,epsilon_N,mu_epsilon,sigma_epsilon,gamma,Tau)){
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
    Strata.model=args$Strata.model
    psi_A<-matrix(args$psi_A,ncol=1)
    psi_C<-matrix(args$psi_C,ncol=1)
    
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
    
    w_A=args$w_A
    w_N=args$w_N
    w_C=args$w_C
    
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
      
      

        f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A)+log(w_A[G=="A"]))+dmvnorm(c(beta_x),mean=rep(0,length(beta_x)),sigma=diag(prior.sd,length(beta_x)),log=T)+
          sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N)+log(w_N[G=="N"]))+
          sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C)+log(w_C[G=="C"]))
        
        
    }
    
    
    
    if(sample.param=="tau"){
      
      tau_x<-param
      
      lp.A <- X_A.mat %*% beta_A+epsilon_A*Z_A*gamma
      lp.N <- X_N.mat %*% beta_N+epsilon_N*Z_N*gamma
      lp.C <- X_C.mat %*% beta_C+epsilon_C*Z_C*gamma+tau_x*Z_C
      
      p.A<-1/(1+exp(-lp.A))
      p.N<-1/(1+exp(-lp.N))
      p.C<-1/(1+exp(-lp.C))

        f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A)+log(w_A[G=="A"]))+
          sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N)+log(w_N[G=="N"]))+
          sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C)+log(w_C[G=="C"]))+dnorm(tau_x,mean=0,sd=prior.sd,log=T)

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
        
 
          f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A)+log(w_A[G=="A"]))+dnorm(epsilon_x,mean=mu_epsilon,sd=sigma_epsilon,log=T)+
            sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N)+log(w_N[G=="N"]))+dnorm(epsilon_N,mean=mu_epsilon,sd=sigma_epsilon,log=T)+
            sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C)+log(w_C[G=="C"]))+dnorm(epsilon_C,mean=mu_epsilon,sd=sigma_epsilon,log=T)
          

      }
      
      if(Grp=="N"){
        
  
          f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A)+log(w_A[G=="A"]))+dnorm(epsilon_A,mean=mu_epsilon,sd=sigma_epsilon,log=T)+
            sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N)+log(w_N[G=="N"]))+dnorm(epsilon_x,mean=mu_epsilon,sd=sigma_epsilon,log=T)+
            sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C)+log(w_C[G=="C"]))+dnorm(epsilon_C,mean=mu_epsilon,sd=sigma_epsilon,log=T)
          

      }
      
      if(Grp=="C"){
        

          f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A)+log(w_A[G=="A"]))+dnorm(epsilon_A,mean=mu_epsilon,sd=sigma_epsilon,log=T)+
            sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N)+log(w_N[G=="N"]))+dnorm(epsilon_N,mean=mu_epsilon,sd=sigma_epsilon,log=T)+
            sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C)+log(w_C[G=="C"]))+dnorm(epsilon_x,mean=mu_epsilon,sd=sigma_epsilon,log=T)
          

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

        f=sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A)+log(w_A[G=="A"]))+invgamma::dinvgamma(gamma_x,shape=prior.gamma.shape,scale=prior.gamma.scale,log=T)+
          sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N)+log(w_N[G=="N"]))+
          sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C)+log(w_C[G=="C"]))

    }
    
    
    return(f)
    
  }
  
  
  MH_mcmc_expand<-function(curr_X,args, accept,proposal_rt_df,proposal_rw_width,proposal_gamma_width,proposal_Hessian=NA){

    
    if(args$sample.param=="beta"){
      
      proposed_X=rmvnorm(1,mean=curr_X,sigma=proposal_Hessian)
      
      proposed_X<-matrix(proposed_X,ncol=1)
      
      a=min(c(1,exp(target_param_expand(proposed_X,args)-target_param_expand(curr_X,args))))#+lp_curr-lp_prop))) #random walk 
      
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
      

      if(args$sample.param=="gamma"){
        proposed_X=rnorm(1,mean=curr_X,sd=proposal_gamma_width)
        
        
      }else{
        proposed_X=rnorm(1,mean=curr_X,sd=proposal_gamma_width)
      }
      lp_curr<-dnorm(curr_X,mean=curr_X,sd=proposal_gamma_width,log=T)
      lp_prop<-dnorm(proposed_X,mean=curr_X,sd=proposal_gamma_width,log=T)  
      
      a=min(c(1,exp(target_param_expand(proposed_X,args)-target_param_expand(curr_X,args)+lp_curr-lp_prop)))
      
      
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
  
      a=min(c(1,exp(target_param_expand(proposed_X,args)-target_param_expand(curr_X,args))))
      
      if(runif(1)<a){
        x = proposed_X      # accept move with probabily min(1,A)
        accept=accept+1
        
        
      } else {
        x = curr_X     # otherwise "reject" move, and stay where we are
        accept=accept+0
        
      }  
    }
    
    return(list(new.param=x,accept=accept))
    
  }
  
  ################# DA ######################################
  MCMC_Cross_temporal_binary<-function(Iteration, Data=list(Y,D,Z,X), random.ini=T,
                                       prior_values=list(prior_A,prior_C,prior_N,prior.sd,
                                                         A, proposal_rt_df,proposal_rw_width,
                                                         prior.gamma.shape,
                                                         prior.gamma.scale,
                                                         prior.epsilon.shape,
                                                         prior.epsilon.scale,
                                                         proposal_gamma_width),
                                       ini_values=list(w_A_ini,w_C_ini,w_N_ini,
                                                       beta_A_ini, beta_N_ini, beta_C_ini, 
                                                       gamma_ini,
                                                       alpha_A_ini, alpha_C_ini, alpha_N_ini,alpha_ini,
                                                       mu_alpha_ini,sigma_alpha_ini,Tau_ini),
                                       burn_in,thin,printYes=F)
  {
    
    t=1
    
    if(T){  
      
      y.mean09<-y.mean04<-rep(NA,Iteration)
      # Data 
      Y=Data$Y
      D=Data$D
      Z=Data$Z
      X=Data$X
      N=length(Y)
      
      # imputations
      # initial values
      w_A_ini=ini_values$w_A_ini
      w_N_ini=ini_values$w_N_ini
      w_C_ini=ini_values$w_C_ini
      
      Tau_ini<-ini_values$Tau_ini
      
      # place holders
      DID<-rep(0,Iteration)
      DID_C04<-rep(0,Iteration)
      DID_C09<-rep(0,Iteration)
      
      accept_beta<-c(0,0,0)
      accept_beta_A<-rep(0,ncol(X)+1)
      accept_beta_N<-rep(0,ncol(X)+1)
      accept_beta_C<-rep(0,ncol(X)+1)
      
      accept_alpha.A<-0
      accept_alpha.N<-0
      accept_alpha.C<-0      
      accept_alpha<-0
      accept_tau<-0
      accept_mu.alpha<-0
      accept_sigma.alpha<-0
      accept_gamma<-0

      G.ini<-array(NA,N)
      G.ini[which(D==1&Z==0)]<-"A"
      G.ini[which(D==0&Z==1)]<-"N"
      G.ini[which(D==1&Z==1)]<-sample(c("A","C"),length(which(D==1&Z==1)),replace = T)
      G.ini[which(D==0&Z==0)]<-sample(c("N","C"),length(which(D==0&Z==0)),replace = T)
      
      
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

      psi_A<-array(0,dim=c(Iteration,ncol(X)+1))
      psi_C<-array(0,dim=c(Iteration,ncol(X)+1))
      
      psi0<-rep(0,ncol(X)+1)   #prior mean for psi
      T0<-diag(.01,ncol(X)+1)  # Prior precision of psi
      
      ## calculate marginal G probability
      D.xmat<-as.matrix(cbind(rep(1),X))
      lpA<-D.xmat%*%psi_A[1,]
      lpC<-D.xmat%*%psi_C[1,]
      
      w_A_ini<-exp(lpA)/(1+exp(lpA)+exp(lpC))
      w_C_ini<-exp(lpC)/(1+exp(lpA)+exp(lpC))
      w_N_ini<-1/(1+exp(lpA)+exp(lpC))
      
      Tau<-rep(0,Iteration)
      
      # Priors 
      prior_A=prior_values$prior_A
      prior_C=prior_values$prior_C
      prior_N=prior_values$prior_N
      prior_A_scale=prior_values$A
      prior.sd=prior_values$prior.sd
      proposal_rw_width=prior_values$proposal_rw_width
      proposal_rt_df=prior_values$proposal_rt_df
      proposal_gamma_width=prior_values$proposal_gamma_width
      
      prior.gamma.shape=prior_values$prior.gamma.shape
      prior.gamma.scale=prior_values$prior.gamma.scale
      prior.epsilon.shape=prior_values$prior.epsilon.shape
      prior.epsilon.scale=prior_values$prior.epsilon.scale
      
      
    }
    
    ## Get initial values for coefficients 
    {
      
      id_A=which(G.ini=="A")
      id_C=which(G.ini=="C")
      id_N=which(G.ini=="N")
      
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
      
      ini_dataA<-data.frame(cbind(Y_A,X_A))
      fit0.A<-glm(Y_A~.,data=ini_dataA,family = binomial)
      
      ini_dataN<-data.frame(cbind(Y_N,X_N))
      fit0.N<-glm(Y_N~.,data=ini_dataN,family = binomial)
      
      ini_dataC<-data.frame(cbind(Y_C,X_C))
      fit0.C<-glm(Y_C~.,data=ini_dataC,family = binomial)
      
      beta_A_ini<-matrix(coef(fit0.A),ncol=1)
      beta_N_ini<-matrix(coef(fit0.N),ncol=1)
      beta_C_ini<-matrix(coef(fit0.C),ncol=1)
      
      X.mat<-cbind(1,X)
      
      mu_A_ini<-X.mat%*%beta_A_ini
      mu_N_ini<-X.mat%*%beta_N_ini
      mu_C_ini<-X.mat%*%beta_C_ini
      
    }
    
    while(t <= Iteration){
      
      if(printYes==T){
        
        if(t%%5000==0){
          print(paste("Iteration=",t,sep=""))
            temp<-data.frame(alpha_A=alpha_A[1:t],alpha_N=alpha_N[1:t],alpha_C=alpha_C[1:t],
                             Tau=Tau[1:t],DID_t2=DID_C04[1:t],DID_t2=DID_C09[1:t],
                             beta_A=beta_A[1:t,],beta_C=beta_C[1:t,],beta_N=beta_N[1:t,], 
                             psi_A=psi_A[1:t,],psi_C=psi_C[1:t,])

          
          save(temp,file=paste0("DA_weak_2t_temp_result_",t,".Rdata"))
          
        }
          
        }
      
      
      if(t==1){
        
        
        ## update Group label ###
     
        G<-Update.G(Y,D,Z,G.ini,w_A_ini, w_N_ini, w_C_ini, 
                      mu_A_ini, mu_C_ini, mu_N_ini, alpha_A_ini,alpha_N_ini, alpha_C_ini,Tau_ini)
          
        G<-c(G)
        
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
        X_N.mat<-cbind(rep(1),X_N)
        X_C.mat<-cbind(rep(1),X_C)
        
        X.mat<-cbind(rep(1),X)
        
   ## Gibbs sampler update psi~ multinomial logistic regression ##
            # Define category-specific binary responses (Note: cat N is reference)
            u.a<-1*(G=="A") 
            u.c<-1*(G=="C")
            
            # Update Always-takers
            c.a<-log(1+exp(X.mat%*%psi_C[1,])) 
            eta.a<-X.mat%*%psi_A[1,]-c.a
            w.a<-rpg(N,1,eta.a) 
            z.a<-(u.a-1/2)/w.a+c.a
            v<-solve(T0+crossprod(X.mat*sqrt(w.a)))
            m<-v%*%(T0%*%psi0+t(w.a*X.mat)%*%z.a) 
            psi_A[t,]<-c(rmvnorm(1,m,v))
            
            # Update Compliers
            c.c<-log(1+exp(X.mat%*%psi_A[1,])) 
            eta.c<-X.mat%*%psi_C[1,]-c.c 
            w.c<-rpg(N,1,eta.c) 
            z.c<-(u.c-1/2)/w.c+c.c 
            v<-solve(T0+crossprod(X.mat*sqrt(w.c)))
            m<-v%*%(T0%*%psi0+t(w.c*X.mat)%*%z.c) 
            psi_C[t,]<-c(rmvnorm(1,m,v))
            
            ## calculate marginal G probability
            D.xmat<-as.matrix(cbind(rep(1),X))
            lpA<-D.xmat%*%psi_A[t,]
            lpC<-D.xmat%*%psi_C[t,]
            
            w_A<-exp(lpA)/(1+exp(lpA)+exp(lpC))
            w_C<-exp(lpC)/(1+exp(lpA)+exp(lpC))
            w_N<-1/(1+exp(lpA)+exp(lpC))
 
        ## get starting values for beta 
        ini_dataA<-data.frame(cbind(Y_A,X_A,Z_A))
        fit0.A<-glm(Y_A~.,data=ini_dataA,family = binomial)
        
        ini_dataN<-data.frame(cbind(Y_N,X_N,Z_N))
        fit0.N<-glm(Y_N~.,data=ini_dataN,family = binomial)
        
        ini_dataC<-data.frame(cbind(Y_C,X_C,Z_C))
        fit0.C<-glm(Y_C~.,data=ini_dataC,family = binomial)
        
        
        beta_A_ini<-matrix(coef(fit0.A)[-length(coef(fit0.A))],ncol=1)
        beta_N_ini<-matrix(coef(fit0.N)[-length(coef(fit0.N))],ncol=1)
        beta_C_ini<-matrix(coef(fit0.C)[-length(coef(fit0.C))],ncol=1)
        

        {
  
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
            
          
              f=(sum(Y_A*log(p.A)+(1-Y_A)*log(1-p.A)+log(w_A[G=="A"]))+dnorm(epsilon_A,mean= mu_epsilon_ini,sd=sqrt(sigma_epsilon_ini),log=T)+
                   sum(Y_N*log(p.N)+(1-Y_N)*log(1-p.N)+log(w_N[G=="N"]))+dnorm(epsilon_N,mean= mu_epsilon_ini,sd=sqrt(sigma_epsilon_ini),log=T)+
                   sum(Y_C*log(p.C)+(1-Y_C)*log(1-p.C)+log(w_C[G=="C"]))+dnorm(epsilon_C,mean= mu_epsilon_ini,sd=sqrt(sigma_epsilon_ini),log=T)
                 + invgamma::dinvgamma(gamma_ini,shape=prior.gamma.shape,scale=prior.gamma.scale,log=T)
                 +dnorm(mu_epsilon_ini,mean=0,sd=prior.sd,log=T)+invgamma::dinvgamma(sigma_epsilon,shape=prior.epsilon.shape,scale=prior.epsilon.scale,log=T))
              
            
            
            return(f)
          }
        }
      
          ########## Parameter expansion ###########
          if(T){
            
            if(random.ini==F){
              ml.A <- maxLik(llikfun,start=beta_A_ini,
                             Grp="A",beta_A=beta_A_ini,
                             beta_C=beta_C_ini,
                             beta_N=beta_N_ini,
                             epsilon_A=epsilon_A_ini,
                             epsilon_C=epsilon_C_ini,epsilon_N=epsilon_N_ini,gamma.val=gamma_ini,Tau=Tau_ini)
              
              Sigma_hes.A<-solve(hessian(ml.A))
              Sigma_hes.A<-make.positive.definite(Sigma_hes.A)
              
            }else{
              
              Sigma_hes.A<-diag(0.000001,ncol(X.mat))
            }
            
            tmp.beta.A<-MH_mcmc_expand(curr_X=beta_A_ini,
                                       args=list(sample.param="beta",Grp="A",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                 prior.gamma.shape=prior.gamma.shape,
                                                 prior.gamma.scale=prior.gamma.scale,
                                                 prior.epsilon.shape=prior.epsilon.shape,
                                                 prior.epsilon.scale=prior.epsilon.scale,
                                                 A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                 beta_A=beta_A_ini,beta_C=beta_C_ini,beta_N=beta_N_ini, 
                                                 psi_A=psi_A[t,],psi_C=psi_C[t,],
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
              
              Sigma_hes.N<--solve(hessian(ml.N))
              Sigma_hes.N<-make.positive.definite(Sigma_hes.N)
              
            }else{
              
              Sigma_hes.N<-diag(0.000001,ncol(X.mat))
            }
            
            tmp.beta.N<-MH_mcmc_expand(curr_X=beta_N_ini,
                                       args=list(sample.param="beta",Grp="N",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                 prior.gamma.shape=prior.gamma.shape,
                                                 prior.gamma.scale=prior.gamma.scale,
                                                 prior.epsilon.shape=prior.epsilon.shape,
                                                 prior.epsilon.scale=prior.epsilon.scale,
                                                 A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                 beta_A=beta_A[t,],beta_C=beta_C_ini,beta_N=beta_N_ini, 
                                                 psi_A=psi_A[t,],psi_C=psi_C[t,],
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
              
              Sigma_hes.C<--solve(hessian(ml.C))
              Sigma_hes.C<-make.positive.definite(Sigma_hes.C)
              
            }else{
              
              Sigma_hes.C<-diag(0.000001,ncol(X.mat))
            }
            
            tmp.beta.C<-MH_mcmc_expand(curr_X=beta_C_ini,
                                       args=list(sample.param="beta",Grp="C",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                 prior.gamma.shape=prior.gamma.shape,
                                                 prior.gamma.scale=prior.gamma.scale,
                                                 prior.epsilon.shape=prior.epsilon.shape,
                                                 prior.epsilon.scale=prior.epsilon.scale,
                                                 A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                 beta_A=beta_A[t,],beta_C=beta_C_ini,beta_N=beta_N[t,], 
                                                 psi_A=psi_A[t,],psi_C=psi_C[t,],
                                                 epsilon_A=epsilon_A_ini,epsilon_C=epsilon_C_ini,epsilon_N=epsilon_N_ini,gamma=gamma_ini,
                                                 mu_epsilon=mu_epsilon_ini,sigma_epsilon=sigma_epsilon_ini,Tau=Tau_ini),
                                       accept = accept_beta[3],proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,
                                       proposal_gamma_width = proposal_gamma_width,proposal_Hessian=Sigma_hes.C)
            
            beta_C[t,]<-tmp.beta.C$new.param
            accept_beta[3]<-tmp.beta.C$accept
            
            
            
            ## Update tau ##
            tmp.tau<-MH_mcmc_expand(curr_X=Tau_ini,
                                    args=list(sample.param="tau",Grp=NA,temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                              prior.gamma.shape=prior.gamma.shape,
                                              prior.gamma.scale=prior.gamma.scale,
                                              prior.epsilon.shape=prior.epsilon.shape,
                                              prior.epsilon.scale=prior.epsilon.scale,
                                              A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                              beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                              psi_A=psi_A[t,],psi_C=psi_C[t,],
                                              epsilon_A=epsilon_A_ini,epsilon_C=epsilon_C_ini,epsilon_N=epsilon_N_ini,
                                              mu_epsilon=mu_epsilon_ini,sigma_epsilon=sigma_epsilon_ini,Tau=Tau_ini,gamma=gamma_ini,
                                              Strata.model=Strata.model),
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
            tmp.epsilon.A<-MH_mcmc_expand(curr_X=epsilon_A_ini,
                                          args=list(sample.param="epsilon",Grp="A",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                    prior.gamma.shape=prior.gamma.shape,
                                                    prior.gamma.scale=prior.gamma.scale,
                                                    prior.epsilon.shape=prior.epsilon.shape,
                                                    prior.epsilon.scale=prior.epsilon.scale,
                                                    A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                    beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                    psi_A=psi_A[t,],psi_C=psi_C[t,],
                                                    epsilon_A=epsilon_A_ini,epsilon_C=epsilon_C_ini,epsilon_N=epsilon_N_ini,
                                                    mu_epsilon=mu_epsilon_ini,sigma_epsilon=sigma_epsilon_ini,Tau=Tau[t],gamma=gamma_ini),
                                          accept = accept_alpha.A,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
            
            epsilon_A[t]<-tmp.epsilon.A$new.param
            alpha_A[t]<-tmp.epsilon.A$new.param*gamma_ini
            accept_alpha.A<-tmp.epsilon.A$accept
            
            tmp.epsilon.N<-MH_mcmc_expand(curr_X=epsilon_N_ini,
                                          args=list(sample.param="epsilon",Grp="N",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                    prior.gamma.shape=prior.gamma.shape,
                                                    prior.gamma.scale=prior.gamma.scale,
                                                    prior.epsilon.shape=prior.epsilon.shape,
                                                    prior.epsilon.scale=prior.epsilon.scale,
                                                    A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                    beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                    psi_A=psi_A[t,],psi_C=psi_C[t,],
                                                    epsilon_A=epsilon_A[t],epsilon_C=epsilon_C_ini,epsilon_N=epsilon_N_ini,gamma=gamma_ini,
                                                    mu_epsilon=mu_epsilon_ini,sigma_epsilon=sigma_epsilon_ini,Tau=Tau[t]),
                                          accept = accept_alpha.N,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
            
            epsilon_N[t]<-tmp.epsilon.N$new.param
            alpha_N[t]<-tmp.epsilon.N$new.param*gamma_ini
            accept_alpha.N<-tmp.epsilon.N$accept
            
            tmp.epsilon.C<-MH_mcmc_expand(curr_X=epsilon_C_ini,
                                          args=list(sample.param="epsilon",Grp="C",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                    prior.gamma.shape=prior.gamma.shape,
                                                    prior.gamma.scale=prior.gamma.scale,
                                                    prior.epsilon.shape=prior.epsilon.shape,
                                                    prior.epsilon.scale=prior.epsilon.scale,
                                                    A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                    beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                    psi_A=psi_A[t,],psi_C=psi_C[t,],
                                                    epsilon_A=epsilon_A[t],epsilon_C=epsilon_C_ini,epsilon_N=epsilon_N[t],gamma=gamma_ini,
                                                    mu_epsilon=mu_epsilon_ini,sigma_epsilon=sigma_epsilon_ini,Tau=Tau[t]),
                                          accept = accept_alpha.C,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
            
            epsilon_C[t]<-tmp.epsilon.C$new.param
            alpha_C[t]<-tmp.epsilon.C$new.param*gamma_ini
            accept_alpha.C<-tmp.epsilon.C$accept
            
            ##update gamma ##
            tmp.gamma<-MH_mcmc_expand(curr_X=gamma_ini,
                                      args=list(sample.param="gamma",Grp=NA,temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                prior.gamma.shape=prior.gamma.shape,
                                                prior.gamma.scale=prior.gamma.scale,
                                                prior.epsilon.shape=prior.epsilon.shape,
                                                prior.epsilon.scale=prior.epsilon.scale,
                                                A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                psi_A=psi_A[t,],psi_C=psi_C[t,],
                                                epsilon_A=epsilon_A[t],epsilon_C=epsilon_C[t],epsilon_N=epsilon_N[t],
                                                mu_epsilon=mu_epsilon_ini,sigma_epsilon=sigma_epsilon_ini,Tau=Tau[t],gamma=gamma_ini,
                                                Strata.model=Strata.model),
                                      accept = accept_gamma,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,
                                      proposal_gamma_width=proposal_gamma_width)
            
            gamma[t]<-tmp.gamma$new.param
            accept_gamma<-tmp.gamma$accept
            
            
            ## re-order alpha_g 
            #alpha_A[t]<-max(alpha_At,alpha_Ct,alpha_Nt)
            #alpha_C[t]<-c(alpha_At,alpha_Ct,alpha_Nt)[which(order(c(alpha_At,alpha_Ct,alpha_Nt))==2)]
            #alpha_N[t]<-min(alpha_At,alpha_Ct,alpha_Nt)
            
            ## update mu_epsilon ##
            tmp.mu.epsilon<-MH_mcmc_expand(curr_X=mu_epsilon_ini,
                                           args=list(sample.param="mu_epsilon",Grp=NA,temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                     prior.gamma.shape=prior.gamma.shape,
                                                     prior.gamma.scale=prior.gamma.scale,
                                                     prior.epsilon.shape=prior.epsilon.shape,
                                                     prior.epsilon.scale=prior.epsilon.scale,
                                                     A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                     beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                     psi_A=psi_A[t,],psi_C=psi_C[t,],
                                                     epsilon_A=epsilon_A[t],epsilon_C=epsilon_C[t],epsilon_N=epsilon_N[t],gamma=gamma[t],
                                                     mu_epsilon=mu_epsilon_ini,sigma_epsilon=sigma_epsilon_ini,Tau=Tau[t]),
                                           accept = accept_mu.alpha,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
            
            mu_epsilon[t]<-tmp.mu.epsilon$new.param
            mu_alpha[t]<-mu_epsilon[t]*gamma[t]
            accept_mu.alpha<-tmp.mu.epsilon$accept
            
            ## update sigma_epsilon ##
            tmp.sigma.epsilon<-MH_mcmc_expand(curr_X=sigma_epsilon_ini,
                                              args=list(sample.param="sigma_epsilon",Grp=NA,temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                        prior.gamma.shape=prior.gamma.shape,
                                                        prior.gamma.scale=prior.gamma.scale,
                                                        prior.epsilon.shape=prior.epsilon.shape,
                                                        prior.epsilon.scale=prior.epsilon.scale,
                                                        A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                        beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                        psi_A=psi_A[t,],psi_C=psi_C[t,],
                                                        epsilon_A=epsilon_A[t],epsilon_C=epsilon_C[t],epsilon_N=epsilon_N[t],gamma=gamma[t],
                                                        mu_epsilon=mu_epsilon[t],sigma_epsilon=sigma_epsilon_ini,Tau=Tau[t]),
                                              accept = accept_sigma.alpha,proposal_rw_width = proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
            
            sigma_epsilon[t]<-tmp.sigma.epsilon$new.param
            sigma_alpha[t]<-sigma_epsilon[t]*(gamma[t]^2)
            accept_sigma.alpha<-tmp.sigma.epsilon$accept
          }
          
          
        X_A.mat<-cbind(1,X_A) 
        X_N.mat<-cbind(1,X_N)
        X_C.mat<-cbind(1,X_C)
        
        mu_A=X.mat%*%beta_A[t,]
        mu_N=X.mat%*%beta_N[t,]
        mu_C=X.mat%*%beta_C[t,]

        }else{
        
         
         G<-Update.G(Y,D,Z,G,w_A, w_N, w_C,
                          mu_A, mu_C, mu_N, 
                          alpha_A[t-1],alpha_N[t-1], alpha_C[t-1],Tau[t-1])

          
          G<-c(G)
          
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
          X_N.mat<-cbind(rep(1),X_N)
          X_C.mat<-cbind(rep(1),X_C)
          
            ##Gibbs sampler  for updating stratification model ####3
            if(T){
              
              ## Gibbs sampler update psi~ multinomial logistic regression ##
              # Define category-specific binary responses (Note: cat N is reference)
              u.a<-1*(G=="A") 
              u.c<-1*(G=="C")
              
              psi_A.t<-psi_A[t-1,]
              psi_C.t<-psi_C[t-1,]

              # Update Always-takers
              c.a<-log(1+exp(X.mat%*%psi_C.t)) 
              eta.a<-X.mat%*%psi_A.t-c.a
              w.a<-rpg(N,1,eta.a) 
              z.a<-(u.a-1/2)/w.a+c.a
              v<-solve(T0+crossprod(X.mat*sqrt(w.a)))
              m<-v%*%(T0%*%psi0+t(w.a*X.mat)%*%z.a) 
              psi_A.t<-c(rmvnorm(1,m,v))
              
              # Update Compliers
              c.c<-log(1+exp(X.mat%*%psi_A.t)) 
              eta.c<-X.mat%*%psi_C.t-c.c 
              w.c<-rpg(N,1,eta.c) 
              z.c<-(u.c-1/2)/w.c+c.c 
              v<-solve(T0+crossprod(X.mat*sqrt(w.c)))
              m<-v%*%(T0%*%psi0+t(w.c*X.mat)%*%z.c) 
              psi_C.t<-c(rmvnorm(1,m,v))
              #}
              
              psi_A[t,]<-psi_A.t
              psi_C[t,]<-psi_C.t
              
              ## calculate marginal G probability
              D.xmat<-as.matrix(cbind(rep(1),X))
              lpA<-D.xmat%*%psi_A[t,]
              lpC<-D.xmat%*%psi_C[t,]
              
              w_A<-exp(lpA)/(1+exp(lpA)+exp(lpC))
              w_C<-exp(lpC)/(1+exp(lpA)+exp(lpC))
              w_N<-1/(1+exp(lpA)+exp(lpC))
              
              
            }
            
              
            ##### Update coefficients #########
              tmp.beta.A<-MH_mcmc_expand(curr_X=beta_A[t-1,],
                                         args=list(sample.param="beta",Grp="A",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                   prior.gamma.shape=prior.gamma.shape,
                                                   prior.gamma.scale=prior.gamma.scale,
                                                   prior.epsilon.shape=prior.epsilon.shape,
                                                   prior.epsilon.scale=prior.epsilon.scale,
                                                   A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                   beta_A=beta_A[t-1,],beta_C=beta_C[t-1,],beta_N=beta_N[t-1,], 
                                                   psi_A=psi_A[t,],psi_C=psi_C[t,],
                                                   epsilon_A=epsilon_A[t-1],epsilon_C=epsilon_C[t-1],epsilon_N=epsilon_N[t-1],gamma=gamma[t-1],
                                                   mu_epsilon=mu_epsilon[t-1],sigma_epsilon=sigma_epsilon[t-1],Tau=Tau[t-1]),
                                         accept = accept_beta[1],proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,
                                         proposal_gamma_width = proposal_gamma_width,proposal_Hessian=Sigma_hes.A)
              beta_A[t,]<-tmp.beta.A$new.param
              accept_beta[1]<-tmp.beta.A$accept
              
              
              tmp.beta.N<-MH_mcmc_expand(curr_X=beta_N[t-1,],
                                         args=list(sample.param="beta",Grp="N",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                   prior.gamma.shape=prior.gamma.shape,
                                                   prior.gamma.scale=prior.gamma.scale,
                                                   prior.epsilon.shape=prior.epsilon.shape,
                                                   prior.epsilon.scale=prior.epsilon.scale,
                                                   A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                   beta_A=beta_A[t,],beta_C=beta_C[t-1,],beta_N=beta_N[t-1,], 
                                                   psi_A=psi_A[t,],psi_C=psi_C[t,],
                                                   epsilon_A=epsilon_A[t-1],epsilon_C=epsilon_C[t-1],epsilon_N=epsilon_N[t-1],gamma=gamma[t-1],
                                                   mu_epsilon=mu_epsilon[t-1],sigma_epsilon=sigma_epsilon[t-1],Tau=Tau[t-1]),
                                         accept = accept_beta[2],proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,
                                         proposal_gamma_width = proposal_gamma_width,proposal_Hessian=Sigma_hes.N)
              beta_N[t,]<-tmp.beta.N$new.param
              accept_beta[2]<-tmp.beta.N$accept
              
              tmp.beta.C<-MH_mcmc_expand(curr_X=beta_C[t-1,],
                                         args=list(sample.param="beta",Grp="C",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                   prior.gamma.shape=prior.gamma.shape,
                                                   prior.gamma.scale=prior.gamma.scale,
                                                   prior.epsilon.shape=prior.epsilon.shape,
                                                   prior.epsilon.scale=prior.epsilon.scale,
                                                   A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                   beta_A=beta_A[t,],beta_C=beta_C[t-1,],beta_N=beta_N[t,], 
                                                   psi_A=psi_A[t,],psi_C=psi_C[t,],
                                                   epsilon_A=epsilon_A[t-1],epsilon_C=epsilon_C[t-1],epsilon_N=epsilon_N[t-1],gamma=gamma[t-1],
                                                   mu_epsilon=mu_epsilon[t-1],sigma_epsilon=sigma_epsilon[t-1],Tau=Tau[t-1]),
                                         accept = accept_beta[3],proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,
                                         proposal_gamma_width = proposal_gamma_width,proposal_Hessian=Sigma_hes.C)
              
              
              beta_C[t,]<-tmp.beta.C$new.param
              accept_beta[3]<-tmp.beta.C$accept  
              
              
              ## Update tau ##
              tmp.tau<-MH_mcmc_expand(curr_X=Tau[t-1],
                                      args=list(sample.param="tau",Grp=NA,temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                prior.gamma.shape=prior.gamma.shape,
                                                prior.gamma.scale=prior.gamma.scale,
                                                prior.epsilon.shape=prior.epsilon.shape,
                                                prior.epsilon.scale=prior.epsilon.scale,
                                                A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                psi_A=psi_A[t,],psi_C=psi_C[t,],
                                                epsilon_A=epsilon_A[t-1],epsilon_C=epsilon_C[t-1],epsilon_N=epsilon_N[t-1],gamma=gamma[t-1],
                                                mu_epsilon=mu_epsilon[t-1],sigma_epsilon=sigma_epsilon[t-1],Tau=Tau[t-1]),
                                      accept = accept_tau,proposal_rt_df = proposal_rt_df, proposal_gamma_width = proposal_gamma_width,proposal_rw_width=proposal_rw_width)
              
              Tau[t]<-tmp.tau$new.param
              accept_tau<-tmp.tau$accept
              
              ## Update epsilon_g ##
              tmp.epsilon.A<-MH_mcmc_expand(curr_X=epsilon_A[t-1],
                                            args=list(sample.param="epsilon",Grp="A",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                      prior.gamma.shape=prior.gamma.shape,
                                                      prior.gamma.scale=prior.gamma.scale,
                                                      prior.epsilon.shape=prior.epsilon.shape,
                                                      prior.epsilon.scale=prior.epsilon.scale,
                                                      A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                      beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                      psi_A=psi_A[t,],psi_C=psi_C[t,],
                                                      epsilon_A=epsilon_A[t-1],epsilon_C=epsilon_C[t-1],epsilon_N=epsilon_N[t-1],gamma=gamma[t-1],
                                                      mu_epsilon=mu_epsilon[t-1],sigma_epsilon=sigma_epsilon[t-1],Tau=Tau[t]),
                                            accept = accept_alpha.A,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
              
              
              epsilon_A[t]<-tmp.epsilon.A$new.param
              alpha_A[t]<-epsilon_A[t]*gamma[t-1]
              accept_alpha.A<-tmp.epsilon.A$accept
              
              tmp.epsilon.N<-MH_mcmc_expand(curr_X=epsilon_N[t-1],
                                            args=list(sample.param="epsilon",Grp="N",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                      prior.gamma.shape=prior.gamma.shape,
                                                      prior.gamma.scale=prior.gamma.scale,
                                                      prior.epsilon.shape=prior.epsilon.shape,
                                                      prior.epsilon.scale=prior.epsilon.scale,
                                                      A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                      beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                      psi_A=psi_A[t,],psi_C=psi_C[t,],
                                                      epsilon_A=epsilon_A[t],epsilon_C=epsilon_C[t-1],epsilon_N=epsilon_N[t-1],gamma=gamma[t-1],
                                                      mu_epsilon=mu_epsilon[t-1],sigma_epsilon=sigma_epsilon[t-1],Tau=Tau[t]),
                                            accept = accept_alpha.N,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
              
              
              epsilon_N[t]<-tmp.epsilon.N$new.param
              alpha_N[t]<-epsilon_N[t]*gamma[t-1]
              accept_alpha.N<-tmp.epsilon.N$accept
              
              tmp.epsilon.C<-MH_mcmc_expand(curr_X=epsilon_C[t-1],
                                            args=list(sample.param="epsilon",Grp="C",temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                      prior.gamma.shape=prior.gamma.shape,
                                                      prior.gamma.scale=prior.gamma.scale,
                                                      prior.epsilon.shape=prior.epsilon.shape,
                                                      prior.epsilon.scale=prior.epsilon.scale,
                                                      A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                      beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                      psi_A=psi_A[t,],psi_C=psi_C[t,],
                                                      epsilon_A=epsilon_A[t],epsilon_C=epsilon_C[t-1],epsilon_N=epsilon_N[t],gamma=gamma[t-1],
                                                      mu_epsilon=mu_epsilon[t-1],sigma_epsilon=sigma_epsilon[t-1],Tau=Tau[t]),
                                            accept = accept_alpha.C,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
              
              epsilon_C[t]<-tmp.epsilon.C$new.param
              alpha_C[t]<-epsilon_C[t]*gamma[t-1]
              accept_alpha.C<-tmp.epsilon.C$accept
              
              ##update gamma ##
              tmp.gamma<-MH_mcmc_expand(curr_X=gamma[t-1],
                                        args=list(sample.param="gamma",Grp=NA,temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                  prior.gamma.shape=prior.gamma.shape,
                                                  prior.gamma.scale=prior.gamma.scale,
                                                  prior.epsilon.shape=prior.epsilon.shape,
                                                  prior.epsilon.scale=prior.epsilon.scale,
                                                  A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                  beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                  psi_A=psi_A[t,],psi_C=psi_C[t,],
                                                  epsilon_A=epsilon_A[t],epsilon_C=epsilon_C[t],epsilon_N=epsilon_N[t],
                                                  mu_epsilon=mu_epsilon_ini,sigma_epsilon=sigma_epsilon_ini,Tau=Tau[t],gamma=gamma[t-1],
                                                  Strata.model=Strata.model),
                                        accept = accept_gamma,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
              
              gamma[t]<-tmp.gamma$new.param
              accept_gamma<-tmp.gamma$accept
              
              
              ## update mu_epsilon ##
              tmp.mu.epsilon<-MH_mcmc_expand(curr_X=mu_epsilon[t-1],
                                             args=list(sample.param="mu_epsilon",Grp=NA,temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                       prior.gamma.shape=prior.gamma.shape,
                                                       prior.gamma.scale=prior.gamma.scale,
                                                       prior.epsilon.shape=prior.epsilon.shape,
                                                       prior.epsilon.scale=prior.epsilon.scale,
                                                       A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                       beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                       psi_A=psi_A[t,],psi_C=psi_C[t,],
                                                       epsilon_A=epsilon_A[t],epsilon_C=epsilon_C[t],epsilon_N=epsilon_N[t],gamma=gamma[t],
                                                       mu_epsilon=mu_epsilon[t-1],sigma_epsilon=sigma_epsilon[t-1],Tau=Tau[t]),
                                             accept = accept_mu.alpha,proposal_rt_df = proposal_rt_df, proposal_rw_width=proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
              
              
              mu_epsilon[t]<-tmp.mu.epsilon$new.param
              mu_alpha[t]<-mu_epsilon[t]*gamma[t]
              accept_mu.alpha<-tmp.mu.epsilon$accept
              
              ## update sigma_epsilon ##
              tmp.sigma.epsilon<-MH_mcmc_expand(curr_X=sigma_epsilon[t-1],
                                                args=list(sample.param="sigma_epsilon",Grp=NA,temporal.assumption=temporal.assumption,prior.sd=prior.sd,
                                                          prior.gamma.shape=prior.gamma.shape,
                                                          prior.gamma.scale=prior.gamma.scale,
                                                          prior.epsilon.shape=prior.epsilon.shape,
                                                          prior.epsilon.scale=prior.epsilon.scale,
                                                          A=prior_A_scale,Z=Z,Y=Y,X=X,G=G,w_A=w_A,w_C= w_C, w_N=w_N,
                                                          beta_A=beta_A[t,],beta_C=beta_C[t,],beta_N=beta_N[t,], 
                                                          psi_A=psi_A[t,],psi_C=psi_C[t,],
                                                          epsilon_A=epsilon_A[t],epsilon_C=epsilon_C[t],epsilon_N=epsilon_N[t],gamma=gamma[t],
                                                          mu_epsilon=mu_epsilon[t],sigma_epsilon=sigma_epsilon[t-1],Tau=Tau[t]),
                                                accept = accept_sigma.alpha,proposal_rw_width = proposal_rw_width,proposal_gamma_width = proposal_gamma_width)
              
              
              sigma_epsilon[t]<-tmp.sigma.epsilon$new.param
              sigma_alpha[t]<-sigma_epsilon[t]*gamma[t]^2
              accept_sigma.alpha<-tmp.sigma.epsilon$accept
             
          
          X_A.mat<-cbind(1,X_A)
          X_N.mat<-cbind(1,X_N)
          X_C.mat<-cbind(1,X_C)
          
          mu_A=X.mat%*%beta_A[t,]
          mu_N=X.mat%*%beta_N[t,]
          mu_C=X.mat%*%beta_C[t,]
          
          
        
          
        }

      ### imputation ###
      X_C0.mat<-cbind(1,X[which(G=="C"),])
      X_C1.mat<-cbind(1,X[which(G=="C"),])
      
      mu_A.pred<-X.mat%*%as.matrix(beta_A[t,],ncol=1)
      mu_N.pred<-X.mat%*%as.matrix(beta_N[t,],ncol=1)
      mu_C.pred<-X.mat%*%as.matrix(beta_C[t,],ncol=1)
      
      mu_C0_04<-X_C.mat%*%as.matrix(beta_C[t,],ncol=1)
      mu_C1_04<-X_C.mat%*%as.matrix(beta_C[t,],ncol=1)+Tau[t]
      
      mu_C0_09<-X_C.mat%*%as.matrix(beta_C[t,],ncol=1)+alpha_C[t]
      mu_C1_09<-X_C.mat%*%as.matrix(beta_C[t,],ncol=1)+alpha_C[t]+Tau[t]
      # 
      mu_pred_04=mu_A.pred*(G=="A")+ mu_C.pred*(G=="C")+ mu_N.pred*(G=="N")  
      mu_pred_09=mu_A.pred*(G=="A")+ (mu_C.pred+Tau[t])*(G=="C")+ mu_N.pred*(G=="N") 
      

      Y_imp_C0_04<-rbinom(N_C,1,invlogit(mu_C0_04))
      Y_imp_C1_04<-rbinom(N_C,1,invlogit(mu_C1_04))
      
      Y_imp_C0_09<-rbinom(N_C,1,invlogit(mu_C0_09))
      Y_imp_C1_09<-rbinom(N_C,1,invlogit(mu_C1_09))

      Y_imp_04<-rbinom(N,1,invlogit(mu_pred_04))
      Y_imp_09<-rbinom(N,1,invlogit(mu_pred_09))
      
      #predicted odds ratio
      y.mean09[t]<-mean(Y_imp_09[which(G=="C")])
      y.mean04[t]<-mean(Y_imp_04[which(G=="C")])

        DID[t]=(y.mean09[t]/(1-y.mean09[t]))/(y.mean04[t]/(1-y.mean04[t]))#/exp(alpha_C[t]) #exp(tau)
        
        DID_C04[t]<-(mean(Y_imp_C1_04)/(1-mean(Y_imp_C1_04)))/(mean(Y_imp_C0_04)/(1-mean(Y_imp_C0_04)))
        DID_C09[t]<-(mean(Y_imp_C1_09)/(1-mean(Y_imp_C1_09)))/(mean(Y_imp_C0_09)/(1-mean(Y_imp_C0_09)))
        

      
      t=t+1
      
      
      
    }
    
   
    ## results
        MCMC_list<-data.frame(alpha_A=alpha_A, alpha_N=alpha_N, alpha_C=alpha_C,
                              mu_alpha=mu_alpha, sigma_alpha=sigma_alpha,
                              gamma=gamma,mu_epsilon=mu_epsilon,sigma_epsilon=sigma_epsilon,
                              Tau=Tau,DID=DID,DID_C04=DID_C04,DID_C09=DID_C09,
                              beta_A=beta_A,beta_C=beta_C,beta_N=beta_N, 
                              psi_A=psi_A,psi_C=psi_C)
        
        accept_rates=list(sigma_alpha=accept_sigma.alpha/Iteration,
                          mu_alpha=accept_mu.alpha/Iteration,
                          beta=accept_beta/Iteration,
                          gamma=accept_gamma/Iteration,
                          alpha_A=accept_alpha.A/Iteration,
                          alpha_N=accept_alpha.N/Iteration,
                          alpha_C=accept_alpha.C/Iteration,
                          tau=accept_tau/Iteration)
        

    return(list(param=MCMC_list, G=G,
                DID=DID,
                DID.04=DID_C04,
                DID.09=DID_C09,
                accept_rates=accept_rates))

  }
  
}  

