###########################################################################
###     Data Augmentation Algorithms for Cross-temporal Design          ###
###########################################################################

library(brms)
library(extraDistr)
library(invgamma)
library(MCMCglmm) 
library(MCMCpack)
library(mvtnorm)
library(maxLik)
library(BayesLogit)


### Data Augmentation with exclusion restrictions assumption ######
{
  ### Update Group Label: A, N, C ######
  Update.G.IV<-function(Y,D,Z,G, w_A, w_N, w_C, mu_C1, mu_C0,mu_A, mu_N, Outcome, sigma2=list(sigma2.A,sigma2.N,sigma2.C1,sigma2.C0)){
    
    N<-length(Y)
    
    ## Draw compliers from (1,1) (mixture of always-takers and compliers)
    id_a_or_c<-which(Z==1&D==1)
    Y_a_or_c<-Y[id_a_or_c]
    
    if(Outcome=="Binary"){
    Prob.muA1=invlogit(mu_A[id_a_or_c])
    Prob.muC1=invlogit(mu_C1[id_a_or_c])
    
    Prob_YA1<-dbinom(Y_a_or_c,size=1,prob=Prob.muA1)
    Prob_YC1<-dbinom(Y_a_or_c,size=1,prob=Prob.muC1)
    }
    
    if(Outcome=="Gaussian"){
      
      Prob_YA1<-dnorm(Y_a_or_c,mean=mu_A[id_a_or_c],sd=sqrt(sigma2$sigma2.A))
      Prob_YC1<-dnorm(Y_a_or_c,mean=mu_C1[id_a_or_c],sd=sqrt(sigma2$sigma2.C1)) 
      
    }
    
    
    Prob_C1<-w_C[id_a_or_c]*Prob_YC1/(w_A[id_a_or_c]*Prob_YA1+w_C[id_a_or_c]*Prob_YC1)
    
    G_a_or_c<-ifelse(rbinom(length(id_a_or_c),1,Prob_C1)==1,"C","A")
    
    ## Draw compliers from (0,0)(mixture of neve-takers and compliers)
    id_n_or_c<-which(Z==0&D==0)
    Y_n_or_c<-Y[id_n_or_c]
    
    if(Outcome=="Gaussian"){
      
      Prob_YN0<-dnorm(Y_n_or_c,mean=mu_N[id_n_or_c],sd=sqrt(sigma2$sigma2.N))
      Prob_YC0<-dnorm(Y_n_or_c,mean=mu_C0[id_n_or_c],sd=sqrt(sigma2$sigma2.C0)) 
      
    }
    
    if(Outcome=="Binary"){
    Prob.muN0=invlogit(mu_N[id_n_or_c])
    Prob.muC0=invlogit(mu_C0[id_n_or_c])
    Prob_YN0<-dbinom(Y_n_or_c,size=1,prob=Prob.muN0)
    Prob_YC0<-dbinom(Y_n_or_c,size=1,prob=Prob.muC0)  
    }
    
    Prob_C0<-w_C[id_n_or_c]*Prob_YC0/(w_N[id_n_or_c]*Prob_YN0+w_C[id_n_or_c]*Prob_YC0)
    
    G_n_or_c<-ifelse(rbinom(length(id_n_or_c),1,Prob_C0),"C","N")
    
    G[id_a_or_c]<-G_a_or_c
    G[id_n_or_c]<-G_n_or_c
    
    return(G)  
    
  }

  ### MCMC algorithm  ######
  # Datalist: list(Y=Y, D=D, Z=Z,X=X.matrix)
  # Y: outcome
  # D: 0/1 exposure status
  # Z: 0/1 baseline time=0, post baseline=1
  # X: matrix of covariates 
  # Outcome: "Gaussian" for continuous outcomes, "Binary" for binary outcomes
  # prior_values: only for continuous outcomes, sigma^2_g follows InvGamma(prior_a, prior_b) 
  DA_IV<-function(Iteration, Data,Outcome,prior_values=list(prior_a,prior_b))
  {
    # Data 
    Y=Data$Y
    D=Data$D
    Z=Data$Z
    X=Data$X
    N=length(Y)
    
    # place holders
    DID<-rep(0,Iteration)
    
    G<-array(NA,N)
    G[which(D==1&Z==0)]<-"A"
    G[which(D==0&Z==1)]<-"N"
    G[which(D==1&Z==1)]<-sample(c("A","C"),length(which(D==1&Z==1)),replace = T)
    G[which(D==0&Z==0)]<-sample(c("N","C"),length(which(D==0&Z==0)),replace = T)
    
    if(Outcome=="Binary"){
      
    ## calculate initial values
    temp.A<-as.data.frame(cbind(Y,X)[which(G=="A"),])
    fit.A<-glm(Y~.,data=temp.A,family = binomial)
    beta_A_ini<-matrix(coef(fit.A)[1:(1+ncol(X))] ,ncol=1) 
    
    
    temp.N<-as.data.frame(cbind(Y,X)[which(G=="N"),])
    fit.N<-glm(Y~.,data=temp.N,family = binomial)
    beta_N_ini<-matrix(coef(fit.N)[1:(1+ncol(X))],ncol=1)  
    
    
    temp.C1<-as.data.frame(cbind(Y,X)[which(G=="C"&Z==1),])
    fit.C1<-glm(Y~.,data=temp.C1,family = binomial)
    beta_C1_ini<-matrix(coef(fit.C1)[1:(1+ncol(X))] ,ncol=1) 
    
    temp.C0<-as.data.frame(cbind(Y,X)[which(G=="C"&Z==1),])
    fit.C0<-glm(Y~.,data=temp.C0,family = binomial)
    beta_C0_ini<-matrix(coef(fit.C0)[1:(1+ncol(X))] ,ncol=1) 
    
    }
    
    if(Outcome=="Gaussian"){
      
      prior_invgamma_a<-prior_values$prior_a
      prior_invgamma_b<-prior_values$prior_b
      
      ## calculate initial values
      temp.A<-as.data.frame(cbind(Y,X)[which(G=="A"),])
      fit.A<-lm(Y~.,data=temp.A)
      beta_A_ini<-matrix(coef(fit.A)[1:(1+ncol(X))] ,ncol=1) 
      
      
      temp.N<-as.data.frame(cbind(Y,X)[which(G=="N"),])
      fit.N<-lm(Y~.,data=temp.N)
      beta_N_ini<-matrix(coef(fit.N)[1:(1+ncol(X))],ncol=1)  
      
      
      temp.C1<-as.data.frame(cbind(Y,X)[which(G=="C"&Z==1),])
      fit.C1<-lm(Y~.,data=temp.C1)
      beta_C1_ini<-matrix(coef(fit.C1)[1:(1+ncol(X))] ,ncol=1) 
      
      temp.C0<-as.data.frame(cbind(Y,X)[which(G=="C"&Z==1),])
      fit.C0<-lm(Y~.,data=temp.C0)
      beta_C0_ini<-matrix(coef(fit.C0)[1:(1+ncol(X))] ,ncol=1) 
      
      sigma2.A<-rep(100,Iteration)
      sigma2.N<-rep(100,Iteration)
      sigma2.C1<-rep(100,Iteration)
      sigma2.C0<-rep(100,Iteration)
      
    }
    
    
    beta_A<-array(0,dim=c(Iteration,ncol(X)+1))
    beta_C1<-array(0,dim=c(Iteration,ncol(X)+1))
    beta_N<-array(0,dim=c(Iteration,ncol(X)+1))
    beta_C0<-array(0,dim=c(Iteration,ncol(X)+1))
    
    beta_A[1,]<-beta_A_ini
    beta_C1[1,]<-beta_C1_ini
    beta_N[1,]<-beta_N_ini
    beta_C0[1,]<-beta_C0_ini
    
    
    psi_A<-array(0,dim=c(Iteration,ncol(X)+1))
    psi_C<-array(0,dim=c(Iteration,ncol(X)+1))
    
    
    psi0<-rep(0,ncol(X)+1)   #prior mean for psi
    T0<-diag(.01,ncol(X)+1)  # Prior precision of psi
    
    ## calculate marginal G probability
    D.xmat<-as.matrix(cbind(rep(1),X))
    lpA<-D.xmat%*%psi_A[1,]
    lpC<-D.xmat%*%psi_C[1,]
    
    w_A<-exp(lpA)/(1+exp(lpA)+exp(lpC))
    w_C<-exp(lpC)/(1+exp(lpA)+exp(lpC))
    w_N<-1/(1+exp(lpA)+exp(lpC))
    
    
    t=2
    
    while(t<=Iteration){
      
      ### Impute G ###
      
      X.mat<-cbind(1,X)
      
      mu_A<-X.mat%*%beta_A[t-1,]
      mu_N<-X.mat%*%beta_N[t-1,]
      mu_C1<-X.mat%*%beta_C1[t-1,]
      mu_C0<-X.mat%*%beta_C0[t-1,]
      
      
      
      G<-Update.G.IV(Y,D,Z,G, w_A, w_N, w_C, mu_C1=mu_C1, mu_C0=mu_C0,  mu_A=mu_A, mu_N=mu_N,Outcome=Outcome,
                     sigma2=list(sigma2.A=sigma2.A[t-1],sigma2.N=sigma2.N[t-1],
                                 sigma2.C0=sigma2.C0[t-1],sigma2.C1=sigma2.C1[t-1]))
      
      
      ### Update psi ###
      ## Gibbs sampler update psi~ multinomial logistic regression ##
      # Define category-specific binary responses (Note: cat N is reference)
      u.a<-1*(G=="A") 
      u.c<-1*(G=="C")
      
      # Update Always-takers
      c.a<-log(1+exp(X.mat%*%psi_C[t-1,]))
      eta.a<-X.mat%*%psi_A[t-1,]-c.a
      w.a<-rpg(N,1,eta.a) 
      z.a<-(u.a-1/2)/w.a+c(c.a)
      v<-solve(T0+crossprod(X.mat*sqrt(w.a)))
      m<-v%*%(T0%*%psi0+t(w.a*X.mat)%*%z.a) 
      psi_A[t,]<-c(rmvnorm(1,m,v))
      
      # Update Compliers
      c.c<-log(1+exp(X.mat%*%psi_A[t,])) 
      eta.c<-X.mat%*%psi_C[t-1,]-c.c 
      w.c<-rpg(N,1,eta.c) 
      z.c<-(u.c-1/2)/w.c+c(c.c) 
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
      
      
      ### Update parameters in the outcomes models #####
      if(Outcome=="Gaussian"){
        
        
        ## Big design matrix ##
        q<-ncol(X)
        Bigmat.A<-X.mat[which(G=="A"),]
        colnames(Bigmat.A)<-c(paste("beta",0:q,"A",sep="_"))
        p<-ncol(Bigmat.A)
        mu_beta<-matrix(rep(0,p),ncol=1) #prior
        T0_beta<-diag(10000,p)
        
        ## update betas ##
        inv.Sigma0<-solve(T0_beta)
        M<-solve(inv.Sigma0+1/sigma2.A[t-1]*t(Bigmat.A)%*%Bigmat.A)%*%(inv.Sigma0%*%mu_beta+1/sigma2.A[t-1]*t(Bigmat.A)%*%Y[which(G=="A")])
        V<-solve(inv.Sigma0+1/sigma2.A[t-1]*t(Bigmat.A)%*%Bigmat.A)
        beta_A[t,]<-c(rmvnorm(1,M,V))
        
        ## update sigma2 ##
        post_a<-sum(G=="A")+prior_invgamma_a
        post_b<-prior_invgamma_b+crossprod(Y[which(G=="A")]-Bigmat.A%*%beta_A[t,])
        sigma2.A[t]<-invgamma::rinvgamma(1,shape=post_a/2,rate=post_b/2)
        
        
        q<-ncol(X)
        Bigmat.N<-X.mat[which(G=="N"),]
        colnames(Bigmat.N)<-c(paste("beta",0:q,"N",sep="_"))
        p<-ncol(Bigmat.N)
        mu_beta<-matrix(rep(0,p),ncol=1) #prior
        T0_beta<-diag(10000,p)
        
        ## update betas ##
        inv.Sigma0<-solve(T0_beta)
        M<-solve(inv.Sigma0+1/sigma2.N[t-1]*t(Bigmat.N)%*%Bigmat.N)%*%(inv.Sigma0%*%mu_beta+1/sigma2.N[t-1]*t(Bigmat.N)%*%Y[which(G=="N")])
        V<-solve(inv.Sigma0+1/sigma2.N[t-1]*t(Bigmat.N)%*%Bigmat.N)
        beta_N[t,]<-c(rmvnorm(1,M,V))
        
        ## update sigma2 ##
        post_a<-sum(G=="N")+prior_invgamma_a
        post_b<-prior_invgamma_b+crossprod(Y[which(G=="N")]-Bigmat.N%*%beta_N[t,])
        sigma2.N[t]<-invgamma::rinvgamma(1,shape=post_a/2,rate=post_b/2)

        q<-ncol(X)
        Bigmat.C1<-X.mat[which(G=="C"&Z==1),]
        colnames(Bigmat.C1)<-c(paste("beta",0:q,"C1",sep="_"))
        p<-ncol(Bigmat.C1)
        mu_beta<-matrix(rep(0,p),ncol=1) #prior
        T0_beta<-diag(10000,p)
        
        ## update betas ##
        inv.Sigma0<-solve(T0_beta)
        M<-solve(inv.Sigma0+1/sigma2.C1[t-1]*t(Bigmat.C1)%*%Bigmat.C1)%*%(inv.Sigma0%*%mu_beta+1/sigma2.C1[t-1]*t(Bigmat.C1)%*%Y[which(G=="C"&Z==1)])
        V<-solve(inv.Sigma0+1/sigma2.C1[t-1]*t(Bigmat.C1)%*%Bigmat.C1)
        beta_C1[t,]<-c(rmvnorm(1,M,V))
        
        ## update sigma2 ##
        post_a<-sum(G=="C"&Z==1)+prior_invgamma_a
        post_b<-prior_invgamma_b+crossprod(Y[which(G=="C"&Z==1)]-Bigmat.C1%*%beta_C1[t,])
        sigma2.C1[t]<-invgamma::rinvgamma(1,shape=post_a/2,rate=post_b/2)
        
        
        q<-ncol(X)
        Bigmat.C0<-X.mat[which(G=="C"&Z==0),]
        colnames(Bigmat.C0)<-c(paste("beta",0:q,"C0",sep="_"))
        p<-ncol(Bigmat.C0)
        mu_beta<-matrix(rep(0,p),ncol=1) #prior
        T0_beta<-diag(10000,p)
        
        ## update betas ##
        inv.Sigma0<-solve(T0_beta)
        M<-solve(inv.Sigma0+1/sigma2.C0[t-1]*t(Bigmat.C0)%*%Bigmat.C0)%*%(inv.Sigma0%*%mu_beta+1/sigma2.C0[t-1]*t(Bigmat.C0)%*%Y[which(G=="C"&Z==0)])
        V<-solve(inv.Sigma0+1/sigma2.C0[t-1]*t(Bigmat.C0)%*%Bigmat.C0)
        beta_C0[t,]<-c(rmvnorm(1,M,V))
        
        ## update sigma2 ##
        post_a<-sum(G=="C"&Z==0)+prior_invgamma_a
        post_b<-prior_invgamma_b+crossprod(Y[which(G=="C"&Z==0)]-Bigmat.C0%*%beta_C0[t,])
        sigma2.C0[t]<-invgamma::rinvgamma(1,shape=post_a/2,rate=post_b/2)
        

      }
      
      if(Outcome=="Binary"){
        
      ## Big design matrix ##
      q<-ncol(X)
      Bigmat.A<-X.mat[which(G=="A"),]
      colnames(Bigmat.A)<-c(paste("beta",0:q,"A",sep="_"))
      p<-ncol(Bigmat.A)
      mu_beta<-rep(0,p)
      T0_beta<-diag(0.01,p)
      
      eta.A<-Bigmat.A%*%beta_A[t-1,]
      w<-rpg(sum(G=="A"),1,eta.A)
      z<-(Y[which(G=="A")]-1/2)/w
      v<-solve(crossprod(Bigmat.A*sqrt(w))+T0_beta) 
      m<-v%*%(T0_beta%*%mu_beta+t(w*Bigmat.A)%*%z) 
      beta_A[t,]<-c(rmvnorm(1,m,v))
      
      
      Bigmat.N<-X.mat[which(G=="N"),]
      colnames(Bigmat.N)<-c(paste("beta",0:q,"N",sep="_"))
      p<-ncol(Bigmat.N)
      mu_beta<-rep(0,p)
      T0_beta<-diag(0.01,p)
      
      eta.N<-Bigmat.N%*%beta_N[t-1,]
      w<-rpg(sum(G=="N"),1,eta.N)
      z<-(Y[which(G=="N")]-1/2)/w
      v<-solve(crossprod(Bigmat.N*sqrt(w))+T0_beta) 
      m<-v%*%(T0_beta%*%mu_beta+t(w*Bigmat.N)%*%z) 
      beta_N[t,]<-c(rmvnorm(1,m,v))
      
      
      
      Bigmat.C1<-X.mat[which(G=="C"&Z==1),]
      colnames(Bigmat.C1)<-c(paste("beta",0:q,"C1",sep="_"))
      p<-ncol(Bigmat.C1)
      mu_beta<-rep(0,p)
      T0_beta<-diag(0.01,p)
      
      eta.C1<-Bigmat.C1%*%beta_C1[t-1,]
      w<-rpg(sum(G=="C"&Z==1),1,eta.C1)
      z<-(Y[which(G=="C"&Z==1)]-1/2)/w
      v<-solve(crossprod(Bigmat.C1*sqrt(w))+T0_beta) 
      m<-v%*%(T0_beta%*%mu_beta+t(w*Bigmat.C1)%*%z) 
      beta_C1[t,]<-c(rmvnorm(1,m,v))
      
      
      Bigmat.C0<-X.mat[which(G=="C"&Z==1),]
      colnames(Bigmat.C0)<-c(paste("beta",0:q,"C0",sep="_"))
      p<-ncol(Bigmat.C0)
      mu_beta<-rep(0,p)
      T0_beta<-diag(0.01,p)
      
      eta.C0<-Bigmat.C0%*%beta_C0[t-1,]
      w<-rpg(sum(G=="C"&Z==1),1,eta.C0)
      z<-(Y[which(G=="C"&Z==1)]-1/2)/w
      v<-solve(crossprod(Bigmat.C0*sqrt(w))+T0_beta) 
      m<-v%*%(T0_beta%*%mu_beta+t(w*Bigmat.C0)%*%z) 
      beta_C0[t,]<-c(rmvnorm(1,m,v))
      
      }
      
      mu_A=X.mat%*%beta_A[t,]
      mu_N=X.mat%*%beta_N[t,]
      mu_C0=X.mat%*%beta_C0[t,]
      mu_C1=X.mat%*%beta_C1[t,]
      
      
      ### imputation of PO ###
      # N_C<-sum(Gt=="C")
      
      mu_C0.pred<-mu_C0[which(G=="C")]
      mu_C1.pred<-mu_C1[which(G=="C")]
      
      if(Outcome=="Gaussian"){
        Y_imp_C0<-rnorm(length(mu_C0.pred),mean=mu_C0.pred,sd=sqrt(sigma2.C0[t]))
        Y_imp_C1<-rnorm(length(mu_C1.pred),mean=mu_C1.pred,sd=sqrt(sigma2.C1[t]))
        
  
        DID[t]=mean(Y_imp_C1-Y_imp_C0)
        
        
        MCMC_list<-data.frame(beta_A=beta_A,beta_C1=beta_C1,beta_C0=beta_C0,beta_N=beta_N, 
                              sigma2.A=sigma2.A, sigma2.N=sigma2.N, sigma2.C1=sigma2.C1, sigma2.C0=sigma2.C0,
                              psi_A=psi_A,psi_C=psi_C)
        
      }
      
      if(Outcome=="Binary"){
      Y_imp_C0<-rbinom(length(mu_C0.pred),1,invlogit(mu_C0.pred))
      Y_imp_C1<-rbinom(length(mu_C1.pred),1,invlogit(mu_C1.pred))
      
      
      #predicted odds ratio
      y.C0<-mean(Y_imp_C0)
      y.C1<-mean(Y_imp_C1)
      
      
      DID[t]=((y.C1/(1-y.C1))/(y.C0/(1-y.C0)))
      
      
      MCMC_list<-data.frame(beta_A=beta_A,beta_C1=beta_C1,beta_C0=beta_C0,beta_N=beta_N, 
                            psi_A=psi_A,psi_C=psi_C)
      }
      t=t+1
      
    }
    
    
    
    return(list(DID=DID,param=MCMC_list))
    
  }
  
  ### Example ####
  model_DA_IV<-DA_IV(Iteration=10000, Data=Datalist,Outcome="Gaussian",prior_values=list(prior_a=1,prior_b=1))
  
}
  

### Data Augmentation with common trends assumption ######
{
  ### Update Group Label: A, N, C ######
  Update.G.nh<-function(Y,D,Z,G,w_A, w_N, w_C, mu_A, mu_C, mu_N, alpha, alpha_C ,sigma_A, sigma_C, sigma_N,Strata.model=F){
    
    N<-length(Y)
    
    ## Draw compliers from (1,1) (mixture of always-takers and compliers)
    id_a_or_c<-which(Z==1&D==1)
    
    Y_a_or_c<-Y[id_a_or_c]
    muA1=mu_A[id_a_or_c]+alpha
    muC1=mu_C[id_a_or_c]+alpha+alpha_C
    Prob_YA1<-dnorm(Y_a_or_c,mean=muA1,sd=sqrt(sigma_A))
    Prob_YC1<-dnorm(Y_a_or_c,mean=muC1,sd=sqrt(sigma_C))
    
    
    if(Strata.model==T){
      Prob_C1<-w_C[id_a_or_c]*Prob_YC1/(w_A[id_a_or_c]*Prob_YA1+w_C[id_a_or_c]*Prob_YC1)
    }else{
      Prob_C1<-w_C*Prob_YC1/(w_A*Prob_YA1+w_C*Prob_YC1)
    }
    
    G_a_or_c<-ifelse(rbinom(length(id_a_or_c),1,Prob_C1),"C","A")
    
    ## Draw compliers from (0,0)(mixture of neve-takers and compliers)
    id_n_or_c<-which(Z==0&D==0)
    
    
    Y_n_or_c<-Y[id_n_or_c]
    muN0=mu_N[id_n_or_c]
    muC0=mu_C[id_n_or_c]
    Prob_YN0<-dnorm(Y_n_or_c,mean=muN0,sd=sqrt(sigma_N))
    Prob_YC0<-dnorm(Y_n_or_c,mean=muC0,sd=sqrt(sigma_C))
    
    
    if(Strata.model==T){
      Prob_C0<-w_C[id_n_or_c]*Prob_YC0/(w_N[id_n_or_c]*Prob_YN0+w_C[id_n_or_c]*Prob_YC0)
    }else{
      Prob_C0<-w_C*Prob_YC0/(w_N*Prob_YN0+w_C*Prob_YC0)
    }
    
    G_n_or_c<-ifelse(rbinom(length(id_n_or_c),1,Prob_C0),"C","N")
    G[id_a_or_c]<-G_a_or_c
    G[id_n_or_c]<-G_n_or_c
    
    return(G)
    
    
    
  }
  
  
  ### MCMC algorithm #
  DA_common<-function(Iteration, Data=list(Y,D,Z),
                      prior_values=list(prior_mu,prior_sd2,prior_beta_mu,prior_beta_sd2,prior_a,prior_b,prior_A,prior_C,prior_N),
                      ini_values=list(w_A_ini, w_N_ini, w_C_ini, 
                                      mu_A_ini, mu_C_ini, mu_N_ini, 
                                      beta0_A_ini, beta0_C_ini, beta0_N_ini, 
                                      beta_A_ini, beta_C_ini, beta_N_ini, 
                                      psi_A_ini,psi_C_ini,
                                      alpha_ini, alpha_C_ini,
                                      sigma_A_ini, sigma_C_ini, sigma_N_ini),
                      burn_in,thin){
    t=1
    # Load Data, initial values, prior parameters #
    {Y=Data$Y
      D=Data$D
      Z=Data$Z
      X=Data$X
      X.mat=cbind(rep(1),X)
      N=length(Y)
      
      # imputations
      Y_imp_04<-matrix(0,nrow=Iteration,ncol=N)
      Y_imp_09<-matrix(0,nrow=Iteration,ncol=N)
      DID<-rep(0,Iteration)
      DID_C04<-DID_C09<-rep(0,Iteration)
      
      # initial values
      mu_A_ini=ini_values$mu_A_ini
      mu_C_ini=ini_values$mu_C_ini
      mu_N_ini=ini_values$mu_N_ini
      
      w_A_ini=ini_values$w_A_ini
      w_N_ini=ini_values$w_N_ini
      w_C_ini=ini_values$w_C_ini
      
      psi_A_ini=ini_values$psi_A_ini
      psi_C_ini=ini_values$psi_C_ini
      
      alpha_ini=ini_values$alpha_ini
      alpha_C_ini=ini_values$alpha_C_ini
      
      sigma_A_ini=ini_values$sigma_A_ini
      sigma_C_ini=ini_values$sigma_C_ini
      sigma_N_ini=ini_values$sigma_N_ini
      
    
      psi_A<-array(0,dim=c(Iteration,ncol(X)+1))
      psi_C<-array(0,dim=c(Iteration,ncol(X)+1))
        

      # place holders
      G<-matrix(rep(NA),ncol=N,nrow=Iteration)
      G.ini<-array(NA,N)
      G.ini[which(D==1&Z==0)]<-"A"
      G.ini[which(D==0&Z==1)]<-"N"
      G.ini[which(D==1&Z==1)]<-sample(c("A","C"),length(which(D==1&Z==1)),replace = T)
      G.ini[which(D==0&Z==0)]<-sample(c("N","C"),length(which(D==1&Z==1)),replace = T)
      
      
      w_A<-rep(0,Iteration)
      w_C<-rep(0,Iteration)
      w_N<-rep(0,Iteration)
      
      mu_A<-rep(0,Iteration)
      mu_C<-rep(0,Iteration)
      mu_N<-rep(0,Iteration)
      
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
      psi0<-rep(0,ncol(X)+1)   #prior mean for psi
      T0<-diag(.01,ncol(X)+1)  # Prior precision of psi
      }
    
    # MCMC #
    while(t <= Iteration){
      
      if(t%%10000==0){
        print(paste("Iteration=",t,sep=""))
      }
      
      
      if(t==1){
        
        
        # sample D_mis #
        G[t,]<-Update.G.nh(Y,D,Z,G.ini,rep(w_A_ini,N), rep(w_N_ini,N), rep(w_C_ini,N), 
                           rep(mu_A_ini,N), rep(mu_C_ini,N), rep(mu_N_ini,N), alpha_ini, alpha_C_ini,
                           sigma_A_ini, sigma_C_ini, sigma_N_ini,Strata.model = Strata.model)
        
        id_A=which(G[t,]=="A")
        id_C=which(G[t,]=="C")
        id_N=which(G[t,]=="N")
        
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
        
        X.mat<-cbind(1,X)
        
        ##update MNL 
          ## Gibbs sampler update psi~ multinomial logistic regression ##
          # Define category-specific binary responses (Note: cat N is reference)
          u.a<-1*(G[t,]=="A") 
          u.c<-1*(G[t,]=="C")
          
          # Update Always-takers
          c.a<-log(1+exp(X.mat%*%psi_C_ini)) 
          eta.a<-X.mat%*%psi_A_ini-c.a
          w.a<-rpg(N,1,eta.a) 
          z.a<-(u.a-1/2)/w.a+c.a
          v<-solve(T0+crossprod(X.mat*sqrt(w.a)))
          m<-v%*%(T0%*%psi0+t(w.a*X.mat)%*%z.a) 
          psi_A[t,]<-c(rmvnorm(1,m,v))
          
          # Update Compliers
          c.c<-log(1+exp(X.mat%*%psi_A_ini)) 
          eta.c<-X.mat%*%psi_C_ini-c.c 
          w.c<-rpg(N,1,eta.c) 
          z.c<-(u.c-1/2)/w.c+c.c 
          v<-solve(T0+crossprod(X.mat*sqrt(w.c)))
          m<-v%*%(T0%*%psi0+t(w.c*X.mat)%*%z.c) 
          psi_C[t,]<-c(rmvnorm(1,m,v))
          
          
          D.xmat<-as.matrix(cbind(rep(1),X))
          lpA<-D.xmat%*%psi_A[t,]
          lpC<-D.xmat%*%psi_C[t,]
          
          w_A<-exp(lpA)/(1+exp(lpA)+exp(lpC))
          w_C<-exp(lpC)/(1+exp(lpA)+exp(lpC))
          w_N<-1-w_A-w_C
          
          ## Update coefficients of the outcome model
          ## Big design matrix ##
          Bigmat<-as.matrix(cbind(1*(G[t,]=="A")*X.mat,1*(G[t,]=="N")*X.mat,1*(G[t,]=="C")*cbind(X.mat,Z),Z))
          colnames(Bigmat)<-c("beta0A","beta1A","beta2A","beta0N","beta1N","beta2N","beta0C","beta1C","beta2C","alpha_C","alpha")
          p<-ncol(Bigmat)
          beta_ini<-rep(0,p)
          mu_beta<-rep(0,p)
          T0_beta<-diag(0.001,p)
          
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
        
        G[t,]<-Update.G.nh(Y,D,Z,G[t-1,],w_A, w_N, w_C,
                           mu_A, mu_C, mu_N, alpha[t-1], alpha_C[t-1],
                           sigma_A[t-1], sigma_C[t-1], sigma_N[t-1],Strata.model = Strata.model)
        
        
        # Update parameters # hierarchical model #conjugate priors
        id_A=which(G[t,]=="A")
        id_C=which(G[t,]=="C")
        id_N=which(G[t,]=="N")
        
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
        

          ## Gibbs sampler update psi~ multinomial logistic regression ##
          # Define category-specific binary responses (Note: cat N is reference)
          u.a<-1*(G[t,]=="A") 
          u.c<-1*(G[t,]=="C")
          
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
          
          psi_A[t,]<-psi_A.t
          psi_C[t,]<-psi_C.t
          
          
          # ## calculate marginal G probability
          D.xmat<-as.matrix(cbind(rep(1),X))
          lpA<-D.xmat%*%psi_A[t,]
          lpC<-D.xmat%*%psi_C[t,]
          
          w_A<-exp(lpA)/(1+exp(lpA)+exp(lpC))
          w_C<-exp(lpC)/(1+exp(lpA)+exp(lpC))
          w_N<-1-w_A-w_C#1/(1+exp(lpA)+exp(lpC))
          
        
          
          X_A<-X[id_A,]
          X_C<-X[id_C,]
          X_N<-X[id_N,]
          
          ## Update coefficients for outcome model 
          ## Big design matrix ##
          Bigmat<-as.matrix(cbind(1*(G[t,]=="A")*X.mat,1*(G[t,]=="N")*X.mat,1*(G[t,]=="C")*cbind(X.mat,Z),Z))
          colnames(Bigmat)<-c("beta0A","beta1A","beta2A","beta0N","beta1N","beta2N","beta0C","beta1C","beta2C","alpha_C","alpha")
          p<-ncol(Bigmat)
          mu_beta<-rep(0,p)
          T0_beta<-diag(0.01,p)
          
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
      
      ### imputation ###
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
      
      
      mu_pred_04=mu_A.p*(G[t,]=="A")+ mu_C.p*(G[t,]=="C")+ mu_N.p*(G[t,]=="N")  
      mu_pred_09=(mu_A.p+alpha[t])*(G[t,]=="A")+ (mu_C.p+alpha[t]+alpha_C[t])*(G[t,]=="C")+ (mu_N.p+alpha[t])*(G[t,]=="N") 
      sigma_pred=sigma_A[t]*(G[t,]=="A")+ sigma_C[t]*(G[t,]=="C")+ sigma_N[t]*(G[t,]=="N") 
      
      Y_imp_04[t,]<-rnorm(N,mu_pred_04,sd=sqrt(sigma_pred))
      Y_imp_09[t,]<-rnorm(N,mu_pred_09,sd=sqrt(sigma_pred))

      
      DID[t]=mean(Y_imp_09[t,which(G[t,]=="C")])-mean(Y_imp_04[t,which(G[t,]=="C")])-(mean(Y_imp_09[t,which(G[t,]!="C")])-mean(Y_imp_04[t,which(G[t,]!="C")]))
      DID_C04[t]<-mean(Y_imp_C1_04)-mean(Y_imp_C0_04)
      DID_C09[t]<-mean(Y_imp_C1_09)-mean(Y_imp_C0_09)
      

      t=t+1
      
   
    }
    
    
    CACE.all<-DID[seq(burn_in,Iteration,thin)]
    CACE<-mean(CACE.all,na.rm=T)
    CI<-quantile(CACE.all,c(0.025,0.975),na.rm=T)
    
    CACE.all.04<-DID_C04[seq(burn_in,Iteration,thin)]
    CACE.04<-mean(CACE.all.04,na.rm=T)
    CI.04<-quantile(CACE.all.04,c(0.025,0.975),na.rm=T)
    
    CACE.all.09<-DID_C09[seq(burn_in,Iteration,thin)]
    CACE.09<-mean(CACE.all.09,na.rm=T)
    CI.09<-quantile(CACE.all.09,c(0.025,0.975),na.rm=T)
    

    MCMC_list<-data.frame(alpha=alpha, alpha_C=alpha_C,
                              psiA=psi_A,psiC=psi_C,
                              beta_A0=beta_A0,beta_C0=beta_C0,beta_N0=beta_N0,
                              beta_A=beta_A,beta_C=beta_C,beta_N=beta_N,                            
                              sigma_A=sigma_A,sigma_N=sigma_N,sigma_C=sigma_C)
        
    
    return(list(param=MCMC_list[seq(burn_in,Iteration,thin),],
                cat.number=cat.number[seq(burn_in,Iteration,thin),],
                DID=CACE.all,CACE=CACE,CACE_CI=CI,
                DID.04=CACE.all.04,CACE.04=CACE.04,CACE_CI.04=CI.04,
                DID.09=CACE.all.09,CACE.09=CACE.09,CACE_CI.09=CI.09))

     }

  ### Example ###
  data_actual<-data.all$Data_2000
  data<-data_actual
  Datalist=list(Y=data$Y, D=data$D, Z=data$Z,X=cbind(data$X1,data$X2),X1=data$X1,X2=data$X2)
  
  prior_val=list(prior_mu=0,prior_sd2=10000,prior_beta_mu=0,prior_beta_sd2=10000,prior_a=1,prior_b=1,prior_A=10,prior_C=10,prior_N=10)
  
  ini_val=list(w_A_ini=1/3,w_C_ini=1/3,w_N_ini=1/3,
               mu_A_ini=mean(data_actual$Y[which(data_actual$Z==0&data_actual$D==1)]), 
               mu_C_ini=mean(data_actual$Y), 
               mu_N_ini=mean(data_actual$Y[which(data_actual$Z==1&data_actual$D==0)]),
               beta0_A_ini=mean(data_actual$Y[which(data_actual$Z==0&data_actual$D==1)]), 
               beta0_C_ini=mean(data_actual$Y), 
               beta0_N_ini=mean(data_actual$Y[which(data_actual$Z==1&data_actual$D==0)]), 
               beta_A_ini=c(0,0),beta_C_ini=c(0,0),beta_N_ini=c(0,0),
               psi_A_ini=c(1,2), psi_C_ini=c(1,2), 
               alpha_ini=2,
               alpha_C_ini=4,
               sigma_A_ini=10, sigma_C_ini=10, sigma_N_ini=10) 
  
  model_common<-DA_common(Iteration, Data=Datalist,Intercept.model=F,
                          prior_values=prior_val,
                          ini_values=ini_val,Strata.model=F,
                          burn_in=burn_in,thin=thin)
}


### Data Augmentation with weaker temporal trends assumption ######
{
  Update.G<-function(Y,D,Z,G, w_A, w_N, w_C, mu_A, mu_C, mu_N, 
                       alpha_A, alpha_N, alpha_C, Tau,sigma_A, sigma_C, sigma_N)
    {
      N<-length(Y)
      
      ## Draw compliers from (1,1) (mixture of always-takers and compliers)
      id_a_or_c<-which(Z==1&D==1)
      
      Y_a_or_c<-Y[id_a_or_c]
      muA1=mu_A[id_a_or_c]+alpha_A
      muC1=mu_C[id_a_or_c]+alpha_C+Tau
      Prob_YA1<-dnorm(Y_a_or_c,mean=muA1,sd=sqrt(sigma_A))
      Prob_YC1<-dnorm(Y_a_or_c,mean=muC1,sd=sqrt(sigma_C))

      Prob_C1<-w_C[id_a_or_c]*Prob_YC1/(w_A[id_a_or_c]*Prob_YA1+w_C[id_a_or_c]*Prob_YC1)
      G_a_or_c<-ifelse(rbinom(length(id_a_or_c),1,Prob_C1),"C","A")
      
      ## Draw compliers from (0,0)(mixture of neve-takers and compliers)
      id_n_or_c<-which(Z==0&D==0)
      
        Y_n_or_c<-Y[id_n_or_c]
        muN0=mu_N[id_n_or_c]
        muC0=mu_C[id_n_or_c]
        Prob_YN0<-dnorm(Y_n_or_c,mean=muN0,sd=sqrt(sigma_N))
        Prob_YC0<-dnorm(Y_n_or_c,mean=muC0,sd=sqrt(sigma_C))

        Prob_C0<-w_C[id_n_or_c]*Prob_YC0/(w_N[id_n_or_c]*Prob_YN0+w_C[id_n_or_c]*Prob_YC0)

      G_n_or_c<-ifelse(rbinom(length(id_n_or_c),1,Prob_C0),"C","N")
      G[id_a_or_c]<-G_a_or_c
      G[id_n_or_c]<-G_n_or_c
      
      return(G)
      
    }
    

    # MCMC algorithm 
    DA_weak<-function(Iteration, Data=list(Y,D,Z,X), 
                      prior_values=list(prior_A=10,prior_C1=10,prior_C2=10,prior_N=10, prior_a=1, prior_b=1,hyperparam=list(a_l=1,b_l=1,a_k=1,b_k=1),
                                        prior.sd,proposal_rt_df,proposal_rw_width),
                      ini_values=list(w_A_ini,w_C_ini,w_N_ini,
                                      beta_A0_ini, beta_N0_ini, beta_C0_ini, 
                                      beta_A_ini, beta_N_ini, beta_C_ini, 
                                      psi_A_ini, psi_C_ini,  
                                      alpha_A_ini, alpha_C_ini, alpha_N_ini,
                                      sigma2_ini,mu_alpha_var_ini,
                                      mu_alpha_ini,sigma_alpha_ini,Tau_ini),
                      burn_in,thin,printYes=F){
      
      t=1
      
      if(T){  
        # Data 
        Y=Data$Y
        D=Data$D
        Z=Data$Z
        X=Data$X
        N=length(Y)
        X.mat<-cbind(rep(1),X)
        
        # imputations
        Y_imp_04<-matrix(0,nrow=Iteration,ncol=N)
        Y_imp_09<-matrix(0,nrow=Iteration,ncol=N)
        
        
        # initial values
        w_A_ini=ini_values$w_A_ini
        w_N_ini=ini_values$w_N_ini
        w_C_ini=ini_values$w_C_ini
        
        beta_A0_ini=ini_values$beta_A0_ini
        beta_C0_ini=ini_values$beta_C0_ini
        beta_N0_ini=ini_values$beta_N0_ini
        
        beta_A_ini=as.matrix(ini_values$beta_A_ini,ncol=1)
        beta_C_ini=as.matrix(ini_values$beta_C_ini,ncol=1)
        beta_N_ini=as.matrix(ini_values$beta_N_ini,ncol=1)
        
        alpha_A_ini=ini_values$alpha_A_ini
        alpha_C_ini=ini_values$alpha_C_ini
        alpha_N_ini=ini_values$alpha_N_ini

        sigma2_ini=ini_values$sigma2_ini
        mu_alpha_ini<-ini_values$mu_alpha_ini
        sigma_alpha_ini<-ini_values$sigma_alpha_ini
        Tau_ini<-ini_values$Tau_ini
        mu_alpha_var_ini=ini_values$mu_alpha_var_ini
        
        ## hyper-parameters 
        a_l=prior_values$hyperparam$a_l #prior for mu.alpha.var
        b_l=prior_values$hyperparam$b_l
        a_k=prior_values$hyperparam$a_k #prior for sigma_alpha
        b_k=prior_values$hyperparam$b_k
        
        
        psi_A_ini<-matrix(ini_values$psi_A_ini,ncol=1)
        psi_C_ini<-matrix(ini_values$psi_C_ini,ncol=1)
        
        # place holders
        DID<-rep(0,Iteration)
        DID_C04<-rep(0,Iteration)
        DID_C09<-rep(0,Iteration)
        
        accept_sigma<-0
        
        G<-matrix(rep(NA),ncol=N,nrow=Iteration)
        G.ini<-array(NA,N)
        G.ini[which(D==1&Z==0)]<-"A"
        G.ini[which(D==0&Z==1)]<-"N"
        

        beta_A0<-rep(0,Iteration)
        beta_C0<-rep(0,Iteration)
        beta_N0<-rep(0,Iteration)
        
        beta_A<-array(0,dim=c(Iteration,ncol(X)))
        beta_C<-array(0,dim=c(Iteration,ncol(X)))
        beta_N<-array(0,dim=c(Iteration,ncol(X)))
        
        sigma2<-rep(1,Iteration)
 
        alpha_A<-rep(0,Iteration)
        alpha_C<-rep(0,Iteration)
        alpha_N<-rep(0,Iteration)
        
        psi_A<-array(0,dim=c(Iteration,ncol(X)+1))
        psi_C<-array(0,dim=c(Iteration,ncol(X)+1))
        
        Tau<-rep(0,Iteration)
        
        mu_alpha<-rep(0,Iteration)
        #sigma_alpha<-rep(sigma_alpha_ini,Iteration)
        sigma_alpha<-rep(1,Iteration)
        
        # Priors 
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
        psi0<-rep(0,ncol(X)+1)   #prior mean for psi
        T0<-diag(.01,ncol(X)+1)  # Prior precision of psi
        
      }
      
      while(t <= Iteration){
        
        if(printYes==T){
          
          if(t%%5000==0){
            print(paste("Iteration=",t,sep=""))
          }
        }
        
        if(t==1){
          
          mu_A_ini<-beta_A0_ini+X%*%beta_A_ini
          mu_N_ini<-beta_N0_ini+X%*%beta_N_ini
          mu_C_ini<-beta_C0_ini+X%*%beta_C_ini
          
          G[t,]<-Update.G(Y,D,Z,G.ini,w_A_ini, w_N_ini, w_C_ini, 
                          mu_A_ini, mu_C_ini, mu_N_ini, alpha_A_ini,alpha_N_ini, alpha_C_ini,Tau_ini,
                          sigma2_ini, sigma2_ini, sigma2_ini)
          
          id_A=which(G[t,]=="A")
          id_C=which(G[t,]=="C")
          id_N=which(G[t,]=="N")
          
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
  
            ## Gibbs sampler update psi~ multinomial logistic regression ##
            # Define category-specific binary responses (Note: cat N is reference)
            u.a<-1*(G[t,]=="A") 
            u.c<-1*(G[t,]=="C")
            
            # Update Always-takers
            c.a<-log(1+exp(X.mat%*%psi_C_ini)) 
            eta.a<-X.mat%*%psi_A_ini-c.a
            w.a<-rpg(N,1,eta.a) 
            z.a<-(u.a-1/2)/w.a+c.a
            v<-solve(T0+crossprod(X.mat*sqrt(w.a)))
            m<-v%*%(T0%*%psi0+t(w.a*X.mat)%*%z.a) 
            psi_A[t,]<-c(rmvnorm(1,m,v))
            
            # Update Compliers
            c.c<-log(1+exp(X.mat%*%psi_A_ini)) 
            eta.c<-X.mat%*%psi_C_ini-c.c 
            w.c<-rpg(N,1,eta.c) 
            z.c<-(u.c-1/2)/w.c+c.c 
            v<-solve(T0+crossprod(X.mat*sqrt(w.c)))
            m<-v%*%(T0%*%psi0+t(w.c*X.mat)%*%z.c) 
            psi_C[t,]<-c(rmvnorm(1,m,v))
            
            
            D.xmat<-as.matrix(cbind(rep(1),X))
            lpA<-D.xmat%*%psi_A[t,]
            lpC<-D.xmat%*%psi_C[t,]
            
            w_A<-exp(lpA)/(1+exp(lpA)+exp(lpC))
            w_C<-exp(lpC)/(1+exp(lpA)+exp(lpC))
            w_N<-1-w_A-w_C
          
          
          ### Update betas, tau ########
          ## matrix without the Z column 
          X.mat<-cbind(rep(1),X)
          q<-ncol(X.mat)
          Bigmat<-as.matrix(cbind(1*(G[t,]=="A")*X.mat,1*(G[t,]=="N")*X.mat,1*(G[t,]=="C")*X.mat,1*(G[t,]%in%c("C"))*D))
          colnames(Bigmat)<-c(paste("beta",c(0:q),"A",sep="_"),paste("beta",c(0:q),"N",sep="_"),paste("beta",c(0:q),"C",sep="_"),"tau")
          p<-ncol(Bigmat)
          beta_ini<-matrix(rep(0,p),ncol=1) #prior mu
          T0_beta<-diag(10000,p) #prior sigma
          
          ## calculate Ystar=Y-Z*temporal coef
          Y_star<-Y
          Y_star[id_A]<-Y[id_A]-Z[id_A]*alpha_A_ini
          Y_star[id_N]<-Y[id_N]-Z[id_N]*alpha_N_ini
          Y_star[id_C]<-Y[id_C]-Z[id_C]*alpha_C_ini
          
          ## update betas ##
          inv.Sigma0<-solve(T0_beta)
          M<-solve(inv.Sigma0+1/sigma2_ini*t(Bigmat)%*%Bigmat)%*%(inv.Sigma0%*%beta_ini+1/sigma2_ini*t(Bigmat)%*%Y_star)
          V<-solve(inv.Sigma0+1/sigma2_ini*t(Bigmat)%*%Bigmat)
          beta<-c(rmvnorm(1,M,V))
          beta_A0[t]<-beta[1]
          beta_A[t,]<-beta[2:3]
          beta_N0[t]<-beta[4]
          beta_N[t,]<-beta[5:6]
          beta_C0[t]<-beta[7]
          beta_C[t,]<-beta[8:9]
          Tau[t]<-beta[10]
          
          ## Update  mu_alpha ## (M=4, four groups)
          inv.sigma_alpha2<-solve(sigma_alpha_ini)
          alpha.sum<-alpha_A_ini+alpha_N_ini+alpha_C_ini
          
          mu.alpha<-solve(3*inv.sigma_alpha2+1/mu_alpha_var_ini)*inv.sigma_alpha2*alpha.sum
          sigma.alpha<-solve(3*inv.sigma_alpha2+1/mu_alpha_var_ini)
          mu_alpha[t]<-rmvnorm(1,mean=mu.alpha,sigma = sigma.alpha)
          
          Y_star<-Y-Bigmat%*%matrix(beta,ncol=1)
          
          temp.A<-solve(t(Z_A)%*%Z_A+inv.sigma_alpha2)%*%(t(Z_A)%*%Y_star[id_A]+inv.sigma_alpha2*mu_alpha[t])
          temp.N<-solve(t(Z_N)%*%Z_N+inv.sigma_alpha2)%*%(t(Z_N)%*%Y_star[id_N]+inv.sigma_alpha2*mu_alpha[t])
          temp.C<-solve(t(Z_C)%*%Z_C+inv.sigma_alpha2)%*%(t(Z_C)%*%Y_star[id_C]+inv.sigma_alpha2*mu_alpha[t])
          
          temp.sigma.A<-1/(t(Z_A)%*%Z_A+inv.sigma_alpha2) #sigma2_ini*
          temp.sigma.N<-1/(t(Z_N)%*%Z_N+inv.sigma_alpha2)
          temp.sigma.C<-1/(t(Z_C)%*%Z_C+inv.sigma_alpha2)
          
          alpha_A[t]<-rnorm(1,mean=temp.A,sd=sqrt(temp.sigma.A))
          alpha_N[t]<-rnorm(1,mean=temp.N,sd=sqrt(temp.sigma.N))
          alpha_C[t]<-rnorm(1,mean=temp.C,sd=sqrt(temp.sigma.C))
          
          ### Update mu.alpha.var
          mu_alpha_var[t]<-invgamma::rinvgamma(1,shape=a_l+1/2,rate=b_l+1/2*mu_alpha[t]^2)
          
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
          
          
          
        }else{
          
     
             G[t,]<-Update.G(Y,D,Z, G[t-1,],w_A, w_N, w_C,
                              mu_A, mu_C, mu_N, alpha_A[t-1], alpha_N[t-1],alpha_C[t-1],Tau[t-1],
                              sigma2[t-1], sigma2[t-1], sigma2[t-1])  
              
            
            
            id_A=which(G[t,]=="A")
            id_C=which(G[t,]=="C")
            id_N=which(G[t,]=="N")
            
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
            N=length(Y)
            
            X_A<-X[id_A,]
            X_C<-X[id_C,]
            X_N<-X[id_N,]
            

              ## Gibbs sampler update psi~ multinomial logistic regression ##
              # Define category-specific binary responses (Note: cat N is reference)
              u.a<-1*(G[t,]=="A") 
              u.c<-1*(G[t,]=="C")
              
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
              
              # ## calculate marginal G probability
              D.xmat<-as.matrix(cbind(rep(1),X))
              lpA<-D.xmat%*%psi_A[t,]
              lpC<-D.xmat%*%psi_C[t,]
              
              w_A<-exp(lpA)/(1+exp(lpA)+exp(lpC))
              w_C<-exp(lpC)/(1+exp(lpA)+exp(lpC))
              w_N<-1-w_A-w_C#1/(1+exp(lpA)+exp(lpC))
              

            X.mat<-cbind(rep(1),X)
            q<-ncol(X.mat)
            Bigmat<-as.matrix(cbind(1*(G[t,]=="A")*X.mat,1*(G[t,]=="N")*X.mat,1*(G[t,]=="C")*X.mat,1*(G[t,]=="C")*D))
            colnames(Bigmat)<-c(paste("beta",c(0:q),"A",sep="_"),paste("beta",c(0:q),"N",sep="_"),paste("beta",c(0:q),"C",sep="_"),"tau")
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
            
            ### Update temporal coefficients ####
            ## Update  mu_alpha ## (M=4, four groups)
            inv.sigma_alpha2<-solve(sigma_alpha[t-1])
            alpha.sum<-alpha_A[t-1]+alpha_N[t-1]+alpha_C[t-1]
            
            mu.alpha<-solve(3*inv.sigma_alpha2+1/mu_alpha_var[t-1])*inv.sigma_alpha2*alpha.sum
            sigma.alpha<-solve(3*inv.sigma_alpha2+1/mu_alpha_var[t-1])
            mu_alpha[t]<-rmvnorm(1,mean=mu.alpha,sigma = sigma.alpha)
            
            Y_star<-Y-Bigmat%*%matrix(beta,ncol=1)
            
            temp.A<-solve(t(Z_A)%*%Z_A+inv.sigma_alpha2)%*%(t(Z_A)%*%Y_star[id_A]+inv.sigma_alpha2*mu_alpha[t])
            temp.N<-solve(t(Z_N)%*%Z_N+inv.sigma_alpha2)%*%(t(Z_N)%*%Y_star[id_N]+inv.sigma_alpha2*mu_alpha[t])
            temp.C<-solve(t(Z_C)%*%Z_C+inv.sigma_alpha2)%*%(t(Z_C)%*%Y_star[id_C]+inv.sigma_alpha2*mu_alpha[t])
            
            temp.sigma.A<-1/(t(Z_A)%*%Z_A+inv.sigma_alpha2) #sigma2[t-1]*
            temp.sigma.N<-1/(t(Z_N)%*%Z_N+inv.sigma_alpha2)
            temp.sigma.C<-1/(t(Z_C)%*%Z_C+inv.sigma_alpha2)
            
            alpha_A[t]<-rnorm(1,mean=temp.A,sd=sqrt(temp.sigma.A))
            alpha_N[t]<-rnorm(1,mean=temp.N,sd=sqrt(temp.sigma.N))
            alpha_C[t]<-rnorm(1,mean=temp.C,sd=sqrt(temp.sigma.C))
            
            ### Update mu.alpha.var
            mu_alpha_var[t]<-invgamma::rinvgamma(1,shape=a_l+1/2,rate=b_l+1/2*mu_alpha[t]^2)
            
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
            

          
        }
        
        
        ### imputation ###
        mu_A<-beta_A0[t]+X%*%as.matrix(beta_A[t,],ncol=1)
        mu_N<-beta_N0[t]+X%*%as.matrix(beta_N[t,],ncol=1)
        mu_C<-beta_C0[t]+X%*%as.matrix(beta_C[t,],ncol=1)
        
        mu_C0_04<-beta_C0[t]+X_C%*%as.matrix(beta_C[t,],ncol=1)
        mu_C1_04<-beta_C0[t]+X_C%*%as.matrix(beta_C[t,],ncol=1)+Tau[t]
        
        mu_C0_09<-beta_C0[t]+X_C%*%as.matrix(beta_C[t,],ncol=1)+alpha_C[t]
        mu_C1_09<-beta_C0[t]+X_C%*%as.matrix(beta_C[t,],ncol=1)+alpha_C[t]+Tau[t]
        

        Y_imp_C0_04<-rnorm(N_C,mu_C0_04,sd=sqrt(sigma2[t]))
        Y_imp_C1_04<-rnorm(N_C,mu_C1_04,sd=sqrt(sigma2[t]))
        
        Y_imp_C0_09<-rnorm(N_C,mu_C0_09,sd=sqrt(sigma2[t]))
        Y_imp_C1_09<-rnorm(N_C,mu_C1_09,sd=sqrt(sigma2[t]))
        
        
        mu_pred_04=mu_A*(G[t,]=="A")+ mu_C*(G[t,]=="C")+ mu_N*(G[t,]=="N")  
        mu_pred_09=(mu_A)*(G[t,]=="A")+ (mu_C+Tau[t])*(G[t,]=="C")+ (mu_N)*(G[t,]=="N") 
        sigma_pred=sigma2[t]#*(G[t,]=="A")+ sigma_C[t]*(G[t,]=="C")+ sigma_N[t]*(G[t,]=="N") 
        
        Y_imp_04[t,]<-rnorm(N,mu_pred_04,sd=sqrt(sigma_pred))
        Y_imp_09[t,]<-rnorm(N,mu_pred_09,sd=sqrt(sigma_pred))
        

        DID[t]=mean(Y_imp_09[t,which(G[t,]=="C")])-mean(Y_imp_04[t,which(G[t,]=="C")])#-alpha_C[t]
        DID_C04[t]<-mean(Y_imp_C1_04)-mean(Y_imp_C0_04)
        DID_C09[t]<-mean(Y_imp_C1_09)-mean(Y_imp_C0_09)
          

        
        t=t+1
        
        
      }
      
      
          MCMC_list<-data.frame(alpha_A=alpha_A, alpha_N=alpha_N, alpha_C=alpha_C,
                                mu_alpha=mu_alpha, sigma_alpha=sigma_alpha,
                                Tau=Tau,
                                beta0_A=beta_A0,beta0_C=beta_C0,beta0_N=beta_N0, 
                                beta_A=beta_A,beta_C=beta_C,beta_N=beta_N, 
                                psi_A=psi_A,psi_C=psi_C, 
                                sigma2=sigma2)  

        CACE.all<-DID[seq(burn_in,Iteration,thin)]
        CACE<-mean(CACE.all,na.rm=T)
        CI<-quantile(CACE.all,c(0.025,0.975),na.rm=T)
        
        CACE.all.04<-DID_C04[seq(burn_in,Iteration,thin)]
        CACE.04<-mean(CACE.all.04,na.rm=T)
        CI.04<-quantile(CACE.all.04,c(0.025,0.975),na.rm=T)
        
        CACE.all.09<-DID_C09[seq(burn_in,Iteration,thin)]
        CACE.09<-mean(CACE.all.09,na.rm=T)
        CI.09<-quantile(CACE.all.09,c(0.025,0.975),na.rm=T)
        

      

      return(list(param=MCMC_list[seq(burn_in,Iteration,thin),],
                  DID=CACE.all,CACE=CACE,CACE_CI=CI,
                  DID.04=CACE.all.04,CACE.04=CACE.04,CACE_CI.04=CI.04,
                  DID.09=CACE.all.09,CACE.09=CACE.09,CACE_CI.09=CI.09,
                  accept_rates=accept_sigma/Iteration))
      
    }
    
  
    ### Example ###
    prior_val=list(prior_beta_mu=0,prior_beta_sd2=10000,
                   prior_A=10,prior_C=10,prior_N=10,
                   prior_tau_mu=0,prior_tau_sigma=10000,
                   prior_a=1, prior_b=1,
                   hyperparam=list(a_l=1,b_l=1,a_k=1,b_k=1),
                   A=1,proposal_width=1)
    
    ## strat value ##
    ini_val=list(w_A_ini=1/3,w_C_ini=1/3,w_N_ini=1/3,
                 beta_A0_ini=mean(data_actual$Y[which(data_actual$Z==0&data_actual$D==1)]), 
                 beta_C0_ini=mean(data_actual$Y), 
                 beta_N0_ini=mean(data_actual$Y[which(data_actual$Z==1&data_actual$D==0)]), 
                 beta_A_ini=c(0,0), beta_C_ini=c(0,0),beta_N_ini=c(0,0),
                 psi_A_ini=c(1,2), psi_C_ini=c(1,2), 
                 alpha_A_ini=4,alpha_C_ini=4,alpha_N_ini=4, 
                 sigma2_ini=10,mu_alpha_var_ini=1,
                 mu_alpha_ini=4,sigma_alpha_ini=1,Tau_ini=1)
    
    model_weak<-DA_weak(Iteration, Data=Datalist,
                        prior_values=prior_val,
                        ini_values=ini_val,
                        burn_in=burn_in,thin=thin,printYes=T)
  

}

