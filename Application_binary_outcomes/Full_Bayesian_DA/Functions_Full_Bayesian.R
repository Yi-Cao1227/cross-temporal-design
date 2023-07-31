###########################################################################
###                Functions for Application                           ###
###########################################################################


############### DA algorithm #########################
library(brms)
library(extraDistr)
library(invgamma)
library(MCMCglmm) 
library(MCMCpack)
library(mvtnorm)
library(maxLik)
library(BayesLogit)
library(mvtnorm)
library(truncnorm)
#library(truncnorm)


######### Functions file ################
# sd_alpha^2/sd_A^2=var.Ratio
invlogit<-function(x){
  
  exp(x)/(1+exp(x))
}

logit<-function(p){
  
  log(p/(1-p))
}

 
 ### Update Group Label with exclusion restriction ######
  Update.G.IV<-function(Y,D,Z,G, w_A, w_N, w_C, mu_C1, mu_C0,mu_A, mu_N){
    
    N<-length(Y)
    
    ## Draw compliers from (1,1) (mixture of always-takers and compliers)
    id_a_or_c<-which(Z==1&D==1)
    Y_a_or_c<-Y[id_a_or_c]
    
    Prob.muA1=invlogit(mu_A[id_a_or_c])
    Prob.muC1=invlogit(mu_C1[id_a_or_c])
    
    Prob_YA1<-dbinom(Y_a_or_c,size=1,prob=Prob.muA1)
    Prob_YC1<-dbinom(Y_a_or_c,size=1,prob=Prob.muC1)
    
    Prob_C1<-w_C[id_a_or_c]*Prob_YC1/(w_A[id_a_or_c]*Prob_YA1+w_C[id_a_or_c]*Prob_YC1)
    
    G_a_or_c<-ifelse(rbinom(length(id_a_or_c),1,Prob_C1)==1,"C","A")
    
    ## Draw compliers from (0,0)(mixture of neve-takers and compliers)
    id_n_or_c<-which(Z==0&D==0)
    Y_n_or_c<-Y[id_n_or_c]
    Prob.muN0=invlogit(mu_N[id_n_or_c])
    Prob.muC0=invlogit(mu_C0[id_n_or_c])
    Prob_YN0<-dbinom(Y_n_or_c,size=1,prob=Prob.muN0)
    Prob_YC0<-dbinom(Y_n_or_c,size=1,prob=Prob.muC0)  
    
    Prob_C0<-w_C[id_n_or_c]*Prob_YC0/(w_N[id_n_or_c]*Prob_YN0+w_C[id_n_or_c]*Prob_YC0)
    
    G_n_or_c<-ifelse(rbinom(length(id_n_or_c),1,Prob_C0),"C","N")
    
    G[id_a_or_c]<-G_a_or_c
    G[id_n_or_c]<-G_n_or_c
    
    return(G)  
    
  }
  
 ### Update Group Label without exclusion restriction ######
  Update.G.2t<-function(Y,D,Z,G, w_A, w_N, w_C, mu_A, mu_C, mu_N, 
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
  

  App_IV_DA<-function(Iteration, Data)
  {

    # Data 
    Y=Data$Y
    D=Data$D
    Z=Data$Z
    X=Data$X
    N=length(Y)
    
    # place holders
    DID<-rep(0,Iteration)
    
    #G<-matrix(rep(NA),ncol=N,nrow=Iteration)
    G<-array(NA,N)
    G[which(D==1&Z==0)]<-"A"
    G[which(D==0&Z==1)]<-"N"
    G[which(D==1&Z==1)]<-sample(c("A","C"),length(which(D==1&Z==1)),replace = T)
    G[which(D==0&Z==0)]<-sample(c("N","C"),length(which(D==0&Z==0)),replace = T)
    
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
      
      
      
      G<-Update.G.IV(Y,D,Z,G, w_A, w_N, w_C, mu_C1=mu_C1, mu_C0=mu_C0,  mu_A=mu_A, mu_N=mu_N)
      
      
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
      
      ## Big design matrix ##
      q<-ncol(X)
      Bigmat.A<-X.mat[which(G=="A"),]
      colnames(Bigmat.A)<-c(paste(colnames(X.mat),"A",sep="_"))
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
      colnames(Bigmat.N)<-c(paste(colnames(X.mat),"N",sep="_"))
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
      colnames(Bigmat.C1)<-c(paste(colnames(X.mat),"C1",sep="_"))
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
      colnames(Bigmat.C0)<-c(paste(colnames(X.mat),"C0",sep="_"))
      p<-ncol(Bigmat.C0)
      mu_beta<-rep(0,p)
      T0_beta<-diag(0.01,p)
      
      eta.C0<-Bigmat.C0%*%beta_C0[t-1,]
      w<-rpg(sum(G=="C"&Z==1),1,eta.C0)
      z<-(Y[which(G=="C"&Z==1)]-1/2)/w
      v<-solve(crossprod(Bigmat.C0*sqrt(w))+T0_beta) 
      m<-v%*%(T0_beta%*%mu_beta+t(w*Bigmat.C0)%*%z) 
      beta_C0[t,]<-c(rmvnorm(1,m,v))
      
      mu_A=X.mat%*%beta_A[t,]
      mu_N=X.mat%*%beta_N[t,]
      mu_C0=X.mat%*%beta_C0[t,]
      mu_C1=X.mat%*%beta_C1[t,]
      
      
      ### imputation of PO ###
      mu_C0.pred<-mu_C0[which(G=="C")]
      mu_C1.pred<-mu_C1[which(G=="C")]
      
      Y_imp_C0<-rbinom(length(mu_C0.pred),1,invlogit(mu_C0.pred))
      Y_imp_C1<-rbinom(length(mu_C1.pred),1,invlogit(mu_C1.pred))
      
      
      #predicted odds ratio
      y.C0<-mean(Y_imp_C0)
      y.C1<-mean(Y_imp_C1)
      
      
      DID[t]=((y.C1/(1-y.C1))/(y.C0/(1-y.C0)))
      
      t=t+1
      
    }
    
    MCMC_list<-data.frame(beta_A=beta_A,beta_C1=beta_C1,beta_C0=beta_C0,beta_N=beta_N, 
                          psi_A=psi_A,psi_C=psi_C)
    
    return(list(DID=DID,param=MCMC_list))
    
  }
  
  App_Cross_temporal_DA<-function(Iteration, Data=list(Y,D,Z,X),
                                  prior_values=list(prior_A=10,prior_C=10,prior_C1=10,prior_C2=10,prior_N=10, prior_a=1, prior_b=1,hyperparam=list(a_l=1,b_l=1,a_k=1,b_k=1)),
                                  ini_values=list(w_A_ini=1/4,w_C_ini=1/3,w_N_ini=1/4,w_C1_ini=1/4,w_C2_ini=1/4,
                                                  psi_A_ini=c(0,0,0),psi_C_ini=c(0,0,0), psi_C1_ini=c(0,0,0),psi_C2_ini=c(0,0,0),  
                                                  alpha_ini=c(1,1),
                                                  Tau_ini=c(5,5)),
                                  burn_in,thin,printYes=T,Model.name){
    
    t=1 
    
      y.mean09<-y.mean04<-rep(NA,Iteration)
      # Data 
      Y=Data$Y
      D=Data$D
      Z=Data$Z
      X=Data$X
      N=length(Y)
      

      # initial values
      w_A_ini=ini_values$w_A_ini
      w_N_ini=ini_values$w_N_ini
      w_C_ini=ini_values$w_C_ini
      alpha_ini=ini_values$alpha_ini
      Tau_ini<-ini_values$Tau_ini
      
      # place holders
      DID<-rep(0,Iteration)
      DID_t1<-rep(0,Iteration)
      DID_t2<-rep(0,Iteration)
      
      
      #G<-matrix(rep(NA),ncol=N,nrow=Iteration)
      G.ini<-array(NA,N)
      G.ini[which(D==1&Z==0)]<-"A"
      G.ini[which(D==0&Z==1)]<-"N"
      G.ini[which(D==1&Z==1)]<-sample(c("A","C"),length(which(D==1&Z==1)),replace = T)
      G.ini[which(D==0&Z==0)]<-sample(c("N","C"),length(which(D==0&Z==0)),replace = T)
      
      ## calculate initial values
      temp.A<-as.data.frame(cbind(Y,X)[which(G.ini=="A"),])
      fit.A<-glm(Y~.,data=temp.A,family = binomial("probit"))
      beta_A_ini<-matrix(coef(fit.A)[1:(1+ncol(X))] ,ncol=1)    
      alpha_ini<-ini_values$alpha_ini#coef(fit.A)[2+ncol(X)]
      
      temp.N<-as.data.frame(cbind(Y,X)[which(G.ini=="N"),])
      fit.N<-glm(Y~.,data=temp.N,family = binomial("probit"))
      beta_N_ini<-matrix(coef(fit.N)[1:(1+ncol(X))],ncol=1)   
      
      temp.C<-as.data.frame(cbind(Y,X)[which(G.ini=="C"),])
      fit.C<-glm(Y~.,data=temp.C,family = binomial("probit"))
      beta_C_ini<-matrix(coef(fit.C)[1:(1+ncol(X))] ,ncol=1) 
      Tau_ini<-ini_values$Tau_ini
      
      beta_A<-array(0,dim=c(Iteration,ncol(X)+1))
      beta_C<-array(0,dim=c(Iteration,ncol(X)+1))
      beta_N<-array(0,dim=c(Iteration,ncol(X)+1))
      
      alpha<-rep(0,Iteration)
      
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

    while(t <= Iteration){
      
      if(printYes==T){
        
        if(t%%5000==0){
          print(paste("Iteration=",t,sep=""))
        }
      }
      
      
      if(t==1){
        
        
        X.mat<-cbind(1,X)
        
        
        ## update Group label ###
          mu_A_ini<-X.mat%*%beta_A_ini
          mu_N_ini<-X.mat%*%beta_N_ini
          mu_C_ini<-X.mat%*%beta_C_ini
          
          
          G<-Update.G.2t(Y,D,Z,G.ini, w_A_ini, w_N_ini, w_C_ini, mu_A_ini, mu_C_ini, mu_N_ini, 
                         alpha_A=alpha_ini, alpha_N=alpha_ini, alpha_C=alpha_ini, Tau=Tau_ini,Strata.model=T) 
          
          id.A<-which(G=="A")
          id.N<-which(G=="N")
          id.C<-which(G=="C")
          
        
        ##### Update parameters for Stratification model #####
          ## Gibbs sampler update psi~ multinomial logistic regression ##
          # Define category-specific binary responses (Note: cat N is reference)
          u.a<-1*(G=="A") 
          u.c<-1*(G=="C")
          
          # Update Always-takers
          c.a<-log(1+exp(X.mat%*%psi_C[1,])) 
          eta.a<-X.mat%*%psi_A[1,]-c(c.a)
          w.a<-rpg(N,1,eta.a) 
          z.a<-(u.a-1/2)/w.a+c(c.a)
          v<-solve(T0+crossprod(X.mat*sqrt(w.a)))
          m<-v%*%(T0%*%psi0+t(w.a*X.mat)%*%z.a) 
          psi_A[t,]<-c(rmvnorm(1,m,v))
          
          # Update Compliers
          c.c<-log(1+exp(X.mat%*%psi_A[t,])) 
          eta.c<-X.mat%*%psi_C[1,]-c(c.c) 
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
          
        #### Update parameters for Outcome model ######
          Gt<-c(G)
          ## Big design matrix ##
          q<-ncol(X)
          Bigmat<-as.matrix(cbind(1*(Gt=="A")*X.mat,1*(Gt=="N")*X.mat,1*(Gt=="C")*cbind(X.mat,Z),Z))
          
          #Bigmat<-as.matrix(cbind(1*(Gt=="A")*X.mat,1*(Gt=="N")*X.mat,1*(Gt=="C")*X.mat,Z,1*(Gt=="C")*D))
          colnames(Bigmat)<-c(paste("beta",c(0:q),"A",sep="_"),paste("beta",c(0:q),"N",sep="_"),paste("beta",c(0:q),"C",sep="_"),"tau","alpha")
          p<-ncol(Bigmat)
          beta_ini<-c(beta_A_ini,beta_N_ini,beta_C_ini,Tau_ini,alpha_ini)
          
          ##logistic regression
          if(F){
          mu_beta<-rep(0,p)
          T0_beta<-diag(0.01,p)
          
          eta<-Bigmat%*%beta_ini
          w<-rpg(N,1,eta)
          z<-(Y-1/2)/w
          v<-solve(crossprod(Bigmat*sqrt(w))+T0_beta) 
          m<-v%*%(T0_beta%*%mu_beta+t(w*Bigmat)%*%z) 
          beta<-c(rmvnorm(1,m,v))
          }
          
          
          ## probit regression ##
          {
            mu_beta<-rep(0,p)
            T0_beta<- diag(0.1, p)
            
            ## augmented latent variable
            lat.y<-rep(0, N)
            V <- solve(T0_beta + crossprod(Bigmat, Bigmat))
            
            # Update Mean of lat.y
            mu_lat.y <- Bigmat %*% beta_ini
            # Draw latent variable z from its full conditional: z | \theta, y, X
            lat.y[Y == 0] <- rtruncnorm(sum(Y==0), mean = mu_lat.y[Y == 0], sd = 1, a = -Inf, b = 0)
            lat.y[Y == 1] <- rtruncnorm(sum(Y==1), mean = mu_lat.y[Y == 1], sd = 1, a = 0, b = Inf)
            
            # Compute posterior mean of beta
            M <- V %*% (T0_beta %*% mu_beta + crossprod(Bigmat, lat.y))
            # Draw variable \theta from its full conditional: \theta | z, X
            beta <- c(rmvnorm(1, M, V))
            
          }
          
          
          
          beta_A[t,]<-beta[1:(q+1)]
          beta_N[t,]<-beta[(q+2):(2*(q+1))]
          beta_C[t,]<-beta[(2*(q+1)+1):(3*(q+1))]
          Tau[t]<-beta[p-1]
          alpha[t]<-beta[p]
          
          mu_A=X.mat%*%beta_A[t,]
          mu_N=X.mat%*%beta_N[t,]
          mu_C=X.mat%*%beta_C[t,]
          

      }else{
        
        ## update Group label ###

          G<-Update.G.2t(Y,D,Z,G,w_A, w_N, w_C, mu_A, mu_C, mu_N,  
                         alpha_A=alpha[t-1],alpha_N=alpha[t-1], alpha_C=alpha[t-1],Tau[t-1],Strata.model = T)
          

          
          id.A<-which(G=="A")
          id.N<-which(G=="N")
          id.C<-which(G=="C")
          
        ## Update stratification model parameters ###
          u.a<-1*(G=="A") 
          u.c<-1*(G=="C")
          
          psi_A.t<-psi_A[t-1,]
          psi_C.t<-psi_C[t-1,]
          
          # Update Always-takers
          c.a<-log(1+exp(X.mat%*%psi_C.t)) 
          eta.a<-X.mat%*%psi_A.t-c(c.a)
          w.a<-rpg(N,1,eta.a) 
          z.a<-(u.a-1/2)/w.a+c(c.a)
          v<-solve(T0+crossprod(X.mat*sqrt(w.a)))
          m<-v%*%(T0%*%psi0+t(w.a*X.mat)%*%z.a) 
          psi_A.t<-c(rmvnorm(1,m,v))
          
          # Update Compliers
          c.c<-log(1+exp(X.mat%*%psi_A.t)) 
          eta.c<-X.mat%*%psi_C.t-c(c.c) 
          w.c<-rpg(N,1,eta.c) 
          z.c<-(u.c-1/2)/w.c+c(c.c) 
          v<-solve(T0+crossprod(X.mat*sqrt(w.c)))
          m<-v%*%(T0%*%psi0+t(w.c*X.mat)%*%z.c) 
          psi_C.t<-c(rmvnorm(1,m,v))

          psi_A[t,]<-psi_A.t
          psi_C[t,]<-psi_C.t
          
          ## calculate marginal G probability
          D.xmat<-as.matrix(cbind(rep(1),X))
          lpA<-D.xmat%*%psi_A[t,]
          lpC<-D.xmat%*%psi_C[t,]
          
          w_A<-exp(lpA)/(1+exp(lpA)+exp(lpC))
          w_C<-exp(lpC)/(1+exp(lpA)+exp(lpC))
          w_N<-1/(1+exp(lpA)+exp(lpC))
          
        ## Update Outcome model parameters ####
        Gt<-c(G)
          
          ## Big design matrix ##
          q<-ncol(X)
          Bigmat<-as.matrix(cbind(1*(Gt=="A")*X.mat,1*(Gt=="N")*X.mat,1*(Gt=="C")*cbind(X.mat,Z),Z))
          
          #Bigmat<-as.matrix(cbind(1*(Gt=="A")*X.mat,1*(Gt=="N")*X.mat,1*(Gt=="C")*X.mat,Z,1*(Gt=="C")*D))
          colnames(Bigmat)<-c(paste("beta",c(0:q),"A",sep="_"),paste("beta",c(0:q),"N",sep="_"),paste("beta",c(0:q),"C",sep="_"),"tau","alpha")
          p<-ncol(Bigmat)
          
          ## logistic regression
          if(F){
          mu_beta<-rep(0,p) #prior
          T0_beta<-diag(0.01,p)
          
          eta<-Bigmat%*%beta
          w<-rpg(N,1,eta)
          z<-(Y-1/2)/w
          v<-solve(crossprod(Bigmat*sqrt(w))+T0_beta) 
          m<-v%*%(T0_beta%*%mu_beta+t(w*Bigmat)%*%z) 
          beta<-c(rmvnorm(1,m,v))
          }
          
          ## probit regression ##
          {
            mu_beta<-rep(0,p)
            T0_beta<- diag(0.1, p)
            
            ## augmented latent variable
            lat.y<-rep(0, N)
            V <- solve(T0_beta + crossprod(Bigmat, Bigmat))
            
            # Update Mean of lat.y
            mu_lat.y <- Bigmat %*% beta_ini
            # Draw latent variable z from its full conditional: z | \theta, y, X
            lat.y[Y == 0] <- rtruncnorm(sum(Y==0), mean = mu_lat.y[Y == 0], sd = 1, a = -Inf, b = 0)
            lat.y[Y == 1] <- rtruncnorm(sum(Y==1), mean = mu_lat.y[Y == 1], sd = 1, a = 0, b = Inf)
            
            # Compute posterior mean of beta
            M <- V %*% (T0_beta %*% mu_beta + crossprod(Bigmat, lat.y))
            # Draw variable \theta from its full conditional: \theta | z, X
            beta <- c(rmvnorm(1, M, V))
            
          }
          
          
          
          beta_A[t,]<-beta[1:(q+1)]
          beta_N[t,]<-beta[(q+2):(2*(q+1))]
          beta_C[t,]<-beta[(2*(q+1)+1):(3*(q+1))]
          alpha[t]<-beta[p]
          Tau[t]<-beta[p-1]
          
          mu_A=X.mat%*%beta_A[t,]
          mu_N=X.mat%*%beta_N[t,]
          mu_C=X.mat%*%beta_C[t,]
          
        
      }
      
      
      ## Imputation ##

        N_C<-sum(Gt=="C")
        
        X_C.mat<-cbind(1,X[which(G=="C"),])

        mu_C.pred<-X.mat%*%as.matrix(beta_C[t,],ncol=1)
      
        mu_C0_04<-X_C.mat%*%as.matrix(beta_C[t,],ncol=1)
        mu_C1_04<-X_C.mat%*%as.matrix(beta_C[t,],ncol=1)+Tau[t]
        
        mu_C0_09<-X_C.mat%*%as.matrix(beta_C[t,],ncol=1)+alpha[t]
        mu_C1_09<-X_C.mat%*%as.matrix(beta_C[t,],ncol=1)+alpha[t]+Tau[t]
        
        Y_imp_C0_04<-rbinom(N_C,1,pnorm(lat.y[which(Gt=="C")],mean=mu_C0_04,sd=1))
        Y_imp_C1_04<-rbinom(N_C,1,pnorm(lat.y[which(Gt=="C")],mean=mu_C1_04,sd=1))
        
        Y_imp_C0_09<-rbinom(N_C,1,pnorm(lat.y[which(Gt=="C")],mean=mu_C0_09,sd=1))
        Y_imp_C1_09<-rbinom(N_C,1,pnorm(lat.y[which(Gt=="C")],mean=mu_C1_09,sd=1))
        
        
        #predicted odds ratio
        DID_t1[t]<-(mean(Y_imp_C1_04)/(1-mean(Y_imp_C1_04)))/(mean(Y_imp_C0_04)/(1-mean(Y_imp_C0_04)))
        DID_t2[t]<-(mean(Y_imp_C1_09)/(1-mean(Y_imp_C1_09)))/(mean(Y_imp_C0_09)/(1-mean(Y_imp_C0_09)))
        
      
      t=t+1
      
      
      if(t%%5000==0){
        
        if(Num.time.pts=="Two"){
          temp<-data.frame(alpha=alpha[1:t],
                           Tau=Tau[1:t],DID_t1=DID_t1[1:t],DID_t2=DID_t2[1:t],
                           beta_A=beta_A[1:t,],beta_C=beta_C[1:t,],beta_N=beta_N[1:t,], 
                           psi_A=psi_A[1:t,],psi_C=psi_C[1:t,])
        }else{
          
          temp<-data.frame(alpha=alpha,
                           Tau=Tau,DID_t2=DID_t2,DID_t3=DID_t3,
                           beta_A=beta_A,beta_C1=beta_C1,beta_C2=beta_C2,beta_N=beta_N,
                           psi_A=psi_A,psi_C1=psi_C1,psi_C2=psi_C2,sigma2=sigma2)
          
          
        }
        
        save(temp,file=paste0(Model.name,"_temp_result_",t,".Rdata"))
        
      }
      
      
    }
    
    
    ## results ####
    if(Num.time.pts=="Two"){
      
      MCMC_list<-data.frame(alpha=alpha,
                            Tau=Tau,DID_t1=DID_t1,DID_t2=DID_t2,
                            beta_A=beta_A,beta_C=beta_C,beta_N=beta_N, 
                            psi_A=psi_A,psi_C=psi_C)
      
    }else{
      
      MCMC_list<-data.frame(alpha=alpha,
                            Tau=Tau,DID_t2=DID_t2,DID_t3=DID_t3,
                            beta_A=beta_A,beta_C1=beta_C1,beta_C2=beta_C2,beta_N=beta_N,
                            psi_A=psi_A,psi_C1=psi_C1,psi_C2=psi_C2,sigma2=sigma2)
      
      
    }
    
    
    
    return(list(MCMC_list,G=G))
    
    
  }
  












