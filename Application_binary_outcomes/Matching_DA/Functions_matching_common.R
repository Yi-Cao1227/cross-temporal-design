
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
library(truncnorm)


######### Functions file ################
# sd_alpha^2/sd_A^2=var.Ratio
invlogit<-function(x){
  
  exp(x)/(1+exp(x))
}

logit<-function(p){
  
  log(p/(1-p))
}



######### Matching Method  ###############

MCMC_binary_Cross_matching<-function(Iteration, Data=list(Y,D,Z,X,G),
                                     temporal.assumption="common",
                                     prior_values=list(prior.sd,A, proposal_rt_df,proposal_rw_width),
                                     ini_values=list(beta_A_ini, beta_N_ini, beta_C_ini, 
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
    
    # imputations
    
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
    
    
    
    beta_A<-array(0,dim=c(Iteration,ncol(X)+1))
    beta_C<-array(0,dim=c(Iteration,ncol(X)+1))
    beta_N<-array(0,dim=c(Iteration,ncol(X)+1))
    
    sigma_A<-rep(1,Iteration)
    sigma_C<-rep(1,Iteration)
    sigma_N<-rep(1,Iteration)
    
      alpha_ini=ini_values$alpha_ini
      alpha<-rep(0,Iteration)
      
    
    
    Tau<-rep(0,Iteration)
    
    # Priors 
    prior_A_scale=prior_values$A
    prior.sd=prior_values$prior.sd
    proposal_rw_width=prior_values$proposal_rw_width
    proposal_rt_df=prior_values$proposal_rt_df
    
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
      fit0.A<-glm(Y_A~.,data=ini_dataA,family = binomial("probit"))
      
      ini_dataN<-data.frame(cbind(Y_N,X_N))
      fit0.N<-glm(Y_N~.,data=ini_dataN,family = binomial("probit"))
      
      ini_dataC<-data.frame(cbind(Y_C,X_C))
      fit0.C<-glm(Y_C~.,data=ini_dataC,family = binomial("probit"))
      
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
      

       
        ## Big design matrix ##
        Bigmat<-as.matrix(cbind(1*(G=="A")*X.mat,1*(G=="N")*X.mat,1*(G=="C")*cbind(X.mat,Z),Z))
        colnames(Bigmat)<-c(paste("beta_A", colnames(X.mat),sep="_"),paste("beta_N", colnames(X.mat),sep="_"),paste("beta_C", colnames(X.mat),sep="_"),"tau","alpha")
        p<-ncol(Bigmat)
        beta_ini<-rep(0,p)
        
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
        

        beta_A[t,]<-beta[1:ncol(X.mat)]
        beta_N[t,]<-beta[(1+ncol(X.mat)):(2*ncol(X.mat))]
        beta_C[t,]<-beta[(1+2*ncol(X.mat)):(3*ncol(X.mat))]
        Tau[t]<-beta[3*ncol(X.mat)+1]
        alpha[t]<-beta[3*ncol(X.mat)+2]

      
      
      
    }else{
      

        X_A.mat<-cbind(rep(1),X_A)
        X_N.mat<-cbind(rep(1),X_N)
        X_C.mat<-cbind(rep(1),X_C)
        
          
          ## Big design matrix ##
          Bigmat<-as.matrix(cbind(1*(G=="A")*X.mat,1*(G=="N")*X.mat,1*(G=="C")*cbind(X.mat,Z),Z))
          colnames(Bigmat)<-c(paste("beta_A", colnames(X.mat),sep="_"),paste("beta_N", colnames(X.mat),sep="_"),paste("beta_C", colnames(X.mat),sep="_"),"tau","alpha")
          p<-ncol(Bigmat)
          
          if(F){
          beta_ini<-rep(0,p)
          mu_beta<-rep(0,p)
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
          
          beta_A[t,]<-beta[1:ncol(X.mat)]
          beta_N[t,]<-beta[(1+ncol(X.mat)):(2*ncol(X.mat))]
          beta_C[t,]<-beta[(1+2*ncol(X.mat)):(3*ncol(X.mat))]
          Tau[t]<-beta[3*ncol(X.mat)+1]
          alpha[t]<-beta[3*ncol(X.mat)+2]
        

      
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
      
      mu_C0_09<-X_C.mat%*%as.matrix(beta_C[t,],ncol=1)+alpha[t]
      mu_C1_09<-X_C.mat%*%as.matrix(beta_C[t,],ncol=1)+alpha[t]+Tau[t]
      
      mu_pred_04=mu_A.pred*(G=="A")+ mu_C.pred*(G=="C")+ mu_N.pred*(G=="N")  
      mu_pred_09=mu_A.pred*(G=="A")+ (mu_C.pred+Tau[t])*(G=="C")+ mu_N.pred*(G=="N") 
      
    
    Y_imp_C0_04<-rbinom(N_C,1,pnorm(lat.y[which(Gt=="C")],mean=mu_C0_04,sd=1))
    Y_imp_C1_04<-rbinom(N_C,1,pnorm(lat.y[which(Gt=="C")],mean=mu_C1_04,sd=1))
    
    Y_imp_C0_09<-rbinom(N_C,1,pnorm(lat.y[which(Gt=="C")],mean=mu_C0_09,sd=1))
    Y_imp_C1_09<-rbinom(N_C,1,pnorm(lat.y[which(Gt=="C")],mean=mu_C1_09,sd=1))
    
    
      DID_C04[t]<-(mean(Y_imp_C1_04)/(1-mean(Y_imp_C1_04)))/(mean(Y_imp_C0_04)/(1-mean(Y_imp_C0_04)))
      DID_C09[t]<-(mean(Y_imp_C1_09)/(1-mean(Y_imp_C1_09)))/(mean(Y_imp_C0_09)/(1-mean(Y_imp_C0_09)))
      
    
    t=t+1
    
    
  }
  
  
  ## results
  
    MCMC_list<-data.frame(alpha=alpha, 
                          Tau=Tau,
                          beta_A=beta_A,beta_C=beta_C,beta_N=beta_N)


  return(list(param=MCMC_list,
              DID=DID,
              DID.04=DID_C04,
              DID.09=DID_C09))
  
  
  
  
}



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



Cross_binary_mathching<-function(Datalist,Potential.outcome=F,Iteration,temporal.assumption="common", binary.model="logit",
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
    
    #data_model$G_pred<-factor(data_model$G_pred,levels = c("Always_taker","Compliers","Never_taker"),labels=c("A","C","N"))
    
    #data_PO<-data_model
    
    #DX.mat<-model.matrix(~.-1,data=data_PO[,which(colnames(data_PO)%in%var.choose)])[,-5]
    
    #Datalist=list(Y=data_PO$Y, D=data_PO$D, Z=data_PO$Z,X=DX.mat,G=data_PO$G_pred)
    
    if(temporal.assumption=="common"){
      
      
      #if(binary.model=="logit"){
      
      #### Method3.1: Potential outcome framework + common temporal assumption #####
      MCMC_model<-MCMC_binary_Cross_matching(Iteration=Iteration, Data=Datalist,
                                             temporal.assumption="common",
                                             prior_values=prior_values,
                                             ini_values=ini_values,
                                             burn_in=burn_in,thin=thin,printYes=T)
      
      
      
      
      
      
    }else if(temporal.assumption=="weak"){
      
      MCMC_model<-MCMC_binary_Cross_matching(Iteration=Iteration, Data=Datalist,
                                             temporal.assumption="weak",
                                             prior_values=prior_values,
                                             ini_values=ini_values,
                                             burn_in=burn_in,thin=thin,printYes=T)
      
    }
    
    
    
   
    
    ##quantile
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














