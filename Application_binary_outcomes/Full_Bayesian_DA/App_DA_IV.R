source("Functions_Full_Bayesian.R")


#load("data_for_run_2017_latest.Rdata")

model_DA_IV_latest<-App_IV_DA(Iteration=20000,Data=Datalist)

#save(model_DA_IV_latest,file="App_DA_iv_result.Rdata")
