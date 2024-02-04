set.seed(1234)
#setwd("/blue/daniels/s.bhandari/ARIC/Code here")
install.packages('GcompBART2_0.1.0.tar.gz', repos=NULL, type='source')
library("GcompBART2")

library(coda)
library(haven)
library(plyr)
library(tidyverse)
library(nlme)
library(lme4)
library(dplyr)
library(data.table)
library(car)
library(ggplot2)
library(tidyr)


# Code for G-Computation

source("pred_comb_mediator.R")
source("pred_surv.R")
source("pwbart_it.R")
source("expit.R")
source("quiet.R")
source("pred_comb_Y.R")
source("pred_comb.R")
source("pred_y.R")




#.....................Bayesian Bootstrap Function......................//
library(gtools)

getSampleByBB = function(input_df , C_star , BBweightsIT){ #C_star = size of pseudo data
  BBdraws = rmultinom(1,C_star,BBweightsIT) # vector of number of balls in each box: put (by replacement) 'C_star' balls in length(BBweightsIT) boxes
  return_df = as.data.frame(lapply(input_df, rep, BBdraws)) #repeat i^th row of input_df, BBdraw[i] times
  return(return_df)
}




#.......................Survival Function............................//
get_S_zzstar = function(P_mat){
  s_temp <- apply(P_mat, 2, function(x){1-x})
  s_mat <- apply(s_temp, 2, cumprod)
  return(s_mat)
}




#.......................Cumulative Incidence Function.......................//
get_F_k_zzstar = function(P1_mat, P2_mat,k){
  l <- 1
  s1_temp <- apply(P1_mat, 2, function(x){1-x})
  s2_temp <- apply(P2_mat, 2, function(x){1-x})
  s_tilde_mat <- apply(s1_temp*s2_temp, 2, cumprod)
  f_k_mat <-matrix(0, nrow = 1, ncol = ncol(s_tilde_mat))
  if(k==1){
    
    f_k_mat <- rbind(f_k_mat,s_tilde_mat[1:nrow(s_tilde_mat)-1,] * P1_mat[2:nrow(P1_mat),])
  }else if(k==2){
    f_k_mat <- rbind(f_k_mat,s_tilde_mat[1:nrow(s_tilde_mat)-1,] * P2_mat[2:nrow(P2_mat),])
  }
  f_k_return <- apply(f_k_mat,2,cumsum)
  return(f_k_return)
}


get_IDE <- function(s_zzstar, s_zstarzstar){  #Both arguments must be matrices of same dimensions
  return(s_zzstar - s_zstarzstar)
}
get_IIE <- function(s_zz,s_zzstar){           #Both arguments must be matrices of same dimensions
  return(s_zz-s_zzstar)
}
get_TotalEffects <- function(IDE_mat, IIE_mat){       #Both arguments must be matrices of same dimensions
  return(IDE_mat + IIE_mat)
}






calculate_95bayesian_credible_interval <- function(data) {
  cat("Posterior Means: \n")
  cred_int = function(data_row){
    posterior_samples <- data_row
    posterior_means <- as.numeric(format(round(mean(posterior_samples) , 4), nsmall = 4))
    print(posterior_means)
    quantiles <- quantile(posterior_samples, c(0.025, 0.975))
    interval <- as.numeric(format(round(quantiles, 4), nsmall = 4))
    return(interval)
  }
  
  return_mat = apply(data,1,cred_int)
  return(t(return_mat))
  
} 




plot_ConfInt <-function(conf_int_mat) { 
  conf_int_df = data.frame(conf_int_mat)
  conf_int_df$CI <- 1:nrow(conf_int_df)
  ci_long <- gather(conf_int_df, key = "Bound", value = "Value", -CI)
  ggplot(ci_long, aes(x = CI, ymin = Value, ymax = Value, color = Bound )) +
    geom_errorbar(width = 0.2)+
    labs(x = "Visits", y = "Confidence/Credible Intervals")+
    scale_color_manual(values = c("lower" = "red", "upper" = "blue")) +
    theme_minimal()
  
}





# Fit BART models

##########################################################################################
#########............Function to fit BART Models (Competing Event).................#######

FitBARTCompEvent <- function(data, # Dataframe in temporal order: L0, Z,L,M,Y,D
                             var.type, # Vector of variable specifications for data.  Fi=fixed (e.g. the exposure), L=confounder, M= mediator, Y=outcome, D = event indicator. 
                             k,
                             J=10000, # Size of pseudo data. Default is set to 2,000.
                             Ndraws=1000, # Number of posterior draws. Default is set to 200.
                             Nskip = 1000, # Number of burn-in samples. Default is set to 100.
                             Ntree = 100, # Number of trees. Default is set to 100.
                             Keepevery = 10, # Keep every k:th draw. Default is set to 1.
                             Suppress = TRUE, # Indicates if the output should be suppressed. Default is TRUE
                             By = Ndraws/10, # If Suppress is set to FALSE, output is provided for the By:th iteration.
                             Sparse = TRUE, # Indicates if the sparse dirichlet hyper prior should be used. Default is TRUE
                             ...
){
  
  
  
  
  ## test if the variables in data are continuous or binary  
  continuous <- apply(data, 2, function(x) !all(na.omit(x) %in% 0:1)) # Continuous + Multinomial
  
  
  n_Y1 <- length(which(var.type=="Y1")) # Number of outcome variables, >1 allowed
  n_Y2 <- length(which(var.type=="Y2")) # Number of outcome variables, >1 allowed
  n_L0 <- length(which(var.type=="L0")) # Number of baseline confounder variables, >1 allowed
  n_Fi <- length(which(var.type=="Fi")) # Number of fixed variables, >1 allowed
  
  
  if (ncol(data) != length(var.type)) stop("The number of columns in data is not equal to the length of var.type ")
  
  
  
  
  ################################################################
  #Fit the BART models and store results as a list
  BModels <- vector("list", ncol(data) - (n_L0+1))
  for (i in (n_L0+1):(ncol(data)-1)) {
    # Identify complete observations
    id <- ifelse(apply(data[, 1:(1 + i)], 1, function(x) all(!is.na(x))), 1, 0)
    
    if(var.type[1+i] == "Fi") {
      BModels[[i]]  <- NA
    } else {
      if (continuous[1+i] == FALSE & Suppress == TRUE & k==1) {
        quiet(BModels[[i]]  <- BART::pbart(data[id == 1, which(var.type[1:i]!="d" & var.type[1:i]!="Y2")], data[id == 1, (1 + i)], ndpost = Ndraws, nkeeptrain = 0, nkeeptest = 0, ntree = Ntree, nskip = Nskip, sparse = Sparse, keepevery = Keepevery))
      } else if (continuous[1+i] == TRUE & Suppress == TRUE & k==1) {
        quiet(BModels[[i]]  <- BART::wbart(data[id == 1, which(var.type[1:i] != "d" & var.type[1:i]!="Y2")], data[id == 1, (1 + i)], ndpost = Ndraws, nkeeptrain = 0, nkeeptest = 0, ntree = Ntree, nskip = Nskip, sparse = Sparse, keepevery = Keepevery))
      } else if (continuous[1+i] == FALSE & Suppress == FALSE & k==1) {
        BModels[[i]]  <- BART::pbart(data[id == 1, which(var.type[1:i] != "d" & var.type[1:i]!="Y2")], data[id == 1, (1 + i)], ndpost = Ndraws, nkeeptrain = 0, nkeeptest = 0, ntree = Ntree, nskip = Nskip, sparse = Sparse, keepevery = Keepevery)
      } else if (continuous[1+i] == TRUE & Suppress == FALSE & k==1) {
        BModels[[i]]  <- BART::wbart(data[id == 1, which(var.type[1:i] != "d" & var.type[1:i]!="Y2")], data[id == 1, (1 + i)], ndpost = Ndraws, nkeeptrain = 0, nkeeptest = 0, ntree = Ntree, nskip = Nskip, sparse = Sparse, keepevery = Keepevery)
      } else if (continuous[1+i] == FALSE & Suppress == TRUE & k==2) {
        quiet(BModels[[i]]  <- BART::pbart(data[id == 1, which(var.type[1:i]!="d" & var.type[1:i]!="Y1")], data[id == 1, (1 + i)], ndpost = Ndraws, nkeeptrain = 0, nkeeptest = 0, ntree = Ntree, nskip = Nskip, sparse = Sparse, keepevery = Keepevery))
      } else if (continuous[1+i] == TRUE & Suppress == TRUE & k==2) {
        quiet(BModels[[i]]  <- BART::wbart(data[id == 1, which(var.type[1:i] != "d" & var.type[1:i]!="Y1")], data[id == 1, (1 + i)], ndpost = Ndraws, nkeeptrain = 0, nkeeptest = 0, ntree = Ntree, nskip = Nskip, sparse = Sparse, keepevery = Keepevery))
      } else if (continuous[1+i] == FALSE & Suppress == FALSE & k==2) {
        BModels[[i]]  <- BART::pbart(data[id == 1, which(var.type[1:i] != "d" & var.type[1:i]!="Y1")], data[id == 1, (1 + i)], ndpost = Ndraws, nkeeptrain = 0, nkeeptest = 0, ntree = Ntree, nskip = Nskip, sparse = Sparse, keepevery = Keepevery)
      } else if (continuous[1+i] == TRUE & Suppress == FALSE & k==2) {
        BModels[[i]]  <- BART::wbart(data[id == 1, which(var.type[1:i] != "d" & var.type[1:i]!="Y1")], data[id == 1, (1 + i)], ndpost = Ndraws, nkeeptrain = 0, nkeeptest = 0, ntree = Ntree, nskip = Nskip, sparse = Sparse, keepevery = Keepevery)
      }
    }
  }
  
  return(BModels) 
}








# G-computation


##########################################################################################
#########............Function for G-computation (Competing Event).................#######
MCIntegCompEvent <- function(data, # Dataframe in temporal order: L0, Y,Z,L,M,....Y,D
                             var.type, # Vector of variable specifications for data.  Fi=fixed (e.g. the exposure), L=confounder, M= mediator, Y=outcome, D = event indicator. 
                             BModels,
                             k, # event typw: 1 (Main Event) vs. 2 (Competing Event)
                             fixed.regime.trt, # A vector specifying the Exposure  regime Z
                             fixed.regime.control, # A vector specifying the Exposure  regime Z_star
                             conditionOnBaselineAge = FALSE, # Indicates if the output should be conditioned on baseline age. Default is FALSE
                             startBaselineAge, 
                             endBaselineAge = startBaselineAge,
                             J=10000, # Size of pseudo data. Default is set to 2,000.
                             Ndraws=200, # Number of posterior draws. Default is set to 200.
                             #Nskip = 100, # Number of burn-in samples. Default is set to 100.
                             #Ntree = 100, # Number of trees. Default is set to 100.
                             #Keepevery = 1, # Keep every k:th draw. Default is set to 1.
                             Suppress = TRUE, # Indicates if the output should be suppressed. Default is TRUE
                             By = Ndraws/10, # If Suppress is set to FALSE, output is provided for the By:th iteration.
                             Sparse = TRUE, # Indicates if the sparse dirichlet hyper prior should be used. Default is TRUE
                             ...
){
  
  
  
  
  n_Y1 <- length(which(var.type=="Y1")) # Number of outcome variables, >1 allowed
  n_Y2 <- length(which(var.type=="Y2")) # Number of outcome variables, >1 allowed
  n_L0 <- length(which(var.type=="L0")) # Number of baseline confounder variables, >1 allowed
  n_Fi <- length(which(var.type=="Fi")) # Number of fixed variables, >1 allowed
  if(!is.null(fixed.regime.trt)) {
    n_Reg_trt <- max(ifelse(!is.null(fixed.regime.trt), length(fixed.regime.trt), 0))
  } else {
    n_Reg_trt <- 1
  }
  
  if(!is.null(fixed.regime.control)) {
    n_Reg_control <- max(ifelse(!is.null(fixed.regime.control), length(fixed.regime.control), 0))
  } else {
    n_Reg_control <- 1
  }
  
  
  if (ncol(data) != length(var.type)) stop("The number of columns in data is not equal to the length of var.type ")
  
  if(!is.null(fixed.regime.trt)){
    if (length(fixed.regime.trt) != n_Fi) stop("Warning: The number of fixed variables is not equal to the length of the fixed regime(s)")
  }
  
  if(!is.null(fixed.regime.control)){
    if (length(fixed.regime.control) != n_Fi) stop("Warning: The number of fixed variables is not equal to the length of the fixed regime(s)")
  }
  
  
  
  
  
  
  
  ## test if the variables in data are continuous or binary  
  continuous <- apply(data, 2, function(x) !all(na.omit(x) %in% 0:1)) # Continuous + Multinomial
  for (i in (n_L0+1):(ncol(data)-1)) {
    # Identify complete observations
    id <- ifelse(apply(data[, 1:(1 + i)], 1, function(x) all(!is.na(x))), 1, 0)
  }
  
  
  
  
  
  
  
  
  
  
  
  ########################################################
  ## Create matrices to store the output. One for the outcome and one for survival probabilities.
  
  if(k==1){
    Pk_zzstar_mat <- matrix(nrow =  n_Y1 , ncol = Ndraws)
  } else if(k==2){
    Pk_zzstar_mat <- matrix(nrow =  n_Y2 , ncol = Ndraws)
  }
  Pk_zzstar_mat[1,] <- 0                              #assuming everyone is alive at visit 1 
  
  
  
  
  if(conditionOnBaselineAge == FALSE){
    L0data <- data[id==1,c(1:n_L0)] #complete L0 data                         #(numL0Vars+1)
  } else if(conditionOnBaselineAge == TRUE){
    subset_data <- data[data$AGE %in% c(startBaselineAge:endBaselineAge),]
    L0data <- subset_data[id==1,c(1:n_L0)]
  }
  
  
  
  
  
  
  for(it in 1:Ndraws) {
    m <- 2
    l <- 1
    s_hat <- NULL
    # sample data for the baseline confounder using BAYESIAN BOOTSTRAP
    BBweights <- rdirichlet(1,rep(1,nrow(L0data))) #need weights at each iteration
    L0samples <- getSampleByBB(L0data,J,BBweights) #returns a df of baseline confounder samples with J rows
    if(var.type[n_L0+3] == "Fi") {
      x_trt <- data.frame(rep(fixed.regime.trt[1], J))
      x_control <- data.frame(rep(fixed.regime.control[1], J))
      l <- l + 1
    } else {
      x_trt <- data.frame(sample(data[, n_L0+3], size = J, replace = TRUE))
      x_control <- data.frame(sample(data[, n_L0+3], size = J, replace = TRUE))
    }
    x_trt <- cbind(L0samples,x_trt) # Merge baseline confounder samples with the first column of the prediction variable
    x_control <- cbind(L0samples,x_control)
    for (j in (n_L0+4):ncol(data)) {
      if(var.type[j] == "Fi") {
        x_trt <- cbind(x_trt,rep(fixed.regime.trt[l],nrow(x_trt)))
        x_control <- cbind(x_control,rep(fixed.regime.control[l],nrow(x_control)))
        l <- l + 1
      }
      else if(var.type[j] == "L") {
        x_trt <- pred_comb(continuous[j], BModels[[j - 1]], x_trt, it)
        x_control <- pred_comb(continuous[j], BModels[[j - 1]], x_control, it)
      }
      else if(var.type[j] == "M") { 
        x_trt <- pred_comb_mediator(continuous[j], BModels[[j - 1]], x_trt, x_control, it)
        x_control <- pred_comb(continuous[j], BModels[[j - 1]], x_control, it)
      }
      else if(var.type[j] == "Y1" & k==1) {
        Ykprob_trt_it <- pred_surv(BModels[[j - 1]], x_trt, it) #failure probablity
        Ykprob_control_it <- pred_surv(BModels[[j - 1]], x_control, it) #failure probablity
        surv_hat_trt <- 1- as.numeric(Ykprob_trt_it)
        surv_hat_control <- 1- as.numeric(Ykprob_control_it)
        if(!is.null(s_hat)){
          s_temp_trt <- s_hat*surv_hat_trt
          s_temp_control <- s_hat*surv_hat_control
        } else{
          s_temp_trt <- surv_hat_trt
          s_temp_control <- surv_hat_control
        }
        
        Pk_zzstar_mat[m,it] <- mean(Ykprob_trt_it)
        m <- m+1
        surv_trt <- rbinom(length(surv_hat_trt),1,prob = surv_hat_trt)
        
        surv_control <- rbinom(length(surv_hat_trt),1,prob = surv_hat_control)
        x_trt <- cbind(x_trt, surv_trt)
        
        x_trt <- x_trt[surv_trt==1,]
        x_control <- cbind(x_control, surv_control)
        x_control <- x_control[surv_trt==1,]
        s_hat <- s_temp_trt[surv_trt==1]
      }
      else if(var.type[j] == "Y2" & k==2) {
        Ykprob_trt_it <- pred_surv(BModels[[j - 1]], x_trt, it) #failure probablity
        Ykprob_control_it <- pred_surv(BModels[[j - 1]], x_control, it) #failure probablity
        surv_hat_trt <- 1- as.numeric(Ykprob_trt_it)
        surv_hat_control <- 1- as.numeric(Ykprob_control_it)
        if(!is.null(s_hat)){
          s_temp_trt <- s_hat*surv_hat_trt
          s_temp_control <- s_hat*surv_hat_control
        } else{
          s_temp_trt <- surv_hat_trt
          s_temp_control <- surv_hat_control
        }
        Pk_zzstar_mat[m,it] <- mean(Ykprob_trt_it)
        m <- m+1
        surv_trt <- rbinom(length(surv_hat_trt),1,prob = surv_hat_trt)
        surv_control <- rbinom(length(surv_hat_trt),1,prob = surv_hat_control)
        x_trt <- cbind(x_trt, surv_trt)
        x_trt <- x_trt[surv_trt==1,]
        x_control <- cbind(x_control, surv_control)
        x_control <- x_control[surv_trt==1,]
        s_hat <- s_temp_trt[surv_trt==1]
      }
    }
  }
  return(Pk_zzstar_mat)
}







