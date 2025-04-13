#setwd("/blue/daniels/s.bhandari/ARIC/Code here")
#install.packages('GcompBART2_0.1.0.tar.gz', repos=NULL, type='source')
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
library(gtools)





ARIC_DATA = read_sas('ARIC.sas7bdat')
ARIC_df = data.frame(ARIC_DATA)
head(ARIC_df)


ARIC_workdf = data.frame(ID_FINAL = ARIC_df$ID_FINAL,  EXAM = ARIC_df$EXAM,RACE = ARIC_df$RACE, BMI = ARIC_df$BMINOW, 
                         AGE = ARIC_df$AGE,SEX = ARIC_df$SEX,EDU = ARIC_df$edu_g, RXHYP = ARIC_df$RXHYP, 
                         SMOKE_STATUS_BIN =ARIC_df$SMOKER, SMOKE_STATUS_MULT =ARIC_df$SMOKE_STATUS, SBP = ARIC_df$SBP,
                         DBP = ARIC_df$DBP, FOLLOWUPYRS = ARIC_df$LENYFL, CVD_DTH = ARIC_df$CVD_DTH, 
                         TOT_DTH =ARIC_df$TOT_DTH, GLUCOSE = ARIC_df$GLUCOSE, HXDIAB = ARIC_df$HXDIAB, 
                         TOT_HDL_CHL_RATIO = ARIC_df$TOTCHL/ARIC_df$HDLCHL)

head(ARIC_workdf)

# Define Mean BP
MeanBP = (ARIC_workdf$SBP + ARIC_workdf$DBP)/2
ARIC_workdf$MeanBP = MeanBP



#head(ARIC_workdf)
setDT(ARIC_workdf)
ARIC_workdf_wide = dcast(ARIC_workdf, ID_FINAL + SEX+ EDU +RACE   ~ EXAM, 
                         value.var = c("AGE", "BMI", "FOLLOWUPYRS", "RXHYP","SMOKE_STATUS_BIN","SMOKE_STATUS_MULT","MeanBP",
                                       "CVD_DTH","TOT_DTH", "GLUCOSE", "HXDIAB","TOT_HDL_CHL_RATIO"))


# Define DIAB at baseline
diab_baseline = ifelse(ARIC_workdf_wide$GLUCOSE_1 >= 126 | ARIC_workdf_wide$HXDIAB_1 ==1, 1,0) 

ARIC_workdf_wide$diab_baseline = diab_baseline
head(ARIC_workdf_wide)





#..............................REMOVE Non-MONOTONE MISSINGNESS.......................//


#We create three dataframes (one each for Exposure, Confounder, and Mediator) where a cell is TRUE if the actual observation is NA and FALSE otherwise.
Exp_NAdf = data.frame(ID_FINAL = ARIC_workdf_wide$ID_FINAL, ISZ1NA = is.na(ARIC_workdf_wide$RXHYP_1),ISZ2NA = is.na(ARIC_workdf_wide$RXHYP_2)
                      ,ISZ3NA =is.na(ARIC_workdf_wide$RXHYP_3),ISZ4NA =is.na(ARIC_workdf_wide$RXHYP_4))
Med_NAdf = data.frame(ID_FINAL = ARIC_workdf_wide$ID_FINAL, ISM1NA = is.na(ARIC_workdf_wide$MeanBP_1),ISM2NA = is.na(ARIC_workdf_wide$MeanBP_2)
                      ,ISM3NA =is.na(ARIC_workdf_wide$MeanBP_3),ISM4NA =is.na(ARIC_workdf_wide$MeanBP_4))
Conf_NAdf = data.frame(ID_FINAL = ARIC_workdf_wide$ID_FINAL, ISL1NA = is.na(ARIC_workdf_wide$SMOKE_STATUS_BIN_1), ISL2NA = is.na(ARIC_workdf_wide$SMOKE_STATUS_BIN_2)
                       ,ISL3NA = is.na(ARIC_workdf_wide$SMOKE_STATUS_BIN_3), ISL4NA = is.na(ARIC_workdf_wide$SMOKE_STATUS_BIN_4))





#input df must be in format [ID, T/F, ..., T/F] where T = NA, F = No NA
#returns all IDs with monotone data
get_monotone_data_ID = function(patternDF){
  pat1 = cbind(FALSE,FALSE,FALSE,FALSE)
  pat2 = cbind(FALSE,FALSE,FALSE,TRUE)
  pat3 = cbind(FALSE,FALSE,TRUE,TRUE)
  pat4 = cbind(FALSE,TRUE,TRUE,TRUE)
  #only keep IDs that match these patterns
  keepId = apply(patternDF[,2:5],1,function(x){all(x==pat1)|all(x==pat2)|all(x==pat3)|all(x==pat4)}) 
  patternDF$keepId = keepId #add a new column in the input dataframe (TRUE = we wish to keep this ID, FALSE = discard)
  returnID = patternDF[patternDF$keepId == TRUE,][,1]
  return(returnID)
}



get_monotone_Exp_ID = get_monotone_data_ID(Exp_NAdf)
monotone_Exp_dataframe = Exp_NAdf[Exp_NAdf$ID_FINAL %in% get_monotone_Exp_ID, ]

get_monotone_Conf_ID = get_monotone_data_ID(Conf_NAdf)
monotone_Conf_dataframe = Conf_NAdf[Conf_NAdf$ID_FINAL %in% get_monotone_Conf_ID, ]

get_monotone_Med_ID = get_monotone_data_ID(Med_NAdf)
monotone_Med_dataframe = Med_NAdf[Med_NAdf$ID_FINAL %in% get_monotone_Med_ID, ]


monotoneDataID = intersect(get_monotone_Exp_ID,get_monotone_Conf_ID)
monotoneDataID = intersect(monotoneDataID,get_monotone_Med_ID )

ARIC_workdf_wide = ARIC_workdf_wide[ARIC_workdf_wide$ID_FINAL %in% monotoneDataID]
head(ARIC_workdf_wide)







#................................CREATE OUTCOME VARIABLES............................//
AllDeathIndicator = ARIC_workdf_wide$TOT_DTH_1
#delta indicates CVD death vs censoring
delta = ARIC_workdf_wide$CVD_DTH_1

# d indicates death by CVD vs death by other causes
d = ifelse(AllDeathIndicator == 1, 2,0) # first set all death types = 2 and censoring to 0
d = ifelse(delta ==1, delta, d) # Next set only CVD deaths (main event) to 1 keeping other deaths (competing event) to 2 and censoring to 0




ageAtDeathOrCensor = ARIC_workdf_wide$AGE_1 + ARIC_workdf_wide$FOLLOWUPYRS_1 #baseline age + follow-up years remaining at visit 1
ARIC_workdf_wide$ageAtDeathOrCensor = ageAtDeathOrCensor

Y1 = rep(0,length(ARIC_workdf_wide$ID_FINAL))   #Assumption that everyone is alive at the first visit

Y2 = ifelse(!(is.na(ARIC_workdf_wide$AGE_1)) & is.na(ARIC_workdf_wide$AGE_2) & is.na(ARIC_workdf_wide$AGE_3) & is.na(ARIC_workdf_wide$AGE_4) & delta ==1 , 1,0) 
AGE_2 = replace(ARIC_workdf_wide$AGE_2, is.na(ARIC_workdf_wide$AGE_2), 0)
Y2 = ifelse((ARIC_workdf_wide$ageAtDeathOrCensor > ARIC_workdf_wide$AGE_1) &
              (ARIC_workdf_wide$ageAtDeathOrCensor <= AGE_2)  & 
              delta ==1,
            1,Y2) 


Y3 = ifelse(!(is.na(ARIC_workdf_wide$AGE_2)) & is.na(ARIC_workdf_wide$AGE_3) & is.na(ARIC_workdf_wide$AGE_4)  & delta ==1 , 1,0) 
AGE_3 = replace(ARIC_workdf_wide$AGE_3, is.na(ARIC_workdf_wide$AGE_3), 0)
Y3 = ifelse((ARIC_workdf_wide$ageAtDeathOrCensor > ARIC_workdf_wide$AGE_2) &
              (ARIC_workdf_wide$ageAtDeathOrCensor <= AGE_3)  & 
              delta ==1,
            1,Y3) 



Y4 = ifelse(!(is.na(ARIC_workdf_wide$AGE_3)) & is.na(ARIC_workdf_wide$AGE_4)   & delta ==1 , 1,0) 
AGE_4 = replace(ARIC_workdf_wide$AGE_4, is.na(ARIC_workdf_wide$AGE_4), 0)
Y4 = ifelse((ARIC_workdf_wide$ageAtDeathOrCensor > ARIC_workdf_wide$AGE_3) &
              (ARIC_workdf_wide$ageAtDeathOrCensor <= AGE_4)  & 
              delta ==1,
            1,Y4) 


# Y5 = ifelse(!(is.na(ARIC_workdf_wide$AGE_4)) & is.na(ARIC_workdf_wide$AGE_5)   & delta ==1 , 1,0) 
# AGE_5 = replace(ARIC_workdf_wide$AGE_5, is.na(ARIC_workdf_wide$AGE_5), 0)
# Y5 = ifelse((ARIC_workdf_wide$ageAtDeathOrCensor > ARIC_workdf_wide$AGE_4) &
#               (ARIC_workdf_wide$ageAtDeathOrCensor <= AGE_5)  & 
#               delta ==1,
#             1,Y5) 



Y5 = ifelse(!(is.na(ARIC_workdf_wide$AGE_4)) & is.na(ARIC_workdf_wide$AGE_5)   & delta ==1 , 1,0) 
AGE_5 = replace(ARIC_workdf_wide$AGE_5, is.na(ARIC_workdf_wide$AGE_5), 0)
Y5 = ifelse((ARIC_workdf_wide$ageAtDeathOrCensor > ARIC_workdf_wide$AGE_4) &
              delta ==1,
            1,Y5) 






ARIC_workdf_wide$Y1 = Y1
ARIC_workdf_wide$Y2 = Y2
ARIC_workdf_wide$Y3 = Y3
ARIC_workdf_wide$Y4 = Y4
ARIC_workdf_wide$delta = delta
ARIC_workdf_wide$d = d
#ARIC_workdf_wide$AllDeathIndicator = AllDeathIndicator
#head(ARIC_workdf_wide)



### Individuals who are hypertensive (at baseline) in ARIC

# Create a new variable hyp_at_baseline
ARIC_workdf <- ARIC_workdf %>%
  group_by(ID_FINAL) %>%
  mutate(hyp_at_baseline = ifelse(row_number() == 1 & (SBP >= 140 | DBP >= 90), 1, 0))
# Extract IDs with hyp_at_baseline == 1


ID_ARIC_hyp_at_baseline <- ARIC_workdf %>%
  filter(hyp_at_baseline == 1) %>%
  select(ID_FINAL) %>%
  unique()
ID_ARIC_hyp_at_baseline = as.vector(unlist(ID_ARIC_hyp_at_baseline))



# Update ARIC
ARIC_workdf = ARIC_workdf[ARIC_workdf$ID_FINAL %in% ID_ARIC_hyp_at_baseline,]
ARIC_workdf_wide = ARIC_workdf_wide[ARIC_workdf_wide$ID_FINAL %in% ID_ARIC_hyp_at_baseline,]


#######################################################################################################
#######################################################################################################


############################################################
####....... INPUT DATA FOR PARAMETRIC REGRESSION.....#########




param_reg_df_Y2 = data.frame(V0 = ARIC_workdf_wide$SMOKE_STATUS_MULT_1, BMI = ARIC_workdf_wide$BMI_1 ,
                             AGE = ARIC_workdf_wide$AGE_1, RACE = ARIC_workdf_wide$RACE, 
                             SEX = ARIC_workdf_wide$SEX, EDU = ARIC_workdf_wide$EDU,
                             Z_1 = ARIC_workdf_wide$RXHYP_1, 
                             L_1 = ARIC_workdf_wide$SMOKE_STATUS_BIN_1, 
                             M_1 = ARIC_workdf_wide$MeanBP_1,  
                             Y_2 = ARIC_workdf_wide$Y2)

param_reg_df_Y2_input_list = list(param_reg_df_Y2)

param_reg_df_Y3 = data.frame(V0 = ARIC_workdf_wide$SMOKE_STATUS_MULT_1, BMI = ARIC_workdf_wide$BMI_1 ,
                             AGE = ARIC_workdf_wide$AGE_1, RACE = ARIC_workdf_wide$RACE, 
                             SEX = ARIC_workdf_wide$SEX, EDU = ARIC_workdf_wide$EDU, 
                             Z_1 = ARIC_workdf_wide$RXHYP_1, 
                             L_1 = ARIC_workdf_wide$SMOKE_STATUS_BIN_1, 
                             M_1 = ARIC_workdf_wide$MeanBP_1,
                             Z_2 = ARIC_workdf_wide$RXHYP_2, 
                             L_2 = ARIC_workdf_wide$SMOKE_STATUS_BIN_2,
                             M_2 = ARIC_workdf_wide$MeanBP_2, 
                             Y_3 = ARIC_workdf_wide$Y3)

param_reg_df_Y3_input_list = list(param_reg_df_Y3)

param_reg_df_ZLM= data.frame(V0 = ARIC_workdf_wide$SMOKE_STATUS_MULT_1, BMI = ARIC_workdf_wide$BMI_1 ,
                             AGE = ARIC_workdf_wide$AGE_1, RACE = ARIC_workdf_wide$RACE, 
                             SEX = ARIC_workdf_wide$SEX, EDU = ARIC_workdf_wide$EDU, 
                             Z_1 = ARIC_workdf_wide$RXHYP_1, 
                             L_1 = ARIC_workdf_wide$SMOKE_STATUS_BIN_1, 
                             M_1 = ARIC_workdf_wide$MeanBP_1,  
                             Z_2 = ARIC_workdf_wide$RXHYP_2, 
                             L_2 = ARIC_workdf_wide$SMOKE_STATUS_BIN_2,
                             M_2 = ARIC_workdf_wide$MeanBP_2)


param_reg_df_ZLM_input_list = list(param_reg_df_ZLM)   




## Ground truth (parameters) using parametric linear models


is_binary = function(x) {
  unique_values = unique(na.omit(x))
  
  # Check if the vector is binary
  if (all(unique_values %in% c(0, 1))) {
    return(TRUE)
  } else  {
    return(FALSE)
  }
}

# This function fits the model and returns entire model outcome 
fit_parametric_glm_models = function(covariates_df, # Dataframe in temporal order: L0,Y, Z,L,M,Y,D
                                     outcome_vec
){
  
  # Combine the outcome and covariates into a single dataframe
  data = data.frame(outcome = outcome_vec, covariates_df)
  # Remove rows with any NA values
  data = na.omit(data)
  
  # Determine if the outcome is binary or continuous
  outcome_type_is_binary = is_binary(outcome_vec)
  # Create a formula for the model
  formula = as.formula(paste("outcome ~", paste(names(covariates_df), collapse = " + ")))
  
  
  # Fit the appropriate linear/logistic model
  if (outcome_type_is_binary == TRUE) {                        
    model_fit = glm(formula, data = data, family = binomial(link="logit"))
    #return(as.numeric(predict(model, type = "response")))
  } else  {
    model_fit = lm(formula, data = data)
    #return_list = vector("list", 2)
    #return_list[[1]] = as.numeric(predict(model, type = "response"))
    #return_list[[2]] = sigma(model)
    #return(return_list)
  } 
  return(model_fit)
}



get_parametric_coeffs_iter = function(param_reg_df_Y2,   #regression data for fitting Y2 model
                                      param_reg_df_Y3,   #regression data for fitting Y2 model
                                      param_reg_df_ZLM #regression data for fitting Z/L/M model
){
  
  
  Y2_cov_df = param_reg_df_Y2 %>% select(1:(which(names(param_reg_df_Y2) == "Y_2") - 1))
  fit_glm_Y2 = fit_parametric_glm_models(Y2_cov_df, param_reg_df_Y2$Y_2)
  if(is_binary(param_reg_df_Y2$Y_2) == TRUE){
    fit_glm_Y2_coeff = as.data.frame(coef(fit_glm_Y2))
  }else{
    fit_glm_Y2_coeff = as.data.frame(coef(fit_glm_Y2))
    fit_glm_Y2_coeff = rbind(fit_glm_Y2_coeff, sigma = sigma(fit_glm_Y2))
  }
  
  
  Y3_cov_df = param_reg_df_Y3 %>% select(1:(which(names(param_reg_df_Y3) == "Y_3") - 1))
  fit_glm_Y3 = fit_parametric_glm_models(Y3_cov_df, param_reg_df_Y3$Y_3)
  if(is_binary(param_reg_df_Y3$Y_3) == TRUE){
    fit_glm_Y3_coeff = as.data.frame(coef(fit_glm_Y3))
  }else{
    fit_glm_Y3_coeff = as.data.frame(coef(fit_glm_Y3))
    fit_glm_Y3_coeff = rbind(fit_glm_Y3_coeff, sigma = sigma(fit_glm_Y3))
  }
  
  
  
  Z1_cov_df = param_reg_df_ZLM %>% select(1:(which(names(param_reg_df_ZLM) == "Z_1") - 1))
  fit_glm_Z1 = fit_parametric_glm_models(Z1_cov_df, param_reg_df_ZLM$Z_1)
  if(is_binary(param_reg_df_ZLM$Z_1) == TRUE){
    fit_glm_Z1_coeff = as.data.frame(coef(fit_glm_Z1))
  }else{
    fit_glm_Z1_coeff = as.data.frame(coef(fit_glm_Z1))
    fit_glm_Z1_coeff = rbind(fit_glm_Z1_coeff, sigma = sigma(fit_glm_Z1))
  }
  
  
  Z2_cov_df = param_reg_df_ZLM %>% select(1:(which(names(param_reg_df_ZLM) == "Z_2") - 1))
  fit_glm_Z2 = fit_parametric_glm_models(Z2_cov_df, param_reg_df_ZLM$Z_2)
  if(is_binary(param_reg_df_ZLM$Z_2) == TRUE){
    fit_glm_Z2_coeff = as.data.frame(coef(fit_glm_Z2))
  }else{
    fit_glm_Z2_coeff = as.data.frame(coef(fit_glm_Z2))
    fit_glm_Z2_coeff = rbind(fit_glm_Z2_coeff, sigma = sigma(fit_glm_Z2))
  }
  
  
  
  
  
  L1_cov_df = param_reg_df_ZLM %>% select(1:(which(names(param_reg_df_ZLM) == "L_1") - 1))
  fit_glm_L1 = fit_parametric_glm_models(L1_cov_df, param_reg_df_ZLM$L_1)
  if(is_binary(param_reg_df_ZLM$L_1) == TRUE){
    fit_glm_L1_coeff = as.data.frame(coef(fit_glm_L1))
  }else{
    fit_glm_L1_coeff = as.data.frame(coef(fit_glm_L1))
    fit_glm_L1_coeff = rbind(fit_glm_L1_coeff, sigma = sigma(fit_glm_L1))
  }
  
  
  L2_cov_df = param_reg_df_ZLM %>% select(1:(which(names(param_reg_df_ZLM) == "L_2") - 1))
  fit_glm_L2 = fit_parametric_glm_models(L2_cov_df, param_reg_df_ZLM$L_2)
  if(is_binary(param_reg_df_ZLM$L_2) == TRUE){
    fit_glm_L2_coeff = as.data.frame(coef(fit_glm_L2))
  }else{
    fit_glm_L2_coeff = as.data.frame(coef(fit_glm_L2))
    fit_glm_L2_coeff = rbind(fit_glm_L2_coeff, sigma = sigma(fit_glm_L2))
  }
  
  
  
  
  
  M1_cov_df = param_reg_df_ZLM %>% select(1:(which(names(param_reg_df_ZLM) == "M_1") - 1))
  fit_glm_M1 = fit_parametric_glm_models(M1_cov_df, param_reg_df_ZLM$M_1)
  if(is_binary(param_reg_df_ZLM$M_1) == TRUE){
    fit_glm_M1_coeff = as.data.frame(coef(fit_glm_M1))
  }else{
    fit_glm_M1_coeff = as.data.frame(coef(fit_glm_M1))
    fit_glm_M1_coeff = rbind(fit_glm_M1_coeff, sigma = sigma(fit_glm_M1))
  }
  
  
  M2_cov_df = param_reg_df_ZLM %>% select(1:(which(names(param_reg_df_ZLM) == "M_2") - 1))
  fit_glm_M2 = fit_parametric_glm_models(M2_cov_df, param_reg_df_ZLM$M_2)
  if(is_binary(param_reg_df_ZLM$M_2) == TRUE){
    fit_glm_M2_coeff = as.data.frame(coef(fit_glm_M2))
  }else{
    fit_glm_M2_coeff = as.data.frame(coef(fit_glm_M2))
    fit_glm_M2_coeff = rbind(fit_glm_M2_coeff, sigma = sigma(fit_glm_M2))
  }
  
  
  
  
  param_list_iter =  list(param_L1 = fit_glm_L1_coeff,
                          param_M1 = fit_glm_M1_coeff, 
                          param_Y2 = fit_glm_Y2_coeff,
                          param_L2 = fit_glm_L2_coeff,
                          param_M2 = fit_glm_M2_coeff, 
                          param_Y3 = fit_glm_Y3_coeff)
  
  
  
  return(param_list_iter)
  
  
}




get_parametric_coeffs = function(param_reg_list_Y2,   #reg data list (of dfs) for fitting Y2 model
                                 param_reg_list_Y3,   #reg data list (of dfs) for fitting Y2 model
                                 param_reg_list_ZLM,  #reg data list (of dfs) for fitting Z/L/M model
                                 num_iters            # number of MCMC iterations Q
){
  param_list = replicate(num_iters, list(), simplify = FALSE)    # this will be a list of Q lists
  
  for(it in 1:num_iters){
    param_list[[it]] = get_parametric_coeffs_iter(param_reg_list_Y2[[it]],  
                                                  param_reg_list_Y3[[it]],
                                                  param_reg_list_ZLM[[it]]
    )
    
  }
  return(param_list)
  
  
}





ground_truth_parametric_coeffs_list = get_parametric_coeffs(param_reg_list_Y2 =param_reg_df_Y2_input_list,
                                                            param_reg_list_Y3 =param_reg_df_Y3_input_list,
                                                            param_reg_list_ZLM = param_reg_df_ZLM_input_list,
                                                            num_iters= 1)

ground_truth_parametric_coeffs_list



#######################################################################################################
#######################################################################################################

## G-Computation

source("pred_comb_mediator.R")
source("pred_surv.R")
source("pwbart_it.R")
source("expit.R")
source("quiet.R")
source("pred_comb_Y.R")
source("pred_comb_mediator.R")
source("pred_comb.R")
source("pred_y.R")


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





library(ggplot2)

calculate_95confidence_interval <- function(data) {
  cat("Posterior Means: \n")
  conf_int <- function(data_row){
    n <- length(data_row)
    mean_value <- as.numeric(format(round(mean(data_row), 4), nsmall = 4))
    print(mean_value)
    std_error <- sd(data_row) / sqrt(n)
    z_value <- qnorm(0.975)  # Z-value for a 95% confidence interval
    lower <- mean_value - z_value * std_error
    lower <- as.numeric(format(round(lower, 4), nsmall = 4))
    upper <- mean_value + z_value * std_error
    upper <- as.numeric(format(round(upper, 4), nsmall = 4))
    interval <- cbind(lower, upper)
    return(interval)
  }
  return_mat = apply(data,1,conf_int)
  return(t(return_mat))
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


library(tidyr)

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

#######################################################################################################
#######################################################################################################

### MC Integration


pred_comb <- function(cont,                           # TRUE if the given var is continuous
                      param_model_coeff,              # this is a df that containes coeff (df is an element of a list)
                      MCdata                          # MC (synthetic) data upto this point
){ 
  param_model_coeff = unlist(param_model_coeff); 
  sd_hat = as.numeric(param_model_coeff[length(param_model_coeff)])
  
  MCdata = as.matrix(MCdata)
  param_model_coeff = matrix(unlist(param_model_coeff), nrow = length(unlist(param_model_coeff)), ncol =1)
  if(cont == FALSE) {
    #cat("reached hereL2\n")
    #linear_predictor = as.numeric(MCdata %*%param_model_coeff)
    #x_hat = exp(linear_predictor) / (1 + exp(linear_predictor))
    x_hat <- plogis(as.numeric(MCdata %*%param_model_coeff))
    new_MCdata <- cbind(MCdata, rbinom(length(x_hat), 1, prob = x_hat))
    #print(new_MCdata)
  } else if(cont == TRUE) {
    param_model_coeff <- param_model_coeff[-nrow(param_model_coeff), , drop = FALSE]
    x_hat <- as.numeric(MCdata %*%param_model_coeff)
    new_MCdata <- cbind(MCdata, rnorm(length(x_hat), 
                                      mean = x_hat, sd = sd_hat))
  }
  return(new_MCdata)
}







pred_surv <- function(param_model_coeff, MCdata){
  #intercept = rep(1, nrow(MCdata))
  MCdata = as.matrix(MCdata)
  #print(head(MCdata))
  param_model_coeff = matrix(unlist(param_model_coeff), nrow = length(unlist(param_model_coeff)), ncol =1)
  #linear_predictor = as.numeric(MCdata %*%param_model_coeff)
  #s_hat = exp(linear_predictor) / (1 + exp(linear_predictor))
  s_hat <- plogis(as.numeric(MCdata %*%param_model_coeff))
  #print(s_hat)
  return(s_hat)
}


pred_comb_mediator <- function(cont, param_model_coeff, MCdata_trt, MCdata_control){
  
  MCdata_control = as.matrix(MCdata_control)
  if(cont == FALSE) {
    param_model_coeff = matrix(unlist(param_model_coeff), nrow = length(unlist(param_model_coeff)), ncol =1)
    x_hat <- plogis(as.numeric(MCdata_control %*%param_model_coeff)) #make predictions conditional on counterfactual data
    new_MCdata <- cbind(MCdata_trt, rbinom(length(x_hat), 1, prob = x_hat))
  } else if(cont == TRUE) {
    param_model_coeff = unlist(param_model_coeff); 
    sd_hat = as.numeric(param_model_coeff[length(param_model_coeff)])
    param_model_coeff = param_model_coeff[1:(length(param_model_coeff) - 1)]
    param_model_coeff = matrix(param_model_coeff, nrow = length(param_model_coeff), ncol =1)
    #cat("Print param_model_coeff:\n")
    #print(param_model_coeff)
    x_hat = as.numeric(MCdata_control %*%param_model_coeff) 
    #cat("Print mean:\n")
    #print(x_hat)
    #cat("Print sd:\n")
    #print(sd_hat)
    new_predicted_M = rnorm(length(x_hat), 
                            mean = x_hat, sd = sd_hat)
    #print(new_predicted_M)
    new_MCdata = cbind(MCdata_trt,new_predicted_M )
  }
  #print(new_MCdata)
  return(new_MCdata)
}







########################################################################################
######.............Function for MC integration (No Competing Event).............########
MCInteg_param_model <- function(data, # Dataframe in temporal order: L0,Y, Z,L,M,...,Y,D
                                var.type, # Vector of variable specifications for data.  Fi=fixed (e.g. the exposure), L=confounder, M= mediator, Y=outcome, D = event indicator. 
                                #var.type.temp.order, #var.type indicating temporal oder to longitudinal variables: L0,Y1, Z1,L1,M1,...,YJ,D
                                PModels, #fitted parametric model: list of coefficients for all covariates (also sig_sq for cont variables)
                                fixed.regime.trt, # A vector specifying the Exposure  regime Z
                                fixed.regime.control, # A vector specifying the Exposure  regime Z_star
                                conditionOnBaselineAge = FALSE, # Indicates if the output should be conditioned on baseline age. Default is FALSE
                                startBaselineAge, 
                                endBaselineAge = startBaselineAge,
                                J=2000, # Size of pseudo data. Default is set to 2,000.
                                Ndraws=200, # Number of posterior draws. Default is set to 200.
                                #Nskip = 100, # Number of burn-in samples. Default is set to 100.
                                #Ntree = 100, # Number of trees. Default is set to 100.
                                #Keepevery = 1, # Keep every k:th draw. Default is set to 1.
                                Suppress = TRUE, # Indicates if the output should be suppressed. Default is TRUE
                                By = Ndraws/10, # If Suppress is set to FALSE, output is provided for the By:th iteration.
                                Sparse = TRUE, # Indicates if the sparse dirichlet hyper prior should be used. Default is TRUE
                                ...
){
  
  n_Y <- length(which(var.type=="Y")) # Number of outcome variables, >1 allowed
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
  for (i in (n_L0+1):(ncol(data)-1)) {  #(numL0Vars+1)
    # Identify complete observations
    id <- ifelse(apply(data[, 1:(1 + i)], 1, function(x) all(!is.na(x))), 1, 0)
  }
  
  
  
  
  
  
  
  
  ########################################################
  ## Create matrices to store the output. 
  
  
  P_zzstar_mat <- matrix(nrow =  n_Y+1 , ncol = Ndraws)  #nrow =  n_Y+1 because we removed Y1 from input data but still have analysis until Y3
  P_zzstar_mat[1,] <- 0                              #assuming everyone is alive at visit 1
  
  
  
  
  if(conditionOnBaselineAge == FALSE){
    L0data <- data[id==1,c(1:n_L0)] #complete L0 data                         #(numL0Vars+1)
  } else if(conditionOnBaselineAge == TRUE){
    subset_data <- data[data$AGE %in% c(startBaselineAge:endBaselineAge),]
    L0data <- subset_data[id==1,c(1:n_L0)]
  }
  for(it in 1:Ndraws) {
    index_coeff_list <- 1
    m <- 2                               # m indexes rows of the outcome matrix P_zzstar
    l <- 1                               # l indexes FIXED longitudinal trt regime in G-comp
    s_hat <- NULL
    # sample data for the baseline confounder using BAYESIAN BOOTSTRAP
    BBweights <- rdirichlet(1,rep(1,nrow(L0data))) #need weights at each iteration
    L0samples <- getSampleByBB(L0data,J,BBweights) #returns a df of baseline confounder samples with J rows
    if(var.type[n_L0+1] == "Fi") {
      x_trt <- data.frame(rep(fixed.regime.trt[1], J))
      x_control <- data.frame(rep(fixed.regime.control[1], J))
      l <- l + 1
    } else {
      x_trt <- data.frame(sample(data[, n_L0+1], size = J, replace = TRUE))
      x_control <- data.frame(sample(data[, n_L0+1], size = J, replace = TRUE))
    }
    x_trt <- cbind(L0samples,x_trt) # Merge baseline confounder samples with the first column of the prediction variable
    x_control <- cbind(L0samples,x_control)
    
    intercept = rep(1, J)
    x_trt <- cbind(intercept,x_trt)        # Merge intercept for parametric regression prediction as the first col
    x_control <- cbind(intercept,x_control)
    for (j in (n_L0+2):ncol(data)) {                                         #(numL0Vars+1), j indexes longitudinal data
      if(var.type[j] == "Fi") {
        x_trt <- cbind(x_trt,rep(fixed.regime.trt[l],nrow(x_trt)))
        x_control <- cbind(x_control,rep(fixed.regime.control[l],nrow(x_control)))
        l <- l + 1
      }
      else if(var.type[j] == "L") {
        #print(PModels[[index_coeff_list]])
        x_trt <- pred_comb(continuous[j], PModels[[it]][index_coeff_list], x_trt)
        x_control <- pred_comb(continuous[j], PModels[[it]][index_coeff_list], x_control)
        index_coeff_list<- index_coeff_list+1
        #cat("reached L end.\n")
      }
      else if(var.type[j] == "M") { 
        #print(PModels[[index_coeff_list]])
        x_trt <- pred_comb_mediator(continuous[j], PModels[[it]][index_coeff_list], x_trt, x_control)
        #print(x_trt)
        x_control <- pred_comb(continuous[j], PModels[[it]][index_coeff_list], x_control)
        #print(x_control)
        index_coeff_list<- index_coeff_list+1
        #cat("reached M end.\n")
      }
      # else if(var.type[j] == "Y") {
      # Yprob_trt_it <- lapply(x_trt, function(x2) {pred_surv(BModels[[j - 1]], data.frame(x2), it)})
      # P_zzstar_mat[m,it] <- mean(unlist(Yprob_trt_it))
      # m <- m+1
      # x_trt <- lapply(x_trt, function(x2) {pred_comb(continuous[j], BModels[[j - 1]], data.frame(x2), it)})
      # x_control <- lapply(x_control, function(x2) {pred_comb(continuous[j], BModels[[j - 1]], data.frame(x2), it)})
      # }
      else if(var.type[j] == "Y") {
        #print(PModels[[index_coeff_list]])
        Yprob_trt_it <- pred_surv(PModels[[it]][index_coeff_list], x_trt) #failure probablity
        Yprob_control_it <- pred_surv(PModels[[it]][index_coeff_list], x_control) #failure probablity
        index_coeff_list<- index_coeff_list+1
        
        surv_hat_trt <- 1- as.numeric(Yprob_trt_it)
        surv_hat_control <- 1- as.numeric(Yprob_control_it)
        if(!is.null(s_hat)){
          s_temp_trt <- s_hat*surv_hat_trt
          s_temp_control <- s_hat*surv_hat_control
        } else{
          s_temp_trt <- surv_hat_trt
          s_temp_control <- surv_hat_control
        }
        #cat("The Monte Carlo Standard Error","for iteration", it ,"is:", sd(Yprob_trt_it)/sqrt(length(Yprob_trt_it)), "\n")
        #cat("The ratio  SE/mean(estimates) is:", (sd(Yprob_trt_it)/sqrt(length(Yprob_trt_it)))/mean(Yprob_trt_it), "\n")
        P_zzstar_mat[m,it] <- mean(Yprob_trt_it)
        m <- m+1
        surv_trt <- rbinom(length(surv_hat_trt),1,prob = surv_hat_trt)
        #cat("Out of", J, "MC samples, the number of people alive at visit",m-1, "of iteration", it,":", sum(surv_trt), "\n")
        surv_control <- rbinom(length(surv_hat_trt),1,prob = surv_hat_control)
        #x_trt <- cbind(x_trt, surv_trt)     # DONOT cbind here for parametric regression
        x_trt <- x_trt[surv_trt==1,]
        #x_control <- cbind(x_control, surv_control)      # DONOT cbind here for parametric regression
        x_control <- x_control[surv_trt==1,]
        s_hat <- s_temp_trt[surv_trt==1]
        #cat("reached Y end.\n")
      }
    }
  }
  
  return(P_zzstar_mat)
}




#######################################################################################################
#######################################################################################################


## Function Calls


############################################################
####....... FUNCTION CALL (NO COMPETING EVENT).....#########

input_data1 = data.frame(V0 = ARIC_workdf_wide$SMOKE_STATUS_MULT_1, BMI = ARIC_workdf_wide$BMI_1 ,
                         AGE = ARIC_workdf_wide$AGE_1, RACE = ARIC_workdf_wide$RACE, 
                         SEX = ARIC_workdf_wide$SEX, EDU = ARIC_workdf_wide$EDU, 
                         Z1 = ARIC_workdf_wide$RXHYP_1, 
                         L1 = ARIC_workdf_wide$SMOKE_STATUS_BIN_1, 
                         M1 = ARIC_workdf_wide$MeanBP_1,  
                         Y2 = ARIC_workdf_wide$Y2,
                         Z2 = ARIC_workdf_wide$RXHYP_2, 
                         L2 = ARIC_workdf_wide$SMOKE_STATUS_BIN_2,
                         M2 = ARIC_workdf_wide$MeanBP_2, 
                         Y3 = ARIC_workdf_wide$Y3, 
                         delta = ARIC_workdf_wide$delta)


input_varType1 = c("L0","L0","L0","L0","L0","L0","Fi", "L", "M","Y", "Fi", "L", "M", "Y","D")
#input_varType2 = c("L0","L0","L0","L0","L0","L0","L0","L0","Y1","Fi1", "L1", "M1","Y2", "Fi2", "L2", "M2", "Y3","D")







#################################################################################
#...................NO CONDITION ON BASELINE AGE................................#
P_zzstar = MCInteg_param_model(data = input_data1,
                               var.type = input_varType1,
                               #var.type.temp.order = input_varType2,
                               PModels = ground_truth_parametric_coeffs_list,
                               fixed.regime.trt = rep(1,2),
                               fixed.regime.control = rep(0,2),
                               #conditionOnBaselineAge = TRUE, 
                               #startBaselineAge = 45, 
                               #endBaselineAge = 48,
                               J=10000,          #10000
                               Ndraws= 1,                  #set to 1 because we are computing the ground truth
                               #Nskip = 1000, 
                               #Ntree = 100, 
                               #Keepevery = 5, 
                               Suppress = TRUE, 
                               By = Ndraws/2, 
                               Sparse = TRUE)







S_zzstar = get_S_zzstar(P_zzstar)



P_zz = MCInteg_param_model(data = input_data1,
                           var.type = input_varType1,
                           #var.type.temp.order = input_varType2,
                           PModels = ground_truth_parametric_coeffs_list,
                           fixed.regime.trt = rep(1,2),
                           fixed.regime.control = rep(1,2),
                           #conditionOnBaselineAge = TRUE, 
                           #startBaselineAge = 45, 
                           #endBaselineAge = 48,
                           J=10000, 
                           Ndraws= 1, 
                           #Nskip = 1000, 
                           #Ntree = 100, 
                           #Keepevery = 5,   
                           Suppress = TRUE, 
                           By = Ndraws/2, 
                           Sparse = TRUE)







S_zz = get_S_zzstar(P_zz)



P_zstarzstar = MCInteg_param_model(data = input_data1,
                                   var.type = input_varType1,
                                   #var.type.temp.order = input_varType2,
                                   PModels = ground_truth_parametric_coeffs_list,
                                   fixed.regime.trt = rep(0,2),
                                   fixed.regime.control = rep(0,2),
                                   #conditionOnBaselineAge = TRUE, 
                                   #startBaselineAge = 45, 
                                   #endBaselineAge = 48,
                                   J=10000, 
                                   Ndraws= 1, 
                                   #Nskip = 1000, 
                                   #Ntree = 100, 
                                   #Keepevery = 5,  
                                   Suppress = TRUE, 
                                   By = Ndraws/2, 
                                   Sparse = TRUE)







S_zstarzstar = get_S_zzstar(P_zstarzstar)


IDE = get_IDE(S_zzstar,S_zstarzstar)


IIE = get_IIE(S_zz,S_zzstar)





Total_Effects = get_TotalEffects(IDE,IIE)





true_causal_estimands = data.frame(IDE = IDE, IIE= IIE, TE= Total_Effects)
true_causal_estimands


###############################################################





#######################################################################################################
#######################Change the local directories accordingly########################################
#######################################################################################################


## Results after running the simulation 1000 times


library(ggplot2)
library(gridExtra)


sim_BART_result_df = data.frame(data.table::fread("~/Documents/My Computer/Project 1-Discrete time Mediation Analysis/Codes and Results/Project Coding/GcompBART2/R/Results_ARIC_BART_Models_result.txt", sep = "\t"))

head(sim_BART_result_df)

sim_Parametric_result_df = data.frame(data.table::fread("~/Documents/My Computer/Project 1-Discrete time Mediation Analysis/Codes and Results/Project Coding/GcompBART2/R/Results_ARIC_Stan_Param_Models_result.txt", sep = "\t"))


head(sim_Parametric_result_df)






#######################################################################################################
#######################################################################################################

## Compute bias, coverage and MSE





calculate_statistics = function(results_df, true_value) {
  # Calculate Bias
  bias = as.numeric(results_df$posterior_mean) - true_value
  
  # Calculate Coverage (1 if true_value is within the CI, 0 otherwise)
  coverage = ifelse(as.numeric(results_df$lower_CI) <= true_value & as.numeric(results_df$upper_CI) >= true_value, 1, 0)
  
  # Calculate MSE (Mean Squared Error)
  mse = (as.numeric(results_df$posterior_mean) - true_value)^2
  
  # Summarize Bias, Coverage, and MSE
  result = data.frame(
    bias = sprintf("%.5e", mean(bias)),
    coverage = mean(coverage),
    mse = sprintf("%.5e", mean(mse))
  )
  
  return(result)
}




results_table_visit = function(IDE_df_visit_j, IIE_df_visit_j, TE_df_visit_j, true_value_df, visit){
  cat("Results at visit ", visit, ":\n")
  result_df = calculate_statistics(IDE_df_visit_j, true_causal_estimands[visit,1])
  result_df = rbind( result_df, calculate_statistics(IIE_df_visit_j, true_causal_estimands[visit,2]))
  result_df = rbind( result_df, calculate_statistics(TE_df_visit_j, true_causal_estimands[visit,3]))
  colnames(result_df) = c("BIAS",  "95% CI coverage rate","MSE")
  rownames(result_df) = c("IDE", "IIE", "TE")
  print(result_df)
}




BART_IDE_result_df_visit2 = data.frame(posterior_mean = sim_BART_result_df$IDE2_post_mean,
                                       lower_CI =  sim_BART_result_df$IDE2_lowerCI,
                                       upper_CI = sim_BART_result_df$IDE2_upperCI)


BART_IIE_result_df_visit2 = data.frame(posterior_mean = sim_BART_result_df$IIE2_post_mean,
                                       lower_CI =  sim_BART_result_df$IIE2_lowerCI,
                                       upper_CI = sim_BART_result_df$IIE2_upperCI)


BART_TE_result_df_visit2 = data.frame(posterior_mean = sim_BART_result_df$TE2_post_mean,
                                      lower_CI =  sim_BART_result_df$TE2_lowerCI,
                                      upper_CI = sim_BART_result_df$TE2_upperCI)




BART_IDE_result_df_visit3 = data.frame(posterior_mean = sim_BART_result_df$IDE3_post_mean,
                                       lower_CI =  sim_BART_result_df$IDE3_lowerCI,
                                       upper_CI = sim_BART_result_df$IDE3_upperCI)

BART_IIE_result_df_visit3 = data.frame(posterior_mean = sim_BART_result_df$IIE3_post_mean,
                                       lower_CI =  sim_BART_result_df$IIE3_lowerCI,
                                       upper_CI = sim_BART_result_df$IIE3_upperCI)

BART_TE_result_df_visit3 = data.frame(posterior_mean = sim_BART_result_df$TE3_post_mean,
                                      lower_CI =  sim_BART_result_df$TE3_lowerCI,
                                      upper_CI = sim_BART_result_df$TE3_upperCI)



Parametric_IDE_result_df_visit2 = data.frame(posterior_mean = sim_Parametric_result_df$IDE2_post_mean,
                                             lower_CI =  sim_Parametric_result_df$IDE2_lowerCI,
                                             upper_CI = sim_Parametric_result_df$IDE2_upperCI)

Parametric_IIE_result_df_visit2 = data.frame(posterior_mean = sim_Parametric_result_df$IIE2_post_mean,
                                             lower_CI =  sim_Parametric_result_df$IIE2_lowerCI,
                                             upper_CI = sim_Parametric_result_df$IIE2_upperCI)

Parametric_TE_result_df_visit2 = data.frame(posterior_mean = sim_Parametric_result_df$TE2_post_mean,
                                            lower_CI =  sim_Parametric_result_df$TE2_lowerCI,
                                            upper_CI = sim_Parametric_result_df$TE2_upperCI)

Parametric_IDE_result_df_visit3 = data.frame(posterior_mean = sim_Parametric_result_df$IDE3_post_mean,
                                             lower_CI =  sim_Parametric_result_df$IDE3_lowerCI,
                                             upper_CI = sim_Parametric_result_df$IDE3_upperCI)

Parametric_IIE_result_df_visit3 = data.frame(posterior_mean = sim_Parametric_result_df$IIE3_post_mean,
                                             lower_CI =  sim_Parametric_result_df$IIE3_lowerCI,
                                             upper_CI = sim_Parametric_result_df$IIE3_upperCI)

Parametric_TE_result_df_visit3 = data.frame(posterior_mean = sim_Parametric_result_df$TE3_post_mean,
                                            lower_CI =  sim_Parametric_result_df$TE3_lowerCI,
                                            upper_CI = sim_Parametric_result_df$TE3_upperCI)






## Print Results at visit 2

results_table_visit(Parametric_IDE_result_df_visit2 , Parametric_IIE_result_df_visit2 , Parametric_TE_result_df_visit2 ,
                    true_causal_estimands, 2)

results_table_visit(BART_IDE_result_df_visit2 , BART_IIE_result_df_visit2 , BART_TE_result_df_visit2 ,
                    true_causal_estimands, 2)



## Print Results at visit 3

results_table_visit(Parametric_IDE_result_df_visit3 , Parametric_IIE_result_df_visit3 , Parametric_TE_result_df_visit3 ,
                    true_causal_estimands, 3)

results_table_visit(BART_IDE_result_df_visit3 , BART_IIE_result_df_visit3 , BART_TE_result_df_visit3 ,
                    true_causal_estimands, 3)



