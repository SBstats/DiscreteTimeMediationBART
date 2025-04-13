#-------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------Prepare for simulation----------------------------------------------#

run_ID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(run_ID)

txt.title = paste0("param_simu_results/Results_ARIC_BART_MissModels_result.txt")
if (run_ID == 1) {
  df = data.frame(matrix(ncol = 19, nrow = 0))                 #each row contains results from each replication/dataset/run
  df_col_names = c("run_ID",
                   "IDE2_post_mean", "IDE2_lowerCI", "IDE2_upperCI",
                   "IIE2_post_mean", "IIE2_lowerCI", "IIE2_upperCI",
                   "TE2_post_mean", "TE2_lowerCI", "TE2_upperCI",
                   "IDE3_post_mean", "IDE3_lowerCI", "IDE3_upperCI",
                   "IIE3_post_mean", "IIE3_lowerCI", "IIE3_upperCI",
                   "TE3_post_mean", "TE3_lowerCI", "TE3_upperCI")
  colnames(df) = df_col_names
  write.table(df, file = txt.title, sep = "\t", row.names = FALSE, col.names = TRUE)
}

AGE_lg_coeff = -1.5; BMI_lg_coeff = 1.5; BMIXSEX_coeff = 0.25;


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ----------LOAD CODE FOR ONE REPLICATION HERE----------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------Fit Models in ARIC-------------------------------------------------#




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

#library(rstan)
#options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)


ARIC_DATA = read_sas('ARIC.sas7bdat')
ARIC_df = data.frame(ARIC_DATA)
#head(ARIC_df)
ARIC_workdf = data.frame(ID_FINAL = ARIC_df$ID_FINAL,  EXAM = ARIC_df$EXAM,RACE = ARIC_df$RACE, BMI = ARIC_df$BMINOW, 
                         AGE = ARIC_df$AGE,SEX = ARIC_df$SEX,EDU = ARIC_df$edu_g, RXHYP = ARIC_df$RXHYP, 
                         SMOKE_STATUS_BIN =ARIC_df$SMOKER, SMOKE_STATUS_MULT =ARIC_df$SMOKE_STATUS, SBP = ARIC_df$SBP,
                         DBP = ARIC_df$DBP, FOLLOWUPYRS = ARIC_df$LENYFL, CVD_DTH = ARIC_df$CVD_DTH, 
                         TOT_DTH =ARIC_df$TOT_DTH, GLUCOSE = ARIC_df$GLUCOSE, HXDIAB = ARIC_df$HXDIAB, 
                         TOT_HDL_CHL_RATIO = ARIC_df$TOTCHL/ARIC_df$HDLCHL)

#head(ARIC_workdf)

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
#head(ARIC_workdf_wide)



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
#head(ARIC_workdf_wide)

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








#Note that the number of people who died between E3 and E4 is sum(Y4) - sum(Y3)
ARIC_workdf_wide$Y1 = Y1
ARIC_workdf_wide$Y2 = Y2
ARIC_workdf_wide$Y3 = Y3
ARIC_workdf_wide$Y4 = Y4
ARIC_workdf_wide$delta = delta
ARIC_workdf_wide$d = d
#ARIC_workdf_wide$AllDeathIndicator = AllDeathIndicator
#head(ARIC_workdf_wide)










#..................................Ykj where k = 1(main event) vs 2(competing event)...........//


Y11 = rep(0,length(ARIC_workdf_wide$ID_FINAL))   #Assumption that everyone is alive at the first visit

Y12 = ifelse(!(is.na(ARIC_workdf_wide$AGE_1)) & is.na(ARIC_workdf_wide$AGE_2) & is.na(ARIC_workdf_wide$AGE_3) & is.na(ARIC_workdf_wide$AGE_4) & d ==1 , 1,0) 
Y12 = ifelse((ARIC_workdf_wide$ageAtDeathOrCensor > ARIC_workdf_wide$AGE_1) &
               (ARIC_workdf_wide$ageAtDeathOrCensor <= AGE_2)  & 
               d ==1,
             1,Y12) 


Y13 = ifelse(!(is.na(ARIC_workdf_wide$AGE_2)) & is.na(ARIC_workdf_wide$AGE_3) & is.na(ARIC_workdf_wide$AGE_4)  & d ==1 , 1,0) 
Y13 = ifelse((ARIC_workdf_wide$ageAtDeathOrCensor > ARIC_workdf_wide$AGE_2) &
               (ARIC_workdf_wide$ageAtDeathOrCensor <= AGE_3) & 
               d ==1,
             1,Y13) 



Y14 = ifelse(!(is.na(ARIC_workdf_wide$AGE_3)) & is.na(ARIC_workdf_wide$AGE_4)   & d ==1 , 1,0) 
Y14 = ifelse((ARIC_workdf_wide$ageAtDeathOrCensor > ARIC_workdf_wide$AGE_3) &
               (ARIC_workdf_wide$ageAtDeathOrCensor <= AGE_4)  & 
               d ==1,
             1,Y14) 


# Y15 = ifelse(!(is.na(ARIC_workdf_wide$AGE_4)) & is.na(ARIC_workdf_wide$AGE_5)   & d ==1 , 1,0) 
# Y15 = ifelse((ARIC_workdf_wide$ageAtDeathOrCensor > ARIC_workdf_wide$AGE_4) &
#               (ARIC_workdf_wide$ageAtDeathOrCensor <= AGE_5)  & 
#               d ==1,
#             1,Y15) 


Y15 = ifelse(!(is.na(ARIC_workdf_wide$AGE_4)) & is.na(ARIC_workdf_wide$AGE_5)   & d ==1 , 1,0) 
Y15 = ifelse((ARIC_workdf_wide$ageAtDeathOrCensor > ARIC_workdf_wide$AGE_4) &
               d ==1,
             1,Y15) 





Y21 = rep(0,length(ARIC_workdf_wide$ID_FINAL))   #Assumption that everyone is alive at the first visit

Y22 = ifelse(!(is.na(ARIC_workdf_wide$AGE_1)) & is.na(ARIC_workdf_wide$AGE_2) & is.na(ARIC_workdf_wide$AGE_3) & is.na(ARIC_workdf_wide$AGE_4) & d ==2 , 1,0) 
Y22 = ifelse((ARIC_workdf_wide$ageAtDeathOrCensor > ARIC_workdf_wide$AGE_1) &
               (ARIC_workdf_wide$ageAtDeathOrCensor <= AGE_2)  & 
               d ==2,
             1,Y22) 


Y23 = ifelse(!(is.na(ARIC_workdf_wide$AGE_2)) & is.na(ARIC_workdf_wide$AGE_3) & is.na(ARIC_workdf_wide$AGE_4)  & d ==2 , 1,0) 
Y23 = ifelse((ARIC_workdf_wide$ageAtDeathOrCensor > ARIC_workdf_wide$AGE_2) &
               (ARIC_workdf_wide$ageAtDeathOrCensor <= AGE_3)  & 
               d ==2,
             1,Y23) 



Y24 = ifelse(!(is.na(ARIC_workdf_wide$AGE_3)) & is.na(ARIC_workdf_wide$AGE_4)   & d ==2 , 1,0) 
Y24 = ifelse((ARIC_workdf_wide$ageAtDeathOrCensor > ARIC_workdf_wide$AGE_3) &
               (ARIC_workdf_wide$ageAtDeathOrCensor <= AGE_4)  & 
               d ==2,
             1,Y24) 


# Y25 = ifelse(!(is.na(ARIC_workdf_wide$AGE_4)) & is.na(ARIC_workdf_wide$AGE_5)   & d ==2 , 1,0) 
# Y25 = ifelse((ARIC_workdf_wide$ageAtDeathOrCensor > ARIC_workdf_wide$AGE_4) &
#                (ARIC_workdf_wide$ageAtDeathOrCensor < AGE_5)  & 
#                d ==2,
#              1,Y25) 



Y25 = ifelse(!(is.na(ARIC_workdf_wide$AGE_4)) & is.na(ARIC_workdf_wide$AGE_5)   & d ==2 , 1,0) 
Y25 = ifelse((ARIC_workdf_wide$ageAtDeathOrCensor > ARIC_workdf_wide$AGE_4) &
               d ==2,
             1,Y25) 













#Note that the number of people who died between E3 and E4 is sum(Y4) - sum(Y3)
ARIC_workdf_wide$Y11 = Y11
ARIC_workdf_wide$Y12 = Y12
ARIC_workdf_wide$Y13 = Y13
ARIC_workdf_wide$Y14 = Y14



ARIC_workdf_wide$Y21 = Y21
ARIC_workdf_wide$Y22 = Y22
ARIC_workdf_wide$Y23 = Y23
ARIC_workdf_wide$Y24 = Y24




#-------------------------------------------------------------------------------------------------------------------#
# Individuals who are hypertensive (at baseline) in ARIC

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

#-------------------------------------------------------------------------------------------------------------------#


















#-------------------------------------------------------------------------------------------------------------------#
#------------------------------Fit Parametric Regressions-----------------------------------------------------------#






############################################################
####....... INPUT DATA FOR PARAMETRIC REGRESSION.....#########


param_reg_df_Y2_true = data.frame(V0 = ARIC_workdf_wide$SMOKE_STATUS_MULT_1, 
                                  BMI = ARIC_workdf_wide$BMI_1 ,
                                  AGE = ARIC_workdf_wide$AGE_1,
                                  RACE = ARIC_workdf_wide$RACE, 
                                  SEX = ARIC_workdf_wide$SEX,
                                  EDU = ARIC_workdf_wide$EDU,
                                  BMI_lg = log(ARIC_workdf_wide$BMI_1),
                                  AGE_lg = log(ARIC_workdf_wide$AGE_1),
                                  BMIXSEX = ARIC_workdf_wide$BMI_1 *ARIC_workdf_wide$SEX,
                                  Z_1 = ARIC_workdf_wide$RXHYP_1, 
                                  L_1 = ARIC_workdf_wide$SMOKE_STATUS_BIN_1,
                                  M_1 = ARIC_workdf_wide$MeanBP_1,
                                  Y_2 = ARIC_workdf_wide$Y2)



param_reg_df_Y3_true = data.frame(V0 = ARIC_workdf_wide$SMOKE_STATUS_MULT_1, 
                                  BMI = ARIC_workdf_wide$BMI_1 ,
                                  AGE = ARIC_workdf_wide$AGE_1,
                                  RACE = ARIC_workdf_wide$RACE, 
                                  SEX = ARIC_workdf_wide$SEX,
                                  EDU = ARIC_workdf_wide$EDU,
                                  BMI_lg = log(ARIC_workdf_wide$BMI_1),
                                  AGE_lg = log(ARIC_workdf_wide$AGE_1),
                                  BMIXSEX = ARIC_workdf_wide$BMI_1 *ARIC_workdf_wide$SEX,
                                  Z_1 = ARIC_workdf_wide$RXHYP_1, 
                                  L_1 = ARIC_workdf_wide$SMOKE_STATUS_BIN_1, 
                                  M_1 = ARIC_workdf_wide$MeanBP_1,
                                  Z_2 = ARIC_workdf_wide$RXHYP_2, 
                                  L_2 = ARIC_workdf_wide$SMOKE_STATUS_BIN_2,
                                  M_2 = ARIC_workdf_wide$MeanBP_2,
                                  Y_3 = ARIC_workdf_wide$Y3)



param_reg_df_ZLM_true= data.frame(V0 = ARIC_workdf_wide$SMOKE_STATUS_MULT_1, 
                                  BMI = ARIC_workdf_wide$BMI_1 ,
                                  AGE = ARIC_workdf_wide$AGE_1,
                                  RACE = ARIC_workdf_wide$RACE, 
                                  SEX = ARIC_workdf_wide$SEX,
                                  EDU = ARIC_workdf_wide$EDU,
                                  BMI_lg = log(ARIC_workdf_wide$BMI_1),
                                  AGE_lg = log(ARIC_workdf_wide$AGE_1),
                                  BMIXSEX = ARIC_workdf_wide$BMI_1 *ARIC_workdf_wide$SEX,
                                  Z_1 = ARIC_workdf_wide$RXHYP_1, 
                                  L_1 = ARIC_workdf_wide$SMOKE_STATUS_BIN_1, 
                                  M_1 = ARIC_workdf_wide$MeanBP_1,
                                  Z_2 = ARIC_workdf_wide$RXHYP_2, 
                                  L_2 = ARIC_workdf_wide$SMOKE_STATUS_BIN_2,
                                  M_2 = ARIC_workdf_wide$MeanBP_2)


# param_reg_df_Y2 = data.frame(V0 = ARIC_workdf_wide$SMOKE_STATUS_MULT_1, BMI = ARIC_workdf_wide$BMI_1 ,
#                              AGE = ARIC_workdf_wide$AGE_1, RACE = ARIC_workdf_wide$RACE, 
#                              SEX = ARIC_workdf_wide$SEX, EDU = ARIC_workdf_wide$EDU, 
#                              Z_1 = ARIC_workdf_wide$RXHYP_1, 
#                              L_1 = ARIC_workdf_wide$SMOKE_STATUS_BIN_1, 
#                              M_1 = ARIC_workdf_wide$MeanBP_1,  
#                              Y_2 = ARIC_workdf_wide$Y2)
# 
# 
# param_reg_df_Y3 = data.frame(V0 = ARIC_workdf_wide$SMOKE_STATUS_MULT_1, BMI = ARIC_workdf_wide$BMI_1 ,
#                              AGE = ARIC_workdf_wide$AGE_1, RACE = ARIC_workdf_wide$RACE, 
#                              SEX = ARIC_workdf_wide$SEX, EDU = ARIC_workdf_wide$EDU, 
#                              Z_1 = ARIC_workdf_wide$RXHYP_1, 
#                              L_1 = ARIC_workdf_wide$SMOKE_STATUS_BIN_1, 
#                              M_1 = ARIC_workdf_wide$MeanBP_1,
#                              Z_2 = ARIC_workdf_wide$RXHYP_2, 
#                              L_2 = ARIC_workdf_wide$SMOKE_STATUS_BIN_2,
#                              M_2 = ARIC_workdf_wide$MeanBP_2, 
#                              Y_3 = ARIC_workdf_wide$Y3)
# 
# 
# 
# param_reg_df_ZLM= data.frame(V0 = ARIC_workdf_wide$SMOKE_STATUS_MULT_1, BMI = ARIC_workdf_wide$BMI_1 ,
#                              AGE = ARIC_workdf_wide$AGE_1, RACE = ARIC_workdf_wide$RACE, 
#                              SEX = ARIC_workdf_wide$SEX, EDU = ARIC_workdf_wide$EDU, 
#                              Z_1 = ARIC_workdf_wide$RXHYP_1, 
#                              L_1 = ARIC_workdf_wide$SMOKE_STATUS_BIN_1, 
#                              M_1 = ARIC_workdf_wide$MeanBP_1,  
#                              Z_2 = ARIC_workdf_wide$RXHYP_2, 
#                              L_2 = ARIC_workdf_wide$SMOKE_STATUS_BIN_2,
#                              M_2 = ARIC_workdf_wide$MeanBP_2)
# 




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




#-------------------------------------------------------------------------------------------------------------------#
#------------------------------Fit Parametric Regressions with real data--------------------------------------------#

if(is_binary(param_reg_df_Y2_true$Y_2) == TRUE){
  fit_glm_Y2 = glm(Y_2 ~  V0 +  RACE + SEX + EDU + Z_1 + L_1 + 
                     M_1, data = param_reg_df_Y2_true, family = binomial(link="logit"))
  fit_glm_Y2_coeff = as.data.frame(coef(fit_glm_Y2))
  print(dim(fit_glm_Y2_coeff))
  fit_glm_Y2_coeff = rbind(fit_glm_Y2_coeff[1:5, , drop = FALSE],
                           BMI_lg_coeff, AGE_lg_coeff, BMIXSEX_coeff,
                           fit_glm_Y2_coeff[6:nrow(fit_glm_Y2_coeff), , drop = FALSE])
}else{
  fit_glm_Y2 = lm(Y_2 ~  V0 +  RACE + SEX + EDU + Z_1 + L_1 + 
                    M_1, data = param_reg_df_Y2_true)
  fit_glm_Y2_coeff = as.data.frame(coef(fit_glm_Y2))
  fit_glm_Y2_coeff = rbind(fit_glm_Y2_coeff[1:5, , drop = FALSE], 
                           BMI_lg_coeff, AGE_lg_coeff, BMIXSEX_coeff,
                           fit_glm_Y2_coeff[6:nrow(fit_glm_Y2_coeff), , drop = FALSE])
}




# formula_Y3= as.formula("Y_3 ~ V0 + BMI + I(BMI^2) + AGE+ I(AGE^2) + V0 +  RACE + SEX + EDU + Z_1 + L_1 
#                 + M_1 + Z_2 + L_2 +M_2")
if(is_binary(param_reg_df_Y3_true$Y_3) == TRUE){
  fit_glm_Y3 = glm(Y_3 ~  V0 +  RACE + SEX + EDU + Z_1 + L_1 + M_1
                   + Z_2 + L_2 +M_2, data = param_reg_df_Y3_true, family = binomial(link="logit"))
  fit_glm_Y3_coeff = as.data.frame(coef(fit_glm_Y3))
  fit_glm_Y3_coeff = rbind(fit_glm_Y3_coeff[1:5, , drop = FALSE], 
                           BMI_lg_coeff, AGE_lg_coeff, BMIXSEX_coeff,
                           fit_glm_Y3_coeff[6:nrow(fit_glm_Y3_coeff), , drop = FALSE])
}else{
  fit_glm_Y3 = lm(Y_3 ~  V0 +  RACE + SEX + EDU + Z_1 + L_1 + M_1
                  + Z_2 + L_2 +M_2, data = param_reg_df_Y3_true)
  fit_glm_Y3_coeff = as.data.frame(coef(fit_glm_Y3))
  fit_glm_Y3_coeff = rbind(fit_glm_Y3_coeff[1:5, , drop = FALSE],
                           BMI_lg_coeff, AGE_lg_coeff, BMIXSEX_coeff,
                           fit_glm_Y3_coeff[6:nrow(fit_glm_Y3_coeff), , drop = FALSE])
}


#formula_Z1= as.formula("Z_1 ~  V0 +  RACE + SEX + EDU")
if(is_binary(param_reg_df_ZLM_true$Z_1) == TRUE){
  fit_glm_Z1 = glm(Z_1 ~  V0 +  RACE + SEX + EDU, 
                   data = param_reg_df_ZLM_true, family = binomial(link="logit"))
  fit_glm_Z1_coeff = as.data.frame(coef(fit_glm_Z1))
  fit_glm_Z1_coeff = rbind(fit_glm_Z1_coeff[1:5, , drop = FALSE], 
                           BMI_lg_coeff, AGE_lg_coeff, BMIXSEX_coeff)
}else{
  fit_glm_Z1 = lm(Z_1 ~  V0 +  RACE + SEX + EDU,
                  data = param_reg_df_ZLM_true)
  fit_glm_Z1_coeff = as.data.frame(coef(fit_glm_Z1))
  fit_glm_Z1_coeff = rbind(fit_glm_Z1_coeff[1:5, , drop = FALSE], 
                           BMI_lg_coeff, AGE_lg_coeff, BMIXSEX_coeff)
}

# formula_Z2= as.formula("Z_2 ~  V0 +  RACE + SEX + EDU + Z_1 + L_1 +   
#                        M_1")
if(is_binary(param_reg_df_ZLM_true$Z_2) == TRUE){
  fit_glm_Z2 = glm(Z_2 ~  V0 +  RACE + SEX + EDU + Z_1 + L_1 +   
                     M_1, data = param_reg_df_ZLM_true, family = binomial(link="logit"))
  fit_glm_Z2_coeff = as.data.frame(coef(fit_glm_Z2))
  fit_glm_Z2_coeff = rbind(fit_glm_Z2_coeff[1:5, , drop = FALSE], 
                           BMI_lg_coeff, AGE_lg_coeff, BMIXSEX_coeff,
                           fit_glm_Z2_coeff[6:nrow(fit_glm_Z2_coeff), , drop = FALSE])
}else{
  fit_glm_Z2 = lm(Z_2 ~  V0 +  RACE + SEX + EDU + Z_1 + L_1 +   
                    M_1, data = param_reg_df_ZLM_true)
  fit_glm_Z2_coeff = as.data.frame(coef(fit_glm_Z2))
  fit_glm_Z2_coeff = rbind(fit_glm_Z2_coeff[1:5, , drop = FALSE], 
                           BMI_lg_coeff, AGE_lg_coeff, BMIXSEX_coeff,
                           fit_glm_Z2_coeff[6:nrow(fit_glm_Z2_coeff), , drop = FALSE])
}





#formula_L1= as.formula("L_1 ~  V0 +  RACE + SEX + EDU + Z_1")
if(is_binary(param_reg_df_ZLM_true$L_1) == TRUE){
  fit_glm_L1 = glm(L_1 ~  V0 +  RACE + SEX + EDU + Z_1, 
                   data = param_reg_df_ZLM_true, family = binomial(link="logit"))
  fit_glm_L1_coeff = as.data.frame(coef(fit_glm_L1))
  fit_glm_L1_coeff = rbind(fit_glm_L1_coeff[1:5, , drop = FALSE],
                           BMI_lg_coeff, AGE_lg_coeff, BMIXSEX_coeff,
                           fit_glm_L1_coeff[6:nrow(fit_glm_L1_coeff), , drop = FALSE])
}else{
  fit_glm_L1 = lm(L_1 ~  V0 +  RACE + SEX + EDU + Z_1,
                  data = param_reg_df_ZLM_true)
  fit_glm_L1_coeff = as.data.frame(coef(fit_glm_L1))
  fit_glm_L1_coeff = rbind(fit_glm_L1_coeff[1:5, , drop = FALSE],
                           BMI_lg_coeff, AGE_lg_coeff, BMIXSEX_coeff,
                           fit_glm_L1_coeff[6:nrow(fit_glm_L1_coeff), , drop = FALSE])
}


# formula_L2= as.formula("L_2 ~  V0 +  RACE + SEX + EDU + Z_1 + L_1 
#             + M_1+ Z_2")
if(is_binary(param_reg_df_ZLM_true$L_2) == TRUE){
  fit_glm_L2 = glm(L_2 ~  V0 +  RACE + SEX + EDU + Z_1 + L_1 + M_1
                   + Z_2, data = param_reg_df_ZLM_true, family = binomial(link="logit"))
  fit_glm_L2_coeff = as.data.frame(coef(fit_glm_L2))
  fit_glm_L2_coeff = rbind(fit_glm_L2_coeff[1:5, , drop = FALSE],
                           BMI_lg_coeff, AGE_lg_coeff, BMIXSEX_coeff,
                           fit_glm_L2_coeff[6:nrow(fit_glm_L2_coeff), , drop = FALSE])
}else{
  fit_glm_L2 = lm(L_2 ~  V0 +  RACE + SEX + EDU + Z_1 + L_1 + M_1
                  + Z_2, data = param_reg_df_ZLM_true)
  fit_glm_L2_coeff = as.data.frame(coef(fit_glm_L2))
  fit_glm_L2_coeff = rbind(fit_glm_L2_coeff[1:5, , drop = FALSE], 
                           BMI_lg_coeff, AGE_lg_coeff, BMIXSEX_coeff,
                           fit_glm_L2_coeff[6:nrow(fit_glm_L2_coeff), , drop = FALSE])
}





#formula_M1= as.formula("M_1 ~  V0 +  RACE + SEX + EDU + Z_1 + L_1")
if(is_binary(param_reg_df_ZLM_true$M_1) == TRUE){
  fit_glm_M1 = glm(M_1 ~  V0 +  RACE + SEX + EDU + Z_1 + L_1,
                   data = param_reg_df_ZLM_true, family = binomial(link="logit"))
  fit_glm_M1_coeff = as.data.frame(coef(fit_glm_M1))
  fit_glm_M1_coeff = rbind(fit_glm_M1_coeff[1:5, , drop = FALSE], 
                           BMI_lg_coeff, AGE_lg_coeff, BMIXSEX_coeff,
                           fit_glm_M1_coeff[6:nrow(fit_glm_M1_coeff), , drop = FALSE])
}else{
  fit_glm_M1 = lm(M_1 ~  V0 +  RACE + SEX + EDU + Z_1 + L_1,
                  data = param_reg_df_ZLM_true)
  fit_glm_M1_coeff = as.data.frame(coef(fit_glm_M1))
  fit_glm_M1_coeff = rbind(fit_glm_M1_coeff[1:5, , drop = FALSE], 
                           BMI_lg_coeff, AGE_lg_coeff, BMIXSEX_coeff,
                           fit_glm_M1_coeff[6:nrow(fit_glm_M1_coeff), , drop = FALSE])
}

# formula_M2= as.formula("M_2 ~  V0 +  RACE + SEX + EDU + Z_1 + L_1 
#             + M_1+ Z_2 + L_2")
if(is_binary(param_reg_df_ZLM_true$M_2) == TRUE){
  fit_glm_M2 = glm(M_2 ~  V0 +  RACE + SEX + EDU + Z_1 + L_1 + M_1
                   + Z_2 + L_2, data = param_reg_df_ZLM_true, family = binomial(link="logit"))
  fit_glm_M2_coeff = as.data.frame(coef(fit_glm_M2))
  fit_glm_M2_coeff = rbind(fit_glm_M2_coeff[1:5, , drop = FALSE], 
                           BMI_lg_coeff, AGE_lg_coeff, BMIXSEX_coeff,
                           fit_glm_M2_coeff[6:nrow(fit_glm_M2_coeff), , drop = FALSE],
                           sigma = sigma(fit_glm_M2))
}else{
  fit_glm_M2 = lm(M_2 ~  V0 +  RACE + SEX + EDU + Z_1 + L_1 + M_1
                  + Z_2 + L_2, data = param_reg_df_ZLM_true)
  fit_glm_M2_coeff = as.data.frame(coef(fit_glm_M2))
  fit_glm_M2_coeff = rbind(fit_glm_M2_coeff[1:5, , drop = FALSE], 
                           BMI_lg_coeff, AGE_lg_coeff, BMIXSEX_coeff,
                           fit_glm_M2_coeff[6:nrow(fit_glm_M2_coeff), , drop = FALSE])
  
}
#------------------------------------------------------------------------------------------------







#-------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------Generate one dataset (1 run in Hipergator)----------------------------#

L0_df_data_gen = data.frame(V0 = ARIC_workdf_wide$SMOKE_STATUS_MULT_1, 
                            RACE = ARIC_workdf_wide$RACE, 
                            SEX = ARIC_workdf_wide$SEX,
                            EDU = ARIC_workdf_wide$EDU,
                            BMI_lg = log(ARIC_workdf_wide$BMI_1),
                            AGE_lg = log(ARIC_workdf_wide$AGE_1),
                            BMIXSEX = ARIC_workdf_wide$BMI_1 *ARIC_workdf_wide$SEX)

L0_df_data_gen = na.omit(L0_df_data_gen)

L0_df = data.frame(V0 = ARIC_workdf_wide$SMOKE_STATUS_MULT_1,
                   BMI = ARIC_workdf_wide$BMI_1,
                   AGE = ARIC_workdf_wide$AGE_1,
                   RACE = ARIC_workdf_wide$RACE, 
                   SEX = ARIC_workdf_wide$SEX,
                   EDU = ARIC_workdf_wide$EDU)

L0_df = na.omit(L0_df)

#L0_df = rbind(L0_df,L0_df,L0_df) #increasing the sample size for each replication
sample_size = nrow(L0_df_data_gen)


generate_longitudinal_data_bin = function(sample_size_for_pred,     
                                          fitted_model_coef,
                                          baseline_L0_df,
                                          longitudinal_covariates_df = NULL){
  
  intercept = rep(1,nrow(baseline_L0_df))
  if(is.null(longitudinal_covariates_df)){
    covariates_df = cbind(intercept, baseline_L0_df)
  } else{
    covariates_df = cbind(intercept,baseline_L0_df, longitudinal_covariates_df)
  }
  covariates_mat = as.matrix(covariates_df, nrow=nrow(covariates_df), ncol = ncol(covariates_df))
  fitted_model_coef_mat = as.matrix(fitted_model_coef, nrow = nrow(fitted_model_coef),
                                    ncol = ncol(fitted_model_coef))
  #print(head(covariates_df))
  prob_vec = as.numeric(plogis(covariates_mat %*% fitted_model_coef_mat))
  return(rbinom(sample_size_for_pred,1,prob_vec))
}


generate_longitudinal_data_cont = function(sample_size_for_pred,     
                                           fitted_model_coef,
                                           sigma_hat,
                                           baseline_L0_df,
                                           longitudinal_covariates_df = NULL
){
  intercept = rep(1,nrow(baseline_L0_df))
  if(is.null(longitudinal_covariates_df)){
    covariates_df = cbind(intercept, baseline_L0_df)
  } else{
    covariates_df = cbind(intercept,baseline_L0_df, longitudinal_covariates_df)
  }
  covariates_mat = as.matrix(covariates_df, nrow=nrow(covariates_df), ncol = ncol(covariates_df))
  fitted_model_coef_mat = as.matrix(fitted_model_coef, nrow = nrow(fitted_model_coef),
                                    ncol = ncol(fitted_model_coef))
  #print(head(covariates_df))
  mean_vec = as.numeric(covariates_mat %*% fitted_model_coef_mat)
  
  return(rnorm(sample_size_for_pred,mean = mean_vec,sd = sigma_hat))
  
}



#-------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------#
generate_longitudinal_data_bin_old = function(sample_size_for_pred,     
                                              fitted_model,
                                              baseline_L0_df,
                                              longitudinal_covariates_df = NULL){
  
  intercept = rep(1,nrow(baseline_L0_df))
  if(is.null(longitudinal_covariates_df)){
    covariates_df = cbind(intercept, baseline_L0_df)
  } else{
    covariates_df = cbind(intercept,baseline_L0_df, longitudinal_covariates_df)
  }
  #print(head(covariates_df))
  prob_vec = as.numeric(predict(object = fitted_model,
                                newdata = covariates_df,
                                type = "response"))
  return(rbinom(sample_size_for_pred,1,prob_vec))
}



generate_longitudinal_data_cont_old = function(sample_size_for_pred,     
                                               fitted_model,
                                               baseline_L0_df,
                                               longitudinal_covariates_df = NULL){
  intercept = rep(1,nrow(baseline_L0_df))
  if(is.null(longitudinal_covariates_df)){
    covariates_df = cbind(intercept, baseline_L0_df)
  } else{
    covariates_df = cbind(intercept,baseline_L0_df, longitudinal_covariates_df)
  }
  #print(head(covariates_df))
  mean_vec = as.numeric(predict(object = fitted_model,
                                newdata = covariates_df,
                                type = "response"))
  sd_val = sigma(fitted_model)
  return(rnorm(sample_size_for_pred,mean = mean_vec,sd = sd_val))
  
}
#-------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------#






#-------------------------------------------------------------------------------------------------------------------#

Y_1 = rep(0,sample_size)



Z_1 = generate_longitudinal_data_bin(sample_size_for_pred = sample_size, 
                                     fitted_model_coef= fit_glm_Z1_coeff,
                                     baseline_L0_df = L0_df_data_gen,
                                     longitudinal_covariates_df = NULL)


L_1 = generate_longitudinal_data_bin(sample_size_for_pred = sample_size, 
                                     fitted_model_coef = fit_glm_L1_coeff,
                                     baseline_L0_df = L0_df_data_gen,
                                     longitudinal_covariates_df = Z_1)

M_1 = generate_longitudinal_data_cont(sample_size_for_pred = sample_size, 
                                      fitted_model_coef = fit_glm_M1_coeff,
                                      sigma_hat = sigma(fit_glm_M1),
                                      baseline_L0_df = L0_df_data_gen,
                                      longitudinal_covariates_df = cbind(Z_1,L_1))



Y_2 = generate_longitudinal_data_bin(sample_size_for_pred = sample_size, 
                                     fitted_model_coef = fit_glm_Y2_coeff,
                                     baseline_L0_df = L0_df_data_gen,
                                     longitudinal_covariates_df = cbind(Z_1,L_1,M_1))



Z_2 = generate_longitudinal_data_bin(sample_size_for_pred = sample_size, 
                                     fitted_model_coef = fit_glm_Z2_coeff,
                                     baseline_L0_df = L0_df_data_gen,
                                     longitudinal_covariates_df = cbind(Z_1,L_1,M_1))

L_2 = generate_longitudinal_data_bin(sample_size_for_pred = sample_size, 
                                     fitted_model_coef = fit_glm_L2_coeff,
                                     baseline_L0_df = L0_df_data_gen,
                                     longitudinal_covariates_df = cbind(Z_1,L_1,M_1,Z_2))

M_2 = generate_longitudinal_data_cont(sample_size_for_pred = sample_size, 
                                      fitted_model_coef = fit_glm_M2_coeff,
                                      sigma_hat = sigma(fit_glm_M2),
                                      baseline_L0_df = L0_df_data_gen,
                                      longitudinal_covariates_df = cbind(Z_1,L_1,M_1,Z_2,L_2))


Y_3 = generate_longitudinal_data_bin(sample_size_for_pred = sample_size, 
                                     fitted_model_coef = fit_glm_Y3_coeff,
                                     baseline_L0_df = L0_df_data_gen,
                                     longitudinal_covariates_df = cbind(Z_1,L_1,M_1,Z_2,L_2,M_2))

delta_generated = ifelse(Y_2 ==1|Y_3 ==1,1,0)















#-------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------STAN CODE---------------------------------------------------------#

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





library(ggplot2)

calculate_95confidence_interval <- function(data) {
  conf_int <- function(data_row){
    n <- length(data_row)
    mean_value <- as.numeric(format(round(mean(data_row), 4), nsmall = 4))
    std_error <- sd(data_row) 
    #z_value <- qnorm(0.975,0,1)  # Z-value for a 95% confidence interval
    t_value <- qt(p = 1 - (0.05/2), df = n-1)
    lower <- mean_value - t_value * std_error
    lower <- as.numeric(format(round(lower, 4), nsmall = 4))
    upper <- mean_value + t_value * std_error
    upper <- as.numeric(format(round(upper, 4), nsmall = 4))
    interval <- cbind(lower, upper)
    return(interval)
  }
  return_mat = apply(data,1,conf_int)
  return(t(return_mat))
}



calculate_95bayesian_credible_interval <- function(data) {
  #cat("Posterior Means: \n")
  cred_int = function(data_row){
    posterior_samples <- data_row
    posterior_means <- as.numeric(format(round(mean(posterior_samples) , 4), nsmall = 4))
    #print(posterior_means)
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




#-------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------#






# Functions needed for MC Integration 


source("pred_surv.R")
source("pwbart_it.R")
source("expit.R")
source("quiet.R")
#source("pred_comb_Y.R")

source("pred_comb.R")
source("pred_y.R")









pred_comb_mediator <- function(cont, BM, MCdata_trt, MCdata_control, IT){
  if(cont == FALSE) {
    x_hat <- pnorm(pwbart_it(MCdata_control, BM$treedraws, mu = BM$binaryOffset, it=IT)) #make predictions conditional on counterfactual data
    new_MCdata <- cbind(MCdata_trt, rbinom(length(x_hat), 1, prob = x_hat))
  } else if(cont == TRUE) {
    x_hat <- pwbart_it(MCdata_control, BM$treedraws, mu = BM$mu, it = IT)
    new_MCdata <- cbind(MCdata_trt, rnorm(length(x_hat), mean = x_hat, sd = mean(BM$sigma)))
  }
  return(new_MCdata)
}






########################################################################################
######.............Function to fit BART Models (No Competing Event).............########


FitBART <- function(data, # Dataframe in temporal order: L0, Z,L,M,Y,D
                    var.type, # Vector of variable specifications for data.  Fi=fixed (e.g. the exposure), L=confounder, M= mediator, Y=outcome, D = event indicator. 
                    J=2000, # Size of pseudo data. Default is set to 2,000.
                    Ndraws=200, # Number of posterior draws. Default is set to 200.
                    Nskip = 100, # Number of burn-in samples. Default is set to 100.
                    Ntree = 100, # Number of trees. Default is set to 100.
                    Keepevery = 1, # Keep every k:th draw. Default is set to 1.
                    Suppress = TRUE, # Indicates if the output should be suppressed. Default is TRUE
                    By = Ndraws/10, # If Suppress is set to FALSE, output is provided for the By:th iteration.
                    Sparse = TRUE, # Indicates if the sparse dirichlet hyper prior should be used. Default is TRUE
                    ...
){
  
  
  
  
  ## test if the variables in data are continuous or binary  
  continuous <- apply(data, 2, function(x) !all(na.omit(x) %in% 0:1)) # Continuous + Multinomial
  
  
  n_Y <- length(which(var.type=="Y")) # Number of outcome variables, >1 allowed
  n_L0 <- length(which(var.type=="L0")) # Number of baseline confounder variables, >1 allowed
  n_Fi <- length(which(var.type=="Fi")) # Number of fixed variables, >1 allowed
  
  
  if (ncol(data) != length(var.type)) stop("The number of columns in data is not equal to the length of var.type ")
  
  
  
  
  
  
  ################################################################
  #Fit the BART models and store results as a list
  BModels <- vector("list", ncol(data) - (n_L0+1))  #(numL0Vars+1)
  for (i in (n_L0+1):(ncol(data)-1)) {  #(numL0Vars+1)
    # Identify complete observations
    id <- ifelse(apply(data[, 1:(1 + i)], 1, function(x) all(!is.na(x))), 1, 0)
    
    if(var.type[1+i] == "Fi") {
      BModels[[i]]  <- NA
    } else {
      if (continuous[1+i] == FALSE & Suppress == TRUE) {
        start_time <- Sys.time()
        quiet(BModels[[i]]  <- BART::pbart(data[id == 1, which(var.type[1:i]!="D")], data[id == 1, (1 + i)], ndpost = Ndraws, nkeeptrain = 0, nkeeptest = 0, ntree = Ntree, nskip = Nskip, sparse = Sparse, keepevery = Keepevery))
        end_time <- Sys.time()
        cat("BART fit time is :",end_time - start_time, ".\n")
      } else if (continuous[1+i] == TRUE & Suppress == TRUE) {
        start_time <- Sys.time()
        quiet(BModels[[i]]  <- BART::wbart(data[id == 1, which( var.type[1:i]!="D")], data[id == 1, (1 + i)], ndpost = Ndraws, nkeeptrain = 0, nkeeptest = 0, ntree = Ntree, nskip = Nskip, sparse = Sparse, keepevery = Keepevery))
        end_time <- Sys.time()
        cat("BART fit time is :",end_time - start_time, ".\n")
      } else if (continuous[1+i] == FALSE & Suppress == FALSE) {
        start_time <- Sys.time()
        BModels[[i]]  <- BART::pbart(data[id == 1, which(var.type[1:i]!="D")], data[id == 1, (1 + i)], ndpost = Ndraws, nkeeptrain = 0, nkeeptest = 0, ntree = Ntree, nskip = Nskip, sparse = Sparse, keepevery = Keepevery)
        end_time <- Sys.time()
        cat("BART fit time is :",end_time - start_time, ".\n")
      } else if (continuous[1+i] == TRUE & Suppress == FALSE) {
        start_time <- Sys.time()
        BModels[[i]]  <- BART::wbart(data[id == 1, which(var.type[1:i]!="D")], data[id == 1, (1 + i)], ndpost = Ndraws, nkeeptrain = 0, nkeeptest = 0, ntree = Ntree, nskip = Nskip, sparse = Sparse, keepevery = Keepevery)
        end_time <- Sys.time()
        cat("BART fit time is :",end_time - start_time, ".\n")
      }
    }
  }
  
  return(BModels) 
}















#-------------------------------------------------------------------------------------------------------------------#
#---------------------------------------G-computation---------------------------------------------------------------#




########################################################################################
######.............Function for MC integration (No Competing Event).............########
MCInteg <- function(data, # Dataframe in temporal order: L0,Y, Z,L,M,...,Y,D
                    var.type, # Vector of variable specifications for data.  Fi=fixed (e.g. the exposure), L=confounder, M= mediator, Y=outcome, D = event indicator. 
                    BModels, #fitted BART model
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
  
  
  P_zzstar_mat <- matrix(nrow =  n_Y , ncol = Ndraws)
  P_zzstar_mat[1,] <- 0                              #assuming everyone is alive at visit 1
  
  
  
  
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
    if(var.type[n_L0+2] == "Fi") {
      x_trt <- data.frame(rep(fixed.regime.trt[1], J))
      x_control <- data.frame(rep(fixed.regime.control[1], J))
      l <- l + 1
    } else {
      x_trt <- data.frame(sample(data[, n_L0+2], size = J, replace = TRUE))
      x_control <- data.frame(sample(data[, n_L0+2], size = J, replace = TRUE))
    }
    x_trt <- cbind(L0samples,x_trt) # Merge baseline confounder samples with the first column of the prediction variable
    x_control <- cbind(L0samples,x_control)
    for (j in (n_L0+3):ncol(data)) {                                         #(numL0Vars+1)
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
      # else if(var.type[j] == "Y") {
      # Yprob_trt_it <- lapply(x_trt, function(x2) {pred_surv(BModels[[j - 1]], data.frame(x2), it)})
      # P_zzstar_mat[m,it] <- mean(unlist(Yprob_trt_it))
      # m <- m+1
      # x_trt <- lapply(x_trt, function(x2) {pred_comb(continuous[j], BModels[[j - 1]], data.frame(x2), it)})
      # x_control <- lapply(x_control, function(x2) {pred_comb(continuous[j], BModels[[j - 1]], data.frame(x2), it)})
      # }
      else if(var.type[j] == "Y") {
        Yprob_trt_it <- pred_surv(BModels[[j - 1]], x_trt, it) #failure probablity
        Yprob_control_it <- pred_surv(BModels[[j - 1]], x_control, it) #failure probablity
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
        x_trt <- cbind(x_trt, surv_trt)
        x_trt <- x_trt[surv_trt==1,]
        x_control <- cbind(x_control, surv_control)
        x_control <- x_control[surv_trt==1,]
        s_hat <- s_temp_trt[surv_trt==1]
      }
    }
  }
  
  return(P_zzstar_mat)
}











#-------------------------------------------------------------------------------------------------------------------#
#--------------------------------------Store Results----------------------------------------------------------------#

############################################################
####....... FUNCTION CALL (NO COMPETING EVENT).....#########







#generated data

input_data1 = data.frame(L0_df, Y1 = Y_1, Z1 = Z_1, L1 = L_1, M1 = M_1,
                         Y2 = Y_2, Z2 = Z_2, L2 = L_2, M2 = M_2, 
                         Y3 = Y_3, delta = delta_generated)

input_varType1 = c("L0","L0","L0","L0","L0","L0","Y","Fi", "L", "M","Y", "Fi", "L", "M", "Y","D")


#################################################################################
#................................FIT BART MODELS................................#
BartModel = FitBART(data = input_data1,
                    var.type = input_varType1,
                    J=10000, 
                    Ndraws= 1000,         # Actually runs 'Ndraws * Keepevery + Nskip' iterations but keeps only 'Ndraws' iterations
                    Nskip = 1000, 
                    Ntree = 100, 
                    Keepevery = 10, 
                    Suppress = TRUE, 
                    By = Ndraws/2, 
                    Sparse = TRUE)


#################################################################################
#...................NO CONDITION ON BASELINE AGE................................#

start_time1 <- Sys.time()
P_zzstar = MCInteg(data = input_data1,
                   var.type = input_varType1,
                   #var.type.temp.order = input_varType2,
                   BModels = BartModel,
                   fixed.regime.trt = rep(1,2),
                   fixed.regime.control = rep(0,2),
                   #conditionOnBaselineAge = TRUE, 
                   #startBaselineAge = 45, 
                   #endBaselineAge = 48,
                   J=10000,          #10000
                   Ndraws= 1000,                  #set to 1 because we are computing the ground truth
                   #Nskip = 1000, 
                   #Ntree = 100, 
                   #Keepevery = 5, 
                   Suppress = TRUE, 
                   By = Ndraws/2, 
                   Sparse = TRUE)

end_time1 <- Sys.time()
cat("G-comp runtime is :",end_time1 - start_time1, ".\n")


#print(P_zzstar)


S_zzstar = get_S_zzstar(P_zzstar)
#print(S_zzstar)
#conf_int1 = calculate_95confidence_interval(S_zzstar) 
#print(conf_int1)
#plot_ConfInt(conf_int1)
#cred_int1 = calculate_95bayesian_credible_interval(S_zzstar) 
#print(cred_int1)
#plot_ConfInt(cred_int1)


P_zz = MCInteg(data = input_data1,
               var.type = input_varType1,
               #var.type.temp.order = input_varType2,
               BModels = BartModel,
               fixed.regime.trt = rep(1,2),
               fixed.regime.control = rep(1,2),
               #conditionOnBaselineAge = TRUE, 
               #startBaselineAge = 45, 
               #endBaselineAge = 48,
               J=10000, 
               Ndraws= 1000, 
               #Nskip = 1000, 
               #Ntree = 100, 
               #Keepevery = 5,   
               Suppress = TRUE, 
               By = Ndraws/2, 
               Sparse = TRUE)




#print(P_zz)


S_zz = get_S_zzstar(P_zz)
#print(S_zz)
# conf_int2 = calculate_95confidence_interval(S_zz) 
# print(conf_int2)
# plot_ConfInt(conf_int2)
# cred_int2 = calculate_95bayesian_credible_interval(S_zz) 
# print(cred_int2)
# plot_ConfInt(cred_int2)


P_zstarzstar = MCInteg(data = input_data1,
                       var.type = input_varType1,
                       #var.type.temp.order = input_varType2,
                       BModels = BartModel,
                       fixed.regime.trt = rep(0,2),
                       fixed.regime.control = rep(0,2),
                       #conditionOnBaselineAge = TRUE, 
                       #startBaselineAge = 45, 
                       #endBaselineAge = 48,
                       J=10000, 
                       Ndraws= 1000, 
                       #Nskip = 1000, 
                       #Ntree = 100, 
                       #Keepevery = 5,  
                       Suppress = TRUE, 
                       By = Ndraws/2, 
                       Sparse = TRUE)




#print(P_zstarzstar)


S_zstarzstar = get_S_zzstar(P_zstarzstar)
#print(S_zstarzstar)
# conf_int3 = calculate_95confidence_interval(S_zstarzstar) 
# print(conf_int3)
# plot_ConfInt(conf_int3)
#cred_int3 = calculate_95bayesian_credible_interval(S_zstarzstar) 
# print(cred_int3)
# plot_ConfInt(cred_int3)

IDE = get_IDE(S_zzstar,S_zstarzstar)
# proportionsDE = apply(IDE > 0, 1, mean)
# cat("P(IDE > 0) at each visit is:", proportionsDE, ".\n")
# print(rowMeans(IDE))
#conf_int_IDE = calculate_95confidence_interval(IDE) 
# print(conf_int4)
# plot_ConfInt(conf_int4)
cred_int_IDE = calculate_95bayesian_credible_interval(IDE) 
# print(cred_int_IDE)
# plot_ConfInt(cred_int_IDE)



IIE = get_IIE(S_zz,S_zzstar)
# proportionsIE = apply(IIE > 0, 1, mean)
# cat("P(IIE > 0) at each visit is:", proportionsIE, ".\n")
# print(rowMeans(IIE))
#conf_int_IIE = calculate_95confidence_interval(IIE) 
# print(conf_int5)
# plot_ConfInt(conf_int5)
cred_int_IIE = calculate_95bayesian_credible_interval(IIE) 
# print(cred_int_IIE)
# plot_ConfInt(cred_int_IIE)




Total_Effects = get_TotalEffects(IDE,IIE)

# proportions <- apply(Total_Effects > 0, 1, mean)
# cat("P(TE > 0) at each visit is:", proportions, ".\n")
# print(rowMeans(Total_Effects))
#conf_int_TE = calculate_95confidence_interval(Total_Effects) 
# print(conf_int6)
# plot_ConfInt(conf_int6)
cred_int_TE = calculate_95bayesian_credible_interval(Total_Effects) 
# print(cred_int_TE)
# plot_ConfInt(cred_int_TE)










#-------------------------------------------------------------------------------------------------------------------#
#--------------------------------WRITE RESULTS ON A DF---------------------__---------------------------------------#


allinfo = data.frame( runID = run_ID,
                      IDE2_post_mean = as.numeric(format(round(rowMeans(IDE)[2] , 4), nsmall = 4)),
                      IDE2_lowerCI = cred_int_IDE[2,1],
                      IDE2_upperCI = cred_int_IDE[2,2],
                      IIE2_post_mean = as.numeric(format(round(rowMeans(IIE)[2] , 4), nsmall = 4)),
                      IIE2_lowerCI = cred_int_IIE[2,1],
                      IIE2_upperCI = cred_int_IIE[2,2],
                      TE2_post_mean = as.numeric(format(round(rowMeans(Total_Effects)[2] , 4), nsmall = 4)),
                      TE2_lowerCI = cred_int_TE[2,1],
                      TE2_upperCI = cred_int_TE[2,2],
                      IDE3_post_mean = as.numeric(format(round(rowMeans(IDE)[3] , 4), nsmall = 4)),
                      IDE3_lowerCI = cred_int_IDE[3,1],
                      IDE3_upperCI = cred_int_IDE[3,2],
                      IIE3_post_mean = as.numeric(format(round(rowMeans(IIE)[3] , 4), nsmall = 4)),
                      IIE3_lowerCI = cred_int_IIE[3,1],
                      IIE3_upperCI = cred_int_IIE[3,2],
                      TE3_post_mean = as.numeric(format(round(rowMeans(Total_Effects)[3] , 4), nsmall = 4)),
                      TE3_lowerCI = cred_int_TE[3,1],
                      TE3_upperCI = cred_int_TE[3,2] )
write.table(allinfo, file = txt.title, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
#save(allinfo, file = txt.title)
rm(list="IDE")
rm(list="IIE")
rm(list="Total_Effects")
rm(list="cred_int_IDE")
rm(list="cred_int_IIE")
rm(list="cred_int_TE")
gc()





#-------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------#













#-------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------#