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




source("pred_comb_mediator.R")
source("pred_surv.R")
source("pwbart_it.R")
source("expit.R")
source("quiet.R")
source("pred_comb_Y.R")
source("pred_comb_mediator.R")
source("pred_comb.R")
source("pred_y.R")
source("DiscTimeMediationBARTCompEvent.R")




############################################################
####....... FUNCTION CALLS (COMPETING EVENT).....#########


# ARIC_workdf_wide is the analysis dataset in the WIDE FORMAT


input_data2 = data.frame(V0 = ARIC_workdf_wide$SMOKE_STATUS_MULT_1, AGE = ARIC_workdf_wide$AGE_1, 
                         BMI = ARIC_workdf_wide$BMI_1, RACE = ARIC_workdf_wide$RACE,
                         SEX = ARIC_workdf_wide$SEX, EDU = ARIC_workdf_wide$EDU,
                         DIAB_BASELINE = ARIC_workdf_wide$diab_baseline,
                         TOT_HDL_CHL_RATIO_base = ARIC_workdf_wide$TOT_HDL_CHL_RATIO_1,
                         Y11 = ARIC_workdf_wide$Y11, Y21 = ARIC_workdf_wide$Y21,
                         Z1 = ARIC_workdf_wide$RXHYP_1, L1 = ARIC_workdf_wide$SMOKE_STATUS_BIN_1,
                         M1 = ARIC_workdf_wide$MeanBP_1, Y12 = ARIC_workdf_wide$Y12, Y22 = ARIC_workdf_wide$Y22,
                         Z2 = ARIC_workdf_wide$RXHYP_2, L2 = ARIC_workdf_wide$SMOKE_STATUS_BIN_2, 
                         M2 = ARIC_workdf_wide$MeanBP_2,Y13 = ARIC_workdf_wide$Y13, Y23 = ARIC_workdf_wide$Y23,
                         Z3 = ARIC_workdf_wide$RXHYP_3, L3 = ARIC_workdf_wide$SMOKE_STATUS_BIN_3, 
                         M3 = ARIC_workdf_wide$MeanBP_3, Y14 = ARIC_workdf_wide$Y14, Y24 = ARIC_workdf_wide$Y24, 
                         d = ARIC_workdf_wide$d)

input_varType2 = c("L0","L0","L0","L0","L0","L0","L0","L0","Y1","Y2","Fi", "L", "M", "Y1","Y2","Fi", "L", "M", 
                   "Y1","Y2","Fi", "L", "M", "Y1","Y2","d")


#################################################################################
#................................FIT BART MODELS................................#


BartModelComp1 = FitBARTCompEvent(data = input_data2,
                                  var.type = input_varType2,
                                  k =1,
                                  J=10000, 
                                  Ndraws= 1000, 
                                  Nskip = 1000, 
                                  Ntree = 100, 
                                  Keepevery = 10, 
                                  Suppress = TRUE, 
                                  By = Ndraws/10, 
                                  Sparse = TRUE) 




BartModelComp2 = FitBARTCompEvent(data = input_data2,
                                  var.type = input_varType2,
                                  k =2,
                                  J=10000, 
                                  Ndraws= 1000, 
                                  Nskip = 1000, 
                                  Ntree = 100, 
                                  Keepevery = 10,  
                                  Suppress = TRUE, 
                                  By = Ndraws/10, 
                                  Sparse = TRUE) 






#################################################################################################
#...................G-Computatation: NO CONDITION ON BASELINE AGE................................#
#################################################################################################


P1_zzstar = MCIntegCompEvent(data = input_data2,
                             var.type = input_varType2,
                             BModels = BartModelComp1, 
                             k = 1, #k = event type (1 vs 2)
                             fixed.regime.trt = rep(1,3),
                             fixed.regime.control = rep(0,3),
                             #conditionOnBaselineAge = TRUE, 
                             #startBaselineAge = 45, 
                             #endBaselineAge = 48,
                             J=10000, 
                             Ndraws=1000, 
                             #Nskip = 1000, 
                             #Ntree = 100, 
                             #Keepevery = 5, 
                             Suppress = TRUE, 
                             By = Ndraws/10, 
                             Sparse = TRUE) 

#print(P1_zzstar)







P2_zzstar = MCIntegCompEvent(data = input_data2,
                             var.type = input_varType2,
                             BModels = BartModelComp2,
                             k = 2, #k = event type (1 vs 2)
                             fixed.regime.trt = rep(1,3),
                             fixed.regime.control = rep(0,3),
                             #conditionOnBaselineAge = TRUE, 
                             #startBaselineAge = 45, 
                             #endBaselineAge = 48,
                             J=10000, 
                             Ndraws= 1000, 
                             #Nskip = 1000, 
                             #Ntree = 100, 
                             #Keepevery = 5, 
                             Suppress = TRUE, 
                             By = Ndraws/10, 
                             Sparse = TRUE) 

#print(P2_zzstar)


F_1_zzstar = get_F_k_zzstar(P1_zzstar,P2_zzstar,1)
#print(F_1_zzstar)
conf_int1 = calculate_95confidence_interval(F_1_zzstar) 
#print(conf_int1)
plot_ConfInt(conf_int1)
cred_int1 = calculate_95bayesian_credible_interval(F_1_zzstar) 
#print(cred_int1)
plot_ConfInt(cred_int1)


F_2_zzstar =  get_F_k_zzstar(P1_zzstar,P2_zzstar,2)
#print(F_2_zzstar)
conf_int2 = calculate_95confidence_interval(F_2_zzstar) 
#print(conf_int2)
plot_ConfInt(conf_int2)
cred_int2 = calculate_95bayesian_credible_interval(F_2_zzstar) 
#print(cred_int2)
plot_ConfInt(cred_int2)







P1_zz = MCIntegCompEvent(data = input_data2,
                         var.type = input_varType2,
                         BModels = BartModelComp1, 
                         k = 1, #k = event type (1 vs 2)
                         fixed.regime.trt = rep(1,3),
                         fixed.regime.control = rep(1,3),
                         #conditionOnBaselineAge = TRUE, 
                         #startBaselineAge = 45, 
                         #endBaselineAge = 48,
                         J=10000, 
                         Ndraws= 1000, 
                         #Nskip = 1000, 
                         #Ntree = 100, 
                         #Keepevery = 5,  
                         Suppress = TRUE, 
                         By = Ndraws/10, 
                         Sparse = TRUE) 

#print(P1_zz)







P2_zz = MCIntegCompEvent(data = input_data2,
                         var.type = input_varType2,
                         BModels = BartModelComp2,
                         k = 2, #k = event type (1 vs 2)
                         fixed.regime.trt = rep(1,3),
                         fixed.regime.control = rep(1,3),
                         #conditionOnBaselineAge = TRUE, 
                         #startBaselineAge = 45, 
                         #endBaselineAge = 48,
                         J=10000, 
                         Ndraws= 1000, 
                         #Nskip = 1000, 
                         #Ntree = 100, 
                         #Keepevery = 5,  
                         Suppress = TRUE, 
                         By = Ndraws/10, 
                         Sparse = TRUE) 

#print(P2_zz)


F_1_zz = get_F_k_zzstar(P1_zz,P2_zz,1)
#print(F_1_zz)
conf_int3 = calculate_95confidence_interval(F_1_zz) 
#print(conf_int3)
plot_ConfInt(conf_int3)
cred_int3 = calculate_95bayesian_credible_interval(F_1_zz) 
#print(cred_int3)
plot_ConfInt(cred_int3)


F_2_zz =  get_F_k_zzstar(P1_zz,P2_zz,2)
#print(F_2_zz)
conf_int4 = calculate_95confidence_interval(F_2_zz) 
#print(conf_int4)
plot_ConfInt(conf_int4)
cred_int4 = calculate_95bayesian_credible_interval(F_2_zz) 
#print(cred_int4)
plot_ConfInt(cred_int4)







P1_zstarzstar = MCIntegCompEvent(data = input_data2,
                                 var.type = input_varType2,
                                 BModels = BartModelComp1, 
                                 k = 1, #k = event type (1 vs 2)
                                 fixed.regime.trt = rep(0,3),
                                 fixed.regime.control = rep(0,3),
                                 #conditionOnBaselineAge = TRUE, 
                                 #startBaselineAge = 45, 
                                 #endBaselineAge = 48,
                                 J=10000, 
                                 Ndraws= 1000, 
                                 #Nskip = 1000, 
                                 #Ntree = 100, 
                                 #Keepevery = 5,   
                                 Suppress = TRUE, 
                                 By = Ndraws/10, 
                                 Sparse = TRUE) 

#print(P1_zstarzstar)







P2_zstarzstar = MCIntegCompEvent(data = input_data2,
                                 var.type = input_varType2,
                                 BModels = BartModelComp2,
                                 k = 2, #k = event type (1 vs 2)
                                 fixed.regime.trt = rep(0,3),
                                 fixed.regime.control = rep(0,3),
                                 #conditionOnBaselineAge = TRUE, 
                                 #startBaselineAge = 45, 
                                 #endBaselineAge = 48,
                                 J=10000, 
                                 Ndraws= 1000, 
                                 #Nskip = 1000, 
                                 #Ntree = 100, 
                                 #Keepevery = 5,  
                                 Suppress = TRUE, 
                                 By = Ndraws/10, 
                                 Sparse = TRUE) 

#print(P2_zstarzstar)


F_1_zstarzstar = get_F_k_zzstar(P1_zstarzstar,P2_zstarzstar,1)
#print(F_1_zstarzstar)
conf_int5 = calculate_95confidence_interval(F_1_zstarzstar) 
#print(conf_int5)
plot_ConfInt(conf_int5)
cred_int5 = calculate_95bayesian_credible_interval(F_1_zstarzstar) 
#print(cred_int5)
plot_ConfInt(cred_int5)


F_2_zstarzstar =  get_F_k_zzstar(P1_zstarzstar,P2_zstarzstar,2)
#print(F_2_zstarzstar)
conf_int6 = calculate_95confidence_interval(F_2_zstarzstar) 
#print(conf_int6)
plot_ConfInt(conf_int6)
cred_int6 = calculate_95bayesian_credible_interval(F_2_zstarzstar) 
#print(cred_int6)
plot_ConfInt(cred_int6)




#########################################################################################################
#...................Compute Causal Effects................................#
#########################################################################################################


IDE1 = get_IDE(F_1_zzstar,F_1_zstarzstar)
#print(rowMeans(IDE1))
conf_int7 = calculate_95confidence_interval(IDE1) 
#print(conf_int7)
plot_ConfInt(conf_int7)
cred_int7 = calculate_95bayesian_credible_interval(IDE1) 
#print(cred_int7)
plot_ConfInt(cred_int7)

IIE1 = get_IIE(F_1_zz,F_1_zzstar)
#print(rowMeans(IIE1))
conf_int8 = calculate_95confidence_interval(IIE1) 
#print(conf_int8)
plot_ConfInt(conf_int8)
cred_int8 = calculate_95bayesian_credible_interval(IIE1) 
#print(cred_int8)
plot_ConfInt(cred_int8)


Total_Effects1 = get_TotalEffects(IDE1,IIE1)
#print(rowMeans(Total_Effects1))
conf_int9 = calculate_95confidence_interval(Total_Effects1) 
#print(conf_int9)
plot_ConfInt(conf_int9)
cred_int9 = calculate_95bayesian_credible_interval(Total_Effects1) 
#print(cred_int9)
plot_ConfInt(cred_int9)



IDE2 = get_IDE(F_2_zzstar,F_2_zstarzstar)
#print(rowMeans(IDE2))
conf_int10 = calculate_95confidence_interval(IDE2) 
#print(conf_int10)
plot_ConfInt(conf_int10)
cred_int10 = calculate_95bayesian_credible_interval(IDE2) 
#print(cred_int10)
plot_ConfInt(cred_int10)

IIE2 = get_IIE(F_2_zz,F_2_zzstar)
#print(rowMeans(IIE2))
conf_int11 = calculate_95confidence_interval(IIE2) 
#print(conf_int11)
plot_ConfInt(conf_int11)
cred_int11 = calculate_95bayesian_credible_interval(IIE2) 
#print(cred_int11)
plot_ConfInt(cred_int11)




Total_Effects2 = get_TotalEffects(IDE2,IIE2)
#print(rowMeans(Total_Effects2))
conf_int12 = calculate_95confidence_interval(Total_Effects2) 
#print(conf_int12)
plot_ConfInt(conf_int12)
cred_int12 = calculate_95bayesian_credible_interval(Total_Effects2) 
#print(cred_int12)
plot_ConfInt(cred_int12)





