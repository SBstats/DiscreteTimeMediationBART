set.seed(1234)

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
source("DiscTimeMediationBARTNoCompEvent.R")




############################################################
####....... FUNCTION CALLS (NO COMPETING EVENT).....#########


# ARIC_workdf_wide is the analysis dataset in the WIDE FORMAT

input_data1 = data.frame(V0 = ARIC_workdf_wide$SMOKE_STATUS_MULT_1, BMI = ARIC_workdf_wide$BMI_1 ,
                         AGE = ARIC_workdf_wide$AGE_1, RACE = ARIC_workdf_wide$RACE, 
                         SEX = ARIC_workdf_wide$SEX, EDU = ARIC_workdf_wide$EDU, 
                         DIAB_BASELINE = ARIC_workdf_wide$diab_baseline,
                         TOT_HDL_CHL_RATIO_base = ARIC_workdf_wide$TOT_HDL_CHL_RATIO_1,
                         Y1 = ARIC_workdf_wide$Y1, Z1 = ARIC_workdf_wide$RXHYP_1, 
                         L1 = ARIC_workdf_wide$SMOKE_STATUS_BIN_1, M1 = ARIC_workdf_wide$MeanBP_1,  
                         Y2 = ARIC_workdf_wide$Y2, Z2 = ARIC_workdf_wide$RXHYP_2, 
                         L2 = ARIC_workdf_wide$SMOKE_STATUS_BIN_2, M2 = ARIC_workdf_wide$MeanBP_2, 
                         Y3 = ARIC_workdf_wide$Y3, Z3 = ARIC_workdf_wide$RXHYP_3, 
                         L3 = ARIC_workdf_wide$SMOKE_STATUS_BIN_3, M3 = ARIC_workdf_wide$MeanBP_3, 
                         Y4 = ARIC_workdf_wide$Y4, delta = ARIC_workdf_wide$delta)


input_varType1 = c("L0","L0","L0","L0","L0","L0","L0","L0","Y","Fi", "L", "M","Y", "Fi", "L", "M", "Y","Fi", "L", "M", "Y","D")



#################################################################################
#................................FIT BART MODELS................................#
#################################################################################
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





#########################################################################################################
#...................G-Computatation: NO CONDITION ON BASELINE AGE................................#
#########################################################################################################

P_zzstar = MCInteg(data = input_data1,
                   var.type = input_varType1,
                   BModels = BartModel,
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
                   By = Ndraws/2, 
                   Sparse = TRUE)




print(P_zzstar)


S_zzstar = get_S_zzstar(P_zzstar)
#print(S_zzstar)
conf_int1 = calculate_95confidence_interval(S_zzstar) 
#print(conf_int1)
plot_ConfInt(conf_int1)
cred_int1 = calculate_95bayesian_credible_interval(S_zzstar) 
#print(cred_int1)
plot_ConfInt(cred_int1)


P_zz = MCInteg(data = input_data1,
               var.type = input_varType1,
               BModels = BartModel,
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
               By = Ndraws/2, 
               Sparse = TRUE)




#print(P_zz)


S_zz = get_S_zzstar(P_zz)
#print(S_zz)
conf_int2 = calculate_95confidence_interval(S_zz) 
#print(conf_int2)
plot_ConfInt(conf_int2)
cred_int2 = calculate_95bayesian_credible_interval(S_zz) 
#print(cred_int2)
plot_ConfInt(cred_int2)


P_zstarzstar = MCInteg(data = input_data1,
                       var.type = input_varType1,
                       BModels = BartModel,
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
                       By = Ndraws/2, 
                       Sparse = TRUE)




#print(P_zstarzstar)


S_zstarzstar = get_S_zzstar(P_zstarzstar)
#print(S_zstarzstar)
conf_int3 = calculate_95confidence_interval(S_zstarzstar) 
#print(conf_int3)
plot_ConfInt(conf_int3)
cred_int3 = calculate_95bayesian_credible_interval(S_zstarzstar) 
#print(cred_int3)
plot_ConfInt(cred_int3)






#########################################################################################################
#...................Compute Causal Effects................................#
#########################################################################################################

IDE = get_IDE(S_zzstar,S_zstarzstar)
proportionsDE = apply(IDE > 0, 1, mean)
cat("P(IDE > 0) at each visit is:", proportionsDE, ".\n")
#print(rowMeans(IDE))
conf_int4 = calculate_95confidence_interval(IDE) 
#print(conf_int4)
plot_ConfInt(conf_int4)
cred_int4 = calculate_95bayesian_credible_interval(IDE) 
#print(cred_int4)
plot_ConfInt(cred_int4)

IIE = get_IIE(S_zz,S_zzstar)
proportionsIE = apply(IIE > 0, 1, mean)
cat("P(IIE > 0) at each visit is:", proportionsIE, ".\n")
#print(rowMeans(IIE))
conf_int5 = calculate_95confidence_interval(IIE) 
#print(conf_int5)
plot_ConfInt(conf_int5)
cred_int5 = calculate_95bayesian_credible_interval(IIE) 
#print(cred_int5)
plot_ConfInt(cred_int5)




Total_Effects = get_TotalEffects(IDE,IIE)
proportions <- apply(Total_Effects > 0, 1, mean)
cat("P(TE > 0) at each visit is:", proportions, ".\n")
#print(rowMeans(Total_Effects))
conf_int6 = calculate_95confidence_interval(Total_Effects) 
#print(conf_int6)
plot_ConfInt(conf_int6)
cred_int6 = calculate_95bayesian_credible_interval(Total_Effects) 
#print(cred_int6)
plot_ConfInt(cred_int6)







###############################################################





