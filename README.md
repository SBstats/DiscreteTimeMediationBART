# DiscreteTimeMediationBART

This repository contains code to implement methods described in the article "A Bayesian semi-parametric approach to causal mediation for longitudinal mediators and time-to-event outcomes with application to a cardiovascular disease cohort study." A brief description of the R files is given below:



## Real data analysis codes

- The files **DiscTimeMediationBARTNoCompEvent.R** and **DiscTimeMediationBARTCompEvent.R** contain code to implement the methods described in the paper for settings without and with a competing event, respectively. Specifically, these files include functions to fit BART models and perform G-computation.

- Similarly, the files **No Comp Event Function Calls.R** and **Comp Event Function Calls.R** contain code that makes function calls to the respective implementation files: **DiscTimeMediationBARTNoCompEvent.R** and **DiscTimeMediationBARTCompEvent.R**. 

To run the code on real data, users must first import their datasets in the function call files and update the variable names at the top of each file accordingly.




## Simulation study codes

- The R files used for simulation have filenames beginning with **Sim-**.

- The files **Sim-Fitted BART Models For True Data Generating Models.R** and **Sim-Fitted Parametric Models For True Data Generating Models.R** contain code for running the simulation on **one data replication** generated using the true data-generating models. To run these files, users must first import their datasets and update the variable names at the top of each file accordingly.

- We executed these simulations on a supercomputer by submitting batch files to conduct the studies in parallel across replicated datasets. The output for all replications is saved in a single `.txt` file after the code is executed.

- This `.txt` file is then imported by **Sim-TRUE MODELS Compute Ground Truth AND PRINT BIAS MSE 95 COVERAGE.R**, which computes the bias, mean squared error (MSE), and 95% coverage.

- For data replications generated under misspecified models, the corresponding files are:
  - **Sim-Fitted BART Models For Misspecified Data Generating Models.R**
  - **Sim-Fitted Parametric Models For Misspecified Data Generating Models.R**
  - **Sim-MISSPECIFIED MODELS Compute Ground Truth AND PRINT BIAS MSE 95 COVERAGE.R**

These files are used to run simulations and compute the bias, MSE, and 95% coverage under the misspecified model setting.
