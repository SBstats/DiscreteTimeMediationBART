# DiscreteTimeMediationBART




## Example sbatch file

\#!/bin/bash
\#SBATCH --job-name=ARICBARTSim                          
\#SBATCH --ntasks=1
\#SBATCH --mem=15GB
\#SBATCH --mail-type=ALL                                 
\#SBATCH --mail-user=s.bhandari\@ufl.edu                 
\#SBATCH --time=50:00:00                                 
\#SBATCH --array=1-1000
date; hostname; pwd

mkdir -p param_simu_results

module load R

Rscript Simulation-BARTModelsRun1000.R