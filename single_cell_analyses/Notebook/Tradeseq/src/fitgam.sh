#!/bin/bash
#SBATCH --time=30:00:00
#SBATCH --mail-user=ambuj.kumar@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE
#SBATCH --partition=himem
#SBATCH --mem=600G    
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --account=gdstantonlab
#SBATCH --job-name=phase1_fitgam
#SBATCH --output=phase2_fitgam.out        
#SBATCH --error=phase2_fitgam.err  


set -e

cd /home/gdwanglab/axk201/BenLab/BenLab_PAX3_FOXO1_SingleCellAnalysis/Notebook/Tradeseq

ml GCC/9.3.0
ml OpenMPI/4.0.3
ml R/4.4.1

Rscript ./src/fitgam_condition_phase2.r