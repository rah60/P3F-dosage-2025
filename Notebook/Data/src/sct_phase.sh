#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mail-user=ambuj.kumar@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=himem 
#SBATCH --mem 600G
#SBATCH --account=gdstantonlab
#SBATCH --job-name=sct
#SBATCH --output=sct.out
#SBATCH --error=sct.err

set -euo pipefail

cd /home/gdwanglab/axk201/BenLab/BenLab_PAX3_FOXO1_SingleCellAnalysis/Notebook/Data

ml GCC/9.3.0
ml OpenMPI/4.0.3
ml R/4.4.1

Rscript src/sct_phase.r