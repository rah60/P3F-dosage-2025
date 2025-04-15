#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=bash

#call motifs for both Figures 6C & Extended Data Figure 5C

ml purge

ml load homer/4.11.1
ml rstudio

file_dir="/home/gdstantonlab/rah013/repliATAC"

output_dir="/home/gdstantonlab/rah013/repliATAC"

prefixList=(cluster_1_repliATAC_chase cluster_2_repliATAC_chase cluster_3_repliATAC_chase)

for(( f=0; f<${#prefixList[@]}; f++ ))
  
  do 
  
  out_name=/home/gdstantonlab/rah013/repliATAC/homer_results/${prefixList[$f]}
  
  findMotifsGenome.pl ${file_dir}/${prefixList[$f]}.bed hg38 $out_name -size 200 -preparsedDir /home/gdstantonlab/rah013/repliATAC/${prefixList[0]}
  
  done
  
  
  ml purge