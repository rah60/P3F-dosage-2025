#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=bash

#call motifs for both Figure 5D & Extended Data Figure 4C

ml purge

ml load homer/4.11.1
ml rstudio

file_dir="/home/gdstantonlab/rah013/repliATAC"

output_dir="/home/gdstantonlab/rah013/repliATAC"

prefixList=(0p_peaks_103024_pulse 75p_peaks_103024_pulse 500p_peaks_103024_pulse 0c_peaks_103024_chase 75c_peaks_103024_chase 500c_peaks_103024_chase 0p_75p_500p_peaks_103024_pulse 0c_75c_500c_peaks_103024_chase)

for(( f=0; f<${#prefixList[@]}; f++ ))
  
  do 
  
  out_name=/home/gdstantonlab/rah013/repliATAC/homer_results/${prefixList[$f]}
  
  findMotifsGenome.pl ${file_dir}/${prefixList[$f]}.bed hg38 $out_name -size 200 -preparsedDir /home/gdstantonlab/rah013/repliATAC/${prefixList[0]}
  
  done
  
  
  ml purge