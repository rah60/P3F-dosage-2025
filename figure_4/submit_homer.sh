#!/bin/bash
  
#SBATCH --job-name=bash

#motif calling from clusters established in rna_chip_figures.R
#motif calls used in motif_freq_plot_kmeans.R and new_motif_plot_kmeans_clusters.R

ml purge

ml load homer/4.11.1 rstudio 

prefixList=(1 2 3 4 5 6)

for(( f=0; f<${#prefixList[@]}; f++ ))
 
	do 

	first_path=/home/gdstantonlab/rah013/dbt_dosage_two_reps/clusters/cluster_${prefixList[$f]}_peaks_102124.bed #bed files established from k-means clusters in rna_chip_figures.R
   
	out_name=/home/gdstantonlab/rah013/dbt_dosage_two_reps/clusters/homer_results/cluster_${prefixList[$f]}
    
    findMotifsGenome.pl $first_path hg38 $out_name -size 200 -preparsedDir /home/gdstantonlab/rah013/dbt_dosage_two_reps/${prefixList2[0]} 

	done


ml purge

