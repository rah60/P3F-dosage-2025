#!/bin/bash
  
#SBATCH --cpus-per-task=1
#SBATCH --job-name=bash

#call motifs on peaks and annotate peaks with nearest gene, genomic features

ml purge

ml load homer/4.11.1 rstudio

#prefixes for categories of peaks to call motifs on, based on categories from peak_set_overlaps.R
prefixList=(250 500 150_500 75_150_250_500 75_150_250_500_1000 150 75_150_500 150_250 75_150_250_500_1000_0 150_250_500 250_500 1000 75_150_250 75 75_150 75_250 75_150_1000 250_500_1000 75_150_250_1000 75_250_500 75_150_250_500_0 150_250_500_1000 500_1000 75_500 75_250_500_1000 250_1000 150_500_1000 75_150_250_0 150_1000 150_250_1000 0 150_250_500_1000_0 75_500_0 75_0 75_150_500_1000 75_150_0 75_1000 75_150_500_0 150_250_500_0)

for(( f=0; f<${#prefixList[@]}; f++ ))

	do 

	first_path=/home/gdstantonlab/rah013/dbt_dosage_two_reps/peak_categories/with_zero/${prefixList[$f]}_peaks_2w_090624.bed

    #add peak id to each peak, needed for homer
    Rscript ~/scripts/add_peak_id_bed_centered.R $first_path
   
	out_name=/home/gdstantonlab/rah013/dbt_dosage_two_reps/peak_categories/with_zero/homer_results/${prefixList[$f]}
    
	#call motifs
    findMotifsGenome.pl /home/gdstantonlab/rah013/dbt_dosage_two_reps/peak_categories/with_zero/homer_${prefixList[$f]}_peaks_2w_090624.bed hg38 $out_name -size 200 -preparsedDir /home/gdstantonlab/rah013/dbt_dosage_two_reps/${prefixList2[0]} 

	done


#total peaks in each dosage 
prefixList=(0_2w_50bp_peaks 75_2w_50bp_peaks 150_2w_50bp_peaks 250_2w_50bp_peaks 500_2w_50bp_peaks 1000_2w_50bp_peaks 0_2w_unique_intersect 75_2w_unique_intersect 150_2w_unique_intersect 250_2w_unique_intersect 500_2w_unique_intersect 1000_2w_unique_intersect 75_150_250_500_common_intersect 75_150_250_500_1000_common_intersect)

for(( f=0; f<${#prefixList[@]}; f++ ))

	do 

	out_name=/home/gdstantonlab/rah013/dbt_dosage_two_reps/homer_results/${prefixList[$f]}_annotated.txt

	#annotate peaks--nearest gene & type of region in genome
	annotatePeaks.pl /home/gdstantonlab/rah013/dbt_dosage_two_reps/${prefixList[$f]}.bed hg38 > $out_name  

	done




ml purge

