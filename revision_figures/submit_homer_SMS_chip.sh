#!/bin/bash
  
#SBATCH --cpus-per-task=1
#SBATCH --job-name=bash

#call motifs on peaks in ATAC-seq peaks from SMS-CTR-iP3F, used in S12C

ml purge

ml load homer/4.11.1 rstudio

#prefixes for categories of peaks to call motifs on, based on categories from peak_set_overlaps.R

prefixList=(SMS SMS-0 SMS-75 SMS-250 SMS-500)

for(( f=0; f<${#prefixList[@]}; f++ ))

	do 

	first_path=/home/gdstantonlab/rah013/RD_SMS_iP3F/ATAC_seq/peak_categories_chip/${prefixList[$f]}_peaks_040926.bed

    #add peak id to each peak, needed for homer
    Rscript ~/scripts/add_peak_id_bed_centered.R $first_path
   
	out_name=/home/gdstantonlab/rah013/RD_SMS_iP3F/ATAC_seq/peak_categories_chip/${prefixList[$f]}
    
	#call motifs
    findMotifsGenome.pl /home/gdstantonlab/rah013/RD_SMS_iP3F/ATAC_seq/peak_categories_chip/homer_${prefixList[$f]}_peaks_040926.bed hg38 $out_name -size 200 -preparsedDir /home/gdstantonlab/rah013/RD_SMS_iP3F/ATAC_seq/peak_categories_chip/${prefixList2[0]} 

	done


ml purge

