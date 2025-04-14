#!/bin/bash
  
#SBATCH --job-name=bash

#peak annotation, for rna_chip_figures.R

ml purge

ml load homer/4.11.1

annotatePeaks.pl /home/gdstantonlab/rah013/dbt_dosage_two_reps/all_peaks_categories.bed hg38 >  /home/gdstantonlab/rah013/dbt_dosage_two_reps/all_peaks_categories_annotated.txt

ml purge


