#!/bin/bash
  
#SBATCH --cpus-per-task=1
#SBATCH --job-name=bash
#SBATCH --partition=himem

#set up annotations for Extended Data figure 4B, peak_annotation_fig.R

ml purge

ml load homer/4.11.1

annotatePeaks.pl /home/gdstantonlab/rah013/repliATAC/pulse_GR_103024.bed hg38 >  /home/gdstantonlab/rah013/repliATAC/pulse_GR_103024_annotated.txt

annotatePeaks.pl /home/gdstantonlab/rah013/repliATAC/chase_GR_103024.bed hg38 >  /home/gdstantonlab/rah013/repliATAC/chase_GR_103024_annotated.txt

ml purge
