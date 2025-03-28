#!/bin/bash

#SBATCH --job-name=bash
#SBATCH --time=01-00:00:00

#takes output from ENCODE ChIP-seq pipeline and extends peak size (prior to looking for overlapping peaks and motif calling downstream)

ml purge 

ml load GCC/10.3.0 BEDTools/2.30.0

#read in bed files

file_dir="/home/gdstantonlab/lab/GSL-RH-3507"

output_dir="/home/gdstantonlab/rah013/dbt_dosage_two_reps"

samples=(DBT-MYCN-IP3F_DOSAGE0_2w DBT-MYCN-IP3F_DOSAGE75_2w DBT-MYCN-IP3F_DOSAGE150_2w DBT-MYCN-IP3F_DOSAGE250_2w DBT-MYCN-IP3F_DOSAGE500_2w DBT-MYCN-IP3F_DOSAGE1000_2w)

samples_less=(0_2w 75_2w 150_2w 250_2w 500_2w 1000_2w)

#for each file, get the file name prefix and add 50 bp to each side, save with prefix as new bed file

for(( f=0; f<${#samples[@]}; f++ ))
	do 

	file_name=${file_dir}/${samples[$f]}/${samples[$f]}_IDR_Conservative.bed

	output_name=${output_dir}/${samples_less[$f]}_50bp_peaks.bed
	
	bedtools slop -i $file_name -g ~/chip_p3f_dosage/genome/genome.fa.fai -b 50 > $output_name 

	done