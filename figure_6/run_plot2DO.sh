#!/bin/bash
  
#SBATCH --cpus-per-task=10
#SBATCH --job-name=bash
#SBATCH --partition=himem

#need plot2DO package to run this script, place script in plot2DO-master directory once downloaded

#Figure 6D & E, Extended Data Figure 5D

ml purge 

ml rstudio

#bed files to use
peakFiles=(cluster_1_repliATAC_pulse.bed cluster_2_repliATAC_pulse.bed cluster_3_repliATAC_pulse.bed)

#max read coverage value
colorVar=(0.15 0.06 0.02)

#bam files for pulse samples
bamFiles=(0_pulse_merged.bam 75_pulse_merged.bam 500_pulse_merged.bam)

#run plot2DO for each bam file and bed file
for(( f=0; f<${#peakFiles[@]}; f++ ))
    do
    var1=$(echo ${peakFiles[$f]} | awk -F '.bed' '{print $1}') #make name of site be the peak file name

    for(( g=0; g<${#bamFiles[@]}; g++ ))
        do
        Rscript ~/plot2DO-master/plot2DO.R  -f ${bamFiles[$g]} \
        -g hg38 -s ${peakFiles[$f]} --siteLabel=$var1 -L 300 -l 30 -u 500 -d 500 --simplifyPlot=on --squeezePlot=on -m ${colorVar[$f]}
        done
    
    done

#repeat above but with bam files for chase samples 
colorVar=(0.35 0.15 0.04)
bamFiles=(0_chase_merged.bam 75_chase_merged.bam 500_chase_merged.bam)

for(( f=0; f<${#peakFiles[@]}; f++ ))
    do
    var1=$(echo ${peakFiles[$f]} | awk -F '.bed' '{print $1}') #make name of site be the peak file name

    for(( g=0; g<${#bamFiles[@]}; g++ ))
        do
        Rscript ~/plot2DO-master/plot2DO.R  -f ${bamFiles[$g]} \
        -g hg38 -s ${peakFiles[$f]} --siteLabel=$var1 -L 300 -l 30 -u 500 -d 500 --simplifyPlot=on --squeezePlot=on -m ${colorVar[$f]}
        done
    
    done

ml purge

