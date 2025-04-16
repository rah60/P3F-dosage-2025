#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --job-name=bash

ml purge 

ml  GCC/9.3.0  OpenMPI/4.0.3 deepTools/3.3.1-Python-3.8.2

#Figure 6A heatmaps

#original heatmap to initially generate k-means clusters in pulse data via heatmap

# computeMatrix reference-point -R ~/dbt_dosage_two_reps/peak_categories_GR_090624.bed \
# -S ~/repliATAC/bigwigs/0_pulse_merged_rpgc.bw ~/repliATAC/bigwigs/75_pulse_merged_rpgc.bw \
# ~/repliATAC/bigwigs/500_pulse_merged_rpgc.bw \
# --referencePoint center -a 1000 -b 1000 -bs 50 \
# -o ~/dosage_manuscript/figure_5/P3F_pulse_110624_rpgc_all_sites.gz \
# -p 10

# plotHeatmap -m ~/dosage_manuscript/figure_5/P3F_pulse_110624_rpgc_all_sites.gz \
# --samplesLabel "0 ng/ml chase" "75 ng/ml chase" "500 ng/ml chase"   --colorMap Blues \
# --missingDataColor "white" \
# --xAxisLabel "distance from peak" \
# --refPointLabel "0" \
# --sortUsingSamples 3 \
# --kmeans 3 \
# --outFileSortedRegions  ~/dosage_manuscript/figure_5/P3F_pulse_heatmap_110624_rpgc_kmeans3_all_sites.bed \
# -o ~/dosage_manuscript/figure_5/P3F_pulse_heatmap_110624_rpgc_kmeans3_all_sites.pdf

#plotting heatmap again with pulse & chase repliATAC data in clusters determined above

computeMatrix reference-point -R ~/repliATAC/cluster_1_repliATAC_pulse.bed ~/repliATAC/cluster_2_repliATAC_pulse.bed ~/repliATAC/cluster_3_repliATAC_pulse.bed \
-S ~/repliATAC/bigwigs/0_pulse_merged_rpgc.bw ~/repliATAC/bigwigs/75_pulse_merged_rpgc.bw \
~/repliATAC/bigwigs/500_pulse_merged_rpgc.bw ~/repliATAC/bigwigs/0_chase_merged_rpgc.bw \
~/repliATAC/bigwigs/75_chase_merged_rpgc.bw ~/repliATAC/bigwigs/500_chase_merged_rpgc.bw \
--referencePoint center -a 1000 -b 1000 -bs 50 \
-o ~/dosage_manuscript/figure_5/P3F_pulse_chase_111824_rpgc_pulse_clusters.gz \
-p 10

plotHeatmap -m ~/dosage_manuscript/figure_5/P3F_pulse_chase_111824_rpgc_pulse_clusters.gz --regionsLabel "cluster 1" "cluster 2" "cluster 3" \
--samplesLabel "0 ng/ml" "75 ng/ml" "500 ng/ml" "0 ng/ml" "75 ng/ml" "500 ng/ml"   --colorMap YlGnBu YlGnBu YlGnBu YlOrRd YlOrRd YlOrRd \
--missingDataColor "white" \
--zMin 0 0 0 0 0 0 \
--zMax 20 20 20 70 70 70 \
--xAxisLabel "distance from peak" \
--whatToShow "heatmap and colorbar" \
--refPointLabel "0" \
--sortUsingSamples 3 \
--heatmapWidth 3 \
--heatmapHeight 15 \
-o ~/dosage_manuscript/figure_5/P3F_pulse_chase_111824_rpgc_pulse_clusters.pdf
