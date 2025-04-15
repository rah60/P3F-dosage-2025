
rm(list = ls())  
graphics.off()

library(GenomicRanges) 
library(tidyverse) 
library(RColorBrewer)
library(pheatmap)
library(ggpubr)
library(aplot)
library(lemon)
library(cowplot)
library(enrichR)
library(ggsignif)

# Figure 5F, Extended Data Figure 4D, set up for Table S2

#function for plotting Enrichr output
plot_gene_sets <- function(gene_sets, num_sets = 10, fill_col="black", plot_title = ""){

    #gene_sets = output from Enrichr for one specific gene set
    #num_sets = number of sets to plot if more than num_sets are significant
    #fill_col = color for bars on barplot
    #plot_title = plot title

    subset_set <- gene_sets[which(gene_sets$Adjusted.P.value < 0.05),]

    if(nrow(subset_set) == 0){
        return(NULL)
    }
    
    subset_set$log10_adj_pval <- -log10(subset_set$Adjusted.P.value)
    subset_set$log10_pval <- -log10(subset_set$Adjusted.P.value)
    
    #select top num_sets terms & plot
    if(nrow(subset_set) >= num_sets){
        select_20 <- subset_set[order(subset_set$log10_adj_pval, decreasing=T),][1:num_sets,]
    }else{
        select_20 <- subset_set
    }
    
    order.v <- select_20[order(select_20$log10_adj_pval),1]
    
    p1 <- ggplot(select_20) +
        geom_col(aes(x=factor(Term, level=order.v), y=log10_adj_pval, fill="blank"))+
        coord_flip()+
        theme_classic(base_size=16)+
        scale_fill_manual(values=fill_col)+
        scale_y_continuous(expand = c(0,0))+
        labs(x= "Gene sets" , y="-Log10 adjusted p-value", title=plot_title)+
        guides(fill="none")

    return(p1)
}

#Enrichr set up for later
#is connection working? T/F
websiteLive <- getOption("enrichR.live")
options(enrichR.base.address="https://amp.pharm.mssm.edu/Enrichr/")
#View(listEnrichrDbs())

#choose gene set for Enrichr
dbs_to_use <- c("MSigDB_Hallmark_2020")

#PULSE, Extended Data Fig 4D

pulse_peaks <- read.table("~/dosage_manuscript/csv/peak_categories_GR_103024_pulse.tsv",sep="\t")

#add peak ids for homer
peak_ids <- c()
for(i in 1:nrow(pulse_peaks)){
    new_name <- paste0("peak_",i)
    peak_ids <- c(peak_ids , new_name)
}

new_pulse <- cbind(pulse_peaks[,c(1:3)], peak_ids,pulse_peaks[,c(4,5)] )
pulse_peaks$peak_ids <- peak_ids

#table with peak ids for homer
#write.table(new_pulse, "~/repliATAC/pulse_GR_103024.bed", sep="\t", quote=F, row.names=F, col.names=F, eol='\r\n')

#after homer, annotated file
new_pulse <- read.table("~/repliATAC/pulse_GR_103024_annotated.txt", sep="\t", fill=T, header=T, quote="")

colnames(new_pulse)[1] <- "peak_ids"

#merge annotated and original dataframes for peaks so nearest gene info is in dataframe
pulse_merged <- merge(pulse_peaks, new_pulse) 

#compare genes to genes expressed in Dbt (see motif plots for fxn to do this)
rna_seq <- read.table("~/dosage_manuscript/csv/DE_result.txt", sep="\t", as.is=T, header=T, quote="\"")

get_dosage_gene_vector <- function(x, genes_rds = "~/dosage_manuscript/rds/gene_expressed_at_dosages.rds"){ #x is one dosage
  genes_df <- readRDS(genes_rds)
  subset_x <- genes_df[which(genes_df$Dosage == x),]
  subset_x$gene_name
}

#id genes that are expressed at chosen dosages
all_exp_genes <- lapply(c("0","75","150", "250", "500","1000"), get_dosage_gene_vector) %>%
                      unlist() %>%
                      unique()

#subset peaks that are near expressed genes
pulse_2 <- pulse_merged[which(pulse_merged$Nearest.Ensembl %in% all_exp_genes) , ] 

#normalize & scale rna-seq data
rna_seq2 <- rna_seq[, c(8:13)] %>% 
    as.matrix() 
rna_seq2 <- log2(rna_seq2 + 1)

rna_seq2 <- t(scale(t(rna_seq2)))

rna_seq[,c(8:13)] <- rna_seq2

#add expression data to peak dataframe based on nearest gene
pulse_2 <- left_join(pulse_2, rna_seq, by=join_by(Nearest.Ensembl == ensembl))

#use categories for peaks to determine if peaks are P3F independent or dependent for their accessibility
P3F_independent <- c("0p","0p_75p","0p_75p_500p","0p_500p")
P3F_dependent <- c("75p","75p_500p","500p")

pulse_2$p3f_status <- NA

pulse_2[which(pulse_2$category_string %in% P3F_dependent), 98 ] <- "P3F_dependent"

pulse_2[which(pulse_2$category_string %in% P3F_independent), 98 ] <- "P3F_independent"

independent <- pulse_2[which(pulse_2$p3f_status == "P3F_independent"),]

dependent <- pulse_2[which(pulse_2$p3f_status == "P3F_dependent"),]

delete_1 <- independent[which( (independent$gene_name %in% dependent$gene_name) ),]$peak_ids #remove overlap only from independent--one dependent peak = dependent gene

keep_these_indices <- which( !(pulse_2$peak_ids %in% delete_1) )

pulse_2 <- pulse_2[keep_these_indices,]

#need to unique gene names...
pulse_2$duplicate_row <- FALSE
for(w in 1:nrow(pulse_2)){

    if(pulse_2[w, ncol(pulse_2) ]){
        next
    }else{

    these_ones <- which( pulse_2$gene_name == pulse_2[w, 89 ])

    minus_that_one <- these_ones[which(these_ones != w)]

    pulse_2[minus_that_one, ncol(pulse_2)] <- TRUE

    }

}

pulse_2 <- pulse_2[which(!pulse_2$duplicate_row),]

#Table S2
#write.csv(pulse_2, "~/dosage_manuscript/figure_5/pulse_inde_dep_sites.csv", row.names=F)

#prepare for plotting
pulse_2_long <- pulse_2 %>%
                pivot_longer(cols=c(Control,X75ng.ml,X500ng.ml), names_to = "Dosage", values_to="Expression")

p3f_de <- subset(pulse_2_long, p3f_status == "P3F_dependent")
p3f_inde <- subset(pulse_2_long, p3f_status == "P3F_independent")

#determine significance
aov_res <- aov(Expression ~ Dosage, data = p3f_de)
summary(aov(Expression ~ Dosage, data = p3f_de))
TukeyHSD(aov_res)
#all significant
#500, 0 ****
#0, 75 **
#75, 500 ****

aov_res <- aov(Expression ~ Dosage, data = p3f_inde)
summary(aov(Expression ~ Dosage, data = p3f_inde))
TukeyHSD(aov_res)
#0, 500 ***

my_comparisons <- list(c("Control", "X75ng.ml"),c("X500ng.ml", "X75ng.ml"),c("Control", "X500ng.ml"))
####### need to use actual anova/Tukey res & p-values...will need to use a different fxn to plot this

#Extended Data Figure 4D, left

p1 <- ggplot(p3f_de, aes(x=Dosage, y=Expression, fill = p3f_status))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size=16)+
    scale_fill_manual(values=c("#BB5566"))+
    labs(title="PAX3::FOXO1 dependent\naccessible sites", x="Doxycycline dosage (ng/ml)", y="Scaled log2 expression")+
    theme(plot.title = element_text( hjust=0.5),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=16))+
    scale_x_discrete(limits =  c("Control","X75ng.ml","X500ng.ml"), labels=c("0","75","500"))+
    geom_hline(yintercept = 0, linewidth = 1, linetype = 2, color="#DDAA33")+
    geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, annotation=c("**","****","****"), y_position=c(2,2.5,3))

my_comparisons2 <- list(c("Control", "X500ng.ml"))

p2 <- ggplot(p3f_inde, aes(x=Dosage, y=Expression, fill = p3f_status))+
    geom_boxplot(show.legend=F, notch=T, linewidth =1.2)+
    theme_classic(base_size=16)+
    scale_fill_manual(values=c("#004488"))+
    labs(title="PAX3::FOXO1 independent\naccessible sites", x="Doxycycline dosage (ng/ml)", y="")+
    theme(plot.title = element_text( hjust=0.5),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=16))+
    scale_x_discrete(limits =  c("Control","X75ng.ml","X500ng.ml"), labels=c("0","75","500"))+
    geom_hline(yintercept = 0, linewidth = 1, linetype = 2, color="#DDAA33")+
    geom_signif(comparisons = my_comparisons2, map_signif_level = TRUE, annotation=c("***"), y_position=c(2))+
    ylim2(p1)

p3 <- p1 %>% insert_right(p2) 


 png(paste0("~/dosage_manuscript/figure_5/pulse_all_repliATAC_sites_in_dep_exp_revised.png" ), width = 7.5, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
  print(p3)
  dev.off()

enriched_pulse1 <- enrichr(pulse_2[which(pulse_2$p3f_status  == "P3F_dependent"),25], dbs_to_use)

clust_1 <- plot_gene_sets(enriched_pulse1[["MSigDB_Hallmark_2020"]],  plot_title = "P3F dependent", num_sets=5,  fill_col="#BB5566")

enriched_pulse2 <- enrichr(pulse_2[which(pulse_2$p3f_status  == "P3F_independent"),25], dbs_to_use)

clust_2 <- plot_gene_sets(enriched_pulse2[["MSigDB_Hallmark_2020"]],  plot_title = "P3F independent", num_sets=5, fill_col="#004488")

#Extended Data Figure 4B, right

png( paste0("~/dosage_manuscript/figure_5/hallmarks_enrichr_P3F_in_dependent_pulse_plot_revised.png"), width = 5.5, height = 7, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(
plot_grid(clust_1, clust_2,ncol=1)
)
dev.off()

##############CHASE##################
#same as pulse analysis above

#chase
chase_peaks <- read.table("~/dosage_manuscript/csv/peak_categories_GR_103024_chase.tsv",sep="\t")

peak_ids <- c()
for(i in 1:nrow(chase_peaks)){
    new_name <- paste0("peak_",i)
    peak_ids <- c(peak_ids , new_name)
}

new_chase <- cbind(chase_peaks[,c(1:3)], peak_ids,chase_peaks[,c(4,5)] )
chase_peaks$peak_ids <- peak_ids

#write.table(new_chase, "~/repliATAC/chase_GR_103024.bed", sep="\t", quote=F, row.names=F, col.names=F, eol='\r\n')

new_chase <- read.table("~/repliATAC/chase_GR_103024_annotated.txt", sep="\t", fill=T, header=T, quote="")

colnames(new_chase)[1] <- "peak_ids"

chase_merged <- merge(chase_peaks, new_chase) 
#write.csv(chase_merged, "~/repliATAC/chase_test.csv")

#compare genes to genes expressed in Dbt (see motif plots for fxn to do this)
rna_seq <- read.table("~/dosage_manuscript/csv/DE_result.txt", sep="\t", as.is=T, header=T, quote="\"")

genes_df <- readRDS("~/dosage_manuscript/rds/gene_expressed_at_dosages.rds")

get_dosage_gene_vector <- function(x){ #x is one dosage
  subset_x <- genes_df[which(genes_df$Dosage == x),]
  subset_x$gene_name
}

all_exp_genes <- lapply(c("0","75","150", "250", "500","1000"), get_dosage_gene_vector) %>%
                      unlist() %>%
                      unique()

chase_2 <- chase_merged[which(chase_merged$Nearest.Ensembl %in% all_exp_genes) , ] 
#write.csv(chase_2, "~/repliATAC/chase_test_2.csv")

rna_seq2 <- rna_seq[, c(8:13)] %>% 
    as.matrix() 
rna_seq2 <- log2(rna_seq2 + 1)

rna_seq2 <- t(scale(t(rna_seq2)))

rna_seq[,c(8:13)] <- rna_seq2

chase_2 <- left_join(chase_2, rna_seq, by=join_by(Nearest.Ensembl == ensembl))

P3F_independent <- c("0c","0c_75c","0c_75c_500c","0c_500c")
P3F_dependent <- c("75c","75c_500c","500c")

chase_2$p3f_status <- NA

chase_2[which(chase_2$category_string %in% P3F_dependent), 98 ] <- "P3F_dependent"

chase_2[which(chase_2$category_string %in% P3F_independent), 98 ] <- "P3F_independent"

independent <- chase_2[which(chase_2$p3f_status == "P3F_independent"),]

dependent <- chase_2[which(chase_2$p3f_status == "P3F_dependent"),]

delete_1 <- independent[which( (independent$gene_name %in% dependent$gene_name) ),]$peak_ids #remove overlap only from independent--one dependent peak = dependent gene

keep_these_indices <- which( !(chase_2$peak_ids %in% delete_1) )

chase_2 <- chase_2[keep_these_indices,]
#write.csv(chase_2, "~/repliATAC/chase_test_2.csv")


#need to unique gene names...
chase_2$duplicate_row <- FALSE
for(w in 1:nrow(chase_2)){

    if(chase_2[w, ncol(chase_2) ]){
        next
    }else{

    these_ones <- which( chase_2$gene_name == chase_2[w, 89 ])

    minus_that_one <- these_ones[which(these_ones != w)]

    chase_2[minus_that_one, ncol(chase_2)] <- TRUE

    }

}

chase_2 <- chase_2[which(!chase_2$duplicate_row),]

#Table S2
#write.csv(chase_2, "~/dosage_manuscript/figure_5/chase_inde_dep_sites.csv", row.names=F)

chase_2_long <- chase_2 %>%
                pivot_longer(cols=c(Control,X75ng.ml,X500ng.ml), names_to = "Dosage", values_to="Expression")

p3f_de <- subset(chase_2_long, p3f_status == "P3F_dependent")
p3f_inde <- subset(chase_2_long, p3f_status == "P3F_independent")

aov_res <- aov(Expression ~ Dosage, data = p3f_de)
summary(aov(Expression ~ Dosage, data = p3f_de))
TukeyHSD(aov_res)
#all significant
#500, 0 ****
#0, 75 **
#75, 500 ****

aov_res <- aov(Expression ~ Dosage, data = p3f_inde)
summary(aov(Expression ~ Dosage, data = p3f_inde))
TukeyHSD(aov_res)
#0, 500 **

#Figure 5F, left

my_comparisons <- list(c("Control", "X75ng.ml"),c("X500ng.ml", "X75ng.ml"),c("Control", "X500ng.ml"))
####### need to use actual anova/Tukey res & p-values...will need to use a different fxn to plot this

p1 <- ggplot(p3f_de, aes(x=Dosage, y=Expression, fill = p3f_status))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size=16)+
    scale_fill_manual(values=c("#BB5566"))+
    labs(title="PAX3::FOXO1 dependent\naccessible sites", x="Doxycycline dosage (ng/ml)", y="Scaled log2 expression")+
    theme(plot.title = element_text( hjust=0.5),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=16))+
    scale_x_discrete(limits =  c("Control","X75ng.ml","X500ng.ml"), labels=c("0","75","500"))+
    geom_hline(yintercept = 0, linewidth = 1, linetype = 2, color="#DDAA33")+
    geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, annotation=c("**","****","****"), y_position=c(2,2.5,3))

my_comparisons2 <- list(c("Control", "X500ng.ml"))

p2 <- ggplot(p3f_inde, aes(x=Dosage, y=Expression, fill = p3f_status))+
    geom_boxplot(show.legend=F, notch=T, linewidth =1.2)+
    theme_classic(base_size=16)+
    scale_fill_manual(values=c("#004488"))+
    labs(title="PAX3::FOXO1 independent\naccessible sites", x="Doxycycline dosage (ng/ml)", y="")+
    theme(plot.title = element_text( hjust=0.5),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=16))+
    scale_x_discrete(limits =  c("Control","X75ng.ml","X500ng.ml"), labels=c("0","75","500"))+
    geom_hline(yintercept = 0, linewidth = 1, linetype = 2, color="#DDAA33")+
    geom_signif(comparisons = my_comparisons2, map_signif_level = TRUE, annotation=c("**"), y_position=c(2,2.5,3))+
    ylim2(p1)

p3 <- p1 %>% insert_right(p2) 


 png(paste0("~/dosage_manuscript/figure_5/chase_all_repliATAC_sites_in_dep_exp_revised.png" ), width = 7.5, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
  print(p3)
  dev.off()

#Figure 5F, right

enriched_chase1 <- enrichr(chase_2[which(chase_2$p3f_status  == "P3F_dependent"),25], dbs_to_use)

clust_1 <- plot_gene_sets(enriched_chase1[["MSigDB_Hallmark_2020"]],  plot_title = "P3F dependent", num_sets=5,  fill_col="#BB5566")

enriched_chase2 <- enrichr(chase_2[which(chase_2$p3f_status  == "P3F_independent"),25], dbs_to_use)

clust_2 <- plot_gene_sets(enriched_chase2[["MSigDB_Hallmark_2020"]],  plot_title = "P3F independent", num_sets=5, fill_col="#004488")

png( paste0("~/dosage_manuscript/figure_5/hallmarks_enrichr_P3F_in_dependent_chase_plot_revised.png"), width = 5.5, height = 7, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(
plot_grid(clust_1, clust_2,ncol=1)
)
dev.off()









