rm(list = ls())
graphics.off()

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(harmony)
library(fgsea)
library(msigdbr)
library(clustree)
library(Cairo)
library(DoubletFinder)
library(stringr)
library(readr)
library(viridis)
library(tidyverse)
library(RColorBrewer)

run_fgsea <- function(markers_list, hallmarks=T, perm_num=1000){
     
    #runs gene set enrichment analysis on marker genes from each single cell cluster 
    #markers_list = output from FindAllMarkers
    #hallmarks = logical for whether or not to use the Hallmarks gene sets (if T), otherwise use use GO BP gene sets
    #perm_num = number of permutations to run in fgsea (nPerm)

      #libraries for GSEA
      library(fgsea)
      library(msigdbr)

      if(hallmarks){ #choose gene set

      h_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") #human hallmark gene sets: H. GO BP: C5

      }else{

      h_gene_sets <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "GO:BP") #human hallmark gene sets: H. GO BP: C5

      }

      h_gene_sets2 <- split(x=h_gene_sets$gene_symbol, f=h_gene_sets$gs_name) #set up gene sets for fgsea
      fgsea.l <- vector("list") #set up list to store output
      cluster_names <- unique(markers_list$cluster) #get vector of cluster names to iterate on
      for(k in 1:length(unique(markers_list$cluster))){ #go through clusters & run fgsea on each
      print(k)
      clust_num <- k - 1

      markers.v <- markers_list[which(markers_list$cluster == cluster_names[k]),2] #avglogFC

      names(markers.v) <- rownames(markers_list[which(markers_list$cluster == cluster_names[k]),])
      
      if(min(markers.v) > 0){ #determine if there are only positive logFC markers
      fgsea.l[[k]] <- fgsea(pathways = h_gene_sets2, stats=markers.v , maxSize=500, minSize=5, nPerm=perm_num, scoreType="pos" ) %>%
                        arrange(-NES) 
                        print("pos") 
                        next }
      if(max(markers.v) < 0){ #negative markers only
      fgsea.l[[k]] <- fgsea(pathways = h_gene_sets2, stats=markers.v , maxSize=500, minSize=5, nPerm=perm_num , scoreType="neg") %>%
                        arrange(-NES)
                        print("neg")
                        next }
      fgsea.l[[k]] <- fgsea(pathways = h_gene_sets2, stats=markers.v , maxSize=500, minSize=5, nPerm=perm_num ) %>% #or both & run correct FGSEA settings
                        arrange(-NES)
                        print("std")
      }
      return(fgsea.l) #returns list of dataframes with fgsea output for each cluster
}

process_scrna_norm <- function(dataset, output_name, num.cores=1, res = 0.2, cell_cycle_file = "~/Dbt_scRNA/cycle.rda", out_path = "~/Dbt_scRNA/"){
        
        #processes Seurat object using NormalizeData, also scales, clusters, and run UMAP
        #dataset = path to merged Seurat object 
        #output_name = output name minus ".rds" and path
        #num.cores = # of cores available for parallelization
        #res = resolution to use for clustering
        #cell_cycle_file = path to cycle.rda from Seurat for calling cell cycle
        #path for output files
       
        scrna_data <- dataset
        
        scrna_data <- NormalizeData(scrna_data, verbose=F) #log norm
        
        #cell cycle annotations
        load(cell_cycle_file)
        scrna_data <- CellCycleScoring(scrna_data, g2m.features = g2m_genes, s.features = s_genes, vars.to.regress="percent.mt") #determine cell cycle phases, regress out mt %

        scrna_data <- FindVariableFeatures(scrna_data, verbose=F)
        scrna_data <- ScaleData(scrna_data, features=rownames(scrna_data), verbose=F)
        
        scrna_data <- RunPCA(scrna_data, verbose=F)

        scrna_data <- FindNeighbors(scrna_data, dims = 1:30, reduction = "pca", verbose=F)
        scrna_data <- FindClusters(scrna_data, resolution = c(res), verbose=F)
        scrna_data <- RunUMAP(scrna_data, dims=1:30, reduction = "pca", verbose=F)

        saveRDS(scrna_data, paste0(out_path,output_name,".rds") ) #save data
}

plot_gene_sets2 <- function(gene_sets, num_sets = 10, fill_col="black", plot_title = ""){

    #plots gene sets from Enrichr output, need to be Hallmark gene sets
    #gene_sets = Enrichr output of Hallmark gene sets
    #num_sets = number of gene sets to plot (starting with most significant)
    #fill_col = color of bars on plot
    #plot_title = title of plot

    subset_set <- gene_sets[which(gene_sets$Adjusted.P.value < 0.05),] #subset only gene sets p < 0.05


    if(nrow(subset_set) == 0){  #stop if no significant gene sets
        return(NULL)
    }
    
    subset_set$log10_adj_pval <- -log10(subset_set$Adjusted.P.value) #calculate -log10 p val
    
    #select top num_sets terms & plot
    if(nrow(subset_set) >= num_sets){
        select_20 <- subset_set[order(subset_set$log10_adj_pval, decreasing=T),][1:num_sets,]
    }else{
        select_20 <- subset_set
    }
    
    order.v <- select_20[order(select_20$log10_adj_pval),1] #order gene sets
    
    #plot
    p1 <- ggplot(select_20) +
        geom_col(aes(x=factor(Term, level=order.v), y=log10_adj_pval, fill="blank"))+
        coord_flip()+
        theme_classic(base_size=16)+
        scale_fill_manual(values=fill_col)+
        scale_y_continuous(expand = c(0,0))+
        labs(x= "Gene sets" , y="-Log10 adjusted p-value", title=plot_title)+
        guides(fill="none")

    return(p1) #return plot
}
