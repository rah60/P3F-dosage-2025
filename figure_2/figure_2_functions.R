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
        theme_classic(base_size=20)+
        scale_fill_manual(values=fill_col)+
        scale_y_continuous(expand = c(0,0))+
        labs(x= "Gene sets" , y="-Log10 adjusted p-value", title=plot_title)+
        guides(fill="none")

    return(p1) #return plot
}

# from Ambuj
plot_slingshot_trajectory_trimmed <- function(
        seurat_obj,
        slingshot_obj,
        annot_column = "cell_label",
        reduction = "umap",
        trim_end_to_cluster_median = TRUE,
        trim_start_to_cluster_median = FALSE,
        point_size = 4,
        point_alpha = 0.6,
        line_size = 1,
        cluster_colors = NULL
    ) {

    cell_embeddings <- as.data.frame(Embeddings(seurat_obj, reduction))
    colnames(cell_embeddings)[1:2] <- c("umap_1", "umap_2")
    cell_embeddings$Cell_Type <- seurat_obj[[annot_column, drop = TRUE]]

    cluster_medians <- cell_embeddings %>%
        dplyr::group_by(Cell_Type) %>%
        dplyr::summarise(
            umap_1 = median(umap_1, na.rm = TRUE),
            umap_2 = median(umap_2, na.rm = TRUE),
            .groups = "drop"
        )

    cluster_levels <- unique(as.character(cell_embeddings$Cell_Type))
    if (is.null(cluster_colors)) {
        n_clusters <- length(cluster_levels)
        auto_cols <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_clusters)
        cluster_colors <- setNames(auto_cols, cluster_levels)
    } else {
        if (is.null(names(cluster_colors))) {
            stop("cluster_colors must be a named vector, with names matching cluster labels.")
        }
        missing_clusters <- setdiff(cluster_levels, names(cluster_colors))
        if (length(missing_clusters) > 0) {
            stop(
                paste(
                    "Missing colors for clusters:",
                    paste(missing_clusters, collapse = ", ")
                )
            )
        }
        cluster_colors <- cluster_colors[cluster_levels]
    }

    curve_list <- slingshot::slingCurves(slingshot_obj)
    lineage_paths <- slingshot::slingLineages(slingshot_obj)
    curve_dfs <- vector("list", length(curve_list))
    for (i in seq_along(curve_list)) {

        crv <- curve_list[[i]]
        s <- as.data.frame(crv$s)
        s <- s[, 1:2, drop = FALSE]
        colnames(s) <- c("umap_1", "umap_2")
        s$Order <- seq_len(nrow(s))
        s$curve_id <- paste0("Lineage_", i)

        lineage_clusters <- lineage_paths[[i]]
        start_cluster <- lineage_clusters[1]
        end_cluster   <- lineage_clusters[length(lineage_clusters)]

        if (trim_start_to_cluster_median) {
            start_med <- cluster_medians %>% dplyr::filter(Cell_Type == start_cluster)
            if (nrow(start_med) == 1) {
                d_start <- (s$umap_1 - start_med$umap_1)^2 + (s$umap_2 - start_med$umap_2)^2
                idx_start <- which.min(d_start)
                s <- s[idx_start:nrow(s), , drop = FALSE]
                s$Order <- seq_len(nrow(s))
                s$umap_1[1] <- start_med$umap_1
                s$umap_2[1] <- start_med$umap_2
            }
        }

        if (trim_end_to_cluster_median) {
            end_med <- cluster_medians %>% dplyr::filter(Cell_Type == end_cluster)
            if (nrow(end_med) == 1) {
                d_end <- (s$umap_1 - end_med$umap_1)^2 + (s$umap_2 - end_med$umap_2)^2
                idx_end <- which.min(d_end)
                s <- s[seq_len(idx_end), , drop = FALSE]
                s$Order <- seq_len(nrow(s))
            }
        }

        curve_dfs[[i]] <- s
    }

    curves_trimmed <- dplyr::bind_rows(curve_dfs)
    centers <- as.data.frame(do.call(rbind, V(slingMST(slingshot_obj))$coordinates))
    colnames(centers)[1:2] <- c("umap_1", "umap_2")
    centers$label <- V(slingMST(slingshot_obj))$name
    cell_plot <- ggplot(cell_embeddings, aes(x = umap_1, y = umap_2, color = Cell_Type)) +
        geom_point(alpha = point_alpha, size = point_size) +
        scale_color_manual(values = cluster_colors) +
        theme_classic() +
        theme(
            axis.title = element_text(size = 20),
            axis.text = element_text(size = 16)
        ) #+
        #ggtitle("Cells with Slingshot MST Overlay")

    cell_plot <- cell_plot +
        geom_path(
            data = curves_trimmed %>% dplyr::arrange(curve_id, Order),
            aes(x = umap_1, y = umap_2, group = curve_id),
            inherit.aes = FALSE,
            color = "black",
            linewidth = line_size
        )


    return(cell_plot)
}
