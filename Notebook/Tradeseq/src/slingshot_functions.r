
library(RColorBrewer)
library(ggplot2)
library(scales)
library(viridis)
library(scales)
library(slingshot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(patchwork)
library(igraph)
library(ggbeeswarm)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(ggrepel)




plot_slingshot_trajectory_2 <- function(seurat_obj, slingshot_obj, embedding = "PCA",
                                      annot_column = "cell_label") {

  # cell embeddings from the SAME space used by Slingshot
  if (embedding %in% SingleCellExperiment::reducedDimNames(slingshot_obj)) {
    emb <- SingleCellExperiment::reducedDim(slingshot_obj, embedding)
  } else {
    stop(sprintf("Embedding '%s' not found in slingshot_obj reducedDims.", embedding))
  }
  cell_df <- as.data.frame(emb[, 1:2, drop = FALSE])
  colnames(cell_df) <- c("x", "y")
  cell_df$Cell_Type <- seurat_obj[[annot_column, drop = TRUE]]

  # curve coordinates (also in the SAME space)
  curves <- slingshot::slingCurves(slingshot_obj, as.df = TRUE)
  # pick first two numeric coordinate columns and rename to x,y
  keep <- names(curves)[sapply(curves, is.numeric)][1:2]
  curves_df <- curves[, c(keep, "Lineage", "Order")]
  colnames(curves_df)[1:2] <- c("x", "y")

  ggplot(cell_df, aes(x = x, y = y, color = Cell_Type)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_path(data = dplyr::arrange(curves_df, Lineage, Order),
              aes(x = x, y = y, group = Lineage),
              inherit.aes = FALSE, color = "black", size = 1) +
    theme_classic() +
    ggtitle(sprintf("Slingshot MST & Curves (%s space)", embedding))
}


plot_slingshot_trajectory <- function(seurat_obj, slingshot_obj, annot_column="cell_label") {

    curves <- slingCurves(slingshot_obj, as.df=T)
    sling_mst_obj <- slingMST(slingshot_obj, as.df=T)
    centers <- as.data.frame(do.call(rbind, V(slingMST(slingshot_obj))$coordinates))
    labels <- V(slingMST(slingshot_obj))$name
    centers$label <- labels

    cell_embeddings <- as.data.frame(Embeddings(seurat_obj, "umap"))
    cell_embeddings$Cell_Type <- seurat_obj[[annot_column, drop = TRUE]]  

    cell_plot <- ggplot(cell_embeddings, aes(x = umap_1, y = umap_2, color = Cell_Type)) +
        geom_point(alpha = 0.6, size = 4)  +  
        theme_classic() +
        theme(
            axis.title = element_text(size = 20),
            axis.text = element_text(size = 16)
        ) +
        ggtitle("Cells with Slingshot MST Overlay")

    cell_plot <- cell_plot +
    geom_path(data = curves %>% arrange(Order),
            aes(x = umap_1, y = umap_2, group = Lineage),
            color = "black", size = 1)

    cell_plot <- cell_plot +
        ggrepel::geom_text_repel(
          data = centers,
          aes(x = umap_1, y = umap_2, label = label),
          color = "black",
          size = 8,
          box.padding = 0.5,
          point.padding = 0.3,
          max.overlaps = Inf
        )


    return(cell_plot)

}



plot_slingshot_lineage <- function(seurat_obj, slingshot_obj, annot_column="cell_label") {
    slingshot_obj_test <- slingshot_obj
    colData(slingshot_obj_test)$selected_ident <- as.character(seurat_obj[[annot_column, drop = TRUE]])
    slinsghot_data <- colData(slingshot_obj_test)
    column_obj <- grep(pattern="slingPseudotime", x=names(slinsghot_data))
    slingshot_subset <- as.data.frame(slinsghot_data[, column_obj, drop=FALSE])
    slingshot_subset$selected_ident <- slinsghot_data$selected_ident

    slingshot_long <- slingshot_subset %>%
        select(selected_ident, starts_with("slingPseudotime_")) %>% 
        pivot_longer(
            cols = starts_with("slingPseudotime_"),
            names_to = "Lineage",
            values_to = "Pseudotime"
        ) %>%
        drop_na(Pseudotime)

    n_values <- nlevels(as.factor(slingshot_long$selected_ident))
    select_colors <- colorRampPalette(brewer.pal(9, "Set1"))(n_values)

    plot_slingshot_all <- ggplot(slingshot_long, aes(x = Pseudotime, y = Lineage)) +
        geom_quasirandom(
            groupOnX = FALSE,
            size=4,
            aes(color = factor(selected_ident)),
            alpha = 1
        ) +
        scale_color_manual(values = select_colors, name = "Clusters") +
        theme_classic() +
        theme(text = element_text(size = 20)) +
        xlab("Pseudotime") +
        ylab("Lineages") +
        ggtitle("Slingshot Pseudotime by Lineage")

    return(plot_slingshot_all)

}



library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(RColorBrewer)
library(SingleCellExperiment)


plot_slingshot_lineage_clean <- function(seurat_obj, slingshot_obj, annot_column="cell_label") {
  slingshot_obj_test <- slingshot_obj
  colData(slingshot_obj_test)$selected_ident <- as.character(seurat_obj[[annot_column, drop = TRUE]])
  slinsghot_data <- colData(slingshot_obj_test)

  # map lineage index -> clusters on that lineage (from MST path)
  lineage_paths <- slingshot::slingLineages(slingshot_obj_test)
  lineage_clusters <- lapply(lineage_paths, function(x) as.character(x))
  names(lineage_clusters) <- paste0("slingPseudotime_", seq_along(lineage_clusters))

  # build long df of pseudotime
  column_obj <- grep(pattern="slingPseudotime", x=names(slinsghot_data))
  slingshot_subset <- as.data.frame(slinsghot_data[, column_obj, drop=FALSE])
  slingshot_subset$selected_ident <- slinsghot_data$selected_ident

  slingshot_long <- slingshot_subset %>%
    dplyr::select(selected_ident, dplyr::starts_with("slingPseudotime_")) %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("slingPseudotime_"),
      names_to = "Lineage",
      values_to = "Pseudotime"
    ) %>%
    tidyr::drop_na(Pseudotime) %>%
    # keep only clusters that belong to that lineage’s MST path
    dplyr::rowwise() %>%
    dplyr::filter(selected_ident %in% lineage_clusters[[Lineage]]) %>%
    dplyr::ungroup()

  # keep cluster color mapping stable across plots
  all_clusters <- levels(as.factor(seurat_obj[[annot_column, drop = TRUE]]))
  n_values <- length(all_clusters)
  select_colors <- stats::setNames(
    colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n_values),
    all_clusters
  )

  plot_slingshot_all <- ggplot2::ggplot(
      slingshot_long,
      ggplot2::aes(x = Pseudotime, y = Lineage)
    ) +
    ggbeeswarm::geom_quasirandom(
      groupOnX = FALSE,
      size = 4,
      ggplot2::aes(color = factor(selected_ident, levels = all_clusters)),
      alpha = 1
    ) +
    ggplot2::scale_color_manual(values = select_colors, name = "Clusters") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::xlab("Pseudotime") +
    ggplot2::ylab("Lineages") +
    ggplot2::ggtitle("Slingshot Pseudotime by Lineage")

  return(plot_slingshot_all)
}



plot_slingshot_trajectory_fixed <- function(
    seurat_obj,
    slingshot_obj,
    embedding = "PCA",
    annot_column = "cell_label",
    pt_size = 1.5,
    pt_alpha = 0.6
) {

  library(ggplot2)
  library(dplyr)
  library(SingleCellExperiment)
  library(slingshot)
  library(ggrepel)

  # --- get embeddings from slingshot object ---
  if (!(embedding %in% SingleCellExperiment::reducedDimNames(slingshot_obj))) {
    stop(sprintf("Embedding '%s' not found in slingshot_obj reducedDims.", embedding))
  }

  emb <- SingleCellExperiment::reducedDim(slingshot_obj, embedding)
  emb <- as.data.frame(emb[, 1:2, drop = FALSE])
  colnames(emb) <- c("x", "y")

  # --- make sure cell labels align with slingshot_obj cell order ---
  sce_cells <- colnames(slingshot_obj)
  seu_cells <- colnames(seurat_obj)

  common_cells <- intersect(sce_cells, seu_cells)
  if (length(common_cells) == 0) {
    stop("No overlapping cell names between seurat_obj and slingshot_obj.")
  }

  # reorder embeddings and labels to same cell order
  emb <- emb[match(common_cells, sce_cells), , drop = FALSE]
  labels <- seurat_obj[[annot_column, drop = TRUE]][match(common_cells, seu_cells)]

  cell_df <- emb
  cell_df$Cell_Type <- labels
  cell_df$cell_id <- common_cells

  # --- slingshot curves ---
  curves <- slingshot::slingCurves(slingshot_obj, as.df = TRUE)

  # keep first 2 numeric columns as coordinates
  numeric_cols <- names(curves)[sapply(curves, is.numeric)]
  coord_cols <- numeric_cols[1:2]

  curves_df <- curves[, c(coord_cols, "Lineage", "Order")]
  colnames(curves_df)[1:2] <- c("x", "y")
  curves_df <- curves_df %>% arrange(Lineage, Order)

  # --- MST cluster centers ---
  mst_graph <- slingMST(slingshot_obj)
  centers <- as.data.frame(do.call(rbind, igraph::V(mst_graph)$coordinates))
  colnames(centers)[1:2] <- c("x", "y")
  centers$label <- igraph::V(mst_graph)$name

  ggplot(cell_df, aes(x = x, y = y, color = Cell_Type)) +
    geom_point(alpha = pt_alpha, size = pt_size) +
    geom_path(
      data = curves_df,
      aes(x = x, y = y, group = Lineage),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 1.1
    ) +
    ggrepel::geom_text_repel(
      data = centers,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      color = "black",
      size = 5
    ) +
    theme_classic() +
    labs(
      title = sprintf("Slingshot trajectories in %s space", embedding),
      x = paste0(embedding, "_1"),
      y = paste0(embedding, "_2"),
      color = "Cell Type"
    )
}


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
        ) +
        ggtitle("Cells with Slingshot MST Overlay")

    cell_plot <- cell_plot +
        geom_path(
            data = curves_trimmed %>% dplyr::arrange(curve_id, Order),
            aes(x = umap_1, y = umap_2, group = curve_id),
            inherit.aes = FALSE,
            color = "black",
            linewidth = line_size
        )

    cell_plot <- cell_plot +
        ggrepel::geom_text_repel(
            data = centers,
            aes(x = umap_1, y = umap_2, label = label),
            inherit.aes = FALSE,
            color = "black",
            size = 8,
            box.padding = 0.5,
            point.padding = 0.3,
            max.overlaps = Inf
        )

    return(cell_plot)
}