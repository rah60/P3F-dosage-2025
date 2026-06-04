suppressPackageStartupMessages({
    library(Seurat)
    library(qs)
    library(dplyr)
    library(openxlsx)
    library(igraph)
    library(dplyr)
    library(tidygraph)
    library(ggraph)
    library(ggrepel)
    library(dorothea)
    library(patchwork)
})


run_DE_and_write <- function(seurat_obj, object_name, outdir,
                             min.pct = 0.10,
                             logfc.threshold = 0.25,
                             test.use = "wilcox",
                             only.pos = FALSE,
                             top_n_per_cluster = 50) {
  message("Running FindAllMarkers for: ", object_name)
  markers <- FindAllMarkers(
    object = seurat_obj,
    assay = DefaultAssay(seurat_obj),
    only.pos = only.pos,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    test.use = test.use
  )
  if (nrow(markers) == 0) {
    warning("No markers found for ", object_name, " with current thresholds.")
  }

  markers <- markers %>%
    arrange(cluster, p_val_adj, desc(avg_log2FC))

  out_file <- file.path(outdir, paste0(object_name, "_DE_results.xlsx"))

  wb <- createWorkbook()
  addWorksheet(wb, "all_markers")
  writeData(wb, "all_markers", markers)
  clusters <- unique(markers$cluster)
  for (cl in clusters) {
    sheet_name <- paste0("cluster_", cl)
    if (nchar(sheet_name) > 31) sheet_name <- substr(sheet_name, 1, 31)
    top_markers <- markers %>%
      filter(cluster == cl) %>%
      arrange(p_val_adj, desc(avg_log2FC)) %>%
      head(top_n_per_cluster)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, top_markers)
  }

  saveWorkbook(wb, out_file, overwrite = TRUE)
  message("Wrote DE results to: ", out_file)
  write.csv(markers, file = file.path(outdir, paste0(object_name, "_DE_results_full.csv")), row.names = FALSE)
  message("Also wrote CSV: ", file.path(outdir, paste0(object_name, "_DE_results_full.csv")))

  return(markers)
}

make_tf_network_plot <- function(seurat_obj = NULL,
                                 markers_df = NULL,
                                 gene_list = NULL,
                                 cluster = NULL,
                                 lfc_cutoff = 1.0,
                                 p_cutoff = 0.05,
                                 dorothea_confidence = c("A"),
                                 layout = "fr",
                                 seed = 1,
                                 top_n_labels = Inf,
                                 save_to = NULL) {

  pkgs <- c("igraph","dplyr","tidygraph","ggraph","ggrepel","dorothea","tibble")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) stop("Please install package: ", p)
  library(dplyr)
  library(tibble)
  library(tidygraph)
  library(ggraph)
  library(ggrepel)
  library(dorothea)
  library(igraph)

  data(dorothea_hs)   

  if (!is.null(markers_df)) {
    if (!("avg_log2FC" %in% colnames(markers_df))) {
      stop("markers_df must contain column 'avg_log2FC' (or change column name accordingly).")
    }
    if (!is.null(cluster)) {
      de_genes <- markers_df %>%
        filter(cluster == !!cluster,
               avg_log2FC > lfc_cutoff,
               p_val_adj <= p_cutoff) %>%
        pull(gene) %>%
        unique()
    } else {
      de_genes <- markers_df %>%
        filter(avg_log2FC > lfc_cutoff,
               p_val_adj <= p_cutoff) %>%
        pull(gene) %>%
        unique()
    }
  } else if (!is.null(gene_list)) {
    de_genes <- unique(gene_list)
  } else if (!is.null(seurat_obj) && !is.null(cluster) && !is.null(assay <- DefaultAssay(seurat_obj))) {
    message("No markers_df or gene_list provided. Attempting quick FindMarkers for cluster vs rest from the supplied Seurat object.")
    if (!("Seurat" %in% class(seurat_obj))) stop("seurat_obj must be a Seurat object.")
    if (is.null(Idents(seurat_obj))) stop("seurat_obj must have identities set (Idents).")
    if (!(cluster %in% levels(Idents(seurat_obj)))) stop("Cluster '", cluster, "' not found in Idents(seurat_obj).")
    fm <- FindMarkers(seurat_obj, ident.1 = cluster, only.pos = TRUE,
                      logfc.threshold = lfc_cutoff, min.pct = 0.1, test.use = "wilcox")
    if (nrow(fm) == 0) stop("FindMarkers returned no hits. Adjust thresholds or supply markers_df/gene_list.")
    de_genes <- rownames(fm)[which(fm$p_val_adj <= p_cutoff)]
  } else {
    stop("Provide either markers_df, gene_list, or (seurat_obj + cluster) to derive de_genes.")
  }

  if (length(de_genes) == 0) {
    stop("No DE genes found with the provided filters.")
  }

  tf_network <- dorothea_hs %>%
    filter(confidence %in% dorothea_confidence) %>%
    filter(target %in% de_genes) %>%
    distinct(tf, target, .keep_all = TRUE)

  if (nrow(tf_network) == 0) {
    stop("No TF-target pairs found in dorothea_hs for the provided DE genes and confidence level.")
  }

  edges <- tf_network %>% select(tf, target) %>% distinct()
  nodes <- tibble(name = unique(c(edges$tf, edges$target))) %>%
    mutate(type = ifelse(name %in% edges$tf, "TF", "Gene"))

  ig <- igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
  graph <- tidygraph::as_tbl_graph(ig)

  graph <- graph %>%
    activate(nodes) %>%
    mutate(deg = centrality_degree(mode = "all"))

  set.seed(seed)
  p <- ggraph(graph, layout = layout) +
    geom_edge_link(alpha = 0.25,
                   arrow = arrow(length = unit(2.5, "mm")),
                   end_cap = circle(2, "mm")) +
    geom_node_point(aes(size = ifelse(type == "TF", deg + 6, deg + 1),
                        shape = type,
                        color = type),
                    alpha = 0.95) +
    geom_node_text(aes(label = ifelse(type == "TF" & (deg >= 1), name, NA_character_)),
                   repel = TRUE, size = 3.5, max.overlaps = Inf) +
    scale_shape_manual(values = c(TF = 17, Gene = 16)) +
    guides(size = "none") +
    theme_void(base_size = 14) +
    theme(legend.position = "bottom")

  if (is.finite(top_n_labels) && top_n_labels > 0 && top_n_labels < Inf) {
    node_df <- as_tibble(graph, active = "nodes")
    top_genes <- node_df %>%
      filter(type == "Gene") %>%
      arrange(desc(deg)) %>%
      slice_head(n = top_n_labels) %>%
      pull(name)
    if (length(top_genes) > 0) {
      p <- p + geom_node_text(data = node_df %>% filter(name %in% top_genes),
                              aes(x = x, y = y, label = name),
                              repel = TRUE, size = 3)
    }
  }

  if (!is.null(save_to)) {
    ggsave(filename = save_to, plot = p, width = 10, height = 10)
    message("Saved plot to: ", save_to)
  }

  return(list(plot = p, graph = graph, de_genes = de_genes, tf_network = tf_network))
}


comparative_dotplot <- function(seu_list,
                                group_names = NULL,
                                genes,
                                cluster_order,
                                cluster_col = "cell_label",
                                title = "Comparative dotplot",
                                color_low = "#2166AC",
                                color_mid = "grey85",
                                color_high = "#B2182B",
                                midpoint = 0,
                                dotplot_width = 10,
                                dotplot_height = 10,
                                save_to = NULL) {

  pkgs <- c("Seurat", "dplyr", "ggplot2", "tibble")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) stop("Please install package: ", p)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tibble)

  if (inherits(seu_list, "Seurat")) {
    seu_list <- list(sample1 = seu_list)
  } else if (is.list(seu_list) && all(sapply(seu_list, inherits, "Seurat"))) {
    # ok
  } else {
    stop("seu_list must be a Seurat object or a list of Seurat objects.")
  }

  n_objs <- length(seu_list)
  if (!is.null(group_names)) {
    if (length(group_names) != n_objs) stop("group_names length must equal number of Seurat objects.")
    names(seu_list) <- group_names
  } else {
    if (is.null(names(seu_list)) || any(names(seu_list) == "")) {
      names(seu_list) <- paste0("Group", seq_len(n_objs))
    }
  }

  make_dot_df_single <- function(seu, group_name, genes, cluster_order, cluster_col) {
    if (!(cluster_col %in% colnames(seu@meta.data))) {
      stop("cluster_col '", cluster_col, "' not found in Seurat object's meta.data.")
    }

    Idents(seu) <- cluster_col
    seu[[cluster_col]] <- factor(seu[[cluster_col]][,1], levels = cluster_order)
    Idents(seu) <- seu[[cluster_col]][,1]
    seu_sub <- subset(seu, idents = cluster_order)

    genes_present <- intersect(genes, rownames(seu_sub))
    if (length(genes_present) == 0) {
      warning(group_name, ": none of the genes were found in the object. Returning NULL for this group.")
      return(NULL)
    }

    dp <- DotPlot(seu_sub, features = genes_present)$data
    dp <- dp %>%
      mutate(
        group = group_name,
        cluster = factor(id, levels = cluster_order),
        gene = factor(features.plot, levels = genes_present)
      )

    dp <- dp %>%
      select(gene, cluster, pct.exp, avg.exp.scaled, group, everything())

    return(dp)
  }

  df_list <- lapply(seq_along(seu_list), function(i) {
    nm <- names(seu_list)[i]
    make_dot_df_single(seu = seu_list[[i]],
                       group_name = nm,
                       genes = genes,
                       cluster_order = cluster_order,
                       cluster_col = cluster_col)
  })

  keep_idx <- !sapply(df_list, is.null)
  if (!any(keep_idx)) stop("No genes present in any supplied Seurat objects.")
  df_list <- df_list[keep_idx]

  df_all <- bind_rows(df_list)

  df_all <- df_all %>%
    mutate(
      gene = factor(as.character(gene), levels = unique(as.character(gene))),
      cluster = factor(as.character(cluster), levels = cluster_order),
      group = factor(group, levels = unique(names(seu_list)[keep_idx]))
    )

  p <- ggplot(df_all, aes(x = gene, y = cluster)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    facet_wrap(~ group, ncol = 1) +
    theme_classic(base_size = 14) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 14),
      axis.text.y = element_text(size = 14),
      panel.grid.major = element_line(linewidth = 0.2),
      strip.text = element_text(face = "bold", size = 14),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    ) +
    labs(
      title = title,
      color = "Avg. expr\n(scaled)",
      size  = "% expressed"
    ) +
    scale_color_gradient2(low = color_low, mid = color_mid, high = color_high, midpoint = midpoint)
  if (!is.null(save_to)) {
    ggsave(filename = save_to, plot = p, width = dotplot_width, height = dotplot_height)
    message("Saved dotplot to: ", save_to)
  }

  return(list(plot = p, df = df_all))
}

subset_tf_graph_plot <- function(graph,
                                 tf_subset,
                                 layout = "fr",
                                 seed = 1) {

  suppressPackageStartupMessages({
    library(tidygraph)
    library(ggraph)
    library(ggrepel)
    library(igraph)
    library(dplyr)
  })

  if (!inherits(graph, "tbl_graph")) {
    stop("graph must be a tidygraph::tbl_graph object")
  }

  nodes_df <- as_tibble(graph, active = "nodes")
  keep_nodes <- nodes_df %>%
    filter(
      (type == "TF" & name %in% tf_subset) |
      name %in% (
        graph %>%
          activate(edges) %>%
          as_tibble() %>%
          filter(from %in% which(nodes_df$name %in% tf_subset)) %>%
          pull(to) %>%
          { nodes_df$name[.] }
      )
    ) %>%
    pull(name)

  sub_graph <- graph %>%
    activate(nodes) %>%
    filter(name %in% keep_nodes) %>%
    activate(edges) %>%
    filter(.N()$name[from] %in% tf_subset)

  sub_graph <- sub_graph %>%
    activate(nodes) %>%
    mutate(deg = centrality_degree(mode = "all"))

  sub_graph <- sub_graph %>%
      activate(nodes) %>%
      mutate(type = ifelse(name %in% tf_subset, "TF", "Gene"))

  set.seed(seed)
  p <- ggraph(sub_graph, layout = layout) +
    geom_edge_link(
      alpha = 0.3,
      arrow = arrow(length = unit(2.5, "mm")),
      end_cap = circle(2, "mm")
    ) +
    geom_node_point(
      aes(
        size = ifelse(type == "TF", deg + 6, deg + 1),
        shape = type,
        color = type
      ),
      alpha = 0.95
    ) +
    geom_node_text(
      aes(label = ifelse(type == "TF", name, NA_character_)),
      repel = TRUE,
      size = 4
    ) +
    scale_shape_manual(values = c(TF = 17, Gene = 16)) +
    guides(size = "none") +
    theme_void(base_size = 14) +
    theme(legend.position = "bottom")

  return(list(graph = sub_graph, plot = p))
}


rank_tfs_by_centrality <- function(graph,
                                   tf_subset = NULL,
                                   mode = "all",
                                   top_n = NULL) {
  suppressPackageStartupMessages({
    library(tidygraph)
    library(dplyr)
  })

  if (!inherits(graph, "tbl_graph")) {
    stop("graph must be a tidygraph::tbl_graph object")
  }

  tf_table <- graph %>%
    activate(nodes) %>%
    mutate(
      degree = centrality_degree(mode = mode),
      betweenness = centrality_betweenness(directed = TRUE, normalized = TRUE),
      closeness = centrality_closeness(mode = "all")
    ) %>%
    as_tibble() %>%
    filter(type == "TF")

  if (!is.null(tf_subset)) {
    tf_table <- tf_table %>%
      filter(name %in% tf_subset)
  }

  tf_table <- tf_table %>%
    arrange(desc(degree))

  if (!is.null(top_n)) {
    tf_table <- tf_table %>% slice_head(n = top_n)
  }

  tf_table
}

