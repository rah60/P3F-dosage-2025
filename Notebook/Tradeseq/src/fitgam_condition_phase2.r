library(slingshot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(BiocParallel)
library(tradeSeq)
library(qs)

BiocParallel::register(BiocParallel::MulticoreParam())
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 16

seurat_obj    <- qread("../Data/out/seurat_phase2_annotated.qs")
slingshot_obj <- qread("./out/slingshot_obj_phase2.qs")

var_genes <- VariableFeatures(seurat_obj, assay = "RNA")
lyr_names <- Layers(seurat_obj[["RNA"]])
lyr_to_use <- if ("counts" %in% lyr_names) "counts" else lyr_names[1]  
mat <- GetAssayData(seurat_obj, assay = "RNA", layer = lyr_to_use)
counts <- mat[var_genes, ]

sds <- SlingshotDataSet(slingshot_obj)

meta <- seurat_obj@meta.data
cond_vec <- dplyr::recode(
    meta$orig.ident,
    "SEG0503" = "2w_dox500",
    "SEG0504" = "noDox_48h",
    "SEG0505" = "noDox_96h",
    .default   = NA_character_
)

condition <- factor(
    cond_vec,
    levels = c("noDox_48h", "noDox_96h", "2w_dox500")
)
names(condition) <- rownames(meta)

condition <- condition[colnames(counts)]
stopifnot(identical(names(condition), colnames(counts)))
if (any(is.na(condition))) {
    stop("Some cells have missing condition labels. Check orig.ident mapping.")
}

qsave(condition, "./out/fitGAM_conditions_by_cell_input_phase2.rds")

set.seed(42)
fitGAM_res <- tradeSeq::fitGAM(
  counts     = counts,
  sds        = sds,
  conditions = condition,   
  parallel   = TRUE,
  BPPARAM    = BPPARAM
)

qsave(fitGAM_res, "./out/fitGAM_conditions_by_cell_result_phase2.qs")
