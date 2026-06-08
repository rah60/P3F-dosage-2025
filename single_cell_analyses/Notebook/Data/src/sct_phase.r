library(Seurat)
library(qs)


library(future)

options(future.globals.maxSize = 550 * 1024^3)  # 550 GB

plan(multicore, workers = 16)

seurat_obj <- qread("./out/seurat_obj_phase2_lowqiality_filtered.qs")

seurat_obj <- SCTransform(
  seurat_obj,
  vars.to.regress = c("percent.mt", "nCount_RNA"),
  vst.flavor = "v2",
  verbose = FALSE
)


qsave(seurat_obj, "./out/seurat_obj_phase2_postsct.qs")