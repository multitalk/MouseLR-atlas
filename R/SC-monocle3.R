library(SeuratWrappers)
library(monocle3)
# for each cell type in scRNA-seq
# e.g., cholagiocytes stored in obj_chol.rds
obj_seurat <- readRDS("/path/to/obj_chol.rds")
obj_cds <- as.cell_data_set(obj_seurat)
obj_cds <- cluster_cells(cds = obj_cds, reduction_method = "UMAP")
cellcluster <- obj_seurat$subtype
obj_cds@clusters@listData[["UMAP"]][["clusters"]] <- factor(obj_mp$subtype)
obj_cds <- learn_graph(obj_cds, use_partition = FALSE)
obj_cds <- order_cells(obj_cds)
plot_cells(
  cds = obj_cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)