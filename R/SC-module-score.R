library(msigdbr)
library(Seurat)

# for each cell type in scRNA-seq or spot data
# e.g., cholagiocytes stored in obj_chol.rds
obj_seurat <- readRDS("/path/to/obj_chol.rds")
# C5 -- GO:BP
gmt <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
gmt <- unique(gmt[,c("gs_name","gene_symbol")])
colnames(gmt) <- c("term","gene")
gmt$term <- tolower(gmt$term)
geneSets <- lapply(unique(gmt$term),function(x){gmt$gene[gmt$term==x]})
names(geneSets) <- unique(gmt$term)
obj_seurat <- AddModuleScore(object = obj_seurat, features = gmt)
