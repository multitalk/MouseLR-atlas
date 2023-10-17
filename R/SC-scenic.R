library(Seurat)
library(SCENIC)
# for each cell type in scRNA-seq
# e.g., cholagiocytes stored in obj_chol.rds
obj_seurat <- readRDS("/path/to/obj_chol.rds")
exprMat <- as.matrix(obj_seurat[["RNA"]]@counts)
cellInfo <- obj_seurat@meta.data
cellInfo$CellType <- cellInfo$subtype
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")
dbDir <- "/path/to/cisTarget_databases"
scenicOptions <- initializeScenic(org="mgi", dbDir=dbDir, nCores=10)
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
saveRDS(aucell_regulonAUC,"aucell_regulonAUC.rds")
regulons <- loadInt(scenicOptions, "regulons")
saveRDS(regulons,"regulons.rds")

# Cell-type specific regulators (RSS): 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )
rssPlot <- plotRSS(rss)
saveRDS(rss,"rss.rds")
saveRDS(rssPlot,"rssPlot.rds")