#---------------------------------Functions----------------------------------
# revise genes
rev_gene_list <- function(rawdata_list, gene_info){
  library(Matrix)
  for (k in 1:length(rawdata_list)) {
    rawdata <- rawdata_list[[k]]
    rawdata <- rawdata[["RNA"]]@counts
    genename<- rownames(rawdata)
    genename1<- genename[genename %in% gene_info$Symbol]
    genename2<- genename[!genename %in% gene_info$Symbol]
    genename3<- genename2[genename2 %in% gene_info$Synonyms]
    genename4<- rep('NA',length(genename3))
    for (i in 1:length(genename3)) {
      d1<- gene_info[gene_info$Synonyms == genename3[i],]$Symbol
      if(length(d1) == 1){
        genename4[i]<- d1
      }
    }
    genename3<- c(genename1,genename3)
    genename4<- c(genename1,genename4)
    genedata<- data.frame(raw_name = genename3,new_name = genename4,stringsAsFactors = F)
    genedata<- genedata[!genedata$new_name == 'NA',]
    genedata1<- as.data.frame(table(genedata$new_name),stringsAsFactors = F)
    genedata1<- genedata1[genedata1$Freq == 1,]
    genedata<- genedata[genedata$new_name %in% genedata1$Var1,]
    rawdata<- rawdata[genedata$raw_name,]
    rownames(rawdata)<- genedata$new_name
    rawdata_list[[k]] <- rawdata
  }
  return(rawdata_list)
}

#---------------------------------Process-----------------------------------
# [1] quality control
# n_Features > 200 ,percent.mt < 25
library(Seurat)
rawdata_list <- readRDS("/path/to/rawdata_list.rds")
rawdata_list <- lapply(X = rawdata_list, FUN = CreateSeuratObject)
rawdata_list <- lapply(X = rawdata_list, FUN = function(seurat_data){
  seurat_data[["percent.mt"]] <- PercentageFeatureSet(seurat_data, pattern = "^mt-")
  seurat_data <- subset(seurat_data, subset = nFeature_RNA > 200 & percent.mt < 25)
  return(seurat_data)
})

# [2] revise genes with NCBI genes updated in 2022.06.19
gene_info <- readRDS("/path/to/mouse_gene_info.rds")
gene_info <- unique(gene_info[,c("Symbol","Synonyms")])
rawdata_list <- rev_gene_list(rawdata_list = rawdata_list, gene_info = gene_info)
rawdata_list <- lapply(X = rawdata_list, FUN = CreateSeuratObject)
rawdata_list <- lapply(X = rawdata_list, FUN = function(x){
  x_name <- colnames(x)[1]
  x_sample <- paste0(substr(x_name,1,5),"-",substr(x_name,8,11))
  x$sample <- x_sample
  x_ct <- substr(x_name,8,11)
  x$condition_time <- x_ct
  return(x)
})

# [3] integration
rawdata_list <- lapply(X = rawdata_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)
})
features <- SelectIntegrationFeatures(object.list = rawdata_list)
rawdata_seurat <- FindIntegrationAnchors(object.list = rawdata_list, anchor.features = features, reference = 1)
rawdata_seurat <- IntegrateData(anchorset = rawdata_seurat)
DefaultAssay(rawdata_seurat) <- "integrated"
rawdata_seurat <- ScaleData(rawdata_seurat, verbose = FALSE)
rawdata_seurat <- RunPCA(rawdata_seurat, verbose = FALSE)
# rawdata_seurat <- JackStraw(rawdata_seurat, num.replicate = 100, dims = 50)
# rawdata_seurat <- ScoreJackStraw(rawdata_seurat, dims = 1:50)
# JackStrawPlot(rawdata_seurat, dims = 1:50)
# ElbowPlot(rawdata_seurat, ndims = 50)
# dim=30
rawdata_seurat <- RunUMAP(rawdata_seurat, reduction = "pca", dims = 1:30)
rawdata_seurat <- FindNeighbors(rawdata_seurat)
rawdata_seurat <- FindClusters(rawdata_seurat)

