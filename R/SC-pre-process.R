#---------------------------------Functions----------------------------------
# read snRNA-seq data
read_huada <- function(cellname, sc_condition, sc_time, filename){
  get_cell_id <- function(x){
    x_char <- as.character(x)
    x_nchar <- nchar(x_char)
    x_zero <- 5-x_nchar
    x_zero_new <- NULL
    for (i in 1:length(x_zero)) {
      x_zero1 <- x_zero[i]
      if (x_zero1 == 0) {
        x_zero_new <- c(x_zero_new, "")
      }
      if (x_zero1 == 1) {
        x_zero_new <- c(x_zero_new, "0")
      }
      if (x_zero1 == 2) {
        x_zero_new <- c(x_zero_new, "00")
      }
      if (x_zero1 == 3) {
        x_zero_new <- c(x_zero_new, "000")
      }
      if (x_zero1 == 4) {
        x_zero_new <- c(x_zero_new, "0000")
      }
    }
    x_char <- paste0(x_zero_new, x_char)
    return(x_char)
  }
  library(Matrix)
  filename <- paste0(filename, "/")
  filename_sub <- dir(path = filename)
  filename_sub1 <- paste0(filename, filename_sub, "/")
  a_data <- list()
  genename <- NULL
  for (i in 1:length(filename_sub1)) {
    a1 <- Seurat::Read10X(filename_sub1[i], gene.column = 1)
    genename <- c(genename, rownames(a1))
    a_data[[i]] <- a1
  }
  names(a_data) <- filename_sub
  genename <- unique(genename)
  filename_sub <- 1:length(filename_sub)
  filename_sub <- paste0("0",filename_sub)
  barcode_meta <- data.frame()
  for (i in 1:length(a_data)) {
    a1 <- a_data[[i]]
    genename1 <- rownames(a1)
    genename1 <- genename[!genename %in% genename1]
    a1_sub <- matrix(0, nrow = length(genename1), ncol = ncol(a1))
    a1_sub <- as(a1_sub, Class = "dgCMatrix")
    colnames(a1_sub) <- colnames(a1)
    rownames(a1_sub) <- genename1
    a1 <- rbind(a1, a1_sub)
    a1 <- a1[genename,]
    barcode_meta_temp <- data.frame(mouse_id = substr(cellname,1,3), 
                                    lib_raw = names(a_data)[i], lib = filename_sub[i], barcode_raw = colnames(a1),
                                    barcode = paste0(cellname,filename_sub[i],sc_condition, sc_time, "D"),stringsAsFactors = F)
    cellid <- 1:ncol(a1)
    cellid <- get_cell_id(cellid)
    barcode_meta_temp$barcode <- paste0(barcode_meta_temp$barcode, cellid)
    colnames(a1) <- barcode_meta_temp$barcode
    a_data[[i]] <- a1
    barcode_meta <- rbind(barcode_meta, barcode_meta_temp)
  }
  rawdata <- a_data[[1]]
  if (length(a_data) > 1) {
    for (i in 2:length(a_data)) {
      rawdata <- cbind(rawdata, a_data[[i]])
    }
  }
  return(list(rawdata = rawdata, barcode_meta = barcode_meta))
}

#---------------------------------Pre-process--------------------------------
# [1] read scRNA-seq data
# download the SC-Raw.7z and decompress it
setwd("/path/to/SC-Raw/")
sc_meta <- readRDS("/path/to/sc_meta_mouse.rds")
rawdata_list <- list()
sc_meta_barcode <- data.frame()
for (i in 1:nrow(sc_meta)) {
  mouse_id <- sc_meta$mouse_id[i]
  sc_condition <- sc_meta$condition[i]
  sc_time <- sc_meta$time[i]
  sc_time <- substr(sc_time, 4, 4)
  cellname <- paste0(mouse_id, "SC")
  sc_name <- sc_meta$sc_name[i]
  rawdata <- read_huada(cellname, sc_condition, sc_time, sc_name)
  rawdata_list[[i]] <- rawdata[[1]]
  names(rawdata_list)[i] <- paste0(cellname, "-", sc_condition, sc_time, "D")
  sc_meta_barcode <- rbind(sc_meta_barcode, rawdata[[2]])
}

# [2] statistics
library(Seurat)
rawdata_meta <- lapply(X = rawdata_list, FUN = function(x){
  x <- CreateSeuratObject(x)
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^mt-")
  x <- x@meta.data
  return(x)
})

sc_meta_mouse_temp <- lapply(X = rawdata_meta, FUN = function(x){
  x_meta <- c(nrow(x), median(x$nCount_RNA), mean(x$nCount_RNA),
              median(x$nFeature_RNA), mean(x$nFeature_RNA),
              median(x$percent.mt), mean(x$percent.mt))
  return(x_meta)
})
sc_meta_mouse_temp <- as.data.frame(sc_meta_mouse_temp,stringsAsFactors = F)
sc_meta_mouse_temp <- as.data.frame(t(sc_meta_mouse_temp))
colnames(sc_meta_mouse_temp) <- c("cell_num","median_nCount_RNA", "mean_nCount_RNA",
                                  "median_nFeature_RNA", "mean_nFeature_RNA",
                                  "median_percent.mt", "mean_percent.mt")

sc_meta <- cbind(sc_meta,sc_meta_mouse_temp)

# save rds data
saveRDS(rawdata_list,file = "rawdata_list.rds")
saveRDS(sc_meta,file = "sc_meta_mouse.rds")




