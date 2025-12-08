## Creating genes adata
.libPaths("/n/groups/cbdm_lab/odc180/Python/conda/R_4.4.2_clean_env/lib/R/library")

library(Seurat)
#library(DropletUtils)
library(Matrix)
library(MASS)
library(SeuratDisk)

args = commandArgs(TRUE)
so_path <- args[1]
so_name <- args[2]
output_dir <- args[3]
h5_name <- args[4]
subgroup <- args[5]

## Filter for common genes
so_integrated <- readRDS("/n/groups/cbdm_lab/odc180/ImmgenT_workshop/ImmgenT_freeze_20250109/igt1_96_withtotalvi20250109.Rds")
## Comment if querying the whole reference
so_integrated <- so_integrated[,so_integrated$annotation_level1 %in% subgroup]
so_integrated
#ADT <- so_integrated
#DefaultAssay(ADT) <- "ADT"
#ADT[['RNA']] <- NULL
#proteins <- rownames(ADT)

query <- readRDS(paste0(so_path, so_name))
query <- DietSeurat(
	query,
	counts = TRUE,
	data = TRUE,
	scale.data = FALSE,
	features = NULL,
	assays = "RNA",
	dimreducs = NULL,
	graphs = NULL)
##query@assays$RNA@scale.data <- query@assays$RNA@counts
query@assays$RNA@data <- query@assays$RNA@counts 
query

## Create anndata object

query <- query[rownames(query@assays$RNA) %in% rownames(so_integrated@assays$RNA),]
so_integrated <- so_integrated[rownames(so_integrated@assays$RNA) %in% rownames(query@assays$RNA) ,]
so_integrated
query

genes <- rownames(query)
write.csv(genes, "genes.csv")
## query ADT data
#query_ADT_matrix <- matrix(0, 180, length(colnames(query)))
#colnames(query_ADT_matrix) <- colnames(query)
#rownames(query_ADT_matrix) <- proteins
#query_ADT_matrix <- as(query_ADT_matrix, Class = "sparseMatrix")
##query[["ADT"]] <- CreateAssayObject(query_ADT_matrix)
##so_integrated[['ADT']] <- integrated_ADT

## save query metadata
ImmgenT_metadata <- so_integrated@meta.data
ImmgenT_IGTHT <- as.data.frame(ImmgenT_metadata[, c("IGTHT", "IGT")])
colnames(ImmgenT_IGTHT) <- c("IGTHT", "IGT")
##query_metadata$IGTHT <- paste0(query_metadata$batch, ".", query_metadata$sample)
write.csv(ImmgenT_IGTHT, "ImmgenT_IGTHT.csv")

## save query metadata
query_metadata <- query@meta.data
#query_IGTHT <- as.data.frame(query_metadata[, c("organ_batch", "batch")])
query_IGTHT <- as.data.frame(query_metadata[, c("IGTHT", "IGT")])
colnames(query_IGTHT) <- c("IGTHT", "IGT")
##query_metadata$IGTHT <- paste0(query_metadata$batch, ".", query_metadata$sample)
write.csv(query_IGTHT, "query_IGTHT.csv")
write.csv(query_metadata, "query_all_metadata.csv")

## Combine into one seurat object, create anndata object
so_integrated@meta.data <- ImmgenT_IGTHT
query@meta.data <- query_IGTHT
merged <- merge(so_integrated, query, merge.data = F) 
merged

saveRDS(merged, "seurat_object_merged_RNA.Rds")
## Save as anndata objecta
SaveH5Seurat(merged, paste0(output_dir, "MERGED_RNA_", h5_name))
Convert(paste0(output_dir, "MERGED_RNA_", h5_name), dest = "h5ad", raw = FALSE)


##ADT <- so_integrated
##DefaultAssay(ADT) <- "ADT"
##ADT[['RNA']] <- NULL
##so_integrated[['ADT']] <- NULL

## Combine into one seurat object, create anndata object
##merged <- merge(so_integrated, query, merge.data = T) 

##SaveH5Seurat(merged, paste0(output_dir, "RNA_", h5_name))
##Convert(paste0(output_dir, "RNA_combined_", h5_name), dest = "h5ad")

## ADT adata
##ADT_merged <- merge(ADT, query_ADT, merge.data = T)
##SaveH5Seurat(ADT, paste0(output_dir, "ADT_", h5_name))
##Convert(paste0(output_dir, "ADT_reference", h5_name), dest = "h5ad")

#query_ADT <- CreateSeuratObject(counts = query_ADT_matrix)
#query_ADT[['ADT']] <- query_ADT[['RNA']]
#DefaultAssay(query_ADT) <- "ADT"
#query_ADT[['RNA']] <- NULL

#ADT_merged <- merge(ADT, query_ADT, merge.data = T)
#saveRDS(ADT_merged, "seuratobject_Merged_ADT.Rds")

## save as Anndata object - broke
##SaveH5Seurat(ADT_merged, paste0(output_dir, "MERGED_ADT_", h5_name))
##Convert(paste0(output_dir, "MERGED_ADT_", h5_name), dest = "h5ad")
## save as 10xCounts and construct mdata object separately
#ADT_matrix <- ADT_merged@assays$ADT@counts
#rownames(ADT_matrix) <- rownames(ADT)
#colnames(ADT_matrix) <- colnames(merged)
#write10xCounts(
#  paste0(so_path, "/ADT_matrix"),
#  ADT_matrix,
#  #query_df,
#  barcodes = colnames(merged),
#  gene.id = rownames(ADT),
#  gene.symbol = rownames(ADT),
#  overwrite = TRUE
#)
#


## ADT adata
##SaveH5Seurat(ADT, paste0(output_dir, "ADT_", h5_name))
##Convert(paste0(output_dir, "ADT_query_", h5_name), dest = "h5ad")
