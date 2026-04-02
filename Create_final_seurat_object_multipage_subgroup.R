library(Seurat)
library(readxl)
library(tidyverse)
library(gridExtra)
library(scattermore)
library(ZemmourLib)
library(RColorBrewer)

args = commandArgs(TRUE)
matrix_path = args[1]
query_all_metadata = args[2]
output_file_query = args[3]
mde_incremental_path = args[4]
mde_plot_path = args[5]
annotation = as.character(args[6])
subgroup = as.character(args[7])
prefix = args[8]

mypal_level1 <- c(ZemmourLib::immgent_colors$level1, "not classified" = "black")
level1_subgroup <- mypal_level1
level1_subgroup[] <- "black"
mypal_level2 <- c(ZemmourLib::immgent_colors$level2, level1_subgroup, "not classified" = "black")
n = 70
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mypal1 = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
mypal1 = mypal1[-4]

library("pals")
parade = function(n) { return(Seurat::DiscretePalette(n, palette = "parade", shuffle = F)) }

length(glasbey())
length(polychrome())
mypal = c(glasbey(), polychrome(), mypal1)
names(mypal) = NULL

#mde <- read.csv("N:/CBDM_Lab/Odhran/scVerse/ImmgenT_freeze_20250109/mde_plot.csv", row.names = 1)
mde <- read.csv(mde_plot_path, row.names = 1)
mde <- mde[mde$level1 == subgroup, ]

#annotation = "celltypes"
#so <- readRDS("N:/CBDM_Lab/Odhran/Data_Integration_Webpage/Plotting/multiple groups/David_Foxp3_creERT2_TdT_seurat_object.Rds")
#so <- readRDS(so_path)
files <- list.files(matrix_path)

#matrix_in <- Read10X(matrix_path)

#matrix_in <- Read10X(
#  data.dir = matrix_path,
#  gene.column = 2,
#  unique.features = TRUE,
#  strip.suffix = FALSE
#)

read10x_flexible <- function(dir) {
  files <- list.files(dir)

  pick <- function(base) {
    if (paste0(base, ".gz") %in% files) return(paste0(base, ".gz"))
    if (base %in% files) return(base)
    stop("Missing ", base, " or ", base, ".gz in ", dir)
  }

  mtx  <- pick("matrix.mtx")
  bc   <- pick("barcodes.tsv")
  feat <- if ("features.tsv.gz" %in% files || "features.tsv" %in% files) {
    pick("features.tsv")
  } else {
    pick("genes.tsv")
  }

  library(Matrix)
  mat <- readMM(file.path(dir, mtx))
  barcodes <- read.delim(file.path(dir, bc), header = FALSE)
  features <- read.delim(file.path(dir, feat), header = FALSE)
  features[, 1] <- make.unique(features[, 1])

  rownames(mat) <- features[, 1]
  colnames(mat) <- barcodes[, 1]

  mat <- as(mat, "dgCMatrix")
  mat
}

files <- list.files(matrix_path)
if (("features.tsv.gz" %in% files)) {
        matrix_in <- Read10X(matrix_path)
} else {
        matrix_in <- read10x_flexible(matrix_path)

}

#matrix_in <- read10x_flexible(matrix_path)
so <- CreateSeuratObject(matrix_in)
so[['RNA']] <- CreateAssayObject(so@assays$RNA$counts)
#colnames(so) <- gsub("-1", "", colnames(so))


query_all_metadata <- read.csv(query_all_metadata, row.names = 1)
#rownames(query_all_metadata) <- gsub("-1", "", rownames(query_all_metadata))
output_file_query <- read.csv(output_file_query, row.names = 1)
#rownames(output_file_query) <- gsub("-1", "", rownames(output_file_query))
mde_incremental <- read.csv(mde_incremental_path, row.names = 1)
#rownames(mde_incremental) <- paste0(rownames(mde_incremental), "-1")

#mypal_annotation <- setNames(c(mypal,mypal)[1:length(unique(query_all_metadata[, annotation]))], unique(query_all_metadata[, annotation]))

joint.bcs = intersect(rownames(query_all_metadata), rownames(output_file_query))
print(length(joint.bcs))
joint.bcs = intersect(joint.bcs, rownames(mde_incremental))
print(length(joint.bcs))
query_all_metadata <- query_all_metadata[joint.bcs, ]
#output_file <- output_file[joint.bcs, ]
output_file_query <- output_file_query[joint.bcs, ]
#mde_incremental <- read.csv(mde_incremental_path, row.names = 1)
#rownames(mde_incremental) <- paste0(rownames(mde_incremental), "-1")
mde_incremental <- mde_incremental[joint.bcs, ]
colnames(mde_incremental) <- c("level2_MDE1", "level2_MDE2")
write.csv(mde_incremental, paste0(prefix, "/", subgroup, "_mde.csv"), row.names = T, quote = F)


#annotation_conversion <- read.csv("/n/groups/cbdm_lab/odc180/ImmgenT_workshop/ImmgenT_freeze_20250109/new_annotations/annotation_table_20260201_conversion_May_to_now.csv")
##View(so@meta.data)
#for (old_cluster in annotation_conversion$level2.old) {
#  new_cluster <- annotation_conversion$level2.new[annotation_conversion$level2.old == old_cluster]
#  new_level1 <- annotation_conversion$level1.new[annotation_conversion$level2.old == old_cluster]
#  output_file_query$level2_final[output_file_query$level2_final == old_cluster] <- new_cluster
#  #so@meta.data$level1_final[output_file_query$level2_final == old_cluster] <- new_level1
#}

#mde_incremental <- read.csv(mde_incremental_path, row.names = 1)
#rownames(mde_incremental) <- paste0(rownames(mde_incremental), "-1")
#mde_incremental <- mde_incremental[rownames(mde_incremental) %in% joint.bcs, ]
#colnames(mde_incremental) <- c("level2_MDE1", "level2_MDE2")
#write.csv(mde_incremental, paste0(prefix, "/", subgroup, "_mde.csv"), row.names = T, quote = F)
#joint.bcs <- intersect(rownames(output_file), rownames(mde_incremental))
#output_file[joint.bcs, "level2_MDE1"] <- mde_incremental[joint.bcs, "level2_MDE1"]
#output_file[joint.bcs, "level2_MDE2"] <- mde_incremental[joint.bcs, "level2_MDE2"]

output_file_query$level2_MDE1 <- "-3"
output_file_query$level2_MDE2 <- "3"
output_file_query[joint.bcs, ]$level2_MDE1 <- mde_incremental[joint.bcs, ]$level2_MDE1
output_file_query[joint.bcs, ]$level2_MDE2 <- mde_incremental[joint.bcs, ]$level2_MDE2
output_file_query$level2_MDE1 <- as.numeric(output_file_query$level2_MDE1)
output_file_query$level2_MDE2 <- as.numeric(output_file_query$level2_MDE2)

metadata <- cbind(query_all_metadata, output_file_query)

so <- so[, joint.bcs]
so <- AddMetaData(so, metadata)

tmp = table(so$level2_final)
write.table(tmp, paste0(prefix, "/level2_final_table.txt"),quote = F, row.names = F, col.names = F, sep = "\t")

# level2
#level2_coords <- output_file_query[, c("level2_MDE1", "level2_MDE2")]
level2_coords <- mde_incremental
level2_coords_matrix <- as.matrix(level2_coords)
rownames(level2_coords_matrix) <- colnames(so)
colnames(level2_coords_matrix) <- c("level2_1", "level2_2")
level2_dimreduc <- CreateDimReducObject(level2_coords_matrix, key = "level2MDE_", assay = "RNA")
so@reductions$mde_incremental_level2 <- level2_dimreduc

## Cluster cells and create metadata columns for Rosetta
DefaultAssay(so) <- "RNA"
so <- so %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 1e4, verbose = T) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = T) %>%
  ScaleData(features = VariableFeatures(.), verbose = T) %>%     # note the "."
  RunPCA(features = VariableFeatures(.), npcs = 50, verbose = T, reduction.name = "pca_rna") %>%
  FindNeighbors(reduction = "pca_rna", dims = 1:30, verbose = T) %>%
  RunUMAP(reduction = "pca_rna", dims = 1:30, seed.use = 123, verbose = T, reduction.name = "umap_rna") %>%
  FindClusters(resolution = c(0.1, 0.25, 0.4, 0.5, 0.75, 0.8, 1, 2, 3))

so$RNA_clusters <- so@meta.data[, annotation]
so$sample_name <- so@meta.data[, "IGT"]
#so$sample_name <- so$batch_sample

#colnames(so) <- paste0(EXP_name, ".", colnames(so))

## Plots for webpage
level2_confidence <- ggplot(so@meta.data, mapping = aes(x = level2_scanvi_confidence)) + geom_density(fill = "red", alpha = 0.5) + theme_minimal() + xlim(0,1) + ggtitle("level2 scanvi score distribution")

## Plot with ggplot
so <- AddMetaData(so, so[['umap_rna']]@cell.embeddings)

metadata_plot <- so@meta.data

## Remove IGT cells from plots
query_cells <- rownames(metadata_plot)[!(grepl("IGT", rownames(metadata_plot)))]
metadata_plot <- metadata_plot[query_cells, ]
metadata_subset <- metadata_plot[metadata_plot$level2_final == "not classified", ]
ncells_plotted <- length(rownames(metadata_plot))
mypal_annotation <- setNames(mypal, unique(metadata_plot[, annotation]))

percentages <- as.data.frame(table(metadata_plot$level2_final))
rownames(percentages) <- percentages$Var1
percentages$Var1 <- NULL

percentages <- (percentages / sum(percentages$Freq)) * 100
percentages_subset <- percentages[!(rownames(percentages) %in% c("not classified", "nonT", "remove", "unclear", "thymocyte")), , drop = F]

level2_not_classified_percentage <- round((length(rownames(metadata_subset)) / length(query_cells)) * 100, 2)
query_only_not_classified <-  ggplot() + geom_scattermore(metadata_plot, mapping = aes(umaprna_1, umaprna_2), colour = "grey", pointsize = 1.5)  +
  geom_point(metadata_subset, mapping = aes(umaprna_1, umaprna_2, colour = "level2_final"), colour = "red", size = 1/log10(ncells_plotted+1)) + theme_classic() +
  theme(
    axis.title = element_blank(), 
    plot.title = element_text(size = 15)) +
  ggtitle(paste0('Query cells: "not classified" (red - ', level2_not_classified_percentage, '% ) overlayed onto classified (grey)')) 
#query_only_not_classified <-  DimPlot(so[, so$level2_final %in% c(subgroups, "not classified")], reduction = "umap_rna", group.by = "level2_final", pt.size = 1/log10(ncells_plotted+1))  

query_only <- DimPlot(so, reduction = "umap_rna", group.by = "level2_final") + theme(
  axis.title = element_blank())

ncells_plotted <- length(query_cells)
#subgroup_level2_final <- ggplot() + geom_scattermore(data = mde, mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "gray", size = 1/log10(ncells_plotted+1))  +
#  geom_point(data = metadata_plot, mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_final), size = 1/log10(ncells_plotted+1)) + 
#  theme_classic() + scale_colour_manual(values = mypal_level1) + guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) +
#  theme(
#    axis.title = element_blank(), 
#    plot.title = element_text(size = 15)) + ggtitle('Query cells (immgenT annotation) overlayed onto immgenT (grey)') + xlim(-3, 3) + ylim(-3, 3)
#allT_level1_final <- DimPlot(so, reduction = "mde_incremental_allT", group.by = "level1_final", pt.size = 1/log10(ncells_plotted+1)) + xlim(-2,2) + ylim(-2,2)

subgroup_level2_final <- ggplot() + geom_scattermore(data = mde, mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "gray", size = 1/log10(ncells_plotted+1))  +
      geom_point(data = metadata_plot, mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_final, shape = level2_final), size = 1/log10(ncells_plotted+1)) +
      theme_classic() + scale_colour_manual(values = mypal_level2) + guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) +
      scale_shape_manual(values = setNames(ifelse(levels(factor(metadata_plot$level2_final)) == "not classified", 17, 16), levels(factor(metadata_plot$level2_final)))) +
      theme(
        axis.title = element_blank(),
        plot.title = element_text(size = 15)) + ggtitle(paste0('Query cells (immgenT annotation) overlayed onto ', subgroup, ' immgenT (grey)')) + NoAxes() +
     xlim(min(mde$level2_MDE1), max(mde$level2_MDE1)) + ylim(min(mde$level2_MDE2), max(mde$level2_MDE2))


subgroup_annotation <- ggplot() + geom_scattermore(data = mde, mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "gray", size = 1/log10(ncells_plotted+1))  +
geom_point(data = metadata_plot, mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = !!sym(annotation)), size = 1/log10(ncells_plotted+1)) + 
theme_classic() + scale_colour_manual(values = mypal_annotation) + guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) +  # + scale_colour_manual(values = mypal_annotation)
theme(
  axis.title = element_blank(), 
  plot.title = element_text(size = 18)) + ggtitle('Query cells (paper annotation) overlayed onto immgenT (grey)') + xlim(min(mde$level2_MDE1), max(mde$level2_MDE1)) + ylim(min(mde$level2_MDE2), max(mde$level2_MDE2))
#allT_annotation <- DimPlot(so, reduction = "mde_incremental_allT", group.by = annotation, pt.size = 1/log10(ncells_plotted+1)) + xlim(-2,2) + ylim(-2,2)

plots <- list(level2_confidence, query_only_not_classified, subgroup_level2_final, subgroup_annotation)

## Discovery Plot and score
## Create Discovery score plot
scores_tbl_list2 = list()
summaries_list2 = list()

library(RANN)
knn_novelty_scores <- function(old_mat, new_mat, k = 10) {
  old_mat <- as.matrix(old_mat)
  new_mat <- as.matrix(new_mat)
  stopifnot(ncol(old_mat) >= 2, ncol(new_mat) >= 2, ncol(old_mat) == ncol(new_mat))
  if (nrow(new_mat) <= k) stop("Need k < number of new points.")

  # Mean distance to k nearest REFERENCE neighbors
  ref_knn <- nn2(data = old_mat, query = new_mat, k = k)
  mean_d_ref <- rowMeans(ref_knn$nn.dists)

  # Mean distance to k nearest QUERY neighbors (leave-one-out: drop the self at distance 0)
  qry_knn <- nn2(data = new_mat, query = new_mat, k = k + 1)
  # First column is self (distance 0). Remove it; take the next k neighbors.
  mean_d_qry <- rowMeans(qry_knn$nn.dists[, -1, drop = FALSE])

  # Scores
  ratio <- mean_d_ref / mean_d_qry
  score_log <- log(ratio)
  score_sym <- (mean_d_ref - mean_d_qry) / (mean_d_ref + mean_d_qry)  # in [-1, 1]

  tibble::tibble(
    idx = seq_len(nrow(new_mat)),
    mean_d_ref = mean_d_ref,
    mean_d_qry = mean_d_qry,
    ratio = ratio,
    score_log = score_log,
    score_sym = score_sym
  )
}

tmp = so
tmp = tmp[,tmp@meta.data %>% filter(!level2_final %in% c("nonT", "unclear")) %>% row.names()]
# tmp = tmp %>%
#     NormalizeData(normalization.method = "LogNormalize", verbose = F) %>%
#     FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>%
#     ScaleData(features = VariableFeatures(.), verbose = F) %>%
#     RunPCA(features = VariableFeatures(.), npcs = 30, verbose = F) %>%
#     FindNeighbors(dims = 1:30, verbose = F) %>%
#     RunUMAP(dims = 1:30, seed.use = 123, verbose = T)
if (!"mde_incremental_level2" %in% names(tmp@reductions)) {
  message("Skipping : no mde_incremental_allT")
  next
}
if (!"umap_rna" %in% names(tmp@reductions)) {
  message("Skipping : no umap_rna")
  next
}
if (!"pca_rna" %in% names(tmp@reductions)) {
  message("Skipping : no pca_rna")
  next
}

tmp$is_noclassl2 = F
tmp$is_noclassl2[tmp$level2_final %in% c(unique(so$level2_C_scANVI))] = T
old_mat = tmp[["pca_rna"]]@cell.embeddings[tmp$is_noclassl2 == F,]
new_mat = tmp[["pca_rna"]]@cell.embeddings[tmp$is_noclassl2 == T,]
scores_tbl = knn_novelty_scores(old_mat, new_mat, k = 10) %>% as.data.frame()
rownames(scores_tbl) = rownames(new_mat)
scores_tbl$dataset = "dataset"
scores_tbl$mde_incremental_allT_dim1 = tmp[["mde_incremental_level2"]]@cell.embeddings[rownames(scores_tbl),1]
scores_tbl$mde_incremental_allT_dim2 = tmp[["mde_incremental_level2"]]@cell.embeddings[rownames(scores_tbl),2]
scores_tbl_list2[["dataset"]] = scores_tbl

tmp$score_log = NA
tmp$score_log = scores_tbl$score_log[match(colnames(tmp), rownames(scores_tbl))]
tmp$discovery = tmp$score_log > log(2)
tmp$ratio = scores_tbl$ratio[match(colnames(tmp), rownames(scores_tbl))]

discovery_score_allT <-  ggplot() + geom_scattermore(data = mde, mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "gray", size = 1/log10(ncells_plotted+1))  +
  geom_point(data = tmp@meta.data[!is.na(tmp@meta.data$ratio), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = ratio)) +
  theme_classic() + scale_colour_gradient2(low = "black", mid = "white",high = "red", midpoint = 5, limits = c(0,10)) +
  theme(
    axis.title = element_blank(),
    plot.title = element_text(size = 15)) + ggtitle('Query cells discovery score overlayed onto immgenT (grey)') + xlim(min(mde$level2_MDE1), max(mde$level2_MDE1)) + ylim(min(mde$level2_MDE2), max(mde$level2_MDE2))

#discovery_score <-  ggplot() + geom_scattermore(data = mde, mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "gray", size = 1/log10(ncells_plotted+1))  +
#  geom_point(data = tmp@meta.data[!is.na(tmp@meta.data$score_log), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = discovery)) +
#  theme_classic() + scale_colour_manual(values = c("FALSE" = "Black", "TRUE" = "Red")) +
#  theme(
#    axis.title = element_blank(),
#    plot.title = element_text(size = 15)) + ggtitle('Query cells discovery score overlayed onto immgenT (grey)') + xlim(min(mde$level2_MDE1), max(mde$level2_MDE1)) + ylim(min(mde$level2_MDE2), max(mde$level2_MDE2))

discovery_score <- ggplot() + theme_void()

pdf(paste0(prefix, "/TRBI.pdf"), height = 15, width = 15)
## page 1 - TRBI results
p <- grid.arrange(grobs = plots[1:4], ncol = 2)

## page 2 - gene expression
DefaultAssay(so) <- "RNA"
genes <- c("Foxp3", "Cd4", "Cd8b1", "Cd8a", "Trdc", "Zbtb16") 
genes <- genes[genes %in% rownames(so)]
plots_feature  <- FeaturePlot(so, reduction = "mde_incremental_level2", features = genes, order = T, combine = F)
plots_feature <- lapply(plots_feature, function(p) {
  p + xlim(-3.01, 3.01) + ylim(-3.01, 3.01)
})

plot_features <- CombinePlots(plots_feature, ncol = 2)
print(plot_features)

## page 3 - not classified UMAP, discovery score allT
p <- grid.arrange(plots[[3]], discovery_score_allT, discovery_score, ncol = 2, nrow = 2, heights = c(0.5, 0.5))

dev.off()

#saveRDS(so, paste0(prefix, "/Final_Seurat_object.Rds"))
write.csv(percentages, paste0(prefix, "/percentages_final.csv"))

metadata_plot$confidence_score <- metadata_plot$level2_scanvi_confidence
overlapping_cells <- intersect(rownames(metadata_plot), rownames(tmp@meta.data))
#metadata_plot$score_log <- NA
#metadata_plot[overlapping_cells, ]$score_log <- tmp@meta.data[overlapping_cells, ]$score_log
#metadata_plot$discovery <- metadata_plot$score_log > log(2)
metadata_plot$score <- NA
metadata_plot[overlapping_cells, ]$score <- tmp@meta.data[overlapping_cells, ]$ratio
metadata_plot$discovery <- metadata_plot$score> 2
metadata_plot$cellID <- rownames(metadata_plot)
final_output_file <- metadata_plot[, c("cellID", "level2_final", "confidence_score", "score", "discovery")]
colnames(final_output_file) <- c("cellID", "level2", "confidence_score", "discovery_score", "discovery")
write.csv(final_output_file, paste0(prefix, "/Final_user_output_file.csv"), row.names = F, col.names = T, quote = F)

