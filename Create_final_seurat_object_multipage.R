library(Seurat)
library(readxl)
library(tidyverse)
library(gridExtra)
library(scattermore)
library(ZemmourLib)
library(RColorBrewer)

args = commandArgs(TRUE)
so_path = args[1]
query_all_metadata = args[2]
output_file_query = args[3]
mde_plot_path = args[4]
annotation = as.character(args[5])
prefix = args[6]
EXP_name = args[7]

mypal_level1 <- c(ZemmourLib::immgent_colors$level1, "not classified" = "black", "nonT" = "darkgray")
level1_subgroup <- mypal_level1
level1_subgroup[] <- "black"
mypal_level2 <- c(ZemmourLib::immgent_colors$level2, level1_subgroup, "not classified" = "black", "nonT" = "darkgray")
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

#annotation = "celltypes"
#so <- readRDS("N:/CBDM_Lab/Odhran/Data_Integration_Webpage/Plotting/multiple groups/David_Foxp3_creERT2_TdT_seurat_object.Rds")
so <- readRDS(so_path)

query_all_metadata <- read.csv(query_all_metadata, row.names = 1)
output_file_query <- read.csv(output_file_query, row.names = 1)
#mypal_annotation <- setNames(c(mypal,mypal)[1:length(unique(query_all_metadata[, annotation]))], unique(query_all_metadata[, annotation]))

joint.bcs = intersect(rownames(query_all_metadata), rownames(output_file_query))
query_all_metadata <- query_all_metadata[joint.bcs, ]
#output_file <- output_file[joint.bcs, ]
output_file_query <- output_file_query[joint.bcs, ]

metadata <- cbind(query_all_metadata, output_file_query)

so <- so[, joint.bcs]
so <- AddMetaData(so, metadata)

tmp = table(so$level1_final)
write.table(tmp, paste0(prefix, "/level1_final_table.txt"),quote = F, row.names = F, col.names = F, sep = "\t")
tmp = table(so$level2_final)
write.table(tmp, paste0(prefix, "/level2_final_table.txt"),quote = F, row.names = F, col.names = F, sep = "\t")


# level1
allT_coords <- output_file_query[, c("allT_MDE1", "allT_MDE2")]
allT_coords_matrix <- as.matrix(allT_coords)
rownames(allT_coords_matrix) <- colnames(so)
colnames(allT_coords_matrix) <- c("allTMDE_1", "allTMDE_2")
allT_dimreduc <- CreateDimReducObject(allT_coords_matrix, key = "allTMDE_", assay = "RNA")
so@reductions$mde_incremental_allT <- allT_dimreduc

# level2
level2_coords <- output_file_query[, c("level2_MDE1", "level2_MDE2")]
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
so$sample_name <- so@meta.data[, "IGTHT"]
#so$sample_name <- so$batch_sample

colnames(so) <- paste0(EXP_name, ".", colnames(so))

## Plots for webpage
level1_confidence <- ggplot(so@meta.data, mapping = aes(x = level1_scanvi_confidence)) + geom_density(fill = "skyblue", alpha = 0.5) + theme_minimal() + xlim(0,1) + ggtitle("level1 scanvi score distribution")
level2_confidence <- ggplot(so@meta.data, mapping = aes(x = level2_scanvi_confidence)) + geom_density(fill = "red", alpha = 0.5) + theme_minimal() + xlim(0,1) + ggtitle("level2 scanvi score distribution")

## Plot with ggplot
so <- AddMetaData(so, so[['umap_rna']]@cell.embeddings)
metadata_plot <- so@meta.data[!(so$level1_C_scANVI %in% c("nonT", "unclear")),]
metadata_plot[, annotation] <- as.character(metadata_plot[, annotation])
subgroups <- unique(so$level1_C_scANVI)
subgroups <- subgroups[!(subgroups %in% c("nonT", "remove", "unclear", "thymocyte"))]
#ncells_plotted <- length(colnames(so))
metadata_subset <- metadata_plot[metadata_plot$level2_final %in% c(subgroups, "not classified"), ]

## Remove IGT cells from plots
query_cells <- rownames(metadata_plot)[!(grepl("IGT", rownames(metadata_plot)))]
metadata_plot <- metadata_plot[query_cells, ]
metadata_subset <- metadata_subset[rownames(metadata_subset) %in% query_cells, ]
ncells_plotted <- length(rownames(metadata_plot))
mypal_annotation <- setNames(mypal, unique(metadata_plot[, annotation]))

percentages <- as.data.frame(table(metadata_plot$level1_final))
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

query_only <- DimPlot(so, reduction = "umap_rna", group.by = "level1_final") + theme(
  axis.title = element_blank())

ncells_plotted <- length(query_cells)
allT_level1_final <- ggplot() + geom_scattermore(data = mde, mapping = aes(x = allT_MDE1, y = allT_MDE2), colour = "gray", size = 1/log10(ncells_plotted+1))  +
  geom_point(data = metadata_plot, mapping = aes(x = allT_MDE1, y = allT_MDE2, colour = level1_final), size = 1/log10(ncells_plotted+1)) + 
  theme_classic() + scale_colour_manual(values = mypal_level1) + guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) +
  theme(
    axis.title = element_blank(), 
    plot.title = element_text(size = 15)) + ggtitle('Query cells (immgenT annotation) overlayed onto immgenT (grey)') + xlim(-2.5, 2.5) + ylim(-2.5, 2.5)
#allT_level1_final <- DimPlot(so, reduction = "mde_incremental_allT", group.by = "level1_final", pt.size = 1/log10(ncells_plotted+1)) + xlim(-2,2) + ylim(-2,2)
allT_annotation <- ggplot() + geom_scattermore(data = mde, mapping = aes(x = allT_MDE1, y = allT_MDE2), colour = "gray", size = 1/log10(ncells_plotted+1))  +
  geom_point(data = metadata_plot, mapping = aes(x = allT_MDE1, y = allT_MDE2, colour = !!sym(annotation)), size = 1/log10(ncells_plotted+1)) + 
  theme_classic() + scale_colour_manual(values = mypal_annotation) + guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) +  # + scale_colour_manual(values = mypal_annotation)
  theme(
    axis.title = element_blank(), 
    plot.title = element_text(size = 18)) + ggtitle('Query cells (paper annotation) overlayed onto immgenT (grey)') + xlim(-2.5, 2.5) + ylim(-2.5, 2.5)
#allT_annotation <- DimPlot(so, reduction = "mde_incremental_allT", group.by = annotation, pt.size = 1/log10(ncells_plotted+1)) + xlim(-2,2) + ylim(-2,2)

plots <- list(level2_confidence, query_only_not_classified, allT_level1_final, allT_annotation)

#percentages <- as.data.frame(table(metadata_plot$level1_final))
#rownames(percentages) <- percentages$Var1
#percentages$Var1 <- NULL

#percentages <- (percentages / sum(percentages$Freq)) * 100
#percentages <- percentages[!(rownames(percentages) %in% c("not classified", "nonT", "remove", "unclear", "thyocyte"), , drop = F]

pdf(paste0(prefix, "/test_show_results_multipage.pdf"), height = 15, width = 15)
p <- grid.arrange(grobs = plots[1:4], ncol = 2)

if (length(rownames(percentages_subset)[which(percentages_subset > 0)]) == 1) {
  i = 5
  group <- rownames(percentages_subset)[which(percentages_subset > 0)]
  
  group_level2_final <- ggplot() + geom_scattermore(data = mde[mde$level1 == group, ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "gray", size = 1/log10(ncells_plotted+1))  +
    geom_point(data = metadata_plot[metadata_plot$level1_final == group | metadata_plot$level2_final == group, ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_final, shape = level2_final), size = 1/log10(ncells_plotted+1)) + 
    theme_classic() + scale_colour_manual(values = mypal_level2) + guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) + 
    scale_shape_manual(values = setNames(ifelse(levels(factor(metadata_plot$level2_final)) == group, 17, 16), levels(factor(metadata_plot$level2_final)))) +
    theme(
      axis.title = element_blank(), 
      plot.title = element_text(size = 15)) + ggtitle(paste0('Query cells (immgenT annotation) overlayed onto ', group, ' immgenT (grey)')) + NoAxes() +
    xlim(min(mde[mde$level1 == group, ]$level2_MDE1) - 0.1, max(mde[mde$level1 == group, ]$level2_MDE1) + 0.1) + ylim(min(mde[mde$level1 == group, ]$level2_MDE2) - 0.1, max(mde[mde$level1 == group, ]$level2_MDE2) + 0.1)
  plots[[i]] <- group_level2_final
  i = i + 1
  #group_level2_final <- DimPlot(so[, so$level1_final == group], reduction = "mde_incremental_level2", group.by = "level2_final", pt.size = 1/log10(ncells_plotted+1)) + guides(color = guide_legend(ncol = 2))
  
  group_annotation <- ggplot() + geom_scattermore(data = mde[mde$level1 == group | metadata_plot$level2_final == group, ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "gray", size = 1/log10(ncells_plotted+1))  +
    geom_point(data = metadata_plot[metadata_plot$level1_final == group, ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = !!sym(annotation)), size = 1/log10(ncells_plotted+1)) + 
    theme_classic() + scale_colour_manual(values = mypal_annotation) + guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) +    # + scale_colour_manual(values = mypal_annotation)
    theme(
      axis.title = element_blank(), 
      plot.title = element_text(size = 15)) + ggtitle(paste0('Query cells (immgenT annotation) overlayed onto ', group, ' immgenT (grey)')) + NoAxes() +
    xlim(min(mde[mde$level1 == group, ]$level2_MDE1) - 0.1, max(mde[mde$level1 == group, ]$level2_MDE1) + 0.1) + ylim(min(mde[mde$level1 == group, ]$level2_MDE2) - 0.1, max(mde[mde$level1 == group, ]$level2_MDE2) + 0.1) # + xlim(-2.5, 2.5) + ylim(-2.5, 2.5)
  plots[[i]] <- group_annotation
  i = i + 1
  #group_annotation <- DimPlot(so[, so$level1_final == group], reduction = "mde_incremental_level2", group.by = annotation, pt.size = 1/log10(ncells_plotted+1))
  
  ## Arrange plots into one pdf 
  library(gridExtra)
  layout <- rbind(
    c(1, 1, 2, 2),  # p1 spans two columns, p2 on the right
    c(0, 0, 0, 0)
    #c(5, 5, 6, 6) # p4 spans two columns
  )
  #pdf("one_group_page_test.pdf", height = 10, width = 10)
  #p <- grid.arrange(group_level2_final,  group_annotation, layout = layout_matrix, ncol = 2)
  #print(p)
  #dev.off()
}


if (length(rownames(percentages_subset)[which(percentages_subset > 0)]) > 1) {
  groups_sig <- rownames(percentages_subset)[which(percentages_subset > 0)]
  size <- length(groups_sig)
  
  i = 5
  for (group in groups_sig) {
    #plots[[i]] <- DimPlot(so[, so$level1_final == group], reduction = "mde_incremental_level2", group.by = "level2_final") + guides(color = guide_legend(ncol = 2))
    #i = i+1
    #plots[[i]] <- DimPlot(so[, so$level1_final == group], reduction = "mde_incremental_level2", group.by = annotation)
    #i = i+1
    ncells_plotted = length(colnames(so[, so$level1_final == group]))
    group_level2_final <- ggplot() + geom_scattermore(data = mde[mde$level1 == group, ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "gray", size = 1/log10(ncells_plotted+1))  +
      geom_point(data = metadata_plot[metadata_plot$level1_final == group | metadata_plot$level2_final == group, ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_final, shape = level2_final), size = 1/log10(ncells_plotted+1)) + 
      theme_classic() + scale_colour_manual(values = mypal_level2) + guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) +
      scale_shape_manual(values = setNames(ifelse(levels(factor(metadata_plot$level2_final)) == group, 17, 16), levels(factor(metadata_plot$level2_final)))) +
      theme(
        axis.title = element_blank(), 
        plot.title = element_text(size = 15)) + ggtitle(paste0('Query cells (immgenT annotation) overlayed onto ', group, ' immgenT (grey)')) + NoAxes() +
    xlim(min(mde[mde$level1 == group, ]$level2_MDE1) - 0.1, max(mde[mde$level1 == group, ]$level2_MDE1) + 0.1) + ylim(min(mde[mde$level1 == group, ]$level2_MDE2) - 0.1, max(mde[mde$level1 == group, ]$level2_MDE2) + 0.1) # + xlim(-2.5, 2.5) + ylim(-2.5, 2.5)
    plots[[i]] <- group_level2_final
    i = i + 1
    #group_level2_final <- DimPlot(so[, so$level1_final == group], reduction = "mde_incremental_level2", group.by = "level2_final", pt.size = 1/log10(ncells_plotted+1)) + guides(color = guide_legend(ncol = 2))
    
    group_annotation <- ggplot() + geom_scattermore(data = mde[mde$level1 == group, ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "gray", size = 1/log10(ncells_plotted+1))  +
      geom_point(data = metadata_plot[metadata_plot$level1_final == group | metadata_plot$level2_final == group, ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = !!sym(annotation)), size = 1/log10(ncells_plotted+1)) + 
      theme_classic() + scale_colour_manual(values = mypal_annotation) + guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) +  # + scale_colour_manual(values = mypal_annotation)
      theme(
        axis.title = element_blank(), 
        plot.title = element_text(size = 15)) + ggtitle(paste0('Query cells (immgenT annotation) overlayed onto ', group, ' immgenT (grey)')) + NoAxes() +
    xlim(min(mde[mde$level1 == group, ]$level2_MDE1) - 0.1, max(mde[mde$level1 == group, ]$level2_MDE1) + 0.1) + ylim(min(mde[mde$level1 == group, ]$level2_MDE2) - 0.1, max(mde[mde$level1 == group, ]$level2_MDE2) + 0.1) # + xlim(-2.5, 2.5) + ylim(-2.5, 2.5)
    plots[[i]] <- group_annotation
    i = i + 1
    #group_annotation <- DimPlot(so[, so$level1_final == group], reduction = "mde_incremental_level2", group.by = annotation, pt.size = 1/log10(ncells_plotted+1))
    
    #p <- grid.arrange(group_level2_final,  group_annotation, ncol = 1)
    #print(p)
  }
  
  
  ## Arrange plots into one pdf 
  #library(gridExtra)
  #layout <- rbind(
  #  c(1, 1, 2, 2),  # p1 spans two columns, p2 on the right
  #  c(3, 3, 4, 4),
  #  c(5, 5, 6, 6) # p4 spans two columns
  #)
  
  #pdf("multiple_groups_page_test.pdf", height = 10, width = 10)
  #p <- grid.arrange(level2_confidence, query_only_not_classified, allT_level1_final, allT_annotation, CD8_level2_final,  CD8_annotation, layout_matrix = layout, ncol = 4)
  #p <- grid.arrange(grobs = plots, ncol = 2)
  #print(p)
  #dev.off()
}

n_plots <- length(plots)
# --- Remaining plots in consecutive pairs (2 plots per half page) ---
if (n_plots > 4) {
  remaining_plots <- plots[-(1:4)]
  
  # Split remaining plots into pairs (2 plots per pair)
  pairs <- split(remaining_plots, ceiling(seq_along(remaining_plots)/2))
  
  # Process every 2 pairs per page
  for (i in seq(1, length(pairs), by = 2)) {
    #top_pair <- pairs[[i]]
    #bottom_pair <- if (i + 1 <= length(pairs)) pairs[[i + 1]] else NULL
    
    # Prepare top half (2 plots side by side)
    #top_grid <- grid.arrange(grobs = top_pair, ncol = 2)
    
    #if (!is.null(bottom_pair)) {
    #  # Prepare bottom half (2 plots side by side)
    #  bottom_grid <- grid.arrange(grobs = bottom_pair, ncol = 2)
    #  # Stack top and bottom halves on same page
    #  p <- grid.arrange(top_grid, bottom_grid, ncol = 1, heights = c(0.5, 0.5))
    #  print(p)
    #} else {
    #  # Only top half exists on this page
    #  p <- grid.arrange(top_grid, ncol = 1, heights = c(0.5, 0.5))
    #  print(p)
    #}
  #}
    top_pair <- pairs[[i]]
    # Bottom pair (blank if none)
    bottom_pair <- if (i + 1 <= length(pairs)) pairs[[i + 1]] else list(ggplot() + theme_void())
    
    # Combine top and bottom pairs in one page
    grobs_page <- c(top_pair, bottom_pair)
    
    grid.arrange(grobs = grobs_page, ncol = 2, nrow = 2, heights = c(0.5, 0.5))
  }
}

DefaultAssay(so) <- "RNA"
genes <- c("Foxp3", "Cd4", "Cd8b1", "Cd8a", "Trdc", "Zbtb16") 
genes <- genes[genes %in% rownames(so)]
plots  <- FeaturePlot(so, reduction = "mde_incremental_allT", features = genes, order = T, combine = F)
plots <- lapply(plots, function(p) {
  p + xlim(-2.5, 2.5) + ylim(-2.5, 2.5)
})

plot_features <- CombinePlots(plots, ncol = 2)
print(plot_features)

dev.off()

saveRDS(so, paste0(prefix, "/Final_Seurat_object.Rds"))
write.csv(percentages, paste0(prefix, "/percentages_final.csv"))

#predictions_output_file_not_classified <- read.csv("N:/CBDM_Lab/Odhran/Data_Integration_Webpage/Plotting/predictions_output_file_not_classified.csv", row.names = 1)
