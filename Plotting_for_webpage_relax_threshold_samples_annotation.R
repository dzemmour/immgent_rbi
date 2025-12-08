## Sample output sheet for ImmgenT integration webpage ##

library(Seurat)

#Gradients
library(RColorBrewer)
ColorRamp = rev(colorRampPalette(brewer.pal(n = 7,name = "RdYlBu"))(100))
ColorRamp = rev(colorRampPalette(c('red','white','blue'))(20))
ColorRamp = rev(rainbow(10, end = 4/6))
library(viridis)
ColorRamp = rev(viridis(100))
ColorRamp = rev(cividis(100))
#image(1:100, 1, as.matrix(1:100), col = ColorRamp, main = "Color Palette", xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")

#mycol = c("black","blue","red","darkgreen", "orange","purple", "cyan", "greenyellow",    "salmon", "magenta","pink", "tan", "brown") # c("magenta", "red",  "darkgreen", "cyan", "blue", "blue", "blue", "black", "blue", "blue", "blue", "blue", "orange")

#Large Color palette
library("pals")
n = 70
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mypal1 = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
mypal1 = mypal1[-4]
mypal = c(glasbey(), polychrome(), mypal1)
names(mypal) = NULL

args = commandArgs(TRUE)
output_file_path = args[1]
user_output_file_path = args[2]
prefix = args[3]
query_IGTHT_path = args[4]
annotation_column = args[5]
mde_plot_path = args[6]
threshold = 0.85

output_file <- read.csv(output_file_path, row.names = 1)
colnames(output_file) <- c("level1_C_scANVI", "level1_scanvi_confidence", "level1_final", "level2_C_scANVI", "level2_scanvi_confidence", "level2_final", "allT_MDE1", "level2_MDE1", "allT_MDE2", "level2_MDE2")
# Wrangle output file
#output_file$level2_final <- output_file$level2_C_scANVI
#output_file$level2_final[output_file$level2_scanvi_confidence < threshold] <- "not classified"
#output_file <- output_file[!(output_file$level1_C_scANVI %in% c("nonT", "unclear", "remove")),]
#output_file$level1_final <- output_file$level1_C_scANVI
#output_file$level1_final[output_file$level1_scanvi_confidence < threshold] <- "not classified"
#output_file$level2_final[output_file$level1_final == "not classified"] <- "not classified"
#output_file <- output_file[!(output_file$level1_final == "unclear" & output_file$level2_final == "unclear"), ]

# Add new columns to outpur file for plotting
#groups <- sort(unique(output_file$level1_C_scANVI))
#groups <- groups[!(groups %in% c("unclear", "nonT", "remove"))]
output_file$level2_plot <- output_file$level2_final
##output_file$level2_plot[output_file$level2_final == "unclear" & output_file$level1_final == "Treg"] <- "Treg"
##output_file$level2_plot[output_file$level2_final == "unclear" & output_file$level1_final == "CD4"] <- "CD4"

#for(annot in groups) {
#output_file$level2_final[output_file$level2_final == "not classified" & output_file$level1_final == annot] <- annot
#output_file$level2_final[!(grepl(paste0(annot, "_"), output_file$level2_final)) & output_file$level1_final == annot] <- "not classified"

#}

output_file$level1_assignment <- output_file$level1_final

#output_file$level1_final[output_file$level2_plot == "not classified"] <- "not classified"


#user_output_file <- read.csv(user_output_file_path, row.names = 1)
#metadata_query <- output_file[rownames(output_file) %in% rownames(user_output_file), ]
#metadata <- output_file[!(rownames(output_file) %in% rownames(user_output_file)), ]

# Create new dfs for plotting
user_output_file <- read.csv(user_output_file_path, row.names = 1)
query_IGTHT <- read.csv(query_IGTHT_path, row.names=1)
#query_IGTHT <- query_IGTHT[rownames(query_IGTHT) %in% rownames(user_output_file), ]
metadata_query <- output_file[rownames(output_file) %in% rownames(user_output_file), ]
query_IGTHT <- query_IGTHT[rownames(query_IGTHT) %in% rownames(metadata_query), ]
metadata_query[rownames(metadata_query), "IGTHT"] <- query_IGTHT[rownames(metadata_query), "IGTHT"]
metadata_query[rownames(metadata_query), "annotation"] <- query_IGTHT[rownames(metadata_query) , c(annotation_column)]
metadata <- output_file[!(rownames(output_file) %in% rownames(user_output_file)), ]
query_cells <- rownames(metadata_query)[!(grepl("IGT", rownames(metadata_query)))]
metadata_query <- metadata_query[query_cells, ]

mypal_groups <- setNames(mypal, sort(unique(metadata_query$level1_final)))
mypal_clusters <- setNames(mypal, sort(unique(metadata_query$level2_final)))

mde_plot <- read.csv(mde_plot_path, row.names = 1)

## Downsample for plotting
#library(tidyverse)
library(tidyr)
library(dplyr)
metadata_downsampled <- metadata %>%
  group_by(level2_C_scANVI) %>%
  #slice_sample(n = 1000) %>%
  #slice_sample(n = min(n(), 1000)) %>%  
  group_modify(~ slice_sample(.x, n = min(1000, nrow(.x)))) %>%
  ungroup()

#metadata_downsampled <- rbind(metadata_downsampled, metadata_query)

#mypal_level2 <- setNames(mypal, sort(unique(metadata$level1_final)))
#print(mypal_level2)
#names(mypal_level2)[13:(length(sort(unique(metadata$level2_C_scANVI)))+12)] <- sort(unique(metadata$level2_C_scANVI))
#print(mypal_level2)

library(ggplot2)
library(scattermore)
## Sample plots for webpage
#for (sample in unique(metadata_query$IGTHT)) {

#metadata_query_sample <- metadata_query[metadata_query$IGTHT == sample ,]
metadata_query_sample <- metadata_query

#pdf(paste0(prefix, "/", sample,"_samples_annotation_scattermore_fix_level2_relax_threshold_Sample_plots_ImmgenT_data_integration_webpage.pdf"), height = 10, width = 10)
pdf(paste0(prefix, "/samples_annotation_level2_ImmgenT_data_integration_webpage.pdf"), height = 10, width = 10)
for (annotation in unique(metadata_query_sample$annotation)) {
	metadata_query_plot <- metadata_query_sample[metadata_query_sample$annotation == annotation, ]

# Whole dataset
#ggplot(metadata_downsampled, aes(x = level1_MDE1, y = level1_MDE2)) + geom_point(aes(colour = is_ref), size = 0.5) + 
#  theme_classic() + ggtitle("Query data mapped onto whole MDE") +
#  xlim(-2.5, 2.5) + ylim(-2.5, 2.5)
#ggplot(metadata_downsampled, aes(x = level1_MDE1, y = level1_MDE2)) + geom_point(aes(colour = level1), size = 0.5) + scale_color_manual(values = mypal_level2) + 
#  theme_classic() + ggtitle("ImmgenT + query Annotation level 1") +
#  xlim(-2.5, 2.5) + ylim(-2.5, 2.5) 
plot <- ggplot() + geom_scattermore(mde_plot, mapping = aes(x = allT_MDE1, y = allT_MDE2), colour = "grey", size = 0.75, pixels = c(512, 512)) + 
  geom_point(metadata_query_plot[metadata_query_plot$level1_final == "not classified", ], mapping = aes(x = allT_MDE1, y = allT_MDE2), colour = "darkgray", size = 0.75) +
  geom_point(metadata_query_plot[metadata_query_plot$level1_final != "not classified", ], mapping = aes(x = allT_MDE1, y = allT_MDE2, colour = level1_final), size = 0.75) +
  #geom_point(metadata_query[metadata_query$level1_final == "Not Classified", ], mapping = aes(x = level1_MDE1, y = level1_MDE2), colour = "darkgrey", size = 0.25, shape = 17) +
  scale_color_manual(values = mypal_groups) + 
  scale_shape_manual(values = c("Not Classified" = 17)) +
   guides(color = guide_legend(override.aes = list(size = 2)), 
          shape = guide_legend(override.aes = list(size = 2))) +
  theme_classic() + ggtitle(paste0(annotation, ": ImmgenT Reference (grey background) with query cells with ImmgenT annotation (colour) - Overall MDE")) +
  xlim(-3, 3) + ylim(-3, 3) 
print(plot)


#ggplot() + geom_point(metadata_query, mapping = aes(x = level1_MDE1, y = level1_MDE2, colour = level1_scanvi_confidence), size = 0.75) +
#  scale_color_gradient(low = "orange", high = "black") +
#  theme_classic() + ggtitle("query Annotation level 1 scanvi confidence") +
#  xlim(-3, 3) + ylim(-3, 3)

# No coord lims
plot <- ggplot() + geom_scattermore(mde_plot, mapping = aes(x = allT_MDE1, y = allT_MDE2), colour = "grey", size = 0.75, pixels = c(512, 512)) +
  geom_point(metadata_query_plot[metadata_query_plot$level1_final == "not classified", ], mapping = aes(x = allT_MDE1, y = allT_MDE2), colour = "darkgray", size = 0.75) +
  geom_point(metadata_query_plot[metadata_query_plot$level1_final != "not classified", ], mapping = aes(x = allT_MDE1, y = allT_MDE2, colour = level1_final), size = 0.75) +
  #geom_point(metadata_query[metadata_query$level1_final == "Not Classified", ], mapping = aes(x = level1_MDE1, y = allT_MDE2), colour = "darkgrey", size = 0.25, shape = 17) +
  scale_color_manual(values = mypal_groups) +
  scale_shape_manual(values = c("not classified" = 17)) +
   guides(color = guide_legend(override.aes = list(size = 2)),
          shape = guide_legend(override.aes = list(size = 2))) +
  theme_classic() + ggtitle(paste0(annotation, ":ImmgenT Reference (grey background) with query cells with ImmgenT annotation (colour) - Overall MDE")) +
  #xlim(-3, 3) + ylim(-3, 3)
print(plot)




# CD4
#ggplot(metadata_downsampled[metadata_downsampled$level1 == "CD4", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = is_ref), size = 0.5) + 
#  theme_classic() + ggtitle("Query data mapped onto CD4 MDE") +
#  xlim(-2.5, 6) + ylim(-2.5,6)
#ggplot(metadata_downsampled[metadata_downsampled$level1 == "CD4", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = level2), size = 0.5) + scale_color_manual(values = mypal_level2) + 
#  theme_classic() + ggtitle("ImmgenT CD4 + query Annotation level 1") +
#  xlim(-2.5, 6) + ylim(-2.5,6)
#ggplot(metadata_query[metadata_query$level1 == "CD4", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = level2), size = 0.5) + scale_color_manual(values = mypal_level2) + 
#  theme_classic() + ggtitle("ImmgenT CD4 + query Annotation level 2") +
#  xlim(-2.5, 6) + ylim(-2.5,6)

if (length(rownames(metadata_query_plot[(metadata_query_plot$level1_final == "CD4" & grepl("CD4_", metadata_query_plot$level2_plot)), ])) > 1) { 
  plot <- ggplot() + geom_scattermore(mde_plot[mde_plot$level1 == "CD4", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "grey", size = 0.75, pixels = c(512, 512)) + 
    geom_point(metadata_query_plot[metadata_query_plot$level2_plot == "CD4", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "darkgray", size = 0.75) +
    geom_point(metadata_query_plot[(metadata_query_plot$level1_final == "CD4" & grepl("CD4_", metadata_query_plot$level2_plot)), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_plot), size = 0.75) +
    scale_color_manual(values = mypal_clusters) +
    scale_shape_manual(values = c("not classified - level2" = 17)) + 
    guides(color = guide_legend(override.aes = list(size = 2)),
           shape = guide_legend(override.aes = list(size = 2))) +
    theme_classic() + ggtitle(paste0(annotation, ": CD4: ImmgenT reference (grey) with query cells with ImmgenT annotation (colour)")) +
    xlim(-3.5, 3.5) + ylim(-3.5, 3.5)
print(plot)

}

#ggplot() +   geom_point(metadata_query[grep("CD4", metadata_query$level2_plot), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_scanvi_confidence), size = 0.75) +
#  scale_color_gradient(low = "orange", high = "black") +
#  theme_classic() + ggtitle("CD4 query Annotation level2 scanvi confidence") +
#  xlim(-3.5, 3.5) + ylim(-3.5, 3.5)


# CD8
#ggplot(metadata_downsampled[metadata_downsampled$level1 == "CD8", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = is_ref), size = 0.5) + 
#  theme_classic() + ggtitle("Query data mapped onto CD8 MDE") +
#  xlim(-3, 3) + ylim(-3,3)
#ggplot(metadata_downsampled[metadata_downsampled$level1 == "CD8", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = level2), size = 0.5) + scale_color_manual(values = mypal_level2) + 
#  theme_classic() + ggtitle("ImmgenT CD8 + query Annotation level 1") +
#  xlim(-3, 3) + ylim(-3,3)
#ggplot(metadata_query[metadata_query$level1 == "CD8", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = level2), size = 0.5) + scale_color_manual(values = mypal_level2) + 
#  theme_classic() + ggtitle("ImmgenT CD8 + query Annotation level 2") +
#  xlim(-3, 3) + ylim(-3,3)

#ggplot() + geom_point(metadata[metadata$level1_C_scANVI == "CD8", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "grey", size = 0.75) + 
#  geom_point(metadata_query[grep("CD8", metadata_query$level2_plot), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_plot), size = 0.75) +
#  scale_color_manual(values = mypal) + 
#  theme_classic() + ggtitle("ImmgenT CD8 + query Annotation level 2") +
#  xlim(-4, 4) + ylim(-4,4)

if (length(rownames(metadata_query_plot[(metadata_query_plot$level1_final == "CD8" & grepl("CD8_", metadata_query_plot$level2_plot)), ])) > 1) {
  plot <- ggplot() + geom_scattermore(mde_plot[mde_plot$level1 == "CD8", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "grey", size = 0.75, pixels = c(512, 512)) +
    geom_point(metadata_query_plot[metadata_query_plot$level2_plot == "CD8", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "darkgray", size = 0.75) +
    geom_point(metadata_query_plot[(metadata_query_plot$level1_final == "CD8" & grepl("CD8_", metadata_query_plot$level2_plot)), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_plot), size = 0.75) +
    scale_color_manual(values = mypal_clusters) +
    scale_shape_manual(values = c("not classified - level2" = 17)) +
    guides(color = guide_legend(override.aes = list(size = 2)),
           shape = guide_legend(override.aes = list(size = 2))) +
    theme_classic() + ggtitle(paste0(annotation, ": CD8: ImmgenT reference (grey) with query cells with ImmgenT annotation (colour)")) +
    xlim(-4, 4) + ylim(-4, 4)
print(plot)

}


#ggplot() +   geom_point(metadata_query[grep("CD8", metadata_query$level2_plot), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_scanvi_confidence), size = 0.75) +
#  scale_color_gradient(low = "orange", high = "black") +
#  theme_classic() + ggtitle("CD8 query Annotation level2 scanvi confidence") +
#  xlim(-4, 4) + ylim(-4, 4)


# Treg
#ggplot(metadata_downsampled[metadata_downsampled$level1 == "Treg", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = is_ref), size = 0.5) + 
#  theme_classic() + ggtitle("Query data mapped onto Treg MDE") +
#  xlim(-2.5, 3.5) + ylim(-2.5,3.5)
#ggplot(metadata_downsampled[metadata_downsampled$level1 == "Treg", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = level2), size = 0.5) + scale_color_manual(values = mypal_level2) + 
#  theme_classic() + ggtitle("ImmgenT Treg + query Annotation level 1") +
#  xlim(-2.5, 3.5) + ylim(-2.5,3.5)
#ggplot(metadata_query[metadata_query$level1 == "Treg", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = level2), size = 0.5) + scale_color_manual(values = mypal_level2) + 
#  theme_classic() + ggtitle("ImmgenT Treg + query Annotation level 2") +
#  xlim(-2.5, 3.5) + ylim(-2.5,3.5)
#ggplot() + geom_point(metadata[metadata$level1_C_scANVI == "Treg", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "grey", size = 0.75) + 
#  geom_point(metadata_query[grep("Treg", metadata_query$level2_plot), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_plot), size = 0.75) +
#  scale_color_manual(values = mypal) + 
#  theme_classic() + ggtitle("ImmgenT Treg + query Annotation level 2") +
#  xlim(-3.5, 4.5) + ylim(-3.5,4.5)

if (length(rownames(metadata_query_plot[(metadata_query_plot$level1_final == "Treg" & grepl("Treg_", metadata_query_plot$level2_plot)), ])) > 1) {
  plot <- ggplot() + geom_scattermore(mde_plot[mde_plot$level1 == "Treg", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "grey", size = 0.75, pixels = c(512, 512)) +
  geom_point(metadata_query_plot[metadata_query_plot$level2_plot == "Treg", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "darkgray", size = 0.75) +
  geom_point(metadata_query_plot[(metadata_query_plot$level1_final == "Treg" & grepl("Treg_", metadata_query_plot$level2_plot)), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_plot), size = 0.75) +
  scale_color_manual(values = mypal_clusters) +
  scale_shape_manual(values = c("not classified - level2" = 17)) +
  guides(color = guide_legend(override.aes = list(size = 2)),
         shape = guide_legend(override.aes = list(size = 2))) +
  theme_classic() + ggtitle(paste0(annotation, ": Treg: ImmgenT reference (grey) with query cells with ImmgenT annotation (colour)")) +
  xlim(-3.5, 4.5) + ylim(-3.5, 4.5)
print(plot)

}


#ggplot() +   geom_point(metadata_query[grep("Treg", metadata_query$level2_plot), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_scanvi_confidence), size = 0.75) +
#  scale_color_gradient(low = "orange", high = "black") +
#  theme_classic() + ggtitle("Treg query Annotation level2 scanvi confidence") +
#  xlim(-3.5, 4.5) + ylim(-3.5, 4.5)


# DN
#ggplot(metadata_downsampled[metadata_downsampled$level1 == "DN", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = is_ref), size = 0.5) + 
#  theme_classic() + ggtitle("Query data mapped onto whole MDE") +
#  xlim(-3, 3) + ylim(-3,3)
#ggplot(metadata_downsampled[metadata_downsampled$level1 == "DN", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = level2), size = 0.5) + scale_color_manual(values = mypal_level2) + 
#  theme_classic() + ggtitle("ImmgenT + query Annotation level 1") +
#  xlim(-3, 3) + ylim(-3,3)
#ggplot(metadata_query[metadata_query$level1 == "DN", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = level2), size = 0.5) + scale_color_manual(values = mypal_level2) + 
#  theme_classic() + ggtitle("ImmgenT + query Annotation level 2") +
#  xlim(-3, 3) + ylim(-3,3)
#ggplot() + geom_point(metadata[metadata$level1_C_scANVI == "DN", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "grey", size = 0.75) + 
#  geom_point(metadata_query[grep("DN", metadata_query$level2_plot), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_plot), size = 0.75) +
#  scale_color_manual(values = mypal) + 
#  theme_classic() + ggtitle("ImmgenT DN + query Annotation level 2") +
#  xlim(-4, 4) + ylim(-4,4)

if (length(rownames(metadata_query_plot[(metadata_query_plot$level1_final == "DN" & grepl("DN_", metadata_query_plot$level2_plot)), ])) > 1) {
  plot <- ggplot() + geom_scattermore(mde_plot[mde_plot$level1 == "DN", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "grey", size = 0.75, pixels = c(512, 512)) +
    geom_point(metadata_query_plot[metadata_query_plot$level2_plot == "DN", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "darkgray", size = 0.75) +
    geom_point(metadata_query_plot[(metadata_query_plot$level1_final == "DN" & grepl("DN_", metadata_query_plot$level2_plot)), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_plot), size = 0.75) +
    scale_color_manual(values = mypal_clusters) +
    scale_shape_manual(values = c("not classified - level2" = 17)) +
    guides(color = guide_legend(override.aes = list(size = 2)),
           shape = guide_legend(override.aes = list(size = 2))) +
    theme_classic() + ggtitle(paste0(annotation, ": DN: ImmgenT reference (grey) with query cells with ImmgenT annotation (colour)")) +
    xlim(-4, 4) + ylim(-4, 4)
print(plot)

}


#ggplot() +   geom_point(metadata_query[grep("DN", metadata_query$level2_plot), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_scanvi_confidence), size = 0.75) +
#  scale_color_gradient(low = "orange", high = "black") +
#  theme_classic() + ggtitle("DN query Annotation level2 scanvi confidence") +
#  xlim(-4, 4) + ylim(-4, 4)


# CD8aa
#ggplot(metadata_downsampled[metadata_downsampled$level1 == "CD8aa", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = is_ref), size = 0.5) + 
#  theme_classic() + ggtitle("Query data mapped onto CD8aa MDE") +
#  xlim(-3, 3) + ylim(-3,3)
#ggplot(metadata_downsampled[metadata_downsampled$level1 == "CD8aa", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = level2), size = 0.5) + scale_color_manual(values = mypal_level2) + 
#  theme_classic() + ggtitle("ImmgenT CD8aa + query Annotation level 1") +
#  xlim(-3, 3) + ylim(-3,3)
#ggplot(metadata_query[metadata_query$level1 == "CD8aa", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = level2), size = 0.5) + scale_color_manual(values = mypal_level2) + 
#  theme_classic() + ggtitle("ImmgenT CD8aa + query Annotation level 2") +
#  xlim(-3, 3) + ylim(-3,3)
#ggplot() + geom_point(metadata[metadata$level1_C_scANVI == "CD8aa", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "grey", size = 0.75) + 
#  geom_point(metadata_query[grep("CD8aa", metadata_query$level2_plot), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_plot), size = 0.75) +
#  scale_color_manual(values = mypal) + 
#  theme_classic() + ggtitle("ImmgenT CD8aa + query Annotation level 2") +
#  xlim(-4, 4) + ylim(-4,4)

if (length(rownames(metadata_query_plot[(metadata_query_plot$level1_final == "CD8aa" & grepl("CD8aa_", metadata_query_plot$level2_plot)), ])) > 1) {
  plot <- ggplot() + geom_scattermore(mde_plot[mde_plot$level1 == "CD8aa", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "grey", size = 0.75, pixels = c(512, 512)) +
    geom_point(metadata_query_plot[metadata_query_plot$level2_plot == "CD8aa", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "darkgray", size = 0.75) +
    geom_point(metadata_query_plot[(metadata_query_plot$level1_final == "CD8aa" & grepl("CD8aa_", metadata_query_plot$level2_plot)), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_plot), size = 0.75) +
    scale_color_manual(values = mypal_clusters) +
    scale_shape_manual(values = c("not classified - level2" = 17)) +
    guides(color = guide_legend(override.aes = list(size = 2)),
           shape = guide_legend(override.aes = list(size = 2))) +
    theme_classic() + ggtitle(paste0(annotation, ": CD8aa: ImmgenT reference (grey) with query cells with ImmgenT annotation (colour)")) +
    xlim(-4, 4) + ylim(-4, 4)
print(plot)

}


#ggplot() +   geom_point(metadata_query[grep("CD8aa", metadata_query$level2_plot), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_scanvi_confidence), size = 0.75) +
#  scale_color_gradient(low = "orange", high = "black") +
#  theme_classic() + ggtitle("CD8aa query Annotation level2 scanvi confidence") +
#  xlim(-4, 4) + ylim(-4, 4)


# nonconv
#ggplot(metadata_downsampled[metadata_downsampled$level1 == "nonconv", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = is_ref), size = 0.5) + 
#  theme_classic() + ggtitle("Query data mapped onto nonconv MDE") +
#  xlim(-3, 3) + ylim(-3,3)
#ggplot(metadata_downsampled[metadata_downsampled$level1 == "nonconv", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = level2), size = 0.5) + scale_color_manual(values = mypal_level2) + 
#  theme_classic() + ggtitle("ImmgenT nonconv + query Annotation level 1") +
#  xlim(-3, 3) + ylim(-3,3)
#ggplot(metadata_query[metadata_query$level1 == "nonconv", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = level2), size = 0.5) + scale_color_manual(values = mypal_level2) + 
#  theme_classic() + ggtitle("ImmgenT nonconv + query Annotation level 2") +
#  xlim(-3, 3) + ylim(-3,3)
#ggplot() + geom_point(metadata[metadata$level1_C_scANVI == "nonconv", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "grey", size = 0.75) + 
#  geom_point(metadata_query[grep("nonconv", metadata_query$level2_plot), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_plot), size = 0.75) +
#  scale_color_manual(values = mypal) + 
#  theme_classic() + ggtitle("ImmgenT nonconv + query Annotation level 2") +
#  xlim(-4, 4) + ylim(-4,4)

if (length(rownames(metadata_query_plot[(metadata_query_plot$level1_final == "nonconv" & grepl("nonconv_", metadata_query_plot$level2_plot)), ])) > 1) {
  plot <- ggplot() + geom_scattermore(mde_plot[mde_plot$level1 == "nonconv", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "grey", size = 0.75, pixels = c(512, 512)) +
    geom_point(metadata_query_plot[metadata_query_plot$level2_plot == "nonconv", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "darkgray", size = 0.75) +
    geom_point(metadata_query_plot[(metadata_query_plot$level1_final == "nonconv" & grepl("nonconv_", metadata_query_plot$level2_plot)), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_plot), size = 0.75) +
    scale_color_manual(values = mypal_clusters) +
    scale_shape_manual(values = c("not classified - level2" = 17)) +
    guides(color = guide_legend(override.aes = list(size = 2)),
           shape = guide_legend(override.aes = list(size = 2))) +
    theme_classic() + ggtitle(paste0(annotation, ": nonconv: ImmgenT reference (grey) with query cells with ImmgenT annotation (colour)")) +
    xlim(-4, 4) + ylim(-4, 4)
print(plot)

}


#ggplot() +   geom_point(metadata_query[grep("nonconv", metadata_query$level2_plot), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_scanvi_confidence), size = 0.75) +
#  scale_color_gradient(low = "orange", high = "black") +
#  theme_classic() + ggtitle("nonconv query Annotation level2 scanvi confidence") +
#  xlim(-4, 4) + ylim(-4, 4)


# gdT
#ggplot(metadata_downsampled[metadata_downsampled$level1 == "gdT", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = is_ref), size = 0.5) + 
#  theme_classic() + ggtitle("Query data mapped onto gdT MDE") +
#  xlim(-2, 3) + ylim(-2,3)
#ggplot(metadata_downsampled[metadata_downsampled$level1 == "gdT", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = level2), size = 0.5) + scale_color_manual(values = mypal_level2) + 
#  theme_classic() + ggtitle("ImmgenT gdT + query Annotation level 1") +
#  xlim(-2, 3) + ylim(-2,3)
#ggplot(metadata_query[metadata_query$level1 == "gdT", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = level2), size = 0.5) + scale_color_manual(values = mypal_level2) + 
#  theme_classic() + ggtitle("ImmgenT gdT + query Annotation level 2") +
#  xlim(-2, 3) + ylim(-2,3)
#ggplot() + geom_point(metadata[metadata$level1_C_scANVI == "gdT", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "grey", size = 0.75) + 
#  geom_point(metadata_query[grep("gdT", metadata_query$level2_plot), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_plot), size = 0.75) +
#  scale_color_manual(values = mypal) + 
#  theme_classic() + ggtitle("ImmgenT gdT + query Annotation level 2") +
#  xlim(-3, 4) + ylim(-3,4)

if (length(rownames(metadata_query_plot[(metadata_query_plot$level1_final == "gdT" & grepl("gdT_", metadata_query_plot$level2_plot)), ])) > 1) {
  plot <- ggplot() + geom_scattermore(mde_plot[mde_plot$level1 == "gdT", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "grey", size = 0.75, pixels = c(512, 512)) +
    geom_point(metadata_query_plot[metadata_query_plot$level2_plot == "gdT", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "darkgray", size = 0.75) +
    geom_point(metadata_query_plot[(metadata_query_plot$level1_final == "gdT" & grepl("gdT_", metadata_query_plot$level2_plot)), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_plot), size = 0.75) +
    scale_color_manual(values = mypal_clusters) +
    scale_shape_manual(values = c("not classified - level2" = 17)) +
    guides(color = guide_legend(override.aes = list(size = 2)),
           shape = guide_legend(override.aes = list(size = 2))) +
    theme_classic() + ggtitle(paste0(annotation, ": gdT: ImmgenT reference (grey) with query cells with ImmgenT annotation (colour)")) +
    xlim(-3, 4) + ylim(-3, 4)
print(plot)

}


#ggplot() +   geom_point(metadata_query[grep("gdT", metadata_query$level2_plot), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_scanvi_confidence), size = 0.75) +
#  scale_color_gradient(low = "orange", high = "black") +
#  theme_classic() + ggtitle("gdT query Annotation level2 scanvi confidence") +
#  xlim(-3, 4) + ylim(-3, 4)


# DP
#ggplot(metadata_downsampled[metadata_downsampled$level1 == "DP", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = is_ref), size = 0.5) + 
#  theme_classic() + ggtitle("Query data mapped onto DP MDE") +
#  xlim(-2.5, 2.5) + ylim(-2.5,2.5)
#ggplot(metadata_downsampled[metadata_downsampled$level1 == "DP", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = level2), size = 0.5) + scale_color_manual(values = mypal_level2) + 
#  theme_classic() + ggtitle("ImmgenT DP + query Annotation level 1") +
#  xlim(-2.5, 2.5) + ylim(-2.5,2.5)
#ggplot(metadata_query[metadata_query$level1 == "DP", ], aes(x = level2_MDE1, y = level2_MDE2)) + geom_point(aes(colour = level2), size = 0.5) + scale_color_manual(values = mypal_level2) + 
#  theme_classic() + ggtitle("ImmgenT DP + query Annotation level 2") +
#  xlim(-2.5, 2.5) + ylim(-2.5,2.5)
#ggplot() + geom_point(metadata[metadata$level1_C_scANVI == "DP", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "grey", size = 0.75) + 
#  geom_point(metadata_query[grep("DP", metadata_query$level2_plot), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_plot), size = 0.75) +
#  scale_color_manual(values = mypal) + 
#  theme_classic() + ggtitle("ImmgenT DP + query Annotation level 2") +
#  xlim(-3.5, 3.5) + ylim(-3.5, 3.5)

if (length(rownames(metadata_query_plot[(metadata_query_plot$level1_final == "DP" & grepl("DP_", metadata_query_plot$level2_plot)), ])) > 1) {
  plot <- ggplot() + geom_scattermore(mde_plot[mde_plot$level1 == "DP", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "grey", size = 0.75, pixels = c(512, 512)) +
    geom_point(metadata_query_plot[metadata_query_plot$level2_plot == "DP", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "darkgray", size = 0.75) +
    geom_point(metadata_query_plot[(metadata_query_plot$level1_final == "DP" & grepl("DP_", metadata_query_plot$level2_plot)), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_plot), size = 0.75) +
    scale_color_manual(values = mypal_clusters) +
    scale_shape_manual(values = c("not classified - level2" = 17)) +
    guides(color = guide_legend(override.aes = list(size = 2)),
           shape = guide_legend(override.aes = list(size = 2))) +
    theme_classic() + ggtitle(paste0(annotation, ": DP: ImmgenT reference (grey) with query cells with ImmgenT annotation (colour)")) +
    xlim(-3.5, 3.5) + ylim(-3.5, 3.5)
print(plot)

}


#ggplot() +   geom_point(metadata_query[grep("DP", metadata_query$level2_plot), ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_scanvi_confidence), size = 0.75) +
#  scale_color_gradient(low = "orange", high = "black") +
#  theme_classic() + ggtitle("DP query Annotation level2 scanvi confidence") +
#  xlim(-3.5, 3.5) + ylim(-3.5, 3.5)
}

dev.off()
#}

write.csv(metadata_query, paste0(prefix, "/metadata_query.csv"), quote = F, row.names = T, col.names = T)

## Save stats for level1 and level2
tmp = table(metadata_query$level1_final)
write.table(tmp, paste0(prefix, "/level1_final_table.txt"),quote = F, row.names = F, col.names = F, sep = "\t")
tmp = table(metadata_query$level1_assignment)
write.table(tmp, paste0(prefix, "/level1_assignment_table.txt"),quote = F, row.names = F, col.names = F, sep = "\t")
tmp = table(metadata_query$level2_plot)
write.table(tmp, paste0(prefix, "/level2_plot_table.txt"),quote = F, row.names = F, col.names = F, sep = "\t")
tmp = table(metadata_query$level2_final)
write.table(tmp, paste0(prefix, "/level2_final_table.txt"),quote = F, row.names = F, col.names = F, sep = "\t")

## Save MDEs
#MDE_save <- sort(unique(output_file$level1_C_scANVI))
#MDE_save <- MDE_save[!(MDE_save %in% c("unclear", "nonT", "remove"))]
#for (group in MDE_save) {
#  print(group)
#  vp_rows <- grep("VP", rownames(output_file))
#  group_rows <- grep(group, output_file$level2_plot)
#  
#  mde <- output_file[intersect(vp_rows,  group_rows), c("level2_plot", "level1_MDE1", "level1_MDE2", "level2_MDE1", "level2_MDE2")]
#  if (length(rownames(mde)) > 20) {
#    write.csv(mde, paste0(prefix, "/", group ,"_mde.csv"), quote = F, row.names = T, col.names = T)
#  }
#}

##CD4_mde <- output_file[grep("CD4", output_file$level2_plot), c("level2_plot", "level1_MDE1", "level1_MDE2", "level2_MDE1", "level2_MDE2")]
##write.csv(CD4_mde, paste0(prefix, "/CD4_mde.csv"), quote = F, row.names = T, col.names = T)

