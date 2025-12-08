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
mde_incremental_path = args[3]
prefix = args[4]
query_IGTHT_path = args[5]
predictions_output_file_not_classified = args[6]
annotation_column = args[7]
threshold = 0.85

output_file <- read.csv(output_file_path, row.names = 1)
mde_incremental <- read.csv(mde_incremental_path, row.names = 1)
print("no")
colnames(mde_incremental) <- c("level2_MDE1", "level2_MDE2")
print("yes")
#joint.bcs <- intersect(rownames(output_file), rownames(mde_incremental)) 
#output_file[joint.bcs, "level2_MDE1"] <- mde_incremental[joint.bcs, "level2_MDE1"]
#output_file[joint.bcs, "level2_MDE2"] <- mde_incremental[joint.bcs, "level2_MDE2"]

output_file$level2_MDE1 <- mde_incremental$level2_MDE1
output_file$level2_MDE2 <- mde_incremental$level2_MDE2


predictions_output_file_not_classified <- read.csv(predictions_output_file_not_classified, row.names = 1)
## Store originals
output_file$level2_C_scANVI_orig <- "NA"
output_file$level2_scanvi_confidence_orig <- "NA"


# Ensure both objects have rownames set
common_ids <- intersect(rownames(output_file), rownames(predictions_output_file_not_classified))

# Assign only to the common cells
output_file$level2_C_scANVI_orig <- output_file$level2_C_scANVI
output_file$level2_scanvi_confidence_orig <- output_file$level2_scanvi_confidence

## Replace with new for plotting
# Ensure both objects have rownames set
common_ids <- intersect(rownames(output_file), rownames(predictions_output_file_not_classified))

# Assign only to the common cells
output_file[common_ids, "level2_C_scANVI"] <- predictions_output_file_not_classified[common_ids, "level2_C_scANVI_not_classified"]
output_file[common_ids, "level2_scanvi_confidence"] <- predictions_output_file_not_classified[common_ids, "level2_scanvi_confidence_not_classified"]

#output_file$level2_C_scANVI[rownames(predictions_output_file_not_classified)] <- predictions_output_file_not_classified$level2_C_scANVI_not_classified
#output_file$level2_scanvi_confidence[rownames(predictions_output_file_not_classified)] <- predictions_output_file_not_classified$level2_scanvi_confidence_not_classified

# Wrangle output file
output_file$level2_final <- output_file$level2_C_scANVI
output_file$level2_final[output_file$level2_scanvi_confidence < threshold] <- "not classified"

# Create new dfs for plotting
#user_output_file <- read.csv(user_output_file_path, row.names = 1)
user_output_file <- output_file[!(grepl("IGT", rownames(output_file))), ]
query_IGTHT <- read.csv(query_IGTHT_path, row.names=1)
query_IGTHT <- query_IGTHT[rownames(query_IGTHT) %in% rownames(user_output_file), ]
metadata_query <- output_file[rownames(output_file) %in% rownames(user_output_file), ]
metadata_query$IGTHT <- query_IGTHT$IGTHT
metadata_query$annotation <- query_IGTHT$annotation
#common_idsx <- intersect(rownames(metadata_query), rownames(query_IGTHT))
#metadata_query$annotation[common_idsx] <- query_IGTHT[common_idsx, c(annotation_column)]
metadata <- output_file[!(rownames(output_file) %in% rownames(user_output_file)), ]

mypal_clusters <- setNames(mypal, sort(unique(metadata_query$level2_final)))

## Downsample for plotting
#library(tidyverse)
library(tidyr)
library(dplyr)
#metadata_downsampled <- metadata %>%
#  group_by(level2_C_scANVI) %>%
#  #slice_sample(n = 1000) %>%
#  #slice_sample(n = min(n(), 1000)) %>%  
#  group_modify(~ slice_sample(.x, n = min(1000, nrow(.x)))) %>%
#  ungroup()

#metadata_downsampled <- rbind(metadata_downsampled, metadata_query)


#mypal_level2 <- setNames(mypal, sort(unique(metadata$level1_final)))
#print(mypal_level2)
#names(mypal_level2)[13:(length(sort(unique(metadata$level2_C_scANVI)))+12)] <- sort(unique(metadata$level2_C_scANVI))

library(ggplot2)
library(scattermore)
## Sample plots for webpage
pdf(paste0(prefix, "/sample_split_not_classified_scattermore_fix_level2_relax_threshold_Sample_plots_ImmgenT_data_integration_webpage.pdf"), height = 10, width = 10)
# Dataset

plot_main <- ggplot() + geom_scattermore(metadata, mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "grey", size = 0.75, pixels = c(512, 512)) +
  geom_point(metadata_query[metadata_query$level2_final == "not classified", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "black", size = 0.75, shape = 17) +
  geom_point(metadata_query[metadata_query$level2_final != "not classified", ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_final), size = 0.75) +
  #geom_point(metadata_query[metadata_query$level1_final == "Not Classified", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "darkgrey", size = 0.25, shape = 17) +
  scale_color_manual(values = mypal_clusters) +
  scale_shape_manual(values = c("not classified" = 17)) +
   guides(color = guide_legend(override.aes = list(size = 2)),
          shape = guide_legend(override.aes = list(size = 2))) +
  theme_classic() + ggtitle("ImmgenT Reference (grey background) with query cells with ImmgenT annotation (colour) - CD8 MDE") +
  xlim(-3, 3) + ylim(-3, 3)

print(plot_main)

# Split by sample
for (annotation in unique(metadata_query$annotation)) {
	plot <- ggplot() + geom_scattermore(metadata, mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "grey", size = 0.75, pixels = c(512, 512)) +
		geom_point(metadata_query[metadata_query$level2_final == "not classified" & metadata_query$annotation == annotation, ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "black", size = 0.75, shape = 17) +
		geom_point(metadata_query[metadata_query$level2_final != "not classified" & metadata_query$annotation == annotation, ], mapping = aes(x = level2_MDE1, y = level2_MDE2, colour = level2_final), size = 0.75) +
		#geom_point(metadata_query[metadata_query$level1_final == "Not Classified", ], mapping = aes(x = level2_MDE1, y = level2_MDE2), colour = "darkgrey", size = 0.25, shape = 17) +
		scale_color_manual(values = mypal_clusters) +
		scale_shape_manual(values = c("not classified" = 17)) +
		guides(color = guide_legend(override.aes = list(size = 2)),
		       shape = guide_legend(override.aes = list(size = 2))) +
theme_classic() + ggtitle(paste0(annotation, ": ImmgenT Reference (grey background) with query cells with ImmgenT annotation (colour) - CD8 MDE")) +
xlim(-3, 3) + ylim(-3, 3)
print(plot)

}

dev.off()

write.csv(metadata_query, paste0(prefix, "/metadata_query_not_classified.csv"), quote = F, row.names = T, col.names = T)

tmp = table(metadata_query$level2_final)
write.table(tmp, paste0(prefix, "/level2_final_not_classified_table.txt"),quote = F, row.names = F, col.names = F, sep = "\t")

