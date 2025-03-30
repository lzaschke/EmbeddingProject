# ------------------ Load libraries ----------------
# Data Manipulation
library(tidyverse)  # Data science toolkit (ggplot2, dplyr, readr, etc.)
library(dplyr)      # Data manipulation
library(magrittr)   # Pipe operator (%>%)
library(readr)      # Read CSV and delimited files
# Visualization
library(ggplot2)    # Plotting
library(patchwork)  # Arrange multiple plots
library(pheatmap)   # Heatmaps
library(viridis)    # Colorblind-friendly palettes
library(cowplot)
# Single-Cell & RNA-Seq Analysis
library(Seurat)     # Single-cell RNA-seq
library(SeuratDisk) # Save Seurat objects
library(DoubletFinder) # remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DESeq2)     # Differential expression (bulk RNA-seq)
library(SingleR)    # Cell type annotation
library(SoupX)      # Remove ambient RNA
library(harmony)    # Data integration
# Genomic Analysis (Bioconductor)
library(BiocGenerics)         # Bioconductor core functions
library(S4Vectors)            # Data structures
library(GenomeInfoDb)         # Genomic coordinate management
library(SummarizedExperiment) # Omics data container
library(infercnv)             # Copy number variation
library(celldex)              # Cell reference datasets

# ------------------ Load postSoupX Data ----------------
O_RTK1_m_2a = readRDS("/data/cephfs-1/work/projects/aging-gbm/analysis/PostSoupX/O_RTK1_m_2a_seurat_postSoupX.rds")
O_RTK1_m_2b = readRDS("/data/cephfs-1/work/projects/aging-gbm/analysis/PostSoupX/O_RTK1_m_2b_seurat_postSoupX.rds")

sobj = merge(x = O_RTK1_m_2a, y = O_RTK1_m_2b, 
             add.cell.ids = c("O_RTK1_m_2a","O_RTK1_m_2b"))

# ------------------ QC ---------------------
sobj <- PercentageFeatureSet(sobj, "^MT-", col.name = "percent_mito") 
sobj <- PercentageFeatureSet(sobj, "^HB[^(P|E|S)]", col.name = "percent_hb")

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_hb")
VlnPlot(sobj, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3)
FeatureScatter(sobj, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = .5)

print(paste0("pre-filterung #cells = ", ncol(sobj))) # 27337
sobj <- subset(sobj, 
               subset = nFeature_RNA > 400 & 
                 nCount_RNA > 450 & 
                 percent_mito < 15 & 
                 percent_hb < 5)
print(paste0("post-filterung #cells = ", ncol(sobj))) # 27169
VlnPlot(sobj, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + NoLegend()

# ------------------ DoubletFinder ------------------
samp_split <- SplitObject(sobj, split.by = "orig.ident")
sample <- sample_list[[1]]

min_pc_cal <- function(sample){
  stdv <- sample[["pca"]]@stdev
  percent_stdv <- (stdv/sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1] # co1 is the first principal component where the cumulative percentage of variance explained is greater than 90%, and the contribution of that individual component is less than 5%.
  co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - # co2 identifies which consecutive principal components have a large drop of more than 10% (or 0.1) in their percentage contribution to variance.
                       percent_stdv[2:length(percent_stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min_pc <- min(co1, co2)
  return(min_pc)
}

# pre-processing
run_doubletfinder_custom <- function(sample){
  sample <- NormalizeData(sample)
  sample <- FindVariableFeatures(sample)
  sample <- ScaleData(sample)
  sample <- RunPCA(sample)
  # find minimum relevant #PC
  stdv <- sample[["pca"]]@stdev
  percent_stdv <- (stdv/sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1] # co1 is the first principal component where the cumulative percentage of variance explained is greater than 90%, and the contribution of that individual component is less than 5%.
  co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - # co2 identifies which consecutive principal components have a large drop of more than 10% (or 0.1) in their percentage contribution to variance.
                       percent_stdv[2:length(percent_stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min_pc <- min(co1, co2)
  # continue pre-processing
  sample <- RunUMAP(sample, dims = 1:min_pc)
  sample <- FindNeighbors(sample, dims = 1:min_pc)              
  sample <- FindClusters(sample, resolution = 0.1)
  
  # pK-Identification
  sweep_list <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE)   
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats) # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
  optimal.pk <- bcmvn %>% # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
  
  # Expected number of doublets
  # Homotypic doublet proportion estimate
  annotations <- sample@meta.data$seurat_clusters # use the clusters as the user-defined cell types
  homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets
  # Multiplet_rate calculation
  # 10X multiplet rates table: https://rpubs.com/kenneditodd/doublet_finder_example
  multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                    'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
                                    'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))
  multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(sample@meta.data)) %>% 
    dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
    dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
  print(paste('Setting multiplet rate to', multiplet_rate))
  nExp.poi <- round(multiplet_rate * nrow(sample@meta.data)) # multiply by number of cells to get the number of expected multiplets
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
  
  # run DoubletFinder
  sample <- doubletFinder(seu = sample, 
                          PCs = 1:min_pc, 
                          pK = optimal.pk, # the neighborhood size used to compute the number of artificial nearest neighbours
                          nExp = nExp.poi.adj) # number of expected real doublets
  colnames(sample@meta.data)[grepl('DF.classifications.*', colnames(sample@meta.data))] <- "doublet_finder" # change name of metadata column with Singlet/Doublet information
  
  double_finder_res <- sample@meta.data['doublet_finder'] # get the metadata column with singlet, doublet info
  double_finder_res <- rownames_to_column(double_finder_res, "row_names") # add the cell IDs as new column to be able to merge correctly
  return(double_finder_res)
}

samp_split <- lapply(samp_split, run_doubletfinder_custom) # get singlet/doublet assigned to each of the cell IDs (each element of the list is a different sample)
sobj_metadata <- data.frame(bind_rows(samp_split)) # merge to a single dataframe
sobj <- AddMetaData(sobj, sobj_metadata, col.name = 'doublet_finder')

# Check how doublets singlets differ in QC measures per sample.
VlnPlot(sobj, group.by = 'orig.ident', split.by = "doublet_finder",
        features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_hb"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')

# Get doublets per sample
doublets_summary <- sobj@meta.data %>% 
  group_by(orig.ident, doublet_finder) %>% 
  summarise(total_count = n(),.groups = 'drop') %>% as.data.frame() %>% ungroup() %>%
  group_by(orig.ident) %>%
  mutate(countT = sum(total_count)) %>%
  group_by(doublet_finder, .add = TRUE) %>%
  mutate(percent = paste0(round(100 * total_count/countT, 2),'%')) %>%
  dplyr::select(-countT)
print(doublets_summary)
write.table(doublets_summary, file = file.path("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/", paste0('doubletfinder_doublets_summary.txt')), quote = FALSE, row.names = FALSE, sep = '\t')

VlnPlot(sample, split.by = "doublet_finder",
        features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo", "percent_hb"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')

sobj_dblt <- subset(sobj, subset = doublet_finder == 'Singlet')

saveRDS(sobj, file = file.path("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/", 'seu_pre__dblt_filt.rds'))
saveRDS(sobj_dblt, file = file.path("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/", 'seu_dblt_filt.rds'))

sobj = readRDS("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/seu_dblt_filt.rds")

# ------------------ Identify malignant cells -------------------
# Filter high-count cells
sobj <- subset(sobj, 
               subset = nFeature_RNA < 4000 &
                 nCount_RNA < 10000)
samples <- SplitObject(sobj, split.by = "orig.ident")

for (sample in samples){
  # Add normal reference (oligodendrocytes)
  mg_od_original = readRDS("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/original_filtered_mg_od_Y_RTK2_m_1b.rds")
  od = subset(mg_od_original, subset = cell_type == "od")
  od_cm = as.matrix(GetAssayData(od))
  od_annotation = data.frame(
    cell_type = od$cell_type 
  )
  
  sobj_cm = as.matrix(GetAssayData(sample, slot = "counts"))
  sobj_annotation = data.frame(cell_type = rep("unknown", ncol(sobj_cm)))
  rownames(sobj_annotation) <- colnames(sobj_cm)
  
  combined_annotation <- rbind(od_annotation, sobj_annotation)
  genes_normal <- rownames(od_cm)
  genes_unknown <- rownames(sobj_cm)
  common_genes <- intersect(genes_normal, genes_unknown)
  normal_cells_data_aligned <- od_cm[common_genes, , drop = FALSE]
  unknown_cells_data_aligned <- sobj_cm[common_genes, , drop = FALSE]
  combined_cm <- cbind(normal_cells_data_aligned, unknown_cells_data_aligned)
  
  # Create OutputFolder
  sample_name = sample@meta.data$orig.ident[[1]]
  main_dir <- "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/inferCNV/"
  subfolder <- file.path(main_dir, sample_name)  # Create full path
  if (!dir.exists(subfolder)) {
    dir.create(subfolder, recursive = TRUE)  # 'recursive = TRUE' ensures parent dirs exist
  }
  
  # Run InferCNV
  file_path_gene_pos <- "/data/cephfs-1/work/projects/aging-gbm/analysis/inferCNV/gencode_v19_gene_pos.txt"
  infercnv_obj_full = CreateInfercnvObject(raw_counts_matrix = combined_cm,
                                           annotations_file = combined_annotation,
                                           gene_order_file= file_path_gene_pos,
                                           delim = "\t",
                                           ref_group_names = c("od"))
  
  infercnv_obj_full_run = infercnv::run(infercnv_obj_full,
                                        out_dir = subfolder,
                                        cutoff=0.1, #minimum reads that a gene must have on average across all cells not to be filtered out # SMARTseq: 1, 10x: 0.1
                                        cluster_by_groups = T,
                                        HMM = T,
                                        leiden_resolution = 0.003,
                                        #up_to_step = 15, #stop analysis after a certain step to allow reevaltuation
                                        num_threads = 20)
}
output_dir_full = "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/inferCNV_O_RTK1_m_2b_od_only_reference"
#output_dir_full = "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/inferCNV_O_RTK1_m_2a_od_only_reference_No.3"
# Load HMM-results from inferCNV folder (calculated on leiden subcluster level)
infercnv_obj <- readRDS(file.path(output_dir_full, "19_HMM_pred.Bayes_NetHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.infercnv_obj"))

# Prepare CNV matrix with prefix
cnv_matrix <- infercnv_obj@expr.data
colnames(cnv_matrix) <- paste0("O_RTK1_m_2b_", colnames(cnv_matrix))

# Load gene position info and filter for chromosomes 7 & 10
gene_info <- read.delim(file_path_gene_pos, header = FALSE, col.names = c("gene", "chr", "start", "end"))
chr_genes <- split(gene_info$gene, gene_info$chr)

chr7_genes <- gene_info$gene[gene_info$chr == "chr7"]
chr10_genes <- gene_info$gene[gene_info$chr == "chr10"]

chr7_genes_in_matrix <- intersect(chr7_genes, rownames(cnv_matrix))
chr10_genes_in_matrix <- intersect(chr10_genes, rownames(cnv_matrix))

chr7_expr <- colMeans(cnv_matrix[intersect(chr_genes$chr7, rownames(cnv_matrix)), ], na.rm = TRUE)
chr10_expr <- colMeans(cnv_matrix[intersect(chr_genes$chr10, rownames(cnv_matrix)), ], na.rm = TRUE)

# Identify malignant cells based on expression thresholds
gain_threshold <- 4  
loss_threshold <- 2  

is_malignant <- chr7_expr > gain_threshold & chr10_expr < loss_threshold
names(is_malignant) <- colnames(cnv_matrix)

# Map to Seurat object
sample = samples[[2]]
common_cells <- intersect(colnames(cnv_matrix), colnames(sample))
sample$is_malignant <- factor(replace(rep(NA, ncol(sample)), common_cells, ifelse(is_malignant[common_cells], "Malignant", "Non-Malignant")),
                              levels = c("Non-Malignant", "Malignant"))

samples[[2]] = sample

view(samples[[2]]@meta.data)

# ------------------ Annotation -------------------
# samples <- SplitObject(sobj, split.by = "orig.ident")

seurat_clustering <- function(sample, resolution_value){
  sample <- NormalizeData(sample)
  sample <- FindVariableFeatures(sample)
  sample <- ScaleData(sample)
  sample <- RunPCA(sample)
  # find minimum relevant #PC
  stdv <- sample[["pca"]]@stdev
  percent_stdv <- (stdv/sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1] # co1 is the first principal component where the cumulative percentage of variance explained is greater than 90%, and the contribution of that individual component is less than 5%.
  co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - # co2 identifies which consecutive principal components have a large drop of more than 10% (or 0.1) in their percentage contribution to variance.
                       percent_stdv[2:length(percent_stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min_pc <- min(co1, co2)
  # continue pre-processing
  sample <- FindNeighbors(sample, dims = 1:min_pc, reduction = "pca")              
  sample <- FindClusters(sample, resolution = resolution_value)
  sample <- RunUMAP(sample, dims = 1:min_pc, reduction = "pca")
  sample <- RunTSNE(sample, dims = 1:min_pc, reduction = "pca")
  return(sample)
}

samples <- lapply(samples, seurat_clustering, resolution_value = 2)

annotation <- function(sample, cell_markers, malignant_threshold, score_threshold){
  sample <- AddModuleScore(
    object = sample,
    features = cell_markers,
    name = "ModuleScore"
  )
  
  old_colnames <- grep("^ModuleScore", colnames(sample@meta.data), value = TRUE)
  names(old_colnames) <- names(cell_markers)
  
  new_colnames <- names(cell_markers)
  colnames(sample@meta.data)[match(old_colnames, colnames(sample@meta.data))] <- new_colnames 
  
  score_names = names(cell_markers)
  malignant_scores = score_names[1:8]
  non_malignant_scores = score_names[9:length(score_names)]
  
  clusters <- levels(sample$seurat_clusters)
  cluster_annotations <- setNames(vector("character", length(clusters)), clusters)
  
  for (cl in clusters){
    cells <- WhichCells(sample, idents = cl)
    valid_cells <- cells[!is.na(sample@meta.data[cells, "is_malignant"])] # remove NA
    if (length(valid_cells) > 0) {
      malignant_ratio <- mean(sample@meta.data[valid_cells, "is_malignant"] == "Malignant")
    } else {
      malignant_ratio <- 0  # Default to 0 if no valid cells
    }
    if (malignant_ratio > malignant_threshold) {
      scores = malignant_scores
    }
    else {
      scores = non_malignant_scores
    }
    avg_mod_scores <- colMeans(sample@meta.data[cells, scores, drop = FALSE], na.rm = TRUE)
    # Determine the max-scoring cell type
    max_score_index <- which.max(avg_mod_scores)
    max_score <- avg_mod_scores[max_score_index]
    if (max_score < score_threshold) {
      cluster_annotations[cl] <- "Unknown"  # If no score is high, assign "Unknown"
    } else {
      cluster_annotations[cl] <- scores[max_score_index]
    }
  }
  
  sample$Cluster_Annotations <- plyr::mapvalues(
    sample$seurat_clusters, 
    from = names(cluster_annotations), 
    to = cluster_annotations
  )
  
  return(sample)
}

cell_markers = read.csv("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/Annotation_modules.csv", sep = ";")
malignant_threshold = 0.7
score_threshold = 0.2

samples <- lapply(samples, annotation, 
               cell_markers = cell_markers, 
               malignant_threshold = malignant_threshold, 
               score_threshold = score_threshold)

sample = samples[[1]]

DimPlot(sample, group.by = "Cluster_Annotations")

saveRDS(samples, "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/samples_processed.rds" )
samples = readRDS("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/samples_processed.rds")
# Back to original sobj
sobj = readRDS("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/seu_dblt_filt.rds")

merged_metadata <- do.call(rbind, lapply(samples, function(x) {
  meta <- x@meta.data[, c("is_malignant","seurat_clusters","Cluster_Annotations"), drop = FALSE]
  meta$cell <- rownames(meta)  # Store cell names explicitly
  meta$is_malignant <- as.character(meta$is_malignant)  # Convert factor to character
  return(meta)
}))

matching_cells <- intersect(rownames(sobj@meta.data), merged_metadata$cell)

sobj@meta.data[matching_cells, "is_malignant"] <- merged_metadata$is_malignant[match(matching_cells, merged_metadata$cell)]
sobj@meta.data[matching_cells, "seurat_clusters"] <- merged_metadata$seurat_clusters[match(matching_cells, merged_metadata$cell)]
sobj@meta.data[matching_cells, "Cluster_Annotations"] <- merged_metadata$Cluster_Annotations[match(matching_cells, merged_metadata$cell)]

sobj <- subset(sobj, cells = rownames(sobj@meta.data)[!is.na(sobj@meta.data$Cluster_Annotations)])

saveRDS(sobj, "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/sobj_annotated.rds")
#sobj = readRDS("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/sobj_annotated.rds")

tester = samples[[1]]
FeaturePlot(tester, features = c("NCBI"), reduction = "umap")
# --------------- Export unintegrated ---------------------
spliter = SplitObject(sobj, split.by = "orig.ident")
export_data = spliter[[2]]
write.csv(GetAssayData(export_data, slot = "counts.O_RTK1_m_2b"), "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/O_RTK1_m_2b_expression_matrix.csv")
write.csv(export_data@meta.data, "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/O_RTK1_m_2b_metadata.csv")

view(export_data@meta.data)

# --------------- Sample Integration ---------------------- 
#saver = samples
samples <- SplitObject(sobj, split.by = "orig.ident")

samples <- lapply(samples, function(sample) {
  sample <- NormalizeData(sample) #log1p
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
  return(sample)
})

anchors <- FindIntegrationAnchors(object.list = samples, dims = 1:30)
integrated_sobj <- IntegrateData(anchorset = anchors, dims = 1:30)

saveRDS(integrated_sobj, "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/integrated_sobj_annotated.rds")
integrated_sobj = readRDS("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/integrated_sobj_annotated.rds")
# ------------------ Transfer to Python -------------- 
# Create unique cluster IDs
integrated_sobj$Sample_Cluster <- paste0(integrated_sobj$orig.ident, "_", integrated_sobj$seurat_clusters)
view(integrated_sobj@meta.data)
unique_clusters <- sort(unique(integrated_sobj$Sample_Cluster))  # Sort for consistency
cluster_map <- setNames(seq_along(unique_clusters), unique_clusters)

integrated_sobj$Unique_Cluster_ID <- as.character(cluster_map[integrated_sobj$Sample_Cluster])
integrated_sobj$Unique_Cluster_ID <- factor(integrated_sobj$Unique_Cluster_ID, levels = as.character(seq_along(unique_clusters)))

integrated_sobj@meta.data$Sample_Cluster <- NULL

view(integrated_sobj@meta.data)
# Transfer
cm = integrated_sobj@assays$integrated$data
meta_data = integrated_sobj@meta.data
view(meta_data)
write.csv(cm, "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/integrated_sobj_expression_matrix.csv")
write.csv(meta_data, "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/integrated_sobj_metadata.csv")

##### Draw Final Figures ##### 
DefaultAssay(integrated_sobj) <- "integrated"
integrated_sobj <- ScaleData(integrated_sobj)
integrated_sobj <- RunPCA(integrated_sobj)
min_pc <- min_pc_cal(integrated_sobj)
integrated_sobj <- RunUMAP(integrated_sobj, dims = 1:min_pc)
integrated_sobj <- FindNeighbors(integrated_sobj, dims = 1:min_pc)
integrated_sobj <- FindClusters(integrated_sobj, resolution = 2)



color_mapping <- c(
  "Astrocyte" = "#1f77b4",            # Deep blue 
  "BDM" = "#ff7f0e",                  # Bright orange 
  "Endothelial" = "#2ca02c",          # Strong green 
  "GB.AC" = "#d62728",                # Deep red
  "GB.G1.S" = "#9467bd",              # Purple
  "GB.G2.M" = "#8c564b",              # Brown 
  "GB.MES2" = "#e377c2",              # Pink
  "GB.NPC1" = "#7f7f7f",              # Gray
  "MG" = "#bcbd22",                   # Olive green 
  "Neurons_excit" = "#17becf",        # Cyan
  "Neurons_inhib" = "#aec7e8",        # Light blue
  "OPC" = "#ff9896",                  # Soft pink
  "Oligodendrocyte" = "#98df8a",      # Light green
  "Pericyte" = "#d62728",             # Red
  "Perivascular_fibroblast" = "#c49c94", # Beige
  "SMC" = "#f7b6d2",                  # Light pin
  "T_cell" = "#ffbb78"                # Light orange 
)

sample = integrated_sobj #samples[[2]] #

sample_filtered <- subset(sample, subset = Cluster_Annotations != "Unknown")

sample_filtered$Cluster_Annotations <- factor(sample_filtered$Cluster_Annotations, levels = names(color_mapping))

umap_plot <- DimPlot(sample_filtered, reduction = "umap", group.by = "Cluster_Annotations", pt.size = 0.01, label = T, label.size = 3, label.box = F) +
  #scale_color_manual(values = color_mapping) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),  # Remove axis text
    axis.ticks = element_blank(), # Remove axis ticks
    axis.title = element_blank(), # Remove axis titles
    panel.grid = element_blank()  # Remove grid lines
  )

print(umap_plot)
ggsave("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/O_RTK1_m_2b_seuratclusters.png", plot = umap_plot, width = 6, height = 5,dpi = 300)
ggsave("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/O_RTK1_m_2b_seuratclusters.pdf", plot = umap_plot, width = 6, height = 5)






clusters <- integrated_sobj$integrated_snn_res.2

annotations <- integrated_sobj$Cluster_Annotations

# Create a data frame with both
annotation_df <- data.frame(
  cluster = clusters,
  cell_type = annotations
)

# Find the most frequent annotation for each cluster
cluster_to_name <- list()
for (cluster_id in unique(annotation_df$cluster)) {
  # Get annotations for this cluster
  cluster_cells <- annotation_df[annotation_df$cluster == cluster_id, "cell_type"]
  
  # Find the most common annotation
  freq_table <- table(cluster_cells)
  most_common <- names(freq_table)[which.max(freq_table)]
  
  # Store in our list
  cluster_to_name[[as.character(cluster_id)]] <- most_common
}

cluster_mapping_df <- data.frame(
  cluster_id = names(cluster_to_name),
  cell_type = unlist(cluster_to_name)
)

# Print the mapping
print(cluster_mapping_df)















sample$Cluster_Annotations <- factor(sample$Cluster_Annotations, levels = names(color_mapping))

sample_filtered = 
  
DimPlot(sample, reduction = "umap", group.by = "Cluster_Annotations") +
scale_color_manual(values = color_mapping) 


DimPlot(sample, reduction = "umap", group.by = "Cluster_Annotations")






test = integrated_sobj

DefaultAssay(integrated_sobj) <- "integrated"

test <- ScaleData(test)
test <- RunPCA(test)
test <- RunUMAP(test, dims = 1:17)
DimPlot(test, group.by = "is_malignant", label = T)

integrated_sobj <- ScaleData(integrated_sobj)
integrated_sobj <- RunPCA(integrated_sobj)
min_pc <- min_pc_cal(integrated_sobj)
integrated_sobj <- RunUMAP(integrated_sobj, dims = 1:min_pc)
integrated_sobj <- FindNeighbors(integrated_sobj, dims = 1:min_pc)
integrated_sobj <- FindClusters(integrated_sobj, resolution = 1)

DimPlot(integrated_sobj, group.by = "orig.ident")
DimPlot(integrated_sobj, group.by = "Cluster_Annotations", label = T)
DimPlot(integrated_sobj, group.by = "seurat_clusters", label = T)

saveRDS(integrated_sobj, "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/integrated_sobj_annotated.rds")
#integrated_sobj  <- readRDS("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/integrated_sobj_annotated.rds")

# Re-annotation of integrated object

test <- subset(integrated_sobj, cells = rownames(integrated_sobj@meta.data)[!(sobj@meta.data$Cluster_Annotations == "Unknown")])

test <- AddModuleScore(
  object = test,
  features = cell_markers,
  name = "ModuleScore"
)
view(integrated_sobj@assays$integrated$data)
head(rownames(integrated_sobj@assays$integrated$data))
test <- annotation(test, 
                  cell_markers = cell_markers, 
                  malignant_threshold = malignant_threshold, 
                  score_threshold = score_threshold)














DimPlot(test, group.by = "Cluster_Annotations")

view(sample@meta.data)

cell_markers = read.csv("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/Annotation_modules.csv", sep = ";")



sample <- AddModuleScore(
  object = sample,
  features = cell_markers,
  name = "ModuleScore"
)

old_colnames <- grep("^ModuleScore", colnames(sample@meta.data), value = TRUE)
names(old_colnames) <- names(cell_markers)

new_colnames <- names(cell_markers)
colnames(sample@meta.data)[match(old_colnames, colnames(sample@meta.data))] <- new_colnames 

score_names = names(cell_markers)
malignant_scores = score_names[1:8]
non_malignant_scores = score_names[9:length(score_names)]

malignant_threshold = 0.6
score_threshold = 0.2

clusters <- levels(sample$seurat_clusters)
cluster_annotations <- setNames(vector("character", length(clusters)), clusters)

for (cl in clusters){
  cells <- WhichCells(sample, idents = cl)
  valid_cells <- cells[!is.na(sample@meta.data[cells, "is_malignant"])] # remove NA
  malignant_ratio <- mean(sample@meta.data[valid_cells, "is_malignant"] == "Malignant")
  if (malignant_ratio > malignant_threshold) {
    scores = malignant_scores
  }
  else {
    scores = non_malignant_scores
  }
  avg_mod_scores <- colMeans(sample@meta.data[cells, scores, drop = FALSE], na.rm = TRUE)
  # Determine the max-scoring cell type
  max_score_index <- which.max(avg_mod_scores)
  max_score <- avg_mod_scores[max_score_index]
  if (max_score < score_threshold) {
    cluster_annotations[cl] <- "Unknown"  # If no score is high, assign "Unknown"
  } else {
    cluster_annotations[cl] <- scores[max_score_index]
  }
}

sample$Cluster_Annotations <- plyr::mapvalues(
  sample$seurat_clusters, 
  from = names(cluster_annotations), 
  to = cluster_annotations
)

DimPlot(sample, group.by ="Cluster_Annotations")


view(sample@meta.data)

print(avg_mod)

print(non_malignant_scores)
FeaturePlot(sample, features = score_names)

marker_list <- list(
  GBM = c("EGFR", "PDGFRA", "CD44", "CHI3L1", "SOX2", "OLIG2", "MKI67"),
  Astrocyte = c("TPD52L1", "ALDOC", "CLDN10", "ETNPPL"),
  Oligodendrocyte = c("PLP1", "MOG", "MOBP", "CLDN11"),
  Inhibitory_Neurons = c("GAD2", "GAD1", "SLC32A1"),
  Excitatory_Neurons = c("SLC17A7", "NRGN", "NEUROD2", "NEUROD6"),
  Microglia = c("CD74", "RGS1", "FCER1G", "LAPTM5", "CCL3", "TMEM119"),
  Macrophage = c("CD163", "MS4A6A", "PLAUR", "CD14", "CCL3"),
  Tcell = c("CD2", "CD3D", "CD3G"),
  Endothelial = c("ESAM", "PECAM1", "TBX3"),
  Pericytes = c("PDGFRB", "MYL9", "RGS5", "KCNJ8", "DES")
)

add_module_score <- function(sample){
  sample <- AddModuleScore(
    object = sample,
    features = marker_list,
    name = "ModuleScore"
  )
  
  # Rename metadata columns
  old_colnames <- grep("^ModuleScore", colnames(sample@meta.data), value = TRUE)
  names(old_colnames) <- names(marker_list) 
  new_colnames <- paste0("Score_", names(marker_list))
  
  colnames(sample@meta.data)[match(old_colnames, colnames(sample@meta.data))] <- new_colnames 
  return(sample)
}

samples <- lapply(samples, add_module_score)

score_names <- c("Score_GBM", "Score_Astrocyte", "Score_Oligodendrocyte", 
                 "Score_Inhibitory_Neurons", "Score_Excitatory_Neurons", 
                 "Score_Microglia", "Score_Macrophage", "Score_Tcell", 
                 "Score_Endothelial", "Score_Pericytes")

sample = samples[[1]]

VlnPlot(sample, features = score_names, pt.size = 0, ncol = 3)
FeaturePlot(sample, features = score_names, reduction = "umap", ncol = 3)

## Semi-automated annotation ## 
automated_annotation <- function(sample, malignant_threshold = 0.7, score_threshold = 0.2) {
  clusters <- levels(sample)  # Get cluster levels
  cluster_annotations <- setNames(vector("character", length(clusters)), clusters)
  
  for (cl in clusters) {
    cells <- WhichCells(sample, idents = cl)
    
    # Remove NA values before calculating the malignant ratio
    valid_cells <- cells[!is.na(sample@meta.data[cells, "is_malignant"])]
    
    if (length(valid_cells) > 0) {
      malignant_ratio <- mean(sample@meta.data[valid_cells, "is_malignant"] == "Malignant")
    } else {
      malignant_ratio <- 0  # Default to 0 if no valid cells
    }
    
    if (malignant_ratio > malignant_threshold) {
      cluster_annotations[cl] <- "Malignant"
      next  # SKIP the rest of the loop for this cluster (prevent overwriting)
    }
  }
  
  # Process non-malignant clusters separately
  for (cl in clusters) {
    if (cluster_annotations[cl] == "Malignant") next  # Ensure malignant clusters remain unchanged
    
    cells <- WhichCells(sample, idents = cl)
    
    # Compute mean module scores
    avg_mod_scores <- colMeans(sample@meta.data[cells, grep("Score", colnames(sample@meta.data))], na.rm = TRUE)
    
    # Determine the max-scoring cell type
    max_score_index <- which.max(avg_mod_scores)
    max_score <- avg_mod_scores[max_score_index]
    
    if (max_score < score_threshold) {
      cluster_annotations[cl] <- "Unknown"  # If no score is high, assign "Unknown"
    } else {
      cluster_annotations[cl] <- names(marker_list)[max_score_index]
    }
  }
  
  print(cluster_annotations)  # Check if plausible
  
  # Assign new annotations to metadata
  sample$Cluster_Annotations <- plyr::mapvalues(
    sample$seurat_clusters, 
    from = names(cluster_annotations), 
    to = cluster_annotations
  )
  
  return(sample)
}

samples <- lapply(samples,automated_annotation)

sample = samples[[1]]

DimPlot(sample, reduction = "umap", group.by = "Cluster_Annotations", label = FALSE)
DimPlot(sample, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(sample, reduction = "umap", group.by = "is_malignant", label = TRUE)

# ------------------ Further annotate malignant cells (Neftel et al.) -------------------
GBM_features = read.csv("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/metamodules.csv", sep = ";")

sub_annotate_malignant <- function(sample, GBM_features, resolution = 0.5, score_threshold = 0.2) {
  # Subset malignant cells
  malignant_cells <- WhichCells(sample, expression = Cluster_Annotations == "Malignant")
  if (length(malignant_cells) == 0) {
    sample$Final_Annotation <- sample$Cluster_Annotations
    return(sample)
  }
  
  malignant_subset <- subset(sample, cells = malignant_cells)
  
  # Re-cluster malignant cells
  malignant_subset <- subset(sample, cells = malignant_cells)
  malignant_subset <- FindVariableFeatures(malignant_subset)
  malignant_subset <- ScaleData(malignant_subset)
  malignant_subset <- RunPCA(malignant_subset)
  malignant_subset <- FindNeighbors(malignant_subset, dims = 1:10)
  malignant_subset <- FindClusters(malignant_subset, resolution = 0.5)
  malignant_subset <- RunUMAP(malignant_subset, dims = 1:10)
  
  # Calculate Module scores per cell with proper column names
  malignant_subset <- AddModuleScore(
    object = malignant_subset,
    features = GBM_features,
    name = "SubtypeScore"
  )
  
  old_colnames <- grep("^SubtypeScore", colnames(malignant_subset@meta.data), value = TRUE)
  names(old_colnames) <- names(GBM_features)
  
  new_colnames <- paste0("SubtypeScore_", names(GBM_features))
  colnames(malignant_subset@meta.data)[match(old_colnames, colnames(malignant_subset@meta.data))] <- new_colnames 
  
  # Assign subtypes per malignant cluster
  cluster_annotations <- setNames(rep("GB.Unknown", length(levels(malignant_subset$seurat_clusters))),
                                  levels(malignant_subset$seurat_clusters))
  
  for (cl in levels(malignant_subset$seurat_clusters)) {
    cells <- WhichCells(malignant_subset, idents = cl)
    
    # Ensure only valid columns are used in colMeans()
    valid_scores <- intersect(new_colnames, colnames(malignant_subset@meta.data))
    
    if (length(valid_scores) > 0) {
      avg_mod_scores <- colMeans(malignant_subset@meta.data[cells, valid_scores, drop = FALSE], na.rm = TRUE)
      
      max_score_index <- which.max(avg_mod_scores)
      max_score <- avg_mod_scores[max_score_index]
      
      if (max_score >= score_threshold) {
        cluster_annotations[cl] <- names(GBM_features)[max_score_index]  # Use clean names (no prefix)
      }
    }
  }
  
  # Create a new metadata column in malignant_subset
  malignant_subset$Malignant_Subtype <- plyr::mapvalues(
    malignant_subset$seurat_clusters,
    from = names(cluster_annotations),
    to = cluster_annotations
  )
  
  # Transfer scores & annotations back to full Seurat object
  sample$Malignant_Subtype <- NA  
  sample$Malignant_Subtype[malignant_cells] <- as.character(malignant_subset$Malignant_Subtype)
  
  for (score in valid_scores) {
    sample@meta.data[[score]] <- NA  # Initialize in main object
    sample@meta.data[malignant_cells, score] <- malignant_subset@meta.data[malignant_cells, score]
  }
  
  sample$Final_Annotation <- ifelse(
    sample$Cluster_Annotations == "Malignant" & !is.na(sample$Malignant_Subtype), 
    sample$Malignant_Subtype, 
    as.character(sample$Cluster_Annotations)
  )
  
  return(sample)
}

samples <- lapply(samples, sub_annotate_malignant, GBM_features = GBM_features)

sample = samples[[2]]

VlnPlot(sample, features = names(GBM_features), pt.size = 0, group.by = "Final_Annotation")
FeaturePlot(sample, features = names(GBM_features), reduction = "umap", ncol = 3)

DimPlot(sample, reduction ="umap", group.by = "Malignant_Subtype")
DimPlot(sample, reduction ="umap", group.by = "Final_Annotation")
DimPlot(sample, reduction ="umap", group.by = "seurat_clusters", label = T)

OPC_markers = c("DSCAM", "LHFPL3", "PCDH15", "TNR", "CSMD1", "OPCML", "KCNIP4", "LRRC4C", 
                 "FGF14", "LRRTM4", "KCND2", "CA10", "SNTG1", "LUZP2", "CSMD3", "NXPH1", 
                 "MMP16", "GRM7", "ATRNL1", "VCAN", "MDGA2", "FGF12", "MEGF11", "XYLT1", 
                 "BRINP3", "GRIK1", "PTPRZ1", "DCC", "GRIK2", "NRXN1", "DPP6", "SEZ6L", 
                 "SMOC1", "SLC35F1", "NLGN1", "SGCZ", "CNTNAP5", "GRM5", "HS6ST3", "NRCAM", 
                 "AC004852.2", "GALNT13", "TMEM132D", "PDZRN4", "PCDH7", "NLGN4X", "CHL1", 
                 "GRID2", "ARPP21", "SOX6", "XKR4", "RALYL", "FAM19A1", "GRIA2", "SORCS3", 
                 "SCN1A", "DGKG", "GRIA4", "COL11A1", "LRP1B", "LRRC7", "LINC01322", "CHST11", 
                 "LRRTM3", "CACNA1A", "SHISA9", "LINC00511", "FAM155A", "PTPRT", "ALK", "CSMD2", 
                 "DLGAP1", "KHDRBS3", "KAZN", "AL445250.1", "HECW1", "SLC8A1", "NOVA1", 
                 "SLC24A3", "CDH13", "FAM19A2", "STK32A", "LINC02588", "OPHN1", "CDH10", 
                 "RGS7", "SORCS1", "TRIM9", "ZFPM2", "NELL1", "ADAMTS17", "CNTN3", "KCNIP1", 
                 "DAB1", "KIF26B", "LRRN1", "KIAA1217", "NLGN4Y", "ADARB2", "UNC80", "SHC3", 
                 "MIR3681HG", "MYT1", "MAP2", "PLPPR1", "SEMA5A", "PLPP4", "PDZD2", "ADGRL3", 
                 "ZNF462", "EPN2", "POU6F2", "NEGR1", "SCN3A", "NTNG1", "RNF150", "KAT2B", 
                 "CHST9", "STK32B", "GRIA3", "SLC44A5", "RIMS2", "ERC2", "GSG1L", "TMEM108", 
                 "TMEM132C", "MEG8", "TMEM132B", "NOL4", "FERMT1", "RIT2", "DGKI", "STXBP5L", 
                 "SCN9A", "AC092691.1", "ASTN1", "RAB31", "COL9A1", "CDH18", "C1orf21", 
                 "CNTN1", "TNK2", "SGCD", "MYT1L", "LRFN5", "TOX", "ANKS1B", "RIMS1", 
                 "ZDHHC14", "MACROD2", "AL392023.2", "SCD5", "SEMA3E", "ESRRG", "ITGA8", 
                 "AL589740.1", "LSAMP", "KCNMA1", "ADGRB3", "TMEM163", "RBFOX1", "ERBB4", 
                 "CSRNP3", "UNC79", "GLCCI1", "ATP13A4", "SLC2A13", "TENM1", "KCTD16", 
                 "MYO16", "CACNG4", "SNAP91", "PRKG2", "IL1RAP", "DNER", "THRB", "GRIN2B", 
                 "PID1", "KIF13A", "BCAS1", "U91319.1", "ZNF804A")

FeaturePlot(sample, features="OPC_Score1")
sample <- AddModuleScore(sample, features = list(OPC_markers), name = "OPC_Score")

sample = samples[[1]]


# Plausibility check
markers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top5_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

view(top5_markers)

# Sub-annotate clusters
## immune cells
FeaturePlot(sample, features= c("CD163","MRC1","MS4A7","TGFBI"))
FeaturePlot(sample, features= c("PTPRC"))
FeaturePlot(sample, features= c("EGFR"))

# Transfer Annotation results to original seurat_obj
sobj = readRDS("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/seu_dblt_filt.rds")

merged_metadata <- do.call(rbind, lapply(samples, function(x) {
  meta <- x@meta.data[, c("Final_Annotation", "is_malignant","seurat_clusters"), drop = FALSE]
  meta$cell <- rownames(meta)  # Store cell names explicitly
  meta$is_malignant <- as.character(meta$is_malignant)  # Convert factor to character
  return(meta)
}))

matching_cells <- intersect(rownames(sobj@meta.data), merged_metadata$cell)


sobj@meta.data[matching_cells, "Final_Annotation"] <- merged_metadata$Final_Annotation[match(matching_cells, merged_metadata$cell)]
sobj@meta.data[matching_cells, "is_malignant"] <- merged_metadata$is_malignant[match(matching_cells, merged_metadata$cell)]
#sobj@meta.data[matching_cells, "seurat_clusters"] <- merged_metadata$seurat_clusters[match(matching_cells, merged_metadata$cell)]


view(sobj@meta.data)
#saver = samples



# ------------------ Data Export ------------------------
spliter = SplitObject(sobj, split.by = "orig.ident")
export_data = spliter[[1]]
write.csv(GetAssayData(export_data, slot = "counts.O_RTK1_m_2a"), "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/O_RTK1_m_2a_expression_matrix_joined.csv")
write.csv(export_data@meta.data, "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/O_RTK1_m_2a_metadata_joined.csv")
view(export_data@meta.data)


view(export_data@assays$RNA$counts.O_RTK1_m_2a)




# ------------------ Harmony data integration ---------------
sobj <- subset(sobj, 
               subset = nFeature_RNA < 4000 &
                 nCount_RNA < 10000)
samples <- SplitObject(sobj, split.by = "orig.ident")

# Embeddings
sobj <- NormalizeData(sobj)
sobj <- FindVariableFeatures(sobj)
sobj <- ScaleData(sobj)
sobj <- RunPCA(sobj)
# find minimum relevant #PC
stdv <- sobj[["pca"]]@stdev
percent_stdv <- (stdv/sum(stdv)) * 100
cumulative <- cumsum(percent_stdv)
co1 <- which(cumulative > 90 & percent_stdv < 5)[1] # co1 is the first principal component where the cumulative percentage of variance explained is greater than 90%, and the contribution of that individual component is less than 5%.
co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - # co2 identifies which consecutive principal components have a large drop of more than 10% (or 0.1) in their percentage contribution to variance.
                     percent_stdv[2:length(percent_stdv)]) > 0.1), 
            decreasing = T)[1] + 1
min_pc <- min(co1, co2)
# continue pre-processing
sobj <- FindNeighbors(sobj, dims = 1:min_pc, reduction = "pca")              
sobj <- FindClusters(sobj, resolution = 0.1, cluster.name = "unintegratedCluster")
sobj <- RunUMAP(sobj, dims = 1:min_pc, reduction = "pca", reduction.name = "umap.unintegrated")

DimPlot(sobj, reduction = "umap.unintegrated", group.by = "orig.ident")

sobj_int <- IntegrateLayers(object=sobj,
                            method = HarmonyIntegration,
                            orig.reduction = "pca",
                            new.reduction = "harmony",
                            verbose = T)
#sobj_int[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
sobj_int <- FindNeighbors(sobj_int, reduction = "harmony", dims = 1:min_pc)
sobj_int <- FindNeighbors(sobj_int, reduction = "harmony", dims = 1:min_pc)
sobj_int <- FindClusters(sobj_int, resolution = 1, cluster.name = "harmonyCluster")
sobj_int <- RunUMAP(sobj_int, dims = 1:min_pc, reduction = "harmony", reduction.name ="umap.harmony")
sobj_int <- RunTSNE(sobj_int, dims = 1:min_pc, reduction = "harmony", reduction.name ="tsne.harmony")
Idents(sobj_int) <- "orig.ident"
view(sobj_int@meta.data)
unique(sobj_int$orig.ident)
DimPlot(sobj_int, reduction = "umap.harmony", group.by = "orig.ident", 
        cells.highlight = list("O_RTK1_m_2a" = WhichCells(sobj_int, ident = "O_RTK1_m_2a")),  
        cols.highlight = c("red"),  
        sizes.highlight = 0.01,  
        label = F
)

DimPlot(sobj_int, reduction = "umap.harmony", group.by = "harmonyCluster", label=T)
DimPlot(sobj_int, reduction = "umap.harmony", group.by = "Final_Annotation", label=T)
DimPlot(sobj_int, reduction = "umap.harmony", group.by = "orig.ident", label=T)

# ------------------ Transfer to Python ----------------- 
#saveRDS(sobj_int, "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/sobj_int_annotated.rds")
#sobj_int = readRDS("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/sobj_int_annotated.rds")

sobj_int_layJoint = sobj_int
sobj_int_layJoint[["RNA"]] <- JoinLayers(sobj_int_layJoint[["RNA"]])

write.csv(GetAssayData(sobj_int_layJoint, slot = "data"), "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/expression_matrix_joined.csv")
write.csv(sobj_int_layJoint@meta.data, "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/metadata_joined.csv")
write.csv(Embeddings(sobj_int_layJoint, reduction = "harmony"), "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/harmony_pca_embeddings.csv")

#write.csv(GetAssayData(sobj_int, slot = "data"), "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/expression_matrix_joint.csv")

##### Alternative export ######
export_data = saver[[1]]
write.csv(GetAssayData(export_data, slot = "data"), "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/expression_matrix_joined.csv")
write.csv(export_data@meta.data, "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/metadata_joined.csv")















## Initially automated 
# Reference to be used 
ref <- celldex::HumanPrimaryCellAtlasData()
# Get counts data
sobj_int_layJoint = sobj_int
sobj_int_layJoint[["RNA"]] <- JoinLayers(sobj_int_layJoint[["RNA"]])
sobj_counts <- GetAssayData(sobj_int_layJoint, layer = "counts")
# SingleR
pred <- SingleR(test = sobj_counts,
                ref = ref,
                labels = ref$label.main)

sobj_int$singleR.labels <- pred$labels[match(rownames(sobj_int@meta.data), rownames(pred))]
DimPlot(sobj_int, reduction = "umap.harmony", group.by = "singleR.labels", label=T)

# inferCNV
mg_od <- readRDS("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/filtered_mg_od_Y_RTK2_m_1b.rds")
mg_od_counts <- as.matrix(GetAssayData(mg_od, slot = "counts"))





# Scevan

library(devtools)
library(SCEVAN)
library(doParallel)
library(foreach)

scevan_results <- pipelineCNA(sobj_counts)
saveRDS(sobj_int_layJoint, "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/sobj_int_layJoint.rds")
sobj_int_layJoint <- AddMetaData(sobj_int_layJoint, metadata = scevan_results$class, col.name = "scevan_results")
DimPlot(sobj_int_layJoint, reduction="umap.harmony", group.by="scevan_results")
view(sobj_int_layJoint@meta.data)
plot(DimPlot(sobj_int_layJoint, reduction = "umap.harmony", group.by = "orig.ident", label=T))

## InferCNV
# Pre-annotation
cell_annotation <- ifelse(
  sobj_int$orig.ident == "O_RTK1_m_2b" & sobj_int$seurat_clusters == 2, 
  "Normal", 
  "Unknown"
)
annotation_df <- data.frame(annotation = cell_annotation)
rownames(annotation_df) <- colnames(sobj_int)
# Gene-order file
file_path_gene_pos <- "/data/cephfs-1/work/projects/aging-gbm/analysis/inferCNV/gencode_v19_gene_pos.txt"
# Run
infercnv_obj_full = CreateInfercnvObject(raw_counts_matrix = as.matrix(GetAssayData(sobj_int_layJoint, layer = "counts")),
                                         annotations_file = annotation_df,
                                         gene_order_file= file_path_gene_pos,
                                         delim = "\t",
                                         ref_group_names = c("Normal"))
output_dir_full = "/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/inferCNV5"
infercnv_obj_full_run = infercnv::run(infercnv_obj_full,
                                      out_dir = output_dir_full,
                                      cutoff=0.1,
                                      window_length = 201,
                                      cluster_by_groups = T,
                                      HMM = T, # Hidden Markov Model: segmentation and prediction analysis that will generate clear cut regions with CNVs and assign them specific states of change
                                      output_format = "pdf",
                                      #denoise = T, #filtering of the residual expression which is the signal at the end of the infercnv analysis outside HMM
                                      #leiden_resolution = 0.00005, # resolution parameter for the leiden clustering, lower = less clusters, std= 0.05
                                      hclust_method='ward.D2',
                                      analysis_mode='subclusters',
                                      write_expr_matrix = T
                                      #up_to_step = 15, #stop analysis after a certain step to allow reevaltuation
                                      #debug = TRUE
                                      #num_threads = 20
)
infercnv::plot_cnv(infercnv_obj_full,
                   out_dir=output_dir_full,
                   cluster_by_groups=T,
                   color_safe_pal=FALSE,
                   x.center=1,
                   x.range=c(0.5,1.5),
                   title="inferCNV",
                   obs_title="Observations (Cells)",
                   ref_title="References (Cells)",
                   output_filename="heatmap",
                   hclust_method = "ward.D2",
                   write_expr_matrix = TRUE
)


a = as.matrix(infercnv_obj_full_run@tumor_subclusters$subclusters$Unknown$Unknown_s1)
a[,1] = "Malignant"
b = as.matrix(infercnv_obj_full_run@tumor_subclusters$subclusters$Unknown$Unknown_s2)
b[,1] = "Non-malignant"
c = as.matrix(infercnv_obj_full_run@tumor_subclusters$subclusters$Unknown$Unknown_s3)
c[,1] = "Non-malignant"
d = as.matrix(infercnv_obj_full_run@tumor_subclusters$subclusters$Unknown$Unknown_s4)
d[,1] = "Non-malignant"
e = as.matrix(infercnv_obj_full_run@tumor_subclusters$subclusters$Unknown$Unknown_s5)
e[,1] = "Non-malignant"
f= rbind(a,b,c,d,e)

sobj_int_layJoint = AddMetaData(sobj_int_layJoint, metadata = f, col.name = "infercnv_malignant")
DimPlot(sobj_int_layJoint, group.by="infercnv_malignant")
sobj_int_layJoint@meta.data$infercnv_malignant <- NULL

sobj_int_layJoint = infercnv::add_to_seurat(infercnv_output_path=output_dir_full,
                                            seurat_obj=sobj_int_layJoint, # optional
                                            top_n=10
)


# Tumor cell calling
## Extract data
observations = read.table(paste0(output_dir_full, "/infercnv.observations.txt"), header = TRUE, sep = " ")
colnames(observations) <- gsub("\\.1$", "-1", colnames(observations))
references = read.table(paste0(output_dir_full, "/infercnv.references.txt"), header = TRUE, sep = " ")
colnames(references) <- gsub("\\.1$", "-1", colnames(references))
all_cnv_values = cbind(observations, references)


head(references)
cnv_scores <- colSums(all_cnv_values^2)
hist(cnv_scores, breaks=200)
?hist
## Labels for obs vs ref
obs_ref = c(rep("obs", dim(observations)[2]), rep("ref", dim(references)[2]))
names(obs_ref) = c(colnames(observations), colnames(references))
## Sort the order of cells in cnv_values and obs_ref to the same as cm
cm = as.matrix(GetAssayData(sobj_int_layJoint, layer = "counts"))
all_cnv_values = all_cnv_values[colnames(cm)]
obs_ref = obs_ref[colnames(cm)]
all_cnv_values = all_cnv_values - 1
## Save
saveRDS(all_cnv_values, paste0(output_dir_full, "all_cnv_values.rds"))
saveRDS(obs_ref, paste0(output_dir_full, "obs_ref.rds"))
## 
samples = sobj_int_layJoint$orig.ident
malig_status = annotation_df
## Clustering
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend)
library(weights)

cnv_hc_each_sample = list()
cnv_cutree_each_sample = list()
weights_list = list()

for (sample in sort(unique(samples))){
  ## Subset each sample + (mg and od)
  to_subset = (samples == sample | samples == "mg" | samples == 'od')
  message("======================================")
  message("Processing sample: ", sample)
  message("Start time: ", Sys.time())
  ## Compute weights
  weights = apply(all_cnv_values[,to_subset], 1, var)
  weights_list[[sample]] = weights
  ## Compute pairwise correlations and distances
  pairwise_cor = wtd.cors(all_cnv_values[,to_subset], weight = weights) 
  pairwise_dist = 1 - pairwise_cor
  ## Hierarchical clustering and cutree 
  hc = hclust(as.dist(pairwise_dist), method="ward.D2")
  cnv_hc_each_sample[[sample]] = hc
  cnv_cutree_each_sample[[sample]] = cutree(hc, 4)
  message("End time: ", Sys.time())
}
`



















## Call functions
samples = sobj_int_layJoint$orig.ident
malig_status = annotation_df
malginant_cells = callTumorCell(all_cnv_values, samples, malig_status)
sobj_int_layJoint = AddMetaData(sobj_int_layJoint, metadata = malginant_cells, col.name = "infercnv_malignant")
DimPlot(sobj_int_layJoint, group.by="infercnv_malignant")

head(malginant_cells)
## Helper functions
callTumorCell <- function(cnv_values, samples, malig_status, cnv_score_cutoff=0.03, cnv_coef_cutoff=0.2){
  ## Calculate CNV score
  message("Calculating CNV scores...")
  cnv_score = colMeans(cnv_values^2)
  ## Calculate CNV correlation 
  message("Calculating CNV correlations...")
  mean_cnv_value = calculateAverageCnvValue(cnv_values, samples, malig_status)
  cor_coef = calculateCnvCor(cnv_values, mean_cnv_value, samples)
  ## Call malignant cells
  message("Call malignant cells...")
  result = ifelse((cnv_score < cnv_score_cutoff & cor_coef < cnv_coef_cutoff), "non-malignant", "malignant")
  names(result) = names(samples) #change??
  return(result)
}
calculateAverageCnvValue <- function(cnv_values, samples, malig_status, min_tumor_cell=3){
  try (if(any(colnames(cnv_values) != names(samples) | any(colnames(cnv_values) != rownames(malig_status))))
    stop("CNV values, sample names, and malignancy status should have identical cell names and orders"))
  
  result = list()
  unique_sample_names = names(table(samples))
  message("Subsetting cnv_values and samples with only expressionally defined tumor cells...")
  putative_malig = sapply(malig_status, function(x) grepl("^Unknown", x))
  cnv_values_tumor = cnv_values[,putative_malig]
  samples = samples[putative_malig]
  
  try(if(dim(cnv_values_tumor)[2] != length(samples)) stop("Subsetted cnv values and sample names should have equal length"))
  
  message("Computing average CNV values for each sample...")  
  for (sample_name in unique_sample_names){
    tmp_cnv_values = cnv_values_tumor[,which(samples == sample_name)]
    message("The number of putative tumor cells in sample ", sample_name, " is ", dim(tmp_cnv_values)[2])
    if (dim(tmp_cnv_values)[2] >= min_tumor_cell){
      result[[sample_name]] = rowMeans(tmp_cnv_values)
    }
    else{
      result[[sample_name]] = NA
    }
  }
  return(result)
}
calculateCnvCor <- function(cnv_values, mean_cnv_values, samples){
  try (if(any(colnames(cnv_values) != names(samples)))
    stop("CNV values and sample names should have identical cell names and orders"))
  
  cor_coef = NULL
  for (i in seq(dim(cnv_values)[2])){
    sample = as.character(samples[i])
    #print(sample)
    #print(head(mean_cnv_values$sample))
    if (sum(is.na(mean_cnv_values[[sample]]))>0){
      cor_coef = c(cor_coef, 0)
    }else{ 
      ##cor_coef = c(cor_coef, cor(cnv_values[,i], mean_cnv_values[[sample]], method = "spearman"))
      cor_obj = cor.test(cnv_values[,i], mean_cnv_values[[sample]], 
                         alternative = "two.sided", method = "spearman", exact = FALSE)
      if (cor_obj$p.value > 0.05){
        cor_coef = c(cor_coef, 0)
      }
      else{
        cor_coef = c(cor_coef, cor_obj$estimate)
      }
    }
  }
  
  names(cor_coef) = names(samples)
  return(cor_coef)
}






# Load CNV predictions per cell
cnv_cells <- read.table("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/saveDir/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat", 
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Load cell cluster assignments
cell_clusters <- read.table("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/saveDir/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings", 
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(cell_clusters) <- c("cluster", "cell_id")  # Assign column names correctly

# Count CNVs per cell using the correct column
cnv_counts_per_cell <- cnv_cells %>%
  group_by(cell_group_name) %>%  # Use correct column for cell IDs
  summarise(cnv_count = n())

# Merge CNV counts with cluster information
cnv_counts_with_clusters <- merge(cnv_counts_per_cell, cell_clusters, by.x = "cell_group_name", by.y = "cell_id")

# Count average CNVs per cluster
cnv_counts_per_cluster <- cnv_counts_with_clusters %>%
  group_by(cluster) %>%
  summarise(avg_cnv_count = mean(cnv_count))

# Define automatic CNV threshold for tumor classification
cnv_threshold <- median(cnv_counts_per_cluster$avg_cnv_count) + sd(cnv_counts_per_cluster$avg_cnv_count)
tumor_clusters <- cnv_counts_per_cluster$cluster[cnv_counts_per_cluster$avg_cnv_count > cnv_threshold]

# Label cells as Tumor/Non-Tumor based on cluster assignment
cell_clusters$tumor_status <- ifelse(cell_clusters$cluster %in% tumor_clusters, "Tumor", "Non-Tumor")

sobj_int$tumor_status <- cell_clusters$tumor_status[match(rownames(sobj_int@meta.data), cell_clusters$cell_id)]
DimPlot(sobj_int, group.by = "tumor_status")

view(sobj_int@meta.data)

DimPlot(sobj_int, group.by = "orig.ident")
DimPlot(sobj_int, group.by = "harmonyCluster", label=T)

plot_cnv(inferCNV_HMM, title = "CNV Heatmap")

head(cnv_genes)
head(cnv_regions)
head(tumor_clusters)
print(inferCNV_HMM)
tumor_clusters <- unique(cnv_regions$cell_group_name)



head(cell_groups)
scores=apply(inferCNVresults@expr.data,2,function(x){ sum(x < 0.95 | x > 1.05)/length(x) })
df <- data.frame(values = scores)
hist(scores, breaks=100)
# Plot histogram
ggplot(df, aes(x = values)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Histogram of Numeric Vector", x = "Values", y = "Frequency")
hist(scores)
view(scores)

inferCNVresults = readRDS("/data/cephfs-1/work/projects/aging-gbm/EmbeddingProject/saveDir/run.final.infercnv_obj")
print(inferCNVresults)

normals_path <- "/data/cephfs-1/work/projects/aging-gbm/analysis/inferCNV/cm_raw_normal_cells_DIPG.rds"
normal_cells_data <- as.matrix(readRDS(normals_path))
normal_cells_annotation <- data.frame(cell_type = rep("Normal", ncol(normal_cells_data)))
rownames(normal_cells_annotation) <- colnames(normal_cells_data)

unknown_cells_data <- as.matrix(GetAssayData(seurat_obj, slot = "counts"))
unknown_cells_annotations <- data.frame(cell_type = rep("Unknown", ncol(unknown_cells_data)))
rownames(unknown_cells_annotations) <- colnames(unknown_cells_data)

combined_cells_annotation <- rbind(normal_cells_annotation, unknown_cells_annotations)
# Get the gene names
genes_normal <- rownames(normal_cells_data)
genes_unknown <- rownames(unknown_cells_data)
# Find the intersection of both gene sets
common_genes <- intersect(genes_normal, genes_unknown)
# Subset both matrices to include only the common genes
normal_cells_data_aligned <- normal_cells_data[common_genes, , drop = FALSE]
unknown_cells_data_aligned <- unknown_cells_data[common_genes, , drop = FALSE]
combined_counts_matrix <- cbind(normal_cells_data_aligned, unknown_cells_data_aligned)

file_path_gene_pos <- "/data/cephfs-1/work/projects/aging-gbm/analysis/inferCNV/gencode_v19_gene_pos.txt"
infercnv_obj_full = CreateInfercnvObject(raw_counts_matrix = combined_counts_matrix,
                                         annotations_file = combined_cells_annotation,
                                         gene_order_file= file_path_gene_pos,
                                         delim = "\t",
                                         ref_group_names = c("Normal"))


output_dir_full = "/data/cephfs-1/work/projects/aging-gbm/analysis/inferCNV/FigureMDPhD"
infercnv_obj_full_run = infercnv::run(infercnv_obj_full,
                                      out_dir = output_dir_full,
                                      cutoff=1, #minimum reads that a gene must have on average across all cells not to be filtered out # SMARTseq: 1, 10x: 0.1
                                      #cluster_by_groups = T,
                                      HMM = T, # Hidden Markov Model: segmentation and prediction analysis that will generate clear cut regions with CNVs and assign them specific states of change
                                      denoise = T, #filtering of the residual expression which is the signal at the end of the infercnv analysis outside HMM
                                      leiden_resolution = 0.00005, # resolution parameter for the leiden clustering, lower = less clusters, std= 0.05
                                      up_to_step = 15, #stop analysis after a certain step to allow reevaltuation
                                      #debug = TRUE)
                                      num_threads = 20)

plot_cnv(infercnv_obj_full_run,
         out_dir = output_dir_full,
         title = "infercnv figure",
         obs_title = "Malginant cells",
         ref_title = "Normal cells",
         cluster_by_groups = T,
         plot_chr_scale = T,
         color_safe_pal = T,
         dynamic_resize = 0.5, # how high is the bottow head map?
         output_filename = "infercnv_scaled_to_chr"
)

# cluster normal cells and tumor cells separately



# ----------------- Move Object to Python ----------------
SaveH5Seurat(sobj_int_layJoint, filename = "seurat_data.h5Seurat", overwrite = TRUE)

