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



# ------------------ PostSoupX Data ----------------
path_output = "/data/cephfs-1/work/projects/aging-gbm/analysis/PostSoupX/"
paths_output_files = list()

for (i in 1:4) {
  for (j in 1:length(samples[[i]])) { 
    sample_name = samples[[i]][j]
    sample_path = paths[i]
    # Load data and initialize SoupChannel
    raw <- Read10X(data.dir = paste0(sample_path,sample_name,"/count/sample_raw_feature_bc_matrix"))
    filtered <- Read10X(data.dir = paste0(sample_path,sample_name,"/count/sample_filtered_feature_bc_matrix"))
    raw=raw[rownames(filtered),]
    sc = SoupChannel(raw, filtered)
    # Add metadata
    cluster.info <- read.csv(paste0(sample_path,sample_name,"/count/analysis/clustering/gene_expression_graphclust/clusters.csv"))
    umap.info <- read.csv(paste0(sample_path,sample_name,"/count/analysis/umap/gene_expression_2_components/projection.csv"))
    MetaData <- merge(umap.info, cluster.info, by = "Barcode")
    sc = setClusters(sc, setNames(MetaData$Cluster, rownames(MetaData)))
    sc = setDR(sc, umap.info[, 2:3])
    # Clean data
    sc = autoEstCont(sc, doPlot=FALSE)
    out = adjustCounts(sc, roundToInt = TRUE)
    # Create & seurat obj
    sample_name = substr(sample_name, nchar(sample_name) - 10, nchar(sample_name))
    sobj <- CreateSeuratObject(out, project=sample_name)
    filename_sample = paste0(path_output,sample_name,"_seurat_postSoupX.rds")
    saveRDS(sobj, filename_sample)
    paths_output_files = append(paths_output_files,list(filename_sample))
    print(filename_sample)
  }
}

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

# ------------------ Draw Pre-processing figures ---------------------
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
