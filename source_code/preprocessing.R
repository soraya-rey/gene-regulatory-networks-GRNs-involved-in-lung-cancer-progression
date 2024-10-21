# Main project file

# 0. Library importations
library(Seurat)
library(GENIE3)
library(igraph)
library(ggplot2)
library(stringr) 

# Sources:
#Article: Single-cell RNA sequencing demonstrates the molecular and cellular reprogramming of metastatic lung adenocarcinoma
#Tutorial for data preprocessing and visualisation: https://holab-hku.github.io/Fundamental-scRNA/downstream.html



# 1. Loading Data
rna_seq_tumour_file <- "./data/Epithelial_cells/Tumour_lung_III_matrix.csv"
rna_seq_normal_file <- "./data/Epithelial_cells/Normal_lung_III_matrix.csv"
name <- "./data/Epithelial_cells/lung_stage_III"
separator=','
rna_seq_tumour <- read.csv(rna_seq_tumour_file, sep=separator,header = TRUE, row.names = 1)
rna_seq_normal <- read.csv(rna_seq_normal_file, sep=separator,header = TRUE, row.names = 1)

# 2. Creating Seurat Objects (contains both data for a single-cell dataset and its analysis)

create_seurat <- function(rna_seq_data, project){
# pre-processing function: creates a normalized Seurat object from RNAseq file
# execute quality control
  
  # 2.1 Minimum cells required per gene calculation (from article 0.1%)
  min_cells <- floor((length(rna_seq_data) - 1) * 0.001)
  
  # 2.2 Minimum read per gene and cell (from article 200)
  min_gene <- 200
  
  # 2.3 Seurat object
  seurat_object <- CreateSeuratObject(rna_seq_data, project=project, min.cells = min_cells, min.features = min_gene)
  
  
  
  # 3. Quality Control
  
  # 3.1 Add number of genes per UMI for each cell to metadata
  # (quantity of genes detected per Unique Molecular Identifier = molecular barcode)
  seurat_object$log10GenesPerUMI <- log10(seurat_object$nFeature_RNA) / log10(seurat_object$nCount_RNA)
  
  # 3.2 Computing mitochondrial percentage and filtering (from article max.mito = 0.2)
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  seurat_object[["percent.mt"]] <- seurat_object[["percent.mt"]]/sum(seurat_object[["percent.mt"]])
  View(seurat_object@meta.data)
  max_mito <- 0.2
  seurat_mito <- subset(x = seurat_object, cells = WhichCells(seurat_object, expression = percent.mt < max_mito))
  View(seurat_mito@meta.data)
  
  # 3.3 Visualisation
  VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0001)
  plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
    theme(legend.position="none")
  plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
    theme(legend.position="none")
  plot1 + plot2
  
  

  # 4. Transforming data
  
  # 4.1 Computing normalized datas
  seurat_object <- NormalizeData(seurat_object)
  
  # 4.2 Visualisation
  set.seed(123)
  par(mfrow=c(1,2))
  # original expression distribution
  raw_geneExp = as.vector(seurat_object[['RNA']]@counts) %>% sample(10000)
  raw_geneExp = raw_geneExp[raw_geneExp != 0]
  hist(raw_geneExp)
  # expression distribution after normalization
  logNorm_geneExp = as.vector(seurat_object[['RNA']]@data) %>% sample(10000)
  logNorm_geneExp = logNorm_geneExp[logNorm_geneExp != 0]
  hist(logNorm_geneExp)

  
  return(seurat_object)
} 

# calling pre-processing for both files
seurat_tumour <- create_seurat(rna_seq_tumour, 'tumour')
seurat_normal <- create_seurat(rna_seq_normal, 'normal')

# seurat object fusion
seurat_object <- merge(seurat_tumour, y = seurat_normal, add.cell.ids = c("tumour", "normal"), project = "lung_stage_I",
                                          merge.data = TRUE)



# 4.3 Computing hightly variable features (best for further analysis)
# 4.3.1 Literature's filtration parameters
min_exp <- 0.0125
max_exp <- 3
min_quantile_normalized_variance <- 0.5
# 4.3.2 Computing variability with parameters
seurat_var <- FindVariableFeatures(seurat_object, mean.cutoff = c(min_exp, max_exp), dispertion.cutoff = c(min_quantile_normalized_variance, Inf))

# 4.3.3 Showing high variability values
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_var), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_var) + 
  theme(legend.position="top")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + 
  theme(legend.position="none")
plot1 + plot2

# 4.4 Scalling Data (standard pre-processing: shifts mean expression to 0, variance to 1)
all.genes <- rownames(seurat_var)
seurat_scaled <- ScaleData(seurat_var, features = all.genes, verbose = FALSE)



# 5. Unsupervised dimensional reduction and clustering (following article's process)

# 5.1 Applying PCA
seurat_object <- RunPCA(seurat_scaled)

# 5.2 Vizualisation
print(x = seurat_object[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)
VizDimLoadings(seurat_object, dims = 1:2, reduction = "pca")
DimPlot(seurat_object, reduction = "pca") + theme(legend.position="none")
DimHeatmap(seurat_object, dims = 1:9, cells = 500, balanced = TRUE)

# 5.2 Selecting PCs with high statistical significance
# 5.2.1 Computing standard deviation of PCs with JackStraw function (from article)
seurat_jack <- JackStraw(seurat_object)
# 5.2.2 Plotting values
plot <- ElbowPlot(seurat_jack)
print(plot)
pc_values <- plot$data$stdev
# 5.2.3 Selecting PC subset (cutting subset where more PC don't explain much more variance)
treshold <- 0.1
last_pc <- length(pc_values)
last_pc_value <- pc_values[1]
for (i in c(2:length(pc_values))){
  current_pc_value <- pc_values[i]
  if (abs(last_pc_value - current_pc_value) < treshold) {
    last_pc <- i - 1
    break
  }
  last_pc_value <- current_pc_value
}

# 5.3 Clustering data (using TSNE according to article)
seurat_neighbors <- FindNeighbors(object = seurat_jack, 
                            dims = 1:last_pc)
seurat_clustered <- FindClusters(seurat_neighbors, dims = 1:last_pc)
seurat_clustered <- RunTSNE(seurat_clustered, dims.use = 1:last_pc, check_duplicates = FALSE)
DimPlot(seurat_clustered, reduction = "tsne")



# 6. Selecting data for marker genes

# 6.1 Finding marker genes matrix (with params from article) and selecting 
markers_matrix <- FindAllMarkers(seurat_clustered, logfc.threshold = 1, min.pct=0.25, test.use = "t")

# 6.2 selecting markers with Student pvalue < 0.01 and Bonferri adjustement pvalue < 0.01
# (from article's parameters)
# Moreover, we are selecting genes over the log2 Fold Change, in order to have an amount of 
# genes close to the expected (67). We choose to keep the 100 most significative genes
# 50 lowest below 0, 50 hightest above 0

markers_matrix <- subset(markers_matrix, subset = p_val < 0.01 & p_val_adj < 0.01)
markers_sorted <- markers_matrix[order(markers_matrix[["avg_log2FC"]]), ]
last_index <- length(markers_sorted[["avg_log2FC"]])
first_index <- last_index - 49
markers_mat <- markers_sorted[c(c(1:50),c(first_index: last_index)),]
View(markers_mat)

# 6.3 Selecting corresponding data in seurat object matrix
regulators <- unique(rownames(markers_mat))
expression_matrix <- as.matrix(seurat_object[["RNA"]]$data)
processed_data <-subset(expression_matrix, subset = (rownames(expression_matrix) %in% regulators))

columns <- colnames(processed_data)
classes <- c(1:length(columns))
classes[1:length(classes)] <- rep('Normal')
for (i in c(1:length(classes))) {
  if (str_sub(columns[i],1,3) == 'tum') {
    classes[i] = 'Tumoral'
  }
}
processed_data <- rbind(processed_data, classes)

# 7. Writing Seurat Object and Processed expression matrix for later use
seurat_file <- paste(name, "_seurat.rds", sep='')
saveRDS(seurat_clustered, file = seurat_file)
matrix_file <-  paste(name, "_processed.csv", sep='')
write.csv(processed_data, file = matrix_file)

