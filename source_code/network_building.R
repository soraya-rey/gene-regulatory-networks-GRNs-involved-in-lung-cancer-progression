# Network building file: 
# From pre-processed scRNAseq data, creates 3 networks using three methods: 
# Bayesian, GENIE3 and Aracne

# 0. Library importations
library(igraph)
library(GENIE3)
library(vegan)
library(Rgraphviz)
library(doParallel)
library(doRNG)
library(Matrix)
library(qgraph)
library(minet)

# 1. Data path
matrix_file <- './data/Epithelial_cells/lung_stage_I_processed.csv'
rds_file <- './data/Epithelial_cells/lung_stage_I_seurat.rds'
name <- './data/Epithelial_cells/lung_stage_I'

# 2. Data importations

# 2.1 Getting Seurat Object
#seurat_object <- readRDS(rds_file)

# 2.2 Getting processed matrix
matrix_data <- read.csv(matrix_file, header = TRUE, row.names = 1)
classes <- matrix_data["classes",]

matrix <- as.matrix(matrix_data[!(row.names(matrix_data) %in% c("classes")),])
expression_matrix <-  matrix(as.numeric(matrix),
                             (nrow(matrix_data)-1), 
                             ncol(matrix_data))
colnames(expression_matrix) <- colnames(matrix)
rownames(expression_matrix) <- rownames(matrix)

binary_classes <- matrix(0,1,length(matrix_data))
binary_classes[which(classes == "Normal")] <- -1.
binary_classes[which(classes == "Tumoral")] <- 1.
colnames(binary_classes) <- colnames(matrix_data)
rownames(binary_classes) <- c("class")

expression_matrix <- as.matrix(expression_matrix[!(row.names(expression_matrix) %in% c("class")),])
expression_matrix <- rbind(expression_matrix, binary_classes)
# 3. Creating Network with GENIE3

# 3.1 Adjacency matrix inference
adj_matrix <- GENIE3(expression_matrix)

# 3.3 Save network
write.csv(adj_matrix, file = paste(name,"genie3_matrix.csv",sep='_'))
pdf(paste(name,"genie3.pdf",sep='_'))
qgraph(adj_matrix, layout="spring", groups=which(colnames(adj_matrix) == 'class'), palette=c('pastel'))
dev.off()


# 4. Creating network with ARACNE

# 4.1  Getting distance matrix with Jaccard distance (from article: as results were not good,
# we used built in build.min method from minet with spearman estimator)
# 4.1.1 Calling Jaccard
# dist <- vegdist(expression_matrix, method = "jaccard")
# 4.1.2 Converting to distance
# distance_matrix <- matrix(0,nrow(expression_matrix),nrow(expression_matrix))
# distance_matrix[lower.tri(distance_matrix)] <- dist
# distance_matrix[upper.tri(distance_matrix)] <- t(distance_matrix)[upper.tri(distance_matrix)]
# colnames(distance_matrix) <- rownames(expression_matrix)
# rownames(distance_matrix) <- rownames(expression_matrix)

# 4.2 Computing network adjacency matrix
mim <- build.mim(t(expression_matrix),estimator="spearman")

# 4.3 Running ARACNE
net <- minet::aracne(mim, eps=0.05)

# 4.3 Save network
write.csv(net, file = paste(name,"aracne_matrix.csv",sep='_'))
pdf(paste(name,"aracne.pdf",sep='_'))
qgraph(net, layout="spring", groups=which(colnames(net) == 'class'), palette=c('pastel'))
dev.off()

