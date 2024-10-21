# librairies
library(igraph)

# imports
file1 <- "./data/Epithelial_cells/lung_stage_I_aracne_matrix.csv"
file2 <- "./data/Epithelial_cells/lung_stage_II_aracne_matrix.csv"
file3 <- "./data/Epithelial_cells/lung_stage_III_aracne_matrix.csv"
name <- ""

graph1 <- read.csv(file1, header = TRUE, row.names = 1)
graph2 <- read.csv(file2, header = TRUE, row.names = 1)
graph3 <- read.csv(file3, header = TRUE, row.names = 1)

# getting genes hightly correlated to class in networks
subgraph1 <- graph1[order(graph1[["class"]], decreasing = TRUE)]
colnames(subgraph1[1:5])
subgraph2 <- graph2[order(graph2[["class"]], decreasing = TRUE)]
colnames(subgraph2[1:5])
subgraph3 <- graph3[order(graph3[["class"]], decreasing = TRUE)]
colnames(subgraph3[1:5])

graph1 <- graph_from_adjacency_matrix(as.matrix(graph1), mode='undirected')
graph2 <- graph_from_adjacency_matrix(as.matrix(graph2), mode='undirected')
graph3 <- graph_from_adjacency_matrix(as.matrix(graph3), mode='undirected')

plt_density <- function(d1, d2, d3, title) {
  plot(d1, col='blue', 
       main=title,
       xlim = c(min(d1$x,d2$x,d3$x), max(d1$x, d2$x, d3$x)),
       ylim = c(min(d1$y,d2$y,d3$y), max(d1$y, d2$y, d3$y)))
  lines(d2, col='green')
  lines(d3, col="red")
  
  legend("topright", lwd=1, legend = c("Stage 1", "Stage 2", "Stage 3"), col = c("blue", "green", "red"))
}

# computing betweeness
betweeness1 <- betweenness(graph1, cutoff = 0)
density1 <- density(betweeness1[betweeness1!=0])
betweeness2 <- betweenness(graph2, cutoff = 0)
density2 <- density(betweeness2[betweeness2!=0])
betweeness3 <- betweenness(graph3, cutoff = 0)
density3 <- density(betweeness3[betweeness3!=0])

plt_density(density1, density2, density3, 'Betweeness centrality distribution (MIIC)')

pdf(paste(name,"betweeness_cent_distri_miic.pdf",sep='_'))
plt_density(density1, density2, density3, 'Betweeness centrality distribution (MIIC)')
dev.off()

# computing degrees
degrees1 <- igraph::degree(graph1)
density1 <- density(degrees1[degrees1!=0])
degrees2 <- igraph::degree(graph2)
density2 <- density(degrees2[degrees2!=0])
degrees3 <- igraph::degree(graph3)
density3 <- density(degrees3[degrees3!=0])

plt_density(density1, density2, density3, 'Degree centrality distribution (MIIC)')

pdf(paste(name,"degree_cent_distri_miic.pdf",sep='_'))
plt_density(density1, density2, density3, 'Degree centrality distribution (MIIC)')
dev.off()

# computing degrees
degrees1 <- igraph::degree(graph1)
density1 <- density(degrees1)
degrees2 <- igraph::degree(graph2)
density2 <- density(degrees2)
degrees3 <- igraph::degree(graph3)
density3 <- density(degrees3)

plt_density(density1, density2, density3, 'Degree centrality distribution (MIIC)')

pdf(paste(name,"degree_cent_distri_miic.pdf",sep='_'))
plt_density(density1, density2, density3, 'Degree centrality distribution (MIIC)')
dev.off()

# computing closeness
clos1 <- igraph::closeness(graph1)
density1 <- density(clos1)
clos2 <- igraph::closeness(graph2)
density2 <- density(clos2)
clos3 <- igraph::closeness(graph3)
density3 <- density(clos3)

plt_density(density1, density2, density3, 'Closeness centrality distribution (MIIC)')

pdf(paste(name,"closeness_cent_distri_miic.pdf",sep='_'))
plt_density(density1, density2, density3, 'Closeness centrality distribution (MIIC)')
dev.off()

