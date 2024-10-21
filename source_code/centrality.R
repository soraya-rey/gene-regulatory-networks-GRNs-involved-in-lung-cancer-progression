if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("bnlearn")
library(igraph)
library(ggplot2)
install.packages("dplyr")
library(dplyr)

stade.1.adj.matrix <- as.matrix(stade.1.adj.matrix)
graph1 <- graph.adjacency(stade.1.adj.matrix)

stade.2.adj.matrix <- as.matrix(stade.2.matrix.adj)
graph2 <- graph.adjacency(stade.2.adj.matrix)

stade.3.adj.matrix <- as.matrix(stade.3.matrix.adj)
graph3 <- graph.adjacency(stade.3.adj.matrix)


# Calcul de la centralité des nœuds
degree_centrality1 <- degree(graph1)
closeness_centrality1 <- closeness(graph1)
betweenness_centrality1 <- betweenness(graph1)

degree_centrality1_ordered <-  degree_centrality1[order(degree_centrality1)]
best_stade_1_dg <- degree_centrality1_ordered[45:50]
print("best degree centrality stade 1 :")
print(best_stade_1_dg)

closeness_centrality1_ordered <-  closeness_centrality1[order(closeness_centrality1)]
best_stade_1_cc <- closeness_centrality1_ordered[45:50]
print("best closeness centrality stade 1 :")
print(best_stade_1_cc)

betweenness_centrality1_ordered <-  betweenness_centrality1[order(betweenness_centrality1)]
best_stade_1_bc <- betweenness_centrality1_ordered[45:50]
print("best betweenness centrality stade 1 :")
print(best_stade_1_bc)



# Calcul de la centralité des nœuds
degree_centrality2 <- degree(graph2)
closeness_centrality2 <- closeness(graph2)
betweenness_centrality2 <- betweenness(graph2)

degree_centrality2_ordered <- degree_centrality2[order(degree_centrality1)]
best_stade_2_dg <- degree_centrality2_ordered[45:50]
print("best degree centrality stade 2:")
print(best_stade_2_dg)

closeness_centrality2_ordered <- closeness_centrality2[order(closeness_centrality2)]
best_stade_2_cc <- closeness_centrality2_ordered[45:50]
print("best closeness centrality stade 2:")
print(best_stade_2_cc)

betweenness_centrality2_ordered <- betweenness_centrality2[order(betweenness_centrality2)]
best_stade_2_bc <- betweenness_centrality2_ordered[45:50]
print("best betweenness centrality stade 2:")
print(best_stade_2_bc)



# Calcul de la centralité des nœuds
degree_centrality3 <- degree(graph3)
closeness_centrality3 <- closeness(graph3)
betweenness_centrality3 <- betweenness(graph3)

degree_centrality3_ordered <- degree_centrality3[order(degree_centrality1)]
best_stade_3_dg <- degree_centrality3_ordered[45:50]
print("best degree centrality stade 3:")
print(best_stade_3_dg)

closeness_centrality3_ordered <- closeness_centrality3[order(closeness_centrality3)]
best_stade_3_cc <- closeness_centrality3_ordered[45:50]
print("best closeness centrality stade 3:")
print(best_stade_3_cc)

betweenness_centrality3_ordered <- betweenness_centrality3[order(betweenness_centrality3)]
best_stade_3_bc <- betweenness_centrality3_ordered[45:50]
print("best betweenness centrality stade 3:")
print(best_stade_3_bc)



degree_centrality3 <- degree_centrality3[1:50]
closeness_centrality3 <- closeness_centrality3[1:50]
betweenness_centrality3 <- betweenness_centrality3[1:50]

degree_centrality1 <- degree_centrality1[1:50]
closeness_centrality1 <- closeness_centrality1[1:50]
betweenness_centrality1 <- betweenness_centrality1[1:50]

degree_centrality2 <- degree_centrality2[1:50]
closeness_centrality2 <- closeness_centrality2[1:50]
betweenness_centrality2 <- betweenness_centrality2[1:50]



# Combine the data for all stages
best_betweenness_data <- data.frame(
  best_Gene = c(names(best_stade_1_bc), names(best_stade_2_bc), names(best_stade_3_bc)),
  Stage = c("Stage 1", "Stage 2", "Stage 3"),
  best_betweenness = c(best_stade_1_bc, best_stade_2_bc, best_stade_3_bc)
)
# Plotting
ggplot(best_betweenness_data, aes(x = best_Gene, y = best_betweenness, fill = Stage)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Stage 1" = "red", "Stage 2" = "green", "Stage 3" = "blue")) +
  labs(title = "betweenness Centrality Across Stages", x = "Gene", y = "Betweenness Centrality") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, margin = margin(t = 10)))


# Combine the data for all stages
best_degree_data <- data.frame(
  best_Gene = c(names(best_stade_1_dg), names(best_stade_2_dg), names(best_stade_3_dg)),
  Stage = c("Stage 1", "Stage 2", "Stage 3"),
  best_degree = c(best_stade_1_dg, best_stade_2_dg, best_stade_3_dg)
)
# Plotting
ggplot(best_degree_data, aes(x = best_Gene, y = best_degree, fill = Stage)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Stage 1" = "red", "Stage 2" = "green", "Stage 3" = "blue")) +
  labs(title = "Degree Centrality Across Stages", x = "Gene", y = "degree Centrality") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, margin = margin(t = 10)))


# Combine the data for all stages
best_closeness_data <- data.frame(
  best_Gene = c(names(best_stade_1_cc), names(best_stade_2_cc), names(best_stade_3_cc)),
  Stage = c("Stage 1", "Stage 2", "Stage 3"),
  best_closeness = c(best_stade_1_cc, best_stade_2_cc, best_stade_3_cc)
)
# Plotting
ggplot(best_closeness_data, aes(x = best_Gene, y = best_closeness, fill = Stage)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Stage 1" = "red", "Stage 2" = "green", "Stage 3" = "blue")) +
  labs(title = " Closeness Centrality Across Stages", x = "Gene", y = "Closeness Centrality") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, margin = margin(t = 10)))
