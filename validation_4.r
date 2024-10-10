source("utilities.r")
library(ranger, lib.loc = "/home/christelsirocchi/R/x86_64-pc-linux-gnu-library/4.3") 

########################################################
############# GENERATE DATA & TRAIN MODELS #############
########################################################

cluster_number = 4

# Combine cluster centers into a matrix
cluster_centers <- get_cluster_centers(cluster_number)
k <- dim(cluster_centers)[1]
n_feat <- dim(cluster_centers)[2]
n_samples_per_cluster <- 50  
n_iter <- 30

# save training data and models
model_list <- list()
train_list <- list()

for (n in 1:n_iter) {
  train <- get_synthetic_dataset(cluster_centers, n_samples_per_cluster)
  train_list[[n]] <- train
  target <- rep(1:k, each = n_samples_per_cluster)
  #ggplot(train, aes(x = V2, y = V3, color = as.factor(target))) + geom_point() 
  
  model = ranger(x=train, num.trees=500, 
                 clustering = TRUE,  probability= FALSE, oob.error = FALSE,
                 mtry = round(sqrt(ncol(train))), min.bucket = 5)
  
  DISTANCE = distance_rf(model, train)
  hc = hclust(as.dist(DISTANCE), method="ward.D2")
  res1 = cutree(hc, k)
  print(ARI(res1, target))
  print(NMI(res1, target))
  model_list[[n]] <- model
}

# in this experiment only consider "sample" criterion
edge_types = c("sample") #c("present", "fixation", "level", "sample")
all_edge_matrix <- array(0, dim = c(n_feat+1, n_feat+1, n_iter, length(edge_types)))

for (e in 1:length(edge_types)){ 
  for (n in 1:n_iter) {
    print(n)
    model <- model_list[[n]]
    train <- train_list[[n]]
    all_edge_matrix[,,n,e] <- compute_edge_matrix(model, train, edge_types[e])
  }
}

# save results
var_names <- model$forest$independent.variable.names #paste("V", seq(n_feat), sep="")
dimnames(all_edge_matrix) <- list(c(var_names,"T"), c(var_names,"T"), paste("iter", seq(n_iter), sep=""), edge_types)
edge_df <- melt(all_edge_matrix)
colnames(edge_df) <- c("node_from", "node_to", "iteration", "edge_type", "value")
file_name_all_e <- file.path("synthetic", paste("v4_all_edges_cluster_number", cluster_number, "repeated_.csv", sep = "_"))
write.csv(edge_df, file = file_name_all_e, row.names = FALSE)

stop("finished computation")

##############################################################
############# READ DATA & ANALYSE FEATURE GRAPHS #############
##############################################################

###### RDUNDANT FEATURES EXPERIMENT - CENTRALITY #############

cluster_number = 4
file_name_all_e <- file.path("synthetic", paste("v4_all_edges_cluster_number", cluster_number, "repeated_.csv", sep = "_"))
edge_df <- read.csv(file = file_name_all_e)
# get parameters  
all_edge_matrix <- acast(edge_df, node_from ~ node_to ~ iteration ~ edge_type)
n_feat <- dim(all_edge_matrix)[1]-1
max_n_feat <- n_feat - 1 
n_iter <- dim(all_edge_matrix)[3]
var_names <- paste("V", seq(n_feat), sep="")
all_edge_matrix <- all_edge_matrix[c(var_names,"T"), c(var_names,"T"),,] #reorder features

# compute out-degree centrality for each edge matrix
all_out_degree_centralities <- array(0, dim = c(length(c(var_names,"T")), n_iter))  
for (n in 1:n_iter) {
  edge_matrix <- all_edge_matrix[,,n]
  all_out_degree_centralities[,n] <- rowSums(edge_matrix)/ sum(rowSums(edge_matrix)) *100
}  
dimnames(all_out_degree_centralities) <- list(c(var_names,"T"), paste("iter",seq(n_iter),sep=""))
all_cent_df <- melt(all_out_degree_centralities)
colnames(all_cent_df) <- c("node", "iteration", "value")

# plot centralities (to show that all relevant features have similar centrality)
ggplot(all_cent_df[all_cent_df$node != "T",], aes(x = node, y = value)) + geom_boxplot() + theme_bw() +
  labs(#title = "node centrality of relevant and irrelevant features",
       y = "out-degree centrality", x = "features") + 
  coord_cartesian(ylim = c(0, 18)) 

# compute the weight of every clique of size 3 
clique_size <- 3 
cliques <- t(combn(var_names, clique_size))
all_clique_weights <- array(0, dim = c(nrow(cliques), n_iter))
for (n in 1:n_iter){
  # normalise edge matrix
  edge_matrix <- all_edge_matrix[var_names,var_names,n] 
  edge_matrix <- edge_matrix/sum(edge_matrix)*100
  # compute weight of each triad
  g <- graph.adjacency(as.matrix(edge_matrix), weighted = TRUE, mode = "directed")
  all_clique_weights[,n] <- apply(cliques, 1, function(clique) sum(E(induced.subgraph(graph=simplify(g),vids=clique))$weight/6))
}
# save weights to dataframe
rownames(all_clique_weights) <- apply( cliques, 1 , paste , collapse = "-" )
colnames(all_clique_weights) <- paste("iter", seq(n_iter), sep="")
all_clique_weights <- all_clique_weights[order(rowMeans(all_clique_weights), decreasing = TRUE),]
all_clique_df <- melt(all_clique_weights)
colnames(all_clique_df) <- c("clique", "iteration", "value")

# label optimal combinations of relevant features (that do not include repeated info)
sel_triples <- expand.grid(c("V1", "V2"), c("V3", "V4"), c("V5", "V6"))
sel_triples_names <- apply( sel_triples, 1, paste, collapse = "-" )
all_clique_df$optimal <- all_clique_df$clique %in% sel_triples_names

# plot triads (to show that optimal combinations of three features correspond to the heaviest triads)
ggplot(all_clique_df, aes(x = clique, y = value, fill = optimal)) +
  geom_boxplot(outlier.size = 0.7) +
  labs(x = "feature triads", y = "average edge weight (AW)") +
  theme_minimal() + guides(fill="none") +
  scale_fill_manual(values=c("white", "#EE4C4C")) + theme_bw() +
  theme( axis.text.y = element_text( size = 9 ),
         axis.text.x = element_text( size = 8, angle = 90),
         axis.title = element_text( size = 10 ),
         strip.text = element_text( size = 11 ))#, face = "bold"

# confirm that in all iterations the heaviest triad correspond to one of the optimal combinations
heavy_triad_df <- all_clique_df %>% group_by(iteration) %>% filter(value == max(value, na.rm=TRUE)) %>% as.data.frame()

