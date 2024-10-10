# load packages
library(reshape2)
library(igraph)
library(pheatmap)
library(MASS)
library(ggplot2)
library(matrixcalc)
library(scales)
library(RColorBrewer)
library(gridExtra)
library(aricode)
library(dplyr)
library(fpc)
library(HDclassif)
library(kmed)
library(mlbench)
library(cluster)
library(tidyr)
library(fpc)
library(clValid)
library(mclust)
library(clusterCrit)

# function to compute diversity in the fixation index
get_diversity <- function(values){
  if (length(values) == 1){
    return(0)
  }
  else{
    dist_matrix <- outer(values, values, FUN = function(x, y) (x - y)^2)
    diversity <- mean(c(as.dist(dist_matrix)))
    return(diversity)
  }
}

# function to compute the fixation index given a vector and a split value 
fixation_index <- function(var_data, splitval){
  # get diversity for left values
  values_l <- var_data[var_data <= splitval]
  diversity1 <- get_diversity(values_l)
  # get diversity for right values
  values_r <- var_data[var_data > splitval]
  diversity2 <- get_diversity(values_r)
  # and average
  diversity = (diversity1 + diversity2)/2
  # compute divergence
  dist_matrix <- outer(values_l, values_r, FUN = function(x, y) (x - y)^2)
  divergence <- mean(dist_matrix)
  return(1 - diversity/divergence)
}

# recursive function to calculate the distance of a node from the root
calculate_level <- function(node_id, tree_data, current_level = 1) {
  if (node_id == 0) {
    return(current_level)
  } else {
    parent_id <- which(tree_data$leftChild == node_id | tree_data$rightChild == node_id) - 1
    calculate_level(parent_id, tree_data, current_level + 1)
  }
}

# recursive function to traverse the tree and update the edge matrix
fixation_traverse_tree <- function(tree_data, nodeID, train, edge_matrix, type, sample_size, cluster_index = NULL) {
  splitvarname <- tree_data$splitvarName[nodeID + 1]
  splitval <- tree_data$splitval[nodeID + 1]
  leftID <- tree_data$leftChild[nodeID + 1]
  rightID <- tree_data$rightChild[nodeID + 1]
  left_name <- tree_data$splitvarName[leftID + 1]
  right_name <- tree_data$splitvarName[rightID + 1]
  # evaluate split if node is not terminal
  if (tree_data[nodeID + 1,]$terminal=="FALSE") {
    # split dataset at the node
    left_df <- train[train[,splitvarname] <= splitval,]
    right_df <- train[train[,splitvarname] > splitval,]
    # define left and right weights based on edge criteria
    if (type == "present") {
      # unitary weight
      left_weight <- right_weight <- 1
    } else if (type == "fixation") {
      # weight equal to the fixation index of the split
      left_weight <- right_weight <- fixation_index(train[,splitvarname], splitval)
    } else if (type == "level") { 
      # weight scaled by distance from the root node
      left_weight <- right_weight <- 1/tree_data$level[nodeID + 1] #0.5^tree_data$level[nodeID + 1] 
    } else if (type == "sample") {
      # weight scaled by the number of samples traversing the edge
      left_weight <- NROW(left_df) / sample_size
      right_weight <- NROW(right_df) / sample_size
    } 
    # if the cluster criteria is applied, scale the weight by the fraction of nodes belonging to a given cluster
    if (!is.null(cluster_index)) { 
      left_factor <- ifelse(NROW(left_df) > 0, length(intersect(as.numeric(rownames(left_df)), cluster_index)) / NROW(left_df), 0)
      right_factor<- ifelse(NROW(right_df) > 0, length(intersect(as.numeric(rownames(right_df)), cluster_index)) / NROW(right_df), 0)
      left_weight <- left_weight*left_factor
      right_weight<- right_weight*right_factor
    } 
    # update edge weight
    edge_matrix[splitvarname, left_name] <- edge_matrix[splitvarname, left_name] + left_weight   
    edge_matrix[splitvarname, right_name] <- edge_matrix[splitvarname, right_name] + right_weight
    # continue traversing the tree
    edge_matrix <- fixation_traverse_tree(tree_data, leftID, train[train[,splitvarname] <= splitval,], edge_matrix, type, sample_size, cluster_index)
    edge_matrix <- fixation_traverse_tree(tree_data, rightID, train[train[,splitvarname] > splitval,], edge_matrix, type, sample_size, cluster_index)
    }
  return(edge_matrix)
}

# function to compute the adjacency edge matrix 
compute_edge_matrix <- function(model, train, edge_type, cluster_index = NULL){
  # initialise edge matrix
  var_names <- model$forest$independent.variable.names
  edge_matrix <- matrix(0, nrow = length(var_names)+1, ncol = length(var_names)+1)
  rownames(edge_matrix) <- colnames(edge_matrix) <- c(var_names,"T")
  # update edge matrix after traversing each tree
  for (i in 1:model$num.trees) {
    # read tree data
    tree_data <- treeInfo(model, i)
    # label terminal nodes
    tree_data$splitvarName[tree_data$terminal=="TRUE"] <- "T"
    # add level to each split
    tree_data$level <- sapply(0:(nrow(tree_data) - 1), function(node_id) { calculate_level(node_id, tree_data) })
    # recursively traverse tree and update edge matrix
    edge_matrix <- fixation_traverse_tree(tree_data, 0, train, edge_matrix, edge_type, nrow(train), cluster_index) 
  }
  return(edge_matrix)
}

# recursive function to traverse the tree selecting edges traversed by a given sample
sample_traverse_tree <- function(tree_data, nodeID, train, sample, edge_matrix, type, sample_size, cluster_index = NULL) {
  splitvarname <- tree_data$splitvarName[nodeID + 1]
  splitval <- tree_data$splitval[nodeID + 1]
  direction <- ifelse(as.vector(sample[splitvarname]) <= splitval, "left", "right")
  childID <- ifelse(direction == "left", tree_data$leftChild[nodeID + 1], tree_data$rightChild[nodeID + 1])
  child_name <- tree_data$splitvarName[childID + 1]
  # evaluate split if node is not terminal
  if (tree_data[nodeID + 1,]$terminal=="FALSE") {
    # split dataset at the node
    if(direction == "left") {
      sub_df <-  train[train[,splitvarname] <= splitval,]
    } else{
      sub_df <- train[train[,splitvarname] > splitval,]
    }
    # define left and right weights based on edge criteria
    if (type == "present") {
      # unitary weight
      weight <- 1
    } else if (type == "fixation") {
      # weight equal to the fixation index of the split
      weight <- fixation_index(train[,splitvarname], splitval)
    } else if (type == "level") { 
      # weight scaled by distance from the root node
      weight <- 1/tree_data$level[nodeID + 1] #0.5^tree_data$level[nodeID + 1] 
    } else if (type == "sample") {
      # weight scaled by the number of samples traversing the edge
      weight <- NROW(sub_df) / sample_size
    } 
    # if the cluster criteria is applied, scale the weight by the fraction of nodes belonging to a given cluster
    if (!is.null(cluster_index)) { 
      factor <- ifelse(NROW(sub_df) > 0, length(intersect(as.numeric(rownames(sub_df)), cluster_index)) / NROW(sub_df), 0)
      weight <- weight*factor
    } 
    # update edge weight
    edge_matrix[splitvarname, child_name] <- edge_matrix[splitvarname, child_name] + weight   
    # continue traversing the tree
    edge_matrix <- sample_traverse_tree(tree_data, childID, sub_df, sample, edge_matrix, type, sample_size, cluster_index)
  }
  return(edge_matrix)
}


# function to compute the adjacency edge matrix for edges traversed by a given sample
sample_edge_matrix <- function(model, train, sample, edge_type, cluster_index = NULL){
  # initialise edge matrix
  var_names <- model$forest$independent.variable.names
  edge_matrix <- matrix(0, nrow = length(var_names)+1, ncol = length(var_names)+1)
  rownames(edge_matrix) <- colnames(edge_matrix) <- c(var_names,"T")
  # update edge matrix after traversing each tree
  for (i in 1:model$num.trees) {
    # read tree data
    tree_data <- treeInfo(model, i)
    # label terminal nodes
    tree_data$splitvarName[tree_data$terminal=="TRUE"] <- "T"
    # add level to each split
    tree_data$level <- sapply(0:(nrow(tree_data) - 1), function(node_id) { calculate_level(node_id, tree_data) })
    # recursively traverse tree and update edge matrix
    edge_matrix <- sample_traverse_tree(tree_data, 0, train, sample, edge_matrix, edge_type, nrow(train), cluster_index) 
  }
  return(edge_matrix)
}


# function to compute a matrix indicating the number of clusters a pair of features can discriminate
get_discr_matrix <- function(cluster_centers){
  n_vars = ncol(cluster_centers)
  discr_matrix <- matrix(0, nrow = n_vars, ncol = n_vars)
  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      discr_matrix[i,j] <- length(unique(apply(cluster_centers[,c(i,j)], 1, function(col) paste(col, collapse = ""))))
    }
  }
  return(discr_matrix)
}

# plot heatmap of adjacency matrix
plot_heatmap <- function(edge_matrix){
  pheatmap(edge_matrix,
           #legend_breaks = c(20, 30, 40, 50, 60, max(edge_matrix)), 
           main = "", #legend_labels = c("20", "30", "40", "50", "60", "edge\nweight"),
           fontsize_row = 10, fontsize_col = 10, #breaks = seq(0, max(all_edge_matrix), length.out = 30),
           color = colorRampPalette(c("white", "darkblue"))(30),
           cluster_cols = FALSE, cluster_rows = FALSE,
           cellheight=15, cellwidth = 15)  
}

# function to compute urf distance from trained model
distance_rf <- function(model, train){
  pred.ranger <- predict(model, train, type = "terminalNodes")
  lyves <- pred.ranger$prediction
  # Get Similarity
  SIMILARITY = matrix(0, dim(lyves)[1], dim(lyves)[1])
  for (xx in 1:(dim(SIMILARITY)[1]-1)){
    for (yy in xx:dim(SIMILARITY)[1]){
      res = lyves[c(xx,yy),]
      hit = sum(apply(res,2,function(x){x[1]==x[2]}))
      SIMILARITY[xx,yy] = hit
      SIMILARITY[yy,xx] = hit
    }
  }
  DISTANCE <- 1 - SIMILARITY/max(SIMILARITY)
  return(DISTANCE)
}

# function to compute weighted monotonicity of a sequence
compute_weighted_monotonicity <- function(sequence) {
  total_decrease <- sum(abs(pmin(0, diff(sequence))))
  total_variation <- sum(abs(diff(sequence)))
  if (total_variation == 0) return(1)  # Handle constant sequence case
  return(1 - total_decrease / total_variation)
}

# function to compute stepmonotonicity of a sequence
compute_step_monotonicity <- function(sequence) {
  n <- length(sequence)
  non_decreasing_steps <- sum(diff(sequence) >= 0)
  return(non_decreasing_steps / (n - 1))
}

# function to visualize the structure of the tree from treeinfo data
# example usage: visualize_tree(treeInfo(model, 1))
visualize_tree <- function(data) {
  # create an empty graph
  g <- graph.empty(directed = TRUE)
  # add vertices to the graph with nodeIDs and labels
  for (i in 1:nrow(data)) {
    vertex_color <- ifelse(data$terminal[i], "green", "orange") 
    vertex_label <- ifelse(is.na(data$splitvarName)[i], "", data$splitvarName[i]) 
    g <- add.vertices(g, n = 1, name = as.character(data$nodeID[i]), label = vertex_label, color = vertex_color)
  }
  # add edges to the graph
  for (i in 1:nrow(data)) {
    if (!is.na(data$leftChild[i])) {
      g <- add_edges(g, edges = c(as.character(data$nodeID[i]), as.character(data$leftChild[i])))
    }
    if (!is.na(data$rightChild[i])) {
      g <- add_edges(g, edges = c(as.character(data$nodeID[i]), as.character(data$rightChild[i])))
    }
  }
  # plot the graph
  plot(g, layout = layout_as_tree(g), vertex.label = V(g)$label, vertex.size = 12)
}

# function to visualize the structure of the graph from treeinfo data
# example usage: visualize_graph(treeInfo(model, 1))
visualize_graph <- function(data) {
  # create an empty graph
  g <- graph.empty(directed = TRUE)
  # add vertices 
  unique_vertices <- unique(na.omit(data$splitvarName))
  g <- add_vertices(g, nv = length(unique_vertices))
  V(g)$name <- unique_vertices
  # add edges to the graph
  for (i in 1:nrow(data)) {
    if (!is.na(data$leftChild[i]) && !is.na(data$splitvarName[data$leftChild[i] + 1])) {
      g <- add_edges(g, edges = c(data$splitvarName[i], data$splitvarName[data$leftChild[i] + 1]))
    }
    if (!is.na(data$rightChild[i]) && !is.na(data$splitvarName[data$rightChild[i] + 1])) {
      g <- add_edges(g, edges = c(data$splitvarName[i], data$splitvarName[data$rightChild[i] + 1]))
    }
  }
  # plot the graph
  plot(g, layout = layout_as_tree(g), vertex.label = V(g)$label, vertex.size = 12)
}

# function to generate synthetic data around predefined cluster centers
get_synthetic_dataset <- function(cluster_centers, n_samples_per_cluster, sd = 0.2){
  k <- dim(cluster_centers)[1]
  n_feat <- dim(cluster_centers)[2]
  train <- data.frame()
  for (i in 1:k) { cluster_data <- matrix( rnorm(n_feat*n_samples_per_cluster, mean = cluster_centers[i,], sd = sd),
                                           nrow = n_samples_per_cluster, ncol = n_feat, byrow=TRUE )
    cluster_data <- as.data.frame(cluster_data)
    train <- rbind(train, cluster_data) 
  }
  return(train)
}

# function to compute the heaviest clique of a given size from adjacency matrix
get_heaviest_clique <- function(edge_matrix, size){
  # get all possible cliques of a given size
  var_names <- rownames(edge_matrix)
  cliques <- t(combn(var_names, size)) 
  # compute the weight of every clique
  g <- graph_from_adjacency_matrix(as.matrix(edge_matrix), weighted = TRUE, mode = "directed")
  clique_weights <- apply(cliques, 1, function(clique) sum(E(induced_subgraph(graph=simplify(g),vids=clique))$weight))
  return(list(vars = cliques[which.max(clique_weights),], weight = max(clique_weights)))
}

# function to add a feature to the feature set in a greedy manner
select_feature_greedy <- function(edge_matrix, var_added, var_to_add){
  # if no features have been selected, choose the two features connected by the heaviest edge
  if(length(var_added)==0){
    selected <- rownames(which(edge_matrix == max(edge_matrix), arr.ind = TRUE))
    weight <- max(edge_matrix)
  }
  # otherwise choose the feature connected to the selected features by the heaviest edges
  else{
    sub_edge <- edge_matrix[unlist(var_to_add), unlist(var_added), drop=FALSE]
    selected <- rownames(sub_edge)[which.max(rowMeans(sub_edge))]
    weight <- mean(sub_edge[selected,])    
  }
  return(list(feature = selected, avg_weight = weight))
}

# function to select all features in a greedy manner (it assumes that the graph is connected)
select_greedy <- function(edge_matrix){
  # remove self loops and make unweighted
  diag(edge_matrix) <- 0
  edge_matrix <- (edge_matrix + t(edge_matrix))/2
  # get parameters
  n_feat <- nrow(edge_matrix)
  var_names <- rownames(edge_matrix)
  var_added <- list()
  var_to_add <- as.list(var_names)
  # keep track of two types of weights
  # 1 - the average edge weight of the current clique; 2 - the average edge weights of the new feature to the previous clique
  res_greedy_clique <- array(0, dim = c(n_feat-1))
  res_greedy_new <- array(0, dim = c(n_feat-1))
  for (i in 2:n_feat){
    selected <- select_feature_greedy(edge_matrix, var_added, var_to_add)
    var_added <- c(var_added, selected$feature)
    var_to_add <- var_to_add[!var_to_add %in% var_added]
    res_greedy_clique[i-1] <- sum(edge_matrix[unlist(var_added), unlist(var_added)])/(i*(i-1)/2)
    res_greedy_new[i-1] <- selected$avg_weight
  }
  return(list(feature = unlist(var_added), weight_clique = res_greedy_clique, weight_new = res_greedy_new))
}

# function to normalise a 3d matrix along the first two dimensions
norm_matrix <- function(all_edge_matrix){
  n_iter <- dim(all_edge_matrix)[3]
  for (n in 1:n_iter){
    edge_matrix <- all_edge_matrix[,,n] 
    all_edge_matrix[,,n] <- all_edge_matrix[,,n]/sum(all_edge_matrix[,,n] )*100
  }
  return(all_edge_matrix)
}  

# function to scale values in a [0,1] range
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# function to generate cluster centers for synthetic datasets
get_cluster_centers <- function(cluster_number){
  # centrality test
  if (cluster_number==1){
    cluster_center_1 <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_2 <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_3 <- c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_4 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  } 
  # discriminating power test
  if (cluster_number==2){
    cluster_center_1 <- c(1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0)
    cluster_center_2 <- c(0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0)
    cluster_center_3 <- c(0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0)
    cluster_center_4 <- c(1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0)
  } 
  # cluster specific test
  if (cluster_number==3){
    cluster_center_1 <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_2 <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_3 <- c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_4 <- c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  } 
  # repeated features test
  if (cluster_number==4){
    cluster_center_1 <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_2 <- c(0, 0, 1, 1, 0, 0, 0, 0, 0, 0)
    cluster_center_3 <- c(0, 0, 0, 0, 1, 1, 0, 0, 0, 0)
    cluster_center_4 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  } 
  # alternative test for redundant features
  if (cluster_number==5){
    cluster_center_1 <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_2 <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_3 <- c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_4 <- c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  } 
  cluster_centers <- rbind(cluster_center_1, cluster_center_2, cluster_center_3, cluster_center_4)
  return(cluster_centers)
}

# function to generate larger cluster centers for synthetic datasets
get_large_cluster_centers <- function(cluster_number){
  # top-n features experiment
  if (cluster_number==1){
    cluster_center_1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_2 <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_3 <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_4 <- c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_5 <- c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_6 <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_7 <- c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_8 <- c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
  }
  # other top-n features experiments
  if (cluster_number==2){
    cluster_center_1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_2 <- c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_3 <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_4 <- c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_5 <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_6 <- c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_7 <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_8 <- c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  }
  if (cluster_number==3){
    cluster_center_1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_2 <- c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_3 <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_4 <- c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_5 <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_6 <- c(1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_7 <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_8 <- c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  }
  if (cluster_number==4){
    cluster_center_1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_2 <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_3 <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_4 <- c(0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_5 <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_6 <- c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_7 <- c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_8 <- c(1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
  }
  if (cluster_number==5){
    cluster_center_1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_2 <- c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
    cluster_center_3 <- c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_4 <- c(0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_5 <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_6 <- c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
    cluster_center_7 <- c(0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0)
    cluster_center_8 <- c(0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0)
  }
  cluster_centers <- rbind(cluster_center_1, cluster_center_2, cluster_center_3, cluster_center_4,
                           cluster_center_5, cluster_center_6, cluster_center_7, cluster_center_8)
  return(cluster_centers)
} 


# get benchmark datasets
get_dataset <- function(DATASET){
  # -------------
  if(DATASET=="STATLOG"){
    base_name <- file.path("clustering-data-v1", "uci", "statlog")
    train    <- as.matrix(read.table(paste0(base_name, ".data.gz")))
    target  <- scan(paste0(base_name, ".labels0.gz"), integer())
    k = length(unique(target))
  }
  if(DATASET=="YEAST"){
    base_name <- file.path("clustering-data-v1", "uci", "yeast")
    train    <- as.matrix(read.table(paste0(base_name, ".data.gz")))
    target  <- scan(paste0(base_name, ".labels0.gz"), integer())
    k = length(unique(target))
  }
  if(DATASET=="WDBC"){
    base_name <- file.path("clustering-data-v1", "uci", "wdbc")
    train    <- as.matrix(read.table(paste0(base_name, ".data.gz")))
    target  <- scan(paste0(base_name, ".labels0.gz"), integer())
    k = length(unique(target))
  }
  if(DATASET=="SONAR"){
    base_name <- file.path("clustering-data-v1", "uci", "sonar")
    train    <- as.matrix(read.table(paste0(base_name, ".data.gz")))
    target  <- scan(paste0(base_name, ".labels0.gz"), integer())
    k = length(unique(target))
  }
  if(DATASET=="ECOLI"){
    base_name <- file.path("clustering-data-v1", "uci", "ecoli")
    train    <- as.matrix(read.table(paste0(base_name, ".data.gz")))
    target  <- scan(paste0(base_name, ".labels0.gz"), integer())
    k = length(unique(target))
  }
  if(DATASET=="IONOSPHERE"){
    base_name <- file.path("clustering-data-v1", "uci", "ionosphere")
    train    <- as.matrix(read.table(paste0(base_name, ".data.gz")))
    target  <- scan(paste0(base_name, ".labels0.gz"), integer())
    k = length(unique(target))
  }
  if(DATASET=="PARKINSON"){# yes - but bad signal
    my_url="https://archive.ics.uci.edu/ml/machine-learning-databases/parkinsons/parkinsons.data"
    dataset2=read.csv(my_url, header=TRUE, sep=",")
    train = dataset2[,-c(1,18)]
    target = dataset2[,18]
    k = length(unique(target))
  }
  if(DATASET=="HEART"){
    library(kmed)
    data(heart)
    train = heart[, -ncol(heart)]
    train = matrix(as.numeric(unlist(train)), dim(train)[1], dim(train)[2])
    colnames(train) = paste("V", 1:ncol(train), sep="")
    target = heart[,ncol(heart)]
    target[target!=0] = 1
    k = length(unique(target))
  }
  if(DATASET=="IRIS"){
    data(iris)
    train = iris[,-ncol(iris)]
    train = as.matrix(train)
    train = matrix(as.numeric(train), dim(train)[1], dim(train)[2])
    colnames(train) = paste("V", 1:ncol(train), sep="")
    target = iris[,ncol(iris)]
    k = length(unique(target))
  }
  if(DATASET=="WINE"){
    library(HDclassif)
    data(wine)
    train = wine[,-1]
    train = as.matrix(train)
    train = matrix(as.numeric(train), dim(train)[1], dim(train)[2])
    colnames(train) = paste("V", 1:ncol(train), sep="")
    target = wine[,1]
    k = length(unique(target))
  }
  if(DATASET=="GLASS"){
    library(mlbench)
    data(Glass)
    train = Glass[,-10]
    train = as.matrix(train)
    train = matrix(as.numeric(train), dim(train)[1], dim(train)[2])
    colnames(train) = paste("V", 1:ncol(train), sep="")
    target = Glass[,10]
    target = as.numeric(target)
    ids = which(target>=5)
    target[ids] =  4
    k = length(unique(target))
  }
  if(DATASET=="BREAST"){
    train = read.csv(file.path("clustering-data-v1", "uci", "breast_tissue.csv"), header=TRUE)
    target = train[,2]
    train = train[,-c(1,2)]
    train = as.matrix(train)
    train = matrix(as.numeric(train), dim(train)[1], dim(train)[2])
    colnames(train) = paste("V", 1:ncol(train), sep="")
    k = length(unique(target))
  }
  if(DATASET=="LIVER"){
    data <- read.csv(file.path("clustering-data-v1", "uci", "liver_disease.csv"))
    train <- data[, -ncol(data)]
    target <- data[, ncol(data)]
    k <- length(unique(target))
  }
  if(DATASET=="LYMPHO"){
    data <- read.csv(file.path("clustering-data-v1", "uci", "lymphography.csv"))
    train <- data[, -ncol(data)]
    target <- data[, ncol(data)]
    k <- length(unique(target))
  } 
  return(list(train = train, target = target, k = k))
}

# function to calculate centralities for weighted directed graphs
get_wd_centralities <- function(edge_matrix) {
  g1 <- graph.adjacency(edge_matrix, mode = "directed", weighted = TRUE)#, diag = TRUE)
  g2 <- graph.adjacency(t(edge_matrix), mode = "directed", weighted = TRUE)#, diag = TRUE)
  g3 <- graph.adjacency(1/(edge_matrix+1), mode = "directed", weighted = TRUE)
  
  centralities <- data.frame(
    Out_Degree = rowSums(edge_matrix), #degree(g, mode = "out"),
    Out_Eigenvector = evcent(g2, directed = TRUE)$vector,
    Out_Closeness = closeness(g3, mode = "out"),
    In_Degree = colSums(edge_matrix), #degree(g, mode = "in"),
    In_Eigenvector = evcent(g1, directed = TRUE)$vector,
    In_Closeness = closeness(g3, mode = "in")
  )
  return(centralities)
} 

# function to compute silhouette index
compute_silhouette <- function(distance_matrix, cluster_labels) {
  # Compute the silhouette widths
  silhouette_values <- silhouette(as.numeric(cluster_labels), distance_matrix)
  # Return the average silhouette width
  return(mean(silhouette_values[, 'sil_width']))
}

# function to compute Separability Index
compute_separability_index <- function(distance_matrix, cluster_labels) {
  n <- nrow(distance_matrix)
  separability_count <- 0
  # fraction of samples whose nearest neighbour is of the same class
  for (i in 1:n) {
    distances <- distance_matrix[i, ]
    sorted_indices <- order(distances)
    nearest_neighbor_index <- sorted_indices[2]  # nearest neighbor index (excluding self)
    if (cluster_labels[i] == cluster_labels[nearest_neighbor_index]) {
      separability_count <- separability_count + 1
    }
  }
  return(separability_count / n)
}

# function to compute hypothesis margin
compute_hypothesis_margin <- function(distance_matrix, cluster_labels) {
  n <- nrow(distance_matrix)
  margins <- numeric(n)
  for (i in 1:n) {
    distances <- distance_matrix[i, ]
    sorted_indices <- order(distances)
    near_hit <- sorted_indices[distances[sorted_indices] > 0 & cluster_labels[sorted_indices] == cluster_labels[i]][1]
    near_miss <- sorted_indices[distances[sorted_indices] > 0 & cluster_labels[sorted_indices] != cluster_labels[i]][1]
    margins[i] <- distances[near_miss] - distances[near_hit] 
  }
  return(mean(margins, na.rm = TRUE))
}


# functions for phylogeny calculation 
convert_to_phylo = function(cl.obj){
  dend = to.dend(cl.obj)
  dend.str = convert_to_par(dend)
  tree = ape::read.tree(text = dend.str)
  return(tree)
}

convert_to_par <- function(dend, first_it = TRUE)
{
  if (first_it == TRUE){
    dist = as.double(attr(dend, "height"))
    dend.object1 = dend[[1]]
    dist1 = as.double(attr(dend.object1, "height"))
    dend.object2 = dend[[2]]
    dist2 = as.double(attr(dend.object2, "height"))
    first_it = FALSE
    if (is.list(dend.object1) == TRUE & is.list(dend.object2) == TRUE){
      return(paste0("((",
                    convert_to_par(dend.object1, first_it),
                    "):",
                    dist - dist1,
                    ",",
                    "(",
                    convert_to_par(dend.object2, first_it),
                    "):",
                    dist - dist2,
                    ");"))}else{
                      if(is.list(dend.object1) == TRUE & is.list(dend.object2) == FALSE){
                        label = attr(dend.object2, "label")
                        return(paste0("((",
                                      convert_to_par(dend.object1, first_it),
                                      "):",
                                      dist - dist1,
                                      ",",
                                      label,
                                      ":",
                                      dist - dist2,
                                      ");"))
                      }else{
                        label = attr(dend.object1, "label")
                        return(paste0("(",
                                      label,
                                      ":",
                                      dist - dist1,
                                      ",",
                                      "(",
                                      convert_to_par(dend.object2, first_it),
                                      "):",
                                      dist - dist2,
                                      ");"))
                      }
                    }
  }else{
    dist = as.double(attr(dend, "height"))
    dend.object1 =  dend[[1]]
    dist1 = as.double(attr(dend.object1, "height"))
    dend.object2 = dend[[2]]
    dist2 = as.double(attr(dend.object2, "height"))
    if(is.list(dend.object1) == TRUE & is.list(dend.object2) == TRUE){
      return(paste0("(",
                    convert_to_par(dend.object1, first_it),
                    "):", dist - dist1,
                    ",",
                    "(",
                    convert_to_par(dend.object2, first_it),
                    "):", dist - dist2))
    }else{
      if(is.list(dend.object1) == FALSE & is.list(dend.object2) == TRUE){
        label1 = attr(dend.object1,"label")
        return(paste0(label1, ":",
                      dist - dist1,
                      ",",
                      "(",
                      convert_to_par(dend.object2, first_it),
                      "):", dist - dist2))
      }else{
        if(is.list(dend.object1) == TRUE & is.list(dend.object2) == FALSE){
          label2 = attr(dend.object2, "label")
          return(paste0("(",
                        convert_to_par(dend.object1, first_it),
                        "):", dist - dist1, ",",
                        label2, ":",
                        dist - dist2))
        }else{
          label1 = attr(dend.object1, "label")
          label2 = attr(dend.object2, "label")
          return(paste0(label1, ":", dist - dist1,
                        ",",
                        label2, ":", dist - dist2))
        }
      }
    }
  }
}

to.dend = function(cl.obj){
  # converting generical clustering object into hclust object first
  cluster.obj = cl.obj %>% stats::as.hclust()
  dend = cluster.obj %>% stats::as.dendrogram(hang = -1, check = TRUE)
  return(dend)
}
