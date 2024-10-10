source("utilities.r")
library(ranger, lib.loc = "/home/christelsirocchi/R/x86_64-pc-linux-gnu-library/4.3") 

########################################################
############# GENERATE DATA & TRAIN MODELS #############
########################################################

cluster_number = 3 # cluster-specific features

# get cluster centers and define experiment parameters
cluster_centers <- get_cluster_centers(cluster_number)
k <- dim(cluster_centers)[1]
n_feat <- dim(cluster_centers)[2]
n_samples_per_cluster <- 50  
n_iter <- 30

# save training data, predictions and models
ress <- array(0, dim = c(n_samples_per_cluster*k, n_iter))
model_list <- list()
train_list <- list()

# generate synthetic data around predefined cluster centers and train uRF
for (n in 1:n_iter) {
  # generate synthetic data 
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
  ress[,n] <- res1
}

# compute edge matrix for each dataset, edge type, and cluster
edge_types = c("present", "fixation", "level", "sample")
all_edge_matrix <- array(0, dim = c(n_feat+1, n_feat+1, k, n_iter, length(edge_types)))

for (e in 1:length(edge_types)){ 
  for (n in 1:n_iter) {
    print(n)
    model <- model_list[[n]]
    train <- train_list[[n]]
    for (cluster in 1:k){ 
      # get indices of samples mapped to a given cluster
      cluster_index <- which(ress[,n] == cluster)
      all_edge_matrix[,,cluster,n,e] <- compute_edge_matrix(model, train, edge_types[e], cluster_index)
    }
  }
}

# save results
var_names <- model$forest$independent.variable.names
dimnames(all_edge_matrix) <- list(c(var_names,"T"), c(var_names,"T"), paste("cluster",seq(k),sep=""), paste("iter",seq(n_iter),sep=""), edge_types)
edge_df <- melt(all_edge_matrix)
colnames(edge_df) <- c("node_from", "node_to", "cluster", "iteration", "edge_type", "value")
file_name_all_e <- file.path("synthetic", paste("v2_all_edges_cluster_number", cluster_number, "cluster_.csv", sep = "_"))
write.csv(edge_df, file = file_name_all_e, row.names = FALSE)

stop("finished computation")


##############################################################
############# READ DATA & ANALYSE FEATURE GRAPHS #############
##############################################################

cluster_number = 3
edge_types = c("present", "fixation", "level", "sample")
file_name_all_e <- file.path("synthetic", paste("v2_all_edges_cluster_number", cluster_number, "cluster_.csv", sep = "_"))
edge_df <- read.csv(file = file_name_all_e)
all_edge_matrix <- acast(edge_df, node_from ~ node_to ~ cluster ~ iteration ~ edge_type)
n_feat <- dim(all_edge_matrix)[1]-1
var_names <- paste("V", seq(n_feat), sep="")
n_clusters <- dim(all_edge_matrix)[3]
n_iter <- dim(all_edge_matrix)[4]
all_edge_matrix <- all_edge_matrix[c(var_names,"T"), c(var_names,"T"),,,]

# compute out-degree centralities for each node in each generated network
all_out_degree_centralities <- array(0, dim = c(length(c(var_names,"T")), n_clusters, n_iter, length(edge_types)))
for (e in 1:length(edge_types)){
  for (n in 1:n_iter) {
    for (c in 1:n_clusters){
      edge_matrix <- all_edge_matrix[,,c,n,e]
      all_out_degree_centralities[,c,n,e] <- rowSums(edge_matrix)/ sum(rowSums(edge_matrix)) *100    
    }
  }
}

dimnames(all_out_degree_centralities) <- list(c(var_names,"T"),paste("cluster",seq(n_clusters),sep=""),paste("iter",seq(n_iter),sep=""),edge_types)
all_cent_df <- melt(all_out_degree_centralities)
colnames(all_cent_df) <- c("node", "cluster", "iteration", "edge_type", "value")
all_cent_df = all_cent_df[all_cent_df$node != "T",]

# group nodes in (a) cluster-specific features (b) sub-relevant features (c) not-relevant features
all_cent_df$relevance <- "sub-relevant"
all_cent_df[all_cent_df$node %in% c("V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "T"),"relevance"] <-"irrelevant"
all_cent_df[(all_cent_df$node == "V1") & (all_cent_df$cluster == "cluster1"),"relevance"] <-"cluster-specific"
all_cent_df[(all_cent_df$node == "V2") & (all_cent_df$cluster == "cluster2"),"relevance"] <-"cluster-specific"
all_cent_df[(all_cent_df$node == "V3") & (all_cent_df$cluster == "cluster3"),"relevance"] <-"cluster-specific"
all_cent_df[(all_cent_df$node == "V4") & (all_cent_df$cluster == "cluster4"),"relevance"] <-"cluster-specific"
all_cent_df$relevance <- factor(all_cent_df$relevance, levels = c("cluster-specific","sub-relevant","irrelevant"))

# plot centralities                                
ggplot(all_cent_df, aes(x = cluster, y = value, fill = relevance)) + 
  geom_boxplot(outlier.size = 0.5) +
  labs(#title = "cluster specific node centrality",
       y = "out-degree centrality", x = "clusters") + theme_bw() +
  coord_cartesian(ylim = c(0, 21)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(size=12)) +
  facet_wrap(~ edge_type, nrow = 1) + theme(legend.position = "bottom", legend.title = element_blank())
