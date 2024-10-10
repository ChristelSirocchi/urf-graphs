source("utilities.r")
library(ranger, lib.loc = "/home/christelsirocchi/R/x86_64-pc-linux-gnu-library/4.3") 

########################################################
############# GENERATE DATA & TRAIN MODELS #############
########################################################

cluster_number = 1 # top-n features

# get cluster centers and define experiment parameters
n_samples_per_cluster <- 50
cluster_centers <- get_large_cluster_centers(cluster_number)
n_feat <- dim(cluster_centers)[2]
n_iter <- 30

# save all models and datasets as list of lists
model_list_all <- list()
train_list_all <- list()

# first generate all synthetic datasets and train models
for (k in 4:8){
  # generate datasets with increasing number of clusters and relevant features (k clusters, k-1 relevant features)
  cluster_centers_sub <- cluster_centers[rownames(cluster_centers)[1:k],]
  model_list <- list()
  train_list <- list()
  # generate synthetic data and run the model
  for (n in 1:n_iter) {
    train <- get_synthetic_dataset(cluster_centers_sub, n_samples_per_cluster)
    train_list[[n]] <- train
    target <- rep(1:k, each = n_samples_per_cluster)
    
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
  model_list_all[[k]] <- model_list
  train_list_all[[k]] <- train_list
}

# then compute all adjacency matrices
all_edge_matrix <- array(0, dim = c(n_feat+1, n_feat+1, n_iter, length(edge_types), 5))
for (k in 4:8){
  # derive the graph
  for (e in 1:length(edge_types)){ 
    for (n in 1:n_iter) {
      model <- model_list_all[[k]][[n]]
      train <- train_list_all[[k]][[n]]
      all_edge_matrix[,,n,e,k-3] <- compute_edge_matrix(model, train, edge_types[e])
    }
  }
}

# save results
var_names <- model$forest$independent.variable.names
dimnames(all_edge_matrix) <- list(c(var_names,"T"), c(var_names,"T"), paste("iter",seq(n_iter),sep=""), edge_types, paste("sub",seq(5),sep=""))
edge_df <- melt(all_edge_matrix)
colnames(edge_df) <- c("node_from", "node_to", "iteration", "edge_type", "subset", "value")
file_name_all_e <- file.path("synthetic", paste("v3_all_edges_large_cluster_number", cluster_number, "sub_separability.csv", sep = "_"))
write.csv(edge_df, file = file_name_all_e, row.names = FALSE)

saveRDS(train_list_all, file = file.path("synthetic", "v3_train_list_all.rds"))
saveRDS(model_list_all, file = file.path("synthetic", "v3_model_list_all.rds"))


##############################################################
############# READ DATA & ANALYSE FEATURE GRAPHS #############
##############################################################

# import edge files
cluster_number = 1
edge_types = c("present", "fixation", "level", "sample")
file_name_all_e <- file.path("synthetic", paste("v3_all_edges_large_cluster_number", cluster_number, "sub_separability.csv", sep = "_"))
edge_df <- read.csv(file = file_name_all_e)
all_edge_matrix <- acast(edge_df, node_from ~ node_to ~ iteration ~ edge_type ~ subset)
n_feat <- dim(all_edge_matrix)[1]-1
n_iter <- dim(all_edge_matrix)[3]
n_subs <- dim(all_edge_matrix)[5]
var_names <- paste("V", seq(n_feat), sep="")
max_n_feat <- n_feat - 1 # consider all possible cliques up to size n_feat-1

# import models
n_samples_per_cluster <- 50
cluster_centers <- get_large_cluster_centers(cluster_number)
n_iter <- 30

train_list_all <- readRDS(file = file.path("synthetic", "v3_train_list_all.rds"))
model_list_all <- readRDS(file = file.path("synthetic", "v3_model_list_all.rds"))


##################################################################################
######################## BRUTE FORCE SELECTION EXPERIMENT ########################
##################################################################################

res_brute_m <- array(0, dim = c(max_n_feat-1, length(edge_types), n_subs))
for (k in 4:8){
  for (e in 1:length(edge_types)){ 
    # get edge matrix (do not include terminal node)
    edge_matrix <- apply(norm_matrix(all_edge_matrix[var_names,var_names,,e,k-3]), c(1, 2), mean)
    # for each matrix and each clique size compute the average weight of the clique
    for (i in 2:max_n_feat){
      max_clique <- get_heaviest_clique(edge_matrix, i)
      res_brute_m[i-1,e,k-3] <- max_clique$weight/(i*(i-1))
    }
  }
}
# save data
dimnames(res_brute_m) <- list(2:max_n_feat, edge_types, paste(3:7, "relevant", sep=" "))
brute_df_m = melt(res_brute_m)
colnames(brute_df_m) <- c("relevant", "edge_type", "selected_features", "value")
file_name_all_e <- file.path("synthetic", paste("v3_all_edges_large_cluster_number", cluster_number, "brute_mean.csv", sep = "_"))
write.csv(brute_df_m, file = file_name_all_e, row.names = FALSE)

######################## SEPARABILITY COMPUTATION on SAMPLE configuration ########

separability <- array(0, dim = c(11,n_iter,3,5))
for (k in 4:8){
  print(k)
  cluster_labels <- rep(1:k, each = n_samples_per_cluster)
  # get edge matrix (do not include terminal node)
  edge_matrix <- apply(norm_matrix(all_edge_matrix[var_names,var_names,,4,k-3]), c(1, 2), mean)
  # for each matrix and each clique size compute the average weight of the clique
  for (i in 2:max_n_feat){
    sel_vars <- get_heaviest_clique(edge_matrix, i)$vars
    for (n in 1:n_iter) {      
      train <- train_list_all[[k]][[n]]
      sub_train <- train[sel_vars]
      # compute distance over the selected features
      model = ranger(x=sub_train, num.trees=500, 
                     clustering = TRUE,  probability= FALSE, oob.error = FALSE,
                     mtry = round(sqrt(ncol(sub_train))), min.bucket = 5)
      distance_matrix <- distance_rf(model, sub_train)
      # compute separability indices
      separability[i-1,n,1,k-3] <- compute_silhouette(distance_matrix, cluster_labels)
      separability[i-1,n,2,k-3] <- compute_separability_index(distance_matrix, cluster_labels)
      separability[i-1,n,3,k-3] <- compute_hypothesis_margin(distance_matrix, cluster_labels)
    }
  }
}
sep_brute <- separability
dimnames(sep_brute) <- list(2:12, paste("iter",seq(n_iter),sep=""), 
                            c("silhouette", "separability index", "hypothesis margin"), paste(3:7, "relevant", sep=" "))
sep_brute_df = melt(sep_brute)
colnames(sep_brute_df) <- c("selected_features", "iteration", "metric", "relevant", "value")
write.csv(sep_brute_df, file = file.path("synthetic", "v3_all_edges_brute_separability.csv"), row.names = FALSE)


#############################################################################
######################## GREEDY SELECTION EXPERIMENT ########################
#############################################################################

res_greedy_m <- array(0, dim = c(max_n_feat-1, length(edge_types), n_subs))
for (k in 4:8){
  #print(k)
  for (e in 1:length(edge_types)){ 
    for (n in 1:n_iter){
      # get edge matrix (do not include terminal node)
      edge_matrix <- apply(norm_matrix(all_edge_matrix[var_names,var_names,,e,k-3]), c(1, 2), mean)
      # consider graph unweighted
      edge_matrix <- (edge_matrix + t(edge_matrix))/2
      # don't consider self-loops
      diag(edge_matrix) <- 0
      # for each matrix select top-n features with a greedy approach
      var_added <- list()
      var_to_add <- as.list(var_names)
      for (i in 2:max_n_feat){
        selected <- select_feature_greedy(edge_matrix, var_added, var_to_add)
        var_added <- c(selected$feature, var_added)
        var_to_add <- var_to_add[!var_to_add %in% var_added]
        # compute average weight of the selected subgraph
        res_greedy_m[i-1,e,k-3] <- sum(edge_matrix[unlist(var_added), unlist(var_added)])/(i*(i-1))
      }
    }
  }
}
# save data
dimnames(res_greedy_m) <- list(2:max_n_feat, edge_types, paste(3:7, "relevant", sep=" "))
greedy_df_m = melt(res_greedy_m)
colnames(greedy_df_m) <- c("relevant", "edge_type", "selected_features", "value")
file_name_all_e <- file.path("synthetic", paste("v3_all_edges_large_cluster_number", cluster_number, "greedy_mean.csv", sep = "_"))
write.csv(greedy_df_m, file = file_name_all_e, row.names = FALSE)

######################## SEPARABILITY COMPUTATION on SAMPLE configuration ########

separability <- array(0, dim = c(11,n_iter,3,5))
for (k in 4:8){
  print(k)
  cluster_labels <- rep(1:k, each = n_samples_per_cluster)
  # get edge matrix (do not include terminal node)
  edge_matrix <- apply(norm_matrix(all_edge_matrix[var_names,var_names,,4,k-3]), c(1, 2), mean)
  # consider graph unweighted
  edge_matrix <- (edge_matrix + t(edge_matrix))/2
  # don't consider self-loops
  diag(edge_matrix) <- 0
  # for each matrix select top-n features with a greedy approach
  var_added <- list()
  var_to_add <- as.list(var_names)
  
  for (i in 2:max_n_feat){
    selected <- select_feature_greedy(edge_matrix, var_added, var_to_add)
    var_added <- c(selected$feature, var_added)
    var_to_add <- var_to_add[!var_to_add %in% var_added]
    sel_vars <- unlist(var_added)
    print(sel_vars)
    for (n in 1:n_iter) {      
      train <- train_list_all[[k]][[n]]
      sub_train <- train[sel_vars]
      model = ranger(x=sub_train, num.trees=500, 
                     clustering = TRUE,  probability= FALSE, oob.error = FALSE,
                     mtry = round(sqrt(ncol(sub_train))), min.bucket = 5)
      distance_matrix <- distance_rf(model, sub_train)
      separability[i-1,n,1,k-3] <- compute_silhouette(distance_matrix, cluster_labels)
      separability[i-1,n,2,k-3] <- compute_separability_index(distance_matrix, cluster_labels)
      separability[i-1,n,3,k-3] <- compute_hypothesis_margin(distance_matrix, cluster_labels)
    }
  }
}
sep_greedy <- separability
dimnames(sep_greedy) <- list(2:12, paste("iter",seq(n_iter),sep=""),
                             c("silhouette", "separability index", "hypothesis margin"), paste(3:7, "relevant", sep=" "))
sep_greedy_df = melt(sep_greedy)
colnames(sep_greedy_df) <- c("selected_features", "iteration", "metric", "relevant", "value")
write.csv(sep_greedy_df, file = file.path("synthetic", "v3_all_edges_greedy_separability.csv"), row.names = FALSE)

#####################################
########### READ AND PLOT ###########
#####################################

cluster_number = 1 

#read data
brute_df <- read.csv(file = file.path("synthetic", paste("v3_all_edges_large_cluster_number", cluster_number, "brute_mean.csv", sep = "_")))
brute_df$method <- "brute-force"
greedy_df <- read.csv(file = file.path("synthetic", paste("v3_all_edges_large_cluster_number", cluster_number, "greedy_mean.csv", sep = "_")))
greedy_df$method <- "greedy"

all_df <- rbind(brute_df,greedy_df)
all_df$edge_type <- factor(all_df$edge_type, levels = edge_types)
ggplot(all_df, aes(x = relevant, y = value, colour = edge_type, group = edge_type)) + 
  geom_point(size = 1.2) + geom_line() + theme_bw() +
  scale_x_continuous(breaks = seq(2, 12, by = 2)) +
  facet_grid(~ method ~ selected_features, scales = "free") +
  labs( y = "average edge weight (AW)", x = "number of selected features") + 
  guides(color=guide_legend("edge type")) +
  theme(legend.position = "bottom", legend.title = element_blank())


######### SEPARABILITY PLOT ##########

sep_brute_df <- read.csv(file = file.path("synthetic", "v3_all_edges_brute_separability.csv"))
separability <- acast(sep_brute_df, selected_features ~ iteration ~ metric ~ relevant)
separability_mean <- apply(separability, c(1, 3, 4), mean)
separability_norm <- apply(separability_mean, c(2, 3), function(x) (x-min(x)) / (max(x)-min(x)))
sep_df1 = melt(separability_norm)
colnames(sep_df1) <- c("selected_features", "metric", "relevant", "value")
sep_df1$metric <- factor(sep_df1$metric, levels = c("silhouette", "separability index", "hypothesis margin"))

sep_greedy_df <- read.csv(file = file.path("synthetic", "v3_all_edges_greedy_separability.csv"))
separability <- acast(sep_greedy_df, selected_features ~ iteration ~ metric ~ relevant)
separability_mean <- apply(separability, c(1, 3, 4), mean)
separability_norm <- apply(separability_mean, c(2, 3), function(x) (x-min(x)) / (max(x)-min(x)))
sep_df2 = melt(separability_norm)
colnames(sep_df2) <- c("selected_features", "metric", "relevant", "value")

sep_df1$method = "brute-graph"
sep_df2$method = "greedy-graph"
sep_df = rbind(sep_df1,sep_df2)

ggplot(sep_df, aes(x = selected_features, y = value, color = metric, shape = metric, group = metric)) +
  geom_point(size = 1.1) + geom_line(size = 0.3) + theme_bw() + scale_x_continuous(breaks = seq(2, 12, by = 2)) +
  facet_grid(method ~ relevant, scales = "free") +
  labs(y = "separability metric value", x = "number of selected features") + 
  theme(legend.position = "bottom", legend.title = element_blank())

