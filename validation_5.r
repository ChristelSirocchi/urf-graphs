#install.packages("remotes")
#remotes::install_github("Monoxido45/PhyloHclust")
library(PhyloHclust)
source("utilities.r")
library(ranger, lib.loc = "/home/christelsirocchi/R/x86_64-pc-linux-gnu-library/4.3") 

# compute feature importance with the proposed methods and three state-of-the-art methods

datasets = c("IRIS", "BREAST", "ECOLI", "GLASS",  "WINE", "PARKINSON", "LIVER", "LYMPHO", "SONAR", "IONOSPHERE")
n_iter <- 30

##################################################
##### train uRF and compute adjacency matrix #####
##################################################

edge_types =  c("sample") #c("present", "fixation", "level", "sample")

# train all models
for (dataset in datasets){
  # get data and parameters
  data <- get_dataset(dataset)
  print(dataset)
  target <- data$target 
  train <- data$train 
  k <- data$k
  n_feat <- dim(train)[2]
  # store trained models and predictions
  model_list <- list()
  ress = array(0, dim = c(nrow(train), n_iter))
  # train models
  for (n in 1:n_iter) {
    model = ranger(x=train, num.trees=1000,
                   clustering = TRUE,  probability= FALSE, oob.error = FALSE,
                   mtry = round(sqrt(ncol(train))), min.bucket = 5)
    DISTANCE = distance_rf(model, train)
    hc = hclust(as.dist(DISTANCE), method="ward.D2")
    res1 = cutree(hc, k)
    print(ARI(res1, target))
    print(NMI(res1, target))
    ress[,n] <- res1
    model_list[[n]] <- model
  }
  # save trained models and predictions
  colnames(ress) <- paste("iter", seq(n_iter),sep="")
  write.csv(ress, file = file.path("benchmarks", dataset, paste("v5_predictions_uRF", dataset, "1000trees.csv", sep = "_")), row.names = FALSE) 
  saveRDS(model_list, file = file.path("benchmarks", dataset,paste("v5_models_uRF", dataset, "1000trees.rds", sep = "_")))

  # compute all adjacency matrices
  all_edge_matrix <- array(0, dim = c(n_feat+1, n_feat+1, n_iter, length(edge_types)))
  for (n in 1:n_iter) {
    model <- model_list[[n]]
    for (e in 1:length(edge_types)){
      all_edge_matrix[,,n,e] <- compute_edge_matrix(model, train, edge_types[e])
    }
  }
  # save results
  var_names <- model$forest$independent.variable.names
  dimnames(all_edge_matrix) <- list(c(var_names,"T"), c(var_names,"T"), paste("iter", seq(n_iter),sep=""), edge_types)
  edge_df <- melt(all_edge_matrix)
  colnames(edge_df) <- c("node_from", "node_to", "iteration", "edge_type", "value")
  file_name_all_e <- file.path("benchmarks", dataset, paste("v5_all_edges_uRF", dataset, "1000trees.csv", sep = "_"))
  write.csv(edge_df, file = file_name_all_e, row.names = FALSE)
}

#################################################################
### compute feature importance by GREEDY graph based approach ###
#################################################################

for (dataset in datasets){
  print(dataset)
  data <- get_dataset(dataset)
  train <- data$train 
  k <- data$k
  n_feat <- dim(train)[2]
  var_names <- colnames(train)
  # get unsupervised feature importances
  file_name_all_e <- file.path("benchmarks", dataset, paste("v5_all_edges_uRF", dataset, "1000trees.csv", sep = "_"))
  edge_df <- read.csv(file = file_name_all_e)
  all_edge_matrix_all <- acast(edge_df, node_from ~ node_to ~ iteration ~ edge_type)
  all_edge_matrix <- all_edge_matrix_all[,,,"sample"]
  
  # normalise each edge matrix, compute average and do greedy selection
  edge_mean <- apply(norm_matrix(all_edge_matrix[var_names,var_names,]), c(1, 2), mean)
  feat_imp <- select_greedy(edge_mean[var_names,var_names])$feature
  feat_imp_df <- as.data.frame(feat_imp)
  file_name_imp <- file.path("benchmarks", dataset, paste("v5_feature_imp_GREEDY", dataset, "1000trees.csv", sep = "_"))
  write.csv(feat_imp_df, file = file_name_imp, row.names = FALSE) 
}


######################################################################
### compute feature importance by BRUTE-FORCE graph based approach ###
######################################################################

for (dataset in datasets){
  print(dataset)
  data <- get_dataset(dataset)
  train <- data$train 
  k <- data$k
  n_feat <- dim(train)[2]
  var_names <- colnames(train)
  # get unsupervised feature importances
  file_name_all_e <- file.path("benchmarks", dataset, paste("v5_all_edges_uRF", dataset, "1000trees.csv", sep = "_"))
  edge_df <- read.csv(file = file_name_all_e)
  all_edge_matrix_all <- acast(edge_df, node_from ~ node_to ~ iteration ~ edge_type)
  all_edge_matrix <- all_edge_matrix_all[,,,"sample"]
  
  # normalise each edge matrix, compute average and do greedy selection
  edge_mean <- apply(norm_matrix(all_edge_matrix[var_names,var_names,]), c(1, 2), mean)
  cliques = list()
  #res_greedy <- select_greedy(edges[var_names,var_names])$feature
  for (i in 2:min(n_feat,12)){
    print(i)
    cliques[[i]] <- tryCatch(get_heaviest_clique(edge_mean, i)$vars, silent = T, error=function(msg){ return(NA) })  
  }
  saveRDS(cliques, file =file.path("benchmarks", dataset, paste("v5_cliques_greedy", dataset, "1000trees.rds", sep = "_")))
}


###################################################################
### compute feature importance by leave-one-out approach (LOVO) ###
###################################################################

for (dataset in datasets) {
  print(dataset)
  data <- get_dataset(dataset)
  train <- data$train 
  k <- data$k
  n_feat <- dim(train)[2]
  n_row <- dim(train)[1]
  feat_imp <- array(0, dim = c(n_feat, n_iter))
  
  for (n in 1:n_iter) {
    lovo_values <- numeric(n_feat)
    print(n)
    for (x in 1:n_feat) {
      new_data <- train[, -x]
      # obtaining partition for leave one out
      model <- ranger(x = new_data, num.trees = 1000,
                      clustering = TRUE, probability = FALSE, oob.error = FALSE,
                      mtry = round(sqrt(ncol(new_data))), min.bucket = 5)
      
      DISTANCE <- distance_rf(model, new_data)
      hc <- hclust(as.dist(DISTANCE), method = "ward.D2")
      part <- cutree(hc, k)
      # within cluster heterogeneity for each variable
      lovo <- numeric(k)
      for (t in 1:k) {
        cluster_data <- new_data[part==t,]
        trace_t <- sum(diag(var(cluster_data)))
        nt <- nrow(cluster_data)
        lovo[t] <- nt / n_row * trace_t
      }
      lovo_values[x] <- mean(lovo)
    }
    feat_imp[, n] <- lovo_values
  }
  dimnames(feat_imp) <- list(colnames(train), paste("iter", seq(n_iter), sep = ""))
  feat_imp_df <- melt(feat_imp)
  colnames(feat_imp_df) <- c("features", "iteration", "value")
  file_name_imp <- file.path("benchmarks", dataset, paste("v5_feature_imp_LOVO", dataset, "_1000trees.csv", sep = "_"))
  write.csv(feat_imp_df, file = file_name_imp, row.names = FALSE) 
}


################################################################
### compute feature importance by phylogeny approach (PHYLO) ###
################################################################

for (dataset in datasets){
  print(dataset)
  data <- get_dataset(dataset)
  train <- data$train 
  k <- data$k
  n_feat <- dim(train)[2]
  model_list <- readRDS(file.path("benchmarks", dataset, paste("v5_models_uRF", dataset, "1000trees.rds", sep = "_")))
  feat_imp <- array(0, dim = c(n_feat, n_iter))
  for (n in 1:n_iter) {
    print(n)
    model <- model_list[[n]]
    DISTANCE = distance_rf(model, train)
    hc = hclust(as.dist(DISTANCE), method="ward.D2")
    tree <- convert_to_phylo(hc)
    #feat_imp[,n] <- L_score(tree, data.frame(train), score = TRUE, ncores = 2)
    results <- numeric(n_feat)
    for (j in 1:n_feat) {
      test_data = as.data.frame(train[, j])
      colnames(test_data) = colnames(train)[j]
      rownames(test_data) = rownames(train)
      results[j] = tryCatch(L_score(tree, test_data, score = TRUE, ncores = 2), silent = T, error=function(msg){return(NA)})   
    }
    feat_imp[,n] <- results
  }
  dimnames(feat_imp) <- list(model$forest$independent.variable.names, paste("iter", seq(n_iter), sep=""))
  feat_imp_df <- melt(feat_imp)
  colnames(feat_imp_df) <- c("features", "iteration", "value")
  file_name_imp <- file.path("benchmarks", dataset, paste("v5_feature_imp_PHYLO", dataset, "_1000trees.csv", sep = "_"))
  write.csv(feat_imp_df, file = file_name_imp, row.names = FALSE) 
}

################################################################
# compute feature importance by classification (supervised RF) #
################################################################

library(ranger, lib.loc = "/home/christelsirocchi/R/x86_64-pc-linux-gnu-library/RF") # import original ranger

importance_types <- c("impurity", "impurity_corrected")

# train with uRF learnt clusters
for (d in datasets){
  data <- get_dataset(d)
  train <- data$train 
  ress <- read.csv(file = file.path("benchmarks", dataset, paste("v5_predictions_uRF", dataset, "1000trees.csv", sep = "_"))) 
  k <- data$k
  n_feat <- dim(train)[2]
  feat_imp <- array(0, dim = c(n_feat, n_iter, length(importance_types)))
  
  for (e in 1:length(importance_types)){
    for (n in 1:n_iter) {
      model = ranger(x = train, y = as.factor(ress[,n]), 
                     num.trees = 1000, mtry = round(sqrt(ncol(train))), 
                     min.bucket = 5, importance = importance_types[e], oob.error = FALSE)
      feat_imp[,n,e] <- as.vector(model$variable.importance)
    }
  }
  dimnames(feat_imp) <- list(model$forest$independent.variable.names, paste("iter", seq(n_iter), sep=""), importance_types)
  feat_imp_df <- melt(feat_imp)
  colnames(feat_imp_df) <- c("features", "iteration", "importance", "value")
  file_name_all_e <- file.path("benchmarks", dataset, paste("v5_feature_imp_RF", d, "clusters_1000trees.csv", sep = "_"))
  write.csv(feat_imp_df, file = file_name_all_e, row.names = FALSE)    
}


