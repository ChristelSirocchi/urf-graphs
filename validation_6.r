source("utilities.r")
library(ranger, lib.loc = "/home/christelsirocchi/R/x86_64-pc-linux-gnu-library/4.3") 

datasets = c("IRIS", "BREAST", "ECOLI", "GLASS",  "WINE", "PARKINSON", "LIVER", "LYMPHO", "SONAR", "IONOSPHERE")
n_iter <- 30

#########################################################
### RETRAIN uRF on TOP k features GREEDY ################
#########################################################

for (dataset in datasets){
  print(dataset)
  data <- get_dataset(dataset)
  train <- data$train
  target <- data$target
  k <- data$k
  var_names <- colnames(train)
  n_rows <- dim(train)[1]
  
  file_name_all_e <- file.path("benchmarks", dataset, paste("v5_feature_imp_GREEDY", dataset, "1000trees.csv", sep = "_"))
  feat_imp_df <- read.csv(file = file_name_all_e, header = TRUE)  
  ordered_feat <- feat_imp_df$feat_imp
  
  max_feat <- min(length(var_names), 12)
  min_feat <- 2
  
  aris <- array(0, dim = c(n_iter, max_feat-min_feat+1, 2))
  for (v in min_feat:max_feat){
    print(v)
    sup_selected <- ordered_feat[1:v]
    print(sup_selected)
    for (n in 1:n_iter) {
      model = ranger(x = train[,sup_selected], num.trees = 500,
                     clustering = TRUE,  probability = FALSE, oob.error = FALSE,
                     mtry = round(sqrt(ncol(train[,sup_selected]))), min.bucket = 5)
      
      DISTANCE = distance_rf(model, train[,sup_selected])
      hc = hclust(as.dist(DISTANCE), method="ward.D2")
      res1 = cutree(hc, k)
      aris[n,v-min_feat+1,1] <- ARI(res1, target)
      aris[n,v-min_feat+1,2] <- NMI(res1, target)
      print(ARI(res1, target))
    }
  }  
  # save to file
  dimnames(aris) <- list(paste("iter", seq(n_iter), sep=""), min_feat:max_feat, c("ARI","NMI"))
  aris_df <- melt(aris)
  aris_df$method <- "greedy" 
  colnames(aris_df) <- c("iteration", "n_feat", "metric", "value", "method")
  aris_df  <-  aris_df[c("iteration", "n_feat", "method", "metric", "value")]
  file_name_all_e <- file.path("benchmarks", dataset, paste("v5_selected_aris_greedy", dataset, "500trees.csv", sep = "_"))
  write.csv(aris_df, file = file_name_all_e, row.names = FALSE)    
}


##############################################################
### RETRAIN uRF on TOP k features BRUTE FORCE ################
##############################################################

#datasets = c("IRIS", "BREAST", "ECOLI", "GLASS",  "LIVER", "WINE", "LYMPHO", "PARKINSON")

for (dataset in datasets){
  print(dataset)
  data <- get_dataset(dataset)
  train <- data$train
  target <- data$target
  k <- data$k
  var_names <- colnames(train)
  n_rows <- dim(train)[1]
  cliques <- readRDS(file.path("benchmarks", dataset, paste("v5_cliques_greedy", dataset, "1000trees.rds", sep = "_")))
  max_feat <- min(length(var_names), 12)
  min_feat <- 2
  
  aris <- array(0, dim = c(n_iter, max_feat-min_feat+1, 2))
  for (v in min_feat:max_feat){
    print(v)
    sup_selected <- cliques[[v]]#ordered_feat[1:v]
    print(sup_selected)
    for (n in 1:n_iter) {
      model = ranger(x = train[,sup_selected], num.trees = 500,
                     clustering = TRUE,  probability = FALSE, oob.error = FALSE,
                     mtry = round(sqrt(ncol(train[,sup_selected]))), min.bucket = 5)
      
      DISTANCE = distance_rf(model, train[,sup_selected])
      hc = hclust(as.dist(DISTANCE), method="ward.D2")
      res1 = cutree(hc, k)
      aris[n,v-min_feat+1,1] <- ARI(res1, target)
      aris[n,v-min_feat+1,2] <- NMI(res1, target)
      print(ARI(res1, target))
    }
  }  
  # save to file
  dimnames(aris) <- list(paste("iter", seq(n_iter), sep=""), min_feat:max_feat, c("ARI","NMI"))
  aris_df <- melt(aris)
  aris_df$method <- "brute" 
  colnames(aris_df) <- c("iteration", "n_feat", "metric", "value", "method")
  aris_df  <-  aris_df[c("iteration", "n_feat", "method", "metric", "value")]
  file_name_all_e <- file.path("benchmarks", dataset, paste("v5_selected_aris_brute", dataset, "500trees.csv", sep = "_"))
  write.csv(aris_df, file = file_name_all_e, row.names = FALSE)    
}


################################################################
### RETRAIN uRF on TOP k features LEAVE-ONE-OUT ################
################################################################

for (dataset in datasets){
  print(dataset)
  data <- get_dataset(dataset)
  train <- data$train
  target <- data$target
  k <- data$k
  var_names <- colnames(train)
  n_rows <- dim(train)[1]
  
  # read lovo feature importances
  file_name_all_e <- file.path("benchmarks", dataset, paste("v5_feature_imp_LOVO", dataset, "_1000trees.csv", sep = "_"))
  feat_imp_df <- read.csv(file = file_name_all_e, header = TRUE)   
  feat_imp <- acast(feat_imp_df, features ~ iteration)
  lovo_imp <- rowMeans(feat_imp)
  sup_lovo_imp <- var_names[order(lovo_imp, decreasing = TRUE)]
  max_feat <- min(length(var_names), 12)
  min_feat <- 2
  
  aris <- array(0, dim = c(n_iter, max_feat-min_feat+1, 2))
  for (v in min_feat:max_feat){
    print(v)
    sup_selected <- sup_lovo_imp[1:v]
    print(sup_selected)
    for (n in 1:n_iter) {
      model = ranger(x = train[,sup_selected], num.trees = 500,
                     clustering = TRUE,  probability = FALSE, oob.error = FALSE,
                     mtry = round(sqrt(ncol(train[,sup_selected]))), min.bucket = 5)
      
      DISTANCE = distance_rf(model, train[,sup_selected])
      hc = hclust(as.dist(DISTANCE), method="ward.D2")
      res1 = cutree(hc, k)
      aris[n,v-min_feat+1,1] <- ARI(res1, target)
      aris[n,v-min_feat+1,2] <- NMI(res1, target)
      print(ARI(res1, target))
    }
  }  
  # save to file
  dimnames(aris) <- list(paste("iter", seq(n_iter), sep=""), min_feat:max_feat, c("ARI","NMI"))
  aris_df <- melt(aris)
  aris_df$method <- "lovo" 
  colnames(aris_df) <- c("iteration", "n_feat", "metric", "value", "method")
  aris_df  <-  aris_df[c("iteration", "n_feat", "method", "metric", "value")]
  file_name_all_e <- file.path("benchmarks", dataset, paste("v5_selected_aris_lovo", dataset, "500trees.csv", sep = "_"))
  write.csv(aris_df, file = file_name_all_e, row.names = FALSE)    
}


###############################################################
### RETRAIN uRF on TOP k features PHYLOGENETIC ################
###############################################################

for (dataset in datasets){
  print(dataset)
  data <- get_dataset(dataset)
  train <- data$train
  target <- data$target
  k <- data$k
  var_names <- colnames(train)
  n_rows <- dim(train)[1]
  
  # read phylogeny feature importances
  file_name_all_e <- file.path("benchmarks", dataset, paste("v5_feature_imp_PHYLO", dataset, "_1000trees.csv", sep = "_"))
  feat_imp_df <- read.csv(file = file_name_all_e, header = TRUE)  
  feat_imp <- acast(feat_imp_df, features ~ iteration)
  phylo_imp <- 1-rowMeans(feat_imp)
  sup_phylo_imp <- var_names[order(phylo_imp, decreasing = TRUE)]

  max_feat <- min(length(var_names), 12)
  min_feat <- 2
  
  aris <- array(0, dim = c(n_iter, max_feat-min_feat+1, 2))
  for (v in min_feat:max_feat){
    print(v)
    sup_selected <- sup_phylo_imp[1:v]
    print(sup_selected)
    for (n in 1:n_iter) {
      model = ranger(x = train[,sup_selected], num.trees = 500,
                     clustering = TRUE,  probability = FALSE, oob.error = FALSE,
                     mtry = round(sqrt(ncol(train[,sup_selected]))), min.bucket = 5)
      
      DISTANCE = distance_rf(model, train[,sup_selected])
      hc = hclust(as.dist(DISTANCE), method="ward.D2")
      res1 = cutree(hc, k)
      aris[n,v-min_feat+1,1] <- ARI(res1, target)
      aris[n,v-min_feat+1,2] <- NMI(res1, target)
      print(ARI(res1, target))
    }
  }  
  # save to file
  dimnames(aris) <- list(paste("iter", seq(n_iter), sep=""), min_feat:max_feat, c("ARI","NMI"))
  aris_df <- melt(aris)
  aris_df$method <- "phylo" 
  colnames(aris_df) <- c("iteration", "n_feat", "metric", "value", "method")
  aris_df  <-  aris_df[c("iteration", "n_feat", "method", "metric", "value")]
  file_name_all_e <- file.path("benchmarks", dataset, paste("v5_selected_aris_phylo_up", dataset, "500trees.csv", sep = "_"))
  write.csv(aris_df, file = file_name_all_e, row.names = FALSE)    
}


#################################################################
### RETRAIN uRF on TOP k features CLASSIFICATION ################
#################################################################

for (dataset in datasets){
  print(dataset)
  data <- get_dataset(dataset)
  train <- data$train
  target <- data$target
  k <- data$k
  var_names <- colnames(train)
  n_rows <- dim(train)[1]
  
  # read phylogeny feature importances
  file_name_all_e <- file.path("benchmarks", dataset, paste("v5_feature_imp_RF", dataset, "clusters_1000trees.csv", sep = "_"))
  feat_imp_df <- read.csv(file = file_name_all_e, header = TRUE)   
  feat_imp <- acast(feat_imp_df, features ~ iteration ~ importance)
  sup_corrected <- rowMeans(feat_imp[,,"impurity_corrected"])[var_names]
  sup_clusters_corrected <- var_names[order(sup_corrected, decreasing = TRUE)]

  max_feat <- min(length(var_names), 12)
  min_feat <- 2
  
  aris <- array(0, dim = c(n_iter, max_feat-min_feat+1, 2))
  for (v in min_feat:max_feat){
    print(v)
    sup_selected <- sup_clusters_corrected[1:v]
    print(sup_selected)
    for (n in 1:n_iter) {
      model = ranger(x = train[,sup_selected], num.trees = 500,
                     clustering = TRUE,  probability = FALSE, oob.error = FALSE,
                     mtry = round(sqrt(ncol(train[,sup_selected]))), min.bucket = 5)
      
      DISTANCE = distance_rf(model, train[,sup_selected])
      hc = hclust(as.dist(DISTANCE), method="ward.D2")
      res1 = cutree(hc, k)
      aris[n,v-min_feat+1,1] <- ARI(res1, target)
      aris[n,v-min_feat+1,2] <- NMI(res1, target)
      print(ARI(res1, target))
    }
  }  
  # save to file
  dimnames(aris) <- list(paste("iter", seq(n_iter), sep=""), min_feat:max_feat, c("ARI","NMI"))
  aris_df <- melt(aris)
  aris_df$method <- "sup_corrected" 
  colnames(aris_df) <- c("iteration", "n_feat", "metric", "value", "method")
  aris_df  <-  aris_df[c("iteration", "n_feat", "method", "metric", "value")]
  file_name_all_e <- file.path("benchmarks", dataset, paste("v5_selected_aris_rf", dataset, "500trees.csv", sep = "_"))
  write.csv(aris_df, file = file_name_all_e, row.names = FALSE)    
}

###############################################################################################
##################################### STATISTICAL TESTING #####################################
###############################################################################################

small_datasets = c("IRIS", "BREAST", "ECOLI", "GLASS", "LIVER")
large_datasets = c("PARKINSON", "SONAR", "LYMPHO", "IONOSPHERE", "WINE")

all_res = data.frame()
for (dataset in datasets){
  data <- get_dataset(dataset)
  n_feat <- dim(data$train)[2]
  file_name_all_e <- file.path("benchmarks", dataset, paste("v5_selected_aris_greedy", dataset, "500trees.csv", sep = "_"))
  aris_df0 <- read.csv(file = file_name_all_e)    
  file_name_all_e <- file.path("benchmarks", dataset, paste("v5_selected_aris_rf", dataset, "500trees.csv", sep = "_"))
  aris_df1 <- read.csv(file = file_name_all_e)  
  file_name_all_e <- file.path("benchmarks", dataset, paste("v5_selected_aris_phylo_up", dataset, "500trees.csv", sep = "_"))
  aris_df2 <- read.csv(file = file_name_all_e)
  file_name_all_e <- file.path("benchmarks", dataset, paste("v5_selected_aris_lovo", dataset, "500trees.csv", sep = "_"))
  aris_df3 <- read.csv(file = file_name_all_e)
  aris_df <- rbind(aris_df0, aris_df1, aris_df2, aris_df3)
  f <- file.path("benchmarks", dataset, paste("v5_selected_aris_brute", dataset, "500trees.csv", sep = "_"))
  if (file.exists(f)){
    aris_df4 <- read.csv(file = f)
    aris_df <- rbind(aris_df, aris_df4)
  }
  max_feat <- min(n_feat-1, 12)
  aris_df$sizedataset <- ifelse(dataset %in% small_datasets, "small", "large")
  aris_df$dataset = dataset 
  aris_df$method <- factor(aris_df$method, levels=c("unsup_network", "brute", "sup_corrected", "phylo", "lovo"))
  levels(aris_df$method) <- c("greedy", "brute", "classification", "phylogenetic", "leave-one-variable-out")
  all_res <- rbind(all_res, aris_df)
}

##################################### PAIRED WILCOXON TEST ####################################

options(warn=-1)
res_table <- NULL
names <- NULL
all_res_m <- all_res[(all_res$metric=="NMI"),] # or ARI

for (dataset in datasets){
  data <- get_dataset(dataset)
  n_feat <- dim(data$train)[2]
  max_feat <- min(n_feat-1, 12)
  aris_sub <- all_res_m[all_res_m$dataset == dataset, ]
  aris_sub <- aris_sub[(aris_sub$n_feat<=max_feat),]
  # Perform pairwise Wilcoxon tests
  pairwise_result <- pairwise.wilcox.test(aris_sub$value, aris_sub$method, p.adjust.method = "bonferroni", paired = TRUE)
  agg <- aggregate(value ~ method, data = aris_sub, FUN = mean)
  means <- agg$value
  names(means) <- agg$method
  row1 <- means[c("greedy", "brute", "classification", "phylogenetic", "leave-one-variable-out")]
  row2 <- pairwise_result$p.value[,1][c("greedy", "brute", "classification", "phylogenetic", "leave-one-variable-out")]
  names <- c(names, paste0(dataset, " mean"), paste0(dataset," p-value"))
  res_table <- rbind(res_table, row1, row2)
}
rownames(res_table) <- names
colnames(res_table) <- c("greedy", "brute", "classification", "phylogenetic", "leave-one-variable-out") 
round(res_table,4)

##################################### MONOTONICITY #####################################

methods = c("greedy","classification", "phylogenetic", "leave-one-variable-out")

all_res_m <- all_res[(all_res$metric=="ARI"),]
all_res_m <- all_res_m[all_res_m$method %in% methods,]
all_mono <- array(0, dim = c(n_iter, length(methods), length(datasets)))

for (d in seq_along(datasets)){
  dataset = datasets[d]
  data <- get_dataset(dataset)
  n_feat <- dim(data$train)[2]
  max_feat <- min(n_feat, 12)
  data_sub <- all_res_m[all_res_m$dataset == dataset, ]
  data_sub <- data_sub[(data_sub$n_feat<=max_feat),]
  m_sub <- acast(data_sub, iteration ~ n_feat ~ method)
  all_mono[,,d] <- apply(m_sub, c(1,3), compute_weighted_monotonicity)
}
dimnames(all_mono) <- list(paste("iter", seq(n_iter), sep=""), methods, datasets) 
round(apply(all_mono, c(2,3), mean),4)
round(apply(all_mono, 2, mean),4)

mono_df <- melt(all_mono)
colnames(mono_df) <- c("iteration", "method", "dataset", "value")
pairwise.wilcox.test(mono_df$value, mono_df$method, p.adjust.method = "bonferroni", paired = TRUE)
