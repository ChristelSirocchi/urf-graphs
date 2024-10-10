source("utilities.r")
library(ranger, lib.loc = "uRF_PATH") # change with the path of the directory where the custom ranger library was installed

########################################################
############# GENERATE DATA & TRAIN MODELS #############
########################################################

cluster_number = 2 #or 1
# 1 for centrality study
# 2 for edge weight study

# get cluster centers and define experiment parameters
cluster_centers <- get_cluster_centers(cluster_number)
k <- dim(cluster_centers)[1]
n_feat <- dim(cluster_centers)[2]
n_samples_per_cluster <- 50  
n_iter <- 30

# ress <- array(0, dim = c(n_samples_per_cluster*k, n_iter)) # save clustering prediction if needed
model_list <- list()
train_list <- list()

# generate synthetic data around predefined cluster centers and train uRF
for (n in 1:n_iter) {
  train <- get_synthetic_dataset(cluster_centers, n_samples_per_cluster)
  train_list[[n]] <- train
  target <- rep(1:k, each = n_samples_per_cluster)
  #ggplot(train, aes(x = V2, y = V3, color = as.factor(target))) + geom_point() # plot data points if needed
  
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

# compute edge matrix for each dataset and edge type
edge_types = c("present", "fixation", "level", "sample")
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
var_names <- model$forest$independent.variable.names 
dimnames(all_edge_matrix) <- list(c(var_names,"T"), c(var_names,"T"), paste("iter", seq(n_iter), sep=""), edge_types)
edge_df <- melt(all_edge_matrix)
colnames(edge_df) <- c("node_from", "node_to", "iteration", "edge_type", "value")
file_name_all_e <- file.path("synthetic", paste("v1_all_edges_cluster_number", cluster_number, ".csv", sep = "_"))
write.csv(edge_df, file = file_name_all_e, row.names = FALSE)

stop("finished computation")


##############################################################
############# READ DATA & ANALYSE FEATURE GRAPHS #############
##############################################################

############## 1. CENTRALITY VALIDATION EXPERIMENT ##############

cluster_number = 1
file_name_all_e <- file.path("synthetic", paste("v1_all_edges_cluster_number", cluster_number, ".csv", sep = "_"))
edge_df <- read.csv(file = file_name_all_e)
all_edge_matrix <- acast(edge_df, node_from ~ node_to ~ iteration ~ edge_type) 
var_names <- paste("V", seq(dim(all_edge_matrix)[1]-1), sep="")
all_edge_matrix <- all_edge_matrix[c(var_names,"T"), c(var_names,"T"), , ] # reorder variables

# compute out-degree centralities (row-sums in the adjacency matrix) for each dataset and edge type
n_iter = 30
edge_types = c("present", "fixation", "level", "sample")
all_out_degree_centralities <- array(0, dim = c(length(c(var_names,"T")), n_iter, length(edge_types)))
for (e in 1:length(edge_types)){
  for (n in 1:n_iter) {
    edge_matrix <- all_edge_matrix[,,n,e]
    # compute normalised out-degree centrality
    all_out_degree_centralities[,n,e] <- rowSums(edge_matrix)/ sum(rowSums(edge_matrix)) *100
  }
}
dimnames(all_out_degree_centralities) <- list(c(var_names,"T"), paste("iter",seq(n_iter),sep=""), edge_types)
all_cent_df <- melt(all_out_degree_centralities)
colnames(all_cent_df) <- c("node", "iteration", "edge_type", "value")
all_cent_df$relevance <- as.factor(all_cent_df$node %in% c("V1","V2","V3"))

# plot boxplots (one for each feature)
ggplot(all_cent_df[all_cent_df$node != "T",], aes(x = node, y = value, fill = relevance)) + geom_boxplot(outlier.size = 1) + theme_bw() +
  labs(y = "out-degree centrality", x = "feature") + 
  theme(plot.title = element_text(size=10)) +
  coord_cartesian(ylim = c(4, 16)) +
  #scale_x_discrete(labels= c("I", "R")) +
  facet_wrap(~ edge_type, nrow = 1) +
  theme(strip.text = element_text(size=14))

# plot boxplots (grouping relevant and irrelevant features)
ggplot(all_cent_df[all_cent_df$node != "T",], aes(x = relevance, y = value)) + geom_boxplot() + theme_bw() +
  labs(y = "out-degree centrality", x = "feature type", title = "Relevant (R) and Irrelevant (I) features") + 
       theme(plot.title = element_text(size=10)) +
  coord_cartesian(ylim = c(0, 18)) +
  scale_x_discrete(labels= c("I", "R")) +
  facet_wrap(~ edge_type, nrow = 1)


############## 2. EDGE WEIGHT VALIDATION EXPERIMENT ##############

cluster_number = 2
file_name_all_e <- file.path("synthetic", paste("v1_all_edges_cluster_number", cluster_number, ".csv", sep = "_"))
all_edge_df <- read.csv(file = file_name_all_e)
all_edge_matrix <- acast(all_edge_df, node_from ~ node_to ~ iteration ~ edge_type) 
var_names <- paste("V", seq(dim(all_edge_matrix)[1]-1), sep="")

# compute matrix indicating the number of clusters a pair of features can discriminate
n_iter = 30
cluster_centers <- get_cluster_centers(cluster_number)
discr_matrix <- get_discr_matrix(cluster_centers)
colnames(discr_matrix) <- rownames(discr_matrix) <- var_names

# only consider edges between variables
all_edge_df <- all_edge_df[(all_edge_df$node_from!="T") & (all_edge_df$node_to!="T"),]
# scale edges by the sum of the adjecency matrix
all_edge_df <- all_edge_df %>%
  group_by(iteration, edge_type) %>%
  mutate(scaled = value/sum(value)*100) %>% as.data.frame()

# add number of discriminated cluster
all_edge_df$discriminant <- sapply(1:nrow(all_edge_df), 
  function(x) { discr_matrix[all_edge_df[x, "node_from"], all_edge_df[x, "node_to"]] })
all_edge_df$discriminant <- as.factor(all_edge_df$discriminant)

# compute correlation coefficient
edge_types = c("present", "fixation", "level", "sample")
labels <- array(0, length(edge_types))
for (e in 1:length(edge_types)){
  edge_df = all_edge_df[all_edge_df$edge_type==edge_types[e],]
  test = cor.test(as.integer(edge_df$discriminant), edge_df$scaled, method="pearson")
  labels[e] = paste(edge_types[e], "\nPearson's corr.\ncoef. :", round(test$estimate,3))
}

#plot boxplots (grouping by discriminated clusters)
all_edge_df$edge_type <- factor(all_edge_df$edge_type, levels = edge_types, labels = labels)
ggplot(all_edge_df, aes(x = discriminant, y = scaled)) + geom_boxplot(outlier.size = 0.5) + theme_bw() +
  labs(#title = "edge weight of feature pair discriminating power",
       y = "norm. edge weight", x = "discriminated clusters", title = "Discriminating power of feature pairs") + 
  theme(plot.title = element_text(size=10, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 2.2)) +
  facet_wrap(~ edge_type, nrow = 1)
