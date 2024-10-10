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
n_cols <- dim(cluster_centers)[2]

# generate synthetic data around predefined cluster centers and train uRF
mtrys = c(1,2,round(sqrt(n_cols)), n_cols)
ntrees = c(10,50,100,500,1000)

# compute edge matrix for each dataset and edge type
edge_types = c("present", "fixation", "level", "sample")
all_edge_matrix <- array(0, dim = c(n_feat+1, n_feat+1, length(mtrys), length(ntrees), n_iter, length(edge_types)))

for (t in 1:length(ntrees)){
  print(ntrees[t])
  for (m in 1:length(mtrys)){
    print(mtrys[m])
    for (n in 1:n_iter) {
      print(n)
      train <- get_synthetic_dataset(cluster_centers, n_samples_per_cluster)
      target <- rep(1:k, each = n_samples_per_cluster)
      #ggplot(train, aes(x = V2, y = V3, color = as.factor(target))) + geom_point() 
      model = ranger(x=train, num.trees=ntrees[t], 
                     clustering = TRUE,  probability= FALSE, oob.error = FALSE,
                     mtry = mtrys[m], min.bucket = 5)
      for (e in 1:length(edge_types)){ 
          all_edge_matrix[,,m,t,n,e] <- compute_edge_matrix(model, train, edge_types[e])
      }
    }
  }
}

# save results
var_names <- model$forest$independent.variable.names
dimnames(all_edge_matrix) <- list(c(var_names,"T"), c(var_names,"T"), mtrys, ntrees, paste("iter", seq(n_iter), sep=""), edge_types)
edge_df <- melt(all_edge_matrix)
colnames(edge_df) <- c("node_from", "node_to", "mtrys", "ntrees","iteration", "edge_type", "value")
file_name_all_e <- file.path("synthetic", paste("v1_hyperparams_all_edges_cluster_number", cluster_number, ".csv", sep = "_"))
write.csv(edge_df, file = file_name_all_e, row.names = FALSE)

stop("finished computation")


##############################################################
############# READ DATA & ANALYSE FEATURE GRAPHS #############
##############################################################

cluster_centers <- get_cluster_centers(cluster_number)
k <- dim(cluster_centers)[1]
n_feat <- dim(cluster_centers)[2]
n_samples_per_cluster <- 50  
n_cols <- dim(cluster_centers)[2]

mtrys = c(1,2,round(sqrt(n_cols)), n_cols)
ntrees = c(10,50,100,500,1000)


############## 1. CENTRALITY VALIDATION EXPERIMENT w HYPERPARAMETERS ##############

cluster_number = 1
file_name_all_e <- file.path("synthetic", paste("v1_hyperparams_all_edges_cluster_number", cluster_number, ".csv", sep = "_"))
edge_df <- read.csv(file = file_name_all_e)
all_edge_matrix <- acast(edge_df, node_from ~ node_to ~ mtrys ~ ntrees ~ iteration ~ edge_type) 
var_names <- paste("V", seq(dim(all_edge_matrix)[1]-1), sep="")
all_edge_matrix <- all_edge_matrix[c(var_names,"T"), c(var_names,"T"), , , , ] # reorder variables

# compute out-degree centralities (row-sums in the adjacency matrix) for each dataset and edge type
n_iter <- 30
edge_types = c("present", "fixation", "level", "sample")
all_out_degree_centralities <- array(0, dim = c(length(c(var_names,"T")), length(mtrys), length(ntrees), n_iter, length(edge_types)))

for (t in 1:length(ntrees)){
  for (m in 1:length(mtrys)){
    for (e in 1:length(edge_types)){
      for (n in 1:n_iter) {
        edge_matrix <- all_edge_matrix[,,m,t,n,e]
        # compute normalised out-degree centrality
        all_out_degree_centralities[,m,t,n,e] <- rowSums(edge_matrix)/ sum(rowSums(edge_matrix)) *100
      }
    }
  }
}
dimnames(all_out_degree_centralities) <- list(c(var_names,"T"), paste(mtrys, "mtry", sep=" "), paste(ntrees, "trees", sep=" "), paste("iter",seq(n_iter),sep=""), edge_types)
all_cent_df <- melt(all_out_degree_centralities)
colnames(all_cent_df) <- c("node", "mtrys", "ntrees", "iteration", "edge_type", "value")
all_cent_df$relevance <- as.factor(all_cent_df$node %in% c("V1","V2","V3"))
levels(all_cent_df$mtrys) <- c("MTRY = 1", "MTRY = 2", "MTRY = sqrt", "MTRY = n")
levels(all_cent_df$relevance) <- c("irrelevant features", "relevant features")

# Perform t-test for each combination of mtrys, ntrees, and edge_type
t_test_results <- all_cent_df[all_cent_df$node != "T",] %>%
  group_by(mtrys, ntrees, edge_type) %>%
  summarise(
    statistic = t.test(value[relevance == "relevant features"], value[relevance == "irrelevant features"])$statistic,
    p_value = t.test(value[relevance == "relevant features"], value[relevance == "irrelevant features"])$p.value
  )
t_test_results["log_p_value"] <- -log10(t_test_results$p_value)


############## 2. EDGE WEIGHT VALIDATION EXPERIMENT w HYPERPARAMETERS ##############

cluster_number = 2
file_name_all_e <- file.path("synthetic", paste("v1_hyperparams_all_edges_cluster_number", cluster_number, ".csv", sep = "_"))
all_edge_df <- read.csv(file = file_name_all_e)
all_edge_matrix <- acast(all_edge_df, node_from ~ node_to ~ mtrys ~ ntrees ~ iteration ~ edge_type) 
var_names <- paste("V", seq(dim(all_edge_matrix)[1]-1), sep="")
n_iter <- 30
edge_types = c("present", "fixation", "level", "sample")

# compute matrix indicating the number of clusters a pair of features can discriminate
cluster_centers <- get_cluster_centers(cluster_number)
discr_matrix <- get_discr_matrix(cluster_centers)
colnames(discr_matrix) <- rownames(discr_matrix) <- var_names

# only consider edges between variables
all_edge_df <- all_edge_df[(all_edge_df$node_from!="T") & (all_edge_df$node_to!="T"),]
# scale edges by the sum of the adjecency matrix
all_edge_df <- all_edge_df %>%
  group_by(iteration, edge_type, mtrys, ntrees) %>%
  mutate(scaled = value/sum(value)*100) %>% as.data.frame()

# add number of discriminated cluster
all_edge_df$discriminant <- sapply(1:nrow(all_edge_df), 
  function(x) { discr_matrix[all_edge_df[x, "node_from"], all_edge_df[x, "node_to"]] })
all_edge_df$discriminant <- as.factor(all_edge_df$discriminant)
all_edge_df$mtrys <- as.factor(all_edge_df$mtrys)
all_edge_df$ntrees <- as.factor(all_edge_df$ntrees)
levels(all_edge_df$mtrys) <- c("MTRY = 1", "MTRY = 2", "MTRY = sqrt", "MTRY = n")
levels(all_edge_df$ntrees) <- c("10 trees", "50 trees", "100 trees", "500 trees", "1000 trees")

# Perform pearson correlation for each combination of mtrys, ntrees, and edge_type
p_test_results <- all_edge_df %>%
  group_by(mtrys, ntrees, edge_type) %>%
  summarise(
    statistic = cor(as.integer(discriminant), scaled, method = "pearson"),
    p_value = cor.test(as.integer(discriminant), scaled, method = "pearson")$p.value
  )

#### combine results from both evaluations ####

t_test_results$study = "t-test\nof node centralities"
p_test_results$study = "Pearson's correlation\nof edge weights"
t_test_results$statistic = t_test_results$log_p_value

all_tests = rbind(t_test_results, p_test_results)

levels(all_tests$ntrees) <- ntrees
all_tests$study <- factor(all_tests$study, c("t-test\nof node centralities","Pearson's correlation\nof edge weights"))
all_tests$edge_type <- factor(all_tests$edge_type, edge_types)
ggplot(all_tests, aes(x = factor(ntrees), y = statistic, color = edge_type, group = edge_type)) + 
  geom_line() + geom_point() + theme_bw() +
  labs( y = "Pearson's corr. coef.           -log10 (p-value)", x = "number of trees") + 
  theme(
    plot.title = element_text(size = 10), 
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "bottom", 
    legend.title = element_blank()
  ) +
  geom_hline(data = all_tests %>% filter(study == "t-test\nof node centralities"),
             aes(yintercept = -log10(0.05)), col = "blue", linetype = 2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  facet_grid(study ~ mtrys,  scales = "free_y")
