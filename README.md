# Feature Graphs for Interpretable Unsupervised Tree Ensembles: Centrality, Interaction, and Application in Disease Subtyping

This repository contains the code and data for the article: **Feature graphs for interpretable unsupervised tree ensembles: centrality, interaction, and application in disease subtyping** by Christel Sirocchi, Martin Urschler, and Bastian Pfeifer.

## Overview

The repository includes the implementation of the proposed graph-building and graph-mining methods, as well as the code to replicate all experiments presented in the paper.

The proposed method constructs feature graphs from unsupervised random forests (uRF), which are trained using a Fixation Index (FST) split rule, as introduced in the article *Federated unsupervised random forest for privacy-preserving patient stratification* by Pfeifer et al. (2024). [DOI: https://doi.org/10.1093/bioinformatics/btae382](https://doi.org/10.1093/bioinformatics/btae382)

In this project, we compile the source code as a modified `ranger` package to access all parameters and configurations, rather than using the uRF library directly. The code will be fully integrated into the uRF library in future work.

## Installation

This project requires two versions of the `ranger` package, installed in different locations:

1. **Original `ranger` package**:  
   Install from CRAN by specifying a library path:  
   ```r
   install.packages("ranger", lib = "RF_PATH")
2. **Modified `ranger` package**:  
   Compile from the source file `uRF_source.tar.gz` found in this repository:  
   ```r
   install.packages("uRF_source.tar.gz", repos = NULL, type = "source", lib = "uRF_PATH")

## Repository Structure

- **`clustering-data-v1/`**: Contains benchmark datasets used in the experiments. These datasets can be retrieved from [this GitHub repository](https://github.com/gagolews/clustering-data-v1) or downloaded directly from UCI.

- **`synthetic/`**: Contains files generated for experiments on synthetic datasets.

- **`benchmark/`**: Contains files generated for experiments on real-world datasets.

- **`validation_1.r`**: Evaluates feature centrality and edge weight in synthetic datasets with relevant and irrelevant features (results in Section 5.1.1).

- **`validation_1_params.r`**: Hyperparameter study, evaluating the effect of the number of trees and the number of features selected at each node on feature centrality and edge weight in synthetic datasets (results in Section 5.1.1).

- **`validation_2.r`**: Evaluates cluster-specific graphs on synthetic datasets, grouping features into cluster-specific, sub-relevant, and irrelevant (results in Section 5.1.2).

- **`validation_3.r`**: Evaluates the proposed greedy and brute-force algorithms on synthetic datasets with varying numbers of relevant features and computes separability scores (results in Section 5.1.3).

- **`validation_4.r`**: Evaluates feature selection with redundant features in synthetic datasets (results in Section 5.1.4).

- **`validation_5.r`**: Performs feature selection using the proposed approaches and three state-of-the-art methods over 10 datasets.

- **`validation_6.r`**: Computes the ARI and NMI scores for clustering solutions obtained on the top k features selected by each method, and conducts statistical analysis on these scores (results in Section 5.2).

- **`utilities.r`**: Contains functions for implementing the proposed graph-building and graph-mining methods, as well as functions for computing metrics (monotonicity, separability) and plotting feature graphs.
   
