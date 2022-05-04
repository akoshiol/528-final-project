## Author: Adalee Koshiol
## Class: BF528
## Assignment: Final Project
## Start Date: 4/26/22

##################
# file description: for programmer and analyst roles from project 1
##################

############## PROGRAMMER ROLE ##############
# if you haven't installed the biocmanager, then do so
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

# install all the necessary packages for the project:
#     affy, affyPLM, sva, AnnotationDbi, and hgu133plus2.db

BiocManager::install(c("affy","affyPLM","sva","AnnotationDbi","hgu133plus2.db"))

# confirm they are installed and load the pacakages
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(ggplot2)
library(tidyverse)


## read all the CEL files
# read in the CEL files (with the symbolic link) with the ReadAffy function
affy1 <- ReadAffy(celfile.path = '/projectnb/bf528/users/frazzled/project_1/samples/CEL_files/')
# read in the CEL file (within our actual project folder) with the ReadAffy function
affy2 <- ReadAffy(celfile.path = '/projectnb/bf528/users/frazzled/project_1/samples/')
# merge the two affys together
affy <- merge(affy1, affy2)

## use the rma function to normalize
normalized_affy <- affy::rma(affy)

## convert normalized affy into PLM set for RLE and NUSE
fit <- fitPLM(affy, normalize=TRUE, background=TRUE)
RLE <- RLE(fit, type='stats')
NUSE <- NUSE(fit, type='stats')

## make histograms for the deliverables using ggplot
# first rows of NUSE and RLE are the medians
df <- tibble(RLE_median = RLE[1,], NUSE_median = NUSE[1,])
# make the RLE histogram
RLE_histogram <- ggplot(df, aes(x=RLE_median)) +
  geom_histogram(bins = 50) +
  xlab("Median") +
  ylab("Count") +
  ggtitle("Counts of the RLE Medians")
# make the NUSE histogram
NUSE_histogram <- ggplot(df, aes(x=NUSE_median)) +
  geom_histogram(bins = 50) +
  xlab("Median") +
  ylab("Count") +
  ggtitle("Counts of the NUSE Medians")

## read in the provided annotation and get out the necessary variables
annotation <- readr::read_csv('/project/bf528/project_1/doc/proj_metadata.csv')
# get the batch effects
batch_effects <- annotation$normalizationcombatbatch
# get the features of interest
features <- model.matrix(~as.factor(normalizationcombatmod), data = annotation)

#transforming normalized data into a matrix
normalized_mat <- exprs(normalized_affy)

#correcting batch effects while preserving features of interest
combat_data <- ComBat(normalized_mat, batch_effects, features)

# perform PCA on normalized data
t_mat <- t(normalized_mat)
s_mat <- scale(t_mat)
pca_mat <- t(s_mat)
pca_results <- prcomp(pca_mat, scale = FALSE, center = FALSE)

# grab the first two principal components
pc1 <- pca_results$rotation[,1]
pc2 <- pca_results$rotation[,2]
pc_df <- tibble(pc1, pc2)

# grab percent variability of the PCs
v <- pca_results$sdev^2
ve <- round((v / sum(v)) * 100, digits = 2)

# plot the principal components
plot <- ggplot(pc_df, aes(x=pc1, y=pc2)) +
  geom_point() +
  xlab(paste('PC1 (Variance Explained: ', as.character(ve[1]), '%)')) +
  ylab(paste('PC2 (Variance Explained: ', as.character(ve[2]), '%)')) +
  ggtitle('PC1 vs. PC2')


############## ANALYST ROLE ##############

filter1 <- function(expression_data){
  #make a copy without the probeids
  copy <- dplyr::select(expression_data, -probeids)
  ## filter the data so that 20% of the gene values are above log2(15)
  # find the number of samples that is 20% of all samples
  num_genes <- .20 * ncol(copy)
  # grab the threshold of log2(15)
  threshold <- log2(15)
  # find which ones in the expression matrix have values above the threshold
  above_threshold <- as_tibble(sapply(copy, function(x){ifelse(x>threshold, 1, 0)}))
  # count the number of samples that are above the threshold
  num_above_threshold <- as_tibble(apply(above_threshold, 1, sum))
  # merge the count of those above the threshold with the original expression data
  expr_data_threshold <- dplyr::mutate(expression_data, threshold = num_above_threshold)
  # filter the merged tibble for those with .20 of values above the theshold
  expr_data_filtered <- dplyr::filter(expr_data_threshold, threshold > num_genes) %>%
    dplyr::select(-threshold)
  return(expr_data_filtered)
}

filter2 <- function(expression_data){
  #make a copy without the probeids
  copy <- dplyr::select(expression_data, -probeids)
  ## have variance significantly different from the median variance
  # find the degrees of freedom
  degrees <- ncol(copy) - 1
  # calculate the chi-squared
  chi_sq = qchisq((1 - 0.01)/2, degrees, lower.tail = FALSE)
  # find the variances of each sample
  variance <- apply(copy, 1, var)
  # compute the test statistic for each gene
  test_stat <- (degrees*variance/median(variance))
  # filter the data for variance
  expr_data_filtered <- dplyr::mutate(expression_data, test_statistic = test_stat) %>% 
    dplyr::filter(test_statistic > chi_sq) %>% 
    dplyr::select(-test_statistic)
  return(expr_data_filtered)
}

filter3 <- function(expression_data){
  #make a copy without the probeids
  copy <- dplyr::select(expression_data, -probeids)
  ## have a coefficient of variation above 0.186
  # calculate the coefficient of variant for each gene (the row)
  cv <- as_tibble(apply(copy, 1, function(x){sd(x)/mean(x)}))
  # filter the data for the coefficient of variation
  expr_data_filtered <- dplyr::mutate(expression_data, cv = cv) %>%
    dplyr::filter(cv > 0.186) %>%
    dplyr::select(-cv)
  return(expr_data_filtered)
}

# expression data comes from the combat data matrix
expr_data <- as_tibble(combat_data)
expr_data$probeids <- rownames(combat_data)

# apply the first filter
expr_fil1 <- filter1(expr_data)
# apply the second filter
expr_fil2 <- filter2(expr_fil1)
# apply the third filter
expr_fil3 <- filter3(expr_fil2)

## perform hierarchical clustering on filtered data matrix
# transpose the data
transposed <- as_tibble(cbind(nms = names(expr_fil3), t(expr_fil3)))
columns <- c("sample", expr_data$probeids)
colnames(transposed) <- columns
expr_t <- transposed[1:nrow(transposed)-1,]

# make a distance matrix
dist_mat <- dist(expr_t, method = "euclidean")
# produce the dendrogram
hclust <- hclust(dist_mat, method = 'average')
# cut the tree into two clusters
cut <- cutree(hclust, k = 2)

## create a heatmap of the gene-expression
# prep the data for the heatmap
copy_expr_t <- dplyr::select(expr_t, -sample)
expr_hm <- as.matrix(sapply(copy_expr_t, as.numeric))
# make colors depending upon if the sample belongs to C3
colors <- transmute(annotation, color = if_else(SixSubtypesClassification == 'C3', 'red', 'blue'))
t_colors <- t(colors)
# make the heatmap
hm <- heatmap(expr_hm)

## cut the data into the clusters
# label which samples are for which cluster
cut_labeled <- copy_expr_t
cut_labeled$cluster <- cut
# put the samples in their respective clusters
cluster1 <- dplyr::filter(cut_labeled, cluster == "1") %>%
  dplyr::select(-cluster) 
cluster1 <- as_tibble(sapply(cluster1, as.numeric))
cluster2 <- dplyr::filter(cut_labeled, cluster == "2") %>%
  dplyr::select(-cluster)
cluster2 <- as-tibble(sapply(cluster2, as.numeric))

## do the t test
# get the number of columns in the clusters
num_cols <- ncol(cluster1)
# make an empty datafram to put the p-values in 
empty <- data.frame(matrix(ncol=num_cols, nrow = 1))
colnames(empty) <- colnames(cluster1)
# iterate through the clusters
for (x in 1:num_cols){
  test <- t.test(cluster1[,x], cluster2[,x])
  empty[1,x] <- test$p.value
}

