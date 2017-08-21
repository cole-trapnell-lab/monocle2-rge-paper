rm(list = ls())

library(densityClust)
library(monocle)
#variance of gene expression increases near the branch point: 

load('./RData/fig_SI3.RData') #note the simulation results 

var(as.vector(exprs(HSMM_myo[1, ])))

# > nao_sim_cds@auxOrderingData[['DDRTree']]$branch_points
# [1] "Cell_430" "Cell_836"

reducedDimS(nao_sim_cds) <- reducedDimK(nao_sim_cds)

adjusted_K <- Matrix::t(reducedDimK(nao_sim_cds))
dp <- as.matrix(dist(adjusted_K))
cellPairwiseDistances(nao_sim_cds) <- dp
gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
dp_mst <- minimum.spanning.tree(gp)
minSpanningTree(nao_sim_cds) <- dp_mst
nao_sim_cds@dim_reduce_type <- "DDRTree"
nao_sim_cds <- findNearestPointOnMST(nao_sim_cds)
nao_sim_cds <- orderCells(nao_sim_cds)
nao_sim_cds <- orderCells(nao_sim_cds, root_state = 4)
plot_cell_trajectory(nao_sim_cds, color_by = 'Pseudotime')

# > nao_sim_cds@auxOrderingData[['DDRTree']]$branch_points
# [1] "Cell_430" "Cell_836"
pData(nao_sim_cds[, c('Cell_430', 'Cell_836')])
range(pData(nao_sim_cds)$Pseudotime)

#test on the real data: 
nao_sim_cds_subset <- nao_sim_cds[, 801:1200]
cluster <- rep(1:10, each = 20)
cor_df <- lapply(unique(cluster), function(x) apply(as.matrix(exprs(nao_sim_cds_subset[, cluster == x])), 1, function(y) cor(y[1:(length(y) - 1)], y[-1])))
cor_mat <- do.call(cbind.data.frame, cor_df)
colnames(cor_mat) <- paste('cluster', unique(cluster), sep = '_')

cor_mat_mlt <- melt(cor_mat)
qplot(variable, value, data = cor_mat_mlt, geom = 'boxplot', color = variable) #show the variance
qplot(value, data = cor_mat_mlt, geom = 'density', color = variable, facets = '~variable') #show the variance

#separate into three different paths: 
#1. 4 -> 3
#2. 4 -> 2 -> 5
#3. 4 -> 2 -> 1 

nao_sim_cds_subset <- nao_sim_cds[, pData(nao_sim_cds)$State %in% c(1, 2, 4)]
cell_order <- order(pData(nao_sim_cds_subset)$Pseudotime)

cluster <- rep(1, ncol(nao_sim_cds_subset)) 
cluster[cell_order] <- rep(1:10, each = 46)

# high variance
var_df <- lapply(unique(cluster), function(x) apply(as.matrix(exprs(nao_sim_cds_subset[, cluster == x])), 1, var))
var_mat <- do.call(cbind.data.frame, var_df)
colnames(var_mat) <- paste('cluster', unique(cluster), sep = '_')

var_mat_mlt <- melt(var_mat)
qplot(variable, value + 1, data = var_mat_mlt, geom = 'boxplot', color = variable, log = 'y') #show the variance

# high correlation
cor_df <- lapply(unique(cluster), function(x) apply(as.matrix(exprs(nao_sim_cds_subset[, cluster == x])), 1, function(y) cor(y[1:(length(y) - 1)], y[-1])))
cor_mat <- do.call(cbind.data.frame, cor_df)
colnames(cor_mat) <- paste('cluster', unique(cluster), sep = '_')

cor_mat_mlt <- melt(cor_mat)
qplot(variable, value, data = cor_mat_mlt, geom = 'boxplot', color = variable) #show the variance

pData(nao_sim_cds_subset)$cluster <- cluster
plot_cell_trajectory(nao_sim_cds_subset, color_by = 'Pseudotime') + facet_wrap(~cluster)
#Cluster pseudotime into bins and then calculate the variance at each time point for the ordering genes (cells near branch time point have high variance)
densityCluster <- function(data, rho=2, delta=2) {
  Dist <- dist(data)
  Clust <- densityClust(Dist, gaussian=TRUE)
  plot(Clust) # Inspect clustering attributes to define thresholds
  
  Clust <- findClusters(Clust, rho=rho, delta=rho)
  # split(iris[,5], irisClust$clusters)
  Clust
}

# df <- data.frame(x = HSMM_myo@reducedDimS[1, ], y = HSMM_myo@reducedDimS[2, ], cluster = as.character(res$clusters))
# qplot(x, y, color = cluster, data = df, facets = "~cluster") 

pData(HSMM_myo)$trajectory_cluster <- as.character(res$clusters)
plot_cell_trajectory(HSMM_myo, color_by = 'trajectory_cluster') + facet_wrap(~trajectory_cluster)
res <- densityCluster(t(HSMM_myo@reducedDimS))

var(as.vector(exprs(HSMM_myo[1, ])))

var_df <- lapply(unique(res$clusters), function(x) apply(as.matrix(exprs(HSMM_myo[fData(HSMM_myo)$use_for_ordering, res$clusters == x])), 1, var))
var_mat <- do.call(cbind.data.frame, var_df)
colnames(var_mat) <- paste('cluster', unique(res$clusters), sep = '_')

var_mat_mlt <- melt(var_mat)
qplot(variable, value + 1, data = var_mat_mlt, geom = 'boxplot', color = variable, log = 'y')

# #implement the algorithm from DNB: 
# #1. DEG changes at each pseudotime interval: 
# #2. how to determine the number of clusters: 
# 
# #function for performing enrichment analysis on this group: 
# 
# #function to find the genes at certain point: 
# extractDBN <- function(cds = AT12_cds_subset_all_gene, pseudotime_interval = 10,  Max_index_all = test$Index) { 
#   max_ind <- which(Max_index_all == max(Max_index_all), arr.ind = T)  
#   
#   grp_num = ncol(Max_index_all) - 1
#   cds_exprs <- exprs(cds)
#   
#   sample_time <- pData(cds)$Pseudotime
#   
#   rng <- range(pData(cds)$Pseudotime)
#   time_stamp <- seq(rng[1], rng[2], by = pseudotime_interval)
#   
#   valid_cells <- pData(cds)$Pseudotime < time_stamp[max_ind[1]] & pData(cds)$Pseudotime >= time_stamp[max_ind[1] - 1]
#   cds_exprs <- cds_exprs[, valid_cells]
#   valid_genes <- which(apply(cds_exprs, 1, sd) != 0)
#   cds_exprs <- cds_exprs[valid_genes, ]
#   
#   correlation <- cor(t(cds_exprs), method = 'pearson')
#   fit <- hclust(as.dist(correlation), method = 'ward.D2')
#   plot(fit)
#   #     grp_num <- as.integer(readline(prompt="Confirm number of clusters (integer): "))
#   rect.hclust(fit, k=grp_num, border="red")
#   
#   grp <- cutree(fit, k = grp_num)
#   
#   return(list(max_ind = max_ind, all_grp_assignment = grp, DNB = names(grp[grp == max_ind[2]])))
# }
# 
# # #grp_num: maximal number of clusters should have
# # #update the function to output all the datasets generate during the calculations 
# findDBN <- function(cds = AT12_cds_subset_all_gene, pseudotime_interval = 10, grp_num = 4, 
#                     grps_vec = rep(4, 5), use_DEG = F, q_threshold = 0.1) {
#   sample_time <- pData(cds)$Pseudotime
#   
#   rng <- range(pData(cds)$Pseudotime)
#   time_stamp <- seq(rng[1], rng[2], by = pseudotime_interval)
#   
#   index_vec <- matrix(rep(0, grp_num * length(time_stamp)), nrow = length(time_stamp))
#   consec_diff <- array(rep(0, length(time_stamp) * 3), dim = c(length(time_stamp), 3)) #columns: max value, grp, all grps
#   consec_diff <- as.list(numeric(length(time_stamp) * 3)) #use the trick of the list ? 
#   dim(consec_diff) <- c(length(time_stamp), 3)
#   
#   #print(consec_diff)
#   grp_df <- matrix(rep(0, nrow(cds) * (length(time_stamp) - 2)), nrow = length(time_stamp) - 2)
#   colnames(grp_df) <- row.names(cds) 
#   
#   #set the pseudotime group for each cell: 
#   pData(cds)$Pseudotime_rgn = ceiling(pData(cds)$Pseudotime / pseudotime_interval)
#   
#   #calcualte the composite index: 
#   
#   if(length(time_stamp) <= 2) stop('Please storten the pseudotime range') 
#   
#   #length - 2 is equal to number of intervals
#   for(interval_index in 1:(length(time_stamp) - 2)){ 
#     cds_exprs <- exprs(cds)
#     
#     #select cells in the range of pseudotime: 
#     if(interval_index < length(time_stamp) - 2)
#       valid_cells <- pData(cds)$Pseudotime < time_stamp[interval_index + 1] & pData(cds)$Pseudotime >= time_stamp[interval_index]
#     else 
#       valid_cells <- pData(cds)$Pseudotime >= time_stamp[interval_index] 
#     
#     message(paste('The number of cells in the time step is', sum(valid_cells)))
#     #1. just use a group of genes (2. calculate the genes changes between groups): 
#     
#     if(use_DEG) {
#       DEG_res <- differentialGeneTest(cds, fullModelFormulaStr = "~Pseudotime_rgn", cores = round(detectCores() / 2))
#       DEG_genes <- row.names(DEG_res[DEG_res$qval < q_threshold, ])
#       cds_exprs <- cds_exprs[DEG_genes, ]
#     }
#     
#     #cluster genes at each sampling time points by correlation: 
#     cds_exprs <- cds_exprs[, valid_cells]
#     valid_genes <- which(apply(cds_exprs, 1, sd) != 0)
#     cds_exprs <- cds_exprs[valid_genes, ]
#     
#     correlation <- cor(t(cds_exprs), method = 'pearson')
#     fit <- hclust(as.dist(correlation), method = 'ward.D2')
#     plot(fit)
#     #     grp_num <- as.integer(readline(prompt="Number of clusters (integer): "))
#     
#     rect.hclust(fit, k=grps_vec[interval_index], border="red")
#     
#     #     grp_num <- as.integer(readline(prompt="Confirm number of clusters (integer): "))
#     rect.hclust(fit, k=grps_vec[interval_index], border="red")
#     grp <- cutree(fit, k = grps_vec[interval_index])
#     grp_df[interval_index, names(grp)] <- grp
#   }
#   
#   #save(cds, grp, grp_num, time_stamp, interval_index, file = 'DNB_bug.RData')
#   
#   #calculate the sum of difference between two consecutive time points: 
#   
#   #3. Iterate for all the time points and all the genes sets:     
#   
#   #calculate the composite index for the cluster which shows the largest consecutive changes between time points: 
#   #at each interval, there are different number of clusters, calculate the composite index for particular clustering across all time points
#   per_interval_grp_index <- apply(grp_df, 1, function(x) { 
#     gene_id <- which(x > 0) 
#     grp <- x[x > 0]
#     names(grp) <- colnames(grp_df)[gene_id]
#     grp_num <- max(x)
#     
#     grp_index <- lapply(1:(length(time_stamp) - 2), function(y) 
#       calCompIndex(cds, grp, grp_num, time_stamp, interval_index = y))
#     tmp <- do.call(rbind.data.frame, grp_index)
#     colnames(tmp) <- paste('g', 1:ncol(tmp), sep = '')
#     tmp
#   })
#   
#   #per_interval_grp_index is a list of length 4 (number of intervals) and each element contains a vector with dataframe with num_interval * cluster_num
#   max_index_diff_per_interval_grp <-lapply(per_interval_grp_index, function(x) { 
#     index_diff_gene <- apply(x, 2, function(y){
#       index_diff <- 2* y[2:(length(y) - 1)] - y[1:(length(y) - 2)]  -  y[3:length(y)]
#       which(index_diff == max(index_diff), arr.ind = T)
#       max(index_diff)  #largest composite index for each clusters at each interval
#     })
#     max(index_diff_gene) #largest composite index for all clusters at each interval
#   }
#   )
#   
#   interval_grp_max_index_diff <- which(unlist(max_index_diff_per_interval_grp) == max(unlist(max_index_diff_per_interval_grp)))
#   Max_index_all <- per_interval_grp_index[[interval_grp_max_index_diff]] #find the interval where we have the highest values
#   #   save(Max_index_all, file = 'Max_index_all')
#   colnames(Max_index_all) <- paste('g', 1:ncol(Max_index_all), sep = '')
#   Max_index_all$Time <- 1:nrow(Max_index_all)
#   mlt_Max_index_all <- melt(Max_index_all, id.vars = 'Time')
#   
#   #make the plot of the composite index along the pseudotime: 
#   p1 <- qplot(Time, value, color = variable, data = mlt_Max_index_all, geom = c('point', 'line'))
#   
#   return(list(Index = Max_index_all, max_index_diff_per_interval_grp = max_index_diff_per_interval_grp, per_interval_grp_index = per_interval_grp_index, cluster_df = grp_df, plot = p1))
# }
# 
# # #calcualte the composite index for clusters at interval_index: 
# # #2. Determine the dominant group or the DNB by significance analysis: calculate all the three components for the model
# # #calculate the three parts of the values: 
# # #SD: is the average SD of the dominant group: 
# # #|PCC|: average PCC of the dominant group in absolute value 
# # #|PCC_o|: average PCC between the dominant group and others in absolute value 
# calCompIndex <- function(cds, grp, grp_num, time_stamp, interval_index){
#   if(interval_index < length(time_stamp))
#     valid_cells <- pData(cds)$Pseudotime < time_stamp[interval_index + 1] & pData(cds)$Pseudotime > time_stamp[interval_index]
#   else 
#     valid_cells <- pData(cds)$Pseudotime >= time_stamp[interval_index] & pData(cds)$Pseudotime > time_stamp[interval_index]
#   
#   cds_exprs <- exprs(cds)[, valid_cells]
#   valid_genes <- which(apply(cds_exprs, 1, sd) != 0) 
#   valid_genes <- intersect(names(valid_genes), names(grp)) #inconsistent for the gene clusters across time points? 
#   cds_exprs <- cds_exprs[valid_genes, ] 
#   grp <- grp[valid_genes]
#   
#   Index <- rep(0, grp_num)
#   for(grp_ind in 1:grp_num) {
#     SD <- mean(apply(cds_exprs[names(grp[grp == grp_ind]), ], 1, sd))
#     PCC <- mean(abs(cor(t(cds_exprs[names(grp[grp == grp_ind]), ]), method = 'pearson')))
#     PCC_O <- mean(abs(cor(t(cds_exprs))[names(grp[grp == grp_ind]), names(grp[grp != grp_ind])]))
#     
#     Index[grp_ind] <- SD * PCC / PCC_O
#   }
#   
#   Index[is.na(Index)] <- 0
#   return(Index)
# }