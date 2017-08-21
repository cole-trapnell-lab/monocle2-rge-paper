# this file is used to combine DPT with Monocle 2 for trajectory reconstruction (not used in the paper)

library(monocle)
library(destiny)
library(mcclust)
library(plyr)

extract_ddrtree_ordering_xj <- function(dp_mst, dp = dp, root_cell, verbose=T) # dp,
{
  nei <- NULL
  type <- NULL
  pseudo_time <- NULL
  
  curr_state <- 1
  
  res <- list(subtree = dp_mst, root = root_cell)
  
  states = rep(1, ncol(dp))
  names(states) <- V(dp_mst)$name
  
  pseudotimes = rep(0, ncol(dp))
  names(pseudotimes) <- V(dp_mst)$name
  
  parents = rep(NA, ncol(dp))
  names(parents) <- V(dp_mst)$name
  
  mst_traversal <- graph.dfs(dp_mst, 
                             root=root_cell, 
                             neimode = "all", 
                             unreachable=FALSE, 
                             father=TRUE)
  mst_traversal$father <- as.numeric(mst_traversal$father)
  curr_state <- 1
  
  state_stat <- table(degree(dp_mst)[degree(dp_mst) > 2])
  state_total_num <- sum(state_stat * 2:(length(state_stat) + 1))
  
  node_num <-  state_total_num + 2
  state_mst <- make_empty_graph(n = node_num)  #number of states
  state_mst <- add_edges(state_mst, c(node_num, 1))
  
  for (i in 1:length(mst_traversal$order)){
    curr_node = mst_traversal$order[i]
    curr_node_name = V(dp_mst)[curr_node]$name
    
    if (is.na(mst_traversal$father[curr_node]) == FALSE){
      parent_node = mst_traversal$father[curr_node]
      parent_node_name = V(dp_mst)[parent_node]$name
      parent_node_pseudotime = pseudotimes[parent_node_name]
      parent_node_state = states[parent_node_name]
      curr_node_pseudotime = parent_node_pseudotime + dp[curr_node_name, parent_node_name]
      if (degree(dp_mst, v=parent_node_name) > 2){
        curr_state <- curr_state + 1
        # if(curr_state >= 1405){
        #   # browser()
        # }
        message('current state is ', curr_state, 'parent state is ', parent_node_state)
        state_mst <- add_edges(state_mst, c(parent_node_state, curr_state))
      }
    }else{
      parent_node = NA
      parent_node_name = NA
      curr_node_pseudotime = 0
    }
    
    curr_node_state = curr_state
    pseudotimes[curr_node_name] <- curr_node_pseudotime
    states[curr_node_name] <- curr_node_state
    parents[curr_node_name] <- parent_node_name
  }
  
  ordering_df <- data.frame(sample_name = names(states),
                            cell_state = factor(states),
                            pseudo_time = as.vector(pseudotimes),
                            parent = parents)
  row.names(ordering_df) <- ordering_df$sample_name
  # ordering_df <- plyr::arrange(ordering_df, pseudo_time)
  E(state_mst)$weight <- c(0.1, table(ordering_df$cell_state))
  V(state_mst)$name <- as.character(1:node_num)
  state_mst <- as.undirected(state_mst)
  return(list(ordering_df = ordering_df, state_mst = state_mst))
  
  state_mst
  
}

##################################################################################################################################################################
# MARS-seq dataset
##################################################################################################################################################################

benchmark_type <- "na_sim_data_wishbone_original_monocle2"
load(paste('./RData/', benchmark_type, '_real_simulation_na.RData', sep = ''))
load('valid_subset_GSE72857_cds2')
# 
# root_cell <- paste('Y_', valid_subset_GSE72857_cds2@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[row.names(subset(pData(valid_subset_GSE72857_cds2), Pseudotime == 0)), 1], sep = '')
# test_res <- extract_ddrtree_ordering_xj(valid_subset_GSE72857_cds2@minSpanningTree, cellPairwiseDistances(valid_subset_GSE72857_cds2),
#                                         root_cell = root_cell, verbose=T)
# plot(test_res$state_mst, layout = layout_as_tree(test_res$state_mst))
# 
# # do the pqtree on this complicated principal tree:
# 
# order_res <- extract_ddrtree_ordering_xj(dp_mst, dp = dp, root_cell = which(degree(dp_mst) == 1)[1])
# 
# #assign the branches as well as the branch time point:
# cds <- cds_downsampled_cells_ordered_update[[36]]
# 
# next_node <<- 0
# 
# dp_mst <- cds@minSpanningTree
# dp <- cellPairwiseDistances(cds)
# res <- pq_helper(dp_mst, use_weights=T, root_node=which(degree(dp_mst) == 1)[1])
# 
next_node <<- 0
dp_mst <- test_res$state_mst
dp_mst <- as.undirected(dp_mst)
E(dp_mst)$weight[11] <- 0.1
dp = distances(dp_mst, v = V(dp_mst), to = V(dp_mst), weights = NULL)


which.min(pData(valid_subset_GSE72857_cds2)$Pseudotime)

res <- monocle:::pq_helper(dp_mst, use_weights=T, root_node=24)
# res <- pq_helper(dp_mst, use_weights=T, root_node=which(degree(dp_mst) == 1)[1])


if(is.null(branch_num))
  branch_num <- sum(degree(dp_mst) > 2) + 1

branch_num <- 7 #6 cell types in the end
order_list <- monocle:::extract_good_branched_ordering(res$subtree, res$root, dp, branch_num, FALSE)

cc_ordering <- order_list$ordering_df
row.names(cc_ordering) <- cc_ordering$sample_name
# 
# test_state_tree <- layout_as_tree(dp_mst); test_state_tree <- as.data.frame(test_state_tree)
# colnames(test_state_tree) <- c('x', 'y')
# test_state_tree$State <- cc_ordering[as.character(1:12), 'cell_state']
# 
# plot(dp_mst, layout = layout_as_tree(dp_mst), vertex.color = cc_ordering[as.character(1:12), 'cell_state'])
# qplot(x, y, color = State, data = test_state_tree)
# 
# qplot(dm$DC1, dm$DC2, color = cc_ordering[row.names(data_ori), 'cell_state'])
# 
# dp <- DPT_res[1:cell_num, 1:cell_num]
# dimnames(dp) <- list(colnames(cds_downsampled_cells_ordered[[36]])[!duplicated(data)], colnames(cds_downsampled_cells_ordered[[36]])[!duplicated(data)])
# gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
# dp_mst <- minimum.spanning.tree(gp)

valid_subset_GSE72857_cds2
data_ori <- as.matrix(t(exprs(valid_subset_GSE72857_cds2)))
data_uniq <- data_ori[!duplicated(data_ori), ]
dm <- DiffusionMap(as.matrix(data_ori))
DPT_res <- DPT(dm)
cell_num <- length(DPT_res$DPT1)

dp <- DPT_res[1:cell_num, 1:cell_num]
dimnames(dp) <- list(colnames(valid_subset_GSE72857_cds2)[!duplicated(data_ori)], colnames(valid_subset_GSE72857_cds2)[!duplicated(data_ori)])
gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
dp_mst <- minimum.spanning.tree(gp)

root_cell <- row.names(subset(pData(valid_subset_GSE72857_cds2), Pseudotime == 0))
test_res <- extract_ddrtree_ordering_xj(dp_mst, dp,
                                        root_cell = root_cell, verbose=T)
plot(test_res$state_mst, layout = layout_as_tree(test_res$state_mst))

qplot(DPT_res@dm$DC1, DPT_res@dm$DC2, color = pData(valid_subset_GSE72857_cds2)$State)

qplot(DPT_res@dm$DC1, DPT_res@dm$DC2, color = cc_ordering[as.character(test_res$ordering_df[colnames(valid_subset_GSE72857_cds2), 'cell_state']), 'cell_state'])

pData(valid_subset_GSE72857_cds)$cell_type2 <- revalue(as.character(pData(valid_subset_GSE72857_cds)$cluster),
  c("1" = 'Ery', "2" = 'Ery', "3" = 'Ery', "4" = 'Ery', "5" = 'Ery', "6" = 'Ery',
  "7" = 'MP/EP', "8" = 'MK', "9" = 'GMP', "10" = 'GMP',
  "11" = 'DC',
  "12" = 'B', "13" = 'B', "14" = 'M', "15" = 'M', "16" = 'N', "17" = 'N', "18" = 'E',
  "19" = 'lymphoid'))

calClusteringMetrics(cc_ordering[as.character(test_res$ordering_df[colnames(valid_subset_GSE72857_cds2), 'cell_state']), 'cell_state'], pData(valid_subset_GSE72857_cds2)$cell_type2)
calClusteringMetrics(pData(valid_subset_GSE72857_cds2)$State, pData(valid_subset_GSE72857_cds2)$cell_type2)

#load the data from the Fabian group's processed results: 
load('./script_for_reproduce/MARSseq_analysis_tutorial.RData')
calClusteringMetrics(pData(valid_subset_GSE72857_cds2)$cluster, cluster.id[cluster.id[, 1] != 19, 1]) #confirm that index matches up 
calClusteringMetrics(branching[cluster.id[, 1] != 19], pData(valid_subset_GSE72857_cds2)$cell_type2) #check the result 

ARI_branches <- data.frame(ARI = c(0.5923145, 0.6566774, 0.7225239), Type = c('DPT (original)', 'DPT + Monocle 2', 'Monocle 2'))
qplot(ARI, geom = 'bar', data = ARI_branches, stat = 'identity')

pdf(paste(SI_fig_dir, 'MARS_seq_ARI_branch_cluster.pdf', sep = ''), width = 2, height = 1)
ggplot(data = ARI_branches, aes(Type, ARI)) + geom_bar(stat = 'identity', aes(fill = Type)) + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab('')
dev.off()

pdf(paste(SI_fig_dir, 'MARS_seq_pseudotime_correspondence.pdf', sep = ''), width = 1, height = 1)
qplot(DPT_res$DPT24, pData(valid_subset_GSE72857_cds2)$Pseudotime, color = pData(valid_subset_GSE72857_cds2)$State, size = I(0.25)) + xlab('DPT diffusion \n pseudotime') + ylab('Monocle 2 \n pseudotime') + nm_theme()
dev.off()

pdf(paste(SI_fig_dir, 'MARS_seq_pseudotime_correspondence_helper.pdf', sep = ''))
qplot(DPT_res$DPT24, pData(valid_subset_GSE72857_cds2)$Pseudotime, color = pData(valid_subset_GSE72857_cds2)$State, size = I(0.25)) + xlab('DPT diffusion pseudotime') + ylab('Monocle 2 pseudotime')
dev.off()

##################################################################################################################################################################
# The full dataset for MARS-seq study
##################################################################################################################################################################
load('/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/RData/fig_SI6.RData')

valid_subset_GSE72857_cds2 <- reduceDimension(valid_subset_GSE72857_cds, norm_method = 'log', verbose = T, max_components = 10, ncenter = 300) 

all_GSE72857_cds <- all_GSE72857_cds[fData(all_GSE72857_cds)$use_for_ordering, ]
plot_pc_variance_explained(all_GSE72857_cds, norm_method = 'log')
all_GSE72857_cds <- reduceDimension(all_GSE72857_cds, norm_method = 'log', verbose = T, max_components = 10, ncenter = 300)
all_GSE72857_cds <- orderCells(all_GSE72857_cds)

data_ori <- as.matrix(t(exprs(all_GSE72857_cds[fData(all_GSE72857_cds)$use_for_ordering, ])))
data_uniq <- data_ori[!duplicated(data_ori), ]
dm <- DiffusionMap(as.matrix(data_ori))
DPT_res <- DPT(dm)
cell_num <- length(DPT_res$DPT1)

dp <- DPT_res[1:cell_num, 1:cell_num]
dimnames(dp) <- list(colnames(all_GSE72857_cds)[!duplicated(data_ori)], colnames(all_GSE72857_cds)[!duplicated(data_ori)])
gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
dp_mst <- minimum.spanning.tree(gp)

root_cell <- row.names(subset(pData(all_GSE72857_cds), Pseudotime == 0))
test_res <- extract_ddrtree_ordering_xj(dp_mst, dp,
                                        root_cell = root_cell, verbose=T)

next_node <<- 0
dp_mst <- test_res$state_mst
dp_mst <- as.undirected(dp_mst)
dp = distances(dp_mst, v = V(dp_mst), to = V(dp_mst), weights = NULL)

which.min(pData(all_GSE72857_cds)$Pseudotime)

res <- monocle:::pq_helper(dp_mst, use_weights=T, root_node=4446)
# res <- pq_helper(dp_mst, use_weights=T, root_node=which(degree(dp_mst) == 1)[1])

if(is.null(branch_num))
  branch_num <- sum(degree(dp_mst) > 2) + 1

branch_num <- 7 #6 cell types in the end
order_list <- monocle:::extract_good_branched_ordering(res$subtree, res$root, dp, branch_num, FALSE)

cc_ordering <- order_list$ordering_df
row.names(cc_ordering) <- cc_ordering$sample_name

plot(test_res$state_mst, layout = layout_as_tree(test_res$state_mst))

qplot(DPT_res@dm$DC1, DPT_res@dm$DC2, color = pData(all_GSE72857_cds)$State)

qplot(DPT_res@dm$DC1, DPT_res@dm$DC2, color = cc_ordering[as.character(test_res$ordering_df[colnames(all_GSE72857_cds), 'cell_state']), 'cell_state'])

cellPairwiseDistances(all_GSE72857_cds) <- DPT_res[1:cell_num, 1:cell_num]
all_GSE72857_cds@minSpanningTree <- dp_mst

all_GSE72857_cds <- diameter_path_ordering(all_GSE72857_cds, num_paths = 7, root_state = Root_state(all_GSE72857_cds)) # 'Cell_2'
plot_cell_trajectory(all_GSE72857_cds)
all_GSE72857_cds <- diameter_path_ordering(all_GSE72857_cds, num_paths = 7, root_state = Root_state(all_GSE72857_cds)) # 'Cell_2'
plot_cell_trajectory(all_GSE72857_cds)

pData(all_GSE72857_cds)$cell_type2 <- revalue(as.character(pData(all_GSE72857_cds)$cluster),
                                               c("1" = 'Ery', "2" = 'Ery', "3" = 'Ery', "4" = 'Ery', "5" = 'Ery', "6" = 'Ery',
                                                 "7" = 'MP/EP', "8" = 'MK', "9" = 'GMP', "10" = 'GMP',
                                                 "11" = 'DC',
                                                 "12" = 'B', "13" = 'B', "14" = 'M', "15" = 'M', "16" = 'N', "17" = 'N', "18" = 'E',
                                                 "19" = 'lymphoid'))

calClusteringMetrics(DPT_res@branch[, 1], pData(all_GSE72857_cds)$cell_type)
calClusteringMetrics(cc_ordering[as.character(test_res$ordering_df[colnames(all_GSE72857_cds), 'cell_state']), 'cell_state'], pData(all_GSE72857_cds)$cell_type)
calClusteringMetrics(pData(all_GSE72857_cds[, colnames(valid_subset_GSE72857_cds2)])$State, pData(all_GSE72857_cds[, colnames(valid_subset_GSE72857_cds2)])$cell_type)

ARI_branches <- data.frame(ARI = c(0.11459862, 0.549, 0.408), Type = c('DPT (original)', 'DPT + Monocle 2', 'Monocle 2'))
qplot(ARI, geom = 'bar', data = ARI_branches, stat = 'identity')

pdf(paste(SI_fig_dir, 'MARSEQ_all_ARI_branch_cluster.pdf', sep = ''), width = 1, height = 1)
ggplot(data = ARI_branches, aes(Type, ARI)) + geom_bar(stat = 'identity', aes(fill = Type)) + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab('')
dev.off()

which(pData(all_GSE72857_cds)$Pseudotime == 0)
pdf(paste(SI_fig_dir, 'MARSEQ_all_pseudotime_correspondence.pdf', sep = ''), width = 1, height = 1)
qplot(DPT_res$DPT77, pData(all_GSE72857_cds)$Pseudotime, color = pData(all_GSE72857_cds)$State, size = I(0.25)) + xlab('DPT diffusion \n pseudotime') + ylab('Monocle 2 \n pseudotime') + nm_theme()
dev.off()

pdf(paste(SI_fig_dir, 'MARSEQ_all_pseudotime_correspondence_helper.pdf', sep = ''))
qplot(DPT_res$DPT77, pData(all_GSE72857_cds)$Pseudotime, color = pData(all_GSE72857_cds)$State, size = I(0.25)) + xlab('DPT diffusion pseudotime') + ylab('Monocle 2 pseudotime')
dev.off()

##################################################################################################################################################################
# Olsson dataset (WT dataset)
##################################################################################################################################################################
# load('./RData/fig5.RData')
load('/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/RData/fig5_4_22.RData')

#do the samething as above but for the Olsson dataset:  
data_ori <- as.matrix(t(exprs(URMM_all_fig1b[fData(URMM_all_fig1b)$use_for_ordering, ])))
data_uniq <- data_ori[!duplicated(data_ori), ]
dm <- DiffusionMap(as.matrix(data_ori))
DPT_res <- DPT(dm)
cell_num <- length(DPT_res$DPT1)

dp <- DPT_res[1:cell_num, 1:cell_num]
dimnames(dp) <- list(colnames(URMM_all_fig1b)[!duplicated(data_ori)], colnames(URMM_all_fig1b)[!duplicated(data_ori)])
gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
dp_mst <- minimum.spanning.tree(gp)

root_cell <- row.names(subset(pData(URMM_all_fig1b), Pseudotime == 0))
test_res <- extract_ddrtree_ordering_xj(dp_mst, dp,
                                        root_cell = root_cell, verbose=T)

next_node <<- 0
dp_mst <- test_res$state_mst
dp_mst <- as.undirected(dp_mst)
dp = distances(dp_mst, v = V(dp_mst), to = V(dp_mst), weights = NULL)

which.min(pData(URMM_all_fig1b)$Pseudotime)

res <- monocle:::pq_helper(dp_mst, use_weights=T, root_node=77)
# res <- pq_helper(dp_mst, use_weights=T, root_node=which(degree(dp_mst) == 1)[1])

if(is.null(branch_num))
  branch_num <- sum(degree(dp_mst) > 2) + 1

branch_num <- 3 #6 cell types in the end
order_list <- monocle:::extract_good_branched_ordering(res$subtree, res$root, dp, branch_num, FALSE)

cc_ordering <- order_list$ordering_df
row.names(cc_ordering) <- cc_ordering$sample_name

plot(test_res$state_mst, layout = layout_as_tree(test_res$state_mst))

qplot(DPT_res@dm$DC1, DPT_res@dm$DC2, color = pData(URMM_all_fig1b)$State)

qplot(DPT_res@dm$DC1, DPT_res@dm$DC2, color = cc_ordering[as.character(test_res$ordering_df[colnames(URMM_all_fig1b), 'cell_state']), 'cell_state'])

pData(URMM_all_fig1b)$cell_type <- revalue(as.character(pData(URMM_all_fig1b)$cluster),
                                                       c("HSCP-1" = 'HSC', "HSCP-2" = 'HSC', "Meg" = 'MEK', "Eryth" = 'MEK', "Multi-Lin" = 'HSC', "MDP" = 'HSC',
                                                         "Mono" = 'Mono', "Gran" = 'Gran', "Myelocyte" = 'Gran'))

calClusteringMetrics(DPT_res@branch[, 1], pData(URMM_all_fig1b)$cell_type)
calClusteringMetrics(cc_ordering[as.character(test_res$ordering_df[colnames(URMM_all_fig1b), 'cell_state']), 'cell_state'], pData(URMM_all_fig1b)$cell_type)
calClusteringMetrics(pData(URMM_all_fig1b)$State, pData(URMM_all_fig1b)$cell_type)

ARI_branches <- data.frame(ARI = c(0.11459862, 0.549, 0.408), Type = c('DPT (original)', 'DPT + Monocle 2', 'Monocle 2'))
qplot(ARI, geom = 'bar', data = ARI_branches, stat = 'identity')

pdf(paste(SI_fig_dir, 'URMM_ARI_branch_cluster.pdf', sep = ''), width = 1, height = 1)
ggplot(data = ARI_branches, aes(Type, ARI)) + geom_bar(stat = 'identity', aes(fill = Type)) + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab('')
dev.off()

which(pData(URMM_all_fig1b)$Pseudotime == 0)
pdf(paste(SI_fig_dir, 'URMM_pseudotime_correspondence.pdf', sep = ''), width = 1, height = 1)
qplot(DPT_res$DPT77, pData(URMM_all_fig1b)$Pseudotime, color = pData(URMM_all_fig1b)$State, size = I(0.25)) + xlab('DPT diffusion \n pseudotime') + ylab('Monocle 2 \n pseudotime') + nm_theme()
dev.off()

pdf(paste(SI_fig_dir, 'URMM_pseudotime_correspondence_helper.pdf', sep = ''))
qplot(DPT_res$DPT77, pData(URMM_all_fig1b)$Pseudotime, color = pData(URMM_all_fig1b)$State, size = I(0.25)) + xlab('DPT diffusion pseudotime') + ylab('Monocle 2 pseudotime')
dev.off()

##################################################################################################################################################################
# Olsson dataset (all dataset)
##################################################################################################################################################################

#do the samething as above but for the Olsson dataset:  
data_ori <- as.matrix(t(exprs(URMM_all_abs[fData(URMM_all_abs)$use_for_ordering, ])))
data_uniq <- data_ori[!duplicated(data_ori), ]
dm <- DiffusionMap(as.matrix(data_ori))
DPT_res <- DPT(dm)
cell_num <- length(DPT_res$DPT1)

dp <- DPT_res[1:cell_num, 1:cell_num]
dimnames(dp) <- list(colnames(URMM_all_abs)[!duplicated(data_ori)], colnames(URMM_all_abs)[!duplicated(data_ori)])
gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
dp_mst <- minimum.spanning.tree(gp)

root_cell <- row.names(subset(pData(URMM_all_abs), Pseudotime == 0))
test_res <- extract_ddrtree_ordering_xj(dp_mst, dp,
                                        root_cell = root_cell, verbose=T)

next_node <<- 0
dp_mst <- test_res$state_mst
dp_mst <- as.undirected(dp_mst)
dp = distances(dp_mst, v = V(dp_mst), to = V(dp_mst), weights = NULL)

which.min(pData(URMM_all_abs)$Pseudotime)

res <- monocle:::pq_helper(dp_mst, use_weights=T, root_node=77)
# res <- pq_helper(dp_mst, use_weights=T, root_node=which(degree(dp_mst) == 1)[1])

if(is.null(branch_num))
  branch_num <- sum(degree(dp_mst) > 2) + 1

branch_num <- 3 #6 cell types in the end
order_list <- monocle:::extract_good_branched_ordering(res$subtree, res$root, dp, branch_num, FALSE)

cc_ordering <- order_list$ordering_df
row.names(cc_ordering) <- cc_ordering$sample_name

plot(test_res$state_mst, layout = layout_as_tree(test_res$state_mst))

qplot(DPT_res@dm$DC1, DPT_res@dm$DC2, color = pData(URMM_all_abs)$State)

qplot(DPT_res@dm$DC1, DPT_res@dm$DC2, color = cc_ordering[as.character(test_res$ordering_df[colnames(URMM_all_abs), 'cell_state']), 'cell_state'])

pData(URMM_all_abs)$cell_type <- revalue(as.character(pData(URMM_all_abs)$paper_cluster),
                                           c("HSCP-1" = 'HSC', "HSCP-2" = 'HSC', "Meg" = 'MEK', "Eryth" = 'MEK', "Multi-Lin" = 'MPP', "MDP" = 'MPP',
                                             "Mono" = 'Mono', "Gran" = 'Gran', "Myelocyte" = 'Gran'))

calClusteringMetrics(DPT_res@branch[colnames(URMM_all_abs) %in% colnames(URMM_all_fig1b), 1], pData(URMM_all_abs[, colnames(URMM_all_fig1b)])$cell_type)
calClusteringMetrics(cc_ordering[as.character(test_res$ordering_df[colnames(URMM_all_fig1b), 'cell_state']), 'cell_state'], pData(URMM_all_abs[, colnames(URMM_all_fig1b)])$cell_type)
calClusteringMetrics(pData(URMM_all_abs[, colnames(URMM_all_fig1b)])$State, pData(URMM_all_abs[, colnames(URMM_all_fig1b)])$cell_type)

ARI_branches <- data.frame(ARI = c(0.1976251, 0.4797827, 0.4311858), Type = c('DPT (original)', 'DPT + Monocle 2', 'Monocle 2'))
qplot(ARI, geom = 'bar', data = ARI_branches, stat = 'identity')

pdf(paste(SI_fig_dir, 'URMM_abs_ARI_branch_cluster.pdf', sep = ''), width = 1, height = 1)
ggplot(data = ARI_branches, aes(Type, ARI)) + geom_bar(stat = 'identity', aes(fill = Type)) + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab('')
dev.off()

which(pData(URMM_all_abs)$Pseudotime == 0)
pdf(paste(SI_fig_dir, 'URMM_abs_pseudotime_correspondence.pdf', sep = ''), width = 1, height = 1)
qplot(DPT_res$DPT77, pData(URMM_all_abs)$Pseudotime, color = pData(URMM_all_abs)$State, size = I(0.25)) + xlab('DPT diffusion \n pseudotime') + ylab('Monocle 2 \n pseudotime') + nm_theme()
dev.off()

pdf(paste(SI_fig_dir, 'URMM_abs_pseudotime_correspondence_helper.pdf', sep = ''))
qplot(DPT_res$DPT77, pData(URMM_all_abs)$Pseudotime, color = pData(URMM_all_abs)$cell_type, size = I(0.25)) + xlab('DPT diffusion pseudotime') + ylab('Monocle 2 pseudotime')
dev.off()

save.image('./RData/revision_1_dpt_multiple_branches.RData')

