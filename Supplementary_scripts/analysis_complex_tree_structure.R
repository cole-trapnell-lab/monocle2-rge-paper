rm(list = ls())

#benchmark with complicated tree structure: 
library(SLICER)
library(destiny)
library(dpt)
library(monocle)
library(simplePPT)

##################################################################################################################################################################
main_fig_dir <- "./Figures/main_figures/"
SI_fig_dir <- "./Figures/supplementary_figures/"

##################################################################################################################################################################
#state1 #E669A7 ; state 2: #D09329 state 3: #B27AB4 state 4: #8FAA3C state 5: #6E94CC state 6: #02B8E3 state 7: #02B8E3
color_scale = c('1' = '#F2756D', '2' = '#D29429', '3' = '#93AC3D', '4' = '#29B34A', '5' = '#2DB99C', '6' = '#0BB9E2', '7' = '#6F94CC', '8' = '#B37BB5', '9' = '#E76BA8')
#run new dpt 
run_new_dpt_exprs <- function(exprs, branching = T, normalize = T, root = NULL, color_by = 'Hours'){
  if(normalize)
    data <- t(log(exprs + 1))
  else data <- t(exprs)
  
  dm <- DiffusionMap(data)
  dpt <- DPT(dm)
  
  ts <- dm@transitions
  M <- destiny:::accumulated_transitions(dm)
  
  if(is.null(root)){
    
    plot(dpt, root = 2, paths_to = c(1,3), col_by = 'branch', pch = 20)
    # plot(dpt, col_by = 'branch', divide = 3, dcs = c(-1,3,-2), pch = 20)
  }
  else{
    dm <- DiffusionMap(data)
    dpt <- DPT(dm, tips = 1)
    plot(dpt, root = root, paths_to = c(1,3), col_by = 'branch', pch = 20)
    plot(dpt, col_by = 'branch', divide = 3, dcs = c(-1,3,-2), pch = 20)
  }
  
  # if('Hours' %in% colnames(pData(cds)))
  #   pData(cds)$Hours <- pData(cds)[, color_by]
  p1 <- 'test'#qplot(DM$DC1, DM$DC2, colour = pData(cds)$Hours)
  
  dp_res <- list(dm = dm, pt = dpt, ts = ts, M = M, ev = dm@eigenvectors, p1 = p1)
  
  return(dp_res)
}

source('./scripts/function.R', echo=TRUE)
library(simplePPT)
library(R.matlab)
library(xacHelper)

params.maxIter = 100;

# load data
tree_300 <- R.matlab::readMat('./data/tree_300.mat')
params.lambda = 0.0464; s = 0.05;
X <- tree_300$X
X = t(X);
D = nrow(X); N <- ncol(X)

params.bandwidth = 2 * s * s;

# X <- rbind(X, colSums(X))
#principal tree method
pt_res <- principal_tree(X, MU = NULL, lambda = params.lambda, bandwidth = params.bandwidth, verbose = T);

# plot the learned principal tree
qplot(X[1, ], X[2, ]) + geom_point(aes( x = pt_res$MU[1, ], y = pt_res$MU[2, ]), color = 'red')

#calculate pseudotime: 
graph_mat <- as(as.matrix(pt_res$stree), "matrix")
row.names(graph_mat) <- rep(1:nrow(graph_mat))
colnames(graph_mat) <- rep(1:ncol(graph_mat))
graph_mat[graph_mat > 0] <- 1
dp_mst <- graph_from_adjacency_matrix(graph_mat > 0, weighted = NULL, mode = 'undirected')

dp <- as.matrix(dist(t(pt_res$MU)))
root_cell <- which(pt_res$MU[2, ] == min(pt_res$MU[2, ]))
test_DDRTree_dpt_ordering <- extract_ddrtree_ordering_xj(dp_mst = dp_mst, dp = dp, root_cell)

ddrtree_state = test_DDRTree_dpt_ordering$cell_state
ddrtree_pseudotime = test_DDRTree_dpt_ordering$pseudo_time

pdf(paste(SI_fig_dir, "SI6b_SimpPPT.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(pt_res$MU[1, ], pt_res$MU[2, ], size = ddrtree_pseudotime, color = ddrtree_state) + scale_size(range = c(0.1, 1)) + 
  xlab('Dimension 1') + ylab('Dimension 2') + geom_point(aes(X[1, ], X[2, ]), size = 0.2, alpha = 0.7, color = 'black') + nm_theme() 
dev.off()

pdf(paste(SI_fig_dir, "SI6b_SimpPPT.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(pt_res$MU[1, ], pt_res$MU[2, ], size = ddrtree_pseudotime, color = ddrtree_state) + scale_size(range = c(0.1, 1)) + 
  xlab('Dimension 1') + ylab('Dimension 2') + geom_point(aes(X[1, ], X[2, ]), size = 0.2, alpha = 0.7, color = 'black') + nm_theme() 
dev.off()

pdf(paste(SI_fig_dir, "SI6b_SimpPPT_helper.pdf", sep = ''))
qplot(pt_res$MU[1, ], pt_res$MU[2, ], size = ddrtree_pseudotime, color = ddrtree_state) + scale_size(range = c(0.1, 1)) + 
  xlab('Dimension 1') + ylab('Dimension 2') + geom_point(aes(X[1, ], X[2, ]), size = 0.2, alpha = 0.7, color = 'black') 
dev.off()

#compare with dpt method: 
tree_dpt_res <- run_new_dpt_exprs(X, normalize = F, color_by = 'Time')

dpt_Branch =  as.character(tree_dpt_res$pt@branch[, 1])
apply(X, 1, min)
which(X[2, ] == 0.0186)

DPT = tree_dpt_res$pt$DPT103
pdf(paste(SI_fig_dir, "SI6b_SimpPPT_dpt.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(X[1, ], X[2, ], size = DPT, color = dpt_Branch)  + scale_size(range = c(0.1, 1)) +
  xlab('Dimension 1') + ylab('Dimension 2') + nm_theme() + scale_color_manual(values = color_scale)
dev.off()

pdf(paste(SI_fig_dir, "SI6b_SimpPPT_dpt.1.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(tree_dpt_res$dm$DC1, tree_dpt_res$dm$DC2, size = DPT, color = dpt_Branch) + scale_size(range = c(0.1, 1)) + 
  xlab('Dimension 1') + ylab('Dimension 2') + nm_theme() + scale_color_manual(values = color_scale) 
dev.off()

#pseudotime distribution: 
qplot(DPT)

#compare with DDRTree method: 
DDRTree_res <- DDRTree(X, dimensions = 2, verbose = T, maxIter = 5, param.gamma = 0.01)# #, maxIter = 5, sig, initial_method = PCAma = 1e-2, lambda = 1, ncenter = 3, param.gamma = 10, tol = 1e-2

graph_mat <- as(as.matrix(DDRTree_res$stree), "matrix")
row.names(graph_mat) <- rep(1:nrow(graph_mat))
colnames(graph_mat) <- rep(1:ncol(graph_mat))
graph_mat[graph_mat > 0] <- 1
ddtree_mst <- graph_from_adjacency_matrix(graph_mat > 0, weighted = NULL, mode = 'undirected')

dp <- as.matrix(dist(t(DDRTree_res$Y)))
root_cell <- which(DDRTree_res$Y[2, ] == min(DDRTree_res$Y[2, ]))
ddrtree_ordering <- extract_ddrtree_ordering_xj(dp_mst = dp_mst, dp = dp, "103")
ddrtree_pseudotime <- ddrtree_ordering$pseudo_time
State <- ddrtree_ordering$cell_state

pdf(paste(SI_fig_dir, "SI6b_DDRTree.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(X[1, ], X[2, ], size = ddrtree_pseudotime, color = as.character(ddrtree_state)) + scale_size(range = c(0.1, 1))  + 
  xlab('Dimension 1') + ylab('Dimension 2') + nm_theme() + scale_color_manual(values = color_scale)
dev.off()

qplot(DDRTree_res$Z[1, ], DDRTree_res$Z[2, ], size = ddrtree_pseudotime, color = ddrtree_state) + 
  geom_point(aes(DDRTree_res$Z[1, ], DDRTree_res$Z[2, ]), color = 'gray') + 
  xlab('Dimension 1') + ylab('Dimension 2') + nm_theme() + scale_size(range = c(0.5, 0.5)) 

pdf(paste(SI_fig_dir, "SI6b_DDRTree.1.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(DDRTree_res$Y[1, ], DDRTree_res$Y[2, ], size = ddrtree_pseudotime, color = State) + scale_size(range = c(0.1, 1)) + 
  xlab('Dimension 1') + ylab('Dimension 2') + nm_theme() 
dev.off()

qplot(DPT, ddrtree_pseudotime)

##################################################################################################################################################################
#run slicer and wishbone: 
run_slicer_exprs <- function(exprs, normalize = F, start = NULL, min_branch_len = 10){
  if(normalize)
    traj <- t(log(exprs + 1))
  else traj <- t(exprs)
  
  k = select_k(traj, kmin=5)
  traj_lle = lle(traj[,], m=2, k)$Y
  traj_graph = conn_knn_graph(traj_lle,5)
  ends = find_extreme_cells(traj_graph, traj_lle)
  
  if(is.null(start))
    start <- ends[1]
  
  cells_ordered = cell_order(traj_graph, start)
  branches = assign_branches(traj_graph,start, min_branch_len = min_branch_len)
  
  return(list(traj_lle = traj_lle, ends = ends, order_df = data.frame(cells_ordered = cells_ordered, branches = branches)))
}

tree_slicer_res <- run_slicer_exprs(X, normalize = F, start = 103)

Branch =  as.character(tree_dpt_res$pt@branch[, 1])
apply(X, 1, min)
which(X[2, ] == 0.0186)

slicer_branch = as.character(tree_slicer_res$order_df$branches)
slicer_time = tree_slicer_res$order_df$cells_ordered

pdf(paste(SI_fig_dir, "SI6b_Slicer.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(X[1, ], X[2, ], size = slicer_time, color = slicer_branch)  + scale_size(range = c(0.1, 1)) +
  xlab('Dimension 1') + ylab('Dimension 2') + nm_theme() + scale_color_manual(values = color_scale) 
dev.off()

##################################################################################################################################################################
#add monocle 1: 
PD <- data.frame(Index = rep(1:300, 1), row.names = paste('Cell_', rep(1:300, 1), sep = ''))
FD <- data.frame(gene_short_names = c('Gene_1', 'Gene_2', 'Gene_3'), row.names = c('Gene_1', 'Gene_2', 'Gene_3'))

X_update <- rbind(X, Gene_3 = colSums(X))
dimnames(X_update) <- list(row.names(FD), row.names(PD))
exprs_mat <- as.matrix(X_update)
sim_tree_cds <-  newCellDataSet(X_update, 
                                phenoData = new("AnnotatedDataFrame", data = PD), 
                                featureData = new("AnnotatedDataFrame", data = FD), 
                                expressionFamily = negbinomial(), 
                                lowerDetectionLimit = 0.1)
sim_tree_cds <- estimateSizeFactors(sim_tree_cds)
# sim_tree_cds <- estimateDispersions(sim_tree_cds)
pData(sim_tree_cds)$Size_Factor <- 1

sim_tree_cds <- setOrderingFilter(sim_tree_cds, ordering_genes = row.names(sim_tree_cds))
sim_tree_cds <- reduceDimension(sim_tree_cds, norm_method = 'none', max_components = 2, verbose = T, reduction_method = 'ICA')
sim_tree_cds <- orderCells(sim_tree_cds, num_paths = 5)
plot_cell_trajectory(sim_tree_cds) + facet_wrap(~State)

sim_tree_cds <- orderCells(sim_tree_cds, num_paths = 5, root_state = 5)
monocle1_time <- pData(sim_tree_cds)$Pseudotime
monocle1_branch <-  pData(sim_tree_cds)$State
pdf(paste(SI_fig_dir, "SI6b_monocle.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(X[1, ], X[2, ], size = monocle1_time, color = monocle1_branch) + scale_size(range = c(0.1, 1)) +
  xlab('Dimension 1') + ylab('Dimension 2') + nm_theme() 
dev.off()

####################################################################################################################################################################################
# test on the simulated neuron datasets: 
####################################################################################################################################################################################
load('./RData/neuron_sim_data.RData')

#neuron + astrocyte + oligodendrocyte
#############################################
#simple PPT
#############################################
neuron_pt_res <- principal_tree(exprs(nao_sim_cds), MU = NULL, lambda = 1, bandwidth = 30, verbose = T);

# plot the learned principal tree
qplot(exprs(nao_sim_cds)['Olig2', ], exprs(nao_sim_cds)['Scl', ])

#calculate pseudotime: 
neuron_graph_mat <- as(as.matrix(neuron_pt_res$stree), "matrix")
row.names(neuron_graph_mat) <- rep(1:nrow(neuron_graph_mat))
colnames(neuron_graph_mat) <- rep(1:ncol(neuron_graph_mat))
neuron_graph_mat[neuron_graph_mat > 0] <- 1
neuron_dp_mst <- graph_from_adjacency_matrix(neuron_graph_mat > 0, weighted = NULL, mode = 'undirected')

neuron_dp <- as.matrix(dist(t(neuron_pt_res$MU)))
root_cell <- which(neuron_pt_res$MU[2, ] == min(neuron_pt_res$MU[2, ]))

dimnames(neuron_dp) <- list(c(1:nrow(neuron_dp)), c(1:ncol(neuron_dp)))
neuron_DDRTree_dpt_ordering <- extract_ddrtree_ordering_xj(dp_mst = neuron_dp_mst, dp = neuron_dp, 1083)

neuron_ddrtree_state = neuron_DDRTree_dpt_ordering$cell_state
neuron_ddrtree_pseudotime = neuron_DDRTree_dpt_ordering$pseudo_time

pdf(paste(SI_fig_dir, "SI6a_SimplePPT.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(exprs(nao_sim_cds)['Mash1', ], exprs(nao_sim_cds)['Hes5', ], size = neuron_ddrtree_pseudotime, color = as.character(neuron_ddrtree_state)) + scale_size(range = c(0.1, 1)) + 
  xlab('Mash1') + ylab('Hes5') + nm_theme() + scale_color_manual(values = color_scale)
dev.off()

pdf(paste(SI_fig_dir, "SI6a_DDRTree.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(exprs(nao_sim_cds)['Olig2', ], exprs(nao_sim_cds)['Scl', ], size = neuron_ddrtree_pseudotime, color = as.character(neuron_ddrtree_state)) + scale_size(range = c(0.1, 1)) + 
  xlab('Olig2') + ylab('Scl') + nm_theme()  + scale_color_manual(values = color_scale)
dev.off()

#############################################
#monocle 1 
#############################################
sim_tree_cds <- orderCells(sim_tree_cds, num_paths = 5, root_state = 5)

nao_sim_cds_subset <- reduceDimension(nao_sim_cds[, sort(sample(1:1200, 800))], reduction_method = 'ICA', norm_method = 'none', verbose = T, scaling = F)
nao_sim_cds_subset <- orderCells(nao_sim_cds_subset, num_paths = 3)
plot_cell_trajectory(sim_tree_cds) + facet_wrap(~State)

plot_cell_trajectory(nao_sim_cds_subset, color_by = 'State') + scale_x_reverse() + scale_y_reverse() + nm_theme() #
pdf(paste(SI_fig_dir, "SI6a_monocle.pdf", sep = ''), height = 1.5, width = 1.5)
plot_cell_trajectory(nao_sim_cds_subset, color_by = 'State') + scale_x_reverse() + scale_y_reverse() + nm_theme() + scale_color_manual(values = color_scale)#
dev.off()

#############################################
#slicer 
#############################################
nao_slicer_res <- run_slicer_exprs(exprs(nao_sim_cds), normalize = F, start = 1)
slicer_branch = as.character(nao_slicer_res$order_df$branches)
slicer_time = nao_slicer_res$order_df$cells_ordered

pdf(paste(SI_fig_dir, "SI6a_slicer.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(nao_slicer_res$traj_lle[, 1], nao_slicer_res$traj_lle[, 2], size = slicer_time, color = slicer_branch)  + scale_size(range = c(0.1, 1)) +
  xlab('Dimension 1') + ylab('Dimension 2') + nm_theme() # + scale_color_manual(values = color_scale)
dev.off()

#############################################
#dpt 
#############################################
nao_dpt_res <- run_new_dpt_exprs(exprs(nao_sim_cds), normalize = F, color_by = 'Time')

dpt_Branch =  as.character(nao_dpt_res$pt@branch[, 1])
# apply(X, 1, min)
# which(X[2, ] == 0.0186)

DPT = nao_dpt_res$pt$DPT793
pdf(paste(SI_fig_dir, "SI6a_DPT.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(nao_dpt_res$dm$DC1, nao_dpt_res$dm$DC2, size = DPT, color = dpt_Branch) + scale_size(range = c(0.1, 1)) +
  xlab('Dimension 1') + ylab('Dimension 2') + nm_theme() + scale_color_manual(values = color_scale)
dev.off()

#############################################
#use dpt result for ordering: 
#############################################

neuron_pt_res_dpt_dm <- principal_tree(t(nao_dpt_res$dm[, 1:3]), MU = NULL, lambda = params.lambda, bandwidth = params.bandwidth, verbose = T);


#############################################
#wishbone 
#############################################
write.csv(file = './csv_data/Wishbone_test_data/neuron_simulation.txt', t(exprs(nao_sim_cds)), quote = F, row.names = T)
write.csv(file = './csv_data/Wishbone_test_data/tree_simulation.txt', t(X_update), quote = F, row.names = T)

neuron_sim_wishbone_res <- read.table("/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/neuron_simulation_wishbone_res.txt", header = T, sep= '\t', row.names=NULL)
complex_tree_wishbone_res <- read.table("/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree//wishbone/complex_tree_simulation_wishbone_res.txt", header = T, sep= '\t', row.names=NULL)

pdf(paste(SI_fig_dir, "SI6a_wishbone.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(neuron_sim_wishbone_res$dm1, neuron_sim_wishbone_res$dm2, size = neuron_sim_wishbone_res$trajectory, color = as.character(neuron_sim_wishbone_res$branch))  + scale_size(range = c(0.1, 1)) +
  xlab('Dimension 1') + ylab('Dimension 2') + nm_theme() + scale_color_manual(values = color_scale)
dev.off()

pdf(paste(SI_fig_dir, "SI6a_wishbone.1.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(complex_tree_wishbone_res$dm1, complex_tree_wishbone_res$dm2)  + scale_size(range = c(0.1, 1)) + #, size = complex_tree_wishbone_res$trajectory, color = as.character(complex_tree_wishbone_res$branch)
  xlab('Dimension 1') + ylab('Dimension 2') + nm_theme() 
dev.off()

pdf(paste(SI_fig_dir, "SI6b_wishbone.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(X[1, ], X[2, ], color = I('black'))  + scale_size(range = c(0.1, 1)) + #, size = complex_tree_wishbone_res$trajectory, color = as.character(complex_tree_wishbone_res$branch)
  xlab('Dimension 1') + ylab('Dimension 2') + nm_theme() 
dev.off()

qplot(X[1, ], X[2, ], size = monocle1_time, color = monocle1_branch) + scale_size(range = c(0.1, 1)) +
  xlab('Dimension 1') + ylab('Dimension 2') + nm_theme() 

##################################################################################################################################################################
# L1 graph (not used in the paper)
##################################################################################################################################################################
load("./RData/simulation_l1_res.RData")
filenames = c('Circle','two_moon','tree_300','Spiral_SUN','three_clusters','DistortedSShape')

tree_300 <- create_cds(data = res_list[[5]]$X)
return_myself <- function(data) {
  return(t(data))
}

tree_300 <- reduceDimension(tree_300, norm_method = 'none', pseudo_expr = 0, scaling = F,
                                   reduction_method = 'L1-span-tree', maxiter = maxiter, initial_method = return_myself,
                                   eps = eps, lambda = 1, gamma = 10,
                                   sigma = 0.01, nn = 5, verbose = T)
tree_300 <- orderCells(tree_300)

pdf(paste(SI_fig_dir, "fig_si2i_l1_tree.pdf", sep = ''), height = 1.5, width = 1.5)
plot_cell_trajectory(tree_300, color_by = 'State', show_branch_points = F, cell_size = 0.4) + nm_theme() + scale_color_manual(values = color_scale)#
dev.off()

tree_300_graph <- reduceDimension(tree_300, norm_method = 'none', pseudo_expr = 0, scaling = F,
                            reduction_method = 'L1-graph', maxiter = maxiter, initial_method = return_myself,
                            eps = eps, lambda = 1, gamma = 10,
                            sigma = 0.01, nn = 5, verbose = T)

pdf(paste(SI_fig_dir, "fig_si2i_l1_graph.pdf", sep = ''), height = 1.5, width = 1.5)
plot_cell_trajectory(tree_300_graph, color_by = 'State', show_branch_points = F, cell_size = 0.4) + nm_theme() + scale_color_manual(values = color_scale)#
dev.off()

##################################################################################################################################################################
#save the data: 
save.image('./RData/analysis_complex_tree_structure.RData')
