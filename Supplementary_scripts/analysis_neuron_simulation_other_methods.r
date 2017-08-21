rm(list = ls())

#1. prove the HSMM branching tree is correct; 
#2. theoretical discussion of the pseudotime  vs real time for different methods (simple and complex trajectory)
#2.1 consider adding branch information into the algorithm 

library(monocle)
library(reshape2)
library(xacHelper)
library(grid)
library(destiny)
library(SLICER)

# load('./RData/fig3.RData')
load('./RData/fig3.RData')

################################################################################################################################################################################################################################################
## Monocle 2
################################################################################################################################################################################################################################################
neuron_dist_mat <- as.matrix(dist(t(exprs(neuron_sim_cds))))
neuron_dist_manifold_dist <- cumsum(c(0, diag(neuron_dist_mat[-1, ])))

low_d_neuron_dist_mat <- as.matrix(dist(t(neuron_sim_cds@reducedDimS))) #reducedDimK
low_d_neuron_dist_manifold_dist <- cumsum(c(0, diag(low_d_neuron_dist_mat[-1, ])))

neuron_sim_cds <- orderCells(neuron_sim_cds, reverse = T)

#do this for the two branchpoints simulation: 
nao_sim_cds <- orderCells(nao_sim_cds, root_state = 4)

plot_cell_trajectory(nao_sim_cds, color_by = 'Pseudotime')
nao_dist_mat <- as.matrix(dist(t(exprs(nao_sim_cds))))
nao_dist_manifold_dist <- cumsum(c(0, diag(nao_dist_mat[-1, ])))

low_d_nao_dist_mat <- as.matrix(dist(t(nao_sim_cds@reducedDimS))) #reducedDimK
low_d_nao_dist_manifold_dist <- cumsum(c(0, diag(low_d_nao_dist_mat[-1, ])))

################################################################################################################################################################################################################################################
## Monocle 1
################################################################################################################################################################################################################################################
# neuron_sim_cds_monocle1 <- reduceDimension(neuron_sim_cds, reduction_method = 'ICA', norm_method = 'none')
# neuron_sim_cds_monocle1 <- orderCells(neuron_sim_cds_monocle1, num_paths = 1)
# plot_cell_trajectory(neuron_sim_cds_monocle1, color_by = 'Pseudotime')

# nao_sim_cds_monocle1 <- reduceDimension(nao_sim_cds, reduction_method = 'ICA', norm_method = 'none')
# nao_sim_cds_monocle1 <- orderCells(nao_sim_cds_monocle1, num_paths = 1)
# plot_cell_trajectory(nao_sim_cds_monocle1, color_by = 'Pseudotime')


################################################################################################################################################################################################################################################
## dpt
################################################################################################################################################################################################################################################
run_new_dpt_exprs <- function (cds, branching = T, root = NULL) 
{
    data <-  t(as.matrix(exprs(cds)[fData(cds)$use_for_ordering, ]))
    dm <- DiffusionMap(data)
    dpt <- DPT(dm, branching = branching, tips = 1)
    ts <- dm@transitions
    M <- destiny:::accumulated_transitions(dm)
    if (is.null(root)) {
        plot(dpt, root = 2, paths_to = c(1, 3), col_by = "branch", 
            pch = 20)
    }
    else {
        dm <- DiffusionMap(data)
        dpt <- DPT(dm, tips = 1)
        plot(dpt, root = root, paths_to = c(1, 3), col_by = "branch", 
            pch = 20)
    }
    dp_res <- list(dm = dm, pt = dpt, ts = ts, M = M, ev = dm@eigenvectors)
    return(dp_res)
}

# avoid the error: stats$g >= gmin is not TRUE
DPT <- function (dm, tips = random_root(dm), ..., w_width = 0.1)
{
  if (!is(dm, "DiffusionMap"))
    stop("dm needs to be of class DiffusionMap, not ", class(dm))
  if (!length(tips) %in% 1:3)
    stop("you need to specify 1-3 tips, got ", length(tips))
  dpt <- destiny:::dummy_dpt(dm)
  all_cells <- seq_len(nrow(dpt))
  stats <- destiny:::tipstats(dpt, all_cells, tips)
  branches <- auto_branch(dpt, all_cells, stats, w_width)
  colnames(branches$branch) <- paste0("Branch", seq_len(ncol(branches$branch)))
  colnames(branches$tips) <- paste0("Tips", seq_len(ncol(branches$tips)))
  dpt@branch <- branches$branch
  dpt@tips <- branches$tips
  dpt
}

auto_branch <- function (dpt, cells, stats, w_width, nmin = 10L, gmin = 1.1)
{
  n <- length(cells)
  stopifnot(n >= nmin)
  if(stats$g < gmin)
    stats$g <- 1.1
  # stopifnot(stats$g >= gmin)
  branches <- destiny:::cut_branches(dpt[cells, stats$tips], cells, w_width)
  branch <- matrix(destiny:::idx_list_to_vec(branches, cells, n), n,
                   1L)
  tips <- matrix(logical(n), n, 1L)
  tips[match(stats$tips, cells), 1L] <- TRUE
  subs <- mapply(function(idx_sub, i) {
    if (length(idx_sub) < nmin || !i %in% idx_sub)
      return(NULL)
    sub_stats <- destiny:::tipstats(dpt, idx_sub, i)
    if (sub_stats$g < gmin)
      return(NULL)
    destiny:::auto_branch(dpt, idx_sub, sub_stats, w_width, nmin, gmin)
  }, branches, stats$tips, SIMPLIFY = FALSE)
  nonnull_subs <- vapply(subs, Negate(is.null), logical(1L))
  if (any(nonnull_subs)) {
    n_sublevels <- do.call(max, lapply(subs[nonnull_subs],
                                       function(s) ncol(s$branch)))
    branch <- cbind(branch, matrix(NA_integer_, n, n_sublevels))
    tips <- cbind(tips, matrix(NA, n, n_sublevels))
    for (s in which(nonnull_subs)) {
      sub <- subs[[s]]
      idx_sub <- branches[[s]]
      idx_newcol <- seq.int(ncol(branch) - n_sublevels +
                              1L, length.out = ncol(sub$branch))
      stopifnot(ncol(sub$branch) == ncol(sub$tips))
      branch_offset <- max(branch, na.rm = TRUE)
      branch[match(idx_sub, cells), idx_newcol] <- sub$branch +
        branch_offset
      tips[match(idx_sub, cells), idx_newcol] <- sub$tips
    }
  }
  stopifnot(ncol(branch) == ncol(tips))
  list(branch = branch, tips = tips)
}

nao_dpt_res <- run_new_dpt_exprs(nao_sim_cds)
neuron_dpt_res <- run_new_dpt_exprs(neuron_sim_cds, branching = F) #stat$g < 1.1; mannually change that 

################################################################################################################################################################################################################################################
## slicer
################################################################################################################################################################################################################################################
run_slicer_exprs <- function (cds, select_genes = F, start = NULL, min_branch_len = 10) 
{
    traj <- t(exprs(cds))
    if (select_genes) 
        genes = select_genes(traj)
    else genes = 1:ncol(traj)
    k = select_k(traj[, genes], kmin = 5)
    traj_lle = lle(traj[, genes], m = 2, k)$Y
    traj_graph = conn_knn_graph(traj_lle, 5)
    ends = find_extreme_cells(traj_graph, traj_lle)
    if (is.null(start)) 
        start <- ends[1]
    cells_ordered = cell_order(traj_graph, start)
    branches = assign_branches(traj_graph, start, min_branch_len = min_branch_len)
    return(list(traj_lle = traj_lle, ends = ends, order_df = data.frame(cells_ordered = cells_ordered, 
        branches = branches)))
}
        
neuron_slicer_res <- run_slicer_exprs(neuron_sim_cds, select_genes = F)
nao_slicer_res <- run_slicer_exprs(nao_sim_cds, select_genes = F)

################################################################################################################################################################################################################################################
# compare all the five different softwares
################################################################################################################################################################################################################################################

nao_dpt_tmp <- nao_dpt_res$pt$DPT1
nao_dpt <- c(nao_dpt_tmp[1:400], 0, nao_dpt_tmp[401:799], 0, nao_dpt_tmp[800:1198])

## raw pseudotime 
all_time_df <- data.frame("real time" = c(1:ncol(neuron_sim_cds)) * 0.05, 
                          'monocle 2' = pData(neuron_sim_cds)$Pseudotime, 
                          'monocle 1' = pData(neuron_sim_cds)$Pseudotime, 
                          'High dimension' = neuron_dist_manifold_dist, 
                          'Reduced dimension' = low_d_neuron_dist_manifold_dist, 
                          dpt = neuron_dpt_res$pt$DPT1, 
                          slicer = neuron_slicer_res$order_df$cells_ordered)
#do this for all the two branchpoints data: 
nao_all_time_df <- data.frame("real time" = rep(1:400, 3) * 0.05, 
                              'monocle 2' = pData(nao_sim_cds)$Pseudotime, 
                              'monocle 1' = pData(nao_sim_cds)$Pseudotime, 
                              'High dimension' = nao_dist_manifold_dist, 
                              'Reduced dimension' = low_d_nao_dist_manifold_dist, 
                              dpt = nao_dpt, 
                              slicer = nao_slicer_res$order_df$cells_ordered)
mlt_all_time_df <- melt(all_time_df, id.var = 'real.time')
colnames(mlt_all_time_df)[2] <- 'Type'
colnames(mlt_all_time_df)[3] <- 'Time.metric'

pdf(paste(SI_fig_dir, "neuron_comparing_pseudotime.pdf", sep = ''), height = 1, width = 6)
qplot(real.time, Time.metric, data = mlt_all_time_df, color = Type, size = 0.5)  + geom_abline()  + 
  scale_size(range=c(0.5)) + facet_wrap(~Type, scale = 'free', nrow = 1) + 
  monocle:::monocle_theme_opts() + nm_theme()
dev.off()

mlt_nao_all_time_df <- melt(nao_all_time_df, id.var = 'real.time')
colnames(mlt_nao_all_time_df)[2] <- 'Type'
colnames(mlt_nao_all_time_df)[3] <- 'Time.metric'

pdf(paste(SI_fig_dir, "nao_comparing_pseudotime.pdf", sep = ''), height = 1, width = 6)
qplot(real.time, Time.metric, data = mlt_nao_all_time_df, color = Type, size = 0.5)  + geom_abline()  + 
  scale_size(range=c(0.5)) + facet_wrap(~Type, scale = 'free', nrow = 1) + 
  monocle:::monocle_theme_opts() + nm_theme()
dev.off()

## pseudotime ordering 
order_all_time_df <- data.frame('real time' = c(1:ncol(neuron_sim_cds)), 
                          'monocle 2' = order(pData(neuron_sim_cds)$Pseudotime), 
                          'monocle 1' = order(pData(neuron_sim_cds)$Pseudotime), 
                          'High dimension' = order(neuron_dist_manifold_dist), 
                          'Reduced dimension' = order(low_d_neuron_dist_manifold_dist), 
                          dpt = neuron_dpt_res$pt$DPT1, 
                          slicer = neuron_slicer_res$order_df$cells_ordered)
#do this for all the two branchpoints data: 
order_nao_all_time_df <- data.frame('real time' = rep(1:400, 3), 
                                    'monocle 2' = order(pData(neuron_sim_cds)$Pseudotime), 
                                    'monocle 1' = order(pData(neuron_sim_cds)$Pseudotime), 
                                    'High dimension' = order(nao_dist_manifold_dist), 
                                    'Reduced dimension' = order(low_d_nao_dist_manifold_dist), 
                              dpt = order(nao_dpt), 
                              slicer = nao_slicer_res$order_df$cells_ordered)
order_mlt_all_time_df <- melt(order_all_time_df, id.var = 'real.time')
colnames(order_mlt_all_time_df)[2] <- 'Type'
colnames(order_mlt_all_time_df)[3] <- 'Time.metric'

pdf(paste(SI_fig_dir, "ordering_neuron_comparing_pseudotime.pdf", sep = ''), height = 1, width = 6)
qplot(real.time, Time.metric, data = order_mlt_all_time_df, color = Type, size = 0.5)  + geom_abline()  + 
  scale_size(range=c(0.5)) + facet_wrap(~Type, scale = 'free', nrow = 1) + 
  monocle:::monocle_theme_opts() + nm_theme()
dev.off()

order_mlt_nao_all_time_df <- melt(order_nao_all_time_df, id.var = 'real.time')
colnames(order_mlt_nao_all_time_df)[2] <- 'Type'
colnames(order_mlt_nao_all_time_df)[3] <- 'Time.metric'

pdf(paste(SI_fig_dir, "ordering_nao_comparing_pseudotime.pdf", sep = ''), height = 1, width = 6)
qplot(real.time, Time.metric, data = order_mlt_nao_all_time_df, color = Type, size = 0.5)  + geom_abline()  + 
  scale_size(range=c(0.5)) + facet_wrap(~Type, scale = 'free', nrow = 1) + 
  monocle:::monocle_theme_opts() + nm_theme()
dev.off()

#################################################################################################################################################################
# run analysis_complex_tree_structure
#################################################################################################################################################################
source('./analysis_complex_tree_structure.R', echo = T)

#################################################################################################################################################################
# add l1-graph and l1-tree: 
color_scale = c('1' = '#F2756D', '2' = '#D29429', '3' = '#93AC3D', '4' = '#29B34A', '5' = '#2DB99C', '6' = '#0BB9E2', '7' = '#6F94CC', '8' = '#B37BB5', '9' = '#E76BA8')
return_myself <- function(data) {
  return(t(data))
}

maxiter = 20
eps = 1e-5
gstruct = "l1-graph"
lambda = 1.0
gamma = 0.5
sigma = 0.01
nn = 5

l1_graph_nao_sim_cds <- reduceDimension(nao_sim_cds, reduction_method = 'L1-graph',  
                                        max_components = nrow(nao_sim_cds) - 1, initial_method = return_myself, 
                                        scaling = F, norm_method = 'none', #C0 = PCA_HSMM_myo_l1_span_tree@reducedDimK, maxiter = 2,
                                         eps = eps, lambda = lambda, gamma = gamma,
                                         sigma = sigma, #lambda = 0.01, gamma = 0.01, sigma = 100, nn = 5, 
                                         verbose = T)

l1_tree_nao_sim_cds <- orderCells(l1_tree_nao_sim_cds)

l1_tree_nao_sim_cds <- reduceDimension(nao_sim_cds, reduction_method = 'L1-span-tree', 
                                       max_components = nrow(nao_sim_cds) - 1, initial_method = return_myself, 
                                       scaling = F, norm_method = 'none', #C0 = PCA_HSMM_myo_l1_span_tree@reducedDimK, maxiter = 2,
                                       # eps = eps, lambda = lambda, gamma = gamma,
                                       # sigma = sigma, #lambda = 0.01, gamma = 0.01, sigma = 100, nn = 5, 
                                       verbose = T)

l1_tree_nao_sim_cds <- orderCells(l1_tree_nao_sim_cds)

pdf(paste(SI_fig_dir, "l1_graph_fig_si3a_1.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(exprs(nao_sim_cds)['Mash1', ], exprs(nao_sim_cds)['Hes5', ], size = neuron_ddrtree_pseudotime, color = as.character(neuron_ddrtree_state)) + scale_size(range = c(0.1, 1)) + 
  xlab('Mash1') + ylab('Hes5') + nm_theme() + scale_color_manual(values = color_scale)
dev.off()

pdf(paste(SI_fig_dir, "l1_tree_fig_si3a_1.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(exprs(l1_tree_nao_sim_cds)['Olig2', ], exprs(l1_tree_nao_sim_cds)['Scl', ], size = pData(l1_tree_nao_sim_cds)$Pseudotime, color = as.character(pData(l1_tree_nao_sim_cds)$State)) + scale_size(range = c(0.1, 1)) + 
  xlab('Olig2') + ylab('Scl') + nm_theme()  #+ scale_color_manual(values = color_scale)
dev.off()

#################################################################################################################################################################
#save the data 
#################################################################################################################################################################
save.image('./RData/fig_SI3.RData')
