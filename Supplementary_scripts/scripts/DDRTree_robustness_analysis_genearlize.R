# ############################################################################################################################################################
# #1. boxplot for the kendall's tau

# #2. consistency of the branch assignment

# #3. boxplot for distribution of pseudotime across different time points (we don't have this...)

# #4. smoothness and lag-1 autocorrelation:

# #5. branch assignment:

# #6. Correlation matrices for genes selected by SLICER. Blue indicates negative correlation and red indicates positive correlation. (a) Genes selected from the distal lung epithelium dataset. (b) Genes selected from the neural stem cell dataset.

# #7. effects on dpt pseudotime calculation from uneven sampling of cells from different branch

# #downsampling cells and test the consistency of branch point, branches and other assignments
# ############################################################################################################################################################
# # Downsampling the number of cells
# ############################################################################################################################################################
# library(devtools)
# load_all('~/Dropbox (Personal)/Projects/monocle-dev')
library(monocle)
library(xacHelper)
library(plyr)
library(stringr)
library(dplyr)
library(grid)
library(gridExtra)
library(mcclust)
library(flexclust)
library(igraph)
library(destiny)
library(dpt)

# source('~/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/scripts/function.R', echo=TRUE)

if(!exists('norm_method_data')) 
  norm_method_data <- 'log'
if(benchmark_type %in% c('na_sim_data', 'na_sim_data_dpt_dim"'))
  norm_method_data <- 'none'
if(!exists("num_dim"))
  num_dim <- 2

############################################################################################################################################################
#using monocle2 to order cells and match the states to the original state from the full cds
order_shalek_cells_by_original_states <- function(cds_subset, root_state, cells_state_2, cells_state_3, 
                                                  auto_param = auto_param_selection, norm_method = norm_method_data, 
                                                  max_components = num_dim, scaling = !(benchmark_type %in% c("na_sim_data_dpt_dim", 'na_sim_data'))) {
  
  if(!(benchmark_type == c('na_sim_data', 'URMM_all_fig1b', 'na_sim_data_dpt_dim')))
    run_new_dpt_function <- NULL
  else
    run_new_dpt_function <- function(data){
      dm <- DiffusionMap(data)
      return(dm@eigenvectors)
    }
    
  cds_subset = reduceDimension(cds_subset, norm_method = norm_method, auto_param_selection = auto_param, max_components = max_components, scaling = scaling, initial_method = run_new_dpt_function, maxIter = 20)
  cds_subset = orderCells(cds_subset)

  if(length(unique(pData(cds_subset)$State)) > 3)
    cds_subset <- trimTree(cds_subset)

  #determine the mapping between original state 1/2/3 and new state 1/2/3:
  overlap_state_1 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 1, ]), root_state))
  overlap_state_2 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 2, ]), root_state))
  overlap_state_3 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 3, ]), root_state))

  #find the state corresponding to the original root state
  overlap_vec <- c(overlap_state_1, overlap_state_2, overlap_state_3)
  max_ind <- which(overlap_vec == max(overlap_vec))

  cds_subset = orderCells(cds_subset, root_state = max_ind)
  if(length(unique(pData(cds_subset)$State)) > 3)
    cds_subset <- trimTree(cds_subset)

  overlap_state_2 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 2, ]), cells_state_2))
  overlap_state_3 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 3, ]), cells_state_2))

  #find the new state corresponding to the original state 2:
  overlap_state_2_vec <- c(overlap_state_2, overlap_state_3)
  state_2_max_ind <- which(overlap_state_2_vec == max(overlap_state_2_vec)) #root state

  if(state_2_max_ind != 1) {
    State_vec <- pData(cds_subset)$State
    pData(cds_subset)$State[State_vec == 2] <- 3
    pData(cds_subset)$State[State_vec == 3] <- 2
  }

  cds_subset
}

#using raw DDRTree algorithm to order cells and match the states to the original state from the full cds
order_shalek_cells_by_original_states_custom_ordering_function <- function(cds_subset, root_state, cells_state_2, cells_state_3, initial_method = PCA, norm_method = norm_method_data, max_components = num_dim) {
  #don't perform the projection from the reduced dimension space to the graph:
  cds_norm_exprs <- convert2DRData(cds_subset, norm_method = norm_method)
  # DDRTree_res <- DDRTree(cds_norm_exprs, dimensions = 2, verbose = F, initial_method = initial_method) #, maxIter = 5, sigma = 1e-2, lambda = 1, ncenter = 3, param.gamma = 10, tol = 1e-2
  DDRTree_res <- DDRTree(as.matrix(log2(exprs(cds_subset[fData(cds_subset)$use_for_ordering, ]) + 1)), dimensions = 2, verbose = F, initial_method = initial_method, max_components = max_components) #, maxIter = 5, sigma = 1e-2, lambda = 1, ncenter = 3, param.gamma = 10, tol = 1e-2
  dpt_ordering <- custom_ordering(DDRTree_res, branch_num = 2)

  #set the name for the ordering dataframe:
  cc_ordering <- dpt_ordering$cc_ordering
  row.names(cc_ordering) <- colnames(cds_norm_exprs)[as.numeric(as.character(cc_ordering$sample_name))]

  #determine the mapping between original state 1/2/3 and new state 1/2/3:
  overlap_state_1 <- length(intersect(row.names(cc_ordering[cc_ordering$cell_state == 1, ]), root_state))
  overlap_state_2 <- length(intersect(row.names(cc_ordering[cc_ordering$cell_state == 2, ]), root_state))
  overlap_state_3 <- length(intersect(row.names(cc_ordering[cc_ordering$cell_state == 3, ]), root_state))

  #find the state corresponding to the original root state
  overlap_vec <- c(overlap_state_1, overlap_state_2, overlap_state_3)
  max_ind <- which(overlap_vec == max(overlap_vec))[1]

  dp <- as.matrix(dist(t(DDRTree_res$Y)))
  gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)
  root_cell <- intersect(which(degree(dp_mst) == 1), cc_ordering[cc_ordering$cell_state == max_ind, "sample_name"]) #find the tip cell of the root state cells
  dpt_ordering <- custom_ordering(DDRTree_res, root_cell = root_cell)

  #set the name for the ordering dataframe:
  cc_ordering <- dpt_ordering$cc_ordering
  row.names(cc_ordering) <- colnames(cds_norm_exprs)[as.numeric(as.character(cc_ordering$sample_name))]

  overlap_state_2 <- length(intersect(row.names(cc_ordering[cc_ordering$cell_state == 2, ]), cells_state_2))
  overlap_state_3 <- length(intersect(row.names(cc_ordering[cc_ordering$cell_state == 3, ]), cells_state_2))

  #find the new state corresponding to the original state 2:
  overlap_state_2_vec <- c(overlap_state_2, overlap_state_3)
  state_2_max_ind <- which(overlap_state_2_vec == max(overlap_state_2_vec)) #root state

  if(state_2_max_ind != 1) {
    State_vec <- cc_ordering$cell_state
    cc_ordering$cell_state[State_vec == 2] <- 3
    cc_ordering$cell_state[State_vec == 3] <- 2
  }

  cc_ordering
}

# #functions to run ica ordering; need to run separately because it is extremely slow or require downsampling whne we deal with big datasets
ICA_order_shalek_cells_by_original_states <- function(cds_subset, root_state, cells_state_2, cells_state_3, norm_method = norm_method_data) {
  cds_subset = reduceDimension(cds_subset, norm_method = norm_method, reduction_method = 'ICA')
  cds_subset = orderCells(cds_subset, num_paths = 2, reverse = T)

  #determine the mapping between original state 1/2/3 and new state 1/2/3:
  overlap_state_1 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 1, ]), root_state))
  overlap_state_2 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 2, ]), root_state))
  overlap_state_3 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 3, ]), root_state))

  #find the state corresponding to the original root state
  overlap_vec <- c(overlap_state_1, overlap_state_2, overlap_state_3)
  max_ind <- which(overlap_vec == max(overlap_vec))

  cds_subset = orderCells(cds_subset, max_ind, num_paths = 2, reverse = T)
  overlap_state_2 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 2, ]), cells_state_2))
  overlap_state_3 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 3, ]), cells_state_2))

  #find the new state corresponding to the original state 2:
  overlap_state_2_vec <- c(overlap_state_2, overlap_state_3)
  state_2_max_ind <- which(overlap_state_2_vec == max(overlap_state_2_vec)) #root state

  if(state_2_max_ind != 1) {
    State_vec <- pData(cds_subset)$State
    pData(cds_subset)$State[State_vec == 2] <- 3
    pData(cds_subset)$State[State_vec == 3] <- 2
  }

  cds_subset
}

plot_monocle_spanning_tree_vectorized = function(cds) {
  monocle::plot_spanning_tree(cds, color_by="Time", cell_size=2) + nm_theme()
}

if(!exists("auto_param_selection"))
  auto_param_selection <- F

# Get the states assignment used in the manuscript
if(!exists("root_state"))
  root_state <- row.names(subset(pData(absolute_cds), State == 1))
if(!exists("AT1_state"))
  AT1_state <- row.names(subset(pData(absolute_cds), State == 2))
if(!exists("AT2_state"))
  AT2_state <- row.names(subset(pData(absolute_cds), State == 3))

############################################################################################################################################################
# Generate downsampled sets of cells
set.seed(20161130)
MIN_PROPORTION = 0.1
MAX_PROPORTION = 1
STEP = 0.1
REPS_PER = 3
EXTRA_PROPORTIONS = c(0.85, 0.95)

if(!exists("repeat_downsampling_num"))
  repeat_downsampling_num <- 100

downsampled_proportions = sort(rep(c(seq(MIN_PROPORTION, MAX_PROPORTION, by=STEP), EXTRA_PROPORTIONS) , REPS_PER))
names(downsampled_proportions) = downsampled_proportions  # will tie CDS objects to proportion for later

set.seed(20161130) #ensure reproducibility

#downsample cells from 0.1 to 1
cds_downsampled_cells_ordered = lapply(downsampled_proportions, function(proportion) {
  message('proportion is ', proportion)
  subset_cds <- absolute_cds[, sort(sample(ncol(absolute_cds), round(ncol(absolute_cds) * proportion)))]  #ensure we get the same values every time we run it (shuffle index will change the kmean values) 

  tryCatch({
    order_shalek_cells_by_original_states(subset_cds, root_state = root_state, cells_state_2 = AT1_state, cells_state_3 = AT2_state) 
  },
  error = function(e) {print(e); subset_cds})
})

#repeated downsample of cells at 0.8
cds_downsampled_cells_ordered_0.8 = lapply(1:repeat_downsampling_num, function(id) {
  message('trial number is ', id)
  subset_cds <- absolute_cds[, sort(sample(ncol(absolute_cds), round(ncol(absolute_cds) * 0.8)))] #ensure we get the same values every time we run it (shuffle index will change the kmean values) 

  tryCatch({
    order_shalek_cells_by_original_states(subset_cds, root_state = root_state, cells_state_2 = AT1_state, cells_state_3 = AT2_state)
  },
  error = function(e) {print(e); NA})
})

save.image(paste('./RData/', benchmark_type, '_DDRTree_robustness_analysis_tmp.RData', sep = ''))

# ################################################################################################################################################
# #use DDRTree only:
# ################################################################################################################################################
# DDRTree_cds_downsampled_cells_ordered = lapply(cds_downsampled_cells_ordered, function(cds) {
#   message('proportion is ', ncol(cds) / ncol(absolute_cds)) #183 full dataset
#   # root_state <- row.names(subset(pData(cds), State == 2))
#   # AT1_state <- row.names(subset(pData(cds), State == 1))
#   # AT2_state <- row.names(subset(pData(cds), State == 3))

#   tryCatch({
#     order_shalek_cells_by_original_states_custom_ordering_function(cds, root_state = root_state, cells_state_2 = AT1_state, cells_state_3 = AT2_state)
#   },
#   error = function(e) {print(e); NA})
# })

# DDRTree_cds_downsampled_cells_ordered_0.8 = lapply(1:repeat_downsampling_num, function(id) {
#   message('trial number is ', id)
#   cds <- cds_downsampled_cells_ordered_0.8[[id]]
#   # root_state <- row.names(subset(pData(cds), State == 2))
#   # AT1_state <- row.names(subset(pData(cds), State == 1))
#   # AT2_state <- row.names(subset(pData(cds), State == 3))

#   tryCatch({
#     order_shalek_cells_by_original_states_custom_ordering_function(cds, root_state = root_state, cells_state_2 = AT1_state, cells_state_3 = AT2_state)
#   },
#   error = function(e) {print(e); NA})
# })

# ################################################################################################################################################
# #use DDRTree only with diffusion map (from destiny package) as initial dimension reduction method:
# ################################################################################################################################################
# destiny_dm_DDRTree_cds_downsampled_cells_ordered = lapply(cds_downsampled_cells_ordered, function(cds) {
#   message('proportion is ', ncol(cds) / ncol(absolute_cds)) #183 full dataset
#   # root_state <- row.names(subset(pData(cds), State == 2))
#   # AT1_state <- row.names(subset(pData(cds), State == 1))
#   # AT2_state <- row.names(subset(pData(cds), State == 3))
#   tryCatch({
#     order_shalek_cells_by_original_states_custom_ordering_function(cds, root_state = root_state, cells_state_2 = AT1_state, cells_state_3 = AT2_state, initial_method = destiny_diffusionMaps)
#   },
#   error = function(e) {print(e); NA})
# })

# destiny_dm_DDRTree_cds_downsampled_cells_ordered_0.8 = lapply(1:repeat_downsampling_num, function(id) {
#   message('trial number is ', id)
#   cds <- cds_downsampled_cells_ordered_0.8[[id]]
#   # root_state <- row.names(subset(pData(cds), State == 2))
#   # AT1_state <- row.names(subset(pData(cds), State == 1))
#   # AT2_state <- row.names(subset(pData(cds), State == 3))
#   tryCatch({
#     order_shalek_cells_by_original_states_custom_ordering_function(cds, root_state = root_state, cells_state_2 = AT1_state, cells_state_3 = AT2_state, initial_method = destiny_diffusionMaps)
#   },
#   error = function(e) {print(e); NA})
# })

# ##########################################################################################################################################################################
# # run on each piece when the code is buggy
# ##########################################################################################################################################################################

# destiny_dm_DDRTree_cds_downsampled_cells_ordered_0.8 = lapply(1:repeat_downsampling_num, function(id) {
#   message('trial number is ', id)
#   cds <- cds_downsampled_cells_ordered_0.8[[id]]
#   # root_state <- row.names(subset(pData(cds), State == 2))
#   # AT1_state <- row.names(subset(pData(cds), State == 1))
#   # AT2_state <- row.names(subset(pData(cds), State == 3))
#   tryCatch({
#     order_shalek_cells_by_original_states_custom_ordering_function(cds, root_state = root_state, cells_state_2 = AT1_state, cells_state_3 = AT2_state, initial_method = destiny_diffusionMaps)
#   },
#   error = function(e) {print(e); NA})
# })

# save(file = 'RData/DDRTree_cds_downsampled_cells_ordered.RData', DDRTree_cds_downsampled_cells_ordered_0.8, DDRTree_cds_downsampled_cells_ordered, destiny_dm_DDRTree_cds_downsampled_cells_ordered_0.8)

###############################################################################################################################################
# Plot the reordered trajectories
###############################################################################################################################################

cds_downsampled_cells_ordered_trajectories = lapply(cds_downsampled_cells_ordered, plot_monocle_spanning_tree_vectorized)

#save the tree structure across different depth:
if(!exists("trajectory_color_by"))
  trajectory_color_by <- "Time"

for(i in 1:length(cds_downsampled_cells_ordered)) {
    p <- plot_cell_trajectory(cds_downsampled_cells_ordered[[i]], color_by=trajectory_color_by, cell_size=1, cell_link_size = 0.1, show_branch_points = F) +  nm_theme() +
      #plot_spanning_tree(cds_downsampled_cells_ordered[[i]], color_by="interaction(experiment_name, time)", cell_size=0.5, cell_link_size = 0.01) +
         #scale_color_manual(values=shalek_custom_color_scale_plus_states) + nm_theme() +
         scale_size(range = c(0.01, 0.5),  limits = c(0.01, 0.5))
    ggsave(p, filename = paste("tmp/", i, "", downsampled_proportions[i], benchmark_type, '.pdf', sep = '_'), height = 1.1, width = 1.1)
}

#get the number of cells:
lapply(cds_downsampled_cells_ordered[c(3, 6, 9, 12, 15, 17, 21, 24, 27, 30, 33, 36)], ncol)
lapply(cds_downsampled_cells_ordered[c(1, 5, 8, 12, 15, 16, 20, 22, 27, 30, 33, 36)], ncol)
lapply(cds_downsampled_cells_ordered[seq(1, 36, by = 3)], ncol)

############################################################################################################################################################
#Robustness of tree construction under downsampling from 0.1 to 1
############################################################################################################################################################

cell_sampling_cmbn_sets <- expand.grid(1:36, 1:36)

cell_sampling_cmbn_sets_list <- split(cell_sampling_cmbn_sets[, ], f = 1:nrow(cell_sampling_cmbn_sets[, ]))
cell_sampling_benchmark_res_list <- mclapply(cell_sampling_cmbn_sets_list, function(x) {
  cds_1 <- cds_downsampled_cells_ordered[[x[[1]]]]
  cds_2 <- cds_downsampled_cells_ordered[[x[[2]]]]

  if(length(unique(pData(cds_1)$State)) > 3){
    cds_1 <- trimTree(cds_1, num_paths = 2, min_branch_thrsld = 0.1)
  }
  if(length(unique(pData(cds_2)$State)) > 3){
    cds_2 <- trimTree(cds_2, num_paths = 2, min_branch_thrsld = 0.1)
  }

  overlap_cells <- intersect(colnames(cds_1), colnames(cds_2))

  if(length(overlap_cells)){
    t_1 <- pData(cds_1[, overlap_cells])$Pseudotime
    t_2 <- pData(cds_2[, overlap_cells])$Pseudotime

    cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    pearson_cor_res <- cor(t_1, t_2)
    
    #branch assignment:
    clusters_1 <- pData(cds_1[, overlap_cells])$State
    clusters_2 <- pData(cds_2[, overlap_cells])$State
    ClusteringMetrics_res <- calClusteringMetrics(as.character(clusters_1), as.character(clusters_2))

    return(list(cor = cor_res, pearson_cor = pearson_cor_res, cluster = ClusteringMetrics_res))
  }
  else
    return(list(cor = NA, pearson_cor = NA, cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
}, mc.cores = detectCores() - 2)

# #DDRTree without projection
# DDRTree_cell_sampling_benchmark_res_list <- mclapply(cell_sampling_cmbn_sets_list, function(x) {
#   cds_1 <- DDRTree_cds_downsampled_cells_ordered[[x[[1]]]]
#   cds_2 <- DDRTree_cds_downsampled_cells_ordered[[x[[2]]]]

#   # if(length(pData(cds_1)$State) > 3){
#   #   cds_1 <- trimTree(cds_1, num_paths = 2, min_branch_thrsld = 0.1)
#   # }
#   # if(length(pData(cds_2)$State) > 3){
#   #   cds_2 <- trimTree(cds_2, num_paths = 2, min_branch_thrsld = 0.1)
#   # }

#   overlpa_cells <- intersect(row.names(cds_1), row.names(cds_2))

#   t_1 <- cds_1[overlpa_cells, ]$pseudo_time
#   t_2 <- cds_2[overlpa_cells, ]$pseudo_time

#   cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')

#   #branch assignment:
#   clusters_1 <- cds_1[overlpa_cells, ]$cell_state
#   clusters_2 <- cds_2[overlpa_cells, ]$cell_state
#   ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)

#   return(list(cor = cor_res, cluster = ClusteringMetrics_res))
# }, mc.cores = detectCores() - 2)

# #DDRTree + destiny DM without projection
# destiny_dm_DDRTree_cell_sampling_benchmark_res_list <- mclapply(cell_sampling_cmbn_sets_list, function(x) {
#   cds_1 <- destiny_dm_DDRTree_cds_downsampled_cells_ordered[[x[[1]]]]
#   cds_2 <- destiny_dm_DDRTree_cds_downsampled_cells_ordered[[x[[2]]]]

#   # if(length(pData(cds_1)$State) > 3){
#   #   cds_1 <- trimTree(cds_1, num_paths = 2, min_branch_thrsld = 0.1)
#   # }
#   # if(length(pData(cds_2)$State) > 3){
#   #   cds_2 <- trimTree(cds_2, num_paths = 2, min_branch_thrsld = 0.1)
#   # }

#   overlpa_cells <- intersect(row.names(cds_1), row.names(cds_2))

#   t_1 <- cds_1[overlpa_cells, ]$pseudo_time
#   t_2 <- cds_2[overlpa_cells, ]$pseudo_time

#   cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')

#   #branch assignment:
#   clusters_1 <- cds_1[overlpa_cells, ]$cell_state
#   clusters_2 <- cds_1[overlpa_cells, ]$cell_state
#   ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)

#   return(list(cor = cor_res, cluster = ClusteringMetrics_res))
# }, mc.cores = detectCores() - 2)

############################################################################################################################################################
#Robustness of tree construction under 0.8 downsampling rate
############################################################################################################################################################

#robustness of the ordering algorithm based on the DPT paper:
cmbn_sets <- expand.grid(1:repeat_downsampling_num, 1:repeat_downsampling_num)

cmbn_sets_list <- split(cmbn_sets[, ], f = 1:nrow(cmbn_sets[, ]))
benchmark_res_list <- mclapply(cmbn_sets_list, function(x) {
  cds_1 <- cds_downsampled_cells_ordered_0.8[[x[[1]]]]
  cds_2 <- cds_downsampled_cells_ordered_0.8[[x[[2]]]]

  if(is.na(cds_1) | is.na(cds_2))
    return(list(cor = NA, pearson_cor = NA, cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
  else {
    if(length(unique(pData(cds_1)$State)) > 3){
      cds_1 <- trimTree(cds_1, num_paths = 2, min_branch_thrsld = 0.1)
    }
    if(length(unique(pData(cds_2)$State)) > 3){
      cds_2 <- trimTree(cds_2, num_paths = 2, min_branch_thrsld = 0.1)
    }

    overlap_cells <- intersect(colnames(cds_1), colnames(cds_2))
  
    if(length(overlap_cells)) {
  	  t_1 <- pData(cds_1[, overlap_cells])$Pseudotime
  	  t_2 <- pData(cds_2[, overlap_cells])$Pseudotime
  
  	  cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
  	  pearson_cor_res <- cor(t_1, t_2)
  	  
  	  #branch assignment:
  	  clusters_1 <- pData(cds_1[, overlap_cells])$State
  	  clusters_2 <- pData(cds_2[, overlap_cells])$State
  	  ClusteringMetrics_res <- calClusteringMetrics(as.character(clusters_1), as.character(clusters_2))
  
  	  return(list(cor = cor_res, pearson_cor = pearson_cor_res, cluster = ClusteringMetrics_res))
    }
    else
      return(list(cor = NA, pearson_cor = NA, cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
  }
}, mc.cores = detectCores() - 2)

# #DDRTree without projection
# DDRTree_benchmark_res_list <- mclapply(cmbn_sets_list, function(x) {
#   cds_1 <- DDRTree_cds_downsampled_cells_ordered_0.8[[x[[1]]]]
#   cds_2 <- DDRTree_cds_downsampled_cells_ordered_0.8[[x[[2]]]]
#
#   # if(length(pData(cds_1)$State) > 3){
#   #   cds_1 <- trimTree(cds_1, num_paths = 2, min_branch_thrsld = 0.1)
#   # }
#   # if(length(pData(cds_2)$State) > 3){
#   #   cds_2 <- trimTree(cds_2, num_paths = 2, min_branch_thrsld = 0.1)
#   # }
#
#   overlpa_cells <- intersect(row.names(cds_1), row.names(cds_2))
#
#   t_1 <- cds_1[overlpa_cells, ]$pseudo_time
#   t_2 <- cds_2[overlpa_cells, ]$pseudo_time
#
#   cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
#
#   #branch assignment:
#   clusters_1 <- cds_1[overlpa_cells, ]$cell_state
#   clusters_2 <- cds_1[overlpa_cells, ]$cell_state
#   ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
#
#   return(list(cor = cor_res, cluster = ClusteringMetrics_res))
# }, mc.cores = detectCores() - 2)
#
# #DDRTree + destiny DM without projection
# destiny_dm_DDRTree_benchmark_res_list <- mclapply(cmbn_sets_list, function(x) {
#   cds_1 <- destiny_dm_DDRTree_cds_downsampled_cells_ordered_0.8[[x[[1]]]]
#   cds_2 <- destiny_dm_DDRTree_cds_downsampled_cells_ordered_0.8[[x[[2]]]]
#
#   # if(length(pData(cds_1)$State) > 3){
#   #   cds_1 <- trimTree(cds_1, num_paths = 2, min_branch_thrsld = 0.1)
#   # }
#   # if(length(pData(cds_2)$State) > 3){
#   #   cds_2 <- trimTree(cds_2, num_paths = 2, min_branch_thrsld = 0.1)
#   # }
#
#   overlpa_cells <- intersect(row.names(cds_1), row.names(cds_2))
#
#   t_1 <- cds_1[overlpa_cells, ]$pseudo_time
#   t_2 <- cds_2[overlpa_cells, ]$pseudo_time
#
#   cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
#
#   #branch assignment:
#   clusters_1 <- cds_1[overlpa_cells, ]$cell_state
#   clusters_2 <- cds_1[overlpa_cells, ]$cell_state
#   ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
#
#   return(list(cor = cor_res, cluster = ClusteringMetrics_res))
# }, mc.cores = detectCores() - 2)
#

#plot the results:
############################################################################################################################################################
sampling_res_df <- data.frame(kendall.tau = unlist(lapply(benchmark_res_list, function(x) x$cor)),
                              pearson_rho = unlist(lapply(benchmark_res_list, function(x) x$pearson_cor)),
                                  rand_ind = unlist(lapply(benchmark_res_list, function(x) x$cluster[1, 1])),
                                  var_inf = unlist(lapply(benchmark_res_list, function(x) x$cluster[2, 1])),
                                  adj_rand = unlist(lapply(benchmark_res_list, function(x) x$cluster[3, 1]))
)

# DDRTree_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(DDRTree_benchmark_res_list, function(x) x$cor)),
#                               rand_ind = unlist(lapply(DDRTree_benchmark_res_list, function(x) x$cluster[1, 1])),
#                               var_inf = unlist(lapply(DDRTree_benchmark_res_list, function(x) x$cluster[2, 1])),
#                               adj_rand = unlist(lapply(DDRTree_benchmark_res_list, function(x) x$cluster[3, 1]))
# )
#
# destiny_dm_DDRTree_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(destiny_dm_DDRTree_benchmark_res_list, function(x) x$cor)),
#                                       rand_ind = unlist(lapply(destiny_dm_DDRTree_benchmark_res_list, function(x) x$cluster[1, 1])),
#                                       var_inf = unlist(lapply(destiny_dm_DDRTree_benchmark_res_list, function(x) x$cluster[2, 1])),
#                                       adj_rand = unlist(lapply(destiny_dm_DDRTree_benchmark_res_list, function(x) x$cluster[3, 1]))
# )
#

qplot(sampling_res_df$kendall.tau)
qplot(sampling_res_df$rand_ind)
qplot(sampling_res_df$var_inf)
qplot(sampling_res_df$adj_rand)

############################################################################################################################################################
cell_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(cell_sampling_benchmark_res_list, function(x) x$cor)),
                              pearson_rho = unlist(lapply(cell_sampling_benchmark_res_list, function(x) x$pearson_cor)),
                              rand_ind = unlist(lapply(cell_sampling_benchmark_res_list, function(x) x$cluster[1, 1])),
                              var_inf = unlist(lapply(cell_sampling_benchmark_res_list, function(x) x$cluster[2, 1])),
                              adj_rand = unlist(lapply(cell_sampling_benchmark_res_list, function(x) x$cluster[3, 1]))
)
#
# DDRTree_cell_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(DDRTree_cell_sampling_benchmark_res_list, function(x) x$cor)),
#                                    rand_ind = unlist(lapply(DDRTree_cell_sampling_benchmark_res_list, function(x) x$cluster[1, 1])),
#                                    var_inf = unlist(lapply(DDRTree_cell_sampling_benchmark_res_list, function(x) x$cluster[2, 1])),
#                                    adj_rand = unlist(lapply(DDRTree_cell_sampling_benchmark_res_list, function(x) x$cluster[3, 1]))
# )
#
# destiny_dm_DDRTree_cell_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(destiny_dm_DDRTree_cell_sampling_benchmark_res_list, function(x) x$cor)),
#                                                  rand_ind = unlist(lapply(destiny_dm_DDRTree_cell_sampling_benchmark_res_list, function(x) x$cluster[1, 1])),
#                                                  var_inf = unlist(lapply(destiny_dm_DDRTree_cell_sampling_benchmark_res_list, function(x) x$cluster[2, 1])),
#                                                  adj_rand = unlist(lapply(destiny_dm_DDRTree_cell_sampling_benchmark_res_list, function(x) x$cluster[3, 1]))
# )
#
#
# ICA_cell_sampling_res_df_update <- data.frame(kendall.tau = unlist(lapply(ICA_cell_sampling_benchmark_res_list_update, function(x) x$cor)),
#                                        rand_ind = unlist(lapply(ICA_cell_sampling_benchmark_res_list_update, function(x) x$cluster[1, 1])),
#                                        var_inf = unlist(lapply(ICA_cell_sampling_benchmark_res_list_update, function(x) x$cluster[2, 1])),
#                                        adj_rand = unlist(lapply(ICA_cell_sampling_benchmark_res_list_update, function(x) x$cluster[3, 1]))
# )

progressive_downsampling_rng <- (nrow(cell_sampling_res_df) - 35):nrow(cell_sampling_res_df)
valid_cell_sampling_res_df <- cell_sampling_res_df[progressive_downsampling_rng, ]
# DDRTree_valid_cell_sampling_res_df <- DDRTree_cell_sampling_res_df[1261:1296, ]
# ICA_valid_cell_sampling_res_df_update <- ICA_cell_sampling_res_df_update[1261:1296, ]

qplot(names(cds_downsampled_cells_ordered), valid_cell_sampling_res_df$kendall.tau) + xlab('Proportion of original cells') + ylab("Kendall's tau")
qplot(names(cds_downsampled_cells_ordered), valid_cell_sampling_res_df$rand_ind) + xlab('Proportion of original cells') + ylab("Rand index")
qplot(names(cds_downsampled_cells_ordered), valid_cell_sampling_res_df$var_inf) + xlab('Proportion of original cells') + ylab("Variation of information")
qplot(names(cds_downsampled_cells_ordered), valid_cell_sampling_res_df$adj_rand) + xlab('Proportion of original cells') + ylab("Adjusted rand index")

# ############################################################################################################################################################
# #consistency of branch point identification: this requres development
# ############################################################################################################################################################
# branching_timepoint <- unlist(lapply(cds_downsampled_cells_ordered_0.8, function(cds) min(subset(pData(cds), State == 1)[, 'Pseudotime'])))

############################################################################################################################################################
#save the data
############################################################################################################################################################
save.image(paste('./RData/', benchmark_type, '_DDRTree_robustness_analysis.RData', sep = ''))
