############################################################################################################################################################
#this needs further development
############################################################################################################################################################
library(dplyr)
library(RColorBrewer)
library(stringr)
library(dpt)
library(destiny)
library(SLICER)
library(monocle)
library(xacHelper)
library(plyr)
library(stringr)
library(dplyr)
library(grid)
library(gridExtra)
library(mcclust)
library(flexclust)
library(R.utils)
library(modeest)

source('~/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/scripts/function.R', echo=TRUE)
main_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/Figures/main_figures/"
SI_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/Figures/supplementary_figures/"

############################################################################################################################################################
#using monocle2 to order cells and match the states to the original state from the full cds
#########################################################################################################################################################################
software_custom_color_scale <- c("reference" = "#F3756C",
                                 "wishbone"="#26B24B",
                                 "monocle1" = "#00BCC3",
                                 "slicer" = "#6E95CD",
                                 "dpt" = "#CC71AD",
                                 "monocle2" = "#B8A131")
# + scale_color_manual(values=software_custom_color_scale)
software_levels <- c('monocle2', 'monocle1', 'dpt', 'wishbone', 'slicer')

order_shalek_cells_by_original_states <- function(cds_subset, root_state, cells_state_2, cells_state_3, 
                                                  auto_param = auto_param_selection, norm_method = norm_method_data, 
                                                  max_components = num_dim, scaling = !(benchmark_type == 'na_sim_data')) {
  
  cds_subset = reduceDimension(cds_subset, norm_method = norm_method, auto_param_selection = auto_param, max_components = max_components, scaling = scaling)
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

#########################################################################################################################################################################
# progressively downsampling a branch and compare the result to the full dataset:
#########################################################################################################################################################################
set.seed(20161130)
MIN_PROPORTION = 0 #including removing the entire branch
MAX_PROPORTION = 1
STEP = 0.1
REPS_PER = 3
EXTRA_PROPORTIONS = c(0.85, 0.95)

downsampled_proportions = sort(rep(c(seq(MIN_PROPORTION, MAX_PROPORTION, by=STEP), EXTRA_PROPORTIONS) , REPS_PER))
names(downsampled_proportions) = downsampled_proportions  # will tie CDS objects to proportion for later

# monocle 2
progressively_branch_downsampling <- lapply(downsampled_proportions, function(proportion, downsampling_state) {
  message('proportion is ', proportion)
  downsampling_state_cells <- row.names(subset(pData(absolute_cds), State == downsampling_state))
  tmp <- sample(downsampling_state_cells, round(length(downsampling_state_cells) * proportion))
  tmp <- downsampling_state_cells[rank(pData(absolute_cds[, downsampling_state_cells])$Pseudotime) < round(length(downsampling_state_cells) * proportion)] # ensure uneven branch downsampling

  valid_cells <- c(tmp, row.names(subset(pData(absolute_cds), State != downsampling_state)))

  subset_cds <- absolute_cds[, valid_cells %in% colnames(absolute_cds)] #ensure the order of cells matches up with the original order 

  tryCatch({
    order_shalek_cells_by_original_states(subset_cds, root_state = root_state, cells_state_2 = AT1_state[AT1_state %in% colnames(subset_cds)], cells_state_3 = AT2_state[AT2_state %in% colnames(subset_cds)])
  },
  error = function(e) {print(e); NA})
}, downsampling_state)

save(file = paste('./RData/', benchmark_type, 'progressively_branch_downsampling', sep = ''), progressively_branch_downsampling)

#dpt
# the following code fixes the stopifnot stats$g >= gmin problem
DPT <- function (dm, tips = random_root(dm), ..., w_width = 0.1) {
  if (!is(dm, "DiffusionMap"))
    stop("dm needs to be of class DiffusionMap, not ", class(dm))
  if (!length(tips) %in% 1:3)
    stop("you need to specify 1-3 tips, got ", length(tips))
  dpt <- destiny:::dummy_dpt(dm)
  all_cells <- seq_len(nrow(dpt))
  stats <- destiny:::tipstats(dpt, all_cells, tips)

  if(stats$g < 1.1)
    stats$g <- 1.1

  branches <- destiny:::auto_branch(dpt, all_cells, stats, w_width)
  colnames(branches$branch) <- paste0("Branch", seq_len(ncol(branches$branch)))
  colnames(branches$tips) <- paste0("Tips", seq_len(ncol(branches$tips)))
  dpt@branch <- branches$branch
  dpt@tips <- branches$tips
  dpt
}

dpt_progressively_branch_downsampling = lapply(progressively_branch_downsampling, function(cds) {
  if(is.na(cds))
    return(res = NA)
  else{
    message('proportion is ', ncol(cds) / ncol(absolute_cds)) #183 full dataset
    # cds <- estimateSizeFactors(cds)
    
    cell_name <- colnames(cds)
    tryCatch({
      res <- run_new_dpt(cds, normalize = normalize_data)
      names(res$pt) <- cell_name[!duplicated(t(as.matrix(exprs(cds))))]
    },
    error = function(e) {
      message('there is an error')
      res <- NA
      message("res is ", res)
      #res <- run_new_dpt(cds, normalize = normalize_data) #weird errors when running DPT; each different run gives different results
      #names(res$pt) <- cell_name[!duplicated(t(as.matrix(exprs(cds))))]
    }
    )
    
    res
  }
})

save(file = paste(benchmark_type, 'progressively_branch_downsampling.RData', sep = '_'), dpt_progressively_branch_downsampling, progressively_branch_downsampling)

#wishbone 
lapply(1:length(progressively_branch_downsampling), function(index) {
  cds <- progressively_branch_downsampling[[index]]
  if(!is.na(progressively_branch_downsampling[[index]])) {
    message(index)
    data <- t(convert2DRData(cds = cds, norm_method = 'log'))
    write.csv(file = paste('./csv_data/Wishbone_test_data/', benchmark_type, 'branch_downsampling_', index, ".txt", sep = ''), data, quote = F, row.names = T)
    dm <-  data.frame(DC1 = dpt_progressively_branch_downsampling[[index]]$dm$DC1, DC2 = dpt_progressively_branch_downsampling[[index]]$dm$DC1)
    row.names(dm) <- names(dpt_progressively_branch_downsampling[[index]]$pt)
    write.csv(file = paste('./csv_data/Wishbone_test_data/', benchmark_type, 'branch_downsampling_', index, ".txt", sep = ''), dm, quote = F, row.names = T)
  }
})

#write the root cells:
root_cells_vec <- unlist(lapply(progressively_branch_downsampling, function(x) {
  if(is.na(x)){
    NA
  } else  colnames(x)[which(pData(x)$Pseudotime == 0)]
}
))

write.csv(file = paste('./csv_data/Wishbone_test_data/', benchmark_type, 'branch_downsampling_root_cells.txt', sep = ''), root_cells_vec, quote = F, row.names = T)

#########################################################################################################################################################################
#return the results from the dpt run 
cell_sampling_cmbn_sets <- expand.grid(1:length(progressively_branch_downsampling), length(progressively_branch_downsampling))

cell_sampling_cmbn_sets_list <- split(cell_sampling_cmbn_sets[, ], f = 1:nrow(cell_sampling_cmbn_sets[, ]))
cell_sampling_benchmark_res_list <- mclapply(cell_sampling_cmbn_sets_list, function(x) {
  message(x)
  cds_1 <- progressively_branch_downsampling[[x[[1]]]]
  cds_2 <- progressively_branch_downsampling[[x[[2]]]]
  
  if(is.na(cds_1) | is.na(cds_2))
    return(list(pearson_cor_res = NA, cor = NA, cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
  else {

    if(length(unique(pData(cds_1)$State)) > 3){
      cds_1 <- trimTree(cds_1, num_paths = 2, min_branch_thrsld = 0.1)
    }
    if(length(unique(pData(cds_2)$State)) > 3){
      cds_2 <- trimTree(cds_2, num_paths = 2, min_branch_thrsld = 0.1)
    }
    
    overlap_cells <- intersect(colnames(cds_1), colnames(cds_2))
  
    t_1 <- pData(cds_1[, overlap_cells])$Pseudotime
    t_2 <- pData(cds_2[, overlap_cells])$Pseudotime
    
    cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    pearson_cor_res <- cor(t_1, t_2)
    
    #branch assignment:
    clusters_1 <- pData(cds_1[, overlap_cells])$State
    clusters_2 <- pData(cds_2[, overlap_cells])$State
    ClusteringMetrics_res <- calClusteringMetrics(as.character(clusters_1), as.character(clusters_2))
    
    return(list(pearson_cor_res = pearson_cor_res, cor = cor_res, cluster = ClusteringMetrics_res))
  }
}, mc.cores = detectCores() - 2)

dpt_benchmark_res_list <- mclapply(cell_sampling_cmbn_sets_list, function(x) {
  message(x)

  cds_1 <- progressively_branch_downsampling[[x[[1]]]]
  cds_2 <- progressively_branch_downsampling[[x[[2]]]]
  overlap_cells <- intersect(colnames(cds_1), colnames(cds_2))
  
  if(is.na(dpt_progressively_branch_downsampling[[x[[1]]]]) | is.na(dpt_progressively_branch_downsampling[[x[[2]]]])){
    return(list(pearson_cor_res = NA, cor = NA, cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
  }
  else{
    t_1 <- dpt_progressively_branch_downsampling[[x[[1]]]]$pt[overlap_cells] #[overlpa_cells, 'DPT']
  
    t_2 <- dpt_progressively_branch_downsampling[[x[[2]]]]$pt[overlap_cells]
    
    cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    pearson_cor_res = cor(t_1, t_2) 
    
    clusters_1 <- as.character(dpt_progressively_branch_downsampling[[x[[1]]]]$branch[overlap_cells, 1])
    clusters_2 <- as.character(dpt_progressively_branch_downsampling[[x[[2]]]]$branch[overlap_cells, 1])
    
    ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
    print(list(pearson_cor_res = pearson_cor_res, cor = cor_res, cluster = ClusteringMetrics_res))
    return(list(pearson_cor_res = pearson_cor_res, cor = cor_res, cluster = ClusteringMetrics_res))
  }
}, mc.cores = detectCores() - 2)

#monocle 1
monocle1_progressively_branch_downsampling = lapply(progressively_branch_downsampling, function(cds) {
  if(exists("sample_cells")){
    overlap_cells <- sample_cells[sample_cells %in% colnames(cds)]
    cds <- cds[, overlap_cells]
  }
  
  message('proportion is ', ncol(cds) / ncol(absolute_cds)) #183 full dataset
  
  tryCatch({
    res <- ICA_order_shalek_cells_by_original_states(cds, root_state = root_state, cells_state_2 = AT1_state[AT1_state %in% colnames(cds)], cells_state_3 = AT2_state[AT2_state %in% colnames(subset_cds)])
    res
  },
  error = function(e) {print(e); NA})
})

save(file = 'lung_data_progressively_branch_downsampling.RData', dpt_progressively_branch_downsampling, progressively_branch_downsampling, monocle1_progressively_branch_downsampling)


#slicer
slicer_progressively_branch_downsampling = lapply(progressively_branch_downsampling, function(cds) {
  if(exists("sample_cells")){
    overlap_cells <- sample_cells[sample_cells %in% colnames(cds)]
    cds <- cds[, overlap_cells]
  }
  
  message('proportion is ', ncol(cds) / ncol(absolute_cds)) #183 full dataset
  cds <- estimateSizeFactors(cds)
  tryCatch(
    expr = {
      evalWithTimeout({res <- run_slicer(cds)},
                      timeout = 900) #don't run the script more than 15 minutes
    },
    TimeoutException = function(ex) cat("Timeout. Skipping.\n")
  )
  row.names(res$order_df) <- colnames(cds)
  
  res
})

ICA_cell_sampling_benchmark_res_list <- mclapply(cell_sampling_cmbn_sets_list, function(x) {
  message(x)
  cds_1 <- monocle1_progressively_branch_downsampling[[x[[1]]]]
  cds_2 <- monocle1_progressively_branch_downsampling[[x[[2]]]]
  
  overlap_cells <- intersect(colnames(cds_1), colnames(cds_2))
  
  if(length(overlap_cells)) {
    t_1 <- pData(cds_1[, overlap_cells])$Pseudotime
    t_2 <- pData(cds_2[, overlap_cells])$Pseudotime
    
    cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    pearson_cor_res <-  cor(t_1, t_2) 
    
    #branch assignment:
    clusters_1 <- pData(cds_1[, overlap_cells])$State
    clusters_2 <- pData(cds_2[, overlap_cells])$State
    ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
    
    return(list(pearson_cor_res = pearson_cor_res, cor = cor_res, cluster = ClusteringMetrics_res))
  }
  else
    return(list(pearson_cor_res = NA, cor = NA, cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
}, mc.cores = detectCores() - 2)

slicer_benchmark_res_list <- mclapply(cell_sampling_cmbn_sets_list, function(x) {
  message(x)
  
  cds_1 <- progressively_branch_downsampling[[x[[1]]]]
  cds_2 <- progressively_branch_downsampling[[x[[2]]]]
  cds_1_overlap_cells <- which(colnames(cds_2) %in% colnames(cds_1))
  cds_2_overlap_cells <- which(colnames(cds_1) %in% colnames(cds_2))
  
  t_1 <- slicer_progressively_branch_downsampling[[x[[1]]]]$order_df$cells_ordered[cds_1_overlap_cells] #[overlpa_cells, 'DPT']
  
  t_2 <- slicer_progressively_branch_downsampling[[x[[2]]]]$order_dfcells_ordered[cds_2_overlap_cells]
  
  cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
  pearson_cor_res <- cor(t_1, t_2)
  
  clusters_1 <- as.character(dpt_progressively_branch_downsampling[[x[[1]]]]$order_df$branch[cds_1_overlap_cells] )
  clusters_2 <- as.character(dpt_progressively_branch_downsampling[[x[[2]]]]$order_df$branch[cds_1_overlap_cells] )
  
  ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
  print(list(pearson_cor_res = pearson_cor_res, cor = cor_res, cluster = ClusteringMetrics_res))
  return(list(pearson_cor_res = pearson_cor_res, cor = cor_res, cluster = ClusteringMetrics_res))
}, mc.cores = detectCores() - 2)

wishbone_benchmark_res_list <- mclapply(wishbone_cmbn_sets_list, function(x) {
  cds_1 <- cds_downsampled_cells_ordered[[x[[1]]]]
  cds_2 <- cds_downsampled_cells_ordered[[x[[2]]]]
  overlpa_cells <- intersect(colnames(cds_1), colnames(cds_2))
  
  fraction_wishbone_res_1 <- subset(fraction_wishbone_res, run == x[[1]])
  fraction_wishbone_res_2 <- subset(fraction_wishbone_res, run == x[[2]])
  row.names(fraction_wishbone_res_1) <- fraction_wishbone_res_1$row.names
  row.names(fraction_wishbone_res_2) <- fraction_wishbone_res_2$row.names
  
  t_1 <- fraction_wishbone_res_1[overlpa_cells, 'trajectory']
  t_2 <- fraction_wishbone_res_2[overlpa_cells, 'trajectory']
  
  cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
  pearson_cor_res <- cor(t_1, t_2)
  
  #branch assignment:
  clusters_1 <- as.character(fraction_wishbone_res_1[overlpa_cells, 'branch'])
  clusters_2 <- as.character(fraction_wishbone_res_2[overlpa_cells, 'branch'])
  ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
  
  return(list(pearson_cor_res = pearson_cor_res, cor = cor_res, cluster = ClusteringMetrics_res))
}, mc.cores = detectCores() - 2)

monocle2_cell_sampling_res_df <- data.frame(pearson.rho = unlist(lapply(cell_sampling_benchmark_res_list, function(x) x$pearson_cor_res)),
                                            kendall.tau = unlist(lapply(cell_sampling_benchmark_res_list, function(x) x$cor)),
                                            rand_ind = unlist(lapply(cell_sampling_benchmark_res_list, function(x) x$cluster[1, 1])),
                                            var_inf = unlist(lapply(cell_sampling_benchmark_res_list, function(x) x$cluster[2, 1])),
                                            adj_rand = unlist(lapply(cell_sampling_benchmark_res_list, function(x) x$cluster[3, 1]))
)

monocle1_cell_sampling_res_df <- data.frame(pearson.rho = unlist(lapply(ICA_cell_sampling_benchmark_res_list, function(x) x$pearson_cor_res)),
                                       kendall.tau = unlist(lapply(ICA_cell_sampling_benchmark_res_list, function(x) x$cor)),
                                       rand_ind = unlist(lapply(ICA_cell_sampling_benchmark_res_list, function(x) x$cluster[1, 1])),
                                       var_inf = unlist(lapply(ICA_cell_sampling_benchmark_res_list, function(x) x$cluster[2, 1])),
                                       adj_rand = unlist(lapply(ICA_cell_sampling_benchmark_res_list, function(x) x$cluster[3, 1]))
)

dpt_cell_sampling_res_df <- data.frame(pearson.rho = unlist(lapply(dpt_benchmark_res_list, function(x) x$pearson_cor_res)),
                                       kendall.tau = unlist(lapply(dpt_benchmark_res_list, function(x) x$cor)),
                                       rand_ind = unlist(lapply(dpt_benchmark_res_list, function(x) x$cluster[1, 1])),
                                       var_inf = unlist(lapply(dpt_benchmark_res_list, function(x) x$cluster[2, 1])),
                                       adj_rand = unlist(lapply(dpt_benchmark_res_list, function(x) x$cluster[3, 1]))
)

wishbone_cell_sampling_res_df <- data.frame(pearson.rho = unlist(lapply(wishbone_benchmark_res_list, function(x) x$pearson_cor_res)),
                                       kendall.tau = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cor)),
                                       rand_ind = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cluster[1, 1])),
                                       var_inf = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cluster[2, 1])),
                                       adj_rand = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cluster[3, 1]))
)

slicer_cell_sampling_res_df <- data.frame(pearson.rho = unlist(lapply(slicer_benchmark_res_list, function(x) x$pearson_cor_res)),
                                       kendall.tau = unlist(lapply(slicer_benchmark_res_list, function(x) x$cor)),
                                       rand_ind = unlist(lapply(slicer_benchmark_res_list, function(x) x$cluster[1, 1])),
                                       var_inf = unlist(lapply(slicer_benchmark_res_list, function(x) x$cluster[2, 1])),
                                       adj_rand = unlist(lapply(slicer_benchmark_res_list, function(x) x$cluster[3, 1]))
)
 
#########################################################################################################################################################################
#how to present the results? 
#compare the results: 
index <- 2 
dpt_progressively_branch_downsampling[[1]]
progressively_branch_downsampling

which(names(dpt_progressively_branch_downsampling) == '0.2')
x <- dpt_progressively_branch_downsampling[[7]]$dm$DC1
y <- dpt_progressively_branch_downsampling[[7]]$dm$DC2
branch <- dpt_progressively_branch_downsampling[[7]]$branch[, 1]
Time <- pData(progressively_branch_downsampling[[7]])$Time
pdf(paste(SI_fig_dir, "example_unbalanced_trajectory_dpt_time_0.2.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(x, y, color = as.character(branch)) + xlab('Balanced trajectory') + ylab('Unbalanced trajectory') + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

which(names(dpt_progressively_branch_downsampling) == '0.5')
x <- dpt_progressively_branch_downsampling[[16]]$dm$DC1
y <- dpt_progressively_branch_downsampling[[16]]$dm$DC2
branch <- dpt_progressively_branch_downsampling[[16]]$branch[, 1]
Time <- pData(progressively_branch_downsampling[[16]])$Time
pdf(paste(SI_fig_dir, "example_unbalanced_trajectory_dpt_time_0.5.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(x, y, color = as.character(branch)) + xlab('Balanced trajectory') + ylab('Unbalanced trajectory') + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

which(names(dpt_progressively_branch_downsampling) == '0.85')
x <- dpt_progressively_branch_downsampling[[29]]$dm$DC1
y <- dpt_progressively_branch_downsampling[[29]]$dm$DC2
branch <- dpt_progressively_branch_downsampling[[29]]$branch[, 1]
Time <- pData(progressively_branch_downsampling[[29]])$Time
pdf(paste(SI_fig_dir, "example_unbalanced_trajectory_dpt_time_0.85.pdf", sep = ''), height = 1.5, width = 1.5)
qplot(x, y, color = as.character(branch)) + xlab('Balanced trajectory') + ylab('Unbalanced trajectory') + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pdf(paste(SI_fig_dir, "monocle2_example_unbalanced_trajectory_dpt_time_0.2.pdf", sep = ''), height = 1.5, width = 1.5)
plot_cell_trajectory(progressively_branch_downsampling[[7]], show_branch_points = F) + 
  nm_theme() + scale_size(range = c(0.5, 0.5))
dev.off()

pdf(paste(SI_fig_dir, "monocle2_example_unbalanced_trajectory_dpt_time_0.5.pdf", sep = ''), height = 1.5, width = 1.5)
plot_cell_trajectory(progressively_branch_downsampling[[16]], show_branch_points = F) + 
  nm_theme() + scale_size(range = c(0.5, 0.5))
dev.off()

pdf(paste(SI_fig_dir, "monocle2_example_unbalanced_trajectory_dpt_time_0.85.pdf", sep = ''), height = 1.5, width = 1.5)
plot_cell_trajectory(progressively_branch_downsampling[[29]], show_branch_points = F) + 
  nm_theme() + scale_size(range = c(0.5, 0.5))
dev.off()

#########################################################################################################################################################################
#plot the result
#########################################################################################################################################################################
all_sampling_res_df <- Reduce(rbind , list(monocle2_cell_sampling_res_df, monocle1_cell_sampling_res_df, dpt_cell_sampling_res_df, wishbone_cell_sampling_res_df, slicer_cell_sampling_res_df)) 
all_sampling_res_df$Type <- rep(c('monocle2', 'monocle1', 'dpt', 'wishbone', 'slicer'), each = nrow(monocle2_cell_sampling_res_df))
all_sampling_res_df$Type <- factor(all_sampling_res_df$Type, levels = c('monocle2', 'monocle1', 'dpt', 'wishbone', 'slicer'))

all_sampling_res_df <- Reduce(rbind , list(monocle2_cell_sampling_res_df, monocle1_cell_sampling_res_df, dpt_cell_sampling_res_df)) 
all_sampling_res_df$Type <- rep(c('monocle2', 'monocle1', 'dpt'), each = nrow(monocle2_cell_sampling_res_df))
all_sampling_res_df$Type <- factor(all_sampling_res_df$Type, levels = c('monocle2', 'monocle1', 'dpt'))

pdf(paste(SI_fig_dir, "uneven_downsampling_pearson_rho_comparison_robustness.pdf", sep = ''), height = 1.2, width = 1.2)
qplot(Type, abs(pearson.rho), data = all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Pearson's rho") + nm_theme()  + xlab('') + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, "uneven_downsampling_kendall_tau_comparison_robustness.pdf", sep = ''), height = 1.2, width = 1.2)
qplot(Type, kendall.tau, data = all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Kendall's tau") + nm_theme()  + xlab('') + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, "uneven_downsampling_rand_index_comparison_robustness_helper.pdf", sep = ''), height = 1.2, width = 1.2)
qplot(Type, adj_rand, data = all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Adjusted Rand index") + ylab('cells') + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, "uneven_downsampling_rand_index_comparison_robustness.pdf", sep = ''), height = 1.2, width = 1.2)
qplot(Type, adj_rand, data = all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Adjusted Rand index") + nm_theme() + xlab('') + ylim(0, 1)
dev.off()

#########################################################################################################################################################################
# plot the result 
#########################################################################################################################################################################
all_sampling_res_df <- Reduce(rbind , list(monocle2_cell_sampling_res_df, monocle1_cell_sampling_res_df, dpt_cell_sampling_res_df, wishbone_cell_sampling_res_df, wishbone_cell_sampling_res_df)) 
all_sampling_res_df <- Reduce(rbind , list(monocle2_cell_sampling_res_df, monocle1_cell_sampling_res_df, dpt_cell_sampling_res_df)) 

progressive_downsampling_rng <- (nrow(monocle2_cell_sampling_res_df) - 38):nrow(monocle2_cell_sampling_res_df)
valid_monocle2_cell_sampling_res_df <- monocle2_cell_sampling_res_df[progressive_downsampling_rng, ]
valid_monocle1_cell_sampling_res_df <- monocle1_cell_sampling_res_df[progressive_downsampling_rng, ]
valid_dpt_cell_sampling_res_df <- dpt_cell_sampling_res_df[progressive_downsampling_rng, ]
valid_wishbone_cell_sampling_res_df <- wishbone_cell_sampling_res_df[progressive_downsampling_rng, ]
valid_slicer_cell_sampling_res_df <- slicer_cell_sampling_res_df[progressive_downsampling_rng, ]

all_sampling_res_df <- Reduce(rbind , list(valid_monocle2_cell_sampling_res_df, valid_monocle1_cell_sampling_res_df, valid_dpt_cell_sampling_res_df, valid_wishbone_cell_sampling_res_df, valid_slicer_cell_sampling_res_df)) 
all_sampling_res_df <- Reduce(rbind , list(valid_monocle2_cell_sampling_res_df, valid_monocle1_cell_sampling_res_df, valid_dpt_cell_sampling_res_df)) 
all_sampling_res_df$Type <- rep(c('monocle2', 'monocle1', 'dpt', 'wishbone', 'slicer'), each = nrow(valid_monocle2_cell_sampling_res_df))
all_sampling_res_df$Type <- factor(all_sampling_res_df$Type, levels = c('monocle2', 'monocle1', 'dpt', 'wishbone', 'slicer'))

all_sampling_res_df$Type <- rep(c('monocle2', 'monocle1', 'dpt'), each = nrow(valid_monocle2_cell_sampling_res_df))
all_sampling_res_df$Type <- factor(all_sampling_res_df$Type, levels = c('monocle2', 'monocle1', 'dpt'))

all_sampling_res_df$fraction <- as.character(rep(as.numeric(names(progressively_branch_downsampling)), time = length(unique(all_sampling_res_df$Type))))

pdf(paste(SI_fig_dir, "uneven_downsampling_pearson_rho_comparison_trend.pdf", sep = ''), height = 1.2, width = 6)
qplot(fraction, pearson.rho, data = all_sampling_res_df, color = Type, size = 1, geom = 'boxplot', facets = '~Type') +
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + nm_theme() + scale_size(range = c(0.1, 1)) + monocle_theme_opts() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) 
dev.off()

pdf(paste(SI_fig_dir, "uneven_downsampling_kendall_tau_comparison_trend.pdf", sep = ''), height = 1.2, width = 6)
qplot(fraction, kendall.tau, data = all_sampling_res_df, color = Type, size = 1, geom = 'boxplot', facets = '~Type') +
  xlab('Proportion of original cells') + ylab("Kendall's tau") + nm_theme() + scale_size(range = c(0.1, 1)) + monocle_theme_opts() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) 
dev.off()

pdf(paste(SI_fig_dir, "uneven_downsampling_rand_index_comparison_trend_helper.pdf", sep = ''))
qplot(fraction, adj_rand, data = all_sampling_res_df, color = Type, size = 1, geom = 'boxplot', facets = '~Type') +
  xlab('Proportion of original cells') + ylab("Adjusted Rand index") + theme(axis.text.x = element_text(angle = 30, hjust = 1)) 
dev.off()

pdf(paste(SI_fig_dir, "uneven_downsampling_rand_index_comparison_trend.pdf", sep = ''), height = 1.2, width = 6)
qplot(fraction, adj_rand, data = all_sampling_res_df, color = Type, size = 1, geom = 'boxplot', facets = '~Type') +
  xlab('Proportion of original cells') + ylab("Adjusted Rand index") + nm_theme() + scale_size(range = c(0.1, 1)) + monocle_theme_opts() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) 
dev.off()

#########################################################################################################################################################################
# use the mean/variance plot: 
#########################################################################################################################################################################

process_cell_sampling_res_df <- ddply(all_sampling_res_df, .(Type, fraction), summarize, 
                                      mean_kendall.tau = mean(kendall.tau, na.rm = T), 
                                      sd_kendall.tau = sd(kendall.tau, na.rm = T), 
                                      
                                      mean_pearson_rho = mean(pearson.rho, na.rm = T), 
                                      sd_pearson_rho = sd(pearson.rho, na.rm = T), 
                                      
                                      mean_adj_rand = mean(adj_rand, na.rm = T),
                                      sd_adj_rand = sd(adj_rand, na.rm = T))

pdf(paste(SI_fig_dir, benchmark_type, 'progressive_branch_rand_index_comparison_cell_downsampling2.pdf', sep = ''), width = 1.5, height = 1.5) 
ggplot(aes(fraction, abs(mean_adj_rand)), data = process_cell_sampling_res_df) + 
  geom_line(aes(color = Type, group = Type), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("Adjusted rand index") + scale_size(range = c(0.1, 1)) + nm_theme() +
  monocle_theme_opts() + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'progressive_branch_pearson_rho_comparison_cell_downsampling2.pdf', sep = ''), width = 1.5, height = 1.5) 
ggplot(aes(fraction, abs(mean_pearson_rho)), data = process_cell_sampling_res_df) + 
  geom_line(aes(color = Type, group = Type), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'progressive_branch_pearson_rho_comparison_cell_downsampling2_helper.pdf', sep = '')) 
ggplot(aes(fraction, mean_pearson_rho), data = process_cell_sampling_res_df) + 
  geom_point(aes(color = Type), position=position_dodge(width=0.9), size = 0.5) + 
  geom_errorbar(aes(color = Type, ymax = mean_pearson_rho + sd_pearson_rho, ymin=mean_pearson_rho - sd_pearson_rho), position=position_dodge(width=0.9), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + scale_size(range = c(0.1, 1)) + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()


#########################################################################################################################################################################
#calculate the difference of M matrix between different methods: 
#########################################################################################################################################################################

difference_M_matrix <- function(M_0, M_alter){
  
 square_dif <- (as.vector(M_0[row.names(M_alter), row.names(M_alter)] - M_alter))^2
 L2_norm_res <- sqrt(square_dif / nrow(M_alter))
 
 res(L2_norm_res)
}

#########################################################################################################################################################################
#calculate the difference of M matrix between different methods: 

############################################################################################################################################################
save.image(paste('./RData/', benchmark_type, '_uneven_downsampling.RData', sep = ''))
############################################################################################################################################################
# 
# #debug the error leading to the error in the monocle2 downsampling: 
# which(unlist(lapply(progressively_branch_downsampling, is.na)))
# plot_cell_trajectory(progressively_branch_downsampling[[17]])
# 
# progressively_branch_downsampling2 <- lapply(downsampled_proportions, function(proportion, downsampling_state) {
#   message('proportion is ', proportion)
#   downsampling_state_cells <- row.names(subset(pData(absolute_cds), State == downsampling_state))
#   tmp <- sample(downsampling_state_cells, round(length(downsampling_state_cells) * proportion))
#   tmp <- downsampling_state_cells[rank(pData(absolute_cds[, downsampling_state_cells])$Pseudotime) < round(length(downsampling_state_cells) * proportion)] # ensure uneven branch downsampling
#   
#   valid_cells <- c(tmp, row.names(subset(pData(absolute_cds), State != downsampling_state)))
#   
#   subset_cds <- absolute_cds[, valid_cells]
#   
#   tryCatch({
#     order_shalek_cells_by_original_states(subset_cds, root_state = root_state, cells_state_2 = AT1_state[AT1_state %in% colnames(subset_cds)], cells_state_3 = AT2_state[AT2_state %in% colnames(subset_cds)])
#   },
#   error = function(e) {print(e); subset_cds})
# }, downsampling_state)
# 
# # [AT1_state %in% colnames(subset_cds)]
# # [AT2_state %in% colnames(subset_cds)]
# 
# subset_cds <- progressively_branch_downsampling_cds[[3]]
# order_shalek_cells_by_original_states(subset_cds, root_state = root_state, cells_state_2 = AT1_state, cells_state_3 = AT2_state)
# 
# 
