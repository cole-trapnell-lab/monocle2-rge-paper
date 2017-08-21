#the script to run ICA or slicer
# 
# # #functions to run ica ordering; need to run separately because it is extremely slow or require downsampling whne we deal with big datasets
# ICA_order_shalek_cells_by_original_states <- function(cds_subset, root_state, cells_state_2, cells_state_3) {
#   cds_subset = reduceDimension(cds_subset, norm_method = 'log', reduction_method = 'ICA')
#   cds_subset = orderCells(cds_subset, num_paths = 2, reverse = T)
# 
#   #determine the mapping between original state 1/2/3 and new state 1/2/3:
#   overlap_state_1 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 1, ]), root_state))
#   overlap_state_2 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 2, ]), root_state))
#   overlap_state_3 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 3, ]), root_state))
# 
#   #find the state corresponding to the original root state
#   overlap_vec <- c(overlap_state_1, overlap_state_2, overlap_state_3)
#   max_ind <- which(overlap_vec == max(overlap_vec))
# 
#   cds_subset = orderCells(cds_subset, max_ind, num_paths = 2, reverse = T)
#   overlap_state_2 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 2, ]), cells_state_2))
#   overlap_state_3 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 3, ]), cells_state_2))
# 
#   #find the new state corresponding to the original state 2:
#   overlap_state_2_vec <- c(overlap_state_2, overlap_state_3)
#   state_2_max_ind <- which(overlap_state_2_vec == max(overlap_state_2_vec)) #root state
# 
#   if(state_2_max_ind != 1) {
#     State_vec <- pData(cds_subset)$State
#     pData(cds_subset)$State[State_vec == 2] <- 3
#     pData(cds_subset)$State[State_vec == 3] <- 2
#   }
# 
#   cds_subset
# }
# 
# ################################################################################################################################################
# #run ICA for downsampling (run this in a separate session because it takes a long time)
# ################################################################################################################################################
# 
# ICA_cds_downsampled_cells_ordered = lapply(cds_downsampled_cells_ordered, function(cds) {
#   if(exists("sample_cells")){
#     overlap_cells <- sample_cells[sample_cells %in% colnames(cds)]
#     cds <- cds[, overlap_cells]
#   }
#   message('proportion is ', ncol(cds) / ncol(absolute_cds)) #183 full dataset
# 
#   tryCatch({
#     res <- ICA_order_shalek_cells_by_original_states(cds, root_state = root_state, cells_state_2 = AT1_state, cells_state_3 = AT2_state)
#     res
#   },
#   error = function(e) {print(e); NA})
# })
# 
# ICA_cds_downsampled_cells_ordered_0.8 = lapply(1:repeat_downsampling_num, function(id) {
#   message('trial number is ', id)
#   if(exists("sample_cells")){
#     sample_cells_0.8 <- sample(sample_cells, round(length(sample_cells) * 0.8))
#     cds <- cds[, sample_cells_0.8]
#   }
#   else
#     cds <- cds_downsampled_cells_ordered_0.8[[id]]
#   tryCatch({
#     ICA_order_shalek_cells_by_original_states(cds, root_state = root_state, cells_state_2 = AT1_state, cells_state_3 = AT2_state)
#   },
#   error = function(e) {print(e); NA})
# })
# 
# save(file = paste('./RData/', benchmark_type, '_ICA_downsampling.RData', sep = ''), ICA_cds_downsampled_cells_ordered, ICA_cds_downsampled_cells_ordered_0.8)
# 
# load(paste('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/RData/', benchmark_type, '_ICA_downsampling.RData', sep = ''))

if(exists('ICA_cds_downsampled_cells_ordered_0.8_update')){
  ICA_cds_downsampled_cells_ordered_0.8 <- ICA_cds_downsampled_cells_ordered_0.8_update
  rm(ICA_cds_downsampled_cells_ordered_0.8_update)
}
if(exists('ICA_cds_downsampled_cells_ordered_update')) {
  ICA_cds_downsampled_cells_ordered <- ICA_cds_downsampled_cells_ordered_update
  rm(ICA_cds_downsampled_cells_ordered_update)
}

# ICA_cds_downsampled_cells_ordered_0.8 <- lapply(ICA_cds_downsampled_cells_ordered_0.8, function(x) {
#   x <- x[1, ]
#   x} )
# ICA_cds_downsampled_cells_ordered <- lapply(ICA_cds_downsampled_cells_ordered, function(x) {
#   if(is.na(x))
#     x
#   else {
#     x <- x[1, ]
#     x
#     }
#   } )
#save(file = '/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/RData/ICA_downsampling_update_processed', ICA_cds_downsampled_cells_ordered_0.8, ICA_cds_downsampled_cells_ordered)

#ICA
if(!exists("cell_sampling_cmbn_sets_list")) {
  cell_sampling_cmbn_sets <- expand.grid(1:36, 1:36)
  cell_sampling_cmbn_sets_list <- split(cell_sampling_cmbn_sets[, ], f = 1:nrow(cell_sampling_cmbn_sets[, ]))
}

if(!exists("cmbn_sets_list")) {
  cmbn_sets <- expand.grid(1:repeat_downsampling_num, 1:repeat_downsampling_num)
  cmbn_sets_list <- split(cmbn_sets[, ], f = 1:nrow(cmbn_sets[, ]))
}

ICA_cell_sampling_benchmark_res_list <- mclapply(cell_sampling_cmbn_sets_list, function(x) {
  message(x)
  cds_1 <- ICA_cds_downsampled_cells_ordered[[x[[1]]]]
  cds_2 <- ICA_cds_downsampled_cells_ordered[[x[[2]]]]
  
  # if(length(pData(cds_1)$State) > 3){
  #   cds_1 <- trimTree(cds_1, num_paths = 2, min_branch_thrsld = 0.1)
  # }
  # if(length(pData(cds_2)$State) > 3){
  #   cds_2 <- trimTree(cds_2, num_paths = 2, min_branch_thrsld = 0.1)
  # }
  
  overlap_cells <- intersect(colnames(cds_1), colnames(cds_2))
  
  if(length(overlap_cells)) {
    t_1 <- pData(cds_1[, overlap_cells])$Pseudotime
    t_2 <- pData(cds_2[, overlap_cells])$Pseudotime
    
    cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    pearson_cor_res <- cor(t_1, t_2)
    # AT1_lineage <- c(root_state, AT1_state)
    # AT2_lineage <- c(root_state, AT2_state)
    # AT1_t_1 <- pData(cds_1[, overlap_cells[overlap_cells %in% AT1_state]])$Pseudotime
    # AT1_t_2 <- pData(cds_2[, overlap_cells[overlap_cells %in% AT1_state]])$Pseudotime
    # AT2_t_1 <- pData(cds_1[, overlap_cells[overlap_cells %in% AT2_state]])$Pseudotime
    # AT2_t_2 <- pData(cds_2[, overlap_cells[overlap_cells %in% AT2_state]])$Pseudotime
    # 
    # AT1_cor_res <- cor(AT1_t_1, AT1_t_2, method = 'kendall', use = 'pairwise.complete.obs')
    # AT1_pearson_cor_res <- cor(AT1_t_1, AT1_t_2)
    # AT2_cor_res <- cor(AT2_t_1, AT2_t_2, method = 'kendall', use = 'pairwise.complete.obs')
    # AT2_pearson_cor_res <- cor(AT2_t_1, AT2_t_2)
    
    #branch assignment:
    clusters_1 <- pData(cds_1[, overlap_cells])$State
    clusters_2 <- pData(cds_2[, overlap_cells])$State
    ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
    
    # return(list(cor = c(cor_res, AT1_cor_res, AT2_cor_res), pearson_cor = c(pearson_cor_res, AT1_pearson_cor_res, AT2_pearson_cor_res), cluster = ClusteringMetrics_res))
    return(list(cor = cor_res, pearson_cor = pearson_cor_res, cluster = ClusteringMetrics_res))
    # return(list(cor = c(cor_res, AT1_cor_res, AT2_cor_res), pearson_cor = c(pearson_cor_res, AT1_pearson_cor_res, AT2_pearson_cor_res), cluster = ClusteringMetrics_res))
  }
  else
    return(list(cor = NA, pearson_cor = NA, cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
  # return(list(cor = rep(NA, 3), pearson_cor = rep(NA, 3), cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
}, mc.cores = detectCores() - 2)

#ICA
ICA_benchmark_res_list <- mclapply(cmbn_sets_list, function(x) {
  cds_1 <- ICA_cds_downsampled_cells_ordered_0.8[[x[[1]]]]
  cds_2 <- ICA_cds_downsampled_cells_ordered_0.8[[x[[2]]]]
  
  # if(length(pData(cds_1)$State) > 3){
  #   cds_1 <- trimTree(cds_1, num_paths = 2, min_branch_thrsld = 0.1)
  # }
  # if(length(pData(cds_2)$State) > 3){
  #   cds_2 <- trimTree(cds_2, num_paths = 2, min_branch_thrsld = 0.1)
  # }
  
  overlap_cells <- intersect(colnames(cds_1), colnames(cds_2))
  
  if(length(overlap_cells)) {
    t_1 <- pData(cds_1[, overlap_cells])$Pseudotime
    t_2 <- pData(cds_2[, overlap_cells])$Pseudotime
    
    cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    pearson_cor_res <- cor(t_1, t_2)
    
    # AT1_lineage <- c(root_state, AT1_state)
    # AT2_lineage <- c(root_state, AT2_state)
    # AT1_t_1 <- pData(cds_1[, overlap_cells[overlap_cells %in% AT1_state]])$Pseudotime
    # AT1_t_2 <- pData(cds_2[, overlap_cells[overlap_cells %in% AT1_state]])$Pseudotime
    # AT2_t_1 <- pData(cds_1[, overlap_cells[overlap_cells %in% AT2_state]])$Pseudotime
    # AT2_t_2 <- pData(cds_2[, overlap_cells[overlap_cells %in% AT2_state]])$Pseudotime
    # 
    # AT1_cor_res <- cor(AT1_t_1, AT1_t_2, method = 'kendall', use = 'pairwise.complete.obs')
    # AT1_pearson_cor_res <- cor(AT1_t_1, AT1_t_2)
    # AT2_cor_res <- cor(AT2_t_1, AT2_t_2, method = 'kendall', use = 'pairwise.complete.obs')
    # AT2_pearson_cor_res <- cor(AT2_t_1, AT2_t_2)
    
    #branch assignment:
    clusters_1 <- pData(cds_1[, overlap_cells])$State
    clusters_2 <- pData(cds_2[, overlap_cells])$State
    ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
    
    return(list(cor = cor_res, pearson_cor = pearson_cor_res, cluster = ClusteringMetrics_res))
    
    # return(list(cor = c(cor_res, AT1_cor_res, AT2_cor_res), pearson_cor = c(pearson_cor_res, AT1_pearson_cor_res, AT2_pearson_cor_res), cluster = ClusteringMetrics_res))
  }
  else
    return(list(cor = NA, pearson_cor = NA, cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
  # return(list(cor = rep(NA, 3), pearson_cor = rep(NA, 3), cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
}, mc.cores = detectCores() - 2)

ICA_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(ICA_benchmark_res_list, function(x) mean(x$cor))),
                                  pearson_rho = unlist(lapply(ICA_benchmark_res_list, function(x) mean(x$pearson_cor))),
                                  rand_ind = unlist(lapply(ICA_benchmark_res_list, function(x) x$cluster[1, 1])),
                                  var_inf = unlist(lapply(ICA_benchmark_res_list, function(x) x$cluster[2, 1])),
                                  adj_rand = unlist(lapply(ICA_benchmark_res_list, function(x) x$cluster[3, 1]))
)

ICA_cell_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(ICA_cell_sampling_benchmark_res_list, function(x) mean(x$cor))),
                                       pearson_rho = unlist(lapply(ICA_cell_sampling_benchmark_res_list, function(x) mean(x$pearson_cor))),
                                       rand_ind = unlist(lapply(ICA_cell_sampling_benchmark_res_list, function(x) x$cluster[1, 1])),
                                       var_inf = unlist(lapply(ICA_cell_sampling_benchmark_res_list, function(x) x$cluster[2, 1])),
                                       adj_rand = unlist(lapply(ICA_cell_sampling_benchmark_res_list, function(x) x$cluster[3, 1]))
)

progressive_downsampling_rng <- (nrow(ICA_cell_sampling_res_df) - 35):nrow(ICA_cell_sampling_res_df)
# DDRTree_valid_cell_sampling_res_df <- DDRTree_cell_sampling_res_df[1261:1296, ]
ICA_valid_cell_sampling_res_df <- ICA_cell_sampling_res_df[progressive_downsampling_rng, ]

save(file = (paste('./RData/', benchmark_type, '_ICA_robustness_analysis.RData', sep = '')), ICA_valid_cell_sampling_res_df, ICA_sampling_res_df)
# ############################################################################################################################################################
# # run SLICER algorithm 
# #run slicer (unfortunately, slicer has some running issues and we cannot finish some of the runs)
# ############################################################################################################################################################
# slicer_cds_downsampled_cells_ordered = lapply(cds_downsampled_cells_ordered, function(cds) {
#   message('proportion is ', ncol(cds) / ncol(absolute_cds)) #183 full dataset
#   cds <- estimateSizeFactors(cds)
#   tryCatch(
#     expr = {
#       evalWithTimeout({res <- run_slicer(cds)},
#                       timeout = 900) #don't run the script more than 15 minutes
#     },
#     TimeoutException = function(ex) cat("Timeout. Skipping.\n")
#   )
#   row.names(res$order_df) <- colnames(cds)
# 
#   res
# })
# 
# save.image(file = './RData/slicer_cds_downsampled_cells_ordered.RData', slicer_cds_downsampled_cells_ordered)
# 
# valid_cds_downsampled_cells_ordered_0.8 <- cds_downsampled_cells_ordered_0.8 #[c(26:100)] #avoid cases where we cannot run the software
# 
# slicer_cds_downsampled_cells_ordered_0.8 = lapply(1:length(valid_cds_downsampled_cells_ordered_0.8), function(index) {
#   message('trial number is ', index + 25)
#   cds <- cds_downsampled_cells_ordered_0.8[[index + 25]]
# 
#   tryCatch(
#     expr = {
#       evalWithTimeout({res <- run_slicer(cds)},
#                       timeout = 900) #don't run the script more than 15 minutes
#     },
#     TimeoutException = function(ex) cat("Timeout. Skipping.\n")
#   )
# 
#   row.names(res$order_df) <- colnames(cds)
# })
# 
# save.image(file = './slicer_cds_downsampled_cells_ordered_0.8.RData', slicer_cds_downsampled_cells_ordered_0.8)

############################################################################################################################################################
# load the results from the downsampling benchmark with ICA and slicer 
############################################################################################################################################################
# save.image(paste('./RData/', benchmark_type, '_ICA_robustness_analysis.RData', sep = ''))

