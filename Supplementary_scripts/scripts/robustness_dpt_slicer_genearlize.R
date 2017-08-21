software_custom_color_scale <- c("reference" = "#F3756C",
                                 "Wishbone"="#26B24B",
                                 "Monocle 1" = "#00BCC3",
                                 "Slicer" = "#6E95CD",
                                 "DPT" = "#CC71AD",
                                 "Monocle 2" = "#B8A131")
# + scale_color_manual(values=software_custom_color_scale)
software_levels <- c('Monocle 2', 'Monocle 1', 'DPT', 'Wishbone', 'Slicer')
# SI_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/Figures/supplementary_figures/"
# ############################################################################################################################################################
# #run the same analysis using dpt and slicer:
# # Downsampling the number of cells
# # load('./RData/DDRTree_robustness_analysis.RData')
# 
# # library(devtools)
# # load_all('~/Dropbox (Personal)/Projects/monocle-dev')
# library(monocle)
# library(dpt)
# library(SLICER)
# library(xacHelper)
# library(plyr)
# library(stringr)
# library(dplyr)
# library(grid)
# library(gridExtra)
# library(mcclust)
# library(flexclust)
# library(R.utils)
# library(destiny)
# source('~/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/scripts/function.R', echo=TRUE)
# 
# Time <- pData(absolute_cds)$Time
# names(Time) <- colnames(absolute_cds)
# 
# if(!exists('normalize_data'))
#   normalize_data <- T
# 
# ############################################################################################################################################################
# #run dpt
# ############################################################################################################################################################
# dpt_cds_downsampled_cells_ordered = lapply(cds_downsampled_cells_ordered, function(cds) {
#   if(is.na(cds))
#     return(res = NA)
#   else{
#     message('proportion is ', ncol(cds) / ncol(absolute_cds)) #183 full dataset
#     # cds <- estimateSizeFactors(cds)
# 
#     cell_name <- colnames(cds)
#     res <- NA
#     tryCatch({
#         res <- run_new_dpt(cds, normalize = normalize_data)
#       names(res$pt) <- cell_name[!duplicated(t(as.matrix(exprs(cds))))]
#     },
#     error = function(e) {
#       message('there is an error')
#       res <- NA
#       message("res is ", res)
#       #res <- run_new_dpt(cds, normalize = normalize_data) #weird errors when running DPT; each different run gives different results
#       #names(res$pt) <- cell_name[!duplicated(t(as.matrix(exprs(cds))))]
#     }
#     )
# 
#     res
#   }
# })
# 
# # one bug happen at organize_branches call
# dpt_cds_downsampled_cells_ordered_0.8 = lapply(1:length(cds_downsampled_cells_ordered_0.8), function(index) {
#   message('trial number is ', index)
#   res <- NA
# 
#   tryCatch({
#     cds <- cds_downsampled_cells_ordered_0.8[[index]]
#     res <- run_new_dpt(cds, normalize = normalize_data)
#     names(res$pt) <- colnames(cds)[!duplicated(t(as.matrix(exprs(cds))))]
#     },
#     error = function(e) {
#       message('there is an error')
#       res <- NA
#       message("res is ", res)
# 
#     }
#   )
# 
#   res
# })
# 
# # slicer cannot be run because it throws numerical errors all the time
# # ############################################################################################################################################################
# # #run slicer (unfortunately, slicer has some running issues and we cannot finish some of the runs)
# # ############################################################################################################################################################
# # slicer_cds_downsampled_cells_ordered = lapply(cds_downsampled_cells_ordered, function(cds) {
# #   message('proportion is ', ncol(cds) / ncol(absolute_cds)) #183 full dataset
# #   cds <- estimateSizeFactors(cds)
# #   tryCatch(
# #     expr = {
# #       evalWithTimeout({res <- run_slicer(cds)},
# #                       timeout = 900) #don't run the script more than 15 minutes
# #     },
# #     TimeoutException = function(ex) cat("Timeout. Skipping.\n")
# #   )
# #   row.names(res$order_df) <- colnames(cds)
# 
# #   res
# # })
# 
# # save.image(file = './RData/slicer_cds_downsampled_cells_ordered.RData', slicer_cds_downsampled_cells_ordered)
# 
# # cds_downsampled_cells_ordered_0.8 <- cds_downsampled_cells_ordered_0.8[c(26:100)] #avoid cases where we cannot run the software
# 
# # slicer_cds_downsampled_cells_ordered_0.8 = lapply(1:length(cds_downsampled_cells_ordered_0.8), function(index) {
# #   message('trial number is ', index + 25)
# #   cds <- cds_downsampled_cells_ordered_0.8[[index + 25]]
# 
# #   tryCatch(
# #     expr = {
# #       evalWithTimeout({res <- run_slicer(cds)},
# #                       timeout = 900) #don't run the script more than 15 minutes
# #     },
# #     TimeoutException = function(ex) cat("Timeout. Skipping.\n")
# #   )
# 
# #   row.names(res$order_df) <- colnames(cds)
# # })
# 
# # save.image(file = './slicer_cds_downsampled_cells_ordered_0.8.RData', slicer_cds_downsampled_cells_ordered_0.8)
# 
# ############################################################################################################################################################
# #generate the data.frame for benchmarking dpt results
# ############################################################################################################################################################
# dpt_cell_sampling_cmbn_sets <- expand.grid(1:36, 1:36)
# dpt_cell_sampling_cmbn_sets_list <- split(dpt_cell_sampling_cmbn_sets[, ], f = 1:nrow(dpt_cell_sampling_cmbn_sets[, ]))
# 
# max_root_ind <- c()
# dpt_cell_sampling_benchmark_res_list <- mclapply(dpt_cell_sampling_cmbn_sets_list, function(x) {
#   cds_1 <- cds_downsampled_cells_ordered[[x[[1]]]]
#   cds_2 <- cds_downsampled_cells_ordered[[x[[2]]]]
#   overlap_cells <- intersect(colnames(cds_1), colnames(cds_2))
# 
#   if(length(overlap_cells)) {
#     t_1 <- dpt_cds_downsampled_cells_ordered[[x[[1]]]]$pt[overlap_cells]#[overlpa_cells, 'DPT']
#     # t_1.1 <- dpt_cds_downsampled_cells_ordered[[x[[1]]]]$pt[overlpa_cells, 'DPT.1']
#     # t_1.2 <- dpt_cds_downsampled_cells_ordered[[x[[1]]]]$pt[overlpa_cells, 'DPT.2']
# 
#     t_2 <- dpt_cds_downsampled_cells_ordered[[x[[2]]]]$pt[overlap_cells]#[overlpa_cells, 'DPT']
# 
#     # cor_res.0 <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
#     # cor_res.1 <- cor(t_1.1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
#     # cor_res.2 <- cor(t_1.2, t_2, method = 'kendall', use = 'pairwise.complete.obs')
#     #
#     # cor_res <- max(c(cor_res.0, cor_res.1, cor_res.2))
#     cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
#     pearson_cor_res <- cor(t_1, t_2)
#     # max_root_ind <- c(which(c(cor_res.0, cor_res.1, cor_res.2) == cor_res))
# 
#     #branch assignment:
#     # overlpa_cells_update <- overlpa_cells[intersect(which(as.character(dpt_cds_downsampled_cells_ordered[[x[[1]]]]$pt[overlpa_cells, 'Branch']) %in% c('branch 1', 'branch 2', 'branch 3')),
#     #   which(as.character(dpt_cds_downsampled_cells_ordered[[x[[2]]]]$pt[overlpa_cells, 'Branch']) %in% c('branch 1', 'branch 2', 'branch 3')))] #remove the unassigned 1,2,3, uncertain 1,2,3, cells
#     clusters_1 <- as.character(dpt_cds_downsampled_cells_ordered[[x[[1]]]]$branch[overlap_cells, 1])
#     clusters_2 <- as.character(dpt_cds_downsampled_cells_ordered[[x[[2]]]]$branch[overlap_cells, 1])
#     ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
# 
#     return(list(cor = cor_res, pearson_cor = pearson_cor_res, cluster = ClusteringMetrics_res, max_root_ind = max_root_ind))
#   }
#   else
#     return(list(cor = NA, pearson_cor = NA, cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
# }, mc.cores = detectCores() - 2)
# 
# # ############################################################################################################################################################
# # #knockout experiment, that is removing part of the branch and re-calculating the dpt again
# # dpt_knockout_cds_downsampled_cells_ordered = lapply(cds_downsampled_cells_ordered, function(cds) {
# #   if(is.na(cds)) {
# #     return(NA)
# #   }
# #   else {
# #     message('proportion is ', ncol(cds) / ncol(absolute_cds)) #183 full dataset
# #     cds <- estimateSizeFactors(cds)
# #     cds_remove_AT1 <- setdiff(colnames(cds), AT1_state)
# #
# #     tryCatch({
# #       res <- run_new_dpt(cds[, cds_remove_AT1], branching = T)
# #       names(res$pt) <- colnames(cds[, cds_remove_AT1])
# #     },
# #     error = function(e) {
# #       message('there is an error', e)
# #       res <- NA
# #       message("res is ", res)
# #
# #     }
# #     )
# #
# #     res
# #   }
# # })
# #
# # dpt_knockout_cds_downsampled_cells_ordered_0.8 = lapply(1:length(cds_downsampled_cells_ordered_0.8), function(index) {
# #   message('trial number is ', index)
# #   cds <- cds_downsampled_cells_ordered_0.8[[index]]
# #   cds_remove_AT1 <- setdiff(colnames(cds), AT1_state)
# #
# #   tryCatch({
# #     res <- run_new_dpt(cds[, cds_remove_AT1], branching = T)
# #     names(res$pt) <- colnames(cds[, cds_remove_AT1])
# #   },
# #   error = function(e) {
# #     message('there is an error', e)
# #     res <- NA
# #     message("length of cds_remove_AT1 is ", length(cds_remove_AT1))
# #
# #   }
# #   )
# #
# #   res
# # })
# #
# # dpt_knockout_cell_sampling_benchmark_res_list <- mclapply(dpt_cell_sampling_cmbn_sets_list, function(x) {
# #   print(x)
# #   cds_1 <- cds_downsampled_cells_ordered[[x[[1]]]]
# #   cds_2 <- cds_downsampled_cells_ordered[[x[[2]]]]
# #
# #   if(any(c(is.na(cds_1), is.na(cds_2)))) {
# #     return(list(cor = NA, cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
# #   }
# #   else {
# #     overlap_cells <- intersect(colnames(cds_1), colnames(cds_2))
# #
# #     t_1 <- dpt_knockout_cds_downsampled_cells_ordered[[x[[1]]]]$pt[overlap_cells] #, 'DPT'
# #     t_2 <- dpt_knockout_cds_downsampled_cells_ordered[[x[[2]]]]$pt[overlap_cells] #, 'DPT'
# #
# #     cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
# #
# #     #branch assignment:
# #     # overlpa_cells_update <- overlpa_cells[intersect(which(as.character(dpt_knockout_cds_downsampled_cells_ordered[[x[[1]]]]$pt[overlpa_cells, 'Branch']) %in% c('branch 1', 'branch 2', 'branch 3')),
# #     #                                                 which(as.character(dpt_knockout_cds_downsampled_cells_ordered[[x[[2]]]]$pt[overlpa_cells, 'Branch']) %in% c('branch 1', 'branch 2', 'branch 3')))] #remove the unassigned 1,2,3, uncertain 1,2,3, cells
# #     clusters_1 <- as.character(dpt_knockout_cds_downsampled_cells_ordered[[x[[1]]]]$branch[overlap_cells, 1])
# #     clusters_2 <- as.character(dpt_knockout_cds_downsampled_cells_ordered[[x[[2]]]]$branch[overlap_cells, 1])
# #
# #     ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
# #
# #     return(list(cor = cor_res, cluster = ClusteringMetrics_res))
# #   }
# # }, mc.cores = detectCores() - 2)
# 
# ############################################################################################################################################################
# dpt_cmbn_sets <- expand.grid(1:length(cds_downsampled_cells_ordered_0.8), 1:length(cds_downsampled_cells_ordered_0.8))
# 
# dpt_cmbn_sets_list <- split(dpt_cmbn_sets[, ], f = 1:nrow(dpt_cmbn_sets[, ]))
# dpt_benchmark_res_list <- mclapply(dpt_cmbn_sets_list, function(x) {
#   message(x)
#   print(x)
#   cds_1 <- cds_downsampled_cells_ordered_0.8[[x[[1]]]]
#   cds_2 <- cds_downsampled_cells_ordered_0.8[[x[[2]]]]
# 
#   if(any(c(is.na(cds_1), is.na(cds_2)))) {
#     return(list(cor = NA, pearson_cor = NA, cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
#   }
#   else {  # if(any(x == 66)){
#     #   return(list(cor = NA, cluster = data.frame(randIndex = rep(NA, 3), Type = c('rand index', 'variation of information', 'adjusted rand index'))))
#     # }
#     overlap_cells <- intersect(colnames(cds_1), colnames(cds_2))
# 
#     t_1 <- dpt_cds_downsampled_cells_ordered_0.8[[x[[1]]]]$pt[overlap_cells] #[overlpa_cells, 'DPT']
#     # t_1.1 <- dpt_cds_downsampled_cells_ordered_0.8[[x[[1]]]]$pt[overlpa_cells, 'DPT.1']
#     # t_1.2 <- dpt_cds_downsampled_cells_ordered_0.8[[x[[1]]]]$pt[overlpa_cells, 'DPT.2']
# 
#     t_2 <- dpt_cds_downsampled_cells_ordered_0.8[[x[[2]]]]$pt[overlap_cells]
# 
#     # cor_res.0 <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
#     # cor_res.1 <- cor(t_1.1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
#     # cor_res.2 <- cor(t_1.2, t_2, method = 'kendall', use = 'pairwise.complete.obs')
#     #
#     # cor_res <- max(c(cor_res.0, cor_res.1, cor_res.2))
#     cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
#     pearson_cor_res <- cor(t_1, t_2)
#     #branch assignment:
#     # overlpa_cells_update <- overlpa_cells[intersect(which(as.character(dpt_cds_downsampled_cells_ordered_0.8[[x[[1]]]]$pt[overlpa_cells, 'Branch']) %in% c('branch 1', 'branch 2', 'branch 3')),
#     #                                                 which(as.character(dpt_cds_downsampled_cells_ordered_0.8[[x[[2]]]]$pt[overlpa_cells, 'Branch']) %in% c('branch 1', 'branch 2', 'branch 3')))] #remove the unassigned 1,2,3, uncertain 1,2,3, cells
#     clusters_1 <- as.character(dpt_cds_downsampled_cells_ordered_0.8[[x[[1]]]]$branch[overlap_cells, 1])
#     clusters_2 <- as.character(dpt_cds_downsampled_cells_ordered_0.8[[x[[2]]]]$branch[overlap_cells, 1])
# 
#     ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
# 
#     return(list(cor = cor_res, pearson_cor = pearson_cor_res, cluster = ClusteringMetrics_res))
#   }
# }, mc.cores = detectCores() - 2)
# 
# # #knockout:
# # dpt_knockout_benchmark_res_list <- mclapply(cmbn_sets_list, function(x) {
# 
# #   cds_1 <- cds_downsampled_cells_ordered_0.8[[x[[1]]]]
# #   cds_2 <- cds_downsampled_cells_ordered_0.8[[x[[2]]]]
# 
# #   if(any(c(is.na(cds_1), is.na(cds_2)))) {
# #     return(list(cor = NA, cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
# #   }
# #   else {
# #     overlap_cells <- intersect(colnames(cds_1), colnames(cds_2))
# 
# #     t_1 <- dpt_knockout_cds_downsampled_cells_ordered_0.8[[x[[1]]]]$pt[overlap_cells] #[overlpa_cells, 'DPT']
# #     t_2 <- dpt_knockout_cds_downsampled_cells_ordered_0.8[[x[[2]]]]$pt[overlap_cells] #[overlpa_cells, 'DPT']
# 
# #     cor_res <- abs(cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs'))
# 
# #     #branch assignment:
# #     # overlpa_cells_update <- overlpa_cells[intersect(which(as.character(dpt_knockout_cds_downsampled_cells_ordered_0.8[[x[[1]]]]$pt[overlpa_cells, 'Branch']) %in% c('branch 1', 'branch 2', 'branch 3')),
# #     #                                                 which(as.character(dpt_knockout_cds_downsampled_cells_ordered_0.8[[x[[2]]]]$pt[overlpa_cells, 'Branch']) %in% c('branch 1', 'branch 2', 'branch 3')))] #remove the unassigned 1,2,3, uncertain 1,2,3, cells
# #     clusters_1 <- as.character(dpt_knockout_cds_downsampled_cells_ordered_0.8[[x[[1]]]]$branch[overlap_cells, 1])
# #     clusters_2 <- as.character(dpt_knockout_cds_downsampled_cells_ordered_0.8[[x[[2]]]]$branch[overlap_cells, 1])
# #     ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
# 
# #     return(list(cor = cor_res, cluster = ClusteringMetrics_res))
# #   }
# # }, mc.cores = detectCores() - 2)
# 
# ############################################################################################################################################################
# dpt_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(dpt_benchmark_res_list, function(x) x$cor)),
#                                   pearson_rho = unlist(lapply(dpt_benchmark_res_list, function(x) x$pearson_cor)),
#                                   rand_ind = unlist(lapply(dpt_benchmark_res_list, function(x) x$cluster[1, 1])),
#                                   var_inf = unlist(lapply(dpt_benchmark_res_list, function(x) x$cluster[2, 1])),
#                                   adj_rand = unlist(lapply(dpt_benchmark_res_list, function(x) x$cluster[3, 1]))
# )
# #knockout
# # dpt_knockout_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(dpt_knockout_benchmark_res_list, function(x) x$cor)),
# #                                   rand_ind = unlist(lapply(dpt_knockout_benchmark_res_list, function(x) x$cluster[1, 1])),
# #                                   var_inf = unlist(lapply(dpt_knockout_benchmark_res_list, function(x) x$cluster[2, 1])),
# #                                   adj_rand = unlist(lapply(dpt_knockout_benchmark_res_list, function(x) x$cluster[3, 1]))
# # )
# 
# qplot(dpt_sampling_res_df$kendall.tau)
# qplot(dpt_sampling_res_df$rand_ind)
# qplot(dpt_sampling_res_df$var_inf)
# qplot(dpt_sampling_res_df$adj_rand)
# 
# ############################################################################################################################################################
# dpt_cell_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(dpt_cell_sampling_benchmark_res_list, function(x) x$cor)),
#                                        pearson_rho = unlist(lapply(dpt_cell_sampling_benchmark_res_list, function(x) x$pearson_cor)),
#                                        rand_ind = unlist(lapply(dpt_cell_sampling_benchmark_res_list, function(x) x$cluster[1, 1])),
#                                        var_inf = unlist(lapply(dpt_cell_sampling_benchmark_res_list, function(x) x$cluster[2, 1])),
#                                        adj_rand = unlist(lapply(dpt_cell_sampling_benchmark_res_list, function(x) x$cluster[3, 1]))
# )
# 
# # #knockout:
# # dpt_knockout_cell_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(dpt_knockout_cell_sampling_benchmark_res_list, function(x) x$cor)),
# #                                        rand_ind = unlist(lapply(dpt_knockout_cell_sampling_benchmark_res_list, function(x) x$cluster[1, 1])),
# #                                        var_inf = unlist(lapply(dpt_knockout_cell_sampling_benchmark_res_list, function(x) x$cluster[2, 1])),
# #                                        adj_rand = unlist(lapply(dpt_knockout_cell_sampling_benchmark_res_list, function(x) x$cluster[3, 1]))
# # )
# 
# progressive_downsampling_rng <- (nrow(dpt_cell_sampling_res_df) - 35):nrow(dpt_cell_sampling_res_df)
# valid_dpt_cell_sampling_res_df <- dpt_cell_sampling_res_df[progressive_downsampling_rng, ]
# 
# # valid_dpt_knockout_cell_sampling_res_df <- dpt_knockout_cell_sampling_res_df[progressive_downsampling_rng, ]
# 
# qplot(downsampled_proportions, abs(valid_dpt_cell_sampling_res_df$kendall.tau)) + xlab('Proportion of original cells') + ylab("Kendall's tau")
# qplot(downsampled_proportions, valid_dpt_cell_sampling_res_df$rand_ind) + xlab('Proportion of original cells') + ylab("Rand index")
# qplot(downsampled_proportions, valid_dpt_cell_sampling_res_df$var_inf) + xlab('Proportion of original cells') + ylab("Variation of information")
# qplot(downsampled_proportions, valid_dpt_cell_sampling_res_df$adj_rand) + xlab('Proportion of original cells') + ylab("Adjusted rand index")
# 
# ############################################################################################################################################################
# #run wishbone
# ############################################################################################################################################################
# #fraction downsampling
# 
# lapply(1:length(cds_downsampled_cells_ordered), function(index) {
#   data <- t(convert2DRData(cds_downsampled_cells_ordered[[index]], norm_method = 'log'))
#   write.csv(file = paste('./csv_data/Wishbone_test_data/', benchmark_type, '_fraction_downsampling_', index, ".txt", sep = ''), data, quote = F, row.names = T)
#   dm <-  data.frame(DC1 = dpt_cds_downsampled_cells_ordered[[index]]$dm$DC1, DC2 = dpt_cds_downsampled_cells_ordered[[index]]$dm$DC1)
#   row.names(dm) <- names(dpt_cds_downsampled_cells_ordered[[index]]$pt)
#   write.csv(file = paste('./csv_data/Wishbone_test_data/', benchmark_type, '_dpt_dm_fraction_downsampling_', index, ".txt", sep = ''), dm, quote = F, row.names = T)
# })
# 
# #write the root cells:
# fraction_root_cells_vec <- unlist(lapply(1:36, function(x) {
#   cds <- cds_downsampled_cells_ordered[[x]]
#   colnames(cds)[which(pData(cds)$Pseudotime == 0)]
# }))
# write.csv(file = paste('./csv_data/Wishbone_test_data/', benchmark_type, '_fraction_root_cells.txt', sep = ''), fraction_root_cells_vec, quote = F, row.names = T)
# 
# cds_downsampled_cells_ordered_0.8 <- cds_downsampled_cells_ordered_0.8
# #cell downsampling
# lapply(1:length(cds_downsampled_cells_ordered_0.8), function(index) {
#   cds <- cds_downsampled_cells_ordered_0.8[[index]]
#   if(!is.na(dpt_cds_downsampled_cells_ordered_0.8[[index]])) {
#     message(index)
#     data <- t(convert2DRData(cds = cds, norm_method = 'log'))
#     write.csv(file = paste('./csv_data/Wishbone_test_data/', benchmark_type, 'cell_downsampling_', index, ".txt", sep = ''), data, quote = F, row.names = T)
#     dm <-  data.frame(DC1 = dpt_cds_downsampled_cells_ordered_0.8[[index]]$dm$DC1, DC2 = dpt_cds_downsampled_cells_ordered_0.8[[index]]$dm$DC1)
#     row.names(dm) <- names(dpt_cds_downsampled_cells_ordered_0.8[[index]]$pt)
#     write.csv(file = paste('./csv_data/Wishbone_test_data/', benchmark_type, 'dpt_dm_cell_downsampling_', index, ".txt", sep = ''), dm, quote = F, row.names = T)
#   }
# })
# 
# #write the root cells:
# root_cells_vec <- unlist(lapply(cds_downsampled_cells_ordered_0.8, function(x) {
#   if(is.na(x)){
#     NA
#   } else  colnames(x)[which(pData(x)$Pseudotime == 0)]
# }
# ))
# 
# write.csv(file = paste('./csv_data/Wishbone_test_data/', benchmark_type, 'root_cells.txt', sep = ''), root_cells_vec, quote = F, row.names = T)
# 
# if(benchmark_type == 'lung_data'){
#   wishbone_res_repeat <- "lung_repeat_wishbone_df_wp_20.txt"
#   wishbone_res_fraction <- "lung_fractioin_wishbone_df_wp_20.txt"
# }
# if(benchmark_type == 'MAR_seq_data'){
#   wishbone_res_repeat <- "MAR_seq_repeat_wishbone_df_wp_20.txt"
#   wishbone_res_fraction <- "MAR_seq_fractioin_wishbone_df_wp_20.txt"
# }
# if(benchmark_type == 'na_sim_data'){
#   wishbone_res_repeat <- "na_sim_repeat_wishbone_df_wp_20.txt"
#   wishbone_res_fraction <- "na_sim_data_fractioin_wishbone_df_wp_20.txt"
# }
# if(benchmark_type == 'HSMM_myo'){
#   wishbone_res_repeat <- "HSMM_myo_repeat_wishbone_df_wp_20.txt"
#   wishbone_res_fraction <- "HSMM_myo_fractioin_wishbone_df_wp_20.txt"
# }
# if(benchmark_type == 'URMM_all_fig1b'){
#   wishbone_res_repeat <- "URMM_all_fig1b_repeat_wishbone_df_wp_20.txt"
#   wishbone_res_fraction <- "URMM_all_fig1b_wishbone_df_wp_20.txt"
# }
# 
# wishbone_res <- read.table(paste("/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/", wishbone_res_repeat, sep = ''), header = T, sep= '\t', row.names=NULL)
# # wishbone_res <- read.table("/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/all_wishbone_res_df_nw_10.txt", header = T, sep= '\t', row.names=NULL)
# # # qplot(wishbone_res$tSNE1, wishbone_res$tSNE2, color = pData(cds_downsampled_cells_ordered_0.8[[index]])$Time)
# # qplot(wishbone_res$dm1, wishbone_res$dm2, color = pData(cds_downsampled_cells_ordered_0.8[[1]])$Time)
# 
# # ############################################################################################################################################################
# # #****************** error happens here, need to re-run wishbone: ******************
# # qplot(pData(cds_downsampled_cells_ordered_0.8[[1]])[, 'Pseudotime'], wishbone_res$trajectory)
# # qplot(wishbone_res$dm1, wishbone_res$dm2, size = wishbone_res$trajectory)
# #
# # qplot(pData(cds_downsampled_cells_ordered_0.8[[1]])[, 'State'], wishbone_res$branch)
# #
# # qplot(dpt_cds_downsampled_cells_ordered_0.8[[1]]$pt$DPT, subset(wishbone_res, run ==1)$trajectory)
# # qplot(dpt_cds_downsampled_cells_ordered_0.8[[1]]$pt$Branch, subset(wishbone_res, run ==1)$branch)
# #
# # save.image('./RData/robustness_dpt_slicer_tmp.RData')
# # ############################################################################################################################################################
# # #show benchmark results for wishbone (cell robustness downsampling)):
# # ############################################################################################################################################################
# wishbone_cmbn_sets <- expand.grid(unique(wishbone_res$run), unique(wishbone_res$run))
# 
# wishbone_cmbn_sets_list <- split(wishbone_cmbn_sets[, ], f = 1:nrow(wishbone_cmbn_sets[, ]))
# 
# wishbone_benchmark_res_list <- mclapply(wishbone_cmbn_sets_list, function(x) {
#   message(x)
#   cds_1 <- cds_downsampled_cells_ordered_0.8[[x[[1]]]]
#   cds_2 <- cds_downsampled_cells_ordered_0.8[[x[[2]]]]
#   overlpa_cells <- intersect(colnames(cds_1), colnames(cds_2))
# 
#   wishbone_res_1 <- subset(wishbone_res, run == x[[1]])
#   wishbone_res_2 <- subset(wishbone_res, run == x[[2]])
#   row.names(wishbone_res_1) <- wishbone_res_1$X
#   row.names(wishbone_res_2) <- wishbone_res_2$X
# 
#   t_1 <- wishbone_res_1[overlpa_cells, 'trajectory']
#   t_2 <- wishbone_res_2[overlpa_cells, 'trajectory']
# 
#   cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
#   pearson_cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
# 
#   print(t_1)
#   print(t_2)
# 
#   if(cor_res < 0 & is.finite(cor_res)){
#     start_cell_id <- which(t_1 == min(t_1, na.rm = T)) #finding the starting cell (cell with smallest pseudotime) in the overlapping set
#     t_2_update <- abs(t_2 - t_2[start_cell_id])
#     cor_res <- cor(t_1, t_2_update, method = 'kendall', use = 'pairwise.complete.obs')
#   }
# 
#   #branch assignment:
#   clusters_1 <- as.character(wishbone_res_1[overlpa_cells, 'branch'])
#   clusters_2 <- as.character(wishbone_res_2[overlpa_cells, 'branch'])
#   ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
# 
#   return(list(cor = abs(cor_res), pearson_cor = pearson_cor_res, cluster = ClusteringMetrics_res))
# }, mc.cores = detectCores() - 2)
# 
# wishbone_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cor)),
#                                        pearson_rho = unlist(lapply(wishbone_benchmark_res_list, function(x) x$pearson_cor)),
#                                        rand_ind = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cluster[1, 1])),
#                                        var_inf = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cluster[2, 1])),
#                                        adj_rand = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cluster[3, 1]))
# )
# 
# qplot(wishbone_sampling_res_df$kendall.tau)
# qplot(wishbone_sampling_res_df$rand_ind)
# qplot(wishbone_sampling_res_df$var_inf)
# qplot(wishbone_sampling_res_df$adj_rand)
# 
# ############################################################################################################################################################
# #show benchmark results for wishbone (cell robustness downsampling)):
# ############################################################################################################################################################
# fraction_wishbone_res <- read.table(paste("/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/", wishbone_res_fraction, sep = ''), header = T, row.names=NULL, sep = '\t')
# # fraction_wishbone_res <- read.table("/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/fraction_all_wishbone_res_df_nw_10.txt", header = T, row.names=NULL)
# 
# wishbone_cmbn_sets <- expand.grid(unique(fraction_wishbone_res$run), unique(fraction_wishbone_res$run))
# 
# wishbone_cmbn_sets_list <- split(wishbone_cmbn_sets[, ], f = 1:nrow(wishbone_cmbn_sets[, ]))
# 
# wishbone_benchmark_res_list <- mclapply(wishbone_cmbn_sets_list, function(x) {
#   cds_1 <- cds_downsampled_cells_ordered[[x[[1]]]]
#   cds_2 <- cds_downsampled_cells_ordered[[x[[2]]]]
#   overlpa_cells <- intersect(colnames(cds_1), colnames(cds_2))
# 
#   fraction_wishbone_res_1 <- subset(fraction_wishbone_res, run == x[[1]])
#   fraction_wishbone_res_2 <- subset(fraction_wishbone_res, run == x[[2]])
#   row.names(fraction_wishbone_res_1) <- fraction_wishbone_res_1[, 1]
#   row.names(fraction_wishbone_res_2) <- fraction_wishbone_res_2[, 1]
# 
#   t_1 <- fraction_wishbone_res_1[overlpa_cells, 'trajectory']
#   t_2 <- fraction_wishbone_res_2[overlpa_cells, 'trajectory']
# 
#   cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
#   pearson_cor_res <- cor(t_1, t_2)
# 
#   #branch assignment:
#   clusters_1 <- as.character(fraction_wishbone_res_1[overlpa_cells, 'branch'])
#   clusters_2 <- as.character(fraction_wishbone_res_2[overlpa_cells, 'branch'])
#   ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
# 
#   return(list(cor = cor_res, pearson_cor = pearson_cor_res, cluster = ClusteringMetrics_res))
# }, mc.cores = detectCores() - 2)
# 
# wishbone_cell_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cor)),
#                                             pearson_rho = unlist(lapply(wishbone_benchmark_res_list, function(x) x$pearson_cor)),
#                                             rand_ind = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cluster[1, 1])),
#                                             var_inf = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cluster[2, 1])),
#                                             adj_rand = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cluster[3, 1]))
# )
# 
# wishbone_res_dim <- dim(wishbone_cell_sampling_res_df)
# valid_wishbone_cell_sampling_res_df <- wishbone_cell_sampling_res_df[(wishbone_res_dim[1] - length(unique(fraction_wishbone_res$run)) + 1):wishbone_res_dim[1], ]

############################################################################################################################################################
#generate supplementary figures for benchmarking results
############################################################################################################################################################

###cell downsampling: color by both of the dpt and ddrtree result together:
all_valid_cell_sampling_res_df <- Reduce(rbind , list(valid_dpt_cell_sampling_res_df, valid_cell_sampling_res_df, valid_wishbone_cell_sampling_res_df, ICA_valid_cell_sampling_res_df))

all_valid_cell_sampling_res_df$proportion <- c(rep(downsampled_proportions, 2), names(cds_downsampled_cells_ordered)[unique(fraction_wishbone_res$run)], downsampled_proportions)
all_valid_cell_sampling_res_df$Type <- c(rep('DPT', 36), rep('Monocle 2', 36), rep('Wishbone', length(unique(fraction_wishbone_res$run))), rep('Monocle 1', 36))

all_valid_cell_sampling_res_df$Type <- factor(all_valid_cell_sampling_res_df$Type, levels = c('Monocle 2', 'Monocle 1', "DPT", "Wishbone"))
all_valid_cell_sampling_res_df$se <- 0.1

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_cell_downsampling.pdf', sep = ''), width = 2.5, height = 2)
qplot(proportion, abs(all_valid_cell_sampling_res_df$pearson_rho), data = all_valid_cell_sampling_res_df, color = Type, size = 1, geom = 'boxplot') +
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + nm_theme() + scale_size(range = c(0.1, 1)) + monocle_theme_opts() + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'rand_index_comparison_cell_downsampling.pdf', sep = ''), width = 2.5, height = 2)
qplot(proportion, abs(all_valid_cell_sampling_res_df$adj_rand), data = all_valid_cell_sampling_res_df, color = Type, size = 1, geom = 'boxplot') +
  xlab('Proportion of original cells') + ylab("Adjusted rand index") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

process_cell_sampling_res_df <- ddply(all_valid_cell_sampling_res_df, .(Type, proportion), summarize,
                                      mean_kendall.tau = mean(abs(kendall.tau), na.rm = T),
                                      sd_kendall.tau = sd(abs(kendall.tau), na.rm = T),

                                      mean_pearson_rho = mean(abs(pearson_rho), na.rm = T),
                                      sd_pearson_rho = sd(abs(pearson_rho), na.rm = T),

                                      mean_adj_rand = mean(abs(adj_rand), na.rm = T),
                                      sd_adj_rand = sd(abs(adj_rand), na.rm = T),

                                      se = mean(se))
limits <- aes(ymax = mean_adj_rand + sd_adj_rand, ymin=mean_adj_rand - sd_adj_rand)

pdf(paste(SI_fig_dir, benchmark_type, 'rand_index_comparison_cell_downsampling2.pdf', sep = ''), width = 2.5, height = 2)
ggplot(aes(proportion, mean_adj_rand), data = process_cell_sampling_res_df) +
  geom_point(aes(color = Type), position=position_dodge(width=0.1), size = 0.5) + facet_wrap(~Type) +
  geom_errorbar(aes(color = Type, ymax = mean_adj_rand + sd_adj_rand, ymin=mean_adj_rand - sd_adj_rand), position=position_dodge(width=0.9), size = 0.5) +
  xlab('Proportion of original cells') + ylab("Adjusted rand index") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_cell_downsampling2.pdf', sep = ''), width = 2.5, height = 2)
ggplot(aes(proportion, mean_pearson_rho), data = process_cell_sampling_res_df) +
  geom_point(aes(color = Type), position=position_dodge(width=0.9), size = 0.5) + facet_wrap(~Type) +
  geom_errorbar(aes(color = Type, ymax = mean_pearson_rho + sd_pearson_rho, ymin=mean_pearson_rho - sd_pearson_rho), position=position_dodge(width=0.1), size = 0.5) +
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_cell_downsampling2_helper.pdf', sep = ''))
ggplot(aes(proportion, mean_pearson_rho), data = process_cell_sampling_res_df) +
  geom_point(aes(color = Type), position=position_dodge(width=0.1), size = 0.5) +
  geom_errorbar(aes(color = Type, ymax = mean_pearson_rho + sd_pearson_rho, ymin=mean_pearson_rho - sd_pearson_rho), position=position_dodge(width=0.1), size = 0.5) +
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + scale_size(range = c(0.1, 1)) + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

######################################################################################################################################################
# Cole's request to remove DPT + Monocle 2
######################################################################################################################################################
process_cell_sampling_res_df$proportion <- as.numeric(as.character(process_cell_sampling_res_df$proportion))

pdf(paste(SI_fig_dir, benchmark_type, 'rand_index_comparison_cell_downsampling2_cole.pdf', sep = ''), width = 2.5, height = 2)
ggplot(aes(proportion, mean_adj_rand), data = process_cell_sampling_res_df) +
  geom_point(aes(color = Type), position=position_dodge(width=0.1), size = 0.5) + facet_wrap(~Type) +
  geom_errorbar(aes(color = Type, ymax = mean_adj_rand + sd_adj_rand, ymin=mean_adj_rand - sd_adj_rand), position=position_dodge(width=0.9), size = 0.5) +
  xlab('Proportion of original cells') + ylab("Adjusted rand index") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale) + ylab('Robustness (Branch)') +
  scale_x_continuous(limits=c(0.1, 1), breaks = seq(0.1, 1, by = 0.2))
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_cell_downsampling2_cole.pdf', sep = ''), width = 2.5, height = 2)
ggplot(aes(proportion, mean_pearson_rho), data = process_cell_sampling_res_df) +
  geom_point(aes(color = Type), position=position_dodge(width=0.9), size = 0.5) + facet_wrap(~Type) +
  geom_errorbar(aes(color = Type, ymax = mean_pearson_rho + sd_pearson_rho, ymin=mean_pearson_rho - sd_pearson_rho), position=position_dodge(width=0.1), size = 0.5) +
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)+ ylab('Robustness (Pseudotime)') +
  scale_x_continuous(limits=c(0.1, 1), breaks = seq(0.1, 1, by = 0.2))
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_cell_downsampling2_helper.pdf', sep = ''))
ggplot(aes(proportion, mean_pearson_rho), data = process_cell_sampling_res_df) +
  geom_point(aes(color = Type), position=position_dodge(width=0.1), size = 0.5) +
  geom_errorbar(aes(color = Type, ymax = mean_pearson_rho + sd_pearson_rho, ymin=mean_pearson_rho - sd_pearson_rho), position=position_dodge(width=0.1), size = 0.5) +
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + scale_size(range = c(0.1, 1)) + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

####

pdf(paste(SI_fig_dir, benchmark_type, 'rand_index_comparison_cell_downsampling2_cole2.pdf', sep = ''), width = 5, height = 1)
ggplot(aes(proportion, mean_adj_rand), data = process_cell_sampling_res_df) +
  geom_point(aes(color = Type), position=position_dodge(width=0.1), size = 0.5) + facet_wrap(~Type, nrow = 1) +
  geom_errorbar(aes(color = Type, ymax = mean_adj_rand + sd_adj_rand, ymin=mean_adj_rand - sd_adj_rand), position=position_dodge(width=0.9), size = 0.5) +
  xlab('Proportion of original cells') + ylab("Adjusted rand index") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale) + ylab('Robustness\n(Branch)') +
  scale_x_continuous(limits=c(0.1, 1), breaks = seq(0.1, 1, by = 0.2))
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_cell_downsampling2_cole2.pdf', sep = ''), width = 5, height = 1)
ggplot(aes(proportion, mean_pearson_rho), data = process_cell_sampling_res_df) +
  geom_point(aes(color = Type), position=position_dodge(width=0.9), size = 0.5) + facet_wrap(~Type, nrow = 1) +
  geom_errorbar(aes(color = Type, ymax = mean_pearson_rho + sd_pearson_rho, ymin=mean_pearson_rho - sd_pearson_rho), position=position_dodge(width=0.1), size = 0.5) +
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)+ ylab('Robustness\n(Pseudotime)') +
  scale_x_continuous(limits=c(0.1, 1), breaks = seq(0.1, 1, by = 0.2))
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_cell_downsampling2_helper2.pdf', sep = ''))
ggplot(aes(proportion, mean_pearson_rho), data = process_cell_sampling_res_df) +
  geom_point(aes(color = Type), position=position_dodge(width=0.1), size = 0.5) +
  geom_errorbar(aes(color = Type, ymax = mean_pearson_rho + sd_pearson_rho, ymin=mean_pearson_rho - sd_pearson_rho), position=position_dodge(width=0.1), size = 0.5) +
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + scale_size(range = c(0.1, 1)) + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()
############################################################################################################################################################
#create mean and variance plot:
mean_df <- ddply(all_valid_cell_sampling_res_df, .(Type, proportion), function(x) data.frame(mean_kt = mean(x$pearson_rho),
                                                                                             mean_rand_ind = mean(x$rand_ind),
                                                                                             mean_var_inf = mean(x$var_inf),
                                                                                             mean_adj_rand = mean(x$adj_rand) ))

pdf(paste(SI_fig_dir, benchmark_type, 'mean_pearson_rho_comparison_cell_downsampling.pdf', sep = ''), width = 2.5, height = 2)
qplot(proportion, abs(mean_kt), data = mean_df, group = Type, color = Type, size = 1, geom = 'line') +
  xlab('Proportion of original cells') + ylab("Mean Pearson's Rho") + nm_theme() + scale_size(range = c(0.1, 1)) + monocle_theme_opts() + nm_theme() + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'mean_adj_rand_index_comparison_cell_downsampling.pdf', sep = ''), width = 2.5, height = 2)
qplot(proportion, abs(mean_df$mean_adj_rand), data = mean_df, color = Type, group = Type, size = 1, geom = 'line') +
  xlab('Proportion of original cells') + ylab("Mean adjusted rand index") + nm_theme() + scale_size(range = c(0.1, 1)) + monocle_theme_opts() + nm_theme() + ylim(0, 1)
dev.off()

variance_df <- ddply(all_valid_cell_sampling_res_df, .(Type, proportion), function(x) data.frame(var_kt = var(x$pearson_rho),
                                                                                             var_rand_ind = var(x$rand_ind),
                                                                                             var_var_inf = var(x$var_inf),
                                                                                             var_adj_rand = var(x$adj_rand) ))

pdf(paste(SI_fig_dir, benchmark_type, 'var_pearson_rho_comparison_cell_downsampling.pdf', sep = ''), width = 2.5, height = 2)
qplot(proportion, var_kt, data = variance_df, group = Type, color = Type, size = 1, geom = 'line') +
  xlab('Proportion of original cells') + ylab("Variance of Pearson's Rho") + nm_theme() + scale_size(range = c(0.1, 1)) + monocle_theme_opts() + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'var_adj_rand_index_comparison_cell_downsampling.pdf', sep = ''), width = 2.5, height = 2)
qplot(proportion, var_adj_rand, data = variance_df, group = Type, color = Type, size = 1, geom = 'line') +
  xlab('Proportion of original cells') + ylab("Varianceadjusted rand index") + nm_theme() + scale_size(range = c(0.1, 1)) + monocle_theme_opts() + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'var_adj_rand_index_comparison_cell_downsampling_helper.pdf', sep = ''))
qplot(proportion, var_adj_rand, data = variance_df, group = Type, color = Type, size = 1, geom = 'line') +
  xlab('Proportion of original cells') + ylab("Varianceadjusted rand index") + scale_size(range = c(0.1, 1)) + ylim(0, 1)
dev.off()

############################################################################################################################################################
#use ICA downsample data:
# load('./RData/ICA_downsampling.RData')

############################################################################################################################################################
###robustness downsampling: color by both of the dpt and ddrtree result together: 
all_sampling_res_df <- Reduce(rbind , list(dpt_sampling_res_df, sampling_res_df,  wishbone_sampling_res_df, ICA_sampling_res_df)) # ICA_sampling_res_df,

# #use result from DDRTree
#all_sampling_res_df <- Reduce(rbind , list(dpt_sampling_res_df, destiny_dm_DDRTree_sampling_res_df,  wishbone_sampling_res_df, ICA_sampling_res_df)) #DDRTree_sampling_res_df ICA_sampling_res_df,
# all_sampling_res_df$Type <- c(rep('dpt', 9604), rep('Monocle2', 10000), rep('Wishbone', 6724)) #,  rep('Monocle1', 10000)
# 
# #use result from DDRTree and dpt knockout
#all_sampling_res_df <- Reduce(rbind , list(dpt_sampling_res_df, DDRTree_sampling_res_df,  wishbone_sampling_res_df, ICA_sampling_res_df)) # ICA_sampling_res_df,
# all_sampling_res_df$Type <- c(rep('dpt', 10000), rep('Monocle2', 10000), rep('Wishbone', 6724)) #,  rep('Monocle1', 10000)

all_sampling_res_df$Type <- c(rep('DPT', nrow(dpt_sampling_res_df)), rep('Monocle 2',  nrow(sampling_res_df)), rep('Wishbone', nrow(wishbone_sampling_res_df)), rep('Monocle 1', nrow(ICA_sampling_res_df)))#,  rep('Monocle1', 10000)
all_sampling_res_df$Type <- revalue(as.character(all_sampling_res_df$Type), 
                                    c("monocle2" = 'Monocle 2', "monocle1" = 'Monocle 1', "dpt" = 'DPT', "wishbone" = 'Wishbone'))
all_sampling_res_df$Type <- factor(all_sampling_res_df$Type, levels = c('Monocle 2', 'Monocle 1', "DPT", "Wishbone")) #dpt (non-uniform branch)

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_robustness.pdf', sep = ''), width = 1, height = 1) 
qplot(Type, pearson_rho, data = all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Pearson's Rho") + monocle_theme_opts()  + xlab('') + 
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_robustness_helper.pdf', sep = '')) 
qplot(Type, adj_rand, data = all_sampling_res_df, color = Type, geom = 'boxplot') + xlab("Pearson's Rho") + ylab('cells') + monocle_theme_opts() + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'rand_index_comparison_robustness.pdf', sep = ''), width = 1, height = 1) 
qplot(Type, adj_rand, data = all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Adjusted Rand index") + nm_theme() + xlab('') + monocle_theme_opts() + nm_theme() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

############################################################################################################################################################
# Cole's request
############################################################################################################################################################
# monocle2 monocle1 dpt wishbone
pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_robustness_cole.pdf', sep = ''), width = 1, height = 1) 
qplot(Type, pearson_rho, data = all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Pearson's Rho") + monocle_theme_opts()  + xlab('') + 
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)+ ylab('Robustness\n(Pseudotime)')
dev.off()
pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_robustness_helper_cole.pdf', sep = '')) 
qplot(Type, adj_rand, data = all_sampling_res_df, color = Type, geom = 'boxplot') + xlab("Pearson's Rho") + ylab('cells') + monocle_theme_opts() + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()
pdf(paste(SI_fig_dir, benchmark_type, 'rand_index_comparison_robustness_cole.pdf', sep = ''), width = 1, height = 1) 
qplot(Type, adj_rand, data = all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Adjusted Rand index") + nm_theme() + xlab('') + monocle_theme_opts() + nm_theme() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale) + ylab('Robustness\n(Branch)')
dev.off()
############################################################################################################################################################
#change to plot the trend of mean as well as the trend of variance for each method: 
#mean vs variance: 
pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_robustness_mean.pdf', sep = ''), width = 1.5, height = 1.5) 
qplot(Type, pearson_rho, data = all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Distribution of Pearson's Rho") + nm_theme()  + xlab('') + monocle_theme_opts() + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_robustness_helper.pdf', sep = ''), width = 1.5, height = 1.5) 
qplot(Type, adj_rand, data = all_sampling_res_df, color = Type, geom = 'boxplot') + xlab("Pearson's Rho") + ylab('cells') + monocle_theme_opts() + nm_theme() + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'rand_index_comparison_robustness_variance.pdf', sep = ''), width = 1.5, height = 1.5) 
qplot(Type, adj_rand, data = all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Adjusted Rand index") + nm_theme() + xlab('')  + monocle_theme_opts() + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

############################################################################################################################################################
# save.image(paste('./RData/', benchmark_type, '_robustness_dpt_slicer_wishbone.RData', sep = ''))
############################################################################################################################################################
