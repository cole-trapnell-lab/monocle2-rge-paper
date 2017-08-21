library(monocle)
library(xacHelper)
library(grid)
library(mcclust)
library(plyr)
# 
# # load('./RData/na_sim_data_robustness_dpt_slicer_wishbone.RData')
# 
##########################################################################################################################################################################
# perform the accuracy analysis with the true pseudotime and branch
##########################################################################################################################################################################
# ICA_cds_downsampled_cells_ordered_0.8
# cds_downsampled_cells_ordered_0.8
neuron_cells <- colnames(absolute_cds)[1:400]
astrocyte_cells <- colnames(absolute_cds)[401:800]

root_state <- row.names(subset(pData(absolute_cds), State == 1))
reorder_cell_cds <- function (cds, root_state) {
  #determine the mapping between original state 1/2/3 and new state 1/2/3:
  overlap_state_1 <- length(intersect(row.names(pData(cds)[pData(cds)$State == 1, ]), root_state))
  overlap_state_2 <- length(intersect(row.names(pData(cds)[pData(cds)$State == 2, ]), root_state))
  overlap_state_3 <- length(intersect(row.names(pData(cds)[pData(cds)$State == 3, ]), root_state))

  #find the state corresponding to the original root state
  overlap_vec <- c(overlap_state_1, overlap_state_2, overlap_state_3)
  max_ind <- which(overlap_vec == max(overlap_vec))

  if(0 %in% pData(cds)[pData(cds)$State == max_ind, 'Pseudotime']) #avoid recalculation of the Pseudotime
    return(cds)

  cds = orderCells(cds, root_state = max_ind)
  if(length(unique(pData(cds)$State)) > 3)
    cds <- trimTree(cds)

  print('pass') #

  return(cds)
}

cds_downsampled_cells_ordered_0.8_update <- lapply(cds_downsampled_cells_ordered_0.8, reorder_cell_cds, root_state)
cds_downsampled_cells_ordered_update <- lapply(cds_downsampled_cells_ordered, reorder_cell_cds, root_state)

pairwise_cal_benchmark_res <- function(cds_1, cds_2) {
  
  overlpa_cells <- intersect(colnames(cds_1), colnames(cds_2))
  
  if(length(unique(pData(cds_1)$State)) > 3){
    cds_1 <- trimTree(cds_1, num_paths = 2, min_branch_thrsld = 0.1)
  }
  if(length(unique(pData(cds_2)$State)) > 3){
    cds_2 <- trimTree(cds_2, num_paths = 2, min_branch_thrsld = 0.1)
  }
  
  overlpa_cells <- intersect(colnames(cds_1), colnames(cds_2))
  
  neuron_t_1 <- pData(cds_1[, intersect(overlpa_cells, neuron_cells)])$Pseudotime
  neuron_t_2 <- pData(cds_2[, intersect(overlpa_cells, neuron_cells)])$Pseudotime
  astrocyte_t_1 <- pData(cds_1[, intersect(overlpa_cells, astrocyte_cells)])$Pseudotime
  astrocyte_t_2 <- pData(cds_2[, intersect(overlpa_cells, astrocyte_cells)])$Pseudotime
  
  cor_res <- c(cor(neuron_t_1, neuron_t_2), cor(astrocyte_t_1, astrocyte_t_2))
  kendall_cor_res <- c(cor(neuron_t_1, neuron_t_2, method = 'kendall', use = 'pairwise.complete.obs'), cor(astrocyte_t_1, astrocyte_t_2, method = 'kendall', use = 'pairwise.complete.obs'))
  
  #branch assignment:
  clusters_1 <- as.character(pData(cds_1[, overlpa_cells])$State)
  clusters_2 <- as.character(pData(cds_2[, overlpa_cells])$State)
  ClusteringMetrics_res <- calClusteringMetrics(as.character(clusters_1), as.character(clusters_2))
  
  return(list(cor = mean(abs(cor_res)), kendall_tau = mean(abs(kendall_cor_res)), cluster = ClusteringMetrics_res, raw_cor = cor_res, raw_kendall_tau = kendall_cor_res))
}

# absolute_cds <- cds_downsampled_cells_ordered_update[[36]]
monocle2_benchmark_res_list <- lapply(cds_downsampled_cells_ordered_0.8_update, function(x) pairwise_cal_benchmark_res(x, absolute_cds)) #
progressive_monocle2_benchmark_res_list <- lapply(cds_downsampled_cells_ordered_update, function(x) pairwise_cal_benchmark_res(x, absolute_cds)) #
monocle_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(monocle2_benchmark_res_list, function(x) x$kendall_tau)),
                                      pearson_rho = unlist(lapply(monocle2_benchmark_res_list, function(x) x$cor)),
                                      rand_ind = unlist(lapply(monocle2_benchmark_res_list, function(x) x$cluster[1, 1])),
                                      var_inf = unlist(lapply(monocle2_benchmark_res_list, function(x) x$cluster[2, 1])),
                                      adj_rand = unlist(lapply(monocle2_benchmark_res_list, function(x) x$cluster[3, 1]))
)
progressive_monocle_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(progressive_monocle2_benchmark_res_list, function(x) x$kendall_tau)),
                                                  pearson_rho = unlist(lapply(progressive_monocle2_benchmark_res_list, function(x) x$cor)),
                                                  rand_ind = unlist(lapply(progressive_monocle2_benchmark_res_list, function(x) x$cluster[1, 1])),
                                                  var_inf = unlist(lapply(progressive_monocle2_benchmark_res_list, function(x) x$cluster[2, 1])),
                                                  adj_rand = unlist(lapply(progressive_monocle2_benchmark_res_list, function(x) x$cluster[3, 1]))
)

#dpt_cds_downsampled_cells_ordered_0.8

dpt_benchmark_res_list <- lapply(1:length(dpt_cds_downsampled_cells_ordered_0.8), function(x) {
  dpt_res <- dpt_cds_downsampled_cells_ordered_0.8[[x]]
  overlap_cells <- intersect(names(dpt_res$pt), colnames(absolute_cds))
  
  if(length(overlap_cells)) {
    neuron_t_1 <- dpt_res$pt[intersect(overlpa_cells, neuron_cells)]#[overlpa_cells, 'DPT']
    neuron_t_2 <- pData(absolute_cds[, intersect(overlpa_cells, neuron_cells)])$Pseudotime#[overlpa_cells, 'DPT']
    
    astrocyte_t_1 <- dpt_res$pt[intersect(overlpa_cells, astrocyte_cells)]#[overlpa_cells, 'DPT']
    astrocyte_t_2 <- pData(absolute_cds[, intersect(overlpa_cells, astrocyte_cells)])$Pseudotime#[overlpa_cells, 'DPT']
    # cor_res.0 <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    # cor_res.1 <- cor(t_1.1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    # cor_res.2 <- cor(t_1.2, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    #
    # cor_res <- max(c(cor_res.0, cor_res.1, cor_res.2))
    cor_res <- c(cor(neuron_t_1, neuron_t_2), cor(astrocyte_t_1, astrocyte_t_2))
    kendall_cor_res <- c(cor(neuron_t_1, neuron_t_2, method = 'kendall', use = 'pairwise.complete.obs'), cor(astrocyte_t_1, astrocyte_t_2, method = 'kendall', use = 'pairwise.complete.obs'))
    
    #branch assignment:
    # overlpa_cells_update <- overlpa_cells[intersect(which(as.character(dpt_cds_downsampled_cells_ordered[[x[[1]]]]$pt[overlpa_cells, 'Branch']) %in% c('branch 1', 'branch 2', 'branch 3')),
    #   which(as.character(dpt_cds_downsampled_cells_ordered[[x[[2]]]]$pt[overlpa_cells, 'Branch']) %in% c('branch 1', 'branch 2', 'branch 3')))] #remove the unassigned 1,2,3, uncertain 1,2,3, cells
    clusters_1 <- as.character(dpt_res$branch[overlap_cells, 1])
    clusters_2 <- as.character(pData(absolute_cds[, overlap_cells])$State)
    ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
    
    return(list(cor = mean(abs(cor_res)), kendall_tau = mean(abs(kendall_cor_res)), cluster = ClusteringMetrics_res, raw_cor = cor_res, raw_kendall_tau = kendall_cor_res))
  }
  else
    return(list(cor = NA, kendall_tau = NA, cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index")), raw_cor = NA, raw_kendall_tau = NA))
})

progressive_dpt_benchmark_res_list <- lapply(1:length(dpt_cds_downsampled_cells_ordered), function(x) {
  dpt_res <- dpt_cds_downsampled_cells_ordered[[x]]
  overlap_cells <- intersect(names(dpt_res$pt), colnames(absolute_cds))
  
  if(length(overlap_cells)) {
    neuron_t_1 <- dpt_res$pt[intersect(overlpa_cells, neuron_cells)]#[overlpa_cells, 'DPT']
    neuron_t_2 <- pData(absolute_cds[, intersect(overlpa_cells, neuron_cells)])$Pseudotime#[overlpa_cells, 'DPT']
    
    astrocyte_t_1 <- dpt_res$pt[intersect(overlpa_cells, astrocyte_cells)]#[overlpa_cells, 'DPT']
    astrocyte_t_2 <- pData(absolute_cds[, intersect(overlpa_cells, astrocyte_cells)])$Pseudotime#[overlpa_cells, 'DPT']
    # cor_res.0 <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    # cor_res.1 <- cor(t_1.1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    # cor_res.2 <- cor(t_1.2, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    #
    # cor_res <- max(c(cor_res.0, cor_res.1, cor_res.2))
    cor_res <- c(cor(neuron_t_1, neuron_t_2, use = "na.or.complete"), cor(astrocyte_t_1, astrocyte_t_2, use = "na.or.complete"))
    kendall_cor_res <- c(cor(neuron_t_1, neuron_t_2, method = 'kendall', use = 'pairwise.complete.obs'), cor(astrocyte_t_1, astrocyte_t_2, method = 'kendall', use = 'pairwise.complete.obs'))
    
    #branch assignment:
    # overlpa_cells_update <- overlpa_cells[intersect(which(as.character(dpt_cds_downsampled_cells_ordered[[x[[1]]]]$pt[overlpa_cells, 'Branch']) %in% c('branch 1', 'branch 2', 'branch 3')),
    #   which(as.character(dpt_cds_downsampled_cells_ordered[[x[[2]]]]$pt[overlpa_cells, 'Branch']) %in% c('branch 1', 'branch 2', 'branch 3')))] #remove the unassigned 1,2,3, uncertain 1,2,3, cells
    clusters_1 <- as.character(dpt_res$branch[overlap_cells, 1])
    clusters_2 <- as.character(pData(absolute_cds[, overlap_cells])$State)
    ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
    
    return(list(cor = mean(abs(cor_res)), kendall_tau = mean(abs(kendall_cor_res)), cluster = ClusteringMetrics_res, raw_cor = cor_res, raw_kendall_tau = kendall_cor_res))
  }
  else
    return(list(cor = NA, kendall_tau = NA, cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index")), raw_cor = NA, raw_kendall_tau = NA))
})

dpt_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(dpt_benchmark_res_list, function(x) x$kendall_tau)),
                                  pearson_rho = unlist(lapply(dpt_benchmark_res_list, function(x) x$cor)),
                                  rand_ind = unlist(lapply(dpt_benchmark_res_list, function(x) x$cluster[1, 1])),
                                  var_inf = unlist(lapply(dpt_benchmark_res_list, function(x) x$cluster[2, 1])),
                                  adj_rand = unlist(lapply(dpt_benchmark_res_list, function(x) x$cluster[3, 1]))
)
progressive_dpt_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(progressive_dpt_benchmark_res_list, function(x) x$kendall_tau)),
                                              pearson_rho = unlist(lapply(progressive_dpt_benchmark_res_list, function(x) x$cor)),
                                              rand_ind = unlist(lapply(progressive_dpt_benchmark_res_list, function(x) x$cluster[1, 1])),
                                              var_inf = unlist(lapply(progressive_dpt_benchmark_res_list, function(x) x$cluster[2, 1])),
                                              adj_rand = unlist(lapply(progressive_dpt_benchmark_res_list, function(x) x$cluster[3, 1]))
)

# wishbone_res
# fraction_wishbone_res
wishbone_benchmark_res_list <- lapply(unique(wishbone_res$run), function(ind) {
  message(ind)
  subset_wishbone_res <- subset(wishbone_res, run == ind)
  row.names(subset_wishbone_res) <- subset_wishbone_res[, 1]
  overlpa_cells <- intersect(colnames(absolute_cds), row.names(subset_wishbone_res) )
  
  neuron_t_1 <- pData(absolute_cds)[intersect(overlpa_cells, neuron_cells), 'Pseudotime']
  neuron_t_2 <- subset_wishbone_res[intersect(overlpa_cells, neuron_cells), 'trajectory']
  astrocyte_t_1 <- pData(absolute_cds)[intersect(overlpa_cells, astrocyte_cells), 'Pseudotime']
  astrocyte_t_2 <- subset_wishbone_res[intersect(overlpa_cells, astrocyte_cells), 'trajectory']
  
  cor_res <- c(cor(neuron_t_1, neuron_t_2), cor(astrocyte_t_1, astrocyte_t_2))
  kendall_cor_res <- c(cor(neuron_t_1, neuron_t_2, method = 'kendall', use = 'pairwise.complete.obs'), cor(astrocyte_t_1, astrocyte_t_2, method = 'kendall', use = 'pairwise.complete.obs'))
  
  # print(t_1)
  # print(t_2)
  
  if(cor_res < 0 & is.finite(cor_res)){
    start_cell_id <- which(neuron_t_1 == min(neuron_t_1, na.rm = T)) #finding the starting cell (cell with smallest pseudotime) in the overlapping set
    neuron_t_2_update <- abs(neuron_t_2 - neuron_t_2[start_cell_id])
    start_cell_id <- which(astrocyte_t_1 == min(astrocyte_t_1, na.rm = T)) #finding the starting cell (cell with smallest pseudotime) in the overlapping set
    astrocyte_t_2_update <- abs(astrocyte_t_2 - astrocyte_t_2[start_cell_id])
    
    cor_res <- c(cor(neuron_t_1, neuron_t_2_update), cor(astrocyte_t_1, astrocyte_t_2_update))
    kendall_cor_res <- c(cor(neuron_t_1, neuron_t_2_update, method = 'kendall', use = 'pairwise.complete.obs'), cor(astrocyte_t_1, astrocyte_t_2_update, method = 'kendall', use = 'pairwise.complete.obs'))
  }
  
  #branch assignment:
  clusters_1 <- as.character(pData(absolute_cds)[overlpa_cells, 'State'])
  clusters_2 <- as.character(subset_wishbone_res[overlpa_cells, 'branch'])
  ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
  
  return(list(cor = mean(abs(cor_res)), kendall_tau = mean(abs(kendall_cor_res)), cluster = ClusteringMetrics_res, raw_cor = cor_res, raw_kendall_tau = kendall_cor_res))
})

progressive_wishbone_benchmark_res_list <- lapply(unique(fraction_wishbone_res$run), function(ind) {
  message(ind)
  subset_wishbone_res <- subset(fraction_wishbone_res, run == ind)
  row.names(subset_wishbone_res) <- subset_wishbone_res[, 1]
  overlpa_cells <- intersect(colnames(absolute_cds), row.names(subset_wishbone_res) )
  
  neuron_t_1 <- pData(absolute_cds)[intersect(overlpa_cells, neuron_cells), 'Pseudotime']
  neuron_t_2 <- subset_wishbone_res[intersect(overlpa_cells, neuron_cells), 'trajectory']
  astrocyte_t_1 <- pData(absolute_cds)[intersect(overlpa_cells, astrocyte_cells), 'Pseudotime']
  astrocyte_t_2 <- subset_wishbone_res[intersect(overlpa_cells, astrocyte_cells), 'trajectory']
  
  cor_res <- c(cor(neuron_t_1, neuron_t_2), cor(astrocyte_t_1, astrocyte_t_2))
  kendall_cor_res <- c(cor(neuron_t_1, neuron_t_2, method = 'kendall', use = 'pairwise.complete.obs'), cor(astrocyte_t_1, astrocyte_t_2, method = 'kendall', use = 'pairwise.complete.obs'))
  
  # print(t_1)
  # print(t_2)
  
  if(any(cor_res < 0) & any(is.finite(cor_res))){
    start_cell_id <- which(neuron_t_1 == min(neuron_t_1, na.rm = T)) #finding the starting cell (cell with smallest pseudotime) in the overlapping set
    neuron_t_2_update <- abs(neuron_t_2 - neuron_t_2[start_cell_id])
    start_cell_id <- which(astrocyte_t_1 == min(astrocyte_t_1, na.rm = T)) #finding the starting cell (cell with smallest pseudotime) in the overlapping set
    astrocyte_t_2_update <- abs(astrocyte_t_2 - astrocyte_t_2[start_cell_id])
    
    cor_res <- c(cor(neuron_t_1, neuron_t_2_update), cor(astrocyte_t_1, astrocyte_t_2_update))
    kendall_cor_res <- c(cor(neuron_t_1, neuron_t_2_update, method = 'kendall', use = 'pairwise.complete.obs'), cor(astrocyte_t_1, astrocyte_t_2_update, method = 'kendall', use = 'pairwise.complete.obs'))
  }
  
  #branch assignment:
  clusters_1 <- as.character(pData(absolute_cds)[overlpa_cells, 'State'])
  clusters_2 <- as.character(subset_wishbone_res[overlpa_cells, 'branch'])
  ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
  
  return(list(cor = mean(abs(cor_res)), kendall_tau = mean(abs(kendall_cor_res)), cluster = ClusteringMetrics_res, raw_cor = cor_res, raw_kendall_tau = kendall_cor_res))
})

wishbone_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(wishbone_benchmark_res_list, function(x) x$kendall_tau)),
                                       pearson_rho = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cor)),
                                       rand_ind = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cluster[1, 1])),
                                       var_inf = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cluster[2, 1])),
                                       adj_rand = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cluster[3, 1]))
)
progressive_wishbone_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(progressive_wishbone_benchmark_res_list, function(x) x$kendall_tau)),
                                                   pearson_rho = unlist(lapply(progressive_wishbone_benchmark_res_list, function(x) x$cor)),
                                                   rand_ind = unlist(lapply(progressive_wishbone_benchmark_res_list, function(x) x$cluster[1, 1])),
                                                   var_inf = unlist(lapply(progressive_wishbone_benchmark_res_list, function(x) x$cluster[2, 1])),
                                                   adj_rand = unlist(lapply(progressive_wishbone_benchmark_res_list, function(x) x$cluster[3, 1]))
)


downsamling_marker_all_sampling_res_df <- Reduce(rbind , list(dpt_sampling_res_df, monocle_sampling_res_df,  wishbone_sampling_res_df, dpt_sampling_res_df)) # ICA_sampling_res_df,

downsamling_marker_all_sampling_res_df$Type <- c(rep('dpt', nrow(dpt_sampling_res_df)), rep('monocle2',  nrow(monocle_sampling_res_df)), rep('wishbone', nrow(wishbone_sampling_res_df)), rep('monocle1', nrow(dpt_sampling_res_df)))#,  rep('Monocle1', 10000)
downsamling_marker_all_sampling_res_df$Type <- factor(downsamling_marker_all_sampling_res_df$Type, levels = c('monocle2', 'monocle1', "dpt", "wishbone")) #dpt (non-uniform branch)

pdf(paste(SI_fig_dir, benchmark_type, 'real_simulation_pearson_rho_comparison_robustness.pdf', sep = ''), width = 1, height = 1.5)
qplot(Type, pearson_rho, data = downsamling_marker_all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Pearson's Rho") + monocle_theme_opts()  + xlab('') +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'real_simulation_all_kendall_tau_comparison_robustness.pdf', sep = ''), width = 1, height = 1.5)
qplot(Type, kendall.tau, data = downsamling_marker_all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Kendall's tau") + monocle_theme_opts()  + xlab('') +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'real_simulation_ArI_comparison_robustness.pdf', sep = ''), width = 1, height = 1.5)
qplot(Type, adj_rand, data = downsamling_marker_all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Adjusted Rand index") + monocle_theme_opts()  + xlab('') +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()


progressive_all_valid_cell_sampling_res_df <- Reduce(rbind , list(progressive_dpt_sampling_res_df, progressive_monocle_sampling_res_df,  progressive_wishbone_sampling_res_df, progressive_dpt_sampling_res_df)) # ICA_sampling_res_df,

progressive_all_valid_cell_sampling_res_df$proportion <- c(rep(downsampled_proportions, 2), names(cds_downsampled_cells_ordered)[unique(fraction_wishbone_res$run)], downsampled_proportions)
progressive_all_valid_cell_sampling_res_df$Type <- c(rep('dpt', 36), rep('monocle2', 36), rep('wishbone', length(unique(fraction_wishbone_res$run))), rep('monocle1', 36))

progressive_all_valid_cell_sampling_res_df$Type <- factor(progressive_all_valid_cell_sampling_res_df$Type, levels = c('monocle2', 'monocle1', "dpt", "wishbone")) 
progressive_all_valid_cell_sampling_res_df$se <- 0.1 

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_real_simulation_rho_comparison_cell_downsampling.pdf', sep = ''), width = 2.5, height = 2) 
qplot(proportion, abs(progressive_all_valid_cell_sampling_res_df$pearson_rho), data = progressive_all_valid_cell_sampling_res_df, color = Type, size = 1, geom = 'boxplot') +
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + nm_theme() + scale_size(range = c(0.1, 1)) + monocle_theme_opts() + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'rand_index_comparison_real_simulation_downsampling.pdf', sep = ''), width = 2.5, height = 2) 
qplot(proportion, abs(progressive_all_valid_cell_sampling_res_df$adj_rand), data = progressive_all_valid_cell_sampling_res_df, color = Type, size = 1, geom = 'boxplot') + 
  xlab('Proportion of original cells') + ylab("Adjusted rand index") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

progressive_process_cell_sampling_res_df <- ddply(progressive_all_valid_cell_sampling_res_df, .(Type, proportion), summarize, 
                                                  mean_kendall.tau = mean(abs(kendall.tau), na.rm = T), 
                                                  sd_kendall.tau = sd(abs(kendall.tau), na.rm = T), 
                                                  
                                                  mean_pearson_rho = mean(abs(pearson_rho), na.rm = T), 
                                                  sd_pearson_rho = sd(abs(pearson_rho), na.rm = T), 
                                                  
                                                  mean_adj_rand = mean(abs(adj_rand), na.rm = T),
                                                  sd_adj_rand = sd(abs(adj_rand), na.rm = T),
                                                  
                                                  se = mean(se))
limits <- aes(ymax = mean_adj_rand + sd_adj_rand, ymin=mean_adj_rand - sd_adj_rand)

pdf(paste(SI_fig_dir, benchmark_type, 'rand_index_real_simulation_comparison_cell_downsampling2.pdf', sep = ''), width = 2.5, height = 2) 
ggplot(aes(proportion, mean_adj_rand), data = progressive_process_cell_sampling_res_df) + 
  geom_point(aes(color = Type), position=position_dodge(width=0.1), size = 0.5) + facet_wrap(~Type) + 
  geom_errorbar(aes(color = Type, ymax = mean_adj_rand + sd_adj_rand, ymin=mean_adj_rand - sd_adj_rand), position=position_dodge(width=0.9), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("Adjusted Rand index") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_real_simulation_comparison_cell_downsampling2.pdf', sep = ''), width = 2.5, height = 2) 
ggplot(aes(proportion, mean_pearson_rho), data = progressive_process_cell_sampling_res_df) + 
  geom_point(aes(color = Type), position=position_dodge(width=0.9), size = 0.5) + facet_wrap(~Type) + 
  geom_errorbar(aes(color = Type, ymax = mean_pearson_rho + sd_pearson_rho, ymin=mean_pearson_rho - sd_pearson_rho), position=position_dodge(width=0.1), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'kendall_tau_real_simulation_comparison_cell_downsampling2.pdf', sep = ''), width = 2.5, height = 2) 
ggplot(aes(proportion, mean_kendall.tau), data = progressive_process_cell_sampling_res_df) + 
  geom_point(aes(color = Type), position=position_dodge(width=0.9), size = 0.5) + facet_wrap(~Type) + 
  geom_errorbar(aes(color = Type, ymax = mean_kendall.tau + sd_kendall.tau, ymin=mean_kendall.tau - sd_kendall.tau), position=position_dodge(width=0.1), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("Kendall's Tau") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'kendall_tau_real_simulation_comparison_cell_downsampling2_cole.pdf', sep = ''), width = 3, height = 1) 
ggplot(aes(proportion, mean_kendall.tau), data = progressive_process_cell_sampling_res_df) + 
  geom_point(aes(color = Type), position=position_dodge(width=0.9), size = 0.5) + facet_wrap(~Type) + 
  geom_errorbar(aes(color = Type, ymax = mean_kendall.tau + sd_kendall.tau, ymin=mean_kendall.tau - sd_kendall.tau), position=position_dodge(width=0.1), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("Kendall's Tau") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

##########################################################################################################################################################################
# save the dataset 
##########################################################################################################################################################################
save.image(paste('./RData/', benchmark_type, '_real_simulation_na.RData', sep = ''))

##########################################################################################################################################################################
# use diffusion distance to perform the trajectory reconstruction
##########################################################################################################################################################################
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
# na simulation dataset (not working below, we really need to have a better method for assigning branches based on some graph operations )
##################################################################################################################################################################
# root_cell <- paste('Y_', test@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[row.names(subset(pData(valid_subset_GSE72857_cds2), Pseudotime == 0)), 1], sep = '')
colnames(test)
test_res <- extract_ddrtree_ordering_xj(test@minSpanningTree, cellPairwiseDistances(test),
                                        root_cell = 1, verbose=T) # 'Cell_2'
next_node <<- 0
dp_mst <- test_res$state_mst
dp_mst <- as.undirected(dp_mst)
dp = distances(dp_mst, v = V(dp_mst), to = V(dp_mst), weights = NULL)

# res <- monocle:::pq_helper(dp_mst, use_weights=T, root_node=2)
res <- pq_helper(dp_mst, use_weights=T, root_node=which(degree(dp_mst) == 1)[1])

if(is.null(branch_num))
  branch_num <- sum(degree(dp_mst) > 2) + 1

branch_num <- 2 #6 cell types in the end
order_list <- monocle:::extract_good_branched_ordering(res$subtree, res$root, dp, branch_num, FALSE)

cc_ordering <- order_list$ordering_df
row.names(cc_ordering) <- cc_ordering$sample_name

data_ori <- as.matrix(t(exprs(test)))
data_uniq <- data_ori[!duplicated(data_ori), ]
dm <- DiffusionMap(as.matrix(data_ori))
DPT_res <- DPT(dm)
cell_num <- length(DPT_res$DPT1)
# 
# dp <- DPT_res[1:cell_num, 1:cell_num]
# dimnames(dp) <- list(colnames(test)[!duplicated(data_ori)], colnames(test)[!duplicated(data_ori)])
# gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
# dp_mst <- minimum.spanning.tree(gp)
# 
# root_cell <- row.names(subset(pData(test), Pseudotime == 0))
# test_res <- extract_ddrtree_ordering_xj(dp_mst, dp,
#                                         root_cell = root_cell, verbose=T)
plot(test_res$state_mst, layout = layout_as_tree(test_res$state_mst))

qplot(DPT_res@dm$DC1, DPT_res@dm$DC2, color = pData(test)$State)

qplot(DPT_res@dm$DC1, DPT_res@dm$DC2, color = cc_ordering[as.character(test_res$ordering_df[colnames(test), 'cell_state']), 'cell_state'])
# 
# calClusteringMetrics(cc_ordering[as.character(test_res$ordering_df[colnames(valid_subset_GSE72857_cds2), 'cell_state']), 'cell_state'], pData(valid_subset_GSE72857_cds2)$cell_type2)
# calClusteringMetrics(pData(valid_subset_GSE72857_cds2)$State, pData(valid_subset_GSE72857_cds2)$cell_type2)
# 
# #load the data from the Fabian group's processed results: 
# load('./script_for_reproduce/MARSseq_analysis_tutorial.RData')
# calClusteringMetrics(pData(valid_subset_GSE72857_cds2)$cluster, cluster.id[cluster.id[, 1] != 19, 1]) #confirm that index matches up 
# calClusteringMetrics(branching[cluster.id[, 1] != 19], pData(valid_subset_GSE72857_cds2)$cell_type2) #check the result 
# 
# ARI_branches <- data.frame(ARI = c(0.5923145, 0.6566774, 0.7225239), Type = c('DPT (original)', 'DPT + Monocle 2', 'Monocle 2'))
# qplot(ARI, geom = 'bar', data = ARI_branches, stat = 'identity')
# 
# pdf(paste(SI_fig_dir, 'MARS_seq_ARI_branch_cluster.pdf', sep = ''), width = 2, height = 1)
# ggplot(data = ARI_branches, aes(Type, ARI)) + geom_bar(stat = 'identity', aes(fill = Type)) + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab('')
# dev.off()
# 
# pdf(paste(SI_fig_dir, 'MARS_seq_pseudotime_correspondence.pdf', sep = ''), width = 1, height = 1)
# qplot(DPT_res$DPT24, pData(valid_subset_GSE72857_cds2)$Pseudotime, color = pData(valid_subset_GSE72857_cds2)$State, size = I(0.25)) + xlab('DPT diffusion \n pseudotime') + ylab('Monocle 2 \n pseudotime') + nm_theme()
# dev.off()
# 
# pdf(paste(SI_fig_dir, 'MARS_seq_pseudotime_correspondence_helper.pdf', sep = ''))
# qplot(DPT_res$DPT24, pData(valid_subset_GSE72857_cds2)$Pseudotime, color = pData(valid_subset_GSE72857_cds2)$State, size = I(0.25)) + xlab('DPT diffusion pseudotime') + ylab('Monocle 2 pseudotime')
# dev.off()

diameter_path_ordering <- function(cds, num_paths, root_state) {
  root_cell <- monocle::select_root_cell(cds, root_state = root_state)
  dp_mst <- cds@minSpanningTree
  dp <- cellPairwiseDistances(cds)
  diamter_path_old <- get_diameter(dp_mst)
  branchpoint <- c()
  dp_mst_tmp <- delete.vertices(dp_mst, diamter_path_old$name)
  
  for(i in 1:(num_paths - 1)) {
    diamter_path_current <- get_diameter(dp_mst_tmp)
    start_cell_distance <- distances(dp_mst, v = diamter_path_current$name[1], to = diamter_path_old$name)
    branchpoint <- c(branchpoint, colnames(start_cell_distance)[which.min(start_cell_distance)]) #find the ce
    
    diamter_path_old <- diamter_path_current
    dp_mst_tmp <- delete.vertices(dp_mst_tmp, diamter_path_old$name)
  }
 
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
      if (degree(dp_mst, v=parent_node_name) > 2 & parent_node_name %in% branchpoint){
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
  # return(list(ordering_df = ordering_df, state_mst = state_mst))
  
  ordering_df
  pData(cds)$State <- ordering_df$cell_state
  pData(cds)$Pseudotime <- ordering_df$pseudo_time
  return(cds)
}

i <- 20
test <- cds_downsampled_cells_ordered_update_dpt[[i]]
test <- reduceDimension(test, max_components = 2, norm_method = 'none', scaling = F, reduction_method = 'DPT')
test <- diameter_path_ordering(test, num_paths = 2, root_state = Root_state(test)) # 'Cell_2'
plot_cell_trajectory(test)

cds_downsampled_cells_ordered_update_diffusion_pseudotime <- lapply(1:length(cds_downsampled_cells_ordered_update_dpt), function(i) {
  message('current id is ', i)
  test <- cds_downsampled_cells_ordered_update[[i]]
  test <- reduceDimension(test, max_components = 2, norm_method = 'none', scaling = F, reduction_method = 'DPT')
  tryCatch({
    test <- diameter_path_ordering(test, num_paths = 2, root_state = Root_state(test)) # 'Cell_2'
    return(test)
  }, error = function(e) {
    warning('Your initial method throws numerical errors!')
    test <- orderCells(test)
    return(test)
  })
})

cds_downsampled_cells_ordered_0.8_diffusion_pseudotime <- lapply(1:length(cds_downsampled_cells_ordered_0.8_update), function(i) {
  message('current id is ', i)
  test <- cds_downsampled_cells_ordered_update[[i]]
  test <- reduceDimension(test, max_components = 2, norm_method = 'none', scaling = F, reduction_method = 'DPT')
  tryCatch({
    test <- diameter_path_ordering(test, num_paths = 2, root_state = Root_state(test)) # 'Cell_2'
    return(test)
  }, error = function(e) {
    warning('Your initial method throws numerical errors!')
    test <- orderCells(test)
    return(test)
  })
})

# show the results for enhanced Monocle 2 + DPT



