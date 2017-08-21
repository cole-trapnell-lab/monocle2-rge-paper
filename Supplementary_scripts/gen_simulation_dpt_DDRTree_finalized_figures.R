# This file is used to calculate the benchmark results 

library(stringr)

load('./RData/na_sim_data_wishbone_original_monocle2_real_simulation_na.RData')

ind <- 36
dpt_cds_downsampled_cells_ordered[[ind]]$branch[sort(row.names(dpt_cds_downsampled_cells_ordered[[ind]]$branch)), ]

cell_id <- as.numeric(str_split_fixed(row.names(dpt_cds_downsampled_cells_ordered[[ind]]$branch), pattern = '_', 2)[, 2])
dpt_cds_downsampled_cells_ordered[[ind]]$branch[sort(cell_id)]

run_dpt <- function(data, max_components = 3, norm_method = 'log', verbose = F){
   data <- t(data)
   data <- data[!duplicated(data), ]
   dm <- DiffusionMap(as.matrix(data))
   return(dm@eigenvectors[, 1:ncol(data)]) #[, 1:max_components]
}

#test why dpt is better in simulation dataset:
load('./RData/na_sim_data_robustness_dpt_slicer_wishbone.RData')
load('./RData/na_sim_data_a')

load('./RData/na_sim_data_wishbone_real_simulation_na.RData')
plot_cell_trajectory(cds_downsampled_cells_ordered_0.8[[1]])
test <- cds_downsampled_cells_ordered_update[[36]]

qplot(test2[, 1], test2[, 2], color = pData(test)$real_state)

test2 <- (run_dpt(exprs(test)))

run_dpt(test)

Root_state <- function(cds){
  root_cells <- row.names(subset(pData(absolute_cds), State == 1))
  T0_counts <- table(pData(cds[, intersect(root_cells, colnames(cds))])$State)
  as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))])
}

Root_state(cds_downsampled_cells_ordered_update[[1]])
cds_downsampled_cells_ordered_update_dpt <- lapply(1:length(cds_downsampled_cells_ordered_update), function(i) {
  message('current id is ', i)
  test <- cds_downsampled_cells_ordered_update[[i]]
  test <- reduceDimension(test, max_components = 2, norm_method = 'none', scaling = F, reduction_method = run_dpt)
  tryCatch({
    test <- orderCells(test, num_paths = 2)
    test <- orderCells(test, num_paths = 2, root_state = Root_state(test))
    return(test)
    }, error = function(e) {
      warning('Your initial method throws numerical errors!')
      test <- orderCells(test)
      return(test)
    })
})


cds <- cds_downsampled_cells_ordered_0.8_update_dpt[[36]]

cds <- cds_downsampled_cells_ordered_update_dpt[[36]]
test <- reduceDimension(cds, max_components = 2, norm_method = 'none', scaling = F, reduction_method = 'DPT')
test <- orderCells(test, num_paths = 2) #, num_paths = 2
plot_cell_trajectory(test, color_by = 'Pseudotime')

# use the state

test <- reduceDimension(cds, max_components = 2, norm_method = 'none', scaling = F, reduction_method = 'DPT')
test <- orderCells(test, num_paths = 2)

data_ori <- as.matrix(t(exprs(cds_downsampled_cells_ordered_0.8_update_dpt[[14]])))
data_uniq <- data_ori[!duplicated(data_ori), ]
dm <- DiffusionMap(as.matrix(data_ori))
DPT_res <- DPT(dm)
cell_num <- length(DPT_res$DPT1)

dp <- DPT_res[1:cell_num, 1:cell_num]
dimnames(dp) <- list(colnames(cds_downsampled_cells_ordered_0.8_update_dpt[[14]])[!duplicated(data_ori)], colnames(cds_downsampled_cells_ordered_0.8_update_dpt[[14]])[!duplicated(data_ori)])

gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
dp_mst <- minimum.spanning.tree(gp)
minSpanningTree(cds) <- dp_mst

extract_ddrtree_ordering_xj(test@minSpanningTree, cellPairwiseDistances(test), root_cell = which(degree(dp_mst) == 1)[1])
order_res <- extract_ddrtree_ordering_xj(dp_mst, dp = dp, root_cell = which(degree(dp_mst) == 1)[1])

#assign the branches as well as the branch time point:
next_node <<- 0

res <- pq_helper(dp_mst, use_weights=FALSE, root_node=which(degree(dp_mst) == 1)[1])

if(is.null(branch_num))
  branch_num <- sum(degree(dp_mst) > 2) + 1

order_list <- extract_good_branched_ordering(res$subtree, res$root, dp, branch_num, FALSE)

cc_ordering <- order_list$ordering_df
row.names(cc_ordering) <- cc_ordering$sample_name

qplot(dm$DC1, dm$DC2, color = cc_ordering[row.names(data_ori), 'cell_state'])

dp <- DPT_res[1:cell_num, 1:cell_num]
dimnames(dp) <- list(colnames(cds_downsampled_cells_ordered[[36]])[!duplicated(data)], colnames(cds_downsampled_cells_ordered[[36]])[!duplicated(data)])
gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
dp_mst <- minimum.spanning.tree(gp)
cds@minSpanningTree <- dp_mst

cds <- orderCells(cds, num_paths = 2)
plot_cell_trajectory(cds)

cds <- trimTree(cds)
plot_cell_trajectory(cds)

cds_downsampled_cells_ordered_0.8_update_dpt <- lapply(1:length(cds_downsampled_cells_ordered_0.8_update), function(i) {
  message('current id is ', i)
  test <- cds_downsampled_cells_ordered_update[[i]]
  test <- reduceDimension(test, max_components = 2, norm_method = 'none', scaling = F, reduction_method = run_dpt)
  tryCatch({
    test <- orderCells(test, num_paths = 2)
    test <- orderCells(test, num_paths = 2, root_state = Root_state(test))
    return(test)
  }, error = function(e) {
    warning('Your initial method throws numerical errors!')
    test <- orderCells(test)
    return(test)
  })
})

cds_downsampled_cells_ordered_0.8_update_SGL <- lapply(1:length(cds_downsampled_cells_ordered_0.8_update), function(i) {
  message('current id is ', i)
  test <- cds_downsampled_cells_ordered_update[[i]]
  test <- reduceDimension(test, max_components = 2, norm_method = 'none', scaling = F, initial_method = run_dpt, reduction_method = 'SGL-tree')
  tryCatch({
    test <- orderCells(test, num_paths = 2)
    test <- orderCells(test, num_paths = 2, root_state = Root_state(test))
    return(test)
  }, error = function(e) {
    warning('Your initial method throws numerical errors!')
    test <- orderCells(test)
    return(test)
  })
})

cds_downsampled_cells_ordered_0.8_update_SimplePPT <- lapply(1:length(cds_downsampled_cells_ordered_0.8_update), function(i) {
  message('current id is ', i)
  test <- cds_downsampled_cells_ordered_update[[i]]
  test <- reduceDimension(test, max_components = 2, norm_method = 'none', scaling = F, initial_method = run_dpt, reduction_method = 'SimplePPT')
  tryCatch({
    test <- orderCells(test, num_paths = 2)
    test <- orderCells(test, num_paths = 2, root_state = Root_state(test))
    return(test)
  }, error = function(e) {
    warning('Your initial method throws numerical errors!')
    test <- orderCells(test)
    return(test)
  })
})


test <- reduceDimension(test, max_components = 2, norm_method = 'none', auto_param_selection = T, initial_method = run_dpt, maxIter = 1, scaling = F, reduction_method = 'SimplePPT')
test <- orderCells(test)
plot_cell_trajectory(test)

calClusteringMetrics(dpt_cds_downsampled_cells_ordered[[ind]]$branch[order(cell_id)], pData(cds_downsampled_cells_ordered_update[[36]])$State)

test <- trimTree(test)
calClusteringMetrics(dpt_cds_downsampled_cells_ordered[[ind]]$branch[order(cell_id)], pData(test)$State)

pData(absolute_cds)$State[c(1:33, 401:437)] <- 1
pData(absolute_cds)$State[c(34:400)] <- 2
pData(absolute_cds)$State[c(438:800)] <- 3


reducedDim <- reduction_method(FM, ...)
colnames(reducedDim) <- colnames(FM)
reducedDimW(cds) <- as.matrix(reducedDim)
reducedDimA(cds) <- as.matrix(reducedDim)
reducedDimS(cds) <- as.matrix(reducedDim)
reducedDimK(cds) <- as.matrix(reducedDim)
dp <- as.matrix(dist(reducedDim))
cellPairwiseDistances(cds) <- dp
gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
dp_mst <- minimum.spanning.tree(gp)
minSpanningTree(cds) <- dp_mst
cds@dim_reduce_type <- "function_passed"

weights <- W
A <- Matrix::t(solve(weights) %*% Matrix::t(init_ICA$K))
colnames(A) <- colnames(weights)
rownames(A) <- rownames(FM)
S <- weights %*% x_pca
rownames(S) <- colnames(weights)
colnames(S) <- colnames(FM)
reducedDimW(cds) <- as.matrix(W)
reducedDimA(cds) <- as.matrix(A)
reducedDimS(cds) <- as.matrix(S)
reducedDimK(cds) <- as.matrix(init_ICA$K)
adjusted_S <- Matrix::t(reducedDimS(cds))
dp <- as.matrix(dist(adjusted_S))
cellPairwiseDistances(cds) <- dp
gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
dp_mst <- minimum.spanning.tree(gp)
minSpanningTree(cds) <- dp_mst
cds@dim_reduce_type <- "ICA"

software_custom_color_scale <- c("Reference" = "#F3756C",
                                 "Wishbone"="#26B24B",
                                 "Monocle 2 (DPT DCs)" = "#00BCC3",
                                 "Monocle 2 (DPT)" = "#6E95CD",
                                 "DPT" = "#CC71AD",
                                 "Monocle 2" = "#B8A131")

dpt_monocle2_benchmark_res_list <- lapply(cds_downsampled_cells_ordered_0.8_update_dpt, function(x) pairwise_cal_benchmark_res(x, absolute_cds)) #
dpt_progressive_monocle2_benchmark_res_list <- lapply(cds_downsampled_cells_ordered_update_dpt, function(x) pairwise_cal_benchmark_res(x, absolute_cds)) #

dpt_monocle_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(dpt_monocle2_benchmark_res_list, function(x) x$kendall_tau)),
                                      pearson_rho = unlist(lapply(dpt_monocle2_benchmark_res_list, function(x) x$cor)),
                                      rand_ind = unlist(lapply(dpt_monocle2_benchmark_res_list, function(x) x$cluster[1, 1])),
                                      var_inf = unlist(lapply(dpt_monocle2_benchmark_res_list, function(x) x$cluster[2, 1])),
                                      adj_rand = unlist(lapply(dpt_monocle2_benchmark_res_list, function(x) x$cluster[3, 1]))
)

dpt_progressive_monocle_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(dpt_progressive_monocle2_benchmark_res_list, function(x) x$kendall_tau)),
                                                  pearson_rho = unlist(lapply(dpt_progressive_monocle2_benchmark_res_list, function(x) x$cor)),
                                                  rand_ind = unlist(lapply(dpt_progressive_monocle2_benchmark_res_list, function(x) x$cluster[1, 1])),
                                                  var_inf = unlist(lapply(dpt_progressive_monocle2_benchmark_res_list, function(x) x$cluster[2, 1])),
                                                  adj_rand = unlist(lapply(dpt_progressive_monocle2_benchmark_res_list, function(x) x$cluster[3, 1]))
)

diffusion_dpt_monocle2_benchmark_res_list <- lapply(cds_downsampled_cells_ordered_0.8_diffusion_pseudotime, function(x) pairwise_cal_benchmark_res(x, absolute_cds)) #
diffusion_dpt_progressive_monocle2_benchmark_res_list <- lapply(cds_downsampled_cells_ordered_update_diffusion_pseudotime, function(x) pairwise_cal_benchmark_res(x, absolute_cds)) #
diffusion_dpt_monocle_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(diffusion_dpt_monocle2_benchmark_res_list, function(x) x$kendall_tau)),
                                                    pearson_rho = unlist(lapply(diffusion_dpt_monocle2_benchmark_res_list, function(x) x$cor)),
                                                    rand_ind = unlist(lapply(diffusion_dpt_monocle2_benchmark_res_list, function(x) x$cluster[1, 1])),
                                                    var_inf = unlist(lapply(diffusion_dpt_monocle2_benchmark_res_list, function(x) x$cluster[2, 1])),
                                                    adj_rand = unlist(lapply(diffusion_dpt_monocle2_benchmark_res_list, function(x) x$cluster[3, 1]))
)
diffusion_dpt_progressive_monocle_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(diffusion_dpt_progressive_monocle2_benchmark_res_list, function(x) x$kendall_tau)),
                                                      pearson_rho = unlist(lapply(diffusion_dpt_progressive_monocle2_benchmark_res_list, function(x) x$cor)),
                                                      rand_ind = unlist(lapply(diffusion_dpt_progressive_monocle2_benchmark_res_list, function(x) x$cluster[1, 1])),
                                                      var_inf = unlist(lapply(diffusion_dpt_progressive_monocle2_benchmark_res_list, function(x) x$cluster[2, 1])),
                                                      adj_rand = unlist(lapply(diffusion_dpt_progressive_monocle2_benchmark_res_list, function(x) x$cluster[3, 1]))
)

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
downsamling_marker_all_sampling_res_df <- Reduce(rbind , list(dpt_sampling_res_df, monocle_sampling_res_df,  wishbone_sampling_res_df, dpt_monocle_sampling_res_df, diffusion_dpt_monocle_sampling_res_df)) # ICA_sampling_res_df,

downsamling_marker_all_sampling_res_df$Type <- c(rep('DPT', nrow(dpt_sampling_res_df)), rep('Monocle 2',  nrow(monocle_sampling_res_df)), 
                                                 rep('Wishbone', nrow(wishbone_sampling_res_df)), rep('Monocle 2 (DPT DCs)', nrow(dpt_monocle_sampling_res_df)), 
                                                 rep('Monocle 2 (DPT)', nrow(diffusion_dpt_monocle_sampling_res_df)))#,  rep('Monocle1', 10000)
downsamling_marker_all_sampling_res_df$Type <- factor(downsamling_marker_all_sampling_res_df$Type, levels = c('Monocle 2 (DPT)','Monocle 2', 'Monocle 2 (DPT DCs)', "DPT", "Wishbone")) #dpt (non-uniform branch)

pdf(paste(SI_fig_dir, benchmark_type, 'real_simulation_pearson_rho_comparison_robustness.pdf', sep = ''), width = 1, height = 1.5)
qplot(Type, pearson_rho, data = downsamling_marker_all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Pearson's Rho") + monocle_theme_opts()  + xlab('') +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)  + scale_color_manual(values = software_custom_color_scale) + ylab("Pseudotime accuracy \n (Pearson's Rho)")
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'real_simulation_all_kendall_tau_comparison_robustness.pdf', sep = ''), width = 1, height = 1.5)
qplot(Type, kendall.tau, data = downsamling_marker_all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Kendall's tau") + monocle_theme_opts()  + xlab('') +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale) + ylab("Pseudotime accuracy \n (Kendall's Tau)")
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'real_simulation_ArI_comparison_robustness.pdf', sep = ''), width = 1, height = 1.5)
qplot(Type, adj_rand, data = downsamling_marker_all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Adjusted Rand index") + monocle_theme_opts()  + xlab('') +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale) + ylab("Branch accuracy \n (ARI)")
dev.off()

######################################################################################################################################################
# Cole's request to remove DPT + Monocle 2
######################################################################################################################################################
pdf(paste(SI_fig_dir, benchmark_type, 'real_simulation_pearson_rho_comparison_robustness_cole.pdf', sep = ''), width = 1, height = 1)
qplot(Type, pearson_rho, data = subset(downsamling_marker_all_sampling_res_df, Type != 'Monocle 2 (DCs)'), color = Type, geom = 'boxplot') + ylab("Pearson's Rho") + monocle_theme_opts()  + xlab('') +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)  + scale_color_manual(values = software_custom_color_scale) + ylab('Accuracy\n(Pseudotime)')
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'real_simulation_all_kendall_tau_comparison_robustness_cole.pdf', sep = ''), width = 1, height = 1)
qplot(Type, kendall.tau, data = subset(downsamling_marker_all_sampling_res_df, Type != 'Monocle 2 (DCs)'), color = Type, geom = 'boxplot') + ylab("Kendall's tau") + monocle_theme_opts()  + xlab('') +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale) + ylab("Pseudotime accuracy \n(Kendall's Tau)")
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'real_simulation_ArI_comparison_robustness_cole.pdf', sep = ''), width = 1, height = 1)
qplot(Type, adj_rand, data = subset(downsamling_marker_all_sampling_res_df, Type != 'Monocle 2 (DCs)'), color = Type, geom = 'boxplot') + ylab("Adjusted Rand index") + monocle_theme_opts()  + xlab('') +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale) + ylab("Accuracy\n(Branch)")
dev.off()

######################################################################################################################################################
# Progressive downsampling 
######################################################################################################################################################

progressive_all_valid_cell_sampling_res_df <- Reduce(rbind , list(progressive_dpt_sampling_res_df, progressive_monocle_sampling_res_df,  progressive_wishbone_sampling_res_df, dpt_progressive_monocle_sampling_res_df, diffusion_dpt_progressive_monocle_sampling_res_df)) # ICA_sampling_res_df,

progressive_all_valid_cell_sampling_res_df$proportion <- c(rep(downsampled_proportions, 2), names(cds_downsampled_cells_ordered)[unique(fraction_wishbone_res$run)], downsampled_proportions, downsampled_proportions)
progressive_all_valid_cell_sampling_res_df$Type <- c(rep('DPT', 36), rep('Monocle 2', 36), rep('Wishbone', length(unique(fraction_wishbone_res$run))), rep('Monocle 2 (DPT DCs)', 36), rep('Monocle 2 (DPT)', 36))

progressive_all_valid_cell_sampling_res_df$Type <- factor(progressive_all_valid_cell_sampling_res_df$Type, levels = c('Monocle 2 (DPT)','Monocle 2', 'Monocle 2 (DPT DCs)', "DPT", "Wishbone"))
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
ggplot(aes(proportion, mean_kendall.tau), data = subset(progressive_process_cell_sampling_res_df, Type != c('Monocle 2 (DCs)'))) +
  geom_point(aes(color = Type), position=position_dodge(width=0.9), size = 0.5) + facet_wrap(~Type) +
  geom_errorbar(aes(color = Type, ymax = mean_kendall.tau + sd_kendall.tau, ymin=mean_kendall.tau - sd_kendall.tau), position=position_dodge(width=0.1), size = 0.5) +
  xlab('Proportion of original cells') + ylab("Kendall's Tau") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

######################################################################################################################################################
# Cole's request to remove DPT + Monocle 2
######################################################################################################################################################
progressive_process_cell_sampling_res_df$proportion <- as.numeric(as.character(progressive_process_cell_sampling_res_df$proportion))
pdf(paste(SI_fig_dir, benchmark_type, 'rand_index_real_simulation_comparison_cell_downsampling2_cole.pdf', sep = ''), width = 3, height = 1)
ggplot(aes(proportion, mean_adj_rand), data = subset(progressive_process_cell_sampling_res_df, Type != c('Monocle 2 (DCs)'))) +
  geom_point(aes(color = Type), position=position_dodge(width=0.1), size = 0.5) + facet_wrap(~Type, nrow = 1) +
  geom_errorbar(aes(color = Type, ymax = mean_adj_rand + sd_adj_rand, ymin=mean_adj_rand - sd_adj_rand), position=position_dodge(width=0.9), size = 0.5) +
  xlab('Proportion of original cells') + ylab("Adjusted Rand index") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale) +
  scale_x_continuous(limits=c(0.1, 1), breaks = seq(0.1, 1, by = 0.2))  + ylab('Accuracy\n(Branch)')
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_real_simulation_comparison_cell_downsampling2_cole.pdf', sep = ''), width = 3, height = 1)
ggplot(aes(proportion, mean_pearson_rho), data = subset(progressive_process_cell_sampling_res_df, Type != c('Monocle 2 (DCs)'))) +
  geom_point(aes(color = Type), position=position_dodge(width=0.9), size = 0.5) + facet_wrap(~Type, nrow = 1) +
  geom_errorbar(aes(color = Type, ymax = mean_pearson_rho + sd_pearson_rho, ymin=mean_pearson_rho - sd_pearson_rho), position=position_dodge(width=0.1), size = 0.5) +
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale) +
  scale_x_continuous(limits=c(0.1, 1), breaks = seq(0.1, 1, by = 0.2))  + ylab('Accuracy\n(Pseudotime)')
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'kendall_tau_real_simulation_comparison_cell_downsampling2_cole.pdf', sep = ''), width = 3, height = 1)
ggplot(aes(proportion, mean_kendall.tau), data = subset(progressive_process_cell_sampling_res_df, Type != c('Monocle 2 (DCs)'))) +
  geom_point(aes(color = Type), position=position_dodge(width=0.9), size = 0.5) + facet_wrap(~Type, nrow = 1) +
  geom_errorbar(aes(color = Type, ymax = mean_kendall.tau + sd_kendall.tau, ymin=mean_kendall.tau - sd_kendall.tau), position=position_dodge(width=0.1), size = 0.5) +
  xlab('Proportion of original cells') + ylab("Kendall's Tau") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale) +
  scale_x_continuous(limits=c(0.1, 1), breaks = seq(0.1, 1, by = 0.2)) + ylab('Accuracy\n(Pseudotime 2)')
dev.off()

