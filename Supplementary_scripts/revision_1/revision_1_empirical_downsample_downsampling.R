# show the robustness based on repeated downsampling:
library(monocle)

software_custom_color_scale <- c("reference" = "#F3756C",
                                 "Wishbone"="#26B24B",
                                 "Monocle 1" = "#00BCC3",
                                 "Slicer" = "#6E95CD",
                                 "DPT" = "#CC71AD",
                                 "Monocle 2" = "#B8A131")
# + scale_color_manual(values=software_custom_color_scale)
software_levels <- c('Monocle 2', 'Monocle 1', 'DPT', 'Wishbone', 'Slicer')

# load('./RData/analysis_score_ordering_MAR_seq.RData') #get the score ordering for all cells:
# all_cell_score_df <- score_df
# load('./RData/marseq_downsamling_empirical_ordering.RData')
# 
# downsampling_score_df <- score_df
# save(file = './RData/downsampling_score_df', downsampling_score_df)
# 
# # load('./RData/MAR_seq_data_robustness_dpt_slicer_wishbone.RData')
# # load('./RData/", benchmark_type, "ICA_downsampling.RData')

ICA_cds_downsampled_cells_ordered_0.8
ICA_cds_downsampled_cells_ordered

#ICA
downsampling_score_df_cgmp <- subset(downsampling_score_df, lineage %in% c('CMP', 'GMP'))
downsampling_score_df_cmep <- subset(downsampling_score_df, lineage %in% c('CMP', 'MEP'))

ICA_marker_reference_list <- lapply(ICA_cds_downsampled_cells_ordered_0.8, function(cds) {
  if(!is.na(cds)) {
    overlap_cgmp <- intersect(colnames(cds), row.names(downsampling_score_df_cgmp))
    overlap_cmep <- intersect(colnames(cds), row.names(downsampling_score_df_cmep))
    overalap_all <- intersect(colnames(cds), row.names(downsampling_score_df))
    print('pass here')

    monocle1_kendall_tau <- c(cor(rank(downsampling_score_df_cgmp[overlap_cgmp, "naive_pseudotime"]), rank(pData(cds[, overlap_cgmp])$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs'),
      cor(rank(downsampling_score_df_cmep[overlap_cmep, "naive_pseudotime"]), rank(pData(cds[, overlap_cmep])$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs'))
    monocle1_pearson_rho <- c(cor(rank(downsampling_score_df_cgmp[overlap_cgmp, "naive_pseudotime"]), rank(pData(cds[, overlap_cgmp])$Pseudotime)),
                              cor(rank(downsampling_score_df_cmep[overlap_cmep, "naive_pseudotime"]), rank(pData(cds[, overlap_cmep])$Pseudotime)))
    ARI <- calClusteringMetrics(downsampling_score_df[overalap_all, "lineage"], pData(cds)[overalap_all, "State"])
    cell_type_ARI <- calClusteringMetrics(pData(absolute_cds)[overalap_all, "cell_type"], pData(cds)[overalap_all, "State"])
    all_kendall <- cor(rank(downsampling_score_df[overalap_all, "naive_pseudotime"]), rank(pData(cds[, overalap_all])$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
    all_pearson <- cor(downsampling_score_df[overalap_all, "naive_pseudotime"], pData(cds[, overalap_all])$Pseudotime)

    data.frame(Person = mean(abs(monocle1_pearson_rho)), Kendall = mean(abs(monocle1_kendall_tau)), ARI = ARI[3, 1], all_kendall = all_kendall, all_pearson = all_pearson, cell_type_ARI = cell_type_ARI[3, 1])
  }
})

ICA_marker_reference_res <- do.call(rbind, ICA_marker_reference_list)

#monocle 2, dpt, wishbone:
all_cell_score_df_cgmp <- subset(all_cell_score_df, lineage %in% c('CMP', 'GMP'))
all_cell_score_df_cmep <- subset(all_cell_score_df, lineage %in% c('CMP', 'MEP'))

monocle2_marker_reference_list <- lapply(cds_downsampled_cells_ordered_0.8, function(cds) {
  if(!is.na(cds)) {
    overlap_cgmp <- intersect(colnames(cds), row.names(all_cell_score_df_cgmp))
    overlap_cmep <- intersect(colnames(cds), row.names(all_cell_score_df_cmep))
    overalap_all <- intersect(colnames(cds), row.names(all_cell_score_df))

    print('pass here')
    monocle1_kendall_tau <- c(cor(rank(all_cell_score_df_cgmp[overlap_cgmp, "naive_pseudotime"]), rank(pData(cds[, overlap_cgmp])$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs'),
                              cor(rank(all_cell_score_df_cmep[overlap_cmep, "naive_pseudotime"]), rank(pData(cds[, overlap_cmep])$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs'))
    monocle1_pearson_rho <- c(cor(rank(all_cell_score_df_cgmp[overlap_cgmp, "naive_pseudotime"]), rank(pData(cds[, overlap_cgmp])$Pseudotime)),
                              cor(rank(all_cell_score_df_cmep[overlap_cmep, "naive_pseudotime"]), rank(pData(cds[, overlap_cmep])$Pseudotime)))
    ARI <- calClusteringMetrics(all_cell_score_df[overalap_all, "lineage"], pData(cds)[overalap_all, "State"])
    cell_type_ARI <- calClusteringMetrics(pData(absolute_cds)[overalap_all, "cell_type"], pData(cds)[overalap_all, "State"])

    all_kendall <- cor(rank(all_cell_score_df[overalap_all, "naive_pseudotime"]), rank(pData(cds[, overalap_all])$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
    all_pearson <- cor(all_cell_score_df[overalap_all, "naive_pseudotime"], pData(cds[, overalap_all])$Pseudotime)

    data.frame(Person = mean(abs(monocle1_pearson_rho)), Kendall = mean(abs(monocle1_kendall_tau)), ARI = ARI[3, 1], all_kendall = all_kendall, all_pearson = all_pearson, cell_type_ARI = cell_type_ARI[3, 1])
  }
})

monocle2_marker_reference_res <- do.call(rbind, monocle2_marker_reference_list)

dpt_marker_reference_list <- lapply(dpt_cds_downsampled_cells_ordered_0.8, function(dpt_res) {
  if(!is.na(dpt_res)) {
    overlap_cgmp <- intersect(row.names(all_cell_score_df_cgmp), names(dpt_res$pt))
    overlap_cmep <- intersect(row.names(all_cell_score_df_cmep), names(dpt_res$pt))
    overalap_all <- intersect(names(dpt_res$pt), row.names(all_cell_score_df))
    print('pass here')
    monocle1_kendall_tau <- c(cor(rank(all_cell_score_df_cgmp[overlap_cgmp, "naive_pseudotime"]), rank(dpt_res$pt[overlap_cgmp]), method = 'kendall', use = 'pairwise.complete.obs'),
                              cor(rank(all_cell_score_df_cmep[overlap_cmep, "naive_pseudotime"]), rank(dpt_res$pt[overlap_cmep]), method = 'kendall', use = 'pairwise.complete.obs'))
    monocle1_pearson_rho <- c(cor(rank(all_cell_score_df_cgmp[overlap_cgmp, "naive_pseudotime"]), rank(dpt_res$pt[overlap_cgmp])),
                              cor(rank(all_cell_score_df_cmep[overlap_cmep, "naive_pseudotime"]), rank(dpt_res$pt[overlap_cmep])))

    ARI <- calClusteringMetrics(all_cell_score_df[overalap_all, "lineage"], dpt_res$branch[overalap_all, 'Branch1'])
    cell_type_ARI <- calClusteringMetrics(pData(absolute_cds)[overalap_all, "cell_type"], dpt_res$branch[overalap_all, 'Branch1'])
    all_kendall <- cor(rank(all_cell_score_df[overalap_all, "naive_pseudotime"]), rank(dpt_res$pt[overalap_all]), method = 'kendall', use = 'pairwise.complete.obs')
    all_pearson <- cor(all_cell_score_df[overalap_all, "naive_pseudotime"], dpt_res$pt[overalap_all])

    data.frame(Person = mean(abs(monocle1_pearson_rho)), Kendall = mean(abs(monocle1_kendall_tau)), ARI = ARI[3, 1], all_kendall = all_kendall, all_pearson = all_pearson, cell_type_ARI = cell_type_ARI[3, 1])
  }
})

dpt_marker_reference_res <- do.call(rbind, dpt_marker_reference_list)

wishbone_marker_reference_list <- lapply(1:25, function(ind) {
  subset_wishbone_res <- subset(wishbone_res, run == ind)
  row.names(subset_wishbone_res) <- subset_wishbone_res[, 1]
  if(nrow(subset_wishbone_res) > 0) {
    overlap_cgmp <- intersect(row.names(all_cell_score_df_cgmp), subset_wishbone_res[, 1])
    overlap_cmep <- intersect(row.names(all_cell_score_df_cmep), subset_wishbone_res[, 1])
    overalap_all <- intersect(subset_wishbone_res[, 1], row.names(all_cell_score_df))
    print('pass here')
    monocle1_kendall_tau <- c(cor(rank(all_cell_score_df_cgmp[overlap_cgmp, "naive_pseudotime"]), rank(subset_wishbone_res[overlap_cgmp, 'trajectory']), method = 'kendall', use = 'pairwise.complete.obs'),
                              cor(rank(all_cell_score_df_cmep[overlap_cmep, "naive_pseudotime"]), rank(subset_wishbone_res[overlap_cmep, 'trajectory']), method = 'kendall', use = 'pairwise.complete.obs'))
    monocle1_pearson_rho <- c(cor(rank(all_cell_score_df_cgmp[overlap_cgmp, "naive_pseudotime"]), rank(subset_wishbone_res[overlap_cgmp, 'trajectory'])),
                              cor(rank(all_cell_score_df_cmep[overlap_cmep, "naive_pseudotime"]), rank(subset_wishbone_res[overlap_cmep, 'trajectory'])))

    ARI <- calClusteringMetrics(all_cell_score_df[overalap_all, "lineage"], subset_wishbone_res[overalap_all, "branch"])
    cell_type_ARI <- calClusteringMetrics(pData(absolute_cds)[overalap_all, "cell_type"], subset_wishbone_res[overalap_all, "branch"])
    all_kendall <- cor(rank(all_cell_score_df[overalap_all, "naive_pseudotime"]), rank(subset_wishbone_res[overalap_all, 'trajectory']), method = 'kendall', use = 'pairwise.complete.obs')
    all_pearson <- cor(all_cell_score_df[overalap_all, "naive_pseudotime"], subset_wishbone_res[overalap_all, 'trajectory'])

    data.frame(Person = mean(abs(monocle1_pearson_rho)), Kendall = mean(abs(monocle1_kendall_tau)), ARI = ARI[3, 1], all_kendall = all_kendall, all_pearson = all_pearson, cell_type_ARI = cell_type_ARI[3, 1])
    }
})

wishbone_marker_reference_res <- do.call(rbind, wishbone_marker_reference_list)

# make the histogram:
###robustness downsampling: color by both of the dpt and ddrtree result together:
downsamling_marker_all_sampling_res_df <- Reduce(rbind , list(dpt_marker_reference_res, monocle2_marker_reference_res,  wishbone_marker_reference_res, ICA_marker_reference_res)) # ICA_sampling_res_df,

downsamling_marker_all_sampling_res_df$Type <- c(rep('DPT', nrow(dpt_marker_reference_res)), rep('Monocle 2',  nrow(monocle2_marker_reference_res)), rep('Wishbone', nrow(wishbone_marker_reference_res)), rep('Monocle 1', nrow(ICA_marker_reference_res)))#,  rep('Monocle1', 10000)
downsamling_marker_all_sampling_res_df$Type <- factor(downsamling_marker_all_sampling_res_df$Type, levels = c('Monocle 2', 'Monocle 1', "DPT", "Wishbone")) #dpt (non-uniform branch)

pdf(paste(SI_fig_dir, benchmark_type, 'downsampling_marker_pearson_rho_comparison_robustness.pdf', sep = ''), width = 1, height = 1)
qplot(Type, Person, data = downsamling_marker_all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Accuracy\n(Pseudotime)") + monocle_theme_opts()  + xlab('') +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'downsampling_marker_all_pearson_rho_comparison_robustness.pdf', sep = ''), width = 1, height = 1)
qplot(Type, all_pearson, data = downsamling_marker_all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Accuracy\n(Pseudotime)") + monocle_theme_opts()  + xlab('') +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'downsampling_marker_all_kendall_tau_comparison_robustness.pdf', sep = ''), width = 1, height = 1)
qplot(Type, all_kendall, data = downsamling_marker_all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Accuracy\n(Pseudotime)") + monocle_theme_opts()  + xlab('') +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'downsampling_marker_ARI_comparison_robustness.pdf', sep = ''), width = 1, height = 1)
qplot(Type, ARI, data = downsamling_marker_all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Accuracy\n(Branch)") + monocle_theme_opts()  + xlab('') +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'downsampling_marker_cell_type_ARI_comparison_robustness.pdf', sep = ''), width = 1, height = 1)
qplot(Type, cell_type_ARI, data = downsamling_marker_all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Accuracy\n(Branch)") + monocle_theme_opts()  + xlab('') +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

################################################################################################################################################################################
# do the same analysis for the progressive downsmpaling 
################################################################################################################################################################################
# show the robustness based on repeated downsampling:

#ICA
progressive_ICA_marker_reference_list <- lapply(ICA_cds_downsampled_cells_ordered, function(cds) {
  if(!is.na(cds)) {
    overlap_cgmp <- intersect(colnames(cds), row.names(downsampling_score_df_cgmp))
    overlap_cmep <- intersect(colnames(cds), row.names(downsampling_score_df_cmep))
    overalap_all <- intersect(colnames(cds), row.names(downsampling_score_df))
    print('pass here')
    
    monocle1_kendall_tau <- c(cor(rank(downsampling_score_df_cgmp[overlap_cgmp, "naive_pseudotime"]), rank(pData(cds[, overlap_cgmp])$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs'),
                              cor(rank(downsampling_score_df_cmep[overlap_cmep, "naive_pseudotime"]), rank(pData(cds[, overlap_cmep])$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs'))
    monocle1_pearson_rho <- c(cor(rank(downsampling_score_df_cgmp[overlap_cgmp, "naive_pseudotime"]), rank(pData(cds[, overlap_cgmp])$Pseudotime)),
                              cor(rank(downsampling_score_df_cmep[overlap_cmep, "naive_pseudotime"]), rank(pData(cds[, overlap_cmep])$Pseudotime)))
    ARI <- calClusteringMetrics(downsampling_score_df[overalap_all, "lineage"], pData(cds)[overalap_all, "State"])
    cell_type_ARI <- calClusteringMetrics(pData(absolute_cds)[overalap_all, "cell_type"], pData(cds)[overalap_all, "State"])
    all_kendall <- cor(rank(downsampling_score_df[overalap_all, "naive_pseudotime"]), rank(pData(cds[, overalap_all])$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
    all_pearson <- cor(downsampling_score_df[overalap_all, "naive_pseudotime"], pData(cds[, overalap_all])$Pseudotime)
    
    data.frame(Person = mean(abs(monocle1_pearson_rho)), Kendall = mean(abs(monocle1_kendall_tau)), ARI = ARI[3, 1], all_kendall = all_kendall, all_pearson = all_pearson, cell_type_ARI = cell_type_ARI[3, 1])
  }
})

progressive_ICA_marker_reference_res <- do.call(rbind, progressive_ICA_marker_reference_list)

#monocle 2, dpt, wishbone:
progressive_monocle2_marker_reference_list <- lapply(cds_downsampled_cells_ordered, function(cds) {
  if(!is.na(cds)) {
    overlap_cgmp <- intersect(colnames(cds), row.names(all_cell_score_df_cgmp))
    overlap_cmep <- intersect(colnames(cds), row.names(all_cell_score_df_cmep))
    overalap_all <- intersect(colnames(cds), row.names(all_cell_score_df))
    
    print('pass here')
    monocle1_kendall_tau <- c(cor(rank(all_cell_score_df_cgmp[overlap_cgmp, "naive_pseudotime"]), rank(pData(cds[, overlap_cgmp])$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs'),
                              cor(rank(all_cell_score_df_cmep[overlap_cmep, "naive_pseudotime"]), rank(pData(cds[, overlap_cmep])$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs'))
    monocle1_pearson_rho <- c(cor(rank(all_cell_score_df_cgmp[overlap_cgmp, "naive_pseudotime"]), rank(pData(cds[, overlap_cgmp])$Pseudotime)),
                              cor(rank(all_cell_score_df_cmep[overlap_cmep, "naive_pseudotime"]), rank(pData(cds[, overlap_cmep])$Pseudotime)))
    ARI <- calClusteringMetrics(all_cell_score_df[overalap_all, "lineage"], pData(cds)[overalap_all, "State"])
    cell_type_ARI <- calClusteringMetrics(pData(absolute_cds)[overalap_all, "cell_type"], pData(cds)[overalap_all, "State"])
    
    all_kendall <- cor(rank(all_cell_score_df[overalap_all, "naive_pseudotime"]), rank(pData(cds[, overalap_all])$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
    all_pearson <- cor(all_cell_score_df[overalap_all, "naive_pseudotime"], pData(cds[, overalap_all])$Pseudotime)
    
    data.frame(Person = mean(abs(monocle1_pearson_rho)), Kendall = mean(abs(monocle1_kendall_tau)), ARI = ARI[3, 1], all_kendall = all_kendall, all_pearson = all_pearson, cell_type_ARI = cell_type_ARI[3, 1])
  }
})

progressive_monocle2_marker_reference_res <- do.call(rbind, progressive_monocle2_marker_reference_list)

progressive_dpt_marker_reference_list <- lapply(dpt_cds_downsampled_cells_ordered, function(dpt_res) {
  if(!is.na(dpt_res)) {
    overlap_cgmp <- intersect(row.names(all_cell_score_df_cgmp), names(dpt_res$pt))
    overlap_cmep <- intersect(row.names(all_cell_score_df_cmep), names(dpt_res$pt))
    overalap_all <- intersect(names(dpt_res$pt), row.names(all_cell_score_df))
    print('pass here')
    monocle1_kendall_tau <- c(cor(rank(all_cell_score_df_cgmp[overlap_cgmp, "naive_pseudotime"]), rank(dpt_res$pt[overlap_cgmp]), method = 'kendall', use = 'pairwise.complete.obs'),
                              cor(rank(all_cell_score_df_cmep[overlap_cmep, "naive_pseudotime"]), rank(dpt_res$pt[overlap_cmep]), method = 'kendall', use = 'pairwise.complete.obs'))
    monocle1_pearson_rho <- c(cor(rank(all_cell_score_df_cgmp[overlap_cgmp, "naive_pseudotime"]), rank(dpt_res$pt[overlap_cgmp])),
                              cor(rank(all_cell_score_df_cmep[overlap_cmep, "naive_pseudotime"]), rank(dpt_res$pt[overlap_cmep])))
    
    ARI <- calClusteringMetrics(all_cell_score_df[overalap_all, "lineage"], dpt_res$branch[overalap_all, 'Branch1'])
    cell_type_ARI <- calClusteringMetrics(pData(absolute_cds)[overalap_all, "cell_type"], dpt_res$branch[overalap_all, 'Branch1'])
    all_kendall <- cor(rank(all_cell_score_df[overalap_all, "naive_pseudotime"]), rank(dpt_res$pt[overalap_all]), method = 'kendall', use = 'pairwise.complete.obs')
    all_pearson <- cor(all_cell_score_df[overalap_all, "naive_pseudotime"], dpt_res$pt[overalap_all])
    
    data.frame(Person = mean(abs(monocle1_pearson_rho)), Kendall = mean(abs(monocle1_kendall_tau)), ARI = ARI[3, 1], all_kendall = all_kendall, all_pearson = all_pearson, cell_type_ARI = cell_type_ARI[3, 1])
  }
})

progressive_dpt_marker_reference_res <- do.call(rbind, progressive_dpt_marker_reference_list)

# test DPT on why downsampling = 1, the result of ARI value is close to 0. 
test_DPT <- dpt_cds_downsampled_cells_ordered$`1`
qplot(test_DPT$dm@eigenvectors[, 1], test_DPT$dm@eigenvectors[, 2], color = test_DPT$branch[, 'Branch1'])
test_DPT <- dpt_cds_downsampled_cells_ordered$`1`

pdf('./Figures/test.pdf')
qplot(test_DPT$dm@eigenvectors[, 1], test_DPT$dm@eigenvectors[, 2], color = test_DPT$branch[, 'Branch1'])
dev.off()

pdf('./Figures/test.pdf')
qplot(test_DPT$dm@eigenvectors[, 1], test_DPT$dm@eigenvectors[, 2], color = all_cell_score_df[, 'lineage'])
dev.off()

# monocle results: 
pdf('./Figures/test2.pdf')
plot_cell_trajectory(cds_downsampled_cells_ordered$`1`)
dev.off()

test <- cds_downsampled_cells_ordered$`1`
pData(test)$marker_branch <- all_cell_score_df[, 'lineage']
pdf('./Figures/test2.pdf')
plot_cell_trajectory(test, color_by = 'marker_branch')
dev.off()

progressive_wishbone_marker_reference_list <- lapply(unique(fraction_wishbone_res$run), function(ind) {
  message('ind is ', ind)
  subset_wishbone_res <- subset(fraction_wishbone_res, run == ind)
  row.names(subset_wishbone_res) <- subset_wishbone_res[, 1]
  if(nrow(subset_wishbone_res) > 0) {
    overlap_cgmp <- intersect(row.names(all_cell_score_df_cgmp), subset_wishbone_res[, 1])
    overlap_cmep <- intersect(row.names(all_cell_score_df_cmep), subset_wishbone_res[, 1])
    overalap_all <- intersect(subset_wishbone_res[, 1], row.names(all_cell_score_df))
    print('pass here')
    monocle1_kendall_tau <- c(cor(rank(all_cell_score_df_cgmp[overlap_cgmp, "naive_pseudotime"]), rank(subset_wishbone_res[overlap_cgmp, 'trajectory']), method = 'kendall', use = 'pairwise.complete.obs'),
                              cor(rank(all_cell_score_df_cmep[overlap_cmep, "naive_pseudotime"]), rank(subset_wishbone_res[overlap_cmep, 'trajectory']), method = 'kendall', use = 'pairwise.complete.obs'))
    monocle1_pearson_rho <- c(cor(rank(all_cell_score_df_cgmp[overlap_cgmp, "naive_pseudotime"]), rank(subset_wishbone_res[overlap_cgmp, 'trajectory'])),
                              cor(rank(all_cell_score_df_cmep[overlap_cmep, "naive_pseudotime"]), rank(subset_wishbone_res[overlap_cmep, 'trajectory'])))
    
    ARI <- calClusteringMetrics(all_cell_score_df[overalap_all, "lineage"], subset_wishbone_res[overalap_all, "branch"])
    cell_type_ARI <- calClusteringMetrics(pData(absolute_cds)[overalap_all, "cell_type"], subset_wishbone_res[overalap_all, "branch"])
    all_kendall <- cor(rank(all_cell_score_df[overalap_all, "naive_pseudotime"]), rank(subset_wishbone_res[overalap_all, 'trajectory']), method = 'kendall', use = 'pairwise.complete.obs')
    all_pearson <- cor(all_cell_score_df[overalap_all, "naive_pseudotime"], subset_wishbone_res[overalap_all, 'trajectory'])
    
    data.frame(Person = mean(abs(monocle1_pearson_rho)), Kendall = mean(abs(monocle1_kendall_tau)), ARI = ARI[3, 1], all_kendall = all_kendall, all_pearson = all_pearson, cell_type_ARI = cell_type_ARI[3, 1])
  }
})

progressive_wishbone_marker_reference_res <- do.call(rbind, progressive_wishbone_marker_reference_list)

# make the histogram:
# ###cell downsampling: color by both of the dpt and ddrtree result together: 
progressive_all_valid_cell_sampling_res_df <- Reduce(rbind , list(progressive_dpt_marker_reference_res, progressive_monocle2_marker_reference_res,  progressive_wishbone_marker_reference_res, progressive_ICA_marker_reference_res)) # ICA_sampling_res_df,

progressive_all_valid_cell_sampling_res_df$proportion <- c(rep(downsampled_proportions, 2), names(cds_downsampled_cells_ordered)[unique(fraction_wishbone_res$run)], downsampled_proportions)
progressive_all_valid_cell_sampling_res_df$Type <- c(rep('DPT', 36), rep('Monocle 2', 36), rep('Wishbone', length(unique(fraction_wishbone_res$run))), rep('Monocle 1', 36))

progressive_all_valid_cell_sampling_res_df$Type <- factor(progressive_all_valid_cell_sampling_res_df$Type, levels = c('Monocle 2', 'Monocle 1', "DPT", "Wishbone")) 
progressive_all_valid_cell_sampling_res_df$se <- 0.1 

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_marker_reference_rho_comparison_cell_downsampling.pdf', sep = ''), width = 2.5, height = 2) 
qplot(proportion, abs(progressive_all_valid_cell_sampling_res_df$pearson_rho), data = progressive_all_valid_cell_sampling_res_df, color = Type, size = 1, geom = 'boxplot') +
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + nm_theme() + scale_size(range = c(0.1, 1)) + monocle_theme_opts() + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'rand_index_comparison_cell_downsampling.pdf', sep = ''), width = 2.5, height = 2) 
qplot(proportion, abs(progressive_all_valid_cell_sampling_res_df$adj_rand), data = progressive_all_valid_cell_sampling_res_df, color = Type, size = 1, geom = 'boxplot') + 
  xlab('Proportion of original cells') + ylab("Adjusted rand index") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

progressive_process_cell_sampling_res_df <- ddply(progressive_all_valid_cell_sampling_res_df, .(Type, proportion), summarize, 
                                      mean_kendall.tau = mean(abs(all_kendall), na.rm = T), 
                                      sd_kendall.tau = sd(abs(all_kendall), na.rm = T), 
                                      
                                      mean_pearson_rho = mean(abs(all_pearson), na.rm = T), 
                                      sd_pearson_rho = sd(abs(all_pearson), na.rm = T), 
                                      
                                      mean_adj_rand = mean(abs(ARI), na.rm = T),
                                      sd_adj_rand = sd(abs(ARI), na.rm = T),
                                      
                                      mean_cell_type_adj_rand = mean(abs(cell_type_ARI), na.rm = T),
                                      sd_cell_type_adj_rand = sd(abs(cell_type_ARI), na.rm = T),
                                      
                                      se = mean(se))
limits <- aes(ymax = mean_adj_rand + sd_adj_rand, ymin=mean_adj_rand - sd_adj_rand)
progressive_process_cell_sampling_res_df$proportion <- as.numeric(as.character(progressive_process_cell_sampling_res_df$proportion))

pdf(paste(SI_fig_dir, benchmark_type, 'rand_index_marker_reference_comparison_cell_downsampling2.pdf', sep = ''), width = 2.5, height = 2) 
ggplot(aes(proportion, mean_adj_rand), data = progressive_process_cell_sampling_res_df) + 
  geom_point(aes(color = Type), position=position_dodge(width=0.1), size = 0.5) + facet_wrap(~Type) + 
  geom_errorbar(aes(color = Type, ymax = mean_adj_rand + sd_adj_rand, ymin=mean_adj_rand - sd_adj_rand), position=position_dodge(width=0.9), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("ARI (Marker ordering))") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_marker_reference_comparison_cell_downsampling2.pdf', sep = ''), width = 2.5, height = 2) 
ggplot(aes(proportion, mean_pearson_rho), data = progressive_process_cell_sampling_res_df) + 
  geom_point(aes(color = Type), position=position_dodge(width=0.9), size = 0.5) + facet_wrap(~Type) + 
  geom_errorbar(aes(color = Type, ymax = mean_pearson_rho + sd_pearson_rho, ymin=mean_pearson_rho - sd_pearson_rho), position=position_dodge(width=0.1), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("Accuracy (Branch)") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale) + 
  scale_x_continuous(limits=c(0.1, 1), breaks = seq(0.1, 1, by = 0.2)) 
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'kendall_tau_marker_reference_comparison_cell_downsampling2.pdf', sep = ''), width = 2.5, height = 2) 
ggplot(aes(proportion, mean_kendall.tau), data = progressive_process_cell_sampling_res_df) + 
  geom_point(aes(color = Type), position=position_dodge(width=0.9), size = 0.5) + facet_wrap(~Type) + 
  geom_errorbar(aes(color = Type, ymax = mean_kendall.tau + sd_kendall.tau, ymin=mean_kendall.tau - sd_kendall.tau), position=position_dodge(width=0.1), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("Kendall's Tau") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale) + 
  scale_x_continuous(limits=c(0.1, 1), breaks = seq(0.1, 1, by = 0.2)) 
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'cell_type_ari_marker_reference_comparison_cell_downsampling2.pdf', sep = ''), width = 2.5, height = 2) 
ggplot(aes(proportion, mean_cell_type_adj_rand), data = progressive_process_cell_sampling_res_df) + 
  geom_point(aes(color = Type), position=position_dodge(width=0.9), size = 0.5) + facet_wrap(~Type) + 
  geom_errorbar(aes(color = Type, ymax = mean_cell_type_adj_rand + sd_cell_type_adj_rand, ymin=mean_cell_type_adj_rand - sd_cell_type_adj_rand), position=position_dodge(width=0.1), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("Accuracy (Branch)") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale) + 
  scale_x_continuous(limits=c(0.1, 1), breaks = seq(0.1, 1, by = 0.2)) 
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_marker_reference_comparison_cell_downsampling2_helper.pdf', sep = '')) 
ggplot(aes(proportion, mean_pearson_rho), data = progressive_process_cell_sampling_res_df) + 
  geom_point(aes(color = Type), position=position_dodge(width=0.1), size = 0.5) + 
  geom_errorbar(aes(color = Type, ymax = mean_pearson_rho + sd_pearson_rho, ymin=mean_pearson_rho - sd_pearson_rho), position=position_dodge(width=0.1), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + scale_size(range = c(0.1, 1)) + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

####################################################################################################################################################################################
# Cole layout: 
####################################################################################################################################################################################
progressive_process_cell_sampling_res_df$proportion <- as.numeric(as.character(progressive_process_cell_sampling_res_df$proportion))
pdf(paste(SI_fig_dir, benchmark_type, 'rand_index_marker_reference_comparison_cell_downsampling2_cole.pdf', sep = ''), width = 5, height = 1) 
ggplot(aes(proportion, mean_adj_rand), data = progressive_process_cell_sampling_res_df) + 
  geom_point(aes(color = Type), position=position_dodge(width=0.1), size = 0.5) + facet_wrap(~Type, nrow = 1) + 
  geom_errorbar(aes(color = Type, ymax = mean_adj_rand + sd_adj_rand, ymin=mean_adj_rand - sd_adj_rand), position=position_dodge(width=0.9), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("ARI (Marker ordering))") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()
pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_marker_reference_comparison_cell_downsampling2_cole.pdf', sep = ''), width = 5, height = 1) 
ggplot(aes(proportion, mean_pearson_rho), data = progressive_process_cell_sampling_res_df) + 
  geom_point(aes(color = Type), position=position_dodge(width=0.9), size = 0.5) + facet_wrap(~Type, nrow = 1) + 
  geom_errorbar(aes(color = Type, ymax = mean_pearson_rho + sd_pearson_rho, ymin=mean_pearson_rho - sd_pearson_rho), position=position_dodge(width=0.1), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("Accuracy (Pseudotime)") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale) + 
  scale_x_continuous(limits=c(0.1, 1), breaks = seq(0.1, 1, by = 0.2)) 
dev.off()
pdf(paste(SI_fig_dir, benchmark_type, 'kendall_tau_marker_reference_comparison_cell_downsampling2_cole.pdf', sep = ''), width = 5, height = 1) 
ggplot(aes(proportion, mean_kendall.tau), data = progressive_process_cell_sampling_res_df) + 
  geom_point(aes(color = Type), position=position_dodge(width=0.9), size = 0.5) + facet_wrap(~Type, nrow = 1) + 
  geom_errorbar(aes(color = Type, ymax = mean_kendall.tau + sd_kendall.tau, ymin=mean_kendall.tau - sd_kendall.tau), position=position_dodge(width=0.1), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("Kendall's Tau") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale) + 
  scale_x_continuous(limits=c(0.1, 1), breaks = seq(0.1, 1, by = 0.2)) 
dev.off()
pdf(paste(SI_fig_dir, benchmark_type, 'cell_type_ari_marker_reference_comparison_cell_downsampling2_cole.pdf', sep = ''), width = 5, height = 1) 
ggplot(aes(proportion, mean_cell_type_adj_rand), data = progressive_process_cell_sampling_res_df) + 
  geom_point(aes(color = Type), position=position_dodge(width=0.9), size = 0.5) + facet_wrap(~Type, nrow = 1) + 
  geom_errorbar(aes(color = Type, ymax = mean_cell_type_adj_rand + sd_cell_type_adj_rand, ymin=mean_cell_type_adj_rand - sd_cell_type_adj_rand), position=position_dodge(width=0.1), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("Accuracy (Branch)") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale) + 
  scale_x_continuous(limits=c(0.1, 1), breaks = seq(0.1, 1, by = 0.2)) 
dev.off()
pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_marker_reference_comparison_cell_downsampling2_helper.pdf', sep = '')) 
ggplot(aes(proportion, mean_pearson_rho), data = progressive_process_cell_sampling_res_df) + 
  geom_point(aes(color = Type), position=position_dodge(width=0.1), size = 0.5) + 
  geom_errorbar(aes(color = Type, ymax = mean_pearson_rho + sd_pearson_rho, ymin=mean_pearson_rho - sd_pearson_rho), position=position_dodge(width=0.1), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + scale_size(range = c(0.1, 1)) + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

save.image('./RData/revision_1_empirical_downsample_downsampling.RData')
