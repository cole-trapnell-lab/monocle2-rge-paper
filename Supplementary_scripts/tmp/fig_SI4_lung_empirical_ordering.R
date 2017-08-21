#empirical ordering for the HSMM data:
rm(list = ls())
####################################################################################################################################################################################
#load all package
software_custom_color_scale <- c("reference" = "#F3756C",
                                 "wishbone"="#26B24B",
                                 "monocle1" = "#00BCC3",
                                 "slicer" = "#6E95CD",
                                 "dpt" = "#CC71AD",
                                 "monocle2" = "#B8A131")
# + scale_color_manual(values=software_custom_color_scale)
software_levels <- c('monocle2', 'monocle1', 'dpt', 'wishbone', 'slicer')

library(destiny)
library(dpt)
library(monocle)
library(SLICER)
library(xacHelper)
library(diffusionMap)
library(simplePPT)

load('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/RData/fig1.RData')
load('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/RData/fig2.RData')

source('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/scripts/function.R')
source('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/scripts/plotting.R')
main_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/Figures/main_figures/"
SI_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/Figures/supplementary_figures/"

################################################################################################################################################################################################################################################
# subset cds to ensure there is only a linear trajectory
################################################################################################################################################################################################################################################

# CDK1, ID1, MYOG, MEF2C, MYH3 from Ouija paper
gene_ids <- row.names(subset(fData(HSMM_myo), gene_short_name %in% c('MYH3', 'MYOG', 'ID1', 'ENO3', 'CCNB2', 'CDK1', 'CCNB1'))) # 'MEF2C', 'TNNT1',

#confirm that state 1/2 is the linear trajectory:
plot_cell_trajectory(HSMM_myo)
plot_cell_trajectory(HSMM_myo, color_by = 'Time', markers = 'MYOG')
plot_genes_branched_pseudotime(HSMM_myo[gene_ids, ], color_by = 'Time')

#order tree based on just the marker genes
HSMM_myo_subset <- HSMM_myo[gene_ids, ]
HSMM_myo_subset <- reduceDimension(HSMM_myo_subset, auto_param_selection = F, scaling = F) #scaling = F creates a linear path
HSMM_myo_subset <- orderCells(HSMM_myo_subset, reverse = F, num_paths = 1)
plot_cell_trajectory(HSMM_myo_subset)
plot_genes_in_pseudotime(HSMM_myo_subset)
plot_genes_branched_pseudotime(HSMM_myo_subset)

#subset the state 1, 2 to get a linear trajectory: (there are still two branches)
HSMM_myo_subset2 <- HSMM_myo_subset[, pData(HSMM_myo)$State %in% c(1, 2)]
HSMM_myo_subset2 <- reduceDimension(HSMM_myo_subset2, auto_param_selection = F)
HSMM_myo_subset2 <- orderCells(HSMM_myo_subset2, reverse = F, num_paths = 1)
plot_cell_trajectory(HSMM_myo_subset2)
plot_genes_in_pseudotime(HSMM_myo_subset2)
plot_genes_branched_pseudotime(HSMM_myo_subset2)

HSMM_myo_subset <- HSMM_myo_subset2

HSMM_myo_subset2 <- HSMM_myo[, pData(HSMM_myo_subset)$State %in% c(1, 2)]
HSMM_myo_subset2 <- reduceDimension(HSMM_myo_subset2, auto_param_selection = F, norm_method = 'log', verbose = T)
HSMM_myo_subset2 <- orderCells(HSMM_myo_subset2, reverse = F, num_paths = 1)
plot_cell_trajectory(HSMM_myo_subset2[, ])
# plot_genes_in_pseudotime(HSMM_myo_subset2[, ])
# plot_genes_branched_pseudotime(HSMM_myo_subset2[, ])

################################################################################################################################################################################################################################################
#prinipal curve orderng: (this creates only a linear trajectory)
################################################################################################################################################################################################################################################

FM_subset <- as.matrix(log2(exprs(HSMM_myo_subset) + 1))
FM_subset_ordering <- DDRTree(FM_subset, dimensions = 2)
FM_subset_ordering_res <- custom_ordering(FM_subset_ordering, branch_num = 1)
qplot(FM_subset_ordering$Y[1, ], FM_subset_ordering$Y[2, ], color = pData(HSMM_myo_subset)$Time, size = as.numeric(FM_subset_ordering_res$Pseudotime))
qplot(FM_subset_ordering$Y[1, ], FM_subset_ordering$Y[2, ], color = FM_subset_ordering_res$cc_ordering[as.character(1:ncol(HSMM_myo_subset)), 'cell_state'], size = as.numeric(FM_subset_ordering_res$Pseudotime))

subset_FM_subset_ordering_res <- FM_subset_ordering_res$cc_ordering# subset(FM_subset_ordering_res$cc_ordering, as.numeric(cell_state) == 3)
root_cell <- subset_FM_subset_ordering_res$sample_name[which(subset_FM_subset_ordering_res$pseudo_time == max(subset_FM_subset_ordering_res$pseudo_time))]
FM_subset_ordering_res <- custom_ordering(FM_subset_ordering, root_cell = root_cell) #, root_cell =  add the argument for root cells

qplot(FM_subset_ordering$Y[1, ], FM_subset_ordering$Y[2, ], color = pData(HSMM_myo_subset)$Time, size = as.numeric(FM_subset_ordering_res$Pseudotime))

################################################################################################################################################################################################################################################
#empirical ordering based on marker genes
################################################################################################################################################################################################################################################

down_list <- c('ENSG00000157456.3', 'ENSG00000134057.10', 'ENSG00000170312.11', 'ENSG00000125968.7')
up_list <- c('ENSG00000108515.13', 'ENSG00000109063.9', 'ENSG00000122180.4') #'ENSG00000081189.9', 'ENSG00000105048.12',
down_df <- colMeans(exprs(HSMM_myo_subset[down_list, ]))
up_df <- colMeans(exprs(HSMM_myo_subset[up_list, ]))

down_up_df <- data.frame(up = - up_df, down = down_df)
ordered_down_up_df <- down_up_df[order(down_up_df$down, down_up_df$up, decreasing = T), ]
ordered_down_up_df$Time <- pData(HSMM_myo[, row.names(ordered_down_up_df)])$Time
qplot(-up, down, data = ordered_down_up_df, log = 'xy',  color = Time) + ggtitle('Average ordering')

norm_down_up_df <- scale(down_up_df)
norm_down_up_df <- as.data.frame(norm_down_up_df)

pdf(paste(SI_fig_dir, 'avg_ordering_mutual_exclusion.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(-up, down, data = norm_down_up_df, color = pData(HSMM_myo_subset)$Time, size = 0.5) + scale_size(range = c(0.25, 0.5)) + nm_theme() +
  xlab('Muscle avg. exp.') + ylab('Cell cycle avg. exp.')
dev.off()

#define cell clusters
cell_cycle <- subset(norm_down_up_df, down > 0 & -up < 0)
muscle <- subset(norm_down_up_df, - up > 0 & down < 0)
both <- norm_down_up_df[setdiff(row.names(norm_down_up_df), c(row.names(cell_cycle), row.names(muscle))), ]

cell_cycle_order <-row.names(cell_cycle)[order(cell_cycle$down, decreasing = T)]
muscle_order <- row.names(muscle)[order(muscle$up, decreasing = F)]
both_order <- row.names(both)[order(both$down + both$up, decreasing = T)]
cell_order_by_category <- c(cell_cycle_order, both_order, muscle_order)

#order by down and up
# norm_ordered_down_up_df <- norm_down_up_df[order(norm_down_up_df$down, norm_down_up_df$up, decreasing = T), ]
norm_ordered_down_up_df <- norm_down_up_df[cell_order_by_category, ] #ordering by category approch

#use original data:
FM_subset <- as.matrix(exprs(HSMM_myo_subset)) #

FM_subset <- FM_subset[, row.names(norm_ordered_down_up_df)]
mlt_FM_subset <- melt(FM_subset)
ncells <- ncol(FM_subset)
mlt_FM_subset$Pseudotime <- rep(1:ncells, each = length(c(up_list, down_list)))
mlt_FM_subset$gene_short_name <- fData(HSMM_myo[as.character(mlt_FM_subset$Var1), ])$gene_short_name
mlt_FM_subset$Time <- pData(HSMM_myo_subset[, mlt_FM_subset$Var2])$Time

pdf(paste(SI_fig_dir, './fig3_HSMM_ordering_gene.pdf', sep = ''), height = 1.5, width = 6)
qplot(Pseudotime, value, color = Time, data = mlt_FM_subset, size = 0.5) + facet_wrap("gene_short_name", scales = 'free', nrow = 1) + scale_size(range = c(0.25, 0.5)) + ylab('Expression') + nm_theme()
dev.off()

#add the Time:
mlt_FM_subset$Time <- pData(HSMM_myo)[as.character(mlt_FM_subset$Var2), 'Time']
qplot(Pseudotime, value, color = Time, data = mlt_FM_subset, facets = ~ gene_short_name) + ggtitle('Ordering on (scale)')

#calculate the average order for each time point:
ddply(mlt_FM_subset, .(Time), function(x) mean(x$Pseudotime))

#plot on original scale:
# ori_mlt_FM_subset <- mlt_FM_subset
# ori_mlt_FM_subset$Var2 <-
#add principal curve pseudotime in:

################################################################################################################################################################################################################################################
# compare between empirical ordering and principal curve ordering:
################################################################################################################################################################################################################################################
qplot(Pseudotime, value, color = Time, data = mlt_FM_subset, facets = ~ gene_short_name) + ggtitle('Ordering on (scale)')

DDRTree_Pseudotime_pc <- as.numeric(FM_subset_ordering_res$Pseudotime)
names(DDRTree_Pseudotime_pc) <- colnames(FM_subset)
mlt_FM_subset$Pseudotime_pc <- DDRTree_Pseudotime_pc[as.character(mlt_FM_subset$Var2)]
qplot(Pseudotime_pc, value, color = Time, data = mlt_FM_subset, facets = ~ gene_short_name)

#PCA
pca_res <- PCA(FM_subset)
qplot(pca_res[, 1], pca_res[, 2], color = pData(HSMM_myo_subset)$Time)

#simple ordering:
Pseudotime_pc <- 1:ncells
names(Pseudotime_pc) <- cell_order_by_category

################################################################################################################################################################################################################################################
#panel e, f: smoothness and lag-1 autocorrelation of the pseudotime ordering
################################################################################################################################################################################################################################################
HSMM_myo_subset <- HSMM_myo[, colnames(HSMM_myo_subset)]
HSMM_myo_subset <- setOrderingFilter(HSMM_myo_subset, gene_ids)
HSMM_myo_subset <- reduceDimension(HSMM_myo_subset, auto_param_selection = F, scaling = F) #scaling = F creates a linear path
HSMM_myo_subset <- orderCells(HSMM_myo_subset, reverse = F, num_paths = 1)

plot_cell_trajectory(HSMM_myo_subset)
plot_cell_trajectory(HSMM_myo_subset, color_by = 'Pseudotime')
plot_cell_trajectory(HSMM_myo_subset, color_by = 'Hours')

HSMM_myo_subset <- orderCells(HSMM_myo_subset, reverse = T)

ncells <- ncol(HSMM_myo_subset)

#order cell by different algorithm
ncells <- ncol(HSMM_myo_subset)

#reference pseudotime:
pData(HSMM_myo_subset)[names(table(mlt_FM_subset$Var2)), "reference_pseudotime"] <- 1:ncol(HSMM_myo_subset)

#monocle 1:
HSMM_myo_subset_monolce1 <- reduceDimension(HSMM_myo_subset, norm_method = 'log', reduction_method = 'ICA', auto_param_selection = F)
HSMM_myo_subset_monolce1 <- orderCells(HSMM_myo_subset_monolce1)
plot_cell_trajectory(HSMM_myo_subset_monolce1, color_by = 'Pseudotime')
plot_cell_trajectory(HSMM_myo_subset_monolce1, color_by = 'Hours')
HSMM_myo_subset_monolce1 <- orderCells(HSMM_myo_subset_monolce1)
monocle1_cds_pseudotime <- rank(pData(HSMM_myo_subset_monolce1)$Pseudotime)
pData(HSMM_myo_subset)$monocle1_pseudotime <- monocle1_cds_pseudotime

#slicer:
slicer_res <- run_slicer(HSMM_myo_subset, min_branch_len = 100)
pData(HSMM_myo_subset)$slicer_pseudotime <- order(slicer_res$order_df$cells_ordered)

qplot(slicer_res$traj_lle[slicer_res$order_df$cells_ordered, 1], slicer_res$traj_lle[slicer_res$order_df$cells_ordered, 2], size = 1:ncells, color = pData(HSMM_myo_subset[, ])$Hours)
qplot(slicer_res$traj_lle[slicer_res$order_df$cells_ordered, 1], slicer_res$traj_lle[slicer_res$order_df$cells_ordered, 2], size = 1:ncells, color = pData(HSMM_myo_subset[, ])$Hours)

#dpt:
HSMM_myo_subset <- orderCells(HSMM_myo_subset, reverse = T)
plot_cell_trajectory(HSMM_myo_subset, color_by = 'Time')

dpt_res <- run_new_dpt(HSMM_myo_subset, root = which.min(pData(HSMM_myo_subset)$Pseudotime))
qplot(dpt_res$dm$DC1, dpt_res$dm$DC2, size = dpt_res$pt, color = pData(HSMM_myo_subset)$Time) #

pData(HSMM_myo_subset)$dpt_seudotime <- dpt_res$pt

#wishbone pseudotime:
data <- t(log2(exprs(HSMM_myo_subset)[fData(HSMM_myo_subset)$use_for_ordering, ] + 1))
write.csv(file = paste('./csv_data/Wishbone_test_data/', 'HSMM_data_subset', ".txt", sep = ''), as.matrix(data), quote = F, row.names = T)
row.names(subset(pData(HSMM_myo_subset), Pseudotime == 0))
colnames(HSMM_myo_subset)[which.max(pData(HSMM_myo_subset)$Pseudotime)]
# [1] "T48_CT_A03"

# wishbone_res_HSMM <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/HSMM_data_subset_wishbone_res.txt', header = T, sep = '\t')
# qplot(dm1, dm2, data = wishbone_res_HSMM)
# qplot(tSNE1, tSNE2, data = wishbone_res_HSMM, color = as.character(pData(HSMM_myo_subset)$Time), size = 0.5)
#
# pdf(paste(main_fig_dir, 'wishbone_res_HSMM_wishbone_res.pdf', sep = ''), height = 1.5, width = 1.5)
# qplot(tSNE1, tSNE2, data = wishbone_res_HSMM, color = pData(HSMM_myo_subset)$Time, size = 0.5) + nm_theme() + scale_size(range = c(0.5, 0.5))
# # plot_cell_trajectory(HSMM_myo_subset, color_by = 'Time', cell_size = 0.5, show_branch_points = F) + nm_theme()#, cell_size = pData(HSMM_myo_subset)$Pseudotime
# dev.off()

# pData(HSMM_myo_subset)$wishbone_pseudotime <- wishbone_res_HSMM$trajectory

#################################################################################################################################################################################################################################################
# empirical ordering for all software (dpt, wishbone, monocle1, monocle2, slicer)
#################################################################################################################################################################################################################################################run slicer and dpt on the data:
empirical_ordering_genes_ids <- c(up_list, down_list)

slicer_res_empirical_ordering_genes <- run_slicer(HSMM_myo_subset[empirical_ordering_genes_ids, ],
                                                  start = which.min(pData(HSMM_myo_subset[empirical_ordering_genes_ids, ])$Pseudotime))
slicer_df_empirical_ordering <- data.frame(LLE1 = slicer_res_empirical_ordering_genes$traj_lle[, 1], LLE2 = slicer_res_empirical_ordering_genes$traj_lle[, 2],
                                           pseudotime = slicer_res_empirical_ordering_genes$order_df$cells_ordered, branch = slicer_res_empirical_ordering_genes$order_df$branches,
                                           Time = pData(HSMM_myo_subset[, ])$Time)
# pdf(paste(main_fig_dir, 'marseq_slicer_tree.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(LLE1, LLE2, color = Time, data = slicer_df_empirical_ordering) + nm_theme() + scale_size(range = c(0.25, 0.25))
# dev.off()

monocle1_cds_empirical_ordering_genes <- reduceDimension(HSMM_myo_subset[empirical_ordering_genes_ids, ], norm_method = 'log', reduction_method = 'ICA')
monocle1_cds_empirical_ordering_genes <- orderCells(monocle1_cds_empirical_ordering_genes, num_paths = 2)
plot_cell_trajectory(monocle1_cds_empirical_ordering_genes, color_by = 'Pseudotime')
plot_cell_trajectory(monocle1_cds_empirical_ordering_genes, color_by = 'Time')
plot_cell_trajectory(monocle1_cds_empirical_ordering_genes)

monocle2_cds_empirical_ordering_genes <- reduceDimension(HSMM_myo_subset[empirical_ordering_genes_ids, ], norm_method = 'log')
monocle2_cds_empirical_ordering_genes <- orderCells(monocle2_cds_empirical_ordering_genes)
plot_cell_trajectory(monocle2_cds_empirical_ordering_genes, color_by = 'Pseudotime')
plot_cell_trajectory(monocle2_cds_empirical_ordering_genes, color_by = 'Time')
plot_cell_trajectory(monocle2_cds_empirical_ordering_genes)
monocle2_cds_empirical_ordering_genes <- orderCells(monocle2_cds_empirical_ordering_genes, root_state = 3)

dpt_res_empirical_ordering_genes <- run_new_dpt(monocle2_cds_empirical_ordering_genes)

#write data for running wishbone:
#empirical for wishbone:
data <- t(log2(exprs(HSMM_myo_subset[empirical_ordering_genes_ids, ]) + 1))
write.csv(file = paste('./csv_data/Wishbone_test_data/', 'HSMM_linear_tree_empirical', ".txt", sep = ''), as.matrix(data), quote = F, row.names = T)
row.names(subset(pData(monocle2_cds_empirical_ordering_genes), Pseudotime == 0))

# wishbone_res_downsampling_empirical <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/HSMM_linear_tree_empirical_wishbone_res.txt', header = T, sep = '\t')
# qplot(dm1, dm2, data = wishbone_res_downsampling_empirical)
# qplot(tSNE1, tSNE2, data = wishbone_res_downsampling_empirical, color = as.character(branch), size = 0.5)

# pdf(paste(main_fig_dir, 'marseq_wishbone_tree_downsampling_empirical.pdf', sep = ''), height = 1.5, width = 1.5)
# qplot(tSNE1, tSNE2, data = wishbone_res_downsampling_empirical, color = pData(HSMM_myo_subset)$Time, size = 0.5) + nm_theme() + scale_size(range = c(0.5, 0.5))
# # plot_cell_trajectory(HSMM_myo_subset, color_by = 'Time', cell_size = 0.5, show_branch_points = F) + nm_theme()#, cell_size = pData(HSMM_myo_subset)$Pseudotime
# dev.off()

#################################################################################################################################################################################################################################################
# set pseudotime for the HSMM cells with 1k genes as ordering genes
#################################################################################################################################################################################################################################################run slicer and dpt on the data:
score_df <- unique(mlt_FM_subset[, c('Var2', 'Pseudotime')])
colnames(score_df) <- c('cell_name', 'naive_pseudotime')
row.names(score_df) <- score_df$cell_name
score_df <- score_df[colnames(HSMM_myo_subset), ]

pData(HSMM_myo_subset)$reference_pseudotime <- score_df[colnames(HSMM_myo_subset), 'naive_pseudotime']
pData(HSMM_myo_subset)$dpt_pseudotime <- rank(dpt_res$pt)
pData(HSMM_myo_subset)$wishbone_pseudotime <- rank(dpt_res$pt)# rank(hsmm_wishbone_res$trajectory)
pData(HSMM_myo_subset)$monocle1_pseudotime <- NA
pData(HSMM_myo_subset)[colnames(HSMM_myo_subset_monolce1), 'monocle1_pseudotime'] <- rank(monocle1_cds_pseudotime)
pData(HSMM_myo_subset)$slicer_pseudotime <- NA
pData(HSMM_myo_subset)[colnames(HSMM_myo_subset_monolce1), 'slicer_pseudotime'] <- slicer_res$order_df$cells_ordered

#################################################################################################################################################################################################################################################
# correlation results for result based on 1k genes ordering
#################################################################################################################################################################################################################################################
dpt_kendall_tau <- cor(score_df$naive_pseudotime, dpt_res$pt, method = 'kendall', use = 'pairwise.complete.obs')
dpt_pearson_rho <- cor(score_df$naive_pseudotime, dpt_res$pt)
wishbone_kendall_tau <- cor(score_df[, 'naive_pseudotime'], dpt_res$pt, method = 'kendall', use = 'pairwise.complete.obs') #wishbone_res_HSMM$trajectory
wishbone_pearson_rho <- cor(score_df[, 'naive_pseudotime'], dpt_res$pt) #wishbone_res_HSMM$trajectory
monocle1_kendall_tau <- cor(rank(pData(HSMM_myo_subset)$reference_pseudotime), rank(pData(HSMM_myo_subset_monolce1)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
monocle1_pearson_rho <- cor(score_df$naive_pseudotime, pData(HSMM_myo_subset_monolce1)$Pseudotime)
monocle2_kendall_tau <- cor(rank(score_df$naive_pseudotime), rank(pData(HSMM_myo_subset)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
monocle2_pearson_rho <- cor(score_df$naive_pseudotime, pData(HSMM_myo_subset)$Pseudotime)
slicer_kendall_tau <- cor(score_df[, 'naive_pseudotime'], slicer_res$order_df$cells_ordered, method = 'kendall', use = 'pairwise.complete.obs')
slicer_pearson_rho <- cor(score_df[, 'naive_pseudotime'], slicer_res$order_df$cells_ordered)

#################################################################################################################################################################################################################################################
# marker genes ordering
#################################################################################################################################################################################################################################################
# monocle1_cds_empirical_ordering_genes monocle2_cds_empirical_ordering_genes dpt_res_empirical_ordering_genes slicer_res_empirical_ordering_genes

emp_dpt_kendall_tau <- cor(rank(pData(HSMM_myo_subset)$reference_pseudotime), rank(dpt_res_empirical_ordering_genes$pt), method = 'kendall', use = 'pairwise.complete.obs')
emp_dpt_pearson_rho <- cor(pData(HSMM_myo_subset)$reference_pseudotime, dpt_res_empirical_ordering_genes$pt)
# emp_wishbone_kendall_tau <- cor(rank(pData(HSMM_myo_subset)$reference_pseudotime), rank(pData(HSMM_myo_subset)$wishbone_pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
# emp_wishbone_pearson_rho <- cor(pData(HSMM_myo_subset)$reference_pseudotime, pData(HSMM_myo_subset)$wishbone_pseudotime)
emp_monocle1_kendall_tau <- cor(rank(pData(HSMM_myo_subset)$reference_pseudotime), rank(pData(monocle1_cds_empirical_ordering_genes)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
emp_monocle1_pearson_rho <- cor(pData(HSMM_myo_subset)$reference_pseudotime, pData(monocle1_cds_empirical_ordering_genes)$Pseudotime)
emp_monocle2_kendall_tau <- cor(rank(pData(HSMM_myo_subset)$reference_pseudotime), rank(pData(monocle2_cds_empirical_ordering_genes)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
emp_monocle2_pearson_rho <- cor(pData(HSMM_myo_subset)$reference_pseudotime, pData(monocle2_cds_empirical_ordering_genes)$Pseudotime)
emp_slicer_kendall_tau <- cor(rank(pData(HSMM_myo_subset)$reference_pseudotime), rank(slicer_res_empirical_ordering_genes$order_df$cells_ordered), method = 'kendall', use = 'pairwise.complete.obs')
emp_slicer_pearson_rho <- cor(pData(HSMM_myo_subset)$reference_pseudotime, slicer_res_empirical_ordering_genes$order_df$cells_ordered)

all_cor_df <- data.frame("Pearson_cor" = c(monocle2_pearson_rho, monocle1_pearson_rho, dpt_pearson_rho, slicer_pearson_rho, slicer_pearson_rho,
                                               emp_monocle2_pearson_rho, emp_monocle1_pearson_rho, emp_dpt_pearson_rho, emp_slicer_pearson_rho, emp_slicer_pearson_rho),
                             "Kendall_cor" = c(monocle2_kendall_tau, monocle1_kendall_tau, dpt_kendall_tau, slicer_kendall_tau, slicer_kendall_tau,
                                               emp_monocle2_kendall_tau, emp_monocle1_kendall_tau, emp_dpt_kendall_tau, emp_slicer_kendall_tau, emp_slicer_kendall_tau),
                             names = c("monocle2 (all genes)", "monocle1 (all genes)", "dpt (all genes)", "slicer (all genes)", "wishbone (all genes)",
                                       "monocle2 (markers)", "monocle1 (markers)", "dpt (markers)", "slicer (markers)", "wishbone (markers)"))
mlt_all_cor_df <- melt(all_cor_df, id.vars = c('names'))
levels(mlt_all_cor_df$names) <- c("monocle2 (all genes)", "monocle2 (markers)", "monocle1 (all genes)", "monocle1 (markers)", "dpt (all genes)", "dpt (markers)",
                                  "slicer (all genes)", "slicer (markers)", "wishbone (all genes)","wishbone (markers)")

mlt_all_cor_df$software <- c('dpt', 'monocle1', 'monocle2', 'slicer', 'wishbone')
pdf(paste(main_fig_dir, 'HSMM_pseudotime_correlation.pdf', sep = ''), height = 1.5, width = 3.2)
ggplot(aes(names, abs(value)), data = mlt_all_cor_df) + geom_bar(stat="identity", aes(fill = software)) + facet_wrap(~variable) + nm_theme() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_size(range = c(0.2, 0.2)) + xlab('') + ylab('') + monocle_theme_opts() + scale_fill_manual(values = software_custom_color_scale)
dev.off()

#################################################################################################################################################################################################################################################run slicer and dpt on the data:
#perform DEG test benchmark: (test only on the ordering genes which make more sense)
#################################################################################################################################################################################################################################################run slicer and dpt on the data:
cal_autocorrelation_gene <- function(cds, order_by = 'Pseudotime'){
  # cds <- cds[, pData(cds)$State %in% States]
  ordered_cds <- cds[, order(pData(cds)[, order_by])]
  # exprs(ordered_cds) <- as.matrix(ordered_cds)
  apply(exprs(ordered_cds), 1, function(x) abs(cor(x[-length(x)], x[-1])))
}

cal_smoothness_gene <- function(cds, order_by = 'Pseudotime'){
  # cds <- cds[, pData(cds)$State %in% States]
  ordered_cds <- cds[, order(pData(cds)[, order_by])]
  # exprs(ordered_cds) <- as.matrix(ordered_cds)
  apply(exprs(ordered_cds), 1, function(x) sqrt(mean((diff(x))^2)) ) #sd(diff(x))/abs(mean(diff(x)))
}

cal_smoothness_autocorrelation <- function(cds, modelFormulaStr="~sm.ns(Pseudotime)", order_by = 'Pseudotime') {
  # pData(cds)[, order_by] <- rank(pData(cds)[, order_by])
  full_model_fits <- fitModel(cds, modelFormulaStr=modelFormulaStr, cores = detectCores() - 2)

  expression_curve_matrix <- responseMatrix(full_model_fits, cores = detectCores() - 2)
  # plot_clusters(cds, clusters)

  #different initial method:
  autocorrelation <- cal_autocorrelation_gene(cds, order_by = order_by)
  smoothness <- cal_smoothness_gene(cds, order_by = order_by)
  #different residuals: resid(full_model_fits[[1]], type = "pearson")
  #unlist(lapply(full_model_fits, function(x) sqrt(mean(residuals(x)^2)) )) #
  resi <- apply(exprs(cds) - expression_curve_matrix, 1, function(x) sqrt(mean((x)^2))  )

  return(data.frame(autocorrelation = autocorrelation, smoothness = smoothness, residuals = resi))
}

#run two analysis at once: (run only on the marker genes)
pData(HSMM_myo_subset)$monocle1_pseudotime <- pData(monocle1_cds_empirical_ordering_genes)$Pseudotime
pData(HSMM_myo_subset)$monocle2_pseudotime <- pData(monocle2_cds_empirical_ordering_genes)$Pseudotime
pData(HSMM_myo_subset)$slicer_pseudotime <- slicer_res_empirical_ordering_genes$order_df$cells_ordered
pData(HSMM_myo_subset)$dpt_pseudotime <- dpt_res_empirical_ordering_genes$pt
pData(HSMM_myo_subset)$wishbone_pseudotime <- dpt_res_empirical_ordering_genes$pt

plot_genes_in_pseudotime(HSMM_myo_subset[c(up_list, down_list), ], trend_formula = "~ sm.ns(reference_pseudotime, df=3)")
plot_genes_in_pseudotime(HSMM_myo_subset[c(up_list, down_list), ], trend_formula = "~ sm.ns(rank(monocle1_pseudotime), df=3)")
plot_genes_in_pseudotime(HSMM_myo_subset[c(up_list, down_list), ], trend_formula = "~ sm.ns(rank(Pseudotime), df=3)")
# plot_genes_in_pseudotime(HSMM_myo_subset[c(up_list, down_list), ], trend_formula = "~ sm.ns(slicer_pseudotime, df=3)")
# plot_genes_in_pseudotime(HSMM_myo_subset[c(up_list, down_list), ], trend_formula = "~ sm.ns(dpt_pseudotime, df=3)")
# plot_genes_in_pseudotime(HSMM_myo_subset[c(up_list, down_list), ], trend_formula = "~ sm.ns(wishbone_pseudotime, df=3)")

reference_HSMM_full_model_fits <- cal_smoothness_autocorrelation(HSMM_myo_subset[c(up_list, down_list), ], #HSMM_ordering_genes
                                                                 modelFormulaStr="~sm.ns(reference_pseudotime)", order_by = 'reference_pseudotime')

reference_HSMM_full_model_fits <- cal_smoothness_autocorrelation(HSMM_myo_subset[c(up_list, down_list), ], #HSMM_ordering_genes
                                                                 modelFormulaStr="~sm.ns(reference_pseudotime)", order_by = 'reference_pseudotime')

monocle1_HSMM_full_model_fits <- cal_smoothness_autocorrelation(HSMM_myo_subset[c(up_list, down_list), ], #HSMM_ordering_genes
                                                                modelFormulaStr="~sm.ns(rank(monocle1_pseudotime))", order_by = 'monocle1_pseudotime')

monocle2_HSMM_full_model_fits <- cal_smoothness_autocorrelation(HSMM_myo_subset[c(up_list, down_list), ], #HSMM_ordering_genes
                                                                modelFormulaStr="~sm.ns(rank(monocle2_pseudotime))", order_by = 'monocle2_pseudotime')

slicer_HSMM_full_model_fits <- cal_smoothness_autocorrelation(HSMM_myo_subset[c(up_list, down_list), ], #HSMM_ordering_genes
                                                              modelFormulaStr="~sm.ns(rank(slicer_pseudotime))", order_by = 'slicer_pseudotime')

dpt_HSMM_model_fits <- cal_smoothness_autocorrelation(HSMM_myo_subset[c(up_list, down_list), ], #HSMM_ordering_genes
                                                      modelFormulaStr="~sm.ns(rank(dpt_pseudotime))", order_by = 'dpt_pseudotime')

wishbone_HSMM_full_model_fits <- cal_smoothness_autocorrelation(HSMM_myo_subset[c(up_list, down_list), ],
                                                                modelFormulaStr="~sm.ns(rank(wishbone_pseudotime))", order_by = 'wishbone_pseudotime')

#################################################################################################################################################################################################################################################run slicer and dpt on the data:
# reference_HSMM_full_model_fits
# reference_MEP_full_model_fits
# monocle2_HSMM_full_model_fits
# monocle2_MEP_full_model_fits
# dpt_HSMM_model_fits
# dpt_MEP_full_model_fits

cluster_autocorrelation_df <- data.frame(reference = c(reference_HSMM_full_model_fits$autocorrelation),
                                         monocle2 = c(monocle2_HSMM_full_model_fits$autocorrelation),
                                         monocle1 = c(monocle1_HSMM_full_model_fits$autocorrelation),
                                         #wishbone = 1/c(wishbone_HSMM_full_model_fits$smoothness),
                                         slicer = c(slicer_HSMM_full_model_fits$autocorrelation),
                                         dpt = c(dpt_HSMM_model_fits$autocorrelation),
                                         type = "cluster_autocorrelation", lineage = rep(c('HSMM'), length(reference_HSMM_full_model_fits$autocorrelation)))

cluster_smoothness_df <- data.frame(reference = c(reference_HSMM_full_model_fits$smoothness),
                                    monocle2 = c(monocle2_HSMM_full_model_fits$smoothness),
                                    monocle1 = c(monocle1_HSMM_full_model_fits$smoothness),
                                    #wishbone = 1/c(wishbone_HSMM_full_model_fits$smoothness),
                                    slicer = c(slicer_HSMM_full_model_fits$smoothness),
                                    dpt = c(dpt_HSMM_model_fits$smoothness),
                                    type = "cluster_smooothness", lineage = rep(c('HSMM'), length(reference_HSMM_full_model_fits$smoothness)))

cluster_residual_df <- data.frame(reference = c(reference_HSMM_full_model_fits$residuals),
                                    monocle2 = c(monocle2_HSMM_full_model_fits$residuals),
                                    monocle1 = c(monocle1_HSMM_full_model_fits$residuals),
                                    #wishbone = 1/c(wishbone_HSMM_full_model_fits$smoothness),
                                    slicer = c(slicer_HSMM_full_model_fits$residuals),
                                    dpt = c(dpt_HSMM_model_fits$residuals),
                                    type = "cluster_residual", lineage = rep(c('HSMM'), length(reference_HSMM_full_model_fits$residuals)))

plot_tightness_smoothness <- function(tightness_df, smoothness_df, residual_df) {
    tightness_smoothness_df <- Reduce(rbind, list(tightness_df, smoothness_df, residual_df))

    all_finite_genes <- apply(apply(tightness_smoothness_df[, 1:(ncol(tightness_smoothness_df) - 2)], 2, function(x) is.finite(x)), 1, all)
    all_finite_genes <- names(all_finite_genes)[all_finite_genes]

    tightness_smoothness_df$gene_ids <- row.names(tightness_smoothness_df)
    mlt_ts_df <- melt(tightness_smoothness_df[all_finite_genes, ], id.var = c('type', 'gene_ids'))

    ddply(mlt_ts_df, .(type, variable), function(x) range(x$value, na.rm = T))

    #density plot:
    p1 <- ggplot(aes(x = value + 1), data = mlt_ts_df) + geom_density(aes(color = variable)) + facet_wrap(~type, scale = 'free') + scale_x_log10()

    #scatterplot:
    p2 <- ggplot(aes(x = gene_ids, y = value + 1), data = mlt_ts_df) + geom_point(aes(color = variable)) + scale_y_log10() + facet_wrap(~type, scale = 'free')

    #barplot:
    p3 <- ggplot(aes(x = variable, y = value + 1), data = mlt_ts_df) + geom_boxplot(aes(color = variable)) + scale_y_log10() + facet_wrap(~type, scale = 'free')

    return(list(mlt_ts_df = mlt_ts_df, p1 = p1, p2 = p2, p3 = p3))
  }

tightness_smoothness_res <- plot_tightness_smoothness(cluster_autocorrelation_df[, -ncol(cluster_autocorrelation_df)],
                                                      cluster_smoothness_df[, -ncol(cluster_autocorrelation_df)],
                                                      cluster_residual_df[, -ncol(cluster_autocorrelation_df)])

pdf(paste(SI_fig_dir, 'smoothness_tightness.pdf', sep = ''), height = 2, width = 3)
ggplot(aes(x = variable, y = value), data = tightness_smoothness_res$mlt_ts_df) + geom_boxplot(aes(color = variable))  + monocle_theme_opts() +
  facet_wrap(~type, scale = 'free') + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab('') + ylab('')
dev.off()

#only look at the lag-1 autocorrelation:
subset_tightness_smoothness_res <- subset(tightness_smoothness_res$mlt_ts_df, type == 'cluster_autocorrelation')
subset_tightness_smoothness_res$variable <- as.character(subset_tightness_smoothness_res$variable)
subset_tightness_smoothness_res$variable <- factor(subset_tightness_smoothness_res$variable, levels = c('reference', software_levels))

pdf(paste(main_fig_dir, 'autocorrelation_HSMM_cells.pdf', sep = ''), height = 1.5, width = 1.5)
ggplot(aes(x = variable, y = abs(value)), data = subset_tightness_smoothness_res) + geom_boxplot(aes(color = variable))  + monocle_theme_opts()  +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab('') + ylab('lag1_autocorrelation') #+ scale_color_manual(values = software_custom_color_scale)
dev.off()

#only look at the residual:
subset_tightness_smoothness_res <- subset(tightness_smoothness_res$mlt_ts_df, type == 'cluster_residual')
levels(subset_tightness_smoothness_res$variable) <- c('reference', software_levels)

pdf(paste(main_fig_dir, 'residual_HSMM_cells.pdf', sep = ''), height = 1.5, width = 1.5)
ggplot(aes(x = variable, y = value + 1), data = subset_tightness_smoothness_res) + geom_boxplot(aes(color = variable))  + monocle_theme_opts() + scale_y_log10() + 
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab('') + ylab('Residual of gene expression')# + scale_color_manual(values=software_custom_color_scale)
dev.off()

####################################################################################################################################################################################
# run algorithm with different initial methods  
####################################################################################################################################################################################

HSMM_ordering_genes <- row.names(HSMM_myo[fData(HSMM_myo)$use_for_ordering, ])
ordering_1 <- order_cell_tree(HSMM_ordering_genes, HSMM_myo[, colnames(HSMM_myo_subset)], initial_method = PCA)
ordering_2 <- order_cell_tree(HSMM_ordering_genes, HSMM_myo[, colnames(HSMM_myo_subset)], initial_method = normalized_PCA)
ordering_3 <- order_cell_tree(HSMM_ordering_genes, HSMM_myo[, colnames(HSMM_myo_subset)], initial_method = ICA)
ordering_4 <- order_cell_tree(HSMM_ordering_genes, HSMM_myo[, colnames(HSMM_myo_subset)], initial_method = LLE) #error here
ordering_5 <- order_cell_tree(HSMM_ordering_genes, HSMM_myo[, colnames(HSMM_myo_subset)], initial_method = DM)
ordering_6 <- order_cell_tree(HSMM_ordering_genes, HSMM_myo[, colnames(HSMM_myo_subset)], initial_method = destiny_diffusionMaps)

pdf(paste(main_fig_dir, 'PCA_HSMM_ini.pdf', sep = ''), height = 1, width = 1)
monocle_theme_opts() + nm_theme() 
dev.off()

pdf(paste(main_fig_dir, 'ICA_HSMM_ini.pdf', sep = ''), height = 1, width = 1)
monocle_theme_opts() + nm_theme() 
dev.off()

pdf(paste(main_fig_dir, 'LLE_HSMM_ini.pdf', sep = ''), height = 1, width = 1)
monocle_theme_opts() + nm_theme() 
dev.off()

pdf(paste(main_fig_dir, 'DM_HSMM_ini.pdf', sep = ''), height = 1, width = 1)
monocle_theme_opts() + nm_theme() 
dev.off()

pdf(paste(main_fig_dir, 'destiny_diffusionMaps_ini.pdf', sep = ''), height = 1, width = 1)
monocle_theme_opts() + nm_theme() 
dev.off()

p1 <- plot_order_cell_tree(HSMM_ordering_genes, HSMM_myo, initial_method = PCA) + ggtitle('PCA (top 100)')
p2 <- plot_order_cell_tree(HSMM_ordering_genes, HSMM_myo, initial_method = normalized_PCA) + ggtitle('normalized PCA (top 100)')
p3 <- plot_order_cell_tree(HSMM_ordering_genes, HSMM_myo, initial_method = ICA) + ggtitle('ICA (top 100)')
p4 <- plot_order_cell_tree(HSMM_ordering_genes, HSMM_myo, initial_method = LLE) + ggtitle('LLE (top 100)')
p5 <- plot_order_cell_tree(HSMM_ordering_genes, HSMM_myo, initial_method = DM) + ggtitle('DM (top 100)')
p6 <- plot_order_cell_tree(HSMM_ordering_genes, HSMM_myo, initial_method = destiny_diffusionMaps) + ggtitle('German group; destiny dm (top 100)')

multiplot(plotlist = list(p1, p2, p3, p4, p5, p6), cols = 2)
plot_spanning_tree(HSMM_myo, cell_size=2)
ordering_pseudotime1 <- pData(ordering_1)
ordering_pseudotime2 <- pData(ordering_2)
ordering_pseudotime3 <- pData(ordering_3)
ordering_pseudotime4 <- pData(ordering_4)
ordering_pseudotime5 <- pData(ordering_5)
ordering_pseudotime6 <- pData(ordering_6)

#Kendall's tau analysis on the data: 
cor(rank(ordering_pseudotime2[names(Pseudotime_pc), 'Pseudotime']), rank(ordering_pseudotime1[names(Pseudotime_pc), 'Time']), method = 'kendall', use = 'pairwise.complete.obs')

ddply(ordering_pseudotime1, .(Time), function(x) mean(x$Pseudotime))
ddply(ordering_pseudotime2, .(Time), function(x) mean(x$Pseudotime))
ddply(ordering_pseudotime3, .(Time), function(x) mean(x$Pseudotime))
ddply(ordering_pseudotime4, .(Time), function(x) mean(x$Pseudotime))
ddply(ordering_pseudotime5, .(Time), function(x) mean(x$Pseudotime))
ddply(ordering_pseudotime6, .(Time), function(x) mean(x$Pseudotime))

# PCA, norm_pca, ICA, LLE, DM, destiny_dm
pcor1 <- cor(rank(Pseudotime_pc), rank(ordering_pseudotime1[names(Pseudotime_pc), 'Pseudotime']))
pcor2 <- cor(Pseudotime_pc, ordering_pseudotime2[names(Pseudotime_pc), 'Pseudotime'])
pcor3 <- cor(Pseudotime_pc, ordering_pseudotime3[names(Pseudotime_pc), 'Pseudotime'])
pcor4 <- cor(Pseudotime_pc, ordering_pseudotime4[names(Pseudotime_pc), 'Pseudotime'])
pcor5 <- cor(Pseudotime_pc, ordering_pseudotime5[names(Pseudotime_pc), 'Pseudotime'])
pcor6 <- cor(Pseudotime_pc, ordering_pseudotime6[names(Pseudotime_pc), 'Pseudotime'])

kcor1 <- cor(rank(Pseudotime_pc), rank(ordering_pseudotime1[names(Pseudotime_pc), 'Pseudotime']), method = 'kendall', use = 'pairwise.complete.obs')
kcor2 <- cor(Pseudotime_pc, ordering_pseudotime2[names(Pseudotime_pc), 'Pseudotime'], method = 'kendall', use = 'pairwise.complete.obs')
kcor3 <- cor(Pseudotime_pc, ordering_pseudotime3[names(Pseudotime_pc), 'Pseudotime'], method = 'kendall', use = 'pairwise.complete.obs')
kcor4 <- cor(Pseudotime_pc, ordering_pseudotime4[names(Pseudotime_pc), 'Pseudotime'], method = 'kendall', use = 'pairwise.complete.obs')
kcor5 <- cor(Pseudotime_pc, ordering_pseudotime5[names(Pseudotime_pc), 'Pseudotime'], method = 'kendall', use = 'pairwise.complete.obs')
kcor6 <- cor(Pseudotime_pc, ordering_pseudotime6[names(Pseudotime_pc), 'Pseudotime'], method = 'kendall', use = 'pairwise.complete.obs')

Time_kcor1 <- cor(rank(ordering_pseudotime1[names(Pseudotime_pc), 'Pseudotime']), rank(ordering_pseudotime1[names(Pseudotime_pc), 'Time']), method = 'kendall', use = 'pairwise.complete.obs')
Time_kcor2 <-cor(rank(ordering_pseudotime2[names(Pseudotime_pc), 'Pseudotime']), rank(ordering_pseudotime2[names(Pseudotime_pc), 'Time']), method = 'kendall', use = 'pairwise.complete.obs')
Time_kcor3 <-cor(rank(ordering_pseudotime3[names(Pseudotime_pc), 'Pseudotime']), rank(ordering_pseudotime3[names(Pseudotime_pc), 'Time']), method = 'kendall', use = 'pairwise.complete.obs')
Time_kcor4 <-cor(rank(ordering_pseudotime4[names(Pseudotime_pc), 'Pseudotime']), rank(ordering_pseudotime4[names(Pseudotime_pc), 'Time']), method = 'kendall', use = 'pairwise.complete.obs')
Time_kcor5 <-cor(rank(ordering_pseudotime5[names(Pseudotime_pc), 'Pseudotime']), rank(ordering_pseudotime5[names(Pseudotime_pc), 'Time']), method = 'kendall', use = 'pairwise.complete.obs')
Time_kcor6 <-cor(rank(ordering_pseudotime6[names(Pseudotime_pc), 'Pseudotime']), rank(ordering_pseudotime6[names(Pseudotime_pc), 'Time']), method = 'kendall', use = 'pairwise.complete.obs')

cor_df <- data.frame("Pearson_cor" = c(pcor1, pcor2, pcor3, pcor4, pcor5, pcor6), 
                     "Kendall_cor" = c(kcor1, kcor2, kcor3, kcor4, kcor5, kcor6), 
                     "Time_kcor" = c(Time_kcor1, Time_kcor2, Time_kcor3, Time_kcor4, Time_kcor5, Time_kcor6),
                     names = c("PCA", "norm_pca", "ICA", "LLE", "DM", "destiny_dm"))

mlt_cor_df <- melt(cor_df, id.vars = c('names'))

pdf(paste(SI_fig_dir, 'fig_si3a_initial_method.pdf', sep = ''), height = 1.5, width = 3)
ggplot(aes(names, abs(value)), data = mlt_cor_df) + geom_bar(aes(fill = names), stat = 'identity') + nm_theme() + facet_wrap(~ variable) + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_size(range = c(0.2, 0.2)) + xlab('') + ylab('') + monocle_theme_opts()
dev.off()

####################################################################################################################################################################################
# perform dimension reduction with differnt methods but just learn the structure with simplePPT
####################################################################################################################################################################################
HSMM_ordering_genes <- row.names(HSMM_myo[fData(HSMM_myo)$use_for_ordering, ])

run_dpt_simplePPT <- function(data, branching = F, norm_method = 'log', root = NULL, verbose = F){
  if(verbose)
    message('root should be the id to the cell not the cell name ....')
  
  data <- t(data[!duplicated(data), ])
  dm <- DiffusionMap(as.matrix(data))
  dpt <- DPT(dm, branching = branching)
  
  ts <- dm@transitions
  M <- destiny:::accumulated_transitions(dm)
  
  branch <- dpt@branch
  row.names(branch) <- row.names(data[!duplicated(data), ])
  
  return(dm@eigenvectors[, 1:2])
}

reduceDimension_initial_simplePPT <- function(cds, 
                            max_components=2, 
                            norm_method = c("vstExprs", "log", "none"), 
                            residualModelFormulaStr=NULL,
                            pseudo_expr=NULL, 
                            relative_expr=TRUE,
                            auto_param_selection = TRUE,
                            verbose=FALSE,
                            scaling = TRUE,
                            initial_method = run_dpt_simplePPT,
                            ...){
  extra_arguments <- list(...)
  FM <- normalize_expr_data(cds, norm_method, pseudo_expr)
  
  #FM <- FM[unlist(sparseApply(FM, 1, sd, convert_to_dense=TRUE)) > 0, ]
  xm <- Matrix::rowMeans(FM)
  xsd <- sqrt(Matrix::rowMeans((FM - xm)^2))
  FM <- FM[xsd > 0,]
  
  if (is.null(residualModelFormulaStr) == FALSE) {
    if (verbose) 
      message("Removing batch effects")
    X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr), 
                                       data = pData(cds), drop.unused.levels = TRUE)
    
    fit <- limma::lmFit(FM, X.model_mat, ...)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
  }else{
    X.model_mat <- NULL
  }
  
  if(scaling){
    FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
    FM <- FM[!is.na(row.names(FM)), ]
  } else FM <- as.matrix(FM)
  
  if (nrow(FM) == 0) {
    stop("Error: all rows have standard deviation zero")
  }
  
  FM <- FM[apply(FM, 1, function(x) all(is.finite(x))), ] #ensure all the expression values are finite values

  if(verbose)
    message('running DPT ...')
  reduced_data <- initial_method(FM)
  
  if(verbose)
    message('running simplePPT ...')
  simplePPT_res <- principal_tree(t(reduced_data[, 1:2]), MU = NULL, lambda = 1, bandwidth = 30, maxIter = 10, verbose = verbose)
  
  colnames(simplePPT_res$MU) <- colnames(FM) #paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
  DCs <- t(reduced_data[, 1:2])
  colnames(DCs) <- colnames(FM) #paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
  
  monocle:::reducedDimW(cds) <- DCs
  monocle:::reducedDimS(cds) <- simplePPT_res$MU
  monocle:::reducedDimK(cds) <- simplePPT_res$MU
  cds@auxOrderingData[["simplePPT"]]$objective_vals <- NA
  
  adjusted_K <- Matrix::t(reducedDimK(cds))
  dp <- as.matrix(dist(adjusted_K))
  cellPairwiseDistances(cds) <- dp
  gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)
  minSpanningTree(cds) <- dp_mst
  cds@dim_reduce_type <- "simplePPT"
  
  cds
}

order_cell_tree2 <- function(order_genes, cds, initial_method = PCA, norm_method = 'log',order_by="Time") {
  cds <- setOrderingFilter(cds, order_genes)
  plot_ordering_genes(cds)
  
  cds <- reduceDimension_initial_simplePPT(cds, max_components = 2, norm_method = norm_method, initial_method = initial_method, verbose = T) #, initial_method = DM , initial_method = DM
  cds <- orderCells(cds, num_paths=1, root_state = NULL)
  PD <- pData(cds)
  
  if(!is.null(order_by)){
    avg_pseudotime <- ddply(PD, .(Time), function(x) mean(x$Pseudotime))
    starting_cell_states <- apply(table(PD[, c("Time", "State")]), 1, function(x) which(x == max(x)))[1]
    cds <- orderCells(cds, num_paths=1, root_state = starting_cell_states)
  }
  # return(pData(HSMM_myo))
  return(cds)
}

ordering_1 <- order_cell_tree2(HSMM_ordering_genes, HSMM_myo[, colnames(HSMM_myo_subset)], initial_method = PCA)
ordering_3 <- order_cell_tree2(HSMM_ordering_genes, HSMM_myo[, colnames(HSMM_myo_subset)], initial_method = ICA)
ordering_4 <- order_cell_tree2(HSMM_ordering_genes, HSMM_myo[, colnames(HSMM_myo_subset)], initial_method = LLE) #error here
ordering_5 <- order_cell_tree2(HSMM_ordering_genes, HSMM_myo[, colnames(HSMM_myo_subset)], initial_method = run_dpt_simplePPT)

pdf(paste(SI_fig_dir, 'ordering_1_simplePPT.pdf', sep = ''), height = 1.27, width = 1.27)
plot_cell_trajectory(ordering_1, color_by = "Time") + nm_theme() 
dev.off()

pdf(paste(SI_fig_dir, 'ordering_3_simplePPT.pdf', sep = ''), height = 1.27, width = 1.27)
plot_cell_trajectory(ordering_3, color_by = "Time") + nm_theme() 
dev.off()

pdf(paste(SI_fig_dir, 'ordering_4_simplePPT.pdf', sep = ''), height = 1.27, width = 1.27)
plot_cell_trajectory(ordering_4, color_by = "Time") + nm_theme() 
dev.off()

pdf(paste(SI_fig_dir, 'ordering_5_simplePPT.pdf', sep = ''), height = 1.27, width = 1.27)
plot_cell_trajectory(ordering_5, color_by = "Time") + nm_theme() 
dev.off()

ordering_pseudotime1 <- pData(ordering_1)
# ordering_pseudotime2 <- pData(ordering_2)
ordering_pseudotime3 <- pData(ordering_3)
ordering_pseudotime4 <- pData(ordering_4)
ordering_pseudotime5 <- pData(ordering_5)
# ordering_pseudotime6 <- pData(ordering_6)

#Kendall's tau analysis on the data: 
cor(rank(ordering_pseudotime2[names(Pseudotime_pc), 'Pseudotime']), rank(ordering_pseudotime1[names(Pseudotime_pc), 'Time']), method = 'kendall', use = 'pairwise.complete.obs')

ddply(ordering_pseudotime1, .(Time), function(x) mean(x$Pseudotime))
# ddply(ordering_pseudotime2, .(Time), function(x) mean(x$Pseudotime))
ddply(ordering_pseudotime3, .(Time), function(x) mean(x$Pseudotime))
ddply(ordering_pseudotime4, .(Time), function(x) mean(x$Pseudotime))
ddply(ordering_pseudotime5, .(Time), function(x) mean(x$Pseudotime))
# ddply(ordering_pseudotime6, .(Time), function(x) mean(x$Pseudotime))

# PCA, norm_pca, ICA, LLE, DM, destiny_dm
pcor1 <- cor(rank(Pseudotime_pc), rank(ordering_pseudotime1[names(Pseudotime_pc), 'Pseudotime']))
pcor3 <- cor(Pseudotime_pc, ordering_pseudotime3[names(Pseudotime_pc), 'Pseudotime'])
pcor4 <- cor(Pseudotime_pc, ordering_pseudotime4[names(Pseudotime_pc), 'Pseudotime'])
pcor5 <- cor(Pseudotime_pc, ordering_pseudotime5[names(Pseudotime_pc), 'Pseudotime'])

kcor1 <- cor(rank(Pseudotime_pc), rank(ordering_pseudotime1[names(Pseudotime_pc), 'Pseudotime']), method = 'kendall', use = 'pairwise.complete.obs')
kcor3 <- cor(Pseudotime_pc, ordering_pseudotime3[names(Pseudotime_pc), 'Pseudotime'], method = 'kendall', use = 'pairwise.complete.obs')
kcor4 <- cor(Pseudotime_pc, ordering_pseudotime4[names(Pseudotime_pc), 'Pseudotime'], method = 'kendall', use = 'pairwise.complete.obs')
kcor5 <- cor(Pseudotime_pc, ordering_pseudotime5[names(Pseudotime_pc), 'Pseudotime'], method = 'kendall', use = 'pairwise.complete.obs')

Time_kcor1 <- cor(rank(ordering_pseudotime1[names(Pseudotime_pc), 'Pseudotime']), rank(ordering_pseudotime1[names(Pseudotime_pc), 'Time']), method = 'kendall', use = 'pairwise.complete.obs')
Time_kcor3 <-cor(rank(ordering_pseudotime3[names(Pseudotime_pc), 'Pseudotime']), rank(ordering_pseudotime3[names(Pseudotime_pc), 'Time']), method = 'kendall', use = 'pairwise.complete.obs')
Time_kcor4 <-cor(rank(ordering_pseudotime4[names(Pseudotime_pc), 'Pseudotime']), rank(ordering_pseudotime4[names(Pseudotime_pc), 'Time']), method = 'kendall', use = 'pairwise.complete.obs')
Time_kcor5 <-cor(rank(ordering_pseudotime5[names(Pseudotime_pc), 'Pseudotime']), rank(ordering_pseudotime5[names(Pseudotime_pc), 'Time']), method = 'kendall', use = 'pairwise.complete.obs')

cor_df <- data.frame("Pearson_cor" = c(pcor1, pcor2, pcor4, pcor5), 
                     "Kendall_cor" = c(kcor1, kcor2, kcor4, kcor5), 
                     "Time_kcor" = c(Time_kcor1, Time_kcor2, Time_kcor4, Time_kcor5),
                     names = c("PCA", "ICA", "LLE", "DPT"))

mlt_cor_df <- melt(cor_df, id.vars = c('names'))

pdf(paste(SI_fig_dir, 'fig_si3a_initial_method_simplePPT.pdf', sep = ''), height = 1.5, width = 3)
ggplot(aes(names, abs(value)), data = mlt_cor_df) + geom_bar(aes(fill = names), stat = 'identity') + nm_theme() + facet_wrap(~ variable) + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_size(range = c(0.2, 0.2)) + xlab('') + ylab('') + monocle_theme_opts()
dev.off()

####################################################################################################################################################################################
# save data 
####################################################################################################################################################################################
save.image('./RData/fig_SI4_HSMM_empirical_ordering.RData')
