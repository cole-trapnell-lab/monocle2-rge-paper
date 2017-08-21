rm(list = ls())

software_custom_color_scale <- c("reference" = "#F3756C",
                             "wishbone"="#26B24B",
                             "monocle1" = "#00BCC3",
                             "slicer" = "#6E95CD",
                             "dpt" = "#CC71AD",
                             "monocle2" = "#B8A131")
# + scale_color_manual(values=software_custom_color_scale)
software_levels <- c('monocle2', 'monocle1', 'dpt', 'wishbone', 'slicer')
Mar_seq_cols <- c("CMP" = "#44AF69", "DC" = "#F8333C", "erythroid" = "#FCAB10", "GMP" = "#2B9EB3")

#perform empirical ordering and corresponding benchmarking on the downsampled cells
library(dpt)
library(SLICER)
library(monocle)
library(mcclust)
library(destiny)
library(xacHelper)
library(reshape2)
library(plyr)
library(stringr)
library(GEOquery)
library(R.matlab)
#################################################################################################################################################################################################################################################
#save the data:
load('./RData/analysis_score_ordering_MAR_seq.RData')
main_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/Figures/main_figures/"
SI_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/Figures/supplementary_figures/"

#################################################################################################################################################################################################################################################
#empirical ordering on the downsampled cells:
#################################################################################################################################################################################################################################################
if(!exists("use_downsampling_cell"))
  use_downsampling_cell <- T # only run on the downsampling cell

#sample cells for monocle1 and slicer
set.seed(20170111)
sample_cells <- colnames(valid_subset_GSE72857_cds)[sample(1:ncol(valid_subset_GSE72857_cds), 300)]

valid_subset_GSE72857_cds <- valid_subset_GSE72857_cds[, sample_cells]
#################################################################################################################################################################################################################################################
# run the results
################################################################################################################################################
#mep score VS gmp lineage score:
valid_genes = row.names(valid_subset_GSE72857_cds[rowSums(exprs(valid_subset_GSE72857_cds)) > 0, ])
total_expression = apply(exprs(valid_subset_GSE72857_cds)[valid_genes, ], 1, mean)
quartiles = cut(total_expression, breaks=quantile(total_expression, probs=seq(0,1, by=0.04))[-(1:2)], include.lowest=TRUE)
names(quartiles) = names(total_expression)

retrieve_control_genes <- function(marker_genes, quartiles){
  quartiles_tab <- table(as.character(quartiles[marker_genes]))
  control_ids <- c()
  for(quartile in names(quartiles_tab)) {
    control_ids <- c(control_ids, sample(which(quartiles == quartile), quartiles_tab[quartile] * 25))
  }

  control_genes <- names(quartiles)[control_ids]
}

set.seed(2016)
mep_control_genes <- retrieve_control_genes(mep_genes, quartiles)
set.seed(2016)
cmp_control_genes <- retrieve_control_genes(cmp_genes, quartiles)
set.seed(2016)
gmp_control_genes <- retrieve_control_genes(gmp_genes, quartiles)

mep_lineage_score <- lineage_score(valid_subset_GSE72857_cds, mep_genes, mep_control_genes)
gmp_lineage_score <- lineage_score(valid_subset_GSE72857_cds, gmp_genes, gmp_control_genes)

lineage_score_df <- data.frame(mep_lineage_score = mep_lineage_score, gmp_lineage_score = gmp_lineage_score)

#stemness score:
# Stem(i) = average[Er(G_{stem})] - average[Er(G_{stem}^{cont})] - LIN(i)

cmp_lineage_score_ori <- lineage_score(valid_subset_GSE72857_cds, cmp_genes, cmp_control_genes)
cmp_lineage_score <- cmp_lineage_score_ori - apply(lineage_score_df, 1, max)

#plot the result:
score_df <- data.frame(stemness = cmp_lineage_score, lineage_score = apply(lineage_score_df, 1, function(x) {
  x <- c(-x[1], x[2])
  x[which.max(abs(x))]
}))

score_df$lineage <- "GMP"
score_df$lineage[score_df$lineage_score < 0] <- 'MEP'
score_df$lineage[score_df$stemness > 0] <- 'CMP'
row.names(score_df) <- names(mep_lineage_score)

#use just the distances to the point (0, 0)
score_df$naive_pseudotime <- apply(score_df, 1, function(x) sqrt(sum(as.numeric(x[c(1, 2)])^2)))
score_df$naive_pseudotime[score_df$lineage == 'CMP'] <- score_df$naive_pseudotime[score_df$stemness == max(score_df$stemness)] -
  score_df$naive_pseudotime[score_df$lineage == 'CMP'] #pseduotime starts from the cell with highest stemness score 
score_df$naive_pseudotime[score_df$lineage != 'CMP'] <- score_df$naive_pseudotime[score_df$stemness == max(score_df$stemness)] +
  score_df$naive_pseudotime[score_df$lineage != 'CMP']

score_df$cell_type <- pData(valid_subset_GSE72857_cds)$cell_type

pdf(paste(main_fig_dir, 'downsampling_marseq_score.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(lineage_score, stemness, color = cell_type, size = naive_pseudotime, data = score_df) + nm_theme() + scale_size(range = c(0.2, 1.2))
dev.off()

pdf(paste(main_fig_dir, 'downsampling_marseq_score_helper.pdf', sep = ''))
qplot(lineage_score, stemness, color = lineage, size = naive_pseudotime, data = score_df) + scale_size(range = c(0.2, 1.2))
dev.off()

table(score_df$lineage)

calClusteringMetrics(score_df$lineage, score_df$cell_type)
# randIndex                     Type
# 1 0.2920114               rand index
# 2 1.6663127 variation of information
# 3 0.2920114      adjusted rand index

################################################################################################################################################
#run the dpt, wishbone, slicer:
################################################################################################################################################
num_cells_expressed_percent <- 0.1
qval_thrsld <- 0.01

if(use_downsampling_cell) {
  valid_subset_GSE72857_cds <- valid_subset_GSE72857_cds[, sample_cells]
  #1. determine how many pca dimension you want:
  valid_subset_GSE72857_cds <- detectGenes(valid_subset_GSE72857_cds)
  fData(valid_subset_GSE72857_cds)$use_for_ordering <- F

  num_cells_expressed <- round(num_cells_expressed_percent * ncol(valid_subset_GSE72857_cds))
  fData(valid_subset_GSE72857_cds)$use_for_ordering[fData(valid_subset_GSE72857_cds)$num_cells_expressed > num_cells_expressed] <- T

  valid_subset_GSE72857_cds@auxClusteringData[["tSNE"]]$variance_explained <- NULL
  MAP_pc_variance <- plot_pc_variance_explained(valid_subset_GSE72857_cds, return_all = T)

  #2. run reduceDimension with tSNE as the reduction_method
  # valid_subset_GSE72857_cds <- setOrderingFilter(valid_subset_GSE72857_cds, quake_id)
  valid_subset_GSE72857_cds <- reduceDimension(valid_subset_GSE72857_cds, max_components=2, norm_method = 'log', reduction_method = 'tSNE', num_dim = 6,  verbose = T) #, residualModelFormulaStr = '~groups.1'

  #check the embedding on each PCA components:

  #3. initial run of clusterCells_Density_Peak
  valid_subset_GSE72857_cds <- clusterCells_Density_Peak(valid_subset_GSE72857_cds, verbose = T)

  #4. check the clusters
  plot_cell_clusters(valid_subset_GSE72857_cds, color_by = 'as.factor(Cluster)', show_density = F)
  plot_cell_clusters(valid_subset_GSE72857_cds, color_by = 'cluster', show_density = F)
  plot_cell_clusters(valid_subset_GSE72857_cds, color_by = 'cell_type', show_density = F)

  #5. also check the decision plot
  plot_rho_delta(valid_subset_GSE72857_cds, rho_threshold = 3, delta_threshold = 15 )
  plot_cell_clusters(valid_subset_GSE72857_cds, color_by = 'cell_type', show_density = F, rho_threshold = 3, delta_threshold = 15)

  #6. re-run cluster and skipping calculating the rho_sigma
  valid_subset_GSE72857_cds <- clusterCells_Density_Peak(valid_subset_GSE72857_cds, verbose = T,  rho_threshold = 3, delta_threshold = 15, skip_rho_sigma = T)

  #7. make the final clustering plot:
  plot_cell_clusters(valid_subset_GSE72857_cds, color_by = 'as.factor(Cluster)', show_density = F)
  plot_cell_clusters(valid_subset_GSE72857_cds, color_by = 'as.factor(cell_type)', show_density = F)

  #perform DEG test across clusters:
  valid_subset_GSE72857_cds@expressionFamily <- negbinomial.size()
  pData(valid_subset_GSE72857_cds)$Cluster <- factor(pData(valid_subset_GSE72857_cds)$Cluster)
  MAP_clustering_DEG_genes <- differentialGeneTest(valid_subset_GSE72857_cds, fullModelFormulaStr = '~Cluster', cores = detectCores() - 2)

  MAP_clustering_DEG_genes_subset <- MAP_clustering_DEG_genes[fData(valid_subset_GSE72857_cds)$num_cells_expressed > num_cells_expressed, ]

  #use all DEG gene from the clusters
  MAP_ordering_genes <- row.names(subset(MAP_clustering_DEG_genes, qval < qval_thrsld))

  #
  MAP_ordering_genes <- row.names(MAP_clustering_DEG_genes_subset)[order(MAP_clustering_DEG_genes_subset$qval)][1:1000]

  valid_subset_GSE72857_cds <- setOrderingFilter(valid_subset_GSE72857_cds, ordering_genes = MAP_ordering_genes)
  valid_subset_GSE72857_cds <- reduceDimension(valid_subset_GSE72857_cds, verbose = T, norm_method = 'log')
  valid_subset_GSE72857_cds <- orderCells(valid_subset_GSE72857_cds)

  plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'cell_type')
  plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'Pseudotime')
  plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'State')
  valid_subset_GSE72857_cds <- orderCells(valid_subset_GSE72857_cds, root_state = 2)

  #check the statistics
  cor(score_df$naive_pseudotime, pData(valid_subset_GSE72857_cds)[row.names(score_df), 'Pseudotime'])
  cor(score_df$naive_pseudotime, pData(valid_subset_GSE72857_cds)[row.names(score_df), 'Pseudotime'], method = 'kendall', use = 'pairwise.complete.obs')
  calClusteringMetrics(score_df$lineage, as.character(pData(valid_subset_GSE72857_cds)[row.names(score_df), 'State']))

  pdf(paste(main_fig_dir, 'monocle2_marseq_tree_sample_cells.pdf', sep = ''), height = 1.5, width = 1.5)
  plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'cell_type', cell_size = 0.5) + nm_theme()#, cell_size = pData(valid_subset_GSE72857_cds)$Pseudotime
  dev.off()

  pdf(paste(main_fig_dir, 'monocle2_marseq_tree_sample_cells_helper.pdf', sep = ''))
  plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'cell_type', cell_size = 0.5)#, cell_size = pData(valid_subset_GSE72857_cds)$Pseudotime
  dev.off()
}

ordering_genes <- row.names(subset(fData(valid_subset_GSE72857_cds), use_for_ordering == T))

dpt_res <- run_new_dpt(valid_subset_GSE72857_cds[ordering_genes, ], normalize = F)
qplot(dpt_res$dm$DC1, dpt_res$dm$DC2, color = pData(valid_subset_GSE72857_cds)$cell_type) #ensure the tree is correct
qplot(dpt_res$dm$DC1, dpt_res$dm$DC2, color = as.character(dpt_res$branch[, 1])) #ensure the tree is correct
qplot(dpt_res$dm$DC1, dpt_res$dm$DC2, color = score_df$lineage) #ensure the tree is correct

#check the statistics
cor(score_df$naive_pseudotime, dpt_res$pt)
cor(score_df$naive_pseudotime, dpt_res$pt, method = 'kendall', use = 'pairwise.complete.obs')
calClusteringMetrics(score_df$lineage, as.character(dpt_res$branch[, 1]))

slicer_res <- run_slicer(valid_subset_GSE72857_cds[ordering_genes, sample_cells])
qplot(slicer_res$traj_lle[, 1], slicer_res$traj_lle[, 2], color = pData(valid_subset_GSE72857_cds)$cell_type) #ensure the tree is correct
qplot(slicer_res$traj_lle[, 1], slicer_res$traj_lle[, 2], color = as.character(slicer_res$order_df$branches)) #ensure the tree is correct

#check the statistics
cor(score_df$naive_pseudotime, slicer_res$order_df$cells_ordered)
cor(score_df$naive_pseudotime, slicer_res$order_df$cells_ordered, method = 'kendall', use = 'pairwise.complete.obs')
calClusteringMetrics(score_df$lineage, slicer_res$order_df$branches)

# there are bad genes where naN generated: fix this (done)
monocle1_cds <- reduceDimension(valid_subset_GSE72857_cds[ordering_genes, sample_cells], norm_method = 'log', reduction_method = 'ICA', verbose = T) #sample(ordering_genes, 100) only sample 500 genes
monocle1_cds <- orderCells(monocle1_cds, num_paths = 2)
plot_cell_trajectory(monocle1_cds, color_by = 'Pseudotime')
plot_cell_trajectory(monocle1_cds, color_by = 'cell_type')
plot_spanning_tree(monocle1_cds)

cor(score_df$naive_pseudotime, pData(monocle1_cds)$Pseudotime)
cor(score_df$naive_pseudotime, pData(monocle1_cds)$Pseudotime, method = 'kendall', use = 'pairwise.complete.obs')
calClusteringMetrics(score_df$lineage, pData(monocle1_cds)$State)

monocle1_cds <- orderCells(monocle1_cds, num_paths = 2, root_state = 3)
monocle1_cds_state <- pData(monocle1_cds)$State
monocle1_cds_pseudotime <- pData(monocle1_cds)$Pseudotime

cor(score_df$naive_pseudotime, pData(monocle1_cds)$Pseudotime)
cor(score_df$naive_pseudotime, pData(monocle1_cds)$Pseudotime, method = 'kendall', use = 'pairwise.complete.obs')
calClusteringMetrics(score_df$lineage, pData(monocle1_cds)$State)

pdf(paste(main_fig_dir, 'monocle1_marseq_tree_sample_cells.pdf', sep = ''), height = 1.5, width = 1.5)
plot_cell_trajectory(monocle1_cds, color_by = 'cell_type', cell_size = 0.5) + nm_theme() + scale_color_manual(values = Mar_seq_cols)#, cell_size = pData(valid_subset_GSE72857_cds)$Pseudotime
dev.off()

# dpt_df <- data.frame(x = as.numeric(eigenvectors(diff.plot)[, 1]), y = as.numeric(eigenvectors(diff.plot)[, 2]), branch = as.character(color.branch[branching]), cell_type = pData(valid_subset_GSE72857_cds)$cell_type)
dpt_df <- data.frame(x = dpt_res$dm$DC1, y = dpt_res$dm$DC2, dpt = dpt_res$pt, branch = dpt_res$branch[, 1], cell_type = pData(valid_subset_GSE72857_cds)$cell_type)

pdf(paste(main_fig_dir, 'marseq_dpt_tree_sample_cells.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(x, y, color = cell_type, data = dpt_df, size = 0.25) + nm_theme() + scale_size(range = c(0.25, 0.25))
dev.off()

slicer_df <- data.frame(LLE1 = slicer_res$traj_lle[, 1], LLE2 = slicer_res$traj_lle[, 2],
                        pseudotime = slicer_res$order_df$cells_ordered, branch = slicer_res$order_df$branches,
                        cell_type = pData(valid_subset_GSE72857_cds[, sample_cells])$cell_type)
pdf(paste(main_fig_dir, 'marseq_slicer_tree_sample_cells.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(LLE1, LLE2, color = cell_type, data = slicer_df) + nm_theme() + scale_size(range = c(0.25, 0.25))
dev.off()

pdf(paste(main_fig_dir, 'marseq_monocle1_tree_sample_cells.pdf', sep = ''))
plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'cell_type', cell_size = 0.5)#, cell_size = pData(valid_subset_GSE72857_cds)$Pseudotime
dev.off()

data <- t(log2(exprs(valid_subset_GSE72857_cds) + 1))
write.csv(file = paste('./csv_data/Wishbone_test_data/', 'marseq_data_sample_cells', ".txt", sep = ''), as.matrix(data), quote = F, row.names = T)
row.names(subset(pData(valid_subset_GSE72857_cds), Pseudotime == 0))
# [1] "W31430"

wishbone_res_downsampling <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/MAR_seq_fractioin_wishbone_df_fig_4_downsampling.txt', header = T, sep = '\t')
qplot(dm1, dm2, data = wishbone_res_downsampling)
qplot(tSNE1, tSNE2, data = wishbone_res_downsampling, color = as.character(branch), size = 0.5)

pdf(paste(main_fig_dir, 'marseq_wishbone_tree_downsampling.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(tSNE1, tSNE2, data = wishbone_res_downsampling, color = pData(valid_subset_GSE72857_cds)$cell_type, size = 0.5) + nm_theme() + scale_size(range = c(0.5, 0.5))
# plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'cell_type', cell_size = 0.5, show_branch_points = F) + nm_theme()#, cell_size = pData(valid_subset_GSE72857_cds)$Pseudotime
dev.off()

#################################################################################################################################################################################################################################################run slicer and dpt on the data:
# use the genes from the empirical ordering
#################################################################################################################################################################################################################################################run slicer and dpt on the data:
empirical_ordering_genes <- c(mep_genes, cmp_genes, gmp_genes)
empirical_ordering_genes_ids <- row.names(subset(fData(valid_subset_GSE72857_cds), gene_short_name %in% empirical_ordering_genes))

slicer_res_empirical_ordering_genes <- run_slicer(valid_subset_GSE72857_cds[empirical_ordering_genes_ids, sample_cells],
                                                  start = which.min(pData(valid_subset_GSE72857_cds[empirical_ordering_genes_ids, sample_cells])$Pseudotime))
slicer_df_empirical_ordering <- data.frame(LLE1 = slicer_res_empirical_ordering_genes$traj_lle[, 1], LLE2 = slicer_res_empirical_ordering_genes$traj_lle[, 2],
                                           pseudotime = slicer_res_empirical_ordering_genes$order_df$cells_ordered, branch = slicer_res_empirical_ordering_genes$order_df$branches,
                                           cell_type = pData(valid_subset_GSE72857_cds[, sample_cells])$cell_type)
# pdf(paste(main_fig_dir, 'marseq_slicer_tree.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(LLE1, LLE2, color = cell_type, data = slicer_df_empirical_ordering) + nm_theme() + scale_size(range = c(0.25, 0.25))
# dev.off()

monocle1_cds_empirical_ordering_genes <- reduceDimension(valid_subset_GSE72857_cds[empirical_ordering_genes_ids, sample_cells], norm_method = 'log', reduction_method = 'ICA')
monocle1_cds_empirical_ordering_genes <- orderCells(monocle1_cds_empirical_ordering_genes, num_paths = 2)
plot_cell_trajectory(monocle1_cds_empirical_ordering_genes, color_by = 'Pseudotime')
plot_cell_trajectory(monocle1_cds_empirical_ordering_genes, color_by = 'cell_type')
plot_cell_trajectory(monocle1_cds_empirical_ordering_genes)
monocle1_cds_empirical_ordering_genes <- orderCells(monocle1_cds_empirical_ordering_genes, root_state = 3)

monocle2_cds_empirical_ordering_genes <- reduceDimension(valid_subset_GSE72857_cds[empirical_ordering_genes_ids, ], norm_method = 'log')
monocle2_cds_empirical_ordering_genes <- orderCells(monocle2_cds_empirical_ordering_genes)
plot_cell_trajectory(monocle2_cds_empirical_ordering_genes, color_by = 'Pseudotime')
plot_cell_trajectory(monocle2_cds_empirical_ordering_genes, color_by = 'cell_type')
plot_cell_trajectory(monocle2_cds_empirical_ordering_genes)
monocle2_cds_empirical_ordering_genes <- orderCells(monocle2_cds_empirical_ordering_genes, root_state = 2)

dpt_res_empirical_ordering_genes <- run_new_dpt(monocle2_cds_empirical_ordering_genes)
#################################################################################################################################################################################################################################################
#auto-correlation and lag-1 correlation, etc.
#################################################################################################################################################################################################################################################run slicer and dpt on the data:
#write data for running wishbone:
#empirical for wishbone:
data <- t(log2(exprs(valid_subset_GSE72857_cds[empirical_ordering_genes_ids, ]) + 1))
write.csv(file = paste('./csv_data/Wishbone_test_data/', 'marseq_data_sample_cells_empirical', ".txt", sep = ''), as.matrix(data), quote = F, row.names = T)
row.names(subset(pData(valid_subset_GSE72857_cds), Pseudotime == 0))

wishbone_res_downsampling_empirical <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/MAR_seq_fractioin_wishbone_df_fig_4_downsampling_empirical.txt', header = T, sep = '\t')
qplot(dm1, dm2, data = wishbone_res_downsampling_empirical)
qplot(tSNE1, tSNE2, data = wishbone_res_downsampling_empirical, color = as.character(branch), size = 0.5)

pdf(paste(main_fig_dir, 'marseq_wishbone_tree_downsampling_empirical.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(tSNE1, tSNE2, data = wishbone_res_downsampling_empirical, color = pData(valid_subset_GSE72857_cds)$cell_type, size = 0.5) + nm_theme() + scale_size(range = c(0.5, 0.5))
# plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'cell_type', cell_size = 0.5, show_branch_points = F) + nm_theme()#, cell_size = pData(valid_subset_GSE72857_cds)$Pseudotime
dev.off()

pData(valid_subset_GSE72857_cds)$reference_pseudotime <- score_df$naive_pseudotime
pData(valid_subset_GSE72857_cds)$dpt_pseudotime <- dpt_res$pt
pData(valid_subset_GSE72857_cds)$wishbone_pseudotime <- wishbone_res_downsampling$trajectory

pData(valid_subset_GSE72857_cds)$reference_state <- score_df$lineage
pData(valid_subset_GSE72857_cds)$monocle1_state <- pData(monocle1_cds)$State
pData(valid_subset_GSE72857_cds)$dpt_state <- dpt_res$branch[, 1]
pData(valid_subset_GSE72857_cds)$wishbone_state <- wishbone_res_downsampling$branch #wishbone_pseudotime$state

#################################################################################################################################################################################################################################################
# set pseudotime/state for the downsampling cells
#################################################################################################################################################################################################################################################
pData(valid_subset_GSE72857_cds)$monocle1_pseudotime <- NA
pData(valid_subset_GSE72857_cds)[colnames(monocle1_cds), 'monocle1_pseudotime'] <- monocle1_cds_pseudotime
pData(valid_subset_GSE72857_cds)$slicer_pseudotime <- NA
pData(valid_subset_GSE72857_cds)[colnames(monocle1_cds), 'slicer_pseudotime'] <- slicer_res$order_df$cells_ordered

pData(valid_subset_GSE72857_cds)$monocle1_state <- NA
pData(valid_subset_GSE72857_cds)[colnames(monocle1_cds), "monocle1_state"] <- monocle1_cds_state
pData(valid_subset_GSE72857_cds)$slicer_state <- NA
pData(valid_subset_GSE72857_cds)[colnames(monocle1_cds), 'slicer_state'] <- slicer_res$order_df$branch

table(pData(valid_subset_GSE72857_cds)[, c('dpt_state', 'cell_type')]) #state 1: CMP; state 2: erythroid 3: GMP
plot_cell_trajectory(valid_subset_GSE72857_cds) #State 4, 5: erythroid: State 1-4: GMP

#################################################################################################################################################################################################################################################
# use the reference branch assignment to check the pseudotime correlation 
#################################################################################################################################################################################################################################################
score_df_cgmp <- subset(score_df, lineage %in% c('CMP', 'GMP'))
score_df_cmep <- subset(score_df, lineage %in% c('CMP', 'MEP'))

dpt_kendall_tau <- c(cor(rank(score_df_cgmp$naive_pseudotime), rank(dpt_res$pt[row.names(score_df) %in% row.names(score_df_cgmp)]), method = 'kendall', use = 'pairwise.complete.obs'), 
                     cor(rank(score_df_cmep$naive_pseudotime), rank(dpt_res$pt[row.names(score_df) %in% row.names(score_df_cmep)]), method = 'kendall', use = 'pairwise.complete.obs'))
dpt_pearson_rho <-  c(cor(rank(score_df_cgmp$naive_pseudotime), rank(dpt_res$pt[row.names(score_df) %in% row.names(score_df_cgmp)])), 
                      cor(rank(score_df_cmep$naive_pseudotime), rank(dpt_res$pt[row.names(score_df) %in% row.names(score_df_cmep)])))

monocle1_kendall_tau <- c(cor(rank(score_df_cgmp$naive_pseudotime), rank(pData(valid_subset_GSE72857_cds[, row.names(score_df_cgmp)])$monocle1_pseudotime), method = 'kendall', use = 'pairwise.complete.obs'), 
                          cor(rank(score_df_cmep$naive_pseudotime), rank(pData(valid_subset_GSE72857_cds[, row.names(score_df_cmep)])$monocle1_pseudotime), method = 'kendall', use = 'pairwise.complete.obs'))
monocle1_pearson_rho <- c(cor(rank(score_df_cgmp$naive_pseudotime), rank(pData(valid_subset_GSE72857_cds[, row.names(score_df_cgmp)])$monocle1_pseudotime)), 
                          cor(rank(score_df_cgmp$naive_pseudotime), rank(pData(valid_subset_GSE72857_cds[, row.names(score_df_cgmp)])$monocle1_pseudotime)))

monocle2_kendall_tau <- c(cor(rank(score_df_cgmp$naive_pseudotime), rank(pData(valid_subset_GSE72857_cds[, row.names(score_df_cgmp)])$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs'), 
                          cor(rank(score_df_cmep$naive_pseudotime), rank(pData(valid_subset_GSE72857_cds[, row.names(score_df_cmep)])$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs'))
monocle2_pearson_rho <- c(cor(rank(score_df_cgmp$naive_pseudotime), rank(pData(valid_subset_GSE72857_cds[, row.names(score_df_cgmp)])$Pseudotime)), 
                          cor(rank(score_df_cmep$naive_pseudotime), rank(pData(valid_subset_GSE72857_cds[, row.names(score_df_cmep)])$Pseudotime)))

slicer_kendall_tau <- c(cor(rank(score_df_cgmp$naive_pseudotime), rank(slicer_res$order_df$cells_ordered[row.names(score_df) %in% row.names(score_df_cgmp)]), method = 'kendall', use = 'pairwise.complete.obs'), 
                        cor(rank(score_df_cmep$naive_pseudotime), rank(slicer_res$order_df$cells_ordered[row.names(score_df) %in% row.names(score_df_cmep)]), method = 'kendall', use = 'pairwise.complete.obs'))
slicer_pearson_rho <-  c(cor(rank(score_df_cgmp$naive_pseudotime), rank(slicer_res$order_df$cells_ordered[row.names(score_df) %in% row.names(score_df_cgmp)])), 
                         cor(rank(score_df_cmep$naive_pseudotime), rank(slicer_res$order_df$cells_ordered[row.names(score_df) %in% row.names(score_df_cmep)])))
# wishbone
pData(valid_subset_GSE72857_cds)$reference_state <- score_df$lineage
pData(valid_subset_GSE72857_cds)$dpt_state <- dpt_res$branch[, 1]
pData(valid_subset_GSE72857_cds)$wishbone_state <- wishbone_res_downsampling$branch #wishbone_pseudotime$state
wishbone_kendall_tau <- c(cor(rank(score_df_cgmp$naive_pseudotime), rank(wishbone_res_downsampling[row.names(score_df_cgmp), 'trajectory']), method = 'kendall', use = 'pairwise.complete.obs'), 
                          cor(rank(score_df_cmep$naive_pseudotime), rank(wishbone_res_downsampling[row.names(score_df_cmep), 'trajectory']), method = 'kendall', use = 'pairwise.complete.obs'))
wishbone_pearson_rho <- c(cor(rank(score_df_cgmp$naive_pseudotime), rank(wishbone_res_downsampling[row.names(score_df_cgmp), 'trajectory']), method = 'kendall', use = 'pairwise.complete.obs'), 
                          cor(rank(score_df_cmep$naive_pseudotime), rank(wishbone_res_downsampling[row.names(score_df_cmep), 'trajectory']), method = 'kendall', use = 'pairwise.complete.obs'))

dpt_adj_rand_ind <- calClusteringMetrics(score_df$lineage, pData(valid_subset_GSE72857_cds)$dpt_state)
monocle1_adj_rand_ind <- calClusteringMetrics(score_df[sample_cells, 'lineage'], pData(valid_subset_GSE72857_cds)$monocle1_state)
monocle2_adj_rand_ind <- calClusteringMetrics(score_df$lineage, pData(valid_subset_GSE72857_cds)$State)
slicer_adj_rand_ind <- calClusteringMetrics(score_df[sample_cells, 'lineage'], slicer_res$order_df$branches)
wishbone_adj_rand_ind <- calClusteringMetrics(score_df[sample_cells, 'lineage'], wishbone_res_downsampling$branch)

#################################################################################################################################################################################################################################################
# empirical
emp_dpt_kendall_tau <- c(cor(rank(score_df_cgmp$naive_pseudotime), rank(dpt_res_empirical_ordering_genes$pt[row.names(score_df) %in% row.names(score_df_cgmp)]), method = 'kendall', use = 'pairwise.complete.obs'), 
                         cor(rank(score_df_cmep$naive_pseudotime), rank(dpt_res_empirical_ordering_genes$pt[row.names(score_df) %in% row.names(score_df_cmep)]), method = 'kendall', use = 'pairwise.complete.obs'))
emp_dpt_pearson_rho <- c(cor(rank(score_df_cgmp$naive_pseudotime), rank(dpt_res_empirical_ordering_genes$pt[row.names(score_df) %in% row.names(score_df_cgmp)])), 
                         cor(rank(score_df_cmep$naive_pseudotime), rank(dpt_res_empirical_ordering_genes$pt[row.names(score_df) %in% row.names(score_df_cmep)])))

emp_monocle1_kendall_tau <- c(cor(rank(score_df_cgmp$naive_pseudotime), rank(pData(monocle1_cds_empirical_ordering_genes[, row.names(score_df_cgmp)])$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs'), 
                              cor(rank(score_df_cmep$naive_pseudotime), rank(pData(monocle1_cds_empirical_ordering_genes[, row.names(score_df_cmep)])$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs'))
emp_monocle1_pearson_rho <- c(cor(rank(score_df_cgmp$naive_pseudotime), rank(pData(monocle1_cds_empirical_ordering_genes[, row.names(score_df_cgmp)])$Pseudotime)), 
                              cor(rank(score_df_cmep$naive_pseudotime), rank(pData(monocle1_cds_empirical_ordering_genes[, row.names(score_df_cmep)])$Pseudotime)))

emp_monocle2_kendall_tau <- c(cor(rank(score_df_cgmp$naive_pseudotime), rank(pData(monocle2_cds_empirical_ordering_genes[, row.names(score_df_cgmp)])$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs'), 
                              cor(rank(score_df_cmep$naive_pseudotime), rank(pData(monocle2_cds_empirical_ordering_genes[, row.names(score_df_cmep)])$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs'))
emp_monocle2_pearson_rho <- c(cor(rank(score_df_cgmp$naive_pseudotime), rank(pData(monocle2_cds_empirical_ordering_genes[, row.names(score_df_cgmp)])$Pseudotime)), 
                              cor(rank(score_df_cmep$naive_pseudotime), rank(pData(monocle2_cds_empirical_ordering_genes[, row.names(score_df_cmep)])$Pseudotime)))

emp_slicer_kendall_tau <- c(cor(rank(score_df_cgmp$naive_pseudotime), rank(slicer_res_empirical_ordering_genes$order_df$cells_ordered[row.names(score_df) %in% row.names(score_df_cgmp)]), method = 'kendall', use = 'pairwise.complete.obs'), 
                            cor(rank(score_df_cmep$naive_pseudotime), rank(slicer_res_empirical_ordering_genes$order_df$cells_ordered[row.names(score_df) %in% row.names(score_df_cmep)]), method = 'kendall', use = 'pairwise.complete.obs'))
emp_slicer_pearson_rho <- c(cor(rank(score_df_cgmp$naive_pseudotime), rank(slicer_res_empirical_ordering_genes$order_df$cells_ordered[row.names(score_df) %in% row.names(score_df_cgmp)])), 
                            cor(rank(score_df_cmep$naive_pseudotime), rank(slicer_res_empirical_ordering_genes$order_df$cells_ordered[row.names(score_df) %in% row.names(score_df_cmep)])))

#wishbone
emp_wishbone_kendall_tau <- c(cor(rank(score_df_cgmp$naive_pseudotime), rank(wishbone_res_downsampling_empirical[row.names(score_df_cgmp), 'trajectory']), method = 'kendall', use = 'pairwise.complete.obs'), 
                              cor(rank(score_df_cmep$naive_pseudotime), rank(wishbone_res_downsampling_empirical[row.names(score_df_cmep), 'trajectory']), method = 'kendall', use = 'pairwise.complete.obs'))
emp_wishbone_pearson_rho <- c(cor(rank(score_df_cgmp$naive_pseudotime), rank(wishbone_res_downsampling_empirical[row.names(score_df_cgmp), 'trajectory']), method = 'kendall', use = 'pairwise.complete.obs'), 
                              cor(rank(score_df_cmep$naive_pseudotime), rank(wishbone_res_downsampling_empirical[row.names(score_df_cmep), 'trajectory']), method = 'kendall', use = 'pairwise.complete.obs'))

emp_dpt_adj_rand_ind <- calClusteringMetrics(score_df$lineage, dpt_res_empirical_ordering_genes$branch[, 1])
emp_monocle1_adj_rand_ind <- calClusteringMetrics(score_df[sample_cells, 'lineage'], pData(monocle1_cds_empirical_ordering_genes)$State)
emp_monocle2_adj_rand_ind <- calClusteringMetrics(score_df$lineage, pData(monocle2_cds_empirical_ordering_genes)$State)
emp_slicer_adj_rand_ind <- calClusteringMetrics(score_df[sample_cells, 'lineage'], slicer_res_empirical_ordering_genes$order_df$branches)
emp_wishbone_adj_rand_ind <- calClusteringMetrics(score_df[sample_cells, 'lineage'], wishbone_res_downsampling_empirical$branch)

cor_df <- data.frame("Pearson_cor" = c(monocle2_pearson_rho, monocle1_pearson_rho, dpt_pearson_rho, slicer_pearson_rho, wishbone_pearson_rho, #
                                       emp_monocle2_pearson_rho, emp_monocle1_pearson_rho, emp_dpt_pearson_rho, emp_slicer_pearson_rho, emp_wishbone_pearson_rho), #emp_monocle1_pearson_rho wishbone
                     "Kendall_cor" = c(monocle2_kendall_tau, monocle1_kendall_tau, dpt_kendall_tau, slicer_kendall_tau, wishbone_kendall_tau,
                                       emp_monocle2_kendall_tau, emp_monocle1_kendall_tau, emp_dpt_kendall_tau, emp_slicer_kendall_tau, emp_wishbone_kendall_tau),
                     # "Time_kcor" = c(Time_kcor1, Time_kcor2, Time_kcor3, Time_kcor4, Time_kcor5, Time_kcor6),
                     names = c("monocle2 (all genes)", "monocle1 (all genes)", "dpt (all genes)", "slicer (all genes)", "wishbone (all genes)",
                               "monocle2 (markers)", "monocle1 (markers)", "dpt (markers)", "slicer (markers)", "wishbone (markers)"))

mlt_cor_df <- melt(cor_df, id.vars = c('names'))

pdf(paste(main_fig_dir, 'mar_seq_pseudotime_correlation_sample_cells.pdf', sep = ''), height = 1.5, width = 3.2)
ggplot(aes(names, abs(value)), data = mlt_cor_df) + geom_bar(stat="identity", aes(fill = variable)) + facet_wrap(~variable) + nm_theme() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_size(range = c(0.2, 0.2)) + xlab('') + ylab('') + monocle_theme_opts()
dev.off()

adj_rand_ind_df <- data.frame("Adjusted Rand index" = c(monocle2_adj_rand_ind$randIndex[3], monocle1_adj_rand_ind$randIndex[3], dpt_adj_rand_ind$randIndex[3], slicer_adj_rand_ind$randIndex[3], wishbone_adj_rand_ind$randIndex[3],
                                                        emp_monocle2_adj_rand_ind$randIndex[3], emp_monocle1_adj_rand_ind$randIndex[3], emp_dpt_adj_rand_ind$randIndex[3], emp_slicer_adj_rand_ind$randIndex[3], emp_wishbone_adj_rand_ind$randIndex[3]), #wishbone ,
                              # "Time_kcor" = c(Time_kcor1, Time_kcor2, Time_kcor3, Time_kcor4, Time_kcor5, Time_kcor6),
                              names = c("monocle2 (all genes)", "monocle1 (all genes)", "dpt (all genes)", "slicer (all genes)", "wishbone (all genes)",
                                        "monocle2 (markers)", "monocle1 (markers)", "dpt (markers)", "slicer (markers)", "wishbone (markers)"))

mlt_adj_rand_ind_df <- melt(adj_rand_ind_df, id.vars = c('names'))

pdf(paste(main_fig_dir, 'mar_seq_adj_rand_index_sample_cells.pdf', sep = ''), height = 1.5, width = 1.6)
ggplot(aes(names, abs(value)), data = mlt_adj_rand_ind_df) + geom_bar(stat="identity", aes(fill = variable)) + facet_wrap(~variable) + nm_theme() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_size(range = c(0.2, 0.2)) + xlab('') + ylab('') + monocle_theme_opts()
dev.off()

#################################################################################################################################################################################################################################################
# only look at the all gene section
#################################################################################################################################################################################################################################################
mlt_cor_df_subset <- mlt_cor_df[c(1:5, 11:15, 21:25, 31:35), ]

mlt_cor_df_subset$names <- str_split_fixed(mlt_cor_df_subset$names, ' ', 2)[, 1]
mlt_cor_df_subset$names <- factor(mlt_cor_df_subset$names, levels = software_levels)

mlt_cor_df_subset2 <-ddply(mlt_cor_df_subset, .(names, variable), function(x) mean(abs(x$value)))

pdf(paste(main_fig_dir, 'mar_seq_pseudotime_correlation_sample_cells_subset.pdf', sep = ''), height = 1.5, width = 3.2)
ggplot(aes(names, abs(V1)), data = mlt_cor_df_subset2) + geom_bar(stat="identity", aes(fill = names)) + facet_wrap(~variable) + nm_theme() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_size(range = c(0.2, 0.2)) + xlab('') + ylab('') + monocle_theme_opts() + scale_fill_manual(values=software_custom_color_scale)
dev.off()

mlt_adj_rand_ind_df_subset <- mlt_adj_rand_ind_df[1:5, ]
mlt_adj_rand_ind_df_subset$names <- str_split_fixed(mlt_adj_rand_ind_df_subset$names, ' ', 2)[, 1]
mlt_adj_rand_ind_df_subset$names <- factor(mlt_adj_rand_ind_df_subset$names, levels = software_levels)

pdf(paste(main_fig_dir, 'mar_seq_adj_rand_index_sample_cells_subset.pdf', sep = ''), height = 1.5, width = 1.6)
ggplot(aes(names, abs(value)), data = mlt_adj_rand_ind_df_subset) + geom_bar(stat="identity", aes(fill = names)) + facet_wrap(~variable) + nm_theme() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_size(range = c(0.2, 0.2)) + xlab('') + ylab('') + monocle_theme_opts() + scale_fill_manual(values=software_custom_color_scale)
dev.off()

#################################################################################################################################################################################################################################################
mlt_cor_df_subset <- mlt_cor_df[c(6:10, 16:20, 26:30), ]

mlt_cor_df_subset$names <- str_split_fixed(mlt_cor_df_subset$names, ' ', 2)[, 1]
mlt_cor_df_subset$names <- factor(mlt_cor_df_subset$names, levels = software_levels)

mlt_cor_df_subset2 <-ddply(mlt_cor_df_subset, .(names, variable), function(x) mean(abs(x$value)))

pdf(paste(main_fig_dir, 'mar_seq_pseudotime_correlation_sample_cells_subset2.pdf', sep = ''), height = 1.5, width = 1.6)
ggplot(aes(names, abs(value)), data = mlt_cor_df_subset) + geom_bar(stat="identity", aes(fill = names)) + facet_wrap(~variable) + nm_theme() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_size(range = c(0.2, 0.2)) + xlab('') + ylab('') + monocle_theme_opts() + scale_fill_manual(values=software_custom_color_scale)
dev.off()

mlt_adj_rand_ind_df_subset <- mlt_adj_rand_ind_df[6:10, ]
mlt_adj_rand_ind_df_subset$names <- str_split_fixed(mlt_adj_rand_ind_df_subset$names, ' ', 2)[, 1]
mlt_adj_rand_ind_df_subset$names <- factor(mlt_adj_rand_ind_df_subset$names, levels = software_levels)

pdf(paste(main_fig_dir, 'mar_seq_adj_rand_index_sample_cells_subset2.pdf', sep = ''), height = 1.5, width = 1)
ggplot(aes(names, abs(value)), data = mlt_adj_rand_ind_df_subset) + geom_bar(stat="identity", aes(fill = names)) + facet_wrap(~variable) + nm_theme() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_size(range = c(0.2, 0.2)) + xlab('') + ylab('') + monocle_theme_opts() + scale_fill_manual(values=software_custom_color_scale)
dev.off()

#################################################################################################################################################################################################################################################
# calculate the smoothness of the data:
# run only the marker genes
################################################################################################################################################################################################################################################run slicer and dpt on the data:

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

reference_GMP_full_model_fits <- cal_smoothness_autocorrelation(valid_subset_GSE72857_cds[c(mep_genes, gmp_genes, cmp_genes), pData(valid_subset_GSE72857_cds)$reference_state %in% c('CMP', 'GMP')],
                                                                modelFormulaStr="~sm.ns(reference_pseudotime)", order_by = 'reference_pseudotime')
reference_MEP_full_model_fits <- cal_smoothness_autocorrelation(valid_subset_GSE72857_cds[c(mep_genes, gmp_genes, cmp_genes), pData(valid_subset_GSE72857_cds)$reference_state %in% c('CMP', 'MEP')],
                                                                modelFormulaStr="~sm.ns(reference_pseudotime)", order_by = 'reference_pseudotime')

monocle2_GMP_full_model_fits <- cal_smoothness_autocorrelation(valid_subset_GSE72857_cds[c(mep_genes, gmp_genes, cmp_genes),
                                                                                         pData(valid_subset_GSE72857_cds)$reference_state %in% c('CMP', 'GMP')],
                                                               modelFormulaStr="~sm.ns(Pseudotime)", order_by = 'Pseudotime')
monocle2_MEP_full_model_fits <- cal_smoothness_autocorrelation(valid_subset_GSE72857_cds[c(mep_genes, gmp_genes, cmp_genes),
                                                                                         pData(valid_subset_GSE72857_cds)$reference_state %in% c('CMP', 'MEP')],
                                                               modelFormulaStr="~sm.ns(Pseudotime)", order_by = 'Pseudotime')

dpt_GMP_full_model_fits <- cal_smoothness_autocorrelation(valid_subset_GSE72857_cds[c(mep_genes, gmp_genes, cmp_genes),
                                                                                    pData(valid_subset_GSE72857_cds)$reference_state %in% c('CMP', 'GMP')],
                                                          modelFormulaStr="~sm.ns(dpt_pseudotime)", order_by = 'dpt_pseudotime')
dpt_MEP_full_model_fits <- cal_smoothness_autocorrelation(valid_subset_GSE72857_cds[c(mep_genes, gmp_genes, cmp_genes),
                                                                                    pData(valid_subset_GSE72857_cds)$reference_state %in% c('CMP', 'MEP')],
                                                          modelFormulaStr="~sm.ns(dpt_pseudotime)", order_by = 'dpt_pseudotime')

wishbone_GMP_full_model_fits <- cal_smoothness_autocorrelation(valid_subset_GSE72857_cds[c(mep_genes, gmp_genes, cmp_genes), 
                                                                                         pData(valid_subset_GSE72857_cds)$reference_state %in% c('CMP', 'GMP')],
                                                               modelFormulaStr="~sm.ns(wishbone_pseudotime)", order_by = 'wishbone_pseudotime')
wishbone_MEP_full_model_fits <- cal_smoothness_autocorrelation(valid_subset_GSE72857_cds[c(mep_genes, gmp_genes, cmp_genes), 
                                                                                         pData(valid_subset_GSE72857_cds)$reference_state %in% c('CMP', 'MEP')],
                                                               modelFormulaStr="~sm.ns(wishbone_pseudotime)", order_by = 'wishbone_pseudotime')

monocle1_GMP_full_model_fits <- cal_smoothness_autocorrelation(valid_subset_GSE72857_cds[c(mep_genes, gmp_genes, cmp_genes),
                                                                                         pData(valid_subset_GSE72857_cds)$reference_state %in% c('CMP', 'GMP')],
                                                               modelFormulaStr="~sm.ns(monocle1_pseudotime)", order_by = 'monocle1_pseudotime')
monocle1_MEP_full_model_fits <- cal_smoothness_autocorrelation(valid_subset_GSE72857_cds[c(mep_genes, gmp_genes, cmp_genes),
                                                                                         pData(valid_subset_GSE72857_cds)$reference_state %in% c('CMP', 'MEP')],
                                                               modelFormulaStr="~sm.ns(monocle1_pseudotime)", order_by = 'monocle1_pseudotime')

slicer_GMP_full_model_fits <- cal_smoothness_autocorrelation(valid_subset_GSE72857_cds[c(mep_genes, gmp_genes, cmp_genes),
                                                                                       pData(valid_subset_GSE72857_cds)$reference_state %in% c('CMP', 'GMP')],
                                                             modelFormulaStr="~sm.ns(slicer_pseudotime)", order_by = 'slicer_pseudotime')
slicer_MEP_full_model_fits <- cal_smoothness_autocorrelation(valid_subset_GSE72857_cds[c(mep_genes, gmp_genes, cmp_genes),
                                                                                       pData(valid_subset_GSE72857_cds)$reference_state %in% c('CMP', 'MEP')],
                                                             modelFormulaStr="~sm.ns(slicer_pseudotime)", order_by = 'slicer_pseudotime')

clusters_genes <- c(rep(1, length(mep_genes)), rep(2, length(gmp_genes)), rep(3, length(cmp_genes)))
cluster_autocorrelation_df <- data.frame(reference = c(reference_GMP_full_model_fits$autocorrelation, reference_MEP_full_model_fits$autocorrelation),
                                         monocle2 = c(monocle2_GMP_full_model_fits$autocorrelation, monocle2_MEP_full_model_fits$autocorrelation),
                                         wishbone = c(wishbone_GMP_full_model_fits$autocorrelation, wishbone_MEP_full_model_fits$autocorrelation),
                                         monocle1 = c(monocle1_GMP_full_model_fits$autocorrelation, monocle1_MEP_full_model_fits$autocorrelation),
                                         slicer = c(slicer_GMP_full_model_fits$autocorrelation, slicer_MEP_full_model_fits$autocorrelation),
                                         dpt = c(dpt_GMP_full_model_fits$autocorrelation, dpt_MEP_full_model_fits$autocorrelation),
                                         type = "cluster_autocorrelation", lineage = rep(c('GMP', 'MEP')))

cluster_smoothness_df <- data.frame(reference = c(reference_GMP_full_model_fits$smoothness, reference_MEP_full_model_fits$smoothness),
                                    monocle2 = c(monocle2_GMP_full_model_fits$smoothness, monocle2_MEP_full_model_fits$smoothness),
                                    wishbone = c(wishbone_GMP_full_model_fits$smoothness, wishbone_MEP_full_model_fits$smoothness),
                                    monocle1 = c(monocle1_GMP_full_model_fits$smoothness, monocle1_MEP_full_model_fits$smoothness),
                                    slicer = c(slicer_GMP_full_model_fits$smoothness, slicer_MEP_full_model_fits$smoothness),
                                    dpt = c(dpt_GMP_full_model_fits$smoothness, dpt_MEP_full_model_fits$smoothness),
                                    type = "cluster_smoothness", lineage = rep(c('GMP', 'MEP')))

cluster_residual_df <- data.frame(reference = c(reference_GMP_full_model_fits$residual, reference_MEP_full_model_fits$residual),
                                  monocle2 = c(monocle2_GMP_full_model_fits$residual, monocle2_MEP_full_model_fits$residual),
                                  wishbone = c(wishbone_GMP_full_model_fits$residual, wishbone_MEP_full_model_fits$residual),
                                  monocle1 = c(monocle1_GMP_full_model_fits$residual, monocle1_MEP_full_model_fits$residual),
                                  slicer = c(slicer_GMP_full_model_fits$residual, slicer_MEP_full_model_fits$residual),
                                  dpt = c(dpt_GMP_full_model_fits$residual, dpt_MEP_full_model_fits$residual),
                                  type = "cluster_residual", lineage = rep(c('GMP', 'MEP')))

cluster_residual_df[, 1:6] <- cluster_residual_df[, 1:6] - cluster_residual_df[, 1] 

tightness_smoothness_res <- plot_tightness_smoothness(cluster_autocorrelation_df[, -ncol(cluster_autocorrelation_df)], 
                                                      cluster_smoothness_df[, -ncol(cluster_autocorrelation_df)], 
                                                      cluster_residual_df[, -ncol(cluster_autocorrelation_df)])
ggplot(aes(x = variable, y = value), data = tightness_smoothness_res$mlt_ts_df) + geom_boxplot(aes(color = variable))  + monocle_theme_opts() +
  facet_wrap(~type, scale = 'free') + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab('') + ylab('')

pdf(paste(main_fig_dir, 'smoothness_tightness_sample_cells.pdf', sep = ''), height = 2, width = 3)
ggplot(aes(x = variable, y = value), data = subset(tightness_smoothness_res$mlt_ts_df, value < 50, variable != 'slicer')) + geom_boxplot(aes(color = variable))  + monocle_theme_opts() +
  facet_wrap(~type, scale = 'free') + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab('') + ylab('')
dev.off()

#only look at the lag-1 autocorrelation:
subset_tightness_smoothness_res <- subset(tightness_smoothness_res$mlt_ts_df, type == 'cluster_autocorrelation')
subset_tightness_smoothness_res$variable <- as.character(subset_tightness_smoothness_res$variable)
subset_tightness_smoothness_res$variable <- factor(subset_tightness_smoothness_res$variable, levels = c('reference', software_levels))

pdf(paste(main_fig_dir, 'autocorrelation_sample_cells.pdf', sep = ''), height = 1.5, width = 1.5)
ggplot(aes(x = variable, y = abs(value) ), data = subset_tightness_smoothness_res) + geom_boxplot(aes(color = variable))  + monocle_theme_opts() +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab('') + ylab('lag1_autocorrelation') + scale_color_manual(values = software_custom_color_scale)
dev.off()

#################################################################################################################################################################################################################################################
# save the data
#################################################################################################################################################################################################################################################
save.image('./RData/marseq_downsamling_empirical_ordering.RData')
