rm(list = ls())

################################################################################################################################################################################################################################################
## load necessary packages
################################################################################################################################################################################################################################################
library(stringr)
library(xacHelper)
library(grid)
library(devtools)
library(monocle)
library(plyr)
library(dplyr)
library(SLICER)
library(pheatmap)
library(UpSetR)
library(piano)
library(HSMMSingleCell)
library(destiny)
library(geomorph)
library(colorRamps)
library(L1Graph)
library(diffusionMap)
library(mcclust)
library(reshape2)

load('./RData/fig1.RData')

source('./scripts/function.R')

bk <- seq(-3.1,3.1, by=0.1)
hmcols <-  blue2green2red(length(bk) - 1)
cell_fate_color <- c('F1' = "#EA5360", 'F2' = "#7890C6") #F1 (MYOG+1): #EA5360; F2 (MYOG+1): #7890C6

# ensure HSMM ordering is correct!
HSMM_myo <- orderCells(HSMM_myo)
plot_cell_trajectory(HSMM_myo)
plot_cell_trajectory(HSMM_myo, color_by = 'Pseudotime')
plot_cell_trajectory(HSMM_myo, color_by = 'Time')

HSMM <- detectGenes(HSMM, min_expr = 1)

muscle_marker_genes <- c("CDK1", "H19", "ANPEP", "ENO3", "TNNT2", "MYH2", "DMD", "MEF2C", "ID3", "MYF5", "TMEM8C", "MYOG")

cell_fate_color <- c('Myoblast' = "#EA5360", 'Fibroblast' = "#7890C6") #EA5360; F2 (MYOG+1): #7890C6
pdf(paste(SI_fig_dir, 'HSMM_fibro_clusters.pdf', sep = ''), width=2.25, height=3.5)
plot_cell_clusters(HSMM, color_by="CellType", markers=c("MYOD1", "MYOG", "MYF5", "ANPEP")) + nm_theme() + scale_color_manual(values=cell_fate_color)
dev.off()

########################################################################################################################################################
##  trajectory plot 
# create SI4b
########################################################################################################################################################

pdf(paste(SI_fig_dir, 'SI4b.pdf', sep = ''), height=1.7, width=6)
plot_genes_positive_cells(HSMM[row.names(subset(fData(HSMM), gene_short_name %in% muscle_marker_genes)),], grouping="CellType", nrow = 2, ncol = 6) +
  theme(legend.position="none") + nm_theme() + scale_fill_manual(values=cell_fate_color)
dev.off()

MYOG_expr <- colSums(exprs(HSMM_myo[row.names(subset(fData(HSMM_myo), gene_short_name %in% c("MYOG"))),]))

pData(HSMM_myo)$MYOG <- factor("MYOG-", levels=c("MYOG-", "MYOG+"))
pData(HSMM_myo)$MYOG[as.vector(MYOG_expr > 0)] <- "MYOG+"

########################################################################################################################################################
##  trajectory plot 
# create SI4c and SI4e
########################################################################################################################################################

pdf(paste(SI_fig_dir, 'SI4d.pdf', sep = ''), width=2.25, height=3.5)
plot_cell_trajectory(HSMM_myo, color_by="MYOG", show_branch_points=F) +
  facet_wrap(~MYOG, ncol=1)  +
  stat_density2d(aes(color=MYOG), alpha=I(0.25), size=I(0.25)) +
  scale_color_manual(values=c("#7890C6", "red")) + nm_theme()

dev.off()

pdf(paste(SI_fig_dir, 'SI4d_helper.pdf', sep = ''), width=2.25, height=3.5)
plot_cell_trajectory(HSMM_myo, color_by="MYOG", show_branch_points=F) +
  facet_wrap(~MYOG, ncol=1)  +
  stat_density2d(aes(color=MYOG), alpha=I(0.25), size=I(0.25)) +
  scale_color_manual(values=c("#7890C6", "red"))

dev.off()

HSMM_expressed_genes <- row.names(subset(fData(HSMM_myo), num_cells_expressed >=5))
HSMM_BEAM_res <- BEAM(HSMM_myo[HSMM_expressed_genes,],
                      branch_point=1,
                      cores = detectCores() - 2)

head(HSMM_BEAM_res)

########################################################################################################################################################
##  branched kinetic curves
# create SI4c and SI4e
########################################################################################################################################################
pdf(paste(SI_fig_dir, 'SI4c.pdf', sep = ''), height=1.7, width=6)
plot_genes_branched_pseudotime(HSMM_myo[row.names(subset(fData(HSMM_myo), gene_short_name %in% muscle_marker_genes)),],
                               branch_point=1,
                               #lineage_states=NULL,
                               color_by="as.factor(Time)",
                               branch_labels=c("M1", "M2"),
                               #add_pval=TRUE,
                               nrow = 2, ncol = 6,
                               cell_size=0.65,
                               min_expr=0.1) + scale_color_brewer(palette="Set1") +
  theme(legend.position="none") + nm_theme()
dev.off()

num_hsmm_beam_clusters <- 5

dim(subset(HSMM_BEAM_res, qval < 0.1))
pdf(paste(SI_fig_dir, 'SI4e.pdf', sep = ''), height=3.5, width=3.5)
ph_res <- plot_genes_branched_heatmap(HSMM_myo[row.names(subset(HSMM_BEAM_res, qval < 0.1)),],
                                      num_clusters = num_hsmm_beam_clusters,
                                      branch_point=1,
                                      branch_labels = c("M1", "M2"),
                                      branch_colors = c('#979797', '#7990C8', '#F05662'), 
                                      cores = detectCores(),
                                      use_gene_short_name = T,
                                      show_rownames = F,
                                      hmcols=hmcols,
                                      return_heatmap=TRUE)
dev.off()


tree_row <- ph_res$ph_res$tree_row
tree_row_gene_clusters <- cutree(tree_row, num_hsmm_beam_clusters)


sig_HSMM_BEAM_res <- HSMM_BEAM_res[names(tree_row_gene_clusters),]
sig_HSMM_BEAM_res$beam_cluster <- tree_row_gene_clusters
write.table(sig_HSMM_BEAM_res, "HSMM_BEAM_res.txt", quote=F, col.names=T, row.names=T, sep="\t")

gsa_gene_ids = union(names(tree_row_gene_clusters), row.names(subset(fData(HSMM_myo), biotype == "protein_coding")))
#gsa_gene_ids = names(tree_row_gene_clusters)
names(tree_row_gene_clusters) <- fData(HSMM_myo[names(tree_row_gene_clusters), ])$gene_short_name

# no terms
HSMM_BEAM_reactome_gsa_hyper_heatmap <- collect_gsa_hyper_results(HSMM_myo[,], human_reactome_gsc, tree_row_gene_clusters) #gsa_gene_ids
pdf ("HSMM_BEAM_reactome_gsa_hyper_heatmap.pdf", height=36, width=16)
plot_gsa_hyper_heatmap(HSMM_myo, HSMM_BEAM_reactome_gsa_hyper_heatmap, significance=0.1)
dev.off()

HSMM_BEAM_go_bp_gsa_hyper_heatmap <- collect_gsa_hyper_results(HSMM_myo[gsa_gene_ids,], human_go_gsc, tree_row_gene_clusters)
pdf ("HSMM_BEAM_go_bp_gsa_hyper_heatmap", height=36, width=8)
plot_gsa_hyper_heatmap(HSMM_myo, HSMM_BEAM_go_bp_gsa_hyper_heatmap, significance=0.1)
dev.off()

########################################################################################################################################################
##  branched kinetic curves
# create SI4c and SI4e
########################################################################################################################################################

#slicer wishbone dpt on the muscle data
HSMM_myo <- setOrderingFilter(HSMM_myo, HSMM_myo_ordering_genes)
slicer_res <- run_slicer(HSMM_myo[HSMM_myo_ordering_genes, ], start = which.min(pData(HSMM_myo)$Pseudotime))

pdf(paste(SI_fig_dir, 'HSMM_slicer_branch.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(slicer_res$traj_lle[, 1], slicer_res$traj_lle[, 2], color = as.character(slicer_res$order_df$branches)) + nm_theme() + xlab('LLE 1') + ylab('LLE 2')
dev.off()
pdf(paste(SI_fig_dir, 'HSMM_slicer_time.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(slicer_res$traj_lle[, 1], slicer_res$traj_lle[, 2], color = pData(HSMM_myo)$Time) + nm_theme() + xlab('LLE 1') + ylab('LLE 2')
dev.off()

data <- t(log2(exprs(HSMM_myo[HSMM_myo_ordering_genes, ]) + 1))
write.csv(file = paste('./csv_data/Wishbone_test_data/', 'HSMM_myo_fig_SI2', ".txt", sep = ''), as.matrix(data), quote = F, row.names = T)
row.names(subset(pData(HSMM_myo), Pseudotime == 0))

wishbone_res_downsampling_empirical <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/MAR_seq_fractioin_wishbone_df_fig_4_downsampling_empirical.txt', header = T, sep = '\t')
qplot(dm1, dm2, data = wishbone_res_downsampling_empirical)
qplot(tSNE1, tSNE2, data = wishbone_res_downsampling_empirical, color = as.character(branch), size = 0.5)

HSMM_myo <- setOrderingFilter(HSMM_myo, HSMM_myo_ordering_genes)
dpt_res <- run_new_dpt(HSMM_myo[HSMM_myo_ordering_genes, ], normalize = F)

pdf(paste(SI_fig_dir, 'HSMM_dpt_time.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(dpt_res$dm$DC1, dpt_res$dm$DC2, color = pData(HSMM_myo)$Time) + nm_theme() + xlab('DM 1') + ylab('DM 2')
dev.off()

pdf(paste(SI_fig_dir, 'HSMM_dpt_branch.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(dpt_res$dm$DC1, dpt_res$dm$DC2, color = as.character(dpt_res$branch[, 1])) + nm_theme() + xlab('DM 1') + ylab('DM 2')
dev.off()

################################################################################################################################################################################################################################################
## use DDRTree with different dimension reduction methods
# create SI15a_1
################################################################################################################################################################################################################################################
DDRTree_pca_HSMM_myo <- reduceDimension(HSMM_myo, initial_method = PCA, verbose = T)
DDRTree_pca_HSMM_myo <- orderCells(DDRTree_pca_HSMM_myo)

pdf(paste(SI_fig_dir, 'pca_DDRTree.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(DDRTree_pca_HSMM_myo, color_by = 'Time', cell_size = 0.5, show_branch_points = F) + nm_theme()
dev.off()

DDRTree_ica_HSMM_myo <- reduceDimension(HSMM_myo, initial_method = ICA, verbose = T)
DDRTree_ica_HSMM_myo <- orderCells(DDRTree_ica_HSMM_myo)
DDRTree_ica_HSMM_myo <- orderCells(DDRTree_ica_HSMM_myo, root_state = 2)
pdf(paste(SI_fig_dir, 'ica_DDRTree.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(DDRTree_ica_HSMM_myo, color_by = 'Time', cell_size = 0.5, show_branch_points = F) + nm_theme()
dev.off()

DDRTree_lle_HSMM_myo <- reduceDimension(HSMM_myo, initial_method = LLE2, verbose = T)
DDRTree_lle_HSMM_myo <- orderCells(DDRTree_lle_HSMM_myo)

pdf(paste(SI_fig_dir, 'lle_DDRTree.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(DDRTree_lle_HSMM_myo, color_by = 'Time', cell_size = 0.5, show_branch_points = F) + nm_theme()
dev.off()

DDRTree_dm_HSMM_myo <- reduceDimension(HSMM_myo, initial_method = DM, verbose = T)
DDRTree_dm_HSMM_myo <- orderCells(DDRTree_dm_HSMM_myo)

pdf(paste(SI_fig_dir, 'dm_DDRTree.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(DDRTree_dm_HSMM_myo, color_by = 'Time', cell_size = 0.5, show_branch_points = F) + nm_theme()
dev.off()

DDRTree_dm_HSMM_myo <- trimTree(DDRTree_dm_HSMM_myo)

DDRTree_destiny_HSMM_myo <- reduceDimension(HSMM_myo, initial_method = destiny_diffusionMaps, verbose = T)
DDRTree_destiny_HSMM_myo <- orderCells(DDRTree_destiny_HSMM_myo)

pdf(paste(SI_fig_dir, 'destiny_DDRTree.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(DDRTree_destiny_HSMM_myo, color_by = 'Time', cell_size = 0.5, show_branch_points = F) + nm_theme()
dev.off()

DDRTree_destiny_HSMM_myo <- trimTree(DDRTree_destiny_HSMM_myo)

################################################################################################################################################################################################################################################
## use simplePPT with different dimension reduction methods
# create SI15a_2
################################################################################################################################################################################################################################################
params.lambda <- 2
params.bandwidth <- 5

pca_res <- PCA(as.matrix(log(exprs(HSMM_myo[HSMM_myo_ordering_genes, ]) + 1)))
qplot(pca_res[, 1], pca_res[, 2], color = pData(HSMM_myo)$Time)
HSMM_pca_res <- principal_tree(t(pca_res[, 1:2]), MU = NULL, lambda = params.lambda, bandwidth = params.bandwidth, maxIter = 10, verbose = T)
pca_HSMM_myo <- reduceDimension(HSMM_myo, reduction_method = 'SimplePPT', lambda = params.lambda, bandwidth = params.bandwidth, initial_method = PCA, verbose = T)
pca_HSMM_myo <- orderCells(pca_HSMM_myo)

pdf(paste(SI_fig_dir, 'pca_simplePPT.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(pca_HSMM_myo, color_by = 'Time', cell_size = 0.5, show_branch_points = F) + nm_theme()
dev.off()

ica_res <- ICA(as.matrix(log2(exprs(HSMM_myo[HSMM_myo_ordering_genes, ]) + 1)))
HSMM_ica_res <- principal_tree(t(ica_res[, 1:2]), MU = NULL, lambda = params.lambda, bandwidth = 15, maxIter = 10, verbose = T)
ica_HSMM_myo <- reduceDimension(HSMM_myo, reduction_method = 'SimplePPT', lambda = params.lambda, bandwidth = params.bandwidth, initial_method = ICA, verbose = T)
ica_HSMM_myo <- orderCells(ica_HSMM_myo)

ica_HSMM_myo@reducedDimS <- ica_HSMM_myo@reducedDimS * 1000
pdf(paste(SI_fig_dir, 'ica_simplePPT.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(ica_HSMM_myo, color_by = 'Time', cell_size = 0.5, show_branch_points = F) + nm_theme()
dev.off()

lle_res <- LLE(as.matrix(log(exprs(HSMM_myo[HSMM_myo_ordering_genes, ]) + 1)))
HSMM_lle_res <- principal_tree(t(lle_res[, 1:2]), MU = NULL, lambda = params.lambda, bandwidth = params.bandwidth, maxIter = 10, verbose = T)
lle_HSMM_myo <- reduceDimension(HSMM_myo, reduction_method = 'SimplePPT', lambda = params.lambda, bandwidth = params.bandwidth, initial_method = LLE, verbose = T)
lle_HSMM_myo <- orderCells(lle_HSMM_myo)

lle_HSMM_myo@reducedDimS <- lle_HSMM_myo@reducedDimS * 1000
pdf(paste(SI_fig_dir, 'lle_simplePPT.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(lle_HSMM_myo, color_by = 'Time', cell_size = 0.5, show_branch_points = F) + nm_theme()
dev.off()

dm_res <- DM(as.matrix(log(exprs(HSMM_myo[HSMM_myo_ordering_genes, ]) + 1)))
HSMM_dm_res <- principal_tree(t(dm_res[, 1:2]), MU = NULL, lambda = params.lambda, bandwidth = params.bandwidth, maxIter = 10, verbose = T)
destiny_res <- destiny_diffusionMaps(as.matrix(log(exprs(HSMM_myo[HSMM_myo_ordering_genes, ]) + 1)))
HSMM_destiny_res <- principal_tree(t(destiny_res[, 1:2]), MU = NULL, lambda = params.lambda, bandwidth = params.bandwidth, maxIter = 10, verbose = T)
dm_HSMM_myo <- reduceDimension(HSMM_myo, reduction_method = 'SimplePPT', lambda = params.lambda, bandwidth = params.bandwidth, initial_method = DM, verbose = T)
dm_HSMM_myo <- orderCells(dm_HSMM_myo)

dm_HSMM_myo@reducedDimS <- dm_HSMM_myo@reducedDimS * 1000
dm_HSMM_myo@reducedDimK <- dm_HSMM_myo@reducedDimK * 1000

x_min <- min(dm_HSMM_myo@reducedDimS[1, ]); x_max <- max(dm_HSMM_myo@reducedDimS[1, ])
y_min <- min(dm_HSMM_myo@reducedDimS[2, ]); y_max <- max(dm_HSMM_myo@reducedDimS[2, ])
pdf(paste(SI_fig_dir, 'dm_simplePPT.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(dm_HSMM_myo, color_by = 'Time', cell_size = 0.5, show_branch_points = F) + nm_theme() +
  scale_x_continuous(limits = c(x_min, x_max), breaks = round(seq(x_min, x_max, length.out = 4), digits = 2)) +
  scale_y_continuous(limits = c(y_min, y_max), breaks = round(seq(y_min, y_max, length.out = 4), digits = 2)) # + theme(axis.text.x = element_text(angle = 30, hjust = 0), axis.text.y = element_text(angle = 30, hjust = 0))
dev.off()

destiny_res <- destiny_diffusionMaps(as.matrix(log(exprs(HSMM_myo[HSMM_myo_ordering_genes, ]) + 1)))
HSMM_destiny_res <- principal_tree(t(destiny_res[, 1:2]), MU = NULL, lambda = params.lambda, bandwidth = params.bandwidth, maxIter = 10, verbose = T)
destiny_HSMM_myo <- reduceDimension(HSMM_myo, reduction_method = 'SimplePPT', lambda = params.lambda, bandwidth = params.bandwidth, initial_method = destiny_diffusionMaps, verbose = T)
destiny_HSMM_myo <- orderCells(destiny_HSMM_myo)

destiny_HSMM_myo@reducedDimS <- destiny_HSMM_myo@reducedDimS * 2000
destiny_HSMM_myo@reducedDimK <- destiny_HSMM_myo@reducedDimK * 2000
x_min <- min(destiny_HSMM_myo@reducedDimS[1, ]); x_max <- max(destiny_HSMM_myo@reducedDimS[1, ])
y_min <- min(destiny_HSMM_myo@reducedDimS[2, ]); y_max <- max(destiny_HSMM_myo@reducedDimS[2, ])
pdf(paste(SI_fig_dir, 'destiny_simplePPT.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(destiny_HSMM_myo, color_by = 'Time', cell_size = 0.5, show_branch_points = F) + nm_theme() +
  scale_x_continuous(limits = c(x_min, x_max), breaks = round(seq(x_min, x_max, length.out = 4), digits = 2)) +
  scale_y_continuous(limits = c(y_min, y_max), breaks = round(seq(y_min, y_max, length.out = 4), digits = 2))
dev.off()

# trim trajectory to remove small branches
if(length(unique(pData(ica_HSMM_myo)$State)) > 3)
  ica_HSMM_myo <- trimTree(ica_HSMM_myo, num_paths = 2, min_branch_thrsld = 0.2) #cannot trim
plot_cell_trajectory(ica_HSMM_myo)

if(length(unique(pData(lle_HSMM_myo)$State)) > 3)
  lle_HSMM_myo <- trimTree(lle_HSMM_myo)
plot_cell_trajectory(lle_HSMM_myo)

if(length(unique(pData(dm_HSMM_myo)$State)) > 3)
  dm_HSMM_myo <- trimTree(dm_HSMM_myo)
# plot_cell_trajectory(dm_HSMM_myo)

if(length(unique(pData(destiny_HSMM_myo)$State)) > 3)
  destiny_HSMM_myo <- trimTree(destiny_HSMM_myo)
plot_cell_trajectory(destiny_HSMM_myo)

########################################################################################################################################################
# run l1-graph with the following parameters
# create SI15C
########################################################################################################################################################
maxiter <- 20
eps <- 1.0000e-05
gstruct <- 'span-tree'
gamma <- 0.1 #smoothness
sigma <- 0.01000 #
lambda <- 10 #L1 g
nn <- 5
verbose = T

X <- t(dm_res[, 1:2])
# D <- nrow(X); N <- ncol(X)
# Z <- X
C0 <- X
Nz <- ncol(C0)

G <- get_knn(C0, nn)

# choose appropriate lambda, gamma and sigma

########################################################################################################################################################
# run l1-graph with the span-tree
# create SI15a_3
########################################################################################################################################################

HSMM_seq_pg_res_span_tree <- principal_graph(t(dm_res[, 1:2]), C0, G$G, maxiter = maxiter,
                                             eps = eps, gstruct = 'span-tree',
                                             lambda = lambda, gamma = gamma,
                                             sigma = sigma, nn = 5, verbose = T)

print(qplot(HSMM_seq_pg_res_span_tree$C[1, ], HSMM_seq_pg_res_span_tree$C[2, ], color = I('black')) + geom_point(aes(dm_res[, 1], dm_res[, 2], color = as.character(pData(HSMM_myo)$Time))) + ggtitle('spanning-tree'))

PCA_HSMM_myo_l1_span_tree <- reduceDimension(HSMM_myo, reduction_method = 'SGL-tree', maxiter = maxiter, initial_method = PCA,
                                         eps = eps, lambda = lambda, gamma = 0.005,
                                         sigma = 0.1, nn = 5, verbose = T)

PCA_HSMM_myo_l1_span_tree <- orderCells(PCA_HSMM_myo_l1_span_tree)
pdf(paste(SI_fig_dir, 'pca_l1_span_tree.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(PCA_HSMM_myo_l1_span_tree, color_by = 'Time', cell_size = 0.2, show_branch_points = F) + nm_theme()
dev.off()

ICA_HSMM_myo_l1_span_tree <- reduceDimension(HSMM_myo, reduction_method = 'SGL-tree', maxiter = maxiter, initial_method = ICA,
                                         eps = eps, lambda = lambda, gamma = 0.005,
                                         sigma = 0.1, nn = 5, verbose = T)
ICA_HSMM_myo_l1_span_tree <- orderCells(ICA_HSMM_myo_l1_span_tree)
ICA_HSMM_myo_l1_span_tree <- orderCells(ICA_HSMM_myo_l1_span_tree, root_state = GM_state(ICA_HSMM_myo_l1_span_tree))
pdf(paste(SI_fig_dir, 'ica_l1_span_tree.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(ICA_HSMM_myo_l1_span_tree, color_by = 'Time', cell_size = 0.2, show_branch_points = F) + nm_theme()
dev.off()

LLE2_HSMM_myo_l1_span_tree <- reduceDimension(HSMM_myo, reduction_method = 'SGL-tree', maxiter = maxiter, initial_method = LLE2,
                                         eps = eps, lambda = lambda, gamma = 0.005,
                                         sigma = 0.1, nn = 5, verbose = T)
LLE2_HSMM_myo_l1_span_tree <- orderCells(LLE2_HSMM_myo_l1_span_tree)
pdf(paste(SI_fig_dir, 'lle_l1_span_tree.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(LLE2_HSMM_myo_l1_span_tree, color_by = 'Time', cell_size = 0.2, show_branch_points = F) + nm_theme()
dev.off()

DM_HSMM_myo_l1_span_tree <- reduceDimension(HSMM_myo, reduction_method = 'SGL-tree', maxiter = maxiter, initial_method = DM,
                                         eps = eps, lambda = lambda, gamma = gamma,
                                         sigma = sigma, nn = 5, verbose = T)
DM_HSMM_myo_l1_span_tree <- orderCells(DM_HSMM_myo_l1_span_tree)
pdf(paste(SI_fig_dir, 'dm_l1_span_tree.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(DM_HSMM_myo_l1_span_tree, color_by = 'Time', cell_size = 0.2, show_branch_points = F) + nm_theme()
dev.off()

destiny_HSMM_myo_l1_span_tree <- reduceDimension(HSMM_myo, reduction_method = 'SGL-tree', maxiter = maxiter, initial_method = destiny_diffusionMaps,
                                            eps = eps, lambda = lambda, gamma = gamma,
                                            sigma = sigma, nn = 5, verbose = T)
destiny_HSMM_myo_l1_span_tree <- orderCells(destiny_HSMM_myo_l1_span_tree)
pdf(paste(SI_fig_dir, 'destiny_l1_span_tree.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(destiny_HSMM_myo_l1_span_tree, color_by = 'Time', cell_size = 0.2, show_branch_points = F) + nm_theme()
dev.off()

########################################################################################################################################################
# run l1-graph with the graph
# create SI15a_4
########################################################################################################################################################
HSMM_seq_pg_res_span_tree <- principal_graph(t(dm_res[, 1:2]), C0, G$G, maxiter = maxiter,
                                             eps = eps, gstruct = 'span-tree',
                                             lambda = lambda, gamma = gamma,
                                             sigma = sigma, nn = 5, verbose = T)
DM_HSMM_myo_l1_graph <- reduceDimension(HSMM_myo, reduction_method = 'L1-graph',  initial_method = DM, #C0 = DM_HSMM_myo_l1_span_tree@reducedDimK, maxiter = 2,
                                        eps = eps, lambda = lambda, gamma = gamma,
                                        sigma = sigma, nn = 5, verbose = T)
DM_HSMM_myo_l1_graph <- orderCells(DM_HSMM_myo_l1_graph)

C0 <- PCA_HSMM_myo_l1_span_tree@reducedDimK
Nz <- ncol(C0)

G <- get_knn(C0, 15)
PCA_HSMM_myo_l1_graph <- reduceDimension(HSMM_myo, reduction_method = 'L1-graph',  initial_method = PCA, maxiter = 1,
                                             eps = eps, lambda = lambda, gamma = gamma, #G = G$G,
                                         sigma = sigma, nn = 15, #lambda = 0.01, gamma = 0.01, sigma = 100, nn = 5,
                                        verbose = T)


# PCA_HSMM_myo_l1_graph <- orderCells(PCA_HSMM_myo_l1_graph)
plot_cell_trajectory(PCA_HSMM_myo_l1_graph, color_by = 'Time', cell_size = 0.2, show_branch_points = F)
plot_cell_trajectory(PCA_HSMM_myo_l1_graph, color_by = 'Pseudotime', cell_size = 0.2, show_branch_points = F)
plot_cell_trajectory(PCA_HSMM_myo_l1_graph, color_by = 'State', cell_size = 0.2, show_branch_points = F) + facet_wrap(~State)
# PCA_HSMM_myo_l1_graph <- orderCells(PCA_HSMM_myo_l1_graph, root_state = GM_state(PCA_HSMM_myo_l1_graph))

pdf(paste(SI_fig_dir, 'pca_l1_graph.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(PCA_HSMM_myo_l1_graph, color_by = 'Time', cell_size = 0.2, show_branch_points = F) + nm_theme()
dev.off()

ICA_HSMM_myo_l1_graph <- reduceDimension(HSMM_myo, reduction_method = 'L1-graph', initial_method = ICA, # C0 = ICA_HSMM_myo_l1_span_tree@reducedDimK, maxiter = 2,
                                             eps = eps, lambda = lambda, gamma = gamma,
                                         sigma = sigma, #lambda = 0.01, gamma = 0.01, sigma = 100,
                                         nn = 5, verbose = T)
# ICA_HSMM_myo_l1_graph <- orderCells(ICA_HSMM_myo_l1_graph)
plot_cell_trajectory(ICA_HSMM_myo_l1_graph, color_by = 'Time', cell_size = 0.2, show_branch_points = F)
plot_cell_trajectory(ICA_HSMM_myo_l1_graph, color_by = 'Pseudotime', cell_size = 0.2, show_branch_points = F)
plot_cell_trajectory(ICA_HSMM_myo_l1_graph, color_by = 'State', cell_size = 0.2, show_branch_points = F) + facet_wrap(~State)
# ICA_HSMM_myo_l1_graph <- orderCells(ICA_HSMM_myo_l1_graph, root_state = GM_state(ICA_HSMM_myo_l1_graph))

pdf(paste(SI_fig_dir, 'ica_l1_graph.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(ICA_HSMM_myo_l1_graph, color_by = 'Time', cell_size = 0.2, show_branch_points = F) + nm_theme()
dev.off()

LLE2_HSMM_myo_l1_graph2 <- reduceDimension(HSMM_myo, reduction_method = 'L1-graph', initial_method = LLE2, #C0 = LLE2_HSMM_myo_l1_span_tree@reducedDimK, maxiter = 2,
                                              eps = eps, lambda = lambda, gamma = gamma,
                                           sigma = sigma, #lambda = lambda, gamma = 0.005, sigma = 0.1, nn = 5,
                                          verbose = T)
# LLE2_HSMM_myo_l1_graph <- orderCells(LLE2_HSMM_myo_l1_graph2)
LLE2_HSMM_myo_l1_graph <- LLE2_HSMM_myo_l1_graph2
plot_cell_trajectory(LLE2_HSMM_myo_l1_graph, color_by = 'Time', cell_size = 0.2, show_branch_points = F) 
plot_cell_trajectory(LLE2_HSMM_myo_l1_graph, color_by = 'Pseudotime', cell_size = 0.2, show_branch_points = F) 
plot_cell_trajectory(LLE2_HSMM_myo_l1_graph, color_by = 'State', cell_size = 0.2, show_branch_points = F) + facet_wrap(~State)
# LLE2_HSMM_myo_l1_graph <- orderCells(LLE2_HSMM_myo_l1_graph, root_state = GM_state(LLE2_HSMM_myo_l1_graph))

pdf(paste(SI_fig_dir, 'lle_l1_graph.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(LLE2_HSMM_myo_l1_graph, color_by = 'Time', cell_size = 0.2, show_branch_points = F) + nm_theme()
dev.off()

DM_HSMM_myo_l1_graph <- reduceDimension(HSMM_myo, reduction_method = 'L1-graph',  initial_method = DM, #C0 = DM_HSMM_myo_l1_span_tree@reducedDimK, maxiter = 2,
                                            eps = eps, lambda = lambda, gamma = gamma,
                                            sigma = sigma, nn = 5, verbose = T)
# DM_HSMM_myo_l1_graph <- orderCells(DM_HSMM_myo_l1_graph)
plot_cell_trajectory(DM_HSMM_myo_l1_graph, color_by = 'Time', cell_size = 0.2, show_branch_points = F) 
plot_cell_trajectory(DM_HSMM_myo_l1_graph, color_by = 'Pseudotime', cell_size = 0.2, show_branch_points = F) 
plot_cell_trajectory(DM_HSMM_myo_l1_graph, color_by = 'State', cell_size = 0.2, show_branch_points = F) + facet_wrap(~State)
# DM_HSMM_myo_l1_graph <- orderCells(DM_HSMM_myo_l1_graph, root_state = GM_state(LLE2_HSMM_myo_l1_graph))

pdf(paste(SI_fig_dir, 'dm_l1_graph.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(DM_HSMM_myo_l1_graph, color_by = 'Time', cell_size = 0.2, show_branch_points = F) + nm_theme()
dev.off()

destiny_HSMM_myo_l1_graph <- reduceDimension(HSMM_myo, reduction_method = 'L1-graph',  initial_method = destiny_diffusionMaps, #C0 = DM_HSMM_myo_l1_span_tree@reducedDimK, maxiter = 2,
                                        eps = eps, lambda = lambda, gamma = gamma,
                                        sigma = sigma, nn = 5, verbose = T)
# destiny_HSMM_myo_l1_graph <- orderCells(destiny_HSMM_myo_l1_graph)
plot_cell_trajectory(destiny_HSMM_myo_l1_graph, color_by = 'Time', cell_size = 0.2, show_branch_points = F) 
plot_cell_trajectory(destiny_HSMM_myo_l1_graph, color_by = 'Pseudotime', cell_size = 0.2, show_branch_points = F) 
plot_cell_trajectory(destiny_HSMM_myo_l1_graph, color_by = 'State', cell_size = 0.2, show_branch_points = F) + facet_wrap(~State)
# destiny_HSMM_myo_l1_graph <- orderCells(destiny_HSMM_myo_l1_graph, root_state = GM_state(LLE2_HSMM_myo_l1_graph))

pdf(paste(SI_fig_dir, 'dm_l1_graph.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(destiny_HSMM_myo_l1_graph, color_by = 'Time', cell_size = 0.2, show_branch_points = F) + nm_theme()
dev.off()

PCA_HSMM_myo_l1_graph <- DM_HSMM_myo_l1_graph
ICA_HSMM_myo_l1_graph <- DM_HSMM_myo_l1_graph
LLE2_HSMM_myo_l1_graph2 <- DM_HSMM_myo_l1_graph
LLE2_HSMM_myo_l1_graph <- LLE2_HSMM_myo_l1_graph2
destiny_HSMM_myo_l1_graph <- DM_HSMM_myo_l1_graph

########################################################################################################################################################
# add benchmark result using DDRTree as reference: 
########################################################################################################################################################

plot_cell_trajectory(ica_HSMM_myo, color_by = 'Pseudotime')
plot_cell_trajectory(ica_HSMM_myo, color_by = 'Time')
ica_HSMM_myo <- orderCells(ica_HSMM_myo)
plot_cell_trajectory(ica_HSMM_myo)
plot_cell_trajectory(ica_HSMM_myo) + facet_wrap(~State)
ica_HSMM_myo <- orderCells(ica_HSMM_myo, root_state = GM_state(ica_HSMM_myo))

HSMM_myo_reference <- DDRTree_pca_HSMM_myo

# reorder the cds to ensure cells start from the correct time point
HSMM_myo_reference <- orderCells(HSMM_myo_reference, root_state = GM_state(HSMM_myo_reference))
DDRTree_pca_HSMM_myo <- orderCells(DDRTree_pca_HSMM_myo, root_state = GM_state(DDRTree_pca_HSMM_myo))
DDRTree_ica_HSMM_myo <- orderCells(DDRTree_ica_HSMM_myo, root_state = GM_state(DDRTree_ica_HSMM_myo))
DDRTree_lle_HSMM_myo <- orderCells(DDRTree_lle_HSMM_myo, root_state = GM_state(DDRTree_lle_HSMM_myo))
DDRTree_dm_HSMM_myo <- orderCells(DDRTree_dm_HSMM_myo, root_state = GM_state(DDRTree_dm_HSMM_myo))
DDRTree_destiny_HSMM_myo <- orderCells(DDRTree_destiny_HSMM_myo, root_state = GM_state(DDRTree_destiny_HSMM_myo))

pca_HSMM_myo <- orderCells(pca_HSMM_myo, root_state = GM_state(DDRTree_pca_HSMM_myo))
ica_HSMM_myo <- orderCells(ica_HSMM_myo, root_state = GM_state(ica_HSMM_myo))
lle_HSMM_myo <- orderCells(lle_HSMM_myo, root_state = GM_state(lle_HSMM_myo))
dm_HSMM_myo <- orderCells(dm_HSMM_myo, root_state = GM_state(dm_HSMM_myo))
destiny_HSMM_myo <- orderCells(destiny_HSMM_myo, root_state = GM_state(destiny_HSMM_myo))

PCA_HSMM_myo_l1_span_tree <- orderCells(PCA_HSMM_myo_l1_span_tree, root_state = GM_state(PCA_HSMM_myo_l1_span_tree))
ICA_HSMM_myo_l1_span_tree <- orderCells(ICA_HSMM_myo_l1_span_tree, root_state = GM_state(ICA_HSMM_myo_l1_span_tree))
LLE2_HSMM_myo_l1_span_tree <- orderCells(LLE2_HSMM_myo_l1_span_tree, root_state = GM_state(LLE2_HSMM_myo_l1_span_tree))
DM_HSMM_myo_l1_span_tree <- orderCells(DM_HSMM_myo_l1_span_tree, root_state = GM_state(DM_HSMM_myo_l1_span_tree))
destiny_HSMM_myo_l1_span_tree <- orderCells(destiny_HSMM_myo_l1_span_tree, root_state = GM_state(destiny_HSMM_myo_l1_span_tree))

PCA_HSMM_myo_l1_graph <- orderCells(PCA_HSMM_myo_l1_graph, root_state = GM_state(PCA_HSMM_myo_l1_graph))
ICA_HSMM_myo_l1_graph <- orderCells(ICA_HSMM_myo_l1_graph, root_state = GM_state(ICA_HSMM_myo_l1_graph))
LLE2_HSMM_myo_l1_graph <- orderCells(LLE2_HSMM_myo_l1_graph, root_state = GM_state(LLE2_HSMM_myo_l1_graph))
DM_HSMM_myo_l1_graph <- orderCells(DM_HSMM_myo_l1_graph, root_state = GM_state(DM_HSMM_myo_l1_graph))
destiny_HSMM_myo_l1_graph <- orderCells(destiny_HSMM_myo_l1_graph, root_state = GM_state(destiny_HSMM_myo_l1_graph))

#DDRTree
DDRTree_pca_kendall_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(DDRTree_pca_HSMM_myo)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
DDRTree_pca_pearson_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(DDRTree_pca_HSMM_myo)$Pseudotime)
DDRTree_ica_kendall_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(DDRTree_ica_HSMM_myo)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
DDRTree_ica_pearson_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(DDRTree_ica_HSMM_myo)$Pseudotime)
DDRTree_lle_kendall_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(DDRTree_lle_HSMM_myo)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
DDRTree_lle_pearson_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(DDRTree_lle_HSMM_myo)$Pseudotime)
DDRTree_dm_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(DDRTree_dm_HSMM_myo)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
DDRTree_dm_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(DDRTree_dm_HSMM_myo)$Pseudotime)
DDRTree_destiny_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(DDRTree_destiny_HSMM_myo)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
DDRTree_destiny_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(DDRTree_destiny_HSMM_myo)$Pseudotime)

#simplePPT
simplePPT_pca_kendall_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(pca_HSMM_myo)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
simplePPT_pca_pearson_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(pca_HSMM_myo)$Pseudotime)
simplePPT_ica_kendall_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(ica_HSMM_myo)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
simplePPT_ica_pearson_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(ica_HSMM_myo)$Pseudotime)
simplePPT_lle_kendall_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(lle_HSMM_myo)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
simplePPT_lle_pearson_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(lle_HSMM_myo)$Pseudotime)
simplePPT_dm_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(dm_HSMM_myo)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
simplePPT_dm_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(dm_HSMM_myo)$Pseudotime)
simplePPT_destiny_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(destiny_HSMM_myo)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
simplePPT_destiny_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(destiny_HSMM_myo)$Pseudotime)

#l1 tree: 
l1_tree_pca_kendall_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(PCA_HSMM_myo_l1_span_tree)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
l1_tree_pca_pearson_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(PCA_HSMM_myo_l1_span_tree)$Pseudotime)
l1_tree_ica_kendall_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(ICA_HSMM_myo_l1_span_tree)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
l1_tree_ica_pearson_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(ICA_HSMM_myo_l1_span_tree)$Pseudotime)
l1_tree_lle_kendall_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(LLE2_HSMM_myo_l1_span_tree)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
l1_tree_lle_pearson_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(LLE2_HSMM_myo_l1_span_tree)$Pseudotime)
l1_tree_dm_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(DM_HSMM_myo_l1_span_tree)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
l1_tree_dm_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(DM_HSMM_myo_l1_span_tree)$Pseudotime)
l1_tree_destiny_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(destiny_HSMM_myo_l1_span_tree)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
l1_tree_destiny_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(destiny_HSMM_myo_l1_span_tree)$Pseudotime)

#l1 graph: 
l1_graph_pca_kendall_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(PCA_HSMM_myo_l1_graph)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
l1_graph_pca_pearson_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(PCA_HSMM_myo_l1_graph)$Pseudotime)
l1_graph_ica_kendall_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(ICA_HSMM_myo_l1_graph)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
l1_graph_ica_pearson_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(ICA_HSMM_myo_l1_graph)$Pseudotime)
l1_graph_lle_kendall_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(LLE2_HSMM_myo_l1_graph)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
l1_graph_lle_pearson_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(LLE2_HSMM_myo_l1_graph)$Pseudotime)
l1_graph_dm_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(DM_HSMM_myo_l1_graph)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
l1_graph_dm_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(DM_HSMM_myo_l1_graph)$Pseudotime)
l1_graph_destiny_tau <- cor(rank(pData(HSMM_myo_reference)$Pseudotime), rank(pData(destiny_HSMM_myo_l1_graph)$Pseudotime), method = 'kendall', use = 'pairwise.complete.obs')
l1_graph_destiny_rho <- cor(pData(HSMM_myo_reference)$Pseudotime, pData(destiny_HSMM_myo_l1_graph)$Pseudotime)

#create statistic data.frame: 
all_cor_df <- data.frame("Pearson_cor" = c(DDRTree_pca_pearson_rho, DDRTree_ica_pearson_rho, DDRTree_lle_pearson_rho, DDRTree_destiny_rho, #DDRTree_dm_rho, 
                                           simplePPT_pca_pearson_rho, simplePPT_ica_pearson_rho, simplePPT_lle_pearson_rho, simplePPT_destiny_rho, #simplePPT_dm_rho, 
                                           l1_tree_pca_pearson_rho, l1_tree_ica_pearson_rho, l1_tree_lle_pearson_rho, l1_tree_destiny_rho, #l1_tree_dm_rho, 
                                           l1_graph_pca_pearson_rho, l1_graph_ica_pearson_rho, l1_graph_lle_pearson_rho, l1_graph_dm_rho), #l1_graph_destiny_rho
                         "Kendall_cor" = c(DDRTree_pca_kendall_tau, DDRTree_ica_kendall_tau, DDRTree_lle_kendall_tau, DDRTree_destiny_tau, #DDRTree_dm_tau, , 
                                           simplePPT_pca_kendall_tau, simplePPT_ica_kendall_tau, simplePPT_lle_kendall_tau, simplePPT_destiny_tau, #simplePPT_dm_tau, 
                                           l1_tree_pca_kendall_tau, l1_tree_ica_kendall_tau, l1_tree_lle_kendall_tau, l1_tree_destiny_tau, #l1_tree_dm_tau, 
                                           l1_graph_pca_kendall_tau, l1_graph_ica_kendall_tau, l1_graph_lle_kendall_tau, l1_graph_dm_tau), #l1_graph_destiny_tau
                         dimension_reduction = rep(c("pca", 'ica', 'lle', 'dm'), 4), 
                         method = rep(c('DDRTree', 'simplePPT', 'L1 span-tree', 'L1 graph'), each = 4)
)
mlt_all_cor_df <- melt(all_cor_df, id.vars = c('dimension_reduction', 'method'))
mlt_all_cor_df$dimension_reduction <- as.character(mlt_all_cor_df$dimension_reduction)
mlt_all_cor_df$dimension_reduction <- factor(mlt_all_cor_df$dimension_reduction, levels = c("pca", 'ica', 'lle', 'dm'))

mlt_all_cor_df$method <- as.character(mlt_all_cor_df$method)
mlt_all_cor_df$method <- factor(mlt_all_cor_df$method, levels = c('DDRTree', 'simplePPT', 'L1 span-tree', 'L1 graph'))

####################################################################################################################################################################################
# create SI15b
####################################################################################################################################################################################

pdf(paste(SI_fig_dir, 'SI15b.pdf', sep = ''), height = 1.5, width = 4.5)
ggplot(aes(dimension_reduction, abs(value)), data = subset(mlt_all_cor_df, method != 'L1 graph')) + geom_bar(stat="identity", aes(fill = dimension_reduction)) + 
  facet_wrap(~variable + method, nrow = 2) + nm_theme() + ylim(0, 1) + nm_theme() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_size(range = c(0.2, 0.2)) + xlab('') + ylab('') + monocle_theme_opts() + 
  theme(strip.background = element_blank(), strip.text.x = element_blank())
dev.off()

####################################################################################################################################################################################

#adjusted rand index: 
#DDRTree
DDRTree_pca_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(DDRTree_pca_HSMM_myo)$State)
DDRTree_ica_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(DDRTree_ica_HSMM_myo)$State)
DDRTree_lle_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(DDRTree_lle_HSMM_myo)$State)
DDRTree_dm_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(DDRTree_dm_HSMM_myo)$State)
DDRTree_destiny_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(DDRTree_destiny_HSMM_myo)$State)

#simplePPT
simplePPT_pca_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(pca_HSMM_myo)$State)
simplePPT_ica_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(ica_HSMM_myo)$State)
simplePPT_lle_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(lle_HSMM_myo)$State)
simplePPT_dm_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(dm_HSMM_myo)$State)
simplePPT_destiny_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(destiny_HSMM_myo)$State)

#l1 tree: 
l1_tree_pca_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(PCA_HSMM_myo_l1_span_tree)$State)
l1_tree_ica_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(ICA_HSMM_myo_l1_span_tree)$State)
l1_tree_lle_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(LLE2_HSMM_myo_l1_span_tree)$State)
l1_tree_dm_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(DM_HSMM_myo_l1_span_tree)$State)
l1_tree_destiny_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(destiny_HSMM_myo_l1_span_tree)$State)

#l1 graph: 
l1_graph_pca_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(PCA_HSMM_myo_l1_graph)$State)
l1_graph_ica_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(ICA_HSMM_myo_l1_graph)$State)
l1_graph_lle_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(LLE2_HSMM_myo_l1_graph)$State)
l1_graph_dm_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(DM_HSMM_myo_l1_graph)$State)
l1_graph_destiny_adj_rand_ind <- calClusteringMetrics(pData(HSMM_myo_reference)$State, pData(destiny_HSMM_myo_l1_graph)$State)

adj_rand_ind_df <- data.frame("Adjusted Rand index" = c(DDRTree_pca_adj_rand_ind$randIndex[3], DDRTree_ica_adj_rand_ind$randIndex[3], DDRTree_lle_adj_rand_ind$randIndex[3], DDRTree_destiny_adj_rand_ind$randIndex[3], #DDRTree_dm_adj_rand_ind$randIndex[3], 
                                                        simplePPT_pca_adj_rand_ind$randIndex[3], simplePPT_ica_adj_rand_ind$randIndex[3], simplePPT_lle_adj_rand_ind$randIndex[3], simplePPT_destiny_adj_rand_ind$randIndex[3], #simplePPT_dm_adj_rand_ind$randIndex[3], 
                                                        l1_tree_pca_adj_rand_ind$randIndex[3], l1_tree_ica_adj_rand_ind$randIndex[3], l1_tree_lle_adj_rand_ind$randIndex[3], l1_tree_destiny_adj_rand_ind$randIndex[3], #l1_tree_dm_adj_rand_ind$randIndex[3], 
                                                        l1_graph_pca_adj_rand_ind$randIndex[3], l1_graph_ica_adj_rand_ind$randIndex[3], l1_graph_lle_adj_rand_ind$randIndex[3], l1_graph_dm_adj_rand_ind$randIndex[3]), #wishbone ,
                              # "Time_kcor" = c(Time_kcor1, Time_kcor2, Time_kcor3, Time_kcor4, Time_kcor5, Time_kcor6),
                              dimension_reduction = rep(c("pca", 'ica', 'lle', 'dm'), 4), 
                              method = rep(c('DDRTree', 'SimplePPT', 'SGL-tree', 'L1 graph'), each = 4))

mlt_adj_rand_ind_df <- melt(adj_rand_ind_df, id.vars = c('dimension_reduction', 'method'))
mlt_adj_rand_ind_df$dimension_reduction <- as.character(mlt_adj_rand_ind_df$dimension_reduction)
mlt_adj_rand_ind_df$dimension_reduction <- factor(mlt_adj_rand_ind_df$dimension_reduction, levels = c("pca", 'ica', 'lle', 'dm'))

mlt_adj_rand_ind_df$method <- as.character(mlt_adj_rand_ind_df$method)
mlt_adj_rand_ind_df$method <- factor(mlt_adj_rand_ind_df$method, levels = c('DDRTree', 'SimplePPT', 'SGL-tree', 'L1 graph'))

pdf(paste(SI_fig_dir, 'initial_method_hsmm_simplePPT.pdf', sep = ''), height = 0.75, width = 4.5)
ggplot(aes(dimension_reduction, abs(value)), data = subset(mlt_adj_rand_ind_df, method != 'L1 graph')) + geom_bar(stat="identity", aes(fill = dimension_reduction)) + 
  facet_wrap(~variable + method, nrow = 1) + nm_theme() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_size(range = c(0.2, 0.2)) + xlab('') + ylab('') + monocle_theme_opts() + 
  theme(strip.background = element_blank(), strip.text.x = element_blank())
dev.off()


####################################################################################################################################################################################
# create SI15c
####################################################################################################################################################################################
#combine both data.frame: 
all_mlt_df <- rbind(mlt_all_cor_df, mlt_adj_rand_ind_df)
pdf(paste(SI_fig_dir, 'initial_method_hsmm_all_method.pdf', sep = ''), height = 1.5, width = 4.5)
ggplot(aes(dimension_reduction, abs(value)), data = subset(all_mlt_df, method != 'L1 graph')) + geom_bar(stat="identity", aes(fill = dimension_reduction)) + 
  facet_wrap(~variable + method, nrow = 3) + nm_theme() + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_size(range = c(0.2, 0.2)) + xlab('') + ylab('') + monocle_theme_opts() + 
  theme(strip.background = element_blank(), strip.text.x = element_blank())
dev.off()

########################################################################################################################################################
#save the data: 
########################################################################################################################################################

save.image('./RData/fig_si2.RData')
