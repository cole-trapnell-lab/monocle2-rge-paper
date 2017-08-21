rm(list = ls())
# 1. The legend for panel A needs to use complete words for the various cell types.
# 2. Panel C should become panel B. Exclude the knockouts and the extra gates from this panel if they’re not already.
# 3. Panel E (or whichever one is the wildtype BEAM heatmap) should become panel C.
# 4. We need selected GO categories for the clusterers in the BEAM heatmap.
# 5. The fonts on the BEAM heatmap are unreadably small. Use the same legends and axis annotations as we have in the Census AI files. Make the heatmap pretty and clear :)
# 6. The Knockout panels should be labeled with standard notation (Gfi1-/-, Irf8-/-, and  Gfi1-/-/Irf8-/- I think)
# 7. Move the transient gate panels to the supplement.
# 8. Panel G is a bit confusing, because if I rememeber correctly, the Irf8 genes and Gfi1 genes are from the diff that’s *conditioned* on their pseudotimes. That’s not obvious from the figure and going to confuse the reader. What are we trying to communicate with this panel? I think it might be better to make a different UpSetR panel: one that shows the differentially between the Irf8-/-, Gfi-/-, and double KO cells and the WT cells collected at the corresponding gates. That is, presumably they used just one of their gates (maybe LKCD34+?) to collect the KO cells. So we should compare to the WT cells in that gate. Then we can show the intersection between those genes and the genes with ChIP peaks at their promoters to give a sense for how many direct targets are actually affected by loss of the regulator. Then we can show the intersection with the BEAM genes to show that BEAM is really picking up a big fraction of those genes, along with some others.
# 9. We need to make it clear what the overlap is between genes that have Gfi1 ChIP peaks and/or Irf8 ChIP peaks and BEAM genes. Maybe add this to panel G?
# 10. Split panel H in to two panels: one for Gfi1 and one for Irf8.
# 11. Panel J is not clear and right now doesn’t add anything. We should drop this or reformat so it says something useful.

get_correct_root_state <- function(cds){
  T0_counts <- table(pData(cds)$State, pData(cds)$Type)[,"Lsk"]
  as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))])
}


####################################################################################################################################################################################
#load all package
####################################################################################################################################################################################
library(RColorBrewer)
library(stringr)
library(UpSetR)
library(devtools)
library(rEDM)
library(xlsx)
library(xacHelper)
library(reshape2)
library(pheatmap)
library(d3Network)
library(netbiov)
library(monocle)
library(plyr)
library(dplyr)
library(venneuler)
library(SLICER)
library(destiny)
library(grid)
library(HSMMSingleCell)
library(piano)

source('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/scripts/function.R')

# fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/Figures/"
# tmp_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/tmp_figures/"
fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/Figures/"
tmp_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/Figures/tmp_figures/"

#figure 4 panels:
# main_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/Figures/main_figures/"
# SI_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/Figures/supplementary_figures/"
main_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/Figures/main_figures/"
SI_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/Figures/supplementary_figures/"

#########################################################################################################################################################################
# prepare the cds for both of the WT dataset and the full dataset
#########################################################################################################################################################################

#reading the exprs data and create a cell dataset:
# source('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/scripts/nature_hta_analysis_tmp.R')
hta_exprs <- read.csv("./Nature_hta_paper_Monocle_run_archive/input/Olsson_RSEM_SingleCellRNASeq.csv",row.names=1)
sample_sheet <- data.frame(groups = str_split_fixed(colnames(hta_exprs), "\\.+", 3), row.names = colnames(hta_exprs))
gene_ann <- data.frame(gene_short_name = row.names(hta_exprs), row.names = row.names(hta_exprs))
pd <- new("AnnotatedDataFrame",data=sample_sheet)
fd <- new("AnnotatedDataFrame",data=gene_ann)

tpm_mat <- hta_exprs
tpm_mat <- apply(tpm_mat, 2, function(x) x / sum(x) * 1e6)

URMM_all_std <- newCellDataSet(as.matrix(tpm_mat),phenoData = pd,featureData =fd,
                               expressionFamily = negbinomial.size(),
                               lowerDetectionLimit=1)
URMM_all_std <- estimateSizeFactors(URMM_all_std)
URMM_all_std <- estimateDispersions(URMM_all_std)

#set up the experimental type for each cell
pData(URMM_all_std)[, 'Type'] <- as.character(pData(URMM_all_std)[, 'groups.1']) #WT cells
pData(URMM_all_std)[453:593, 'Type'] <- paste(as.character(pData(URMM_all_std)[453:593, 'groups.1']), '_knockout', sep = '') #KO cells
pData(URMM_all_std)[594:640, 'Type'] <- paste(pData(URMM_all_std)[594:640, 'groups.1'], pData(URMM_all_std)[594:640, 'groups.2'], 'knockout', sep = '_') #double KO cells

#run Census to get the transcript counts
URMM_all_abs_list <- relative2abs(URMM_all_std, t_estimate = estimate_t(URMM_all_std), return_all = T)
URMM_all_abs <- newCellDataSet(as.matrix(URMM_all_abs_list$norm_cds),
                               phenoData = new("AnnotatedDataFrame",data=pData(URMM_all_std)),
                               featureData = new("AnnotatedDataFrame",data=fData(URMM_all_std)),
                               expressionFamily = negbinomial.size(),
                               lowerDetectionLimit=1)
URMM_all_abs <- estimateSizeFactors(URMM_all_abs)
URMM_all_abs <- estimateDispersions(URMM_all_abs)

URMM_all_abs <- setOrderingFilter(URMM_all_abs, row.names(fData(URMM_all_abs)))

#########################################################################################################################################################################
#read data from figure 1b
fig1b <- read.csv("./Nature_hta_paper_Monocle_run_archive/input/fig1b.txt",row.names=1, sep = '\t')

#match up the column name in fig1b to the colnames in URMM_all_fig1b
#note that you should not run this mutliple times
URMM_all_fig1b <- URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Lsk', 'Cmp', 'Gmp', 'LK')]

fig1b_names <- colnames(fig1b)
match_id <- which(str_split_fixed(colnames(fig1b), "\\.+", 2)[, 2] %in% colnames(URMM_all_fig1b) == T)
fig1b_names[match_id] <- str_split_fixed(colnames(fig1b), "\\.+", 2)[match_id, 2]
no_match_id <- which(str_split_fixed(colnames(fig1b), "\\.+", 2)[, 2] %in% colnames(URMM_all_fig1b) == F)
fig1b_names[no_match_id] <- str_split_fixed(colnames(fig1b), "\\.\\.", 2)[no_match_id, 2]
colnames(fig1b)[2:383] <- fig1b_names[2:383]

#learn trajectory using genes from figure 1b:
URMM_all_fig1b <- order_cell_tree(row.names(fig1b)[2:533], cds = URMM_all_fig1b,
                                  initial_method = PCA, norm_method = 'vstExprs', order_by = NULL)

#set up the color for each experiment
cols <- c("Lsk" = "#edf8fb", "Cmp" = "#ccece6", "Gmp" = "#99d8c9", "GG1" = "#66c2a4", "IG2" = "#41ae76", "Irf8" = "#238b45", "LK" = "#005824",
          "Irf8_knockout" = "#fc8d59", "Gfi1_Irf8_knockout" = "#636363", "Gfi1_knockout" = "#dd1c77")
plot_spanning_tree(URMM_all_fig1b, color_by = 'Type', cell_size = 2) + scale_color_manual(values = cols, name = "Type") #+ nm_theme()

#########################################################################################################################################################################
#assign clusters to each cell based on the clustering in the original study
pData(URMM_all_fig1b)$cluster <- 0
cluster_assignments <- as.numeric(fig1b[1, 2:383])
cluster_assignments <- revalue(as.factor(cluster_assignments), c("1" = "HSCP-1", "2" = "HSCP-2", "3" = "Meg", "4" = "Eryth",
                                                                 "5" = "Multi-Lin", "6" = "MDP", "7" = "Mono", "8" = "Gran", "9" = "Myelocyte"))

names(cluster_assignments) <- colnames(fig1b[1, 2:383])
pData(URMM_all_fig1b)$cluster <- cluster_assignments[row.names(pData(URMM_all_fig1b))]
#########################################################################################################################################################################
URMM_all_fig1b <- recreate_cds(URMM_all_fig1b)

URMM_all_fig1b <- setOrderingFilter(URMM_all_fig1b, ordering_genes = row.names(fig1b))
URMM_pc_variance <- plot_pc_variance_explained(URMM_all_fig1b, return_all = T, norm_method = 'log')

#2. run reduceDimension with tSNE as the reduction_method
set.seed(2017)
URMM_all_fig1b <- reduceDimension(URMM_all_fig1b, max_components=2, norm_method = 'log', reduction_method = 'tSNE', num_dim = 12,  verbose = T)

#3. initial run of clusterCells_Density_Peak
URMM_all_fig1b <- clusterCells_Density_Peak(URMM_all_fig1b, verbose = T, num_clusters = 5)

# plot_cell_clusters(URMM_all_fig1b, color_by = 'as.factor(cluster)', show_density_peak = T, show_density = T) # show_density = F,

#4. check the clusters
plot_cell_clusters(URMM_all_fig1b, color_by = 'as.factor(Cluster)', show_density_peak = T) # show_density = F,
plot_cell_clusters(URMM_all_fig1b, color_by = 'cluster')
plot_cell_clusters(URMM_all_fig1b, color_by = 'Cluster') + facet_wrap(~cluster) #, show_density = F
plot_cell_clusters(URMM_all_fig1b, color_by = 'Type') #, show_density = F
# # plot_cell_clusters(URMM_all_fig1b, color_by = 'State', show_density = F) + scale_color_brewer(palette = 'Set1') + facet_wrap(~paper_cluster)
# # plot_cell_clusters(URMM_all_fig1b, color_by = 'Cluster', show_density = F) + scale_color_brewer(palette = 'Set1') + facet_wrap(~paper_cluster)
#
# #5. also check the decision plot
plot_rho_delta(URMM_all_fig1b) #, rho_threshold = 6, delta_threshold = 7)
#
# #6. re-run cluster and skipping calculating the rho_sigma
#URMM_all_fig1b <- clusterCells_Density_Peak(URMM_all_fig1b, verbose = T,  rho_threshold = 6, delta_threshold = 7, skip_rho_sigma = T)
# # URMM_all_fig1b <- clusterCells_Density_Peak(URMM_all_fig1b, verbose = T,  num_clusters = 4)

#7. make the final clustering plot:
plot_cell_clusters(URMM_all_fig1b, color_by = 'as.factor(Cluster)') #, rho_threshold = 6, delta_threshold = 2, show_density = F
plot_cell_clusters(URMM_all_fig1b, color_by = 'as.factor(cluster)') #, rho_threshold = 6, delta_threshold = 2, show_density = F
plot_cell_clusters(URMM_all_fig1b, color_by = 'as.factor(Type)') #, rho_threshold = 6, delta_threshold = 2, show_density = F

#perform DEG test across clusters:
URMM_all_fig1b@expressionFamily <- negbinomial.size()
pData(URMM_all_fig1b)$Cluster <- factor(pData(URMM_all_fig1b)$Cluster)
URMM_clustering_DEG_genes <- differentialGeneTest(URMM_all_fig1b, fullModelFormulaStr = '~Cluster', cores = detectCores() - 2)

#use all DEG gene from the clusters
URMM_ordering_genes <- row.names(URMM_clustering_DEG_genes)[order(URMM_clustering_DEG_genes$qval)][1:1000]
#
# URMM_ordering_genes <- row.names(URMM_clustering_DEG_genes_subset)[order(URMM_clustering_DEG_genes_subset$qval)][1:1000]

URMM_all_fig1b <- setOrderingFilter(URMM_all_fig1b, ordering_genes = c(URMM_ordering_genes))
URMM_all_fig1b <- reduceDimension(URMM_all_fig1b, verbose = T, scaling = T, max_components = 4, maxIter = 100, norm_method = 'log',  lambda = 20 * ncol(URMM_all_fig1b)) #, maxIter = 100, initial_method = DM, R_tSNE, tSNE, destiny_diffusionMaps, maxIter = 100 , param.gamma = 100, ncenter = 100
URMM_all_fig1b <- orderCells(URMM_all_fig1b)
plot_cell_trajectory(URMM_all_fig1b, color_by = 'Type')
plot_cell_trajectory(URMM_all_fig1b, color_by = 'cluster') + facet_wrap(~cluster)
plot_cell_trajectory(URMM_all_fig1b, color_by = 'cluster', x = 1, y = 3) + facet_wrap(~cluster)

pData(URMM_all_fig1b)$paper_cluster <- pData(URMM_all_fig1b)[colnames(URMM_all_fig1b), 'cluster']

pdf(paste(SI_fig_dir, "URMM_cluster_DEG_qval_0.01.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(URMM_all_fig1b, color_by = 'paper_cluster', show_branch_points = F) + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

pdf(paste(SI_fig_dir, "URMM_cluster_DEG_qval_0.01_helper.pdf", sep = ''))
plot_cell_trajectory(URMM_all_fig1b, color_by = 'paper_cluster')
dev.off()

########################################################################################################################################################
#run on the full dataset
URMM_all_abs <- recreate_cds(URMM_all_abs)
pData(URMM_all_abs)$paper_cluster <- NA
pData(URMM_all_abs)[colnames(URMM_all_fig1b), 'paper_cluster'] <- as.character(pData(URMM_all_fig1b)[, 'paper_cluster'])

URMM_all_abs <- setOrderingFilter(URMM_all_abs, ordering_genes = URMM_ordering_genes)
URMM_all_abs <- reduceDimension(URMM_all_abs, verbose = T, scaling = T, maxIter = 100, norm_method = 'log', max_components = 4, param.gamma = 100, lambda = 15 * ncol(URMM_all_fig1b)) # norm_method = 'log',, param.gamma = 100
URMM_all_abs <- orderCells(URMM_all_abs)
plot_cell_trajectory(URMM_all_abs, color_by = 'Type')
plot_cell_trajectory(URMM_all_abs, color_by = 'Type', x = 1, y = 3) + facet_wrap(~paper_cluster)

URMM_all_abs <- orderCells(URMM_all_abs, root_state = get_correct_root_state(URMM_all_abs))

# if(as.numeric(unique(pData(URMM_all_abs)$State)) > 3)
#   URMM_all_abs <- trimTree(URMM_all_abs)
plot_cell_trajectory(URMM_all_abs, color_by = 'State') + facet_wrap(~Type)

pdf(paste(SI_fig_dir, "all_data_URMM_cluster_DEG_qval_0.01.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(URMM_all_abs, color_by = 'Type') + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

pdf(paste(SI_fig_dir, "all_data_URMM_cluster_DEG_qval_0.01_helper.pdf", sep = ''))
plot_cell_trajectory(URMM_all_abs, color_by = 'paper_cluster')
dev.off()

########################################################################################################################################################
#add dpt, slicer, wishbone analysis for this dataset
########################################################################################################################################################
# URMM_all_fig1b

URMM_all_fig1b <- setOrderingFilter(URMM_all_fig1b, URMM_ordering_genes)
slicer_res <- run_slicer(URMM_all_fig1b[URMM_ordering_genes, ], start = which.min(pData(URMM_all_fig1b)$Pseudotime))

pdf(paste(SI_fig_dir, 'URMM_all_fig1b_slicer_branch.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(slicer_res$traj_lle[, 1], slicer_res$traj_lle[, 2], color = as.character(slicer_res$order_df$branches)) + nm_theme() + xlab('LLE 1') + ylab('LLE 2')
dev.off()
pdf(paste(SI_fig_dir, 'URMM_all_fig1b_slicer_time.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(slicer_res$traj_lle[, 1], slicer_res$traj_lle[, 2], color = pData(URMM_all_fig1b)$paper_cluster) + nm_theme() + xlab('LLE 1') + ylab('LLE 2')
dev.off()

data <- t(log2(exprs(URMM_all_fig1b[URMM_ordering_genes, ]) + 1))
write.csv(file = paste('./csv_data/Wishbone_test_data/', 'URMM_all_abs_fig_SI2', ".txt", sep = ''), as.matrix(data), quote = F, row.names = T)
row.names(subset(pData(URMM_all_fig1b), Pseudotime == 0))

URMM_ordering_geness_wishbone_res_downsampling_empirical <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/MAR_seq_fractioin_wishbone_df_fig_4_downsampling_empirical.txt', header = T, sep = '\t')
# qplot(dm1, dm2, data = wishbone_res_downsampling_empirical)
qplot(tSNE1, tSNE2, data = URMM_ordering_geness_wishbone_res_downsampling_empirical, color = as.character(branch), size = 0.5)

benchmark_type <- "URMM_all_fig1b"
dpt_res <- run_new_dpt(URMM_all_fig1b[URMM_ordering_genes, ], normalize = F)

pdf(paste(SI_fig_dir, 'URMM_all_fig1b_dpt_time.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(dpt_res$dm$DC1, dpt_res$dm$DC2, color = pData(URMM_all_fig1b)$paper_cluster) + nm_theme() + xlab('DM 1') + ylab('DM 2')
dev.off()

pdf(paste(SI_fig_dir, 'URMM_all_fig1b_dpt_time.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(dpt_res$dm$DC1, dpt_res$dm$DC2, color = as.factor(dpt_res$branch[, 1]))
dev.off()

pdf(paste(SI_fig_dir, 'URMM_all_fig1b_dpt_branch.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(dpt_res$dm$DC1, dpt_res$dm$DC2, color = as.character(dpt_res$branch[, 1])) + nm_theme() + xlab('DM 1') + ylab('DM 2')
dev.off()

# URMM_all_abs

URMM_all_abs <- setOrderingFilter(URMM_all_abs, URMM_ordering_genes)
slicer_res <- run_slicer(URMM_all_abs[URMM_ordering_genes, ], start = which.min(pData(URMM_all_abs)$Pseudotime))

pdf(paste(SI_fig_dir, 'URMM_all_abs_slicer_branch.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(slicer_res$traj_lle[, 1], slicer_res$traj_lle[, 2], color = as.character(slicer_res$order_df$branches)) + nm_theme() + xlab('LLE 1') + ylab('LLE 2')
dev.off()
pdf(paste(SI_fig_dir, 'URMM_all_abs_slicer_time.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(slicer_res$traj_lle[, 1], slicer_res$traj_lle[, 2], color = pData(URMM_all_abs)$paper_cluster) + nm_theme() + xlab('LLE 1') + ylab('LLE 2')
dev.off()

data <- t(log2(exprs(URMM_all_abs[URMM_ordering_genes, ]) + 1))
write.csv(file = paste('./csv_data/Wishbone_test_data/', 'URMM_all_abs_fig_SI2', ".txt", sep = ''), as.matrix(data), quote = F, row.names = T)
row.names(subset(pData(URMM_all_abs), Pseudotime == 0))

URMM_all_abs_wishbone_res_downsampling_empirical <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/MAR_seq_fractioin_wishbone_df_fig_4_downsampling_empirical.txt', header = T, sep = '\t')
# qplot(dm1, dm2, data = wishbone_res_downsampling_empirical)
# qplot(tSNE1, tSNE2, data = wishbone_res_downsampling_empirical, color = as.character(branch), size = 0.5)

dpt_res <- run_new_dpt(URMM_all_abs[URMM_ordering_genes, ], normalize = F)

pdf(paste(SI_fig_dir, 'URMM_dpt_time.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(dpt_res$dm$DC1, dpt_res$dm$DC2, color = pData(URMM_all_abs)$paper_cluster) + nm_theme() + xlab('DM 1') + ylab('DM 2')
dev.off()

pdf(paste(SI_fig_dir, 'URMM_dpt_branch.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(dpt_res$dm$DC1, dpt_res$dm$DC2, color = as.character(dpt_res$branch[, 1])) + nm_theme() + xlab('DM 1') + ylab('DM 2')
dev.off()

########################################################################################################################################################
#panel 4a
########################################################################################################################################################

return_rotation_mat <- function(theta) {
  theta <- theta / 180 * pi
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
}

rotation_matrix <- matrix(rep(0, 4), nrow = 2)

#set the color:
type_vec <- unique(pData(URMM_all_abs)$Type)
type_cols <- brewer.pal(9, name = 'Set1')
type_cols[6] <- "#6A3D9A"
names(type_cols) <- type_vec

cluster_vec <- unique(pData(URMM_all_fig1b)$cluster)
cluster_cols <- type_cols
cluster_cols[10] <- "#0600FC"

names(cluster_cols) <- cluster_vec

pData(URMM_all_fig1b)$Type <- factor(pData(URMM_all_fig1b[, colnames(URMM_all_fig1b)])$Type, levels = c('Lsk', 'Cmp', 'Gmp', 'LK'))
pdf(paste(main_fig_dir, 'fig4a.pdf', sep = ''), width = 6, height = 2.75)
plot_cell_trajectory(URMM_all_fig1b[, ], color_by = 'cluster', show_branch_points = T, theta = 120, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~Type, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = cluster_cols, name = "cluster") +theme (legend.position="right", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) + theme (legend.position="top", legend.title=element_blank()) #theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pdf(paste(main_fig_dir, 'fig4a_helper.pdf', sep = ''), width = 6, height = 2.75)
plot_cell_trajectory(URMM_all_fig1b[, ], color_by = 'cluster', show_branch_points = T, theta = 120, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~Type, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = cluster_cols, name = "cluster") +theme (legend.position="right", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) + theme (legend.position="top", legend.title=element_blank()) #theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pdf(paste(main_fig_dir, 'fig4a.1.pdf', sep = ''), width = 6, height = 2.75)
plot_complex_cell_trajectory(URMM_all_fig1b[, ], color_by = 'cluster', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~Type, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = cluster_cols, name = "cluster") +theme (legend.position="right", legend.title=element_blank()) +
  theme (legend.position="top", legend.title=element_blank()) #theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pdf(paste(main_fig_dir, 'fig4a.1.1.pdf', sep = ''), width = 2.85, height = 1.5)
plot_complex_cell_trajectory(URMM_all_fig1b[, ], color_by = 'cluster', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~Type, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = cluster_cols, name = "cluster") +theme (legend.position="right", legend.title=element_blank()) +
  theme (legend.position="none", legend.title=element_blank()) + theme_void() + theme (legend.position="none", legend.title=element_blank()) #theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pdf(paste(main_fig_dir, 'fig4a_helper.pdf', sep = ''))
plot_cell_trajectory(URMM_all_fig1b[, ], color_by = 'cluster', show_branch_points = F) + facet_wrap(~Type, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = cluster_cols, name = "cluster") +theme (legend.position="right", legend.title=element_blank()) +
  stat_density2d(aes(color=cluster), alpha=I(0.25), size=I(0.25)) #theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pdf(paste(main_fig_dir, 'fig4a.1.pdf', sep = ''), height = 5, width = 4)
plot_cell_trajectory(URMM_all_fig1b[, ], color_by = 'cluster', show_branch_points = F, theta = 90, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~cluster, nrow = 3, ncol = 3) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = cluster_cols, name = "cluster") +theme (legend.position="right", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) + theme (legend.position="top", legend.title=element_blank()) #theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pdf(paste(main_fig_dir, 'fig4a.1_helper.pdf', sep = ''))
plot_cell_trajectory(URMM_all_fig1b[, ], color_by = 'cluster', show_branch_points = F, theta = 90, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~cluster, nrow = 3, ncol = 3) + scale_size(range = c(0.2, 0.2)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = cluster_cols, name = "cluster") +theme (legend.position="right", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) + theme (legend.position="top", legend.title=element_blank()) #theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

plot_cell_trajectory(URMM_all_fig1b, color_by = 'State')
plot_cell_trajectory(URMM_all_fig1b, color_by = 'Pseudotime')
URMM_all_fig1b <- orderCells(URMM_all_fig1b, root_state = get_correct_root_state(URMM_all_abs))
########################################################################################################################################################
#panel 4b
########################################################################################################################################################
important_genes <- c('Gfi1', 'Irf8', 'Klf4', "Zeb2", "Per3", "S100a8", "Cebpe")

# branchpoints <- URMM_all_fig1b@auxOrderingData$DDRTree$branch_points
# pr_graph_cell_proj_closest_vertex <- URMM_all_fig1b@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex
# pr_graph_cell_proj_closest_vertex[, 1] <- paste('Y_', pr_graph_cell_proj_closest_vertex[, 1], sep = '')
# branchpoints_cluster <- pData(URMM_all_fig1b)$cluster[which(pr_graph_cell_proj_closest_vertex == branchpoints)]
#
# fig1b_gm_branch_point <- which(branchpoints_cluster %in% c('Mono', 'Gran')) # 1
# fig1b_ery_branch_point <- which(!(branchpoints_cluster %in% c('Mono', 'Gran'))) # 2

fig1b_gm_branch_point <- 2
fig1b_ery_branch_point <- 1

# branchpoints <- URMM_all_abs@auxOrderingData$DDRTree$branch_points
# pr_graph_cell_proj_closest_vertex <- URMM_all_abs@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex
# pr_graph_cell_proj_closest_vertex[, 1] <- paste('Y_', pr_graph_cell_proj_closest_vertex[, 1], sep = '')
# branchpoints_cluster <- pData(URMM_all_abs)$paper_cluster[which(pr_graph_cell_proj_closest_vertex == branchpoints)]
#
# fig1b_gm_branch_point <- which(branchpoints_cluster %in% c('Mono', 'Gran'))
# fig1b_ery_branch_point <- which(!(branchpoints_cluster %in% c('Mono', 'Gran')))

#fig1c table
fig1c <- read.csv("./Nature_hta_paper_Monocle_run_archive/input/fig1c.txt",row.names=1, sep = '\t')
fig1c_lineage_genes <- row.names(fig1c)[which(row.names(fig1c) == 'Cebpa'):nrow(fig1c)]
pData(URMM_all_abs)$paper_cluster <- NA
pData(URMM_all_abs)[colnames(URMM_all_fig1b), 'paper_cluster'] <- as.character(pData(URMM_all_fig1b)[, 'paper_cluster'])

pData(URMM_all_abs)$paper_cluster <- factor(pData(URMM_all_abs)$paper_cluster, levels = c("HSCP-1", "HSCP-2", "Multi-Lin", "MDP", "Eryth", "Meg", "Mono", "Gran", "Myelocyte"))
pData(URMM_all_abs)$Type <- factor(pData(URMM_all_abs)$Type, levels = c("Lsk", 'Cmp','Gmp', 'LK', "GG1", 'IG2', 'Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout'))
pData(URMM_all_fig1b)$Type <- factor(pData(URMM_all_fig1b)$Type, levels = c('Lsk', 'Cmp', 'Gmp', 'LK'))

pdf(paste(main_fig_dir, 'hta_cluster_dispersion_genes_fig1b_data_branch_curves.pdf', sep = ''), height = 13, width = 5)
plot_genes_branched_pseudotime(URMM_all_fig1b[fig1c_lineage_genes, ], color_by = 'cluster', branch_point=fig1b_gm_branch_point,  cell_size = 0.2, nrow = 3, ncol = 5) + scale_colour_brewer(palette = "Set1") +
  scale_size(range = c(0.2, 0.2)) #+ nm_theme()
dev.off()

pdf(paste(main_fig_dir, 'fig4b.pdf', sep = ''), height = 1.5, width = 4)
plot_cell_trajectory(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout')], color_by = 'Type', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3, theta = 90) + facet_wrap(~Type, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = type_cols, name = "Type")  + theme (legend.position="none", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) #theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pdf(paste(main_fig_dir, 'fig4b_helper.pdf', sep = ''), height = 1.5, width = 4)
plot_cell_trajectory(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout')], color_by = 'Type', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3, theta = 90) + facet_wrap(~Type, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = type_cols, name = "Type")  + theme (legend.position="none", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) #theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pdf(paste(main_fig_dir, 'fig4b.1.pdf', sep = ''), width = 4, height = 1.5)
plot_complex_cell_trajectory(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout')], color_by = 'Type', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~Type, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = type_cols, name = "Type")  + theme (legend.position="none", legend.title=element_blank())
#theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

# fake the root cell as Irf8_knockout for just making the figure
pData(URMM_all_abs)[which(pData(URMM_all_abs)$Pseudotime == 0), 'Type'] <- 'Irf8_knockout'
pdf(paste(main_fig_dir, 'fig4b.1.1.pdf', sep = ''), width = 2.7, height = 1.1)
plot_complex_cell_trajectory(URMM_all_abs[, which(pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout'))], color_by = 'Type', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~Type, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = type_cols, name = "Type")  + theme (legend.position="none", legend.title=element_blank())
#theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pData(URMM_all_abs)[which(pData(URMM_all_abs)$Pseudotime == 0), 'Type'] <- 'Irf8_knockout'
pdf(paste(main_fig_dir, 'fig4b.1.1_cole.pdf', sep = ''), width = 2.7, height = 1.1)
plot_complex_cell_trajectory(URMM_all_abs[, which(pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gmp'))], color_by = 'Type', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~Type, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = type_cols, name = "Type")  + theme (legend.position="none", legend.title=element_blank())
#theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pData(URMM_all_abs)[which(pData(URMM_all_abs)$Pseudotime == 0), 'Type'] <- 'Lsk'

pdf(paste(main_fig_dir, 'fig4c.pdf', sep = ''), height = 1.5, width = 2.667)
plot_cell_trajectory(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('GG1', 'IG2')], color_by = 'Type', show_branch_points = F, cell_size = 0.5, cell_link_size = 0.3, theta = 90) + facet_wrap(~Type) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_color_manual(values = type_cols, name = "paper_cluster") + theme (legend.position="none", legend.title=element_blank()) +
  stat_density2d(color='black', alpha=I(0.25), size=I(0.25))
#theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pdf(paste(main_fig_dir, 'fig4d.pdf', sep = ''), height = 2, width = 3)
plot_genes_branched_pseudotime(URMM_all_fig1b[important_genes[1:6], ], color_by = 'Type', cell_size = 0.2, nrow = 2, branch_point = fig1b_gm_branch_point,  ncol = 3) + scale_colour_brewer(palette = "Set1")  +
  scale_size(range = c(0.2, 0.2)) + nm_theme()
dev.off()

pdf(paste(main_fig_dir, 'fig4d.1.pdf', sep = ''), height = 2, width = 3)
plot_genes_branched_pseudotime(URMM_all_fig1b[c('Gfi1', "Irf8", "Ly86", "Csf1r", "Cebpe", "S100a8"), ], color_by = 'Type', branch_point = fig1b_gm_branch_point, cell_size = 0.2, nrow = 2, ncol = 3) + scale_colour_brewer(palette = "Set1")  +
  scale_size(range = c(0.2, 0.2)) + nm_theme()
dev.off()
# pdf(paste(SI_fig_dir, 'fig3a.pdf', sep = ''), height = 2, width = 3)
# plot_genes_branched_pseudotime(URMM_all_fig1b[c('Gata2', 'Flt3', 'Cebpa', 'Spi1'), ], color_by = 'cluster',
#                                branch_point=1,  cell_size = 1, nrow = 2,  ncol = 2, panel_order = c('Gata2', 'Flt3','Cebpa', 'Spi1')) + scale_colour_brewer(palette = "Set1") +
#   scale_size(range = c(0.2, 0.2)) + scale_color_manual(values = cluster_cols, name = "cluster") + nm_theme()
# dev.off()

# pdf(paste(SI_fig_dir, 'fig3b.pdf', sep = ''), height = 2, width = 3)
# plot_genes_branched_pseudotime(URMM_all_fig1b[c('Irf8', 'Csf1r', 'Gfi1', 'Cebpe'), ], color_by = 'cluster',
#                                branch_point=1,  cell_size = 1, nrow = 2,  ncol = 2, panel_order = c('Irf8', 'Csf1r', 'Gfi1', 'Cebpe')) + scale_colour_brewer(palette = "Set1") +
#   scale_size(range = c(0.2, 0.2)) +  scale_color_manual(values = cluster_cols, name = "cluster") + nm_theme()
# dev.off()

###################################################################################################################################################################################################
#panel 4c: we may need to add more annotation again here
########################################################################################################################################################
#order the cells in the correct ordering:
plot_cell_trajectory(URMM_all_fig1b)
plot_cell_trajectory(URMM_all_fig1b, color_by = 'Type')
plot_cell_trajectory(URMM_all_fig1b, color_by = 'Pseudotime')
URMM_all_fig1b <- orderCells(URMM_all_fig1b, root_state = get_correct_root_state(URMM_all_fig1b))
plot_cell_trajectory(URMM_all_fig1b, color_by = 'Pseudotime')

plot_cell_trajectory(URMM_all_abs)
plot_cell_trajectory(URMM_all_abs, color_by = 'Type')
plot_cell_trajectory(URMM_all_abs, color_by = 'Pseudotime')
URMM_all_abs <- orderCells(URMM_all_abs, root_state = get_correct_root_state(URMM_all_abs))
plot_cell_trajectory(URMM_all_abs, color_by = 'Pseudotime')

#run BEAM:
fig1b_beam_genes <- BEAM(URMM_all_fig1b, branch_point = fig1b_gm_branch_point, verbose = T, cores = detectCores() - 2)
URMM_all_abs_beam_genes <- BEAM(URMM_all_abs, branch_point = fig1b_gm_branch_point, verbose = T, cores = detectCores() - 2)

fig1b_beam_genes_proj_dup <- BEAM(URMM_all_fig1b, branch_point = fig1b_gm_branch_point, verbose = T, progenitor_method = 'duplicate', cores = detectCores() - 2)
URMM_all_abs_beam_genes_proj_dup <- BEAM(URMM_all_abs, branch_point = fig1b_gm_branch_point, verbose = T, progenitor_method = 'duplicate', cores = detectCores() - 2)

#make the branch heatmaps for the wiltype dataset:
pdf(paste(main_fig_dir, 'fig4e.pdf', sep = ''))
plot_genes_branched_heatmap(URMM_all_fig1b[row.names(fig1b_beam_genes)[fig1b_beam_genes$qval < 1e-1], ],
                            num_clusters = 8, show_rownames = F, branch_point = fig1b_gm_branch_point, cores = detectCores() - 2,
                            branch_labels = c("Granulocyte", "Monocyte"))# + nm_theme()
dev.off()

#perform GO term enrichment:
ph <- plot_genes_branched_heatmap(URMM_all_fig1b[row.names(fig1b_beam_genes)[fig1b_beam_genes$qval < 1e-1], ],
                                  num_clusters = 8, show_rownames = F, branch_point = fig1b_gm_branch_point, cores = detectCores() - 2,
                                  branch_labels = c("Granulocyte", "Monocyte"), return_heatmap = T)# + nm_theme()

#run the enrichment on a cluster of genes from the heatmap:
URMM_fig1b_clustering <- data.frame(Cluster=factor(cutree(ph$ph_res$tree_row, 8)))
clusters <- as.numeric(URMM_fig1b_clustering$Cluster)
names(clusters) <- fData(URMM_all_fig1b[row.names(URMM_fig1b_clustering), ])$gene_short_name

root_directory <- "/Users/xqiu/Dropbox (Cole Trapnell's Lab)/Shared Data/"
mouse_go_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/GO/MOUSE_GO_bp_with_GO_iea_symbol.gmt", sep=""), encoding="latin1")
names(mouse_go_gsc$gsc) <- str_split_fixed(names(mouse_go_gsc$gsc), "%", 2)[,1]

URMM_fig1_gsa_results_dp_genes <- collect_gsa_hyper_results(URMM_all_fig1b[, ], mouse_go_gsc, clusters)

pdf(paste(main_fig_dir, "URMM_fig1_gsa_results_dp_genes_go_enrichment.pdf", sep = ''), height=100, width=15)
plot_gsa_hyper_heatmap(URMM_all_fig1b, URMM_fig1_gsa_results_dp_genes, significance = 1e-1)
dev.off()

# pdf(paste(SI_fig_dir, 'hta_type_dispersion_genes_fig1b_data_branch_heatmap_fig1c.pdf', sep = ''))
# plot_genes_branched_heatmap(URMM_all_fig1b[row.names(fig1c), ],
#                             num_clusters = 4, show_rownames = T, branch_point = 2, cores = detectCores() - 2,
#                             branch_labels = c("Monocyte", "Granulocyte"))# + nm_theme()
# dev.off()

# #make the branch heatmaps for the full dataset:
# pdf(paste(SI_fig_dir, 'hta_type_dispersion_genes_all_data_branch_heatmap_all_significant.pdf', sep = ''))
# plot_genes_branched_heatmap(URMM_all_abs[row.names(URMM_all_abs_beam_genes)[URMM_all_abs_beam_genes$qval < 1e-1], ],
#                             num_clusters = 8, show_rownames = F, branch_point = 1, cores = detectCores() - 2,
#                             branch_labels = c("Monocyte", "Granulocyte"))# + nm_theme()
# dev.off()

# pdf(paste(SI_fig_dir, 'hta_type_dispersion_genes_fig1b_data_branch_heatmap_fig1c.pdf', sep = ''))
# plot_genes_branched_heatmap(URMM_all_abs[row.names(fig1c), ],
#                             num_clusters = 4, show_rownames = T, branch_point = 1, cores = detectCores() - 2,
#                             branch_labels = c("Monocyte", "Granulocyte"))# + nm_theme()
# dev.off()

###################################################################################################################################################################################################
# #save the data:
# URMM_all_fig1b_wt_exprs <- exprs(URMM_all_fig1b)[, order(pData(URMM_all_fig1b)$Pseudotime)]
# write.table(file = 'blood_wt_data.txt', as.matrix(URMM_all_fig1b_wt_exprs), col.names = T, row.names = T, quote = F)
#
# pData_wt <- pData(URMM_all_fig1b[, order(pData(URMM_all_fig1b)$Pseudotime)])
# fData_wt <- fData(URMM_all_fig1b[, order(pData(URMM_all_fig1b)$Pseudotime)])
# write.table(file = 'blood_wt_data_phenotype_data.txt', pData_wt, col.names = T, row.names = T, quote = F)
# write.table(file = 'blood_wt_data_genotype_data.txt', fData_wt, col.names = T, row.names = T, quote = F)
#
# #knockout only data:
# URMM_all_fig1b_knockout_exprs <- exprs(URMM_all_abs)[, pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout')]
# write.table(file = 'blood_knockout_data.txt', as.matrix(URMM_all_fig1b_knockout_exprs), col.names = T, row.names = T, quote = F)
#
# pData_knockout <- pData(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout')])
# fData_knockout <- fData(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout')])
# write.table(file = 'blood_knockout_data_phenotype_data.txt', pData_knockout, col.names = T, row.names = T, quote = F)
# write.table(file = 'blood_knockout_data_genotype_data.txt', fData_knockout, col.names = T, row.names = T, quote = F)
#
# #all data:
# URMM_all_abs_exprs <- exprs(URMM_all_abs)[, order(pData(URMM_all_abs)$Pseudotime)]
# write.table(file = 'blood_all_data.txt', as.matrix(URMM_all_abs_exprs), col.names = T, row.names = T, quote = F)
#
# pData_all <- pData(URMM_all_abs[, order(pData(URMM_all_abs)$Pseudotime)])
# fData_all <- fData(URMM_all_abs)
# write.table(file = 'blood_all_data_phenotype_data.txt', pData_all, col.names = T, row.names = T, quote = F)
# write.table(file = 'blood_all_data_genotype_data.txt', fData_all, col.names = T, row.names = T, quote = F)
#
# # > dim(URMM_all_abs_all_dispersion[row.names(fig1b_beam_genes)[fig1b_beam_genes$qval < 1e-1], ])
# # Features  Samples
# # 1101      640

###################################################################################################################################################################################################
#panel 4d: show only the KO cells on the tree
#######################################################################################################################################################
pData(URMM_all_abs)$paper_cluster <- factor(pData(URMM_all_abs)$paper_cluster, levels = c("HSCP-1", "HSCP-2", "Multi-Lin", "MDP", "Eryth", "Meg", "Mono", "Gran", "Myelocyte"))
pData(URMM_all_abs)$Type <- factor(pData(URMM_all_abs)$Type, levels = c("Lsk", 'Cmp','Gmp', 'LK', "GG1", 'IG2', 'Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout'))
pData(URMM_all_fig1b)$Type <- factor(pData(URMM_all_fig1b)$Type, levels = c('Lsk', 'Cmp', 'Gmp', 'LK'))

# plot_cell_trajectory(URMM_all_abs, color_by = 'paper_cluster', show_branch_points = F, theta = 90)

# pdf(paste(main_fig_dir, 'hta_type_dispersion_genes_all_data_branch_curves.pdf', sep = ''), height = 3, width = 1.5)
# plot_genes_branched_pseudotime(URMM_all_abs[important_genes, ], color_by = 'Type', cell_size = 1) + scale_colour_brewer(palette = "Set1")  +
#   scale_size(range = c(0.2, 0.2)) + nm_theme()
# dev.off()

# pdf(paste(main_fig_dir, 'all_data_same_wt_genes_tree.pdf', sep = ''), height = 1.5, width = 5)
# plot_cell_trajectory(URMM_all_abs, color_by = 'Type', show_branch_points = F, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~Type, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
#   nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = type_cols, name = "Type")  +theme (legend.position="right", legend.title=element_blank()) +
#   stat_density2d(aes(color=Type), alpha=I(0.25), size=I(0.25)) #theme(axis.text.x = element_text(angle = 30, hjust = 1))
# dev.off()

pdf(paste(main_fig_dir, 'fig4b.pdf', sep = ''), height = 1.5, width = 4)
plot_cell_trajectory(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout')], color_by = 'Type', show_branch_points = F, cell_size = 0.5, cell_link_size = 0.3, theta = 90) + facet_wrap(~Type, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = type_cols, name = "Type")  + theme (legend.position="none", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) #theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pdf(paste(SI_fig_dir, 'fig_si4a.pdf', sep = ''), height = 1.5, width = 2.667)
plot_cell_trajectory(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('GG1', 'IG2')], color_by = 'Type', show_branch_points = F, cell_size = 0.5, cell_link_size = 0.3, theta = 90) + facet_wrap(~Type) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_color_manual(values = type_cols, name = "paper_cluster") + theme (legend.position="none", legend.title=element_blank()) +
  stat_density2d(color='black', alpha=I(0.25), size=I(0.25))
#theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

# pdf(paste(main_fig_dir, 'hta_cluster_dispersion_genes_all_data.pdf', sep = ''), height = 2.5, width = 8)
# plot_cell_trajectory(URMM_all_abs, color_by = 'paper_cluster', show_branch_points = F, theta = 90) + facet_wrap(~Type, nrow = 2) + scale_size(range = c(0.2, 0.2)) +
#   nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
#   scale_color_manual(values = type_cols, name = "paper_cluster") +theme (legend.position="right", legend.title=element_blank()) +
#   stat_density2d(aes(color=paper_cluster), alpha=I(0.25), size=I(0.25))
# #theme(axis.text.x = element_text(angle = 30, hjust = 1))
# dev.off()

# pdf(paste(main_fig_dir, 'hta_type_dispersion_genes_all_data.pdf', sep = ''), height = 3, width = 5)
# plot_cell_trajectory(URMM_all_abs, color_by = 'Type', show_branch_points = F, theta = 90) + facet_wrap(~Type, nrow = 3) + scale_size(range = c(0.2, 0.2)) +
#   nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
#   scale_color_manual(values = type_cols, name = "Type") +theme (legend.position="right", legend.title=element_blank()) +
#   stat_density2d(aes(color=Type), alpha=I(0.25), size=I(0.25)) #theme(axis.text.x = element_text(angle = 30, hjust = 1))
# dev.off()

# plot_cell_trajectory(URMM_all_abs, color_by = 'State')
# plot_cell_trajectory(URMM_all_abs, color_by = 'cluster')
# plot_cell_trajectory(URMM_all_abs, color_by = 'Pseudotime')
# #URMM_all_abs <- orderCells(URMM_all_abs, root_state = 3)

# pdf(paste(main_fig_dir, 'hta_cluster_dispersion_genes_all_data_branch_curves.pdf', sep = ''), height = 2, width = 3)
# plot_genes_branched_pseudotime(URMM_all_abs[important_genes[1:6], ], color_by = 'paper_cluster', cell_size = 0.2, branch_point = 1, nrow = 2, ncol = 3) + scale_colour_brewer(palette = "Set1") +
#   scale_size(range = c(0.2, 0.2)) + nm_theme()
# dev.off()

# pdf(paste(main_fig_dir, 'hta_type_dispersion_genes_all_data_branch_curves.pdf', sep = ''), height = 3, width = 1.5)
# plot_genes_branched_pseudotime(URMM_all_abs[important_genes, ], color_by = 'Type', cell_size = 0.2, branch_point = 1) + scale_colour_brewer(palette = "Set1")  +
#   scale_size(range = c(0.2, 0.2)) + nm_theme()
# dev.off()

# #
# pdf(paste(main_fig_dir, 'hta_type_dispersion_genes_all_data_branch_heatmap.pdf', sep = ''))
# plot_genes_branched_heatmap(URMM_all_abs[row.names(fig1c), ], num_clusters = 4, show_rownames = T, branch_point = 1)# + nm_theme()
# dev.off()
###################################################################################################################################################################################################
#how to test the IG1, GG2 as well as the three knockouts:
#test on extra genotype column:
plot_cell_trajectory(URMM_all_abs) #branch 2 is the major branch
pData(URMM_all_abs)$Genotype <- 'WildType'
pData(URMM_all_abs)$Genotype[pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout')] <- as.character(pData(URMM_all_abs)$Type[pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout')])
URMM_all_abs_genotype_beam_genes <- BEAM(URMM_all_abs, branch_point = fig1b_gm_branch_point, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch + Genotype", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch")

pData(URMM_all_abs)$Genotype2 <- 'WildType'
pData(URMM_all_abs)$Genotype2[pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout', 'GG1', 'IG2')] <- as.character(pData(URMM_all_abs)$Type[pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout', 'GG1', 'IG2')])
URMM_all_abs_genotype2_beam_genes <- BEAM(URMM_all_abs, branch_point = fig1b_gm_branch_point, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch + Genotype2", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch")

pData(URMM_all_abs)$Genotype3 <- 'WildType'
pData(URMM_all_abs)$Genotype3[pData(URMM_all_abs)$Type %in% c('Gfi1_Irf8_knockout', 'GG1', 'IG2')] <- as.character(pData(URMM_all_abs)$Type[pData(URMM_all_abs)$Type %in% c('Gfi1_Irf8_knockout', 'GG1', 'IG2')])
URMM_all_abs_genotype3_beam_genes <- BEAM(URMM_all_abs, branch_point = fig1b_gm_branch_point, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch + Genotype3", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch")

pData(URMM_all_abs)$Genotype3.1 <- 'WildType'
pData(URMM_all_abs)$Genotype3.1[pData(URMM_all_abs)$Type %in% c('GG1', 'IG2')] <- as.character(pData(URMM_all_abs)$Type[pData(URMM_all_abs)$Type %in% c('GG1', 'IG2')])
URMM_all_abs_genotype3.1_beam_genes <- BEAM(URMM_all_abs, branch_point = fig1b_gm_branch_point, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch + Genotype3.1", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch")

pData(URMM_all_abs)$Genotype4 <- 'WildType'
pData(URMM_all_abs)$Genotype4[pData(URMM_all_abs)$Type %in% c('Irf8_knockout')] <- as.character(pData(URMM_all_abs)$Type[pData(URMM_all_abs)$Type %in% c('Irf8_knockout')])
URMM_all_abs_genotype4_beam_genes <- BEAM(URMM_all_abs, branch_point = fig1b_gm_branch_point, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch + Genotype4", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch")

###################################################################################################################################################################################################
#directly test on the mainfold using straight DEG test instead of BEAM
#identify knockout cells:
# range(pData(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Irf8_knockout')])$Pseudotime)
# [1] 19.03098 26.73462
rng <- range(pData(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Irf8_knockout')])$Pseudotime)

Irf8_ko_mainfold_cells <- URMM_all_abs[, pData(URMM_all_abs)$Pseudotime > floor(rng[1]) & pData(URMM_all_abs)$Pseudotime < ceiling(rng[2]) & !(pData(URMM_all_abs)$Type %in% c('Gfi1_knockout', 'Gfi1_Irf8_knockout'))]
pData(Irf8_ko_mainfold_cells)$Genotype4 <- 'WildType'
pData(Irf8_ko_mainfold_cells)$Genotype4[pData(Irf8_ko_mainfold_cells)$Type %in% c('Irf8_knockout')] <- as.character(pData(Irf8_ko_mainfold_cells)$Type[pData(Irf8_ko_mainfold_cells)$Type %in% c('Irf8_knockout')])
URMM_all_abs_genotype4.1_beam_genes <- differentialGeneTest(Irf8_ko_mainfold_cells, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*State + Genotype4", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*State")
###################################################################################################################################################################################################

pData(URMM_all_abs)$Genotype5 <- 'WildType'
pData(URMM_all_abs)$Genotype5[pData(URMM_all_abs)$Type %in% c('Gfi1_knockout')] <- as.character(pData(URMM_all_abs)$Type[pData(URMM_all_abs)$Type %in% c('Gfi1_knockout')])
URMM_all_abs_genotype5_beam_genes <- BEAM(URMM_all_abs, branch_point = fig1b_gm_branch_point, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch + Genotype5", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch")

###################################################################################################################################################################################################
#directly test on the mainfold using straight DEG test instead of BEAM
#identify knockout cells:
# range(pData(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Gfi1_knockout')])$Pseudotime)
# [1]  7.763848 23.435222
rng <- range(pData(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Gfi1_knockout')])$Pseudotime)

Gfi1_ko_mainfold_cells <- URMM_all_abs[, pData(URMM_all_abs)$Pseudotime > floor(rng[1]) & pData(URMM_all_abs)$Pseudotime < ceiling(rng[2]) & !(pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_Irf8_knockout'))]
pData(Gfi1_ko_mainfold_cells)$Genotype5 <- 'WildType'
pData(Gfi1_ko_mainfold_cells)$Genotype5[pData(Gfi1_ko_mainfold_cells)$Type %in% c('Gfi1_knockout')] <- as.character(pData(Gfi1_ko_mainfold_cells)$Type[pData(Gfi1_ko_mainfold_cells)$Type %in% c('Gfi1_knockout')])
URMM_all_abs_genotype5.1_beam_genes <- differentialGeneTest(Gfi1_ko_mainfold_cells, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*State + Genotype5", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*State")
###################################################################################################################################################################################################

pData(URMM_all_abs)$Genotype6 <- 'WildType'
pData(URMM_all_abs)$Genotype6[pData(URMM_all_abs)$Type %in% c('Gfi1_Irf8_knockout')] <- as.character(pData(URMM_all_abs)$Type[pData(URMM_all_abs)$Type %in% c('Gfi1_Irf8_knockout')])
URMM_all_abs_genotype6_beam_genes <- BEAM(URMM_all_abs, branch_point = fig1b_gm_branch_point, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch + Genotype6", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch")

###################################################################################################################################################################################################
#directly test on the mainfold using straight DEG test instead of BEAM
#identify knockout cells:
# range(pData(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Gfi1_Irf8_knockout')])$Pseudotime)
# [1] 18.23346 25.85384
rng <- range(pData(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Gfi1_Irf8_knockout')])$Pseudotime)

Gfi1_Irf8_ko_mainfold_cells <- URMM_all_abs[, pData(URMM_all_abs)$Pseudotime > floor(rng[1]) & pData(URMM_all_abs)$Pseudotime < ceiling(rng[2]) & !(pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout'))]
pData(Gfi1_Irf8_ko_mainfold_cells)$Genotype6 <- 'WildType'
pData(Gfi1_Irf8_ko_mainfold_cells)$Genotype6[pData(Gfi1_Irf8_ko_mainfold_cells)$Type %in% c('Gfi1_Irf8_knockout')] <- as.character(pData(Gfi1_Irf8_ko_mainfold_cells)$Type[pData(Gfi1_Irf8_ko_mainfold_cells)$Type %in% c('Gfi1_Irf8_knockout')])
URMM_all_abs_genotype6.1_beam_genes <- differentialGeneTest(Gfi1_Irf8_ko_mainfold_cells, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*State + Genotype6", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*State")
###################################################################################################################################################################################################

pData(URMM_all_abs)$Genotype7 <- 'WildType'
pData(URMM_all_abs)$Genotype7[pData(URMM_all_abs)$Type %in% c('GG1')] <- as.character(pData(URMM_all_abs)$Type[pData(URMM_all_abs)$Type %in% c('GG1')])
URMM_all_abs_genotype7_beam_genes <- BEAM(URMM_all_abs, branch_point = fig1b_gm_branch_point, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch + Genotype7", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch")

###################################################################################################################################################################################################
#directly test on the mainfold using straight DEG test instead of BEAM
#identify knockout cells:
# range(pData(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('GG1')])$Pseudotime)
# [1] 0.9855298 20.7861759
rng <- range(pData(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('GG1')])$Pseudotime)

GG1_ko_mainfold_cells <- URMM_all_abs[, pData(URMM_all_abs)$Pseudotime > floor(rng[1]) & pData(URMM_all_abs)$Pseudotime < ceiling(rng[2]) & !(pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout'))]
pData(GG1_ko_mainfold_cells)$Genotype7 <- 'WildType'
pData(GG1_ko_mainfold_cells)$Genotype7[pData(GG1_ko_mainfold_cells)$Type %in% c('GG1')] <- as.character(pData(GG1_ko_mainfold_cells)$Type[pData(GG1_ko_mainfold_cells)$Type %in% c('GG1')])
URMM_all_abs_genotype7.1_beam_genes <- differentialGeneTest(GG1_ko_mainfold_cells, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*State + Genotype7", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*State")
###################################################################################################################################################################################################

pData(URMM_all_abs)$Genotype8 <- 'WildType'
pData(URMM_all_abs)$Genotype8[pData(URMM_all_abs)$Type %in% c('IG2')] <- as.character(pData(URMM_all_abs)$Type[pData(URMM_all_abs)$Type %in% c('IG2')])
URMM_all_abs_genotype8_beam_genes <- BEAM(URMM_all_abs, branch_point = fig1b_gm_branch_point, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch + Genotype8", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch")

###################################################################################################################################################################################################
#directly test on the mainfold using straight DEG test instead of BEAM
#identify knockout cells:
range(pData(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('IG2')])$Pseudotime)
# [1]  6.964654 24.455794
rng <- range(pData(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('IG2')])$Pseudotime)

IG2_ko_mainfold_cells <- URMM_all_abs[, pData(URMM_all_abs)$Pseudotime > floor(rng[1]) & pData(URMM_all_abs)$Pseudotime < ceiling(rng[2])  & !(pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout'))]
pData(IG2_ko_mainfold_cells)$Genotype8 <- 'WildType'
pData(IG2_ko_mainfold_cells)$Genotype8[pData(IG2_ko_mainfold_cells)$Type %in% c('IG2')] <- as.character(pData(IG2_ko_mainfold_cells)$Type[pData(IG2_ko_mainfold_cells)$Type %in% c('IG2')])
URMM_all_abs_genotype8.1_beam_genes <- differentialGeneTest(IG2_ko_mainfold_cells, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*State + Genotype8", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*State")
###################################################################################################################################################################################################

# #do each test separetely: (cannot subset cds when ncenter is used)
# URMM_all_abs_subset <- SubSet_cds(URMM_all_abs, cells = colnames(URMM_all_abs)[pData(URMM_all_abs)$Genotype2 == 'Irf8_knockout'])
# URMM_all_abs_Irf8_knockout_beam_genes <- BEAM(, branch_point = 1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch")
# URMM_all_abs_Gfi1_knockout_beam_genes <- BEAM(URMM_all_abs[, pData(URMM_all_abs)$Genotype2 == 'Gfi1_knockout'], branch_point = 1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch")
# URMM_all_abs_Gfi1_Irf8_knockout_beam_genes <- BEAM(URMM_all_abs[, pData(URMM_all_abs)$Genotype2 == 'Gfi1_Irf8_knockout'], branch_point = 1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch")
# URMM_all_abs_GG1_beam_genes <- BEAM(URMM_all_abs[, pData(URMM_all_abs)$Genotype2 == 'GG1'], branch_point = 1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch")
# URMM_all_abs_IG2_beam_genes <- BEAM(URMM_all_abs[, pData(URMM_all_abs)$Genotype2 == 'IG2'], branch_point = 1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch")
#
# URMM_all_abs_Irf8_knockout_beam_genes <- BEAM(URMM_all_abs[, pData(URMM_all_abs)$Genotype2 %in% c('WildType', "Irf8_knockout")], branch_point = 1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch + Genotype2")
# URMM_all_abs_Gfi1_knockout_beam_genes <- BEAM(URMM_all_abs[, pData(URMM_all_abs)$Genotype2 == c('WildType', 'Gfi1_knockout')], branch_point = 1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch + Genotype2")
# URMM_all_abs_Gfi1_Irf8_knockout_beam_genes <- BEAM(URMM_all_abs[, pData(URMM_all_abs)$Genotype2 == c('WildType', 'Gfi1_Irf8_knockout')], branch_point = 1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch + Genotype2")
# URMM_all_abs_GG1_beam_genes <- BEAM(URMM_all_abs[, pData(URMM_all_abs)$Genotype2 == c('WildType', 'GG1')], branch_point = 1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch + Genotype2")
# URMM_all_abs_IG2_beam_genes <- BEAM(URMM_all_abs[, pData(URMM_all_abs)$Genotype2 == c('WildType', 'IG2')], branch_point = 1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch + Genotype2")

#test on just the three knockouts:
knockout_data <- URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout')]
knockout_genes <- differentialGeneTest(knockout_data, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~Type")

gate_cell_data <- URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('GG1', 'IG2')]
gate_cell_genes <- differentialGeneTest(gate_cell_data, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~Type")

double_knockout_gate_cell_data <- URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('GG1', 'IG2', 'Gfi1_Irf8_knockout')]
double_knockout_gate_cell_data_genes <- differentialGeneTest(double_knockout_gate_cell_data, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~Type")

#create a upset plot?
listInput <- list('knockout' = row.names(knockout_genes)[knockout_genes$qval <0.1],
                  'URMM_all_abs_genotype_beam_genes' = row.names(URMM_all_abs_genotype_beam_genes)[URMM_all_abs_genotype_beam_genes$qval <0.1],
                  'URMM_all_abs_genotype2_beam_genes' = row.names(URMM_all_abs_genotype2_beam_genes)[URMM_all_abs_genotype2_beam_genes$qval <0.1],
                  'gate_cell_data' =  row.names(gate_cell_genes)[gate_cell_genes$qval <0.1],
                  'double_knockout_gate_cell_data' = row.names(double_knockout_gate_cell_data_genes)[double_knockout_gate_cell_data_genes$qval <0.1])

pdf(paste(main_fig_dir, "knockout_gates_upset_groups_all.pdf", sep = ''), height = 3, width = 6)
upset(fromList(listInput), nsets = 5, order.by = "freq") + nm_theme()
dev.off()

listInput_types <- list('knockout cells' = row.names(URMM_all_abs_genotype_beam_genes)[URMM_all_abs_genotype_beam_genes$qval <0.1], 'gates' = row.names(URMM_all_abs_genotype3.1_beam_genes)[URMM_all_abs_genotype3.1_beam_genes$qval <0.1],
                        'knockout cells + gates' = row.names(URMM_all_abs_genotype2_beam_genes)[URMM_all_abs_genotype2_beam_genes$qval <0.1],
                        'double knockout + gates' = row.names(URMM_all_abs_genotype3_beam_genes)[URMM_all_abs_genotype3_beam_genes$qval <0.1],
                        'Irf8 knockout' = row.names(URMM_all_abs_genotype4_beam_genes)[URMM_all_abs_genotype4_beam_genes$qval <0.1],
                        'Gif1 knockout' = row.names(URMM_all_abs_genotype5_beam_genes)[URMM_all_abs_genotype5_beam_genes$qval <0.1],
                        'double knockout' = row.names(URMM_all_abs_genotype6_beam_genes)[URMM_all_abs_genotype6_beam_genes$qval <0.1],
                        'GG1' = row.names(URMM_all_abs_genotype7_beam_genes)[URMM_all_abs_genotype7_beam_genes$qval <0.1],
                        'IG2' = row.names(URMM_all_abs_genotype8_beam_genes)[URMM_all_abs_genotype8_beam_genes$qval <0.1])

listInput_types.1 <- list('knockout cells' = row.names(URMM_all_abs_genotype_beam_genes)[URMM_all_abs_genotype_beam_genes$qval <0.1], 'gates' = row.names(URMM_all_abs_genotype3.1_beam_genes)[URMM_all_abs_genotype3.1_beam_genes$qval <0.1],
                          'knockout cells + gates' = row.names(URMM_all_abs_genotype2_beam_genes)[URMM_all_abs_genotype2_beam_genes$qval <0.1],
                          'double knockout + gates' = row.names(URMM_all_abs_genotype3_beam_genes)[URMM_all_abs_genotype3_beam_genes$qval <0.1],
                          'Irf8 knockout' = row.names(URMM_all_abs_genotype4.1_beam_genes)[URMM_all_abs_genotype4.1_beam_genes$qval <0.1],
                          'Gif1 knockout' = row.names(URMM_all_abs_genotype5.1_beam_genes)[URMM_all_abs_genotype5.1_beam_genes$qval <0.1],
                          'double knockout' = row.names(URMM_all_abs_genotype6.1_beam_genes)[URMM_all_abs_genotype6.1_beam_genes$qval <0.1],
                          'GG1' = row.names(URMM_all_abs_genotype7.1_beam_genes)[URMM_all_abs_genotype7.1_beam_genes$qval <0.1],
                          'IG2' = row.names(URMM_all_abs_genotype8.1_beam_genes)[URMM_all_abs_genotype8.1_beam_genes$qval <0.1])

pdf(paste(main_fig_dir, "knockout_gates_upset_individual_test_all_types.pdf", sep = ''), height = 4, width = 8)
upset(fromList(listInput_types), nsets = 9, order.by = "freq") + nm_theme()
dev.off()

# save.image('./RData/fig5_tmp.RData') #change to fig5
save.image('/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/RData/fig5_tmp.RData')
#check the setdiff of BEAM genes from wild type vs that from all the data
###################################################################################################################################################################################################
#test on carefully on the set of genes from the manifold test:
intersect(listInput_types$GG1, listInput_types$IG2) #0
intersect(listInput_types$`Irf8 knockout`, listInput_types$`Gif1 knockout`) #32; each single knockout: 195 genes
intersect(listInput_types$`double knockout`, listInput_types$`Irf8 knockout`) #39; double knockout: 226 genes
intersect(listInput_types$`double knockout`, listInput_types$`Gif1 knockout`) #31
intersect(listInput_types$`double knockout`, union(listInput_types$`Irf8 knockout`, listInput_types$`Gif1 knockout`)) #60
intersect(listInput_types$`double knockout`, intersect(listInput_types$`Irf8 knockout`, listInput_types$`Gif1 knockout`)) #10

#target genes:
target_vec <- c('Irf8', 'Gfi1', 'Klf4', 'Per3', 'Zeb2', 'Ets1', 'Irf5', 'Csf1r')
intersect(listInput_types$`double knockout`, target_vec) # "Gfi1" "Irf8"
intersect(listInput_types$`Irf8 knockout`, target_vec) # "Irf8"
intersect(listInput_types$`Gif1 knockout`, target_vec) #  "Ets1" "Gfi1"

write.table(file = 'Irf8_knockout.txt', listInput_types$`Irf8 knockout`, quote = F, row.names = F)
write.table(file = 'Gif1_knockout.txt', listInput_types$`Gif1 knockout`, quote = F, row.names = F)

###################################################################################################################################################################################################
#branch heatmap:
cole_markers <- c("Cebpe", "Fosb", "Fos", "Gata2", "Ikzf2", "Klf2", "Gfi1", "Egr1", "Jun", "Zeb1", "Cebpe", "Crebbp", "Gcfc2", "Dach1", "Etv5", "Irf8", "Mitf", "Mycn", "Zeb1", "Csrnp2", "Egr1", "Fosl2", "Irf7", "Irf8", "Satb1", "Tcf4", "Tcf7l2", "Tada2b", "Zscan26", "Zfp60")
valid_cole_markers <- unique(cole_markers[cole_markers %in% row.names(URMM_all_abs)])
plot_genes_in_pseudotime(URMM_all_abs[cole_markers, ], nrow = 5, ncol = 6)
plot_genes_branched_pseudotime(URMM_all_abs[valid_cole_markers, ], branch_point = fig1b_gm_branch_point, color_by = 'Genotype', nrow = 5, ncol = 6)
plot_genes_branched_pseudotime(URMM_all_abs[valid_cole_markers[c(21, 22)], ], branch_point = fig1b_gm_branch_point, color_by = 'Genotype')

plot_genes_branched_pseudotime(URMM_all_abs[c('Gfi1', 'Irf8'), ], branch_point = fig1b_gm_branch_point, color_by = 'Genotype', nrow = 5, ncol = 6)

###################################################################################################################################################################################################
################################################################################################################################################################################################
#function to generate new cds:
make_cds <- function (exprs_matrix, pd, fd, expressionFamily) {
  cds <- newCellDataSet(exprs_matrix,
                        phenoData = new("AnnotatedDataFrame", data = pd),
                        featureData = new("AnnotatedDataFrame", data = fd),
                        expressionFamily = expressionFamily,
                        lowerDetectionLimit = 0.1)

  if(identical(expressionFamily, tobit())) {
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
  }
  return(cds)
}

################################################################################################################################################################################################
# gene_ids <- unique(c(listInput_types$`double knockout`, listInput_types$`Irf8 knockout`, listInput_types$`Gif1 knockout`, 'Irf8', 'Gfi1'))
# abs_URMM_mat <- as.matrix(t(exprs(URMM_all_abs)[gene_ids, order(pData(URMM_all_abs)$Pseudotime)]))
#
# abs_parallel_res_lineage1 <- parallelCCM(ordered_exprs_mat = abs_URMM_mat[pData(URMM_all_abs)$State %in% c(2, 3), ], cores = detectCores() - 2)
# abs_parallel_res_lineage2 <- parallelCCM(ordered_exprs_mat = abs_URMM_mat[pData(URMM_all_abs)$State %in% c(1, 2), ], cores = detectCores() - 2)
#
# abs_parallel_res_mat1 <- prepare_ccm_res(abs_parallel_res_lineage1, gene_names = gene_ids)
# abs_parallel_res_mat2 <- prepare_ccm_res(abs_parallel_res_lineage2)
#
# pdf('abs_ccm.pdf', width = 30, height = 30)
# pheatmap::pheatmap(abs_parallel_res_mat1[, ], useRaster = T, cluster_cols = T, cluster_rows = T) #, annotation_col = F, annotation_row = F
# dev.off()
#
# abs_parallel_res_mat1 <- prepare_ccm_res(abs_parallel_res)
#
# std_parallel_res <- parallelCCM(ordered_exprs_mat = std_lung_mat[, gene_ids], cores = detectCores())
# std_parallel_res_mat <- prepare_ccm_res(std_parallel_res)
#
# pdf('std_ccm.pdf', width = 30, height = 30)
# pheatmap::pheatmap(std_parallel_res_mat1[, ], useRaster = T, cluster_cols = T, cluster_rows = T, annotation_col = F, annotation_row = F)
# dev.off()

################################################################################################################################################################################################
#ensure we use correct states for doing BEAM test

plot_cell_trajectory(URMM_all_abs)
plot_cell_trajectory(URMM_all_fig1b)
plot_cell_trajectory(URMM_all_fig1b, color_by = 'cluster')

fig1b_beam_states <- as.numeric(apply(table(pData(URMM_all_fig1b)[, c('cluster', 'State')])[c('Mono', 'Myelocyte'), ], 1, function(x) which.max(as.numeric(x))))
all_beam_states <- as.numeric(apply(table(pData(URMM_all_abs)[, c('paper_cluster', 'State')])[c('Mono', 'Myelocyte'), ], 1, function(x) which.max(as.numeric(x))))
#calculate branch time point:
#perform BEAM branchtimepoint calculation:
all_gene_ILRs_list <- calILRs(cds = URMM_all_abs[, ], branch_point = fig1b_gm_branch_point, stretch = T, cores = detectCores() - 2, #beam_genes
                              trend_formula = "~sm.ns(Pseudotime, df = 3) * Branch", ILRs_limit = 3, relative_expr = T, label_by_short_name = F,
                              useVST = F, round_exprs = FALSE, output_type = "all", file = "all_gene_ILRs_list", return_all = T)

calILRs(cds = URMM_all_abs[c('Irf8', 'Gfi1'), ], branch_point = fig1b_gm_branch_point, stretch = T, cores = detectCores() - 2, #beam_genes
        trend_formula = "~sm.ns(Pseudotime, df = 3) * Branch", ILRs_limit = 3, relative_expr = T, label_by_short_name = F,
        useVST = T, round_exprs = FALSE, output_type = "all", file = "all_gene_ILRs_list", return_all = T)

all_gene_ILRs_list_fig1b <- calILRs(cds = URMM_all_fig1b[, ], branch_point = fig1b_gm_branch_point, stretch = T, cores = detectCores() - 2, #beam_genes
                                    trend_formula = "~sm.ns(Pseudotime, df = 3) * Branch", ILRs_limit = 3, relative_expr = T, label_by_short_name = F,
                                    useVST = F, round_exprs = FALSE, output_type = "all", file = "all_gene_ILRs_list", return_all = T)

#get the corresponding branch time for the branchpoint of the MST tree:
# fig1b_beam_genes URMM_all_abs_beam_genes
beam_genes <- row.names(subset(URMM_all_abs_beam_genes, qval < 0.1))
cds_full <- buildBranchCellDataSet(URMM_all_abs[unique(beam_genes), ], progenitor_method = 'duplicate', branch_states = all_beam_states, stretch = T, branch_labels = NULL)
cds_fig1b <- buildBranchCellDataSet(URMM_all_fig1b[unique(beam_genes), ], progenitor_method = 'duplicate', branch_states = fig1b_beam_states, stretch = T, branch_labels = NULL)

#check the longest pseudotime for common progenitors:
range(pData(cds_full[, pData(cds_full)$State == 1])$Pseudotime)

abs_bifurcation_time <- abs(detectBifurcationPoint(all_gene_ILRs_list$norm_str_logfc_df[, 60:100], ILRs_threshold = 0.1, return_cross_point = F)) + 60#
names(abs_bifurcation_time) <- fData(URMM_all_abs[names(abs_bifurcation_time), ])$gene_short_name

range(pData(cds_fig1b[, pData(cds_fig1b)$State == 1])$Pseudotime)
abs_bifurcation_time_fig1b <- abs(detectBifurcationPoint(all_gene_ILRs_list_fig1b$norm_str_logfc_df[, 50:100], ILRs_threshold = 0.1, return_cross_point = F)) + 50 #
names(abs_bifurcation_time_fig1b) <- fData(URMM_all_fig1b[names(abs_bifurcation_time_fig1b), ])$gene_short_name

#don't subset cds
abs_bifurcation_time2 <- abs(detectBifurcationPoint(all_gene_ILRs_list$norm_str_logfc_df[, 0:100], ILRs_threshold = 0.1, return_cross_point = F))#
names(abs_bifurcation_time2) <- fData(URMM_all_abs[names(abs_bifurcation_time2), ])$gene_short_name

abs_bifurcation_time_fig1b2 <- abs(detectBifurcationPoint(all_gene_ILRs_list_fig1b$norm_str_logfc_df[, 0:100], ILRs_threshold = 0.1, return_cross_point = F)) #
names(abs_bifurcation_time_fig1b2) <- fData(URMM_all_fig1b[names(abs_bifurcation_time_fig1b2), ])$gene_short_name

################################################################################################################################################################################################
#make the barplot:
#get the Irf8/Gfi1 targets:
Irf8_targets <- read.table('./csv_data/OG_cancer/Irf8_targets.txt') #2290
Gfi1_targets <- read.table('./csv_data/OG_cancer/Gfi_targets.txt') #9423

load('./RData/mouse_fd')
row.names(mouse_fd) <- str_split_fixed(row.names(mouse_fd), '\\.', 2)[, 1]

Irf8_trim_ensembl_id <- str_split_fixed(Irf8_targets$V1[-1], '\\.', 2)[, 1]
Gfi1_trim_ensembl_id <- str_split_fixed(Gfi1_targets$V1[-1], '\\.', 2)[, 1]

Irf8_targets_gene_short_name <- as.character(mouse_fd[Irf8_trim_ensembl_id, 'gene_short_name'])
Gfi1_targets_gene_short_name <- as.character(mouse_fd[Gfi1_trim_ensembl_id, 'gene_short_name'])

Bcell_TF_5k_enrichment_gsc <- loadGSCSafe('./csv_data/OG_cancer/Bcell_Intersect_PKs.bed_JASPAR_5kb_hits_olap.gmt', encoding="latin1")
CD4cell_TF_5k_enrichment_gsc <- loadGSCSafe('./csv_data/OG_cancer/CD4cell_Intersect_PKs.bed_JASPAR_5kb_hits_olap.gmt', encoding="latin1")
CD8_TF_5k_enrichment_gsc <- loadGSCSafe('./csv_data/OG_cancer/CD8_PKs.bed_JASPAR_5kb_hits_olap.gmt', encoding="latin1")
CMP_TF_5k_enrichment_gsc <- loadGSCSafe('./csv_data/OG_cancer/CMP_Intersect_PKs.bed_JASPAR_5kb_hits_olap.gmt', encoding="latin1")
EryA_TF_5k_enrichment_gsc <- loadGSCSafe('./csv_data/OG_cancer/EryA_PKs.bed_JASPAR_5kb_hits_olap.gmt', encoding="latin1")
GMP_TF_5k_enrichment_gsc <- loadGSCSafe('./csv_data/OG_cancer/GMP_Intersect_PKs.bed_JASPAR_5kb_hits_olap.gmt', encoding="latin1")
Granulocytes_TF_5k_enrichment_gsc <- loadGSCSafe('./csv_data/OG_cancer/Granulocytes_Intersect_PKs.bed_JASPAR_5kb_hits_olap.gmt', encoding="latin1")
Lsk_TF_5k_enrichment_gsc <- loadGSCSafe('./csv_data/OG_cancer/Lsk_Intersect_PKs.bed_JASPAR_5kb_hits_olap.gmt', encoding="latin1")
MEP_TF_5k_enrichment_gsc <- loadGSCSafe('./csv_data/OG_cancer/MEP_Intersect_PKs.bed_JASPAR_5kb_hits_olap.gmt', encoding="latin1")
Monocytes_TF_5k_enrichment_gsc <- loadGSCSafe('./csv_data/OG_cancer/Monocytes_PKs.bed_JASPAR_5kb_hits_olap.gmt', encoding="latin1")
NK_TF_5k_enrichment_gsc <- loadGSCSafe('./csv_data/OG_cancer/NK_Intersect_PKs.bed_JASPAR_5kb_hits_olap.gmt', encoding="latin1")

#intersect Irf8 / Gfi1 target genes from chip-seq with BEAM genes:
URMM_all_abs_beam_genes_vec <- row.names(URMM_all_abs_beam_genes)[URMM_all_abs_beam_genes$qval < 1e-1]

intersect(Irf8_targets_gene_short_name, URMM_all_abs_beam_genes_vec)
intersect(Gfi1_targets_gene_short_name, URMM_all_abs_beam_genes_vec)

length(intersect(c(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name), URMM_all_abs_beam_genes_vec)) / length(URMM_all_abs_beam_genes_vec)
#0.5805981

tmp <- str_split_fixed(names(GMP_TF_5k_enrichment_gsc$gsc), "\\(", 2)[, 1]
GMP_TF_names <- unique(c(str_split_fixed(tmp, "::", 3)))[-1]
tmp <- str_split_fixed(names(Monocytes_TF_5k_enrichment_gsc$gsc), "\\(", 2)[, 1]
Monocytes_TF_names <- unique(c(str_split_fixed(tmp, "::", 3)))[-1]
tmp <- str_split_fixed(names(Granulocytes_TF_5k_enrichment_gsc$gsc), "\\(", 2)[, 1]
Granulocytes_TF_names <- unique(c(str_split_fixed(tmp, "::", 3)))[-1]

c(Monocytes_TF_names, Granulocytes_TF_names, GMP_TF_names)

################################################################################################################################################################################################
return_category <- function(checked_TFs_names, checked_TFs_names2,  Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name, Gfi1_targets_names, Irf8_targets_names){
  other_genes <- setdiff(checked_TFs_names, c(Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name))
  both_genes <- intersect(Gfi1_targets_names, Irf8_targets_names)
  Gfi1_genes <- setdiff(Gfi1_targets_names, Irf8_targets_names)
  Irf8_genes <- setdiff(Irf8_targets_names, Gfi1_targets_names)

  other_genes_ids <- which(toupper(checked_TFs_names2) %in% toupper(other_genes))
  both_genes_ids <- which(toupper(checked_TFs_names2) %in% toupper(both_genes))
  Gfi1_genes_ids <- which(toupper(checked_TFs_names2) %in% c(toupper(Gfi1_genes)))
  Irf8_genes_ids <- which(toupper(checked_TFs_names2) %in% c(toupper(Irf8_genes), "IRF8"))

  Gfi1_self_id <- which(toupper(checked_TFs_names2) %in% c( "GFI1"))
  Irf8_self_id <- which(toupper(checked_TFs_names2) %in% c( "IRF8"))

  category <- rep("both", length(checked_TFs_names2))
  category[other_genes_ids] <- 'other'
  category[Gfi1_genes_ids] <- 'Gfi1'
  category[Irf8_genes_ids] <- 'Irf8'

  category[Gfi1_self_id] <- 'Gfi1_self'
  category[Irf8_self_id] <- 'Irf8_self'

  return(category)
}

#gene name based on WT cells tree only
beam_genes <- row.names(subset(fig1b_beam_genes_proj_dup, qval < 1e-1)) # fig1b_beam_genes
beam_genes_all <- row.names(subset(URMM_all_abs_beam_genes_proj_dup, qval < 1e-1)) # fig1b_beam_genes

#Irf8/Gfi1's direct targets:
intersect(Irf8_targets_gene_short_name, URMM_all_abs_beam_genes_vec)
intersect(Gfi1_targets_gene_short_name, URMM_all_abs_beam_genes_vec)

checked_TFs_jaspar_fig1b <- GMP_TF_names[toupper(GMP_TF_names) %in% toupper(beam_genes)]
checked_TFs_names_fig1b <- beam_genes[toupper(beam_genes) %in% toupper(GMP_TF_names)]

#remove 'Irf8', 'Gfi1', "Klf4"
checked_TFs_jaspar_fig1b <- setdiff(checked_TFs_jaspar_fig1b, c("Klf4")) #'IRF8', 'Gfi1',
checked_TFs_names_fig1b <- setdiff(checked_TFs_names_fig1b, c("Klf4")) #'Irf8', 'Gfi1',

#
# checked_TFs_jaspar_fig1b <- GMP_TF_names[toupper(GMP_TF_names) %in% toupper(WT_beam_genes)]
# checked_TFs_names_fig1b <- c("Bcl2", "Bcl2a1b", "Bcl2a1d", "Cebpe", "Elk3", "Ets1", "Hoxa9", "Id2", "Gfi1", "Irf5", "Irf9", "Kif3b", "Klf4", "Mef2d", # "Irf8",
#                        "Meis1", "Pou2f2", "Prdm5", "Rhox8", "Sox12", "Sox4", "Sp140", "Tcf7l2", "Unc119", "Unc119b", "Zbed3",  "Zbtb48", "Zeb1",
#                        "Zeb2") #"Tcf19", "Zbtb2",
#
# #remove outlier genes with much late branch time points
# checked_TFs_jaspar_fig1b_valid_ids <- which(!(checked_TFs_jaspar_fig1b %in% c('Hoxa', 'ZBTB7', 'Rhox1', 'ELK', 'BCL6', 'PRDM', 'MEIS')))        #Rhox1      ELK     BCL6     PRDM     MEIS
# checked_TFs_names_fig1b_valid_ids <- which(!(checked_TFs_jaspar_fig1b %in% c('Unc119', 'Mef2d', 'MEIS', 'Rhox8', 'Elk', 'Prdm5')))
#
# checked_TFs_jaspar_fig1b <- checked_TFs_jaspar_fig1b[checked_TFs_jaspar_fig1b_valid_ids]
# checked_TFs_names_fig1b <- checked_TFs_names_fig1b[checked_TFs_names_fig1b_valid_ids]

#Irf8/Gfi1's secondary targets:
cmp_sets <- unlist(lapply(checked_TFs_jaspar_fig1b, function(x) grep(x, names(GMP_TF_5k_enrichment_gsc$gsc))))
target_genes <- setdiff(intersect(unlist(GMP_TF_5k_enrichment_gsc$gsc[cmp_sets]), beam_genes), checked_TFs_names_fig1b) #URMM_all_abs_beam_genes_vec

################################################################################################################################
#get the mean branching time point for each direct targets' targets (secondary targets)
cmp_sets_time <- unlist(lapply(checked_TFs_jaspar_fig1b, function(x) {
  cmp_sets <- grep(x, names(GMP_TF_5k_enrichment_gsc$gsc))
  target_genes <- setdiff(intersect(unlist(GMP_TF_5k_enrichment_gsc$gsc[cmp_sets]), beam_genes), checked_TFs_names_fig1b)
  mean_bif_time <- mean(abs_bifurcation_time_fig1b[target_genes], na.rm = T)
}))
names(cmp_sets_time) <- checked_TFs_jaspar_fig1b
################################################################################################################################

Monocytes_sets <- unlist(lapply(checked_TFs_jaspar_fig1b, function(x) grep(x, names(Monocytes_TF_5k_enrichment_gsc$gsc))))
Monocytes_target_genes <- setdiff(intersect(unlist(Monocytes_TF_5k_enrichment_gsc$gsc[Monocytes_sets]), beam_genes), checked_TFs_names_fig1b) #URMM_all_abs_beam_genes_vec
Granulocytes_sets <- unlist(lapply(checked_TFs_jaspar_fig1b, function(x) grep(x, names(Granulocytes_TF_5k_enrichment_gsc$gsc))))
Granulocytess_target_genes <- setdiff(intersect(unlist(Granulocytes_TF_5k_enrichment_gsc$gsc[Granulocytes_sets]), beam_genes), checked_TFs_names_fig1b) #URMM_all_abs_beam_genes_vec

#perform the enrichment analysis:
#add Irf8 / Gfi targets:
Irf8_targets_names <- intersect(checked_TFs_names_fig1b, Irf8_targets_gene_short_name)
Gfi1_targets_names <- intersect(checked_TFs_names_fig1b, Gfi1_targets_gene_short_name)

#genes not binding by Irf8 / Gfi1:
setdiff(checked_TFs_names_fig1b, c(Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name))
# [1]  "Egr3" "Klf4"
intersect(Gfi1_targets_names, Irf8_targets_names)
# [1]  "Tcf7l2"
setdiff(Gfi1_targets_names, Irf8_targets_names)
# "Cebpe" "Ets1"  "Irf8"  "Sox4"
setdiff(Irf8_targets_names, Gfi1_targets_names)
# NULL

checked_TFs_jaspar_fig1b_category <- return_category(checked_TFs_names_fig1b, checked_TFs_jaspar_fig1b,  Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name, Gfi1_targets_names, Irf8_targets_names)
checked_TFs_names_fig1b_category <- return_category(checked_TFs_names_fig1b, checked_TFs_names_fig1b,  Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name, Gfi1_targets_names, Irf8_targets_names)

# checked_TFs_jaspar_fig1b_category <- checked_TFs_jaspar_fig1b_category[checked_TFs_jaspar_fig1b_valid_ids]
# checked_TFs_names_fig1b_category <- checked_TFs_names_fig1b_category[checked_TFs_names_fig1b_valid_ids]

Gfi1_cmp_sets <- unlist(lapply(checked_TFs_jaspar_fig1b[checked_TFs_jaspar_fig1b_category == 'Gfi1'], function(x) grep(x, names(GMP_TF_5k_enrichment_gsc$gsc))))
Gfi1_secondary_target_genes <- setdiff(intersect(unlist(GMP_TF_5k_enrichment_gsc$gsc[Gfi1_cmp_sets]), beam_genes), checked_TFs_names_fig1b) #URMM_all_abs_beam_genes_vec

Irf8_cmp_sets <- unlist(lapply(checked_TFs_jaspar_fig1b[checked_TFs_jaspar_fig1b_category == 'Irf8'], function(x) grep(x, names(GMP_TF_5k_enrichment_gsc$gsc))))
Irf8_secondary_target_genes <- setdiff(intersect(unlist(GMP_TF_5k_enrichment_gsc$gsc[Irf8_cmp_sets]), beam_genes), checked_TFs_names_fig1b) #URMM_all_abs_beam_genes_vec

both_cmp_sets <- unlist(lapply(checked_TFs_jaspar_fig1b[checked_TFs_jaspar_fig1b_category == 'both'], function(x) grep(x, names(GMP_TF_5k_enrichment_gsc$gsc))))
both_secondary_target_genes <- setdiff(intersect(unlist(GMP_TF_5k_enrichment_gsc$gsc[both_cmp_sets]), beam_genes), checked_TFs_names_fig1b) #URMM_all_abs_beam_genes_vec

other_cmp_sets <- unlist(lapply(checked_TFs_jaspar_fig1b[checked_TFs_jaspar_fig1b_category == 'other'], function(x) grep(x, names(GMP_TF_5k_enrichment_gsc$gsc))))
other_secondary_target_genes <- setdiff(intersect(unlist(GMP_TF_5k_enrichment_gsc$gsc[other_cmp_sets]), beam_genes), checked_TFs_names_fig1b) #URMM_all_abs_beam_genes_vec

Gfi1_self_cmp_sets <- unlist(lapply(checked_TFs_jaspar_fig1b[checked_TFs_jaspar_fig1b_category == 'Gfi1_self'], function(x) grep(x, names(GMP_TF_5k_enrichment_gsc$gsc))))
Gfi1_self_secondary_target_genes <- setdiff(intersect(unlist(GMP_TF_5k_enrichment_gsc$gsc[Gfi1_self_cmp_sets]), beam_genes),
                                            c(checked_TFs_names_fig1b, Gfi1_secondary_target_genes, Irf8_secondary_target_genes, both_secondary_target_genes, other_secondary_target_genes) ) #URMM_all_abs_beam_genes_vec

Irf8_self_cmp_sets <- unlist(lapply(checked_TFs_jaspar_fig1b[checked_TFs_jaspar_fig1b_category == 'Irf8_self'], function(x) grep(x, names(GMP_TF_5k_enrichment_gsc$gsc))))
Irf8_self_secondary_target_genes <- setdiff(intersect(unlist(GMP_TF_5k_enrichment_gsc$gsc[Irf8_self_cmp_sets]), beam_genes),
                                            c(checked_TFs_names_fig1b, Irf8_secondary_target_genes, both_secondary_target_genes, other_secondary_target_genes)) #URMM_all_abs_beam_genes_vec

#ignore the targets in the Mono / Gran lineages:
Monocytes_target_genes <- NULL
Granulocytess_target_genes <- NULL

secondary_target_genes <- c(Gfi1_secondary_target_genes, Irf8_secondary_target_genes, both_secondary_target_genes,
                            Monocytes_target_genes, Granulocytess_target_genes)

#direct target branch time point
top_group_branch_time <- abs_bifurcation_time_fig1b2[c(checked_TFs_names_fig1b, Irf8_self_secondary_target_genes, Gfi1_self_secondary_target_genes, other_secondary_target_genes)] #abs_bifurcation_time
#secondary target branch time point
bottom_group_branch_time <- abs_bifurcation_time_fig1b2[secondary_target_genes] #abs_bifurcation_time_fig1b

df <- data.frame(labels = c(c(checked_TFs_names_fig1b, Irf8_self_secondary_target_genes, Gfi1_self_secondary_target_genes, other_secondary_target_genes), secondary_target_genes),
                 Gfi1_Irf8_type = c(checked_TFs_names_fig1b_category,
                                    rep('Gfi1_self', length(Irf8_self_secondary_target_genes)),
                                    rep('Irf8_self', length(Gfi1_self_secondary_target_genes)),
                                    rep('other', length(other_secondary_target_genes)),
                                    rep('Gfi1', length(Gfi1_secondary_target_genes)),
                                    rep('Irf8', length(Irf8_secondary_target_genes)),
                                    rep('both', length(both_secondary_target_genes)),
                                    rep('Gfi1', length(Monocytes_target_genes)),
                                    rep('Irf8', length(Granulocytess_target_genes))
                 ),
                 bifurcation_time = c(top_group_branch_time, bottom_group_branch_time),
                 type = c(rep("Direct targets", length(top_group_branch_time)), rep("Secondary targets", length(bottom_group_branch_time))))

df[, 'type'] <- as.character(df[, 'type'])
df[which(df$labels %in% c('Irf8', 'Gfi1')), 'type'] <- 'Master regulator'
df[which(df$Gfi1_Irf8_type %in% 'other'), 'type'] <- 'Other BEAM genes'
df_all <- rbind(df, data.frame(labels = setdiff(beam_genes, df$labels),
                               Gfi1_Irf8_type = 'other',
                               bifurcation_time = abs_bifurcation_time_fig1b2[setdiff(beam_genes, df$labels)],
                               type = 'Other BEAM genes'))

df_all[, 'type'] <- factor(df_all[, 'type'], levels = rev(c('Master regulator', 'Direct targets', 'Secondary targets', 'Other BEAM genes')))

qplot(type, abs(bifurcation_time), color = type, geom = c('jitter'), data = subset(df_all, abs(bifurcation_time) > 60), alpha = I(0.7), log = 'y', size = 1) +  #, 'boxplot'
  xlab('') + ylab('bifurcation time point') + coord_flip() + scale_size(range = c(0.01, 1), limits = c(0.1, 1)) +
  geom_boxplot(aes(color = type), fatten = 0.5, lwd = 0.5, outlier.shape=NA, alpha = 0.5) +
  scale_y_continuous(breaks = round(seq(min(abs(df$bifurcation_time), na.rm = T), max(abs(df$bifurcation_time) + 5, na.rm = T), by = 9),1)) + #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`))
  nm_theme()

subset_df <- subset(df_all, abs(bifurcation_time) > 42 & !(labels %in% c('Gfi1', 'Irf8')))
subset_df <- rbind(subset_df, subset(df_all, labels %in% c('Gfi1', 'Irf8')))
pdf(paste(main_fig_dir, 'fig4g_new.pdf', sep = ''), width = 2.2, height = 1.5)
qplot(type, abs(bifurcation_time), color = type, geom = c('jitter'), data = subset_df, alpha = I(0.7), log = 'y', size = 1) +  #, 'boxplot'
  xlab('') + ylab('bifurcation time point') + coord_flip() + scale_size(range = c(0.01, 1), limits = c(0.1, 1)) +
  geom_boxplot(aes(color = type), fatten = 0.5, lwd = 0.5, outlier.shape=NA, alpha = 0.5) +
  scale_y_continuous(breaks = round(seq(min(abs(df$bifurcation_time), na.rm = T), max(abs(df$bifurcation_time) + 5, na.rm = T), by = 9),1)) + #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`))
  nm_theme()
dev.off()

################################################################################################################################################################################################
# #do this for all the dataset:
#
# #overlapping of BEAM genes, enriched TFs:
# checked_TFs_jaspar <- GMP_TF_names[toupper(GMP_TF_names) %in% toupper(beam_genes_all)]
# checked_TFs_names <- beam_genes_all[toupper(beam_genes_all) %in% toupper(GMP_TF_names)]
#
# #get the list of targets:
# cmp_sets <- unlist(lapply(checked_TFs_jaspar, function(x) grep(x, names(GMP_TF_5k_enrichment_gsc$gsc))))
# target_genes <- setdiff(intersect(unlist(GMP_TF_5k_enrichment_gsc$gsc[cmp_sets]), beam_genes), checked_TFs_names) #URMM_all_abs_beam_genes_vec
#
# Monocytes_sets <- unlist(lapply(checked_TFs_jaspar, function(x) grep(x, names(Monocytes_TF_5k_enrichment_gsc$gsc))))
# Monocytes_target_genes <- setdiff(intersect(unlist(Monocytes_TF_5k_enrichment_gsc$gsc[cmp_sets]), beam_genes), checked_TFs_names) #URMM_all_abs_beam_genes_vec
# Granulocytes_sets <- unlist(lapply(checked_TFs_jaspar, function(x) grep(x, names(Granulocytes_TF_5k_enrichment_gsc$gsc))))
# Granulocytess_target_genes <- setdiff(intersect(unlist(Granulocytes_TF_5k_enrichment_gsc$gsc[cmp_sets]), beam_genes), checked_TFs_names) #URMM_all_abs_beam_genes_vec
#
# #perform the enrichment analysis:
# #add Irf8 / Gfi targets:
# Irf8_targets_names <- intersect(checked_TFs_names, Irf8_targets_gene_short_name)
# Gfi1_targets_names <- intersect(checked_TFs_names, Gfi1_targets_gene_short_name)
#
# #genes not binding by Irf8 / Gfi1:
# setdiff(checked_TFs_names, c(Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name))
# #  "Atf3"  "Hoxa5" "Hoxa9" "Klf4"  "Mef2d" "Vdr"
# intersect(Gfi1_targets_names, Irf8_targets_names)
# # [1]  "Id2"    "Mef2c"  "Pou2f2" "Tcf4"   "Tcf7l2"
# setdiff(Gfi1_targets_names, Irf8_targets_names)
# # "Cebpe"  "Elk3"   "Ets1"   "Gata2"  "Hoxa10" "Irf8"   "Mecom"  "Sox4"   "Tal1"
# setdiff(Irf8_targets_names, Gfi1_targets_names)
# # NULL
#
# #
#
# # target TFs:
# checked_TFs_jaspar_category <- return_category(checked_TFs_names, checked_TFs_jaspar,  Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name, Gfi1_targets_names, Irf8_targets_names)
# checked_TFs_names_category <- return_category(checked_TFs_names, checked_TFs_names,  Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name, Gfi1_targets_names, Irf8_targets_names)
#
# Gfi1_cmp_sets <- unlist(lapply(checked_TFs_jaspar[checked_TFs_jaspar_category == 'Gfi1'], function(x) grep(x, names(GMP_TF_5k_enrichment_gsc$gsc))))
# Gfi1_secondary_target_genes <- setdiff(intersect(unlist(GMP_TF_5k_enrichment_gsc$gsc[Gfi1_cmp_sets]), beam_genes), checked_TFs_names) #URMM_all_abs_beam_genes_vec
#
# Irf8_cmp_sets <- unlist(lapply(checked_TFs_jaspar[checked_TFs_jaspar_category == 'Irf8'], function(x) grep(x, names(GMP_TF_5k_enrichment_gsc$gsc))))
# Irf8_secondary_target_genes <- setdiff(intersect(unlist(GMP_TF_5k_enrichment_gsc$gsc[Irf8_cmp_sets]), beam_genes), checked_TFs_names) #URMM_all_abs_beam_genes_vec
#
# both_cmp_sets <- unlist(lapply(checked_TFs_jaspar[checked_TFs_jaspar_category == 'both'], function(x) grep(x, names(GMP_TF_5k_enrichment_gsc$gsc))))
# both_secondary_target_genes <- setdiff(intersect(unlist(GMP_TF_5k_enrichment_gsc$gsc[both_cmp_sets]), beam_genes), checked_TFs_names) #URMM_all_abs_beam_genes_vec
#
# other_cmp_sets <- unlist(lapply(checked_TFs_jaspar[checked_TFs_jaspar_category == 'other'], function(x) grep(x, names(GMP_TF_5k_enrichment_gsc$gsc))))
# other_secondary_target_genes <- setdiff(intersect(unlist(GMP_TF_5k_enrichment_gsc$gsc[other_cmp_sets]), beam_genes), checked_TFs_names) #URMM_all_abs_beam_genes_vec
#
# Monocytes_target_genes <- NULL
# Granulocytess_target_genes <- NULL
#
# secondary_target_genes <- (c(Gfi1_secondary_target_genes, Irf8_secondary_target_genes, both_secondary_target_genes, other_secondary_target_genes, Monocytes_target_genes, Granulocytess_target_genes))
#
# top_group_branch_time <- abs_bifurcation_time2[c(checked_TFs_names)] #abs_bifurcation_time_fig1b
# bottom_group_branch_time <- abs_bifurcation_time2[secondary_target_genes]
# df <- data.frame(labels = c(c(checked_TFs_names), secondary_target_genes),
#                  Gfi1_Irf8_type = c(checked_TFs_names_category,
#                                     rep('Gfi1', length(Gfi1_secondary_target_genes)),
#                                     rep('Irf8', length(Irf8_secondary_target_genes)),
#                                     rep('both', length(both_secondary_target_genes)),
#                                     rep('other', length(other_secondary_target_genes)),
#                                     rep('Gfi1', length(Monocytes_target_genes)),
#                                     rep('Irf8', length(Granulocytess_target_genes))
#                  ),
#                  bifurcation_time = c(top_group_branch_time, bottom_group_branch_time),
#                  type = c(rep("Direct targets", length(top_group_branch_time)), rep("Secondary targets", length(bottom_group_branch_time))))
#
# df[, 'type'] <- as.character(df[, 'type'])
# df[which(df$labels %in% c('Irf8', 'Gfi1')), 'type'] <- 'Master regulator'
# df[, 'type'] <- factor(df[, 'type'], levels = rev(c('Master regulator', 'Direct targets', 'Secondary targets')))
#
# ################################################################################################################################################################################################
# pdf('main_figures//fig4g.pdf', width = 2, height = 1)
# qplot(type, abs(bifurcation_time), color = type, geom = c('jitter'), data = subset(df, bifurcation_time > 50), alpha = I(0.7), log = 'y', size = 1) +  #, 'boxplot'
#   xlab('') + ylab('bifurcation time point') + coord_flip() + scale_size(range = c(0.01, 1), limits = c(0.1, 1)) +
#   geom_boxplot(aes(color = type), fatten = 0.5, lwd = 0.5, outlier.shape=NA, alpha = 0.5) +
#   scale_y_continuous(breaks = round(seq(min(abs(df$bifurcation_time), na.rm = T), max(abs(df$bifurcation_time) + 5, na.rm = T), by = 9),1)) + #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`))
#   nm_theme()
# dev.off()
#
# pdf('main_figures//fig4g_helper.pdf')
# qplot(type, abs(bifurcation_time), color = Gfi1_Irf8_type, geom = c('jitter'), data = df, alpha = I(0.7), log = 'y', size = 1) +  #, 'boxplot'
#   xlab('') + ylab('bifurcation time point') + coord_flip() + scale_size(range = c(0.01, 1), limits = c(0.1, 1)) +
#   geom_boxplot(aes(color = Gfi1_Irf8_type), fatten = 0.5, lwd = 0.5, outlier.shape=NA, alpha = 0.5) +
#   scale_y_continuous(breaks = round(seq(min(abs(df$bifurcation_time), na.rm = T), max(abs(df$bifurcation_time) + 5, na.rm = T), by = 9),1))  #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`))
# dev.off()
#
# pdf('main_figures//fig4g2.pdf', width = 3.5, height = 2)
# qplot(type, abs(bifurcation_time), color = Gfi1_Irf8_type, geom = c('jitter'), data = df, alpha = I(0.7), log = 'y', size = 1) +  #, 'boxplot'
#   xlab('') + ylab('bifurcation time point') + coord_flip() + scale_size(range = c(0.01, 1), limits = c(0.1, 1)) +
#   geom_boxplot(aes(color = Gfi1_Irf8_type), fatten = 0.5, lwd = 0.5, outlier.shape=NA, alpha = 0.5) +
#   scale_y_continuous(breaks = round(seq(min(abs(df$bifurcation_time), na.rm = T), max(abs(df$bifurcation_time) + 5, na.rm = T), by = 9),1)) + nm_theme()  #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`))
# dev.off()
#
# ############################################################################################################################################################
# #use all the WT cells for the analysis
# ############################################################################################################################################################
# # source('./scripts/analysis_URMM_all_WT_cells.R', echo = T)
#
# ############################################################################################################################################################
# #create the venn diagram
# ############################################################################################################################################################
#
# intersect(URMM_all_abs_beam_genes_vec, row.names(subset(fig1b_beam_genes, qval < 0.1))) #852, 543, 290
# beam_genes <- row.names(subset(fig1b_beam_genes, qval < 0.1)) #URMM_all_abs_beam_genes_vec
# Irf8_ko <- listInput_types$`Irf8 knockout`
# Gif1_ko <- listInput_types$`Gif1 knockout`
# double_ko <- listInput_types$`double knockout`
#
# element_all <- c(beam_genes,
#                  Irf8_ko,
#                  Gif1_ko,
#                  double_ko)
# sets_all <- c(rep(paste('BEAM genes', sep = ''), length(beam_genes)),
#               rep(paste('Irf8 knockout genes', sep = ''), length(Irf8_ko)),
#               rep(paste('Gfi1 knockout genes', sep = ''), length(Gif1_ko)),
#               rep(paste('Double knockout genes', sep = ''), length(double_ko))
# )
# table(sets_all) #number of genes
#
# intersect(c(Irf8_ko, Gif1_ko), double_ko)
# intersect(c(Irf8_ko, Gif1_ko, double_ko), beam_genes)
# intersect(c(Irf8_ko, Gif1_ko, double_ko), URMM_all_abs_beam_genes_vec)
#
# BEAM_cmpr <- list('BEAM genes' = beam_genes, #URMM_all_abs_beam_genes_vec
#                   'Irf8 knockout' = Irf8_ko,
#                   'Gif1 knockout' = Gif1_ko,
#                   'Double knockout genes' = double_ko
# )
#
# pdf(paste(main_fig_dir, "knockout_gates_upset_individual_test_ko_beam.pdf", sep = ''), height = 4, width = 8)
# upset(fromList(BEAM_cmpr), nsets = 4, order.by = "freq") #+ nm_theme()
# dev.off()
#
# listInput_types <- list('knockout cells' = row.names(URMM_all_abs_genotype_beam_genes)[URMM_all_abs_genotype_beam_genes$qval <0.1],
#                         'gates' = row.names(URMM_all_abs_genotype3.1_beam_genes)[URMM_all_abs_genotype3.1_beam_genes$qval <0.1],
#                         'knockout cells + gates' = row.names(URMM_all_abs_genotype2_beam_genes)[URMM_all_abs_genotype2_beam_genes$qval <0.1],
#                         'double knockout + gates' = row.names(URMM_all_abs_genotype3_beam_genes)[URMM_all_abs_genotype3_beam_genes$qval <0.1],
#                         'Irf8 knockout' = row.names(URMM_all_abs_genotype4_beam_genes)[URMM_all_abs_genotype4_beam_genes$qval <0.1],
#                         'Gif1 knockout' = row.names(URMM_all_abs_genotype5_beam_genes)[URMM_all_abs_genotype5_beam_genes$qval <0.1],
#                         'double knockout' = row.names(URMM_all_abs_genotype6_beam_genes)[URMM_all_abs_genotype6_beam_genes$qval <0.1],
#                         'GG1' = row.names(URMM_all_abs_genotype7_beam_genes)[URMM_all_abs_genotype7_beam_genes$qval <0.1],
#                         'IG2' = row.names(URMM_all_abs_genotype8_beam_genes)[URMM_all_abs_genotype8_beam_genes$qval <0.1])
#
# pdf(paste(main_fig_dir, "knockout_gates_upset_individual_test_all_types.pdf", sep = ''), height = 4, width = 8)
# upset(fromList(listInput_types), nsets = 9, order.by = "freq") #+ nm_theme()
# dev.off()
#
# listInput_types <- list('Irf8 knockout' = row.names(URMM_all_abs_genotype4_beam_genes)[URMM_all_abs_genotype4_beam_genes$qval <0.1],
#                         'Gif1 knockout' = row.names(URMM_all_abs_genotype5_beam_genes)[URMM_all_abs_genotype5_beam_genes$qval <0.1],
#                         'double knockout' = row.names(URMM_all_abs_genotype6_beam_genes)[URMM_all_abs_genotype6_beam_genes$qval <0.1],
#                         'GG1' = row.names(URMM_all_abs_genotype7_beam_genes)[URMM_all_abs_genotype7_beam_genes$qval <0.1],
#                         'IG2' = row.names(URMM_all_abs_genotype8_beam_genes)[URMM_all_abs_genotype8_beam_genes$qval <0.1])
#
# pdf(paste(main_fig_dir, "knockout_gates_upset_individual_test_ko_gates_beam.pdf", sep = ''), height = 4, width = 8)
# upset(fromList(listInput_types), nsets = 5, order.by = "freq") #+ nm_theme()
# dev.off()
#
# pdf(paste(main_fig_dir, 'knockout_gates_upset_individual_test_all_types_venn_diagram.pdf', sep = ''))
# venneuler_venn(element_all, sets_all)
# dev.off()
# table(sets_all) #number of genes
#
# ############################################################################################################################################################
# #create the venn diagram
# ############################################################################################################################################################
# element_all <- c(Gfi1_targets_names,
#                  Irf8_targets_names)
# sets_all <- c(rep(paste('Regulator (Gfi1)', sep = ''), length(Gfi1_targets_names)),
#               rep(paste('Regulator (Irf8)', sep = ''), length(Irf8_targets_names)))
#
# table(sets_all) #number of genes
# intersect(Gfi1_targets_names, Irf8_targets_names)
# setdiff(Gfi1_targets_names, Irf8_targets_names)
# setdiff(Irf8_targets_names, Gfi1_targets_names)
#
# pdf(paste(main_fig_dir, 'fig4h.1.pdf', sep = ''))
# venneuler_venn(element_all, sets_all)
# dev.off()
# table(sets_all) #number of genes
#
# ############################################################################################################################################################
# #create the venn diagram for secondary targets
# ############################################################################################################################################################
# element_all <- c(Irf8_secondary_target_genes, both_secondary_target_genes,
#                  Gfi1_secondary_target_genes)
# sets_all <- c(rep(paste('Regulator (Gfi1)', sep = ''), length(Gfi1_secondary_target_genes)),
#               rep(paste('Regulator (both)', sep = ''), length(both_secondary_target_genes)),
#               rep(paste('Regulator (Irf8)', sep = ''), length(Irf8_secondary_target_genes)))
#
# table(sets_all) #number of genes
#
# pdf(paste(main_fig_dir, 'fig4h.2.pdf', sep = ''))
# venneuler_venn(element_all, sets_all)
# dev.off()
# table(sets_all) #number of genes
#
# ############################################################################################################################################################
# #create the venn diagram for secondary targets
# ############################################################################################################################################################
#make the network:
#1. intersect the chip-seq targets with the knockout manifold test
#2. identifying TFs from the intersection set
#3. find the targets of those TFs

#link Gfi1, Irf8 to the TFs
adj_df <- data.frame(source = c('Gfi1'), target = c('Irf8'),  target_type = c('Gfi1'))
adj_df <- data.frame(source = c('Gfi1'), target = c('Gfi1'),  target_type = c('Gfi1'))
for(jaspar_ind in 1:length(checked_TFs_jaspar_fig1b)){
  jaspar_name <- checked_TFs_jaspar_fig1b[jaspar_ind]
  official_gene_name <- checked_TFs_names_fig1b[toupper(checked_TFs_names_fig1b) %in% toupper(jaspar_name)]

  if(checked_TFs_jaspar_fig1b_category[jaspar_ind] == 'both')
    adj_df = rbind(adj_df, data.frame(source = c("Gfi1", "Irf8"), target = c(official_gene_name), target_type = 'both'))
  else if(checked_TFs_jaspar_fig1b_category[jaspar_ind] == 'Gfi1')
    adj_df = rbind(adj_df, data.frame(source = c("Gfi1"), target = c(official_gene_name), target_type = 'Gfi1'))
  else if(checked_TFs_jaspar_fig1b_category[jaspar_ind] == 'Irf8')
    adj_df = rbind(adj_df, data.frame(source = c("Irf8"), target = c(official_gene_name), target_type = 'Irf8'))
}

#link TFs to their targets (only use Gfi / Irf8 targets)
valid_checked_TFs_jaspar_fig1b <- checked_TFs_jaspar_fig1b[checked_TFs_jaspar_fig1b_category != 'other'] #checked_TFs_jaspar_category

for(jaspar_ind in 1:length(valid_checked_TFs_jaspar_fig1b)){
  message(jaspar_ind)
  tmp_cmp_sets <- unlist(lapply(valid_checked_TFs_jaspar_fig1b[jaspar_ind], function(x) grep(x, names(GMP_TF_5k_enrichment_gsc$gsc))))
  inter_cmp <- intersect(unlist(GMP_TF_5k_enrichment_gsc$gsc[tmp_cmp_sets]), beam_genes)

  if('Tmcc2' %in% inter_cmp){
    print('Tmcc2 is here')
    print(inter_cmp)
  }

  if(length(inter_cmp) > 0){
    tmp_secondary_target_genes <- setdiff(inter_cmp, checked_TFs_names_fig1b) #select genes target by direct targets but ont include direct targets (select all secondary targets)
    jaspar_name <- valid_checked_TFs_jaspar_fig1b[jaspar_ind]
    official_gene_name <- checked_TFs_names_fig1b[toupper(checked_TFs_names_fig1b) %in% toupper(jaspar_name)]
    adj_df = rbind(adj_df, data.frame(source = official_gene_name, target = tmp_secondary_target_genes, target_type = checked_TFs_jaspar_fig1b_category[jaspar_ind])) #checked_TFs_jaspar_category
  }
}

uniq_genes <- setdiff(unique(c(as.character(adj_df[, 1]), as.character(adj_df[, 2]))), c('Klf4', 'IRF8')) #, 'Vdr', 'Hoxa9'
hbp_network_mat <- as.data.frame(matrix(rep(0, length(uniq_genes)^2), nrow = length(uniq_genes), dimnames = list(uniq_genes,  uniq_genes)))
for(ind in 1:nrow(adj_df)){
  message(ind)
  if(all(c(as.character(adj_df[ind, 1]), as.character(adj_df[ind, 2])) %in% uniq_genes))
    hbp_network_mat[as.character(adj_df[ind, 1]), as.character(adj_df[ind, 2])] <- 1
}

g1 <- graph_from_adjacency_matrix(as.matrix(hbp_network_mat), mode = "directed")

setdiff(V(g1)$name, as.character(df$labels))
setdiff(uniq_genes, as.character(df$labels))

# #
# load('./RData/ccm_res')
#
# abs_parallel_res_mat1 <- prepare_ccm_res(abs_parallel_res_lineage1, gene_names = colnames(abs_parallel_res_mat1))
# abs_parallel_res_mat2 <- prepare_ccm_res(abs_parallel_res_lineage2, gene_names = colnames(abs_parallel_res_mat1))
#
# #create a graph and then plot the graph as a hiearchical fashion:
# abs_parallel_res_mat1_g <- abs_parallel_res_mat1
# abs_parallel_res_mat1_g[which(abs(abs_parallel_res_mat1_g) > 0.1, arr.ind = T)] <- abs_parallel_res_mat1_g[which(abs(abs_parallel_res_mat1_g) > 0.1, arr.ind = T)]
# abs_parallel_res_mat1_g[which(abs(abs_parallel_res_mat1_g) < 0.1, arr.ind = T)] <- 0
#
# #netbiov, igraph,
# #build a graph:
# layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1]
# g1 <- graph_from_adjacency_matrix(as.matrix(abs_parallel_res_mat1_g), mode = "directed", weighted = T)
#
# plot(g1, layout = layout.reingold.tilford(g1, root = c('Irf8','Gfi1')), vertex.size = 6)
# # plot(g1, layout = layout_nicely(g1, ))
#
# g5 <- graph.adjacency(m5$adjacency, mode="undirected")
# plot(g5, layout = layout_as_tree(g5) )
# tkplot(g1)
#
# # create data frames for vertices and edges with the right variable names
# MapperNodes <- mapperVertices(mapper_res, 1:length(mapper_res$points_in_vertex) )
# dev.off()
#
# MapperLinks <- mapperEdges(mapper_res)
#
# # interactive plot
# forceNetwork(Nodes = MapperNodes, Links = MapperLinks,
#              Source = "Linksource", Target = "Linktarget",
#              Value = "Linkvalue", NodeID = "Nodename",
#              Group = "Nodegroup", opacity = 0.8,
#              linkDistance = 10, charge = -400)
#
# mst.plot.mod(g1, v.size=1.5,e.size=.25,
#              colors=c("red",   "orange",   "yellow",   "green"),
#              mst.e.size=1.2,expression=abs(runif(vcount(g1),
#                                                  max=5,   min=1)),   sf=1,   v.sf=1,
#              mst.edge.col="white",   layout.function=layout.fruchterman.reingold)
#
#
# plot.abstract.nodes(g1,
#                     lab.cex=1, lab.color="white", v.sf=0, e.sf = 0,
#                     layout.function=layout.fruchterman.reingold)

# level.plot(g1, layout.function=NULL, type=1, initial_nodes=c('Gfi1', 'Irf8'),
#            init_nodes=0, order_degree = "in", plotsteps = FALSE,
#            saveplots=FALSE, dirname=NULL, vertex.colors=NULL,
#            edge.col=NULL, tkplot=F, nodeset=NULL,path.col="green",
#            col.s1="red", col.s2="yellow", nodes.on.path=TRUE,v.size=2,
#            e.size=.5, v.lab=T, bg="black", v.lab.cex=0.5,
#            v.lab.col="skyblue",sf=4,e.path.width=1,e.curve=.5,
#            level.spread=FALSE)

# pdf('./main_figures/fig4i.pdf')
# level.plot(g1, layout.function=NULL, type=1, initial_nodes=c(1, 2),
#            init_nodes=0, order_degree = "in", plotsteps = FALSE,
#            saveplots=FALSE, dirname=NULL, vertex.colors=NULL,
#            edge.col=NULL, tkplot=F, nodeset=NULL,path.col="green",
#            col.s1="red", col.s2="yellow", nodes.on.path=TRUE,v.size=2,
#            e.size=.5, v.lab=T, bg="black", v.lab.cex=0.5,
#            v.lab.col="skyblue",sf=4,e.path.width=1,e.curve=.5,
#            level.spread=FALSE)
# dev.off()
V(g1)$label <- V(g1)$name

pdf(paste(main_fig_dir, 'fig4i.pdf', sep = ''))
res <- level.plot(g1)
dev.off()

save(file = '/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/RData/res', res)
#save the data:
################################################################################################################################################################################################################################################
# save.image('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/RData/fig5.RData')
save.image('/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/RData/fig5.RData')
################################################################################################################################################################################################################################################

# #convert to jason: only restricted to tree structure not others
# library(data.tree)
# tree <- FromDataFrameNetwork(as.data.frame(as.matrix(hbp_network_mat)))

#create the tree with just igraph:
# V(g1)$label <- V(g1)$name
# V(g1)$name <- 1:length(V(g1)$name)
# res <- level.plot(g1 , initial_nodes=1:2, v.lab = T) # , initial_nodes=c('Gfi1', 'Irf8'
#
# # res <- level.plot(g1 , initial_nodes=c('Gfi1', 'Irf8')) # , initial_nodes=c('Gfi1', 'Irf8'
# # res <- level.plot(g1, v.lab = T, init_nodes = 2) # , initial_nodes=c('Gfi1', 'Irf8'
#
# res <- level.plot(g1) # , initial_nodes=c('Gfi1', 'Irf8'

# load("/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/res") #load the res (we need to mannually set the initial_nodes 1:2 in .process_graph1 function)
load("/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/RData/res") #load the res (we need to mannually set the initial_nodes 1:2 in .process_graph1 function)

cus_layout <- res$layout
master_regulator_id <- 1:2 #cus_layout[, 2] %in% min(unique(cus_layout[, 2])) #1:2 #
direct_target_id <- cus_layout[, 2] %in% unique(cus_layout[, 2])[order(unique(cus_layout[, 2])) == 2]
secondary_target_id <- cus_layout[, 2] %in% max(unique(cus_layout[, 2]))

cus_layout[master_regulator_id, 1]
cus_layout[master_regulator_id, 1] <- c(-10, 10)
cus_layout[direct_target_id, 1] <- cus_layout[direct_target_id, 1] * 1.5 #order(cus_layout[direct_target_id, 1]) * 2  - length(cus_layout[direct_target_id, 1]) * 20
cus_layout[secondary_target_id, 1] <- cus_layout[secondary_target_id, 1] * 1#order(cus_layout[secondary_target_id, 1]) * 10 - length(cus_layout[secondary_target_id, 1])  * 5

v_size <- rep(res$vertex.size, nrow(cus_layout))
v_size[master_regulator_id] <- 12
v_size[direct_target_id] <- 4#v_size[direct_target_id] * 1.5
v_size[secondary_target_id] <- 2.5

cus_layout[master_regulator_id, 2] <- 0
cus_layout[direct_target_id, 2] <- -1
cus_layout[secondary_target_id, 2] <- -2

# [1] "g"                  "layout"             "vertex.color"
# [4] "edge.color"         "edge.arrow.size"    "vertex.label.cex"
# [7] "vertex.label"       "lab.color"          "lab.dist"
# [10] "vertex.frame.color" "edge.width"         "bg"
# [13] "edge.curved"
res$vertex.label.cex <- rep(0.25, length(v_size))
res$vertex.label.cex[1:2] <- 1
res$layout <- cus_layout
res$bg <- 'white'
res$vertex.size <- v_size

res$vertex.label <- V(res$g)$name
res$lab.color <- 'black'

res$vertex.color[1:2] <- c('#BCA0CC')
res$vertex.color[direct_target_id] <- '#77CCD2'
res$vertex.color[res$vertex.color == 'orange2'] <- '#7EB044'

#label the second layer genes without outdegree
secondary_layer_degree <- degree(res$g, v = V(res$g)$name[cus_layout[, 2] == -1], mode = c("out"),
                                 loops = TRUE, normalized = FALSE)
res$vertex.color[which(cus_layout[, 2] == -1)[secondary_layer_degree == 0]] <- '#F3756C'
# secondary_layer <- degree(res$g, v = V(res$g)$name[cus_layout[, 2] == -1], mode = c("out"),
#                           loops = TRUE, normalized = FALSE)

res$vertex.frame.color <- res$vertex.color
res$vertex.label <- c(res$vertex.label[1:2], rep('', length(res$vertex.label) - 2))

pdf(paste(main_fig_dir, 'fig5i.pdf', sep = ''), width = 12, height = 7)
plot.netbiov(res)
dev.off()

# plot.modules(g1, layout.function = layout.reingold.tilford,
#              + col.grad=list(color.list$citynight), tkplot=FALSE)
# res <- level.plot(g1, tkplot=FALSE, level.spread=FALSE, initial_nodes = c(1, 2),
#            layout.function=layout.lgl) #, v.lab = T
#
# res <- level.plot(g1, tkplot=FALSE, level.spread=FALSE,
#                   layout.function=layout.fruchterman.reingold, v.lab = T) #
#
# V(g1)$label %in% checked_TFs_jaspar_fig1b[which(checked_TFs_jaspar_fig1b_category == 'other')]

# plot(res$g, layout = cus_layout, vertex.size = v_size, vertex.label.cex = res$vertex.label.cex,
#      vertex.color = res$vertex.color, edge.color = res$vertex.color, edge.arrow.size = res$edge.arrow.size,
#      lab.color = res$lab.color, lab.dist = res$lab.dist, vertex.frame.color = res$vertex.frame.color,
#      edge.width = res$edge.width, edge.curved = res$edge.curved)
# ################################################################################################################################################################################################################################################
# # check the branch time points for other regulators:
# ################################################################################################################################################################################################################################################
# # Hoxa9 Vdr Atf3
#
# # df[df$labels %in% c('Hoxa9', 'Vdr', 'Atf3', 'Gfi1', 'Irf8'), ]
# #
# # df[df$labels %in% c('Hoxa9', 'Vdr', 'Atf3'), 'type'] <- 'Master regulator'
# #
# # #add all other BEAM genes:
# #
# # pdf('main_figures//fig4g_new.pdf', width = 2.2, height = 1.5)
# # qplot(type, abs(bifurcation_time), color = type, geom = c('jitter'), data = subset(df_all, abs(bifurcation_time) > 60), alpha = I(0.7), log = 'y', size = 1) +  #, 'boxplot'
# #   xlab('') + ylab('bifurcation time point') + coord_flip() + scale_size(range = c(0.01, 1), limits = c(0.1, 1)) +
# #   geom_boxplot(aes(color = type), fatten = 0.5, lwd = 0.5, outlier.shape=NA, alpha = 0.5) +
# #   scale_y_continuous(breaks = round(seq(min(abs(df_all$bifurcation_time), na.rm = T), max(abs(df_all$bifurcation_time) + 5, na.rm = T), by = 9),1)) + #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`))
# #   nm_theme()
# # dev.off()
# #
#
#
# ################################################################################################################################################################################################################################################
# # show cluster distribution in states
# ################################################################################################################################################################################################################################################
# state_cluster_stat <- table(pData(URMM_all_fig1b)[, c('State', 'cluster')])
#
# state_cluster_stat <- apply(state_cluster_stat, 2, function(x) x / sum(x))
# state_cluster_stat_ordered <- t(state_cluster_stat)
#
# pdf(paste(main_fig_dir, "URMM_all_fig1b_state_branch_distribution.pdf", sep = ''), height = 6, width = 6)
# pheatmap::pheatmap(state_cluster_stat_ordered, cluster_cols = F, cluster_rows = F, color = colorRampPalette(RColorBrewer::brewer.pal(n=9, name='Greys'))(10))
# dev.off()
#
# state_cluster_stat <- table(pData(URMM_all_abs)[, c('State', 'paper_cluster')])
#
# state_cluster_stat <- apply(state_cluster_stat, 2, function(x) x / sum(x))
# state_cluster_stat_ordered <- t(state_cluster_stat)
#
# pdf(paste(main_fig_dir, "URMM_all_abs_state_branch_distribution.pdf", sep = ''), height = 6, width = 6)
# pheatmap::pheatmap(state_cluster_stat_ordered, cluster_cols = F, cluster_rows = F, color = colorRampPalette(RColorBrewer::brewer.pal(n=9, name='Greys'))(10))
# dev.off()
#
################################################################################################################################################################################################################################################
# BEAM test for the Ery/Meg cells
################################################################################################################################################################################################################################################
#run BEAM:
fig1b_beam_genes_ery_meg <- BEAM(URMM_all_fig1b, branch_point = 1, verbose = T, cores = detectCores() - 2)
URMM_all_abs_beam_genes_ery_meg <- BEAM(URMM_all_abs, branch_point = 1, verbose = T, cores = detectCores() - 2)

fig1b_beam_genes_proj_dup_ery_meg <- BEAM(URMM_all_fig1b, branch_point = 1, verbose = T, progenitor_method = 'duplicate', cores = detectCores() - 2)
URMM_all_abs_beam_genes_proj_dup_ery_meg <- BEAM(URMM_all_abs, branch_point = 1, verbose = T, progenitor_method = 'duplicate', cores = detectCores() - 2)

fig1b_beam_genes_ery_meg_ILRs <- calILRs(URMM_all_fig1b, branch_point = 1, verbose = T, cores = detectCores() - 2)
URMM_all_absbeam_genes_ery_meg_ILRs <- calILRs(URMM_all_abs, branch_point = 1, verbose = T, cores = detectCores() - 2)

fig1b_beam_genes_ery_meg_ILRs_all <- calILRs(URMM_all_fig1b, branch_point = 1, verbose = T, cores = detectCores() - 2, return_all = T)
URMM_all_absbeam_genes_ery_meg_ILRs_all <- calILRs(URMM_all_abs, branch_point = 1, verbose = T, cores = detectCores() - 2, return_all = T)

# run pseudotime test:
fig1b_beam_genes_pseudotime <- differentialGeneTest(URMM_all_fig1b, verbose = T, cores = detectCores() - 2)
URMM_all_abs_pseudotime <- differentialGeneTest(URMM_all_abs, verbose = T, cores = detectCores() - 2)

# run pseudotime test across two lineages: (IMPORTANT: make sure the cell state match the cell type assignment!!!)
URMM_all_fig1b_MEP <- URMM_all_fig1b[, pData(URMM_all_fig1b)$State %in% c(1, 5)]
URMM_all_fig1b_GMP <- URMM_all_fig1b[, pData(URMM_all_fig1b)$State %in% c(1:4)]
URMM_all_abs_MEP <- URMM_all_abs[, pData(URMM_all_abs)$State %in% c(1, 2)]
URMM_all_abs_GMP <- URMM_all_abs[, pData(URMM_all_abs)$State %in% c(1, 3:5)]

URMM_all_fig1b_progenitor_GM <- URMM_all_fig1b[, pData(URMM_all_fig1b)$State %in% c(1, 2)]
URMM_all_abs_progenitor_GM <- URMM_all_abs[, pData(URMM_all_abs)$State %in% c(1, 2)]

fig1b_pseudotime_MEP <- differentialGeneTest(URMM_all_fig1b_MEP, verbose = T, cores = detectCores() - 2)
fig1b_pseudotime_GMP <- differentialGeneTest(URMM_all_fig1b_GMP, verbose = T, cores = detectCores() - 2)
URMM_all_abs_pseudotime_MEP <- differentialGeneTest(URMM_all_abs_MEP, verbose = T, cores = detectCores() - 2)
URMM_all_abs_pseudotime_GMP <- differentialGeneTest(URMM_all_abs_GMP, verbose = T, cores = detectCores() - 2)

fig1b_pseudotime_progenitor_GM <- differentialGeneTest(URMM_all_fig1b_progenitor_GM, verbose = T, cores = detectCores() - 2)
URMM_all_abs_pseudotime_progenitor_GM <- differentialGeneTest(URMM_all_abs_progenitor_GM, verbose = T, cores = detectCores() - 2)

################################################################################################################################################################################################################################################
# Make lineage-score tree plot and related kinetic curves as well as multi-way heatmaps
################################################################################################################################################################################################################################################
# load MARS-seq dataset:
load('./valid_subset_GSE72857_cds2') #126 nodes in 10 dimensions
# score the erythroid / gmp scores:
valid_subset_GSE72857_cds2 <- orderCells(valid_subset_GSE72857_cds2, root_state = 10) # ensure pseudotime starts from State
# 1. identify the lineage for each branch gene
rowMeans(fig1b_beam_genes_ery_meg_ILRs)[row.names(subset(fig1b_beam_genes_ery_meg, qval <0.01))]
qplot(rowMeans(fig1b_beam_genes_ery_meg_ILRs)[row.names(subset(fig1b_beam_genes_ery_meg, qval <0.01))])

# 2. overlap positive/negative lineage score genes with get MARS-seq dataset:
ery_meg_lineage_score <- rowMeans(fig1b_beam_genes_ery_meg_ILRs)[row.names(subset(fig1b_beam_genes_ery_meg, qval <0.01))]
qplot(ery_meg_lineage_score)

qplot(ery_meg_lineage_score[intersect(row.names(valid_subset_GSE72857_cds2), names(ery_meg_lineage_score))])

positive_score_genes <- intersect(row.names(valid_subset_GSE72857_cds2), names(ery_meg_lineage_score[ery_meg_lineage_score > 0.5])) #positive genes (Ery/Meg lineage)
negtive_score_genes <- intersect(row.names(valid_subset_GSE72857_cds2), names(ery_meg_lineage_score[ery_meg_lineage_score < -0.5])) #negative genes (GMP lineage)

cell_ery_meg_lineage_score <- esApply(URMM_all_fig1b[c(positive_score_genes, negtive_score_genes), ], 2, function(x) mean(x[1:length(positive_score_genes)]) - mean(x[length(positive_score_genes):length(x)]))
pData(URMM_all_fig1b)$ery_meg_lineage_score <- cell_ery_meg_lineage_score

plot_complex_cell_trajectory(URMM_all_fig1b, color_by = 'ery_meg_lineage_score') +  scale_colour_gradient2()

pdf(paste(main_fig_dir, 'urmm_complex_tree_ery_meg_lineage_score.pdf', sep = ''), width = 1, height = 1.5)
plot_complex_cell_trajectory(URMM_all_fig1b, color_by = 'ery_meg_lineage_score', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001) +
  scale_colour_gradient2() + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
  theme(legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank())
dev.off()

pData(URMM_all_fig1b)$ery_meg_lineage_score <- - pData(URMM_all_fig1b)$ery_meg_lineage_score
pdf(paste(main_fig_dir, 'urmm_complex_tree_ery_meg_lineage_score2.pdf', sep = ''), width = 1, height = 1.5)
plot_complex_cell_trajectory(URMM_all_fig1b, color_by = 'ery_meg_lineage_score', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001) +
  scale_colour_gradient2() + scale_size(range = c(0.2, 0.2)) + scale_x_reverse() +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
  theme(legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank()) + theme_void() + theme (legend.position="none", legend.title=element_blank())
dev.off()

cell_ery_meg_lineage_score <- esApply(valid_subset_GSE72857_cds2[c(positive_score_genes, negtive_score_genes), ], 2, function(x) mean(x[1:length(positive_score_genes)]) - mean(x[length(positive_score_genes):length(x)]))

pData(valid_subset_GSE72857_cds2)$ery_meg_lineage_score <- cell_ery_meg_lineage_score
plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'ery_meg_lineage_score') +  scale_colour_gradient2()

pdf(paste(main_fig_dir, 'marseq_complex_tree_ery_meg_lineage_score.pdf', sep = ''), width = 1, height = 1.5)
plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'ery_meg_lineage_score', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001) +
  scale_colour_gradient2() + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
  theme(legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank())
dev.off()

pData(valid_subset_GSE72857_cds2)$ery_meg_lineage_score <- - pData(valid_subset_GSE72857_cds2)$ery_meg_lineage_score # again, accounting the switching of step 
pdf(paste(main_fig_dir, 'marseq_complex_tree_ery_meg_lineage_score2.pdf', sep = ''), width = 1, height = 1.5)
plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'ery_meg_lineage_score', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001, root_states = 10) +
  scale_colour_gradient2() + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
  theme(legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank()) + theme_void() + theme (legend.position="none", legend.title=element_blank())
dev.off()

pdf(paste(main_fig_dir, 'marseq_complex_tree_ery_meg_lineage_score2_helper.pdf', sep = ''))
plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'ery_meg_lineage_score', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001) +
  scale_colour_gradient2() + scale_size(range = c(0.2, 0.2))
dev.off()

pdf(paste(main_fig_dir, 'marseq_complex_tree_ery_meg_lineage_score_helper.pdf', sep = ''), width = 4, height = 1.5)
plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'ery_meg_lineage_score', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001) +
  scale_colour_gradient2() + scale_size(range = c(0.2, 0.2))
dev.off()

# make the heatmap and kinetic plot for selected genes:
sort(esApply(valid_subset_GSE72857_cds2[c(positive_score_genes, negtive_score_genes)], 1, mean))
c('Car1', 'Elane', 'Car2', 'Prtn3') %in% positive_score_genes #choose the top 4 genes

marseq_complex_tree_ery_meg_lineage_score_kinetic_curves_cols <- c("Neu" = "#6BBE44", "DC" = "#31C5F4", "Baso" = "#D493C0", "Mono" = "#F58F89", "MK" = "#4CBD8E", "Ery" = "#BCC132")
pdf(paste(main_fig_dir, 'marseq_complex_tree_ery_meg_lineage_score_kinetic_curves.pdf', sep = ''), width = 3, height = 1)
plot_multiple_branches_pseudotime(valid_subset_GSE72857_cds2[c('Car1', 'Elane', 'Car2', 'Prtn3'),],
                                  branches=c(1, 3, 4, 6, 11, 9),
                                  branches_name=c("Neu", "DC", "Baso", "Mono", "MK", "Ery"), color_by = 'Branch',
                                  nrow = 1, ncol = 4) + nm_theme() + scale_color_manual(values = marseq_complex_tree_ery_meg_lineage_score_kinetic_curves_cols)
dev.off()

pdf(paste(main_fig_dir, 'marseq_complex_tree_ery_meg_lineage_score_kinetic_curves_helper.pdf', sep = ''), width = 2, height = 2)
plot_multiple_branches_pseudotime(valid_subset_GSE72857_cds2[c('Car1', 'Elane', 'Car2', 'Prtn3'),],
                                  branches=c(1, 3, 4, 6, 11, 9),
                                  branches_name=c("Neu", "DC", "Baso", "Mono", "MK", "Ery"), color_by = 'Branch',
                                  nrow = 1, ncol = 4) + scale_color_manual(values = marseq_complex_tree_ery_meg_lineage_score_kinetic_curves_cols)
dev.off()

pdf(paste(main_fig_dir, 'marseq_complex_tree_ery_meg_lineage_score_multiway_heatmap.pdf', sep = ''), width = 10, height = 8)
plot_multiple_branches_heatmap(valid_subset_GSE72857_cds2[c(positive_score_genes, negtive_score_genes),],
                               branches=c(1, 3, 4, 6, 11, 9),
                               branches_name=c("Neu", "DC", "Baso", "Mono", "MK", "Ery"),
                               show_rownames=T,
                               num_clusters=4)
dev.off()
#B6A5D0 #76C143 #3DC6F2
URMM_complex_tree_ery_meg_lineage_score_kinetic_curves_cols <- c("Ery/Meg" = "#B6A5D0", "Monocyte" = "#76C143", "Granulocyte" = "#3DC6F2")
plot_spanning_tree(URMM_all_fig1b, color_by = 'Type', cell_size = 2) + scale_color_manual(values = cols, name = "Type")
pdf(paste(main_fig_dir, 'URMM_complex_tree_ery_meg_lineage_score_kinetic_curves.pdf', sep = ''), width = 2, height = 2)
plot_multiple_branches_pseudotime(URMM_all_fig1b[c('Car1', 'Elane', 'Car2', 'Prtn3'),],
                                  branches=c(3, 4, 5), color_by = 'Branch',
                                  branches_name=c("Granulocyte", "Monocyte", "Ery/Meg"), nrow = 2, ncol = 2) +
  scale_color_manual(values = URMM_complex_tree_ery_meg_lineage_score_kinetic_curves_cols) + nm_theme()
dev.off()

pdf(paste(main_fig_dir, 'URMM_complex_tree_ery_meg_lineage_score_kinetic_curves_helper.pdf', sep = ''), width = 2, height = 2)
plot_multiple_branches_pseudotime(URMM_all_fig1b[c('Car1', 'Elane', 'Car2', 'Prtn3'),],
                                  branches=c(3, 4, 5), color_by = 'Branch',
                                  branches_name=c("Granulocyte", "Monocyte", "Ery/Meg"), nrow = 2, ncol = 2) +
  scale_color_manual(values = URMM_complex_tree_ery_meg_lineage_score_kinetic_curves_cols)
dev.off()

# need to switch the color because of the ordering mismatching 
pdf(paste(main_fig_dir, 'URMM_complex_tree_ery_meg_lineage_score_multiway_heatmap.pdf', sep = ''), width = 10, height = 8) 
plot_multiple_branches_heatmap(URMM_all_fig1b[unique(c(positive_score_genes, negtive_score_genes)),],
                               branches=c(3, 4, 5),
                               branches_name=c("Granulocyte", "Monocyte", "Ery/Meg"),
                               show_rownames=T,
                               num_clusters=4)
dev.off()

################################################################################################################################################################################################
# show the stemness goes down in MARS-seq using genes from URMM dataset
# (only to the GMP branch not to the end of monocyte/granulocyte branch: Paul dataset only includes progenitor cells mostly):
################################################################################################################################################################################################
# identify list of pluripotent genes:
# > mep_genes
# [1] "Hba-a2" "Car2"   "Cited4" "Klf1"
# > cmp_genes
# [1] "Pf4"  "Apoe" "Flt3" "Cd74"
# > gmp_genes
# [1] "Mpo"   "Prg2"  "Prtn3" "Ctsg"
#
# #identify all pseudotime significant genes go down from the HSCP:
# # URMM_all_abs_pseudotime_progenitor_GM fig1b_pseudotime_progenitor_GM
# duplicated_cds <- buildBranchCellDataSet(URMM_all_fig1b, progenitor_method = 'duplicate', branch_point = 1)
# range(subset(pData(duplicated_cds), State == 1)[, 'Pseudotime']) #end of HSC
# range(subset(pData(duplicated_cds), State == 3)[, 'Pseudotime']) #end of GMP
#
# initial_reference_val_ery <- rowMeans(fig1b_beam_genes_ery_meg_ILRs_all$str_branchA_expression_curve_matrix[, 25:30])
# initial_reference_val_gmp <- rowMeans(fig1b_beam_genes_ery_meg_ILRs_all$str_branchB_expression_curve_matrix[, 25:30])
#
# diff_reference_ery <- fig1b_beam_genes_ery_meg_ILRs_all$str_branchA_expression_curve_matrix[, 31:100] - matrix(rep(initial_reference_val_ery, 70), ncol = 70)
# diff_reference_gmp <- fig1b_beam_genes_ery_meg_ILRs_all$str_branchB_expression_curve_matrix[, 31:55] - matrix(rep(initial_reference_val_gmp, 25), ncol = 25)
#
# all_down <- apply(diff_reference_ery, 1, function(x) all(x < 0)) & apply(diff_reference_gmp, 1, function(x) all(x < 0))
#
# pseudotime_degs <- intersect(row.names(subset(fig1b_pseudotime_MEP, qval < 0.01)), row.names(subset(fig1b_pseudotime_GMP, qval < 0.01)))
# all_down_valid <- Reduce(intersect, list(row.names(valid_subset_GSE72857_cds2), names(all_down[all_down]), pseudotime_degs))
#
# cell_stemness_score <- esApply(URMM_all_fig1b[all_down_valid, ], 2, function(x) mean(x))
# pData(URMM_all_fig1b)$cell_stemness_score <- cell_stemness_score
# plot_complex_cell_trajectory(URMM_all_fig1b, color_by = 'cell_stemness_score') +  scale_colour_gradientn(colours = terrain.colors(10))
#
# pdf(paste(main_fig_dir, 'urmm_complex_tree_stemness_score_pseudotime.pdf', sep = ''), width = 1, height = 1.5)
# plot_complex_cell_trajectory(URMM_all_fig1b, color_by = 'cell_stemness_score', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001) +
#   scale_colour_gradientn(colours = terrain.colors(10)) + scale_size(range = c(0.2, 0.2)) +
#   nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
#   theme(legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank())
# dev.off()
#
# pdf(paste(main_fig_dir, 'urmm_complex_tree_stemness_score_pseudotime2.pdf', sep = ''), width = 1, height = 1.5)
# plot_complex_cell_trajectory(URMM_all_fig1b, color_by = 'cell_stemness_score', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001) +
#   scale_colour_gradientn(colours = terrain.colors(10)) + scale_size(range = c(0.2, 0.2)) +
#   nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
#   theme(legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank())+ theme_void() + theme (legend.position="none", legend.title=element_blank())
# dev.off()
#
# #normalize UMI to TPM:
# tpm_exprs_valid_subset_GSE72857_cds2 <- esApply(valid_subset_GSE72857_cds2, 2, function(x) x / sum(x) * 1e6)
# cell_stemness_score <- apply(tpm_exprs_valid_subset_GSE72857_cds2[all_down_valid, ], 2, function(x) mean(x)) #look for relative expression changes
# # cell_stemness_score <- esApply(valid_subset_GSE72857_cds2[all_down_valid, ], 2, function(x) mean(x))
#
# pData(valid_subset_GSE72857_cds2)$cell_stemness_score <- log2(cell_stemness_score)
# pData(valid_subset_GSE72857_cds2)$no_expression <- cell_stemness_score == 0
# #plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'cell_stemness_score') +  scale_colour_gradient()
#
# plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'cell_stemness_score') +  scale_colour_gradientn(colours = terrain.colors(10)) + facet_wrap(~no_expression)
#
# pdf(paste(main_fig_dir, 'marseq_complex_tree_stemness_score_pseudotime.pdf', sep = ''), width = 2, height = 1.5) # #use this because those cells at the end of Paul data are still progenitors
# plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'cell_stemness_score', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001) +
#   scale_colour_gradientn(colours = terrain.colors(10)) + facet_wrap(~no_expression) + scale_size(range = c(0.2, 0.2)) +
#   nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
#   theme(legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank())
# dev.off()
#
# pdf(paste(main_fig_dir, 'marseq_complex_tree_stemness_score_pseudotime2.pdf', sep = ''), width = 2, height = 1.5) # #use this because those cells at the end of Paul data are still progenitors
# plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'cell_stemness_score', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001) +
#   scale_colour_gradientn(colours = terrain.colors(10)) + facet_wrap(~no_expression) + scale_size(range = c(0.2, 0.2)) +
#   nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
#   theme(legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank()) + theme_void() + theme (legend.position="none", legend.title=element_blank())
# dev.off()
#
# pdf(paste(main_fig_dir, 'marseq_complex_tree_stemness_score_pseudotime2_helper.pdf', sep = '')) # #use this because those cells at the end of Paul data are still progenitors
# plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'cell_stemness_score', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001) +
#   scale_colour_gradientn(colours = terrain.colors(10)) + facet_wrap(~no_expression) + scale_size(range = c(0.2, 0.2))
# dev.off()
#
# plot_multiple_branches_pseudotime(valid_subset_GSE72857_cds2[all_down_valid, ],
#                                   branches=c(1, 3, 4, 6, 11, 9), color_by = 'Branch',
#                                   branches_name=c("Neu", "DC", "Baso", "Mono", "MK", "Ery"), TPM = F,
#                                   nrow = 6, ncol = 6) + scale_y_log10()
#
# plot_multiple_branches_heatmap(valid_subset_GSE72857_cds2[all_down_valid,],
#                                branches=c(1, 3, 4, 6, 11, 9),
#                                branches_name=c("Neu", "DC", "Baso", "Mono", "MK", "Ery"),
#                                show_rownames=T,
#                                num_clusters=4)
#
# plot_multiple_branches_pseudotime(URMM_all_fig1b[c(all_down_valid),],
#                                   branches=c(2, 4, 5), color_by = 'Branch',
#                                   branches_name=c("Ery/Meg", "Monocyte", "Granulocyte"), nrow = 6, ncol = 6)
#
# plot_multiple_branches_heatmap(URMM_all_fig1b[unique(all_down_valid),],
#                                branches=c(2, 4, 5),
#                                branches_name=c("Ery/Meg", "Monocyte", "Granulocyte"),
#                                show_rownames=T,
#                                num_clusters=4)
#
################################################################################################################################################################################################
# show the stemness goes down in URMM dataset using genes calculated using all cells:
################################################################################################################################################################################################
duplicated_cds <- buildBranchCellDataSet(URMM_all_fig1b, progenitor_method = 'duplicate', branch_point = 1)

initial_reference_val_ery <- rowMeans(fig1b_beam_genes_ery_meg_ILRs_all$str_branchA_expression_curve_matrix[, 25:30])
initial_reference_val_gmp <- rowMeans(fig1b_beam_genes_ery_meg_ILRs_all$str_branchB_expression_curve_matrix[, 25:30])

diff_reference_ery <- fig1b_beam_genes_ery_meg_ILRs_all$str_branchA_expression_curve_matrix[, 31:100] - matrix(rep(initial_reference_val_ery, 70), ncol = 70)
diff_reference_gmp <- fig1b_beam_genes_ery_meg_ILRs_all$str_branchB_expression_curve_matrix[, 31:100] - matrix(rep(initial_reference_val_gmp, 70), ncol = 70)

all_down <- apply(diff_reference_ery, 1, function(x) all(x < 0)) & apply(diff_reference_gmp, 1, function(x) all(x < 0))

pseudotime_degs <- row.names(subset(fig1b_beam_genes_pseudotime, qval < 0.01))
all_down_valid <- Reduce(intersect, list(row.names(valid_subset_GSE72857_cds2), names(all_down[all_down]), pseudotime_degs))

cell_stemness_score <- esApply(URMM_all_fig1b[all_down_valid, ], 2, function(x) mean(x))
pData(URMM_all_fig1b)$cell_stemness_score <- cell_stemness_score
plot_complex_cell_trajectory(URMM_all_fig1b, color_by = 'cell_stemness_score') +  scale_colour_gradientn(colours = terrain.colors(10))

pdf(paste(main_fig_dir, 'urmm_complex_tree_stemness_score.pdf', sep = ''), width = 1, height = 1.5) #use this because we should look for genes has expression less than progenitor until terminal fates
plot_complex_cell_trajectory(URMM_all_fig1b, color_by = 'cell_stemness_score', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001) +
  scale_colour_gradientn(colours = terrain.colors(10)) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
  theme(legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank())
dev.off()

pdf(paste(main_fig_dir, 'urmm_complex_tree_stemness_score2.pdf', sep = ''), width = 1, height = 1.5) #use this because we should look for genes has expression less than progenitor until terminal fates
plot_complex_cell_trajectory(URMM_all_fig1b, color_by = 'cell_stemness_score', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001) +
  scale_colour_gradientn(colours = terrain.colors(10)) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) + scale_x_reverse() + 
  theme(legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank()) + theme_void() + theme (legend.position="none", legend.title=element_blank())
dev.off()

pdf(paste(main_fig_dir, 'urmm_complex_tree_stemness_score_helper.pdf', sep = ''), width = 4, height = 1.5)
plot_complex_cell_trajectory(URMM_all_fig1b, color_by = 'cell_stemness_score', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001) +
  scale_colour_gradientn(colours = terrain.colors(10)) + scale_size(range = c(0.2, 0.2))
dev.off()

#normalize UMI to TPM:
tpm_exprs_valid_subset_GSE72857_cds2 <- esApply(valid_subset_GSE72857_cds2, 2, function(x) x / sum(x) * 1e6)
cell_stemness_score <- apply(tpm_exprs_valid_subset_GSE72857_cds2[all_down_valid, ], 2, function(x) mean(x))
# cell_stemness_score <- esApply(valid_subset_GSE72857_cds2[all_down_valid, ], 2, function(x) mean(x))

pData(valid_subset_GSE72857_cds2)$cell_stemness_score <- log2(cell_stemness_score)
pData(valid_subset_GSE72857_cds2)$no_expression <- cell_stemness_score == 0
#plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'cell_stemness_score') +  scale_colour_gradient()

plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'cell_stemness_score') +  scale_colour_gradientn(colours = terrain.colors(10)) + facet_wrap(~no_expression)

pdf(paste(main_fig_dir, 'marseq_complex_tree_stemness_score.pdf', sep = ''), width = 2, height = 1.5)
plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'cell_stemness_score', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001) +
  scale_colour_gradientn(colours = terrain.colors(10)) + facet_wrap(~no_expression) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
  theme(legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank())
dev.off()

pdf(paste(main_fig_dir, 'marseq_complex_tree_stemness_score2.pdf', sep = ''), width = 2, height = 1.5)
plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'cell_stemness_score', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001, root_states = 10) +
  scale_colour_gradientn(colours = terrain.colors(10)) + facet_wrap(~no_expression) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
  theme(legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank()) + theme_void() + theme (legend.position="none", legend.title=element_blank())
dev.off()

plot_multiple_branches_heatmap(URMM_all_fig1b[unique(cole_markers),],
                               branches=c(3, 4, 5),
                               #branches_name=c("Ery/Meg", "Monocyte", "Granulocyte"),
                               show_rownames=T,
                               num_clusters=4)

plot_multiple_branches_heatmap(valid_subset_GSE72857_cds2[c(mep_genes, cmp_genes, gmp_genes),],
                               branches=c(1, 3, 4, 6, 11, 9),
                               branches_name=c("Neu", "DC", "Baso", "Mono", "MK", "Ery"),
                               show_rownames=T,
                               num_clusters=4)

plot_multiple_branches_pseudotime(valid_subset_GSE72857_cds2[c(mep_genes, cmp_genes, gmp_genes),],
                                  branches=c(1, 3, 4, 6, 11, 9), color_by = 'Branch',
                                  branches_name=c("Neu", "DC", "Baso", "Mono", "MK", "Ery"),
                                  nrow = 3, ncol = 4)

plot_multiple_branches_pseudotime(URMM_all_fig1b[c(mep_genes, cmp_genes, gmp_genes),],
                                  branches=c(3, 4, 5), color_by = 'Branch',
                                  branches_name=c("Ery/Meg", "Monocyte", "Granulocyte"), nrow = 3, ncol = 4)

plot_multiple_branches_pseudotime(URMM_all_fig1b[c(mep_genes, cmp_genes, gmp_genes),],
                                  branches=c(3, 4, 5), color_by = 'Branch',
                                  branches_name=c("Ery/Meg", "Monocyte", "Granulocyte"), nrow = 3, ncol = 4)

plot_multiple_branches_pseudotime(URMM_all_fig1b[c(cole_markers),],
                                  branches=c(3, 4, 5), color_by = 'Branch',
                                  branches_name=c("Ery/Meg", "Monocyte", "Granulocyte"), nrow = 5, ncol = 6)

intersect(c(cole_markers), row.names(valid_subset_GSE72857_cds2))
plot_multiple_branches_pseudotime(valid_subset_GSE72857_cds2[intersect(c(cole_markers), row.names(valid_subset_GSE72857_cds2)), ],
                                  branches=c(1, 3, 4, 6, 11, 9), color_by = 'Branch',
                                  branches_name=c("Neu", "DC", "Baso", "Mono", "MK", "Ery"),
                                  nrow = 4, ncol = 4)

HSC_genes <- c("Dusp2", "Egr1", "Fos", "Zfp36", "Dusp1", "Jun", "Btg2", "Fosb", "Cpa3", "Csrp3", "Pbx1", "Itga2b", "F2r", "Treml1", "F2rl2", "Sdpr", "Pf4", "Rab27b", "Mmrn1", "Mecom", "Gucy1a3", "Gpc6", "Tmem181c-ps", "Gimap6", "Ctla2b", "Gcnt2", "Calcrl", "Gimap1", "Tgtp2", "B2m", "Zfp608", "Pde4b", "Ikzf2", "Cbfa2t3", "Tmem176b", "Abcg1", "Serpina3g", "Myct1", "Mpl", "Tsc22d1", "Msi2", "Hlf", "Rgs1", "Angpt1", "Gpr56", "Car2", "Ctla2a", "Meis1", "Lat", "Nrgn", "Fyb", "Dpysl2", "Ifitm1", "Samsn1", "Spns2", "Ifitm3", "Pdzk1ip1", "Pitpnc1", "Rbp1", "Arglu1", "Rbm26", "Zfp36l2", "Gnb1", "Gata2", "Vamp5", "Slc22a3", "Vezf1", "Muc13", "Gnb4", "Ptrf", "Rtp4", "Sox4", "Cd34", "Eltd1", "Pan3", "Tspan13", "Flt3", "Mn1", "Dntt", "Wfdc17", "Il1r1", "Satb1", "Egfl7", "Tcf7l2", "Notch2", "Kdm6b", "Gvin1", "Gm4070", "Gm17757", "Gm1966", "Pigr", "G630071F17Rik")
Multi_lineage <- c("Gpc6", "Tmem181c-ps", "Gimap6", "Ctla2b", "Gcnt2", "Calcrl", "Gimap1", "Tgtp2", "B2m", "Zfp608", "Pde4b", "Ikzf2", "Cbfa2t3", "Tmem176b", "Abcg1", "Serpina3g", "Myct1", "Mpl", "Tsc22d1", "Msi2", "Hlf", "Rgs1", "Angpt1", "Gpr56", "Car2", "Ctla2a", "Meis1", "Lat", "Nrgn", "Fyb", "Dpysl2", "Ifitm1", "Samsn1", "Spns2", "Ifitm3", "Pdzk1ip1", "Pitpnc1", "Rbp1", "Arglu1", "Rbm26", "Zfp36l2", "Gnb1", "Gata2", "Vamp5", "Slc22a3", "Vezf1", "Muc13", "Gnb4", "Ptrf", "Rtp4", "Sox4", "Cd34", "Eltd1", "Pan3", "Tspan13", "Flt3", "Mn1", "Dntt", "Wfdc17", "Il1r1", "Satb1", "Egfl7", "Tcf7l2", "Notch2", "Kdm6b", "Gvin1", "Gm4070", "Gm17757", "Gm1966", "Pigr", "G630071F17Rik")
genes_goes_down_from_HSC <- intersect(row.names(subset(fig1b_beam_genes_pseudotime, qval < 0.01)), Multi_lineage)

fig1d <- read.csv("./Nature_hta_paper_Monocle_run_archive/input/fig1d.txt",row.names=1, sep = '\t')
Multi_lineage <- row.names(fig1d)[fig1d$row_clusters.flat == 1][-1]

plot_multiple_branches_pseudotime(valid_subset_GSE72857_cds2[all_down_valid, ],
                                  branches=c(1, 3, 4, 6, 11, 9), color_by = 'Branch',
                                  branches_name=c("Neu", "DC", "Baso", "Mono", "MK", "Ery"), TPM = T,
                                  nrow = 6, ncol = 6)

plot_multiple_branches_pseudotime(URMM_all_fig1b[c(all_down_valid),],
                                  branches=c(3, 4, 5), color_by = 'Branch',
                                  branches_name=c("Ery/Meg", "Monocyte", "Granulocyte"), nrow = 6, ncol = 6)

plot_multiple_branches_heatmap(URMM_all_fig1b[unique(all_down_valid),],
                               branches=c(3, 4, 5),
                               branches_name=c("Ery/Meg", "Monocyte", "Granulocyte"),
                               show_rownames=T,
                               num_clusters=4)

plot_multiple_branches_heatmap(valid_subset_GSE72857_cds2[all_down_valid,],
                               branches=c(1, 3, 4, 6, 11, 9),
                               branches_name=c("Neu", "DC", "Baso", "Mono", "MK", "Ery"),
                               show_rownames=T,
                               num_clusters=4)

################################################################################################################################################################################################
# save the dataset 
################################################################################################################################################################################################

save.image('./RData/fig5.RData')

