# #placeholder for fig SI5 

#### currently it is just used for debugging...  

# 
# rm(list = ls())
# # 1. The legend for panel A needs to use complete words for the various cell types.
# # 2. Panel C should become panel B. Exclude the knockouts and the extra gates from this panel if they’re not already.
# # 3. Panel E (or whichever one is the wildtype BEAM heatmap) should become panel C.
# # 4. We need selected GO categories for the clusterers in the BEAM heatmap.
# # 5. The fonts on the BEAM heatmap are unreadably small. Use the same legends and axis annotations as we have in the Census AI files. Make the heatmap pretty and clear :)
# # 6. The Knockout panels should be labeled with standard notation (Gfi1-/-, Irf8-/-, and  Gfi1-/-/Irf8-/- I think)
# # 7. Move the transient gate panels to the supplement.
# # 8. Panel G is a bit confusing, because if I rememeber correctly, the Irf8 genes and Gfi1 genes are from the diff that’s *conditioned* on their pseudotimes. That’s not obvious from the figure and going to confuse the reader. What are we trying to communicate with this panel? I think it might be better to make a different UpSetR panel: one that shows the differentially between the Irf8-/-, Gfi-/-, and double KO cells and the WT cells collected at the corresponding gates. That is, presumably they used just one of their gates (maybe LKCD34+?) to collect the KO cells. So we should compare to the WT cells in that gate. Then we can show the intersection between those genes and the genes with ChIP peaks at their promoters to give a sense for how many direct targets are actually affected by loss of the regulator. Then we can show the intersection with the BEAM genes to show that BEAM is really picking up a big fraction of those genes, along with some others.
# # 9. We need to make it clear what the overlap is between genes that have Gfi1 ChIP peaks and/or Irf8 ChIP peaks and BEAM genes. Maybe add this to panel G?
# # 10. Split panel H in to two panels: one for Gfi1 and one for Irf8.
# # 11. Panel J is not clear and right now doesn’t add anything. We should drop this or reformat so it says something useful.
# 
# get_correct_root_state <- function(cds){
#   T0_counts <- table(pData(cds)$State, pData(cds)$Type)[,"Lsk"]
#   as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))])
# }
# ####################################################################################################################################################################################
# #load all package
# ####################################################################################################################################################################################
# library(dplyr)
# library(RColorBrewer)
# library(stringr)
# library(UpSetR)
# library(devtools)
# library(rEDM)
# library(xlsx)
# library(xacHelper)
# library(reshape2)
# library(pheatmap)
# library(d3Network)
# library(netbiov)
# library(monocle)
# library(plyr)
# library(dplyr)
# library(venneuler)
# library(SLICER)
# library(destiny)
# library(grid)
# library(HSMMSingleCell)
# 
# source('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/scripts/function.R')
# 
# # fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/Figures/"
# # tmp_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/tmp_figures/"
# fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/Figures/"
# tmp_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/Figures/tmp_figures/"
# 
# #figure 4 panels:
# # main_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/Figures/main_figures/"
# # SI_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/Figures/supplementary_figures/"
# main_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/Figures/main_figures/"
# SI_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/Figures/supplementary_figures/"
# 
# #########################################################################################################################################################################
# # prepare the cds for both of the WT dataset and the full dataset
# #########################################################################################################################################################################
# 
# #reading the exprs data and create a cell dataset:
# # source('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/scripts/nature_hta_analysis_tmp.R')
# hta_exprs <- read.csv("./Nature_hta_paper_Monocle_run_archive/input/Olsson_RSEM_SingleCellRNASeq.csv",row.names=1)
# sample_sheet <- data.frame(groups = str_split_fixed(colnames(hta_exprs), "\\.+", 3), row.names = colnames(hta_exprs))
# gene_ann <- data.frame(gene_short_name = row.names(hta_exprs), row.names = row.names(hta_exprs))
# pd <- new("AnnotatedDataFrame",data=sample_sheet)
# fd <- new("AnnotatedDataFrame",data=gene_ann)
# 
# tpm_mat <- hta_exprs
# tpm_mat <- apply(tpm_mat, 2, function(x) x / sum(x) * 1e6)
# 
# URMM_all_std <- newCellDataSet(as.matrix(tpm_mat),phenoData = pd,featureData =fd,
#                                expressionFamily = negbinomial.size(),
#                                lowerDetectionLimit=1)
# URMM_all_std <- estimateSizeFactors(URMM_all_std)
# URMM_all_std <- estimateDispersions(URMM_all_std)
# 
# #set up the experimental type for each cell
# pData(URMM_all_std)[, 'Type'] <- as.character(pData(URMM_all_std)[, 'groups.1']) #WT cells
# pData(URMM_all_std)[453:593, 'Type'] <- paste(as.character(pData(URMM_all_std)[453:593, 'groups.1']), '_knockout', sep = '') #KO cells
# pData(URMM_all_std)[594:640, 'Type'] <- paste(pData(URMM_all_std)[594:640, 'groups.1'], pData(URMM_all_std)[594:640, 'groups.2'], 'knockout', sep = '_') #double KO cells
# 
# #run Census to get the transcript counts
# URMM_all_abs_list <- relative2abs(URMM_all_std, t_estimate = estimate_t(URMM_all_std), return_all = T)
# URMM_all_abs <- newCellDataSet(as.matrix(URMM_all_abs_list$norm_cds),
#                                phenoData = new("AnnotatedDataFrame",data=pData(URMM_all_std)),
#                                featureData = new("AnnotatedDataFrame",data=fData(URMM_all_std)),
#                                expressionFamily = negbinomial.size(),
#                                lowerDetectionLimit=1)
# URMM_all_abs <- estimateSizeFactors(URMM_all_abs)
# URMM_all_abs <- estimateDispersions(URMM_all_abs)
# 
# URMM_all_abs <- setOrderingFilter(URMM_all_abs, row.names(fData(URMM_all_abs)))
# 
# #########################################################################################################################################################################
# #read data from figure 1b
# fig1b <- read.csv("./Nature_hta_paper_Monocle_run_archive/input/fig1b.txt",row.names=1, sep = '\t')
# 
# #match up the column name in fig1b to the colnames in URMM_all_fig1b
# #note that you should not run this mutliple times
# URMM_all_fig1b <- URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Lsk', 'Cmp', 'Gmp', 'LK')]
# 
# fig1b_names <- colnames(fig1b)
# match_id <- which(str_split_fixed(colnames(fig1b), "\\.+", 2)[, 2] %in% colnames(URMM_all_fig1b) == T)
# fig1b_names[match_id] <- str_split_fixed(colnames(fig1b), "\\.+", 2)[match_id, 2]
# no_match_id <- which(str_split_fixed(colnames(fig1b), "\\.+", 2)[, 2] %in% colnames(URMM_all_fig1b) == F)
# fig1b_names[no_match_id] <- str_split_fixed(colnames(fig1b), "\\.\\.", 2)[no_match_id, 2]
# colnames(fig1b)[2:383] <- fig1b_names[2:383]
# 
# #learn trajectory using genes from figure 1b:
# URMM_all_fig1b <- order_cell_tree(row.names(fig1b)[2:533], cds = URMM_all_fig1b,
#                                   initial_method = PCA, norm_method = 'vstExprs', order_by = NULL)
# 
# #set up the color for each experiment
# cols <- c("Lsk" = "#edf8fb", "Cmp" = "#ccece6", "Gmp" = "#99d8c9", "GG1" = "#66c2a4", "IG2" = "#41ae76", "Irf8" = "#238b45", "LK" = "#005824",
#           "Irf8_knockout" = "#fc8d59", "Gfi1_Irf8_knockout" = "#636363", "Gfi1_knockout" = "#dd1c77")
# plot_spanning_tree(URMM_all_fig1b, color_by = 'Type', cell_size = 2) + scale_color_manual(values = cols, name = "Type") #+ nm_theme()
# 
# #########################################################################################################################################################################
# #assign clusters to each cell based on the clustering in the original study
# pData(URMM_all_fig1b)$cluster <- 0
# cluster_assignments <- as.numeric(fig1b[1, 2:383])
# cluster_assignments <- revalue(as.factor(cluster_assignments), c("1" = "HSCP-1", "2" = "HSCP-2", "3" = "Meg", "4" = "Eryth",
#                                                                  "5" = "Multi-Lin", "6" = "MDP", "7" = "Mono", "8" = "Gran", "9" = "Myelocyte"))
# 
# names(cluster_assignments) <- colnames(fig1b[1, 2:383])
# pData(URMM_all_fig1b)$cluster <- cluster_assignments[row.names(pData(URMM_all_fig1b))]
# #########################################################################################################################################################################
# URMM_all_fig1b <- recreate_cds(URMM_all_fig1b)
# 
# num_cells_expressed <- round(ncol(URMM_all_fig1b) * 0.05)
# #1. determine how many pca dimension you want:
# URMM_all_fig1b <- detectGenes(URMM_all_fig1b)
# fData(URMM_all_fig1b)$use_for_ordering[fData(URMM_all_fig1b)$num_cells_expressed > num_cells_expressed] <- T
# URMM_pc_variance <- plot_pc_variance_explained(URMM_all_fig1b, return_all = T)
# 
# #2. run reduceDimension with tSNE as the reduction_method
# URMM_all_fig1b <- reduceDimension(URMM_all_fig1b, max_components=2, norm_method = 'log', reduction_method = 'tSNE', num_dim = 5,  verbose = T)
# 
# #3. initial run of clusterCells_Density_Peak
# URMM_all_fig1b <- clusterCells_Density_Peak(URMM_all_fig1b, verbose = T, num_clusters = 3)
# 
# #4. check the clusters
# plot_cell_clusters(URMM_all_fig1b, color_by = 'as.factor(Cluster)', show_density = F, show_density_peak = T)
# 
# plot_cell_clusters(URMM_all_fig1b, color_by = 'cluster', show_density = F)
# plot_cell_clusters(URMM_all_fig1b, color_by = 'Type', show_density = F)
# # plot_cell_clusters(URMM_all_fig1b, color_by = 'State', show_density = F) + scale_color_brewer(palette = 'Set1') + facet_wrap(~paper_cluster)
# # plot_cell_clusters(URMM_all_fig1b, color_by = 'Cluster', show_density = F) + scale_color_brewer(palette = 'Set1') + facet_wrap(~paper_cluster)
# 
# #5. also check the decision plot
# plot_rho_delta(URMM_all_fig1b, rho_threshold = 6, delta_threshold = 16)
# 
# #6. re-run cluster and skipping calculating the rho_sigma
# URMM_all_fig1b <- clusterCells_Density_Peak(URMM_all_fig1b, verbose = T,  rho_threshold = 6, delta_threshold = 15, skip_rho_sigma = T)
# 
# #7. make the final clustering plot:
# plot_cell_clusters(URMM_all_fig1b, color_by = 'as.factor(Cluster)', show_density = F)
# plot_cell_clusters(URMM_all_fig1b, color_by = 'as.factor(Type)', show_density = F)
# plot_cell_clusters(URMM_all_fig1b, color_by = 'as.factor(cluster)', show_density = F)
# 
# #perform DEG test across clusters:
# URMM_all_fig1b@expressionFamily <- negbinomial.size()
# pData(URMM_all_fig1b)$Cluster <- factor(pData(URMM_all_fig1b)$Cluster)
# URMM_clustering_DEG_genes <- differentialGeneTest(URMM_all_fig1b, fullModelFormulaStr = '~Cluster', cores = detectCores() - 2)
# # URMM_clustering_DEG_genes_subset <- URMM_clustering_DEG_genes[fData(URMM_all_fig1b)$num_cells_expressed > num_cells_expressed, ]
# 
# #use all DEG gene from the clusters
# URMM_ordering_genes <- row.names(URMM_clustering_DEG_genes)[order(URMM_clustering_DEG_genes$qval)][1:1000]
# #
# #URMM_ordering_genes <- row.names(URMM_clustering_DEG_genes_subset)[order(URMM_clustering_DEG_genes_subset$qval)][1:1000]
# URMM_ordering_genes
# 
# URMM_all_fig1b <- setOrderingFilter(URMM_all_fig1b, ordering_genes = URMM_ordering_genes)
# URMM_all_fig1b <- reduceDimension(URMM_all_fig1b, verbose = T, scaling = T, max_components = 2) #, maxIter = 100, initial_method = DM, R_tSNE, tSNE, destiny_diffusionMaps, maxIter = 100 , param.gamma = 100, ncenter = 100
# URMM_all_fig1b <- orderCells(URMM_all_fig1b)
# plot_cell_trajectory(URMM_all_fig1b, color_by = 'Type')
# plot_cell_trajectory(URMM_all_fig1b, color_by = 'Type') + facet_wrap(~cluster)
# 
# pData(URMM_all_fig1b)$paper_cluster <- pData(URMM_all_fig1b)[colnames(URMM_all_fig1b), 'cluster']
# 
# pdf(paste(SI_fig_dir, "URMM_cluster_DEG_qval_0.01.pdf", sep = ''), height = 2, width = 2)
# plot_cell_trajectory(URMM_all_fig1b, color_by = 'paper_cluster', show_branch_points = F) + monocle:::monocle_theme_opts() + nm_theme()
# dev.off()
# 
# pdf(paste(SI_fig_dir, "URMM_cluster_DEG_qval_0.01_helper.pdf", sep = ''))
# plot_cell_trajectory(URMM_all_fig1b, color_by = 'paper_cluster')
# dev.off()
# 
# ########################################################################################################################################################
# #run on the full dataset
# URMM_all_abs <- recreate_cds(URMM_all_abs)
# pData(URMM_all_abs)$paper_cluster <- NA
# pData(URMM_all_abs)[colnames(URMM_all_fig1b), 'paper_cluster'] <- as.character(pData(URMM_all_fig1b)[, 'paper_cluster'])
# 
# URMM_all_abs <- setOrderingFilter(URMM_all_abs, ordering_genes = URMM_ordering_genes)
# URMM_all_abs <- reduceDimension(URMM_all_abs, verbose = T, norm_method = 'log', scaling = T, maxIter = 100) # norm_method = 'log',, param.gamma = 100
# URMM_all_abs <- orderCells(URMM_all_abs)
# plot_cell_trajectory(URMM_all_abs, color_by = 'Type') + facet_wrap(~paper_cluster)
# 
# if(as.numeric(unique(pData(URMM_all_abs)$State)) > 3)
#   URMM_all_abs <- trimTree(URMM_all_abs)
# plot_cell_trajectory(URMM_all_abs, color_by = 'State') + facet_wrap(~Type)
# 
# pdf(paste(SI_fig_dir, "all_data_URMM_cluster_DEG_qval_0.01.pdf", sep = ''), height = 2, width = 2)
# plot_cell_trajectory(URMM_all_abs, color_by = 'Type') + monocle:::monocle_theme_opts() + nm_theme()
# dev.off()
# 
# pdf(paste(SI_fig_dir, "all_data_URMM_cluster_DEG_qval_0.01_helper.pdf", sep = ''))
# plot_cell_trajectory(URMM_all_abs, color_by = 'paper_cluster')
# dev.off()

########################################################################################################################################################
#add dpt, slicer, wishbone analysis for this dataset
########################################################################################################################################################
# URMM_all_fig1b

# URMM_all_fig1b <- setOrderingFilter(URMM_all_fig1b, URMM_ordering_genes)
# slicer_res <- run_slicer(URMM_all_fig1b[URMM_ordering_genes, ], start = which.min(pData(URMM_all_fig1b)$Pseudotime))

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
qplot(dm1, dm2, data = URMM_ordering_geness_wishbone_res_downsampling_empirical)
qplot(tSNE1, tSNE2, data = URMM_ordering_geness_wishbone_res_downsampling_empirical, color = as.character(branch), size = 0.5)

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

URMM_all_abs <- setOrderingFilter(URMM_all_fig1b, URMM_ordering_genes)
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
qplot(dm1, dm2, data = wishbone_res_downsampling_empirical)
qplot(tSNE1, tSNE2, data = wishbone_res_downsampling_empirical, color = as.character(branch), size = 0.5)

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
plot_cell_trajectory(URMM_all_fig1b[, ], color_by = 'cluster', show_branch_points = F, theta = 120, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~Type, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = cluster_cols, name = "cluster") +theme (legend.position="right", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) + theme (legend.position="top", legend.title=element_blank()) #theme(axis.text.x = element_text(angle = 30, hjust = 1))
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
#fig1c table
fig1c <- read.csv("./Nature_hta_paper_Monocle_run_archive/input/fig1c.txt",row.names=1, sep = '\t')
fig1c_lineage_genes <- row.names(fig1c)[which(row.names(fig1c) == 'Cebpa'):nrow(fig1c)]
pData(URMM_all_abs)$paper_cluster <- NA
pData(URMM_all_abs)[colnames(URMM_all_fig1b), 'paper_cluster'] <- as.character(pData(URMM_all_fig1b)[, 'paper_cluster'])

pData(URMM_all_abs)$paper_cluster <- factor(pData(URMM_all_abs)$paper_cluster, levels = c("HSCP-1", "HSCP-2", "Multi-Lin", "MDP", "Eryth", "Meg", "Mono", "Gran", "Myelocyte"))
pData(URMM_all_abs)$Type <- factor(pData(URMM_all_abs)$Type, levels = c("Lsk", 'Cmp','Gmp', 'LK', "GG1", 'IG2', 'Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout'))
pData(URMM_all_fig1b)$Type <- factor(pData(URMM_all_fig1b)$Type, levels = c('Lsk', 'Cmp', 'Gmp', 'LK'))

pdf(paste(main_fig_dir, 'hta_cluster_dispersion_genes_fig1b_data_branch_curves.pdf', sep = ''), height = 13, width = 5)
plot_genes_branched_pseudotime(URMM_all_fig1b[fig1c_lineage_genes, ], color_by = 'cluster', branch_point=1,  cell_size = 0.2, nrow = 3, ncol = 5) + scale_colour_brewer(palette = "Set1") +
  scale_size(range = c(0.2, 0.2)) #+ nm_theme()
dev.off()

pdf(paste(main_fig_dir, 'fig4b.pdf', sep = ''), height = 1.5, width = 4)
plot_cell_trajectory(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout')], color_by = 'Type', show_branch_points = F, cell_size = 0.5, cell_link_size = 0.3, theta = 90) + facet_wrap(~Type, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = type_cols, name = "Type")  + theme (legend.position="none", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) #theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pdf(paste(main_fig_dir, 'fig4c.pdf', sep = ''), height = 1.5, width = 2.667)
plot_cell_trajectory(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('GG1', 'IG2')], color_by = 'Type', show_branch_points = F, cell_size = 0.5, cell_link_size = 0.3, theta = 90) + facet_wrap(~Type) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_color_manual(values = type_cols, name = "paper_cluster") + theme (legend.position="none", legend.title=element_blank()) +
  stat_density2d(color='black', alpha=I(0.25), size=I(0.25))
#theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pdf(paste(main_fig_dir, 'fig4d.pdf', sep = ''), height = 2, width = 3)
plot_genes_branched_pseudotime(URMM_all_fig1b[important_genes[1:6], ], color_by = 'Type', cell_size = 0.2, nrow = 2, ncol = 3) + scale_colour_brewer(palette = "Set1")  +
  scale_size(range = c(0.2, 0.2)) + nm_theme()
dev.off()

pdf(paste(main_fig_dir, 'fig4d.1.pdf', sep = ''), height = 2, width = 3)
plot_genes_branched_pseudotime(URMM_all_fig1b[c('Gfi1', "Irf8", "Ly86", "Csf1r", "Cebpe", "S100a8"), ], color_by = 'Type', cell_size = 0.2, nrow = 2, ncol = 3) + scale_colour_brewer(palette = "Set1")  +
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
