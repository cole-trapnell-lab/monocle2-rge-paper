#1. lung, blood, (maybe mar-seq data) including:
#a. dp clustering results
#b. upset plot
#c. heatmap vs annotation plot
#
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
library(piano)

qval_thrsld <- 0.1
hclust_method <-"ward.D2"
################################################################################################################################################################################################################################################
## set up directories
################################################################################################################################################################################################################################################

main_fig_dir <- "./Figures/main_figures/"
SI_fig_dir <- "./Figures/supplementary_figures/"
source('./scripts/function.R', echo = T)
load('./RData/fig1.RData')

################################################################################################################################################
# lung data
################################################################################################################################################
load('./RData/prepare_lung_data.RData')

#1. determine how many pca dimension you want:
absolute_cds <- recreate_cds(absolute_cds)

absolute_cds <- detectGenes(absolute_cds)
fData(absolute_cds)$use_for_ordering <- F

num_cells_expressed <- round(0.05 * ncol(absolute_cds))
fData(absolute_cds)$use_for_ordering[fData(absolute_cds)$num_cells_expressed > num_cells_expressed] <- T

absolute_cds@auxClusteringData[["tSNE"]]$variance_explained <- NULL
lung_pc_variance <- plot_pc_variance_explained(absolute_cds, return_all = T)

#2. run reduceDimension with tSNE as the reduction_method
# absolute_cds <- setOrderingFilter(absolute_cds, quake_id)
absolute_cds <- reduceDimension(absolute_cds, max_components=2, norm_method = 'log', reduction_method = 'tSNE', num_dim = 5,  verbose = T)

#3. initial run of clusterCells_Density_Peak
absolute_cds <- clusterCells_Density_Peak(absolute_cds, verbose = T)

#4. check the clusters
plot_cell_clusters(absolute_cds, color_by = 'as.factor(Cluster)', show_density = F)
plot_cell_clusters(absolute_cds, color_by = 'Time', show_density = F)
plot_cell_clusters(absolute_cds, color_by = 'State', show_density = F)

#5. also check the decision plot
plot_rho_delta(absolute_cds, rho_threshold = 3, delta_threshold = 9)
plot_cell_clusters(absolute_cds, color_by = 'Time', show_density = F, rho_threshold = 3, delta_threshold = 9)

#6. re-run cluster and skipping calculating the rho_sigma
absolute_cds <- clusterCells_Density_Peak(absolute_cds, verbose = T,  rho_threshold = 3, delta_threshold = 9, skip_rho_sigma = T)

#7. make the final clustering plot:
plot_cell_clusters(absolute_cds, color_by = 'as.factor(Cluster)', show_density = F)
plot_cell_clusters(absolute_cds, color_by = 'as.factor(Time)', show_density = F)

pData(absolute_cds)$Pseudotime <- pData(AT12_cds_subset_all_gene)$Pseudotime
pData(absolute_cds)$State <- pData(AT12_cds_subset_all_gene)$State
plot_cell_clusters(absolute_cds, color_by = 'Pseudotime', show_density = F)
plot_cell_clusters(absolute_cds, color_by = 'State', show_density = F)

#perform DEG test across clusters:
absolute_cds@expressionFamily <- negbinomial.size()
pData(absolute_cds)$Cluster <- factor(pData(absolute_cds)$Cluster)
lung_clustering_DEG_genes <- differentialGeneTest(absolute_cds, fullModelFormulaStr = '~Cluster', cores = detectCores() - 2)
lung_clustering_DEG_genes_subset <- lung_clustering_DEG_genes[fData(absolute_cds)$num_cells_expressed > num_cells_expressed, ]

#use all DEG gene from the clusters
lung_ordering_genes <- row.names(subset(lung_clustering_DEG_genes, qval < qval_thrsld))

lung_ordering_genes <- row.names(lung_clustering_DEG_genes_subset)[order(lung_clustering_DEG_genes_subset$qval)][1:1000] #1971

absolute_cds <- setOrderingFilter(absolute_cds, ordering_genes = lung_ordering_genes)
absolute_cds <- reduceDimension(absolute_cds, norm_method = 'log', verbose = T)
absolute_cds <- orderCells(absolute_cds)

################################################################################################################################################
## feature selection by marker genes
########################################################################################################################################################

#use all marker gene for making the trajectory
absolute_cds <- setOrderingFilter(absolute_cds, ordering_genes = quake_id)
absolute_cds <- reduceDimension(absolute_cds, norm_method = 'log', verbose = T)

pdf(paste(SI_fig_dir, "lung_marker.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(absolute_cds, color_by = 'Time', show_branch_points = F) + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

pdf(paste(SI_fig_dir, "lung_myo_marker_time_dp.pdf", sep = ''), height = 2, width = 2)
plot_cell_clusters(absolute_cds, color_by = 'Time', show_branch_points = F) + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

################################################################################################################################################
## feature selection by PCA
# create SI2c (PCA)
########################################################################################################################################################
##Lung data
lung_pca_res <- selectOrderingGenes2(absolute_cds, dim = c(1, 2), top_gene_num = 500, method = c('PCA'), verbose = T)

#use all pca gene for making the trajectory
lung <- setOrderingFilter(absolute_cds, ordering_genes = lung_pca_res$gene_vec)
lung <- reduceDimension(lung, norm_method = 'log', verbose = T)

pdf(paste(SI_fig_dir, "lung_pca.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(lung, color_by = 'Time', show_branch_points = F)  + monocle:::monocle_theme_opts() + nm_theme()
dev.off()
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 #use the fstree top 500 genes for making the trajectory
# lung_fstree <- setOrderingFilter(absolute_cds, ordering_genes = lung_pca_res$gene_vec[order(weight, decreasing = T)][1:500])
# lung_fstree <- reduceDimension(lung_fstree, norm_method = 'log', verbose = T)

# pdf(paste(SI_fig_dir, "lung_fstree.pdf", sep = ''), height = 2, width = 2)
# plot_cell_trajectory(lung_fstree, color_by = 'Time') + monocle:::monocle_theme_opts() + nm_theme()
# dev.off()

################################################################################################################################################
##  feature selection by SLICER
# create SI2c (SLICER)
########################################################################################################################################################

lung_slicer_genes <- slicer_feature_genes(absolute_cds)
lung_mat <- log(exprs(absolute_cds[lung_slicer_genes, ]) + 1)
pheatmap(lung_mat, cluster_rows = T, cluster_cols = T, scale = 'row')

plot_pseudotime_heatmap(absolute_cds[lung_slicer_genes, ], use_gene_short_name = F)

#use all slicer gene for making the trajectory
lung_slicer <- setOrderingFilter(absolute_cds, ordering_genes = row.names(absolute_cds)[lung_slicer_genes])
lung_slicer <- reduceDimension(lung_slicer, norm_method = 'log', verbose = T)

pdf(paste(SI_fig_dir, "lung_slicer.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(lung_slicer, color_by = 'Time', show_branch_points = F)  + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

########################################################################################################################################################
##  feature selection by Over-dispersion
# create SI2c (Over-dispersion)
########################################################################################################################################################
disp_table <- dispersionTable(absolute_cds)
lung_overdispersion_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)
absolute_cds <- setOrderingFilter(absolute_cds, lung_overdispersion_genes$gene_id)
lung_over_dispersion <- setOrderingFilter(absolute_cds, lung_overdispersion_genes$gene_id)
lung_over_dispersion <- reduceDimension(lung_over_dispersion, verbose = T)

pdf(paste(SI_fig_dir, "lung_over_dispersion.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(lung_over_dispersion, color_by = 'Time', show_branch_points = F)  + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

########################################################################################################################################################
##  feature selection by differential gene expression test
# create SI2c (DEG)
########################################################################################################################################################
lung_time_dependent_genes <- differentialGeneTest(absolute_cds[,], fullModelFormulaStr="~Time", reducedModelFormulaStr="~1", cores=detectCores())
lung_deg_ordering_genes <- row.names(subset(lung_time_dependent_genes, qval < 1e-2))
lung_deg_ordering_genes <- row.names(head(lung_time_dependent_genes[order(lung_time_dependent_genes$qval),], n=1000))

################################################################################################################################################################################################################################################
# create SI2b
################################################################################################################################################################################################################################################
listInput <- list(Overdispersion = lung_overdispersion_genes$gene_id, Slicer = row.names(absolute_cds)[lung_slicer_genes],
                  PCA = lung_pca_res$gene_vec,'DP feature' = lung_ordering_genes, 'Time DEG' = lung_deg_ordering_genes)

pdf(paste(SI_fig_dir, "SI2b.pdf", sep = ''), height = 2.5, width = 5)
upset(fromList(listInput), nsets = 6, order.by = "freq")
dev.off()

################################################################################################################################################################################################################################################
## Panel D (heatmaps are not used in the paper)
################################################################################################################################################################################################################################################

#lung data:
lung_heatmap_matrix <- log(exprs(absolute_cds)[lung_ordering_genes, ] + 1)

#scale and trim the data:
scale_max=3
scale_min=-3

lung_heatmap_matrix=lung_heatmap_matrix[!apply(lung_heatmap_matrix, 1, sd)==0,]
lung_heatmap_matrix=Matrix::t(scale(Matrix::t(lung_heatmap_matrix),center=TRUE))
lung_heatmap_matrix=lung_heatmap_matrix[is.na(row.names(lung_heatmap_matrix)) == FALSE,]
lung_heatmap_matrix[is.nan(lung_heatmap_matrix)] = 0
lung_heatmap_matrix[lung_heatmap_matrix>scale_max] = scale_max
lung_heatmap_matrix[lung_heatmap_matrix<scale_min] = scale_min

lung_heatmap_matrix_ori <- lung_heatmap_matrix
lung_heatmap_matrix <- lung_heatmap_matrix[is.finite(lung_heatmap_matrix[, 1]), ] #remove the NA fitting failure genes for each branch

cols_dist <- as.dist((1 - cor(as.matrix(lung_heatmap_matrix)))/2)
cols_dist[is.na(cols_dist)] <- 1

ph <- pheatmap(lung_heatmap_matrix,
               useRaster = T,
               cluster_cols=T,
               cluster_rows=T,
               show_rownames=F,
               show_colnames=F,
               #scale="row",
               clustering_distance_cols=cols_dist,
               clustering_method = hclust_method,
               cutree_cols = 3,
               silent=F
               # filename=NA,
               # color=hmcols
               #color=hmcols#,
               # filename="expression_pseudotime_pheatmap.pdf",
)

annotation_col <- data.frame(State=factor(pData(absolute_cds[, ph$tree_col$labels])$State),
                             Time = factor(pData(absolute_cds[, ph$tree_col$labels])$Time),
                             'Density peak clusters'=factor(pData(absolute_cds[, ph$tree_col$labels])$Cluster),
                             row.names = ph$tree_col$labels)
ph_res <- pheatmap(lung_heatmap_matrix,
               useRaster = T,
               cluster_cols=T,
               cluster_rows=T,
               show_rownames=F,
               show_colnames=F,
               #scale="row",
               clustering_distance_cols=cols_dist,
               clustering_method = hclust_method,
               cutree_cols = 4,
               cutree_rows = 8,
               annotation_col=annotation_col,
               width = 2.5,
               height = 5,
               silent=F
               # filename=NA,
               # color=hmcols
               #color=hmcols#,
               # filename="expression_pseudotime_pheatmap.pdf",
)

pdf(paste(SI_fig_dir, "fig_si1d.pdf", sep = ''))
grid::grid.rect(gp = grid::gpar("fill", col = NA))
grid::grid.draw(ph_res$gtable)
dev.off()

#run the GO enrichment analysis for the following gene lists in the datasets:
# HSMM_myo_ordering_genes
# lung_ordering_genes
# URMM_ordering_genes
# URMM_all_ordering_genes

lung_clustering <- data.frame(Cluster=factor(cutree(ph_res$tree_row, 8)))
clusters <- as.numeric(lung_clustering$Cluster)
names(clusters) <- fData(absolute_cds[row.names(lung_clustering), ])$gene_short_name

lung_gsa_results_dp_genes <- collect_gsa_hyper_results(absolute_cds[, ], mouse_go_gsc, clusters)

pdf(paste(SI_fig_dir, "lung_gsa_results_dp_genes_go_enrichment.pdf", sep = ''), height=100, width=15)
plot_gsa_hyper_heatmap(absolute_cds, lung_gsa_results_dp_genes, significance = 1e-5)
dev.off()

dp_cluster_group <- rep(1, length(lung_ordering_genes))
names(dp_cluster_group) <- as.character(fData(absolute_cds)[lung_ordering_genes, 'gene_short_name'])
lung_cluster_hyper_geometric_results_go <- collect_gsa_hyper_results(absolute_cds, gsc = mouse_go_gsc, dp_cluster_group)
plot_gsa_hyper_heatmap(absolute_cds, lung_cluster_hyper_geometric_results_go, significance = 1e-50)

pdf(paste(SI_fig_dir, "fig_si1d_2.pdf", sep = ''), height = 5, width = 4)
plot_gsa_hyper_heatmap(absolute_cds, lung_cluster_hyper_geometric_results_go, significance = 1e-80) + monocle:::monocle_theme_opts() + nm_theme()
dev.off()
################################################################################################################################################
# create supplementary figures for the lung data
# create SI2a
################################################################################################################################################
lung_cols <- c("E14.5" = "#bdc131", "E16.5" = "#afa9d3", "E18.5" = "#35c5f4", "Adult" = "#d593c0")

pdf(paste(SI_fig_dir, "lung_cluster_DEG_qval_0.01.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(absolute_cds, color_by = 'Time', show_branch_points = F) + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

pdf(paste(SI_fig_dir, "lung_cluster_DEG_qval_0.01_helper.pdf", sep = ''))
plot_cell_trajectory(absolute_cds, color_by = 'Time')
dev.off()

pdf(paste(SI_fig_dir, "SI2a.pdf", sep = ''), height = 2, width = 2)
plot_cell_clusters(absolute_cds, color_by = 'Time', cell_size = 0.5) + geom_point(aes(shape = Cluster, color = factor(Time))) +
  monocle:::monocle_theme_opts() + scale_size(range = c(0.5, 0.5)) + 
  nm_theme() +
  scale_color_manual(values = lung_cols)
dev.off()

pdf(paste(SI_fig_dir, "SI2a_helper.pdf", sep = ''))
plot_cell_clusters(absolute_cds, color_by = 'Time', cell_size = 0.5) + geom_point(aes(shape = Cluster, color = factor(Time))) +
  monocle:::monocle_theme_opts() + scale_size(range = c(0.5, 0.5)) + 
  scale_color_manual(values = lung_cols)
dev.off()

################################################################################################################################################################################################################################################
## nature blood data
# create SI2d 
################################################################################################################################################################################################################################################
load('./RData/fig5.RData')

pdf(paste(SI_fig_dir, "URMM_fig_si1b.pdf", sep = ''), height = 2, width = 2)
plot_cell_clusters(URMM_all_fig1b, color_by = 'paper_cluster', cell_size = 0.5) + geom_point(aes(shape = Cluster, color = factor(paper_cluster))) + 
 monocle:::monocle_theme_opts() + scale_size(range = c(0.5, 0.5)) + nm_theme()
dev.off()

pdf(paste(SI_fig_dir, "URMM_fig_si1b_helper.pdf", sep = ''))
plot_cell_clusters(URMM_all_fig1b, color_by = 'paper_cluster', cell_size = 0.5) + geom_point(aes(shape = Cluster, color = factor(paper_cluster))) 
dev.off()

URMM_heatmap_matrix <- log(exprs(URMM_all_fig1b)[URMM_ordering_genes, ] + 1)

#scale and trim the data:
scale_max=3
scale_min=-3

URMM_heatmap_matrix=URMM_heatmap_matrix[!apply(URMM_heatmap_matrix, 1, sd)==0,]
URMM_heatmap_matrix=Matrix::t(scale(Matrix::t(URMM_heatmap_matrix),center=TRUE))
URMM_heatmap_matrix=URMM_heatmap_matrix[is.na(row.names(URMM_heatmap_matrix)) == FALSE,]
URMM_heatmap_matrix[is.nan(URMM_heatmap_matrix)] = 0
URMM_heatmap_matrix[URMM_heatmap_matrix>scale_max] = scale_max
URMM_heatmap_matrix[URMM_heatmap_matrix<scale_min] = scale_min

URMM_heatmap_matrix_ori <- URMM_heatmap_matrix
URMM_heatmap_matrix <- URMM_heatmap_matrix[is.finite(URMM_heatmap_matrix[, 1]), ] #remove the NA fitting failure genes for each branch

hclust_method <-"ward.D2"

cols_dist <- as.dist((1 - cor(as.matrix(URMM_heatmap_matrix)))/2)
cols_dist[is.na(cols_dist)] <- 1

ph <- pheatmap(URMM_heatmap_matrix,
               useRaster = T,
               cluster_cols=T,
               cluster_rows=T,
               show_rownames=F,
               show_colnames=F,
               #scale="row",
               clustering_distance_cols=cols_dist,
               clustering_method = hclust_method,
               cutree_cols = 4,
               silent=F
               # filename=NA,
               # color=hmcols
               #color=hmcols#,
               # filename="expression_pseudotime_pheatmap.pdf",
)

annotation_col <- data.frame(State=factor(pData(URMM_all_fig1b[, ph$tree_col$labels])$State),
                             'Cell Type' = factor(pData(URMM_all_fig1b[, ph$tree_col$labels])$paper_cluster),
                             'Density peak clusters'=factor(pData(URMM_all_fig1b[, ph$tree_col$labels])$Cluster),
                             row.names = ph$tree_col$labels)

ph_res <- pheatmap(URMM_heatmap_matrix,
                   useRaster = T,
                   cluster_cols=T,
                   cluster_rows=T,
                   show_rownames=F,
                   show_colnames=F,
                   #scale="row",
                   clustering_distance_cols=cols_dist,
                   clustering_method = hclust_method,
                   cutree_cols = 4,
                   cutree_rows = 8,
                   annotation_col=annotation_col,
                   width = 2.5,
                   height = 5,
                   silent=F
                   # filename=NA,
                   # color=hmcols
                   #color=hmcols#,
                   # filename="expression_pseudotime_pheatmap.pdf",
)
pdf(paste(SI_fig_dir, "URMM_fig_si1d.pdf", sep = ''))
grid::grid.rect(gp = grid::gpar("fill", col = NA))
grid::grid.draw(ph_res$gtable)
dev.off()

URMM_clustering <- data.frame(Cluster=factor(cutree(ph_res$tree_row, 8)))
clusters <- as.numeric(URMM_clustering$Cluster)
names(clusters) <- fData(URMM_all_fig1b[row.names(URMM_clustering), ])$gene_short_name

URMM_gsa_results_dp_genes <- collect_gsa_hyper_results(URMM_all_fig1b[, ], mouse_go_gsc, clusters)

pdf(paste(SI_fig_dir, "URMM_gsa_results_dp_genes_go_enrichment.pdf", sep = ''), height=10, width=15)
plot_gsa_hyper_heatmap(URMM_all_fig1b, URMM_gsa_results_dp_genes, significance = 1e-5)
dev.off()

dp_cluster_group <- rep(1, length(URMM_ordering_genes))
names(dp_cluster_group) <- as.character(fData(URMM_all_fig1b)[URMM_ordering_genes, 'gene_short_name'])
URMM_cluster_hyper_geometric_results_go <- collect_gsa_hyper_results(URMM_all_fig1b, gsc = mouse_go_gsc, dp_cluster_group)
plot_gsa_hyper_heatmap(URMM_all_fig1b, URMM_cluster_hyper_geometric_results_go, significance = 1e-5)

pdf(paste(main_fig_dir, "URMM_fig_si1d.1.pdf", sep = ''), height = 5, width = 4)
plot_gsa_hyper_heatmap(URMM_all_fig1b, URMM_cluster_hyper_geometric_results_go, significance = 1e-5) + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

dp_cluster_group <- rep(1, length(URMM_ordering_genes))
names(dp_cluster_group) <- as.character(fData(URMM_all_fig1b)[URMM_ordering_genes, 'gene_short_name'])
URMM_cluster_hyper_geometric_results_go <- collect_gsa_hyper_results(URMM_all_fig1b, gsc = mouse_go_gsc, dp_cluster_group)
plot_gsa_hyper_heatmap(URMM_all_fig1b, URMM_cluster_hyper_geometric_results_go, significance = 1e-5)

pdf(paste(SI_fig_dir, "URMM_fig_si1d_enrichment.pdf", sep = ''), height = 5, width = 4)
plot_gsa_hyper_heatmap(URMM_all_fig1b, URMM_cluster_hyper_geometric_results_go, significance = 1e-5) + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

################################################################################################################################################
## feature selection by PCA
# create SI2f (PCA)
########################################################################################################################################################
URMM_pca_res <- selectOrderingGenes2(URMM_all_fig1b, dim = c(1, 2), top_gene_num = 1000, method = c('PCA'), verbose = T)

#use all pca gene for making the trajectory
URMM <- setOrderingFilter(URMM_all_fig1b, ordering_genes = URMM_pca_res$gene_vec)
URMM <- reduceDimension(URMM, norm_method = 'log', verbose = T)
URMM <- orderCells(URMM)

pdf(paste(SI_fig_dir, "URMM_pca.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(URMM, color_by = 'paper_cluster', show_branch_points = F)  + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

################################################################################################################################################
##  feature selection by SLICER
# create SI2d (SLICER)
########################################################################################################################################################
URMM_slicer_genes <- slicer_feature_genes(URMM_all_fig1b)
URMM_mat <- log(exprs(URMM_all_fig1b[URMM_slicer_genes, ]) + 1)
pheatmap(URMM_mat, cluster_rows = T, cluster_cols = T, scale = 'row')

plot_pseudotime_heatmap(URMM_all_fig1b[URMM_slicer_genes, ], use_gene_short_name = F)

#use all slicer gene for making the trajectory
URMM_slicer <- setOrderingFilter(URMM_all_fig1b, ordering_genes = row.names(URMM_all_fig1b)[URMM_slicer_genes])
URMM_slicer <- reduceDimension(URMM_slicer, norm_method = 'log', verbose = T)
URMM_slicer <- orderCells(URMM_slicer)

pdf(paste(SI_fig_dir, "URMM_slicer.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(URMM_slicer, color_by = 'paper_cluster', show_branch_points = F)  + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

########################################################################################################################################################
##  feature selection by Over-dispersion
# create SI2d (Over-dispersion)
########################################################################################################################################################
#Over-dispersion genes:
disp_table <- dispersionTable(URMM_all_fig1b)
URMM_overdispersion_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)
URMM_dispersion <- setOrderingFilter(URMM_all_fig1b, URMM_overdispersion_genes$gene_id)
plot_ordering_genes(URMM_dispersion)

URMM_dispersion <- reduceDimension(URMM_dispersion, norm_method = 'log', verbose = T)
URMM_dispersion <- orderCells(URMM_dispersion)

pdf(paste(SI_fig_dir, "URMM_dispersion.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(URMM_dispersion, color_by = 'paper_cluster', show_branch_points = F)  + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

########################################################################################################################################################
##  feature selection by differential gene expression test
# create SI2d (DEG)
########################################################################################################################################################
URMM_time_dependent_genes <- differentialGeneTest(URMM_all_fig1b[,], fullModelFormulaStr="~paper_cluster", reducedModelFormulaStr="~1", cores=detectCores())
URMM_deg_ordering_genes <- row.names(subset(URMM_time_dependent_genes, qval < 1e-2))
URMM_deg_ordering_genes <- row.names(head(URMM_time_dependent_genes[order(URMM_time_dependent_genes$qval),], n=1000))

URMM_DEG <- setOrderingFilter(URMM_all_fig1b, URMM_deg_ordering_genes)
URMM_DEG <- reduceDimension(URMM_DEG, norm_method = 'log', verbose = T)
URMM_DEG <- orderCells(URMM_DEG)

pdf(paste(SI_fig_dir, "URMM_DEG.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(URMM_DEG, color_by = 'paper_cluster', show_branch_points = F)  + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

################################################################################################################################################################################################################################################
## Panel C
# create SI2d (Over-dispersion)
################################################################################################################################################################################################################################################
#use upset to show the overlapping for different feature selection methods: (pca, high-dispersion, time-dependent de, dp feature genes, slicer feature selection)
# example of list input (list of named vectors)
URMM_listInput <- list('High dispersion' = HSMM_overdispersion_genes$gene_id, Slicer = row.names(HSMM_myo)[HSMM_slicer_genes],
                  PCA = HSMM_pca_res$gene_vec, 'DP genes' = HSMM_myo_ordering_genes, 'Time DEG' = HSMM_deg_ordering_genes)

pdf(paste(SI_fig_dir, "SI2e.pdf", sep = ''), height = 2.5, width = 5)
upset(fromList(URMM_listInput), nsets = 6, order.by = "freq") + nm_theme()
dev.off()

################################################################################################################################################################################################################################################
## mar-seq data
# create SI2d (Over-dispersion)
################################################################################################################################################################################################################################################
source('./script_for_reproduce/dpFeature_for_all_datasets.R', echo = T)

# reupdate the directory inform again
main_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/Figures/main_figures/"
SI_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/Figures/supplementary_figures/"

pdf(paste(SI_fig_dir, "SI2g.pdf", sep = ''), height = 2, width = 2)
plot_cell_clusters(valid_subset_GSE72857_cds, color_by = 'cell_type', cell_size = 0.2) + geom_point(aes(shape = Cluster, color = cell_type)) +
  monocle:::monocle_theme_opts() + nm_theme() + scale_size(range = c(0.2, 0.2))
dev.off()

pdf(paste(SI_fig_dir, "SI2g_helper.pdf", sep = ''))
plot_cell_clusters(valid_subset_GSE72857_cds, color_by = 'cell_type', cell_size = 0.2) + geom_point(aes(shape = Cluster, color = cell_type)) 
dev.off()

MAP_heatmap_matrix <- log(exprs(valid_subset_GSE72857_cds)[MAP_ordering_genes, ] + 1)

#scale and trim the data:
scale_max=3
scale_min=-3

MAP_heatmap_matrix=MAP_heatmap_matrix[!apply(MAP_heatmap_matrix, 1, sd)==0,]
MAP_heatmap_matrix=Matrix::t(scale(Matrix::t(MAP_heatmap_matrix),center=TRUE))
MAP_heatmap_matrix=MAP_heatmap_matrix[is.na(row.names(MAP_heatmap_matrix)) == FALSE,]
MAP_heatmap_matrix[is.nan(MAP_heatmap_matrix)] = 0
MAP_heatmap_matrix[MAP_heatmap_matrix>scale_max] = scale_max
MAP_heatmap_matrix[MAP_heatmap_matrix<scale_min] = scale_min

MAP_heatmap_matrix_ori <- MAP_heatmap_matrix
MAP_heatmap_matrix <- MAP_heatmap_matrix[is.finite(MAP_heatmap_matrix[, 1]), ] #remove the NA fitting failure genes for each branch

cols_dist <- as.dist((1 - cor(as.matrix(MAP_heatmap_matrix)))/2)
cols_dist[is.na(cols_dist)] <- 1

hclust_method <-"ward.D2"

ph <- pheatmap(MAP_heatmap_matrix,
               useRaster = T,
               cluster_cols=T,
               cluster_rows=T,
               show_rownames=F,
               show_colnames=F,
               #scale="row",
               clustering_distance_cols=cols_dist,
               clustering_method = hclust_method,
               cutree_cols = 3,
               silent=F
               # filename=NA,
               # color=hmcols
               #color=hmcols#,
               # filename="expression_pseudotime_pheatmap.pdf",
)

annotation_col <- data.frame(State=factor(pData(valid_subset_GSE72857_cds[, ph$tree_col$labels])$State),
                             cell_type = factor(pData(valid_subset_GSE72857_cds[, ph$tree_col$labels])$cell_type),
                             'Density peak clusters'=factor(pData(valid_subset_GSE72857_cds[, ph$tree_col$labels])$Cluster),
                             row.names = ph$tree_col$labels)

ph_res <- pheatmap(MAP_heatmap_matrix,
               useRaster = T,
               cluster_cols=T,
               cluster_rows=T,
               show_rownames=F,
               show_colnames=F,
               #scale="row",
               clustering_distance_cols=cols_dist,
               clustering_method = hclust_method,
               cutree_cols = 4,
               cutree_rows = 8,
               annotation_col=annotation_col,
               width = 2.5,
               height = 5,
               silent=F
               # filename=NA,
               # color=hmcols
               #color=hmcols#,
               # filename="expression_pseudotime_pheatmap.pdf",
)
pdf(paste(SI_fig_dir, "MAP_fig_si1d.pdf", sep = ''))
grid::grid.rect(gp = grid::gpar("fill", col = NA))
grid::grid.draw(ph_res$gtable)
dev.off()

#run the GO enrichment analysis for the following gene lists in the datasets:
MAP_clustering <- data.frame(Cluster=factor(cutree(ph_res$tree_row, 8)))
clusters <- as.numeric(MAP_clustering$Cluster)
names(clusters) <- fData(valid_subset_GSE72857_cds[row.names(MAP_clustering), ])$gene_short_name

MAP_gsa_results_dp_genes <- collect_gsa_hyper_results(valid_subset_GSE72857_cds[, ], mouse_go_gsc, clusters)

pdf(paste(SI_fig_dir, "URMM_gsa_results_dp_genes_go_enrichment.pdf", sep = ''), height=100, width=15)
plot_gsa_hyper_heatmap(valid_subset_GSE72857_cds, MAP_gsa_results_dp_genes, significance = 1e-5)
dev.off()

dp_cluster_group <- rep(1, length(MAP_ordering_genes))
names(dp_cluster_group) <- as.character(fData(valid_subset_GSE72857_cds)[MAP_ordering_genes, 'gene_short_name'])
MAP_cluster_hyper_geometric_results_go <- collect_gsa_hyper_results(valid_subset_GSE72857_cds, gsc = mouse_go_gsc, dp_cluster_group)
plot_gsa_hyper_heatmap(valid_subset_GSE72857_cds, MAP_cluster_hyper_geometric_results_go, significance = 1e-5)

pdf(paste(main_fig_dir, "MAP_fig_si1d.1.pdf", sep = ''), height = 5, width = 4)
plot_gsa_hyper_heatmap(valid_subset_GSE72857_cds, MAP_cluster_hyper_geometric_results_go, significance = 1e-5) + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

dp_cluster_group <- rep(1, length(MAP_ordering_genes))
names(dp_cluster_group) <- as.character(fData(valid_subset_GSE72857_cds)[MAP_ordering_genes, 'gene_short_name'])
MAP_cluster_hyper_geometric_results_go <- collect_gsa_hyper_results(valid_subset_GSE72857_cds, gsc = mouse_go_gsc, dp_cluster_group)
plot_gsa_hyper_heatmap(valid_subset_GSE72857_cds, MAP_cluster_hyper_geometric_results_go, significance = 1e-5)

pdf(paste(SI_fig_dir, "MAP_fig_si1d_heatmap.pdf", sep = ''), height = 5, width = 4)
plot_gsa_hyper_heatmap(valid_subset_GSE72857_cds, MAP_cluster_hyper_geometric_results_go, significance = 1e-5) + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

################################################################################################################################################
## feature selection by PCA
# create SI2i (PCA)
########################################################################################################################################################
MAP_pca_res <- selectOrderingGenes2(valid_subset_GSE72857_cds, dim = c(1, 2), top_gene_num = 1000, method = c('PCA'), verbose = T)

#use all pca gene for making the trajectory
MAP <- setOrderingFilter(valid_subset_GSE72857_cds, ordering_genes = MAP_pca_res$gene_vec)
MAP <- reduceDimension(MAP, norm_method = 'log', verbose = T)
MAP <- orderCells(MAP)

pdf(paste(SI_fig_dir, "MAP_pca.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(MAP, color_by = 'cell_type', show_branch_points = F)  + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

##########################################################################################################################
##  feature selection by SLICER
# create SI2i (SLICER)
########################################################################################################################################################
MAP_slicer_genes <- slicer_feature_genes(valid_subset_GSE72857_cds)
MAP_mat <- log(exprs(valid_subset_GSE72857_cds[MAP_slicer_genes, ]) + 1)
pheatmap(MAP_mat, cluster_rows = T, cluster_cols = T, scale = 'row')

plot_pseudotime_heatmap(valid_subset_GSE72857_cds[MAP_slicer_genes, ], use_gene_short_name = F)

#use all slicer gene for making the trajectory 
MAP_slicer <- setOrderingFilter(valid_subset_GSE72857_cds, ordering_genes = row.names(valid_subset_GSE72857_cds)[MAP_slicer_genes])
MAP_slicer <- reduceDimension(MAP_slicer, norm_method = 'log', verbose = T)

pdf(paste(SI_fig_dir, "MAP_slicer.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(MAP_slicer, color_by = 'cell_type', show_branch_points = F)  + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

########################################################################################################################################################
##  feature selection by Over-dispersion
# create SI2i (Over-dispersion)
########################################################################################################################################################
#Over-dispersion genes: 
disp_table <- dispersionTable(valid_subset_GSE72857_cds)
MAP_overdispersion_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)
MAP_dispersion <- setOrderingFilter(valid_subset_GSE72857_cds, MAP_overdispersion_genes$gene_id)
plot_ordering_genes(MAP_dispersion)

MAP_dispersion <- reduceDimension(MAP_dispersion, norm_method = 'log', verbose = T)
MAP_dispersion <- orderCells(MAP_dispersion)

pdf(paste(SI_fig_dir, "MAP_dispersion.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(MAP_dispersion, color_by = 'cell_type', show_branch_points = F)  + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

########################################################################################################################################################
##  feature selection by differential gene expression test
# create SI2i (DEG)
########################################################################################################################################################  
MAP_time_dependent_genes <- differentialGeneTest(valid_subset_GSE72857_cds[,], fullModelFormulaStr="~cell_type", reducedModelFormulaStr="~1", cores=detectCores())
MAP_deg_ordering_genes <- row.names(subset(MAP_time_dependent_genes, qval < 1e-2))
MAP_deg_ordering_genes <- row.names(head(MAP_time_dependent_genes[order(MAP_time_dependent_genes$qval),], n=1000))

MAP_DEG <- setOrderingFilter(valid_subset_GSE72857_cds, MAP_deg_ordering_genes)
MAP_DEG <- reduceDimension(MAP_DEG, norm_method = 'log', verbose = T)
MAP_DEG <- orderCells(MAP_DEG)

pdf(paste(SI_fig_dir, "MAP_DEG.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(MAP_DEG, color_by = 'cell_type', show_branch_points = F)  + monocle:::monocle_theme_opts() + nm_theme()
dev.off()
################################################################################################################################################################################################################################################
# create SI2h (Over-dispersion)
################################################################################################################################################################################################################################################
#use upset to show the overlapping for different feature selection methods: (pca, high-dispersion, time-dependent de, dp feature genes, slicer feature selection)
# example of list input (list of named vectors)
MAP_listInput <- list('High dispersion' = HSMM_overdispersion_genes$gene_id, Slicer = row.names(HSMM_myo)[HSMM_slicer_genes], 
                  PCA = HSMM_pca_res$gene_vec, 'DP genes' = HSMM_myo_ordering_genes, 'Time DEG' = HSMM_deg_ordering_genes)

pdf(paste(SI_fig_dir, "MAP_fig1c.pdf", sep = ''), height = 2.5, width = 5)
upset(fromList(MAP_listInput), nsets = 6, order.by = "freq") + nm_theme()
dev.off()

########################################################################################################################################################
#save the data: 
########################################################################################################################################################
save.image('./RData/fig_si1.RData')

