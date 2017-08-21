################################################################################################################################################
#run densityPeak (or other clustering method) for clustering
#use the selectOrderingGene function?
################################################################################################################################################

#source function
source('./scripts/function.R', echo = T)
source('./scripts/plotting.R', echo = T)

#load libraries
library(monocle)
library(destiny)
library(xacHelper)
library(grid)

gene_num <- 1000
num_cells_expressed_percent <- 0.1
qval_thrsld <- 1e-2

main_fig_dir <- "./Figures/main_figures/"
SI_fig_dir <- "./Figures/supplementary_figures/"

################################################################################################################################################
#lung data 
load("./RData/prepare_lung_data.RData")
################################################################################################################################################
#1. determine how many pca dimension you want: 
absolute_cds <- recreate_cds(absolute_cds)

absolute_cds <- detectGenes(absolute_cds)
fData(absolute_cds)$use_for_ordering <- F

num_cells_expressed <- round(num_cells_expressed_percent * ncol(absolute_cds))
fData(absolute_cds)$use_for_ordering[fData(absolute_cds)$num_cells_expressed > num_cells_expressed] <- T
lung_pc_variance <- plot_pc_variance_explained(absolute_cds, return_all = T)

#2. run reduceDimension with tSNE as the reduction_method 
# absolute_cds <- setOrderingFilter(absolute_cds, quake_id)
absolute_cds <- reduceDimension(absolute_cds, max_components=2, norm_method = 'log', reduction_method = 'tSNE', num_dim = 3,  verbose = T)

#3. initial run of clusterCells_Density_Peak
absolute_cds <- clusterCells_Density_Peak(absolute_cds, verbose = T)

#4. check the clusters 
plot_cell_clusters(absolute_cds, color_by = 'as.factor(Cluster)', show_density = F)
plot_cell_clusters(absolute_cds, color_by = 'Time', show_density = F)
plot_cell_clusters(absolute_cds, color_by = 'State', show_density = F)

#5. also check the decision plot 
plot_rho_delta(absolute_cds, rho_threshold = 3, delta_threshold = 5 )
plot_cell_clusters(absolute_cds, color_by = 'Time', show_density = F, rho_threshold = 3, delta_threshold = 5)

#6. re-run cluster and skipping calculating the rho_sigma 
absolute_cds <- clusterCells_Density_Peak(absolute_cds, verbose = T,  rho_threshold = 3, delta_threshold = 5, skip_rho_sigma = T)

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

# lung_ordering_genes <- row.names(lung_clustering_DEG_genes_subset)[order(lung_clustering_DEG_genes_subset$qval)][1:500] #1971

absolute_cds <- setOrderingFilter(absolute_cds, ordering_genes = lung_ordering_genes)
absolute_cds <- reduceDimension(absolute_cds, norm_method = 'log', verbose = T)
absolute_cds <- orderCells(absolute_cds)

pdf(paste(SI_fig_dir, "lung_cluster_DEG_qval_0.01.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(absolute_cds, color_by = 'Time') + monocle::monocle_theme_opts() + nm_theme()
dev.off()

pdf(paste(SI_fig_dir, "lung_cluster_DEG_qval_0.01_helper.pdf", sep = ''))
plot_cell_trajectory(absolute_cds, color_by = 'Time') 
dev.off()
#save the data for HSMM and lung 
save.image('./RData/clustering_cells_select_features_tmp.RData')

################################################################################################################################################
#blood data 
load('./RData/fig4.RData')
################################################################################################################################################
URMM_all_fig1b <- recreate_cds(URMM_all_fig1b)

num_cells_expressed <- round(ncol(URMM_all_fig1b) * num_cells_expressed_percent)
#1. determine how many pca dimension you want: 
URMM_all_fig1b <- detectGenes(URMM_all_fig1b)
fData(URMM_all_fig1b)$use_for_ordering[fData(URMM_all_fig1b)$num_cells_expressed > num_cells_expressed] <- T
URMM_pc_variance <- plot_pc_variance_explained(URMM_all_fig1b, return_all = T)

#2. run reduceDimension with tSNE as the reduction_method 
URMM_all_fig1b <- reduceDimension(URMM_all_fig1b, max_components=2, norm_method = 'log', reduction_method = 'tSNE', num_dim = 7,  verbose = T)

#3. initial run of clusterCells_Density_Peak
URMM_all_fig1b <- clusterCells_Density_Peak(URMM_all_fig1b, verbose = T)

#4. check the clusters 
plot_cell_clusters(URMM_all_fig1b, color_by = 'as.factor(Cluster)', show_density = T, show_density_peak = T)

plot_cell_clusters(URMM_all_fig1b, color_by = 'Type', show_density = F)
# plot_cell_clusters(URMM_all_fig1b, color_by = 'State', show_density = F) + scale_color_brewer(palette = 'Set1') + facet_wrap(~paper_cluster)
# plot_cell_clusters(URMM_all_fig1b, color_by = 'Cluster', show_density = F) + scale_color_brewer(palette = 'Set1') + facet_wrap(~paper_cluster)

#5. also check the decision plot 
plot_rho_delta(URMM_all_fig1b, rho_threshold = 6, delta_threshold = 10 )

#6. re-run cluster and skipping calculating the rho_sigma 
URMM_all_fig1b <- clusterCells_Density_Peak(URMM_all_fig1b, verbose = T,  rho_threshold = 6, delta_threshold = 10, skip_rho_sigma = T)

#7. make the final clustering plot: 
plot_cell_clusters(URMM_all_fig1b, color_by = 'as.factor(Cluster)', show_density = F)
plot_cell_clusters(URMM_all_fig1b, color_by = 'as.factor(Type)', show_density = F)
plot_cell_clusters(URMM_all_fig1b, color_by = 'as.factor(cluster)', show_density = F)

#perform DEG test across clusters: 
URMM_all_fig1b@expressionFamily <- negbinomial.size()
pData(URMM_all_fig1b)$Cluster <- factor(pData(URMM_all_fig1b)$Cluster)
URMM_clustering_DEG_genes <- differentialGeneTest(URMM_all_fig1b, fullModelFormulaStr = '~Cluster', cores = detectCores() - 2)
URMM_clustering_DEG_genes_subset <- URMM_clustering_DEG_genes[fData(URMM_all_fig1b)$num_cells_expressed > num_cells_expressed, ]

#use all DEG gene from the clusters
#URMM_ordering_genes <- row.names(subset(URMM_clustering_DEG_genes, qval < qval_thrsld))
URMM_ordering_genes <- row.names(URMM_clustering_DEG_genes)[order(URMM_clustering_DEG_genes$qval)][1:1000]
URMM_ordering_genes

URMM_all_fig1b <- setOrderingFilter(URMM_all_fig1b, ordering_genes = URMM_ordering_genes)
URMM_all_fig1b <- reduceDimension(URMM_all_fig1b, norm_method = 'log', verbose = T, scaling = T) #, param.gamma = 100, ncenter = 100
URMM_all_fig1b <- orderCells(URMM_all_fig1b)
plot_cell_trajectory(URMM_all_fig1b, color_by = 'Type')

pData(URMM_all_fig1b)$paper_cluster <- pData(URMM_all_abs_all_dispersion)[colnames(URMM_all_fig1b), 'cluster']

pdf(paste(SI_fig_dir, "URMM_cluster_DEG_qval_0.01.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(URMM_all_fig1b, color_by = 'paper_cluster') + monocle::monocle_theme_opts() + nm_theme()
dev.off()

pdf(paste(SI_fig_dir, "URMM_cluster_DEG_qval_0.01_helper.pdf", sep = ''))
plot_cell_trajectory(URMM_all_fig1b, color_by = 'paper_cluster') 
dev.off()
########################################################################################################################################################
#run on the full dataset
URMM_all_abs <- recreate_cds(URMM_all_abs)
URMM_all_abs <- setOrderingFilter(URMM_all_abs, ordering_genes = URMM_ordering_genes)
URMM_all_abs <- reduceDimension(URMM_all_abs, norm_method = 'log', verbose = T) # 
URMM_all_abs <- orderCells(URMM_all_abs)
plot_cell_trajectory(URMM_all_abs, color_by = 'Type')

plot_cell_trajectory(URMM_all_abs, color_by = 'Type') + facet_wrap(~Type)

pData(URMM_all_abs)$paper_cluster <- NA 
pData(URMM_all_abs)[colnames(URMM_all_fig1b), 'paper_cluster'] <- as.character(pData(URMM_all_fig1b)[, 'paper_cluster'])

pdf(paste(SI_fig_dir, "all_data_URMM_cluster_DEG_qval_0.01.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(URMM_all_abs, color_by = 'paper_cluster') + monocle::monocle_theme_opts() + nm_theme()
dev.off()

pdf(paste(SI_fig_dir, "all_data_URMM_cluster_DEG_qval_0.01_helper.pdf", sep = ''))
plot_cell_trajectory(URMM_all_abs, color_by = 'paper_cluster') 
dev.off()
########################################################################################################################################################
# pData(URMM_all_abs)$paper_cluster <- pData(URMM_all_abs_all_dispersion)[colnames(URMM_all_abs), 'cluster']
#1. determine how many pca dimension you want: 
URMM_all_abs <- detectGenes(URMM_all_abs)
fData(URMM_all_abs)$use_for_ordering[fData(URMM_all_abs)$num_cells_expressed > ncol(URMM_all_abs) * num_cells_expressed_percent] <- T
URMM_pc_variance <- plot_pc_variance_explained(URMM_all_abs, return_all = T)

#2. run reduceDimension with tSNE as the reduction_method 
URMM_all_abs <- reduceDimension(URMM_all_abs, max_components=2, norm_method = 'log', reduction_method = 'tSNE', num_dim = 4,  verbose = T)

#3. initial run of clusterCells_Density_Peak
URMM_all_abs <- clusterCells_Density_Peak(URMM_all_abs, verbose = T)

#4. check the clusters 
plot_cell_clusters(URMM_all_abs, color_by = 'as.factor(Cluster)', show_density = F)

plot_cell_clusters(URMM_all_abs, color_by = 'Type', show_density = F)
plot_cell_clusters(URMM_all_abs, color_by = 'Cluster', show_density = F) #+ scale_color_brewer(palette = 'Set1') + facet_wrap(~paper_cluster)
plot_cell_clusters(URMM_all_abs, color_by = 'paper_cluster', show_density = F) #+ scale_color_brewer(palette = 'Set1') + facet_wrap(~paper_cluster)

#5. also check the decision plot 
plot_rho_delta(URMM_all_abs, rho_threshold = 5, delta_threshold = 10 )

#6. re-run cluster and skipping calculating the rho_sigma 
URMM_all_abs <- clusterCells_Density_Peak(URMM_all_abs, verbose = T,  rho_threshold = 5, delta_threshold = 10, skip_rho_sigma = T)

#7. make the final clustering plot: 
plot_cell_clusters(URMM_all_abs, color_by = 'as.factor(Cluster)', show_density = F)
plot_cell_clusters(URMM_all_abs, color_by = 'as.factor(Type)', show_density = F)

#perform DEG test across clusters: 
URMM_all_abs@expressionFamily <- negbinomial.size()
pData(URMM_all_abs)$Cluster <- factor(pData(URMM_all_abs)$Cluster)
URMM_all_clustering_DEG_genes <- differentialGeneTest(URMM_all_abs, fullModelFormulaStr = '~Cluster', cores = detectCores() - 2)
URMM_all_clustering_DEG_genes_subset <- URMM_all_clustering_DEG_genes[fData(URMM_all_abs)$num_cells_expressed > ncol(URMM_all_abs) * 0.05, ]

#use all DEG gene from the clusters
URMM_all_ordering_genes <- row.names(URMM_all_clustering_DEG_genes_subset)[order(URMM_all_clustering_DEG_genes_subset$qval)][1:1000]

URMM_all_ordering_genes <- row.names(URMM_all_clustering_DEG_genes)[URMM_all_clustering_DEG_genes$qval < 0.01]
URMM_all_ordering_genes

URMM_all_abs <- setOrderingFilter(URMM_all_abs, ordering_genes = URMM_all_ordering_genes)
URMM_all_abs <- reduceDimension(URMM_all_abs, norm_method = 'log', verbose = T) #
plot_cell_trajectory(URMM_all_abs, color_by = 'Type')

pdf(paste(SI_fig_dir, "URMM_all_cluster_DEG_qval_0.01.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(URMM_all_abs, color_by = 'Type') + monocle::monocle_theme_opts() + nm_theme()
dev.off()

pdf(paste(SI_fig_dir, "URMM_all_cluster_DEG_qval_0.01_helper.pdf", sep = ''))
plot_cell_trajectory(URMM_all_abs, color_by = 'Type') 
dev.off()

URMM_clustering_DEG_genes_subset <- URMM_clustering_DEG_genes[fData(URMM_all_fig1b)$num_cells_expressed > num_cells_expressed, ]

########################################################################################################################################################
# create the MAR dataset with all genes
################################################################################################################################################
#read valid_subset_GSE72857
valid_subset_GSE72857_exprs <- read.table('./data/GSE72857_umitab.txt', header = T, row.names = 1)

#filtering cells to include only the ones which were assigned a cluster id: 
MAP_cells_clusters <- read.csv('./data/MAP.csv', header = F)
row.names(MAP_cells_clusters) <- MAP_cells_clusters$V1
design_mat <- read.table('./csv_data/GSE72857_experimental_design.txt', header = T, row.names = 1, skip = 19, sep = '\t')
design_mat$cluster <- MAP_cells_clusters[row.names(design_mat), 'V2']
valid_design_mat <- subset(design_mat, !is.na(cluster))

setdiff(colnames(valid_subset_GSE72857_exprs[, row.names(subset(design_mat, Batch_desc %in% c( 'Unsorted myeloid')) )])[(apply(valid_subset_GSE72857_exprs[, row.names(subset(design_mat, Batch_desc %in% c( 'Unsorted myeloid')) )], 2, sum) > 500)], row.names(MAP_cells_clusters))
colSums(valid_subset_GSE72857_exprs[, c('W31450', 'W36988')]) #both cells have UMI counts 501 

fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = row.names(valid_subset_GSE72857_exprs), row.names = row.names(valid_subset_GSE72857_exprs)))
pd <- new("AnnotatedDataFrame", data = valid_design_mat)

# Now, make a new CellDataSet using the RNA counts
valid_subset_GSE72857_cds <- newCellDataSet(as(as.matrix(valid_subset_GSE72857_exprs[, row.names(valid_design_mat)]), 'sparseMatrix'), 
                                            phenoData = pd, 
                                            featureData = fd,
                                            lowerDetectionLimit=1,
                                            expressionFamily=negbinomial.size())

pData(valid_subset_GSE72857_cds)$cell_type <- revalue(as.character(pData(valid_subset_GSE72857_cds)$cluster), 
                                                      c("1" = 'erythroid', "2" = 'erythroid', "3" = 'erythroid', "4" = 'erythroid', "5" = 'erythroid', "6" = 'erythroid', 
                                                        "7" = 'CMP', "8" = 'CMP', "9" = 'CMP', "10" = 'CMP',
                                                        "11" = 'DC', 
                                                        "12" = 'GMP', "13" = 'GMP', "14" = 'GMP', "15" = 'GMP', "16" = 'GMP', "17" = 'GMP', "18" = 'GMP', 
                                                        "19" = 'lymphoid'))

#remove all lymphoid cells
valid_subset_GSE72857_cds <- valid_subset_GSE72857_cds[, pData(valid_subset_GSE72857_cds)$cell_type != 'lymphoid']
valid_subset_GSE72857_cds <- estimateSizeFactors(valid_subset_GSE72857_cds)
valid_subset_GSE72857_cds <- estimateDispersions(valid_subset_GSE72857_cds)

################################################################################################################################################

#1. determine how many pca dimension you want: 
valid_subset_GSE72857_cds <- detectGenes(valid_subset_GSE72857_cds)
fData(valid_subset_GSE72857_cds)$use_for_ordering <- F

num_cells_expressed <- round(num_cells_expressed_percent * ncol(valid_subset_GSE72857_cds)) #
fData(valid_subset_GSE72857_cds)$use_for_ordering[fData(valid_subset_GSE72857_cds)$num_cells_expressed > num_cells_expressed] <- T

valid_subset_GSE72857_cds@auxClusteringData[["tSNE"]]$variance_explained <- NULL
MAP_pc_variance <- plot_pc_variance_explained(valid_subset_GSE72857_cds, return_all = T)

#2. run reduceDimension with tSNE as the reduction_method 
# valid_subset_GSE72857_cds <- setOrderingFilter(valid_subset_GSE72857_cds, quake_id)
valid_subset_GSE72857_cds <- reduceDimension(valid_subset_GSE72857_cds, max_components=2, norm_method = 'log', reduction_method = 'tSNE', num_dim = 5,  verbose = T) #, residualModelFormulaStr = '~groups.1'

#check the embedding on each PCA components: 

#3. initial run of clusterCells_Density_Peak
valid_subset_GSE72857_cds <- clusterCells_Density_Peak(valid_subset_GSE72857_cds, verbose = T)

#4. check the clusters 
plot_cell_clusters(valid_subset_GSE72857_cds, color_by = 'as.factor(Cluster)', show_density = F)
plot_cell_clusters(valid_subset_GSE72857_cds, color_by = 'cluster', show_density = F)
plot_cell_clusters(valid_subset_GSE72857_cds, color_by = 'cell_type', show_density = F)

#5. also check the decision plot 
plot_rho_delta(valid_subset_GSE72857_cds, rho_threshold = 3, delta_threshold = 40)
plot_cell_clusters(valid_subset_GSE72857_cds, color_by = 'cell_type', rho_threshold = 3, delta_threshold = 40)

#6. re-run cluster and skipping calculating the rho_sigma 
valid_subset_GSE72857_cds <- clusterCells_Density_Peak(valid_subset_GSE72857_cds, verbose = T,  rho_threshold = 3, delta_threshold = 40, skip_rho_sigma = T)

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
valid_subset_GSE72857_cds <- reduceDimension(valid_subset_GSE72857_cds, verbose = T)
valid_subset_GSE72857_cds <- orderCells(valid_subset_GSE72857_cds)

plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'cell_type') 

pdf(paste(SI_fig_dir, "MAP_cluster_DEG_qval_1k.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'cell_type') + monocle_theme_opts() + nm_theme()
dev.off()

#run dpt and wishbone: 
MAP_dpt_res <- run_new_dpt(valid_subset_GSE72857_cds, normalize = F)

qplot(MAP_dpt_res$dm$DC1, MAP_dpt_res$dm$DC2, color = pData(valid_subset_GSE72857_cds)$cell_type)

#comparing the pseudotime: 
qplot(MAP_dpt_res$pt, pData(valid_subset_GSE72857_cds)$Pseudotime)
qplot(pData(valid_subset_GSE72857_cds)$Pseudotime)
qplot(MAP_dpt_res$pt)

########################################################################################################################################################
#save the data 
########################################################################################################################################################
save.image('./RData/dpFeature_for_all_datasets.RData')
