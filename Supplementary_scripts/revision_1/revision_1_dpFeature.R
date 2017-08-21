rm(list = ls())
#analyze the effects of each parameters to the software on the HSMM dataset

# number of PCs, perplexity, number of genes
library(xacHelper)
library(grid)
# library(devtools)
library(monocle)
clusterCells_Density_Peak <- clusterCells
revision_1_fig_dir <- "./Figures/First_revision/"

# load_all('/Users/xqiu/Dropbox (Personal)/Projects/monocle-dev')

load('./RData/fig1.RData')

#1. determine how many pca dimension you want:
fData(HSMM_myo)$use_for_ordering <- fData(HSMM_myo)$num_cells_expressed > round(ncol(HSMM_myo) / 10)
HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(HSMM_myo, return_all = F)  + geom_vline(xintercept = 15)
pdf(paste(revision_1_fig_dir, "HSMM_tSNE_PC_variance_explained.pdf", sep = ''), height = 1, width = 1)
plot_pc_variance_explained(HSMM_myo, return_all = F)  + geom_vline(xintercept = 15) + nm_theme()
dev.off()

# HSMM@auxClusteringData[["tSNE"]]$variance

HSMM_expressed_genes <- row.names(subset(fData(HSMM_myo),
                                         num_cells_expressed >= 15 &
                                           biotype %in% c("protein_coding", "lincRNA")))

##########################################################################################################################################################################
# A function to do all the process
##########################################################################################################################################################################
for(num_pc in c(2, 4, 6, 8, 10, 12, 15)) { #
  #1. determine how many pca dimension you want:

  fData(HSMM_myo)$use_for_ordering <- fData(HSMM_myo)$num_cells_expressed > round(ncol(HSMM_myo) / 10)
  HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL
  plot_pc_variance_explained(HSMM_myo, return_all = F)

  #2. run reduceDimension with tSNE as the reduction_method
  HSMM_myo <- reduceDimension(HSMM_myo, max_components=2, norm_method = 'log', num_dim = num_pc, reduction_method = 'tSNE', verbose = T)

  p <- plot_cell_clusters(HSMM_myo, color_by = 'Hours', cell_size = 0.25, show_density_peak = F) + nm_theme()
  pdf(paste(revision_1_fig_dir, "HSMM_tSNE_PC_", num_pc, ".pdf", sep = ''), height = 1, width = 1)
  print(p)
  dev.off()

  #3. initial run of clusterCells
  HSMM_myo <- clusterCells(HSMM_myo, verbose = T)

  plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)', show_density = F, rho_threshold = 3.5, delta_threshold = 1.2)
  plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Hours)', show_density = F, rho_threshold = 3.5, delta_threshold = 1.2)

  #6. re-run cluster and skipping calculating the rho_sigma
  HSMM_myo <- clusterCells(HSMM_myo, verbose = T,  rho_threshold = 3.5, delta_threshold = 1.2, skip_rho_sigma = T)

  #7. make the final clustering plot:
  plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)', show_density = F)
  plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Hours)', show_density = F)

  #find the important genes and use that for lineage tree construction
  #perform DEG test across clusters:
  HSMM_myo@expressionFamily <- negbinomial.size()
  pData(HSMM_myo)$Cluster <- as.character(pData(HSMM_myo)$Cluster)
  HSMM_myo_DEG_genes <- differentialGeneTest(HSMM_myo, fullModelFormulaStr = '~Cluster', cores = detectCores() - 2)
  # HSMM_myo_DEG_genes_subset <- HSMM_myo[fData(HSMM_myo)$num_cells_expressed > round(ncol(HSMM_myo) / 10, ]

  #use all DEG gene from the clusters
  HSMM_myo_ordering_genes <- row.names(HSMM_myo_DEG_genes)[order(HSMM_myo_DEG_genes$qval)][1:1000] #row.names(subset(HSMM_myo_DEG_genes, qval < qval_thrsld))

  # HSMM_myo_ordering_genes <- row.names(HSMM_myo_clustering_DEG_genes_subset)[order(HSMM_myo_clustering_DEG_genes_subset$qval)][1:1000]
  HSMM_myo_ordering_genes

  HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes = HSMM_myo_ordering_genes)
  HSMM_myo <- reduceDimension(HSMM_myo, verbose = T, auto_param_selection = F)
  HSMM_myo <- orderCells(HSMM_myo)
  plot_cell_trajectory(HSMM_myo, color_by = 'Hours')

  p <- plot_spanning_tree(HSMM_myo, color_by="Hours", show_branch_points = F, cell_size = 0.25) + nm_theme()
  pdf(paste(revision_1_fig_dir, "HSMM_DM", num_pc, ".pdf", sep = ''), height = 1, width = 1)
  print(p)
  dev.off()

}

##########################################################################################################################################################################
# show the effect of number of genes:
##########################################################################################################################################################################
#1. determine how many pca dimension you want:
fData(HSMM_myo)$use_for_ordering <- fData(HSMM_myo)$num_cells_expressed > round(ncol(HSMM_myo) / 10)
HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(HSMM_myo, return_all = F)

#2. run reduceDimension with tSNE as the reduction_method
HSMM_myo <- reduceDimension(HSMM_myo, max_components=2, norm_method = 'log', num_dim = 4, reduction_method = 'tSNE', verbose = T)

#3. initial run of clusterCells
HSMM_myo <- clusterCells(HSMM_myo, verbose = T)

#4. check the clusters (there are three clusters)
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)', show_density = F)
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Hours)', show_density = F)

#5. also check the decision plot
plot_rho_delta(HSMM_myo, rho_threshold = 3, delta_threshold = 1)
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)', show_density = F, rho_threshold = 3.8, delta_threshold = 4)
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Hours)', show_density = F, rho_threshold = 3.8, delta_threshold = 5)

#6. re-run cluster and skipping calculating the rho_sigma
HSMM_myo <- clusterCells(HSMM_myo, verbose = T,  rho_threshold = 3.8, delta_threshold = 5, skip_rho_sigma = T)

#7. make the final clustering plot:
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)', show_density = F)

#find the important genes and use that for lineage tree construction
#perform DEG test across clusters:
HSMM_myo@expressionFamily <- negbinomial.size()
pData(HSMM_myo)$Cluster <- as.character(pData(HSMM_myo)$Cluster)
HSMM_myo_DEG_genes <- differentialGeneTest(HSMM_myo, fullModelFormulaStr = '~Cluster', cores = detectCores() - 2)
# HSMM_myo_DEG_genes_subset <- HSMM_myo[fData(HSMM_myo)$num_cells_expressed > round(ncol(HSMM_myo) / 10, ]

for(gene_num in c(100, 200, 500, 800, 1200, 1400, 1500, 1800, 2000, 2500, 3000, 4000, 5000, 10000)) {
  #use all DEG gene from the clusters
  HSMM_myo_ordering_genes <- row.names(HSMM_myo_DEG_genes)[order(HSMM_myo_DEG_genes$qval)][1:gene_num] #row.names(subset(HSMM_myo_DEG_genes, qval < qval_thrsld))

  # HSMM_myo_ordering_genes <- row.names(HSMM_myo_clustering_DEG_genes_subset)[order(HSMM_myo_clustering_DEG_genes_subset$qval)][1:1000]
  HSMM_myo_ordering_genes

  HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes = HSMM_myo_ordering_genes)
  HSMM_myo <- reduceDimension(HSMM_myo, verbose = T, auto_param_selection = F)
  HSMM_myo <- orderCells(HSMM_myo)
  plot_cell_trajectory(HSMM_myo, color_by = 'Hours')

  p <- plot_spanning_tree(HSMM_myo, color_by="Hours", show_branch_points = F, cell_size = 0.25) + nm_theme()

  pdf(paste(revision_1_fig_dir, "HSMM_DM_gene_numer_", gene_num, ".pdf", sep = ''), height = 1, width = 1)
  print(p)
  dev.off()
}

##########################################################################################################################################################################
# show the effect of perplexity:
##########################################################################################################################################################################
for(perplexity in c(5, 10, 15, 20, 25, 35, 40, 45)) {
  #1. determine how many pca dimension you want:
  fData(HSMM_myo)$use_for_ordering <- fData(HSMM_myo)$num_cells_expressed > round(ncol(HSMM_myo) / 10)
  HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL
  plot_pc_variance_explained(HSMM_myo, return_all = F)

  #2. run reduceDimension with tSNE as the reduction_method
  HSMM_myo <- reduceDimension(HSMM_myo, max_components=2, norm_method = 'log', num_dim = 4, reduction_method = 'tSNE', verbose = T, perplexity = perplexity)

  p <- plot_cell_clusters(HSMM_myo, color_by="Hours", show_branch_points = F, cell_size = 0.25, show_density_peak = F) + nm_theme()

  pdf(paste(revision_1_fig_dir, "HSMM_DM_perplexity_", perplexity, ".pdf", sep = ''), height = 1, width = 1)
  print(p)
  dev.off()

  #3. initial run of clusterCells
  HSMM_myo <- clusterCells(HSMM_myo, verbose = T)

  #4. check the clusters (there are three clusters)
  plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)', show_density = F)
  plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Hours)', show_density = F)

  #5. also check the decision plot
  plot_rho_delta(HSMM_myo, rho_threshold = 3.8, delta_threshold = 4 )
  plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)', show_density = F, rho_threshold = 3.8, delta_threshold = 4)
  plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Hours)', show_density = F, rho_threshold = 3.8, delta_threshold = 4)

  #6. re-run cluster and skipping calculating the rho_sigma
  HSMM_myo <- clusterCells(HSMM_myo, verbose = T,  rho_threshold = 3.8, delta_threshold = 4, skip_rho_sigma = T)

  #7. make the final clustering plot:
  plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)', show_density = F)

  #find the important genes and use that for lineage tree construction
  #perform DEG test across clusters:
  HSMM_myo@expressionFamily <- negbinomial.size()
  pData(HSMM_myo)$Cluster <- as.character(pData(HSMM_myo)$Cluster)
  HSMM_myo_DEG_genes <- differentialGeneTest(HSMM_myo, fullModelFormulaStr = '~Cluster', cores = detectCores() - 2)
  # HSMM_myo_DEG_genes_subset <- HSMM_myo[fData(HSMM_myo)$num_cells_expressed > round(ncol(HSMM_myo) / 10, ]

  #use all DEG gene from the clusters
  HSMM_myo_ordering_genes <- row.names(HSMM_myo_DEG_genes)[order(HSMM_myo_DEG_genes$qval)][1:1000] #row.names(subset(HSMM_myo_DEG_genes, qval < qval_thrsld))

  # HSMM_myo_ordering_genes <- row.names(HSMM_myo_clustering_DEG_genes_subset)[order(HSMM_myo_clustering_DEG_genes_subset$qval)][1:1000]
  HSMM_myo_ordering_genes

  HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes = HSMM_myo_ordering_genes)
  HSMM_myo <- reduceDimension(HSMM_myo, verbose = T, auto_param_selection = F)
  HSMM_myo <- orderCells(HSMM_myo)
  plot_cell_trajectory(HSMM_myo, color_by = 'Hours')

  p <- plot_spanning_tree(HSMM_myo, color_by="Hours", show_branch_points = F, cell_size = 0.25) + nm_theme()

  pdf(paste(revision_1_fig_dir, "HSMM_DM_perplexity_tree", perplexity, ".pdf", sep = ''), height = 1, width = 1)
  print(p)
  dev.off()
}

##########################################################################################################################################################################
# show the effect of delta: 
##########################################################################################################################################################################
# change the delta and rho
for(delta in c(1, 2, 3, 4, 5)) {
  #1. determine how many pca dimension you want:
  fData(HSMM_myo)$use_for_ordering <- fData(HSMM_myo)$num_cells_expressed > round(ncol(HSMM_myo) / 10)
  HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL
  plot_pc_variance_explained(HSMM_myo, return_all = F)
  
  #2. run reduceDimension with tSNE as the reduction_method
  HSMM_myo <- reduceDimension(HSMM_myo, max_components=2, norm_method = 'log', num_dim = 4, reduction_method = 'tSNE', verbose = T)
  
  #3. initial run of clusterCells
  HSMM_myo <- clusterCells(HSMM_myo, verbose = T)
  
  #4. check the clusters (there are three clusters)
  plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)', show_density = F)
  plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Hours)', show_density = F)
  
  #5. also check the decision plot

  p <- plot_rho_delta(HSMM_myo, rho_threshold = 2, delta_threshold = delta) + nm_theme()
  
  pdf(paste(revision_1_fig_dir, "HSMM_DM_delta_", delta, ".pdf", sep = ''), height = 1, width = 1)
  print(p)
  dev.off()
  # 
  #6. re-run cluster and skipping calculating the rho_sigma
  HSMM_myo <- clusterCells(HSMM_myo, verbose = T,  rho_threshold = 2, delta_threshold = delta, skip_rho_sigma = T)

  #7. make the final clustering plot:
  plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)', show_density = F)

  #find the important genes and use that for lineage tree construction
  #perform DEG test across clusters:
  HSMM_myo@expressionFamily <- negbinomial.size()
  pData(HSMM_myo)$Cluster <- as.character(pData(HSMM_myo)$Cluster)
  HSMM_myo_DEG_genes <- differentialGeneTest(HSMM_myo, fullModelFormulaStr = '~Cluster', cores = detectCores() - 2)
  # HSMM_myo_DEG_genes_subset <- HSMM_myo[fData(HSMM_myo)$num_cells_expressed > round(ncol(HSMM_myo) / 10, ]

  #find the important genes and use that for lineage tree construction
  #perform DEG test across clusters:
  HSMM_myo@expressionFamily <- negbinomial.size()
  pData(HSMM_myo)$Cluster <- as.character(pData(HSMM_myo)$Cluster)
  HSMM_myo_DEG_genes <- differentialGeneTest(HSMM_myo, fullModelFormulaStr = '~Cluster', cores = detectCores() - 2)
  # HSMM_myo_DEG_genes_subset <- HSMM_myo[fData(HSMM_myo)$num_cells_expressed > round(ncol(HSMM_myo) / 10, ]

  #use all DEG gene from the clusters
  HSMM_myo_ordering_genes <- row.names(HSMM_myo_DEG_genes)[order(HSMM_myo_DEG_genes$qval)][1:1000] #row.names(subset(HSMM_myo_DEG_genes, qval < qval_thrsld))

  # HSMM_myo_ordering_genes <- row.names(HSMM_myo_clustering_DEG_genes_subset)[order(HSMM_myo_clustering_DEG_genes_subset$qval)][1:1000]
  HSMM_myo_ordering_genes

  HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes = HSMM_myo_ordering_genes)
  HSMM_myo <- reduceDimension(HSMM_myo, verbose = T, auto_param_selection = F)
  HSMM_myo <- orderCells(HSMM_myo)
  plot_cell_trajectory(HSMM_myo, color_by = 'Hours')

  p <- plot_spanning_tree(HSMM_myo, color_by="Hours", show_branch_points = F, cell_size = 0.25) + nm_theme()

  pdf(paste(revision_1_fig_dir, "HSMM_DM_delta_tree", delta, ".pdf", sep = ''), height = 1, width = 1)
  print(p)
  dev.off()
}

##########################################################################################################################################################################
# run analysis on the blood dataset (not very useful here)
##########################################################################################################################################################################
save.image('./RData/revision_1_dpFeature.RData')

