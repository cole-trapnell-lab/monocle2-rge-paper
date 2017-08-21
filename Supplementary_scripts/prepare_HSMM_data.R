rm(list = ls())

# Chunk 1: package_loads
library(Biobase)
library(knitr)
library(reshape2)
library(ggplot2)

knitr::opts_chunk$set(autodep=TRUE, cache=FALSE, warning=FALSE)
set.seed(0)

# Chunk 2: init_monocle
library(HSMMSingleCell)
library(monocle)
data(HSMM_expr_matrix)
data(HSMM_gene_annotation)
data(HSMM_sample_sheet)

# Chunk 7: build_cell_data_Set_RPC
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)

# First create a CellDataSet from the relative expression levels
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),   
                       phenoData = pd, 
                       featureData = fd)

# Next, use it to estimate RNA counts
rpc_matrix <- relative2abs(HSMM)

# Now, make a new CellDataSet using the RNA counts
HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                       phenoData = pd, 
                       featureData = fd,
                       lowerDetectionLimit=0.5,
                       expressionFamily=negbinomial.size())


# Chunk 8: estimate_size_and_dispersion
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

# Chunk 9: detect_genes
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))

# Chunk 10: show_pData
print(head(pData(HSMM)))

# Chunk 12: show_mRNA_totals
pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))


HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) + 2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) - 2*sd(log10(pData(HSMM)$Total_mRNAs)))
qplot(Total_mRNAs, data=pData(HSMM), color=Hours, geom="density") + 
  geom_vline(xintercept=lower_bound) + 
  geom_vline(xintercept=upper_bound)

HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound & 
               pData(HSMM)$Total_mRNAs < upper_bound]								  
HSMM <- detectGenes(HSMM, min_expr = 0.1)

# Chunk 13: lognormal_plot
# Log-transform each value in the expression matrix.
L <- log(exprs(HSMM[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily"
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
qplot(value, geom="density", data=melted_dens_df) +  stat_function(fun = dnorm, size=0.5, color='red') + 
  xlab("Standardized log(FPKM)") +
  ylab("Density")

# Chunk 14: setup_gates
MYF5_id <- row.names(subset(fData(HSMM), gene_short_name == "MYF5"))
ANPEP_id <- row.names(subset(fData(HSMM), gene_short_name == "ANPEP"))

cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Myoblast", classify_func=function(x) {x[MYF5_id,] >= 1})
cth <- addCellType(cth, "Fibroblast", classify_func=function(x)
{x[MYF5_id,] < 1 & x[ANPEP_id,] > 1})

# Chunk 15: count_cells_unsup
HSMM <- classifyCells(HSMM, cth, 0.1)

# Chunk 16: count_cells_unsup_readout
table(pData(HSMM)$CellType)

pie <- ggplot(pData(HSMM), aes(x = factor(1), fill = factor(CellType))) +
  geom_bar(width = 1)
pie + coord_polar(theta = "y") + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

# Chunk 17: cluster_cells_unsup_gene_pick
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)

# Chunk 18: cluster_cells_unsup_no_covariate
# HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(HSMM, return_all = F) # norm_method = 'log',
HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 6, 
                        reduction_method = 'tSNE', verbose = T) 
HSMM <- clusterCells(HSMM,
                     num_clusters=2)
plot_cell_clusters(HSMM, 1, 2, color="CellType", markers=c("MYF5", "ANPEP"))

# Chunk 19: cluster_cells_unsup_plot_by_media
plot_cell_clusters(HSMM, 1, 2, color="Media")

# Chunk 20: cluster_cells_unsup_control_for_media
HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 2, reduction_method = 'tSNE', 
                        residualModelFormulaStr="~Media + num_genes_expressed", verbose = T) #
HSMM <- clusterCells(HSMM, num_clusters=2)
plot_cell_clusters(HSMM, 1, 2, color="CellType")

# Chunk 21: cluster_cells_unsup_plot_by_cell_type
HSMM <- clusterCells(HSMM, num_clusters=2)
plot_cell_clusters(HSMM, 1, 2, color="Cluster") + facet_wrap(~CellType)

# Chunk 22: cluster_cells_diff_table
marker_diff <- markerDiffTable(HSMM[expressed_genes,], 
                               cth, 
                               residualModelFormulaStr="~Media + num_genes_expressed",
                               cores=1)

# Chunk 23: cluster_cells_semisup_show_marker_spec
candidate_clustering_genes <- row.names(subset(marker_diff, qval < 0.05))
marker_spec <- calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth)
head(selectTopMarkers(marker_spec, 3))

# Chunk 24: cluster_cells_semisup_pick_genes
semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes)
plot_ordering_genes(HSMM)

# Chunk 25: cluster_cells_semisup_clustering_no_impute
plot_pc_variance_explained(HSMM, return_all = F) # norm_method = 'log',
HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 2, reduction_method = 'tSNE', 
                        residualModelFormulaStr="~Media + num_genes_expressed", verbose = T) 
HSMM <- clusterCells(HSMM, num_clusters=2,
                     clustering_genes=semisup_clustering_genes) 

plot_cell_clusters(HSMM, 1, 2, color="CellType")

# Chunk 26: cluster_cells_semisup_clustering_with_impute
HSMM <- clusterCells(HSMM,
                     num_clusters=2, 
                     frequency_thresh=0.1,
                     cell_type_hierarchy=cth,
                     clustering_genes=row.names(subset(marker_diff, qval < 0.05)))
plot_cell_clusters(HSMM, 1, 2, color="CellType", markers = c("MYF5", "ANPEP"))

# Chunk 27: count_cells_semisup_pie
pie <- ggplot(pData(HSMM), aes(x = factor(1), fill = factor(CellType))) +
  geom_bar(width = 1)
pie + coord_polar(theta = "y") + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

# Chunk 28: select_myoblasts
HSMM_myo <- HSMM[,pData(HSMM)$CellType == "Myoblast"]	
HSMM_myo <- estimateDispersions(HSMM_myo)

################################################################################################################################################################################################################################################
# run DPFeature to select ordering genes 
################################################################################################################################################################################################################################################
HSMM_expressed_genes <- row.names(subset(fData(HSMM_myo),
                                         num_cells_expressed >= 15 &
                                           biotype %in% c("protein_coding", "lincRNA")))

GM_state <- function(cds){
  T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
  as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))])
}

HSMM_myo <- detectGenes(HSMM_myo, min_expr=0.1)
fData(HSMM_myo)$use_for_ordering <- fData(HSMM_myo)$num_cells_expressed > 0.05 * ncol(HSMM_myo)
plot_pc_variance_explained(HSMM_myo, return_all = F) #look at the plot and decide how many dimensions you need. It is determined by a huge drop of variance at that dimension. pass that number to num_dim in the next function.

HSMM_myo <- reduceDimension(HSMM_myo, max_components=2, norm_method = 'log', num_dim = 3, 
                              reduction_method = 'tSNE', verbose = T)
HSMM_myo <- clusterCells(HSMM_myo, verbose = F)

plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)')
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Hours)')
plot_rho_delta(HSMM_myo, rho_threshold = 2.6, delta_threshold = 4)

HSMM_myo <- clusterCells(HSMM_myo,  
  rho_threshold = 2.6, 
  delta_threshold = 4, 
  skip_rho_sigma = T, 
  verbose = F)

plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)')
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Hours)')

clustering_DEG_genes <- differentialGeneTest(HSMM_myo[HSMM_expressed_genes,], 
  fullModelFormulaStr = '~Cluster', 
  cores = detectCores())

HSMM_myo_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000] 
HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes = HSMM_myo_ordering_genes)
HSMM_myo <- reduceDimension(HSMM_myo)
HSMM_myo <- orderCells(HSMM_myo)
HSMM_myo <- orderCells(HSMM_myo, root_state=GM_state(HSMM_myo))
plot_cell_trajectory(HSMM_myo, color_by="Hours")

################################################################################################################################################################################################################################################
# save.image('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/RData/fig1_run_monocle_vignette_cluster_cells.RData')
save.image('./RData/fig1_run_monocle_vignette_cluster_cells.RData')
################################################################################################################################################################################################################################################
