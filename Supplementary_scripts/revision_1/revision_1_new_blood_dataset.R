library(monocle)
library(xacHelper)

# analyze the new blood dataset
res_a <- read.csv('../csv_data/STEM_CLOUD/GSE75478_transcriptomics_facs_indeces_filtered_I1 (1).csv', row.names = 1)
res_a.1 <- read.csv('../csv_data/STEM_CLOUD/GSE75478_transcriptomics_facs_indeces_filtered_I2 (1).csv.gz', row.names = 1)
res_b <- read.csv('../csv_data/STEM_CLOUD/GSE75478_transcriptomics_normalized_filtered_I1 (2).csv.gz', row.names = 1)
res_b.1 <- read.csv('../csv_data/STEM_CLOUD/GSE75478_transcriptomics_normalized_filtered_I2 (1).csv.gz', row.names = 1)
# res_c <- read.csv('./csv_data/STEM_CLOUD/GSE75478_transcriptomics_raw_filtered_I2 (1).csv')
# res_d <- read.csv('./csv_data/STEM_CLOUD/GSE75478_RAW (1)/GSM1955701_counts_plate3_A_1.csv.gz')

setdiff(c(colnames(res_b), colnames(res_b.1)), c(colnames(res_a), colnames(res_a.1)))
intersect_facs <- intersect(row.names(res_a), row.names(res_a.1))
intersect_genes <- intersect(row.names(res_b), row.names(res_b.1))

fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = row.names(res_b[intersect_genes, ]), row.names = row.names(res_b[intersect_genes, , ])))
pd <- new("AnnotatedDataFrame", data = as.data.frame(t(cbind(res_a[intersect_facs, ], res_a.1[intersect_facs, ]))))

# Now, make a new CellDataSet using the RNA counts
res_all <- cbind(res_b[intersect_genes, ], res_b.1[intersect_genes, ])
CLOUD_cds <- newCellDataSet(as(as.matrix(res_all[, row.names(pd)]), 'sparseMatrix'), 
                                            phenoData = pd, 
                                            featureData = fd,
                                            lowerDetectionLimit=1,
                                            expressionFamily=gaussianff())

#1. determine how many pca dimension you want:
fData(CLOUD_cds)$use_for_ordering <- T
CLOUD_cds@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(CLOUD_cds, return_all = F, norm_method = 'none')

#2. run reduceDimension with tSNE as the reduction_method
CLOUD_cds <- reduceDimension(CLOUD_cds, max_components=2, norm_method = 'none', num_dim = 15, reduction_method = 'tSNE', verbose = T)

p <- plot_cell_clusters(CLOUD_cds, color_by = 'Hours', cell_size = 0.25) + nm_theme()
pdf(paste(revision_1_fig_dir, "CLOUD_tSNE_PC_", num_pc, ".pdf", sep = ''), height = 1, width = 1)
p
dev.off()

#3. initial run of clusterCells_Density_Peak
CLOUD_cds <- clusterCells(CLOUD_cds, verbose = T)

plot_cell_clusters(CLOUD_cds, color_by = 'as.factor(Cluster)', show_density = F, rho_threshold = 3.5, delta_threshold = 1.2)
plot_cell_clusters(CLOUD_cds, color_by = 'as.factor(Hours)', show_density = F, rho_threshold = 3.5, delta_threshold = 1.2)

#6. re-run cluster and skipping calculating the rho_sigma
CLOUD_cds <- clusterCells_Density_Peak(CLOUD_cds, verbose = T,  rho_threshold = 3.5, delta_threshold = 1.2, skip_rho_sigma = T)

#7. make the final clustering plot:
plot_cell_clusters(CLOUD_cds, color_by = 'as.factor(Cluster)', show_density = F)
plot_cell_clusters(CLOUD_cds, color_by = 'as.factor(Hours)', show_density = F)

#find the important genes and use that for lineage tree construction
#perform DEG test across clusters:
pData(CLOUD_cds)$Cluster <- as.character(pData(CLOUD_cds)$Cluster)
CLOUD_cds_DEG_genes <- differentialGeneTest(CLOUD_cds, fullModelFormulaStr = '~Cluster', cores = detectCores() - 2)
# CLOUD_cds_DEG_genes_subset <- CLOUD_cds[fData(CLOUD_cds)$num_cells_expressed > round(ncol(CLOUD_cds) / 10, ]

#use all DEG gene from the clusters
CLOUD_cds_ordering_genes <- row.names(CLOUD_cds_DEG_genes)[order(CLOUD_cds_DEG_genes$qval)][1:1000] #row.names(subset(CLOUD_cds_DEG_genes, qval < qval_thrsld))

# CLOUD_cds_ordering_genes <- row.names(CLOUD_cds_clustering_DEG_genes_subset)[order(CLOUD_cds_clustering_DEG_genes_subset$qval)][1:1000]
CLOUD_cds_ordering_genes

CLOUD_cds <- setOrderingFilter(CLOUD_cds, ordering_genes = CLOUD_cds_ordering_genes)

plot_pc_variance_explained(CLOUD_cds, norm_method = 'none')
CLOUD_cds <- reduceDimension(CLOUD_cds, verbose = T, auto_param_selection = T, max_components = 10, norm_method = 'none')
CLOUD_cds <- orderCells(CLOUD_cds)
plot_cell_trajectory(CLOUD_cds, color_by = 'FACS_cd34')

plot_complex_cell_trajectory(CLOUD_cds, color_by = 'FACS_cd34')

##########################################################################################################################################################################
pd <- new("AnnotatedDataFrame", data = as.data.frame(t(cbind(res_a[intersect_facs, ], res_a.1[intersect_facs, ]))))

# Now, make a new CellDataSet using the RNA counts
res_all <- as.data.frame(cbind(res_a[intersect_facs, ], res_a.1[intersect_facs, ]))[-(1:4), ]
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = row.names(res_all), row.names = row.names(res_all))) #+ min(res_all)

protein_CLOUD_cds <- newCellDataSet(as(as.matrix(res_all[, row.names(pd)]), 'sparseMatrix'), 
                            phenoData = pd, 
                            featureData = fd,
                            lowerDetectionLimit=-5772,
                            expressionFamily=gaussianff())

protein_CLOUD_cds <- setOrderingFilter(protein_CLOUD_cds, ordering_genes = row.names(protein_CLOUD_cds))

protein_CLOUD_cds <- reduceDimension(protein_CLOUD_cds, verbose = T, auto_param_selection = T, max_components = 8, norm_method = 'none')
CLOUD_cds <- orderCells(protein_CLOUD_cds)
plot_cell_trajectory(protein_CLOUD_cds, color_by = 'FACS_cd34')

##########################################################################################################################################################################
# save the data here
##########################################################################################################################################################################
save.image('../RData/revision_1_new_blood_dataset.RData')