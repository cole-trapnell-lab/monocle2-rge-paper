library(monocle)
library(R.matlab)
library(plyr)

#test with Wishbone result: 
scRNA_seq_wb_res <- read.csv("../wishbone/new_version/wishbone/notebooks/scRNA_seq_wb_res.txt", sep = '\t', row.names = 1)
scRNA_seq_exprs_res <- read.csv("../wishbone/new_version/wishbone/notebooks/scRNA_seq_exprs_res.txt", sep = '\t', row.names = 1)
scRNA_seq_exprs_res <- t(scRNA_seq_exprs_res)

gene_ann <- data.frame(gene_short_name = row.names(scRNA_seq_exprs_res), row.names = row.names(scRNA_seq_exprs_res))
pd <- new("AnnotatedDataFrame",data=scRNA_seq_wb_res)
fd <- new("AnnotatedDataFrame",data=gene_ann)

wishbone_scRNA_seq <- newCellDataSet(as(scRNA_seq_exprs_res, "sparseMatrix"),phenoData = pd,featureData =fd,
                               expressionFamily = negbinomial.size(),
                               lowerDetectionLimit=1)

wishbone_scRNA_seq <- estimateSizeFactors(wishbone_scRNA_seq)
wishbone_scRNA_seq <- estimateDispersions(wishbone_scRNA_seq)

wishbone_scRNA_seq <- reduceDimension(wishbone_scRNA_seq, verbose = T, norm_method = 'log') 
wishbone_scRNA_seq <- orderCells(wishbone_scRNA_seq)
plot_cell_trajectory(wishbone_scRNA_seq)
plot_cell_trajectory(wishbone_scRNA_seq, color_by = 'as.factor(branch)')

#confirm with the cluster information from the original paper 
MAP_cells_clusters <- readMat('../data/MAP_cells_clusters.mat')
MAP_cells_clusters <- MAP_cells_clusters$MAP.cells.clusters
storage.mode(MAP_cells_clusters) <- 'numeric'
MAP_cells_clusters <- read.csv('../csv_data/MAP.csv', header = F)
row.names(MAP_cells_clusters) <- MAP_cells_clusters$V1

pData(wishbone_scRNA_seq)$cluster <- NA
pData(wishbone_scRNA_seq)[row.names(MAP_cells_clusters), 'cluster'] <- as.character(MAP_cells_clusters$V2)

pData(wishbone_scRNA_seq)$cluster <- revalue(as.character(pData(wishbone_scRNA_seq)$cluster), 
                                             c("1" = 'erythroid', "2" = 'erythroid', "3" = 'erythroid', "4" = 'erythroid', "5" = 'erythroid', "6" = 'erythroid', 
                                               "7" = 'CMP', "8" = 'CMP', "9" = 'CMP', "10" = 'CMP',
                                               "11" = 'DC', 
                                               "12" = 'GMP', "13" = 'GMP', "14" = 'GMP', "15" = 'GMP', "16" = 'GMP', "17" = 'GMP', "18" = 'GMP', 
                                               "19" = 'lymphoid'))

cor(pData(wishbone_scRNA_seq)$trajectory, pData(wishbone_scRNA_seq)$Pseudotime)
qplot(pData(wishbone_scRNA_seq)$trajectory, pData(wishbone_scRNA_seq)$Pseudotime)

calClusteringMetrics(pData(wishbone_scRNA_seq)$cluster[!is.na(pData(wishbone_scRNA_seq)$cluster)], 
  pData(wishbone_scRNA_seq)$branch[!is.na(pData(wishbone_scRNA_seq)$cluster)])

calClusteringMetrics(pData(wishbone_scRNA_seq)$cluster[!is.na(pData(wishbone_scRNA_seq)$cluster)], 
                     pData(wishbone_scRNA_seq)$State[!is.na(pData(wishbone_scRNA_seq)$cluster)])

plot_cell_trajectory(wishbone_scRNA_seq, color_by = 'as.factor(cluster)')

# wb_subset_all_GSE72857_cds <- as.matrix(exprs(all_GSE72857_cds)[, colnames(wishbone_scRNA_seq)])
# 
# wishbone_scRNA_seq <- estimateSizeFactors(wishbone_scRNA_seq)
# wishbone_scRNA_seq <- estimateDispersions(wishbone_scRNA_seq)
# all_GSE72857_cds <- reduceDimension(all_GSE72857_cds, norm_method = 'log', verbose = T)
# all_GSE72857_cds <- orderCells(all_GSE72857_cds)
                               
# test the MASS cytometry dataset: 
mass_wb_res <- read.csv("../wishbone/new_version/wishbone/notebooks/MASS_cytometry_wb_res.txt", sep = '\t', row.names = 1)
mass_exprs_res <- read.csv("../wishbone/new_version/wishbone/notebooks/MASS_cytometry_exprs_res.txt", sep = '\t', row.names = 1)
mass_exprs_res <- t(mass_exprs_res)

gene_ann <- data.frame(gene_short_name = row.names(mass_exprs_res), row.names = row.names(mass_exprs_res))
pd <- new("AnnotatedDataFrame",data=mass_wb_res)
fd <- new("AnnotatedDataFrame",data=gene_ann)

wishbone_mass <- newCellDataSet(as(mass_exprs_res, "sparseMatrix"),phenoData = pd,featureData =fd,
                                     expressionFamily = gaussianff(),
                                     lowerDetectionLimit=1)

wishbone_mass <- reduceDimension(wishbone_mass, norm_method = 'none', verbose = T) 

# 
# ###################################################################################################################################################################################################################
# # save the data
# ###################################################################################################################################################################################################################
# 
save.image('./RData/revision_1_test_wishbone_res.RData')

