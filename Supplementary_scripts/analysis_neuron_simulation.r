rm(list = ls())

####################################################################################################################################################################################
# Please note that you may need to use the develop version of igraph to run this script; otherwise your R session will be crashed when creating a CDS with Monocle  
####################################################################################################################################################################################

####################################################################################################################################################################################
#load all package
####################################################################################################################################################################################
library(monocle)
library(R.matlab)

source('./scripts/function.R')
source('./scripts/plotting.R')

main_fig_dir <- "./Figures/main_figures/"
SI_fig_dir <- "./Figures/supplementary_figures/"

####################################################################################################################################################################################
# load all required simulation data 
####################################################################################################################################################################################
#CNS simulation data:
cell_simulate <- readMat('./mat_data/cell_simulate.mat')
all_cell_simulation <- cell_simulate$cell.simulate[, 1:400, ] #time 0-20 are the period two branches appear
row.names(all_cell_simulation) <-  c('Pax6', 'Mash1', 'Brn2', 'Zic1', 'Tuj1', 'Hes5', 'Scl', 'Olig2', 'Stat3', 'Myt1L', 'Aldh1L', 'Sox8', 'Mature')

#obtain the corresponding lineage for each simulation run:
neuron_cell_ids <- which(all_cell_simulation['Mash1', 400, ] > 3)
astrocyte_cell_ids <- which(all_cell_simulation['Scl', 400, ] > 2)
oligodendrocyte_cell_ids <- which(all_cell_simulation['Olig2', 400, ] > 2)

cols <- c("lap" = "black", "original" = "gray", "pt" = "red")
cell_type_cols <- c("neuron" = "#F2746B", "astrocyte" = "#29B14A", "oligodendrocyte" = "#6E95CD")
  
####################################################################################################################################################################################
# create the necessary cds for making all the following figures
####################################################################################################################################################################################
gene_names <- c('Pax6', 'Mash1', 'Brn2', 'Zic1', 'Tuj1', 'Hes5', 'Scl', 'Olig2', 'Stat3', 'Myt1L', 'Aldh1L', 'Sox8', 'Mature')

neuron_exprs_mat <- as.matrix(all_cell_simulation[, , neuron_cell_ids[1]])
dimnames(neuron_exprs_mat) = list(gene_names, paste('Cell_', rep(1:400, 1), sep = ''))

astrocyte_exprs_mat <- as.matrix(all_cell_simulation[, , astrocyte_cell_ids[1]])
dimnames(astrocyte_exprs_mat) = list(gene_names, paste('Cell_', rep(1:400, 1), sep = ''))

oligodendrocyte_exprs_mat <- as.matrix(all_cell_simulation[, , oligodendrocyte_cell_ids[1]])
dimnames(oligodendrocyte_exprs_mat) = list(gene_names, paste('Cell_', rep(1:400, 1), sep = ''))

PD <- data.frame(Time = rep(1:400, 1), row.names = paste('Cell_', rep(1:400, 1), sep = ''))
FD <- data.frame(gene_short_names = gene_names, row.names = gene_names)

pd <- new("AnnotatedDataFrame",data=PD)
fd <- new("AnnotatedDataFrame",data=FD)

make_cds <- function (exprs_matrix, pd, fd, lowerDetectionLimit = 1, expressionFamily) {
  cds <- newCellDataSet(exprs_matrix, phenoData = pd, featureData = fd, lowerDetectionLimit = 1, expressionFamily)
  return(cds)
}

neuron_sim_cds <- newCellDataSet(neuron_exprs_mat, pd, fd, expressionFamily = negbinomial())
pData(neuron_sim_cds)$cell_type <- 'neuron'
astrocyte_sim_cds <- newCellDataSet(astrocyte_exprs_mat, pd, fd, expressionFamily = negbinomial())
pData(neuron_sim_cds)$cell_type <- 'astrocyte'
oligodendrocyte_sim_cds <- newCellDataSet(oligodendrocyte_exprs_mat, pd, fd, expressionFamily = negbinomial())
pData(neuron_sim_cds)$cell_type <- 'oligodendrocyte'

na_exprs_mat <- as.matrix(all_cell_simulation[, , c(neuron_cell_ids[1], astrocyte_cell_ids[1])])
dim(na_exprs_mat) <- c(13, 400 * 2)

no_exprs_mat <- as.matrix(all_cell_simulation[, , c(neuron_cell_ids[1], oligodendrocyte_cell_ids[1])])
dim(no_exprs_mat) <- c(13, 400 * 2)

ao_exprs_mat <- as.matrix(all_cell_simulation[, , c(astrocyte_cell_ids[1], oligodendrocyte_cell_ids[2])])
dim(ao_exprs_mat) <- c(13, 400 * 2)

nao_exprs_mat <- as.matrix(all_cell_simulation[, , c(neuron_cell_ids[1], astrocyte_cell_ids[1], oligodendrocyte_cell_ids[1])])
dim(nao_exprs_mat) <- c(13, 400 * 3)

PD2 <- data.frame(Time = rep(1:400, 2), row.names = paste('Cell_', 1:800, sep = ''))
PD3 <- data.frame(Time = rep(1:400, 3), row.names = paste('Cell_', 1:1200, sep = ''))

pd2 <- new("AnnotatedDataFrame",data=PD2)
pd3 <- new("AnnotatedDataFrame",data=PD3)

FD <- data.frame(gene_short_names = gene_names, row.names = gene_names)
fd <- new("AnnotatedDataFrame",data=FD)

na_sim_cds <- make_cds(na_exprs_mat, pd2, fd, expressionFamily = negbinomial())
pData(na_sim_cds)$cell_type <- rep(c('neuron', 'astrocyte'), each = 400)
no_sim_cds <- make_cds(no_exprs_mat, pd2, fd, expressionFamily = negbinomial())
pData(no_sim_cds)$cell_type <- rep(c('neuron', 'oligodendrocyte'), each = 400)
ao_sim_cds <- make_cds(ao_exprs_mat, pd2, fd, expressionFamily = negbinomial())
pData(ao_sim_cds)$cell_type <- rep(c('astrocyte', 'oligodendrocyte'), each = 400)
nao_sim_cds <- make_cds(nao_exprs_mat, pd3, fd, expressionFamily = negbinomial())
pData(nao_sim_cds)$cell_type <- rep(c('neuron', 'astrocyte', 'oligodendrocyte'), each = 400)

neuron_sim_cds <- estimateSizeFactors(neuron_sim_cds)
pData(neuron_sim_cds)$Size_Factor <- 1
astrocyte_sim_cds <- estimateSizeFactors(astrocyte_sim_cds)
pData(astrocyte_sim_cds)$Size_Factor <- 1
oligodendrocyte_sim_cds <- estimateSizeFactors(oligodendrocyte_sim_cds)
pData(oligodendrocyte_sim_cds)$Size_Factor <- 1

na_sim_cds <- estimateSizeFactors(na_sim_cds)
pData(na_sim_cds)$Size_Factor <- 1

no_sim_cds <- estimateSizeFactors(no_sim_cds)
pData(no_sim_cds)$Size_Factor <- 1

ao_sim_cds <- estimateSizeFactors(ao_sim_cds)
pData(ao_sim_cds)$Size_Factor <- 1

nao_sim_cds <- estimateSizeFactors(nao_sim_cds)
pData(nao_sim_cds)$Size_Factor <- 1

####################################################################################################################################################################################
#create figure SI5b (figSI5b.1.pdf is the file used in the paper, same pattern as below)
####################################################################################################################################################################################
neuron_sim_cds <- setOrderingFilter(neuron_sim_cds, ordering_genes = row.names(neuron_sim_cds))
neuron_sim_cds <- reduceDimension(neuron_sim_cds, norm_method = 'none', verbose = T, auto_param_selection = F, maxIter = 100,  pseudo_expr=0, scaling = F)
neuron_sim_cds <- orderCells(neuron_sim_cds)

pData(neuron_sim_cds)$Cell_type <- 'neuron'
pdf(paste(main_fig_dir, "figSI5b.pdf", sep = ''), height = 1.1, width = 1.1)
plot_cell_trajectory(neuron_sim_cds, color_by = 'Cell_type', cell_size = 0.2) +  #geom_point(aes(neuron_DDRTree_res$Z[1, ], neuron_DDRTree_res$Z[2, ]))
  scale_x_reverse() + scale_y_reverse() + nm_theme() + scale_color_manual(values = cell_type_cols)#
dev.off()

pdf(paste(main_fig_dir, "figSI5b.1.pdf", sep = ''), height = 1.1, width = 1.1)
qplot(neuron_sim_cds@reducedDimS[1, ], neuron_sim_cds@reducedDimS[2, ], alpha = 0.5, color = I('gray'), size = 0.2) + 
  geom_point(aes(neuron_sim_cds@reducedDimK[1, ], neuron_sim_cds@reducedDimK[2, ], size = 3 * c(1:400) / 1000, color = I(cell_type_cols[1]), alpha = 0.5)) + 
  scale_size(range = range(3 * c(1:400) / 1000), limits = range(3 * c(1:400) / 1000)) + 
  scale_x_reverse() + scale_y_reverse() + nm_theme() + xlab('Component 1') + ylab('Component 2') + monocle:::monocle_theme_opts()
dev.off()

pdf(paste(main_fig_dir, "figSI5b.1.pdf", sep = ''), height = 1.1, width = 1.1)
qplot(neuron_sim_cds@reducedDimK[1, ], neuron_sim_cds@reducedDimK[2, ], color = I('black'), size = 0.1) + 
  geom_point(aes(neuron_sim_cds@reducedDimS[1, ], neuron_sim_cds@reducedDimS[2, ], 
                 size = 3 *  rep(1:400, 1) / 1000, alpha = 0.5, color = I(rep(cell_type_cols[c(1)], each = 400)))) + 
  scale_size(range = range(3 *  rep(1:400, 1) / 1000), limits = range(3 *  rep(1:400, 1) / 1000)) + 
  scale_x_reverse() + scale_y_reverse() + nm_theme() + xlab('Component 1') + ylab('Component 2') + monocle:::monocle_theme_opts()
dev.off()

####################################################################################################################################################################################
#create figure SI5c (figSI5c.1)
####################################################################################################################################################################################
astrocyte_sim_cds <- setOrderingFilter(astrocyte_sim_cds, ordering_genes = row.names(neuron_sim_cds))
astrocyte_sim_cds <- reduceDimension(astrocyte_sim_cds, norm_method = 'none', verbose = T, auto_param_selection = F, maxIter = 100,  pseudo_expr=0, scaling = F)
astrocyte_sim_cds <- orderCells(astrocyte_sim_cds)
plot_cell_trajectory(astrocyte_sim_cds, color_by = 'Time')

oligodendrocyte_sim_cds <- setOrderingFilter(oligodendrocyte_sim_cds, ordering_genes = row.names(neuron_sim_cds))
oligodendrocyte_sim_cds <- reduceDimension(oligodendrocyte_sim_cds, norm_method = 'none', verbose = T, auto_param_selection = F, maxIter = 100,  pseudo_expr=0, scaling = F)
oligodendrocyte_sim_cds <- orderCells(oligodendrocyte_sim_cds)
plot_cell_trajectory(oligodendrocyte_sim_cds, color_by = 'Time')

na_sim_cds <- setOrderingFilter(na_sim_cds, ordering_genes = row.names(neuron_sim_cds))
na_sim_cds <- reduceDimension(na_sim_cds, norm_method = 'none', verbose = T, auto_param_selection = F, maxIter = 100, pseudo_expr=0, scaling = F)
na_sim_cds <- orderCells(na_sim_cds)

plot_cell_trajectory(na_sim_cds, color_by = 'Time')
pData(na_sim_cds)$Cell_type <- rep(c('neuron', 'astrocyte'), each = 400)
pdf(paste(main_fig_dir, "figSI5c.pdf", sep = ''), height = 1.1, width = 1.1)
plot_cell_trajectory(na_sim_cds, color_by = 'Cell_type', cell_size = 0.2, show_branch_points = F) +  #geom_point(aes(neuron_DDRTree_res$Z[1, ], neuron_DDRTree_res$Z[2, ]))
  scale_x_reverse() + scale_y_reverse() + nm_theme() + scale_color_manual(values = cell_type_cols)#
dev.off()

pdf(paste(main_fig_dir, "figSI5c.1.pdf", sep = ''), height = 1.1, width = 1.1)
qplot(na_sim_cds@reducedDimK[1, ], na_sim_cds@reducedDimK[2, ], color = I('black'), size = 0.1) + 
  geom_point(aes(na_sim_cds@reducedDimS[1, ], na_sim_cds@reducedDimS[2, ], 
                 size = 3 *  rep(1:400, 2) / 1000, alpha = 0.5, color = I(rep(cell_type_cols[c(1:2)], each = 400)))) + 
  scale_size(range = range(3 *  rep(1:400, 2) / 1000), limits = range(3 *  rep(1:400, 2) / 1000)) + 
  scale_x_reverse() + scale_y_reverse() + nm_theme() + xlab('Component 1') + ylab('Component 2') + monocle:::monocle_theme_opts()
dev.off()

####################################################################################################################################################################################
#create figure SI5d
####################################################################################################################################################################################
ao_sim_cds <- setOrderingFilter(ao_sim_cds, ordering_genes = row.names(neuron_sim_cds))
ao_sim_cds <- reduceDimension(ao_sim_cds, norm_method = 'none', verbose = T, auto_param_selection = F, maxIter = 10,  pseudo_expr=0, scaling = F)
ao_sim_cds <- orderCells(ao_sim_cds)
plot_cell_trajectory(ao_sim_cds, color_by = 'Time')

plot_cell_trajectory(ao_sim_cds, color_by = 'Time')
pData(ao_sim_cds)$Cell_type <- rep(c('astrocyte', 'oligodendrocyte'), each = 400)
pdf(paste(main_fig_dir, "figSI5d.pdf", sep = ''), height = 1.1, width = 1.1)
plot_cell_trajectory(ao_sim_cds, color_by = 'Cell_type', cell_size = 0.2, show_branch_points = F) +  #geom_point(aes(neuron_DDRTree_res$Z[1, ], neuron_DDRTree_res$Z[2, ]))
  scale_x_reverse() + scale_y_reverse() + nm_theme() + scale_color_manual(values = cell_type_cols)#
dev.off()

pdf(paste(main_fig_dir, "figSId.1.pdf", sep = ''), height = 1.1, width = 1.1)
qplot(ao_sim_cds@reducedDimK[1, ], ao_sim_cds@reducedDimK[2, ], color = I('black'), size = 0.1) + 
  geom_point(aes(ao_sim_cds@reducedDimS[1, ], ao_sim_cds@reducedDimS[2, ], 
                 size = 3 *  rep(1:400, 2) / 1000, alpha = 0.5, color = I(rep(cell_type_cols[c(2, 3)], each = 400)))) + 
  scale_size(range = range(3 *  rep(1:400, 2) / 1000), limits = range(3 *  rep(1:400, 2) / 1000)) + 
  scale_x_reverse() + scale_y_reverse() + nm_theme() + xlab('Component 1') + ylab('Component 2') + monocle:::monocle_theme_opts()
dev.off()

pdf(paste(main_fig_dir, "figSI5d_helper.pdf", sep = ''))
qplot(ao_sim_cds@reducedDimS[1, ], ao_sim_cds@reducedDimS[2, ], alpha = 0.5, color = I('gray'), size = 0.2) + 
  geom_point(aes(ao_sim_cds@reducedDimK[1, ], ao_sim_cds@reducedDimK[2, ], size = 3 *  rep(1:400, 2) / 1000, color = I(rep(cell_type_cols[c(1, 3)], each = 400)), alpha = 0.5)) + 
  scale_size(range = range(3 *  rep(1:400, 2) / 1000), limits = range(3 *  rep(1:400, 2) / 1000)) + 
  scale_x_reverse() + scale_y_reverse() + nm_theme() + xlab('Component 1') + ylab('Component 2')
dev.off()

####################################################################################################################################################################################
#create figure SI5e
####################################################################################################################################################################################
nao_sim_cds <- setOrderingFilter(nao_sim_cds, ordering_genes = row.names(neuron_sim_cds))
nao_sim_cds <- reduceDimension(nao_sim_cds, norm_method = 'none', verbose = T, auto_param_selection = F, maxIter = 100, max_components = 3, pseudo_expr=0, scaling = F)
nao_sim_cds <- orderCells(nao_sim_cds)
plot_cell_trajectory(nao_sim_cds, color_by = 'Time')

nao_sim_cds_iter_1 <- reduceDimension(nao_sim_cds, norm_method = 'none', verbose = T, auto_param_selection = F, maxIter = 1, max_components = 3, pseudo_expr=0, scaling = F)
nao_sim_cds_iter_1 <- orderCells(nao_sim_cds_iter_1)

plot_cell_trajectory(nao_sim_cds_iter_1, color_by = 'Time')
pData(nao_sim_cds_iter_1)$Cell_type <- rep(c('neuron', 'astrocyte', 'oligodendrocyte'), each = 400)
pdf(paste(main_fig_dir, "SI5e.pdf", sep = ''), height = 1.1, width = 1.1)
plot_cell_trajectory(nao_sim_cds_iter_1, color_by = 'Cell_type', cell_size = 0.2, show_branch_points = T) +  #"Time"
  #scale_size(range = range(3 *  rep(1:400, 3) / 1000), limits = range(3 *  rep(1:400, 3) / 1000)) + #geom_point(aes(neuron_DDRTree_res$Z[1, ], neuron_DDRTree_res$Z[2, ]))
  scale_x_reverse() + scale_y_reverse() + nm_theme() + scale_color_manual(values = cell_type_cols)#
dev.off()

pdf(paste(main_fig_dir, "SI5e.1.pdf", sep = ''), height = 1.1, width = 1.1)
qplot(nao_sim_cds@reducedDimK[1, ], nao_sim_cds@reducedDimK[2, ], color = I('black'), size = 0.1) + 
  geom_point(aes(nao_sim_cds_iter_1@reducedDimS[1, ], nao_sim_cds_iter_1@reducedDimS[2, ], 
             size = 3 *  rep(1:400, 3) / 1000, alpha = 0.5, color = as.factor(rep(1:3, each = 400)))) + 
  scale_size(range = range(3 *  rep(1:400, 3) / 1000), limits = range(3 *  rep(1:400, 3) / 1000)) + 
  scale_x_reverse() + scale_y_reverse() + nm_theme() + xlab('Component 1') + ylab('Component 2') + monocle:::monocle_theme_opts()
dev.off()

####################################################################################################################################################################################
#save the simulation data cds for future usage: 
####################################################################################################################################################################################
save(file = './RData/neuron_sim_data.RData', nao_sim_cds, ao_sim_cds, na_sim_cds, neuron_sim_cds, astrocyte_sim_cds, oligodendrocyte_sim_cds)

####################################################################################################################################################################################
# run simplePPT algorithm on the simulated data (Figures from here are not used in the paper)
####################################################################################################################################################################################
#run matlab code and get the best lambda based on the gap statistics

#principal tree method
# bandwidth: 3
# lambda: 0.0036
params.lambda <- 0.0036
params.bandwidth <- 3

num_steps <- 200
num_cells <- 3
subset_data_two_cells <- all_cell_simulation[, 1:num_steps, c(neuron_cell_ids[1], astrocyte_cell_ids[1], oligodendrocyte_cell_ids[1])]
dim(subset_data_two_cells) <- c(13, num_steps * num_cells)

pt_res <- principal_tree(subset_data_two_cells, MU = NULL, lambda = params.lambda, bandwidth = params.bandwidth, verbose = T);

pdf('./Figures/supplementary_figures/principal_path_neuron_branch.pdf', height = 1.5, width = 1.5)
qplot(subset_data_two_cells[2, ], subset_data_two_cells[6, ]) + geom_point(aes( x = pt_res$MU[2, ], y = pt_res$MU[6, ]), color = 'red', size = 0.3) + 
  xlab('Mash1') + ylab('Hes5') + monocle:::monocle_theme_opts() + 
  nm_theme() + scale_size(range = c(0.2, 0.3), limits = c(0.25, 0.25))
dev.off()

pdf('./Figures/supplementary_figures/principal_path_ao_branch.pdf', height = 1.5, width = 1.5)
qplot(subset_data_two_cells[7, ], subset_data_two_cells[8, ]) + geom_point(aes( x = pt_res$MU[7, ], y = pt_res$MU[8, ]), color = 'red', size = 0.3) + 
  xlab('Scl') + ylab('Olig2') + 
  nm_theme() + scale_size(range = c(0.3, 0.3), limits = c(0.25, 0.25)) + monocle:::monocle_theme_opts()
dev.off()

qplot(subset_data_two_cells[2, ], subset_data_two_cells[6, ]) + geom_point(aes( x = pt_res$MU[2, ], y = pt_res$MU[6, ]), color = 'red')
qplot(subset_data_two_cells[7, ], subset_data_two_cells[8, ]) + geom_point(aes( x = pt_res$MU[7, ], y = pt_res$MU[8, ]), color = 'red')

#write the data in matlab and run it in matlab: 
# writeMat('./mat_data/subset_data_two_cells.mat', X = subset_data_two_cells)

row.names(subset_data_two_cells) <-  c('Pax6', 'Mash1', 'Brn2', 'Zic1', 'Tuj1', 'Hes5', 'Scl', 'Olig2', 'Stat3', 'Myt1L', 'Aldh1L', 'Sox8', 'Mature')
colnames(subset_data_two_cells) <- paste('cell_', 1:(num_steps * num_cells), sep = '')
mlt_subset_data_two_cells <- melt(subset_data_two_cells)
mlt_subset_data_two_cells$Var2 <- 'original'

row.names(all_data_combined) <-  c('Pax6', 'Mash1', 'Brn2', 'Zic1', 'Tuj1', 'Hes5', 'Scl', 'Olig2', 'Stat3', 'Myt1L', 'Aldh1L', 'Sox8', 'Mature')
colnames(all_data_combined) <- paste('cell_', 1:ncol(all_data_combined), sep = '')
mlt_all_data_combined <- melt(all_data_combined)
mlt_all_data_combined$Var2 <- 'original'

pt_path <- pt_res$MU
row.names(pt_path) <-  c('Pax6', 'Mash1', 'Brn2', 'Zic1', 'Tuj1', 'Hes5', 'Scl', 'Olig2', 'Stat3', 'Myt1L', 'Aldh1L', 'Sox8', 'Mature')
colnames(pt_path) <- paste('cell_', 1:(num_steps * num_cells), sep = '')
mlt_pt_path <- melt(pt_path)
mlt_pt_path$Var2 <- 'pt'

####################################################################################################################################################################################
# use Tang Ying'data to generate the result for the least action path 
# Panel H
####################################################################################################################################################################################
Tangying_lap <- readMat('/Users/xqiu/Downloads/ansT10N150New.mat')

Tangying_lap$ans[1, 1, 1]$Path
valid_ids <- which(unlist(lapply(Tangying_lap$ans[1, 1, 1]$Path, length)) != 0)
valid_ids <- c(18, 24, 30)
lap_df <- do.call(cbind, lapply(Tangying_lap$ans[1, 1, 1]$Path[valid_ids], function(x) x[[1]]))

row.names(lap_df) <-  c('Pax6', 'Mash1', 'Brn2', 'Zic1', 'Tuj1', 'Hes5', 'Scl', 'Olig2', 'Stat3', 'Myt1L', 'Aldh1L', 'Sox8')
colnames(lap_df) <- paste('cell_', 1:ncol(lap_df), sep = '')
mlt_lap_df <- melt(lap_df)
mlt_lap_df$Var2 <- 'lap'

mlt_all_df <- Reduce(rbind, list(mlt_subset_data_two_cells, mlt_pt_path, mlt_lap_df)) #mlt_all_data_combined

pdf(paste(main_fig_dir, "SI5h.pdf", sep = ''), height = 1.1, width = 1.1)
qplot(subset(mlt_all_df)[mlt_all_df$Var1 == 'Mash1', 'value'], subset(mlt_all_df)[mlt_all_df$Var1 == 'Hes5', 'value'], color = subset(mlt_all_df)[mlt_all_df$Var1 == 'Hes5', 'Var2'], size = 0.25) + 
  #geom_point(size = 0.2) + #, alpha = 0.3 
  xlab('Mash1') + ylab('Hes5') + #stat_density2d() + 
  scale_size(range = c(0.25, 0.25), limits = c(0.25, 0.25)) + nm_theme() + scale_color_manual(values = cols, name = "Type") #+ monocle:::monocle_theme_opts()
dev.off()

pdf(paste(main_fig_dir, "SI5h_helper.pdf", sep = ''))
qplot(subset(mlt_all_df)[mlt_all_df$Var1 == 'Mash1', 'value'], subset(mlt_all_df)[mlt_all_df$Var1 == 'Hes5', 'value'], color = subset(mlt_all_df)[mlt_all_df$Var1 == 'Hes5', 'Var2'], size = 0.25) + 
  geom_point(size = 0.3) + #, alpha = 0.3 
  xlab('Mash1') + ylab('Hes5') + #stat_density2d() + 
  scale_size(range = c(0.2, 0.2),  limits = c(0.25, 0.25))
dev.off()

# this doesn't work, convert into png with pdf 
# png(paste(main_fig_dir, "fig2l.png", sep = ''), height = 1.1, width = 1.1)
# qplot(subset(mlt_all_df)[mlt_all_df$Var1 == 'Mash1', 'value'], subset(mlt_all_df)[mlt_all_df$Var1 == 'Hes5', 'value'], color = subset(mlt_all_df)[mlt_all_df$Var1 == 'Hes5', 'Var2'], size = 0.25) +  
#   xlab('Mash1') + ylab('Hes5') + geom_point(size = 0.3) + 
#   nm_theme() + scale_size(range = c(0.1, 0.3),  limits = c(0.25, 0.25)) + monocle:::monocle_theme_opts()
# dev.off()

####################################################################################################################################################################################
#panel I
####################################################################################################################################################################################
pdf(paste(main_fig_dir, "SI5i.pdf", sep = ''), height = 1.1, width = 1.1)
qplot(subset(mlt_all_df)[mlt_all_df$Var1 == 'Scl', 'value'], subset(mlt_all_df)[mlt_all_df$Var1 == 'Olig2', 'value'], size = 0.25, color = subset(mlt_all_df)[mlt_all_df$Var1 == 'Scl', 'Var2']) +
  xlab('Scl') + ylab('Olig2') + scale_size(range = c(0.2, 0.2),  limits = c(0.25, 0.25)) + monocle:::monocle_theme_opts() + 
  nm_theme()  + scale_color_manual(values = cols, name = "Type")
dev.off()

# png(paste(main_fig_dir, "fig2m_simplePPT.png", sep = ''), height = 1.1, width = 1.1)
# qplot(subset(mlt_all_df)[mlt_all_df$Var1 == 'Scl', 'value'], subset(mlt_all_df)[mlt_all_df$Var1 == 'Olig2', 'value'], color = subset(mlt_all_df)[mlt_all_df$Var1 == 'Scl', 'Var2'], size = 0.3) +
#   xlab('Scl') + ylab('Olig2') + 
#   nm_theme() + monocle:::monocle_theme_opts() + scale_size(range = c(0.2, 0.2),  limits = c(0.25, 0.25)) + scale_color_manual(values = cols, name = "Type")
# dev.off()

pdf(paste(main_fig_dir, "SI5i_helper.pdf", sep = ''))
qplot(subset(mlt_all_df)[mlt_all_df$Var1 == 'Scl', 'value'], subset(mlt_all_df)[mlt_all_df$Var1 == 'Olig2', 'value'], color = subset(mlt_all_df)[mlt_all_df$Var1 == 'Scl', 'Var2'], size = 0.3) +
  xlab('Scl') + ylab('Olig2') + scale_size(range = c(0.2, 0.2), limits = c(0.25, 0.25)) + scale_color_manual(values = cols, name = "Type")
dev.off()

####################################################################################################################################################################################
# Use the reverse embedding points to create SI5F, G: 
####################################################################################################################################################################################

# nao_sim_cds_reverse <- reverseEnbeddingCDS(nao_sim_cds)
# reverse_exprs <- reducedDimW(nao_sim_cds) %*% reducedDimS(nao_sim_cds)

# mlt_pt_path <- melt(pt_path)
row.names(reverse_exprs) <- row.names(pt_path)[1:nrow(reverse_exprs)]
mlt_pt_path <- melt(reverse_exprs)

mlt_pt_path$Var2 <- 'pt'

mlt_all_df <- Reduce(rbind, list(mlt_subset_data_two_cells, mlt_pt_path, mlt_lap_df)) #mlt_all_data_combined

pdf(paste(main_fig_dir, "SI5f.pdf", sep = ''), height = 1.1, width = 1.1)
qplot(subset(mlt_all_df)[mlt_all_df$Var1 == 'Mash1', 'value'], subset(mlt_all_df)[mlt_all_df$Var1 == 'Hes5', 'value'], color = subset(mlt_all_df)[mlt_all_df$Var1 == 'Hes5', 'Var2'], size = 0.25) + 
  #geom_point(size = 0.2) + #, alpha = 0.3 
  xlab('Mash1') + ylab('Hes5') + #stat_density2d() + 
  scale_size(range = c(0.25, 0.25), limits = c(0.25, 0.25)) + nm_theme() + scale_color_manual(values = cols, name = "Type") #+ monocle:::monocle_theme_opts()
dev.off()

pdf(paste(main_fig_dir, "SI5f_helper.pdf", sep = ''))
qplot(subset(mlt_all_df)[mlt_all_df$Var1 == 'Mash1', 'value'], subset(mlt_all_df)[mlt_all_df$Var1 == 'Hes5', 'value'], color = subset(mlt_all_df)[mlt_all_df$Var1 == 'Hes5', 'Var2'], size = 0.25) + 
  geom_point(size = 0.3) + #, alpha = 0.3 
  xlab('Mash1') + ylab('Hes5') + #stat_density2d() + 
  scale_size(range = c(0.2, 0.2),  limits = c(0.25, 0.25))
dev.off()

# this doesn't work, convert into png with pdf 
# png(paste(main_fig_dir, "fig2l.png", sep = ''), height = 1.1, width = 1.1)
# qplot(subset(mlt_all_df)[mlt_all_df$Var1 == 'Mash1', 'value'], subset(mlt_all_df)[mlt_all_df$Var1 == 'Hes5', 'value'], color = subset(mlt_all_df)[mlt_all_df$Var1 == 'Hes5', 'Var2'], size = 0.25) +  
#   xlab('Mash1') + ylab('Hes5') + geom_point(size = 0.3) + 
#   nm_theme() + scale_size(range = c(0.1, 0.3),  limits = c(0.25, 0.25)) + monocle:::monocle_theme_opts()
# dev.off()

####################################################################################################################################################################################
#panel M
####################################################################################################################################################################################
pdf(paste(main_fig_dir, "SI5g.pdf", sep = ''), height = 1.1, width = 1.1)
qplot(subset(mlt_all_df)[mlt_all_df$Var1 == 'Scl', 'value'], subset(mlt_all_df)[mlt_all_df$Var1 == 'Olig2', 'value'], size = 0.25, color = subset(mlt_all_df)[mlt_all_df$Var1 == 'Scl', 'Var2']) +
  xlab('Scl') + ylab('Olig2') + scale_size(range = c(0.2, 0.2),  limits = c(0.25, 0.25)) + monocle:::monocle_theme_opts() + 
  nm_theme()  + scale_color_manual(values = cols, name = "Type")
dev.off()

png(paste(main_fig_dir, "SI5g.png", sep = ''), height = 1.1, width = 1.1)
qplot(subset(mlt_all_df)[mlt_all_df$Var1 == 'Scl', 'value'], subset(mlt_all_df)[mlt_all_df$Var1 == 'Olig2', 'value'], color = subset(mlt_all_df)[mlt_all_df$Var1 == 'Scl', 'Var2'], size = 0.3) +
  xlab('Scl') + ylab('Olig2') + 
  nm_theme() + monocle:::monocle_theme_opts() + scale_size(range = c(0.2, 0.2),  limits = c(0.25, 0.25)) + scale_color_manual(values = cols, name = "Type")
dev.off()

pdf(paste(main_fig_dir, "SI5g_helper.pdf", sep = ''))
qplot(subset(mlt_all_df)[mlt_all_df$Var1 == 'Scl', 'value'], subset(mlt_all_df)[mlt_all_df$Var1 == 'Olig2', 'value'], color = subset(mlt_all_df)[mlt_all_df$Var1 == 'Scl', 'Var2'], size = 0.3) +
  xlab('Scl') + ylab('Olig2') + scale_size(range = c(0.2, 0.2), limits = c(0.25, 0.25)) + scale_color_manual(values = cols, name = "Type")
dev.off()

####################################################################################################################################################################################
# the figures below are not used in the paper
####################################################################################################################################################################################
# use the W matrix to convert the theoretical to low dimension 
# This mat has least action paths in ans.Path. It's a 6*6 matrix,
# the 1 is the intermediate between neuron and astrocyte / oligodendrocyte,
# the 2 is  intermediate states of the astrocyte and oligodendrocyte,
# the 4 is  neuron,
# and 3/5 is astrocyte / oligodendrocyte.
# the 6 is your initial point with x_{1}=2, others=0. However, I am not sure whether least action path starting from this 6 can corresponds to ODE.

#identify cell types: 
Tangying_lap <- readMat('/Users/xqiu/Downloads/ansT10N150New.mat')

pData(nao_sim_cds)$cell_type <- rep(c('neuron', 'astrocyte', 'oligodendrocyte'), each = 400)
plot_cell_trajectory(nao_sim_cds, color_by = 'cell_type')

nao_sim_cds_reverse <- reverseEnbeddingCDS(nao_sim_cds)

valid_ids <- c(18, 24, 30) #astrocyte, neuron, oligodendrocyte
neuron_lap_df <- do.call(cbind, lapply(Tangying_lap$ans[1, 1, 1]$Path[valid_ids[3]], function(x) x[[1]])) #set the W matrix as the time steps

neuron_low_dim <- ginv(neuron_sim_cds@reducedDimW) %*% neuron_lap_df
pdf(paste(main_fig_dir, "fig2h_lap.pdf", sep = ''), height = 1.1, width = 1.1)
qplot(neuron_low_dim[1, ], neuron_low_dim[2, ], alpha = 0.5, color = I('gray'), size = 0.2) + 
  geom_point(aes(neuron_sim_cds@reducedDimK[1, ], neuron_sim_cds@reducedDimK[2, ], size = 3 * c(1:400) / 1000, color = I('red'), alpha = 0.5)) + 
  scale_size(range = range(3 * c(1:400) / 1000), limits = range(3 * c(1:400) / 1000)) + 
  scale_x_reverse() + scale_y_reverse() + nm_theme() + xlab('Component 1') + ylab('Component 2') + monocle:::monocle_theme_opts()
dev.off()

na_lap_df <- do.call(cbind, lapply(Tangying_lap$ans[1, 1, 1]$Path[valid_ids[2:3]], function(x) x[[1]])) #set the W matrix as the time steps
na_low_dim <- ginv(na_sim_cds@reducedDimW) %*% na_lap_df
pdf(paste(main_fig_dir, "fig2i_lap.pdf", sep = ''), height = 1.1, width = 1.1)
qplot(na_low_dim[1, ], na_low_dim[2, ], alpha = 0.5, color = I('gray')) + 
  geom_point(aes(na_sim_cds@reducedDimK[1, ], na_sim_cds@reducedDimK[2, ], size = 3 * rep(1:400, 2) / 1000, color = as.factor(rep(1:2, each = 400)), alpha = 0.5)) + 
  scale_size(range = c(0.003, 1.200), limits = c(0.003, 1.200)) + 
  scale_x_reverse() + scale_y_reverse() + nm_theme() + xlab('Component 1') + ylab('Component 2') + monocle:::monocle_theme_opts()
dev.off()

ao_lap_df <- do.call(cbind, lapply(Tangying_lap$ans[1, 1, 1]$Path[valid_ids[c(3, 5)]], function(x) x[[1]])) #set the W matrix as the time steps
ao_low_dim <- ginv(ao_sim_cds@reducedDimW) %*% ao_lap_df
pdf(paste(main_fig_dir, "fig2j_lap.pdf", sep = ''), height = 1.1, width = 1.1)
qplot(ao_low_dim[1, ], ao_low_dim[2, ], alpha = 0.5, color = I('gray'), size = 0.2) + 
  geom_point(aes(ao_sim_cds@reducedDimK[1, ], ao_sim_cds@reducedDimK[2, ], size = 3 *  rep(1:400, 2) / 1000, color = as.factor(rep(1:2, each = 400)), alpha = 0.5)) + 
  scale_size(range = range(3 *  rep(1:400, 2) / 1000), limits = range(3 *  rep(1:400, 2) / 1000)) + 
  scale_x_reverse() + scale_y_reverse() + nm_theme() + xlab('Component 1') + ylab('Component 2') + monocle:::monocle_theme_opts()
dev.off()

nao_lap_df <- do.call(cbind, lapply(Tangying_lap$ans[1, 1, 1]$Path[valid_ids[1]], function(x) x[[1]])) #set the W matrix as the time steps

nao_lap_df <- do.call(cbind, lapply(Tangying_lap$ans[1, 1, 1]$Path[valid_ids[3]], function(x) x[[1]])) #set the W matrix as the time steps
nao_low_dim <- ginv(nao_sim_cds_iter_1@reducedDimW) %*% nao_lap_df
pdf(paste(main_fig_dir, "fig2k_lap.pdf", sep = ''), height = 1.1, width = 1.1)
qplot(nao_sim_cds_iter_1@reducedDimS[1, ], nao_sim_cds_iter_1@reducedDimS[2, ], alpha = 0.5, color = I('gray'), size = 0.2) + 
  geom_point(aes(nao_sim_cds@reducedDimK[1, ], nao_sim_cds@reducedDimK[2, ], size = 3 *  rep(1:400, 3) / 1000, color = as.factor(rep(1:3, each = 400)), alpha = 0.5)) + 
  scale_size(range = range(3 *  rep(1:400, 3) / 1000), limits = range(3 *  rep(1:400, 3) / 1000)) + 
  scale_x_reverse() + scale_y_reverse() + nm_theme() + xlab('Component 1') + ylab('Component 2') + monocle:::monocle_theme_opts()
dev.off()

####################################################################################################################################################################################
# save data 
####################################################################################################################################################################################
save.image('./RData/fig3.RData')

####################################################################################################################################################################################
# add results with sample different cells to reconstruct a developmental trajectory (Note included ) 
####################################################################################################################################################################################
