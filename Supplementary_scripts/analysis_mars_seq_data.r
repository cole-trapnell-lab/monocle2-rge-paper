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

####################################################################################################################################################################################
#load all package
####################################################################################################################################################################################
library(dpt)
library(SLICER)
library(monocle)
library(mcclust)
library(destiny)
library(xacHelper)
library(reshape2)
library(plyr)
library(stringr)
library(GEOquery)
library(R.matlab)

load('./script_for_reproduce/MARSseq_analysis_tutorial.RData') # load the RData from Maren Büttner (https://github.com/theislab/scAnalysisTutorial), shared by Fabian 
source('~/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/scripts/function.R', echo=TRUE)
source('~/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/scripts/plotting.R', echo=TRUE)

main_fig_dir <- "./Figures/main_figures/"
SI_fig_dir <- "./Figures/supplementary_figures/"

#define a score for the erythroid or the GMP branch:
mep_genes <- c('Hba-a2', 'Car2', 'Cited4', 'Klf1')
cmp_genes <- c('Pf4', 'Apoe', 'Flt3', 'Cd74')
gmp_genes <- c('Mpo', 'Prg2', 'Prtn3', 'Ctsg')

Mar_seq_cols <- c("CMP" = "#44AF69", "DC" = "#F8333C", "erythroid" = "#FCAB10", "GMP" = "#2B9EB3")

################################################################################################################################################
# create the CDS which will be used for benchmark and learn the complicate trajectory for the Paul dataset         
################################################################################################################################################
MAP_cells_clusters <- readMat('/Users/xqiu/Dropbox (Personal)/bifurcation_path/simplePPT/data/MAP_cells_clusters.mat')
MAP_cells_clusters <- MAP_cells_clusters$MAP.cells.clusters
storage.mode(MAP_cells_clusters) <- 'numeric'
MAP_cells_clusters <- read.csv('/Users/xqiu/Downloads/MAP.csv', header = F)
row.names(MAP_cells_clusters) <- MAP_cells_clusters$V1

#use GEOquery to read related annotation information
GSE72857_parse <- getGEO(filename = './csv_data/GSE72857_family.soft.gz')

#valid_subset_GSE72857
valid_subset_GSE72857_exprs <- read.table('/Users/xqiu/Downloads/GSE72857_umitab.txt', header = T, row.names = 1)
writeMat('valid_subset_GSE72857_exprs', valid_subset_GSE72857_exprs = as.matrix(valid_subset_GSE72857_exprs))
# rownames(valid_subset_GSE72857_exprs) <- paste('g', 1:nrow(valid_subset_GSE72857_exprs), sep = '')
# colnames(valid_subset_GSE72857_exprs) <- MAP_cells_clusters$V1

#filtering cells to include only the ones which were assigned a cluster id: 
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

common_genes <- rownames(valid_subset_GSE72857_cds)[rownames(valid_subset_GSE72857_cds) %in% info.genes]
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = common_genes, row.names = common_genes))
valid_subset_GSE72857_cds <- newCellDataSet(as(as.matrix(data.info.genes[common_genes, ]), 'sparseMatrix'), 
                                            phenoData = pd, 
                                            featureData = fd,
                                            lowerDetectionLimit=1,
                                            expressionFamily=negbinomial.size())
valid_subset_GSE72857_cds <- estimateSizeFactors(valid_subset_GSE72857_cds)
valid_subset_GSE72857_cds <- estimateDispersions(valid_subset_GSE72857_cds)

pData(valid_subset_GSE72857_cds)$cell_type <- revalue(as.character(pData(valid_subset_GSE72857_cds)$cluster), 
                                                      c("1" = 'erythroid', "2" = 'erythroid', "3" = 'erythroid', "4" = 'erythroid', "5" = 'erythroid', "6" = 'erythroid', 
                                                        "7" = 'CMP', "8" = 'CMP', "9" = 'CMP', "10" = 'CMP',
                                                        "11" = 'DC', 
                                                        "12" = 'GMP', "13" = 'GMP', "14" = 'GMP', "15" = 'GMP', "16" = 'GMP', "17" = 'GMP', "18" = 'GMP', 
                                                        "19" = 'lymphoid'))

#remove all lymphoid cells because they don't belong to the meyloid lineage
valid_subset_GSE72857_cds <- valid_subset_GSE72857_cds[, pData(valid_subset_GSE72857_cds)$cell_type != 'lymphoid']

################################################################################################################################################
# calculate the lineage and stemness score
# create the SI9 
################################################################################################################################################
# mep score VS gmp lineage score: 
# Lin_ij = average[Er(G_j, i)] - average[Er(G_j^cont, i)]]
valid_genes = row.names(valid_subset_GSE72857_cds[rowSums(exprs(valid_subset_GSE72857_cds)) > 0, ])
total_expression = apply(exprs(valid_subset_GSE72857_cds)[valid_genes, ], 1, mean)
quartiles = cut(total_expression, breaks=quantile(total_expression, probs=seq(0,1, by=0.04))[-(1:2)], include.lowest=TRUE)
names(quartiles) = names(total_expression)

retrieve_control_genes <- function(marker_genes, quartiles){
  quartiles_tab <- table(as.character(quartiles[marker_genes]))
  control_ids <- c()
  for(quartile in names(quartiles_tab)) {
    control_ids <- c(control_ids, sample(which(quartiles == quartile), quartiles_tab[quartile] * 25)) #100
  }

  control_genes <- names(quartiles)[control_ids]
}

set.seed(2016)
mep_control_genes <- retrieve_control_genes(mep_genes, quartiles)
cmp_control_genes <- retrieve_control_genes(cmp_genes, quartiles)
gmp_control_genes <- retrieve_control_genes(gmp_genes, quartiles)

#
lineage_score <- function(cds, marker_genes, marker_control_genes){
  esApply(cds[c(marker_genes, marker_control_genes), ], 2, function(x) {
    print(length(x))
    marker_genes_mean_exp <- mean(x[1:length(marker_genes)])
    marker_control_genes_mean_exp <- mean(x[(length(marker_genes) + 1):(length(x))])
    
    marker_genes_mean_exp - marker_control_genes_mean_exp
  })
}

mep_lineage_score <- lineage_score(valid_subset_GSE72857_cds, mep_genes, mep_control_genes)
gmp_lineage_score <- lineage_score(valid_subset_GSE72857_cds, gmp_genes, gmp_control_genes)

lineage_score_df <- data.frame(mep_lineage_score = mep_lineage_score, gmp_lineage_score = gmp_lineage_score)

#stemness score: 
# Stem(i) = average[Er(G_{stem})] - average[Er(G_{stem}^{cont})] - LIN(i)

cmp_lineage_score_ori <- lineage_score(valid_subset_GSE72857_cds, cmp_genes, cmp_control_genes)
cmp_lineage_score <- cmp_lineage_score_ori - apply(lineage_score_df, 1, max)

#plot the result: 
score_df <- data.frame(stemness = cmp_lineage_score, lineage_score = apply(lineage_score_df, 1, function(x) {
  x <- c(-x[1], x[2])
  x[which.max(abs(x))]
}))

score_df$lineage <- "GMP"
score_df$lineage[score_df$lineage_score < 0] <- 'MEP'
score_df$lineage[score_df$stemness > 0] <- 'CMP'
row.names(score_df) <- names(mep_lineage_score)

#use just the distances to the point (0, 0)
score_df$naive_pseudotime <- apply(score_df, 1, function(x) sqrt(sum(as.numeric(x[c(1, 2)])^2)))
score_df$naive_pseudotime[score_df$lineage == 'CMP'] <- score_df$naive_pseudotime[score_df$stemness == max(score_df$stemness)] - 
  score_df$naive_pseudotime[score_df$lineage == 'CMP']
score_df$naive_pseudotime[score_df$lineage != 'CMP'] <- score_df$naive_pseudotime[score_df$stemness == max(score_df$stemness)] + 
  score_df$naive_pseudotime[score_df$lineage != 'CMP']

score_df$cell_type <- pData(valid_subset_GSE72857_cds)$cell_type
  
pdf(paste(main_fig_dir, 'marseq_score.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(lineage_score, stemness, color = cell_type, size = naive_pseudotime, data = score_df) + nm_theme() + scale_size(range = c(0.2, 1.2)) + scale_color_manual(values = Mar_seq_cols)
dev.off()

pdf(paste(main_fig_dir, 'marseq_score_helper.pdf', sep = ''))
qplot(lineage_score, stemness, color = lineage, size = naive_pseudotime, data = score_df) + scale_size(range = c(0.2, 1.2))
dev.off()

table(score_df$lineage)

#order by monocle2
valid_subset_GSE72857_cds <- setOrderingFilter(valid_subset_GSE72857_cds, ordering_genes = row.names(data.info.genes))
valid_subset_GSE72857_cds <- reduceDimension(valid_subset_GSE72857_cds, norm_method = 'log', verbose = T) #
valid_subset_GSE72857_cds <- orderCells(valid_subset_GSE72857_cds)

plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'Pseudotime')
plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'cell_type')
plot_cell_trajectory(valid_subset_GSE72857_cds)
valid_subset_GSE72857_cds <- orderCells(valid_subset_GSE72857_cds, root_state = 2)

pdf(paste(main_fig_dir, 'SI9a.pdf', sep = ''), height = 1.5, width = 1.5)
plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'cell_type', cell_size = 0.5, show_branch_points = F) + nm_theme() + 
  scale_color_manual(values = Mar_seq_cols)#, cell_size = pData(valid_subset_GSE72857_cds)$Pseudotime
dev.off()

# pdf(paste(main_fig_dir, 'monocle2_marseq_tree.pdf', sep = ''), height = 1.5, width = 1.5)
# plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'cell_type', cell_size = 0.5, show_branch_points = F) + nm_theme() 
# dev.off()

pdf(paste(main_fig_dir, 'SI9a_helper.pdf', sep = ''))
plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'cell_type', cell_size = 0.5)#, cell_size = pData(valid_subset_GSE72857_cds)$Pseudotime
dev.off()

#compare the kendall'tau and adjusted rand index
calClusteringMetrics(score_df$lineage, pData(valid_subset_GSE72857_cds[, row.names(score_df)])$State)

#################################################################################################################################################################################################################################################
#save the data: 
save.image('./RData/analysis_score_ordering_MAR_seq.RData')

#################################################################################################################################################################################################################################################
# run dpt, wishbone, slicer algorithms and create all other panels for the SI9a
# sample cells for monocle 1 and slicer 
#################################################################################################################################################################################################################################################
set.seed(2017)
sample_cells <- colnames(valid_subset_GSE72857_cds)[sample(1:ncol(valid_subset_GSE72857_cds), 300)]
# colnames(valid_subset_GSE72857_cds)[c(sample(which(pData(valid_subset_GSE72857_cds)$cell_type == 'CMP'), 100),
#                                                       sample(which(pData(valid_subset_GSE72857_cds)$cell_type == 'GMP'), 100),
#                                                       sample(which(pData(valid_subset_GSE72857_cds)$cell_type == 'erythroid'), 100))]
ordering_genes <- row.names(subset(fData(valid_subset_GSE72857_cds), use_for_ordering == T))

dpt_res <- run_new_dpt(valid_subset_GSE72857_cds[ordering_genes, ], normalize = T)
slicer_res <- run_slicer(valid_subset_GSE72857_cds[ordering_genes, sample_cells])
monocle1_cds <- reduceDimension(valid_subset_GSE72857_cds[ordering_genes, sample_cells], reduction_method = 'ICA', verbose = T)
monocle1_cds <- orderCells(monocle1_cds, num_paths = 2)
plot_cell_trajectory(monocle1_cds, color_by = 'Pseudotime')
plot_cell_trajectory(monocle1_cds, color_by = 'cell_type')
plot_cell_trajectory(monocle1_cds)

monocle1_cds <- orderCells(monocle1_cds, root_state = 3)

# dpt_df <- data.frame(x = as.numeric(eigenvectors(diff.plot)[, 1]), y = as.numeric(eigenvectors(diff.plot)[, 2]), branch = as.character(color.branch[branching]), cell_type = pData(valid_subset_GSE72857_cds)$cell_type)
dpt_df <- data.frame(x = dpt_res$dm$DC1, y = dpt_res$dm$DC2, dpt = dpt_res$pt, branch = dpt_res$branch[, 1], cell_type = pData(valid_subset_GSE72857_cds)$cell_type)

pdf(paste(main_fig_dir, 'SI9a_dpt_tree.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(x, y, color = cell_type, data = dpt_df, size = 0.25) + nm_theme() + scale_size(range = c(0.25, 0.25)) + scale_color_manual(values = Mar_seq_cols)
dev.off()

# plot(eigenvectors(dpt_res$dm)[, 1:2], col=color.branch[branching], pch=20, xlab='', ylab='', 
#      axes=FALSE)
# par(cex.lab=1.8)
# axis(side=1, at=c(-0.05, 0.1), labels=NA,  col.ticks = 'white')
# axis(side=2, at=c(-0.11, 0.11), labels=NA, col.ticks = 'white')

slicer_df <- data.frame(LLE1 = slicer_res$traj_lle[, 1], LLE2 = slicer_res$traj_lle[, 2], 
                        pseudotime = slicer_res$order_df$cells_ordered, branch = slicer_res$order_df$branches, 
                        cell_type = pData(valid_subset_GSE72857_cds[, sample_cells])$cell_type)
pdf(paste(main_fig_dir, 'SI9a_slicer_tree.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(LLE1, LLE2, color = cell_type, data = slicer_df) + nm_theme() + scale_size(range = c(0.25, 0.25)) + scale_color_manual(values = Mar_seq_cols)
dev.off()

pdf(paste(main_fig_dir, 'SI9a_monocle1_tree.pdf', sep = ''), height = 1.5, width = 1.5)
plot_cell_trajectory(monocle1_cds, color_by = 'cell_type', cell_size = 0.5) + scale_color_manual(values = Mar_seq_cols) + nm_theme()#, cell_size = pData(valid_subset_GSE72857_cds)$Pseudotime
dev.off()

pdf(paste(main_fig_dir, 'mSI9a_monocle1_tree_helper.pdf', sep = ''))
plot_cell_trajectory(monocle1_cds, color_by = 'cell_type', cell_size = 0.5) + scale_color_manual(values = Mar_seq_cols) #, cell_size = pData(valid_subset_GSE72857_cds)$Pseudotime
dev.off()

data <- t(log(exprs(valid_subset_GSE72857_cds) + 1))
write.csv(file = paste('./csv_data/Wishbone_test_data/', 'marseq_data', ".txt", sep = ''), as.matrix(data), quote = F, row.names = T)

row.names(subset(pData(valid_subset_GSE72857_cds), Pseudotime == 0))
# [1] "W31587"

wishbone_res <- read.table('./wishbone/MAR_seq_fractioin_wishbone_df_fig_4.txt', header = T, sep = '\t')
qplot(dm1, dm2, data = wishbone_res)
qplot(tSNE1, tSNE2, data = wishbone_res, color = as.character(branch), size = 0.5) 

pdf(paste(main_fig_dir, 'SI9a_wishbone_tree.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(tSNE1, tSNE2, data = wishbone_res, color = pData(valid_subset_GSE72857_cds)$cell_type, size = 0.5) + nm_theme() + scale_size(range = c(0.5, 0.5)) + scale_color_manual(values = Mar_seq_cols)
# plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'cell_type', cell_size = 0.5, show_branch_points = F) + nm_theme()#, cell_size = pData(valid_subset_GSE72857_cds)$Pseudotime
dev.off()

########################################################################################################################################################
#run the benchmark pipeline: 
# run benchmark pipeline with DPT, Monocle 1, Slicer and Wishbone (note that Monocle 1 takes a lot time to run)
########################################################################################################################################################
absolute_cds <- valid_subset_GSE72857_cds
auto_param_selection <- T
repeat_downsampling_num <- 25
benchmark_type <- 'MAR_seq_data' #lung_data
pData(absolute_cds)$Time <- pData(absolute_cds)$cell_type
root_state <- 2
AT1_state <- 1
AT2_state <- 3
source('./scripts/DDRTree_robustness_analysis_genearlize.R', echo = T)

# be careful to run this file, it takes a long time! 
source('./script_for_reproduce/ICA_robustness_analysis_genearlize.R', echo = T) 

#use monocle2 for both monocle 1 and wishbone for downstream analysis: 
ICA_sampling_res_df <- sampling_res_df
ICA_valid_cell_sampling_res_df <- valid_cell_sampling_res_df
source('./scripts/robustness_dpt_slicer_genearlize.R', echo = T)

downsampling_state <- 2 #check this 
source('./scripts/progressive_branch_downsampling_generalize.R', echo = T)
########################################################################################################################################################

########################################################################################################################################################
save.image('./RData/fig4_marseq.RData')
########################################################################################################################################################
