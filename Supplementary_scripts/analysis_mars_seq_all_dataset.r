#MAR-seq ko vs wt cells
#tree for wt, all wt and corresponding branched heatmap for the blood data 

rm(list = ls())
########################################################################################################################################################
#load all the necessary libraries 
library(monocle)
library(plyr)
library(dplyr)
library(destiny)
library(xacHelper)
library(grid)

source('./scripts/function.R')

main_fig_dir <- "./Figures/main_figures/"
SI_fig_dir <- "./Figures/supplementary_figures/"

load('./RData/analysis_score_ordering_MAR_seq.RData')
# load('./RData/analysis_score_ordering_MAR_seq.RData')
# load('./RData/dpFeature_for_all_datasets.RData')

########################################################################################################################################################
# prepare cds for all datasets: 
########################################################################################################################################################

all_valid_cells <- colnames(valid_subset_GSE72857_exprs)[apply(valid_subset_GSE72857_exprs, 2, sum) > 501] #remove cells with very low capture rate 
all_valid_design_mat <- design_mat[all_valid_cells, ]
pd <- new("AnnotatedDataFrame", data = all_valid_design_mat)
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = row.names(valid_subset_GSE72857_exprs), row.names = row.names(valid_subset_GSE72857_exprs)))

all_GSE72857_cds <- newCellDataSet(as(as.matrix(valid_subset_GSE72857_exprs[, all_valid_cells]), 'sparseMatrix'), 
                                   phenoData = pd, 
                                   featureData = fd,
                                   lowerDetectionLimit=1,
                                   expressionFamily=negbinomial.size())

rm(valid_subset_GSE72857_exprs) #remove the gigantic data 

all_GSE72857_cds <- estimateSizeFactors(all_GSE72857_cds)
all_GSE72857_cds <- estimateDispersions(all_GSE72857_cds)

pData(all_GSE72857_cds)$cluster <- as.character(pData(all_GSE72857_cds)$cluster)
pData(all_GSE72857_cds)$cell_type <- revalue(pData(all_GSE72857_cds)$cluster, 
                                             c("1" = 'erythroid', "2" = 'erythroid', "3" = 'erythroid', "4" = 'erythroid', "5" = 'erythroid', "6" = 'erythroid', 
                                               "7" = 'CMP', "8" = 'CMP', "9" = 'CMP', "10" = 'CMP',
                                               "11" = 'DC', 
                                               "12" = 'GMP', "13" = 'GMP', "14" = 'GMP', "15" = 'GMP', "16" = 'GMP', "17" = 'GMP', "18" = 'GMP', 
                                               "19" = 'lymphoid'))

########################################################################################################################################################
# remove lymphoid cells and redo the ordering 
########################################################################################################################################################
all_GSE72857_cds <- all_GSE72857_cds[, setdiff(1:ncol(all_GSE72857_cds), which(pData(all_GSE72857_cds)$cell_type == 'lymphoid'))]

MAP_ordering_genes <- row.names(subset(fData(valid_subset_GSE72857_cds), use_for_ordering == T)) #row.names(MAP_clustering_DEG_genes_subset)[order(MAP_clustering_DEG_genes_subset$qval)][1:1000] #1971
all_GSE72857_cds <- setOrderingFilter(all_GSE72857_cds, ordering_genes = MAP_ordering_genes)
all_GSE72857_cds <- reduceDimension(all_GSE72857_cds, norm_method = 'log', verbose = T)
all_GSE72857_cds <- orderCells(all_GSE72857_cds)

plot_cell_trajectory(all_GSE72857_cds, color_by = 'cell_type') + facet_wrap(~cell_type)
plot_cell_trajectory(all_GSE72857_cds, color_by = 'Batch_desc') + facet_wrap(~Batch_desc)
plot_cell_trajectory(all_GSE72857_cds)
plot_cell_trajectory(all_GSE72857_cds, color_by = 'Pseudotime')

all_GSE72857_cds <- orderCells(all_GSE72857_cds, root_state = 3)

#run dpt: 
all_MAP_dpt_res <- run_new_dpt(all_GSE72857_cds, normalize = F)
qplot(all_MAP_dpt_res$dm$DC1, all_MAP_dpt_res$dm$DC2, color = pData(all_GSE72857_cds)$cell_type)
qplot(all_MAP_dpt_res$dm$DC2, all_MAP_dpt_res$dm$DC3, color = pData(all_GSE72857_cds)$cell_type)

########################################################################################################################################################
# only look at the wild-type cells: 
########################################################################################################################################################

table(pData(all_GSE72857_cds)$Batch_desc)
all_wt_valid_cells <- row.names(subset(pData(all_GSE72857_cds), !(Batch_desc %in% c('Cebpa KO', 'Cebpe KO', 'Cebpa control', 'Cebpe control'))))
all_wt_GSE72857_cds <- setOrderingFilter(all_GSE72857_cds[, all_wt_valid_cells], ordering_genes = MAP_ordering_genes)
all_wt_GSE72857_cds <- reduceDimension(all_wt_GSE72857_cds, norm_method = 'log', verbose = T)
all_wt_GSE72857_cds <- orderCells(all_wt_GSE72857_cds)

plot_cell_trajectory(all_wt_GSE72857_cds)
plot_cell_trajectory(all_wt_GSE72857_cds, color_by = 'cell_type') + facet_wrap(~cell_type) + stat_density_2d()
plot_cell_trajectory(all_wt_GSE72857_cds, color_by = 'Pseudotime') 
all_wt_GSE72857_cds <- orderCells(all_wt_GSE72857_cds, root_state = 3)

pdf(paste(SI_fig_dir, "all_WT_cells_mar_seq_data.pdf", sep = ''), height = 3, width = 3)
plot_cell_trajectory(all_wt_GSE72857_cds, color_by = 'cell_type', theta = 30, show_branch_points = F) + scale_y_reverse() + facet_wrap(~Batch_desc, nrow = 2) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="none", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25))
dev.off()

pdf(paste(SI_fig_dir, "all_WT_cells_mar_seq_data_helper.pdf", sep = ''))
plot_cell_trajectory(all_wt_GSE72857_cds, color_by = 'cell_type', theta = 30) + scale_y_reverse() + facet_wrap(~Batch_desc, nrow = 2) + scale_size(range = c(0.2, 0.2)) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25))
dev.off()

# #only look at the GMP cells
# table(pData(valid_subset_GSE72857_cds)$cell_type)
# GMP_valid_cells <- row.names(subset(pData(valid_subset_GSE72857_cds), (cell_type %in% c('GMP'))))
# GMP_valid_subset_GSE72857_cds <- setOrderingFilter(valid_subset_GSE72857_cds[, myeloid_valid_cells], ordering_genes = blood_beam_genes$V1)
# GMP_valid_subset_GSE72857_cds <- reduceDimension(GMP_valid_subset_GSE72857_cds, norm_method = 'log', verbose = T)
# GMP_valid_subset_GSE72857_cds <- orderCells(GMP_valid_subset_GSE72857_cds)
# 
# #we may detect the monocyte / granuocyte branch
# plot_cell_trajectory(GMP_valid_subset_GSE72857_cds)
# plot_genes_branched_pseudotime(GMP_valid_subset_GSE72857_cds[c('Gfi1', 'Irf8'), ])
#  

########################################################################################################################################################
# perform some DEG test between the KO cells and the WT cells on the same location of the manifold
# create Figure SI19G
########################################################################################################################################################

plot_cell_trajectory(all_GSE72857_cds, color_by = 'cell_type') + facet_wrap(~cell_type)
plot_cell_trajectory(all_GSE72857_cds, color_by = 'Batch_desc') + facet_wrap(~Batch_desc)
plot_cell_trajectory(all_GSE72857_cds)
plot_cell_trajectory(all_GSE72857_cds, color_by = 'Pseudotime')

# + scale_color_manual(values = type_cols, name = "Type")
pdf(paste(SI_fig_dir, 'SI19G_old.pdf', sep = ''), height = 3, width = 4)
plot_cell_trajectory(all_GSE72857_cds, color_by = 'Batch_desc', show_branch_points = F, cell_size = 0.5, cell_link_size = 0.3, theta = 30) + scale_y_reverse() + facet_wrap(~Batch_desc, nrow = 2) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="none", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) #theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

# add this to the figure 6
pdf(paste(main_fig_dir, 'SI19G.pdf', sep = ''), height = 1.7, width = 5.5)
plot_cell_trajectory(all_GSE72857_cds[, pData(all_GSE72857_cds)$Batch_desc %in% c('Cebpa control', 'Cebpa KO', 'Cebpe control', 'Cebpe KO')], 
                     color_by = 'Batch_desc', show_branch_points = F, cell_size = 0.5, cell_link_size = 0.3, theta = 30) + scale_y_reverse() + facet_wrap(~Batch_desc, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="none", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) #theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

plot_cell_trajectory(all_GSE72857_cds, color_by = 'Batch_desc', show_branch_points = F, cell_size = 0.5, cell_link_size = 0.3) 

#######################################################################################################################################################################################################
# Cebpa leaking cells (Cebpa leaking cells)
#######################################################################################################################################################################################################
# leaking cells: 
range(all_GSE72857_cds[, pData(all_GSE72857_cds)$State == 2 & pData(all_GSE72857_cds)$Batch_desc %in% c('Cebpa KO')]$Pseudotime)
Cebpa_ko_leaking_cells <- all_GSE72857_cds[, pData(all_GSE72857_cds)$State == 2 & pData(all_GSE72857_cds)$Batch_desc %in% c('Cebpa KO')]
#[1]  3.662649 20.724515

Cebpa_ko_mainfold_cells <- all_GSE72857_cds[, pData(all_GSE72857_cds)$Pseudotime > 3 & pData(all_GSE72857_cds)$Pseudotime < 21 & pData(all_GSE72857_cds)$State == 2 & !(pData(all_GSE72857_cds)$Batch_desc %in% c('Cebpa KO'))]
pData(Cebpa_ko_mainfold_cells)$figSI4_Cebpa_ko <- 'WildType'
pData(Cebpa_ko_mainfold_cells)$figSI4_Cebpa_ko[pData(Cebpa_ko_mainfold_cells)$Batch_desc %in% c('Cebpa KO')] <- "leaking_cells"
all_GSE72857_cds_genotype_figSI41_beam_genes <- differentialGeneTest(Cebpa_ko_mainfold_cells, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3) + figSI4_Cebpa_ko", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)")
all_GSE72857_cds_genotype_figSI41_group_deg <- differentialGeneTest(Cebpa_ko_mainfold_cells, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~figSI4_Cebpa_ko", reducedModelFormulaStr = "~1")

###########################################################################################################################################################################
# test on the supposed branch: 
###########################################################################################################################################################################
range(all_GSE72857_cds[, pData(all_GSE72857_cds)$State == 1 & pData(all_GSE72857_cds)$Batch_desc %in% c('Cebpa KO')]$Pseudotime)
# [1]  3.662649 20.724515
Cebpa_ko_leaking_cells <- all_GSE72857_cds[, pData(all_GSE72857_cds)$State == 1 & pData(all_GSE72857_cds)$Batch_desc %in% c('Cebpa KO')]

Cebpa_ko_mainfold_cells.1 <- all_GSE72857_cds[, pData(all_GSE72857_cds)$Pseudotime > 3 & pData(all_GSE72857_cds)$Pseudotime < 21 & pData(all_GSE72857_cds)$State == 1 & !(pData(all_GSE72857_cds)$Batch_desc %in% c('Cebpe KO'))]
pData(Cebpa_ko_mainfold_cells.1)$figSI4_Cebpa_ko.1 <- 'WildType'
pData(Cebpa_ko_mainfold_cells.1)$figSI4_Cebpa_ko.1[pData(Cebpa_ko_mainfold_cells.1)$Batch_desc %in% c('Cebpa KO')] <- "Cebpa KO"
# all_GSE72857_cds_genotype_figSI4_Cebpa_ko.1_beam_genes.1 <- differentialGeneTest(Cebpa_ko_mainfold_cells.1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3) + figSI4_Cebpa_ko.1", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)")
all_GSE72857_cds_genotype_figSI4_Cebpa_ko.1_group_deg.1 <- differentialGeneTest(Cebpa_ko_mainfold_cells.1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~figSI4_Cebpa_ko.1", reducedModelFormulaStr = "~1")

###########################################################################################################################################################################
# test KO cells from two branches: 
###########################################################################################################################################################################
Cebpa_ko_leaking_cells.2 <- all_GSE72857_cds[, pData(all_GSE72857_cds)$State %in% c(1, 2) & pData(all_GSE72857_cds)$Batch_desc %in% c('Cebpa KO')]

all_GSE72857_cds_genotype_fig6_Cebpa_ko_group_deg.2 <- differentialGeneTest(Cebpa_ko_leaking_cells.2, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~State", reducedModelFormulaStr = "~1")

########################################################################################################################################################################################################
# test Cebpe KO cells
########################################################################################################################################################################################################

# leaking cells: 
range(all_GSE72857_cds[, pData(all_GSE72857_cds)$State == 1 & pData(all_GSE72857_cds)$Batch_desc %in% c('Cebpe KO')]$Pseudotime)
Cebpe_ko_leaking_cells <- all_GSE72857_cds[, pData(all_GSE72857_cds)$State == 1 & pData(all_GSE72857_cds)$Batch_desc %in% c('Cebpe KO')]
#[1]  3.500532 20.479481

Cebpe_ko_mainfold_cells <- all_GSE72857_cds[, pData(all_GSE72857_cds)$Pseudotime > 3 & pData(all_GSE72857_cds)$Pseudotime < 21 & pData(all_GSE72857_cds)$State == 1 & !(pData(all_GSE72857_cds)$Batch_desc %in% c('Cebpa KO'))]
pData(Cebpe_ko_mainfold_cells)$figSI4_Cebpe_ko <- 'WildType'
pData(Cebpe_ko_mainfold_cells)$figSI4_Cebpe_ko[pData(Cebpe_ko_mainfold_cells)$Batch_desc %in% c('Cebpe KO')] <- "Cebpe KO"
# all_GSE72857_cds_genotype_figSI41_beam_genes <- differentialGeneTest(Cebpe_ko_mainfold_cells, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3) + figSI4_Cebpe_ko", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)")
all_GSE72857_cds_genotype_figSI41_group_deg_cebpe <- differentialGeneTest(Cebpe_ko_mainfold_cells, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~figSI4_Cebpe_ko", reducedModelFormulaStr = "~1")

###########################################################################################################################################################################
# test on the supposed branch: 
###########################################################################################################################################################################
range(all_GSE72857_cds[, pData(all_GSE72857_cds)$State == 2 & pData(all_GSE72857_cds)$Batch_desc %in% c('Cebpe KO')]$Pseudotime)
# [1]  4.459144 13.817094
Cebpe_ko_leaking_cells <- all_GSE72857_cds[, pData(all_GSE72857_cds)$State == 2 & pData(all_GSE72857_cds)$Batch_desc %in% c('Cebpe KO')]

Cebpe_ko_mainfold_cells.1 <- all_GSE72857_cds[, pData(all_GSE72857_cds)$Pseudotime > 4 & pData(all_GSE72857_cds)$Pseudotime < 14 & pData(all_GSE72857_cds)$State == 2 & !(pData(all_GSE72857_cds)$Batch_desc %in% c('Cebpa KO'))]
pData(Cebpe_ko_mainfold_cells.1)$figSI4_Cebpe_ko.1 <- 'WildType'
pData(Cebpe_ko_mainfold_cells.1)$figSI4_Cebpe_ko.1[pData(Cebpe_ko_mainfold_cells.1)$Batch_desc %in% c('Cebpe KO')] <- "Cebpe KO"
# all_GSE72857_cds_genotype_figSI41_beam_genes.1 <- differentialGeneTest(Cebpe_ko_mainfold_cells.1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3) + figSI4_Cebpe_ko.1", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)")
all_GSE72857_cds_genotype_figSI41_group_deg.1 <- differentialGeneTest(Cebpe_ko_mainfold_cells.1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~figSI4_Cebpe_ko.1", reducedModelFormulaStr = "~1")

###########################################################################################################################################################################
# test KO cells from two branches: 
###########################################################################################################################################################################
Cebpe_ko_leaking_cells.2 <- all_GSE72857_cds[, pData(all_GSE72857_cds)$State %in% c(1, 2) & pData(all_GSE72857_cds)$Batch_desc %in% c('Cebpe KO')]

all_GSE72857_cds_genotype_fig6_Cebpe_ko_group_deg.2 <- differentialGeneTest(Cebpe_ko_leaking_cells.2, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~State", reducedModelFormulaStr = "~1")

########################################################################################################################################################
# BEAM test on the WT cells and branchTimePoint analysis 
########################################################################################################################################################
valid_subset_GSE72857_cds_beam_genes <- BEAM(valid_subset_GSE72857_cds, branch_point = 1, verbose = T, cores = detectCores() - 2)
valid_subset_GSE72857_cds_beam_genes_proj_dup <- BEAM(valid_subset_GSE72857_cds, branch_point = 1, verbose = T, progenitor_method = 'duplicate', cores = detectCores() - 2)

all_GSE72857_cds <- orderCells(all_GSE72857_cds)
all_GSE72857_cds <- orderCells(all_GSE72857_cds, root_state = 3)

all_GSE72857_cds_beam_genes <- BEAM(all_GSE72857_cds, branch_point = 1, verbose = T, cores = detectCores() - 2)
all_GSE72857_cds_beam_genes_proj_dup <- BEAM(all_GSE72857_cds, branch_point = 1, verbose = T, progenitor_method = 'duplicate', cores = detectCores() - 2)

plot_cell_trajectory(valid_subset_GSE72857_cds)
plot_cell_trajectory(all_GSE72857_cds)
branch_states <- c(1, 2) 
branch_states_all <- c(1,2) ####check this

valid_subset_GSE72857_cds_all_gene_ILRs_list <- calILRs(cds = valid_subset_GSE72857_cds[, ],  branch_point = 1, stretch = T, cores = 1, #beam_genes
                              trend_formula = "~sm.ns(Pseudotime, df = 3) * Branch", ILRs_limit = 3, relative_expr = T, label_by_short_name = F,
                              useVST = F, round_exprs = FALSE, output_type = "all", file = "all_gene_ILRs_list", return_all = T)
all_GSE72857_cds_all_gene_ILRs_list <- calILRs(cds = all_GSE72857_cds[, ], branch_point = 1, stretch = T, cores = 1, #beam_genes
                                                        trend_formula = "~sm.ns(Pseudotime, df = 3) * Branch", ILRs_limit = 3, relative_expr = T, label_by_short_name = F,
                                                        useVST = F, round_exprs = FALSE, output_type = "all", file = "all_gene_ILRs_list", return_all = T)
#get the corresponding branch time for the branchpoint of the MST tree:
beam_genes <- row.names(subset(valid_subset_GSE72857_cds_beam_genes, qval < 0.1)) #valid_subset_GSE72857_cds_beam_genes_proj_dup all_GSE72857_cds_beam_genes all_GSE72857_cds_beam_genes_proj_dup

#cds_full <- buildBranchCellDataSet(valid_subset_GSE72857_cds[unique(beam_genes), ], progenitor_method = 'duplicate', branch_states = branch_states, stretch = T, branch_labels = NULL)
#cds_fig1b <- buildBranchCellDataSet(valid_subset_GSE72857_cds[unique(beam_genes), ], progenitor_method = 'duplicate', branch_states = branch_states, stretch = T, branch_labels = NULL)

#check the longest pseudotime for common progenitors:
#range(pData(cds_full[, pData(cds_full)$State == 1])$Pseudotime)

# abs_bifurcation_time <- abs(detectBifurcationPoint(valid_subset_GSE72857_cds_all_gene_ILRs_list$norm_str_logfc_df[, 60:100], ILRs_threshold = 0.1, return_cross_point = F)) + 60#

#don't subset cds
abs_bifurcation_time <- abs(detectBifurcationPoint(valid_subset_GSE72857_cds_all_gene_ILRs_list$norm_str_logfc_df[, 0:100], ILRs_threshold = 0.1, return_cross_point = F))#

########################################################################################################################################################
# Cebpa and Irf8's Chip-seq targets:  (no Cebpe Chip-seq data)
########################################################################################################################################################
Cebpa_target_gene_names <- as.character(read.table('./generate_clustering/tf_targets/Cebpa_23403033.txt')$V1)

Irf8_21731497_target_gene_names <- as.character(read.table('./tf_targets/Irf8_21731497.txt')$V1) # Irf8_22096565.txt_targets.txt
Irf8_22096565_target_gene_names <- as.character(read.table('./generate_clustering/tf_targets/Irf8_22096565.txt')$V1) #

########################################################################################################################################################
# test on the supposed branch for Cebpa KO (just check the changes of Gfi1' targets on the correctly positioned cells)
########################################################################################################################################################
cebpa_ko_erythroid_branch_genes <- intersect(beam_genes, row.names(subset(all_GSE72857_cds_genotype_figSI4_Cebpa_ko.1_group_deg.1, qval < 0.1)))
# cebpa_ko_erythroid_branch_genes <- row.names(subset(all_GSE72857_cds_genotype_figSI4_Cebpa_ko.1_group_deg.1, qval < 0.1))

mean_exprs_cebpa_ko_erythroid_branch <- rowMeans(exprs(Cebpa_ko_mainfold_cells.1[cebpa_ko_erythroid_branch_genes, ]))
figSI4_Cebpa_ko.1 <- pData(Cebpa_ko_mainfold_cells.1)$figSI4_Cebpa_ko.1

lfc_cebpa_ko_erythroid_branch <- log2(rowMeans(exprs(Cebpa_ko_mainfold_cells.1[cebpa_ko_erythroid_branch_genes, figSI4_Cebpa_ko.1 == "Cebpa KO"]))/ 
                                 rowMeans(exprs(Cebpa_ko_mainfold_cells.1[cebpa_ko_erythroid_branch_genes, figSI4_Cebpa_ko.1 == "WildType"])))

lfc_cebpa_ko_erythroid_branch_Cebpa_targets <- intersect(toupper(cebpa_ko_erythroid_branch_genes), Cebpa_target_gene_names)
# lfc_double_ko_monocyte_Gfi1_targets <- intersect(cebpa_ko_erythroid_branch_genes, Irf8_21731497_target_gene_names)
cebpa_ko_erythroid_Type <- rep(NA, length(cebpa_ko_erythroid_branch_genes))

# double_ko_Type[double_ko_monocyte_branch_genes %in% lfc_double_ko_monocyte_Irf8_targets] <- 'Irf8'
cebpa_ko_erythroid_Type[toupper(cebpa_ko_erythroid_branch_genes) %in% lfc_cebpa_ko_erythroid_branch_Cebpa_targets] <- 'Cebpa'
qplot(mean_exprs_cebpa_ko_erythroid_branch, lfc_cebpa_ko_erythroid_branch, color = cebpa_ko_erythroid_Type)
# qplot(mean_exprs_cebpa_ko_erythroid_branch, lfc_cebpa_ko_erythroid_branch, color = cebpa_ko_erythroid_Type, facets = '~cebpa_ko_erythroid_Type') + xlab('Mean expression') + ylab('Fold changes')
# 
# qplot(mean_exprs_cebpa_ko_erythroid_branch, lfc_cebpa_ko_erythroid_branch, color = cebpa_ko_erythroid_Type) + 
#   xlab('Mean expression') + ylab('Fold changes') + facet_wrap(~cebpa_ko_erythroid_Type)

########################################################################################################################################################
# check whether or not the leaking cells for Cebpe have overexpression of Cebpa targets 
########################################################################################################################################################
cebpe_ko_erythroid_branch_genes <- intersect(beam_genes, row.names(subset(all_GSE72857_cds_genotype_figSI41_group_deg_cebpe, qval < 0.1)))
# cebpe_ko_erythroid_branch_genes <- row.names(subset(all_GSE72857_cds_genotype_figSI41_group_deg_cebpe, qval < 0.1))

mean_exprs_cebpa_ko_erythroid_branch <- rowMeans(exprs(Cebpe_ko_mainfold_cells[cebpe_ko_erythroid_branch_genes, ]))
figSI4_Cebpe_ko <- pData(Cebpe_ko_mainfold_cells)$figSI4_Cebpe_ko

lfc_cebpe_ko_erythroid_branch <- log2(rowMeans(exprs(Cebpe_ko_mainfold_cells[cebpe_ko_erythroid_branch_genes, figSI4_Cebpe_ko == "leaking_cells"]))/ 
                                        rowMeans(exprs(Cebpe_ko_mainfold_cells[cebpe_ko_erythroid_branch_genes, figSI4_Cebpe_ko == "WildType"])))

lfc_cebpe_ko_erythroid_branch_Cebpa_targets <- intersect(toupper(cebpe_ko_erythroid_branch_genes), Cebpa_target_gene_names)
# lfc_double_ko_monocyte_Gfi1_targets <- intersect(cebpa_ko_erythroid_branch_genes, Irf8_21731497_target_gene_names)
cebpe_ko_erythroid_Type <- rep(NA, length(cebpe_ko_erythroid_branch_genes))

# double_ko_Type[double_ko_monocyte_branch_genes %in% lfc_double_ko_monocyte_Irf8_targets] <- 'Irf8'
cebpe_ko_erythroid_Type[toupper(cebpe_ko_erythroid_branch_genes) %in% lfc_cebpe_ko_erythroid_branch_Cebpa_targets] <- 'Cebpa'

#annotate the genes: 
cebpe_ko_erythroid_df <- data.frame(Type = cebpe_ko_erythroid_Type[!is.na(cebpe_ko_erythroid_Type)], 
                                    mean_expression = mean_exprs_cebpa_ko_erythroid_branch[!is.na(cebpe_ko_erythroid_Type)],
                                    lfc = lfc_cebpe_ko_erythroid_branch[!is.na(cebpe_ko_erythroid_Type)],
                                    row.names = lfc_cebpe_ko_erythroid_branch_Cebpa_targets)
  
pdf(paste(SI_fig_dir, "fig6_cebpe_leaking_example.1.pdf", sep = ''), width = 1, height = 1)
ggplot(aes(mean_expression, lfc, label = rownames(cebpe_ko_erythroid_df)), data = cebpe_ko_erythroid_df) + geom_point(size = 1) + nm_theme() + scale_size() #+ geom_text(size=2) 
dev.off()

#show one gene in WT: 
pdf(paste(SI_fig_dir, "fig6_cebpe_leaking_example.2.pdf", sep = ''), width = 1, height = 1)
plot_genes_branched_pseudotime(all_wt_GSE72857_cds[c("Lcn2", "0610007L01Rik"), ]) + nm_theme()
dev.off()

#show the same gene as scatter plot: 
pdf(paste(SI_fig_dir, "fig6_cebpe_leaking_example.3.pdf", sep = ''), width = 1, height = 1)
plot_genes_jitter(Cebpe_ko_mainfold_cells["Lcn2", ], grouping = "figSI4_Cebpe_ko") + nm_theme()
dev.off()

########################################################################################################################################################
# add Cebpa / Cebpe's branch time point plots: 
########################################################################################################################################################
pdf(paste(SI_fig_dir, "fig6_cebpae_branching_plot.pdf", sep = ''), width = 2.7, height = 3)
plot_genes_branched_pseudotime(all_wt_GSE72857_cds[c('Cebpa', 'Cebpe'), ], branch_labels = c('Erythroid', 'GMP'), min_expr = 0.01) + nm_theme()
dev.off()

########################################################################################################################################################################################################
# add the scatterplot similar to the blood KO datasets 
# cebpe KO cells on the MEP branch (erythroid branch)
########################################################################################################################################################################################################
# valid_subset_GSE72857_cds_beam_genes_proj_dup
# beam_genes <- row.names(subset(valid_subset_GSE72857_cds_beam_genes, qval < 0.01)) #valid_subset_GSE72857_cds_beam_genes_proj_dup all_GSE72857_cds_beam_genes all_GSE72857_cds_beam_genes_proj_dup

mep_beam_manifold_genes <- beam_genes #row.names(subset(valid_subset_GSE72857_cds_beam_genes_proj_dup, qval < 0.1)) #beam_genes #intersect(beam_genes, row.names(subset(URMM_all_abs_genotype_fig6_ko_test2_group_deg.1, qval < 0.1)))
# fold-change for the WT vs KO genes: 
# Gfi1_ko_fold_change <- log2(rowMeans(exprs(Gif1_ko_leaking_cells.1[monocyte_beam_manifold_genes, fig6_ko_test2.1 == "Gfi1 KO"]))/ 
#                                               rowMeans(exprs(Gif1_ko_leaking_cells.1[monocyte_beam_manifold_genes, fig6_ko_test2.1 == "WildType"])))
Cebpe_ko_fold_change <- log2(rowMeans(exprs(all_GSE72857_cds[mep_beam_manifold_genes, pData(all_GSE72857_cds)$State == 2 &  pData(all_GSE72857_cds)$Batch_desc == 'Cebpe KO']))/ 
                               #or pData(all_GSE72857_cds)$Batch_desc != 'Cebpa control'
                              rowMeans(exprs(all_GSE72857_cds[mep_beam_manifold_genes, pData(all_GSE72857_cds)$State == 2 &  pData(all_GSE72857_cds)$Batch_desc == 'Cebpe control'])))

gmp_vs_ery_fold_change <- log2(rowMeans(exprs(valid_subset_GSE72857_cds[mep_beam_manifold_genes, pData(valid_subset_GSE72857_cds)$State == 3]))/ 
                                  rowMeans(exprs(valid_subset_GSE72857_cds[mep_beam_manifold_genes, pData(valid_subset_GSE72857_cds)$State == 1])))

# create the plot:  
Cebpe_Targets <- mep_beam_manifold_genes %in% row.names(subset(all_GSE72857_cds_genotype_figSI41_group_deg.1, qval < 0.1)) # monocyte_beam_manifold_genes %in% Gfi1_targets_gene_short_name

Cebpe_fc_df <- data.frame('GMP Vs. MEP (WT tree)' = gmp_vs_ery_fold_change, 'Cebpe KO Vs. WT cells (MEP)' = Cebpe_ko_fold_change, 'Cebpe targets' = Cebpe_Targets)

pdf(paste(main_fig_dir, "SI19F.pdf", sep = ''), width = 1.5, height = 1.5)
qplot(GMP.Vs..MEP..WT.tree., Cebpe.KO.Vs..WT.cells..MEP., color = Cebpe.targets, data = Cebpe_fc_df, size = 0.5) + #geom_abline() +  
  xlab('GMP Vs. MEP (WT tree)') + ylab('Cebpe KO Vs. WT cells (MEP)') + scale_size(range = c(0.5, 0.5)) + 
  geom_vline(xintercept = 0, linetype = "longdash") + geom_hline(yintercept = 0, linetype = "longdash") + nm_theme() + xlim(c(-10, 10)) + ylim(c(-10, 10))
dev.off()

########################################################################################################################################################
# downsample the erythroid branch to check whether or not the small number of DEGs from the leaking cells to WT are not a power thing: 
########################################################################################################################################################
table(pData(Cebpe_ko_mainfold_cells)$figSI4_Cebpe_ko)
set.seed(20170207)
sample_wt_id <- sample(which(pData(Cebpe_ko_mainfold_cells)$figSI4_Cebpe_ko == 'WildType'), 81)
Cebpe_ko_mainfold_cells_subset <- Cebpe_ko_mainfold_cells[, c(which(pData(Cebpe_ko_mainfold_cells)$figSI4_Cebpe_ko == 'Cebpe KO'), sample_wt_id)]
all_GSE72857_cds_genotype_figSI41_group_deg_cebpe_subset <- differentialGeneTest(Cebpe_ko_mainfold_cells_subset, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~figSI4_Cebpe_ko", reducedModelFormulaStr = "~1")
sum(all_GSE72857_cds_genotype_figSI41_group_deg_cebpe_subset$qval < 0.1)

set.seed(20170207)
sample_wt_id2 <- sample(which(pData(Cebpe_ko_mainfold_cells.1)$figSI4_Cebpe_ko.1 == 'WildType'), 81)
sample_ko_id2 <- sample(which(pData(Cebpe_ko_mainfold_cells.1)$figSI4_Cebpe_ko.1 == 'Cebpe KO'), 81)
Cebpe_ko_mainfold_cells.1_subset <- Cebpe_ko_mainfold_cells.1[, c(sample_wt_id2, sample_ko_id2)]
all_GSE72857_cds_genotype_figSI41_group_deg.1_subset <- differentialGeneTest(Cebpe_ko_mainfold_cells.1_subset, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~figSI4_Cebpe_ko.1", reducedModelFormulaStr = "~1")

sum(all_GSE72857_cds_genotype_figSI41_group_deg.1_subset$qval < 0.1)

########################################################################################################################################################
# save the data 
########################################################################################################################################################
save.image('./RData/fig_SI6.RData')

