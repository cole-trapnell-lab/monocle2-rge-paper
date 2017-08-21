#1. BEAM genes for the HSMM and the URMM 
load('/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/RData/fig_si2.RData')
dim(subset(HSMM_BEAM_res, qval < 0.1))
# [1] 887   8

#2. 
load('/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/RData/fig_SI6.RData')
length(row.names(subset(all_GSE72857_cds_genotype_figSI41_group_deg_cebpe, qval < 0.1))) #305 
length(intersect(row.names(subset(all_GSE72857_cds_genotype_figSI41_group_deg_cebpe, qval < 0.1)), beam_genes)) #overlap genes; #133
length(row.names(subset(all_GSE72857_cds_genotype_figSI41_group_deg.1, qval < 0.1))) #3509

#3.
load('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/RData/analysis_score_ordering_MAR_seq.RData')
nrow(MAP_cells_clusters) - length(pData(valid_subset_GSE72857_cds)$cell_type)
# [1] 31
table(pData(valid_subset_GSE72857_cds)$cell_type)
# CMP        DC erythroid       GMP
# 451        30      1095      1123

#4. 
