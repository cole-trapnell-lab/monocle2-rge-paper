#1. lung data benchmark (fig3.r)
#2. progressive downsampling for the mar-seq data 
#3. complex tree benchmark 
#4. unbalanced branch downsampling benchmark 
#5. all parameters in DDRTree

################################################################################################################################################################################################################################################
# run the complex tree benchmark
################################################################################################################################################################################################################################################
source('./script_for_reproduce/analysis_complex_tree_structure.R', echo = T)

################################################################################################################################################################################################################################################
#fig_si panel c: benchmark on testing different parameters: ncenter, param.gamma, and the number of dimensions used in the reduction
################################################################################################################################################################################################################################################
load('./RData/fig1.RData')
save(file = 'lung_absolute_cds', absolute_cds)
rm(list = ls()); load('')

default_ncenter <- NULL
benchmark_type <- 'lung_data' #
source('./script_for_reproduce/analysis_DDRTree_params.R', echo = T)
save.image('./RData/lung_fig4_all_params.RData')

rm(list = ls())

################################################################################################################################################################################################################################################
#fig_si panel c: benchmark on testing different parameters: ncenter, param.gamma, and the number of dimensions used in the reduction
################################################################################################################################################################################################################################################
load('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/RData/analysis_score_ordering_MAR_seq.RData')

benchmark_type <- 'Marseq_data' 
absolute_cds <- valid_subset_GSE72857_cds
pData(absolute_cds)$Time <- pData(absolute_cds)$cell_type #cell type as Time; Time matches that of cell_type
source('./script_for_reproduce/analysis_DDRTree_params.R', echo = T)
save.image('./RData/fig4_marseq_all_params.RData')

################################################################################################################################################################################################################################################
# unbalanced branch downsampling benchmark (not included in the manuscript)
################################################################################################################################################################################################################################################
