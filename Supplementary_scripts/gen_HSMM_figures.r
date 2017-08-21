rm(list = ls())
####################################################################################################################################################################################
#load all package
####################################################################################################################################################################################
library(monocle)
library(R.matlab)
library(igraph)
library(MASS)
library(simplePPT)
library(xacHelper)
library(reshape2)

source('./scripts/function.R')
source('./scripts/plotting.R')

main_fig_dir <- "./Figures/main_figures/"
SI_fig_dir <- "./Figures/supplementary_figures/"

GM_state <- function(cds){
  T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
  as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))])
}

####################################################################################################################################################################################
load('./RData/analysis_HSMM_data.RData')
marker_ids <- row.names(subset(fData(HSMM), gene_short_name %in% c("CCNB2", "DMD", "MYOG", "TNNT1")))

#CNS simulation data:
cell_simulate <- readMat('./mat_data/cell_simulate.mat')
all_cell_simulation <- cell_simulate$cell.simulate[, 1:400, ] #time 0-20 are the period two branches appear
row.names(all_cell_simulation) <-  c('Pax6', 'Mash1', 'Brn2', 'Zic1', 'Tuj1', 'Hes5', 'Scl', 'Olig2', 'Stat3', 'Myt1L', 'Aldh1L', 'Sox8', 'Mature')

#obtain the corresponding lineage for each simulation run:
neuron_cell_ids <- which(all_cell_simulation['Mash1', 400, ] > 3)
astrocyte_cell_ids <- which(all_cell_simulation['Scl', 400, ] > 2)
oligodendrocyte_cell_ids <- which(all_cell_simulation['Olig2', 400, ] > 2)

####################################################################################################################################################################################
#create muscle tree with monocle 1/2
HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes = HSMM_myo_ordering_genes)
HSMM_myo <- reduceDimension(HSMM_myo, verbose = T, maxIter = 100)
HSMM_myo <- orderCells(HSMM_myo)
HSMM_myo <- orderCells(HSMM_myo, root_state=GM_state(HSMM_myo))

pData(HSMM_myo)$Hours <- pData(HSMM)[colnames(HSMM_myo), 'Hours']
plot_cell_trajectory(HSMM_myo, color_by = 'Hours')
plot_spanning_tree(HSMM_myo, color_by="Hours", markers=c("MYOG", "CCNB2"))

#reorder the tree using ICA with the same genes:
HSMM_myo_ICA <- reduceDimension(HSMM_myo, norm_method = 'log', reduction_method = 'ICA', verbose = T)
HSMM_myo_ICA <- orderCells(HSMM_myo_ICA, num_paths = 2) #origin/nbt_3nd_submission do on nbt 3
HSMM_myo_ICA <- orderCells(HSMM_myo_ICA, root_state = GM_state(HSMM_myo_ICA))
plot_cell_trajectory(HSMM_myo_ICA, color_by = 'Hours')
plot_spanning_tree(HSMM_myo_ICA, color_by="Hours", markers=c("MYOG", "CCNB2"))

####################################################################################################################################################################################
#Create Figure 1b
####################################################################################################################################################################################
pdf(paste(main_fig_dir, "fig1b.pdf", sep = ''), height = 1, width = 1)
plot_cell_trajectory(HSMM_myo_ICA, color_by = 'Hours', cell_size=0.5, theta =180) + nm_theme() +  scale_color_manual(values = HSMM_cols)
dev.off()

####################################################################################################################################################################################
#Create Figure 1c
####################################################################################################################################################################################
pdf(paste(main_fig_dir, "fig1c.pdf", sep = ''), height = 1, width = 1)
plot_cell_trajectory(HSMM_myo, color_by="Hours", show_branch_points = F, cell_size=0.5, theta =30) + nm_theme() +  scale_color_manual(values = HSMM_cols)
dev.off()

####################################################################################################################################################################################
#panel E (Not used in the final version)
####################################################################################################################################################################################
pdf(paste(main_fig_dir, "fig2e_HSMM.pdf", sep = ''), height = 1, width = 3)
plot_genes_in_pseudotime(HSMM_myo_ICA[marker_ids,], color_by="Hours", cell_size=0.5, ncol = 4, min_expr = 0.5) + nm_theme() +  scale_color_manual(values = HSMM_cols) +
  ylab('Expression') + xlab('Pseudotime (stretched)')
dev.off()

pdf(paste(main_fig_dir, "fig2e.1_HSMM.pdf", sep = ''), height = 1.5, width = 1.5)
plot_genes_in_pseudotime(HSMM_myo_ICA[marker_ids,], color_by="Hours", cell_size=0.5, ncol = 2) + nm_theme() +  scale_color_manual(values = HSMM_cols) +
  ylab('Expression') + xlab('Pseudotime (stretched)')
dev.off()

####################################################################################################################################################################################
#panel F (Not used in the final version)
####################################################################################################################################################################################
pData(HSMM_myo)$
pdf(paste(main_fig_dir, "fig2f_HSMM.pdf", sep = ''), height = 1, width = 3)
plot_genes_branched_pseudotime(HSMM_myo[marker_ids,], color_by="Hours", cell_size=0.5, ncol = 4, min_expr = 0.5) + nm_theme() +  scale_color_manual(values = HSMM_cols)
dev.off()

pdf(paste(main_fig_dir, "fig2f_HSMM_helper.pdf", sep = ''))
plot_genes_branched_pseudotime(HSMM_myo[marker_ids,]) 
dev.off()
####################################################################################################################################################################################
#supplmentary figure 4b (percentage of cells)
####################################################################################################################################################################################
c('ANPEP', 'CDK1', 'DMD', 'ENO3', 'H19', 'ID3', 'MEF2C', 'MYF5', 'MYH2', 'MYOG', 'TMEMBC', 'TNNT2')

####################################################################################################################################################################################
# save data 
####################################################################################################################################################################################
save.image('./RData/fig2.RData')

