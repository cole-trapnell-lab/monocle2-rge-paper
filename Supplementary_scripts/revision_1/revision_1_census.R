library(monocle)
library(xacHelper)

load('./RData/fig1.RData')

##########################################################################################################################################################################
# A function to do all the process (Census Vs. not census)
##########################################################################################################################################################################
revision_1_fig_dir <- "./Figures/First_revision/"

#1. determine how many pca dimension you want:
HSMM_myo_std <- HSMM[, colnames(HSMM_myo)]

fData(HSMM_myo_std)$use_for_ordering <- fData(HSMM_myo_std)$num_cells_expressed > round(ncol(HSMM_myo_std) / 10)
HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(HSMM_myo_std, return_all = F)

#2. run reduceDimension with tSNE as the reduction_method
HSMM_myo_std <- reduceDimension(HSMM_myo_std, max_components=2, norm_method = 'log', num_dim = 6, reduction_method = 'tSNE', verbose = T)

#3. initial run of clusterCells_Density_Peak
HSMM_myo_std <- clusterCells_Density_Peak(HSMM_myo_std, verbose = T)

#4. check the clusters (there are three clusters)
plot_cell_clusters(HSMM_myo_std, color_by = 'as.factor(Cluster)', show_density = F)
plot_cell_clusters(HSMM_myo_std, color_by = 'as.factor(Hours)', show_density = F)

#5. also check the decision plot
plot_rho_delta(HSMM_myo_std, rho_threshold = 3.8, delta_threshold = 4 )
plot_cell_clusters(HSMM_myo_std, color_by = 'as.factor(Cluster)', show_density = F, rho_threshold = 3.8, delta_threshold = 5)
plot_cell_clusters(HSMM_myo_std, color_by = 'as.factor(Hours)', show_density = F, rho_threshold = 3.8, delta_threshold = 5)

#6. re-run cluster and skipping calculating the rho_sigma
HSMM_myo_std <- clusterCells_Density_Peak(HSMM_myo_std, verbose = T,  rho_threshold = 3.8, delta_threshold = 5, skip_rho_sigma = T)

#7. make the final clustering plot:
plot_cell_clusters(HSMM_myo_std, color_by = 'as.factor(Cluster)', show_density = F)

#find the important genes and use that for lineage tree construction
#perform DEG test across clusters:
pData(HSMM_myo_std)$Cluster <- as.character(pData(HSMM_myo_std)$Cluster)
HSMM_myo_std_DEG_genes <- differentialGeneTest(HSMM_myo_std, fullModelFormulaStr = '~Cluster', cores = detectCores() - 2)
# HSMM_myo_std_DEG_genes_subset <- HSMM_myo_std[fData(HSMM_myo_std)$num_cells_expressed > round(ncol(HSMM_myo_std) / 10, ]

#use all DEG gene from the clusters
HSMM_myo_std_ordering_genes <- row.names(HSMM_myo_std_DEG_genes)[order(HSMM_myo_std_DEG_genes$qval)][1:1000] #row.names(subset(HSMM_myo_std_DEG_genes, qval < qval_thrsld))

# HSMM_myo_std_ordering_genes <- row.names(HSMM_myo_std_clustering_DEG_genes_subset)[order(HSMM_myo_std_clustering_DEG_genes_subset$qval)][1:1000]
HSMM_myo_std_ordering_genes

HSMM_myo_std <- setOrderingFilter(HSMM_myo_std, ordering_genes = HSMM_myo_std_ordering_genes)
HSMM_myo_std <- reduceDimension(HSMM_myo_std, verbose = T, auto_param_selection = F)
HSMM_myo_std <- orderCells(HSMM_myo_std)
plot_cell_trajectory(HSMM_myo_std, color_by = 'Hours')

plot_spanning_tree(HSMM_myo_std, color_by="Hours", markers=c("MYOG", "CCNB2"))

# HSMM <- setOrderingFilter(HSMM, ordering_genes = HSMM_myo_std_ordering_genes)
# HSMM <- reduceDimension(HSMM, norm_method = 'log', verbose = T, maxIter = 100)
# HSMM <- orderCells(HSMM)
# plot_cell_trajectory(HSMM, color_by = 'Hours')

pdf(paste(revision_1_fig_dir, "HSMM_myo_std_trajectory.pdf", sep = ''), height = 1, width = 1)
plot_cell_trajectory(HSMM_myo_std, color_by = 'Time', show_branch_points = T, cell_size = 0.5) + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

plot_cell_trajectory(HSMM_myo_std, color_by = 'Pseudotime')
plot_cell_trajectory(HSMM_myo_std, color_by = 'State')

HSMM_myo_std <- orderCells(HSMM_myo_std, root_state = GM_state(HSMM_myo_std))
#show BEAM results: 
HSMM_std_expressed_genes <- row.names(subset(fData(HSMM_myo_std), num_cells_expressed >=5))
HSMM_BEAM_std_res <- BEAM(HSMM_myo_std[HSMM_expressed_genes,], 
                      branch_point=2,
                      cores = detectCores() - 2)

head(HSMM_BEAM_std_res) #show the result

load('./RData/fig_si2.RData')

length(intersect(row.names(subset(HSMM_BEAM_res, qval < 0.1)), row.names(subset(HSMM_BEAM_std_res, qval < 0.1)))) / nrow(subset(HSMM_BEAM_res, qval < 0.1))

element_all <- c(row.names(subset(HSMM_BEAM_res, qval < 0.1)), 
                 row.names(subset(HSMM_BEAM_std_res, qval < 0.1)))
sets_all <- c(rep(paste('Transcript counts (Size + VST)', sep = ''), nrow(subset(HSMM_BEAM_res, qval < 0.1))), 
              rep(paste('FPKM', sep = ''), nrow(subset(HSMM_BEAM_std_res, qval < 0.1))))
pdf(paste(revision_1_fig_dir, "census_URMM_BEAM_res.pdf", sep = ''))
venneuler_venn(element_all, sets_all)
dev.off()

# load('./RData/fig_si2_1_30.RData')
save.image('./RData/fig_revision_1_census_HSMM.RData')
  
#clean the data: 
rm(list = ls())
##########################################################################################################################################################################
# run analysis on the blood dataset (not very useful here)
##########################################################################################################################################################################
load('./RData/fig5.RData')

URMM_all_std <- order_cell_tree(row.names(fig1b)[2:533], cds = URMM_all_std[, colnames(URMM_all_fig1b)],
                                 initial_method = NULL, norm_method = 'log', order_by = NULL)
plot_spanning_tree(URMM_all_std, color_by = 'Type', cell_size = 2) + scale_color_manual(values = cols, name = "Type") #+ nm_theme()

#show the trajectory here as well as BEAM results: 
order_cell_tree <- function(order_genes, cds, initial_method = PCA, norm_method = 'log',order_by="Time", max_component = 2) {
  cds <- setOrderingFilter(cds, order_genes)
  plot_ordering_genes(cds)
  
  cds <- reduceDimension(cds, max_components = max_component, norm_method = norm_method, initial_method = initial_method, verbose = T, auto_param_selection = F) #, initial_method = DM , initial_method = DM
  cds <- orderCells(cds, num_paths=1, root_state = NULL)
  PD <- pData(cds)
  
  if(!is.null(order_by)){
    avg_pseudotime <- ddply(PD, .(Time), function(x) mean(x$Pseudotime))
    starting_cell_states <- apply(table(PD[, c("Time", "State")]), 1, function(x) which(x == max(x)))[1]
    cds <- orderCells(cds, num_paths=1, root_state = starting_cell_states)
  }
  # return(pData(HSMM_myo))
  return(cds)
}

URMM_all_std <- setOrderingFilter(URMM_all_std, ordering_genes = c(URMM_ordering_genes))
URMM_all_std <- reduceDimension(URMM_all_std, verbose = T, scaling = T, max_components = 4, maxIter = 100, lambda = 20 * ncol(URMM_all_fig1b)) #, maxIter = 100, initial_method = DM, R_tSNE, tSNE, destiny_diffusionMaps, maxIter = 100 , param.gamma = 100, ncenter = 100
#
URMM_all_std <- orderCells(URMM_all_std)
plot_cell_trajectory(URMM_all_fig1b, color_by = 'Type')

pData(URMM_all_std)$cluster <- pData(URMM_all_fig1b[, colnames(URMM_all_std)])$cluster
pdf(paste(revision_1_fig_dir, "URMM_all_fig1b_std_trajectory.pdf", sep = ''), height = 1, width = 1)
plot_cell_trajectory(URMM_all_std, color_by = 'cluster', cell_size = 0.5) + scale_color_manual(values = cluster_cols, name = "Type") + nm_theme()
dev.off()

# length(intersect(row.names(subset(HSMM_BEAM_res, qval < 0.1)), row.names(subset(HSMM_BEAM_std_res, qval < 0.1)))) / nrow(subset(HSMM_BEAM_res, qval < 0.1))
# 
# element_all <- c(row.names(subset(HSMM_BEAM_res, qval < 0.1)), 
#                  row.names(subset(HSMM_BEAM_std_res, qval < 0.1)))
# sets_all <- c(rep(paste('Transcript counts (Size + VST)', sep = ''), nrow(subset(HSMM_BEAM_res, qval < 0.1))), 
#               rep(paste('FPKM', sep = ''), nrow(subset(HSMM_BEAM_std_res, qval < 0.1))))
# pdf(paste(revision_1_fig_dir, "census_URMM_BEAM_res.pdf", sep = ''))
# venneuler_venn(element_all, sets_all)
# dev.off()

##########################################################################################################################################################################
# run analysis on the blood dataset (not very useful here)
##########################################################################################################################################################################

# save.image('./RData/revision_1_census.RData')
save.image('./RData/fig_revision_1_census_URMM.RData')

