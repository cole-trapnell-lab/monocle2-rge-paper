library(monocle)
library(xacHelper)
library(grid)
library(mcclust)

revision_1_fig_dir <- "../Figures/First_revision/"
#update analysis for the simulation data 
load('../RData/fig3.RData')
cols <- c("lap" = "black", "original" = "gray", "pt" = "red")
cell_type_cols <- c("neuron" = "#F2746B", "astrocyte" = "#29B14A", "oligodendrocyte" = "#6E95CD")

##########################################################################################################################################################################
# branch time point
##########################################################################################################################################################################
nao_sim_cds <- setOrderingFilter(nao_sim_cds, ordering_genes = row.names(neuron_sim_cds))
nao_sim_cds <- reduceDimension(nao_sim_cds, norm_method = 'none', verbose = T, auto_param_selection = F, max_components = 3, pseudo_expr=0, scaling = F)
nao_sim_cds <- orderCells(nao_sim_cds)
plot_cell_trajectory(nao_sim_cds, color_by = 'State')
plot_cell_trajectory(nao_sim_cds, color_by = 'Pseudotime')

pdf(paste(revision_1_fig_dir, "fig2k.1.pdf", sep = ''), height = 2, width = 2)
qplot(nao_sim_cds@reducedDimK[1, ], nao_sim_cds@reducedDimK[2, ], color = I('black'), size = 0.1) + 
  geom_point(aes(nao_sim_cds_iter_1@reducedDimS[1, ], nao_sim_cds_iter_1@reducedDimS[2, ], 
                 size = 3 *  rep(1:400, 3) / 1000, alpha = 0.5, color = as.factor(rep(1:3, each = 400)))) + 
  scale_size(range = range(3 *  rep(1:400, 3) / 1000), limits = range(3 *  rep(1:400, 3) / 1000)) + 
  scale_x_reverse() + scale_y_reverse() + nm_theme() + xlab('Component 1') + ylab('Component 2') + monocle:::monocle_theme_opts()
dev.off()
nao_sim_cds <- orderCells(nao_sim_cds, root_state = 4)

max(pData(nao_sim_cds)[pData(nao_sim_cds)$State == 4, 'Pseudotime'])
min(pData(nao_sim_cds)[pData(nao_sim_cds)$State %in% 3, 'Pseudotime'])

max(pData(nao_sim_cds)[pData(nao_sim_cds)$State == 4, 'Pseudotime'])
min(pData(nao_sim_cds)[pData(nao_sim_cds)$State %in% c(1, 5), 'Pseudotime'])

##########################################################################################################################################################################
# robustness of Monocle 2 VS noise level 
##########################################################################################################################################################################
Noisify <- function(data, sd = 0.001) {
  
  if (is.vector(data)) {
    noise <- rnorm(length(data), 0, sd = sd)
    noisified <- data + noise
  } else {
    length <- dim(data)[1] * dim(data)[2]
    noise <- matrix(rnorm(length, 0, sd), dim(data)[1])
    noisified <- data + noise
  }
  return(noisified)
}

noise_cds_list <- list()
for(sd in c(0.0001, 0.01, 0.1, 0.5, 1, 1.2, 1.5, 1.8, 2, 3, 4, 5, 6)) {
  print(paste('current sd is', sd))
  nao_exprs_mat <- exprs(nao_sim_cds)
  nao_exprs_mat <- Noisify(nao_exprs_mat, sd = sd) 
  nao_sim_cds_noise <- make_cds(nao_exprs_mat, pd3, fd, expressionFamily = gaussianff())
  nao_sim_cds_noise <- reduceDimension(nao_sim_cds_noise, norm_method = 'none', verbose = T, auto_param_selection = F, max_components = 3, pseudo_expr=0, scaling = F)
  nao_sim_cds_noise <- orderCells(nao_sim_cds_noise)
  nao_sim_cds_noise <- orderCells(nao_sim_cds_noise, root_state = pData(nao_sim_cds_noise)$State[1])
  noise_cds_list <- c(noise_cds_list, nao_sim_cds_noise)
  # print(plot_cell_trajectory(nao_sim_cds_noise, color_by = 'Time'))
}
##########################################################################################################################################################################
# Benchmark for the branch time point under different noise level 
##########################################################################################################################################################################
#we can easily determine the branches by looking at the curves: 
marker_gene_df <- data.frame(Time = 1:400, 
                             expression = c(exprs(nao_sim_cds)['Mash1', ], 
                                            exprs(nao_sim_cds)['Hes5', ],
                                            exprs(nao_sim_cds)['Scl', ],
                                            exprs(nao_sim_cds)['Olig2', ]),
                             Gene = c(rep(c('Mash1', 'Hes5', 'Scl', 'Olig2'), each = 1200)),
                             Cell = c(rep(c('Neuron', 'Astrocyte', 'Oligodendrocyte'), each = 400)), 
                             Group = c(rep(c("Branch 1", "Branch 2"), each = 2400)))

pdf(paste(revision_1_fig_dir, "monocle_maker_gene_for_branch_time_point.pdf", sep = ''), height = 2, width = 2)
qplot(Time, expression, color = Cell, data = marker_gene_df, size = I(0.25)) + facet_wrap(~Gene + Group) + geom_vline(xintercept = c(20, 65)) + nm_theme()
dev.off()

pData(nao_sim_cds)$branch <- NA
pData(nao_sim_cds)$branch[c(1:20, 401:420, 801:820)] <- 1
pData(nao_sim_cds)$branch[c(21:400)] <- 2
pData(nao_sim_cds)$branch[c(421:465, 821, 825)] <- 3
pData(nao_sim_cds)$branch[c(466:800)] <- 4
pData(nao_sim_cds)$branch[c(825:1200)] <-5

cor_res <- c()
ari_res <- c()

for(i in 1:length(noise_cds_list)) {
  # show the pseudotime correlation as well as ARI values 
  cor_res <- c(cor_res, cor(pData(noise_cds_list[[i]])$Pseudotime, pData(nao_sim_cds)$Time))
  
  #calculate the branch time point to decide the branches: 
  
  ari_res <- c(ari_res, calClusteringMetrics(pData(noise_cds_list[[i]])$State, pData(nao_sim_cds)$branch)$randIndex[3])
}

noise_benchmark_res <- data.frame(noise = rep(c(0.0001, 0.01, 0.1, 0.5, 1, 1.2, 1.5, 1.8, 2, 3, 4, 5, 6), times = 2), 
           value = c(cor_res, ari_res), 
           Type = c(rep(c('Pseudotime correlation', 'Branch consistency (ARI)'), each = length(cor_res))))

pdf(paste(revision_1_fig_dir, "monocle_robustness_2_noise.pdf", sep = ''), height = 1, width = 1)
qplot(noise, value, data = noise_benchmark_res, color = Type, size = I(0.5)) + xlab('Guassian noise \n (s.d.)') + ylab('') + nm_theme()
dev.off()


pdf(paste(revision_1_fig_dir, "monocle_robustness_2_noise_helper.pdf", sep = ''))
qplot(noise, value, data = noise_benchmark_res, color = Type, size = I(0.5)) + xlab('Guassian noise \n (s.d.)') + ylab('') 
dev.off()
##########################################################################################################################################################################
# noisy expression & noisy tree
##########################################################################################################################################################################

noise_marker_gene_df <- data.frame(Time = 1:400, 
                             expression = c(exprs(noise_cds_list[[6]])['Mash1', ], 
                                            exprs(noise_cds_list[[6]])['Hes5', ],
                                            exprs(noise_cds_list[[6]])['Scl', ],
                                            exprs(noise_cds_list[[6]])['Olig2', ]),
                             Gene = c(rep(c('Mash1', 'Hes5', 'Scl', 'Olig2'), each = 1200)),
                             Cell = c(rep(c('Neuron', 'Astrocyte', 'Oligodendrocyte'), each = 400)), 
                             Group = c(rep(c("Branch 1", "Branch 2"), each = 2400)))

pdf(paste(revision_1_fig_dir, "monocle_maker_gene_for_branch_time_point_noise.pdf", sep = ''), height = 2, width = 2)
qplot(Time, expression, color = Cell, data = noise_marker_gene_df, size = I(0.025), alpha = 0.5) + facet_wrap(~Gene + Group) + geom_vline(xintercept = c(20, 65)) + geom_smooth(size = I(0.3)) + nm_theme()
dev.off()

pdf(paste(revision_1_fig_dir, "fig2k.1.noise_trajectory.pdf", sep = ''), height = 1, width = 1)
qplot(noise_cds_list[[6]]@reducedDimK[1, ], noise_cds_list[[6]]@reducedDimK[2, ], color = I('black'), size = I(0.1)) + 
  geom_point(aes(noise_cds_list[[6]]@reducedDimS[1, ], noise_cds_list[[6]]@reducedDimS[2, ], 
                 size = I(3 *  rep(1:400, 3) / 1000), alpha = 0.5, color = as.factor(rep(1:3, each = 400)))) + 
  # scale_size(range = range(3 *  rep(1:400, 3) / 1000), limits = range(3 *  rep(1:400, 3) / 1000)) + 
  scale_x_reverse() + scale_y_reverse() + nm_theme() + xlab('Component 1') + ylab('Component 2') + monocle:::monocle_theme_opts()
dev.off()

pdf(paste(revision_1_fig_dir, "fig2k.1.noise_trajectory_helper.pdf", sep = ''))
qplot(noise_cds_list[[6]]@reducedDimK[1, ], noise_cds_list[[6]]@reducedDimK[2, ], color = I('black'), size = I(0.1)) + 
  geom_point(aes(noise_cds_list[[6]]@reducedDimS[1, ], noise_cds_list[[6]]@reducedDimS[2, ], 
                 size = I(3 *  rep(1:400, 3) / 1000), alpha = 0.5, color = as.factor(rep(1:3, each = 400)))) #+ 
  # scale_size(range = range(3 *  rep(1:400, 3) / 1000), limits = range(3 *  rep(1:400, 3) / 1000)) + 
dev.off()


##########################################################################################################################################################################
# update the DPT algorithm
##########################################################################################################################################################################

DPT <- function (dm, tips = destiny::random_root(dm), ..., w_width = 0.1) {
  if (!is(dm, "DiffusionMap"))
    stop("dm needs to be of class DiffusionMap, not ", class(dm))
  if (!length(tips) %in% 1:3)
    stop("you need to specify 1-3 tips, got ", length(tips))
  dpt <- destiny:::dummy_dpt(dm)
  all_cells <- seq_len(nrow(dpt))
  stats <- destiny:::tipstats(dpt, all_cells, tips)
  
  if(stats$g < 1.1)
    stats$g <- 1.1
  
  branches <- destiny:::auto_branch(dpt, all_cells, stats, w_width)
  colnames(branches$branch) <- paste0("Branch", seq_len(ncol(branches$branch)))
  colnames(branches$tips) <- paste0("Tips", seq_len(ncol(branches$tips)))
  dpt@branch <- branches$branch
  dpt@tips <- branches$tips
  dpt
}

##########################################################################################################################################################################
# use the ground truth to benchmark our algorithm (update the pseudotime!!!!)
##########################################################################################################################################################################
auto_param_selection <- T
repeat_downsampling_num <- 25
benchmark_type <- 'na_sim_data' #
na_sim_cds <- orderCells(na_sim_cds, root_state = 2)
root_state <- row.names(subset(pData(na_sim_cds), State == 2))
AT1_state <- row.names(subset(pData(na_sim_cds), State == 1))
AT2_state <- row.names(subset(pData(na_sim_cds), State == 3))

pData(absolute_cds)$Pseudotime <- pData(absolute_cds)$Time
marker_gene_df <- data.frame(Time = 1:400, 
                             expression = c(exprs(na_sim_cds)['Mash1', ], 
                                            exprs(na_sim_cds)['Hes5', ]),
                             Gene = c(rep(c('Mash1', 'Hes5'), each = 400)),
                             Cell = c(rep(c('Neuron', 'Astrocyte'), each = 400)))

pdf(paste(revision_1_fig_dir, "na_marker_gene_branch.pdf", sep = ''), height = 1, width = 1)
qplot(Time, expression, color = Cell, data = marker_gene_df, size = I(0.25))  + geom_vline(xintercept = c(20, 65)) + nm_theme()
dev.off()

pData(absolute_cds)$State[c(1:25, 401:425)] <- 1
pData(absolute_cds)$State[c(26:400)] <- 2
pData(absolute_cds)$State[c(426:800)] <- 3

absolute_cds@expressionFamily <- gaussianff()

norm_method_data <- 'None' #for monocle 2 
normalize <- T #for DDRTree 

# DDRTree
source('../scripts/DDRTree_robustness_analysis_genearlize.R', echo = T)

# dpt 
exprs(absolute_cds)[, 401] <- exprs(absolute_cds)[, 401] + 0.001 #avoid removing the duplicated column 
normalize_data <- 'none'
source('./scripts/robustness_dpt_slicer_genearlize.R', echo = T)

# ICA 
#downsampling 300 cells for running ICA function: 
set.seed(20170404)
absolute_cds <- absolute_cds[, sample(1:ncol(absolute_cds), 300)]
source('./scripts/ICA_robustness_analysis_genearlize.R', echo = T)  #not finish the second 0.4 downsampling proportion (taking too long)

#use monocle2 for both monocle 1 and wishbone for downstream analysis: 
# ICA_sampling_res_df <- sampling_res_df
# ICA_valid_cell_sampling_res_df <- valid_cell_sampling_res_df
#load('lung_ICA_downsampling_res') #save the data 

downsampling_state <- 2 #check this 
source('./scripts/progressive_branch_downsampling_generalize.R', echo = T)

source('./scripts/revision_1_accuracy_na_simulation_res.R', echo = T)
##########################################################################################################################################################################
# do the same analysis for the HSMM data 
##########################################################################################################################################################################
load('../RData/fig1.RData')

HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes = HSMM_myo_ordering_genes)
HSMM_myo <- reduceDimension(HSMM_myo, verbose = T, auto_param_selection = T, norm_method = 'log')
HSMM_myo <- orderCells(HSMM_myo)
plot_cell_trajectory(HSMM_myo, color_by = 'Hours')
plot_cell_trajectory(HSMM_myo, color_by = 'State')
plot_cell_trajectory(HSMM_myo, color_by = 'Pseudotime')
HSMM_myo <- orderCells(HSMM_myo, root_state = GM_state(HSMM_myo))

# cds_subset = reduceDimension(cds_subset, norm_method = norm_method, auto_param_selection = auto_param, max_components = max_components, scaling = scaling, initial_method = run_new_dpt_function, maxIter = 20)

absolute_cds <- HSMM_myo
save(file = 'HSMM_myo_absolute_cds', absolute_cds)

auto_param_selection <- T
repeat_downsampling_num <- 25
benchmark_type <- 'HSMM_myo' #
root_state <- row.names(subset(pData(HSMM_myo), State == GM_state(HSMM_myo)))
AT1_state <- row.names(subset(pData(HSMM_myo), State == setdiff(1:3, GM_state(HSMM_myo))[1]))
AT2_state <- row.names(subset(pData(HSMM_myo), State == setdiff(1:3, GM_state(HSMM_myo))[2]))

#norm_method_data <- 'None' #for monocle 2 
#normalize <- T #for DDRTree 

# absolute_cds <- reduceDimension(absolute_cds, auto_param_selection = T)
# absolute_cds <- orderCells(absolute_cds)
# absolute_cds <- orderCells(absolute_cds, root_state = 2)

cor(pData(absolute_cds)$Pseudotime, pData(HSMM_myo)$Pseudotime)
identical(pData(absolute_cds)$State, pData(HSMM_myo)$State)

source('../scripts/DDRTree_robustness_analysis_genearlize.R', echo = T)
source('../scripts/ICA_robustness_analysis_genearlize.R', echo = T)  #not finish the second 0.4 downsampling proportion (taking too long)

#use monocle2 for both monocle 1 and wishbone for downstream analysis: 
# ICA_sampling_res_df <- sampling_res_df
# ICA_valid_cell_sampling_res_df <- valid_cell_sampling_res_df
#load('lung_ICA_downsampling_res') #save the data 
source('../scripts/robustness_dpt_slicer_genearlize.R', echo = T)

downsampling_state <- 2 #check this 
source('../scripts/progressive_branch_downsampling_generalize.R', echo = T)

##########################################################################################################################################################################
# do the same analysis for the Nature blood data 
##########################################################################################################################################################################
# load('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/RData/fig5.RData')
load('../RData/fig5_4_22.RData')
absolute_cds <- URMM_all_fig1b

#URMM_all_fig1b <- reduceDimension(URMM_all_fig1b, verbose = T, scaling = T, max_components = 4, maxIter = 100, lambda = 20 * ncol(URMM_all_fig1b)) #, maxIter = 100, initial_method = DM, R_tSNE, tSNE, destiny_diffusionMaps, maxIter = 100 , param.gamma = 100, ncenter = 100
# 
# absolute_cds <- reduceDimension(absolute_cds, max_components = 4, norm_method = 'log', verbose = T, auto_param_selection = F) # 
# absolute_cds <- orderCells(absolute_cds)
absolute_cds <- trimTree(absolute_cds)

monocle:::normalize_expr_data(absolute_cds, norm_method = 'log', pseudo_expr = NULL)

plot_cell_trajectory(absolute_cds)
#initial_method = PCA, norm_method = 'vstExprs', order_by = NULL

auto_param_selection <- T
repeat_downsampling_num <- 25
benchmark_type <- 'URMM_all_fig1b_centroid' # URMM_all_fig1b
num_dim = 4
root_state <- row.names(subset(pData(absolute_cds), State == 1))
AT1_state <- row.names(subset(pData(absolute_cds), State == 2))
AT2_state <- row.names(subset(pData(absolute_cds), State == 3))

trajectory_color_by <- 'cluster'
#norm_method_data <- 'None' #for monocle 2 
#normalize <- T #for DDRTree 

source('../scripts/DDRTree_robustness_analysis_genearlize.R', echo = T)
source('../scripts/ICA_robustness_analysis_genearlize.R', echo = T)  #not finish the second 0.4 downsampling proportion (taking too long)

#use monocle2 for both monocle 1 and wishbone for downstream analysis: 
# ICA_sampling_res_df <- sampling_res_df
# ICA_valid_cell_sampling_res_df <- valid_cell_sampling_res_df
#load('lung_ICA_downsampling_res') #save the data 
pData(absolute_cds)$Time <- pData(absolute_cds)[, trajectory_color_by]
normalize <- F
source('../scripts/robustness_dpt_slicer_genearlize.R', echo = T)

downsampling_state <- 2 #check this 
source('../scripts/progressive_branch_downsampling_generalize.R', echo = T)

##########################################################################################################################################################################
# save data
##########################################################################################################################################################################
save.image('../RData/revision_1_simulation_data.RData')
