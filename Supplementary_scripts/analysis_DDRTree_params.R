#script to demonstrate the effects of parameters on the tree structure:

# dimensions = 2, initial_method = NULL, maxIter = 20,
# sigma = 0.001, lambda = NULL, ncenter = NULL, param.gamma = 10,
# tol = 0.001
library(monocle)
library(SLICER)
library(dpt)
library(destiny)
library(mcclust)
library(plyr)
library(reshape2)
library(xacHelper)
library(diffusionMap)

SI_fig_dir <- "./Figures/supplementary_figures/"
source('./scripts/function.R', echo = T)

#function for running DDRTree with different parameters:
root_cell_name <- row.names(subset(pData(absolute_cds), Pseudotime == 0))
order_cell_tree_all_params <- function (order_genes, cds, order_by = "Time",
                                        dimensions = 2, initial_method = NULL, maxIter = 20,
                                        sigma = 0.001, lambda = NULL, ncenter = NULL, param.gamma = 10,
                                        tol = 0.001) {
  cds <- setOrderingFilter(cds, order_genes)
  plot_ordering_genes(cds)
  cds <- reduceDimension(cds, max_components = 2, norm_method = "log",
                         initial_method = initial_method, verbose = T,
                         auto_param_selection = F, dimensions = dimensions, maxIter = maxIter,
                         sigma = sigma, lambda = lambda, ncenter = ncenter, param.gamma = param.gamma,
                         tol = tol)
  cds <- orderCells(cds, root_state = NULL)
  PD <- pData(cds)
  if (!is.null(order_by)) {
    avg_pseudotime <- ddply(PD, order_by, function(x) mean(x$Pseudotime))
    #note that for lung or MAR_seq data: the first order_by character alphabetically is #E14.5 or CMP so we can just use the [1] to get the first element
    starting_cell_states <- apply(table(PD[, c(order_by, "State")]),
                                  1, function(x) which(x == max(x)))[1]
    root_state <- pData(cds[, root_cell_name])$State
    if(root_state != 1)
      cds <- orderCells(cds, root_state = starting_cell_states)
  }
  return(cds)
}

# sample_DDRTree_parms_cmbn <- expand.grid(dimensions = seq(2, 10, 1), initial_method = c(ICA, PCA, LLE, DM, destiny_diffusionMaps),
#                                maxIter= seq(10, 100, 5), sigma = seq(0.001, 10, length.out = 20),
#                                lambda = seq(0.5, 50, 20) *  ncol(absolute_cds), ncenter = c(NULL, seq(0.01, 0.99, length.out = 20) * ncol(absolute_cds)),
#                                param.gamma = seq(0.1, 1000, length.out = 20), tol = 0.001)

#run ordering with different parameters to benchmark the results:
if(!exists("auto_param_selection"))
  ordering_genes <- row.names(absolute_cds)[fData(absolute_cds)$use_for_ordering]
if(!exists("default_ncenter"))
  default_ncenter <- monocle:::cal_ncenter(ncol(absolute_cds))

dimensions <- seq(2, 10, 1)
initial_method <- c(ICA, PCA, LLE2, destiny_diffusionMaps)
maxIter <- seq(10, 100, 5)
sigma <-  seq(0.001, 1, length.out = 10)
lambda <- round(seq(0.5, 50, 5) *  ncol(absolute_cds))
ncenter <- round(c(NULL, seq(0.01, 0.99, length.out = 20) * ncol(absolute_cds))[-1])
param.gamma <- c(2, 5, 8, 10, 15, 20, 25, 30, 40, 50)
tol <- 0.001

if(!exists('order_by_type'))
  order_by_type <- 'Time'

sample_dimension_cds_list <- lapply(dimensions,
                                function(x) {
                                  message('Current dimension is ', x)
                                  order_cell_tree_all_params(ordering_genes, absolute_cds, order_by = order_by_type,
                                                             dimensions = x, initial_method = NULL, maxIter = 20,
                                                             sigma = 0.001, lambda = NULL, ncenter = default_ncenter, param.gamma = 10,
                                                             tol = 0.001)
                                })

sample_initial_method_cds_list <- lapply(initial_method,
                                   function(x) {
                                     # message('Current initial_method is ', x)
                                     order_cell_tree_all_params(ordering_genes, absolute_cds, order_by = order_by_type,
                                                                dimensions = 2, initial_method = x, maxIter = 20,
                                                                sigma = 0.001, lambda = NULL, ncenter = default_ncenter, param.gamma = 10,
                                                                tol = 0.001)
                                   })

sample_maxIter_cds_list <- lapply(maxIter,
                                   function(x) {
                                     message('Current maxIter is ', x)
                                     order_cell_tree_all_params(ordering_genes, absolute_cds, order_by = order_by_type,
                                                                dimensions = 2, initial_method = NULL, maxIter = x,
                                                                sigma = 0.001, lambda = NULL, ncenter = default_ncenter, param.gamma = 10,
                                                                tol = 0.001)
                                   })

sample_sigma_cds_list <- lapply(sigma,
                                   function(x) {
                                     message('Current sigma is ', x)
                                     order_cell_tree_all_params(ordering_genes, absolute_cds, order_by = order_by_type,
                                                                dimensions = 2, initial_method = NULL, maxIter = 20,
                                                                sigma = x, lambda = NULL, ncenter = default_ncenter, param.gamma = 10,
                                                                tol = 0.001)
                                   })

sample_lambda_cds_list <- lapply(lambda,
                                   function(x) {
                                     message('Current lambda is ', x)
                                     order_cell_tree_all_params(ordering_genes, absolute_cds, order_by = order_by_type,
                                                                dimensions = 2, initial_method = NULL, maxIter = 20,
                                                                sigma = 0.001, lambda = x, ncenter = default_ncenter, param.gamma = 10,
                                                                tol = 0.001)
                                   })

sample_param.gamma_cds_list <- lapply(param.gamma,
                                function(x) {
                                  message('Current param.gamma is ', x)
                                  order_cell_tree_all_params(ordering_genes, absolute_cds, order_by = order_by_type,
                                                             dimensions = 2, initial_method = NULL, maxIter = 20,
                                                             sigma = 0.001, lambda = NULL, ncenter = default_ncenter, param.gamma = x,
                                                             tol = 0.001)
                                })

order_cell_tree_all_params <- function (order_genes, cds, order_by = order_by_type,
                                        dimensions = 2, initial_method = NULL, maxIter = 20,
                                        sigma = 0.001, lambda = NULL, ncenter = NULL, param.gamma = 10,
                                        tol = 0.001) {
  cds <- setOrderingFilter(cds, order_genes)
  plot_ordering_genes(cds)
  cds <- reduceDimension(cds, max_components = 2, norm_method = "log",
                         initial_method = initial_method, verbose = F,
                         auto_param_selection = F, dimensions = dimensions, maxIter = maxIter,
                         sigma = sigma, lambda = lambda, ncenter = ncenter, param.gamma = param.gamma,
                         tol = tol)
  cds <- orderCells(cds, root_state = NULL)
  PD <- pData(cds)
  if (!is.null(order_by)) {
    avg_pseudotime <- ddply(PD, order_by, function(x) mean(x$Pseudotime))
    #note that for lung or MAR_seq data: the first order_by character alphabetically is #E14.5 or CMP so we can just use the [1] to get the first element
    starting_cell_states <- apply(table(PD[, c(order_by, "State")]),
                                  1, function(x) which(x == max(x)))[1]
    cds <- orderCells(cds, root_state = starting_cell_states)
  }
  return(cds)
}

sample_ncenter_cds_list <- lapply(ncenter,
                                function(x) {
                                  message('Current ncenter is ', x)
                                  tryCatch({
                                    order_cell_tree_all_params(ordering_genes, absolute_cds, order_by = order_by_type,
                                                               dimensions = 2, initial_method = NULL, maxIter = 20,
                                                               sigma = 0.001, lambda = NULL, ncenter = x, param.gamma = 10,
                                                               tol = 0.001)
                                  }, error = function(e){
                                    print(e)
                                    return('there is an error')
                                  })
                                })

################################################################################################################################################################################
# this is used to reorder the cells so that they start from the correct root
################################################################################################################################################################################
reorder_cell_cds <- function (cds, root_cell_name) {
  starting_cell_states <- pData(cds[, root_cell_name])$State
  if(starting_cell_states != 1){
    tryCatch({
      cds <- orderCells(cds, root_state = starting_cell_states)
    }, error = function(e){
      print(e)
      return(cds)
    })
  }
  return(cds)
}

#reorder some of the cds:
sample_dimension_cds_list2 <- lapply(sample_dimension_cds_list, reorder_cell_cds, root_cell_name)
sample_initial_method_cds_list2 <- lapply(sample_initial_method_cds_list, reorder_cell_cds, root_cell_name)
sample_maxIter_cds_list2 <- lapply(sample_maxIter_cds_list, reorder_cell_cds, root_cell_name)
sample_sigma_cds_list2 <- lapply(sample_sigma_cds_list, reorder_cell_cds, root_cell_name)
sample_lambda_cds_list2 <- lapply(sample_lambda_cds_list, reorder_cell_cds, root_cell_name)
sample_ncenter_cds_list2 <- lapply(sample_ncenter_cds_list, reorder_cell_cds, root_cell_name)

################################################################################################################################################################################
#check pseudotime for parameters lambad:
################################################################################################################################################################################
plot_cell_trajectory(sample_lambda_cds_list2[[1]])
plot_cell_trajectory(sample_lambda_cds_list2[[1]], color_by = 'Pseudotime')
sample_lambda_cds_list2[[1]] <- orderCells(sample_lambda_cds_list2[[1]], root_state = 2) #

plot_cell_trajectory(sample_lambda_cds_list2[[5]])
plot_cell_trajectory(sample_lambda_cds_list2[[5]], color_by = 'Pseudotime')
sample_lambda_cds_list2[[5]] <- orderCells(sample_lambda_cds_list2[[5]], root_state = 2)

plot_cell_trajectory(sample_lambda_cds_list2[[6]])
plot_cell_trajectory(sample_lambda_cds_list2[[6]], color_by = 'Pseudotime')
sample_lambda_cds_list2[[6]] <- orderCells(sample_lambda_cds_list2[[6]], root_state = 2)

plot_cell_trajectory(sample_lambda_cds_list2[[7]])
plot_cell_trajectory(sample_lambda_cds_list2[[7]], color_by = 'Pseudotime')
sample_lambda_cds_list2[[7]] <- orderCells(sample_lambda_cds_list2[[7]], root_state = 2)

plot_cell_trajectory(sample_lambda_cds_list2[[8]])
plot_cell_trajectory(sample_lambda_cds_list2[[8]], color_by = 'Pseudotime')
sample_lambda_cds_list2[[8]] <- orderCells(sample_lambda_cds_list2[[8]], root_state = 2)

plot_cell_trajectory(sample_lambda_cds_list2[[9]])
plot_cell_trajectory(sample_lambda_cds_list2[[9]], color_by = 'Pseudotime')
sample_lambda_cds_list2[[9]] <- orderCells(sample_lambda_cds_list2[[9]], root_state = 2)

plot_cell_trajectory(sample_lambda_cds_list2[[10]])
plot_cell_trajectory(sample_lambda_cds_list2[[10]], color_by = 'Pseudotime')
sample_lambda_cds_list2[[10]] <- orderCells(sample_lambda_cds_list2[[10]], root_state = 2)

################################################################################################################################################################################
#check pseudotime for parameters param.gamma:
################################################################################################################################################################################
plot_cell_trajectory(sample_param.gamma_cds_list[[2]])
plot_cell_trajectory(sample_param.gamma_cds_list[[2]], color_by = 'Pseudotime')
sample_param.gamma_cds_list[[2]] <- orderCells(sample_param.gamma_cds_list[[2]], root_state = 2)

plot_cell_trajectory(sample_param.gamma_cds_list[[5]])
plot_cell_trajectory(sample_param.gamma_cds_list[[5]], color_by = 'Pseudotime')
sample_param.gamma_cds_list[[5]] <- orderCells(sample_param.gamma_cds_list[[5]], root_state = 2)

plot_cell_trajectory(sample_param.gamma_cds_list[[6]])
plot_cell_trajectory(sample_param.gamma_cds_list[[6]], color_by = 'Pseudotime')
sample_param.gamma_cds_list[[6]] <- orderCells(sample_param.gamma_cds_list[[6]], root_state = 2)

plot_cell_trajectory(sample_param.gamma_cds_list[[7]])
plot_cell_trajectory(sample_param.gamma_cds_list[[7]], color_by = 'Pseudotime')
sample_param.gamma_cds_list[[7]] <- orderCells(sample_param.gamma_cds_list[[7]], root_state = 2)

plot_cell_trajectory(sample_param.gamma_cds_list[[8]])
plot_cell_trajectory(sample_param.gamma_cds_list[[8]], color_by = 'Pseudotime')
sample_param.gamma_cds_list[[8]] <- orderCells(sample_param.gamma_cds_list[[8]], root_state = 2)

plot_cell_trajectory(sample_param.gamma_cds_list[[9]])
plot_cell_trajectory(sample_param.gamma_cds_list[[9]], color_by = 'Pseudotime')
sample_param.gamma_cds_list[[2]] <- orderCells(sample_param.gamma_cds_list[[9]], root_state = 2)

plot_cell_trajectory(sample_param.gamma_cds_list[[10]])
plot_cell_trajectory(sample_param.gamma_cds_list[[10]], color_by = 'Pseudotime')
sample_param.gamma_cds_list[[10]] <- orderCells(sample_param.gamma_cds_list[[10]], root_state = 2)

################################################################################################################################################################################
#check pseudotime for parameters sigma:
################################################################################################################################################################################
plot_cell_trajectory(sample_sigma_cds_list[[2]])
plot_cell_trajectory(sample_sigma_cds_list[[2]], color_by = 'Pseudotime')
sample_sigma_cds_list[[2]] <- orderCells(sample_sigma_cds_list[[2]], root_state = 2)

plot_cell_trajectory(sample_sigma_cds_list[[3]])
plot_cell_trajectory(sample_sigma_cds_list[[3]], color_by = 'Pseudotime')
sample_sigma_cds_list[[3]] <- orderCells(sample_sigma_cds_list[[3]], root_state = 2)

plot_cell_trajectory(sample_sigma_cds_list2[[7]])
plot_cell_trajectory(sample_sigma_cds_list2[[7]], color_by = 'Pseudotime')
sample_sigma_cds_list2[[7]] <- orderCells(sample_sigma_cds_list2[[7]], root_state = 2)

plot_cell_trajectory(sample_sigma_cds_list2[[8]])
plot_cell_trajectory(sample_sigma_cds_list2[[8]], color_by = 'Pseudotime')
sample_sigma_cds_list2[[8]] <- orderCells(sample_sigma_cds_list2[[8]], root_state = 2)

plot_cell_trajectory(sample_sigma_cds_list2[[9]])
plot_cell_trajectory(sample_sigma_cds_list2[[9]], color_by = 'Pseudotime')
sample_sigma_cds_list2[[9]] <- orderCells(sample_sigma_cds_list2[[9]], root_state = 2)

plot_cell_trajectory(sample_sigma_cds_list2[[10]])
plot_cell_trajectory(sample_sigma_cds_list2[[10]], color_by = 'Pseudotime')
sample_sigma_cds_list2[[10]] <- orderCells(sample_sigma_cds_list2[[10]], root_state = 2)
################################################################################################################################################################################
sample_dimension_cds_list <- sample_dimension_cds_list2
sample_initial_method_cds_list <- sample_initial_method_cds_list2
sample_maxIter_cds_list <- sample_maxIter_cds_list2
sample_sigma_cds_list <- sample_sigma_cds_list2
sample_lambda_cds_list <- sample_lambda_cds_list2
sample_ncenter_cds_list <- sample_ncenter_cds_list2

###############################################################################################################################################################################
# show some statistics to those parameters give pretty consistent results
# pseudotime consistency and branch consistency:

cal_benchmark_res <- function(cmbn_sets_list, cds = sample_dimension_cds_list) {
  mclapply(cmbn_sets_list, function(x, cds) {
    cds_1 <- cds[[x[[1]]]]
    cds_2 <- cds[[x[[2]]]]
    overlpa_cells <- intersect(colnames(cds_1), colnames(cds_2))

    if(length(unique(pData(cds_1)$State)) > 3){
      # cds_1 <- trimTree(cds_1, num_paths = 2, min_branch_thrsld = 0.1)
    }
    if(length(unique(pData(cds_2)$State)) > 3){
      # cds_2 <- trimTree(cds_2, num_paths = 2, min_branch_thrsld = 0.1)
    }

    overlpa_cells <- intersect(colnames(cds_1), colnames(cds_2))

    t_1 <- pData(cds_1[, overlpa_cells])$Pseudotime
    t_2 <- pData(cds_2[, overlpa_cells])$Pseudotime

    kendall_tau_cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    cor_res <- cor(t_1, t_2)

    #branch assignment:
    clusters_1 <- as.character(pData(cds_1[, overlpa_cells])$State)
    clusters_2 <- as.character(pData(cds_2[, overlpa_cells])$State)
    ClusteringMetrics_res <- calClusteringMetrics(as.character(clusters_1), as.character(clusters_2))

    return(list(cor = cor_res, kendall_tau = kendall_tau_cor_res, cluster = ClusteringMetrics_res))
  }, cds = cds, mc.cores = detectCores() - 2)
}

########################################################################################################################################
cmbn_sets <- expand.grid(A = 1:length(sample_dimension_cds_list), B = 1:length(sample_dimension_cds_list))
# cmbn_sets <- expand.grid(A = c(2, 3, 5, 9), B = c(2, 3, 5, 9))

cmbn_sets <- subset(cmbn_sets, A != B)
cmbn_sets_list <- split(cmbn_sets[, ], f = 1:nrow(cmbn_sets[, ]))
dimension_benchmark_res_list <- cal_benchmark_res(cmbn_sets_list, cds = sample_dimension_cds_list) #

########################################################################################################################################
cmbn_sets <- expand.grid(A = 1:length(sample_initial_method_cds_list), B = 1:length(sample_initial_method_cds_list))
# cmbn_sets <- expand.grid(A = c(2, 3, 5, 9), B = c(2, 3, 5, 9))

cmbn_sets <- subset(cmbn_sets, A != B)
cmbn_sets_list <- split(cmbn_sets[, ], f = 1:nrow(cmbn_sets[, ]))
initial_method_benchmark_res_list <- cal_benchmark_res(cmbn_sets_list, cds = sample_initial_method_cds_list)

########################################################################################################################################
maxIter_benchmark_res_list <- cal_benchmark_res(cmbn_sets_list, cds = sample_maxIter_cds_list)

########################################################################################################################################
cmbn_sets <- expand.grid(A = 1:length(sample_sigma_cds_list), B = 1:length(sample_sigma_cds_list))

cmbn_sets <- subset(cmbn_sets, A != B)
cmbn_sets_list <- split(cmbn_sets[, ], f = 1:nrow(cmbn_sets[, ]))
sigma_benchmark_res_list <- cal_benchmark_res(cmbn_sets_list, cds = sample_sigma_cds_list) #

########################################################################################################################################
cmbn_sets <- expand.grid(A = 1:length(sample_lambda_cds_list), B = 1:length(sample_lambda_cds_list))

cmbn_sets <- subset(cmbn_sets, A != B)
cmbn_sets_list <- split(cmbn_sets[, ], f = 1:nrow(cmbn_sets[, ]))
lambda_benchmark_res_list <- cal_benchmark_res(cmbn_sets_list, cds = sample_lambda_cds_list)

########################################################################################################################################
cmbn_sets <- expand.grid(A = 1:length(sample_ncenter_cds_list), B = 1:length(sample_ncenter_cds_list))

cmbn_sets <- subset(cmbn_sets, A != B)
cmbn_sets_list <- split(cmbn_sets[, ], f = 1:nrow(cmbn_sets[, ]))
ncenter_benchmark_res_list <- cal_benchmark_res(cmbn_sets_list, cds = sample_ncenter_cds_list)

########################################################################################################################################
cmbn_sets <- expand.grid(A = 1:length(sample_param.gamma_cds_list), B = 1:length(sample_param.gamma_cds_list))

cmbn_sets <- subset(cmbn_sets, A != B)
cmbn_sets_list <- split(cmbn_sets[, ], f = 1:nrow(cmbn_sets[, ]))
param.gamma_benchmark_res_list <- cal_benchmark_res(cmbn_sets_list, cds = sample_param.gamma_cds_list) #

#make the performance plot:
dimension_res_df <- data.frame(pearson.rho = unlist(lapply(dimension_benchmark_res_list, function(x) x$cor)),
                               kendall.tau = unlist(lapply(dimension_benchmark_res_list, function(x) x$kendall_tau)),
                              rand_ind = unlist(lapply(dimension_benchmark_res_list, function(x) x$cluster[1, 1])),
                              var_inf = unlist(lapply(dimension_benchmark_res_list, function(x) x$cluster[2, 1])),
                              adj_rand = unlist(lapply(dimension_benchmark_res_list, function(x) x$cluster[3, 1]))
)
initial_method_res_df <- data.frame(pearson.rho = unlist(lapply(initial_method_benchmark_res_list, function(x) x$cor)),
                                    kendall.tau = unlist(lapply(initial_method_benchmark_res_list, function(x) x$kendall_tau)),
                               rand_ind = unlist(lapply(initial_method_benchmark_res_list, function(x) x$cluster[1, 1])),
                               var_inf = unlist(lapply(initial_method_benchmark_res_list, function(x) x$cluster[2, 1])),
                               adj_rand = unlist(lapply(initial_method_benchmark_res_list, function(x) x$cluster[3, 1]))
)
maxIter_res_df <- data.frame(pearson.rho = unlist(lapply(maxIter_benchmark_res_list, function(x) x$cor)),
                             kendall.tau = unlist(lapply(maxIter_benchmark_res_list, function(x) x$kendall_tau)),
                               rand_ind = unlist(lapply(maxIter_benchmark_res_list, function(x) x$cluster[1, 1])),
                               var_inf = unlist(lapply(maxIter_benchmark_res_list, function(x) x$cluster[2, 1])),
                               adj_rand = unlist(lapply(maxIter_benchmark_res_list, function(x) x$cluster[3, 1]))
)
sigma_res_df <- data.frame(pearson.rho = unlist(lapply(sigma_benchmark_res_list, function(x) x$cor)),
                           kendall.tau = unlist(lapply(sigma_benchmark_res_list, function(x) x$kendall_tau)),
                               rand_ind = unlist(lapply(sigma_benchmark_res_list, function(x) x$cluster[1, 1])),
                               var_inf = unlist(lapply(sigma_benchmark_res_list, function(x) x$cluster[2, 1])),
                               adj_rand = unlist(lapply(sigma_benchmark_res_list, function(x) x$cluster[3, 1]))
)
lambda_res_df <- data.frame(pearson.rho = unlist(lapply(lambda_benchmark_res_list, function(x) x$cor)),
                            kendall.tau = unlist(lapply(lambda_benchmark_res_list, function(x) x$kendall_tau)),
                               rand_ind = unlist(lapply(lambda_benchmark_res_list, function(x) x$cluster[1, 1])),
                               var_inf = unlist(lapply(lambda_benchmark_res_list, function(x) x$cluster[2, 1])),
                               adj_rand = unlist(lapply(lambda_benchmark_res_list, function(x) x$cluster[3, 1]))
)
ncenter_res_df <- data.frame(pearson.rho = unlist(lapply(ncenter_benchmark_res_list, function(x) x$cor)),
                             kendall.tau = unlist(lapply(ncenter_benchmark_res_list, function(x) x$kendall_tau)),
                               rand_ind = unlist(lapply(ncenter_benchmark_res_list, function(x) x$cluster[1, 1])),
                               var_inf = unlist(lapply(ncenter_benchmark_res_list, function(x) x$cluster[2, 1])),
                               adj_rand = unlist(lapply(ncenter_benchmark_res_list, function(x) x$cluster[3, 1]))
)
param.gamma_res_df <- data.frame(pearson.rho = unlist(lapply(param.gamma_benchmark_res_list, function(x) x$cor)),
                                 kendall.tau = unlist(lapply(param.gamma_benchmark_res_list, function(x) x$kendall_tau)),
                               rand_ind = unlist(lapply(param.gamma_benchmark_res_list, function(x) x$cluster[1, 1])),
                               var_inf = unlist(lapply(param.gamma_benchmark_res_list, function(x) x$cluster[2, 1])),
                               adj_rand = unlist(lapply(param.gamma_benchmark_res_list, function(x) x$cluster[3, 1]))
)

all_res_df <- Reduce(rbind, list(dimension_res_df, initial_method_res_df, maxIter_res_df, sigma_res_df, lambda_res_df, ncenter_res_df, param.gamma_res_df)) #
# all_res_df <- Reduce(rbind, list(dimension_res_df, sigma_res_df, param.gamma_res_df))
all_res_df$Type <- c(rep('Dimension', nrow(dimension_res_df)),
                     rep('initial method', nrow(initial_method_res_df)),
                     rep('maxIter', nrow(maxIter_res_df)),
                     rep('sigma', nrow(sigma_res_df)),
                     rep('lambda', nrow(lambda_res_df)),
                      rep('ncenter', nrow(ncenter_res_df)),
                     rep('param.gamma', nrow(param.gamma_res_df)))

mlt_all_res_df <- melt(all_res_df, id=c('Type'))
pdf(paste(SI_fig_dir, 'mar_seq_all_DDRTree_params_perform.pdf', sep = ''), height = 1.5, width = 2.3)
qplot(Type, abs(value), facets = "~variable", data = subset(mlt_all_res_df, variable %in% c("kendall.tau", "adj_rand")), size = 0.1,
      color = Type, geom = c("jitter", "boxplot"), alpha = I(0.7)) +  scale_size(range = c(0.1, 0.5)) + nm_theme() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + monocle:::monocle_theme_opts() + ylim(0, 1)
dev.off()

########################################################################################################################################################
# plot progressive figures: 
########################################################################################################################################################
pairwise_cal_benchmark_res <- function(cds_1, cds_2) {

  overlpa_cells <- intersect(colnames(cds_1), colnames(cds_2))

  if(length(unique(pData(cds_1)$State)) > 3){
    #cds_1 <- trimTree(cds_1, num_paths = 2, min_branch_thrsld = 0.1)
  }
  if(length(unique(pData(cds_2)$State)) > 3){
    #cds_2 <- trimTree(cds_2, num_paths = 2, min_branch_thrsld = 0.1)
  }

  overlpa_cells <- intersect(colnames(cds_1), colnames(cds_2))

  t_1 <- pData(cds_1[, overlpa_cells])$Pseudotime
  t_2 <- pData(cds_2[, overlpa_cells])$Pseudotime

  cor_res <- cor(t_1, t_2)
  kendall_cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
  
  #branch assignment:
  clusters_1 <- as.character(pData(cds_1[, overlpa_cells])$State)
  clusters_2 <- as.character(pData(cds_2[, overlpa_cells])$State)
  ClusteringMetrics_res <- calClusteringMetrics(as.character(clusters_1), as.character(clusters_2))

  return(list(cor = cor_res, kendall_tau = kendall_cor_res, cluster = ClusteringMetrics_res))
}

dimension_benchmark_res <- lapply(sample_dimension_cds_list, function(x) pairwise_cal_benchmark_res(x, absolute_cds)) #
initial_method_benchmark_res <- lapply(sample_initial_method_cds_list, function(x) pairwise_cal_benchmark_res(x, absolute_cds))
maxIter_benchmark_res <- lapply(sample_maxIter_cds_list, function(x) pairwise_cal_benchmark_res(x, absolute_cds))
sigma_benchmark_res <- lapply(sample_sigma_cds_list, function(x) pairwise_cal_benchmark_res(x, absolute_cds)) #
lambda_benchmark_res <- lapply(sample_lambda_cds_list, function(x) pairwise_cal_benchmark_res(x, absolute_cds))
ncenter_benchmark_res <- lapply(sample_ncenter_cds_list, function(x) pairwise_cal_benchmark_res(x, absolute_cds))
param.gamma_benchmark_res <- lapply(sample_param.gamma_cds_list, function(x) pairwise_cal_benchmark_res(x, absolute_cds)) #
# 
# #make the performance plot: 
dimension_res_df <- data.frame(pearson.rho = unlist(lapply(dimension_benchmark_res, function(x) x$cor)),
                               kendall.tau = unlist(lapply(dimension_benchmark_res, function(x) x$kendall_tau)),
                               rand_ind = unlist(lapply(dimension_benchmark_res, function(x) x$cluster[1, 1])),
                               var_inf = unlist(lapply(dimension_benchmark_res, function(x) x$cluster[2, 1])),
                               adj_rand = unlist(lapply(dimension_benchmark_res, function(x) x$cluster[3, 1]))
)
initial_method_res_df <- data.frame(pearson.rho = unlist(lapply(initial_method_benchmark_res, function(x) x$cor)),
                                    kendall.tau = unlist(lapply(initial_method_benchmark_res, function(x) x$kendall_tau)),
                               rand_ind = unlist(lapply(initial_method_benchmark_res, function(x) x$cluster[1, 1])),
                               var_inf = unlist(lapply(initial_method_benchmark_res, function(x) x$cluster[2, 1])),
                               adj_rand = unlist(lapply(initial_method_benchmark_res, function(x) x$cluster[3, 1]))
)
maxIter_res_df <- data.frame(pearson.rho = unlist(lapply(maxIter_benchmark_res, function(x) x$cor)),
                             kendall.tau = unlist(lapply(maxIter_benchmark_res, function(x) x$kendall_tau)),
                             rand_ind = unlist(lapply(maxIter_benchmark_res, function(x) x$cluster[1, 1])),
                             var_inf = unlist(lapply(maxIter_benchmark_res, function(x) x$cluster[2, 1])),
                             adj_rand = unlist(lapply(maxIter_benchmark_res, function(x) x$cluster[3, 1]))
)
sigma_res_df <- data.frame(pearson.rho = unlist(lapply(sigma_benchmark_res, function(x) x$cor)),
                           kendall.tau = unlist(lapply(sigma_benchmark_res, function(x) x$kendall_tau)),
                           rand_ind = unlist(lapply(sigma_benchmark_res, function(x) x$cluster[1, 1])),
                           var_inf = unlist(lapply(sigma_benchmark_res, function(x) x$cluster[2, 1])),
                           adj_rand = unlist(lapply(sigma_benchmark_res, function(x) x$cluster[3, 1]))
)
lambda_res_df <- data.frame(pearson.rho = unlist(lapply(lambda_benchmark_res, function(x) x$cor)),
                            kendall.tau = unlist(lapply(lambda_benchmark_res, function(x) x$kendall_tau)),
                            rand_ind = unlist(lapply(lambda_benchmark_res, function(x) x$cluster[1, 1])),
                            var_inf = unlist(lapply(lambda_benchmark_res, function(x) x$cluster[2, 1])),
                            adj_rand = unlist(lapply(lambda_benchmark_res, function(x) x$cluster[3, 1]))
)
ncenter_res_df <- data.frame(pearson.rho = unlist(lapply(ncenter_benchmark_res, function(x) x$cor)),
                             kendall.tau = unlist(lapply(ncenter_benchmark_res, function(x) x$kendall_tau)),
                             rand_ind = unlist(lapply(ncenter_benchmark_res, function(x) x$cluster[1, 1])),
                             var_inf = unlist(lapply(ncenter_benchmark_res, function(x) x$cluster[2, 1])),
                             adj_rand = unlist(lapply(ncenter_benchmark_res, function(x) x$cluster[3, 1]))
)
param.gamma_res_df <- data.frame(pearson.rho = unlist(lapply(param.gamma_benchmark_res, function(x) x$cor)),
                                 kendall.tau = unlist(lapply(param.gamma_benchmark_res, function(x) x$kendall_tau)),
                                 rand_ind = unlist(lapply(param.gamma_benchmark_res, function(x) x$cluster[1, 1])),
                                 var_inf = unlist(lapply(param.gamma_benchmark_res, function(x) x$cluster[2, 1])),
                                 adj_rand = unlist(lapply(param.gamma_benchmark_res, function(x) x$cluster[3, 1]))
)

all_res_df <- Reduce(rbind, list(dimension_res_df, initial_method_res_df, lambda_res_df, maxIter_res_df, param.gamma_res_df, sigma_res_df, ncenter_res_df)) #
# all_res_df <- Reduce(rbind, list(dimension_res_df, sigma_res_df, param.gamma_res_df))
all_res_df$Type <- c(rep('Dimension', nrow(dimension_res_df)), #9
                     rep('initial method', nrow(initial_method_res_df)), #5
                     rep('lambda', nrow(lambda_res_df)), #10
                     rep('maxIter', nrow(maxIter_res_df)), #19
                     rep('param.gamma', nrow(param.gamma_res_df)),  #20
                     rep('sigma', nrow(sigma_res_df)), #20
                     rep('ncenter', nrow(ncenter_res_df))) #17

all_res_df$param_vals <- c(dimensions,
                       1:length(initial_method_res_df), #c("ICA", "PCA", "LLE", "DM", "destiny_diffusionMaps"),
                       round(lambda, 1),
                       maxIter,
                       round(param.gamma, 1),
                       round(sigma, 1),
                       ncenter
                       )
#mlt_all_res_df <- melt(all_res_df, id=c('Type', 'values'))

mlt_all_res_df <- melt(all_res_df, id=c('Type', 'param_vals'))

mlt_all_res_df$variable <- revalue(as.character(mlt_all_res_df$variable), 
        c("kendall.tau" = "Kendall'tau", "adj_rand" = 'Adjusted Rand index', "pearson.rho" = "Pearson's rho"))

pdf(paste(SI_fig_dir, benchmark_type, '_progressive_parameters.pdf', sep = ''), height = 3, width = 6)
qplot(param_vals, abs(value), data = subset(mlt_all_res_df, variable %in% c("Kendall'tau", "Adjusted Rand index")), size = 0.1, 
      color = Type, geom = c("point", 'line'), alpha = I(0.7)) + facet_wrap(c("variable", "Type"), scales = 'free', nrow = 2) +  scale_size(range = c(0.1, 0.5)) + nm_theme() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + monocle:::monocle_theme_opts() + ylim(0, 1) + xlab('Parameters') + ylab('')
dev.off()

# pdf(paste(SI_fig_dir, benchmark_type, '_progressive_parameters2.pdf', sep = ''), height = 3, width = 12)
# qplot(param_vals, abs(value), data = subset(mlt_all_res_df, variable %in% c("kendall.tau", "adj_rand")), size = 0.1, 
#       color = Type, geom = c("point", 'line'), alpha = I(0.7)) + facet_grid(c("variable", "Type"), scales = 'free', nrow = 2) +  scale_size(range = c(0.1, 0.5)) + nm_theme() +
#   theme(axis.text.x = element_text(angle = 30, hjust = 1)) + monocle:::monocle_theme_opts() + ylim(0, 1)
# dev.off()

# mlt_all_res_df$Type <- factor(mlt_all_res_df$Type, levels = c())
pdf(paste(SI_fig_dir, benchmark_type, '_progressive_parameters3.pdf', sep = ''), height = 3, width = 6)
qplot(param_vals, abs(value), data = subset(mlt_all_res_df, variable %in% c("Kendall'tau", "Adjusted Rand index") & Type != 'initial method' ), size = 0.1, 
      color = Type, geom = c("point", 'line'), alpha = I(0.7)) + facet_grid(variable ~ Type, scales = 'free') +  scale_size(range = c(0.1, 0.5)) + nm_theme() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + monocle:::monocle_theme_opts() + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, '_progressive_parameters3.1.pdf', sep = ''), height = 3, width = 6)
qplot(param_vals, abs(value), data = subset(mlt_all_res_df, variable %in% c("Pearson's rho", "Adjusted Rand index") & Type != 'initial method' ), size = 0.1, 
      color = Type, geom = c("point", 'line'), alpha = I(0.7)) + facet_grid(variable ~ Type, scales = 'free') +  scale_size(range = c(0.1, 0.5)) + nm_theme() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + monocle:::monocle_theme_opts() + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, '_progressive_parameters3.1.pdf', sep = ''), height = 1, width = 6)
qplot(param_vals, abs(value), data = subset(mlt_all_res_df, variable %in% c("Pearson's rho", "Adjusted Rand index") & Type != 'initial method' ), size = 0.1, 
      color = variable, geom = c("point", 'line'), alpha = I(0.7)) + facet_grid(~Type, scales = 'free') +  scale_size(range = c(0.1, 0.5)) + nm_theme() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + monocle:::monocle_theme_opts() + ylim(0, 1) + xlab('Parameters value')+ ylab('')
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, '_progressive_parameters4.pdf', sep = ''), height = 3, width = 6)
qplot(param_vals, abs(value), data = subset(mlt_all_res_df, variable %in% c("Kendall'tau", "Adjusted Rand index")), size = 0.1, 
      color = Type, geom = c("point", 'line'), alpha = I(0.7)) + facet_grid(variable ~ Type, scales = 'free') +  scale_size(range = c(0.1, 0.5)) + nm_theme() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + monocle:::monocle_theme_opts() + ylim(0, 1)
dev.off()

#######################################################################################################################################################
# save the data 
#######################################################################################################################################################
save.image(paste('./RData/', benchmark_type, 'analysis_DDRTree_param.RData', sep = ''))

