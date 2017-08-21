############################################################################################################################################################
# determine the root_state, AT1_state, AT2_state for each software based on the ordering from the full dataset
############################################################################################################################################################
# monocle 2 state assignment
monocle2_p1 <- plot_cell_trajectory(cds_downsampled_cells_ordered[[36]])
monocle2_p2 <- plot_cell_trajectory(cds_downsampled_cells_ordered[[36]], color_by = 'Time')
monocle2_p3 <- plot_cell_trajectory(cds_downsampled_cells_ordered[[36]], color_by = 'Pseudotime')
multiplot(monocle2_p1, monocle2_p2, monocle2_p3)
monocle2_root_state <- row.names(subset(pData(cds_downsampled_cells_ordered[[36]]), State == 3))
monocle2_AT1_state <- row.names(subset(pData(cds_downsampled_cells_ordered[[36]]), State == 1)) #this should be actually state 2 for history reasons
monocle2_AT2_state <- row.names(subset(pData(cds_downsampled_cells_ordered[[36]]), State == 2))

# monocle 1 state assignment
monocle1_p1 <- plot_cell_trajectory(ICA_cds_downsampled_cells_ordered[[36]])
monocle1_p2 <- plot_cell_trajectory(ICA_cds_downsampled_cells_ordered[[36]], color_by = 'Time')
monocle1_p3 <- plot_cell_trajectory(ICA_cds_downsampled_cells_ordered[[36]], color_by = 'Pseudotime')
multiplot(monocle1_p1, monocle1_p2, monocle1_p3)
monocle1_root_state <- row.names(subset(pData(ICA_cds_downsampled_cells_ordered[[36]]), State == 1))
monocle1_AT1_state <- row.names(subset(pData(ICA_cds_downsampled_cells_ordered[[36]]), State == 2))
monocle1_AT2_state <- row.names(subset(pData(ICA_cds_downsampled_cells_ordered[[36]]), State == 3))

# dpt state assignment
dpt_p1 <- qplot(dpt_cds_downsampled_cells_ordered[[36]]$dm$DC1, dpt_cds_downsampled_cells_ordered[[36]]$dm$DC2, color = as.character(dpt_cds_downsampled_cells_ordered[[36]]$branch[, 1]))
dpt_p2 <- qplot(dpt_cds_downsampled_cells_ordered[[36]]$dm$DC1, dpt_cds_downsampled_cells_ordered[[36]]$dm$DC2, color = pData(cds_downsampled_cells_ordered[[36]])$Time)
dpt_p3 <- qplot(dpt_cds_downsampled_cells_ordered[[36]]$dm$DC1, dpt_cds_downsampled_cells_ordered[[36]]$dm$DC2, color = dpt_cds_downsampled_cells_ordered[[36]]$pt)
multiplot(dpt_p1, dpt_p2, dpt_p3)
dpt_root_state <- row.names(dpt_cds_downsampled_cells_ordered[[36]]$branch)[which(dpt_cds_downsampled_cells_ordered[[36]]$branch[, 1] == 2)]
dpt_AT1_state <- row.names(dpt_cds_downsampled_cells_ordered[[36]]$branch)[which(dpt_cds_downsampled_cells_ordered[[36]]$branch[, 1] == 3)]
dpt_AT2_state <- row.names(dpt_cds_downsampled_cells_ordered[[36]]$branch)[which(dpt_cds_downsampled_cells_ordered[[36]]$branch[, 1] == 1)]

# wishbone state assignment (this part requires re-accessment)
wishbone_p1 <- qplot(dpt_cds_downsampled_cells_ordered[[36]]$dm$DC1, dpt_cds_downsampled_cells_ordered[[36]]$dm$DC2, color = as.character(subset(fraction_wishbone_res, run == 36)[, 'branch']))
wishbone_p2 <- qplot(subset(wishbone_res, run == 36)[, 'dm1'], subset(wishbone_res, run == 36)[, 'dm2'], color = pData(cds_downsampled_cells_ordered[[36]])$Time)
wishbone_p3 <- qplot(subset(wishbone_res, run == 36)[, 'dm1'], subset(wishbone_res, run == 36)[, 'dm2'], color = subset(wishbone_res, run)[, 'trajectory'])
multiplot(wishbone_p1, wishbone_p2, wishbone_p3)
wishbone_root_state <- row.names(subset(ICA_cds_downsampled_cells_ordered[[36]], State == 1))
wishbone_AT1_state <- row.names(subset(ICA_cds_downsampled_cells_ordered[[36]], State == 2))
wishbone_AT2_state <- row.names(subset(ICA_cds_downsampled_cells_ordered[[36]], State == 3))

############################################################################################################################################################
#Robustness of tree construction under downsampling from 0.1 to 1
############################################################################################################################################################

cell_sampling_cmbn_sets <- expand.grid(1:36, 36)

cell_sampling_cmbn_sets_list <- split(cell_sampling_cmbn_sets[, ], f = 1:nrow(cell_sampling_cmbn_sets[, ]))
cell_sampling_benchmark_res_list <- mclapply(cell_sampling_cmbn_sets_list, function(x) {
  cds_1 <- cds_downsampled_cells_ordered[[x[[1]]]]
  cds_2 <- cds_downsampled_cells_ordered[[x[[2]]]]

  if(length(unique(pData(cds_1)$State)) > 3){
    cds_1 <- trimTree(cds_1, num_paths = 2, min_branch_thrsld = 0.1)
  }
  if(length(unique(pData(cds_2)$State)) > 3){
    cds_2 <- trimTree(cds_2, num_paths = 2, min_branch_thrsld = 0.1)
  }

  overlap_cells <- intersect(colnames(cds_1), colnames(cds_2))

  if(length(overlap_cells)){
    t_1 <- pData(cds_1[, overlap_cells])$Pseudotime
    t_2 <- pData(cds_2[, overlap_cells])$Pseudotime

    cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    pearson_cor_res <- cor(t_1, t_2)
    print(AT1_state)
    AT1_lineage <- c(monocle2_root_state, monocle2_AT1_state)
    AT2_lineage <- c(monocle2_root_state, monocle2_AT2_state)
    AT1_t_1 <- pData(cds_1[, overlap_cells[overlap_cells %in% AT1_state]])$Pseudotime
    AT1_t_2 <- pData(cds_2[, overlap_cells[overlap_cells %in% AT1_state]])$Pseudotime
    AT2_t_1 <- pData(cds_1[, overlap_cells[overlap_cells %in% AT2_state]])$Pseudotime
    AT2_t_2 <- pData(cds_2[, overlap_cells[overlap_cells %in% AT2_state]])$Pseudotime

    AT1_cor_res <- cor(AT1_t_1, AT1_t_2, method = 'kendall', use = 'pairwise.complete.obs')
    AT1_pearson_cor_res <- cor(AT1_t_1, AT1_t_2)
    AT2_cor_res <- cor(AT2_t_1, AT2_t_2, method = 'kendall', use = 'pairwise.complete.obs')
    AT2_pearson_cor_res <- cor(AT2_t_1, AT2_t_2)
    
    #branch assignment:
    clusters_1 <- pData(cds_1[, overlap_cells])$State
    clusters_2 <- pData(cds_2[, overlap_cells])$State
    ClusteringMetrics_res <- calClusteringMetrics(as.character(clusters_1), as.character(clusters_2))

    return(list(cor = c(cor_res, AT1_cor_res, AT2_cor_res), pearson_cor = c(pearson_cor_res, AT1_pearson_cor_res, AT2_pearson_cor_res), cluster = ClusteringMetrics_res))
  }
  else
    return(list(cor = rep(NA, 3), pearson_cor = rep(NA, 3), cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
}, mc.cores = detectCores() - 2)

############################################################################################################################################################
#Robustness of tree construction under 0.8 downsampling rate
############################################################################################################################################################

#robustness of the ordering algorithm based on the DPT paper:
cmbn_sets <- expand.grid(1:repeat_downsampling_num, 1:repeat_downsampling_num)

cmbn_sets_list <- split(cmbn_sets[, ], f = 1:nrow(cmbn_sets[, ]))
benchmark_res_list <- mclapply(cmbn_sets_list, function(x) {
  cds_1 <- cds_downsampled_cells_ordered_0.8[[x[[1]]]]
  cds_2 <- cds_downsampled_cells_ordered_0.8[[x[[2]]]]

  if(length(unique(pData(cds_1)$State)) > 3){
    cds_1 <- trimTree(cds_1, num_paths = 2, min_branch_thrsld = 0.1)
  }
  if(length(unique(pData(cds_2)$State)) > 3){
    cds_2 <- trimTree(cds_2, num_paths = 2, min_branch_thrsld = 0.1)
  }

  overlap_cells <- intersect(colnames(cds_1), colnames(cds_2))

  if(length(overlap_cells)) {
    t_1 <- pData(cds_1[, overlap_cells])$Pseudotime
    t_2 <- pData(cds_2[, overlap_cells])$Pseudotime

    cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    pearson_cor_res <- cor(t_1, t_2)
    
    AT1_lineage <- c(monocle2_root_state, monocle2_AT1_state)
    AT2_lineage <- c(monocle2_root_state, monocle2_AT2_state)
    AT1_t_1 <- pData(cds_1[, overlap_cells[overlap_cells %in% AT1_state]])$Pseudotime
    AT1_t_2 <- pData(cds_2[, overlap_cells[overlap_cells %in% AT1_state]])$Pseudotime
    AT2_t_1 <- pData(cds_1[, overlap_cells[overlap_cells %in% AT2_state]])$Pseudotime
    AT2_t_2 <- pData(cds_2[, overlap_cells[overlap_cells %in% AT2_state]])$Pseudotime

    AT1_cor_res <- cor(AT1_t_1, AT1_t_2, method = 'kendall', use = 'pairwise.complete.obs')
    AT1_pearson_cor_res <- cor(AT1_t_1, AT1_t_2)
    AT2_cor_res <- cor(AT2_t_1, AT2_t_2, method = 'kendall', use = 'pairwise.complete.obs')
    AT2_pearson_cor_res <- cor(AT2_t_1, AT2_t_2)
    #branch assignment:
    clusters_1 <- pData(cds_1[, overlap_cells])$State
    clusters_2 <- pData(cds_2[, overlap_cells])$State
    ClusteringMetrics_res <- calClusteringMetrics(as.character(clusters_1), as.character(clusters_2))

    return(list(cor = c(cor_res, AT1_cor_res, AT2_cor_res), pearson_cor = c(pearson_cor_res, AT1_pearson_cor_res, AT2_pearson_cor_res), cluster = ClusteringMetrics_res))
  }
  else
    return(list(cor = rep(NA, 3), pearson_cor = rep(NA, 3), cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
}, mc.cores = detectCores() - 2)

#plot the results:
############################################################################################################################################################
sampling_res_df <- data.frame(kendall.tau = unlist(lapply(benchmark_res_list, function(x) mean(x$cor[2:3]))),
                              pearson_rho = unlist(lapply(benchmark_res_list, function(x) mean(x$pearson_cor[2:3]))),
                                  rand_ind = unlist(lapply(benchmark_res_list, function(x) x$cluster[1, 1])),
                                  var_inf = unlist(lapply(benchmark_res_list, function(x) x$cluster[2, 1])),
                                  adj_rand = unlist(lapply(benchmark_res_list, function(x) x$cluster[3, 1]))
)

############################################################################################################################################################
cell_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(cell_sampling_benchmark_res_list, function(x) mean(x$cor[2:3]))),
                              pearson_rho = unlist(lapply(cell_sampling_benchmark_res_list, function(x) mean(x$pearson_cor[2:3]))),
                              rand_ind = unlist(lapply(cell_sampling_benchmark_res_list, function(x) x$cluster[1, 1])),
                              var_inf = unlist(lapply(cell_sampling_benchmark_res_list, function(x) x$cluster[2, 1])),
                              adj_rand = unlist(lapply(cell_sampling_benchmark_res_list, function(x) x$cluster[3, 1]))
)

progressive_downsampling_rng <- (nrow(cell_sampling_res_df) - 35):nrow(cell_sampling_res_df)
valid_cell_sampling_res_df <- cell_sampling_res_df[progressive_downsampling_rng, ]


#ICA
ICA_cell_sampling_benchmark_res_list <- mclapply(cell_sampling_cmbn_sets_list, function(x) {
  message(x)
  cds_1 <- ICA_cds_downsampled_cells_ordered[[x[[1]]]]
  cds_2 <- ICA_cds_downsampled_cells_ordered[[x[[2]]]]
  
  # if(length(pData(cds_1)$State) > 3){
  #   cds_1 <- trimTree(cds_1, num_paths = 2, min_branch_thrsld = 0.1)
  # }
  # if(length(pData(cds_2)$State) > 3){
  #   cds_2 <- trimTree(cds_2, num_paths = 2, min_branch_thrsld = 0.1)
  # }
  
  overlap_cells <- intersect(colnames(cds_1), colnames(cds_2))
  
  if(length(overlap_cells)) {
    t_1 <- pData(cds_1[, overlap_cells])$Pseudotime
    t_2 <- pData(cds_2[, overlap_cells])$Pseudotime
    
    cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    pearson_cor_res <- cor(t_1, t_2)
    AT1_lineage <- c(monocle1_root_state, monocle1_AT1_state)
    AT2_lineage <- c(monocle1_root_state, monocle1_AT2_state)
    AT1_t_1 <- pData(cds_1[, overlap_cells[overlap_cells %in% AT1_state]])$Pseudotime
    AT1_t_2 <- pData(cds_2[, overlap_cells[overlap_cells %in% AT1_state]])$Pseudotime
    AT2_t_1 <- pData(cds_1[, overlap_cells[overlap_cells %in% AT2_state]])$Pseudotime
    AT2_t_2 <- pData(cds_2[, overlap_cells[overlap_cells %in% AT2_state]])$Pseudotime

    AT1_cor_res <- cor(AT1_t_1, AT1_t_2, method = 'kendall', use = 'pairwise.complete.obs')
    AT1_pearson_cor_res <- cor(AT1_t_1, AT1_t_2)
    AT2_cor_res <- cor(AT2_t_1, AT2_t_2, method = 'kendall', use = 'pairwise.complete.obs')
    AT2_pearson_cor_res <- cor(AT2_t_1, AT2_t_2)

    #branch assignment:
    clusters_1 <- pData(cds_1[, overlap_cells])$State
    clusters_2 <- pData(cds_2[, overlap_cells])$State
    ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
    
    return(list(cor = c(cor_res, AT1_cor_res, AT2_cor_res), pearson_cor = c(pearson_cor_res, AT1_pearson_cor_res, AT2_pearson_cor_res), cluster = ClusteringMetrics_res))
  }
  else
    return(list(cor = rep(NA, 3), pearson_cor = rep(NA, 3), cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
}, mc.cores = detectCores() - 2)

#ICA
ICA_benchmark_res_list <- mclapply(cmbn_sets_list, function(x) {
  cds_1 <- ICA_cds_downsampled_cells_ordered_0.8[[x[[1]]]]
  cds_2 <- ICA_cds_downsampled_cells_ordered_0.8[[x[[2]]]]
  
  # if(length(pData(cds_1)$State) > 3){
  #   cds_1 <- trimTree(cds_1, num_paths = 2, min_branch_thrsld = 0.1)
  # }
  # if(length(pData(cds_2)$State) > 3){
  #   cds_2 <- trimTree(cds_2, num_paths = 2, min_branch_thrsld = 0.1)
  # }
  
  overlap_cells <- intersect(colnames(cds_1), colnames(cds_2))
  
  if(length(overlap_cells)) {
    t_1 <- pData(cds_1[, overlap_cells])$Pseudotime
    t_2 <- pData(cds_2[, overlap_cells])$Pseudotime
    
    cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    pearson_cor_res <- cor(t_1, t_2)
    
    AT1_lineage <- c(monocle1_root_state, monocle1_AT1_state)
    AT2_lineage <- c(monocle1_root_state, monocle1_AT2_state)
    AT1_t_1 <- pData(cds_1[, overlap_cells[overlap_cells %in% AT1_state]])$Pseudotime
    AT1_t_2 <- pData(cds_2[, overlap_cells[overlap_cells %in% AT1_state]])$Pseudotime
    AT2_t_1 <- pData(cds_1[, overlap_cells[overlap_cells %in% AT2_state]])$Pseudotime
    AT2_t_2 <- pData(cds_2[, overlap_cells[overlap_cells %in% AT2_state]])$Pseudotime

    AT1_cor_res <- cor(AT1_t_1, AT1_t_2, method = 'kendall', use = 'pairwise.complete.obs')
    AT1_pearson_cor_res <- cor(AT1_t_1, AT1_t_2)
    AT2_cor_res <- cor(AT2_t_1, AT2_t_2, method = 'kendall', use = 'pairwise.complete.obs')
    AT2_pearson_cor_res <- cor(AT2_t_1, AT2_t_2)

    #branch assignment:
    clusters_1 <- pData(cds_1[, overlap_cells])$State
    clusters_2 <- pData(cds_2[, overlap_cells])$State
    ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)
    
    return(list(cor = c(cor_res, AT1_cor_res, AT2_cor_res), pearson_cor = c(pearson_cor_res, AT1_pearson_cor_res, AT2_pearson_cor_res), cluster = ClusteringMetrics_res))
  }
  else
    return(list(cor = rep(NA, 3), pearson_cor = rep(NA, 3), cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
}, mc.cores = detectCores() - 2)

ICA_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(ICA_benchmark_res_list, function(x) mean(x$cor[2:3]))),
                                  pearson_rho = unlist(lapply(ICA_benchmark_res_list, function(x) mean(x$pearson_cor[2:3]))),
                                  rand_ind = unlist(lapply(ICA_benchmark_res_list, function(x) x$cluster[1, 1])),
                                  var_inf = unlist(lapply(ICA_benchmark_res_list, function(x) x$cluster[2, 1])),
                                  adj_rand = unlist(lapply(ICA_benchmark_res_list, function(x) x$cluster[3, 1]))
)

ICA_cell_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(ICA_cell_sampling_benchmark_res_list, function(x) mean(x$cor[2:3]))),
                                       pearson_rho = unlist(lapply(ICA_cell_sampling_benchmark_res_list, function(x) mean(x$pearson_cor[2:3]))),
                                       rand_ind = unlist(lapply(ICA_cell_sampling_benchmark_res_list, function(x) x$cluster[1, 1])),
                                       var_inf = unlist(lapply(ICA_cell_sampling_benchmark_res_list, function(x) x$cluster[2, 1])),
                                       adj_rand = unlist(lapply(ICA_cell_sampling_benchmark_res_list, function(x) x$cluster[3, 1]))
)

progressive_downsampling_rng <- (nrow(cell_sampling_res_df) - 35):nrow(cell_sampling_res_df)
# DDRTree_valid_cell_sampling_res_df <- DDRTree_cell_sampling_res_df[1261:1296, ]
ICA_valid_cell_sampling_res_df <- ICA_cell_sampling_res_df[progressive_downsampling_rng, ]


############################################################################################################################################################
#generate the data.frame for benchmarking dpt results
############################################################################################################################################################
dpt_cell_sampling_cmbn_sets <- expand.grid(1:36, 36)
dpt_cell_sampling_cmbn_sets_list <- split(dpt_cell_sampling_cmbn_sets[, ], f = 1:nrow(dpt_cell_sampling_cmbn_sets[, ]))

max_root_ind <- c()
dpt_cell_sampling_benchmark_res_list <- mclapply(dpt_cell_sampling_cmbn_sets_list, function(x) {
  cds_1 <- cds_downsampled_cells_ordered[[x[[1]]]]
  cds_2 <- cds_downsampled_cells_ordered[[x[[2]]]]
  overlap_cells <- intersect(colnames(cds_1), colnames(cds_2))

  if(length(overlap_cells)) {
    t_1 <- dpt_cds_downsampled_cells_ordered[[x[[1]]]]$pt[overlap_cells]#[overlpa_cells, 'DPT']
    t_2 <- dpt_cds_downsampled_cells_ordered[[x[[2]]]]$pt[overlap_cells]#[overlpa_cells, 'DPT']

    cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    pearson_cor_res <- cor(t_1, t_2)

    AT1_lineage <- c(dpt_root_state, dpt_AT1_state)
    AT2_lineage <- c(dpt_root_state, dpt_AT2_state)
    AT1_t_1 <- dpt_cds_downsampled_cells_ordered[[x[[1]]]]$pt[overlap_cells[overlap_cells %in% AT1_state]]
    AT1_t_2 <- dpt_cds_downsampled_cells_ordered[[x[[2]]]]$pt[overlap_cells[overlap_cells %in% AT1_state]]
    AT2_t_1 <- dpt_cds_downsampled_cells_ordered[[x[[1]]]]$pt[overlap_cells[overlap_cells %in% AT2_state]]
    AT2_t_2 <- dpt_cds_downsampled_cells_ordered[[x[[2]]]]$pt[overlap_cells[overlap_cells %in% AT2_state]]

    AT1_cor_res <- cor(AT1_t_1, AT1_t_2, method = 'kendall', use = 'pairwise.complete.obs')
    AT1_pearson_cor_res <- cor(AT1_t_1, AT1_t_2)
    AT2_cor_res <- cor(AT2_t_1, AT2_t_2, method = 'kendall', use = 'pairwise.complete.obs')
    AT2_pearson_cor_res <- cor(AT2_t_1, AT2_t_2)
    
    #branch assignment:
    clusters_1 <- as.character(dpt_cds_downsampled_cells_ordered[[x[[1]]]]$branch[overlap_cells, 1])
    clusters_2 <- as.character(dpt_cds_downsampled_cells_ordered[[x[[2]]]]$branch[overlap_cells, 1])
    ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)

    return(list(cor = c(cor_res, AT1_cor_res, AT2_cor_res), pearson_cor = c(pearson_cor_res, AT1_pearson_cor_res, AT2_pearson_cor_res), cluster = ClusteringMetrics_res, max_root_ind = max_root_ind))
  }
  else
    return(list(cor = rep(NA, 3), pearson_cor = rep(NA, 3), cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
}, mc.cores = detectCores() - 2)

############################################################################################################################################################
dpt_cmbn_sets <- expand.grid(1:length(cds_downsampled_cells_ordered_0.8), 1:length(cds_downsampled_cells_ordered_0.8))

dpt_cmbn_sets_list <- split(dpt_cmbn_sets[, ], f = 1:nrow(dpt_cmbn_sets[, ]))
dpt_benchmark_res_list <- mclapply(dpt_cmbn_sets_list, function(x) {
  message(x)
  print(x)
  cds_1 <- cds_downsampled_cells_ordered_0.8[[x[[1]]]]
  cds_2 <- cds_downsampled_cells_ordered_0.8[[x[[2]]]]

  if(any(c(is.na(cds_1), is.na(cds_2)))) {
    return(list(cor = rep(NA, 3), pearson_cor = rep(NA, 3), cluster = data.frame(randIndex = c(NA, NA, NA), Type = c("rand index", "variation of information", "adjusted rand index"))))
  }
  else {  # if(any(x == 66)){
    #   return(list(cor = NA, cluster = data.frame(randIndex = rep(NA, 3), Type = c('rand index', 'variation of information', 'adjusted rand index'))))
    # }
    overlap_cells <- intersect(colnames(cds_1), colnames(cds_2))

    t_1 <- dpt_cds_downsampled_cells_ordered_0.8[[x[[1]]]]$pt[overlap_cells] #[overlpa_cells, 'DPT']
    t_2 <- dpt_cds_downsampled_cells_ordered_0.8[[x[[2]]]]$pt[overlap_cells]

    cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
    pearson_cor_res <- cor(t_1, t_2)

    AT1_lineage <- c(dpt_root_state, dpt_AT1_state)
    AT2_lineage <- c(dpt_root_state, dpt_AT2_state)
    AT1_t_1 <- dpt_cds_downsampled_cells_ordered_0.8[[x[[1]]]]$pt[overlap_cells[overlap_cells %in% AT1_state]]
    AT1_t_2 <- dpt_cds_downsampled_cells_ordered_0.8[[x[[2]]]]$pt[overlap_cells[overlap_cells %in% AT1_state]]
    AT2_t_1 <- dpt_cds_downsampled_cells_ordered_0.8[[x[[1]]]]$pt[overlap_cells[overlap_cells %in% AT2_state]]
    AT2_t_2 <- dpt_cds_downsampled_cells_ordered_0.8[[x[[2]]]]$pt[overlap_cells[overlap_cells %in% AT2_state]]

    AT1_cor_res <- cor(AT1_t_1, AT1_t_2, method = 'kendall', use = 'pairwise.complete.obs')
    AT1_pearson_cor_res <- cor(AT1_t_1, AT1_t_2)
    AT2_cor_res <- cor(AT2_t_1, AT2_t_2, method = 'kendall', use = 'pairwise.complete.obs')
    AT2_pearson_cor_res <- cor(AT2_t_1, AT2_t_2)
    #branch assignment:
    # overlpa_cells_update <- overlpa_cells[intersect(which(as.character(dpt_cds_downsampled_cells_ordered_0.8[[x[[1]]]]$pt[overlpa_cells, 'Branch']) %in% c('branch 1', 'branch 2', 'branch 3')),
    #                                                 which(as.character(dpt_cds_downsampled_cells_ordered_0.8[[x[[2]]]]$pt[overlpa_cells, 'Branch']) %in% c('branch 1', 'branch 2', 'branch 3')))] #remove the unassigned 1,2,3, uncertain 1,2,3, cells
    clusters_1 <- as.character(dpt_cds_downsampled_cells_ordered_0.8[[x[[1]]]]$branch[overlap_cells, 1])
    clusters_2 <- as.character(dpt_cds_downsampled_cells_ordered_0.8[[x[[2]]]]$branch[overlap_cells, 1])

    ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)

    return(list(cor = c(cor_res, AT1_cor_res, AT2_cor_res), pearson_cor = c(pearson_cor_res, AT1_pearson_cor_res, AT2_pearson_cor_res), cluster = ClusteringMetrics_res))
  }
}, mc.cores = detectCores() - 2)

############################################################################################################################################################
dpt_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(dpt_benchmark_res_list, function(x) mean(x$cor[2:3]))),
                                  pearson_rho = unlist(lapply(dpt_benchmark_res_list, function(x) mean(x$pearson_cor[2:3]))),
                                  rand_ind = unlist(lapply(dpt_benchmark_res_list, function(x) x$cluster[1, 1])),
                                  var_inf = unlist(lapply(dpt_benchmark_res_list, function(x) x$cluster[2, 1])),
                                  adj_rand = unlist(lapply(dpt_benchmark_res_list, function(x) x$cluster[3, 1]))
)
############################################################################################################################################################
dpt_cell_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(dpt_cell_sampling_benchmark_res_list, function(x) mean(x$cor[2:3]))),
                                       pearson_rho = unlist(lapply(dpt_cell_sampling_benchmark_res_list, function(x) mean(x$pearson_cor[2:3]))),
                                       rand_ind = unlist(lapply(dpt_cell_sampling_benchmark_res_list, function(x) x$cluster[1, 1])),
                                       var_inf = unlist(lapply(dpt_cell_sampling_benchmark_res_list, function(x) x$cluster[2, 1])),
                                       adj_rand = unlist(lapply(dpt_cell_sampling_benchmark_res_list, function(x) x$cluster[3, 1]))
)

progressive_downsampling_rng <- (nrow(dpt_cell_sampling_res_df) - 35):nrow(dpt_cell_sampling_res_df)
valid_dpt_cell_sampling_res_df <- dpt_cell_sampling_res_df[progressive_downsampling_rng, ]

# ############################################################################################################################################################
# #show benchmark results for wishbone (cell robustness downsampling)):
# ############################################################################################################################################################
wishbone_cmbn_sets <- expand.grid(unique(wishbone_res$run), unique(wishbone_res$run))

wishbone_cmbn_sets_list <- split(wishbone_cmbn_sets[, ], f = 1:nrow(wishbone_cmbn_sets[, ]))

wishbone_benchmark_res_list <- mclapply(wishbone_cmbn_sets_list, function(x) {
  message(x)
  cds_1 <- cds_downsampled_cells_ordered_0.8[[x[[1]]]]
  cds_2 <- cds_downsampled_cells_ordered_0.8[[x[[2]]]]
  overlpa_cells <- intersect(colnames(cds_1), colnames(cds_2))

  wishbone_res_1 <- subset(wishbone_res, run == x[[1]])
  wishbone_res_2 <- subset(wishbone_res, run == x[[2]])
  row.names(wishbone_res_1) <- wishbone_res_1$X
  row.names(wishbone_res_2) <- wishbone_res_2$X

  t_1 <- wishbone_res_1[overlpa_cells, 'trajectory']
  t_2 <- wishbone_res_2[overlpa_cells, 'trajectory']

  cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
  pearson_cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
  
  AT1_lineage <- c(root_state, AT1_state)
  AT2_lineage <- c(root_state, AT2_state)
  AT1_t_1 <- wishbone_res_1[overlpa_cells[overlap_cells %in% AT1_state], 'trajectory']
  AT1_t_2 <- wishbone_res_2[overlpa_cells[overlap_cells %in% AT1_state], 'trajectory']
  AT2_t_1 <- wishbone_res_1[overlpa_cells[overlap_cells %in% AT2_state], 'trajectory']
  AT2_t_2 <- wishbone_res_2[overlpa_cells[overlap_cells %in% AT2_state], 'trajectory']

  AT1_cor_res <- cor(AT1_t_1, AT1_t_2, method = 'kendall', use = 'pairwise.complete.obs')
  AT1_pearson_cor_res <- cor(AT1_t_1, AT1_t_2)
  AT2_cor_res <- cor(AT2_t_1, AT2_t_2, method = 'kendall', use = 'pairwise.complete.obs')
  AT2_pearson_cor_res <- cor(AT2_t_1, AT2_t_2)

  if(cor_res < 0 & is.finite(cor_res)){
    start_cell_id <- which(t_1 == min(t_1, na.rm = T)) #finding the starting cell (cell with smallest pseudotime) in the overlapping set
    t_2_update <- abs(t_2 - t_2[start_cell_id])
    cor_res <- cor(t_1, t_2_update, method = 'kendall', use = 'pairwise.complete.obs')
  }

  #branch assignment:
  clusters_1 <- as.character(wishbone_res_1[overlpa_cells, 'branch'])
  clusters_2 <- as.character(wishbone_res_2[overlpa_cells, 'branch'])
  ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)

  return(list(cor = c(cor_res, AT1_cor_res, AT2_cor_res), pearson_cor = c(pearson_cor_res, AT1_pearson_cor_res, AT2_pearson_cor_res), cluster = ClusteringMetrics_res))
}, mc.cores = detectCores() - 2)

wishbone_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(wishbone_benchmark_res_list, function(x) mean(x$cor))),
                                       pearson_rho = unlist(lapply(wishbone_benchmark_res_list, function(x) mean(x$pearson_cor))),
                                       rand_ind = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cluster[1, 1])),
                                       var_inf = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cluster[2, 1])),
                                       adj_rand = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cluster[3, 1]))
)

############################################################################################################################################################
#show benchmark results for wishbone (cell robustness downsampling)):
############################################################################################################################################################
fraction_wishbone_res <- read.table(paste("/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/", wishbone_res_fraction, sep = ''), header = T, row.names=NULL)
# fraction_wishbone_res <- read.table("/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/fraction_all_wishbone_res_df_nw_10.txt", header = T, row.names=NULL)

wishbone_cmbn_sets <- expand.grid(unique(fraction_wishbone_res$run), unique(fraction_wishbone_res$run))

wishbone_cmbn_sets_list <- split(wishbone_cmbn_sets[, ], f = 1:nrow(wishbone_cmbn_sets[, ]))

wishbone_benchmark_res_list <- mclapply(wishbone_cmbn_sets_list, function(x) {
  cds_1 <- cds_downsampled_cells_ordered[[x[[1]]]]
  cds_2 <- cds_downsampled_cells_ordered[[x[[2]]]]
  overlpa_cells <- intersect(colnames(cds_1), colnames(cds_2))

  fraction_wishbone_res_1 <- subset(fraction_wishbone_res, run == x[[1]])
  fraction_wishbone_res_2 <- subset(fraction_wishbone_res, run == x[[2]])
  row.names(fraction_wishbone_res_1) <- fraction_wishbone_res_1$row.names
  row.names(fraction_wishbone_res_2) <- fraction_wishbone_res_2$row.names

  t_1 <- fraction_wishbone_res_1[overlpa_cells, 'trajectory']
  t_2 <- fraction_wishbone_res_2[overlpa_cells, 'trajectory']

  cor_res <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
  pearson_cor_res <- cor(t_1, t_2)
  
  AT1_lineage <- c(root_state, AT1_state)
  AT2_lineage <- c(root_state, AT2_state)
  AT1_t_1 <- fraction_wishbone_res_1[overlpa_cells[overlap_cells %in% AT1_state], 'trajectory']
  AT1_t_2 <- fraction_wishbone_res_2[overlpa_cells[overlap_cells %in% AT1_state], 'trajectory']
  AT2_t_1 <- fraction_wishbone_res_1[overlpa_cells[overlap_cells %in% AT2_state], 'trajectory']
  AT2_t_2 <- fraction_wishbone_res_2[overlpa_cells[overlap_cells %in% AT2_state], 'trajectory']

  AT1_cor_res <- cor(AT1_t_1, AT1_t_2, method = 'kendall', use = 'pairwise.complete.obs')
  AT1_pearson_cor_res <- cor(AT1_t_1, AT1_t_2)
  AT2_cor_res <- cor(AT2_t_1, AT2_t_2, method = 'kendall', use = 'pairwise.complete.obs')
  AT2_pearson_cor_res <- cor(AT2_t_1, AT2_t_2)

  #branch assignment:
  clusters_1 <- as.character(fraction_wishbone_res_1[overlpa_cells, 'branch'])
  clusters_2 <- as.character(fraction_wishbone_res_2[overlpa_cells, 'branch'])
  ClusteringMetrics_res <- calClusteringMetrics(clusters_1, clusters_2)

  return(list(cor = c(cor_res, AT1_cor_res, AT2_cor_res), pearson_cor = c(pearson_cor_res, AT1_pearson_cor_res, AT2_pearson_cor_res), cluster = ClusteringMetrics_res))
}, mc.cores = detectCores() - 2)

wishbone_cell_sampling_res_df <- data.frame(kendall.tau = unlist(lapply(wishbone_benchmark_res_list, function(x) mean(x$cor))),
                                            pearson_rho = unlist(lapply(wishbone_benchmark_res_list, function(x) mean(x$pearson_cor))),
                                            rand_ind = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cluster[1, 1])),
                                            var_inf = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cluster[2, 1])),
                                            adj_rand = unlist(lapply(wishbone_benchmark_res_list, function(x) x$cluster[3, 1]))
)

wishbone_res_dim <- dim(wishbone_cell_sampling_res_df)
valid_wishbone_cell_sampling_res_df <- wishbone_cell_sampling_res_df[(wishbone_res_dim[1] - length(unique(fraction_wishbone_res$run)) + 1):wishbone_res_dim[1], ]

############################################################################################################################################################
#generate supplementary figures for benchmarking results
############################################################################################################################################################

# ###cell downsampling: color by both of the dpt and ddrtree result together: 
all_valid_cell_sampling_res_df <- Reduce(rbind , list(valid_dpt_cell_sampling_res_df, valid_cell_sampling_res_df, valid_wishbone_cell_sampling_res_df, ICA_valid_cell_sampling_res_df))

all_valid_cell_sampling_res_df$proportion <- c(rep(downsampled_proportions, 2), names(cds_downsampled_cells_ordered)[unique(fraction_wishbone_res$run)], downsampled_proportions)
all_valid_cell_sampling_res_df$Type <- c(rep('dpt', 36), rep('monocle2', 36), rep('wishbone', length(unique(fraction_wishbone_res$run))), rep('monocle1', 36))

all_valid_cell_sampling_res_df$Type <- factor(all_valid_cell_sampling_res_df$Type, levels = c('monocle2', 'monocle1', "dpt", "wishbone")) 
all_valid_cell_sampling_res_df$se <- 0.1 

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_cell_downsampling.pdf', sep = ''), width = 2.5, height = 2) 
qplot(proportion, abs(all_valid_cell_sampling_res_df$pearson_rho), data = all_valid_cell_sampling_res_df, color = Type, size = 1, geom = 'boxplot') +
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + nm_theme() + scale_size(range = c(0.1, 1)) + monocle_theme_opts() + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'rand_index_comparison_cell_downsampling.pdf', sep = ''), width = 2.5, height = 2) 
qplot(proportion, abs(all_valid_cell_sampling_res_df$adj_rand), data = all_valid_cell_sampling_res_df, color = Type, size = 1, geom = 'boxplot') + 
  xlab('Proportion of original cells') + ylab("Adjusted rand index") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

process_cell_sampling_res_df <- ddply(all_valid_cell_sampling_res_df, .(Type, proportion), summarize, 
                                      mean_kendall.tau = mean(abs(kendall.tau), na.rm = T), 
                                      sd_kendall.tau = sd(abs(kendall.tau), na.rm = T), 
                                      
                                      mean_pearson_rho = mean(abs(pearson_rho), na.rm = T), 
                                      sd_pearson_rho = sd(abs(pearson_rho), na.rm = T), 
                                      
                                      mean_adj_rand = mean(abs(adj_rand), na.rm = T),
                                      sd_adj_rand = sd(abs(adj_rand), na.rm = T),
                                      
                                      se = mean(se))
limits <- aes(ymax = mean_adj_rand + sd_adj_rand, ymin=mean_adj_rand - sd_adj_rand)

pdf(paste(SI_fig_dir, benchmark_type, 'rand_index_comparison_cell_downsampling2.pdf', sep = ''), width = 3, height = 2) 
ggplot(aes(proportion, mean_adj_rand), data = process_cell_sampling_res_df) + 
  geom_point(aes(color = Type), position=position_dodge(width=0.1), size = 0.5) + facet_wrap(~Type) + 
  geom_errorbar(aes(color = Type, ymax = mean_adj_rand + sd_adj_rand, ymin=mean_adj_rand - sd_adj_rand), position=position_dodge(width=0.9), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("Adjusted rand index") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_cell_downsampling2.pdf', sep = ''), width = 2.5, height = 2) 
ggplot(aes(proportion, mean_pearson_rho), data = process_cell_sampling_res_df) + 
  geom_point(aes(color = Type), position=position_dodge(width=0.9), size = 0.5) + facet_wrap(~Type) + 
  geom_errorbar(aes(color = Type, ymax = mean_pearson_rho + sd_pearson_rho, ymin=mean_pearson_rho - sd_pearson_rho), position=position_dodge(width=0.1), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + scale_size(range = c(0.1, 1)) + nm_theme() + monocle_theme_opts() + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_cell_downsampling2_helper.pdf', sep = '')) 
ggplot(aes(proportion, mean_pearson_rho), data = process_cell_sampling_res_df) + 
  geom_point(aes(color = Type), position=position_dodge(width=0.1), size = 0.5) + 
  geom_errorbar(aes(color = Type, ymax = mean_pearson_rho + sd_pearson_rho, ymin=mean_pearson_rho - sd_pearson_rho), position=position_dodge(width=0.1), size = 0.5) + 
  xlab('Proportion of original cells') + ylab("Pearson's Rho") + scale_size(range = c(0.1, 1)) + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

############################################################################################################################################################
#create mean and variance plot: 
mean_df <- ddply(all_valid_cell_sampling_res_df, .(Type, proportion), function(x) data.frame(mean_kt = mean(x$pearson_rho), 
                                                                                             mean_rand_ind = mean(x$rand_ind), 
                                                                                             mean_var_inf = mean(x$var_inf), 
                                                                                             mean_adj_rand = mean(x$adj_rand) ))

pdf(paste(SI_fig_dir, benchmark_type, 'mean_pearson_rho_comparison_cell_downsampling.pdf', sep = ''), width = 2.5, height = 2) 
qplot(proportion, abs(mean_kt), data = mean_df, group = Type, color = Type, size = 1, geom = 'line') +
  xlab('Proportion of original cells') + ylab("Mean Pearson's Rho") + nm_theme() + scale_size(range = c(0.1, 1)) + monocle_theme_opts() + nm_theme() + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'mean_adj_rand_index_comparison_cell_downsampling.pdf', sep = ''), width = 2.5, height = 2) 
qplot(proportion, abs(mean_df$mean_adj_rand), data = mean_df, color = Type, group = Type, size = 1, geom = 'line') +
  xlab('Proportion of original cells') + ylab("Mean adjusted rand index") + nm_theme() + scale_size(range = c(0.1, 1)) + monocle_theme_opts() + nm_theme() + ylim(0, 1)
dev.off()

variance_df <- ddply(all_valid_cell_sampling_res_df, .(Type, proportion), function(x) data.frame(var_kt = var(x$pearson_rho), 
                                                                                             var_rand_ind = var(x$rand_ind), 
                                                                                             var_var_inf = var(x$var_inf), 
                                                                                             var_adj_rand = var(x$adj_rand) ))

pdf(paste(SI_fig_dir, benchmark_type, 'var_pearson_rho_comparison_cell_downsampling.pdf', sep = ''), width = 2.5, height = 2) 
qplot(proportion, var_kt, data = variance_df, group = Type, color = Type, size = 1, geom = 'line') +
  xlab('Proportion of original cells') + ylab("Variance of Pearson's Rho") + nm_theme() + scale_size(range = c(0.1, 1)) + monocle_theme_opts() + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'var_adj_rand_index_comparison_cell_downsampling.pdf', sep = ''), width = 2.5, height = 2) 
qplot(proportion, var_adj_rand, data = variance_df, group = Type, color = Type, size = 1, geom = 'line') +
  xlab('Proportion of original cells') + ylab("Varianceadjusted rand index") + nm_theme() + scale_size(range = c(0.1, 1)) + monocle_theme_opts() + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'var_adj_rand_index_comparison_cell_downsampling_helper.pdf', sep = '')) 
qplot(proportion, var_adj_rand, data = variance_df, group = Type, color = Type, size = 1, geom = 'line') + 
  xlab('Proportion of original cells') + ylab("Varianceadjusted rand index") + scale_size(range = c(0.1, 1)) + ylim(0, 1)
dev.off()

############################################################################################################################################################
#use ICA downsample data: 
# load('./RData/ICA_downsampling.RData')

############################################################################################################################################################
###robustness downsampling: color by both of the dpt and ddrtree result together: 
all_sampling_res_df <- Reduce(rbind , list(dpt_sampling_res_df, sampling_res_df,  wishbone_sampling_res_df, ICA_sampling_res_df)) # ICA_sampling_res_df,

# #use result from DDRTree
#all_sampling_res_df <- Reduce(rbind , list(dpt_sampling_res_df, destiny_dm_DDRTree_sampling_res_df,  wishbone_sampling_res_df, ICA_sampling_res_df)) #DDRTree_sampling_res_df ICA_sampling_res_df,
# all_sampling_res_df$Type <- c(rep('dpt', 9604), rep('Monocle2', 10000), rep('Wishbone', 6724)) #,  rep('Monocle1', 10000)
# 
# #use result from DDRTree and dpt knockout
#all_sampling_res_df <- Reduce(rbind , list(dpt_sampling_res_df, DDRTree_sampling_res_df,  wishbone_sampling_res_df, ICA_sampling_res_df)) # ICA_sampling_res_df,
# all_sampling_res_df$Type <- c(rep('dpt', 10000), rep('Monocle2', 10000), rep('Wishbone', 6724)) #,  rep('Monocle1', 10000)

all_sampling_res_df$Type <- c(rep('dpt', nrow(dpt_sampling_res_df)), rep('monocle2',  nrow(sampling_res_df)), rep('wishbone', nrow(wishbone_sampling_res_df)), rep('monocle1', nrow(ICA_sampling_res_df)))#,  rep('Monocle1', 10000)
all_sampling_res_df$Type <- factor(all_sampling_res_df$Type, levels = c('monocle2', 'monocle1', "dpt", "wishbone")) #dpt (non-uniform branch)

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_robustness.pdf', sep = ''), width = 1, height = 1.5) 
qplot(Type, pearson_rho, data = all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Pearson's Rho") + monocle_theme_opts()  + xlab('') + 
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_robustness_helper.pdf', sep = '')) 
qplot(Type, adj_rand, data = all_sampling_res_df, color = Type, geom = 'boxplot') + xlab("Pearson's Rho") + ylab('cells') + monocle_theme_opts() + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'rand_index_comparison_robustness.pdf', sep = ''), width = 1, height = 1.5) 
qplot(Type, adj_rand, data = all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Adjusted Rand index") + nm_theme() + xlab('') + monocle_theme_opts() + nm_theme() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1) + scale_color_manual(values = software_custom_color_scale)
dev.off()
############################################################################################################################################################
#change to plot the trend of mean as well as the trend of variance for each method: 
#mean vs variance: 
pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_robustness_mean.pdf', sep = ''), width = 1.5, height = 1.5) 
qplot(Type, pearson_rho, data = all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Distribution of Pearson's Rho") + nm_theme()  + xlab('') + monocle_theme_opts() + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'pearson_rho_comparison_robustness_helper.pdf', sep = ''), width = 1.5, height = 1.5) 
qplot(Type, adj_rand, data = all_sampling_res_df, color = Type, geom = 'boxplot') + xlab("Pearson's Rho") + ylab('cells') + monocle_theme_opts() + nm_theme() + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

pdf(paste(SI_fig_dir, benchmark_type, 'rand_index_comparison_robustness_variance.pdf', sep = ''), width = 1.5, height = 1.5) 
qplot(Type, adj_rand, data = all_sampling_res_df, color = Type, geom = 'boxplot') + ylab("Adjusted Rand index") + nm_theme() + xlab('')  + monocle_theme_opts() + ylim(0, 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylim(0, 1)
dev.off()

############################################################################################################################################################
save.image(paste('./RData/', benchmark_type, '_calculate_monocle12_slicer_wishbone_dpt_downsampling.RData', sep = ''))
############################################################################################################################################################
