#plotting functions: 
plot_pseudotime_realtime <- function(Pseudotime, Realtime, State = NULL
                                     ){
  cor.coeff <- cor(as.vector(Pseudotime), Realtime, method = "spearman")
  if(cor.coeff < 0)
    Realtime <- rev(Realtime)
  
  maturation_df <- data.frame(cell = rep(paste('cell', names(Pseudotime), sep = '_'), 
                                         2), original_pseudotime = 100 * c(Pseudotime/max(Pseudotime), 
                                                                        Realtime/max(Realtime)), 
                              Type = rep(c("Pseudotime", "Realtime"), each = length(Pseudotime)), rownames = names(Pseudotime))
  
  if(is.null(State))
    maturation_df$State <- 1 
  else 
    maturation_df$State <- State
  
  scaling_df <- ddply(maturation_df, .(State, Type), function(x) {
    rng <- range(as.numeric(x$original_pseudotime))
    scaling_factor <- 100 / diff(rng)
    
    data.frame(min = rng[1], max = rng[2], scaling_factor = scaling_factor)
  })
  
  res <- apply(maturation_df, 1, function(x, scaling_df) {
    # print(x)
    subset_scaling_df <- subset(scaling_df, as.character(Type) == as.character(x[3]) & as.character(State) == as.character(x[5]))
    # cbind(x , scaled_psedutoime = as.vector((as.numeric(x[2]) - subset_scaling_df[, 3]) * subset_scaling_df[, 5]))
    scaled_psedutoime = as.vector((as.numeric(x[2]) - subset_scaling_df[, 3]) * subset_scaling_df[, 5])
  }, scaling_df)
  
  maturation_df <- cbind(maturation_df, maturation_level = res)

  p1 <- ggplot(aes(x = maturation_level, y = Type, group = cell), 
              data = maturation_df) + geom_point(size = 1) + geom_line(aes(color = State), alpha = .3) + 
    xlab("Maturation level") + ylab("Type") + facet_wrap(~State) + 
    annotate("text", x = 80, y = 2.2, label = paste("Spearman correlation:", 
                                                    round(abs(cor.coeff), 2))) + monocle_theme_opts()
  p2 <- ggplot(aes(x = original_pseudotime, y = Type, group = cell), 
               data = maturation_df) + geom_point(size = 1) + geom_line(aes(color = State), alpha = .3) + 
    xlab("Maturation level") + ylab("Type") + facet_wrap(~State) + 
    annotate("text", x = 80, y = 2.2, label = paste("Spearman correlation:", 
                                                    round(abs(cor.coeff), 2))) + monocle_theme_opts()
  
  return(list(p1 = p1, p2 = p2))
}

plot_cross_map <- function(lib_xmap_target_means, target_xmap_lib_means, lib_name, target_name){
  legend_names <- c(paste(lib_name, 'xmap', target_name), paste(target_name, 'xmap', lib_name))
  
  xmap_all <- rbind(a_xmap_t_means, t_xmap_a_means)
  xmap_all$type <- c(rep('a_xmap_t_means', nrow(a_xmap_t_means)), rep('t_xmap_a_means', nrow(t_xmap_a_means)))
  y_max <- max(xmap_all$rho, na.rm = T) + 0.1
  
  lib_rng <- range(xmap_all$lib_size)
  p1 <- ggplot(aes(lib_size, pmax(0, rho)), data = xmap_all) + geom_line(aes(color = type)) + xlim(lib_rng) + 
    xlab("Library Size") + ylab("Cross Map Skill (rho)") + scale_color_discrete(labels=legend_names) + 
    scale_x_discrete(breaks = unique(xmap_all$lib_size)) + monocle_theme_opts()
  
  return(p1)
}

plot_DDTree_complex_tree <- function(DDRTree_res) {
  tree_df <- DDRTree_res$stree[1:ncol(DDRTree_res$Y), 1:ncol(DDRTree_res$Y)]
  stree <- graph.adjacency(tree_df, mode=c("undirected"), weighted=T)
  ordering_res <- custom_ordering(DDRTree_res)
  
  edge_list <- as.data.frame(get.edgelist(stree))
  colnames(edge_list) <- c("source", "target")
  
  root_cell <- subset(ordering_res$cc_ordering, pseudo_time == 0)[, "sample_name"]
  layout_coord <- layout.reingold.tilford(stree, root= as.character(root_cell))
  colnames(layout_coord) <- c('x', 'y')
  row.names(layout_coord) <- as.character(1:nrow(layout_coord))
  
  edge_df <- merge(layout_coord, edge_list, by.x="row.names", by.y="source", all=TRUE)
  
  edge_df <- plyr::rename(edge_df, c("x"="source_x", "y"="source_y"))
  edge_df <- merge(edge_df, layout_coord, by.x="target", by.y="row.names", all=TRUE)
  edge_df <- plyr::rename(edge_df, c("x"="target_x", "y"="target_y"))
  
  g <- ggplot(data=edge_df) + geom_segment(aes_string(x="source_x", y="source_y", xend="target_x", yend="target_y"), 
        size=.3, linetype="solid", na.rm=TRUE, data=edge_df, color = ordering_res$cc_ordering$cell_state) + 
        label() + #color by the cell types 
    
  
  return(list(g = g, edge_df = edge_df, ordering_res = ordering_res))
}

plot_DDTree_complex_tree_3d <- function(DDRTree_res) {
  tree_df <- DDRTree_res$stree[1:ncol(DDRTree_res$Y), 1:ncol(DDRTree_res$Y)]
  stree <- graph.adjacency(tree_df, mode=c("undirected"), weighted=T)
  ordering_res <- custom_ordering(DDRTree_res)
  
  edge_list <- as.data.frame(get.edgelist(stree))
  colnames(edge_list) <- c("source", "target")
  
  root_cell <- subset(ordering_res$cc_ordering, pseudo_time == 0)[, "sample_name"]
  layout_coord <- layout.reingold.tilford(stree, root= as.character(root_cell))
  colnames(layout_coord) <- c('x', 'y')
  row.names(layout_coord) <- as.character(1:nrow(layout_coord))
  
  edge_df <- merge(layout_coord, edge_list, by.x="row.names", by.y="source", all=TRUE)
  
  edge_df <- plyr::rename(edge_df, c("x"="source_x", "y"="source_y"))
  edge_df <- merge(edge_df, layout_coord, by.x="target", by.y="row.names", all=TRUE)
  edge_df <- plyr::rename(edge_df, c("x"="target_x", "y"="target_y"))
  
  plot_ly(all_cell_simulation_df, x = Mash1, y = Scl, z = Olig2, type = "scatter3d", mode = "markers", color = 'red')
  
  g <- plot_ly(data=edge_df) + geom_segment(aes_string(x="source_x", y="source_y", xend="target_x", yend="target_y"), 
        size=.3, linetype="solid", na.rm=TRUE, data=edge_df, color = ordering_res$cc_ordering$cell_state) + 
        label() + #color by the cell types 
    
  
  return(list(g = g, edge_df = edge_df, ordering_res = ordering_res))
}

plot_order_cell_tree <- function(order_genes, HSMM_myo, initial_method = PCA, color_by="Time") {
  HSMM_myo <- setOrderingFilter(HSMM_myo, order_genes)
  plot_ordering_genes(HSMM_myo)
  
  HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, norm_method = 'log', initial_method = initial_method) #, initial_method = DM , initial_method = DM
  HSMM_myo <- orderCells(HSMM_myo, num_paths=1, root_state = NULL)
  plot_spanning_tree(HSMM_myo, color_by=color_by, cell_size=2)
}

#plotting function to show the smoothness and tightness of the fitted trajectories: 
plot_tightness_smoothness <- function(tightness_df, smoothness_df) {
  tightness_smoothness_df <- rbind(tightness_df, smoothness_df)

  all_finite_genes <- apply(apply(tightness_smoothness_df[, 1:(ncol(tightness_smoothness_df) - 2)], 2, function(x) is.finite(x)), 1, all)
  all_finite_genes <- names(all_finite_genes)[all_finite_genes]
  
  tightness_smoothness_df$gene_ids <- row.names(tightness_smoothness_df)
  mlt_ts_df <- melt(tightness_smoothness_df[all_finite_genes, ], id.var = c('type', 'gene_ids'))
  
  ddply(mlt_ts_df, .(type, variable), function(x) range(x$value, na.rm = T))
  
  #density plot: 
  p1 <- ggplot(aes(x = value + 1), data = mlt_ts_df) + geom_density(aes(color = variable)) + facet_wrap(~type, scale = 'free') + scale_x_log10()

  #scatterplot: 
  p2 <- ggplot(aes(x = gene_ids, y = value + 1), data = mlt_ts_df) + geom_point(aes(color = variable)) + scale_y_log10() + facet_wrap(~type, scale = 'free') 
  
  #barplot: 
  p3 <- ggplot(aes(x = variable, y = value + 1), data = mlt_ts_df) + geom_boxplot(aes(color = variable)) + scale_y_log10() + facet_wrap(~type, scale = 'free') 
  
  return(list(mlt_ts_df = mlt_ts_df, p1 = p1, p2 = p2, p3 = p3))
}

#plot the kinetic curve defined by reference ordering: 
plot_gene_clusters <- function(cds, clusters){
  ordered_cds <- cds[names(clusters), order(pData(cds)$Pseudotime)]
  autocorrelation_list <- lapply(unique(clusters), function(x){
    avg_cluster_genes <- esApply(ordered_cds[clusters == x, ], 2, mean)
  })
  res <- do.call(rbind.data.frame, autocorrelation_list)
  colnames(res) <- names(autocorrelation_list[[1]]); row.names(res) <- paste('Cluster', unique(clusters), sep = '_')
  mlt_res <- melt(res)
  mlt_res$cluster <- rep(paste('Cluster', unique(clusters), sep = '_'), ncol(ordered_cds))
  mlt_res$Order <- rep(1:ncol(ordered_cds), each = 4)
  
  qplot(Order, value, color = cluster, data = mlt_res)
  qplot(Order, value, color = cluster, data = mlt_res, geom = 'line')
}

# Modified function: Plot heatmap of 3 branches with the same coloring. Each CDS subset has to have the same set of genes.
plot_multiple_branches_heatmap <- function(cds, 
                                           branches, 
                                           branches_name = NULL, 
                                           cluster_rows = TRUE,
                                           hclust_method = "ward.D2", 
                                           num_clusters = 6,
                                           
                                           hmcols = NULL, 
                                           
                                           add_annotation_row = NULL,
                                           add_annotation_col = NULL,
                                           show_rownames = FALSE, 
                                           use_gene_short_name = TRUE,
                                           
                                           norm_method = c("vstExprs", "log"), 
                                           scale_max=3, 
                                           scale_min=-3, 
                                           
                                           trend_formula = '~sm.ns(Pseudotime, df=3)',
                                           
                                           return_heatmap=FALSE,
                                           cores=1){
  if(!(all(branches %in% pData(cds)$State)) & length(branches) == 1){
    stop('This function only allows to make multiple branch plots where branches is included in the pData')
  }
  
  branch_label <- branches
  if(!is.null(branches_name)){
    if(length(branches) != length(branches_name)){
      stop('branches_name should have the same length as branches')
    }
    branch_label <- branches_name
  }
  
  #test whether or not the states passed to branches are true branches (not truncks) or there are terminal cells 
  g <- cds@minSpanningTree
  m <- NULL
  # branche_cell_num <- c()
  for(branch_in in branches) {
    branches_cells <- row.names(subset(pData(cds), State == branches))
    root_state <- subset(pData(cds), Pseudotime == 0)[, 'State']
    root_state_cells <- row.names(subset(pData(cds), State == root_state))
    
    root_cell <- root_state_cells[which(degree(g, v = root_state_cells) == 1)]
    tip_cell <- branches_cells[which(degree(g, v = branches_cells) == 1)]
    
    traverse_res <- traverseTree(g, root_cell, tip_cell)
    path_cells <- names(traverse_res$shortest_path[[1]])
    cds_subset <- cds[, path_cells]
    
    newdata <- data.frame(Pseudotime = seq(0, max(pData(cds_subset)$Pseudotime),length.out = 100)) 
  
    tmp <- genSmoothCurves(cds_subset, cores=cores, trend_formula = trend_formula,  
                          relative_expr = T, new_data = newdata)
    if(is.null(m))
      m <- tmp
    else
      m <- cbind(m, tmp)
  }
  
  #remove genes with no expression in any condition
  m=m[!apply(m,1,sum)==0,]
  
  norm_method <- match.arg(norm_method)
  
  # FIXME: this needs to check that vst values can even be computed. (They can only be if we're using NB as the expressionFamily)
  if(norm_method == 'vstExprs' && is.null(cds@dispFitInfo[["blind"]]$disp_func) == FALSE) {
    m = vstExprs(cds, expr_matrix=m)
  }     
  else if(norm_method == 'log') {
    m = log10(m+pseudocount)
  }
  
  # Row-center the data.
  m=m[!apply(m,1,sd)==0,]
  m=Matrix::t(scale(Matrix::t(m),center=TRUE))
  m=m[is.na(row.names(m)) == FALSE,]
  m[is.nan(m)] = 0
  m[m>scale_max] = scale_max
  m[m<scale_min] = scale_min
  
  heatmap_matrix <- m
  
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  
  if(is.null(hmcols)) {
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- blue2green2red(length(bks) - 1)
  }
  else {
    bks <- seq(-3.1,3.1, length.out = length(hmcols))
  } 
  
  ph <- pheatmap(heatmap_matrix, 
                 useRaster = T,
                 cluster_cols=FALSE, 
                 cluster_rows=cluster_rows, 
                 show_rownames=F, 
                 show_colnames=F, 
                 clustering_distance_rows=row_dist,
                 clustering_method = hclust_method,
                 cutree_rows=num_clusters,
                 silent=TRUE,
                 filename=NA,
                 breaks=bks,
                 color=hmcols)
  
  annotation_col <- data.frame(Branch=factor(rep(rep(branch_label, each = 100))))
  annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
  gaps_col <- c(1:length(branches) - 1) * 100
  
  if(!is.null(add_annotation_row)) {
    old_colnames_length <- ncol(annotation_row)
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])  
    colnames(annotation_row)[(old_colnames_length+1):ncol(annotation_row)] <- colnames(add_annotation_row)
    # annotation_row$bif_time <- add_annotation_row[as.character(fData(absolute_cds[row.names(annotation_row), ])$gene_short_name), 1]
  }
  
  
  if (use_gene_short_name == TRUE) {
    if (is.null(fData(cds)$gene_short_name) == FALSE) {
      feature_label <- as.character(fData(cds)[row.names(heatmap_matrix), 'gene_short_name'])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
      
      row_ann_labels <- as.character(fData(cds)[row.names(annotation_row), 'gene_short_name'])
      row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
    }
    else {
      feature_label <- row.names(heatmap_matrix)
      row_ann_labels <- row.names(annotation_row)
    }
  }
  else {
    feature_label <- row.names(heatmap_matrix)
    row_ann_labels <- row.names(annotation_row)
  }
  
  row.names(heatmap_matrix) <- feature_label
  row.names(annotation_row) <- row_ann_labels
  
  
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  
  ph_res <- pheatmap(heatmap_matrix[, ], #ph$tree_row$order
                     useRaster = T,
                     cluster_cols = FALSE, 
                     cluster_rows = cluster_rows, 
                     show_rownames=show_rownames, 
                     show_colnames=F, 
                     #scale="row",
                     clustering_distance_rows=row_dist, #row_dist
                     clustering_method = hclust_method, #ward.D2
                     cutree_rows=num_clusters,
                     # cutree_cols = 2,
                     annotation_row=annotation_row,
                     annotation_col=annotation_col,
                     gaps_col = c(100, 200),
                     treeheight_row = 20, 
                     breaks=bks,
                     fontsize = 12,
                     color=hmcols, 
                     silent=TRUE,
                     filename=NA
  )
  
  grid::grid.rect(gp=grid::gpar("fill", col=NA))
  grid::grid.draw(ph_res$gtable)
  if (return_heatmap){
    return(ph_res)
  }
}

scale_pseudotime_for_heatmap <- function(cds, verbose = F) {
  pd <- pData(cds)
  pd$Cell_name <- row.names(pd)
  range_df <- plyr::ddply(pd, .(State), function(x) {
    min_max <- range(x$Pseudotime)
    min_cell <- subset(x, Pseudotime %in% min_max[1])
    max_cell <- subset(x, Pseudotime %in% min_max[2])
    min_cell$fate_type <- 'Start'
    max_cell$fate_type <- 'End'
    rbind(min_cell, max_cell)
  }) #pseudotime range for each state
  
  #1. construct a tree of the selected cells
  #create a cell state tree: 
  cds@auxOrderingData$DDRTree$Y <- cds@reducedDimK
  res <- custom_ordering(cds@auxOrderingData$DDRTree)

  pd <- res$cc_ordering
  
  adj_list <- data.frame(Source = pd$parent[-1] , Target = pd$sample_name[-1], stringsAsFactors = F)
  colnames(adj_list) <- c('Source', 'Target')
  adj_list <- as.data.frame(apply(adj_list, 2, function(x) mapvalues(x, as.character(pd$sample_name), as.numeric(pd$cell_state), warn_missing = F)), stringsAsFactors = F)
  adj_list$Source <- as.numeric(adj_list$Source)
  adj_list$Target <- as.numeric(adj_list$Target)
  
  uniq_cell_list <- as.character(sort(unique(c(adj_list$Source, adj_list$Target))))
  adj_mat <- matrix(rep(0, length(uniq_cell_list)^2), nrow = length(uniq_cell_list), ncol = length(uniq_cell_list), dimnames = list(uniq_cell_list, uniq_cell_list))
  adj_mat[as.matrix(adj_list)] <- 1
  net <- graph.adjacency(as.matrix(adj_mat), mode = 'directed', weighted = NULL, diag = FALSE)
  v_size <- 20 * table(pd$cell_state) / (max(table(pd$cell_state)) - min(table(pd$cell_state)))
  
  plot(net, layout = layout.fruchterman.reingold,
         vertex.size = 20 * table(pd$cell_state) / max(table(pd$cell_state)),
         vertex.color="red",
         vertex.frame.color= "white",
         vertex.label.color = "white",
         vertex.label.cex = 1,
         vertex.label.family = "sans",
         edge.width=2,
         edge.color="black")
  
  #dfs_net <- graph.dfs(net, root = 1, order.out = T) #DFS search for the cell fate tree
  #net_order_out <- as.vector(dfs_net$order.out)
  
  net_leaves <- which(degree(net, v = V(net), mode = "out")==0, useNames = T)
  
  pd$scale_pseudotime <- NA
  
  #2. scale the psedutime:
  for(i in net_leaves) {
    path_vertex <- as.vector(get.all.shortest.paths(net, from = 1, to = i, mode="out")$res[[1]])
    pd_subset <- subset(pd, cell_state %in% path_vertex & is.na(scale_pseudotime))
    
    #scale the pseudotime between the parent cell of the first cell to the last cell on the remaing branch
    # print(path_vertex)
    min_cell_name <- row.names(subset(pd_subset, pseudo_time == min(pseudo_time)))
    # print(min_cell_name)
    
    if(!is.na(pd[min_cell_name, 'parent'])) {
      parent_min_cell <- as.character(pd[min_cell_name, 'parent'])
      subset_min_pseudo <- pd[parent_min_cell, 'pseudo_time']
      scale_pseudotime_ini <- pd[parent_min_cell, 'scale_pseudotime']
      scaling_factor <- (100 - pd[parent_min_cell, 'scale_pseudotime']) / c(max(pd_subset$pseudo_time) - subset_min_pseudo)
    }
    else {
      subset_min_pseudo <- min(pd[, 'pseudo_time'])
      scale_pseudotime_ini <- 0
      scaling_factor <- 100 / c(max(pd_subset$pseudo_time) - min(pd_subset$pseudo_time))
    }
    
    pseudotime_scaled <- (pd_subset$pseudo_time - subset_min_pseudo) * scaling_factor + scale_pseudotime_ini
    
    if(verbose)
      message(i, '\t', range(pseudotime_scaled)[1],'\t', range(pseudotime_scaled)[2])
    
    pd[row.names(pd_subset), 'scale_pseudotime'] <- pseudotime_scaled
    
    if(verbose)
      message(row.names(pd_subset))
  }	

  pData(cds) <- cbind(pData(cds), pd)
  
  return(cds)
}

#plot cellTypeHiearchy: 
plot_cellTypeHiearchy <- function(cth) {
  plot(cth@classificationTree)
}


plot_multiple_branches_pseudotime <- function(cds, 
                                           branches, 
                                           branches_name = NULL, 
                                           cluster_rows = TRUE,
                                           hclust_method = "ward.D2", 
                                           num_clusters = 6,
                                           
                                           hmcols = NULL, 
                                           
                                           add_annotation_row = NULL,
                                           add_annotation_col = NULL,
                                           show_rownames = FALSE, 
                                           use_gene_short_name = TRUE,
                                           
                                           norm_method = c("vstExprs", "log"), 
                                           scale_max=3, 
                                           scale_min=-3, 
                                           
                                           trend_formula = '~sm.ns(Pseudotime, df=3)',
                                           
                                           return_heatmap=FALSE,
                                           cores=1){
  if(!(all(branches %in% pData(cds)$State)) & length(branches) == 1){
    stop('This function only allows to make multiple branch plots where branches is included in the pData')
  }
  
  branch_label <- branches
  if(!is.null(branches_name)){
    if(length(branches) != length(branches_name)){
      stop('branches_name should have the same length as branches')
    }
    branch_label <- branches_name
  }
  
  #test whether or not the states passed to branches are true branches (not truncks) or there are terminal cells 
  g <- cds@minSpanningTree
  m <- NULL
  # branche_cell_num <- c()
  for(branch_in in branches) {
    branches_cells <- row.names(subset(pData(cds), State == branch_in))
    root_state <- subset(pData(cds), Pseudotime == 0)[, 'State']
    root_state_cells <- row.names(subset(pData(cds), State == root_state))
    
    if(cds@dim_reduce_type != 'ICA') {
      root_state_cells <- unique(paste('Y_', cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[root_state_cells, ], sep = ''))
      branches_cells <- unique(paste('Y_', cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[branches_cells, ], sep = ''))
    }
    root_cell <- root_state_cells[which(degree(g, v = root_state_cells) == 1)]
    tip_cell <- branches_cells[which(degree(g, v = branches_cells) == 1)]
    
    traverse_res <- traverseTree(g, root_cell, tip_cell)
    path_cells <- names(traverse_res$shortest_path[[1]])
    
    if(cds@dim_reduce_type != 'ICA') {
      pc_ind <- cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex 
      path_cells <- row.names(pc_ind)[paste('Y_', pc_ind[, 1], sep = '') %in% path_cells]
    }
    
    cds_subset <- cds[, path_cells]
    
    newdata <- data.frame(Pseudotime = seq(0, max(pData(cds_subset)$Pseudotime),length.out = 100)) 
    
    tmp <- genSmoothCurves(cds_subset, cores=cores, trend_formula = trend_formula,  
                           relative_expr = T, new_data = newdata)
    if(is.null(m))
      m <- tmp
    else
      m <- cbind(m, tmp)
  }
  
  #remove genes with no expression in any condition
  m=m[!apply(m,1,sum)==0,]
  
  norm_method <- match.arg(norm_method)
  
  # FIXME: this needs to check that vst values can even be computed. (They can only be if we're using NB as the expressionFamily)
  if(norm_method == 'vstExprs' && is.null(cds@dispFitInfo[["blind"]]$disp_func) == FALSE) {
    m = vstExprs(cds, expr_matrix=m)
  }     
  else if(norm_method == 'log') {
    m = log10(m+pseudocount)
  }
  
  
}