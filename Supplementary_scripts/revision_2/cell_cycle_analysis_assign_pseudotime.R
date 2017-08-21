#assign pseudotime to cells: 
which(degree(dp_mst) == 1)
Mouse_dm_df$tips <- F

qplot(mouse_high_c_result$C[1, ], mouse_high_c_result$C[2, ], color = degree(gp) == 1) 
qplot(mouse_high_c_result$C[1, ], mouse_high_c_result$C[2, ], color = degree(gp) == 1) 

mouse_high_c_result$C[, which(degree(gp) == 1)]

exprs_new_cds@reducedDimS <- t(Mouse_dm_df[, 1:2])
qplot(exprs_new_cds@reducedDimS[1, pData(exprs_new_cds)$cluster_id == 4], exprs_new_cds@reducedDimS[2, pData(exprs_new_cds)$cluster_id == 4])
which.min(exprs_new_cds@reducedDimS[1, pData(exprs_new_cds)$cluster_id == 4]) #179
which(pData(exprs_new_cds)$cluster_id == 4)[179] #cell 5440
exprs_new_cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex['cell_5440', ] #182

colnames(mouse_high_c_result$C) <- paste("Y_", 1:ncol(mouse_high_c_result$C), sep = "")
DCs <- t(Mouse_dm_df[, 1:2])
colnames(DCs) <- paste("cell_", 1:ncol(DCs), sep = "")

monocle:::reducedDimW(exprs_new_cds) <- DCs
monocle:::reducedDimS(exprs_new_cds) <- DCs
monocle:::reducedDimK(exprs_new_cds) <- mouse_high_c_result$C
exprs_new_cds@auxOrderingData[["DDRTree"]]$objective_vals <- tail(mouse_high_c_result$objs, 1)
exprs_new_cds@auxOrderingData[["DDRTree"]]$W <- mouse_high_c_result$W
exprs_new_cds@auxOrderingData[["DDRTree"]]$P <- mouse_high_c_result$P

adjusted_K <- Matrix::t(reducedDimK(exprs_new_cds))
dp <- as.matrix(dist(adjusted_K))
cellPairwiseDistances(exprs_new_cds) <- dp

dimnames(mouse_high_c_result$W) <- list(paste('Y_', 1:ncol(mouse_high_c_result$W), sep = ''), paste('Y_', 1:ncol(mouse_high_c_result$W), sep = ''))
gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)

dp_mst <- minimum.spanning.tree(gp)
neighborhood(dp_mst, nodes = 'Y_182')
intersect( exprs_new_cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex[which(pData(exprs_new_cds)$cluster_id == 4), 1], c(164, 177))
#break 182 -> 177
dp_mst <- delete_edges(dp_mst, 'Y_182|Y_177')
#connect 142 -> 278: 
dp_mst <- add_edges(dp_mst, c('Y_142', 'Y_278'))

minSpanningTree(exprs_new_cds) <- dp_mst
exprs_new_cds@dim_reduce_type <- "DDRTree"

exprs_new_cds <- monocle:::findNearestPointOnMST(exprs_new_cds)

exprs_new_cds <- orderCells(exprs_new_cds)
principal_tree_pseudo <- extract_ddrtree_ordering(exprs_new_cds, 'Y_182')

pData(exprs_new_cds)$Pseudotime <- principal_tree_pseudo[paste('Y_', exprs_new_cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex, sep = ''), 'pseudo_time']
pData(exprs_new_cds)$State <- principal_tree_pseudo[paste('Y_', exprs_new_cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex, sep = ''), 'cell_state']

plot_cell_trajectory(exprs_new_cds, color_by = 'Pseudotime')
plot_cell_trajectory(exprs_new_cds, color_by = 'State')

data_expr <- data.frame(Pseudotime = pData(exprs_new_cds)$Pseudotime, exprs = as.matrix(exprs(exprs_new_cds)[2, ]), cluster_id = pData(exprs_new_cds)$cluster_id, group = 1)

data_expr <- data.frame(Pseudotime = pData(exprs_new_cds)$Pseudotime, exprs = colSums(as.matrix(exprs(exprs_new_cds))), cluster_id = pData(exprs_new_cds)$cluster_id, group = 1)
qplot(Pseudotime, exprs, data = data_expr, color = cluster_id) + geom_smooth(aes(group = 1), span = 0.8, method = "loess")

data_expr <- data.frame(Pseudotime = pData(exprs_new_cds)$Pseudotime, exprs = as.matrix(colSums(exprs(exprs_new_cds)[50:70, ])), cluster_id = pData(exprs_new_cds)$cluster_id, group = 1)
qplot(Pseudotime, exprs, data = data_expr, color = cluster_id) + geom_smooth(aes(group = 1), span = 0.8, method = "loess")

save(file = './RData/exprs_new_cds', exprs_new_cds)

# Y_5440 
# 182 

root_cell <- 'Y_182'

exprs_new_cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex

dp <- as.matrix(dist(t(exprs_new_cds@reducedDimK))) #cellPairwiseDistances(exprs_new_cds)
dimnames(dp) <- dimnames(mouse_high_c_result$W)
dp_mst <- minSpanningTree(exprs_new_cds)

#delete links: 
neighborhood(dp_mst, nodes = 'Y_182')
#dp_mst_ori <- dp_mst
# dp_mst <- delete_edges(dp_mst, 'Y_172|Y_73')
# dp_mst <- delete_edges(dp_mst, 'Y_172|Y_286')
# dp_mst <- delete_edges(dp_mst, 'Y_172|Y_34')

curr_state <- 1

res <- list(subtree = dp_mst, root = root_cell)

states = rep(1, ncol(dp))
names(states) <- V(dp_mst)$name

pseudotimes = rep(0, ncol(dp))
names(pseudotimes) <- V(dp_mst)$name

parents = rep(NA, ncol(dp))
names(parents) <- V(dp_mst)$name

mst_traversal <- graph.dfs(gp,
                           root=root_cell,
                           neimode = "all",
                           unreachable=T,
                           father=TRUE)
mst_traversal$father <- as.numeric(mst_traversal$father)
curr_state <- 1

for (i in 1:length(mst_traversal$order)){
  print(i)
  curr_node = mst_traversal$order[i]
  curr_node_name = V(dp_mst)[curr_node]$name
  
  if (is.na(mst_traversal$father[curr_node]) == FALSE){
    parent_node = mst_traversal$father[curr_node]
    parent_node_name = V(dp_mst)[parent_node]$name
    parent_node_pseudotime = pseudotimes[parent_node_name]
    parent_node_state = states[parent_node_name]
    curr_node_pseudotime = parent_node_pseudotime + dp[curr_node_name, parent_node_name]
    if (degree(dp_mst, v=parent_node_name) > 2){
      curr_state <- curr_state + 1
    }
  }else{
    parent_node = NA
    parent_node_name = NA
    curr_node_pseudotime = 0
  }
  
  curr_node_state = curr_state
  pseudotimes[curr_node_name] <- curr_node_pseudotime
  states[curr_node_name] <- curr_node_state
  parents[curr_node_name] <- parent_node_name
}

ordering_df <- data.frame(sample_name = names(states),
                          cell_state = factor(states),
                          pseudo_time = as.vector(pseudotimes),
                          parent = parents)
row.names(ordering_df) <- ordering_df$sample_name

low_dim_df <- data.frame(x = exprs_new_cds@reducedDimK[1, ], y =  exprs_new_cds@reducedDimK[2, ], color = 1:300 %in% c(34, 73, 172, 286))
mouse_flatten_edge_df$Pseudotime <- ordering_df[as.character(mouse_flatten_edge_df$cell_name), 'pseudo_time']
ggplot(edge_df) + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), #color = "graph_cluster"), 
                               size=1, linetype="solid", na.rm=TRUE)  + #facet_wrap(~graph_cluster) +
  geom_point(aes(x = x, y = y, color = color), data = low_dim_df)

ggplot(edge_df) + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), #color = "graph_cluster"), 
                               size=1, linetype="solid", na.rm=TRUE)  + #facet_wrap(~graph_cluster) +
  geom_point(aes(x, y, size = Pseudotime), data = mouse_flatten_edge_df) + # + geom_label(aes(x, y,label = cell_name), data = flatten_edge_df)
  geom_point(aes(DC1, DC2, color = color), data = Mouse_dm_df)

ggplot(edge_df) + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", 
                                          xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), #, color = "graph_cluster"), 
                               size=1, linetype="solid", na.rm=TRUE) + geom_label(aes(x, y,label = cell_name), data = mouse_flatten_edge_df)

pData(exprs_new_cds)$Pseudotime <- pseudotimes[exprs_new_cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex[, 1]]

ggplot(edge_df) + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2", color = "graph_cluster"), 
                               size=1, linetype="solid", na.rm=TRUE)  + #facet_wrap(~graph_cluster) + 
  geom_point(aes(x, y, size = Pseudotime, color = as.factor(graph_cluster)), data = flatten_edge_df) # + geom_label(aes(x, y,label = cell_name), data = flatten_edge_df)

data_expr <- data.frame(Pseudotime = pData(exprs_new_cds)$Pseudotime, exprs = as.matrix(exprs(exprs_new_cds)[1, ]), cluster_id = pData(exprs_new_cds)$cluster_id, group = 1)
qplot(Pseudotime, exprs, data = data_expr, color = cluster_id) + geom_smooth(aes(group = 1), span = 0.8, method = "loess")

######################################################################################################################################################################################
# human dataset
######################################################################################################################################################################################
load('./RData/Vijay_cds.RData')

#calculate local density
human_hc_Dist <- dist(t(cds@reducedDimS))
human_hc_Clust <- densityClust(human_hc_Dist, gaussian=TRUE)
plot(human_hc_Clust) # Inspect clustering attributes to define thresholds

human_hc_Clust <- findClusters(human_hc_Clust, rho=2, delta=0.002)

plotMDS(human_hc_Clust)

human_dm_df <- data.frame(DC1 = cds@reducedDimS[1, ], DC2 = cds@reducedDimS[2, ], color = as.factor(pData(cds)$cluster_id))
human_dm_df$rho <- human_hc_Clust$rho
human_dm_df$landmark <- F
#downsampling 
human_dm_df$landmark[sample(1:nrow(human_dm_df), size = 400, prob = exp(-(human_dm_df$rho))^2 )] <- T

qplot(DC1, DC2, color = landmark, data = human_dm_df) + facet_wrap(~landmark)
# writeMat('./mat_data/human_landmarks', human_landmarks = as.matrix(subset(human_dm_df, landmark == T )[, 1:2] ))
#run SGL-tree 

#no mannual adjustment
human_high_c_result <- readMat('/Users/xqiu/Dropbox (Personal)/Projects/SPL+L1graph/tree-l1-code/human_high_c_result.mat') #
#project cells to the principal graph 
cds@reducedDimS <- t(human_dm_df[, 1:2])
# qplot(cds@reducedDimS[1, pData(cds)$cluster_id == 4], cds@reducedDimS[2, pData(cds)$cluster_id == 4])
# which.min(cds@reducedDimS[1, pData(cds)$cluster_id == 4]) #179
# which(pData(cds)$cluster_id == 4)[179] #cell 5440
# cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex['cell_5440', ] #182
# 
colnames(human_high_c_result$C) <- paste("Y_", 1:ncol(human_high_c_result$C), sep = "")
DCs <- t(human_dm_df[, 1:2])
colnames(DCs) <- paste("cell_", 1:ncol(DCs), sep = "")

monocle:::reducedDimW(cds) <- DCs
monocle:::reducedDimS(cds) <- DCs
monocle:::reducedDimK(cds) <- human_high_c_result$C
cds@auxOrderingData[["DDRTree"]]$objective_vals <- tail(human_high_c_result$objs, 1)
cds@auxOrderingData[["DDRTree"]]$W <- human_high_c_result$W
cds@auxOrderingData[["DDRTree"]]$P <- human_high_c_result$P

adjusted_K <- Matrix::t(reducedDimK(cds))
dp <- as.matrix(dist(adjusted_K))
cellPairwiseDistances(cds) <- dp

dimnames(human_high_c_result$W) <- list(paste('Y_', 1:ncol(human_high_c_result$W), sep = ''), paste('Y_', 1:ncol(human_high_c_result$W), sep = ''))
gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)

dp_mst <- minimum.spanning.tree(gp)
# neighborhood(dp_mst, nodes = 'Y_182')
# intersect( cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex[which(pData(cds)$cluster_id == 4), 1], c(164, 177))
# #break 182 -> 177
# dp_mst <- delete_edges(dp_mst, 'Y_182|Y_177')
# #connect 142 -> 278: 
# dp_mst <- add_edges(dp_mst, c('Y_142', 'Y_278'))

minSpanningTree(cds) <- dp_mst
cds@dim_reduce_type <- "DDRTree"

cds <- monocle:::findNearestPointOnMST(cds)

cds <- orderCells(cds)

#identify the initial cell: 

principal_tree_pseudo <- extract_ddrtree_ordering(cds, 'Y_182')
pData(cds)$Pseudotime <- principal_tree_pseudo[paste('Y_', cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex, sep = ''), 'pseudo_time']
pData(cds)$State <- principal_tree_pseudo[paste('Y_', cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex, sep = ''), 'cell_state']

plot_cell_trajectory(cds, color_by = 'Pseudotime')
plot_cell_trajectory(cds, color_by = 'State')

data_expr <- data.frame(Pseudotime = pData(cds)$Pseudotime, exprs = as.matrix(exprs(cds)[2, ]), cluster_id = pData(cds)$cluster_id, group = 1)

data_expr <- data.frame(Pseudotime = pData(cds)$Pseudotime, exprs = colSums(as.matrix(exprs(cds))), cluster_id = pData(cds)$cluster_id, group = 1)
qplot(Pseudotime, exprs, data = data_expr, color = cluster_id) + geom_smooth(aes(group = 1), span = 0.8, method = "loess")

data_expr <- data.frame(Pseudotime = pData(cds)$Pseudotime, exprs = as.matrix(colSums(exprs(cds)[50:70, ])), cluster_id = pData(cds)$cluster_id, group = 1)
qplot(Pseudotime, exprs, data = data_expr, color = cluster_id) + geom_smooth(aes(group = 1), span = 0.8, method = "loess")
#calculate pseudotime 

hunman_cds <- cds
save(file = './RData/human_cds', hunman_cds)
