################################################################################################################################################
# load all the necessary packages
################################################################################################################################################
library(monocle)
library(dpt)
library(destiny)
library(diffusionMap)
library(igraph)
library(L1Graph)
library(R.matlab)

################################################################################################################################################
# analyze mouse data 
################################################################################################################################################
cell_cycle_high_c <- t(read.csv('/Users/xqiu/Downloads/Mouse_scaling_coefficients (1).txt.gz', sep = '\t', header = F))
source('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/scripts/function.R', echo = T)

run_dpt <- function(data, branching = T, normalize = T, root = NULL){
  ts <- Transitions(data)
  M <- dpt:::propagation_matrix(ts)
  
  if(is.null(root))
    pt <- dpt(ts, branching = branching)
  else
    pt <- dpt(ts, branching = branching, root = root)
  
  ev <- eigen(as.matrix(ts@transitions), TRUE)$vectors
  dm <- as.data.frame(ev[, -1])
  colnames(dm) <- paste0('DC', seq_len(ncol(dm)))
  
  p1 <- qplot(DC1, DC2, data = dm, colour = pData(cds)$Hours)
  # plot_dpt(ts, pt) # DPT and average path , 1:2
  
  dp_res <- list(dm = dm, pt = pt, ts = ts, M = M, ev = ev, p1 = p1)
  
  return(dp_res)
}

cell_cycle_high_c_center <- scale(t(cell_cycle_high_c), center = TRUE, scale = F)

cell_cycle_high_c_center_dpt <- run_dpt(cell_cycle_high_c_center, normalize = F)
cell_cycle_high_c_center_dm <- DiffusionMap(cell_cycle_high_c_center)

Mouse_cluster_ids <- t(read.csv('/Users/xqiu/Downloads/Mouse_cluster_ids (1).txt', sep = '\t', header = F))

Mouse_dpt_2d_coord <- data.frame(DC1 = cell_cycle_high_c_center_dm$DC1[2:7344], DC2 = cell_cycle_high_c_center_dm$DC2[2:7344])
write.table(file = './csv_data/mouse_dpt_2d_coord.txt', Mouse_dpt_2d_coord, sep = '\t', row.names = F, col.names = F, quote = F)

Mouse_dm_df <- data.frame(DC1 = cell_cycle_high_c_center_dm$DC1, DC2 = cell_cycle_high_c_center_dm$DC2, color = as.factor(Mouse_cluster_ids[1, ]))
qplot(DC1, DC2, data = Mouse_dm_df, color = color) + facet_wrap(~Mouse_dm_df$color)

qplot(DC1, DC2, data = Mouse_dm_df, color = color) #+ facet_wrap(~Mouse_dm_df$color)

qplot(DC1, DC2, data = Mouse_dm_df, color = pData(exprs_new_cds)$Cluster) #+ facet_wrap(~Mouse_dm_df$color)

################################################################################################################################################
# calculate the pseudotime based on our result 
################################################################################################################################################
mouse_high_c_result <- readMat('/Users/xqiu/Dropbox (Personal)/Projects/SPL+L1graph/tree-l1-code/mouse_high_c_result.mat')
names(mouse_high_c_result)

W <- mouse_high_c_result$W
mouse_high_c_result$W[lower.tri(mouse_high_c_result$W)] <- 0
exprs_new_cds@reducedDimK <- mouse_high_c_result$C
cds <- findNearestPointOnMST(exprs_new_cds)

W_weighted <- which(as.matrix(mouse_high_c_result$W) > 0, arr.ind = T)
edge_df <- data.frame(source_prin_graph_dim_1 = mouse_high_c_result$C[1, W_weighted[, 1]],
                      source_prin_graph_dim_2 = mouse_high_c_result$C[2, W_weighted[, 1]],
                      target_prin_graph_dim_1 = mouse_high_c_result$C[1, W_weighted[, 2]],
                      target_prin_graph_dim_2 = mouse_high_c_result$C[2, W_weighted[, 2]] 
)

mouse_flatten_cell_names <- paste('Y_', as.vector(W_weighted), sep = '')
mouse_flatten_edge_df <- data.frame(x = c(source_prin_graph_dim_1 = mouse_high_c_result$C[1, W_weighted[, 1]], source_prin_graph_dim_2 = mouse_high_c_result$C[1, W_weighted[, 2]]),
                              y = c(target_prin_graph_dim_1 = mouse_high_c_result$C[2, W_weighted[, 1]], target_prin_graph_dim_2 = mouse_high_c_result$C[2, W_weighted[, 2]]),
                              cell_name = mouse_flatten_cell_names)

mouse_flatten_edge_df$Pseudotime <- pseudotime[1, as.character(mouse_flatten_edge_df$cell_name)]
ggplot(edge_df[-c(187, 95), ]) + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", 
                                          xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), #, color = "graph_cluster"), 
                               size=1, linetype="solid", na.rm=TRUE)  + #facet_wrap(~graph_cluster) + 
  geom_point(aes(DC1, DC2, color = color), data = Mouse_dm_df) # + geom_label(aes(x, y,label = cell_name), data = flatten_edge_df)

#remove two extra edges:
#unnecessary edges: Y_73 <-> Y_117; Y_137 <-> 
subset(edge_df, (source_prin_graph_dim_1 < -0.025 & target_prin_graph_dim_1 > 0 & source_prin_graph_dim_2 < -0.025 & target_prin_graph_dim_2 > -0.025))
subset(edge_df, (source_prin_graph_dim_1 < -0.125 & target_prin_graph_dim_1 > -0.125 & source_prin_graph_dim_2 < 0.075 & target_prin_graph_dim_2 > 0.075))

gp <- graph.adjacency(mouse_high_c_result$W, mode = "directed", weighted = TRUE)

neighborhood(gp, nodes = '137')

ggplot(edge_df) + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", 
                                          xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), #, color = "graph_cluster"), 
                               size=1, linetype="solid", na.rm=TRUE) + geom_label(aes(x, y,label = cell_name), data = mouse_flatten_edge_df)

pData(exprs_new_cds)$cluster_id <- Mouse_dm_df$color

# state 4 is interphase G1 stage and the first cell can be defined as pseudo-cell cycle 0
exprs_new_cds@reducedDimS <- t(Mouse_dm_df[, 1:2])
qplot(exprs_new_cds@reducedDimS[1, pData(exprs_new_cds)$cluster_id == 4], exprs_new_cds@reducedDimS[2, pData(exprs_new_cds)$cluster_id == 4])
which.min(exprs_new_cds@reducedDimS[1, pData(exprs_new_cds)$cluster_id == 4]) #179

exprs_new_cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex['cell_4096', ]
# cell_4096 
# 181 

root_cell <- 'Y_181'

exprs_new_cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex

dp <- as.matrix(dist(t(exprs_new_cds@reducedDimK))) #cellPairwiseDistances(exprs_new_cds)
dimnames(dp) <- dimnames(W)
dp_mst <- minSpanningTree(exprs_new_cds)

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

flatten_edge_df$Pseudotime <- ordering_df[as.character(flatten_edge_df$cell_name), 'pseudo_time']
ggplot(edge_df) + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2", color = "graph_cluster"), 
                               size=1, linetype="solid", na.rm=TRUE)  + #facet_wrap(~graph_cluster) + 
  geom_point(aes(x, y, size = Pseudotime, color = as.factor(graph_cluster)), data = flatten_edge_df) # + geom_label(aes(x, y,label = cell_name), data = flatten_edge_df)


################################################################################################################################################
# try density peaks: 
################################################################################################################################################

mouse_hc_Dist <- dist(cell_cycle_high_c_center_dpt$dm[, 1:2])
mouse_hc_Clust <- densityClust(mouse_hc_Dist, gaussian=TRUE)
plot(mouse_hc_Clust) # Inspect clustering attributes to define thresholds

mouse_hc_Clust <- findClusters(mouse_hc_Clust, rho=2, delta=0.002)

qplot(mouse_hc_Clust$rho, mouse_hc_Clust$delta) + geom_vline(aes(xintercept = 5)) + geom_hline(aes(yintercept = 0.02))  + geom_hline(aes(yintercept = 0.002))
qplot(DC1, DC2, data = Mouse_dm_df, color = factor(mouse_hc_Clust$clusters)) #+ facet_wrap(~Mouse_dm_df$color)

# select the peaks: 
mouse_hc_Clust$rho[mouse_hc_Clust$de]

rho_delta_data <- rep(NA, length(mouse_hc_Clust$rho))
rho_delta_data[mouse_hc_Clust$delta > 0.002 &  mouse_hc_Clust$rho > 5] <- 'peaks'
rho_delta_data[mouse_hc_Clust$delta > 0.02 &  mouse_hc_Clust$rho < 5] <- 'peaks'

qplot(mouse_hc_Clust$rho, mouse_hc_Clust$delta, color = rho_delta_data) #+ geom_vline(aes(xintercept = 5)) + geom_hline(aes(yintercept = 0.02))  + geom_hline(aes(yintercept = 0.002))

Mouse_dm_df$rho_delta_data <- rho_delta_data
qplot(DC1, DC2, data = Mouse_dm_df, color = rho_delta_data, facets = "~rho_delta_data") #+ facet_wrap(~rho_delta_data)

C0 <- t(subset(Mouse_dm_df, rho_delta_data == "peaks")[, 1:2])
Nz <- ncol(C0)

G <- get_knn(C0, 3)

human_cluster_ids <- read.csv('/Users/xqiu/Downloads/Human_cluster_ids.txt', header = F)

maxiter <- 20
eps <- 1.0000e-05
gstruct <- 'l1-graph' #'span-tree'
gamma <- 1 #smoothness
sigma <- 0.001 #
lambda <- 100 #L1 g
nn <- 3
verbose = T

human_cell_cycle_high_c_center_principal_graph <- principal_graph(t(X[, 1:2]), C0, G$G, maxiter = 100, #t(human_cell_cycle_high_c_center_dpt@eigenvectors[, 1:2])
                                                                  eps = eps, gstruct = gstruct,
                                                                  lambda = lambda, gamma = gamma,
                                                                  sigma = sigma, nn = 3, verbose = T)

tail(human_cluster_ids, 10) #miss the last 1 and have an extra 0 at the beginning 
tail(baz, 10)

dm_df <- data.frame(x = cell_cycle_high_c_center_dpt$dm[, 1], y = cell_cycle_high_c_center_dpt$dm[, 1])
ggplot(aes(x, y), data = dm_df, color = 'black') 

################################################################################################################################################
# get the pseudotime based on the result
################################################################################################################################################

PCA <- function(data, ...) {
  res <- prcomp(t(data), center = T, scale = F)
  res$x
}

cell_cycle_high_c_center_ICA <- ICA(cell_cycle_high_c_center)
cell_cycle_high_c_center_PCA <- PCA(cell_cycle_high_c_center)
cell_cycle_high_c_center_LLE <- lle::lle(cell_cycle_high_c_center, m = 2, k = 5)

colnames(cell_cycle_high_c) <- paste('cell', 1:ncol(cell_cycle_high_c), sep = '_')
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = row.names(cell_cycle_high_c), row.names = row.names(cell_cycle_high_c)))
pd <- new("AnnotatedDataFrame", data = data.frame(cell = colnames(cell_cycle_high_c),
                                                  row.names = colnames(cell_cycle_high_c)))

exprs_new_cds <- newCellDataSet(as(as.matrix(t(cell_cycle_high_c_center)), "sparseMatrix"), #he cofactor parameter is used for arcsinh
                                phenoData = pd, 
                                featureData = fd, 
                                expressionFamily=gaussianff(), 
                                lowerDetectionLimit=1)

################################################################################################################################################
# run tSNE
################################################################################################################################################
#1. determine how many pca dimension you want:
exprs_new_cds <- detectGenes(exprs_new_cds)
exprs_new_cds <- estimateSizeFactors(exprs_new_cds)
exprs_new_cds <- estimateDispersions(exprs_new_cds)

fData(exprs_new_cds)$use_for_ordering <- T#fData(exprs_new_cds)$num_cells_expressed > round(ncol(exprs_new_cds) / 10)
exprs_new_cds@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(exprs_new_cds, return_all = F, norm_method = 'none', pseudo_expr = 1)

#2. run reduceDimension with tSNE as the reduction_method
exprs_new_cds <- reduceDimension(exprs_new_cds, max_components=2, norm_method = 'none', pseudo_expr = 1, num_dim = 15, reduction_method = 'tSNE', verbose = T)

#3. initial run of clusterCells_Density_Peak
exprs_new_cds <- clusterCells_Density_Peak(exprs_new_cds, verbose = T)

#4. check the clusters (there are three clusters)
plot_cell_clusters(exprs_new_cds, color_by = 'as.factor(Cluster)', show_density = F)
plot_cell_clusters(exprs_new_cds, color_by = 'as.factor(Time)', show_density = F)

#5. also check the decision plot
plot_rho_delta(exprs_new_cds, rho_threshold = 3.8, delta_threshold = 12)
plot_cell_clusters(exprs_new_cds, color_by = 'as.factor(Cluster)', show_density = F, rho_threshold = 3.8, delta_threshold = 12)
plot_cell_clusters(exprs_new_cds, color_by = 'as.factor(Time)', show_density = F, rho_threshold = 3.8, delta_threshold = 12)

#6. re-run cluster and skipping calculating the rho_sigma
exprs_new_cds <- clusterCells_Density_Peak(exprs_new_cds, verbose = T,  rho_threshold = 3.8, delta_threshold = 12, skip_rho_sigma = T)

#7. make the final clustering plot:
plot_cell_clusters(exprs_new_cds, color_by = 'as.factor(Cluster)', show_density = F)
plot_cell_clusters(exprs_new_cds, color_by = 'as.factor(Time)', show_density = F)

################################################################################################################################################
# run l1 graph
################################################################################################################################################

exprs_new_cds <- reduceDimension(exprs_new_cds, reduction_method = 'L1-span-tree', norm_method = 'none', verbose = T, pseudo_expr = 0)
exprs_new_cds <- orderCells(exprs_new_cds)

plot_cell_trajectory(exprs_new_cds)

################################################################################################################################################
# analyze human data 
################################################################################################################################################
human_cell_cycle_high_c <- t(read.csv('./csv_data/Human_scaling_coefficients.txt.gz', sep = '\t', header = F))
human_cell_cycle_high_c <- human_cell_cycle_high_c[1:83, ]
human_cell_cycle_high_c_label <- t(read.csv('./csv_data/Human_labels.txt', sep = '\t', header = F))

human_cell_cycle_high_c_center <- scale(t(human_cell_cycle_high_c), center = TRUE, scale = TRUE)

human_cell_cycle_high_c_center_dpt <- DiffusionMap(human_cell_cycle_high_c_center[1:4271, 1:83])
human_cell_cycle_high_c_center_dm <- DM(t(human_cell_cycle_high_c_center))
qplot(human_cell_cycle_high_c_center_dpt$DC1, human_cell_cycle_high_c_center_dpt$DC2, color = human_cell_cycle_high_c_label[1, ])

dpt_2d_coord <- data.frame(dm1 = human_cell_cycle_high_c_center_dpt$DC1, dm2 = human_cell_cycle_high_c_center_dpt$DC2)

write.table(file = './csv_data/dpt_2d_coord.txt', dpt_2d_coord, sep = '\t', row.names = F, col.names = F, quote = F)
writeMat('./mat_data/dpt_2d_coord', dpt_2d_coord = dpt_2d_coord)
# run l1 graph just on the reduced dimension: 
ncenter <- 500
maxiter <- 20
eps <- 1.0000e-05
gstruct <- 'l1-graph' #'span-tree'
gamma <- 1 #smoothness
sigma <- 0.001 #
lambda <- 100 #L1 g
nn <- 3
verbose = T

X <- t(human_cell_cycle_high_c_center_dpt@eigenvectors[, 1:2])
# D <- nrow(X); N <- ncol(X)
# Z <- X

#kmeans to find the ncenters: 
load('/Users/xqiu/Downloads/for_xiaojie.dat')
X <- baz
kmenas_res <- kmeans(X[, 1:2], ncenter, nstart = 20)

C0 <- t(kmenas_res$centers)
Nz <- ncol(C0)

G <- get_knn(C0, 3)

human_cluster_ids <- read.csv('/Users/xqiu/Downloads/Human_cluster_ids.txt', header = F)

human_cell_cycle_high_c_center_principal_graph <- principal_graph(t(X[, 1:2]), C0, G$G, maxiter = 100, #t(human_cell_cycle_high_c_center_dpt@eigenvectors[, 1:2])
                                             eps = eps, gstruct = gstruct,
                                             lambda = lambda, gamma = gamma,
                                             sigma = sigma, nn = 3, verbose = T)

tail(human_cluster_ids, 10) #miss the last 1 and have an extra 0 at the beginning 
tail(baz, 10)

dm_df <- data.frame(x = human_cell_cycle_high_c_center_dpt@eigenvectors[, 1], y = human_cell_cycle_high_c_center_dpt@eigenvectors[, 2], cluster  = factor(as.character(baz$colour)))
ggplot(aes(x, y), data = dm_df) + geom_point(aes(color = cluster)) + facet_wrap(~cluster)

qplot(human_cell_cycle_high_c_center_dpt@eigenvectors[, 1], human_cell_cycle_high_c_center_dpt@eigenvectors[, 2], color = as.factor(human_cluster_ids[, 1])) + 
  facet_wrap(~ factor(as.character(human_cluster_ids[, 1])))

qplot(human_cell_cycle_high_c_center_principal_graph$C[1, ], human_cell_cycle_high_c_center_principal_graph$C[2, ], color = I('black')) + 
         # geom_point(aes(human_cell_cycle_high_c_center_dpt@eigenvectors[, 1], human_cell_cycle_high_c_center_dpt@eigenvectors[, 2], 
         #                color = baz$colour)) + ggtitle('spanning-tree') + 
  geom_point(aes(kmenas_res$centers[, 1], kmenas_res$centers[, 2]), color = I('yellow'))

G <- get_knn(C0, 5)
W <- human_cell_cycle_high_c_center_principal_graph$W #reducedDimW(human_exprs_new_cds)

mat_res <- readMat('/Users/xqiu/Dropbox (Personal)/Projects/SPL+L1graph/tree-l1-code/matlab.mat')

W <- as.matrix(mat_res$W)
dimnames(W) <- list(paste('Y_', 1:nrow(W), sep = ''), paste('Y_', 1:nrow(W), sep = ''))
W[W < 1e-5] <- 0
gp <- graph.adjacency(W, mode = "undirected", weighted = TRUE)
# dp_mst <- minimum.spanning.tree(gp)
minSpanningTree(human_exprs_new_cds) <- gp
human_exprs_new_cds@dim_reduce_type <- "DDRTree"
human_exprs_new_cds <- findNearestPointOnMST(human_exprs_new_cds)
plot_cell_trajectory(human_exprs_new_cds, show_tree = T)

W <- as.matrix(mat_res$W)
dimnames(W) <- list(paste('Y_', 1:nrow(W), sep = ''), paste('Y_', 1:nrow(W), sep = ''))
cds <- human_exprs_new_cds 
l1_graph_res <- mat_res
colnames(l1_graph_res$C) <- paste("Y_", 1:ncol(l1_graph_res$C), sep = "")
DCs <- l1_graph_res$X
row.names(DCs) <- colnames(cds)
reducedDimW(cds) <- DCs
reducedDimS(cds) <- t(DCs) #1_graph_res$X
reducedDimK(cds) <- l1_graph_res$C
cds@auxOrderingData[["DDRTree"]]$objective_vals <- tail(l1_graph_res$objs, 1)
cds@auxOrderingData[["DDRTree"]]$W <- as.matrix(l1_graph_res$W)
cds@auxOrderingData[["DDRTree"]]$P <- l1_graph_res$P
adjusted_K <- Matrix::t(reducedDimK(cds))

# dp <- as.matrix(dist(adjusted_K))
cellPairwiseDistances(cds) <- as.matrix(mat_res$W)# dp
# gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
# dp_mst <- minimum.spanning.tree(gp)
# minSpanningTree(cds) <- dp_mst
W[W < 1e-5] <- 0
gp <- graph.adjacency(W, mode = "undirected", weighted = TRUE)
# dp_mst <- minimum.spanning.tree(gp)
minSpanningTree(cds) <- gp

cds@dim_reduce_type <- "DDRTree"
cds <- findNearestPointOnMST(cds)
cds <- orderCells(cds)
pData(cds)$cluster_id <- baz$colour
plot_cell_trajectory(cds, show_tree = T)

W_weighted <- which(W > 0, arr.ind = T)
edge_df <- data.frame(source_prin_graph_dim_1 = mat_res$C[1, W_weighted[, 1]],
                      source_prin_graph_dim_2 = mat_res$C[2, W_weighted[, 1]],
                      target_prin_graph_dim_1 = mat_res$C[1, W_weighted[, 2]],
                      target_prin_graph_dim_2 = mat_res$C[2, W_weighted[, 2]] 
                      )
flatten_cell_names <- paste('Y_', as.vector(W_weighted), sep = '')
flatten_edge_df <- data.frame(x = c(source_prin_graph_dim_1 = mat_res$C[1, W_weighted[, 1]], source_prin_graph_dim_2 = mat_res$C[1, W_weighted[, 2]]),
                              y = c(target_prin_graph_dim_1 = mat_res$C[2, W_weighted[, 1]], target_prin_graph_dim_2 = mat_res$C[2, W_weighted[, 2]]),
                              cell_name = flatten_cell_names)
dp_mst <- minSpanningTree(cds)
g_clust_res <- clusters(dp_mst)
g_clust_res <- clusters(gp)
flatten_edge_df$graph_cluster <- g_clust_res$membership[flatten_cell_names]

edge_df$graph_cluster <- as.factor(g_clust_res$membership[row.names(W_weighted)])
ggplot(edge_df) + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2", color = "graph_cluster"), 
                        size=1, linetype="solid", na.rm=TRUE) + facet_wrap(~graph_cluster) + 
  geom_point(aes(x, y, size = 2, color = as.factor(graph_cluster)), data = flatten_edge_df) # + geom_label(aes(x, y,label = cell_name), data = flatten_edge_df)

######################################################################################################################################################################################
# mannually calculate the pseudotime 
######################################################################################################################################################################################
# identify the start point: 
# which(V(dp_mst)[degree(dp_mst) == 3] %in%  subset(flatten_edge_df, (x > 0.075 & x < 0.10) & (y > -0.1 & y < -0.05))$cell_name)
# degree(dp_mst)[as.character(subset(flatten_edge_df, (x > 0.075 & x < 0.10) & (y > -0.1 & y < -0.05))$cell_name)]
# 
# degree(dp_mst)[as.character(subset(flatten_edge_df, (x > 0.075 & x < 0.10) & (y > -0.1 & y < -0.05))$cell_name)]
# 
# ggplot(edge_df) + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), 
#                                size=1, linetype="solid", na.rm=TRUE) + 
#   geom_text(aes(x, y,label = cell_name), size = 3, 
#             data = subset(flatten_edge_df, (x > 0.075 & x < 0.10) & (y > -0.1 & y < -0.05)), hjust = 0, nudge_x = 0, angle = 45)
# 
# degree(dp_mst)[as.character(subset(flatten_edge_df, (x > 0.075 & x < 0.10) & (y > -0.1 & y < -0.05))$cell_name)]

############################################################################################################################################
# use the simple to calculate pseudotime for a complicated graph:
# 1. connect nearest points between two connected components in the graph 
############################################################################################################################################

dp <- as.matrix(dist(t(cds@reducedDimK)))
g_clust_res <- clusters(gp)
for(cells in V(gp)[degree(gp)==1]$name) {
  min_val <- min(dp[cells, g_clust_res$membership != g_clust_res$membership[cells]])
  if(min_val < mean(dp)) {
    min_cell <- which.min(dp[cells, g_clust_res$membership != g_clust_res$membership[cells]])
    gp <- add_edges(gp, c(cells, names(min_cell)) )
  }
}

g_clust_res$membership[V(gp)[degree(gp)==1]$name]

for(cells in V(gp)[degree(gp)==1]$name) {
  dp_subset <- dp[V(gp)[degree(gp)==1]$name, g_clust_res$membership != 1] 
  min_val <- min(dp_subset)
  min_cell <- which(dp_subset == min_val, arr.ind = T)
  gp <- add_edges(gp, c(row.names(dp_subset)[min_cell[1]], colnames(dp_subset)[min_cell[2]]) )
}

# state 6 is G1 stage and the first cell can be defined as pseudo-cell cycle 0
qplot(cds@reducedDimS[1, pData(cds)$cluster_id == 6], cds@reducedDimS[2, pData(cds)$cluster_id == 6])
which.min(cds@reducedDimS[2, pData(cds)$cluster_id == 6])
# cell_4096 
# 90  
cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex['cell_4096', ]
# cell_4096 
# 181 

root_cell <- 'Y_181'

cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex

dp <- as.matrix(dist(t(cds@reducedDimK))) #cellPairwiseDistances(cds)
dimnames(dp) <- dimnames(W)
dp_mst <- minSpanningTree(cds)

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

flatten_edge_df$Pseudotime <- ordering_df[as.character(flatten_edge_df$cell_name), 'pseudo_time']
ggplot(edge_df) + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2", color = "graph_cluster"), 
                               size=1, linetype="solid", na.rm=TRUE)  + #facet_wrap(~graph_cluster) + 
  geom_point(aes(x, y, size = Pseudotime, color = as.factor(graph_cluster)), data = flatten_edge_df) # + geom_label(aes(x, y,label = cell_name), data = flatten_edge_df)

######################################################################################################################################################################################
# use the shortest path for pseudotime ordering 
######################################################################################################################################################################################
#break the link from M phase cell to G1: 
cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex[row.names(subset(pData(cds), cluster_id == 2)), ]

neighborhood(gp, nodes = root_cell, order = 1)[[1]]
neighborhood(gp, nodes = 'Y_430', order = 1)[[1]]

gp2 <- delete_edges(gp, 'Y_430|Y_93')
gp2 <- delete_edges(gp, 'Y_181|Y_450')
pseudotime <- distances(gp2, root_cell, V(gp2))

cell_order <- 1:length(V(gp)) 
names(cell_order) <- names(mst_traversal$order)

######################################################################################################################################################################################
# assign all centers with a pseudotime as well as each cell a pseudotime
######################################################################################################################################################################################
dp <- as.matrix(dist(t(cds@reducedDimK)))
for(cells in colnames(pseudotime)[is.infinite(pseudotime)]) {
  min_val <- min(dp[cells, colnames(pseudotime)[is.finite(pseudotime)]])
  min_cell <- which.min(dp[cells, colnames(pseudotime)[is.finite(pseudotime)]])
  pseudotime[, cells] <- pseudotime[, names(min_cell)]
}

flatten_edge_df$Pseudotime <- pseudotime[1, as.character(flatten_edge_df$cell_name)]
ggplot(edge_df) + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2", color = "graph_cluster"), 
                               size=1, linetype="solid", na.rm=TRUE)  + #facet_wrap(~graph_cluster) + 
  geom_point(aes(x, y, size = Pseudotime, color = as.factor(graph_cluster)), data = flatten_edge_df) # + geom_label(aes(x, y,label = cell_name), data = flatten_edge_df)

# each cell: 
pseudotime <- as.numeric(pseudotime)

pData(cds)$Pseudotime <- pseudotime[cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex[, 1]]

data_expr <- data.frame(Pseudotime = pData(cds)$Pseudotime, exprs = as.matrix(exprs(cds)[1, ]), cluster_id = pData(cds)$cluster_id, group = 1)
qplot(Pseudotime, exprs, data = data_expr, color = cluster_id) + geom_smooth(aes(group = 1), span = 0.8, method = "loess")

save(file = './RData/Vijay_cds.RData', cds)
# 
# # use span-tree to get the pseudotime? 
# dp <- as.matrix(dist(adjusted_K))
# cellPairwiseDistances(cds) <- as.matrix(mat_res$W)# dp
# gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
# dp_mst <- minimum.spanning.tree(gp)
# minSpanningTree(cds) <- dp_mst

# ordering_df <- plyr::arrange(ordering_df, pseudo_time)
# irisDist <- dist(human_cell_cycle_high_c_center_dpt@eigenvectors[, 1:2])
# irisClust <- densityClust(irisDist, gaussian=TRUE)
# plot(irisClust) # Inspect clustering attributes to define thresholds
# 
# irisClust <- findClusters(irisClust, rho=2, delta=2)

######################################################################################################################################################################################

PCA <- function(data, ...) {
  res <- prcomp(t(data), center = T, scale = F)
  res$x
}

human_cell_cycle_high_c_center_ICA <- ICA(human_cell_cycle_high_c_center)
human_cell_cycle_high_c_center_PCA <- PCA(human_cell_cycle_high_c)
human_cell_cycle_high_c_center_LLE <- lle::lle(human_cell_cycle_high_c_center, m = 2, k = 5)

dm_df <- data.frame(x = human_cell_cycle_high_c_center_PCA[, 1], y = human_cell_cycle_high_c_center_PCA[, 2], cluster  = as.factor(rev(human_cluster_ids[, 1])))
dm_df$cluster <- as.character(dm_df$cluster)
ggplot(aes(x, y), data = dm_df) + geom_point(aes(color = cluster)) + facet_wrap(~cluster)

colnames(human_cell_cycle_high_c) <- paste('cell', 1:ncol(human_cell_cycle_high_c), sep = '_')
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = row.names(human_cell_cycle_high_c), row.names = row.names(human_cell_cycle_high_c)))
pd <- new("AnnotatedDataFrame", data = data.frame(cell = colnames(human_cell_cycle_high_c),
                                                  row.names = colnames(human_cell_cycle_high_c)))

human_exprs_new_cds <- newCellDataSet(as(as.matrix(t(human_cell_cycle_high_c_center)), "sparseMatrix"), #he cofactor parameter is used for arcsinh
                                phenoData = pd, 
                                featureData = fd, 
                                expressionFamily=gaussianff(), 
                                lowerDetectionLimit=1)
human_exprs_new_cds <- reduceDimension(human_exprs_new_cds, reduction_method = 'L1-span-tree', norm_method = 'none', verbose = T, pseudo_expr = 0)
human_exprs_new_cds <- orderCells(human_exprs_new_cds)

plot_cell_trajectory(human_exprs_new_cds)


