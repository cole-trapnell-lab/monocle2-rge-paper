library(devtools); load_all('/Users/xqiu/Dropbox (Personal)/Projects/monocle-dev');

maxiter <- 20
eps <- 1.0000e-05
gstruct <- 'l1-graph' #'span-tree'
gamma <- 1000 #smoothness
sigma <- 0.0001000 #
lambda <- 10 #L1 g
nn <- 10
verbose = T

human_cell_cycle_high_c_center_principal_graph <- principal_graph(t(human_cell_cycle_high_c_center_dpt@eigenvectors[, 1:2]), C0, G$G, maxiter = maxiter,
                                                                  eps = eps, gstruct = gstruct,
                                                                  lambda = lambda, gamma = gamma,
                                                                  sigma = sigma, nn = 5, verbose = T)

dm_df <- data.frame(x = human_cell_cycle_high_c_center_dpt@eigenvectors[, 1], y = human_cell_cycle_high_c_center_dpt@eigenvectors[, 2], cluster  = factor(as.character(baz$colour)))

qplot(human_cell_cycle_high_c_center_principal_graph$C[1, ], human_cell_cycle_high_c_center_principal_graph$C[2, ], color = I('black')) + 
  # geom_point(aes(human_cell_cycle_high_c_center_dpt@eigenvectors[, 1], human_cell_cycle_high_c_center_dpt@eigenvectors[, 2], 
  #                color = baz$colour)) + ggtitle('spanning-tree') + 
  geom_point(aes(kmenas_res$centers[, 1], kmenas_res$centers[, 2]), color = I('yellow'))

colnames(human_cell_cycle_high_c_center) <- paste('gene', 1:ncol(human_cell_cycle_high_c_center), sep = '_')
row.names(human_cell_cycle_high_c_center) <- paste('cell', 1:nrow(human_cell_cycle_high_c_center), sep = '_')
pd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = row.names(human_cell_cycle_high_c_center), row.names = row.names(human_cell_cycle_high_c_center)))
fd <- new("AnnotatedDataFrame", data = data.frame(cell = colnames(human_cell_cycle_high_c_center),
                                                  row.names = colnames(human_cell_cycle_high_c_center)))

return_dm <- function(data){
  res <- DiffusionMap(t(data))
  res@eigenvectors
}
human_exprs_new_cds <- newCellDataSet(as(as.matrix(t(human_cell_cycle_high_c_center)), "sparseMatrix"), #he cofactor parameter is used for arcsinh
                                      phenoData = pd, 
                                      featureData = fd, 
                                      expressionFamily=gaussianff(), 
                                      lowerDetectionLimit=1)
human_exprs_new_cds <- reduceDimension(human_exprs_new_cds, reduction_method = 'L1-graph', norm_method = 'none', scaling = F, verbose = T, pseudo_expr = 0,
                                       eps = eps, initial_method = return_dm, 
                                       C0 = C0, G = G$G, maxiter = maxiter,
                                       lambda = lambda, gamma = gamma,
                                       sigma = sigma, nn = 10)

cds <- human_exprs_new_cds 
l1_graph_res <- human_cell_cycle_high_c_center_principal_graph
colnames(l1_graph_res$C) <- paste("Y_", 1:ncol(l1_graph_res$C), sep = "")
DCs <- l1_graph_res$X
colnames(DCs) <- paste("Y_", 1:ncol(l1_graph_res$X), sep = "")
reducedDimW(cds) <- DCs
reducedDimS(cds) <- DCs #1_graph_res$X
reducedDimK(cds) <- l1_graph_res$C
cds@auxOrderingData[["DDRTree"]]$objective_vals <- tail(l1_graph_res$objs, 1)
cds@auxOrderingData[["DDRTree"]]$W <- l1_graph_res$W
cds@auxOrderingData[["DDRTree"]]$P <- l1_graph_res$P
adjusted_K <- Matrix::t(reducedDimK(cds))
dp <- as.matrix(dist(adjusted_K))
cellPairwiseDistances(cds) <- dp
gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
dp_mst <- minimum.spanning.tree(gp)
minSpanningTree(cds) <- dp_mst
cds@dim_reduce_type <- "DDRTree"
cds <- findNearestPointOnMST(cds)
plot_cell_trajectory(human_exprs_new_cds, show_tree = T)

human_exprs_new_cds <- orderCells(human_exprs_new_cds)

adjusted_K <- Matrix::t(reducedDimK(human_exprs_new_cds))
dp <- as.matrix(dist(adjusted_K))

cellPairwiseDistances(human_exprs_new_cds) <- dp
gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
dp_mst <- minimum.spanning.tree(gp)
minSpanningTree(human_exprs_new_cds) <- dp_mst
human_exprs_new_cds@dim_reduce_type <- "DDRTree"
human_exprs_new_cds <- findNearestPointOnMST(human_exprs_new_cds)

plot_cell_trajectory(human_exprs_new_cds, color_by = 'Size_Factor') 

save.image('./RData/cell_cycle_analysis_tmp.RData')
