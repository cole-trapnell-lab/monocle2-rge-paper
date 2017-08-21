library(monocle)
library(GEOquery)
library(stringr)
library(xacHelper)
library(destiny)

# run dpt but build the kNN graph and visualize with forced layout: 
run_dpt <- function(data, branching = T, norm_method = 'log', root = NULL, verbose = F){
  if(verbose)
    message('root should be the id to the cell not the cell name ....')
  
  data_ori <- t(data)
  duplicated_cells <- duplicated(data_ori)
  
  # avoid dpt to remove duplicated cells
  data_ori[duplicated_cells, 1] <-  rnorm(sum(duplicated_cells), sd = min(data_ori[data_ori > 0]))
  dm <- DiffusionMap(as.matrix(data_ori))
  
  return(dm@eigenvectors)
}

# run the motor neuron dataset with DDRTree or L1 graph 
revision_1_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/Figures/First_revision/"

# get the annotation data
GSE94883_parse <- getGEO(filename = '/Users/xqiu/Downloads/GSE94883_series_matrix.txt.gz')

# get the expression matrix 
root_dir <- './csv_data/scTDA//' 
motor_neuron_files <- dir(root_dir)

valid_gene_vec <- read.table('./csv_data/scTDA/all_gene_list', stringsAsFactors = F)

umi_matrix <- NULL
cell_num <- 0
for(file in motor_neuron_files) {
  tmp <- read.table(gzfile(paste(root_dir, file, sep = '')), header = F, sep = '\t')
  sample <- str_split_fixed(file, "_", 3)[, 2] 
  colnames(tmp) <- paste(sample, (cell_num + 1):(cell_num + ncol(tmp)), sep = '_c')
  message('processing file: ', file)
  
  tmp <- tmp[tmp[, 1] %in% valid_gene_vec$V1, ]
  row.names(tmp) <- tmp[, 1]
  tmp <- tmp[, -1]
  print(dim(tmp))
  if(is.null(umi_matrix)) {
    umi_matrix <- tmp
  }  
  else {
    sample <- str_split_fixed(file, "_", 3)[, 2] 
    tmp <- tmp[row.names(umi_matrix), ]
    umi_matrix <- cbind(umi_matrix, tmp)
  }
  print(dim(umi_matrix))
  cell_num <- ncol(umi_matrix)
}

# create cds and run Monocle 2 (DDRTree)
pData(GSE94883_parse)[1, c('characteristics_ch1', 'source_name_ch1')]

pData(GSE94883_parse)$title_trim <- str_split_fixed(pData(GSE94883_parse)$title, "_", 2)[, 1]
table(pData(GSE94883_parse)[, c('title_trim', 'characteristics_ch1')])

batch_time <- table(pData(GSE94883_parse)[, c('title_trim', 'characteristics_ch1')])
group <- colnames(batch_time)[apply(batch_time, 1, which.max)]
names(group) <- row.names(batch_time)

cell_batch <- str_split_fixed(colnames(umi_matrix), '_', 2)[, 1]

pd <- data.frame(Cell = group[cell_batch], row.names = colnames(umi_matrix))
fd <- data.frame(gene = row.names(umi_matrix), row.names = row.names(umi_matrix))

motor_neuron_umi_cds <-  newCellDataSet(as.matrix(umi_matrix), 
                               phenoData = new("AnnotatedDataFrame", data = pd), 
                               featureData = new("AnnotatedDataFrame", data = fd), 
                               expressionFamily=negbinomial.size(), 
                               lowerDetectionLimit=1)
motor_neuron_umi_cds <- estimateSizeFactors(motor_neuron_umi_cds)
motor_neuron_umi_cds <- estimateDispersions(motor_neuron_umi_cds)

plot_pc_variance_explained(motor_neuron_umi_cds)

fData(motor_neuron_umi_cds)$use_for_ordering <- T

motor_neuron_umi_cds <- reduceDimension(motor_neuron_umi_cds, max_components=2, norm_method = 'log', num_dim = 6, 
                            reduction_method = 'tSNE', verbose = T)
motor_neuron_umi_cds <- clusterCells(motor_neuron_umi_cds, verbose = F)

plot_cell_clusters(motor_neuron_umi_cds, color_by = 'as.factor(Cluster)')
# plot_cell_clusters(motor_neuron_umi_cds, color_by = 'as.factor(Hours)')
plot_rho_delta(motor_neuron_umi_cds, rho_threshold = 20, delta_threshold = 4)

motor_neuron_umi_cds <- clusterCells(motor_neuron_umi_cds,  
                         rho_threshold = 20, 
                         delta_threshold = 4, 
                         skip_rho_sigma = T, 
                         verbose = F)

plot_cell_clusters(motor_neuron_umi_cds, color_by = 'as.factor(Cluster)')
#plot_cell_clusters(motor_neuron_umi_cds, color_by = 'as.factor(Hours)')

clustering_DEG_genes <- differentialGeneTest(motor_neuron_umi_cds[ ,], 
                                             fullModelFormulaStr = '~Cluster', 
                                             cores = detectCores())

motor_neuron_umi_cds_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000] 
motor_neuron_umi_cds <- setOrderingFilter(motor_neuron_umi_cds, ordering_genes = motor_neuron_umi_cds_ordering_genes)
motor_neuron_umi_cds <- reduceDimension(motor_neuron_umi_cds, norm_method = 'log', max_components = 6, verbose = T)
motor_neuron_umi_cds <- orderCells(motor_neuron_umi_cds)
# motor_neuron_umi_cds <- orderCells(motor_neuron_umi_cds, root_state=GM_state(motor_neuron_umi_cds))
plot_cell_trajectory(motor_neuron_umi_cds, color_by="Cluster")
plot_complex_cell_trajectory(motor_neuron_umi_cds, color_by = 'Cell')
plot_cell_trajectory(motor_neuron_umi_cds, color_by="Cell") + facet_wrap(~Cell)

motor_neuron_umi_cds <- reduceDimension(motor_neuron_umi_cds, norm_method = 'log', max_components = 6, verbose = T, initial_method = run_dpt, reduction_method = 'L1-graph')
motor_neuron_umi_cds <- orderCells(motor_neuron_umi_cds)

Root_state <- function(cds){
  T0_counts <- table(pData(cds)$State, pData(cds)$Cell)[, 1]
  as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))])
}

motor_neuron_umi_cds <- orderCells(motor_neuron_umi_cds, root_state = Root_state(motor_neuron_umi_cds))
plot_cell_trajectory(motor_neuron_umi_cds, color_by="Cell") + facet_wrap(~Cell)
plot_complex_cell_trajectory(motor_neuron_umi_cds, color_by = 'Cell')

pdf(paste(revision_1_fig_dir, 'scTDA_ddrtree_tree.pdf', sep = ''), width = 1.2, height = 1.2) 
plot_complex_cell_trajectory(motor_neuron_umi_cds, color_by = 'Cell', cell_size = 0.5) + nm_theme()
dev.off()

# calculate the Pseudotime correspondence between our method and their method 



# create cds and run L1 graph in high dimension 
motor_neuron_umi_cds_dpt <- run_dpt(exprs(motor_neuron_umi_cds[fData(motor_neuron_umi_cds)$use_for_ordering, ]))

X <- t(motor_neuron_umi_cds_dpt[, 1:6])
D <- nrow(X); N <- ncol(X)
Z <- X
C0 <- Z
Nz <- ncol(C0)

G <- get_knn(C0, 30)

# choose appropriate lambda, gamma and sigma

maxiter <- 20
eps <- 1.0000e-05
gamma <- 1 #smoothness
sigma <- 0.01000 #
lambda <- 50 #L1 g

load(filename);

ncenter = 500;

kmean_res <- kmeans(t(X), centers = ncenter)

C0 <- t(kmean_res$centers)

G <- get_knn(C0, K = 15)

scTDA_pg_res <- principal_graph(t(motor_neuron_umi_cds_dpt[, 1:6]), C0, G$G, maxiter = maxiter,
                                  eps = eps, gstruct = 'l1-graph',
                                  lambda = lambda, gamma = gamma,
                                  sigma = sigma, nn = 5, verbose = T)

print(qplot(scTDA_pg_res$C[1, ], scTDA_pg_res$C[2, ], color = 'red') + geom_point(aes(X[1, ], X[2, ]), color = 'black') + ggtitle('L1-graph'))
print(qplot(scTDA_pg_res$C[1, ], scTDA_pg_res$C[2, ], color = 'red') + geom_point(aes(X[1, ], X[2, ]), color = 'black') + ggtitle('L1-graph'))

# visualize as a graph: 
L1_graph <- igraph::graph_from_adjacency_matrix(scTDA_pg_res$W, mode = "undirected", weighted = T)

pData(motor_neuron_umi_cds)$cell_day <- as.integer(revalue(pData(motor_neuron_umi_cds)$Cell, c("sampling day: Day 2" = 2, "sampling day: Day 3" = 3, 
                                                                   "sampling day: Day 4" = 4, "sampling day: Day 5" = 5,
                                                                   "sampling day: Day 2" = 6)))

pData(motor_neuron_umi_cds)$cell_day[1:80] <- 6
apply(table(apply(scTDA_pg_res$P, 1, which.max)))
V(L1_graph)$color <- color[pData(motor_neuron_umi_cds)$Cell]

L1_layout_coord = layout_with_drl(L1_graph)

pData(motor_neuron_umi_cds)$max_cluster <- apply(P_res, 2, which.max)
P_res <- t(scTDA_pg_res$P)
P_res <- apply(P_res, 2, function(x) {
  index <- x < quantile(x, probs = 0.95)
  x[index] <- 0
  x[x >= quantile(x, probs = 0.95)] <- 1
  x
})

P_res %*% pData(motor_neuron_umi_cds)$cell_day
color <- t(scTDA_pg_res$P) %*% pData(motor_neuron_umi_cds)$cell_day

pdf(paste(revision_1_fig_dir, 'Time_L1_graph_color_state.pdf', sep = '')) 
plot(L1_graph, layout = L1_layout_coord, vertex.size=2, vertex.label=NA, vertex.color = color[pData(motor_neuron_umi_cds)$Cell])
dev.off()

plot(L1_graph, layout = t(scTDA_pg_res$C)[, 1:2], vertex.size=4, vertex.label=NA, vertex.color = color)
ddply(pData(motor_neuron_umi_cds), .(max_cluster), summarize, mean = mean(cell_day))

##################################################################################################################################################################
dx <- FNN::get.knn(motor_neuron_umi_cds_dpt, k = 30)
nn.index <- dx$nn.index
nn.dist <- dx$nn.dist
N <- nrow(nn.index)

color <- c("sampling day: Day 2" =  "blue", "sampling day: Day 3" =  "pink", "sampling day: Day 4" =  "red", "sampling day: Day 5" =  "yellow", "sampling day: Day 6" =  "green")

knn_graph <- NULL
edges <- reshape2::melt(t(nn.index)); colnames(edges) <- c("B", "A", "C"); edges <- edges[,c("A","B","C")]
edges_weight <- reshape2::melt(t(nn.dist)); #colnames(edges_weight) = c("B", "A", "C"); edges_weight = edges_weight[,c("A","B","C")]
edges$B <- edges$C

# if(use_dist)
#   edges$C <- 1 #edges_weight$value
# else
#   edges$C <- 1

#Remove repetitions
edges = unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))

Adj <- Matrix::sparseMatrix(i = c(edges$A, edges$B), j = c(edges$B, edges$A), x = c(edges$C, edges$C), dims = c(N, N))

knn_graph <- igraph::graph_from_adjacency_matrix(Adj, mode = "undirected", weighted = T)

V(knn_graph)$color <- color[pData(motor_neuron_umi_cds)$Cell]

layout_coord = layout_with_drl(knn_graph)

pdf(paste(revision_1_fig_dir, 'Time_knn_graph_color_state.pdf', sep = '')) 
plot(knn_graph, layout = layout_coord, vertex.size=2, vertex.label=NA, vertex.color = color[pData(motor_neuron_umi_cds)$Cell])
dev.off()

pdf(paste(revision_1_fig_dir, 'DDRTree_state_knn_graph_color_state.pdf', sep = '')) 
plot(knn_graph, layout = layout_coord, vertex.size=2, vertex.label=NA, vertex.color = pData(motor_neuron_umi_cds)$State)
dev.off()

################################################################################################################################################################################
# save the data 
################################################################################################################################################################################
save.image('/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision//RData/reversion_1_analysis_motor_neuron_scTDA.RData')

################################################################################################################################################################################
# save the data 
################################################################################################################################################################################

gene_exprs <- exprs(motor_neuron_umi_cds)
dim(gene_exprs)
fData <- fData(motor_neuron_umi_cds)
pData <- pData(motor_neuron_umi_cds)
save(file = 'gene_exprs.txt', gene_exprs)
save(file = 'pData.txt', pData)
save(file = 'fData.txt', fData)

L1graph_res <- readMat('/Users/xqiu/Dropbox (Personal)/Projects/SPL+L1graph/tree-l1-code/mouse_high_c_result.mat')

L1graph <- igraph::graph_from_adjacency_matrix(L1graph_res$W, mode = "undirected", weighted = T)
qplot(L1graph_res$C[1, ],  L1graph_res$C[2, ])

#V(knn_graph)$color <- color[pData(motor_neuron_umi_cds)$Cell]

layout_coord = L1graph_res$centers

plot(L1graph, layout = layout_coord[, 1:2], vertex.size=2, vertex.label=NA)

plot(L1graph, layout = layout.drl(L1graph), vertex.size=2, vertex.label=NA)

qplot(c(L1graph_res$C0[1, ], L1graph_res$centers[, 1]), c(L1graph_res$C0[2, ], L1graph_res$centers[, 2]), color = c(rep('red', ncol(L1graph_res$C0)), c(rep('black', nrow(L1graph_res$centers)))) )

qplot(L1graph_res$centers[, 1],  L1graph_res$centers[, 2])

Cell <- pData(motor_neuron_umi_cds)$Cell
qplot(motor_neuron_umi_cds_dpt[, 2], motor_neuron_umi_cds_dpt[, 3], color = Cell, facets = ~ as.character(Cell))

################################################################################################################################################################################
# L1 graph on forced directed graph coordinates
################################################################################################################################################################################
dpt_knn_layout_coord = layout_with_drl(knn_graph, dim = 2)

X <- t(layout_coord )
D <- nrow(X); N <- ncol(X)
Z <- X
C0 <- Z
Nz <- ncol(C0)

G <- get_knn(C0, 30)

# choose appropriate lambda, gamma and sigma

maxiter <- 20
eps <- 1.0000e-05
gamma <- 1 #smoothness
sigma <- 1 #
lambda <- 65 #L1 g

# load(filename);

ncenter = 300;

kmean_res <- kmeans(t(X), centers = ncenter)

C0 <- t(kmean_res$centers)

G <- get_knn(C0, K = 15)

scTDA_pg_res <- principal_graph(t(layout_coord[, 1:2]), C0, G$G, maxiter = maxiter,
                                eps = eps, gstruct = 'l1-graph',
                                lambda = lambda, gamma = gamma,
                                sigma = sigma, nn = 15, verbose = T)

print(qplot(scTDA_pg_res$C[1, ], scTDA_pg_res$C[2, ], color = 'red') + geom_point(aes(X[1, ], X[2, ]), color = 'black') + ggtitle('L1-graph'))
print(qplot(scTDA_pg_res$C[1, ], scTDA_pg_res$C[2, ], color = 'red') + geom_point(aes(X[1, ], X[2, ]), color = 'black') + ggtitle('L1-graph'))

# visualize as a graph: 
L1_graph <- igraph::graph_from_adjacency_matrix(scTDA_pg_res$W, mode = "undirected", weighted = T)

pData(motor_neuron_umi_cds)$cell_day <- as.integer(revalue(pData(motor_neuron_umi_cds)$Cell, c("sampling day: Day 2" = 2, "sampling day: Day 3" = 3, 
                                                                                               "sampling day: Day 4" = 4, "sampling day: Day 5" = 5,
                                                                                               "sampling day: Day 2" = 6)))

pData(motor_neuron_umi_cds)$cell_day[1:80] <- 6
pData(motor_neuron_umi_cds)$kmean_cluster <- kmean_res$cluster
avg_day <- ddply(pData(motor_neuron_umi_cds), .(kmean_cluster), summarize, mean = mean(cell_day))

V(L1_graph)$color <- avg_day[, 2]

L1_layout_coord = layout_with_drl(L1_graph)

pdf(paste(revision_1_fig_dir, 'Time_L1_graph_color_state.pdf', sep = '')) 
plot(L1_graph, layout = t(scTDA_pg_res$C), vertex.size=4, vertex.label=NA)
dev.off()

qplot(scTDA_pg_res$C[1, ], scTDA_pg_res$C[2, ], color = scTDA_pg_res$C[, 1])

########################################################################################################################################################
# show the result from scTDA software
########################################################################################################################################################

test_res <- read_graph('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/test.graphml', format = 'graphml')
gl_res <- read_graph('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/scTDA_gl.graphml', format = 'graphml')
g_res <- read_graph('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/scTDA_g.graphml', format = 'graphml')
gl_res <- read_graph('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/scTDA_gl.gml', format = 'gml')

node_time_csv <- read.csv('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/csv_data/scTDA/node_time_csv', sep = '\t', row.names = 1, header = F)
node_coord <- read.csv('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/csv_data/scTDA/node_coords', sep = '\t', row.names = 1, header = F)

node_time_csv <- read.csv('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/csv_data/scTDA/node_time_csv2', sep = '\t', stringsAsFactors = F, header = F)
V(res)$color <- node_time_csv[, 1] * 10
plot(res, layout = layout.drl(res), vertex.size=4, vertex.label=NA, vertex.color = node_time_csv[, 1])

res <- delete_vertices(gl_res, which(degree(g_res) == 0))
res <- delete_vertices(g_res, which(degree(g_res) == 0))

V(res)$name <- 1:408
plot(res, layout = as.matrix(node_coord[as.character(sort(as.numeric(row.names(node_time_csv)))), ]), vertex.size=4, vertex.label=NA, vertex.color = node_time_csv[, 1])

scTDA_layout <- layout.drl(res)
qplot(scTDA_layout[, 1], scTDA_layout[, 2], color = node_time_csv[as.character(sort(as.numeric(row.names(node_time_csv)))), 1]) + scale_color_gradient2()
qplot(scTDA_layout[, 1], scTDA_layout[, 2], color = node_time_csv[, 1])

qplot(scTDA_layout[, 1], scTDA_layout[, 2], color =  node_time_csv[, 1])

scTDA_layout <- layout.drl(res)

#order(as.numeric(row.names(node_coord))) 
qplot(node_coord[as.character(sort(as.numeric(row.names(node_coord)))), 1], node_coord[as.character(sort(as.numeric(row.names(node_coord)))), 2], 
      color =  node_time_csv[, 1]) + scale_color_gradient2()

qplot(node_coord[, 1], node_coord[, 2], color =  node_time_csv[, 1]) + scale_color_gradient2()

exp_rng <- range(node_time_csv[, 1])
bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, length.out = 10)
if(is.null(hmcols)) {
  hmcols <- blue2green2red(length(bks) - 1)
}

test_post <- readMat('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/test_post.mat')
dic <- readMat('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/dic.mat')
posgl <- readMat('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/posgl.mat')
pel <- readMat('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/pel.mat')
dr <- readMat('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/dr.mat')

posgl_coord <- do.call(rbind.data.frame, posgl)

posg <- readMat('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/posg.mat')
posg_coord <- do.call(rbind.data.frame, posgl)

num_cells <- unlist(lapply(dic, length)) / 2
qplot(posgl_coord[, 1], posgl_coord[, 2], color = node_time_csv[row.names(posgl_coord), 1])

valid_num_cells <-  num_cells[row.names(posgl_coord)] #num_cells[!names(num_cells) %in% as.character(which(degree(g_res) == 0))]

V(res)$name <- row.names(posgl_coord)
plot(res, layout = as.matrix(posgl_coord), vertex.size = valid_num_cells, vertex.label = NA)
plot(res, layout = as.matrix(posg_coord), vertex.size = valid_num_cells, vertex.label = NA)

plot(res, layout = layout.drl(res), vertex.size = valid_num_cells, vertex.label = NA)

qplot(posgl_coord[, 1], posgl_coord[, 2], color = pel$pel[1, ])
qplot(posgl_coord[, 1], posgl_coord[, 2], color = dr$dr[1, ])

################################################################################################################################################################ 
# calculate the average time 
################################################################################################################################################################ 
table_tpm <- read.table('/Users/xqiu/Downloads/SCTDA_Tutorial/data/table_tpm.tsv', sep = '\t', header = T, row.names = 1)

timepoint <- table_tpm$timepoint
avg_time <- unlist(lapply(dic, function(x) mean(timepoint[x])))
qplot(posgl_coord[, 1], posgl_coord[, 2], color = avg_time[row.names(posgl_coord)])

qplot(scTDA_layout[, 1], scTDA_layout[, 2], color = avg_time[row.names(posgl_coord)])

################################################################################################################################################################ 
# run DDRTree to get the tree structure: 
################################################################################################################################################################ 

cell_batch <- str_split_fixed(colnames(umi_matrix), '_', 2)[, 1]

pd <- data.frame(timepoint = table_tpm$timepoint, lib = table_tpm$lib, row.names = row.names(table_tpm))
fd <- data.frame(gene = colnames(table_tpm)[-c(1:2)], row.names = colnames(table_tpm)[-c(1:2)])

scTDA_exp_1_cds <-  newCellDataSet(t(as.matrix(table_tpm[, -c(1:2)])), 
                                        phenoData = new("AnnotatedDataFrame", data = pd), 
                                        featureData = new("AnnotatedDataFrame", data = fd), 
                                        expressionFamily=negbinomial.size(), 
                                        lowerDetectionLimit=1)
scTDA_exp_1_cds <- estimateSizeFactors(scTDA_exp_1_cds)
scTDA_exp_1_cds <- estimateDispersions(scTDA_exp_1_cds)

scTDA_exp_1_cds <- detectGenes(scTDA_exp_1_cds, min_expr = 0.1)

fData(scTDA_exp_1_cds)$use_for_ordering <- fData(scTDA_exp_1_cds)$num_cells_expressed > round(ncol(scTDA_exp_1_cds) / 50)
plot_pc_variance_explained(scTDA_exp_1_cds)

scTDA_exp_1_cds <- reduceDimension(scTDA_exp_1_cds, max_components=2, norm_method = 'log', num_dim = 6, 
                                        reduction_method = 'tSNE', verbose = T)
scTDA_exp_1_cds <- clusterCells(scTDA_exp_1_cds, verbose = F)

plot_cell_clusters(scTDA_exp_1_cds, color_by = 'as.factor(Cluster)')
# plot_cell_clusters(scTDA_exp_1_cds, color_by = 'as.factor(Hours)')
plot_rho_delta(scTDA_exp_1_cds, rho_threshold = 8, delta_threshold = 8)

scTDA_exp_1_cds <- clusterCells(scTDA_exp_1_cds,  
                                     rho_threshold = 8, 
                                     delta_threshold = 8, 
                                     skip_rho_sigma = T, 
                                     verbose = F)

plot_cell_clusters(scTDA_exp_1_cds, color_by = 'as.factor(Cluster)')
#plot_cell_clusters(scTDA_exp_1_cds, color_by = 'as.factor(Hours)')

clustering_DEG_genes <- differentialGeneTest(scTDA_exp_1_cds[ ,], 
                                             fullModelFormulaStr = '~Cluster', 
                                             cores = detectCores())

scTDA_exp_1_cds_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000] 
scTDA_exp_1_cds <- setOrderingFilter(scTDA_exp_1_cds, ordering_genes = scTDA_exp_1_cds_ordering_genes)
scTDA_exp_1_cds <- reduceDimension(scTDA_exp_1_cds, norm_method = 'log', max_components = 6, verbose = T)
scTDA_exp_1_cds <- orderCells(scTDA_exp_1_cds)
# scTDA_exp_1_cds <- orderCells(scTDA_exp_1_cds, root_state=GM_state(scTDA_exp_1_cds))
plot_cell_trajectory(scTDA_exp_1_cds, color_by="Cluster")
plot_cell_trajectory(scTDA_exp_1_cds, color_by="timepoint")
plot_complex_cell_trajectory(scTDA_exp_1_cds, color_by = 'timepoint')
plot_cell_trajectory(scTDA_exp_1_cds, color_by="Cell") + facet_wrap(~Cell)

scTDA_exp_1_cds <- reduceDimension(scTDA_exp_1_cds, norm_method = 'log', max_components = 6, verbose = T, initial_method = run_dpt, reduction_method = 'L1-graph')
scTDA_exp_1_cds <- orderCells(scTDA_exp_1_cds)

Root_state <- function(cds){
  T0_counts <- table(pData(cds)$State, pData(cds)[, 1])[, 1]
  as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))])
}

scTDA_exp_1_cds <- orderCells(scTDA_exp_1_cds, root_state = Root_state(scTDA_exp_1_cds))
plot_cell_trajectory(scTDA_exp_1_cds, color_by="timepoint") + facet_wrap(~timepoint)
plot_complex_cell_trajectory(scTDA_exp_1_cds, color_by = 'timepoint')

pdf(paste(revision_1_fig_dir, 'exp_1_scTDA_ddrtree_trajectory.pdf', sep = ''), width = 1.2, height = 1.2) 
plot_cell_trajectory(scTDA_exp_1_cds, color_by = 'as.character(timepoint)', cell_size = 0.5) + nm_theme()
dev.off()

pdf(paste(revision_1_fig_dir, 'exp_1_scTDA_ddrtree_tree.pdf', sep = ''), width = 1.2, height = 1.2) 
plot_complex_cell_trajectory(scTDA_exp_1_cds, color_by = 'Cell', cell_size = 0.5) + nm_theme()
dev.off()

##################################################################################################################################################################
scTDA_exp_1_cds_dpt <- run_dpt(exprs(scTDA_exp_1_cds[fData(scTDA_exp_1_cds)$use_for_ordering, ]))
dx <- FNN::get.knn(scTDA_exp_1_cds_dpt, k = 15)
nn.index <- dx$nn.index
nn.dist <- dx$nn.dist
N <- nrow(nn.index)

color <- c("sampling day: Day 2" =  "blue", "sampling day: Day 3" =  "pink", "sampling day: Day 4" =  "red", "sampling day: Day 5" =  "yellow", "sampling day: Day 6" =  "green")

knn_graph <- NULL
edges <- reshape2::melt(t(nn.index)); colnames(edges) <- c("B", "A", "C"); edges <- edges[,c("A","B","C")]
edges_weight <- reshape2::melt(t(nn.dist)); #colnames(edges_weight) = c("B", "A", "C"); edges_weight = edges_weight[,c("A","B","C")]
edges$B <- edges$C

# if(use_dist)
#   edges$C <- 1 #edges_weight$value
# else
#   edges$C <- 1

#Remove repetitions
edges = unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))

Adj <- Matrix::sparseMatrix(i = c(edges$A, edges$B), j = c(edges$B, edges$A), x = c(edges$C, edges$C), dims = c(N, N))

exp1_knn_graph <- igraph::graph_from_adjacency_matrix(Adj, mode = "undirected", weighted = T)

# V(exp1_knn_graph)$color <- color[pData(motor_neuron_umi_cds)$Cell]

layout_coord = layout_with_drl(exp1_knn_graph)

timepoint_color <- c("2" =  "blue", "3" =  "pink", "4" =  "red", "5" =  "yellow", "6" =  "green")

pdf(paste(revision_1_fig_dir, 'exp_1_Time_L1_graph_color_state.pdf', sep = '')) 
plot(exp1_knn_graph, layout = layout_coord, vertex.size = 2, vertex.color = timepoint_color[as.character(pData(scTDA_exp_1_cds)$timepoint)], vertex.label = NA)
dev.off()

qplot(pData(scTDA_exp_1_cds)$timepoint, pData(scTDA_exp_1_cds)$Pseudotime)
qplot(pel$pel[1, ], pData(scTDA_exp_1_cds)$Pseudotime)

avg_pseudo_time <- unlist(lapply(dic, function(x) mean(pData(scTDA_exp_1_cds)$Pseudotime[x])))
qplot(avg_time[row.names(posg_coord)], pel$pel[1, ])

qplot(avg_pseudo_time[row.names(posg_coord)], dr$dr[1, ]) + xlab('Pseudotime (DDRTree)') + ylab('Pseudotime (scTDA)') # 0.902

##################################################################################################################################################################
# analyze the Pancreas dataset with scTDA using data from the scTDA tutorial: 
##################################################################################################################################################################
# This command generates three tab-separated files, named Embryo.all.tsv, Embryo.no_subsampling.tsv, and Embryo.mapper.tsv, 
# where the rows correspond to selected cells. The first column of Embryo.all.tsv and Embryo.no_subsampling.tsv contains a
# unique identifier of the cell, the second column contains the sampling timepoint, the third column contains the library ID,
# and the remaining columns contain  log2(1+TPM)log2(1+TPM)  expression values for all the genes. Embryo.all.tsv contains 
# subsampled data, whereas Embryo.no_subsampling.tsv constains non-subsampled data. Embryo.mapper.tsv contains  log2(1+TPM)
# log2(1+TPM)  expression values for the genes that were selected using the scTDA.Preprocess.select_genes().

Embryo.all <- read.table('/Users/xqiu/Downloads/scTDA Tutorial 2/Embryo.all.tsv', sep = '\t', header = T, row.names = 1)
Embryo.no_subsampling <- read.table('/Users/xqiu/Downloads/scTDA Tutorial 2/Embryo.no_subsampling.tsv', sep = '\t', header = T, row.names = 1)
Embryo.mapper <- read.table('/Users/xqiu/Downloads/scTDA Tutorial 2/Embryo.mapper.tsv', sep = '\t', header = T, row.names = 1)

# create CDS 
pd <- data.frame(timepoint = Embryo.all$timepoint, lib = Embryo.all$lib, row.names = row.names(Embryo.all))
fd <- data.frame(gene = colnames(Embryo.all)[-c(1:2)], row.names = colnames(Embryo.all)[-c(1:2)])

Embryo.all.matrix <- apply(Embryo.all[, -c(1:2)], 1, function(x) as.numeric(x))
row.names(Embryo.all.matrix) <- colnames(Embryo.all)[-c(1:2)]

Embryo.all_cds <-  newCellDataSet(Embryo.all.matrix, 
                                        phenoData = new("AnnotatedDataFrame", data = pd), 
                                        featureData = new("AnnotatedDataFrame", data = fd), 
                                        expressionFamily=gaussianff(), 
                                        lowerDetectionLimit=1)

Embryo.all_cds <- setOrderingFilter(Embryo.all_cds, ordering_genes = c('ABCG2' , 'ACAT2' , 'ACTL8' , 'ADM' , 'AHNAK' , 'ALPL' , 'ALPP' , 'ANP32E' , 'ANPEP', 
                                                     'ANXA6' , 'APOA1' , 'AQP3' , 'ARGFX' , 'ARL4D' , 'ASAH1' , 'ATG9A' , 'ATP1B1' , 'ATP6V0A4', 
                                                     'B3GNT2' , 'BAG4' , 'BIK' , 'BRDT' , 'CCKBR' , 'CCNA1' , 'CD24' , 'CD53' , 'CD81' , 'CD9', 
                                                     'CGA' , 'CLDN10' , 'CLDN4' , 'CLIC1' , 'CRIP1' , 'CUL2' , 'CYP2S1' , 'CYP51A1' , 'DDIT4', 
                                                     'DHCR24' , 'DHCR7' , 'DLG5' , 'DLX3' , 'DNMT1' , 'DPPA5' , 'DSC2' , 'DUXA' , 'EFNA1', 
                                                     'EGLN3' , 'EMP2' , 'EPCAM' , 'ERRFI1' , 'FABP3' , 'FAM151A' , 'FASN' , 'FBXL20' , 'FBXO5', 
                                                     'FDFT1' , 'FOLR1' , 'FOXR1' , 'FRAT2' , 'FSCN1' , 'FZD5' , 'GABARAPL1' , 'GATA2' , 'GATA3', 
                                                     'GATA6' , 'GCM1' , 'GDF9' , 'GINS3' , 'GLIPR2' , 'GPATCH3' , 'GPRC5A' , 'GPX2' , 'GRN', 
                                                     'H1F0' , 'HERPUD1' , 'HINT1' , 'HIPK3' , 'HJURP' , 'HMGCS1' , 'HMOX1' , 'IFI30' , 'IFITM1', 
                                                     'IFITM3' , 'ITM2C' , 'JUP' , 'KLF17' , 'KLF4' , 'KPNA7' , 'KRT18' , 'KRT19' , 'KRT8', 
                                                     'LCP1' , 'LGMN' , 'LMBR1L' , 'LPCAT3' , 'LRP2' , 'LYN' , 'MAGEA4' , 'MFN1' , 'MPP1', 
                                                     'MSMO1' , 'MTMR14' , 'MTRNR2L1' , 'MVD' , 'MX1' , 'NFKBIA' , 'NID1' , 'NLRP4' , 'OOEP', 
                                                     'OVOL1' , 'PDCL3' , 'PDGFA' , 'PDXK' , 'PEBP1' , 'PGF' , 'PHOSPHO1' , 'PIM1' , 'PLEC', 
                                                     'PLTP' , 'PNP' , 'PPP1R14A' , 'PPP1R15A' , 'PRAP1' , 'PRSS23' , 'PRSS8' , 'PSAT1', 
                                                     'PTGES' , 'PTN' , 'PYGB' , 'RAB25' , 'RABGEF1' , 'REEP1' , 'RGS16' , 'RGS2' , 'RNF11', 
                                                     'RNF168' , 'RPS4Y1' , 'RSRP1' , 'S100A13' , 'S100A14' , 'S100A16' , 'S100A6' , 'SEMA6A', 
                                                     'SEPHS2' , 'SERPINB9' , 'SERTAD1' , 'SERTAD3' , 'SGK1' , 'SH2D4A' , 'SHISA5' , 'SLC12A3', 
                                                     'SLC1A3' , 'SLC1A5' , 'SLC25A1' , 'SLC34A2' , 'SLC35F6' , 'SLC38A1' , 'SLC38A2', 
                                                     'SLC4A11' , 'SLC6A19' , 'SLC6A8' , 'SLC7A2' , 'SLC7A4' , 'SLC7A5' , 'SLC7A8' , 'SNAR-E', 
                                                     'SNHG15' , 'SOX15' , 'SPARC' , 'SPTLC2' , 'SQLE' , 'STS' , 'SUSD2' , 'TACSTD2' , 'TBCC', 
                                                     'TBPL1' , 'TCL1A' , 'TKTL1' , 'TMBIM1' , 'TMEM97' , 'TOPORS' , 'TPM4' , 'TPRX1' , 'TRIM60', 
                                                     'TUBB4A' , 'UACA' , 'VSIG10' , 'WDR45' , 'WDR47' , 'WEE2' , 'YIPF4' , 'YPEL2' , 'YPEL5', 
                                                     'ZFAND5' , 'ZNF280A' , 'ZNF296' , 'ZSCAN4'))

Embryo.all_cds <- reduceDimension(Embryo.all_cds, norm_method = 'none', max_components = 6, verbose = T)
Embryo.all_cds <- orderCells(Embryo.all_cds)
root_state <- function(cds){
  T0_counts <- table(pData(cds)$State, pData(cds)$timepoint)[,"3"]
  as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))])
}

# motor_neuron_umi_cds <- orderCells(motor_neuron_umi_cds, root_state=GM_state(motor_neuron_umi_cds))
plot_cell_trajectory(Embryo.all_cds, color_by = 'timepoint')
plot_cell_trajectory(Embryo.all_cds, color_by = 'lib')

# compare with the pseudotime from scTDA

# make the KNN graph 
Embryo.all_scTDA_exp_1_cds_dpt <- run_dpt(exprs(Embryo.all_cds[fData(Embryo.all_cds)$use_for_ordering, ]))
dx <- FNN::get.knn(Embryo.all_scTDA_exp_1_cds_dpt, k = 15)
nn.index <- dx$nn.index
nn.dist <- dx$nn.dist
N <- nrow(nn.index)

knn_graph <- NULL
edges <- reshape2::melt(t(nn.index)); colnames(edges) <- c("B", "A", "C"); edges <- edges[,c("A","B","C")]
edges_weight <- reshape2::melt(t(nn.dist)); #colnames(edges_weight) = c("B", "A", "C"); edges_weight = edges_weight[,c("A","B","C")]
edges$B <- edges$C

# if(use_dist)
#   edges$C <- 1 #edges_weight$value
# else
#   edges$C <- 1

#Remove repetitions
edges = unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))

Adj <- Matrix::sparseMatrix(i = c(edges$A, edges$B), j = c(edges$B, edges$A), x = c(edges$C, edges$C), dims = c(N, N))

Embryo.all_exp1_knn_graph <- igraph::graph_from_adjacency_matrix(Adj, mode = "undirected", weighted = T)

# V(exp1_knn_graph)$color <- color[pData(motor_neuron_umi_cds)$Cell]

Embryo.all_layout_coord = layout_with_drl(Embryo.all_exp1_knn_graph)

# timepoint_color <- c("2" =  "blue", "3" =  "pink", "4" =  "red", "5" =  "yellow", "6" =  "green")
timepoint_color <- c("3" =  "blue", "4" =  "pink", "5" =  "red", "6" =  "yellow", "7" =  "green")

pdf(paste(revision_1_fig_dir, 'Embryo.all_exp_1_Time_L1_graph_color_state.pdf', sep = '')) 
plot(Embryo.all_exp1_knn_graph, layout = Embryo.all_layout_coord, vertex.size = 4, vertex.color = timepoint_color[as.character(pData(Embryo.all_cds)$timepoint)], vertex.label = NA)
dev.off()

pdf(paste(revision_1_fig_dir, 'Embryo.all_exp_1_Time_L1_graph_color_pseudotime.pdf', sep = ''), width = 1, height = 1) 
qplot(Embryo.all_layout_coord[, 1], Embryo.all_layout_coord[, 2], color = pData(Embryo.all_cds)$Pseudotime, size = I(0.25)) + xlab('') + ylab('') + 
  scale_colour_gradient2() + nm_theme()# 0.902
dev.off()

pData(Embryo.all_cds)$Days <- as.character(pData(Embryo.all_cds)$timepoint)
pdf(paste(revision_1_fig_dir, 'Embryo.all_exp_1_Time_L1_graph_color_pseudotime.pdf', sep = ''), width = 1, height = 1) 
plot_cell_trajectory(Embryo.all_cds, color_by = 'Days', show_branch_points = F, cell_size = I(0.5)) + xlab('') + ylab('') +  nm_theme() + scale_color_manual(values = timepoint_color) # 0.902
dev.off()


# compare the KNN graph with their netwrk graph 


# trajectory reconstruction, KNN graph, pseudotime comparision

##################################################################################################################################################################
# save result: 
##################################################################################################################################################################
save.image('./RData/revision_1_analysis_motor_neuron_scTDA.RData')
