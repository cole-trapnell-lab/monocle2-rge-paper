####################################################################################################################################################################################
# run simplePPT on all biological dataset used in this study: 
####################################################################################################################################################################################
#run simplePPT algorithm on the simulated data: 
#run matlab code and get the best lambda based on the gap statistics
#principal tree method
# bandwidth: 0.0050
# lambda: 0.0036
params.lambda <- 1
params.bandwidth <- 30

#run HSMM dataset after using diffusion map dimension reduction:
load('./RData/fig1.RData')
params.lambda <- 1
params.bandwidth <- 30
dpt_res <- run_new_dpt(HSMM_myo[, ], normalize = F)
qplot(dpt_res$dm@eigenvectors[, 1], dpt_res$dm@eigenvectors[, 2])
HSMM_seq_pt_res <- principal_tree(t(dpt_res$dm@eigenvectors[, 1:2]), MU = NULL, lambda = params.lambda, bandwidth = params.bandwidth, maxIter = 10, verbose = T)
pdf(paste(SI_fig_dir, 'simplePPT_HSMM.pdf', sep = ''), height = 1, width = 1)
qplot(HSMM_seq_pt_res$MU[1, ], HSMM_seq_pt_res$MU[2, ], color = pData(HSMM_myo)$Time) + nm_theme() + xlab('DM1') + ylab('DM2')
dev.off()

absolute_cds <- setOrderingFilter(absolute_cds, quake_id)
dpt_res <- run_new_dpt(absolute_cds[, ], normalize = F)
qplot(dpt_res$dm@eigenvectors[, 1], dpt_res$dm@eigenvectors[, 2])
lung_seq_pt_res <- principal_tree(t(dpt_res$dm@eigenvectors[, 1:2]), MU = NULL, lambda = params.lambda, bandwidth = params.bandwidth, maxIter = 10, verbose = T)
pdf(paste(SI_fig_dir, 'simplePPT_lung.pdf', sep = ''), height = 1, width = 1)
qplot(lung_seq_pt_res$MU[1, ], lung_seq_pt_res$MU[2, ], color = pData(absolute_cds)$Time) + nm_theme() + xlab('DM1') + ylab('DM2')
dev.off()

rm(list = ls())
#run blood dataset after using diffusion map dimension reduction:
load('./RData/fig4.RData')
params.lambda <- 1
params.bandwidth <- 30
dpt_res <- run_new_dpt(URMM_all_fig1b[, ], normalize = F)
qplot(dpt_res$dm@eigenvectors[, 1], dpt_res$dm@eigenvectors[, 2], color = pData(URMM_all_fig1b)$cluster)
blood_seq_pt_res <- principal_tree(t(dpt_res$dm@eigenvectors[, 1:2]), MU = NULL, lambda = params.lambda, bandwidth = params.bandwidth, maxIter = 10, verbose = T)
pdf(paste(SI_fig_dir, 'simplePPT_URMM.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(blood_seq_pt_res$MU[1, ], blood_seq_pt_res$MU[2, ], color = pData(URMM_all_fig1b)$cluster) + nm_theme() + xlab('DM1') + ylab('DM2')
dev.off()

rm(list = ls())
#run mar-seq dataset after using diffusion map dimension reduction:
load('./RData/analysis_score_ordering_MAR_seq.RData')
dpt_res <- run_new_dpt(valid_subset_GSE72857_cds[, ], normalize = F)
mar_seq_pt_res <- principal_tree(t(dpt_res$dm@eigenvectors[, 1:2]), MU = NULL, lambda = params.lambda, bandwidth = params.bandwidth, maxIter = 10, verbose = T)

pdf(paste(SI_fig_dir, 'simplePPT_marseq.pdf', sep = ''), height = 1, width = 1)
qplot(mar_seq_pt_res$MU[1, ], mar_seq_pt_res$MU[2, ], color = pData(valid_subset_GSE72857_cds)$cell_type) + nm_theme() + xlab('DM1') + ylab('DM2')
dev.off()

pdf(paste(SI_fig_dir, 'simplePPT_marseq.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(valid_subset_GSE72857_cds2, color_by = "cell_type") + nm_theme() + xlab('DM1') + ylab('DM2')
dev.off()

rm(list = ls())

HSMM_myo2 <- reduceDimension(HSMM_myo[, ], norm_method = 'log', reduction_method = 'simplePPT')
HSMM_myo2 <- orderCells(HSMM_myo2)
plot_cell_trajectory(HSMM_myo2)

absolute_cds <- setOrderingFilter(absolute_cds, quake_id)
absolute_cds2 <- reduceDimension(absolute_cds[, ], norm_method = 'log', reduction_method = 'simplePPT')
absolute_cds2 <- orderCells(absolute_cds2)
plot_cell_trajectory(absolute_cds2)

absolute_cds <- setOrderingFilter(absolute_cds, quake_id)
valid_subset_GSE72857_cds2 <- reduceDimension(valid_subset_GSE72857_cds[, ], norm_method = 'log', reduction_method = 'simplePPT', verbose = T)
valid_subset_GSE72857_cds2 <- orderCells(valid_subset_GSE72857_cds2)
plot_cell_trajectory(valid_subset_GSE72857_cds2)

##################################################################################################################################################################
# create animations: 
##################################################################################################################################################################
# fraction_wishbone_res <- read.table(paste("/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/wishbone/", wishbone_res_fraction, sep = ''), header = T, row.names=NULL)
# 
# p <- ggplot(fraction_wishbone_res, aes(dm1, dm2, size = trajectory, color = branch, frame = run)) +
#   geom_point() +
#   scale_x_log10()

DDRTree_R <- function(X,  dimensions = 2,
                      maxIter = 20,
                      sigma = 1e-3,
                      lambda = NULL,
                      ncenter = NULL,
                      param.gamma = 10,
                      tol = 1e-3,
                      verbose = F) {
  
  D <- nrow(X)
  N <- ncol(X)
  
  #initialization
  W <- pca_projection_R(X %*% t(X), dimensions)
  Z <- t(W) %*% X
  
  if(is.null(ncenter)) {
    K <- N
    Y <- Z[, 1:K]
  }
  else {
    message("running k-means clustering")
    K <- ncenter
    kmean_res <- kmeans(t(Z), K)
    Y <- kmean_res$centers
    Y <- t(Y)
  }
  if (is.null(lambda)){
    lambda = 5 * ncol(X)
  }
  
  #main loop:
  objs <- c()
  history <- list()
  for(iter in 1:maxIter) {
    
    # 		#Kruskal method to find optimal B (use RBGL algorithm: http://stackoverflow.com/questions/16605825/minimum-spaning-tree-with-kruskal-algorithm)
    distsqMU <- sqdist_R(Y, Y)
    # 		#convert with graph packagege to BAM class of graph an calculate mst
    # 		mstKruskalBAM <- mstree.kruskal(graphBAM(as.data.frame(distsqMU)))
    # 		#build new data frame with resut
    # 		stree <- data.frame(cbind(t(mstKruskalBAM$edgeList),
    # 		                                 t(mstKruskalBAM$weight)))
    
    ##########################use mst from igraph: ##########################
    g <- graph.adjacency(distsqMU, mode = 'lower', diag = T, weighted = T)
    g_mst <- mst(g)
    stree <- get.adjacency(g_mst, attr = 'weight', type = 'lower')
    stree_ori <- stree
    
    #convert to matrix:
    stree <- as.matrix(stree)
    stree <- stree + t(stree)
    B_tmp <- stree != 0
    B <- B_tmp
    B[B_tmp == FALSE] <- 0
    B[B_tmp == TRUE] <- 1
    L <- diag(colSums(B)) - B
    
    # #convert back to igraph package
    # stree <- graph.data.frame(mstKruskalDF, directed=FALSE)
    
    #compute R usingmean-shift update rule
    distZY <- sqdist_R(Z, Y)
    min_dist <- matrix(rep(apply(distZY, 1, min), times = K), ncol = K, byrow = F)
    tmp_distZY <- distZY - min_dist
    tmp_R <- exp(-tmp_distZY / sigma)
    #print(tmp_R)
    R <- tmp_R / matrix(rep(rowSums(tmp_R), times = K), byrow = F, ncol = K)
    #print(R)
    Gamma_mat <- matrix(rep(0, ncol(R) ^ 2), nrow = ncol(R))
    diag(Gamma_mat) <- colSums(R)
    
    #termination condition
    obj1 <- - sigma * sum(log(rowSums(exp(-tmp_distZY / sigma)))
                          - min_dist[, 1] /sigma)
    objs[iter] <- (base::norm(X - W %*% Z, '2'))^2 + lambda * sum(diag(Y %*% L %*% t(Y))) + param.gamma * obj1 #sum(diag(A))
    
    if(verbose)
      message('iter = ', iter, ' ', objs[iter])
    
    history$W[[iter]] <- W
    history$Z[[iter]] <- Z
    history$Y[[iter]] <- Y
    history$stree[[iter]] <- stree
    history$R[[iter]] <- R
    
    if(iter > 1) {
      if(abs(objs[iter] - objs[iter - 1]) / abs(objs[iter - 1]) < tol) {
        break
      }
      
    }
    
    #compute low dimension projection matrix
    tmp <- t(solve((((param.gamma + 1) / param.gamma) * ((lambda / param.gamma) * L + Gamma_mat) - t(R) %*% R), t(R)))
    Q <- 1 / (param.gamma + 1) * (diag(1, N) + tmp %*% t(R))
    C <- X %*% Q
    tmp1 <- C %*% t(X)
    W <- pca_projection_R((tmp1 + t(tmp1)) / 2, dimensions)
    Z <- t(W) %*% C
    Y <- t(solve((lambda / param.gamma * L + Gamma_mat), t(Z %*% R)))
    #print (Y)
  }
  
  history$objs <- objs
  
  return(list(W = W, Z = Z, stree = stree_ori, Y = Y, history = history))
}


library(gganimate)
library(ggplot2)
test <- DDRTree_R(log(as.matrix(exprs(absolute_cds[quake_id, ])) + 1), verbose = T)
Y_df_combine <- do.call(rbind, lapply(test$history$Y, t))
row.names(Y_df_combine) <- 1:nrow(Y_df_combine)
Y_df_combine <- as.data.frame(Y_df_combine)
Y_df_combine$iteration <- rep(1:19, each = ncol(absolute_cds))
Y_df_combine$Type <- 'principal graph'

Z_df_combine <- do.call(rbind, lapply(test$history$Z, t))
row.names(Z_df_combine) <- 1:nrow(Z_df_combine)
Z_df_combine <- as.data.frame(Z_df_combine)
Z_df_combine$iteration <- rep(1:19, each = ncol(absolute_cds))
Z_df_combine$Type <- paste('reduced space (', rep(pData(absolute_cds)$Time, 19), ')', sep = '')

YZ_df <- rbind(Y_df_combine, Z_df_combine)
YZ_df$res <- factor(rep(c('principal graph', 'reduced space'), each = nrow(Y_df_combine)), levels = c('reduced space', 'principal graph'))

theme_set(theme_bw())

p <- ggplot(YZ_df, aes(V1, V2, size = 0.5, color = Type, frame = iteration)) + #
  geom_point() + nm_theme() + xlab('Component 1') + ylab('Component 2') + 
  scale_color_manual(values = c("reduced space (E14.5)" = "#bdc131", "reduced space (E16.5)" = "#afa9d3", "reduced space (E18.5)" = "#35c5f4", "reduced space (Adult)" = "#d593c0", 'principal graph' = '#020202'))
gganimate(p, "DDRTree_movie.gif")

p <- ggplot(Y_df_combine, aes(V1, V2, frame = iteration)) +
  geom_point()

gganimate(p)

ica_HSMM_myo <- reduceDimension(HSMM_myo, reduction_method = 'simplePPT', initial_method = ICA, verbose = T)
ica_HSMM_myo <- orderCells(ica_HSMM_myo)
plot_cell_trajectory(ica_HSMM_myo)

