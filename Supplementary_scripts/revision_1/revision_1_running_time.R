#load all the necessary libraries 
library(monocle)
library(plyr)
library(dplyr)
library(destiny)
library(xacHelper)
library(grid)

revision_1_fig_dir <- "../Figures/First_revision/"
################################################################################################################################################################################################
#demonstrate the computational time and memory requirement for monocle 2: 
################################################################################################################################################################################################
#load('./RData/fig_SI6.RData')
load('../RData/fig_SI6.RData')

MAP_ordering_genes <- row.names(subset(fData(valid_subset_GSE72857_cds), use_for_ordering == T)) #row.names(MAP_clustering_DEG_genes_subset)[order(MAP_clustering_DEG_genes_subset$qval)][1:1000] #1971
all_GSE72857_cds <- setOrderingFilter(all_GSE72857_cds, ordering_genes = MAP_ordering_genes)
current <- Sys.time()
all_GSE72857_cds <- reduceDimension(all_GSE72857_cds, norm_method = 'log', verbose = F)
all_GSE72857_cds <- orderCells(all_GSE72857_cds)
monocle2_running_time <- Sys.time() - current   # 9.030406

current <- Sys.time()
all_GSE72857_cds <- reduceDimension(all_GSE72857_cds, norm_method = 'log', verbose = F, maxIter = 5)
all_GSE72857_cds <- orderCells(all_GSE72857_cds)
monocle2_running_time_iter <- Sys.time() - current

current <- Sys.time()
data <- t(convert2DRData(all_GSE72857_cds[, ], norm_method = 'log')) #this is also used in monocle for preprocessing 
dm <- DiffusionMap(data)
dpt <- DPT(dm)
DPT_running_time <- Sys.time() - current # 9.740878

current <- Sys.time()
all_GSE72857_cds <- reduceDimension(all_GSE72857_cds, norm_method = 'log', verbose = F, max_component = 10)
all_GSE72857_cds <- orderCells(all_GSE72857_cds)
monocle2_running_time_10_dims <- Sys.time() - current   # 9.030406

#Show the convergence of the values: 
DDRTree_res <- DDRTree(t(data), maxIter = 20, ncenter = monocle:::cal_ncenter(ncol(all_GSE72857_cds)), verbose = T)

pdf(paste(SI_fig_dir, 'objective_iterations.pdf', sep = ''), width = 1, height = 1)
qplot(x = 1:length(DDRTree_res$objective_vals), y = DDRTree_res$objective_vals, log = 'y', size = I(0.5)) + xlab('Iteration') + ylab('Objective values') + nm_theme()
dev.off()

#show the running as a function cell numbers: 

monocle2_running_time_vec <- rep(0, 14)
i <- 1 
cell_num_vec <-  c(300, 500, 800, 1000, 1200, 1500, 1800, 2000, 2500, 3000, 4000, 5000, 6000, 8000)
for(cell_nums in cell_num_vec) {
  message('current cell nums is: ', cell_nums)
  set.seed(2017)
  subset_all_GSE72857_cds <- all_GSE72857_cds[, sample(1:ncol(all_GSE72857_cds), cell_nums)]
  current <- Sys.time()
  subset_all_GSE72857_cds <- reduceDimension(subset_all_GSE72857_cds, norm_method = 'log', verbose = F, max_components = 10)
  subset_all_GSE72857_cds <- orderCells(subset_all_GSE72857_cds)
  monocle2_running_time_vec[i] <- Sys.time() - current
  i <- i + 1
}

qplot(x = cell_num_vec, y = monocle2_running_time_vec)

dpt_running_time_vec <- rep(0, 14)
i <- 1 
for(cell_nums in cell_num_vec) {
  message('current cell nums is: ', cell_nums)
  set.seed(2017)
  subset_all_GSE72857_cds <- all_GSE72857_cds[, sample(1:ncol(all_GSE72857_cds), cell_nums)]
  current <- Sys.time()
  data <- t(convert2DRData(all_GSE72857_cds[, ], norm_method = 'log')) #this is also used in monocle for preprocessing 
  dm <- DiffusionMap(data)
  dpt <- DPT(dm)
  dpt_running_time_vec[i] <- Sys.time() - current
  i <- i + 1
}

qplot(x = cell_num_vec, y = dpt_running_time_vec)

################################################################################################################################################################################################
# save the data for running Wishbone: 
################################################################################################################################################################################################
for(cell_nums in cell_num_vec) {
  message('current cell nums is: ', cell_nums)
  set.seed(2017)
  subset_all_GSE72857_cds <- all_GSE72857_cds[, sample(1:ncol(all_GSE72857_cds), cell_nums)]
  data <- t(convert2DRData(cds = subset_all_GSE72857_cds, norm_method = 'log'))
  write.csv(file = paste('../csv_data/Wishbone_test_data/running_time_test_', cell_nums, ".txt", sep = ''), data, quote = F, row.names = T)
}

data <- t(convert2DRData(cds = all_GSE72857_cds, norm_method = 'log'))
write.csv(file = paste('../csv_data/Wishbone_test_data/running_time_test_', 8365, ".txt", sep = ''), data, quote = F, row.names = T)

running_time <- data.frame('Monocle 2' = c(monocle2_running_time_vec, monocle2_running_time),
                           DPT = c(dpt_running_time_vec, DPT_running_time), 
                           Wishbone = c(0.08918226957321167, 0.16816035111745198, 0.20987714926401774, 0.39565401474634804, 0.3235554655392965, 0.5412967165311178, 0.7563043316205342, 0.6350314656893412, 1.1284896334012349, 1.6187012314796447, 1.977494215965271, 1.9608156005541484, 2.8247681498527526, 2.6412346839904783, 3.1584882299105326), 
                           cell_num = c(cell_num_vec, ncol(all_GSE72857_cds)))
mlt_running_time <- melt(running_time, id.var = 'cell_num')

pdf(paste(SI_fig_dir, 'running_time.pdf', sep = ''), width = 1, height = 1)
qplot(cell_num, value, color = variable, data = mlt_running_time, size = I(0.5)) + geom_smooth(size = I(0.5), method = 'rlm') + nm_theme() + xlab('Cell number') + ylab('Time (min)') + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pdf(paste(SI_fig_dir, 'running_time_helper.pdf', sep = ''))
qplot(cell_num, value, color = variable, data = mlt_running_time, size = I(0.5)) + geom_smooth(size = I(0.5)) 
dev.off()

qplot(cell_num, value, color = variable, data = mlt_running_time[30:1, ], size = I(0.5)) + geom_smooth(size = I(0.5)) 

# memory: 
all_GSE72857_cds <- reduceDimension(all_GSE72857_cds, norm_method = 'log', verbose = F)
all_GSE72857_cds <- orderCells(all_GSE72857_cds)

data <- t(convert2DRData(all_GSE72857_cds[, ], norm_method = 'log')) #this is also used in monocle for preprocessing 
dm <- DiffusionMap(data)
dpt <- DPT(dm)

# memory_df <- data.frame(Method = c('Monocle 2', 'DPT', 'Wishbone'), Memory = c())
memory_df <- data.frame(Method = factor(c('reduceDimension', 'orderCells', 'DiffusionMap', 'DPT', 'run_pca', 'run_tsne', 'run_diffusion_map', 'Wishbone', 'run_wishbone'), 
                                        levels = c('reduceDimension', 'orderCells', 'DiffusionMap', 'DPT', 'run_pca', 'run_tsne', 'run_diffusion_map', 'Wishbone', 'run_wishbone')), 
                        Memory = c(12001, 36446, 13328.5, 30007, 756.645, 702.047, 403.520, 403.520, 277.312), 
                        Type = factor(c('Monocle 2', 'Monocle 2', 'DPT', 'DPT', 'Wishbone', 'Wishbone', 'Wishbone', 'Wishbone', 'Wishbone'), levels = c('Monocle 2', 'DPT', 'Wishbone')))
pdf(paste(SI_fig_dir, 'running_memory.pdf', sep = ''), width = 2, height = 1.2)
ggplot(aes(Method, Memory), data = memory_df[, ]) + geom_bar(aes(fill = Type), stat = 'identity') + ylab('Memory \n (Mb)') + xlab('') + nm_theme() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

# convergence of Monocle 2: 
qplot(1:iteration, converge_score)

################################################################################################################################################################################################
# learn the structure in high dimension to see whether or not we can get branches associates with other cell types: 
###############################################################################################################################################################################################
plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'cell_type', cell_size = 0.5, show_branch_points = F) + nm_theme() + 
  scale_color_manual(values = Mar_seq_cols)#, cell_size = pData(valid_subset_GSE72857_cds)$Pseudotime

pc_variance_res <- plot_pc_variance_explained(valid_subset_GSE72857_cds, return_all = T) #check how many PCs we should use to recover the structure 

valid_subset_GSE72857_cds <- setOrderingFilter(valid_subset_GSE72857_cds, ordering_genes = row.names(data.info.genes))
valid_subset_GSE72857_cds <- reduceDimension(valid_subset_GSE72857_cds, norm_method = 'log', verbose = T, max_components = 8) 
valid_subset_GSE72857_cds <- orderCells(valid_subset_GSE72857_cds)

valid_subset_GSE72857_cds <- orderCells(valid_subset_GSE72857_cds, root_state=10)
#use the plotting function I made before to show the structure of the data: 
Cell_custom_color_scale <- c("CMP"="#760F23",
                             "DC" = "#E91C26",
                             "erythroid" = "#F49223",
                             "GMP" = "#FAED34")

pdf(paste(main_fig_dir, 'marseq_complex_tree.pdf', sep = ''), width = 2.5, height = 3)
plot_complex_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'as.factor(cluster)', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~cell_type, nrow = 2) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
  theme (legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank()) 
dev.off()

pdf(paste(main_fig_dir, 'marseq_complex_tree2.pdf', sep = ''), width = 2.5, height = 3)
plot_complex_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'as.factor(cluster)', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001) + facet_wrap(~cell_type, nrow = 2) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
  theme (legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank()) 
dev.off()

plot_complex_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'as.factor(cluster)')

# show the cell types: 
plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'Pseudotime') + facet_wrap(~State)

pdf(paste(revision_1_fig_dir, "monocle_marseq_high_dimension12.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'cluster', x = 1, y = 2, cell_size =  I(0.5)) + nm_theme() + xlab('Component 1') + ylab('Component 2')
dev.off()

pdf(paste(revision_1_fig_dir, "monocle_marseq_high_dimension13.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'cluster', x = 1, y = 3, cell_size =  I(0.5)) + nm_theme() + xlab('Component 1') + ylab('Component 3')
dev.off()

pdf(paste(revision_1_fig_dir, "monocle_marseq_high_dimension23.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'cluster', x = 2, y = 3, cell_size =  I(0.5)) + nm_theme() + xlab('Component 2') + ylab('Component 3')
dev.off()

pdf(paste(revision_1_fig_dir, "monocle_marseq_high_dimension24.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'cluster', x = 2, y = 4, cell_size =  I(0.5)) + nm_theme() + xlab('Component 2') + ylab('Component 4')
dev.off()

#plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'cluster', x = 5, y = 8)

#make the heatmap: 
state_cluster_stat <- table(pData(valid_subset_GSE72857_cds)[, c('State', 'cluster')])

state_cluster_stat <- apply(state_cluster_stat, 2, function(x) x / sum(x))
state_cluster_stat_ordered <- t(state_cluster_stat)

annotation_row = data.frame(`Lineage` = c("1" = 'Ery', "2" = 'Ery', "3" = 'Ery', "4" = 'Ery', "5" = 'Ery', "6" = 'Ery', 
                                          "7" = 'MP/EP', "8" = 'MK', "9" = 'GMP', "10" = 'GMP',
                                          "11" = 'DC', 
                                          "12" = 'B', "13" = 'B', "14" = 'M', "15" = 'M', "16" = 'N', "17" = 'N', "18" = 'E', 
                                          "19" = 'lymphoid')) #, 
                            # `cell type` = c("1" = 'erythroid', "2" = 'erythroid', "3" = 'erythroid', "4" = 'erythroid', "5" = 'erythroid', "6" = 'erythroid', 
                            #                 "7" = 'CMP', "8" = 'CMP', "9" = 'CMP', "10" = 'CMP',
                            #                 "11" = 'DC', 
                            #                 "12" = 'GMP', "13" = 'GMP', "14" = 'GMP', "15" = 'GMP', "16" = 'GMP', "17" = 'GMP', "18" = 'GMP', 
                            #                 "19" = 'lymphoid')) 
annotation_colors = c("GMP" = "#EF5B5B ", "CMP" = "#FFBA49 ", "erythroid" = "#8ABF69 ", "DC" = "#0FA3B1 ")

pdf(paste(revision_1_fig_dir, "monocle_marseq_state_branch_distribution.pdf", sep = ''), height = 6, width = 6)
pheatmap::pheatmap(state_cluster_stat_ordered, cluster_cols = F, cluster_rows = F, color = colorRampPalette(RColorBrewer::brewer.pal(n=9, name='Greys'))(10), 
                   annotation_row = annotation_row) + scale_color_manual(values = Mar_seq_cols)
dev.off()


pdf(paste(revision_1_fig_dir, "monocle_marseq_state_branch_distribution.pdf", sep = ''), height = 6, width = 6)
pheatmap::pheatmap(state_cluster_stat_ordered, cluster_cols = F, cluster_rows = F, color = colorRampPalette(RColorBrewer::brewer.pal(n=9, name='Greys'))(10), 
                   annotation_row = annotation_row) + scale_color_manual(values = Mar_seq_cols)
dev.off()

#add the contour plot: 
#color node by the fraction of cells from a particular type of the dataset 
#or just a graph for the states and the hiearchay

plot_ly(x = valid_subset_GSE72857_cds@reducedDimS[1, ], y = valid_subset_GSE72857_cds@reducedDimS[2, ], z = valid_subset_GSE72857_cds@reducedDimS[3, ], type = "scatter3d", mode = "markers", 
        color = Cell_custom_color_scale[pData(valid_subset_GSE72857_cds)$cell_type]) 

##########################################################################################################################################################################
# compare with DPT's result using ARI 
##########################################################################################################################################################################
# show a few kinetic plots for marker genes: 
#make the heatmap when set max_component for reduceDimension as 10: 
valid_subset_GSE72857_cds2 <- reduceDimension(valid_subset_GSE72857_cds, norm_method = 'log', verbose = T, max_components = 10, ncenter = 200) 
valid_subset_GSE72857_cds2 <- orderCells(valid_subset_GSE72857_cds2)
plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'Pseudotime') + facet_wrap(~cell_type)
plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'cell_type', root_states = c(7, 10)) + facet_wrap(~State)
# root_cell <- which(pData(valid_subset_GSE72857_cds)[, 'Pseudotime'] == 0)
root_state <- pData(valid_subset_GSE72857_cds2)[root_cell, 'State']
valid_subset_GSE72857_cds2 <- orderCells(valid_subset_GSE72857_cds2, root_state = root_state)

state_cluster_stat <- table(pData(valid_subset_GSE72857_cds2)[, c('State', 'cluster')])

state_cluster_stat <- apply(state_cluster_stat, 2, function(x) x / sum(x))
state_cluster_stat_ordered <- t(state_cluster_stat)

pdf(paste(revision_1_fig_dir, "monocle_marseq_state_branch_distribution2.pdf", sep = ''), height = 6, width = 6)
pheatmap::pheatmap(state_cluster_stat_ordered, cluster_cols = F, cluster_rows = F, color = colorRampPalette(RColorBrewer::brewer.pal(n=9, name='Greys'))(10), 
                   annotation_row = annotation_row)
dev.off()

detailed_cell_type_color <- c("B" = "#E088B8", "DC" = "#46C7EF", "E" = "#EFAD1E", "Ery" = "#8CB3DF", "M" = "#53C0AD", "MP/EP" = "#4EB859", "GMP" = "#D097C4", "MK" = "#ACC436", "N" = "#F5918A")

pData(valid_subset_GSE72857_cds2)$cell_type2 <- revalue(as.character(pData(valid_subset_GSE72857_cds2)$cluster), 
                                                      c("1" = 'Ery', "2" = 'Ery', "3" = 'Ery', "4" = 'Ery', "5" = 'Ery', "6" = 'Ery', 
                                                        "7" = 'MP/EP', "8" = 'MK', "9" = 'GMP', "10" = 'GMP',
                                                        "11" = 'DC', 
                                                        "12" = 'B', "13" = 'B', "14" = 'M', "15" = 'M', "16" = 'N', "17" = 'N', "18" = 'E', 
                                                        "19" = 'lymphoid'))

pdf(paste(main_fig_dir, 'marseq_complex_tree.1.1.pdf', sep = ''), width = 1.2, height = 1.3)
plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'as.factor(cell_type2)', show_branch_points = T, 
                             cell_size = 0.5, cell_link_size = 0.3, root_states = c(7, 10)) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
  theme (legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank()) + scale_color_manual(values = detailed_cell_type_color) + theme_void() + theme (legend.position="none", legend.title=element_blank())
dev.off()

pdf(paste(main_fig_dir, 'marseq_complex_tree.1.pdf', sep = ''), width = 2.5, height = 3)
plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'as.factor(cluster)', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3, root_states = c(7, 10)) + facet_wrap(~cell_type, nrow = 2) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
  theme (legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank()) 
dev.off()

pdf(paste(main_fig_dir, 'marseq_complex_tree2.1.pdf', sep = ''), width = 2.5, height = 3)
plot_complex_cell_trajectory(valid_subset_GSE72857_cds2, color_by = 'as.factor(cluster)', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001) + facet_wrap(~cell_type, nrow = 2) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
  theme (legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank()) 
dev.off()

valid_subset_GSE72857_cds3 <- reduceDimension(valid_subset_GSE72857_cds, norm_method = 'log', verbose = T, max_components = 10, ncenter = 500) 
valid_subset_GSE72857_cds3 <- orderCells(valid_subset_GSE72857_cds3)
plot_complex_cell_trajectory(valid_subset_GSE72857_cds3, color_by = 'Pseudotime') + facet_wrap(~cell_type)
plot_complex_cell_trajectory(valid_subset_GSE72857_cds3, color_by = 'cell_type') + facet_wrap(~State)
root_cell <- "W37034"
root_state <- pData(valid_subset_GSE72857_cds3)[root_cell, 'State']
valid_subset_GSE72857_cds3 <- orderCells(valid_subset_GSE72857_cds3, root_state = root_state)

state_cluster_stat <- table(pData(valid_subset_GSE72857_cds2)[, c('State', 'cluster')])

state_cluster_stat <- apply(state_cluster_stat, 2, function(x) x / sum(x))
state_cluster_stat_ordered <- t(state_cluster_stat)

pdf(paste(revision_1_fig_dir, "monocle_marseq_state_branch_distribution3.pdf", sep = ''), height = 6, width = 6)
pheatmap::pheatmap(state_cluster_stat_ordered, cluster_cols = F, cluster_rows = F, color = colorRampPalette(RColorBrewer::brewer.pal(n=9, name='Greys'))(10), 
                   annotation_row = annotation_row)
dev.off()

pdf(paste(main_fig_dir, 'marseq_complex_tree.2.pdf', sep = ''), width = 2.5, height = 3)
plot_complex_cell_trajectory(valid_subset_GSE72857_cds3, color_by = 'as.factor(cluster)', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~cell_type, nrow = 2) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
  theme (legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank()) 
dev.off()

pdf(paste(main_fig_dir, 'marseq_complex_tree2.2.pdf', sep = ''), width = 2.5, height = 3)
plot_complex_cell_trajectory(valid_subset_GSE72857_cds3, color_by = 'as.factor(cluster)', show_branch_points = T, cell_link_size = 0.3, cell_size = 0.0001) + facet_wrap(~cell_type, nrow = 2) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
  theme (legend.position="top", legend.title=element_blank()) + theme (legend.position="none", legend.title=element_blank()) 
dev.off()

##########################################################################################################################################################################
# compare with DPT's result using ARI 
##########################################################################################################################################################################
#load dpt result from Fabian group's nice tutorial 
load('../RData/MARSseq_analysis_tutorial.RData')
dpt_state; 
dpt_time; 
cor(pData(valid_subset_GSE72857_cds2)$Pseudotime, dpt_time, method = 'kendall', use = 'pairwise.complete.obs') 
calClusteringMetrics(pData(valid_subset_GSE72857_cds2)$State, dpt_states) #show dpt's states 

##########################################################################################################################################################################
# distribution matrix betwee Monocle 2 and DPT
##########################################################################################################################################################################
calClusteringMetrics <- function(cl1, cl2){
  tab_1 <- table(cl1, cl2)
  randi <- flexclust::randIndex(tab_1)
  
  vi <- vi.dist(cl1, cl2, parts = FALSE, base = 2)
  
  #arandi
  arandi <- arandi(cl1, cl2, adjust = TRUE)
  
  #add all result into a data frame
  randIndex_df <- data.frame(randIndex = c(randi, vi, arandi),
                             Type = c("rand index", "variation of information", "adjusted rand index"))
  
  return(randIndex_df)
}
calClusteringMetrics(pData(valid_subset_GSE72857_cds2)$State, branching[cluster.id[, 1] != 19]) 
ARI_branches <- data.frame(ARI = c(0.5222408, 0.5923145,  0.7225239), Type = c('DPT (original) vs Monocle 2', 'DPT (original) vs Clusters', 'Monocle 2 vs Clusters')) #0.6566774, , 'DPT + Monocle 2'
qplot(ARI, geom = 'bar', data = ARI_branches, stat = 'identity')

pdf(paste(SI_fig_dir, 'MARS_seq_ARI_branch_cluster_cole.pdf', sep = ''), width = 2, height = 2)
ggplot(data = ARI_branches, aes(Type, ARI)) + geom_bar(stat = 'identity', aes(fill = Type)) + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab('')
dev.off()

confusion.matrix2 <- table(branching[cluster.id[, 1] != 19],pData(valid_subset_GSE72857_cds2)$State)
pdf(paste(SI_fig_dir, 'MARS_seq_DPT_monocle_results.pdf', sep = ''), width = 6, height = 4)
pheatmap(t(t(confusion.matrix2) /colSums(confusion.matrix2)*100), 
         border_color = 'grey90', annotation_legend = FALSE,
         breaks = c(0,1,5,seq(10,100, by=5)),
         color = colorRampPalette(c('white','black'))(25),         
         cellwidth = 20, cellheight = 20, annotation=annotation, 
         annotation_colors = list(cluster=color.pal),
         cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

##########################################################################################################################################################################
# run analysis on the blood dataset (not very useful here)
##########################################################################################################################################################################
save.image('../RData/revision_1_running_time.RData')

