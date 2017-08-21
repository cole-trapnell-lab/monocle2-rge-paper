rm(list = ls())

library(monocle)
library(mcclust)
library(xacHelper)
library(reshape2)
library(destiny)
revision_1_fig_dir <- "./Figures/First_revision/"

main_fig_dir <- "./Figures/main_figures/"
SI_fig_dir <- "./Figures/supplementary_figures/"

benchmark_type <- 'HSMM_myo'
##########################################################################################################################################
#run dpt for 1k genes vs all genes: 
##########################################################################################################################################
# HSMM dataset
##########################################################################################################################################
load('./RData/analysis_HSMM_data.RData')

dpt_HSMM_res <- run_new_dpt(HSMM_myo)
HSMM_myo_all_genes <- setOrderingFilter(HSMM_myo, row.names(fData(HSMM_myo)))
all_dpt_HSMM_res <- run_new_dpt(HSMM_myo_all_genes) #run all genes 

qplot(all_dpt_HSMM_res$pt, dpt_HSMM_res$pt)
qplot(all_dpt_HSMM_res$branch[, 1], dpt_HSMM_res$branch[, 2])

#calculate the pearson correlation and the ARI: 
t_1 <- all_dpt_HSMM_res$pt; t_2 <- dpt_HSMM_res$pt; 
HSMM_dpt_feature_cor_tau <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
HSMM_dpt_feature_cor_rho <- cor(t_1, t_2)

clusters_1 <- all_dpt_HSMM_res$branch[, 1]; clusters_2 <- dpt_HSMM_res$branch[, 1]
HSMM_dpt_feature_ARI <- calClusteringMetrics(as.character(clusters_1), as.character(clusters_2))

##########################################################################################################################################
# blood dataset
##########################################################################################################################################
# rm(list = ls())
# load('./RData/fig5.RData')

load('./RData/fig5.RData')

source('./scripts/function.R')
dpt_URMM_all_fig1b_res <- run_new_dpt(URMM_all_fig1b)
URMM_all_fig1b_all_fig1b_all_genes <- setOrderingFilter(URMM_all_fig1b, row.names(fData(URMM_all_fig1b)))
all_dpt_URMM_all_fig1b_res <- run_new_dpt(URMM_all_fig1b_all_fig1b_all_genes) #run all genes 

qplot(all_dpt_URMM_all_fig1b_res$pt, dpt_URMM_all_fig1b_res$pt)
qplot(all_dpt_URMM_all_fig1b_res$branch[, 1], dpt_URMM_all_fig1b_res$branch[, 2])

t_1 <- dpt_URMM_all_fig1b_res$pt; t_2 <- all_dpt_URMM_all_fig1b_res$pt; 
URMM_all_fig1b_dpt_feature_cor_tau <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
URMM_all_fig1b_dpt_feature_cor_rho <- cor(t_1, t_2)

clusters_1 <- dpt_URMM_all_fig1b_res$branch[, 1]; clusters_2 <- all_dpt_URMM_all_fig1b_res$branch[, 1]
URMM_all_fig1b_dpt_feature_ARI <- calClusteringMetrics(as.character(clusters_1), as.character(clusters_2))

##########################################################################################################################################
# lung dataset
##########################################################################################################################################
# rm(list = ls())
# load('./RData/fig_si1.RData')

load('./RData/fig_si1.RData')

source('./scripts/function.R')
absolute_cds <- setOrderingFilter(absolute_cds, quake_id)
dpt_lung_res <- run_new_dpt(absolute_cds)
lung_all_fig1b_all_genes <- setOrderingFilter(absolute_cds, row.names(fData(absolute_cds)))
all_dpt_lung_res <- run_new_dpt(lung_all_fig1b_all_genes) #run all genes 

qplot(all_dpt_lung_res$pt, dpt_lung_res$pt)
qplot(all_dpt_lung_res$branch[, 1], dpt_lung_res$branch[, 2])

t_1 <- dpt_lung_res$pt; t_2 <- all_dpt_lung_res$pt; 
lung_dpt_feature_cor_tau <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
lung_dpt_feature_cor_rho <- cor(t_1, t_2)

clusters_1 <- dpt_lung_res$branch[, 1]; clusters_2 <- all_dpt_lung_res$branch[, 1]
lung_dpt_feature_ARI <- calClusteringMetrics(as.character(clusters_1), as.character(clusters_2))

##########################################################################################################################################
# MAR-seq dataset
##########################################################################################################################################
#rm(list = ls())
# load('./RData/fig_si6.RData')

load('./RData/fig_si6.RData')
absolute_cds <- all_GSE72857_cds[, colnames(valid_subset_GSE72857_cds)]
absolute_cds <- setOrderingFilter(absolute_cds, MAP_ordering_genes)
dpt_MARS_seq_res <- run_new_dpt(absolute_cds[MAP_ordering_genes, ])
MARS_seq_all_genes <- setOrderingFilter(absolute_cds, row.names(fData(absolute_cds)))
MARS_seq_all_res <- run_new_dpt(MARS_seq_all_genes) #run all genes 

qplot(dpt_MARS_seq_res$pt, MARS_seq_all_res$pt)
qplot(dpt_MARS_seq_res$branch[, 1], MARS_seq_all_res$branch[, 2])

t_1 <- dpt_MARS_seq_res$pt; t_2 <- MARS_seq_all_res$pt; 
MARS_seq_dpt_feature_cor_tau <- cor(t_1, t_2, method = 'kendall', use = 'pairwise.complete.obs')
MARS_seq_dpt_feature_cor_rho <- cor(t_1, t_2)

clusters_1 <- dpt_MARS_seq_res$branch[, 1]; clusters_2 <- MARS_seq_all_res$branch[, 1]
MARS_seq_dpt_feature_ARI <- calClusteringMetrics(as.character(clusters_1), as.character(clusters_2))

dpFeature_dpt_res_df <- data.frame("Pearson correlation" = c(HSMM_dpt_feature_cor_rho, URMM_all_fig1b_dpt_feature_cor_rho, NA, MARS_seq_dpt_feature_cor_tau), 
                                   "Kendall's tau" = c(HSMM_dpt_feature_cor_tau, URMM_all_fig1b_dpt_feature_cor_tau, NA, MARS_seq_dpt_feature_cor_rho),
                                   ARI = c(HSMM_dpt_feature_ARI$randIndex[3], URMM_all_fig1b_dpt_feature_ARI$randIndex[3], NA, MARS_seq_dpt_feature_ARI$randIndex[3]),
                                   data = c("HSMM", 'Olsson', 'lung', 'MARS-seq'))
dpFeature_dpt_res_df_mlt <- melt(dpFeature_dpt_res_df)

pdf(paste(revision_1_fig_dir, "dpt_different_feature.pdf", sep = ''), width = 4, height = 2)
ggplot(aes(x=data, y = value), data = dpFeature_dpt_res_df_mlt) + geom_bar(aes(fill = variable), stat = "identity", position = 'dodge') + 
  facet_grid(~variable) + ylim(c(min(dpFeature_dpt_res_df$ARI, na.rm = T), 1)) + ylab('') + nm_theme() + xlab('') +  
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) #theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dev.off()

##########################################################################################################################################
# save data
##########################################################################################################################################
save(file = './RData/dpFeature_dpt_res.RData', dpFeature_dpt_res_df, dpFeature_dpt_res_df_mlt)

