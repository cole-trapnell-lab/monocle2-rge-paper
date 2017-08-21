library(mcclust)
library(pheatmap)

revision_1_fig_dir <- "./Figures/First_revision/"
# compare between branch and pseudotime assignment consistency between different software
load('./RData/marseq_downsamling_empirical_ordering.RData')

pdata <- pData(valid_subset_GSE72857_cds)
all_software_pseudotime_df <- data.frame(DPT = pdata$dpt_pseudotime,
                                         Wishbone = pdata$wishbone_pseudotime, 
                                         Monocle2 = pdata$Pseudotime,
                                         Monocle1 = pdata$monocle1_pseudotime,
                                         SLICER = pdata$slicer_pseudotime)

all_software_state_df <- data.frame(DPT = pdata$dpt_state,
                                         Wishbone = pdata$wishbone_state, 
                                         Monocle2 = pdata$State,
                                         Monocle1 = pdata$monocle1_state,
                                         SLICER = pdata$slicer_state)

software_types <- colnames(all_software_pseudotime_df)
all_software_cor_mat <- matrix(rep(NA, length(software_types)^2), nrow = 5, dimnames = list(software_types, software_types))
all_software_kendall_tau_mat <- matrix(rep(NA, length(software_types)^2), nrow = 5, dimnames = list(software_types, software_types))
all_software_ari_mat <- matrix(rep(NA, length(software_types)^2), nrow = 5, dimnames = list(software_types, software_types))

for(software1 in software_types) {
  for(software2 in setdiff(software_types, software1)) {
    all_software_cor_mat[software1, software2] <- cor(all_software_pseudotime_df[, software1], all_software_pseudotime_df[, software2])
    all_software_kendall_tau_mat[software1, software2] <- cor(all_software_pseudotime_df[, software1], all_software_pseudotime_df[, software2], method = 'kendall', use = 'pairwise.complete.obs')
    all_software_ari_mat[software1, software2] <- calClusteringMetrics(all_software_state_df[, software1], all_software_state_df[, software2])$randIndex[3]
  }
}

pdf(paste(revision_1_fig_dir, 'all_software_cor_mat.pdf'), width = 6, height = 6)
pheatmap(all_software_cor_mat, cluster_rows = F, cluster_cols = F)
dev.off()

pdf(paste(revision_1_fig_dir, 'all_software_kendall_tau_mat.pdf'), width = 6, height = 6)
pheatmap(all_software_kendall_tau_mat, cluster_rows = F, cluster_cols = F)
dev.off()

pdf(paste(revision_1_fig_dir, 'all_software_ari_mat.pdf'), width = 6, height = 6)
pheatmap(all_software_ari_mat, cluster_rows = F, cluster_cols = F)
dev.off()


