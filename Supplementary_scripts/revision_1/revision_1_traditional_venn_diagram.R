library(monocle)
library(piano)
library(UpSetR)
library(xacHelper)
library(venneuler)

revision_1_fig_dir <- "../Figures/First_revision/"

# this script is used to analyze minor analysis asked by the reviewers: 
########################################################################################################################################################
# 1. make a traditional venn diagram 
########################################################################################################################################################
load('../RData/fig1.RData')

listInput <- list('High dispersion' = HSMM_overdispersion_genes$gene_id, Slicer = row.names(HSMM_myo)[HSMM_slicer_genes], 
                  PCA = HSMM_pca_res$gene_vec, 'DP genes' = HSMM_myo_ordering_genes, 'Time DEG' = HSMM_deg_ordering_genes)

upset(fromList(listInput), nsets = 6, order.by = "freq") + nm_theme()
length(HSMM_slicer_genes)

length(intersect(as.character(HSMM_overdispersion_genes$gene_id), row.names(subset(HSMM_BEAM_std_res, qval < 0.1)))) / nrow(subset(HSMM_BEAM_res, qval < 0.1))

element_all <- c(as.character(HSMM_overdispersion_genes$gene_id), 
                 row.names(HSMM_myo)[HSMM_slicer_genes],
                 HSMM_pca_res$gene_vec,
                 HSMM_myo_ordering_genes, 
                 HSMM_deg_ordering_genes
                 )
sets_all <- c(rep(paste('High dispersion', sep = ''), length(HSMM_overdispersion_genes$gene_id)), 
              rep(paste('Slicer', sep = ''), length(HSMM_slicer_genes)),
              rep(paste('PCA', sep = ''), length(HSMM_pca_res$gene_vec)),
              rep(paste('DP genes'), length(HSMM_myo_ordering_genes)),
              rep(paste('Time DEG'), length(HSMM_deg_ordering_genes)))
              
pdf(paste(revision_1_fig_dir, "conventional_venn_diagram.pdf", sep = ''))
venneuler_venn(element_all, sets_all)
dev.off()

########################################################################################################################################################
# 2. split panel B into two subplots: 
########################################################################################################################################################

pdf(paste(main_fig_dir, "fig1b.1.pdf", sep = ''), height = 1, width = 2)
plot_cell_clusters(HSMM_myo, color_by = 'Hours', cell_size = 1) + #geom_point(aes(color = factor(Hours))) + #shape = Cluster, 
  monocle:::monocle_theme_opts() +
  nm_theme() +
  scale_color_manual(values = HSMM_cols)
dev.off()

pdf(paste(main_fig_dir, "fig1b.2.pdf", sep = ''), height = 1, width = 2)
plot_cell_clusters(HSMM_myo, cell_size = 1) + geom_point(aes(color = factor(Cluster))) + #shape = , 
  monocle:::monocle_theme_opts() +
  nm_theme() #+ scale_color_manual(values = HSMM_cols)
dev.off()


########################################################################################################################################################
# 3. change the FDR threshold to see what is downstream effects 
########################################################################################################################################################
load('../RData/fig_si2.RData')
dim(subset(HSMM_BEAM_res, qval < 0.1)) #confirm the BEAM gene list 
dim(subset(HSMM_BEAM_res, qval < 0.05)) #confirm the BEAM gene list 
num_hsmm_beam_clusters <- 5

pdf(paste(SI_fig_dir, 'HSMM_BEAM_res_test.pdf', sep = ''), height=3.5, width=3.5) 
ph_res <- plot_genes_branched_heatmap(HSMM_myo[row.names(subset(HSMM_BEAM_res, qval < 0.1)),], 
                                      num_clusters = num_hsmm_beam_clusters, 
                                      branch_point=1,
                                      branch_labels = c("M1", "M2"),
                                      cores = detectCores(), 
                                      use_gene_short_name = T, 
                                      show_rownames = F,
                                      hmcols=hmcols,
                                      return_heatmap=TRUE)
dev.off()


State <- pData(HSMM_myo)$State
pData(HSMM_myo)$State[State == 3] <- 2
pData(HSMM_myo)$State[State == 2] <- 3

qval_vec <- c(0.05, 0.01, 0.1) 
for(qval_threshold in qval_vec) {
  pdf(paste(SI_fig_dir, "HSMM_BEAM_res_", qval_threshold, ".pdf", sep = ''), height=3.5, width=3.5) 
  ph_res <- plot_genes_branched_heatmap(HSMM_myo[row.names(subset(HSMM_BEAM_res, qval < qval_threshold)),], 
                                        num_clusters = num_hsmm_beam_clusters, 
                                        branch_point=1,
                                        branch_labels = c("M1", "M2"),
                                        branch_colors = c('#979797', '#7990C8', '#F05662'), 
                                        cores = detectCores(), 
                                        use_gene_short_name = T, 
                                        show_rownames = F,
                                        hmcols=hmcols,
                                        return_heatmap=TRUE)
  dev.off()
}
