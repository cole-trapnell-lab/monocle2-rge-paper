rm(list = ls())

source('./Scripts/fig1_run_monocle_vignette_cluster_cells.R', echo = T) #obtain the same clustering result as in vignette

HSMM_myo <- detectGenes(HSMM_myo, min_expr = 0.1) #detect genes for HSMM_myo data
pData(HSMM_myo)$Time <- pData(HSMM_myo)$Hour

################################################################################################################################################################################################################################################
## load necessary packages
################################################################################################################################################################################################################################################
library(stringr)
library(xacHelper)
library(grid)
library(devtools)
library(monocle)
library(plyr)
library(dplyr)
library(SLICER)
library(pheatmap)
library(UpSetR)
library(piano)
library(HSMMSingleCell)

#helper function used for making the marker gradient plot:
plot_cell_clusters_all <- function (cds, x = 1, y = 2, color_by = "Cluster", markers = NULL,
    show_cell_names = FALSE, cell_size = 1.5, cell_name_size = 2, return_all = T,
    ...)
{
    if (is.null(cds@reducedDimA) | length(pData(cds)$Cluster) ==
        0) {
        stop("Error: Clustering is not performed yet. Please call clusterCells() before calling this function.")
    }
    gene_short_name <- NULL
    sample_name <- NULL
    lib_info <- pData(cds)
    tSNE_dim_coords <- reducedDimA(cds)
    data_df <- data.frame(t(tSNE_dim_coords[c(x, y), ]))
    colnames(data_df) <- c("data_dim_1", "data_dim_2")
    data_df$sample_name <- colnames(cds)
    data_df <- merge(data_df, lib_info, by.x = "sample_name",
        by.y = "row.names")
    markers_exprs <- NULL
    if (is.null(markers) == FALSE) {
        markers_fData <- subset(fData(cds), gene_short_name %in%
            markers)
        if (nrow(markers_fData) >= 1) {
            markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData),
                ])))
            colnames(markers_exprs)[1:2] <- c("feature_id", "cell_id")
            markers_exprs <- merge(markers_exprs, markers_fData,
                by.x = "feature_id", by.y = "row.names")
            markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
            markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
        }
    }
    if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) >
        0) {
        data_df <- merge(data_df, markers_exprs, by.x = "sample_name",
            by.y = "cell_id")
        g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2,
            size = log10(value + 0.1))) + facet_wrap(~feature_label)
    }
    else {
        g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
    }
    if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) >
        0) {
        g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE)
    }
    else {
        g <- g + geom_point(aes_string(color = color_by), size = I(cell_size),
            na.rm = TRUE)
    }
    g <- g + monocle_theme_opts() + xlab(paste("Component", x)) +
        ylab(paste("Component", y)) + theme(legend.position = "top",
        legend.key.height = grid::unit(0.35, "in")) + theme(legend.key = element_blank()) +
        theme(panel.background = element_rect(fill = "white"))
    list(g = g, data_df = data_df)
}

################################################################################################################################################################################################################################################
## set up directories and load all GSC file for enrichment analysis
################################################################################################################################################################################################################################################
main_fig_dir <- "./Figures/main_figures/"
SI_fig_dir <- "./Figures/supplementary_figures/"

source('./scripts/function.R', echo = T)

#load the GSC file from human / mouse:
root_directory <- "./data/"

human_go_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Human/symbol/Human_GO_AllPathways_with_GO_iea_June_20_2014_symbol.gmt", sep=""), encoding="latin1")
names(human_go_gsc$gsc) <- str_split_fixed(names(human_go_gsc$gsc), "%", 2)[,1]

human_reactome_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Human/symbol/Pathways/Human_Reactome_June_20_2014_symbol.gmt", sep=""), encoding="latin1")
names(human_reactome_gsc$gsc) <- str_split_fixed(names(human_reactome_gsc$gsc), "%", 2)[,1]

mouse_go_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/GO/MOUSE_GO_bp_with_GO_iea_symbol.gmt", sep=""), encoding="latin1")
names(mouse_go_gsc$gsc) <- str_split_fixed(names(mouse_go_gsc$gsc), "%", 2)[,1]

mouse_reactome_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/Pathways/Mouse_Reactome_June_20_2014_symbol.gmt", sep=""), encoding="latin1")
names(mouse_reactome_gsc$gsc) <- str_split_fixed(names(mouse_reactome_gsc$gsc), "%", 2)[,1]

mouse_kegg_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/Pathways/Mouse_Human_KEGG_June_20_2014_symbol.gmt", sep=""), encoding="latin1")
names(mouse_kegg_gsc$gsc) <- str_split_fixed(names(mouse_kegg_gsc$gsc), "%", 2)[,1]

mouse_go_gsc_cc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/GO/MOUSE_GO_cc_with_GO_iea_symbol.gmt", sep=""), encoding="latin1")
names(mouse_go_gsc_cc$gsc) <- str_split_fixed(names(mouse_go_gsc_cc$gsc), "%", 2)[,1]
mouse_go_gsc_mf <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/GO/MOUSE_GO_mf_with_GO_iea_symbol.gmt", sep=""), encoding="latin1")
names(mouse_go_gsc_mf$gsc) <- str_split_fixed(names(mouse_go_gsc_mf$gsc), "%", 2)[,1]

mouse_reactome_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/Pathways/Mouse_Reactome_June_20_2014_symbol.gmt", sep=""), encoding="latin1")
names(mouse_reactome_gsc$gsc) <- str_split_fixed(names(mouse_reactome_gsc$gsc), "%", 2)[,1]

mouse_kegg_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/Pathways/Mouse_Human_KEGG_June_20_2014_symbol.gmt", sep=""), encoding="latin1")
names(mouse_kegg_gsc$gsc) <- str_split_fixed(names(mouse_kegg_gsc$gsc), "%", 2)[,1]

################################################################################################################################################################################################################################################
# generate supplementary figure SI4a (Note that we need to mannaully adjust the color gradient) 
# figure used in the paper is: HSMM_fibro_clusters_manual.pdf 
################################################################################################################################################################################################################################################

res_all <- plot_cell_clusters_all(HSMM, color="CellType", markers = c("MYF5", "MYOG", "MYOD1", "ANPEP"))

pdf(paste(SI_fig_dir, "SI4a_ori.pdf", sep = ''), height = 2, width = 2)
res_all$g + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

#use color gradient to show gene expression:
# load('./RData/data_df')
data_df <- res_all$data_df
limits = c(-1, 5)
data_df$value_process <- log10(data_df$value + 0.1)
data_df$value_process[data_df$value_process > limits[2]] <- 0.5

g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + facet_wrap(~feature_label)
g <- g + geom_point(aes_string(color = "value_process"), na.rm = TRUE, size = 0.5) + scale_colour_gradient(low = "grey", high = "red", name = "log10(Exp + 0.1)")

g <- g +
  #scale_color_brewer(palette="Set1") +
  monocle:::monocle_theme_opts() +
  xlab(paste("Component", 1)) +
  ylab(paste("Component", 2)) +
  theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
  #guides(color = guide_legend(label.position = "top")) +
  theme(legend.key = element_blank()) +
  theme(panel.background = element_rect(fill='white'))

pdf(paste(SI_fig_dir, 'SI4a.pdf', sep = ''), width=2.25, height=3.1)
g + nm_theme()
dev.off()

pdf(paste(SI_fig_dir, 'SI4a_helper.pdf', sep = ''))
g
dev.off()

################################################################################################################################################################################################################################################
# generate supplementary figure SI1b (two sub-panels: SI1b.1 is on the top and SI1b.2 is on the bottom)
################################################################################################################################################################################################################################################

c("#EF5B5B", "#FFBA49",  "#8ABF69", "#0FA3B1")
HSMM_cols <- c("0" = "#bdc131", "24" = "#afa9d3", "48" = "#35c5f4", "72" = "#d593c0")
HSMM_cols <- c("0" = "#EF5B5B", "24" = "#FFBA49", "48" = "#8ABF69", "72" = "#0FA3B1")

Hours <- pData(HSMM_myo)$Time
pdf(paste(main_fig_dir, "fig1b.pdf", sep = ''), height = 2, width = 2)
plot_cell_clusters(HSMM_myo, color_by = 'Hours', cell_size = 1) + geom_point(aes(shape = Cluster, color = factor(Hours))) +
  monocle:::monocle_theme_opts() +
  nm_theme() +
  scale_color_manual(values = HSMM_cols)
dev.off()

pdf(paste(main_fig_dir, "SI1b.1.pdf", sep = ''), height = 1, width = 2)
plot_cell_clusters(HSMM_myo, color_by = 'Hours', cell_size = 1) + #geom_point(aes(color = factor(Hours))) + #shape = Cluster, 
  monocle:::monocle_theme_opts() +
  nm_theme() +
  scale_color_manual(values = HSMM_cols)
dev.off()

pdf(paste(main_fig_dir, "SI1b.2.pdf", sep = ''), height = 1, width = 2)
plot_cell_clusters(HSMM_myo, cell_size = 1) + geom_point(aes(color = factor(Cluster))) + #shape = , 
  monocle:::monocle_theme_opts() +
  nm_theme() #+ scale_color_manual(values = HSMM_cols)
dev.off()

pdf(paste(SI_fig_dir, "SI1b_helper.pdf", sep = ''))
plot_cell_clusters(HSMM_myo, color_by = 'Hours') + geom_point(aes(shape = Cluster, color = factor(Hours)), size = 2) + scale_color_manual(values = HSMM_cols)
dev.off()

save.image('./RData/analysis_HSMM_data.RData')

################################################################################################################################################
## feature selection by marker genes (not used in the paper)
################################################################################################################################################
marker_gene <- c('CDK1', 'ID1', 'MYOG', 'MEF2C', 'MYH3')
marker_gene_ids <- row.names(subset(fData(HSMM_myo), gene_short_name %in% marker_gene))

#use all marker genes for making the trajectory
HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes = marker_gene_ids)
HSMM_myo <- reduceDimension(HSMM_myo, norm_method = 'log', verbose = T)

pdf(paste(SI_fig_dir, "HSMM_myo_marker.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(HSMM_myo, color_by = 'Time') + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

pdf(paste(SI_fig_dir, "HSMM_myo_marker_time_dp.pdf", sep = ''), height = 2, width = 2)
plot_cell_clusters(HSMM_myo, color_by = 'Time') + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

pdf(paste(SI_fig_dir, "HSMM_myo_marker_clusters.pdf", sep = ''), height = 2, width = 2)
plot_cell_clusters(HSMM_myo, color_by = 'Cluster') + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

################################################################################################################################################
## feature selection by PCA (not used in the paper)
########################################################################################################################################################
##HSMM data
HSMM_pca_res <- selectOrderingGenes2(HSMM_myo, dim = c(1, 2), top_gene_num = 1000, method = c('PCA'), verbose = T)

#use all pca gene for making the trajectory
HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes = HSMM_pca_res$gene_vec)
HSMM_myo <- reduceDimension(HSMM_myo, norm_method = 'log', verbose = T)

pdf(paste(SI_fig_dir, "HSMM_myo_pca.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(HSMM_myo, color_by = 'Time') + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

################################################################################################################################################
##  feature selection by SLICER (not used in the paper)
########################################################################################################################################################
HSMM_slicer_genes <- slicer_feature_genes(HSMM_myo)
HSMM_mat <- log(exprs(HSMM_myo[HSMM_slicer_genes, ]) + 1)
pheatmap(HSMM_mat, cluster_rows = T, cluster_cols = T)
#use all slicer gene for making the trajectory
HSMM_slicer <- setOrderingFilter(HSMM_myo, ordering_genes = row.names(HSMM_myo)[HSMM_slicer_genes])
HSMM_slicer <- reduceDimension(HSMM_slicer, norm_method = 'log', verbose = T)

pdf(paste(SI_fig_dir, "HSMM_slicer.pdf", sep = ''), height = 2, width = 2)
plot_cell_trajectory(HSMM_slicer, color_by = 'Time')  + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

########################################################################################################################################################
##  feature selection by Over-dispersion  (not used in the paper)
########################################################################################################################################################
#Over-dispersion genes: 
disp_table <- dispersionTable(HSMM_myo)
HSMM_overdispersion_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)
HSMM_myo <- setOrderingFilter(HSMM_myo, HSMM_overdispersion_genes$gene_id)
plot_ordering_genes(HSMM_myo)

########################################################################################################################################################
##  feature selection by differential gene expression test
########################################################################################################################################################  
HSMM_myo_time_dependent_genes <- differentialGeneTest(HSMM_myo[HSMM_expressed_genes,], fullModelFormulaStr="~Time", reducedModelFormulaStr="~1", cores=detectCores())
HSMM_deg_ordering_genes <- row.names(subset(HSMM_myo_time_dependent_genes, qval < 1e-2))
HSMM_deg_ordering_genes <- row.names(head(HSMM_myo_time_dependent_genes[order(HSMM_myo_time_dependent_genes$qval),], n=1000))

################################################################################################################################################################################################################################################
## Panel C
################################################################################################################################################################################################################################################
#use upset to show the overlapping for different feature selection methods: (pca, high-dispersion, time-dependent de, dp feature genes, slicer feature selection)
# example of list input (list of named vectors)
listInput <- list('High dispersion' = HSMM_overdispersion_genes$gene_id, Slicer = row.names(HSMM_myo)[HSMM_slicer_genes], 
                  PCA = HSMM_pca_res$gene_vec, 'DP genes' = HSMM_myo_ordering_genes, 'Time DEG' = HSMM_deg_ordering_genes)

pdf(paste(main_fig_dir, "fig1c.pdf", sep = ''), height = 2.5, width = 4)
upset(fromList(listInput), nsets = 6, order.by = "freq") + nm_theme()
dev.off()

################################################################################################################################################################################################################################################
## Heatmap for DPfeature selected genes 
################################################################################################################################################################################################################################################

HSMM_exprs_mat <- log(exprs(HSMM_myo[HSMM_myo_ordering_genes, ]) + 1)

#scale and trim the data: 
scale_max=3 
scale_min=-3

HSMM_exprs_mat=HSMM_exprs_mat[!apply(HSMM_exprs_mat, 1, sd)==0,]
HSMM_exprs_mat=Matrix::t(scale(Matrix::t(HSMM_exprs_mat),center=TRUE))
HSMM_exprs_mat=HSMM_exprs_mat[is.na(row.names(HSMM_exprs_mat)) == FALSE,]
HSMM_exprs_mat[is.nan(HSMM_exprs_mat)] = 0
HSMM_exprs_mat[HSMM_exprs_mat>scale_max] = scale_max
HSMM_exprs_mat[HSMM_exprs_mat<scale_min] = scale_min

HSMM_exprs_mat_ori <- HSMM_exprs_mat
HSMM_exprs_mat <- HSMM_exprs_mat[is.finite(HSMM_exprs_mat[, 1]), ] #remove the NA fitting failure genes for each branch 


cols_dist <- as.dist((1 - cor(as.matrix(HSMM_exprs_mat)))/2)
cols_dist[is.na(cols_dist)] <- 1

hclust_method <-"ward.D2"

ph <- pheatmap(HSMM_exprs_mat, 
               useRaster = T,
               cluster_cols=T, 
               cluster_rows=T, 
               show_rownames=F, 
               show_colnames=F, 
               #scale="row",
               clustering_distance_cols=cols_dist,
               clustering_method = hclust_method,
               cutree_cols = 4,
               silent=F
               # filename=NA,
               # color=hmcols
               #color=hmcols#,
               # filename="expression_pseudotime_pheatmap.pdf",
)

pData(HSMM_myo)$Time <- pData(HSMM[, colnames(HSMM_myo)])$Hours
annotation_col <- data.frame(#State=factor(pData(HSMM_myo[, ph$tree_col$labels])$State),
                             Time = factor(pData(HSMM_myo[, ph$tree_col$labels])$Time), 
                             'Density peak clusters'=factor(pData(HSMM_myo[, ph$tree_col$labels])$Cluster), 
                             row.names = ph$tree_col$labels)

annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, 8)))

annotation_col_time <- annotation_col$Time
names(annotation_col_time) <- unique(annotation_col$Time)
annotation_colors=list("Time"=HSMM_cols, "Density.peak.clusters" = c('1' = "#F3756C", '2' = "#26B24B", '3' = '#6E95CD', '4' = '#B27AB4'))

#F3756C

ph_res <- pheatmap(HSMM_exprs_mat, 
               useRaster = T,
               cluster_cols=T, 
               cluster_rows=T, 
               show_rownames=F, 
               show_colnames=F, 
               #scale="row",
               clustering_distance_cols=cols_dist,
               clustering_method = hclust_method,
               cutree_cols = 4,
               cutree_rows = 8, 
               annotation_col=annotation_col,
               annotation_row=annotation_row,
               annotation_colors=annotation_colors,
               width = 2.5,
               height = 5,
               silent=F
               # filename=NA,
               # color=hmcols
               #color=hmcols#,
               # filename="expression_pseudotime_pheatmap.pdf",
)

pdf(paste(main_fig_dir, "SI1d.pdf", sep = ''))
grid::grid.rect(gp = grid::gpar("fill", col = NA))
grid::grid.draw(ph_res$gtable)
dev.off()


pdf(paste(main_fig_dir, "SI1d_helper.pdf", sep = ''))
grid::grid.rect(gp = grid::gpar("fill", col = NA))
grid::grid.draw(ph_res$gtable)
dev.off()

#run the enrichment on a cluster of genes from the heatmap: 
HSMM_clustering <- data.frame(Cluster=factor(cutree(ph$tree_row, 8)))
clusters <- as.numeric(HSMM_clustering$Cluster)
names(clusters) <- fData(HSMM_myo[row.names(HSMM_clustering), ])$gene_short_name

hsmm_gsa_results_dp_genes <- collect_gsa_hyper_results(HSMM_myo[, ], human_go_gsc, clusters)

pdf(paste(main_fig_dir, "hsmm_gsa_results_dp_genes_go_enrichment.pdf", sep = ''), height=100, width=15)
plot_gsa_hyper_heatmap(HSMM_myo, hsmm_gsa_results_dp_genes, significance = 1e-5)
dev.off()

dp_cluster_group <- rep(1, length(HSMM_myo_ordering_genes))
names(dp_cluster_group) <- as.character(fData(HSMM_myo)[HSMM_myo_ordering_genes, 'gene_short_name'])
HSMM_cluster_hyper_geometric_results_go <- collect_gsa_hyper_results(HSMM_myo, gsc = human_go_gsc, dp_cluster_group)
plot_gsa_hyper_heatmap(HSMM_myo, HSMM_cluster_hyper_geometric_results_go, significance = 1e-50)

pdf(paste(main_fig_dir, "fig_si1d.1.pdf", sep = ''), height = 5, width = 4)
plot_gsa_hyper_heatmap(HSMM_myo, HSMM_cluster_hyper_geometric_results_go, significance = 1e-50) + monocle:::monocle_theme_opts() + nm_theme()
dev.off()

################################################################################################################################################################################################################################################
save.image('./RData/fig1.RData')
################################################################################################################################################################################################################################################
