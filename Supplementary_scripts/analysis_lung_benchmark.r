rm(list = ls())
####################################################################################################################################################################################
#load all package
####################################################################################################################################################################################
# Things need to do: 
# a. finish the benchmark for the MAR-seq data 
# b. get the correct branching ordering based on markers (difficult to do)
# c. make the trajectory plots for all different software

# 1. Include the very nice table you had on your slides yesterday comparing the different features of the programs. You will need to reformat it so that it looks nice. I would suggest making it in Keynote and then taking a high-resolution screenshot for the AI file
# 2. Should show example trajectory plots for dpt, wishbone, and slicer, on the muscle data as Panel A. This will serve both to show qualitative differences between the programs and be a visual cue to the reader that we are comparing to other programs.
# 3. The second panel should be a performance comparison (like the current panels A-B), but with Monocle 1, Monocle 2, Wishbone, and dpt. We need to be able to say who is most accurate with respect to the manually ordered cells.
# 4. Move panels C and D to the supplement.
# 5. Make panels E and F smaller - they are much too large right now.
# 6. Use a solid color fill on the boxplots in panels E and F, with black outlines and outlier dots. This is the more standard way to color boxplots.
# 7. Move the PCA/LLE/ICA comparisons to the supplement or to the end of this figure. I think you already have example trajectories that come out of DDRTree that show how they are different for each initialization method, but they need to be clearly labeled so the reader can see that in most cases, you really do get the same result from DDRTree no matter how you initialize (except with LLE).
# 8. Make panels that show the metrics in panels C and D for a range of values of ncenter, param.gamma, and the number of dimensions used in the reduction. 
# 9. We are going to need to include some discussions of additional datasets.
################################################################################################################################################################################################################################################
library(destiny)
library(dpt)
library(monocle)
library(SLICER)
library(xacHelper)
library(diffusionMap)

load('./RData/fig1.RData')
load('./RData/fig2.RData')

source('./scripts/function.R')
source('./scripts/plotting.R')
main_fig_dir <- "./Figures/main_figures/"
SI_fig_dir <- "./Figures/supplementary_figures/"

################################################################################################################################################################################################################################################
# compare the trajectory plot for all five softwares on lung data
# generate panels a - e 
################################################################################################################################################################################################################################################

# monocle2
absolute_cds <- recreate_cds(absolute_cds)
absolute_cds <- setOrderingFilter(absolute_cds, quake_id)
absolute_cds <- reduceDimension(absolute_cds, norm_method = 'log', auto_param_selection = T)
absolute_cds <- orderCells(absolute_cds)

plot_cell_trajectory(absolute_cds, color_by="Pseudotime")
plot_cell_trajectory(absolute_cds, color_by="Time")
plot_cell_trajectory(absolute_cds, color_by="State")
absolute_cds <- orderCells(absolute_cds, root_state = 3)

pdf(paste(main_fig_dir, 'fig3a.pdf', sep = ''), height = 1.5, width = 1.5)
plot_cell_trajectory(absolute_cds, color_by="Time", show_branch_points = F) + xacHelper::nm_theme() 
dev.off()

# monocle1
absolute_cds_monocle1 <- reduceDimension(absolute_cds, norm_method = 'log', reduction_method = 'ICA')
absolute_cds_monocle1 <- orderCells(absolute_cds_monocle1, num_paths = 2)

plot_cell_trajectory(absolute_cds_monocle1, color_by="Pseudotime")
plot_cell_trajectory(absolute_cds_monocle1, color_by="Time")
plot_cell_trajectory(absolute_cds_monocle1, color_by="State")
absolute_cds_monocle1 <- orderCells(absolute_cds_monocle1, root_state = 3)

pdf(paste(main_fig_dir, 'fig3b.pdf', sep = ''), height = 1.5, width = 1.5)
plot_cell_trajectory(absolute_cds_monocle1, color_by="Time") + xacHelper::nm_theme() 
dev.off()

# dpt
lung_dpt_res <- run_dpt(absolute_cds[quake_id, ], color_by = 'Time')
pdf(paste(main_fig_dir, 'fig3c.pdf', sep = ''), height = 1.5, width = 1.5)
lung_dpt_res$p1 + xacHelper::nm_theme() 
dev.off()

plot.DPT <- function (x, root = NULL, paths_to = integer(0L), dcs = 1:2, 
                      divide = integer(0L), w_width = 0.1, col_by = "dpt", col_path = rev(palette()), 
                      col_tip = "red", ..., col = NULL, legend_main = col_by) 
{
  dpt <- x
  dpt_flat <- branch_divide(dpt, divide)
  if (!is.null(root) && length(root) < 1L) 
    stop("root needs to be specified")
  root <- if (is.null(root)) 
    min(dpt_flat@branch[, 1], na.rm = TRUE)
  else as.integer(root)
  paths_to <- as.integer(paths_to)
  if (length(root) > 1L && length(paths_to) > 0L) 
    stop("(length(root), length(paths_to)) needs to be (1, 0-n) or (2-n, 0), but is (", #what is n there? 
         length(root), ", ", length(paths_to), ")")
  stopifnot(length(dcs) %in% 2:3)
  if (length(root) > 1L && length(paths_to) == 0L) {
    paths_to <- root[-1]
    root <- root[[1]]
  }
  pt_vec <- dpt_for_branch(dpt_flat, root)
  evs <- flipped_dcs(dpt@dm, dcs)
  plot_paths <- function(p, ..., rescale) {
    plot_points <- get_plot_fn(p)
    rescale_fun <- if (is.null(rescale)) 
      identity
    else function(x) rescale_mat(x, from = rescale$from, 
                                 to = rescale$to)
    for (b in seq_along(paths_to)) {
      idx <- dpt@branch[, 1] %in% c(root, paths_to[[b]])
      path <- average_path(pt_vec[idx], evs[idx, ], w_width)
      plot_points(rescale_fun(path), type = "l", col = col_path[[b]], 
                  ...)
    }
    tips <- evs[dpt_flat@tips[, 1], ]
    plot_points(rescale_fun(tips), col = col_tip, ...)
  }
  col <- if (!is.null(col)) 
    col
  else switch(col_by, dpt = pt_vec, branch = , Branch = dpt_flat@branch[, 
                                                                        1], dpt[[col_by]])
  legend_main <- switch(legend_main, dpt = "DPT", branch = "Branch", 
                        legend_main)
  args <- list(dpt@dm, dcs, plot_more = plot_paths, legend_main = legend_main, 
               col = col, ...)
  if (!identical(Sys.getenv("LOG_LEVEL"), "")) 
    message("Args:\n", paste(capture.output(print(args)), 
                             collapse = "\n"))
  do.call(plot, args)
}

pdf(paste(main_fig_dir, 'fig3c.2.pdf', sep = ''))
plot(lung_dpt_res$DPT, root = 2, paths_to = c(1,3), col_by = 'branch', pch = 20)
dev.off()

# slicer 
lung_slicer_res <- run_slicer(abs_AT12_cds_subset_all_gene[quake_id, ], start = which(pData(absolute_cds)$Pseudotime == 0)) 
pdf(paste(main_fig_dir, 'fig3d.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(lung_slicer_res$traj_lle[, 1], lung_slicer_res$traj_lle[, 2], color = as.character(pData(abs_AT12_cds_subset_all_gene)$Time) ) + xacHelper::nm_theme() +  #as.character(lung_slicer_res$order_df$branches)
  xlab('LLE dim 1') + ylab('LLE dim 2') 
dev.off()

# wishbone (generate from the wishbone software)
# check analysis_wishbone.py
data <- t(log2(exprs(absolute_cds)[fData(absolute_cds)$use_for_ordering, ] + 1))
write.csv(file = paste('./csv_data/Wishbone_test_data/', 'lung_data_full', ".txt", sep = ''), as.matrix(data), quote = F, row.names = T)
row.names(subset(pData(absolute_cds), Pseudotime == 0))
# [1] "SRR1034041_0"

wishbone_res_lung <- read.table('./wishbone/lung_data_full_wishbone_res.txt', header = T, sep = '\t')
qplot(dm1, dm2, data = wishbone_res_lung)
qplot(tSNE1, tSNE2, data = wishbone_res_lung, color = as.character(pData(absolute_cds_monocle1)$Time), size = 0.5) 

pdf(paste(main_fig_dir, 'lung_data_full_wishbone_res.pdf', sep = ''), height = 1.5, width = 1.5)
qplot(tSNE1, tSNE2, data = wishbone_res_lung, color = pData(absolute_cds_monocle1)$Time, size = 0.5) + nm_theme() + scale_size(range = c(0.5, 0.5))
# plot_cell_trajectory(valid_subset_GSE72857_cds, color_by = 'cell_type', cell_size = 0.5, show_branch_points = F) + nm_theme()#, cell_size = pData(valid_subset_GSE72857_cds)$Pseudotime
dev.off()

################################################################################################################################################################################################################################################
#panel g, h: pseudotime and branch robustness for 0.8 downsampling as well as progressive downsampling 
################################################################################################################################################################################################################################################
auto_param_selection <- F
repeat_downsampling_num <- 100
benchmark_type <- 'lung_data' #

root_state <- row.names(subset(pData(absolute_cds), State == 1))
AT1_state <- row.names(subset(pData(absolute_cds), State == 2))
AT2_state <- row.names(subset(pData(absolute_cds), State == 3))

source('./script/DDRTree_robustness_analysis_genearlize.R', echo = T)
# be careful to run this file, it takes a long time! 
source('./script_for_reproduce/ICA_robustness_analysis_genearlize.R', echo = T) 
load('lung_ICA_downsampling_res') #save the data 
source('./script/robustness_dpt_slicer_genearlize.R', echo = T)

downsampling_state <- 2 #check this 
source('./script/progressive_branch_downsampling_generalize.R', echo = T)

################################################################################################################################################################################################################################################
#fig_si panel a: different initialization methods (not used in the paper)
################################################################################################################################################################################################################################################

################################################################################################################################################################################################################################################
#fig_si panel c: benchmark on testing different parameters: ncenter, param.gamma, and the number of dimensions used in the reduction
################################################################################################################################################################################################################################################

################################################################################################################################################################################################################################################
save.image('./RData/fig_SI4_lung.RData')
################################################################################################################################################################################################################################################

