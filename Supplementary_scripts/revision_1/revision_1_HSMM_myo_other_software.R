library(monocle)
library(destiny)

# run dpt and wishbone on HSMM_myo dataset: 
load('../RData/fig2.RData')
benchmark_type <- 'HSMM_myo'
dpt_HSMM_myo_res <- run_dpt(HSMM_myo[fData(HSMM_myo)$use_for_ordering, ])

dpt_res <- data.frame(DC1 = dpt_HSMM_myo_res$dm$DC1, DC2 = dpt_HSMM_myo_res$dm$DC2, branch = as.factor(dpt_HSMM_myo_res$pt$Branch))

revision_1_fig_dir <- "../Figures/First_revision/"

pdf(paste(revision_1_fig_dir, benchmark_type, '_dpt.pdf', sep = ''), width = 1.2, height = 1) 
qplot(DC1, DC2, color = branch, data = dpt_res, size = I(0.5)) + nm_theme() + theme(axis.text.x = element_text(angle = 30))
dev.off()

pdf(paste(revision_1_fig_dir, benchmark_type, '_dpt_helper.pdf', sep = '')) 
qplot(DC1, DC2, color = branch, data = dpt_res, size = I(0.5))  + theme(axis.text.x = element_text(angle = 30))
dev.off()

# save data for running wishbone
write.csv(file = paste('./csv_data/Wishbone_test_data/', benchmark_type, '_data.txt', sep = ''), as.matrix(exprs(HSMM_myo[fData(HSMM_myo)$use_for_ordering, ])), quote = F, row.names = T)
# run the same thing with wishbone: 
wishbone_HSMM_res <- read.table('./wishbone/wishbone_HSMM.txt', row.names = 1, header = T, sep = '\t')
wishbone_HSMM_res$branch <- as.character(wishbone_HSMM_res$branch)

pdf(paste(revision_1_fig_dir, benchmark_type, '_dpt_wishbone.pdf', sep = ''), width = 1.2, height = 1) 
qplot(tSNE1, tSNE2, color = branch, data = wishbone_HSMM_res, size = I(0.5)) + nm_theme()
dev.off()

pdf(paste(revision_1_fig_dir, benchmark_type, '_dpt_wishbone_helper.pdf', sep = '')) 
qplot(tSNE1, tSNE2, color = branch, data = wishbone_HSMM_res, size = I(0.5))
dev.off()

calClusteringMetrics(wishbone_HSMM_res$branch, pData(HSMM_myo)$State)
calClusteringMetrics(dpt_HSMM_myo_res$pt$Branch, pData(HSMM_myo)$State)

cor(wishbone_HSMM_res$trajectory, pData(HSMM_myo)$Pseudotime)
cor(dpt_HSMM_myo_res$pt$DPT, pData(HSMM_myo)$Pseudotime)

# make the plot here: 
HSMM_res <- data.frame(Type = c("ARI", "ARI", "Pearson \n  correlation", "Pearson \n  correlation"), Software = c('DPT', "Wishbone", "DPT", "Wishbone"), Value = c(0.80, 0.49, 0.94, 0.76))

pdf(paste(revision_1_fig_dir, benchmark_type, '_dpt_wishbone_comparison.pdf', sep = ''), width = 1.2, height = 1) 
ggplot(aes(Software, Value), data = HSMM_res) + geom_bar(aes(fill = Software), stat = 'identity') + facet_wrap(~Type) + nm_theme() + 
  xlab('') + ylab('') + theme(axis.text.x = element_text(angle = 30)) + ylim(c(0,1))
dev.off()

save.image('../RData/revision_1_HSMM_myo_other_software.RData')

##################################################################################################################################################################
# analyze the URMM dataset (not required )
##################################################################################################################################################################
