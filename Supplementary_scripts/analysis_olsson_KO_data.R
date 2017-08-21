#analysis for the ko dataset
rm(list = ls())
########################################################################################################################################################################################################
# load data
########################################################################################################################################################################################################
library(destiny)
library(RColorBrewer)
library(stringr)
library(UpSetR)
library(devtools)
library(rEDM)
library(xlsx)
library(xacHelper)
library(reshape2)
library(pheatmap)
library(d3Network)
library(netbiov)
library(monocle)
library(plyr)
library(dplyr)
library(venneuler)

load('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/RData/fig4.RData')
# load('/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/RData/fig4.RData')
# load('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/RData/fig5.RData')
load('/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/RData/fig5_7_8.RData')

main_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/Figures/main_figures/"
SI_fig_dir <- "/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/Figures/supplementary_figures/"

########################################################################################################################################################################################################
plot_cell_trajectory(URMM_all_abs)
plot_cell_trajectory(URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout')], color_by = 'Type', show_branch_points = F, cell_size = 0.5, cell_link_size = 0.3, theta = 90) + facet_wrap(~Type, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = type_cols, name = "Type")  + theme (legend.position="none", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) #theme(axis.text.x = element_text(angle = 30, hjust = 1))

granulocyte_state <- as.numeric(names(which.max(table(pData(URMM_all_abs)$State, pData(URMM_all_abs)$paper_cluster)[, 'Myelocyte'])))
monocyte_state <- as.numeric(names(which.max(table(pData(URMM_all_abs)$State, pData(URMM_all_abs)$paper_cluster)[, 'Mono'])))
fig1_granulocyte_state <- as.numeric(names(which.max(table(pData(URMM_all_fig1b)$State, pData(URMM_all_fig1b)$paper_cluster)[, 'Myelocyte'])))
fig1_monocyte_state <- as.numeric(names(which.max(table(pData(URMM_all_fig1b)$State, pData(URMM_all_fig1b)$paper_cluster)[, 'Mono'])))

#WT: state 3: granulocyte; state 4: monocyte; all: State 4 (G) and State 5 (Mono)
########################################################################################################################################################################################################
# test on leaking cells
########################################################################################################################################################################################################
# Irf8 leaking cells
########################################################################################################################################################################################################
rng <- range(URMM_all_abs[, pData(URMM_all_abs)$State == monocyte_state & pData(URMM_all_abs)$Type %in% c('Irf8_knockout')]$Pseudotime)
Irf8_ko_leaking_cells <- URMM_all_abs[, pData(URMM_all_abs)$State == monocyte_state & pData(URMM_all_abs)$Type %in% c('Irf8_knockout')]
# [1] 14.14190 14.26459

Irf8_ko_mainfold_cells <- URMM_all_abs[, pData(URMM_all_abs)$Pseudotime > floor(rng[1]) & pData(URMM_all_abs)$Pseudotime < ceiling(rng[2]) & pData(URMM_all_abs)$State == monocyte_state & !(pData(URMM_all_abs)$Type %in% c('Gfi1_knockout', 'Gfi1_Irf8_knockout'))]
pData(Irf8_ko_mainfold_cells)$fig6_ko_test1 <- 'WildType'
pData(Irf8_ko_mainfold_cells)$fig6_ko_test1[pData(Irf8_ko_mainfold_cells)$Type %in% c('Irf8_knockout')] <- "leaking_cells"
# URMM_all_abs_genotype_fig6_ko_test1_beam_genes <- differentialGeneTest(Irf8_ko_mainfold_cells, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3) + fig6_ko_test1", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)")
URMM_all_abs_genotype_fig6_ko_test1_group_deg <- differentialGeneTest(Irf8_ko_mainfold_cells, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~fig6_ko_test1", reducedModelFormulaStr = "~1")

# > sum(URMM_all_abs_genotype_fig6_ko_test1_group_deg$qval < 0.1)
# [1] 1
# > URMM_all_abs_genotype_fig6_ko_test1_beam_genes[URMM_all_abs_genotype_fig6_ko_test1_beam_genes$qval < 0.1, ]
# status           family         pval       qval gene_short_name use_for_ordering
# Ly86     OK negbinomial.size 5.571019e-07 0.01328521            Ly86             TRUE
# > URMM_all_abs_genotype_fig6_ko_test1_group_deg[URMM_all_abs_genotype_fig6_ko_test1_group_deg$qval < 0.1, ]
# status           family         pval        qval gene_short_name use_for_ordering
# Ly86     OK negbinomial.size 1.219926e-07 0.002922333            Ly86             TRUE
#########################################################
# test on the supposed branch:
#########################################################
rng <- range(URMM_all_abs[, pData(URMM_all_abs)$State == granulocyte_state & pData(URMM_all_abs)$Type %in% c('Irf8_knockout')]$Pseudotime)
# [1] 13.80390 17.10217
Irf8_ko_mainfold_cells.1 <- URMM_all_abs[, pData(URMM_all_abs)$State == granulocyte_state & pData(URMM_all_abs)$Type %in% c('Irf8_knockout')]

Irf8_ko_mainfold_cells.1 <- URMM_all_abs[, pData(URMM_all_abs)$Pseudotime > floor(rng[1]) & pData(URMM_all_abs)$Pseudotime < ceiling(rng[2]) & pData(URMM_all_abs)$State == granulocyte_state & !(pData(URMM_all_abs)$Type %in% c('Gfi1_knockout', 'Gfi1_Irf8_knockout'))]
pData(Irf8_ko_mainfold_cells.1)$fig6_ko_test1.1 <- 'WildType'
pData(Irf8_ko_mainfold_cells.1)$fig6_ko_test1.1[pData(Irf8_ko_mainfold_cells.1)$Type %in% c('Irf8_knockout')] <- "Irf8 KO"
# URMM_all_abs_genotype_fig6_ko_test1_beam_genes.1 <- differentialGeneTest(Irf8_ko_mainfold_cells.1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3) + fig6_ko_test1.1", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)")
URMM_all_abs_genotype_fig6_ko_test1_group_deg.1 <- differentialGeneTest(Irf8_ko_mainfold_cells.1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~fig6_ko_test1.1", reducedModelFormulaStr = "~1")

#########################################################
# test KO cells from two branches:
#########################################################
Irf8_ko_leaking_cells.2 <- URMM_all_abs[, pData(URMM_all_abs)$State %in% c(granulocyte_state, monocyte_state) & pData(URMM_all_abs)$Type %in% c('Irf8_knockout')]

URMM_all_abs_genotype_fig6_ko_test1_group_deg.2 <- differentialGeneTest(Irf8_ko_leaking_cells.2, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~State", reducedModelFormulaStr = "~1")

########################################################################################################################################################################################################
# Gfi1 leaking cells
########################################################################################################################################################################################################
rng <- range(URMM_all_abs[, pData(URMM_all_abs)$State == granulocyte_state & pData(URMM_all_abs)$Type %in% c('Gfi1_knockout')]$Pseudotime)
#[1] 13.80905 14.35581
Gif1_ko_leaking_cells <- URMM_all_abs[, pData(URMM_all_abs)$State == granulocyte_state & pData(URMM_all_abs)$Type %in% c('Gfi1_knockout')]

Gfi1_ko_mainfold_cells <- URMM_all_abs[, pData(URMM_all_abs)$Pseudotime > floor(rng[1]) & pData(URMM_all_abs)$Pseudotime < ceiling(rng[2]) & pData(URMM_all_abs)$State == granulocyte_state & !(pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_Irf8_knockout'))]
pData(Gfi1_ko_mainfold_cells)$fig6_ko_test2 <- 'WildType'
pData(Gfi1_ko_mainfold_cells)$fig6_ko_test2[pData(Gfi1_ko_mainfold_cells)$Type %in% c('Gfi1_knockout')] <- "leaking_cells"
# URMM_all_abs_genotype_fig6_ko_test2_beam_genes <- differentialGeneTest(Gfi1_ko_mainfold_cells, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3) + fig6_ko_test2", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)")
URMM_all_abs_genotype_fig6_ko_test2_group_deg <- differentialGeneTest(Gfi1_ko_mainfold_cells, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~fig6_ko_test2", reducedModelFormulaStr = "~1")

#########################################################
# test on the supposed branch:
#########################################################
rng <- range(URMM_all_abs[, pData(URMM_all_abs)$State == monocyte_state & pData(URMM_all_abs)$Type %in% c('Gfi1_knockout')]$Pseudotime)
# [1] 13.75441 16.13564
Gif1_ko_leaking_cells.1 <- URMM_all_abs[, pData(URMM_all_abs)$State == monocyte_state & pData(URMM_all_abs)$Type %in% c('Gfi1_knockout')]

Gif1_ko_leaking_cells.1 <- URMM_all_abs[, pData(URMM_all_abs)$Pseudotime > floor(rng[1]) & pData(URMM_all_abs)$Pseudotime < ceiling(rng[2]) & pData(URMM_all_abs)$State == monocyte_state & !(pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_Irf8_knockout'))]
pData(Gif1_ko_leaking_cells.1)$fig6_ko_test2.1 <- 'WildType'
pData(Gif1_ko_leaking_cells.1)$fig6_ko_test2.1[pData(Gif1_ko_leaking_cells.1)$Type %in% c('Gfi1_knockout')] <- "Gfi1 KO"
# URMM_all_abs_genotype_fig6_ko_test2_beam_genes.1 <- differentialGeneTest(Gif1_ko_leaking_cells.1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3) + fig6_ko_test2.1", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)")
URMM_all_abs_genotype_fig6_ko_test2_group_deg.1 <- differentialGeneTest(Gif1_ko_leaking_cells.1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~fig6_ko_test2.1", reducedModelFormulaStr = "~1")

#########################################################
# test KO cells from two branches:
#########################################################
Gif1_ko_leaking_cells.2 <- URMM_all_abs[, pData(URMM_all_abs)$State %in% c(monocyte_state, granulocyte_state) & pData(URMM_all_abs)$Type %in% c('Gfi1_knockout')]

URMM_all_abs_genotype_fig6_ko_test2_group_deg.2 <- differentialGeneTest(Gif1_ko_leaking_cells.2, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~State", reducedModelFormulaStr = "~1")

########################################################################################################################################################################################################
# double ko leaking cells (monocyte_state)
########################################################################################################################################################################################################
rng <- range(URMM_all_abs[, pData(URMM_all_abs)$State == monocyte_state & pData(URMM_all_abs)$Type %in% c('Gfi1_Irf8_knockout')]$Pseudotime)
double_ko_leaking_cells_state2 <- URMM_all_abs[, pData(URMM_all_abs)$State == monocyte_state & pData(URMM_all_abs)$Type %in% c('Gfi1_Irf8_knockout')]
# [1] 13.61449 14.28285

double_ko_leaking_cells_state2 <- URMM_all_abs[, pData(URMM_all_abs)$Pseudotime > floor(rng[1]) - 2 & pData(URMM_all_abs)$Pseudotime < ceiling(rng[2]) + 2 & pData(URMM_all_abs)$State == monocyte_state & !(pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout'))]
pData(double_ko_leaking_cells_state2)$fig6_ko_test3 <- 'WildType'
pData(double_ko_leaking_cells_state2)$fig6_ko_test3[pData(double_ko_leaking_cells_state2)$Type %in% c('Gfi1_Irf8_knockout')] <- "double KO"
# URMM_all_abs_genotype_fig6_ko_test3_beam_genes <- differentialGeneTest(double_ko_leaking_cells_state2, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3) + fig6_ko_test3", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)")
URMM_all_abs_genotype_fig6_ko_test3_group_deg <- differentialGeneTest(double_ko_leaking_cells_state2, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~fig6_ko_test3", reducedModelFormulaStr = "~1")

#########################################################
# test with Irf8 ko genes on granulocyte_state
#########################################################
rng <- range(URMM_all_abs[, pData(URMM_all_abs)$State == granulocyte_state & pData(URMM_all_abs)$Type %in% c('Gfi1_Irf8_knockout')]$Pseudotime)
double_ko_leaking_cells_state2.1 <- URMM_all_abs[, pData(URMM_all_abs)$State == granulocyte_state & pData(URMM_all_abs)$Type %in% c('Gfi1_Irf8_knockout')]
# [1] 13.81074 15.98827

double_ko_leaking_cells_state2.1 <- URMM_all_abs[, pData(URMM_all_abs)$Pseudotime > floor(rng[1]) & pData(URMM_all_abs)$Pseudotime < ceiling(rng[2]) & pData(URMM_all_abs)$State == granulocyte_state & pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_Irf8_knockout')]
pData(double_ko_leaking_cells_state2.1)$fig6_ko_test3.1 <- 'Irf8 KO'
pData(double_ko_leaking_cells_state2.1)$fig6_ko_test3.1[pData(double_ko_leaking_cells_state2.1)$Type %in% c('Gfi1_Irf8_knockout')] <- "Double KO"
# URMM_all_abs_genotype_fig6_ko_test3.1_beam_genes <- differentialGeneTest(double_ko_leaking_cells_state2, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3) + fig6_ko_test3.1", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)")
URMM_all_abs_genotype_fig6_ko_test3.1_group_deg <- differentialGeneTest(double_ko_leaking_cells_state2.1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~fig6_ko_test3.1", reducedModelFormulaStr = "~1")

########################################################################################################################################################################################################
# double ko leaking cells (granulocyte_state)
########################################################################################################################################################################################################
rng <- range(URMM_all_abs[, pData(URMM_all_abs)$State == granulocyte_state & pData(URMM_all_abs)$Type %in% c('Gfi1_Irf8_knockout')]$Pseudotime)
double_ko_leaking_cells_state3 <- URMM_all_abs[, pData(URMM_all_abs)$State == granulocyte_state & pData(URMM_all_abs)$Type %in% c('Gfi1_Irf8_knockout')]
# [1] 13.81074 15.98827

double_ko_leaking_cells_state3 <- URMM_all_abs[, pData(URMM_all_abs)$Pseudotime > floor(rng[1]) - 2 & pData(URMM_all_abs)$Pseudotime < ceiling(rng[2]) + 2 & pData(URMM_all_abs)$State == granulocyte_state & !(pData(URMM_all_abs)$Type %in% c('Gfi1_knockout', 'Irf8_knockout'))]
pData(double_ko_leaking_cells_state3)$fig6_ko_test4 <- 'WildType'
pData(double_ko_leaking_cells_state3)$fig6_ko_test4[pData(double_ko_leaking_cells_state3)$Type %in% c('Gfi1_Irf8_knockout')] <- "Double KO"
# URMM_all_abs_genotype_fig6_ko_test4_beam_genes <- differentialGeneTest(double_ko_leaking_cells_state3, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3) + fig6_ko_test4", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)")
URMM_all_abs_genotype_fig6_ko_test4_group_deg <- differentialGeneTest(double_ko_leaking_cells_state3, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~fig6_ko_test4", reducedModelFormulaStr = "~1")

#########################################################
# test with Gfi1 ko genes on monocyte_state
#########################################################
rng <- range(URMM_all_abs[, pData(URMM_all_abs)$State == monocyte_state & pData(URMM_all_abs)$Type %in% c('Gfi1_Irf8_knockout')]$Pseudotime)
# double_ko_leaking_cells_state2 <- URMM_all_abs[, pData(URMM_all_abs)$State == granulocyte_state & pData(URMM_all_abs)$Type %in% c('Gfi1_knockout')]
# [1] 13.61449 14.28285

double_ko_leaking_cells_state3.1 <- URMM_all_abs[, pData(URMM_all_abs)$Pseudotime > floor(rng[1]) & pData(URMM_all_abs)$Pseudotime < ceiling(rng[2]) & pData(URMM_all_abs)$State == monocyte_state & (pData(URMM_all_abs)$Type %in% c('Gfi1_knockout', 'Gfi1_Irf8_knockout'))]
pData(double_ko_leaking_cells_state3.1)$fig6_ko_test4.1 <- 'Gfi1 KO'
pData(double_ko_leaking_cells_state3.1)$fig6_ko_test4.1[pData(double_ko_leaking_cells_state3.1)$Type %in% c('Gfi1_Irf8_knockout')] <- "Double KO"
# URMM_all_abs_genotype_fig6_ko_test4.1_beam_genes <- differentialGeneTest(double_ko_leaking_cells_state3, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3) + fig6_ko_test4.1", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)")
URMM_all_abs_genotype_fig6_ko_test4.1_group_deg <- differentialGeneTest(double_ko_leaking_cells_state3.1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~fig6_ko_test4.1", reducedModelFormulaStr = "~1")

#########################################################
# test double KO cells from two branches:
#########################################################
double_ko_leaking_cells <- URMM_all_abs[, pData(URMM_all_abs)$State %in% c(monocyte_state, granulocyte_state) & pData(URMM_all_abs)$Type %in% c('Gfi1_Irf8_knockout')]

URMM_all_abs_genotype_fig6_double_ko_deg <- differentialGeneTest(double_ko_leaking_cells, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~State", reducedModelFormulaStr = "~1")

########################################################################################################################################################################################################
# create the barplot to demonstrate the following points:
# 1. leaking cells are not different comparing to the WT cells
# 2. double ko cells have a little more DEG comparing to the single ko cells
########################################################################################################################################################################################################
Irf8_ko_deg_df <- data.frame(leaking_cells = sum(URMM_all_abs_genotype_fig6_ko_test1_group_deg$qval < 0.1),
                             proper_positioned_cells = sum(URMM_all_abs_genotype_fig6_ko_test1_group_deg.1$qval < 0.1),
                             cells_between_two_branches = sum(URMM_all_abs_genotype_fig6_ko_test1_group_deg.2$qval < 0.1))
Gfi1_ko_deg_df <- data.frame(leaking_cells = sum(URMM_all_abs_genotype_fig6_ko_test2_group_deg$qval < 0.1),
                             proper_positioned_cells = sum(URMM_all_abs_genotype_fig6_ko_test2_group_deg.1$qval < 0.1),
                             cells_between_two_branches = sum(URMM_all_abs_genotype_fig6_ko_test2_group_deg.2$qval < 0.1))
double_ko_deg_df <- data.frame(state2vsWT = sum(URMM_all_abs_genotype_fig6_ko_test3_group_deg$qval < 0.1),
                               state2vsIrf8_ko = sum(URMM_all_abs_genotype_fig6_ko_test3.1_group_deg$qval < 0.1),
                               state3vsWT = sum(URMM_all_abs_genotype_fig6_ko_test4_group_deg$qval < 0.1),
                               state3vsGfi1_ko = sum(URMM_all_abs_genotype_fig6_ko_test4.1_group_deg$qval < 0.1),
                               cells_between_two_branches = sum(URMM_all_abs_genotype_fig6_double_ko_deg$qval < 0.1))
###### create the three bar plots:
mlt_Irf8_ko_deg_df <- melt(Irf8_ko_deg_df)
pdf(paste(main_fig_dir, 'fig6e.1.pdf', sep = ''), width = 1.5, height = 2)
ggplot(aes(variable, value), data = mlt_Irf8_ko_deg_df) + geom_bar(aes(fill = variable), stat = 'identity') + xlab('DEG test type') + ylab('Number of DEGs') + nm_theme() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

mlt_Gfi1_ko_deg_df <- melt(Gfi1_ko_deg_df)
pdf(paste(main_fig_dir, 'fig6e.2.pdf', sep = ''), width = 1.5, height = 2)
ggplot(aes(variable, value), data = mlt_Gfi1_ko_deg_df) + geom_bar(aes(fill = variable), stat = 'identity') + xlab('DEG test type') + ylab('Number of DEGs') + nm_theme() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

ggplot(aes(variable, value), data = mlt_Gfi1_ko_deg_df) + geom_bar(aes(fill = variable), stat = 'identity') + xlab('DEG test type') + ylab('Number of DEGs')

mlt_double_ko_deg_df <- melt(double_ko_deg_df)
pdf(paste(main_fig_dir, 'fig6e.3.pdf', sep = ''), width = 1.5, height = 2)
ggplot(aes(variable, value), data = mlt_double_ko_deg_df) + geom_bar(aes(fill = variable), stat = 'identity') + xlab('DEG test type') + ylab('Number of DEGs') + nm_theme() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

########################################################################################################################################################################################################
# create a upSetR plot for the double ko genes:
########################################################################################################################################################################################################
listInput <- list(double_ko_MONOCYTEvsWT = row.names(subset(URMM_all_abs_genotype_fig6_ko_test3_group_deg, qval < 0.1)),
                  double_ko_MONOCYTEvsIrf8_ko = row.names(subset(URMM_all_abs_genotype_fig6_ko_test3.1_group_deg, qval < 0.1)),
                  double_ko_GRANUOCYTEvsWT = row.names(subset(URMM_all_abs_genotype_fig6_ko_test4_group_deg, qval < 0.1)),
                  double_ko_GRANUOCYTEvsGfi1_ko = row.names(subset(URMM_all_abs_genotype_fig6_ko_test4.1_group_deg, qval < 0.1)),
                  MONOCYTE_proper_positioned_cells = row.names(subset(URMM_all_abs_genotype_fig6_ko_test1_group_deg.1, qval < 0.1)),
                  GRANUOCYTE_proper_positioned_cells = row.names(subset(URMM_all_abs_genotype_fig6_ko_test2_group_deg.1, qval < 0.1)))

pdf(paste(SI_fig_dir, "branch_knockout_gates_upset_groups_all.pdf", sep = ''), height = 3, width = 7)
upset(fromList(listInput), nsets = 6, order.by = "freq") + nm_theme()
dev.off()

########################################################################################################################################################################################################
# overlap with Chip-seq data
########################################################################################################################################################################################################
Irf8_targets_gene_short_name
Gfi1_targets_gene_short_name

listInput <- list(Irf8_targets = Irf8_targets_gene_short_name[!is.na(Irf8_targets_gene_short_name)],
                  Gfi1_targets = Gfi1_targets_gene_short_name[!is.na(Gfi1_targets_gene_short_name)],
                  double_ko_MONOCYTEvsWT = row.names(subset(URMM_all_abs_genotype_fig6_ko_test3_group_deg, qval < 0.1)),
                  double_ko_MONOCYTEvsIrf8_ko = row.names(subset(URMM_all_abs_genotype_fig6_ko_test3.1_group_deg, qval < 0.1)),
                  double_ko_GRANUOCYTEvsWT = row.names(subset(URMM_all_abs_genotype_fig6_ko_test4_group_deg, qval < 0.1)),
                  double_ko_GRANUOCYTEvsGfi1_ko = row.names(subset(URMM_all_abs_genotype_fig6_ko_test4.1_group_deg, qval < 0.1)),
                  MONOCYTE_proper_positioned_cells = row.names(subset(URMM_all_abs_genotype_fig6_ko_test1_group_deg.1, qval < 0.1)),
                  GRANUOCYTE_proper_positioned_cells = row.names(subset(URMM_all_abs_genotype_fig6_ko_test2_group_deg.1, qval < 0.1)))

pdf(paste(SI_fig_dir, "chip_seq_data_knockout_gates_upset_groups_all.pdf", sep = ''), height = 4, width = 7)
upset(fromList(listInput), nsets = 8, order.by = "freq") + nm_theme()
dev.off()

listInput <- list(Irf8_targets = Irf8_targets_gene_short_name[!is.na(Irf8_targets_gene_short_name)],
                 Gfi1_targets = Gfi1_targets_gene_short_name[!is.na(Gfi1_targets_gene_short_name)],
                 beam_genes = row.names(subset(fig1b_beam_genes_proj_dup, qval < 1e-1)))

pdf(paste(SI_fig_dir, "chip_seq_data_beam_genes.pdf", sep = ''), height = 4, width = 7)
upset(fromList(listInput), nsets = 3, order.by = "freq") + nm_theme()
dev.off()

double_ko_MONOCYTEvsWT <- row.names(subset(URMM_all_abs_genotype_fig6_ko_test3_group_deg, qval < 0.1))
double_ko_GRANUOCYTEvsWT <- row.names(subset(URMM_all_abs_genotype_fig6_ko_test4_group_deg, qval < 0.1))

#double ko VS single-ko reveals the targets from the other genes
double_ko_MONOCYTEvsIrf8_ko <- row.names(subset(URMM_all_abs_genotype_fig6_ko_test3.1_group_deg, qval < 0.1))
double_ko_GRANUOCYTEvsGfi1_ko <- row.names(subset(URMM_all_abs_genotype_fig6_ko_test4.1_group_deg, qval < 0.1))
intersect(double_ko_MONOCYTEvsIrf8_ko, Gfi1_targets_gene_short_name)
intersect(double_ko_GRANUOCYTEvsGfi1_ko, Irf8_targets_gene_short_name)

MONOCYTE_proper_positioned_cells <- row.names(subset(URMM_all_abs_genotype_fig6_ko_test1_group_deg.1, qval < 0.1))
GRANUOCYTE_proper_positioned_cells <- row.names(subset(URMM_all_abs_genotype_fig6_ko_test2_group_deg.1, qval < 0.1))

intersect(MONOCYTE_proper_positioned_cells, Irf8_targets_gene_short_name)
intersect(GRANUOCYTE_proper_positioned_cells, Gfi1_targets_gene_short_name)

###### create two bar plots with the chip-seq data:
single_ko_deg_df <- data.frame(Type = rep(c('Irf8_ko', 'Gfi1_ko'), time = 2),
                               Number = c(Irf8_ko_deg_df$proper_positioned_cells,
                                          Gfi1_ko = Gfi1_ko_deg_df$proper_positioned_cells,
                                          length(intersect(MONOCYTE_proper_positioned_cells, Irf8_targets_gene_short_name)),
                                          length(intersect(GRANUOCYTE_proper_positioned_cells, Gfi1_targets_gene_short_name))),
                              subclass = rep(c('ko', 'chip'), each = 2))
single_ko_deg_df$Number[1:2] <- single_ko_deg_df$Number[1:2] - single_ko_deg_df$Number[3:4]
pdf(paste(main_fig_dir, 'fig6e.1.1.pdf', sep = ''), width = 1, height = 1)
ggplot(aes(Type, Number), data = single_ko_deg_df) + geom_bar(aes(fill = subclass), stat = 'identity') + xlab('DEG test type') + ylab('DEGs') + nm_theme() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

###### create two bar plots with the chip-seq data:
double_ko_deg_df <- data.frame(Type = rep(c('Monocyte (double KO vs WT)', 'Monocyte (double KO vs Irf8 KO)',
                                            'Granuocyte (double KO vs WT)', 'Granuocyte (double KO vs Gfi1 KO)'), time = 2),
                               Number = c(double_ko_deg_df[[1]],
                                          double_ko_deg_df[[2]],
                                          double_ko_deg_df[[3]],
                                          double_ko_deg_df[[4]],
                                          length(intersect(double_ko_MONOCYTEvsWT, c(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),
                                          length(intersect(double_ko_MONOCYTEvsIrf8_ko, Gfi1_targets_gene_short_name)),
                                          length(intersect(double_ko_GRANUOCYTEvsWT, c(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),
                                          length(intersect(double_ko_GRANUOCYTEvsGfi1_ko, Irf8_targets_gene_short_name))
                                          ),
                               subclass = rep(c('ko', 'chip'), each = 4))

double_ko_deg_df$Number[1:4] <- double_ko_deg_df$Number[1:4] - double_ko_deg_df$Number[5:8]
pdf(paste(main_fig_dir, 'fig6e.2.1.pdf', sep = ''), width = 2.5, height = 2)
ggplot(aes(Type, Number), data = double_ko_deg_df) + geom_bar(aes(fill = subclass), stat = 'identity') + xlab('DEG test type') + ylab('Number of DEGs') + nm_theme() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pdf(paste(main_fig_dir, 'fig6e.2.1_helper.pdf', sep = ''))
ggplot(aes(Type, Number), data = single_ko_deg_df) + geom_bar(aes(fill = subclass), stat = 'identity') + xlab('DEG test type') + ylab('Number of DEGs') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

########################################################################################################################################################################################################
# downsampling the correctly positioned cells for deg test: (Irf8)
########################################################################################################################################################################################################
rng <- range(URMM_all_abs[, pData(URMM_all_abs)$State == granulocyte_state & pData(URMM_all_abs)$Type %in% c('Irf8_knockout')]$Pseudotime)
# [1] 13.80390 17.10217
Irf8_ko_leaking_cells <- URMM_all_abs[, pData(URMM_all_abs)$State == granulocyte_state & pData(URMM_all_abs)$Type %in% c('Irf8_knockout')]

Irf8_ko_mainfold_cells.1.1 <- URMM_all_abs[, pData(URMM_all_abs)$Pseudotime > floor(rng[1]) & pData(URMM_all_abs)$Pseudotime < ceiling(rng[2]) & pData(URMM_all_abs)$State == granulocyte_state & !(pData(URMM_all_abs)$Type %in% c('Gfi1_knockout', 'Gfi1_Irf8_knockout'))]
pData(Irf8_ko_mainfold_cells.1.1)$fig6_ko_test1.1 <- 'WildType'
pData(Irf8_ko_mainfold_cells.1.1)$fig6_ko_test1.1[pData(Irf8_ko_mainfold_cells.1.1)$Type %in% c('Irf8_knockout')] <- "Irf8 KO"
Irf8_ko_mainfold_cells.1.1 <- Irf8_ko_mainfold_cells.1.1[, c(sample(colnames(Irf8_ko_mainfold_cells.1.1)[pData(Irf8_ko_mainfold_cells.1.1)$fig6_ko_test1.1 == 'Irf8 KO'], 7),
                                                                                    colnames(Irf8_ko_mainfold_cells.1.1)[pData(Irf8_ko_mainfold_cells.1.1)$fig6_ko_test1.1 != 'Irf8 KO'])]
URMM_all_abs_genotype_fig6_ko_test1_beam_genes.1.1 <- differentialGeneTest(Irf8_ko_mainfold_cells.1.1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3) + fig6_ko_test1.1", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)")
URMM_all_abs_genotype_fig6_ko_test1_group_deg.1.1 <- differentialGeneTest(Irf8_ko_mainfold_cells.1.1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~fig6_ko_test1.1", reducedModelFormulaStr = "~1")

# > sum( URMM_all_abs_genotype_fig6_ko_test1_beam_genes.1.1$qval < 0.1)
# [1] 8
# > sum( URMM_all_abs_genotype_fig6_ko_test1_group_deg.1.1$qval < 0.1)
# [1] 12

########################################################################################################################################################################################################
# downsampling the correctly positioned cells for deg test: (Gfi1)
########################################################################################################################################################################################################
rng <- range(URMM_all_abs[, pData(URMM_all_abs)$State == monocyte_state & pData(URMM_all_abs)$Type %in% c('Gfi1_knockout')]$Pseudotime)
# [1] 13.75441 16.13564
Gif1_ko_leaking_cells.1.1 <- URMM_all_abs[, pData(URMM_all_abs)$State == monocyte_state & pData(URMM_all_abs)$Type %in% c('Gfi1_knockout')]

Gif1_ko_leaking_cells.1.1 <- URMM_all_abs[, pData(URMM_all_abs)$Pseudotime > floor(rng[1]) & pData(URMM_all_abs)$Pseudotime < ceiling(rng[2]) & pData(URMM_all_abs)$State == monocyte_state & !(pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_Irf8_knockout'))]
pData(Gif1_ko_leaking_cells.1.1)$fig6_ko_test2.1 <- 'WildType'
pData(Gif1_ko_leaking_cells.1.1)$fig6_ko_test2.1[pData(Gif1_ko_leaking_cells.1.1)$Type %in% c('Gfi1_knockout')] <- "Gfi1 KO"
Gif1_ko_leaking_cells.1.1 <- Gif1_ko_leaking_cells.1.1[, c(sample(colnames(Gif1_ko_leaking_cells.1.1)[pData(Gif1_ko_leaking_cells.1.1)$fig6_ko_test2.1 == 'Gfi1 KO'], 6),
                                                             colnames(Gif1_ko_leaking_cells.1.1)[pData(Gif1_ko_leaking_cells.1.1)$fig6_ko_test2.1 != 'Gfi1 KO'])]
URMM_all_abs_genotype_fig6_ko_test2_beam_genes.1.1 <- differentialGeneTest(Gif1_ko_leaking_cells.1.1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3) + fig6_ko_test2.1", reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)")
URMM_all_abs_genotype_fig6_ko_test2_group_deg.1.1 <- differentialGeneTest(Gif1_ko_leaking_cells.1.1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~fig6_ko_test2.1", reducedModelFormulaStr = "~1")

# > sum( URMM_all_abs_genotype_fig6_ko_test2_beam_genes.1.1$qval < 0.1)
# [1] 1
# > sum(URMM_all_abs_genotype_fig6_ko_test2_group_deg.1.1$qval < 0.1)
# [1] 1

########################################################################################################################################################################################################
# downsampling the correctly positioned cells for deg test: (comparing cells from two branches)
########################################################################################################################################################################################################
all_state23_cells <- URMM_all_abs[, pData(URMM_all_abs)$State %in% c(monocyte_state, granulocyte_state)]
all_state23_cells_ko_deg <- differentialGeneTest(all_state23_cells, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~State", reducedModelFormulaStr = "~1")

sum(all_state23_cells_ko_deg$qval < 0.1)

## where did we get the number 27 and 16 cells?
all_state23_WT_cells <- URMM_all_abs[, pData(URMM_all_abs)$State %in% c(monocyte_state, granulocyte_state) & !(pData(URMM_all_abs)$Type %in% c('Irf8_knockout', 'Gfi1_knockout', 'Gfi1_Irf8_knockout'))]
sample_cells <- c(sample(colnames(all_state23_WT_cells)[pData(all_state23_WT_cells)$State == monocyte_state], 27),
                  sample(colnames(all_state23_WT_cells)[pData(all_state23_WT_cells)$State == granulocyte_state], 16))
all_state23_WT_cells <- all_state23_WT_cells[, sample_cells]
sample_state23_cells_ko_deg <- differentialGeneTest(all_state23_WT_cells, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~State", reducedModelFormulaStr = "~1")

# > sum(sample_state23_cells_ko_deg$qval < 0.1)
# [1] 70

all_state23_cells_WT_ko_deg <- differentialGeneTest(all_state23_WT_cells, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~State", reducedModelFormulaStr = "~1")

#######################################################################################################################################################################################################
# show the deg vs targets (Gfi1 correct positioned cells)
########################################################################################################################################################################################################
Gfi1_ko_proper_positioned_cells <- intersect(beam_genes, row.names(subset(URMM_all_abs_genotype_fig6_ko_test2_group_deg.1, qval < 0.1)))
# Gfi1_ko_proper_positioned_cells <- row.names(subset(URMM_all_abs_genotype_fig6_ko_test2_group_deg.1, qval < 0.1))

mean_exprs_Gfi1_ko_proper_positioned_cells <- rowMeans(exprs(Gif1_ko_leaking_cells.1[Gfi1_ko_proper_positioned_cells, ]))
fig6_ko_test2.1 <- pData(Gif1_ko_leaking_cells.1)$fig6_ko_test2.1

lfc_Gfi1_ko_proper_positioned_cells <- log2(rowMeans(exprs(Gif1_ko_leaking_cells.1[Gfi1_ko_proper_positioned_cells, fig6_ko_test2.1 == "Gfi1 KO"]))/
                                              rowMeans(exprs(Gif1_ko_leaking_cells.1[Gfi1_ko_proper_positioned_cells, fig6_ko_test2.1 == "WildType"])))

lfc_Gfi1_ko_proper_positioned_cells_Irf8_targets <- intersect(Gfi1_ko_proper_positioned_cells, Irf8_targets_gene_short_name)
lfc_Gfi1_ko_proper_positioned_cells_Gfi1_targets <- intersect(Gfi1_ko_proper_positioned_cells, Gfi1_targets_gene_short_name)
Gfi1_Type <- rep(NA, length(Gfi1_ko_proper_positioned_cells))
# Gfi1_Type[Gfi1_ko_proper_positioned_cells %in% lfc_Gfi1_ko_proper_positioned_cells_Irf8_targets] <- 'Irf8'
Gfi1_Type[Gfi1_ko_proper_positioned_cells %in% lfc_Gfi1_ko_proper_positioned_cells_Gfi1_targets] <- 'Gfi1'

########################################################################################################################################################################################################
Gfi1_ko_monocle_df <- data.frame(Type = Gfi1_Type[!is.na(Gfi1_Type)],
                                 mean_expression = mean_exprs_Gfi1_ko_proper_positioned_cells[!is.na(Gfi1_Type)],
                                 lfc = lfc_Gfi1_ko_proper_positioned_cells[!is.na(Gfi1_Type)],
                                 row.names = lfc_Gfi1_ko_proper_positioned_cells_Gfi1_targets)

pdf(paste(main_fig_dir, "fig6_Gfi1_ko_example.1.pdf", sep = ''), width = 1.2, height = 1.5)
ggplot(aes(mean_expression, lfc, label = rownames(Gfi1_ko_monocle_df)), data = Gfi1_ko_monocle_df) + geom_point(size = 1) + nm_theme() + scale_size() #+ geom_text(size=2)
dev.off()
ggplot(aes(mean_expression, lfc, label = rownames(Gfi1_ko_monocle_df)), data = Gfi1_ko_monocle_df) + geom_point(size = 1) + geom_text(size=6)

#show one gene in WT:
pdf(paste(main_fig_dir, "fig6_Gfi1_leaking_example.2.pdf", sep = ''), width = 1.2, height = 1.5)
plot_genes_branched_pseudotime(URMM_all_fig1b[c("Lcn2", "Mreg"), ], min_expr = 0.1) + nm_theme()
dev.off()

#show the same gene as scatter plot:
pData(Gif1_ko_leaking_cells.1)$fig6_ko_test2.1[pData(Gif1_ko_leaking_cells.1)$fig6_ko_test2.1 == 'leaking_cells'] <- 'Gfi1 KO'

pdf(paste(main_fig_dir, "fig6_Gfi1_leaking_example.3.pdf", sep = ''), width = 1.2, height = 1.5)
plot_genes_jitter(Gif1_ko_leaking_cells.1[c("Lcn2", "Mreg"), ], grouping = "fig6_ko_test2.1") + nm_theme()+ xlab('')
dev.off()

plot_genes_jitter(Gif1_ko_leaking_cells.1[c("Lcn2", "Mreg"), ], grouping = "fig6_ko_test2.1") + nm_theme()

########################################################################################################################################################################################################
# show the deg vs targets (Irf8 correct positioned cells)
########################################################################################################################################################################################################
Irf8_ko_proper_positioned_cells <- intersect(beam_genes, row.names(subset(URMM_all_abs_genotype_fig6_ko_test1_group_deg.1, qval < 0.1)))
# Irf8_ko_proper_positioned_cells <- row.names(subset(URMM_all_abs_genotype_fig6_ko_test1_group_deg.1, qval < 0.1))

mean_exprs_Irf8_ko_proper_positioned_cells <- rowMeans(exprs(Irf8_ko_mainfold_cells.1[Irf8_ko_proper_positioned_cells, ]))
fig6_ko_test1.1 <- pData(Irf8_ko_mainfold_cells.1)$fig6_ko_test1.1

lfc_Irf8_ko_proper_positioned_cells <- log2(rowMeans(exprs(Irf8_ko_mainfold_cells.1[Irf8_ko_proper_positioned_cells, fig6_ko_test1.1 == "Irf8 KO"]))/
                                              rowMeans(exprs(Irf8_ko_mainfold_cells.1[Irf8_ko_proper_positioned_cells, fig6_ko_test1.1 == "WildType"])))

lfc_Irf8_ko_proper_positioned_cells_Irf8_targets <- intersect(Irf8_ko_proper_positioned_cells, Irf8_targets_gene_short_name)
lfc_Irf8_ko_proper_positioned_cells_Gfi1_targets <- intersect(Irf8_ko_proper_positioned_cells, Gfi1_targets_gene_short_name)
Irf8_Type <- rep(NA, length(Irf8_ko_proper_positioned_cells))
Irf8_Type[Irf8_ko_proper_positioned_cells %in% lfc_Irf8_ko_proper_positioned_cells_Irf8_targets] <- 'Irf8'
# Irf8_Type[Irf8_ko_proper_positioned_cells %in% lfc_Irf8_ko_proper_positioned_cells_Gfi1_targets] <- 'Gfi1'

qplot(mean_exprs_Irf8_ko_proper_positioned_cells, lfc_Irf8_ko_proper_positioned_cells, color = Irf8_Type) + xlab('Mean expression') + ylab('Fold changes')

# URMM_all_abs_genotype_fig6_ko_test2_group_deg.1 <- differentialGeneTest(Gif1_ko_leaking_cells.1, verbose = T, cores = detectCores() - 2, fullModelFormulaStr = "~fig6_ko_test2.1", reducedModelFormulaStr = "~1")
########################################################################################################################################################################################################
lfc_Irf8_ko_proper_positioned_cells_Irf8_targets
# [1] "Acpp"   "Ccr2"   "Cfp"    "Chn2"   "Dok2"   "Ltf"    "Ly86"   "Pilra"  "Tcf7l2" "Tespa1" "Trps1"
Irf8_ko_monocle_df <- data.frame(Type = Irf8_Type[!is.na(Irf8_Type)],
                                 mean_expression = mean_exprs_Irf8_ko_proper_positioned_cells[!is.na(Irf8_Type)],
                                 lfc = lfc_Irf8_ko_proper_positioned_cells[!is.na(Irf8_Type)],
                                 row.names = lfc_Irf8_ko_proper_positioned_cells_Irf8_targets)

pdf(paste(main_fig_dir, "fig6_Irf8_ko_example.1.pdf", sep = ''), width = 1.2, height = 1.5)
ggplot(aes(mean_expression, lfc, label = rownames(Irf8_ko_monocle_df)), data = Irf8_ko_monocle_df) + geom_point(size = 1) + nm_theme() + scale_size() #+ geom_text(size=2)
dev.off()

ggplot(aes(mean_expression, lfc, label = rownames(Irf8_ko_monocle_df)), data = Irf8_ko_monocle_df) + geom_point(size = 1) + geom_text(size=6)

#show one gene in WT:
pdf(paste(main_fig_dir, "fig6_Irf8_leaking_example.2.pdf", sep = ''), width = 1.2, height = 1.5)
plot_genes_branched_pseudotime(URMM_all_fig1b[c("Ccr2", "Tcf7l2"), ], min_expr = 0.1) + nm_theme()
dev.off()

#show the same gene as scatter plot:
pdf(paste(main_fig_dir, "fig6_Irf8_leaking_example.3.pdf", sep = ''), width = 1.2, height = 1.5)
plot_genes_jitter(Irf8_ko_mainfold_cells.1[c("Ccr2", "Tcf7l2"), ], grouping = "fig6_ko_test1.1") + nm_theme()+ xlab('')
dev.off()

########################################################################################################################################################################################################
# look at the double KO genes
########################################################################################################################################################################################################
# test on the granulocyte branch:
double_ko_granulocyte_branch_genes <- intersect(beam_genes, row.names(subset(URMM_all_abs_genotype_fig6_ko_test3_group_deg, qval < 0.1)))
# double_ko_granulocyte_branch_genes <- row.names(subset(URMM_all_abs_genotype_fig6_ko_test3_group_deg, qval < 0.1))

mean_exprs_double_ko_granulocyte <- rowMeans(exprs(double_ko_leaking_cells_state2[double_ko_granulocyte_branch_genes, ]))
fig6_ko_test3 <- pData(double_ko_leaking_cells_state2)$fig6_ko_test3

lfc_double_ko_granulocyte <- log2(rowMeans(exprs(double_ko_leaking_cells_state2[double_ko_granulocyte_branch_genes, fig6_ko_test3 == "leaking_cells"]))/
                                              rowMeans(exprs(double_ko_leaking_cells_state2[double_ko_granulocyte_branch_genes, fig6_ko_test3 == "WildType"])))

lfc_double_ko_granulocyte_Irf8_targets <- intersect(double_ko_granulocyte_branch_genes, Irf8_targets_gene_short_name)
lfc_double_ko_granulocyte_Gfi1_targets <- intersect(double_ko_granulocyte_branch_genes, Gfi1_targets_gene_short_name)
double_ko_Type <- rep(NA, length(double_ko_granulocyte_branch_genes))
double_ko_Type[double_ko_granulocyte_branch_genes %in% lfc_double_ko_granulocyte_Irf8_targets] <- 'Irf8'
#double_ko_Type[double_ko_granulocyte_branch_genes %in% lfc_double_ko_granulocyte_Gfi1_targets] <- 'Gfi1'

qplot(mean_exprs_double_ko_granulocyte, lfc_double_ko_granulocyte, color = double_ko_Type) + xlab('Mean expression') + ylab('Fold changes')

plot_genes_branched_heatmap(URMM_all_fig1b[intersect(Gfi1_ko_proper_positioned_cells, lfc_Gfi1_ko_proper_positioned_cells_Gfi1_targets), ])

########################################################################################################################################################################################################
# test on the monocyte branch:
double_ko_monocyte_branch_genes <- intersect(beam_genes, row.names(subset(URMM_all_abs_genotype_fig6_ko_test4_group_deg, qval < 0.1)))
# double_ko_monocyte_branch_genes <- row.names(subset(URMM_all_abs_genotype_fig6_ko_test4_group_deg, qval < 0.1))

mean_exprs_double_ko_monocyte <- rowMeans(exprs(double_ko_leaking_cells_state3[double_ko_monocyte_branch_genes, ]))
fig6_ko_test4 <- pData(double_ko_leaking_cells_state3)$fig6_ko_test4

lfc_double_ko_monocyte <- log2(rowMeans(exprs(double_ko_leaking_cells_state3[double_ko_monocyte_branch_genes, fig6_ko_test4 == "leaking_cells"]))/
                                    rowMeans(exprs(double_ko_leaking_cells_state3[double_ko_monocyte_branch_genes, fig6_ko_test4 == "WildType"])))

lfc_double_ko_monocyte_Irf8_targets <- intersect(double_ko_monocyte_branch_genes, Irf8_targets_gene_short_name)
lfc_double_ko_monocyte_Gfi1_targets <- intersect(double_ko_monocyte_branch_genes, Gfi1_targets_gene_short_name)
double_ko_Type <- rep(NA, length(double_ko_monocyte_branch_genes))

# double_ko_Type[double_ko_monocyte_branch_genes %in% lfc_double_ko_monocyte_Irf8_targets] <- 'Irf8'
double_ko_Type[double_ko_monocyte_branch_genes %in% lfc_double_ko_monocyte_Gfi1_targets] <- 'Gfi1'

qplot(mean_exprs_double_ko_monocyte, lfc_double_ko_monocyte, color = double_ko_Type) + xlab('Mean expression') + ylab('Fold changes')

plot_genes_branched_heatmap(URMM_all_fig1b[intersect(Gfi1_ko_proper_positioned_cells, lfc_Gfi1_ko_proper_positioned_cells_Gfi1_targets), ])

plot_genes_branched_heatmap(URMM_all_fig1b[intersect(Gfi1_ko_proper_positioned_cells, lfc_Gfi1_ko_proper_positioned_cells_Gfi1_targets), ])

########################################################################################################################################################################################################
# save the data
########################################################################################################################################################################################################
# save.image('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/RData/fig6_tmp.RData')
save.image('/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision//RData/fig6_tmp.RData')

########################################################################################################################################################################################################
# draw the scatter plot with the x-axis as BEAM gene fold-change and the y-axis as the WT vs KO genes  fold-change genes
########################################################################################################################################################################################################
beam_genes = row.names(subset(fig1b_beam_genes_proj_dup, qval < 1e-1))
########################################################################################################################################################################################################
# fold-change for the Gfi1 KO cells on monocyte branch:
########################################################################################################################################################################################################
monocyte_beam_manifold_genes <- beam_genes #intersect(beam_genes, row.names(subset(URMM_all_abs_genotype_fig6_ko_test2_group_deg.1, qval < 0.1)))  
# fold-change for the WT vs KO genes:
# Gfi1_ko_fold_change <- log2(rowMeans(exprs(Gif1_ko_leaking_cells.1[monocyte_beam_manifold_genes, fig6_ko_test2.1 == "Gfi1 KO"]))/
#                                               rowMeans(exprs(Gif1_ko_leaking_cells.1[monocyte_beam_manifold_genes, fig6_ko_test2.1 == "WildType"])))
Gfi1_ko_fold_change <- log2(rowMeans(exprs(URMM_all_abs[monocyte_beam_manifold_genes, pData(URMM_all_abs)$State == monocyte_state &  pData(URMM_all_abs)$Genotype == 'Gfi1_knockout']))/
                                  rowMeans(exprs(URMM_all_abs[monocyte_beam_manifold_genes, pData(URMM_all_abs)$State == monocyte_state &  pData(URMM_all_abs)$Genotype == 'WildType'])))

gran_vs_mon_fold_change <- log2(rowMeans(exprs(URMM_all_fig1b[monocyte_beam_manifold_genes, pData(URMM_all_fig1b)$State == fig1_granulocyte_state]))/
                                              rowMeans(exprs(URMM_all_fig1b[monocyte_beam_manifold_genes, pData(URMM_all_fig1b)$State == fig1_monocyte_state])))

# create the plot:
Gfi1_Targets <- monocyte_beam_manifold_genes %in% row.names(subset(URMM_all_abs_genotype_fig6_ko_test2_group_deg.1, qval < 0.1)) # monocyte_beam_manifold_genes %in% Gfi1_targets_gene_short_name

Gfi1_fc_df <- data.frame('Granulocyte Vs. Monocyte (WT tree)' = gran_vs_mon_fold_change, 'Gfi1 KO Vs. WT cells (monocyte)' = Gfi1_ko_fold_change, 'Gfi1 targets' = Gfi1_Targets)

pdf(paste(main_fig_dir, "fig6_scatterplot_gfi1_ko.pdf", sep = ''), width = 1.5, height = 1.5)
qplot(Granulocyte.Vs..Monocyte..WT.tree., Gfi1.KO.Vs..WT.cells..monocyte., color = Gfi1.targets, data = Gfi1_fc_df, size = 0.5) + #geom_abline() +
  xlab('Gran Vs. Mon (WT tree)') + ylab('Gfi1 KO Vs. WT (Monocyte)') + scale_size(range = c(0.5, 0.5)) +
  geom_vline(xintercept = 0, linetype = "longdash") + geom_hline(yintercept = 0, linetype = "longdash") + nm_theme() + xlim(c(-10, 10)) + ylim(c(-10, 10))
dev.off()

# check whether or not they are the first / secondary targets
Gfi1_Targets %in% df$labels

########################################################################################################################################################################################################
# fold-change for the Irf8 KO cells on granulocyte branch:
########################################################################################################################################################################################################
#Irf8_ko_proper_positioned_cells
granulocyte_beam_manifold_genes <- beam_genes #intersect(beam_genes, row.names(subset(URMM_all_abs_genotype_fig6_ko_test1_group_deg.1, qval < 0.1)))

# fold-change for the WT vs KO genes:

# mainfold test fold-change
# Irf8_ko_fold_change <- log2(rowMeans(exprs(Irf8_ko_mainfold_cells.1[granulocyte_beam_manifold_genes, fig6_ko_test1.1 == "Irf8 KO"]))/
#                               rowMeans(exprs(Irf8_ko_mainfold_cells.1[granulocyte_beam_manifold_genes, fig6_ko_test1.1 == "WildType"])))
Irf8_ko_fold_change <- log2(rowMeans(exprs(URMM_all_abs[granulocyte_beam_manifold_genes, pData(URMM_all_abs)$State == granulocyte_state &  pData(URMM_all_abs)$Genotype == 'Irf8_knockout']))/
       rowMeans(exprs(URMM_all_abs[granulocyte_beam_manifold_genes, pData(URMM_all_abs)$State == granulocyte_state &  pData(URMM_all_abs)$Genotype == 'WildType'])))

mon_vs_gran_fold_change <- log2(rowMeans(exprs(URMM_all_fig1b[granulocyte_beam_manifold_genes, pData(URMM_all_fig1b)$State == fig1_monocyte_state]))/
                                  rowMeans(exprs(URMM_all_fig1b[granulocyte_beam_manifold_genes, pData(URMM_all_fig1b)$State == fig1_granulocyte_state])))

# create the plot:
Irf8_Targets <- monocyte_beam_manifold_genes %in% row.names(subset(URMM_all_abs_genotype_fig6_ko_test1_group_deg.1, qval < 0.1)) # granulocyte_beam_manifold_genes %in% Irf8_targets_gene_short_name

Irf8_Targets %in% as.character(df$labels)

pdf(paste(main_fig_dir, "fig6_scatterplot_irf8_ko.pdf", sep = ''), width = 1.5, height = 1.5)
qplot(mon_vs_gran_fold_change, Irf8_ko_fold_change, color = Irf8_Targets, size = 0.5) +
  xlab('Mon Vs. Gran (WT tree)') + ylab('Irf8 KO Vs. WT (Granulocyte)') + scale_size(range = c(0.5, 0.5)) + #geom_abline() +
  geom_vline(xintercept = 0, linetype = "longdash") + geom_hline(yintercept = 0, linetype = "longdash") + nm_theme() + xlim(c(-10, 10)) + ylim(c(-10, 10))
dev.off()

pdf(paste(main_fig_dir, "fig6_scatterplot_irf8_ko_helper.pdf", sep = ''))
qplot(mon_vs_gran_fold_change, Irf8_ko_fold_change, color = Irf8_Targets, size = 0.5) +
  xlab('Mon Vs. Gran (WT tree)') + ylab('Irf8 KO Vs. WT (Granulocyte)') + scale_size(range = c(0.5, 0.5)) + #geom_abline() +
  geom_vline(xintercept = 0, linetype = "longdash") + geom_hline(yintercept = 0, linetype = "longdash") # + nm_theme() + xlim(c(-10, 10)) + ylim(c(-10, 10))
dev.off()

########################################################################################################################################################################################################
# fold-change for the double KO cells on granulocyte branch: (aberrant monocyte genes)
########################################################################################################################################################################################################
granulocyte_beam_manifold_genes <- beam_genes #intersect(beam_genes, row.names(subset(URMM_all_abs_genotype_fig6_ko_test4_group_deg, qval < 0.1)))

# fold-change for the WT vs KO genes:
fig6_ko_test4 <- pData(double_ko_leaking_cells_state3)$fig6_ko_test4
# Irf8_ko_fold_change <- log2(rowMeans(exprs(double_ko_leaking_cells_state3[granulocyte_beam_manifold_genes, fig6_ko_test4 == "Double KO"]))/
#                               rowMeans(exprs(double_ko_leaking_cells_state3[granulocyte_beam_manifold_genes, fig6_ko_test4 == "WildType"])))

granulocyte_double_ko_fold_change <- log2(rowMeans(exprs(URMM_all_abs[granulocyte_beam_manifold_genes, pData(URMM_all_abs)$State == granulocyte_state &  pData(URMM_all_abs)$Genotype == 'Gfi1_Irf8_knockout']))/
                                         rowMeans(exprs(URMM_all_abs[granulocyte_beam_manifold_genes, pData(URMM_all_abs)$State == granulocyte_state &  pData(URMM_all_abs)$Genotype == 'WildType'])))

gran_mon_vs_fold_change <- log2(rowMeans(exprs(URMM_all_fig1b[granulocyte_beam_manifold_genes, pData(URMM_all_fig1b)$State == fig1_monocyte_state]))/
                                  rowMeans(exprs(URMM_all_fig1b[granulocyte_beam_manifold_genes, pData(URMM_all_fig1b)$State == fig1_granulocyte_state])))
# create the plot:
Irf8_Targets <- granulocyte_beam_manifold_genes %in% row.names(subset(URMM_all_abs_genotype_fig6_ko_test4_group_deg, qval < 0.1)) #granulocyte_beam_manifold_genes %in% Irf8_targets_gene_short_name

Irf8_Targets %in% as.character(df$labels)

pdf(paste(main_fig_dir, "fig6_scatterplot_double_ko_granulocyte.pdf", sep = ''), width = 1.5, height = 1.5)
qplot(gran_mon_vs_fold_change, granulocyte_double_ko_fold_change, color = Irf8_Targets, size = 0.5) +
  xlab('Mon Vs. Gran (WT tree)') + ylab('Irf8 KO Vs. WT (Granulocyte)') + scale_size(range = c(0.5, 0.5)) + #geom_abline() +
  geom_vline(xintercept = 0, linetype = "longdash") + geom_hline(yintercept = 0, linetype = "longdash") + nm_theme()
dev.off()

########################################################################################################################################################################################################
# fold-change for the double KO cells on monocyte branch:
########################################################################################################################################################################################################

monocyte_double_ko_beam_manifold_genes <- beam_genes #intersect(beam_genes, row.names(subset(URMM_all_abs_genotype_fig6_ko_test3_group_deg, qval < 0.1)))

# fold-change for the WT vs KO genes:
# monocyte_double_ko_fold_change <- log2(rowMeans(exprs(double_ko_leaking_cells_state2[monocyte_double_ko_beam_manifold_genes, fig6_ko_test3 == "double KO"]))/
#                                          rowMeans(exprs(double_ko_leaking_cells_state2[monocyte_double_ko_beam_manifold_genes, fig6_ko_test3 == "WildType"])))

monocyte_double_ko_fold_change <- log2(rowMeans(exprs(URMM_all_abs[monocyte_double_ko_beam_manifold_genes, pData(URMM_all_abs)$State == monocyte_state &  pData(URMM_all_abs)$Genotype == 'Gfi1_Irf8_knockout']))/
                                         rowMeans(exprs(URMM_all_abs[monocyte_double_ko_beam_manifold_genes, pData(URMM_all_abs)$State == monocyte_state &  pData(URMM_all_abs)$Genotype == 'WildType'])))

mon_gran_vs_fold_change <- log2(rowMeans(exprs(URMM_all_fig1b[monocyte_double_ko_beam_manifold_genes, pData(URMM_all_fig1b)$State == fig1_granulocyte_state]))/ #granulocyte
                                  rowMeans(exprs(URMM_all_fig1b[monocyte_double_ko_beam_manifold_genes, pData(URMM_all_fig1b)$State == fig1_monocyte_state]))) #monocyte
# create the plot:
Gfi1_Targets <- monocyte_double_ko_beam_manifold_genes %in% row.names(subset(URMM_all_abs_genotype_fig6_ko_test3_group_deg, qval < 0.1)) #Gfi1_Targets <- monocyte_double_ko_beam_manifold_genes %in% c(Gfi1_targets_gene_short_name) #Irf8_targets_gene_short_name

Gfi1_Targets %in% as.character(df$labels)

pdf(paste(main_fig_dir, "fig6_scatterplot_double_ko_monocyte.pdf", sep = ''), width = 1.5, height = 1.5)
qplot(mon_gran_vs_fold_change, monocyte_double_ko_fold_change, color = Gfi1_Targets, size = 0.5) +
  xlab('Gran Vs. Mon (WT tree)') + ylab('Double KO Vs. WT (Monocyte)') + scale_size(range = c(0.5, 0.5)) + #geom_abline() +
  geom_vline(xintercept = 0, linetype = "longdash") + geom_hline(yintercept = 0, linetype = "longdash") + nm_theme()
dev.off()

########################################################################################################################################################################################################
# facet by different branch
########################################################################################################################################################################################################
single_ko_deg_df2 <- data.frame(Type = rep(c('Irf8_ko', 'Gfi1_ko'), time = 2),
                               Number = c(Irf8_ko_deg_df$proper_positioned_cells,
                                          Gfi1_ko = Gfi1_ko_deg_df$proper_positioned_cells,
                                          length(intersect(MONOCYTE_proper_positioned_cells, Irf8_targets_gene_short_name)),
                                          length(intersect(GRANUOCYTE_proper_positioned_cells, Gfi1_targets_gene_short_name))),
                               subclass = rep(c('ko', 'chip'), each = 2))
single_ko_deg_df2$Number[1:2] <- single_ko_deg_df2$Number[1:2] - single_ko_deg_df2$Number[3:4]
#pdf('main_figures//fig6e.1.1.pdf', width = 1, height = 1)
ggplot(aes(Type, Number), data = single_ko_deg_df2) + geom_bar(aes(fill = subclass), stat = 'identity') + xlab('DEG test type') + ylab('DEGs') + nm_theme() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
#dev.off()

###### create two bar plots with the chip-seq data:
double_ko_deg_vec <- data.frame(state2vsWT = sum(URMM_all_abs_genotype_fig6_ko_test3_group_deg$qval < 0.1),
                               state2vsIrf8_ko = sum(URMM_all_abs_genotype_fig6_ko_test3.1_group_deg$qval < 0.1),
                               state3vsWT = sum(URMM_all_abs_genotype_fig6_ko_test4_group_deg$qval < 0.1),
                               state3vsGfi1_ko = sum(URMM_all_abs_genotype_fig6_ko_test4.1_group_deg$qval < 0.1),
                               cells_between_two_branches = sum(URMM_all_abs_genotype_fig6_double_ko_deg$qval < 0.1))

double_ko_deg_df2 <- data.frame(Type = rep(c('Monocyte (double KO vs WT)', 'Monocyte (double KO vs Irf8 KO)',
                                            'Granuocyte (double KO vs WT)', 'Granuocyte (double KO vs Gfi1 KO)'), time = 2),
                               Number = c(double_ko_deg_vec[[1]],
                                          double_ko_deg_vec[[2]],
                                          double_ko_deg_vec[[3]],
                                          double_ko_deg_vec[[4]],
                                          length(intersect(double_ko_MONOCYTEvsWT, c(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),
                                          length(intersect(double_ko_MONOCYTEvsIrf8_ko, Gfi1_targets_gene_short_name)),
                                          length(intersect(double_ko_GRANUOCYTEvsWT, c(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),
                                          length(intersect(double_ko_GRANUOCYTEvsGfi1_ko, Irf8_targets_gene_short_name))
                               ),
                               subclass = rep(c('ko', 'chip'), each = 4))

double_ko_deg_df2$Number[1:4] <- double_ko_deg_df2$Number[1:4] - double_ko_deg_df2$Number[5:8]
# pdf('main_figures//fig6e.2.1.pdf', width = 2.5, height = 2)
ggplot(aes(Type, Number), data = double_ko_deg_df2) + geom_bar(aes(fill = subclass), stat = 'identity') + xlab('DEG test type') + ylab('Number of DEGs') + nm_theme() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
# dev.off()

# combine them into a single data.frame to do the faceting:
single_ko_deg_df2 <- data.frame(Type = rep(c('Irf8_ko', 'Gfi1_ko'), each = 4),
                                Number = c(Irf8_ko_deg_df$proper_positioned_cells,
                                           0, #length(intersect(MONOCYTE_proper_positioned_cells, setdiff(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),
                                           length(intersect(MONOCYTE_proper_positioned_cells, setdiff(Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name))),
                                           length(intersect(MONOCYTE_proper_positioned_cells, intersect(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),

                                           Gfi1_ko_deg_df$proper_positioned_cells,
                                           length(intersect(GRANUOCYTE_proper_positioned_cells, setdiff(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),
                                           0,#length(intersect(GRANUOCYTE_proper_positioned_cells, setdiff(Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name))),
                                           length(intersect(GRANUOCYTE_proper_positioned_cells, intersect(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name)))
                                           ),
								branch = rep(c('Granulocyte', 'Monocyte'), each = 4),
                                subclass = rep(c('ko', 'Irf8 only', 'Gfi1 only', 'all'), time = 2))

double_ko_deg_df2 <- data.frame(Type = rep(c('Monocyte (double KO vs WT)', #'Monocyte (double KO vs Irf8 KO)',
                                             'Granulocyte (double KO vs WT)'), #, 'Granulocyte (double KO vs Gfi1 KO)'),
                                           each = 4),
                                Number = c(double_ko_deg_vec[[1]],
                                           length(intersect(double_ko_MONOCYTEvsWT, setdiff(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),
                                           0, #length(intersect(double_ko_MONOCYTEvsWT, setdiff(Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name))),
                                           length(intersect(double_ko_MONOCYTEvsWT, intersect(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),

                                           # double_ko_deg_vec[[2]],
                                           # length(intersect(double_ko_MONOCYTEvsIrf8_ko, setdiff(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),
                                           # length(intersect(double_ko_MONOCYTEvsIrf8_ko, setdiff(Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name))),
                                           # length(intersect(double_ko_MONOCYTEvsIrf8_ko, intersect(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),

                                           double_ko_deg_vec[[3]],
                                           0, #length(intersect(double_ko_GRANUOCYTEvsWT, setdiff(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),
                                           length(intersect(double_ko_GRANUOCYTEvsWT, setdiff(Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name))),
                                           length(intersect(double_ko_GRANUOCYTEvsWT, intersect(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name)))#,

                                           # double_ko_deg_vec[[4]],
                                           # length(intersect(double_ko_GRANUOCYTEvsGfi1_ko, setdiff(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),
                                           # length(intersect(double_ko_GRANUOCYTEvsGfi1_ko, setdiff(Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name))),
                                           # length(intersect(double_ko_GRANUOCYTEvsGfi1_ko, intersect(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))
                                          ),
								branch = rep(c('Monocyte', 'Granulocyte'), each = 4), #8
                                subclass = rep(c('ko', 'Irf8 only', 'Gfi1 only', 'all'), time = 2))

single_double_ko_df <- rbind(single_ko_deg_df2, double_ko_deg_df2)

single_double_ko_df$Number[c(1, 1:3 * 4 + 1)] <- single_double_ko_df$Number[c(1, 1:3 * 4 + 1)] -
													c(sum(single_double_ko_df$Number[2:4]),
													sum(single_double_ko_df$Number[6:8]),
													sum(single_double_ko_df$Number[10:12]),
													sum(single_double_ko_df$Number[14:16]))
pdf(paste(main_fig_dir, '/fig6e.gm_facet.pdf', sep = ''), width = 2.5, height = 2)
ggplot(aes(Type, Number), data = single_double_ko_df) + geom_bar(aes(fill = subclass), stat = 'identity') + xlab('DEG test type') +
ylab('Number of DEGs')  + facet_wrap(~branch, scale = 'free') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + nm_theme()
dev.off()

pdf(paste(main_fig_dir, '/fig6e.gm_facet_helper.pdf', sep = ''), width = 2.5, height = 2)
ggplot(aes(Type, Number), data = single_double_ko_df) + geom_bar(aes(fill = subclass), stat = 'identity') + xlab('DEG test type') +
  ylab('Number of DEGs')  + facet_wrap(~branch, scale = 'free') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) #+ nm_theme()
dev.off()

single_double_ko_df$fraction <- single_double_ko_df$Number
single_double_ko_df$fraction[1:4] <-  single_double_ko_df$Number[1:4] / sum(single_double_ko_df$Number[1:4])
single_double_ko_df$fraction[5:8] <-  single_double_ko_df$Number[5:8] / sum(single_double_ko_df$Number[5:8])
single_double_ko_df$fraction[9:12] <-  single_double_ko_df$Number[9:12] / sum(single_double_ko_df$Number[9:12])
single_double_ko_df$fraction[13:16] <-  single_double_ko_df$Number[13:16] / sum(single_double_ko_df$Number[13:16])

ggplot(aes(Type, fraction), data = single_double_ko_df) + geom_bar(aes(fill = subclass), stat = 'identity') + xlab('DEG test type') +
  ylab('Fraction dof DEGs')  + facet_wrap(~branch, scale = 'free') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) #+ nm_theme()

########################################################################################################################################################################################################
# fraction of enr
########################################################################################################################################################################################################
# combine them into a single data.frame to do the faceting:

granulocyte_exprs_gene <- row.names(URMM_all_abs)[esApply(URMM_all_abs[, pData(URMM_all_abs)$State == 2], 1, function(x) sum(x > 0.1) > 5 )]
monocyte_exprs_gene <- row.names(URMM_all_abs)[esApply(URMM_all_abs[, pData(URMM_all_abs)$State == 3], 1, function(x) sum(x > 0.1) > 5 )]

bg_single_ko_deg_df2 <- data.frame(Type = rep(c('Irf8_ko', 'Gfi1_ko'), each = 4),
                                Number = c(length(granulocyte_exprs_gene),
                                           0, #length(intersect(MONOCYTE_proper_positioned_cells, setdiff(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),
                                           length(intersect(granulocyte_exprs_gene, setdiff(Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name))),
                                           length(intersect(granulocyte_exprs_gene, intersect(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),

                                           length(monocyte_exprs_gene),
                                           length(intersect(monocyte_exprs_gene, setdiff(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),
                                           0,#length(intersect(GRANUOCYTE_proper_positioned_cells, setdiff(Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name))),
                                           length(intersect(monocyte_exprs_gene, intersect(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name)))
                                ),
                                branch = rep(c('Granulocyte', 'Monocyte'), each = 4),
                                subclass = rep(c('ko', 'Irf8 only', 'Gfi1 only', 'all'), time = 2))

bg_double_ko_deg_df2 <- data.frame(Type = rep(c('Monocyte (double KO vs WT)', #'Monocyte (double KO vs Irf8 KO)',
                                             'Granulocyte (double KO vs WT)'), #, 'Granulocyte (double KO vs Gfi1 KO)'),
                                           each = 4),
                                Number = c(length(granulocyte_exprs_gene),
                                           length(intersect(granulocyte_exprs_gene, setdiff(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),
                                           0, #length(intersect(double_ko_MONOCYTEvsWT, setdiff(Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name))),
                                           length(intersect(granulocyte_exprs_gene, intersect(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),

                                           # double_ko_deg_vec[[2]],
                                           # length(intersect(double_ko_MONOCYTEvsIrf8_ko, setdiff(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),
                                           # length(intersect(double_ko_MONOCYTEvsIrf8_ko, setdiff(Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name))),
                                           # length(intersect(double_ko_MONOCYTEvsIrf8_ko, intersect(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),

                                           length(monocyte_exprs_gene),
                                           0, #length(intersect(double_ko_GRANUOCYTEvsWT, setdiff(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),
                                           length(intersect(monocyte_exprs_gene, setdiff(Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name))),
                                           length(intersect(monocyte_exprs_gene, intersect(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name)))#,

                                           # double_ko_deg_vec[[4]],
                                           # length(intersect(double_ko_GRANUOCYTEvsGfi1_ko, setdiff(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))),
                                           # length(intersect(double_ko_GRANUOCYTEvsGfi1_ko, setdiff(Gfi1_targets_gene_short_name, Irf8_targets_gene_short_name))),
                                           # length(intersect(double_ko_GRANUOCYTEvsGfi1_ko, intersect(Irf8_targets_gene_short_name, Gfi1_targets_gene_short_name))
                                ),
                                branch = rep(c('Monocyte', 'Granulocyte'), each = 4), #8
                                subclass = rep(c('ko', 'Irf8 only', 'Gfi1 only', 'all'), time = 2))

bg_single_double_ko_df <- rbind(bg_single_ko_deg_df2, bg_double_ko_deg_df2)

bg_single_double_ko_df$Number[c(1, 1:3 * 4 + 1)] <- bg_single_double_ko_df$Number[c(1, 1:3 * 4 + 1)] -
  c(sum(bg_single_double_ko_df$Number[2:4]),
    sum(bg_single_double_ko_df$Number[6:8]),
    sum(bg_single_double_ko_df$Number[10:12]),
    sum(bg_single_double_ko_df$Number[14:16]))

bg_single_double_ko_df$fraction <- bg_single_double_ko_df$Number
bg_single_double_ko_df$fraction[1:4] <-  bg_single_double_ko_df$Number[1:4] / sum(bg_single_double_ko_df$Number[1:4])
bg_single_double_ko_df$fraction[5:8] <-  bg_single_double_ko_df$Number[5:8] / sum(bg_single_double_ko_df$Number[5:8])
bg_single_double_ko_df$fraction[9:12] <-  bg_single_double_ko_df$Number[9:12] / sum(bg_single_double_ko_df$Number[9:12])
bg_single_double_ko_df$fraction[13:16] <-  bg_single_double_ko_df$Number[13:16] / sum(bg_single_double_ko_df$Number[13:16])

ggplot(aes(Type, fraction), data = bg_single_double_ko_df) + geom_bar(aes(fill = subclass), stat = 'identity') + xlab('DEG test type') +
  ylab('Fraction dof DEGs')  + facet_wrap(~branch, scale = 'free') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) #+ nm_theme()

bg_test_single_double_ko_df <- rbind(single_double_ko_df, bg_single_double_ko_df)
bg_test_single_double_ko_df$background <- c(rep('NO', nrow(single_double_ko_df)), rep('YES', nrow(bg_single_double_ko_df)) )

ggplot(aes(Type, fraction), data = bg_test_single_double_ko_df) + geom_bar(aes(fill = subclass), stat = 'identity') + xlab('DEG test type') +
  ylab('Fraction dof DEGs')  + facet_wrap(~branch + background, scale = 'free') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) #+ nm_theme()

########################################################################################################################################################################################################
# save the data
########################################################################################################################################################################################################
save.image('./RData/fig6.RData')
