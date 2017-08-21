library(dpath)

# loading an expression profile matrix of 548 genes and 100 cells
data(etv2)

# run monocle 2: 
etv2_dimnames <- dimnames(etv2)
colnames(etv2) <- paste('cell', 1:ncol(etv2), sep = '_')
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = row.names(etv2), row.names = row.names(etv2)))
pd <- new("AnnotatedDataFrame", data = data.frame(Time = etv2_dimnames[[2]],
                                                  row.names = colnames(etv2)))
etv2_cds <- newCellDataSet(as(etv2, "sparseMatrix"), #he cofactor parameter is used for arcsinh
                                      phenoData = pd, 
                                      featureData = fd, 
                                      expressionFamily=tobit(), 
                                      lowerDetectionLimit=1)
etv2_cds <- reduceDimension(etv2_cds, norm_method = 'log', pseudo_expr = 1, verbose = T)
etv2_cds <- orderCells(etv2_cds)
plot_cell_trajectory(etv2_cds, color_by = 'Time')

cell.group <- factor(colnames(etv2), c('E7.25', 'E7.75', 'E8.25'))

library(dpath)
library(monocle)
################################################################################################################################################################################################
# run on the lung data 
################################################################################################################################################################################################
lung_data <- exprs(absolute_cds)[quake_id, ]
colnames(lung_data) <- pData(absolute_cds)$Time

time.group <- factor(colnames(lung_data), c('E14.5', 'E16.5', 'E18.5', 'Adult')) #"E18.5" "E14.5" "Adult" "E16.5"
# run wp-NMF with 4 metagenes, using cells from E7.75 and E8.25 to initialize
# the factorization, repeat the factorization for 50 times and utilize 8 CPU
# cores
lung_dp <- dpath(lung_data, K = 4, subset.cell = time.group %in% c('E16.5', 'E18.5'),
            repeat.mf = 50, mc.cores = 2)

lung_gene.list <- quake_id[5:10]	
plot(lung_dp, type = 'markers', genes = lung_gene.list, reorder.genes = FALSE)   

lung_dp <- fitsom(lung_dp, xdim = 15, ydim = 15, n.min = 15)
plot(lung_dp, type = 'cell.cluster', cell.group = cell.group)
plot(lung_dp, type = 'cell.cluster', cell.group = cell.group)
plot(lung_dp, type = 'metagene.entropy', cell.group = cell.group)
plot(lung_dp, type = 'gene.expression', genes = c('Runx1', 'Gata1', 'Etv2', 
                                             'Plvap', 'T', 'Tbx20'))
################################################################################################################################################################################################
# run on the HSMM data 
################################################################################################################################################################################################
load('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/RData/fig1.RData')

HSMM_data <- as.matrix(exprs(HSMM_myo)[fData(HSMM_myo)$use_for_ordering, ])
colnames(HSMM_data) <- pData(HSMM_myo)$Time

time.group <- factor(colnames(HSMM_data), c('0', '24', '48', '72')) #"E18.5" "E14.5" "Adult" "E16.5"
# run wp-NMF with 4 metagenes, using cells from E7.75 and E8.25 to initialize
# the factorization, repeat the factorization for 50 times and utilize 8 CPU
# cores
HSMM_dp <- dpath(HSMM_data, K = 4, subset.cell = time.group %in% c('48', '72'),
                 repeat.mf = 50, mc.cores = 2)

HSMM_gene.list <- row.names(HSMM_myo)[fData(HSMM_myo)$use_for_ordering][5:10]	
plot(HSMM_dp, type = 'markers', genes = HSMM_gene.list, reorder.genes = FALSE)   

HSMM_dp <- fitsom(HSMM_dp, xdim = 15, ydim = 15, n.min = 15)
plot(HSMM_dp, type = 'cell.cluster', cell.group = cell.group)
plot(HSMM_dp, type = 'cell.cluster', cell.group = cell.group)
plot(HSMM_dp, type = 'metagene.entropy', cell.group = cell.group)
plot(HSMM_dp, type = 'gene.expression', genes = c('Runx1', 'Gata1', 'Etv2', 
                                                  'Plvap', 'T', 'Tbx20'))
################################################################################################################################################################################################
data(etv2)
cell.group <- factor(colnames(etv2), c('E7.25', 'E7.75', 'E8.25'))
# run wp-NMF with 4 metagenes, using cells from E7.75 and E8.25 to initialize
# the factorization, repeat the factorization for 50 times and utilize 8 CPU
# cores
dp <- dpath(etv2, K = 4, subset.cell = cell.group %in% c('E7.75', 'E8.25'),
repeat.mf = 50, mc.cores = 4)

# reorder the metagenes from different repetitive runs
#dp <- reorder(dp)

# load the finished dp file
data(dp)

# visualizing metagene coefficients, metagene basis and observed expression
# levels for selected marker genes
gene.list <- c('Gata1', 'Ikzf1', 'Itga2b', 'Hba-a1', 'Runx1', 'Gata4', 
               'Smarcd3', 'Tbx20', 'Alcam', 'Cgnl1', 'Dok4', 'Plvap', 'Emcn', 'Pecam1', 
               'Cd34', 'Cdh5', 'Icam1', 'T', 'Kdr', 'Pdgfra', 'Gli2', 'Pou5f1', 'Nanog')	
dev.new(width = 22, height = 15)
plot(dp, type = 'markers', genes = gene.list, reorder.genes = FALSE)

# fitting a self-organizing map (SOM) using bootstrapped cells and clustering
# metacells by partitinoning the metacell landscape
set.seed(6580)
dp <- fitsom(dp, xdim = 15, ydim = 15, n.min = 15)

# plot the cell clustering results and the mean metagene coefficients of each
# cluster
dev.new(width = 9, height = 6)
plot(dp, type = 'cell.cluster', cell.group = cell.group)

# visualizing metagene entropy on the SOM
dev.new(width = 15, height = 15)
plot(dp, type = 'metagene.entropy', cell.group = cell.group)

# visualizing expression pattern of selected genes on the SOM
dev.new(width = 15, height = 10)
par(mfrow = c(2, 3), mar = c(3, 3, 5, 1))
plot(dp, type = 'gene.expression', genes = c('Runx1', 'Gata1', 'Etv2', 
                                             'Plvap', 'T', 'Tbx20'))

# prioritizing metacells with respect to the committed state of the 1st metagene
# (endothelial metagene)
pg.endothelial <- prioritize(dp, metagene = c(1, 0, 0, 0), direction = 'committed')
library(gplots)
col.metacell <- colorpanel(100, low = 'black', mid = 'white', high = 'purple')

# visualizing the prioritization score of committed endothelial lineages (the
# 1st metagene)
dev.new(width = 7, height = 7)
plot(dp, type = 'metacell.landscape', property = pg.endothelial$metacell, 
     col.metacell = col.metacell)

# show the top 100 most enriched genes in the committed endothelial lineage (the 1st
# metagene)
dev.new(width = 22, height = 15)
plot(dp, type = 'prioritization', score.gene = pg.endothelial[['gene']], 
     score.metacell = pg.endothelial [['metacell']], top = 100)

# prioritizing metacells as the progenitor state for the metagenes 1-3
pg.mpc <- prioritize(dp, metagene = c(1, 1, 1, 0), direction = 'progenitor')

# find the differentiation paths from the progenitor toward the committed
# state of metagene 1 - 3
p2c <- differentiation.path.p2c(dp, metagene = c(1, 2, 3))

# visualizing the differentiation paths toward the committed states of metagene
# 1 - 3
dev.new(width = 7, height = 7)
plot(dp, type = 'metacell.landscape', paths = p2c$paths, 
     property = pg.mpc$metacell, col.metacell = col.metacell)
