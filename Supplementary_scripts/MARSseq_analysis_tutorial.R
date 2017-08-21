# packages = c("ggplot2", "VGAM", "igraph", "pRlyr", "combinat", "fastICA", "irlba", "matrixStats", "reshape2", "stringr", "modeest", "Hmisc", "dplyr", "mixsmsn","ggplot2", "VGAM", "igraph", "pRlyr", "combinat", "fastICA", "irlba", "matrixStats", "reshape2", "stringr","modeest","Hmisc", "DESeq", "DESeq2", "doMC", "data.table", "fitdistrplus", "ggdendro", "gplots", "lmtest", "mixsmsn", "monocle", "pheatmap", "piano", "pscl", "VennDiagram", "venneuler", "zoo")
# install.packages(packages)
#              
# bio_packages = c("Biobase", "BiocGenerics",  "limma", "ComplexHeatmap", "HSMMSingleCell", 'GEOquery')
# source("http://bioconductor.org/biocLite.R")
# biocLite(c('DESeq', 'DESeq2', 'piano', 'edgeR',  'MAST'))
# biocLite(bio_packages)
#              
#install.packages(c('scatterplot3d', 'colorRamps', 'TTR', 'ComplexHeatmap', 'circlize', 'Rtsne', 'FNN', 'reshape', 'reshape2', 'qlcMatrix', 'proxy', 'slam', 'R.matlab', 'flexclust', 'mcclust))

# install.packages(c('raster', 'princurve', 'snow', 'lle', 'tsne', 'roxygen2', 'testthat'))
#R CMD install DDRTree_0.1.5.tar.gz
# library("devtools")
# install_github("jw156605/SLICER")
# install.packages('/Users/xqiu/Downloads/dpt_0.6.0 (2).tar.gz')

# biocLite("ComplexHeatmap")
# install.packages('ComplexHeatmap')
# install_github('Xiaojieqiu/densityClust')
# library(ggplot2)
# library(gplots)
# library(colorRamps)
# library(RColorBrewer)
# library(scatterplot3d)
# library(TTR)
# library(ComplexHeatmap)
# library(circlize)
# library(Rtsne)
# library(destiny)
# library(reshape)
# library(pheatmap)
# 
# options(jupyter.plot_mimetypes = 'image/png')
# 
# 
# load('../RData/Paul_Cell_MARSseq_GSE72857.RData')

ls()

dim(data.debatched)

table(cluster.id)

rownames(data.debatched)[1:10]
info.genes[1:10]

gene.names <-sapply(strsplit(rownames(data.debatched), ";"), "[", 1)
is.informative <- gene.names %in% info.genes[order(info.genes)]
data.info.genes <- data.debatched[is.informative,]
rownames(data.info.genes) <- gene.names[is.informative]

dim(data.info.genes)

info.genes[!info.genes %in% gene.names[is.informative]]

sum(rownames(data.debatched) %in% info.genes)

quantile.cell <- quantile(colSums(data), 0.2)

quantile.cell

total.counts <- colSums(data)
total.genes <- apply(data, 2, function(x) {sum(x>0)})
mean.counts  <- rowMeans(data)
var.counts <- apply(data, 1, var)
cv2.counts <- apply(data, 1, function(x) {var(x)/mean(x)^2})
drop.out <- apply(data, 1, function(x) {sum(x==0)/length(x)})

addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

color.pal <- colorRampPalette(c('blue','red'))(20) 
color.pal2 <- addalpha(colorRampPalette(rev(brewer.pal(8, 'RdYlBu')))(20), alpha=0.3)
transparent.grey <- rgb(0.2,0.2,0.2, 0.3)


par(mfrow=c(2,2)) #smaller plots, 2 per row, 2 per column
plot(total.genes, total.counts, 
     pch=20,col=transparent.grey,
     xlab= 'Number of detected genes',
     ylab= 'Number of UMI filtered reads',
     main='Number of detected genes vs\n Number of Reads')
plot(mean.counts, var.counts, log='xy', 
     pch=20,col=transparent.grey,
     xlab= 'Mean expression',
     ylab= 'Variance in expression',
     main='Variance to mean')
plot(mean.counts, cv2.counts, log='xy', 
     pch=20, col=transparent.grey,
     xlab= 'Mean expression',
     ylab= 'Squared coefficient of variation',
     main='Mean expression vs\n Squared coefficient of variation')
plot(mean.counts, drop.out, log='x', 
     pch=20, col=transparent.grey,
     xlab= 'Mean expression',
     ylab= 'Drop out rate',
     main='Mean expression vs\n Drop out rate')

plot(cv2.counts, drop.out, log='x', 
     pch=20,col=transparent.grey,
     xlab= 'Squared coefficient of variation',
     ylab= 'Drop out rate',
     main='Squared coefficient of variation vs\n Drop out rate')



boxplot(total.genes ~ batch.names,  col=color.pal, main='Number of genes by Batch')
boxplot(total.counts ~ batch.names,  col=color.pal, main='Number of UMI filtered\n reads by Batch')

total.count.debatch <- colSums(data.debatched)
total.genes.debatch <- apply(data.debatched, 2, function(x) {sum(x>0)})
mean.counts.debatch  <- rowMeans(data.debatched)
cv2.counts.debatch <- apply(data.debatched, 1, function(x) {var(x)/mean(x)^2})
var.counts.debatch <- apply(data.debatched, 1, var)
drop.out.debatch <- apply(data.debatched, 1, function(x) {sum(x==0)/length(x)})

par(mfrow=c(2,2)) #smaller plots, 2 per row, 2 per column
plot(total.genes.debatch, total.count.debatch, 
     pch=20,col=transparent.grey,
     xlab= 'Number of detected genes',
     ylab= 'Number of UMI filtered reads',
     main='Number of detected genes vs\n Number of Reads')
plot(mean.counts.debatch, var.counts.debatch, log='xy', 
     pch=20,col=transparent.grey,
     xlab= 'Mean expression',
     ylab= 'Variance in expression',
     main='Variance to mean')
plot(mean.counts.debatch, cv2.counts.debatch, log='xy', 
     pch=20, col=transparent.grey,
     xlab= 'Mean expression',
     ylab= 'Coefficient of variation',
     main='Mean expression vs\n Coefficient of Variation')
points(mean.counts.debatch[is.informative], cv2.counts.debatch[is.informative], col='gold', pch='.')
plot(mean.counts.debatch, drop.out.debatch, log='x', 
     pch=20, col=transparent.grey,
     xlab= 'Mean expression',
     ylab= 'Drop out rate',
     main='Mean expression vs\n Drop out rate')

plot(cv2.counts.debatch, drop.out.debatch, log='x', 
     pch=20,col=transparent.grey,
     xlab= 'Coefficient of variation',
     ylab= 'Drop out rate',
     main='Coefficient of Variation vs\n Drop out rate')


boxplot(total.genes.debatch ~ batch.names,  col=color.pal, main='Number of genes by Batch')
boxplot(total.count.debatch ~ batch.names,  col=color.pal, main='Number of UMI filtered\n reads by Batch')



par(mfrow=c(2,2))
plot(mean.counts , mean.counts.debatch, pch=20,
     col = transparent.grey,
     xlab = 'processed', 
     ylab = 'debatched',
     main = 'Mean expression')
plot(cv2.counts , cv2.counts.debatch, pch=20,
     col= transparent.grey,
     xlab = 'processed', 
     ylab = 'debatched',
     main = 'Squared coefficient of variation')
plot(drop.out, drop.out.debatch, pch=20,
     col = transparent.grey,
     xlab = 'processed', 
     ylab = 'debatched',
     main = 'Drop out rate')

total.count.info <- colSums(data.info.genes)
total.genes.info <- apply(data.info.genes, 2, function(x) {sum(x>0)})
mean.counts.info  <- rowMeans(data.info.genes)
cv2.counts.info <- apply(data.info.genes, 1, function(x) {var(x)/mean(x)^2})
var.counts.info <- apply(data.info.genes, 1, var)
drop.out.info <- apply(data.info.genes, 1, function(x) {sum(x==0)/length(x)})

par(mfrow=c(2,2)) #smaller plots, 2 per row, 2 per column
plot(total.genes.info, total.count.info, 
     pch=20,col=transparent.grey,
     xlab= 'Number of detected genes',
     ylab= 'Number of UMI filtered reads',
     main='Number of detected genes vs\n Number of Reads')
plot(mean.counts.info, var.counts.info, log='xy', 
     pch=20,col=transparent.grey,
     xlab= 'Mean expression',
     ylab= 'Variance in expression',
     main='Variance to mean')
plot(mean.counts.info, cv2.counts.info, log='xy', 
     pch=20, col=transparent.grey,
     xlab= 'Mean expression',
     ylab= 'Coefficient of variation',
     main='Mean expression vs\n Coefficient of Variation')
plot(mean.counts.info, drop.out.info, log='x', 
     pch=20, col=transparent.grey,
     xlab= 'Mean expression',
     ylab= 'Drop out rate',
     main='Mean expression vs\n Drop out rate')

plot(cv2.counts.info, drop.out.info, log='x', 
     pch=20,col=transparent.grey,
     xlab= 'Coefficient of variation',
     ylab= 'Drop out rate',
     main='Coefficient of Variation vs\n Drop out rate')


boxplot(total.genes.info ~ batch.names,  col=color.pal, main='Number of genes by Batch')
boxplot(total.count.info ~ batch.names,  col=color.pal, main='Number of UMI filtered\n reads by Batch')


pca.data <- prcomp(t(data.info.genes), center = TRUE, scale=TRUE)

str(pca.data)

normal <- sum(pca.data$sdev^2)
var <- round((pca.data$sdev)^2 / normal *100,1)
ord <- sample(dim(pca.data$x)[1])

par(mfrow=c(2,2))
plot(pca.data$x[ord,1], pca.data$x[ord,2], 
     col=rgb(0.2,0.2,0.2,0.3), 
     xlab = paste0("PC1, ", var[1], " % variance"),
     ylab = paste0("PC2, ", var[2], " % variance"),
     pch= 20)

plot(pca.data$x[ord,1], pca.data$x[ord,2], 
     col=color.pal2[(1+data.info.genes['Hbb-b1',ord])/quantile(data.info.genes['Hbb-b1',], 0.98)*20], 
     xlab = paste0("PC1, ", var[1], " % variance"),
     ylab = paste0("PC2, ", var[2], " % variance"),
     main='Hbb-b1 expression',
     pch= 20)

plot(pca.data$x[ord,1], pca.data$x[ord,2], 
     col=color.pal2[(1+data.info.genes['Sfpi1',ord])/quantile(data.info.genes['Sfpi1',], 0.98)*20], 
     xlab = paste0("PC1, ", var[1], " % variance"),
     ylab = paste0("PC2, ", var[2], " % variance"),
     main='Sfpi1 expression',
     pch= 20)

plot(pca.data$x[ord,1], pca.data$x[ord,2], 
     col=color.pal2[(1+data.info.genes['Apoe',ord])/quantile(data.info.genes['Apoe',], 0.98)*20], 
     xlab = paste0("PC1, ", var[1], " % variance"),
     ylab = paste0("PC2, ", var[2], " % variance"),
     main='Apoe expression',
     pch= 20)

tsne.done <- Rtsne(t(data.info.genes), perplexity = 40, max_iter = 2000)


tsne.wo.pca <- Rtsne(t(data.info.genes), perplexity = 40, max_iter = 4000, pca = FALSE)

str(tsne.done)

par(mfrow=c(2,2))
plot(tsne.done$Y, col=transparent.grey, 
     xlab= 't-SNE 1',
     ylab= 't-SNE 2',
     pch=20)
plot(tsne.done$Y[ord,], 
     col=color.pal2[(1+data.info.genes['Hbb-b1',ord])/quantile(data.info.genes['Hbb-b1',], 0.98)*20], 
     pch=20, 
     xlab= 't-SNE 1',
     ylab= 't-SNE 2',
     main='Hbb-b1 expression')
plot(tsne.done$Y[ord,], 
     col=color.pal2[(1+data.info.genes['Sfpi1',ord])/quantile(data.info.genes['Sfpi1',], 0.98)*20],  
     pch=20, 
     xlab= 't-SNE 1',
     ylab= 't-SNE 2',
     main='Sfpi1 expression')
plot(tsne.done$Y[ord,], 
     col=color.pal2[(1+data.info.genes['Apoe',ord])/quantile(data.info.genes['Apoe',], 0.98)*20],  
     pch=20, 
     xlab= 't-SNE 1',
     ylab= 't-SNE 2',
     main='Apoe expression')   

plot(tsne.wo.pca$Y, col=transparent.grey, 
     xlab= 't-SNE 1',
     ylab= 't-SNE 2',
     pch=20)
plot(tsne.wo.pca$Y[ord,], 
     col=color.pal2[(1+data.info.genes['Hbb-b1',ord])/quantile(data.info.genes['Hbb-b1',], 0.98)*20], 
     pch=20, 
     xlab= 't-SNE 1',
     ylab= 't-SNE 2',
     main='Hbb-b1 expression')
plot(tsne.wo.pca$Y[ord,], 
     col=color.pal2[(1+data.info.genes['Sfpi1',ord])/quantile(data.info.genes['Sfpi1',], 0.98)*20], 
     pch=20, 
     xlab= 't-SNE 1',
     ylab= 't-SNE 2',
     main='Sfpi1 expression')
plot(tsne.wo.pca$Y[ord,], 
     col=color.pal2[(1+data.info.genes['Apoe',ord])/quantile(data.info.genes['Apoe',], 0.98)*20], 
     pch=20, 
     xlab= 't-SNE 1',
     ylab= 't-SNE 2',
     main='Apoe expression')  

#extra legend (use the grDevices package)
legend_image <- as.raster(matrix(color.pal2, ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'normalised\ngene expression')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, 1, 1, 0,0)


diff.plot <- DiffusionMap(log(t(data.info.genes)+1), k=20, sigma='local')


dim(diff.plot@eigenvectors)


par(mfrow=c(2,2))
plot(diff.plot, 1:2, 
     col=transparent.grey, 
     pch=20, 
     axes=TRUE)

plot(diff.plot, 1:3, angle=60,
     col=transparent.grey, 
     pch=20, 
     axes=TRUE)

plot(diff.plot, 1:2,
     col = color.pal2[(1+data.info.genes['Hbb-b1',])/quantile(data.info.genes['Hbb-b1',], 0.98)*20],
     pch=20, 
     main='Hbb-b1',
     axes=TRUE)

plot(diff.plot, 1:2,
     col = color.pal2[(1+data.info.genes['Sfpi1',])/quantile(data.info.genes['Sfpi1',], 0.98)*20],
     pch=20, 
     main='Sfpi1',
     axes=TRUE)
plot(diff.plot, 1:2,
     col = color.pal2[(1+data.info.genes['Apoe',])/quantile(data.info.genes['Apoe',], 0.98)*20],
     pch=20, 
     main = 'Apoe',
     axes=TRUE)
#extra legend (use the grDevices package)
legend_image <- as.raster(matrix(color.pal2, ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'normalised\ngene expression')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, 1, 1, 0,0)


dpt <- DPT(diff.plot, tips = c(841, 1240, 878))


color.branch <- addalpha(brewer.pal(8,'Paired'), 0.3)

plot(dpt, 1:2, 
     col_by = 'dpt',  root = 1, paths_to = c(2,3),
     pch=20)
plot(dpt, 1:2, 
     col_by = 'Branch',  root = 1, paths_to = c(2,3), pal=color.branch,
     pch=20)

head(dpt@branch)

table(dpt@branch[,1])
sum(is.na(dpt@branch[,1]))


data.info.zoom <- data.info.genes[,dpt@branch[,1] %in% c(NA, 1)] 


dim(data.info.zoom)


new.start <- which(colnames(data.info.zoom)==colnames(data.info.genes)[841])


diff.plot.zoom <- DiffusionMap(log(t(data.info.zoom)+1), k=20, sigma='local', n_local = 10)


dpt.zoom <- DPT(diff.plot.zoom, tips = new.start) 


plot(dpt.zoom, col_by = 'branch',  pch = 20,dcs=c(1:3),angle=135, pal=color.branch[c(1,5,6)])



table(dpt.zoom@branch[,1])


branching <- numeric(dim(data.info.genes)[2])
branching[!is.na(dpt@branch[,1])] <- dpt@branch[!is.na(dpt@branch[,1]),1]
branching[dpt@branch[,1] %in% c(NA,1)][dpt.zoom@branch[,1]==1] <- 1


dpt.zoom@branch[dpt.zoom@branch[,1] %in% 2:3,1] <- dpt.zoom@branch[dpt.zoom@branch[,1] %in% 2:3,1] + 3
dpt.zoom@branch[is.na(dpt.zoom@branch[,1]),1] <- 4
branching[branching==0] <- dpt.zoom@branch[dpt.zoom@branch[,1] !=1,1]


plot(eigenvectors(diff.plot)[,c(1,2)], col=color.branch[branching], pch=20, xlab='', ylab='', 
     axes=FALSE)
par(cex.lab=1.8)
axis(side=1, at=c(-0.05, 0.1), labels=NA,  col.ticks = 'white')
axis(side=2, at=c(-0.11, 0.11), labels=NA, col.ticks = 'white')


boxplot(dpt$DPT1 ~cluster.id, col=color.pal, cex.axis=0.8)


confusion.matrix <- table(branching,cluster.id)


annotation <- data.frame(cluster = 1:19)

pheatmap(t(t(confusion.matrix) /colSums(confusion.matrix)*100), 
         border_color = 'grey90', annotation_legend = FALSE,
         breaks = c(0,1,5,seq(10,100, by=5)),
         color = colorRampPalette(c('white','black'))(25),         
         cellwidth = 20, cellheight = 20, annotation=annotation, 
         annotation_colors = list(cluster=color.pal),
         cluster_rows = FALSE, cluster_cols = FALSE)


ha_clust = HeatmapAnnotation(df = data.frame(cluster = 1:19), 
                             col = list(cluster = c("1" =color.pal[1],"2" =color.pal[2], "3" =color.pal[3],
                                                    "4" =color.pal[4], "5" =color.pal[5], "6" =color.pal[6],
                                                    "7" =color.pal[7], "8" =color.pal[8], "9" =color.pal[9],
                                                    "10" =color.pal[10],"11" =color.pal[11],"12" =color.pal[12],
                                                    "13" =color.pal[13],"14" =color.pal[14],"15" =color.pal[15],
                                                    "16" =color.pal[16],"17" =color.pal[17],"18" =color.pal[18],
                                                    "19" =color.pal[19])), 
                             width = unit(0.5, "cm"), which="row", show_legend = FALSE)


res <-t(matrix(as.numeric(confusion.matrix), ncol=length(unique(cluster.id)))) /colSums(confusion.matrix)*100


dim(res)

rownames(res) <- 1:length(unique(cluster.id))
colnames(res) <- 1:ncol(res)


ha_column = HeatmapAnnotation(cn = function(index) {
  branch = 1:ncol(res)
  grid.text(branch, (branch-0.5)/length(branch), 1, just = c("center", "top"))
})


ht_cm <- Heatmap(res,       
                 show_column_names=FALSE, show_row_names =TRUE,  
                 column_title= 'Branch', row_title = 'Cluster',
                 column_title_side ='bottom', row_title_side='right',
                 heatmap_legend_param = list(color_bar = "continuous"),
                 cluster_columns = FALSE, cluster_rows =FALSE,
                 name = 'Percentage',#show_heatmap_legend = FALSE,
                 row_names_gp=gpar(las=2),
                 col = colorRampPalette(c('white','black'))(50),
                 bottom_annotation = ha_column
)


ha_clust+ht_cm 
ht_global_opt("heatmap_row_names_gp" = gpar(fontsize = 20), "heatmap_column_names_gp" = gpar(fontsize = 20))

decorate_heatmap_body("Percentage", {
  i = 1:ncol(res)
  x = i/ncol(res)
  for (k in i){
    grid.lines(c(x[k], x[k]), c(0, 1), gp = gpar(lwd = 1, col='grey90'))
  }
  grid.lines(c(0, 0), c(0, 1), gp = gpar(lwd = 1, col='grey90'))
  grid.lines(c(0, 1), c(0, 0),  gp = gpar(lwd = 1, col='grey90'))
  i = 1:nrow(res)
  x = i/nrow(res)
  for (k in i){
    grid.lines(c(0, 1),c(x[k], x[k]),  gp = gpar(lwd = 1, col='grey90'))
  }
  
})


ery.cluster <- cluster.id[branching==2][order(cluster.id[branching==2])]
ery.cluster <- ery.cluster[ery.cluster !=9] #exclude the only cell originating from cluster 9 
gmp.cluster <- cluster.id[branching==3][order(cluster.id[branching==3])]
gmp.cluster <- gmp.cluster[3:length(gmp.cluster)] #exclude two single cells of cluster 3 and 4


par(mfrow=c(2,1))

plot(NA, NA, xlim=c(min(dpt$DPT1),max(dpt$DPT1[branching==2])), ylim=c(0,2),
     cex.lab=1.5, cex.axis=1.5,
     xlab='DPT', ylab='density')
for (i in ery.cluster){
  dens.1 <- density(dpt$DPT1[branching==2 & cluster.id==i], adjust = TRUE)
  lines(dens.1, col=color.pal[i], lwd=2)
}

plot(NA, NA, xlim=c(min(dpt$DPT1),max(dpt$DPT1[branching==3])), ylim=c(0,2),
     cex.lab=1.5, cex.axis=1.5,
     xlab='DPT', ylab='density')
for (i in gmp.cluster){
  dens.1 <- density(dpt$DPT1[branching==3 & cluster.id==i], adjust = TRUE)
  lines(dens.1, col=color.pal[i], lwd=2)
}


plot(pca.data$x[ord,1], pca.data$x[ord,2], 
     col=color.branch[branching[ord]], 
     xlab = paste0("PC1, ", var[1], " % variance"),
     ylab = paste0("PC2, ", var[2], " % variance"),
     pch= 20)
scatterplot3d(pca.data$x[ord,1], pca.data$x[ord,2], -pca.data$x[ord,3], 
              color=color.branch[branching[ord]], angle=230,
              col.axis='grey75', col.grid='white', 
              xlab = paste0("PC1, ", var[1], " % variance"),
              ylab = paste0("PC2, ", var[2], " % variance"),
              zlab = paste0("PC3, ", var[3], " % variance"),
              pch= 20)


plot(tsne.done$Y, col=color.branch[branching], pch=20,
     xlab='t-SNE 1', ylab='t-SNE 2')


favour.genes <- c('Gata1', 'Phf10', 'Zfpm1', 'Gfi1b', 'Cited4', 'Klf1', 'Mbd2', 'E2f4', 'Tcf3', 'Phb2', 'Hmgb3',
                  'Cited2',  'Pbx1', 'Mef2c', 'Coro1a', 'Lyz1', 'Ccnb2', 'H2-DMa', 'Itga2b', 'Gata2', 'Sfpi1',
                  'Lmo4', 'Runx1', 'Cebpe', 'Gfi1', 'Irf8', 'Sfpi1', 'Fcgr3',  'Stat3', 'Etv6',
                  'Cebpa', 'Hbb-b1', 'Hba-a2', 'Car1', 'Car2', 'Apoe', 'Prss34', 'Mpo', 'Pf4',
                  'Serpina3f', 'Cpox', 'Cd34')


which(!favour.genes %in% rownames(data.info.genes))
favour.genes[!favour.genes %in% rownames(data.info.genes)]


order.cells <- order(dpt$DPT1)


par(mfrow=c(2,2))
for (i in which(rownames(data.info.genes)%in%favour.genes))
{
  const=i
  plot(dpt$DPT1[order.cells[branching[order.cells]==3]],
       data.info.genes[const,order.cells[branching[order.cells]==3]], 
       main=rownames(data.info.genes)[const], cex=0.6,
       col=addalpha(color.pal,0.3)
       [cluster.id[order.cells[branching[order.cells]==3]]],
       pch=20, ylab='Gene expression (UMI counts)', xlab='DPT')
}



pie(rep(1,length(unique(cluster.id))), col=color.pal, 
    labels = 1:length(unique(cluster.id)), main='Color key (Clusters)')


var.log.data <-apply(X = data.info.genes, MARGIN = 1, FUN = var)


cor.data <- cor(t(data.info.genes[var.log.data!=0,]))


dim(cor.data)


dd.data <- hclust(as.dist(1-cor.data), method="ward.D2")


cluster.hierarch <- cutree(dd.data, k=4)


data.smooth <- data.info.genes


perc.98 <- apply(data.smooth[,branching %in% 1:3], MARGIN=1, quantile, c(0.98,1))
perc.98[1,min(perc.98)==0] <- perc.98[2, min(perc.98)==0] 


data.smooth <- data.smooth/perc.98[1,]


data.smooth.branch1 <- na.omit(apply(log(data.smooth[,order.cells[branching[order.cells]==1]]+0.01),1 , SMA, 20))
data.smooth.branch2 <- na.omit(apply(log(data.smooth[,order.cells[branching[order.cells]==2]]+0.01),1 , SMA, 20))
data.smooth.branch3 <- na.omit(apply(log(data.smooth[,order.cells[branching[order.cells]==3]]+0.01),1 , SMA, 20))


ha.1 = HeatmapAnnotation(
  df = data.frame(cluster=cluster.id[order.cells[branching[order.cells]==1][20:sum(branching==1)]]), 
  col=list(cluster=c("1" =color.pal[1],"2" =color.pal[2], "3" =color.pal[3],
                     "4" =color.pal[4], "5" =color.pal[5], "6" =color.pal[6],
                     "7" =color.pal[7], "8" =color.pal[8], "9" =color.pal[9],
                     "10" =color.pal[10],"11" =color.pal[11],"12" =color.pal[12],
                     "13" =color.pal[13],"14" =color.pal[14],"15" =color.pal[15],
                     "16" =color.pal[16],"17" =color.pal[17],"18" =color.pal[18],
                     "19" =color.pal[19]
  )), 
  show_legend = FALSE, 
  gap = unit(4, "mm"))
ha.2 = HeatmapAnnotation(
  df = data.frame(type=cluster.id[order.cells[branching[order.cells]==2]][20:sum(branching==2)]), 
  col=list(type=c("1" =color.pal[1],"2" =color.pal[2], "3" =color.pal[3],
                  "4" =color.pal[4], "5" =color.pal[5], "6" =color.pal[6],
                  "7" =color.pal[7], "8" =color.pal[8], "9" =color.pal[9],
                  "10" =color.pal[10],"11" =color.pal[11],"12" =color.pal[12],
                  "13" =color.pal[13],"14" =color.pal[14],"15" =color.pal[15],
                  "16" =color.pal[16],"17" =color.pal[17],"18" =color.pal[18],
                  "19" =color.pal[19]
  )), 
  show_legend = FALSE,
  gap = unit(4, "mm"))
ha.3 = HeatmapAnnotation(
  df = data.frame(type=cluster.id[order.cells[branching[order.cells]==3]][20:sum(branching==3)]),
  col=list(type=c("1" =color.pal[1],"2" =color.pal[2], "3" =color.pal[3],
                  "4" =color.pal[4], "5" =color.pal[5], "6" =color.pal[6],
                  "7" =color.pal[7], "8" =color.pal[8], "9" =color.pal[9],
                  "10" =color.pal[10],"11" =color.pal[11],"12" =color.pal[12],
                  "13" =color.pal[13],"14" =color.pal[14],"15" =color.pal[15],
                  "16" =color.pal[16],"17" =color.pal[17],"18" =color.pal[18],
                  "19" =color.pal[19]
  )), 
  show_legend = FALSE, 
  gap = unit(4, "mm"))

cluster.color <- brewer.pal(8, 'Set1')


ha.row = HeatmapAnnotation(df = data.frame(cluster = cluster.hierarch), 
                           col = list(cluster = c("1" = cluster.color[1], 
                                                  "2" = cluster.color[2], 
                                                  "3" = cluster.color[3],
                                                  "4" = cluster.color[4])),
                           width = unit(0.5, "cm"), which="row", show_legend = FALSE)

ht1 <- Heatmap( t(data.smooth.branch1), 
                row_dend_width = unit(3, "cm"),              
                show_row_names = FALSE, 
                show_heatmap_legend = FALSE,
                column_title = 'trunk',
                heatmap_legend_param = list(color_bar = "continuous"),
                cluster_columns = FALSE, 
                cluster_rows =as.dendrogram(dd.data),
                name = 'expression',
                col=colorRamp2(seq(-4.6,-0.2, length.out = 13), c('grey75', rev(brewer.pal(11, "RdYlBu")), 'black')),
                top_annotation=ha.1
)
ht2 <- Heatmap( t(data.smooth.branch2), 
                row_dend_width = unit(3, "cm"),              
                show_row_names = FALSE, 
                show_heatmap_legend = FALSE,column_title = 'branch 2',
                heatmap_legend_param = list(color_bar = "continuous"),
                cluster_columns = FALSE, cluster_rows =as.dendrogram(dd.data),
                name = 'expression',
                col=colorRamp2(seq(-4.6,-0.2, length.out = 13), c('grey75', rev(brewer.pal(11, "RdYlBu")), 'black')),
                top_annotation=ha.2
)
ht3 <- Heatmap( t(data.smooth.branch3), 
                row_dend_width = unit(3, "cm"),              
                show_row_names = FALSE, 
                show_heatmap_legend = TRUE,
                column_title = 'branch 3',
                heatmap_legend_param = list(color_bar = "continuous"),
                cluster_columns = FALSE, cluster_rows =as.dendrogram(dd.data),
                name = 'expression',
                col=colorRamp2(seq(-4.6,-0.2, length.out = 13), c('grey75', rev(brewer.pal(11, "RdYlBu")), 'black')),
                top_annotation=ha.3
)

ht1 +ht2 + ht3 + ha.row

cor.data.fav <- cor(t(data.info.genes[var.log.data!=0 & rownames(data.info.genes) %in% favour.genes,]))

dd.data.fav <- hclust(as.dist(1-cor.data.fav), method="ward.D2")

ht1.fav <- Heatmap( t(data.smooth.branch1[,rownames(data.info.genes) %in% favour.genes]), 
                    row_dend_width = unit(3, "cm"),              
                    show_row_names = FALSE, #width=unit(3, "cm"),
                    show_heatmap_legend = FALSE, column_title = 'trunk',
                    heatmap_legend_param = list(color_bar = "continuous"),
                    cluster_columns = FALSE, cluster_rows =as.dendrogram(dd.data.fav),
                    name = 'expression',#show_heatmap_legend = FALSE,
                    col=colorRamp2(seq(-4.6,-0.2, length.out = 13), c('grey75', rev(brewer.pal(11, "RdYlBu")), 'black')),
                    top_annotation=ha.1#, top_annotation_height = unit(2, "cm")
)
ht2.fav <- Heatmap( t(data.smooth.branch2[,rownames(data.info.genes) %in% favour.genes]),  
                    row_dend_width = unit(3, "cm"),              
                    show_row_names = FALSE, #width=unit(3, "cm"),
                    show_heatmap_legend = FALSE,column_title = 'branch 2',
                    heatmap_legend_param = list(color_bar = "continuous"),
                    cluster_columns = FALSE, cluster_rows =as.dendrogram(dd.data.fav),
                    name = 'log(expression)',#show_heatmap_legend = FALSE,
                    col=colorRamp2(seq(-4.6,-0.2, length.out = 13), c('grey75', rev(brewer.pal(11, "RdYlBu")), 'black')),
                    top_annotation=ha.2#, top_annotation_height = unit(2, "cm")
)
ht3.fav <- Heatmap( t(data.smooth.branch3[,rownames(data.info.genes) %in% favour.genes]),  
                    row_dend_width = unit(3, "cm"),              
                    show_row_names = TRUE, #width=unit(3, "cm"),
                    row_names_gp = gpar(fontsize = 10),
                    show_heatmap_legend = TRUE,
                    column_title = 'branch 3',
                    heatmap_legend_param = list(color_bar = "continuous"),
                    cluster_columns = FALSE, cluster_rows =as.dendrogram(dd.data.fav),
                    name = 'expression',#show_heatmap_legend = FALSE,
                    col=colorRamp2(seq(-4.6,-0.2, length.out = 13), c('grey75', rev(brewer.pal(11, "RdYlBu")), 'black')),
                    top_annotation=ha.3#, top_annotation_height = unit(2, "cm")
)

ha.row.fav = HeatmapAnnotation(df = data.frame(cluster = cluster.hierarch
                                               [rownames(data.info.genes) %in% favour.genes]), 
                               col = list(cluster = c("1" =  cluster.color[1], 
                                                      "2" = cluster.color[2], 
                                                      "3"=cluster.color[3],
                                                      "4"=cluster.color[4])), 
                               width = unit(0.5, "cm"), which="row", show_legend = FALSE)

ht1.fav + ht2.fav + ht3.fav + ha.row.fav

sessionInfo()
save.image('MARSseq_analysis_tutorial.RData')