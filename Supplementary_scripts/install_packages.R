#install this particular version of ggplot2 for making the exact figures as in our manuscript: 
install.packages('./ggplot2_1.0.1.tar.gz', dependencies = TRUE)

#install all packages: 
packages = c("ggplot2", "VGAM", "igraph", "pRlyr", "combinat", "fastICA", "irlba", "matrixStats", "reshape2", "R.utils", "snow", "tsne", "lle", "DDRTree",
		"tidyr",  "dplyr", "stringr", "modeest", "Hmisc", "boot", "doMC", "data.table", "fitdistrplus", "ggdendro", "gplots", "princurve", "sp", "devtools", "gridExtra",
        "lmtest", "MASS", "mixsmsn", "pheatmap", "plyr", "pscl", "RColorBrewer", "VennDiagram", "zoo", "raster", "colorRamps", "grid", "ROCR", "lle", "UpSetR")

install.packages(packages, repo = 'http://cran.fhcrc.org/')

# install bioconductor packages 
bio_packages = c("Biobase", "BiocGenerics", "cummeRbund", "limma", "edgeR", "DESeq", "DESeq2", "piano", "MAST", "Destiny")
source("http://bioconductor.org/biocLite.R")
biocLite(bio_packages)

# install packages I wrote for this project locally (Note that L1Graph and SimplePPT will be on CRAN soon, DDRTree and densityClust will updated in 
# CRAN while Monocle 2 will be updated in bioconductor. The newest DDRTree from CRAN is the same as the one we provided here)
install.packages('./DDRTree_0.1.5.tar.gz', dependencies = TRUE)
install.packages('./L1Graph_0.1.0.tar.gz', dependencies = TRUE) 
install.packages('./densityClust_0.3.tar.gz', dependencies = TRUE)
install.packages('./simplePPT_0.1.0.tar.gz', dependencies = TRUE)
install.packages('./xacHelper_0.0.1.0000.tar.gz', dependencies = TRUE)
install.packages('./monocle_2.2.0.tar.gz', dependencies = TRUE)

# install SLICER 
library("devtools")
install_github("jw156605/SLICER")