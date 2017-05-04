library(cellrangerRkit)

gbmAgg_path<-"/Volumes/oncogxA/Projects/CPTRES/RNAExpression/10x/10x_170206/CPTRES_agg"
gbmAgg<-load_cellranger_matrix(gbmAgg_path)

use_genes <- get_nonzero_genes(gbmAgg)
gbmSub = gbmAgg[use_genes,]

geneBC = as.data.frame(as.matrix(exprs(gbmSub)))

gbmAgg_bcnorm <- normalize_barcode_sums_to_median(gbmAgg[use_genes,])
gbmAgg_log <- log_gene_bc_matrix(gbmAgg_bcnorm,base=10)
print(dim(gbmAgg_log))

tsne_proj<-gbmAgg_results$tsne

genes<-c("ANXA3","BASP1","SPP1","SFTA1P","EMP2","DPYSL4", "IGF2","SPDL1","S100A16")
visualize_gene_markers(gbmAgg_log,genes,tsne_proj[c("TSNE.1","TSNE.2")],limits=c(0,1.5))

#source("http://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
#library(rhdf5)
#
#h5_path<-"/Volumes/oncogxA/Projects/CPTRES/RNAExpression/10x/10x_170206/CPTRES_agg/outs/filtered_gene_bc_matrices_h5.h5"
#h5ls(h5_path)
#filteredGeneBC = h5read(h5_path, "/GRCh38/data")
