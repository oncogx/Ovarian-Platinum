function (gbm, n_pcs = 10, logscale = FALSE, ...) 
{
    if (logscale) {
        gbm_log <- gbm
    }
    else {
        use_genes <- get_nonzero_genes(gbm)
        gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes, 
            ])
        gbm_log <- log_gene_bc_matrix(gbm_bcnorm, ...)
    }
    pca <- sparse_pca(t(exprs(gbm_log)), n_pcs)
    pc_names <- sprintf("PC%d", 1:ncol(pca$x))
    rownames(pca$x) <- colnames(gbm)
    colnames(pca$x) <- pc_names
    rownames(pca$rotation) <- rownames(gbm)[use_genes]
    colnames(pca$rotation) <- pc_names
    var_explained <- sum(pca$var_pcs)
    cat("Variance explained by PCs:", var_explained, "\n")
    return(list(x = pca$x, rotation = pca$rotation, sdev = pca$sdev, 
        tot_var = pca$tot_var, var_pcs = pca$var_pcs, use_genes = use_genes, 
        normalized_mat = exprs(gbm_log), var_explained = var_explained))
}
