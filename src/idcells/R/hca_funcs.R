library(Matrix)

#' Get proportion of counts for a cell assigned to a gene set
#' 
#' @param object an SCESet object
#' @param geneset character vector of HGNC symbols for the given gene set
#' 
#' @export
#' 
get_geneset_score <- function(object, geneset) {
    set_total <- colSums(counts(object)[fData(object)$hgnc_symbol %in% geneset,])
    set_total / object$total_counts
}

#' Get proportion of counts for a cell assigned to a gene set
#' 
#' @param object an SCESet object
#' @param geneset_list a GSEABase GeneSetCollection object
#' 
#' @export
#' 
add_geneset_scores <- function(object, geneset_list) {
    for (i in seq_along(geneset_list))
        pData(object)[, names(geneset_list)[i]] <- 
            get_geneset_score(object, geneset_list[[1]]@geneIds)
    object
}


#' Function to convert dataframe of (gene,barcode,umi_count) to sparse Matrix
#' @export
df2sparse <- function(d) {
    stopifnot(all(c("gene","barcode","count") %in% colnames(d)))
    g <- as.factor(d$gene)
    b <- as.factor(d$barcode)
    m <- Matrix::sparseMatrix(i=as.numeric(g), j=as.numeric(b), x=d$count)
    colnames(m) <- levels(b)
    rownames(m) <- levels(g)
    m
}

#' Profile ambient genes
#' 
#' Use barcodes with low total UMI counts to obtain aggregate profile
#' of gene expression.
#' 
#' @param object an SCESet object
#' @param max_umis numeric, maximum number of total UMIs to define a "low total"
#' barcode
#' 
#' @return list with aggregated counts and aggregated counts-per-million
#' 
profile_ambient <- function(object, max_umis = 1000) {
    agg_counts <- rowSums(counts(object)[, object$total_counts < max_umis])
    agg_cpm <- agg_counts / sum(agg_counts)
    data_frame(agg_counts, agg_cpm)
}

#' Get barcode stats
#' 
#' @param mat numeric matrix, rows represent genes, columns represent cells
#' 
#' @return a data_frame with barcode-level statistics
#' 
#' @export
get_barcode_stats <- function(mat) {
    dplyr::data_frame(
        total_counts = colSums(mat),
        mean_counts = colMeans(mat),
        var_counts = matrixStats::colVars(mat),
        genes_detected = colSums(mat > 0),
        prop_detected = colMeans(mat > 0)
    )
}

#' Get gene-level stats
#' 
#' @param mat numeric matrix, rows represent genes, columns represent cells
#' 
#' @return a data_frame with gene-level statistic
#' 
#' @export
get_gene_stats <- function(mat) {
    dplyr::data_frame(
        total_counts = rowSums(mat),
        mean_counts = rowMeans(mat),
        var_counts = matrixStats::rowVars(mat),
        cells_detected = rowSums(mat > 0),
        prop_detected = rowMeans(mat > 0)
    )
}

#' Compute Shannon information entropy for a count vector
#' 
#' @param x a vector of counts
#' 
#' @return the entropy of the corresponding multinomial distribution
#' 
#' @export
entropy<-function(x){
  lfrq<-log(x)-log(sum(x))
  sum(exp(lfrq + log(-lfrq)))
}

