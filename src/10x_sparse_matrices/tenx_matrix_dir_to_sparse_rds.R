library(argparse)
library(Matrix)

.tenx_to_sparse = function(matrix_dir, genome="GRCh38") {
	  base_path = file.path(matrix_dir, genome)
  base_path.unaggregated = file.path(matrix_dir, genome)

    # Get the file path
    if(file.exists(base_path)) {
	          base_path = base_path
    } else if(file.exists(base_path.unaggregated)) {
	          # This is an aggregated run, change the base path
	          base_path = base_path.unaggregated
      } else {
	            # No expected directories were found
	            stop(paste("Specified genome does not appear in 10X output:", base_path, ' or ', base_path.unaggregated))
        }

    matrix_path = file.path(base_path, "matrix.mtx")
      genes_path = file.path(base_path, "genes.tsv")
      barcodes_path = file.path(base_path, "barcodes.tsv")

        if( ! file.exists(matrix_path) ) { stop(paste("Expression matrix not found in 10X output:", matrix_path)) }
        if( ! file.exists(genes_path) ) { stop(paste("Genes file not found in 10X output:", genes_path)) }
	  if( ! file.exists(barcodes_path) ) { stop(paste("Barcodes file not found in 10X output:", barcodes_path)) }

	  # All files exist, read them in
	  matrix = Matrix::readMM(matrix_path)
	    barcodes = read.table(barcodes_path, header=F, as.is=T)[,1]
	    gene_table = read.table(genes_path, header=F, as.is=T) ## saves for later

	      genes = gene_table[, 1]

	      # Add gene and sample names to expression matrix (adding dataset post-fix in case barcodes appear in multiple samples)
	      row.names(matrix) = genes
	        colnames(matrix) = barcodes
	        matrix = matrix[gene_table[, 1], ] ## ensures order of genes matches

		  row.names(matrix) = gene_table[, 1]

		  return(matrix)
}


parser = argparse::ArgumentParser(description="Script to convert downloaded 10X data to sparse matrix RDS")
parser$add_argument('input_dir', help='Parent directory of matrices to load')
parser$add_argument('sparse_rds', help='Sparse RDS file of matrix.')
parser$add_argument('--genome', default="GRCh38", help='Genome to use.')
args = parser$parse_args()


sparse_matrix = .tenx_to_sparse(args$input_dir, genome=args$genome)
saveRDS(sparse_matrix, file=args$sparse_rds)
