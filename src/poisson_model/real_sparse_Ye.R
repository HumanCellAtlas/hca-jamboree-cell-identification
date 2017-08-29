if (!interactive()) { 
  cmd <- commandArgs(trailing=TRUE)
  in.path <- cmd[1]
  outname <- cmd[2]
}
library(Matrix)
mat <- readMM(file.path(in.path, "Ye2_sparse_molecule_counts.mtx"))

bpath <- file.path(in.path, "Ye2_barcode_id.csv")
b.in <- read.table(bpath, stringsAsFactors=FALSE, sep = ",", row.names = 1)

gpath <- file.path(in.path, "Ye2_gene_id.csv")
g.in <- read.table(gpath, stringsAsFactors=FALSE, sep = ",", row.names = 1)

colnames(mat) <- g.in[,1]
rownames(mat) <- b.in[,1]
saveRDS(file=outname, t(mat))