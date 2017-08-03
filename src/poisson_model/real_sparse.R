if (!interactive()) { 
    cmd <- commandArgs(trailing=TRUE)
    in.path <- cmd[1]
    outname <- cmd[2]
}

library(Matrix)
mat <- readMM(file.path(in.path, "matrix.mtx"))

bpath <- file.path(in.path, "barcodes.tsv")
b.in <- read.table(bpath, stringsAsFactors=FALSE)

gpath <- file.path(in.path, "genes.tsv")
g.in <- read.table(gpath, stringsAsFactors=FALSE)

colnames(mat) <- b.in[,1]
rownames(mat) <- g.in[,1]
saveRDS(file=outname, mat)

