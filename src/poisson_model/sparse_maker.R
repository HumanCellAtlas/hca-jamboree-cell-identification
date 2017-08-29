# Loading into a sparse matrix.
library(Matrix)
f=commandArgs(trailingOnly=T)
print(paste0("Reading file ", f))
df=read.table(f, header=T, sep=",")
fg <- factor(df$row)
fc <- factor(df$column)
m <- sparseMatrix(i=as.integer(fg), j=as.integer(fc), x=df$value)
colnames(m) <- levels(fc)
rownames(m) <- levels(fg)
fout=gsub(".gz$", ".mat.RDS", f)
print(paste0("Saving matrix to file ", fout))
saveRDS(m, file=fout)
# Computing the average profile.
