arg = commandArgs(trailingOnly=TRUE)[1]
if (!is.na(arg)) {
	sim_num=arg
}
calls1 = read.table(paste("~/shared_scratch/group4/betterfit/simu_",sim_num,".tsv", sep=""), header=TRUE)
fdr = calls1[,2];
nNA = sum(is.na(fdr))
fdr[is.na(fdr)] = 1;
calls1 <- calls1[fdr < 0.05,1]
truth = read.table(paste("~/shared_scratch/group4/simdata/real_",sim_num,".gz", sep=""), sep=",", header=TRUE)

truth_1 <- truth[truth[,2]==1,]
truth_2 <- truth[truth[,2]==2,]

TP1 = sum(as.character(calls1) %in% truth_1[,1])
TP2 = sum(as.character(calls1) %in% truth_2[,1])
FP = length(calls1) - TP1 - TP2
FN1 = sum(!(truth_1[,1] %in% calls1))
FN2 = sum(!(truth_2[,1] %in% calls1))
out1 <- data.frame(Group1=c(TP1, FP, FN1, nNA), Group2=c(TP2,FP, FN2, nNA))

out <- cbind(out1)

rownames(out) <- c("TP", "FP", "FN", "nNA");
print(out)

#write.table(as.matrix(calls), file=paste("final_cell_barcodes_sim", sim_num,".txt", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
