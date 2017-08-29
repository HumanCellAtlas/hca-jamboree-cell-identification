library(EmptyDrops)
m <- readRDS("simu_1.mat.RDS")
out <- detectCells(m, npts=20000, BPPARAM=MulticoreParam(3))

sig <- out$FDR <= 0.01 & !is.na(out$FDR)
table(sig, out$Limited) # Any FALSE sig, TRUE limited means npts was too low.

png("simu_1.png", width=12, height=12, units="in", pointsize=12, res=120)
plot(out$Total, out$LR, log="x", col=ifelse(sig, "red", "black"))
o <- order(out$Total)
lines(out$Total[o], out$Expected[o], col="dodgerblue")
legend("topleft", sprintf("%.2f", sum(sig)/length(sig) * 100), bty="n")
dev.off()

truth <- read.csv("real_1.gz", header=TRUE) # Assessment.
table(truth$barc %in% rownames(out)[sig], truth$mode) # false negatives
summary(!rownames(out)[sig] %in% truth$barc) # false positives
