source("model_fitter.R")
m <- readRDS("simu_1.mat.RDS")
out <- detectCells(m, npts=1000)

sig <- out$FDR <= 0.001 & !is.na(out$FDR)
png("simu_1.png", width=12, height=12, units="in", pointsize=12, res=120)
plot(out$Total, out$LR, log="x", col=ifelse(sig, "red", "black"))
o <- order(out$Total)
lines(out$Total[o], out$Expected[o], col="dodgerblue")
legend("topleft", sprintf("%.2f", sum(sig)/length(sig) * 100), bty="n")
dev.off()

