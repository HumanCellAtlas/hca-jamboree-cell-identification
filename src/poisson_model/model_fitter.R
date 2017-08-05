if (!interactive()) {
    cmd.opts <- commandArgs(trailingOnly=TRUE)
    inname <- cmd.opts[1]
    outfix <- cmd.opts[2]
}

library(Matrix)
m <- readRDS(inname)
discard <- rowSums(m) == 0
m <- m[!discard,]

# Computing the average profile.

tail.prob <- 0.85
umi.sum <- colSums(m)
ambient <- (umi.sum <= quantile(umi.sum, tail.prob))
ambient.cells <- m[,ambient]
ambient.prof <- rowSums(ambient.cells)

library(edgeR)
ambient.prop <- goodTuringProportions(ambient.prof)

# Removing cells from the ambient.

m <- m[,!ambient]
gc()

# Calculating the likelihood ratio.

library(methods)
LRFUN <- function(observed, prop) { 
   if (is(observed, "dgCMatrix")) {
       i <- observed@i + 1L
       j <- rep(seq_len(ncol(observed)), diff(observed@p)) 
       x <- observed@x
   } else if (is(observed, "dgTMatrix")) {
       i <- observed@i + 1L
       j <- observed@j + 1L
       x <- observed@x
   } else {
       stop("unsupported matrix type")
   }
  
   total <- colSums(observed)
   p.n0 <- x * log(x/(prop[i]*total[j])) - x # Poisson deviance for observed count - deviance for a zero count.
   by.col <- aggregate(p.n0, list(Col=j), sum)
   p.n0.sum <- numeric(length(total))
   p.n0.sum[by.col$Col] <- by.col$x
   p.0.sum <- total # colsum deviance for zero counts for all entries => colsum of means => totals.
   return(p.n0.sum + p.0.sum)
}

obs <- m
obs.totals <- colSums(obs)
obs.LR <- LRFUN(obs, prop=ambient.prop)

# Computing a simulation with 100000 entries for any "tol"-fold interval around the interrogation point..

tol <- 0.5
lower.pt <- log2(min(obs.totals))-tol
upper.pt <- log2(max(obs.totals))+tol
N <- round(1000 * (upper.pt - lower.pt)/(2*tol))
sim.totals <- 2^seq(from=lower.pt, to=upper.pt, length.out=N)

sim.LR <- numeric(N)
for (x in seq_along(sim.LR)) {
    cur.means <- ambient.prop*sim.totals[x]
    current <- rpois(length(ambient.prop), lambda=cur.means)
    sim.LR[x] <- sum(current * log(current/cur.means), na.rm=TRUE) + sum(cur.means - current)
}

# Modelling the total-dependent trend in the simulated LR (and the variance around the trend).

log.totals <- log(sim.totals)
trend.fit <- lowess(x=log.totals, y=sim.LR, f=0.2)
expected.LR <- spline(trend.fit$x, trend.fit$y, xout=log(obs.totals))$y

spread <- sim.LR/trend.fit$y
plot(sim.totals, spread, log="x")
abline(h=1, col="red")

# Computing a p-value for each observed value.

plot(obs.totals, obs.LR, log="x")
points(obs.totals, expected.LR, col="red", pch=16, cex=0.2)

obs.spread <- obs.LR/expected.LR
p <- numeric(length(obs.LR))
for (x in seq_along(p)) {
    current <- obs.totals[x]*sqrt(2) >= sim.totals & obs.totals[x]/sqrt(2) <= sim.totals
    p[x] <- (sum(spread[current] > obs.spread[x]) + 1)/(sum(current)+1)
}

fdr <- p.adjust(p, method="BH")
all.fdr <- rep(1, nrow(original))
all.fdr[!ambient] <- fdr

write.table(data.frame(Cell=colnames(original), FDR=all.fdr), file=paste0(outfix, ".tsv"),
	    quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# Making a pretty plot.

sig <- fdr <= 0.05
png(paste0(outfix, ".png"), width=12, height=12, units="in", pointsize=12, res=120)
plot(obs.totals, obs.LR, log="x", col=ifelse(sig, "red", "black"))
curve(trend.FUN(x), col="dodgerblue", add=TRUE)
legend("topleft", sprintf("%.2f", sum(sig)/length(sig) * 100), bty="n")
dev.off()
