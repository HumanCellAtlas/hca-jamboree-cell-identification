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

original <- m
m <- m[,!ambient]

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
   p.n0 <- dpois(x, lambda=prop[i] * total[j], log=TRUE) - dpois(0, lambda=prop[i] * total[j], log=TRUE) 
   p.n0.mat <- sparseMatrix(i=i, j=j, x=p.n0, dims=c(nrow(observed), ncol(observed)))
   p.n0.sum <- colSums(p.n0.mat)

   all.totals <- unique(total)
   p.0 <- matrix(dpois(0, lambda=outer(prop, all.totals), log=TRUE), nrow=length(prop), ncol=length(all.totals))
   p.0.sum_x <- colSums(p.0)
   p.0.sum <- p.0.sum_x[match(total, all.totals)]
   
   # Likelihood of a perfit fit (no need to deal with zeroes, as these have likelihood=1).
   p.perfect <- dpois(x, lambda=x, log=TRUE) 
   p.perfect.mat <- sparseMatrix(i=i, j=j, x=p.perfect, dims=c(nrow(observed), ncol(observed)))
   p.perfect.sum <- colSums(p.perfect.mat)

   return(p.perfect.sum - (p.n0.sum + p.0.sum))
}

obs <- m
obs.totals <- colSums(obs)
obs.LR <- LRFUN(obs, prop=ambient.prop)

N <- 20000
sim.totals <- 2^seq(from=log2(min(obs.totals))-1, to=log2(max(obs.totals))+1, length.out=N)
sim <- as(matrix(rpois(length(ambient.prop)*N, outer(ambient.prop, sim.totals)), ncol=N), "dgCMatrix")
sim.LR <- LRFUN(sim, prop=ambient.prop)

# Modelling the total-dependent trend in the LR (and the variance around the trend).

log.totals <- log(sim.totals)
trend.fit <- loess(sim.LR ~ log.totals, span=0.2, degree=1)
trend.FUN <- function(total) {
    ltotal <- pmin(log(total), max(log.totals))
    predict(trend.fit, data.frame(log.totals=ltotal))
}

spread <- sim.LR/fitted(trend.fit)
f <- cut(log.totals, 50)
by.total <- split(spread, f)

bin.var <- unlist(lapply(by.total, var))
bin.mean <- unlist(lapply(by.total, mean))
bin.x <- unlist(lapply(split(log.totals, f), mean))
var.FUN <- splinefun(y=bin.var, x=bin.x)
mean.FUN <- splinefun(y=bin.mean, x=bin.x)

# Computing a p-value for each observed value.

getPValue <- function(total, LR) {
    # Getting the expected value.
    fitted <- trend.FUN(total)

    # Getting the parameters of the "bin" around it. 
    ltotal <- log(total)
    cur.var <- var.FUN(ltotal)
    cur.mean <- mean.FUN(ltotal)

    # Modelling with a gamma distribution
    rate <- cur.mean/cur.var
    shape <- rate * cur.mean

    pgamma(LR/fitted, rate=rate, shape=shape, lower.tail=FALSE)
}

p <- getPValue(obs.totals, obs.LR)
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
