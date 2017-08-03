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

#   grouping <- cut(seq_len(ncol(observed)), max(2, ceiling(ncol(observed)/100)))
#   output <- numeric(ncol(observed))
#   for (g in levels(grouping)) { 
#       chosen <- grouping == g
#       M <- as.matrix(observed[,chosen,drop=FALSE])
#       expected <- outer(ambient.prop, colSums(M))
#       output[chosen] <- colSums(dpois(M, lambda=M, log=TRUE)) - colSums(dpois(M, lambda=expected, log=TRUE))
#   }
#   return(output)
}

N <- 10000
test.totals <- 2^seq(from=0, to=log2(max(colSums(m))), length.out=N)
sim <- as(matrix(rpois(length(ambient.prop)*N, outer(ambient.prop, test.totals)), ncol=N), "dgCMatrix")
sim.LR <- LRFUN(sim, prop=ambient.prop)

obs <- m
obs.LR <- LRFUN(obs, prop=ambient.prop)

# Modelling the total-dependent trend in the LR (and the variance around the trend).

lmeans <- log(colSums(sim))
keep <- is.finite(lmeans)
lmeans.keep <- lmeans[keep]
sim.LR.keep <- sim.LR[keep]
trend.fit <- loess(sim.LR.keep ~ lmeans.keep, span=0.3, degree=1)
mean.FUN <- function(total) {
    ltotal <- pmin(log(total), max(lmeans.keep))
    predict(trend.fit, data.frame(lmeans.keep=ltotal))
}

spread <- log(sim.LR.keep/fitted(trend.fit))
by.mean <- cut(lmeans.keep, 50)
bin.var <- unlist(lapply(split(spread, by.mean), var))
bin.mean <- unlist(lapply(split(lmeans.keep, by.mean), mean))
var.FUN <- approxfun(bin.mean, bin.var, rule=2)

getPValue <- function(total, LR) {
    dev <- log(LR/mean.FUN(total))
    pnorm(dev, mean=0, sd=sqrt(var.FUN(log(total))), lower=FALSE)
}


#y <- numeric(nlevels(group))
#for (i in seq_len(nlevels(group))) {
#    p <- sum(dnbinom(x[[i]], lambda=ambient.prop * sum(x[[i]], size=5), log=TRUE))
#    y[i] <- p
#}

obs.totals <- colSums(obs)
p <- getPValue(obs.totals, obs.LR)
p[obs.totals < 250] <- NA # Killing these guys.
fdr <- p.adjust(p, method="BH")
fdr[is.na(fdr)] <- 1

write.table(data.frame(Cell=colnames(m), FDR=fdr), file=paste0(outfix, ".tsv"),
	    quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# Making a pretty plot.

sig <- fdr <= 0.05
png(paste0(outfix, ".png"), width=12, height=12, units="in", pointsize=12, res=120)
plot(colSums(obs), obs.LR, log="x", col=ifelse(sig, "red", "black"))
curve(mean.FUN(x), col="dodgerblue", add=TRUE)
#points(colSums(sim), sim.LR, col="dodgerblue")
legend("topleft", sprintf("%.2f", sum(sig)/length(sig) * 100), bty="n")
dev.off()
