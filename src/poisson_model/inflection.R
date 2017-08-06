if (!interactive()) {
    cmd.opts <- commandArgs(trailingOnly=TRUE)
    inname <- cmd.opts[1]
    outname <- cmd.opts[2]
}

library(Matrix)
m <- readRDS(inname)
totals <- colSums(m)
totals <- totals[totals > 250]

stuff <- rle(sort(totals, decreasing=TRUE))
y <- log(stuff$values)
x <- log(cumsum(stuff$lengths))

grad <- diff(y)/diff(x)
threshold <- y[which.min(grad)+1]
print(unname(exp(threshold)))

pdf(outname) 
plot(x, y, xlab="Log-rank", ylab="Log-total")
abline(h=threshold, col="red")
text(1, threshold, pos=3, exp(threshold), col="red")
dev.off()
