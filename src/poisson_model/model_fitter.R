library(Matrix)

computePValue <- function(m, lower=250, tol=0.5, npts=10000) 
# A function to compute a non-ambient p-value for each barcode.
# 
# written by Aaron Lun
# created 1 August 2017    
{
    discard <- rowSums(m) == 0
    m <- m[!discard,]

    # Computing the average profile.
    umi.sum <- colSums(m)
    ambient <- (umi.sum <= lower) # lower => "T" in the text.
    ambient.cells <- m[,ambient]
    ambient.prof <- rowSums(ambient.cells)
    ambient.prop <- edgeR::goodTuringProportions(ambient.prof)

    # Removing cells from the ambient.
    obs <- m[,!ambient]
    obs.totals <- colSums(obs)

    # Calculating the likelihood ratio.
    if (is(obs, "dgCMatrix")) {
        i <- obs@i + 1L
        j <- rep(seq_len(ncol(obs)), diff(obs@p)) 
        x <- obs@x
    } else if (is(obs, "dgTMatrix")) {
        i <- obs@i + 1L
        j <- obs@j + 1L
        x <- obs@x
    } else {
        stop("unsupported matrix type")
    }
    
    p.n0 <- x * log(x/(ambient.prop[i]*obs.totals[j])) 
    by.col <- aggregate(p.n0, list(Col=j), sum)
    obs.LR <- numeric(length(obs.totals))
    obs.LR[by.col$Col] <- by.col$x

    # Computing a simulation with "npts" entries for any "tol"-fold interval around the interrogation point..
    tol <- 0.5
    lower.pt <- log2(min(obs.totals))-tol
    upper.pt <- log2(max(obs.totals))+tol
    S <- round(npts * (upper.pt - lower.pt)/(2*tol)) # npts => R in the text.
    sim.totals <- 2^seq(from=lower.pt, to=upper.pt, length.out=S)

    # Computing the deviance estimate for simulated runs.
    sim.LR <- numeric(S)
    for (x in seq_along(sim.LR)) {
        cur.means <- ambient.prop*sim.totals[x]
        current <- rpois(length(ambient.prop), lambda=cur.means)
        sim.LR[x] <- sum(current * log(current/cur.means), na.rm=TRUE) + sum(cur.means - current)
    }

    # Modelling the total-dependent trend in the simulated LR (and the variance around the trend).
    log.totals <- log(sim.totals)
    trend.fit <- lowess(x=log.totals, y=sim.LR, f=0.2)
    spread <- sim.LR/trend.fit$y
    expected.LR <- spline(trend.fit$x, trend.fit$y, xout=log(obs.totals))$y  

    # Computing a p-value for each observed value.
    obs.spread <- obs.LR/expected.LR
    p <- numeric(length(obs.LR))
    limited <- logical(length(obs.LR))
    for (x in seq_along(p)) {
        current <- obs.totals[x]*sqrt(2) >= sim.totals & obs.totals[x]/sqrt(2) <= sim.totals
        nexceed <- sum(spread[current] >= obs.spread[x])
        limited[x] <- nexceed==0L
        p[x] <- (nexceed + 1)/(sum(current)+1)
    }

    all.p <- all.lr <- all.exp <- rep(NA_real_, length(ambient))
    all.lim <- rep(NA, length(ambient))
    all.p[!ambient] <- p
    all.lr[!ambient] <- obs.LR
    all.exp[!ambient] <- expected.LR
    all.lim[!ambient] <- limited
    return(data.frame(Total=umi.sum, LR=all.lr, Expected=all.exp, PValue=all.p, 
                      Limited=all.lim, row.names=colnames(m)))
}

findInflectionPoint <- function(m, lower=250) 
# A function to identify the inflection point from a reverse cumulative distribution.
# 
# written by Aaron Lun
# created 7 August 2017    
{
    totals <- colSums(m)
    totals <- totals[totals >= lower]
    
    stuff <- rle(sort(totals, decreasing=TRUE))
    y <- log(stuff$values)
    x <- log(cumsum(stuff$lengths))
    
    grad <- diff(y)/diff(x)
    threshold <- y[which.min(grad)+1]
    return(unname(exp(threshold)))
}

detectCells <- function(m, lower=250, scale=2, ...) 
# Combined function that puts these all together, always keeping cells above the inflection
# point (they are given p-values of 0, as they are always rejected). 
# 
# written by Aaron Lun
# created 7 August 2017
{
    stats <- computePValue(m, lower=lower, ...)
    inflection <- findInflectionPoint(m, lower=lower)
    always <- stats$Total >= inflection*scale
    tmp <- stats$PValue
    tmp[always] <- 0
    stats$FDR <- p.adjust(tmp, method="BH")
    return(stats)
}
