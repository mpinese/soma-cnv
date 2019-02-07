#!/usr/bin/env Rscript

VERSION = "v0.0.3  7 Feb 2019"

suppressPackageStartupMessages(library(ggplot2))    # For plotting only
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(tvd))        # For plotting only


# Perform per-sample depth anomaly fitting to realise the model:
#    lambda_ij = D_i*(k_ij*a_j*ac_ij)
# => log(lambda_ij) = log(D_i) + log(k_ij / 2) + log(a_j) + log(ac_ij)
# Where
#   lambda_ij is the expected depth in sample i at locus j
#   D_i       is the global depth of sample i
#   k_ij      is the copy number (diploid = 2) of sample i at locus j (might be fractional due to subclonality)
#   a_j       is the affinity anomaly at locus j (equal to affinity$log2_reldp)
#   ac_ij     is the depth anomaly at locus j for sample i (to be fit in this code)
#
# ac_ij will combine sample-specific GC signals as well as sample-specifc
# affinity signal not captured by the first-order global correction in a_j.
#
# So the goal is to fit ac_ij given d_ij ~ Pois(lambda_ij), D_i, and a_j,
# and assuming k_ij = 2.
#
# From the above, assuming k_ij = 2:
#     log(lambda_ij) = [log(D_i) + log(a_j)] + log(ac_ij)
#  => log(lambda_ij) - [log(D_i) + log(a_j)] ~ s(a_j) + s(gc_j)
# Where
#   gc_j      is the gc content at locus j.
#
# We further split the s(gc_j) into multiple terms based on differing
# window sizes for GC content averaging.  PCA is used to remove 
# co-linearity.
#
# The final output from the fit is a sample-specific estimate of
# log(lambda_ij) - [log(D_i) + log(a_j)] for each locus j.  To this
# estimate we add log(D_i) and log(a_j) to yield the final estimate
# of log(lambda_ij) if k_ij = 2.  This is the output stored in 
# lambda.
#
# Perform a basic fit for now, setting lambda_ij = d_ij.  The
# large data size should make up for this liberty with the error model.
fitDepthAnomaly = function(sample_data, affinity, gc, trim = 0)
{
    message("  Merging data...")
    data = merge(sample_data[,c("chrom", "pos", "dp", "ad")], affinity, by = c("chrom", "pos"), all = FALSE)
    data = merge(data, gc, by = c("chrom", "pos"), all = FALSE)
    # message(sprintf(" %d shared loci remain.", nrow(data)))

    logD = log(mean(data$dp, trim = trim))
    logA = log(data$affinity)

    message("  Rotating GC...")
    gc_pca = prcomp(as.matrix(data[,grepl("^gc[0-9]+$", colnames(data))]), center = TRUE, scale = TRUE)

    anomaly = log(data$dp) - logD - logA
    message("  Modelling anomaly...")
    fit = gam(anomaly ~ s(logA) + s(gc_pca$x[,1]) + s(gc_pca$x[,2]) + s(gc_pca$x[,3]) + s(gc_pca$x[,4]) + s(gc_pca$x[,5]))
    # print(anova(fit))
    anomaly_pred = predict(fit)

    log_lambda_hat = logD + logA + anomaly_pred

    # Correct for a slight bias introduced by the move between log and
    # linear space
    message("  Correcting bias...")
    log_lambda_hat = log_lambda_hat + log(mean(data$dp / exp(log_lambda_hat), trim = trim))

    # Overdispersion is observed in some cases.  Address by additionally
    # fitting a negative binomial.
    # fit.nb = gamlss(formula = dp ~ log_lambda_hat, sigma.formula = ~ 1, family = NBI(), data = data)
    # print(summary(fit.nb))

    # Save fit values
    data$pois.lambda = exp(log_lambda_hat)
    # data$nb.mu = exp(predict(fit.nb, "mu"))
    # data$nb.size = 1/exp(predict(fit.nb, "sigma"))

    data
}


testWindowDifference = function(data1, data2)
{
    ad1 = data1$ad
    ad2 = data2$ad
    rd1 = data1$dp - ad1
    rd2 = data2$dp - ad2
    anomaly1 = log(c(rd1, ad1) / rep(data1$pois.lambda, 2))
    anomaly2 = log(c(rd2, ad2) / rep(data2$pois.lambda, 2))
    suppressWarnings(ks.test(anomaly1, anomaly2)$p.value)
}


finetuneWindows = function(windows, data, search_size)
{
    # Fune-tune windows identified by findWindows: for every
    # intra-chromosomal break, find the optimal breakpoint
    # within a +/- search_size interval.
    ddply(windows, .(chrom), finetuneWindows.chrom, data = data, search_size = search_size)
}


finetuneWindows.chrom = function(chrom_windows, data, search_size, window_size = 100)
{
    # Note an issue with this algorithm: it only optimises a single breakpoint, which
    # can be error-prone if multiple adjacent breakpoints are incorrect.  More accurately
    # it should optimise all breakpoints simultaneously.

    if (nrow(chrom_windows) == 1)
        return(chrom_windows)       # No breakpoints so nothing to do

    for (i in 1:(nrow(chrom_windows)-1))
    {
        # Work on the breakpoint between chrom_windows[i,] and chrom_windows[i+1,]
        # The breakpoint position is currently between chrom_windows$end_index[i] and
        # chrom_windows$start_index[i+1].  Search in [chrom_windows$end_index[i]-search_size+1,
        # chrom_windows$end_index[i]+search_size] for a better breakpoint.
        left_start = chrom_windows$start_index[i]
        right_end = chrom_windows$end_index[i+1]
        break_search_start = chrom_windows$end_index[i]-search_size
        break_search_end = chrom_windows$start_index[i+1]+search_size

        # Search for the optimal breakpoint
        candidate_breaks = break_search_start:break_search_end
        candidate_break_pvals = sapply(candidate_breaks, function(j) testWindowDifference(data[(j-window_size+1):j,,drop=FALSE], data[(j+1):(j+window_size),,drop=FALSE]))
        best_break = candidate_breaks[which.min(candidate_break_pvals)]

        # Update the windows to match
        chrom_windows$end_index[i] = best_break
        chrom_windows$start_index[i+1] = best_break+1
        chrom_windows$pr_vs_prev[i+1] = min(candidate_break_pvals, na.rm = TRUE)
    }

    chrom_windows
}


# First split xsomes into small windows.
# Progressively merge by KS test.
# Calculate fits based on windowed regions -- direct optim.
# Requires data to be in chrom, pos order
findWindows = function(data, initial_size = 100, alpha = 0.01)
{
    # Define windows
    data$original_index = 1:nrow(data)
    windows = ddply(data, .(chrom), function(cd) {
        cd$window_id = paste(cd$chrom, floor(cd$original_index/initial_size)+1, sep = ":")
        windows = data.frame(window_id = unique(cd$window_id))
        windows$start_index = tapply(cd$original_index, cd$window_id, min)[windows$window_id]
        windows$end_index = tapply(cd$original_index, cd$window_id, max)[windows$window_id]
        windows
    })

    # Calculate initial split/merge P-values
    windows$pr_vs_prev = c(-1, mapply(function(start1, end1, start2, end2) testWindowDifference(data[start1:end1,], data[start2:end2,]),
        windows$start[-nrow(windows)], windows$end[-nrow(windows)], windows$start[-1], windows$end[-1]))
    # Chromosomes are never merged
    windows$pr_vs_prev[(c(windows$chrom, "") != c("", windows$chrom))[1:nrow(windows)]] = -1

    # Progressive merging of windows
    ntests = sum(windows$pr_vs_prev >= 0)
    while (TRUE)
    {
        # Find the highest P-value
        pmax_index = which.max(windows$pr_vs_prev)
        pmax = windows$pr_vs_prev[pmax_index]

        # Check for termination
        if (pmax * ntests < alpha)
            break

        # Merge windows pmax_index and pmax_index-1
        # message(sprintf("Merging indices %d (%s:%d-%d) and %d (%s:%d-%d), P-value %f", 
        #     pmax_index-1, windows$chrom[pmax_index - 1], windows$start_index[pmax_index - 1], windows$end_index[pmax_index - 1], 
        #     pmax_index, windows$chrom[pmax_index - 1], windows$start_index[pmax_index], windows$end_index[pmax_index], 
        #     pmax))
        windows$end_index[pmax_index - 1] = windows$end_index[pmax_index]
        windows = windows[-pmax_index,]

        # Recalculate the split/merge P-values for either side of the new window
        if (pmax_index - 2 > 0)
        {
            # Upstream is a valid window.
            if (windows$chrom[pmax_index - 2] != windows$chrom[pmax_index - 1])
                windows$pr_vs_prev[pmax_index - 1] = -1    # Upstream is a different chrom
            else
                windows$pr_vs_prev[pmax_index - 1] = testWindowDifference(data[windows$start[pmax_index-2]:windows$end[pmax_index-2],], data[windows$start[pmax_index-1]:windows$end[pmax_index-1],])
            ntests = ntests + 1
        }
        if (pmax_index - 1 != nrow(windows))
        {
            # Downstream is a valid window
            if (windows$chrom[pmax_index] != windows$chrom[pmax_index - 1])
                windows$pr_vs_prev[pmax_index] = -1        # Downstream is a different chrom
            else
                windows$pr_vs_prev[pmax_index] = testWindowDifference(data[windows$start[pmax_index-1]:windows$end[pmax_index-1],], data[windows$start[pmax_index]:windows$end[pmax_index],])
            ntests = ntests + 1
        }
    }

    # Fine tune breakpoints
    windows = finetuneWindows(windows, data, initial_size)

    # Add coordinates and return
    windows$start_pos = data$pos[windows$start_index]
    windows$end_pos = data$pos[windows$end_index]
    return(windows)
}


llik.pois = function(dp1, dp2, f, k1, k2, lambda, w)
{
    if (is.null(w))
        w = rep(1, length(lambda))

    lambda1 = (k1*f + (1-f))/2 * lambda
    lambda2 = (k2*f + (1-f))/2 * lambda
    d12 = dpois(dp1, lambda = lambda1)*dpois(dp2, lambda = lambda2)
    d21 = dpois(dp1, lambda = lambda2)*dpois(dp2, lambda = lambda1)
    d = 0.5*d12 + 0.5*d21
    sum(log(d)*w)
}


llik.nb = function(dp1, dp2, f, k1, k2, mu, size, w)
{
    if (is.null(w))
        w = rep(1, length(mu))

    mu1 = (k1*f + (1-f))/2 * mu
    mu2 = (k2*f + (1-f))/2 * mu
    d12 = dnbinom(dp1, mu = mu1, size = size/2)*dnbinom(dp2, mu = mu2, size = size/2)
    d21 = dnbinom(dp1, mu = mu2, size = size/2)*dnbinom(dp2, mu = mu1, size = size/2)
    d = 0.5*d12 + 0.5*d21
    sum(log(d)*w)
}


fitWindow.kfs = function(window_data, ks, f, isize, w)
{
    # Find the most likely combination of k in the ks, given f and isize,
    # for the data in window_data.  ks is a 2-column matrix with
    # integer columns k1, k2.  Return as list(k1, k2, llik)
    dp1 = window_data$dp - window_data$ad
    dp2 = window_data$ad
    if (isize == 0)
    {
        lliks = apply(ks, 1, function(k) llik.pois(dp1, dp2, f, k[1], k[2], window_data$pois.lambda, w))
        llik_gldup = llik.pois(dp1, dp2, f=1, 1, 2, window_data$pois.lambda, w)
    }
    else
    {
        lliks = apply(ks, 1, function(k) llik.nb(dp1, dp2, f, k[1], k[2], window_data$pois.lambda, 1/isize, w))
        llik_gldup = llik.nb(dp1, dp2, f=1, 1, 2, window_data$pois.lambda, 1/isize, w)
    }

    if (llik_gldup >= max(lliks))
    {
        # Most likely a germline duplication
        return(list(k1 = 1, k2 = 2, f = 1, isize = isize, llik = llik_gldup, type = "gldup"))
    }

    besti = which.max(lliks)
    list(k1 = ks[besti,1], k2 = ks[besti,2], f = f, isize = isize, llik = lliks[besti], type = "som")
}


fitWindows.maxkfs = function(data, windows, maxk, f, isize, w, mink)
{
    # For each window find the most likely combination of 
    # k1, k2 in mink..maxk, given f and isize.  Returns a 
    # data frame of window fits.  Note that 0 <= mink <= maxk, 
    # maxk >= 1.
    stopifnot(maxk >= 1 && mink <= maxk)
    ks = expand.grid(k1 = mink:maxk, k2 = mink:maxk)
    ks = ks[ks$k1 >= ks$k2,]
    window_fits = ddply(windows, colnames(windows), function(d) cbind(d, fit = fitWindow.kfs(data[d$start_index:d$end_index,], ks, f, isize, w)))
    window_fits
}


findLocalMaxima = function(x)
{
    which(diff(diff(c(-Inf, x, -Inf)) < 0) == 1)
}


polishLocalMaximum.poisson = function(objective, fmin, fmax)
{
    opt = optimise(objective, lower = fmin, upper = fmax, maximum = TRUE)
    c(f = opt$maximum, isize = 0, llik = opt$objective)
}


polishLocalMaximum.nb = function(objective, fmin, fmax)
{
    opt = optim(
        par = c(f = (fmin+fmax)/2, isize = 1/500), 
        fn = function(par) objective(par[1], par[2]), 
        method = "L-BFGS-B", 
        lower = c(0.01, 1/1000), 
        upper = c(0.99, 1/1),
        control = list(fnscale = -1, parscale = c(0.01, 1)))
    if (opt$convergence != 0)
        warning(sprintf("Negative binomial optimization failed to converge, code %d, message\"%s\"", opt$convergence, opt$message))
    c(f = opt$par[[1]], isize = opt$par[[2]], llik = opt$value)
}


fitWindows.maxkf = function(data, windows, maxk = 1, f = 0, w = w, mink = 1)
{
    # Given maxk and f, search over s in [1, 1000] to find the most likely
    # value.
    objective = function(isize) sum(fitWindows.maxkfs(data, windows, maxk, f, isize, w, 0)$fit.llik)

    opt = optimise(objective, lower = 1/1000, upper = 1, maximum = TRUE)
    fitWindows.maxkfs(data, windows, maxk = maxk, f = f, isize = opt$maximum, w = w, mink = mink)
}


fitWindows.maxk = function(data, windows, maxk, family, w = NULL, plot = FALSE)
{
    # Given maxk, search over f in (0, 1) and for family = "nb", s in (1, 1000],
    # to find the most likely value.
    objective.rough = function(f) sum(fitWindows.maxkfs(data, windows, maxk, f, 0, w, 0)$fit.llik)
    objective.fine = function(f, isize) sum(fitWindows.maxkfs(data, windows, maxk, f, isize, w, 0)$fit.llik)

    # Assess f at fairly dense equally-spaced points.  Place a slight bias
    # towards low f.  Regardless of family, use a Poisson model here for 
    # speed (objective.rough).  Use the correct family (objective.fine) 
    # in the subsequent polishing step.
    grid_pow = 0.05
    test.f = rev(seq(0.99^(-grid_pow), 0.05^(-grid_pow), length.out = maxk*10))^(-1/grid_pow)
    test.ll = laply(test.f, objective.rough)
    local_maxima.i = findLocalMaxima(test.ll)

    # Polish each local maximum
    polish.objective = c("poisson" = objective.rough, "nb" = objective.fine)[[family]]
    polish.function = c("poisson" = polishLocalMaximum.poisson, "nb" = polishLocalMaximum.nb)[[family]]
    polished_maxima = sapply(local_maxima.i, function(i) polish.function(polish.objective, fmin = test.f[max(1, i - 1)], fmax = test.f[min(length(test.f), i + 1)]))

    # Find the global optimum.  In cases where there are multiple nearly-equivalent
    # optima, choose the one corresponding to the highest value of f.  Define near-equivalence
    # between optima i and j here as |LL_i - LL_j| <= 1.
    max_ll = max(polished_maxima["llik",], na.rm = TRUE)
    delta_ll = max_ll - polished_maxima["llik",]
    near_optimal_maxima = polished_maxima[,delta_ll < 1,drop=FALSE]
    global_optimum = as.list(near_optimal_maxima[,which.max(near_optimal_maxima["f",])])

    if (plot)
    {
        plot(test.ll ~ test.f, xlim = c(0, 1), xlab = "f", ylab = "Log likelihood", main = "Purity scan")
        points(polished_maxima["f",], polished_maxima["llik",], col = "red", pch = 19)
        points(global_optimum$f, global_optimum$llik, col = "red", pch = 1, cex = 2)
    }

    message(sprintf("    mink=0, maxk=%d, optf=%.2f, optisize=%.4f, LL=%.2f", maxk, global_optimum$f, global_optimum$isize, global_optimum$llik))

    list(
        optimum = data.frame(f = global_optimum$f, isize = global_optimum$isize, maxk = maxk, ll = global_optimum$llik), 
        search = data.frame(f = c(test.f, unlist(polished_maxima["f",])), isize = c(rep(0, length(test.f)), unlist(polished_maxima["isize",])), maxk = maxk, ll = c(test.ll, unlist(polished_maxima["llik",])))
    )
}


mergeWindowFits = function(window_fits)
{
    # Merge consecutive window fits if they are on the same chromosome
    # and have matching k1, k2, and f values.  Return the merged
    # fits as a data frame.
    if (nrow(window_fits) <= 1)
        return(window_fits)
    window_fits = window_fits[order(window_fits$chrom, window_fits$start_index),]

    merged = window_fits[FALSE,]
    merged_window.key = c(window_fits$chrom[1], window_fits$fit.k1[1], window_fits$fit.k2[1], window_fits$fit.f[1])
    merged_window.start_i = 1

    for (i in 2:nrow(window_fits))
    {
        current_window.key = c(window_fits$chrom[i], window_fits$fit.k1[i], window_fits$fit.k2[i], window_fits$fit.f[i])

        if (any(current_window.key != merged_window.key))
        {
            # A new window.  Add the old one to the merged data.
            merged_window.end_i = i - 1
            merged_window = window_fits[merged_window.start_i,]
            merged_window$end_index = window_fits$end_index[merged_window.end_i]
            merged_window$fit.llik = sum(window_fits$fit.llik[merged_window.start_i:merged_window.end_i])
            merged = rbind(merged, merged_window)

            # Start the new window
            merged_window.start_i = i
            merged_window.key = current_window.key
        }
    }

    # Add the last merged window
    merged_window.end_i = i
    merged_window = window_fits[merged_window.start_i,]
    merged_window$end_index = window_fits$end_index[merged_window.end_i]
    merged_window$fit.llik = sum(window_fits$fit.llik[merged_window.start_i:merged_window.end_i])
    merged = rbind(merged, merged_window)

    # Return without pr_vs_prev, which is invalidated by the merging
    merged[,colnames(merged) != "pr_vs_prev"]
}


fitWindows = function(data, max_maxk, family = c("poisson", "nb"), w = NULL, window.initial_size = 100, window.alpha = 0.01)
{
    family = match.arg(family)

    message("  Windowing genome")
    windows = findWindows(data, initial_size = window.initial_size, alpha = window.alpha)

    message("  Fitting null")
    if (family == "poisson")
        null_fit = fitWindows.maxkfs(data, windows, maxk = 1, f = 0, isize = 0, w = w, mink = 1)
    else
        null_fit = fitWindows.maxkf(data, windows, maxk = 1, f = 0, w = w, mink = 1)        
    null_loglik = sum(null_fit$fit.llik)
    message(sprintf("    mink=1, maxk=1, optf=0, optisize=%.4f, LL=%.2f", null_fit$fit.isize[1], null_loglik))

    message("  Fitting candidate models")
    model_search = llply(1:max_maxk, function(maxk) fitWindows.maxk(data, windows, maxk, family, w))
    models = ldply(model_search, function(model_search_result) model_search_result$optimum)
    models$p = log((models$maxk+1)*(models$maxk+2)/2)*nrow(windows) + c("poisson" = 0, "nb" = 1)[family]
    # Note log term, which is to account for the discrete nature of the ploidy estimates.
    # There are s=(maxk+1)*(maxk+2)/2 possible ploidy states for each window, with penalty
    # approximately log(s).  Multiplication by nrow(windows) accounts for each window having
    # its own independent ploidy estimate.  The constant term is for the global isize 
    # parameter in nb models.

    models = rbind(models, c(f = 0, isize = null_fit$fit.isize[1], maxk = NA, ll = null_loglik, p = c("poisson" = 0, "nb" = 1)[family]))
    models$BIC = -2*models$ll + log(nrow(data))*models$p
    besti = which.min(models$BIC)

    if (is.na(models$maxk[besti]))
    {
        best.mink = 1
        best.maxk = 1
        best.f = 0
        best.isize = null_fit$fit.isize[1]
    }
    else
    {
        best.mink = 0
        best.maxk = models$maxk[besti]
        best.f = models$f[besti]
        best.isize = models$isize[besti]
    }
    message(sprintf("  Optimal model: f=%.2f, kmin=%d, kmax=%d, isize=%.4f, BIC=%.2f", best.f, best.mink, best.maxk, best.isize, min(models$BIC)))

    fit = fitWindows.maxkfs(data, windows, maxk = best.maxk, f = best.f, isize = best.isize, w = w, mink = best.mink)
    fit = fit[order(fit$chrom, fit$start_index),]

    list(models = models, model_search = ldply(model_search, function(model_search_result) model_search_result$search), fit.orig = fit, fit = mergeWindowFits(fit))
}


plotWindowFit = function(data, fit)
{
    fit$kf1 = fit$fit.f*fit$fit.k1/2 + (1-fit$fit.f)/2
    fit$kf2 = fit$fit.f*fit$fit.k2/2 + (1-fit$fit.f)/2
    fit$d = log2(fit$kf1 + fit$kf2)
    fit$v1 = fit$kf1 / (fit$kf1 + fit$kf2)
    fit$v2 = fit$kf2 / (fit$kf1 + fit$kf2)
    fit$colus = fit$fit.llik / (fit$end_index - fit$start_index + 1)
    fit$col = pmax(0, pmin(1, (fit$colus - (-10)) / 5))
    fit$coldisc = as.integer(ceiling(fit$col * 100))
    fit$coldisc[fit$coldisc == 0] = 1
    pal = rainbow(100, end = 0.33)

    chrom_boundaries = which(data$chrom[-1] != data$chrom[-nrow(data)])
    chrom_midpoints = (c(0, chrom_boundaries) + c(chrom_boundaries, nrow(data))) / 2
    names(chrom_midpoints) = data$chrom[!duplicated(data$chrom)]

    par(mfrow = c(2, 1))
    plot(log2(data$dp) - log2(data$pois.lambda), pch = ".", col = rgb(0, 0, 0, 0.25), xlab = "Genomic position", ylab = "Depth anomaly (log2)", main = "Depth", xaxt = "n")
    axis(1, at = chrom_midpoints, labels = names(chrom_midpoints))
#    lines(runmed(log2(data$dp) - log2(data$pois.lambda), k = 151), col = "red", lwd = 1)
    lines(tvd1d(log2(data$dp) - log2(data$pois.lambda), lambda = 5), col = "red", lwd = 2)
    abline(v = chrom_boundaries, col = "blue", lwd = 2)
    abline(v = fit$start_index, col = "red", lty = "dotted")
    segments(fit$start_index, fit$d, fit$end_index, fit$d, col = pal[fit$coldisc], lwd = 3)

    plot(data$ad / jitter(data$dp, amount = 0.5), pch = ".", ylim = c(0, 1), col = rgb(0, 0, 0, 0.25), xlab = "Genomic position", ylab = "VAF", main = "Allele frequency", xaxt = "n")
    axis(1, at = chrom_midpoints, labels = names(chrom_midpoints))
    abline(v = chrom_boundaries, col = "blue", lwd = 2)
    abline(v = fit$start_index, col = "red", lty = "dotted")
    segments(fit$start_index, fit$v1, fit$end_index, fit$v1, col = pal[fit$coldisc], lwd = 3)
    segments(fit$start_index, fit$v2, fit$end_index, fit$v2, col = pal[fit$coldisc], lwd = 3)
    par(mfrow = c(1, 1))
}


plotRawData = function(data)
{
    chrom_boundaries = which(data$chrom[-1] != data$chrom[-nrow(data)])
    chrom_midpoints = (c(0, chrom_boundaries) + c(chrom_boundaries, nrow(data))) / 2
    names(chrom_midpoints) = data$chrom[!duplicated(data$chrom)]

    par(mfrow = c(2, 1))
    plot(jitter(data$dp, amount = 0.5), pch = ".", col = rgb(0, 0, 0, 0.25), xlab = "Genomic position", ylab = "Depth", main = "Depth", xaxt = "n")
    axis(1, at = chrom_midpoints, labels = names(chrom_midpoints))
    abline(v = chrom_boundaries, col = "blue", lwd = 2)

    plot(data$ad / jitter(data$dp, amount = 0.5), pch = ".", ylim = c(0, 1), col = rgb(0, 0, 0, 0.25), xlab = "Genomic position", ylab = "VAF", main = "Allele frequency", xaxt = "n")
    axis(1, at = chrom_midpoints, labels = names(chrom_midpoints))
    abline(v = chrom_boundaries, col = "blue", lwd = 2)
    par(mfrow = c(1, 1))
}


markFixedLoci = function(data, epsilon, alpha, plot = FALSE)
{
    # Test loci for consistency under the model of being 
    # generated from a fixed locus, with any minor reads
    # being solely due to error.  Return data augmented with
    # the following columns:
    #   prRR: P-value against this locus being derived from
    #         a hom ref (RR) locus.
    #   prAA: P-value against this locus being derived from
    #         a hom alt (AA or AB) locus.
    #   het:  boolean, TRUE if no evidence of this locus being
    #         fixed, else FALSE.
    # 
    # Specifically test for deviation from the two hypotheses:
    #   1. Locus is hom ref
    #   2. Locus is hom alt or het alt
    # 
    # The two components of hypothesis 2 are equivalent,
    # as the counts do not distinguish between the alt alleles.
    #
    # Hypothesis 1 test:
    #   Reject if Pr(ad, rd | G=RR) < alpha
    # Hypothesis 2 test:
    #   Reject if Pr(ad, rd | G=AA) < alpha
    #
    # Overall test:
    #   Reject if BOTH H1 and H2 are rejected.
    # 
    # Note that the overall test size is *not* equal to alpha;
    # it will be greater.

    # Test H1.  In the context of G=RR, any error (at rate epsilon)
    # will lead to an alt read, therefore the rate is epsilon.
    data$prRR = pbinom(data$ad, data$dp, epsilon, lower.tail = FALSE)

    # Test H2.  In the context of G=AA, 2/3 errors yield another
    # alt allele, so wouldn't appear as a ref allele.  The rate is
    # therefore epsilon/3, accounting for the 1/3 cases that would
    # be visible as an A -> R change.
    data$prAA = pbinom(data$dp - data$ad, data$dp, epsilon/3, lower.tail = FALSE)

    homRR = data$prRR >= alpha
    homAA = data$prAA >= alpha
    data$het = !(homRR | homAA)

    if (plot)
    {
        plot_data = data.frame(vaf = jitter(data$ad, amount = 0.5) / data$dp, class = NA)
        plot_data$class[data$het] = "Het"
        plot_data$class[homRR] = "RR"
        plot_data$class[homAA] = "AA"
        print(ggplot(plot_data, aes(x = vaf, fill = class)) + geom_histogram(binwidth = 0.01) + xlab("VAF (jittered)") + xlim(0, 1))
    }

    data
}



main = function()
{
    sprintf(
'Find regions of subclonal aneuploidy in massively-parallel sequencing data.

Usage:
  soma-cnv.R [options] <affinity> <gc> <infile> <outfile>
  soma-cnv.R -h | --help
  soma-cnv.R --version

Parameters:
  <affinity>       Path to input affinity calibration file, tsv format with header, columns chrom, pos, affinity.
  <gc>             Path to input GC file, tsv format with header, columns chrom, pos, gc100, gc200, gc400, gc600, gc800.
  <infile>         Path to input allele depth file, tsv format with header, columns chrom, pos, dp, ad.
  <outfile>        Output RDS.

Options:
  -h --help        Show this message.
  --version        Show version.
  --epsilon=<F>    Per-base error rate [default: 0.02].
  --fixpval=<F>    P-value threshold to exclude fixed loci [default: 0.01].
  --segsize=<I>    Initial genome segment size, in loci [default: 100].
  --segpval=<F>    Segment merge P-value threshold [default: 0.01].
  --model=<S>      Error model to use, either poisson or nb [default: nb].
  --maxploidy=<I>  Maximum allele copy number to consider [default: 2].
  --diag=<P>       Path to PDF of diagnostic plots (if not specified, plots are not generated).

%s
Mark Pinese  <m.pinese@garvan.org.au>', VERSION) -> doc

    suppressPackageStartupMessages(library(docopt))

    opts = docopt(doc)

    if (opts$version)
    {
        cat(VERSION)
        cat("\n")
        quit(save = "no")
    }

    diag_mode = FALSE
    if (!is.null(opts$diag))
    {
        diag_mode = TRUE
        pdf(opts$diag, height = 8, width = 8)
    }

    opts$epsilon = as.numeric(opts$epsilon)
    opts$fixpval = as.numeric(opts$fixpval)
    opts$maxploidy = as.integer(opts$maxploidy)
    opts$model = match.arg(opts$model, choices = c("nb", "poisson"))
    opts$segsize = as.integer(opts$segsize)
    opts$segpval = as.numeric(opts$segpval)

    stopifnot(opts$epsilon >= 0 && opts$epsilon < 2/3)
    stopifnot(opts$fixpval > 0 && opts$fixpval < 1)
    stopifnot(opts$lvt < opts$uvt)
    stopifnot(opts$maxploidy >= 1)
    stopifnot(opts$segsize > 0)
    stopifnot(opts$segpval > 0 && opts$segpval < 1)

    message("Loading affinities... ", appendLF = FALSE)
    affinity = read.table(opts$affinity, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    affinity$affinity = affinity$affinity / mean(affinity$affinity)
    message(sprintf("%d loci read.", nrow(affinity)))
    message("Loading GC... ", appendLF = FALSE)
    gc = read.table(opts$gc, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    message(sprintf("%d loci read.", nrow(gc)))
    message("Loading AFs... ", appendLF = FALSE)
    data = read.table(opts$infile, header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = c("chrom", "pos", "dp", "ad"))
    n = nrow(data)
    message(sprintf("%d loci read.", n))

    message("Retaining only autosomes...", appendLF = FALSE)
    data$chrom = gsub("^chr", "", as.character(data$chrom))
    data = data[grepl("^[0-9]+$", data$chrom),]
    data$chrom = as.integer(data$chrom)
    message(sprintf("%d/%d (%.2f%%) loci remain", nrow(data), n, nrow(data)/n*100))
    n = nrow(data)

    message("Fitting depth anomaly...")
    data = fitDepthAnomaly(data, affinity, gc)
    message(sprintf("%d/%d (%.2f%%) loci remain", nrow(data), n, nrow(data)/n*100))
    n = nrow(data)

    message("Identifying probable fixed loci... ", appendLF = FALSE)
    data = markFixedLoci(data, epsilon = opts$epsilon, alpha = opts$fixpval, plot = diag_mode)
    message(sprintf("%d/%d (%.2f%%) loci likely not fixed", sum(data$het), n, mean(data$het)*100))

    data = data[order(data$chrom, data$pos),]       # Ensure correct ordering

    message("Fitting aneuploidy...")
    fit = fitWindows(data[data$het,], max_maxk = opts$maxploidy, family = opts$model, w = NULL, window.initial_size = opts$segsize, window.alpha = opts$segpval)

    if (diag_mode)
    {
        message("Generating plots")
        plotRawData(data)
        plotWindowFit(data, fit$fit)
        dev.off()
    }

    # Save results
    message("Writing results...")
    saveRDS(list(data = data, fit = fit, opts = opts), opts$outfile)

    message("Done.")
}


main()
