library(tidyverse)
library(survival)
library(goftest)
library(FMStable)

listing_data <- read_csv("TPLTLEx all codes/dataset/listing_data_tidy.csv")
midline_data <- read_csv('TPLTLEx all codes/dataset/midline_data_tidy.csv')
endline_data <- read_csv('TPLTLEx all codes/dataset/endline_data_tidy.csv')

listing_income = listing_data$lis_e34

list_estim <- EmAlgo(listing_income,verbose = T)
list_estim

mid_inc = midline_data$mid_g17
mid_inc <- mid_inc[mid_inc > 0]

mid_estim <- EmAlgo(mid_inc_no_zero,verbose = T)
mid_estim

end_inc = endline_data$end_G17
end_inc <- end_inc[end_inc > 0]

end_estim <- EmAlgo(end_inc,verbose = T)
end_estim


#==========================================================================
#Performance=========================================================
#==========================================================================
ptlestable <- function(x,ustar,p,m,theta,alpha,gamma,delta,
                       lower.tail = TRUE, log.p = FALSE){
  FB <- function(u, m, theta)
    (1 - exp(-theta*u))^m
  FT <- function(u, alpha, gam, zeta)
    FMStable::pEstable(u, setParam(alpha=alpha, location=zeta,
                                   logscale=log(gam), pm=0))
  FB_u <- pmax(FB(ustar, m, theta), 1e-300)
  ST_u <- pmax(1 - FT(ustar, alpha, gamma, delta), 1e-300)
  
  x    <- as.numeric(x)
  cdf  <- numeric(length(x))
  idx1 <- which(x <= ustar)
  idx2 <- which(x >  ustar)
  if (length(idx1) > 0)
    cdf[idx1] <- p * FB(x[idx1], m, theta) / FB_u
  
  if (length(idx2) > 0) {
    FT_x <- FT(x[idx2], alpha, gamma, delta)
    FT_u <- 1 - ST_u
    cdf[idx2] <- p + (1 - p) * (FT_x - FT_u) / ST_u
  }
  
  if (!lower.tail) cdf <- 1 - cdf
  if (log.p)       cdf <- log(cdf)
  
  cdf
}

compute_criteria <- function(logLik, b, n) {
  # b = nbre de paramètres
  AIC  <- -2 * logLik + 2 * b
  BIC  <- -2 * logLik + log(n) * b
  AICc <- AIC + (2*b*(b+1))/(n-b-1)
  HQIC <- -2 * logLik + 2 * log(log(n)) * b
  
  return(c(AIC,BIC,AICc,HQIC))
}

size_list = length(listing_income)
size_mid = length(mid_inc)
size_end = length(end_inc)

list_perf <- compute_criteria(list_estim$CL_final,6,size_list)
list_perf
mid_perf <- compute_criteria(mid_estim$CL_final,6,size_mid)
mid_perf
end_perf <- compute_criteria(end_estim$CL_final,6,size_list)
end_perf

ks_list = ks.test(listing_income,'ptlestable',ustar=6.114245e+04,p = 0.7972747,m = 3.906735,
                           theta=6.372449e-05, alpha=1.5149,gamma=2.806893e+04,delta=4.577066e+04,simulate.p.value = T,B=1000)
ad_list = goftest::ad.test(listing_income,'ptlestable',ustar=6.114245e+04,p = 0.7972747,m = 3.906735,
                           theta=6.372449e-05, alpha=1.5149,gamma=2.806893e+04,delta=4.577066e+04,estimated = T)


cat(sprintf("KS : D = %.5f, p-value = %.4f\n", ks_list$statistic, ks_list$p.value))
cat(sprintf("AD : A = %.5f, p-value = %.4f\n", ad_list$statistic, ad_list$p.value))

ks_mid = ks.test(mid_inc,'ptlestable',ustar=5.955026e+04,p = 0.6831607,m = 1.543060,
                  theta=2.932062e-05, alpha=1.006398,gamma=5.069012e-02,delta=1.834412e+04,simulate.p.value = T,B=1000)
ad_mid = goftest::ad.test(mid_inc,'ptlestable',ustar=5.955026e+04,p = 0.6831607,m = 1.543060,
                          theta=2.932062e-05, alpha=1.006,gamma=5.069012e-02,delta=1.834412e+04,estimated = T)


cat(sprintf("KS : D = %.5f, p-value = %.4f\n", ks_mid$statistic, ks_mid$p.value))
cat(sprintf("AD : A = %.5f, p-value = %.4f\n", ad_mid$statistic, ad_mid$p.value))

ks_end = ks.test(end_inc,'ptlestable',ustar=5.539819e+04,p = 0.652234,m = 1.446717,
                          theta=2.239215e-05, alpha=1.001,gamma=1.000000e-03,delta=2.136853e+04 ,simulate.p.value = T,B=1000)
ad_end = goftest::ad.test(end_inc,'ptlestable',ustar=5.539819e+04,p = 0.652234,m = 1.446717,
                          theta=2.239215e-05, alpha=1.001,gamma=1.000000e-03,delta=2.136853e+04 ,estimated = T)

cat(sprintf("AD : KS = %.5f, p-value = %.4f\n", ks_end$statistic, ks_end$p.value))
cat(sprintf("AD : A = %.5f, p-value = %.4f\n", ad_end$statistic, ad_end$p.value))


#==================================================================================
#Diagnostic graphique==============================================================
#==================================================================================

qtlestable <- function(probs, ustar, p, m, theta, alpha, gamma, delta) {
  
  FB_u <- (1 - exp(-theta * ustar))^m
  FT_u <- FMStable::pEstable(ustar,
                             FMStable::setParam(alpha=alpha, location=delta,
                                                logscale=log(gamma), pm=0))
  ST_u <- 1 - FT_u
  
  sapply(probs, function(u) {
    if (u <= 0)   return(0)
    if (u >= 1)   return(Inf)
    if (u <= p) {
      # Branche corps : inversion analytique de QB
      v <- u * FB_u / p
      v <- pmin(pmax(v, 1e-15), 1 - 1e-15)
      -(1/theta) * log(1 - v^(1/m))          # qB fermé pour TLE
    } else {
      # Branche queue : inversion de QT via qEstable
      v <- FT_u + (u - p) * ST_u / (1 - p)
      v <- pmin(pmax(v, 1e-15), 1 - 1e-15)
      suppressWarnings(
        FMStable::qEstable(v,
                           FMStable::setParam(alpha=alpha, location=delta,
                                              logscale=log(gamma), pm=0))
      )
    }
  })
}

diagnostic_plots <- function(x, pars,
                             conf_level = 0.95) {
  
  n  <- length(x)
  xs <- sort(x)
  
  # Probabilités empiriques (Hazen, évite 0 et 1)
  p_emp <- (seq_len(n) - 0.5) / n
  
  # CDF théorique aux points empiriques (pour PP-plot) 
  p_theo <- ptlestable(xs,
                       ustar = pars$ustar, p     = pars$p,
                       m     = pars$m,     theta = pars$theta,
                       alpha = pars$alpha, gamma = pars$gamma,
                       delta = pars$delta)
  
  # ── Quantiles théoriques (pour QQ-plot) ───────────────────────
  q_theo <- qtlestable(p_emp,
                       ustar = pars$ustar, p     = pars$p,
                       m     = pars$m,     theta = pars$theta,
                       alpha = pars$alpha, gamma = pars$gamma,
                       delta = pars$delta)
  
  # ── Bandes de confiance Kolmogorov-Smirnov (niveau conf_level) ─
  # Basées sur l'approximation de la loi limite de l'ECDF
  epsilon   <- sqrt(-log((1 - conf_level) / 2) / (2 * n))
  p_lo      <- pmax(p_emp - epsilon, 0)
  p_hi      <- pmin(p_emp + epsilon, 1)
  
  # Bandes QQ : inverser les bornes KS
  q_lo <- qtlestable(p_lo,
                     ustar=pars$ustar, p=pars$p, m=pars$m,
                     theta=pars$theta, alpha=pars$alpha,
                     gamma=pars$gamma, delta=pars$delta)
  q_hi <- qtlestable(p_hi,
                     ustar=pars$ustar, p=pars$p, m=pars$m,
                     theta=pars$theta, alpha=pars$alpha,
                     gamma=pars$gamma, delta=pars$delta)
  
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5))
  
  # ── PP-plot ────────────────────────────────────────────────────
  plot(p_emp, p_theo,
       xlab = "Empirical probability",
       ylab = "Theoretical probability",
       pch  = 16, cex = 0.5, col = "steelblue",
       xlim = c(0, 1), ylim = c(0, 1), asp = 1)
  polygon(c(0, 1, 1, 0),
          c(epsilon, 1, 1 - epsilon, 0),
          col = adjustcolor("gray70", alpha.f = 0.25),
          border = NA)
  abline(0, 1, col = "red", lwd = 1.5)
  # ── QQ-plot ────────────────────────────────────────────────────
  ok <- is.finite(q_theo) & is.finite(q_lo) & is.finite(q_hi)
  
  plot(q_theo[ok], xs[ok],
       xlab = "Theoretical quantiles",
       ylab = "Empirical quantiles",
       pch  = 16, cex = 0.5, col = "steelblue")
  polygon(c(q_lo[ok], rev(q_hi[ok])),
          c(xs[ok],   rev(xs[ok])),
          col    = adjustcolor("gray70", alpha.f = 0.25),
          border = NA)
  abline(0, 1, col = "red", lwd = 1.5)
  par(mfrow = c(1, 1))
}


params <- list(
  listing = list(ustar=61142.45, p=0.7972747, m=3.906735,
                 theta=6.372449e-05, alpha=1.5149, gamma=28068.93, delta=45770.66),
  mid     = list(ustar=5.955026e+04,p = 0.6831607,m = 1.543060,
                 theta=2.932062e-05, alpha=1.006,gamma=5.069012e-02,delta=1.834412e+04),
  end     = list(ustar=5.539819e+04,p = 0.652234,m = 1.446717,
                 theta=2.239215e-05, alpha=1.001,gamma=1.000000e-03,delta=2.136853e+04)
)

datasets <- list(
  listing = listing_income,
  mid     = mid_inc,
  end     = end_inc
)

## Enregistrement des figures
for (nm in names(datasets)) {
  
  png(
    filename = paste0("diagnostic_", nm, ".png"),
    width = 1800,
    height = 900,
    res = 300
  )
  
  diagnostic_plots(
    x     = datasets[[nm]],
    pars  = params[[nm]]
  )
  
  dev.off()
}
