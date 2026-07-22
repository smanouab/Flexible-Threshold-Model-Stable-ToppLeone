library(FMStable)
library(statmod)

fB <- function(x, m, theta)
  m * theta * exp(-theta*x) * (1 - exp(-theta*x))^(m-1)
FB <- function(x, m, theta)
  (1 - exp(-theta*x))^m
qB <- function(u, m, theta)
  -(1/theta) * log(1 - u^(1/m))
fT <- function(x, al, ga, ze) 
  FMStable::dEstable(x, setParam(al,location=ze,logscale=log(ga),pm=0))

FT <- function(x, al, ga, ze)
  FMStable::pEstable(x, setParam(al,location=ze,logscale=log(ga),pm=0))

simulate_spliced <- function(n, params) {
  p=params[1]; m=params[2]; theta=params[3]
  al=params[4]; ga=params[5]; ze=params[6]
  g <- function(u){
    hB <- fB(u,m,theta)/FB(u,m,theta) 
    hT <- fT(u,al,ga,ze)/(1-FT(u,al,ga,ze))
    hB_safe <- ifelse(is.nan(hB),0,hB)
    p* hB_safe - (1-p)*hT
  }
  
  mode_B <- if(m > 1) log(m)/theta else 1/theta
  upper_u <- max(50, mode_B * 20, ze)
  uc <- uniroot(g, lower=0.01, upper=upper_u, tol=1e-10)$root
  U0 <- runif(n); U1 <- runif(n); X <- numeric(n)
  for(i in 1:n) {
    if(U0[i] <= p) {
      X[i] <- qB(U1[i]*FB(uc,m,theta), m, theta)
    } else {
      v    <- FT(uc,al,ga,ze) + U1[i]*(1-FT(uc,al,ga,ze))
      X[i] <- qEstable(v, setParam(al,location=ze,logscale=log(ga),pm=0))
    }
  }
  X
}

g_minus_1 <- function(x,alpha, beta){
  scalar_integrand <- function(t){
    # Changement de variable
    r <- t^(1/alpha)
    dr_dt <- (1/alpha) * t^(1/alpha - 1)
    
    trig_part <- cos(x*r - tan(pi*alpha/2)*beta*r^alpha) - 1
    
    val <- trig_part * r^(-2) * dr_dt
    
    return(val)
  }
  
  n <- 400
  gl <- statmod::gauss.quad(n, kind = "laguerre",alpha = 0)
  t <- gl$nodes
  w <- gl$weights
  
  sum(w * sapply(t, scalar_integrand))  
}

# Fonction de superquantile empirique


superquantile_emp <- function(x, p){
  q <- quantile(x, probs = p)
  mean(x[x > q])
} 

quantile_superquantile_estim <- function(Qs, SQs, verbose=F){
  
  objective_function <- function(params, Qs, SQs){
    al <- params[1]
    ga <- params[2]
    ze <- params[3]
    
    # On évite alpha <= 0 pour stabilité
    if(al <= 0 || ga <= 0 ) return(1e10)
    
    # éviter singularités
    if(abs(al - 1) < 1e-2) return(1e10)
    
    # éviter gamma poles
    if(abs(1/al - round(1/al)) < 1e-2) return(1e10)
    
    
    err <- 0
    for(i in 1:length(Qs)){
      Q <- Qs[i]
      P_T_Q <- suppressWarnings(
        pEstable(Q, setParam(alpha=al, location=ze, logscale=log(ga), pm=0),
                 lower.tail=FALSE))
      
      z <- (Q - ze + ga*tan(pi*al/2))/ga
      g1 <- g_minus_1(x=z, alpha=al, beta=0.999)
      
      H_val <- (ze - ga*tan(pi*al/2) - Q)/(2*P_T_Q) +
        (ga/(pi*P_T_Q)) * (gamma(1-1/al) - g1)
      
      err <- err + (SQs[i] - Q - H_val)^2
      
    }
    if (any(!is.finite(err))) return(1e10)
    
    return(err)
  }
  
  init_params <- c(al=1.5, ga= diff(range(Qs))/4, ze=median(Qs))
  
  res <- optim(par = init_params,
               fn = objective_function,
               Qs = Qs,
               SQs = SQs,
               method = "L-BFGS-B",
               lower = c(1.01, 0.01, 0),
               upper = c(1.999, Inf, Inf),
               control=list(maxit=500,factr   = 1e7,pgtol   = 1e-8))
  
  if(verbose){
    if(res$convergence==0) cat("Convergence atteinte à l'initialisation \n")
    else warning("Convergence non atteinte à l'initialisation\n")
  }
  
  return(res$par)
}

incomplete_moment_estim <- function(x_body,
                                    q_values = c(0.5, 1, 1.5, 2),
                                    prob_body = seq(0.1, 0.9, by=0.1),
                                    init_pars = NULL,
                                    verbose   = FALSE) {
  
  IM_q_numint <- function(Qp, m, theta, q) {
    if(Qp <= 0 || !is.finite(Qp)) return(NA_real_)
    integrand <- function(x) {
      vals <- x^q * fB(x, m, theta)
      vals[!is.finite(vals)] <- 0
      vals
    }
    lower_lim <- if(m < 1) 1e-10 else 0
    tryCatch({
      num   <- integrate(integrand, lower=lower_lim, upper=Qp,
                         subdivisions=500, rel.tol=1e-6)$value
      denom <- FB(Qp, m, theta)
      if(!is.finite(num) || denom < 1e-300) NA_real_ else num/denom
    }, error=function(e) NA_real_)
  }
  
  # Quantiles des observations corps directement
  Q_body_emp <- as.numeric(quantile(x_body, probs=prob_body))
  
  # Moments incomplets empiriques sur observations corps
  Imq_emp <- unlist(lapply(q_values, function(q)
    sapply(Q_body_emp, function(Q)
      mean(x_body[x_body <= Q]^q)
    )
  ))
  
  # Vérification dimensions
  expected <- length(q_values) * length(Q_body_emp)
  stopifnot(length(Imq_emp) == expected)
  
  # Objectif : erreur relative au carré (invariance d'échelle)
  objective <- function(pars) {
    m_    <- pars[1]; theta_ <- pars[2]
    if(m_ <= 0 || theta_ <= 0) return(1e10)
    
    err <- 0; idx <- 1
    for(q in q_values) {
      for(Q in Q_body_emp) {
        th_val  <- IM_q_numint(Q, m_, theta_, q)
        emp_val <- Imq_emp[idx]
        if(!is.na(th_val) && is.finite(th_val) &&
           is.finite(emp_val) && abs(emp_val) > 1e-300) {
          err <- err + ((emp_val - th_val) / emp_val)^2
        }
        idx <- idx + 1
      }
    }
    if(!is.finite(err)) 1e10 else err
  }
  
  # Init data-driven si non fournie
  if(is.null(init_pars))
    init_pars <- c(2,1.5)
  
  if(verbose)
    cat(sprintf("Init IM-matching : m=%.3f, theta=%.6f\n",
                init_pars["m"], init_pars["theta"]))
  
  opt <- tryCatch(
    optim(init_pars, objective, method="L-BFGS-B",
          lower=c(0.01, 1e-7), upper=c(Inf, Inf),
          control=list(maxit=2000, factr=1e5, pgtol=1e-9)),
    error=function(e) list(par=init_pars, value=Inf)
  )
  
  names(opt$par) <- c("m","theta")
  opt$par
}

init_body <- function(x, p_init) {
  
  x_body <- sort(x[x <= quantile(x, p_init)])
  
  theta_B <- incomplete_moment_estim(
    x_body  = x_body,
    q_values = c(0.5, 1, 1.5, 2),
    prob_body = seq(0.1, 0.9, by=0.1),
  )
  
  theta_B
}

Em_init <- function(x, p_init = 0.5,verbose=F){
  p <- c(p = p_init)
  
  prob_levels <- p_init + (1 - p_init) * seq(0.1, 0.9, by = 0.01)  # q > p_init
  Q_hat  <- sapply(prob_levels, function(pp) quantile(x, pp))
  SQ_hat <- sapply(prob_levels, function(pp) superquantile_emp(x, pp))
  
  theta_T <- quantile_superquantile_estim(Q_hat, SQ_hat,verbose = verbose)
  theta_B <- init_body(x, p_init)
  if(verbose) cat("Paramètres initiales estimés: ", c(p, theta_B, theta_T),"\n")
  c(p, theta_B, theta_T)
}

# ------------------------------------------------------------------
# EM Algorithme
# ------------------------------------------------------------------

EmAlgo <- function(x, max_iter = 1500, tol = 1e-6, verbose = FALSE) {
  
  n <- length(x)
  
  # ── Fonctions de base ──────────────────────────────────────────
  fB_loc <- function(u, m, theta)
    m * theta * exp(-theta*u) * (1 - exp(-theta*u))^(m-1)
  FB_loc <- function(u, m, theta)
    (1 - exp(-theta*u))^m
  fT_loc <- function(u, alpha, gam, zeta)
    FMStable::dEstable(u, setParam(alpha=alpha, location=zeta,
                                   logscale=log(gam), pm=0))
  FT_loc <- function(u, alpha, gam, zeta)
    FMStable::pEstable(u, setParam(alpha=alpha, location=zeta,
                                   logscale=log(gam), pm=0))
  
  hB_fn <- function(u, m, theta)
    fB_loc(u,m,theta) / pmax(FB_loc(u,m,theta), 1e-300)
  hT_fn <- function(u, alpha, gam, zeta)
    fT_loc(u,alpha,gam,zeta) / pmax(1 - FT_loc(u,alpha,gam,zeta), 1e-300)
  
  profile_p <- function(u, m, theta, alpha, gam, zeta) {
    hb <- tryCatch(hB_fn(u,m,theta),        error=function(e) NA_real_)
    ht <- tryCatch(hT_fn(u,alpha,gam,zeta), error=function(e) NA_real_)
    if(!is.finite(hb) || !is.finite(ht) || (hb+ht) < 1e-300) return(NA_real_)
    ht / (hb + ht)
  }
  
  loglik_spliced_true <- function(pars) {
    u_ <- pars[1]; m_ <- pars[2]; th_ <- pars[3]
    al_ <- pars[4]; ga_ <- pars[5]; ze_ <- pars[6]
    
    if(u_<=0 || m_<=0 || th_<=0 || al_<=1 || al_>=2 || ga_<=0)
      return(-Inf)
    
    p_ <- mean(x <=u_)
    if(is.na(p_) || p_<=0 || p_>=1) return(-Inf)
    
    FB_u <- max(FB_loc(u_, m_, th_), 1e-300)
    ST_u <- max(1 - FT_loc(u_, al_, ga_, ze_), 1e-300)
    
    f_x  <- numeric(n)
    idx1 <- which(x <=  u_)
    idx2 <- which(x > u_)
    
    if(length(idx1) > 0) {
      fv <- fB_loc(x[idx1], m_, th_)
      f_x[idx1] <- p_ * fv / FB_u
    }
    if(length(idx2) > 0) {
      fv <- fT_loc(x[idx2], al_, ga_, ze_)
      f_x[idx2] <- (1-p_) * fv / ST_u
    }
    
    ll <- sum(log(pmax(f_x, 1e-300)))
    if(!is.finite(ll)) -Inf else ll
  }
  
  loglik_spliced <- function(pars, eps=1e-9) {
    u_ <- pars[1]; m_ <- pars[2]; th_ <- pars[3]
    al_ <- pars[4]; ga_ <- pars[5]; ze_ <- pars[6]
    
    if(u_<=0 || m_<=0 || th_<=0 || al_<=1 || al_>=2 || ga_<=0)
      return(-Inf)
    
    p_ <- profile_p(u_, m_, th_, al_, ga_, ze_)
    if (is.na(p_)) return(-Inf)
    if (p_ <= eps || p_ >= 1 - eps) {
      penalty <- -1e8 * (1 + abs(log(pmin(pmax(p_, eps/10), 1-eps/10))))
      return(penalty)
    }
    
    FB_u <- max(FB_loc(u_, m_, th_), 1e-300)
    ST_u <- max(1 - FT_loc(u_, al_, ga_, ze_), 1e-300)
    
    f_x  <- numeric(n)
    idx1 <- which(x <=  u_)
    idx2 <- which(x > u_)
    
    if(length(idx1) > 0) {
      fv <- fB_loc(x[idx1], m_, th_)
      f_x[idx1] <- p_ * fv / FB_u
    }
    if(length(idx2) > 0) {
      fv <- fT_loc(x[idx2], al_, ga_, ze_)
      f_x[idx2] <- (1-p_) * fv / ST_u
    }
    
    ll <- sum(log(pmax(f_x, 1e-300)))
    if(!is.finite(ll)) -Inf else ll
  }
  
  neg_Q <- function(pars, x_body, x_tail) {
    u_ <- pars[1]; m_ <- pars[2]; th_ <- pars[3]
    al_ <- pars[4]; ga_ <- pars[5]; ze_ <- pars[6]
    
    if(u_<=0 || m_<=0 || th_<=0 || al_<=1 || al_>=2 || ga_<=0)
      return(1e10)
    
    p_ <- profile_p(u_, m_, th_, al_, ga_, ze_)
    if(is.na(p_) || p_<=0 || p_>=1) return(1e10)
    
    n1 <- length(x_body); n2 <- length(x_tail)
    FB_u <- max(FB_loc(u_, m_, th_), 1e-300)
    ST_u <- max(1 - FT_loc(u_, al_, ga_, ze_), 1e-300)
    
    Q1 <- n1*log(p_) + sum(log(pmax(fB_loc(x_body,m_,th_),1e-300))) - n1*log(FB_u)
    Q2 <- n2*log(1-p_) + sum(log(pmax(fT_loc(x_tail,al_,ga_,ze_),1e-300))) - n2*log(ST_u)
    
    val <- Q1 + Q2
    if(!is.finite(val)) 1e10 else -val
  }
  
  # ── Initialisation ─────────────
  find_ustar_local <- function(p, m, theta, alpha, gam, zeta, x_data=x) {
    g <- function(u) {
      hb <- tryCatch(hB_fn(u,m,theta),        error=function(e) NA_real_)
      ht <- tryCatch(hT_fn(u,alpha,gam,zeta), error=function(e) NA_real_)
      if(!is.finite(hb) || !is.finite(ht)) return(NA_real_)
      p*hb - (1-p)*ht
    }
    
    lo <- as.numeric(quantile(x_data, 0.01))
    hi <- as.numeric(quantile(x_data, 0.99))
    
    g_lo <- tryCatch(g(lo), error=function(e) NA_real_)
    g_hi <- tryCatch(g(hi), error=function(e) NA_real_)
    
    if (!is.finite(g_lo) || !is.finite(g_hi) || sign(g_lo) == sign(g_hi)) {
      return(NA_real_)
    }
    
    tryCatch(
      uniroot(g, lower=lo, upper=hi, tol=1e-10)$root,
      error = function(e) NA_real_
    )
  }
  
  p_grids <- seq(0.1, 0.9, 0.2)
  init_log_lik <- numeric(length(p_grids))
  init_pars <- vector("list", length(p_grids))
  
  for(i in seq_along(p_grids)){
    if(verbose) cat(sprintf("Initialisation pour p=%.4f\n", p_grids[i]))
    init <- Em_init(x, p_init = p_grids[i], verbose=verbose)
    
    ustar <- find_ustar_local(init[1], init[2], init[3],
                              init[4], init[5], init[6], x_data = x)
    
    if (is.na(ustar)) {
      if(verbose) cat(sprintf("  -> p_init=%.2f rejeté : pas de u* valide\n", p_grids[i]))
      init_log_lik[i] <- -Inf
      init_pars[[i]] <- NULL
      next
    }
    
    init_pars[[i]] <- c(ustar, init[2:6])
    init_log_lik[i] <- loglik_spliced(init_pars[[i]])
    
    if(verbose)
      cat(sprintf("  -> u*=%.4f | CL(profile_p)=%.4f\n", ustar, init_log_lik[i]))
  }
  
  if (all(!is.finite(init_log_lik))) {
    stop("Aucun candidat d'initialisation n'admet de seuil de continuité valide. ",
         "Vérifier l'adéquation du modèle aux données (cf. diagnostic Cause B).")
  }
  
  best_init <- init_pars[[which.max(init_log_lik)]]
  pars <- best_init
  names(pars) <- c("ustar","m","theta","alpha","gamma","zeta")
  
  u_lo    <- as.numeric(quantile(x, 0.01))
  u_hi    <- as.numeric(quantile(x, 0.99))
  lower_b <- c(u_lo, 1e-8, 1e-8, 1.001, 1e-3, 0)
  upper_b <- c(u_hi, 1e9,   1e9,   1.999, 1e9,   1e9)
  
  # ── Historiques ───────────────────────────────────────────────
  CL_history   <- numeric(max_iter)  # CL(η^(t)) = ℓ(η^(t);x)
  z_history    <- vector("list", max_iter)
  
  # CL initiale avant toute itération
  CL_prev <- loglik_spliced(pars)
  CL_best <- CL_prev
  pars_best <- pars
  
  if(verbose)
    cat(sprintf("Init | ustar=%.4f p=%.4f | CL0=%.4f\n",
                pars["ustar"],
                mean(x <=pars["ustar"]),
                CL_prev))
  
  converged <- FALSE
  
  for(iter in 1:max_iter) {
    
    # ── Étape E : assignation hard ─────────────────────────────
    u_cur  <- pars["ustar"]
    x_body <- x[x <=  u_cur]
    x_tail <- x[x > u_cur]
    z_cur  <- x <= u_cur          # partition courante 
    n1 <- length(x_body); n2 <- length(x_tail)
    
    if(n1 < 5 || n2 < 5) {
      warning(sprintf("Iter %d : partition dégénérée (n1=%d,n2=%d).", iter,n1,n2))
      break
    }
    
    # ── Étape M : maximisation jointe, partition FIXÉE ─────────
    opt <- tryCatch(
      nlminb(pars, neg_Q, x_body=x_body, x_tail=x_tail,
             lower=lower_b, upper=upper_b,
             control=list(iter.max=1000, rel.tol=1e-10)),
      error = function(e) {
        if(verbose) cat(sprintf("  [iter %d] nlminb -> Nelder-Mead\n", iter))
        r <- optim(pars, neg_Q, x_body=x_body, x_tail=x_tail,
                   method="Nelder-Mead", control=list(maxit=3000, reltol=1e-8))
        list(par=r$par, objective=r$value, convergence=r$convergence)
      }
    )
    
    pars <- opt$par
    names(pars) <- c("ustar","m","theta","alpha","gamma","zeta")
    
    # On recalcule la partition induite par le NOUVEAU η,
    # puis on évalue la log-vraisemblance composite — c'est CL(η^(t+1)).
    # C'est la quantité dont la différence avec CL(η^(t)) détermine la convergence.
    CL_new <- loglik_spliced(pars)
    CL_history[iter] <- CL_new
    z_new <- x <= pars["ustar"]     # nouvelle partition
    z_history[[iter]] <- z_new
    
    # Mise à jour du meilleur point visité
    if(is.finite(CL_new) && CL_new > CL_best) {
      CL_best   <- CL_new
      pars_best <- pars
    }
    
    p_cur <- profile_p(pars["ustar"],pars["m"],pars["theta"],
                       pars["alpha"],pars["gamma"],pars["zeta"])
    #p_cur <- mean(x<=pars["ustar"])
    
    if(verbose)
      cat(sprintf(
        "Iter %3d | u*=%8.3f p=%.4f | CL(t)=%12.4f CL(t+1)=%12.4f | ΔCL=%.2e\n",
        iter, pars["ustar"], p_cur, CL_prev, CL_new, CL_new - CL_prev))
    
    # ── Critère d'arrêt ────────────
    
    if(abs(CL_new - CL_prev) < tol) {
      cat(sprintf("Convergence (point fixe) — iter %d | CL=%.6f\n",
                  iter, CL_new))
      converged <- TRUE
      break
    }
    
    if(iter == max_iter) warning("max_iter atteint sans convergence.")
    
    CL_prev <- CL_new
  }
  
  p_final <- mean(x < pars_best["ustar"])
  CL_final <- loglik_spliced_true(pars_best)
  
  list(
    estimates    = pars_best,
    p        = as.numeric(p_final),
    CL_final     = CL_final,
    n_iter       = iter,
    CL_history   = CL_history[1:iter],
    converged    = converged
  )
}

#=================================================================================
#========================= Simulation ============================================
#=================================================================================
simulation <- function(n, pars, replications=50){
  ## params est la liste de paramètres nommés p,m,theta,alpha,gamma,delta
  results <- matrix(NA, nrow = replications, ncol = 6)
  
  for (i in 1:replications) {
    sample <- simulate_spliced(n,pars)
    estimation <- EmAlgo(sample,tol=1e-4,verbose=F)
    if (!is.null(estimation)) {
      pars_estimated <- estimation$estimates
      p_estimated <- estimation$p
      results[i, ] <- c(p_estimated,pars_estimated["m"],
                        pars_estimated["theta"],pars_estimated["alpha"],
                        pars_estimated["gamma"], pars_estimated["zeta"])
    }
    else results[i, ] <- c(NA,NA,NA,NA,NA,NA)
  }
  colnames(results) <- c( "p", "m", "theta","alpha","gamma", "zeta")
  # Calcul des métriques
  true_params <- pars[c("p","m","theta","alpha","gamma","zeta")]
  
  true_matrix <- matrix(true_params,nrow = replications,ncol = 6,byrow = TRUE)
  
  mean_estimates <- colMeans(results, na.rm = TRUE)
  biases <- abs(mean_estimates - true_params)
  mse <- colMeans((results - true_matrix)^2, na.rm = TRUE)
  variances <- apply(results, 2, var, na.rm = TRUE)
  
  return(data.frame(Mean_Estimate = mean_estimates,
                    Bias = biases,
                    MSE = mse,
                    Variance = variances))
  
}


# Application of simulation

pars_1 = c(p=0.3, m=10,theta=4,alpha=1.5,gamma=1,zeta = 100)
mc1 <- simulation(250,pars_1)
mc1
mc2 <- simulation(1e3,pars_1)
mc2

pars_2 = c(p=0.8, m=2,theta=0.2,alpha=1.8,gamma=100,zeta = 1000)
mc3 <- simulation(250,pars_2)
mc3


mc4 <- simulation(1e3,pars_2)
mc4

