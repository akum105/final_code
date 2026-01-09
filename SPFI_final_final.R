## =========================================================
##  SPFI (= SFI here) with GMM  --- Simulation 1
##  Output:
##    Figure 1 (M1/M2): Full vs SPFI boxplots
##    Figure 2 (M3/M4): Full vs SPFI boxplots
##    Figure 3: Histogram of selected G (FULL vs SPFI)
##
##  CHANGE (requested):
##    gmm_em_missing() initialization no longer uses mean-impute + kmeans.
##    Instead, it uses COMPLETE-CASE kmeans to initialize mu_g (and alpha_g),
##    with a safe fallback when complete cases are insufficient.
## =========================================================

suppressPackageStartupMessages({
  library(mvtnorm)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(parallel)
})

## -------------------------
## helpers (numerical)
## -------------------------
expit <- function(x) 1 / (1 + exp(-x))

make_pd <- function(S, eps = 1e-10, max_tries = 12) {
  S <- (S + t(S)) / 2
  if (nrow(S) == 1L) {
    if (!is.finite(S[1,1])) S[1,1] <- eps
    if (S[1,1] <= 0) S[1,1] <- eps
    return(S)
  }
  jitter <- eps
  for (k in 0:max_tries) {
    out <- tryCatch(chol(S + diag(jitter, nrow(S))), error = function(e) NULL)
    if (!is.null(out)) return(S + diag(jitter, nrow(S)))
    jitter <- jitter * 10
  }
  ev <- eigen(S, symmetric = TRUE)
  vals <- pmax(ev$values, eps)
  S2 <- ev$vectors %*% diag(vals) %*% t(ev$vectors)
  (S2 + t(S2)) / 2
}

safe_solve <- function(A) {
  A <- make_pd(A)
  solve(A)
}

## MVN observed-marginal loglik with missing
log_mvn_marginal <- function(y, mu, Sigma) {
  obs <- which(!is.na(y))
  yobs <- y[obs]
  muobs <- mu[obs]
  S_oo <- Sigma[obs, obs, drop = FALSE]
  S_oo <- make_pd(S_oo)
  dmvnorm(yobs, mean = muobs, sigma = S_oo, log = TRUE)
}

## -------------------------
## Data generation (Simulation 1)
## -------------------------
Sigma_AR <- function(rho) {
  matrix(c(1, rho, rho^2,
           rho, 1, rho,
           rho^2, rho, 1), nrow = 3, byrow = TRUE)
}

gen_data <- function(model = c("M1","M2","M3","M4"), n = 500) {
  model <- match.arg(model)
  
  if (model %in% c("M1","M2")) {
    alpha <- c(0.3, 0.3, 0.4)
    mu_list <- list(c(-3,-3,-3), c(1,1,1), c(5,5,5))
    Sigma <- Sigma_AR(0.7)
    
    z <- sample(1:3, size = n, replace = TRUE, prob = alpha)
    Y <- matrix(NA_real_, n, 3)
    
    for (g in 1:3) {
      idx <- which(z == g)
      if (length(idx) == 0) next
      if (model == "M2" && g == 2) {
        Y[idx, ] <- matrix(rexp(length(idx) * 3, rate = 1), ncol = 3)
      } else {
        Y[idx, ] <- rmvnorm(length(idx), mean = mu_list[[g]], sigma = Sigma)
      }
    }
    colnames(Y) <- c("Y1","Y2","Y3")
    return(Y)
  }
  
  if (model == "M3") {
    e1 <- rnorm(n, mean = 0, sd = 1)
    e2 <- rgamma(n, shape = 1, rate = 1)  # Exp(1)
    e3 <- rchisq(n, df = 1)
    Y1 <- 1 + e1
    Y2 <- 0.5 * Y1 + e2
    Y3 <- Y2 + e3
    Y <- cbind(Y1,Y2,Y3)
    colnames(Y) <- c("Y1","Y2","Y3")
    return(Y)
  }
  
  if (model == "M4") {
    mu12 <- c(1,2)
    Sig12 <- matrix(c(1,0.5,0.5,1), 2, 2)
    Y12 <- rmvnorm(n, mean = mu12, sigma = Sig12)
    Y1 <- Y12[,1]
    Y2 <- Y12[,2]
    Y3 <- Y2^2 + rnorm(n, mean = 0, sd = 1)
    Y <- cbind(Y1,Y2,Y3)
    colnames(Y) <- c("Y1","Y2","Y3")
    return(Y)
  }
  
  stop("Unknown model")
}

impose_missing <- function(Y) {
  n <- nrow(Y)
  pi2 <- expit(-0.8 + 0.4 * Y[, "Y1"])
  pi3 <- expit( 0.4 - 0.8 * Y[, "Y1"])
  m2 <- rbinom(n, size = 1, prob = pi2)
  m3 <- rbinom(n, size = 1, prob = pi3)
  
  Ymis <- Y
  Ymis[m2 == 1, "Y2"] <- NA_real_
  Ymis[m3 == 1, "Y3"] <- NA_real_
  Ymis
}

get_cuts <- function(model) {
  if (model %in% c("M1","M2")) return(c(c2 = -2, c3 = -2))
  if (model == "M3") return(c(c2 = 2, c3 = 3))
  if (model == "M4") return(c(c2 = 2, c3 = 5))
  stop("unknown model")
}

## -------------------------
## TRUE PARAMETERS (Exact / 1D quadrature)
## -------------------------
pexgauss <- function(x, mu, sigma, tau) {
  if (tau <= 0) stop("tau must be > 0")
  z  <- (x - mu) / sigma
  z2 <- z - sigma / tau
  logPhi_z2 <- pnorm(z2, log.p = TRUE)
  a <- -(x - mu) / tau + (sigma^2) / (2 * tau^2)
  term2 <- exp(a + logPhi_z2)
  F <- pnorm(z) - term2
  F <- pmin(pmax(F, 0), 1)
  bad <- !is.finite(F)
  if (any(bad)) F[bad] <- ifelse(x[bad] < mu, 0, 1)
  F
}

true_params_exact <- function(model, rel.tol = 1e-10, abs.tol = 0) {
  model <- match.arg(model, c("M1","M2","M3","M4"))
  cuts <- get_cuts(model)
  c2 <- as.numeric(cuts["c2"]); c3 <- as.numeric(cuts["c3"])
  
  if (model == "M1") {
    alpha <- c(0.3, 0.3, 0.4)
    mu <- c(-3, 1, 5)
    theta <- sum(alpha * mu)
    P2 <- sum(alpha * pnorm(c2, mean = mu, sd = 1))
    P3 <- sum(alpha * pnorm(c3, mean = mu, sd = 1))
    return(c(theta2 = theta, theta3 = theta, P2 = P2, P3 = P3))
  }
  
  if (model == "M2") {
    alpha <- c(0.3, 0.3, 0.4)
    mu_norm <- c(-3, 5)
    theta <- sum(alpha * c(-3, 1, 5))  # = 1.4
    
    P2 <- alpha[1] * pnorm(c2, mean = mu_norm[1], sd = 1) +
      alpha[2] * pexp(c2, rate = 1) +
      alpha[3] * pnorm(c2, mean = mu_norm[2], sd = 1)
    P3 <- alpha[1] * pnorm(c3, mean = mu_norm[1], sd = 1) +
      alpha[2] * pexp(c3, rate = 1) +
      alpha[3] * pnorm(c3, mean = mu_norm[2], sd = 1)
    
    return(c(theta2 = theta, theta3 = theta, P2 = P2, P3 = P3))
  }
  
  if (model == "M3") {
    mu <- 0.5; sigma <- 0.5; tau <- 1
    theta2 <- mu + tau
    theta3 <- theta2 + 1
    P2 <- pexgauss(c2, mu, sigma, tau)
    
    integrand <- function(t) {
      val <- pexgauss(c3 - t, mu, sigma, tau) * dchisq(t, df = 1)
      val[!is.finite(val)] <- 0
      val
    }
    P3 <- integrate(integrand, lower = 0, upper = Inf,
                    rel.tol = rel.tol, abs.tol = abs.tol)$value
    
    return(c(theta2 = theta2, theta3 = theta3, P2 = P2, P3 = P3))
  }
  
  if (model == "M4") {
    theta2 <- 2
    theta3 <- 1 + 2^2
    P2 <- pnorm(c2, mean = 2, sd = 1)
    
    integrand <- function(y) {
      val <- pnorm(c3 - y^2) * dnorm(y, mean = 2, sd = 1)
      val[!is.finite(val)] <- 0
      val
    }
    P3 <- integrate(integrand, lower = -Inf, upper = Inf,
                    rel.tol = rel.tol, abs.tol = abs.tol)$value
    
    return(c(theta2 = theta2, theta3 = theta3, P2 = P2, P3 = P3))
  }
  
  stop("unknown model")
}

## -------------------------
## GMM-EM with missing (shared covariance)
##   INIT MODIFIED: complete-case kmeans initialization
## -------------------------
gmm_em_missing <- function(Y, G,
                           maxit = 200, tol = 1e-6,
                           ridge = 1e-6, nstart_kmeans = 20,
                           verbose = FALSE) {
  n <- nrow(Y); p <- ncol(Y)
  stopifnot(p == 3)
  
  ## missingness patterns (also used for EM)
  m2 <- is.na(Y[,2])
  m3 <- is.na(Y[,3])
  
  idx_cc  <- which(!m2 & !m3)   # obs: 1,2,3
  idx_m2  <- which( m2 & !m3)   # obs: 1,3 ; mis: 2
  idx_m3  <- which(!m2 &  m3)   # obs: 1,2 ; mis: 3
  idx_m23 <- which( m2 &  m3)   # obs: 1   ; mis: 2,3
  
  Y_cc  <- if (length(idx_cc))  Y[idx_cc,  c(1,2,3), drop = FALSE] else NULL
  Y_13  <- if (length(idx_m2))  Y[idx_m2,  c(1,3),   drop = FALSE] else NULL
  Y_12  <- if (length(idx_m3))  Y[idx_m3,  c(1,2),   drop = FALSE] else NULL
  Y_1   <- if (length(idx_m23)) matrix(Y[idx_m23, 1], ncol = 1) else NULL
  
  ## ---------- INIT (complete-case kmeans) ----------
  ## Fallback: if too few complete cases, use mean-impute+kmeans
  use_cc_init <- (length(idx_cc) >= max(G, 5))  # heuristic safety
  
  if (use_cc_init) {
    km <- kmeans(Y_cc, centers = G, nstart = nstart_kmeans)
    cl <- km$cluster
    
    ## alpha from complete-case proportions
    alpha <- as.numeric(tabulate(cl, nbins = G))
    alpha <- alpha / sum(alpha)
    alpha[alpha <= 0] <- 1 / G
    alpha <- alpha / sum(alpha)
    
    ## mu from kmeans centers
    mu <- km$centers
    mu <- matrix(mu, nrow = G, ncol = p)
    
    ## Sigma from complete cases
    Sigma <- cov(Y_cc)
    Sigma <- (Sigma + t(Sigma)) / 2
    Sigma <- Sigma + diag(ridge, p)
    Sigma <- make_pd(Sigma)
  } else {
    ## fallback: mean-impute + kmeans (only when CC insufficient)
    Yimp <- Y
    for (j in 1:p) Yimp[is.na(Yimp[,j]), j] <- mean(Yimp[,j], na.rm = TRUE)
    
    km <- kmeans(Yimp, centers = G, nstart = nstart_kmeans)
    cl <- km$cluster
    
    alpha <- as.numeric(tabulate(cl, nbins = G)) / n
    mu <- matrix(0, nrow = G, ncol = p)
    for (g in 1:G) {
      idx <- which(cl == g)
      if (length(idx) == 0) {
        mu[g,] <- colMeans(Yimp)
        alpha[g] <- 1 / G
      } else {
        mu[g,] <- colMeans(Yimp[idx, , drop = FALSE])
      }
    }
    
    Sigma <- cov(Yimp)
    Sigma <- (Sigma + t(Sigma)) / 2
    Sigma <- Sigma + diag(ridge, p)
    Sigma <- make_pd(Sigma)
  }
  ## ---------- END INIT ----------
  
  loglik_prev <- NA_real_
  post <- matrix(1/G, nrow = n, ncol = G)
  
  for (it in 1:maxit) {
    
    ## E-step
    logw <- matrix(-Inf, nrow = n, ncol = G)
    
    if (length(idx_cc)) {
      S_oo <- make_pd(Sigma[c(1,2,3), c(1,2,3), drop = FALSE])
      for (g in 1:G) {
        ll <- dmvnorm(Y_cc, mean = mu[g, c(1,2,3)], sigma = S_oo, log = TRUE)
        logw[idx_cc, g] <- log(alpha[g]) + ll
      }
    }
    
    if (length(idx_m2)) {
      S_oo <- make_pd(Sigma[c(1,3), c(1,3), drop = FALSE])
      for (g in 1:G) {
        ll <- dmvnorm(Y_13, mean = mu[g, c(1,3)], sigma = S_oo, log = TRUE)
        logw[idx_m2, g] <- log(alpha[g]) + ll
      }
    }
    
    if (length(idx_m3)) {
      S_oo <- make_pd(Sigma[c(1,2), c(1,2), drop = FALSE])
      for (g in 1:G) {
        ll <- dmvnorm(Y_12, mean = mu[g, c(1,2)], sigma = S_oo, log = TRUE)
        logw[idx_m3, g] <- log(alpha[g]) + ll
      }
    }
    
    if (length(idx_m23)) {
      S_oo <- make_pd(matrix(Sigma[1,1], 1, 1))
      for (g in 1:G) {
        ll <- dmvnorm(Y_1, mean = mu[g, 1], sigma = S_oo, log = TRUE)
        logw[idx_m23, g] <- log(alpha[g]) + ll
      }
    }
    
    mmax <- apply(logw, 1, max)
    w <- exp(logw - mmax)
    rs <- rowSums(w)
    
    bad <- (!is.finite(rs)) | (rs <= 0)
    if (any(bad)) {
      post[bad, ] <- 1/G
      rs[bad] <- 1
    }
    post[!bad, ] <- w[!bad, , drop = FALSE] / rs[!bad]
    
    loglik <- sum(mmax + log(rs))
    
    ## M-step
    alpha_new <- colMeans(post)
    
    S1 <- matrix(0, nrow = G, ncol = p)
    S2 <- array(0, dim = c(p, p, G))
    
    if (length(idx_m2)) {
      obs <- c(1,3); mis <- 2
      S_oo <- make_pd(Sigma[obs, obs, drop = FALSE])
      S_mo <- Sigma[mis, obs, drop = FALSE]
      S_om <- Sigma[obs, mis, drop = FALSE]
      S_mm <- matrix(Sigma[mis, mis], 1, 1)
      S_oo_inv <- safe_solve(S_oo)
      A_m2 <- S_mo %*% S_oo_inv
      cond_var_m2 <- as.numeric(S_mm - S_mo %*% S_oo_inv %*% S_om)
      if (!is.finite(cond_var_m2) || cond_var_m2 < 0) cond_var_m2 <- 0
    }
    
    if (length(idx_m3)) {
      obs <- c(1,2); mis <- 3
      S_oo <- make_pd(Sigma[obs, obs, drop = FALSE])
      S_mo <- Sigma[mis, obs, drop = FALSE]
      S_om <- Sigma[obs, mis, drop = FALSE]
      S_mm <- matrix(Sigma[mis, mis], 1, 1)
      S_oo_inv <- safe_solve(S_oo)
      A_m3 <- S_mo %*% S_oo_inv
      cond_var_m3 <- as.numeric(S_mm - S_mo %*% S_oo_inv %*% S_om)
      if (!is.finite(cond_var_m3) || cond_var_m3 < 0) cond_var_m3 <- 0
    }
    
    if (length(idx_m23)) {
      obs <- 1; mis <- c(2,3)
      S_oo <- make_pd(matrix(Sigma[obs, obs], 1, 1))
      S_mo <- Sigma[mis, obs, drop = FALSE]
      S_om <- Sigma[obs, mis, drop = FALSE]
      S_mm <- Sigma[mis, mis, drop = FALSE]
      S_oo_inv <- safe_solve(S_oo)
      A_m23 <- S_mo %*% S_oo_inv
      cond_cov_m23 <- make_pd(S_mm - S_mo %*% S_oo_inv %*% S_om)
    }
    
    for (g in 1:G) {
      wg_all <- post[, g]
      
      if (length(idx_cc)) {
        w0 <- wg_all[idx_cc]
        sw <- sqrt(w0)
        Ey <- Y_cc
        S1[g, ] <- S1[g, ] + as.vector(crossprod(Ey, w0))
        S2[, , g] <- S2[, , g] + crossprod(Ey * sw)
      }
      
      if (length(idx_m2)) {
        w0 <- wg_all[idx_m2]
        sw <- sqrt(w0)
        delta <- Y_13 - matrix(mu[g, c(1,3)], nrow = nrow(Y_13), ncol = 2, byrow = TRUE)
        Ey2 <- as.numeric(mu[g, 2] + delta %*% t(A_m2))
        Ey <- cbind(Y_13[,1], Ey2, Y_13[,2])
        S1[g, ] <- S1[g, ] + as.vector(crossprod(Ey, w0))
        S2[, , g] <- S2[, , g] + crossprod(Ey * sw)
        S2[2,2,g] <- S2[2,2,g] + sum(w0) * cond_var_m2
      }
      
      if (length(idx_m3)) {
        w0 <- wg_all[idx_m3]
        sw <- sqrt(w0)
        delta <- Y_12 - matrix(mu[g, c(1,2)], nrow = nrow(Y_12), ncol = 2, byrow = TRUE)
        Ey3 <- as.numeric(mu[g, 3] + delta %*% t(A_m3))
        Ey <- cbind(Y_12[,1], Y_12[,2], Ey3)
        S1[g, ] <- S1[g, ] + as.vector(crossprod(Ey, w0))
        S2[, , g] <- S2[, , g] + crossprod(Ey * sw)
        S2[3,3,g] <- S2[3,3,g] + sum(w0) * cond_var_m3
      }
      
      if (length(idx_m23)) {
        w0 <- wg_all[idx_m23]
        sw <- sqrt(w0)
        delta <- Y_1 - matrix(mu[g, 1], nrow = nrow(Y_1), ncol = 1, byrow = TRUE)
        Ey23 <- delta %*% t(A_m23) + matrix(mu[g, c(2,3)], nrow = nrow(Y_1), ncol = 2, byrow = TRUE)
        Ey <- cbind(Y_1[,1], Ey23[,1], Ey23[,2])
        S1[g, ] <- S1[g, ] + as.vector(crossprod(Ey, w0))
        S2[, , g] <- S2[, , g] + crossprod(Ey * sw)
        S2[2:3, 2:3, g] <- S2[2:3, 2:3, g] + sum(w0) * cond_cov_m23
      }
    }
    
    mu_new <- mu
    for (g in 1:G) {
      denom <- sum(post[, g])
      if (denom <= 0) {
        ## robust fallback
        if (length(idx_cc) >= 2) mu_new[g, ] <- colMeans(Y_cc) else mu_new[g, ] <- rep(0, p)
      } else {
        mu_new[g, ] <- S1[g, ] / denom
      }
    }
    
    Sigma_num <- matrix(0, p, p)
    for (g in 1:G) {
      denom <- sum(post[, g])
      mvec <- mu_new[g, ]
      S1g <- S1[g, ]
      S2g <- S2[, , g]
      Sigma_num <- Sigma_num +
        (S2g - tcrossprod(mvec, S1g) - tcrossprod(S1g, mvec) + denom * tcrossprod(mvec))
    }
    
    Sigma_new <- Sigma_num / n
    Sigma_new <- (Sigma_new + t(Sigma_new)) / 2
    Sigma_new <- Sigma_new + diag(ridge, p)
    Sigma_new <- make_pd(Sigma_new)
    
    if (verbose) message(sprintf("G=%d iter=%d loglik=%.6f", G, it, loglik))
    
    converged <- FALSE
    if (is.finite(loglik_prev) && is.finite(loglik)) {
      relchg <- abs(loglik - loglik_prev) / (1 + abs(loglik_prev))
      converged <- is.finite(relchg) && (relchg < tol)
    }
    if (converged) {
      return(list(alpha = alpha_new, mu = mu_new, Sigma = Sigma_new,
                  loglik = loglik, post = post, iter = it, converged = TRUE))
    }
    
    alpha <- alpha_new
    mu <- mu_new
    Sigma <- Sigma_new
    loglik_prev <- loglik
  }
  
  list(alpha = alpha, mu = mu, Sigma = Sigma,
       loglik = loglik_prev, post = post, iter = maxit, converged = FALSE)
}

select_g_by_bic <- function(Y, Gmax = 10, ...) {
  n <- nrow(Y); p <- ncol(Y)
  fits <- vector("list", Gmax)
  bics <- rep(Inf, Gmax)
  
  for (G in 1:Gmax) {
    fit <- gmm_em_missing(Y, G, ...)
    fits[[G]] <- fit
    k <- (G - 1) + G * p
    bics[G] <- -2 * fit$loglik + log(n) * k
  }
  
  Ghat <- which.min(bics)
  list(Ghat = Ghat, fit = fits[[Ghat]], bics = bics, fits = fits)
}

## -------------------------
## SPFI point estimator
## -------------------------
spfi_point <- function(Ymis, c2, c3, M = 2000, Gmax = 10,
                       em_args = list(maxit = 200, tol = 1e-6, ridge = 1e-6),
                       seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(Ymis); p <- ncol(Ymis)
  stopifnot(p == 3)
  
  sel <- do.call(select_g_by_bic, c(list(Y = Ymis, Gmax = Gmax), em_args))
  fit <- sel$fit
  G <- sel$Ghat
  alpha <- fit$alpha
  mu <- fit$mu
  Sigma <- fit$Sigma
  
  post_one <- function(y) {
    logw <- rep(NA_real_, G)
    for (g in 1:G) logw[g] <- log(alpha[g]) + log_mvn_marginal(y, mu[g,], Sigma)
    mmax <- max(logw)
    ww <- exp(logw - mmax)
    s <- sum(ww)
    if (!is.finite(s) || s <= 0) return(rep(1/G, G))
    out <- ww / s
    if (any(!is.finite(out))) rep(1/G, G) else out
  }
  
  t2_sum <- 0; t3_sum <- 0; p2_sum <- 0; p3_sum <- 0
  
  for (i in 1:n) {
    y <- Ymis[i,]
    if (all(!is.na(y))) {
      t2_sum <- t2_sum + unname(y[2])
      t3_sum <- t3_sum + unname(y[3])
      p2_sum <- p2_sum + as.numeric(unname(y[2]) < c2)
      p3_sum <- p3_sum + as.numeric(unname(y[3]) < c3)
      next
    }
    
    pi_g <- post_one(y)
    counts <- as.vector(rmultinom(1, size = M, prob = pi_g))
    obs <- which(!is.na(y))
    mis <- which(is.na(y))
    
    t2_i <- 0; t3_i <- 0; p2_i <- 0; p3_i <- 0
    
    for (g in 1:G) {
      mg <- counts[g]
      if (mg == 0) next
      w_ig <- pi_g[g] / mg
      
      mu_g <- mu[g,]
      yobs <- y[obs]
      muobs <- mu_g[obs]
      mumis <- mu_g[mis]
      
      S_oo <- Sigma[obs, obs, drop = FALSE]
      S_mo <- Sigma[mis, obs, drop = FALSE]
      S_om <- Sigma[obs, mis, drop = FALSE]
      S_mm <- Sigma[mis, mis, drop = FALSE]
      
      S_oo_inv <- safe_solve(S_oo)
      cond_mean_m <- as.vector(mumis + S_mo %*% S_oo_inv %*% (yobs - muobs))
      cond_cov_m  <- S_mm - S_mo %*% S_oo_inv %*% S_om
      cond_cov_m  <- make_pd(cond_cov_m)
      
      ymis_draw <- rmvnorm(mg, mean = cond_mean_m, sigma = cond_cov_m)
      
      for (j in 1:mg) {
        ystar <- y
        ystar[mis] <- ymis_draw[j,]
        
        v2 <- unname(ystar[2])
        v3 <- unname(ystar[3])
        
        t2_i <- t2_i + w_ig * v2
        t3_i <- t3_i + w_ig * v3
        p2_i <- p2_i + w_ig * as.numeric(v2 < c2)
        p3_i <- p3_i + w_ig * as.numeric(v3 < c3)
      }
    }
    
    t2_sum <- t2_sum + t2_i
    t3_sum <- t3_sum + t3_i
    p2_sum <- p2_sum + p2_i
    p3_sum <- p3_sum + p3_i
  }
  
  est <- c(theta2 = t2_sum / n,
           theta3 = t3_sum / n,
           P2 = p2_sum / n,
           P3 = p3_sum / n)
  
  list(est = est, Ghat = G, bic = sel$bics)
}

## -------------------------
## One replication (Full vs SPFI) + Ghat for FULL too
## -------------------------
one_rep <- function(model, n = 500, M = 2000, Gmax = 10,
                    em_args = list(maxit = 200, tol = 1e-6, ridge = 1e-6),
                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cuts <- get_cuts(model); c2 <- cuts["c2"]; c3 <- cuts["c3"]
  
  Yfull <- gen_data(model, n)
  
  full_est <- c(theta2 = mean(Yfull[,2]),
                theta3 = mean(Yfull[,3]),
                P2 = mean(Yfull[,2] < c2),
                P3 = mean(Yfull[,3] < c3))
  
  sel_full <- do.call(select_g_by_bic, c(list(Y = Yfull, Gmax = Gmax), em_args))
  Ghat_full <- sel_full$Ghat
  
  Ymis <- impose_missing(Yfull)
  spfi <- spfi_point(Ymis, c2, c3, M = M, Gmax = Gmax, em_args = em_args)
  
  list(full_est = full_est,
       spfi_est = spfi$est,
       Ghat_full = Ghat_full,
       Ghat_spfi = spfi$Ghat)
}

## -------------------------
## Simulation driver
## -------------------------
run_simulation <- function(models = c("M1","M2","M3","M4"),
                           B = 2000, n = 500, M = 2000, Gmax = 10,
                           ncores = max(1, detectCores() - 1),
                           seed = 123,
                           em_args = list(maxit = 200, tol = 1e-6, ridge = 1e-6),
                           truth_rel_tol = 1e-10,
                           truth_abs_tol = 0,
                           progress = TRUE,
                           chunk_size = 25) {
  set.seed(seed)
  
  if (progress) {
    message(sprintf("[run_simulation] start  (models=%d, B=%d, n=%d, M=%d, Gmax=%d, ncores=%d)  time=%s",
                    length(models), B, n, M, Gmax, ncores, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    flush.console()
  }
  
  if (progress) {
    message("[run_simulation] computing true parameters ...")
    flush.console()
  }
  true_list <- setNames(vector("list", length(models)), models)
  for (m in models) {
    true_list[[m]] <- true_params_exact(m, rel.tol = truth_rel_tol, abs.tol = truth_abs_tol)
  }
  
  use_cluster <- (ncores >= 2)
  if (use_cluster) {
    cl <- makeCluster(ncores)
    on.exit(stopCluster(cl), add = TRUE)
    clusterEvalQ(cl, { library(mvtnorm) })
    
    to_export <- c(
      "one_rep","gen_data","impose_missing","get_cuts",
      "spfi_point","select_g_by_bic","gmm_em_missing",
      "log_mvn_marginal","expit","Sigma_AR","make_pd","safe_solve",
      "pexgauss","true_params_exact"
    )
    clusterExport(cl, varlist = to_export, envir = .GlobalEnv)
  }
  
  out <- list()
  for (mi in seq_along(models)) {
    m <- models[mi]
    seeds <- sample.int(1e9, B)
    
    if (progress) {
      message(sprintf("[Model %s] start (%d/%d)  B=%d  time=%s",
                      m, mi, length(models), B, format(Sys.time(), "%H:%M:%S")))
      flush.console()
    }
    
    if (use_cluster) {
      reps <- vector("list", B)
      idx_all <- 1:B
      chunks <- split(idx_all, ceiling(seq_along(idx_all) / chunk_size))
      
      done <- 0L
      for (ck in seq_along(chunks)) {
        ids <- chunks[[ck]]
        tmp <- parLapply(cl, ids, function(b) {
          one_rep(m, n = n, M = M, Gmax = Gmax, em_args = em_args, seed = seeds[b])
        })
        reps[ids] <- tmp
        
        done <- done + length(ids)
        if (progress) {
          message(sprintf("[Model %s] progress: %d/%d (%.1f%%)  time=%s",
                          m, done, B, 100 * done / B, format(Sys.time(), "%H:%M:%S")))
          flush.console()
        }
      }
    } else {
      reps <- vector("list", B)
      if (progress) {
        pb <- txtProgressBar(min = 0, max = B, style = 3)
        on.exit(try(close(pb), silent = TRUE), add = TRUE)
      }
      for (b in 1:B) {
        reps[[b]] <- one_rep(m, n = n, M = M, Gmax = Gmax, em_args = em_args, seed = seeds[b])
        if (progress) setTxtProgressBar(pb, b)
      }
      if (progress) {
        close(pb)
        flush.console()
      }
    }
    
    full_mat <- do.call(rbind, lapply(reps, `[[`, "full_est"))
    spfi_mat <- do.call(rbind, lapply(reps, `[[`, "spfi_est"))
    
    Ghat_full_vec <- sapply(reps, `[[`, "Ghat_full")
    Ghat_spfi_vec <- sapply(reps, `[[`, "Ghat_spfi")
    
    colnames(full_mat) <- c("theta2","theta3","P2","P3")
    colnames(spfi_mat) <- c("theta2","theta3","P2","P3")
    
    out[[m]] <- list(true = true_list[[m]],
                     full = full_mat,
                     spfi = spfi_mat,
                     Ghat_full = Ghat_full_vec,
                     Ghat_spfi = Ghat_spfi_vec)
    
    if (progress) {
      message(sprintf("[Model %s] done  time=%s", m, format(Sys.time(), "%H:%M:%S")))
      flush.console()
    }
  }
  
  if (progress) {
    message(sprintf("[run_simulation] done  time=%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    flush.console()
  }
  
  out
}


## -------------------------
## Figure 1/2: boxplots (Full vs SPFI)
## -------------------------
make_figure12_full_vs_spfi <- function(sim_out, models, outfile, title) {
  models <- intersect(models, names(sim_out))
  if (!length(models)) stop("No requested models in sim_out.")
  
  df <- bind_rows(lapply(models, function(m) {
    full <- as.data.frame(sim_out[[m]]$full); full$Method <- "Full"
    spfi <- as.data.frame(sim_out[[m]]$spfi); spfi$Method <- "SPFI"
    full$Model <- m; spfi$Model <- m
    bind_rows(full, spfi)
  }))
  
  df$Method <- factor(df$Method, levels = c("Full","SPFI"))
  
  df_long <- df %>%
    pivot_longer(cols = c(theta2, theta3, P2, P3),
                 names_to = "Parameter", values_to = "Estimate")
  
  df_long$Parameter <- factor(
    df_long$Parameter,
    levels = c("theta2","theta3","P2","P3"),
    labels = c("Theta 2","Theta 3","P2","P3")
  )
  
  df_long$Model <- factor(df_long$Model, levels = models)
  df_long$Panel <- interaction(df_long$Model, df_long$Parameter, sep = " : ", lex.order = TRUE)
  
  true_df <- bind_rows(lapply(models, function(m) {
    tr <- sim_out[[m]]$true
    data.frame(
      Model = factor(m, levels = models),
      Parameter = factor(names(tr),
                         levels = c("theta2","theta3","P2","P3"),
                         labels = c("Theta 2","Theta 3","P2","P3")),
      True = as.numeric(tr)
    )
  }))
  true_df$Panel <- interaction(true_df$Model, true_df$Parameter, sep = " : ", lex.order = TRUE)
  
  p <- ggplot(df_long, aes(x = Method, y = Estimate)) +
    geom_boxplot(outlier.size = 0.6) +
    geom_hline(data = true_df, aes(yintercept = True), linetype = "dashed") +
    facet_wrap(~ Panel, ncol = 4, scales = "free_y") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.background = element_rect(fill = "grey90", color = "grey50")
    ) +
    labs(x = NULL, y = NULL, title = title, subtitle = "Dashed lines = true values")
  
  ggsave(outfile, p, width = 11, height = 6, dpi = 250)
  invisible(p)
}

## -------------------------
## Figure 3: histogram of selected G (paper-like)
## -------------------------
make_figure3_hist_G <- function(sim_out,
                                models = c("M1","M2","M3","M4"),
                                outfile = "Figure3_hist_G.png",
                                Gmax = 10) {
  models <- intersect(models, names(sim_out))
  if (!length(models)) stop("No requested models in sim_out.")
  
  dfG <- dplyr::bind_rows(lapply(models, function(m) {
    data.frame(Model = m, Ghat = sim_out[[m]]$Ghat_spfi)
  }))
  
  dfG$Model <- factor(dfG$Model, levels = models)
  
  ## 離散軸にして「棒の幅/隙間」を制御しやすくする（論文の雰囲気）
  dfG$Ghat <- factor(dfG$Ghat, levels = 1:Gmax)
  
  ## ggh4x が入っていれば、論文っぽく “margins だけ軸を表示”
  facet_fun <- if (requireNamespace("ggh4x", quietly = TRUE)) {
    function(...) ggh4x::facet_wrap2(..., axes = "margins")
  } else {
    ggplot2::facet_wrap
  }
  
  p <- ggplot2::ggplot(dfG, ggplot2::aes(x = Ghat)) +
    ggplot2::geom_bar(
      width = 0.55,          # ← 棒を細くして隙間を広げる（重要）
      colour = "black",
      fill   = "black"
    ) +
    facet_fun(~ Model, nrow = 2, ncol = 2) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::labs(x = "G", y = "Count") +     # ← 論文はタイトル無し＆x=G
    ggplot2::theme_gray(base_size = 12) +
    ggplot2::theme(
      ## 外側の余白（ここで「余白が少ない」を改善）
      plot.margin  = grid::unit(c(8, 12, 8, 10), "mm"),
      
      ## パネル間の余白（論文はわりと間がある）
      panel.spacing = grid::unit(1.0, "lines"),
      
      ## facet strip / panel border を論文寄りに
      strip.background = ggplot2::element_rect(fill = "grey80", colour = "grey50"),
      strip.text       = ggplot2::element_text(size = 11),
      panel.border     = ggplot2::element_rect(colour = "grey50", fill = NA, linewidth = 0.6),
      
      ## タイトルは消す（論文はキャプション扱い）
      plot.title = ggplot2::element_blank(),
      
      ## 軸タイトルの間隔
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 6))
    )
  
  ggplot2::ggsave(outfile, p, width = 7.2, height = 4.8, dpi = 300)
  invisible(p)
}



## =========================================================
## RUN
## =========================================================
## (Paper: B=2000, M=2000. まずは B を小さくして動作確認推奨)
sim_all <- run_simulation(models = c("M1","M2","M3","M4"),
                          B = 2000, n = 500, M = 2000, Gmax = 10,
                          ncores = 16, seed = 123)

make_figure12_full_vs_spfi(sim_all, models = c("M1","M2"),
                           outfile = "Figure1_Full_vs_SPFI_M1_M2.png",
                           title = "Figure 1 style: M1 & M2 (Full vs SPFI)")

make_figure12_full_vs_spfi(sim_all, models = c("M3","M4"),
                           outfile = "Figure2_Full_vs_SPFI_M3_M4.png",
                           title = "Figure 2 style: M3 & M4 (Full vs SPFI)")

make_figure3_hist_G(sim_all, outfile = "Figure3_hist_G.png")

saveRDS(sim_all, file = "sim_spfi.rds")








## =========================================================
##  Simulation 1 (paper) : FULL vs PFI  [ALL-IN-ONE SCRIPT]
##  - PFI here = Parametric Fractional Imputation under SINGLE MVN for (Y1,Y2,Y3)
##  - Missingness: select exactly floor(0.25*n) for Y2 and for Y3 independently
##                 with PPS probs proportional to pi_ij = expit(logit using Y1)
##  - Check items: theta2 = E(Y2), theta3 = E(Y3),  P2 = Pr(Y2 < c2), P3 = Pr(Y3 < c3)
##  - Output:
##      Figure1-like: M1 & M2 (theta2, theta3, P2, P3)
##      Figure2-like: M3 & M4 (theta2, theta3, P2, P3)
## =========================================================

suppressPackageStartupMessages({
  library(mvtnorm)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(parallel)
})

## -------------------------
## helpers (numerical)
## -------------------------
expit <- function(x) 1 / (1 + exp(-x))

make_pd <- function(S, eps = 1e-10, max_tries = 12) {
  S <- (S + t(S)) / 2
  if (nrow(S) == 1L) {
    if (!is.finite(S[1,1])) S[1,1] <- eps
    if (S[1,1] <= 0) S[1,1] <- eps
    return(S)
  }
  jitter <- eps
  for (k in 0:max_tries) {
    out <- tryCatch(chol(S + diag(jitter, nrow(S))), error = function(e) NULL)
    if (!is.null(out)) return(S + diag(jitter, nrow(S)))
    jitter <- jitter * 10
  }
  ev <- eigen(S, symmetric = TRUE)
  vals <- pmax(ev$values, eps)
  S2 <- ev$vectors %*% diag(vals) %*% t(ev$vectors)
  (S2 + t(S2)) / 2
}

safe_solve <- function(A) {
  A <- make_pd(A)
  solve(A)
}

log_mvn_marginal <- function(y, mu, Sigma) {
  obs <- which(!is.na(y))
  yobs <- y[obs]
  muobs <- mu[obs]
  S_oo <- Sigma[obs, obs, drop = FALSE]
  S_oo <- make_pd(S_oo)
  dmvnorm(yobs, mean = muobs, sigma = S_oo, log = TRUE)
}

## -------------------------
## Data generation (Simulation 1)
## -------------------------
Sigma_AR <- function(rho) {
  matrix(c(1, rho, rho^2,
           rho, 1, rho,
           rho^2, rho, 1), nrow = 3, byrow = TRUE)
}

gen_data <- function(model = c("M1","M2","M3","M4"), n = 500) {
  model <- match.arg(model)
  
  if (model %in% c("M1","M2")) {
    alpha <- c(0.3, 0.3, 0.4)
    mu_list <- list(c(-3,-3,-3), c(1,1,1), c(5,5,5))
    Sigma <- Sigma_AR(0.7)
    
    z <- sample(1:3, size = n, replace = TRUE, prob = alpha)
    Y <- matrix(NA_real_, n, 3)
    
    for (g in 1:3) {
      idx <- which(z == g)
      if (length(idx) == 0) next
      if (model == "M2" && g == 2) {
        Y[idx, ] <- matrix(rexp(length(idx) * 3, rate = 1), ncol = 3)
      } else {
        Y[idx, ] <- rmvnorm(length(idx), mean = mu_list[[g]], sigma = Sigma)
      }
    }
    colnames(Y) <- c("Y1","Y2","Y3")
    return(Y)
  }
  
  if (model == "M3") {
    e1 <- rnorm(n, mean = 0, sd = 1)
    e2 <- rgamma(n, shape = 1, rate = 1)  # Exp(1)
    e3 <- rchisq(n, df = 1)
    Y1 <- 1 + e1
    Y2 <- 0.5 * Y1 + e2
    Y3 <- Y2 + e3
    Y <- cbind(Y1,Y2,Y3)
    colnames(Y) <- c("Y1","Y2","Y3")
    return(Y)
  }
  
  if (model == "M4") {
    mu12 <- c(1,2)
    Sig12 <- matrix(c(1,0.5,0.5,1), 2, 2)
    Y12 <- rmvnorm(n, mean = mu12, sigma = Sig12)
    Y1 <- Y12[,1]
    Y2 <- Y12[,2]
    Y3 <- Y2^2 + rnorm(n, mean = 0, sd = 1)
    Y <- cbind(Y1,Y2,Y3)
    colnames(Y) <- c("Y1","Y2","Y3")
    return(Y)
  }
  
  stop("Unknown model")
}

## =========================================================
## Missingness (REVISED to match paper description)
##  - For j=2,3 independently:
##      compute pi_ij
##      select exactly floor(miss_frac*n) indices without replacement
##      with probs proportional to pi_ij
## =========================================================
impose_missing <- function(Y, miss_frac = 0.25) {
  n <- nrow(Y)
  k <- as.integer(floor(miss_frac * n))
  if (k <= 0 || k >= n) stop("miss_frac must satisfy 0 < floor(miss_frac*n) < n")
  
  ## selection probabilities based on Y1
  pi2 <- expit(-0.8 + 0.4 * Y[, "Y1"])
  pi3 <- expit( 0.4 - 0.8 * Y[, "Y1"])
  
  ## numerical safeguard: ensure strictly positive weights
  w2 <- pmax(pi2, 1e-12)
  w3 <- pmax(pi3, 1e-12)
  
  ## independently select 25% for Y2 and 25% for Y3
  idx2 <- sample.int(n, size = k, replace = FALSE, prob = w2)
  idx3 <- sample.int(n, size = k, replace = FALSE, prob = w3)
  
  Ymis <- Y
  Ymis[idx2, "Y2"] <- NA_real_
  Ymis[idx3, "Y3"] <- NA_real_
  Ymis
}

## (optional) quick check
missing_summary <- function(Ymis) {
  n <- nrow(Ymis)
  pY2 <- mean(is.na(Ymis[, "Y2"]))
  pY3 <- mean(is.na(Ymis[, "Y3"]))
  poverall <- mean(is.na(Ymis[, "Y2"]) | is.na(Ymis[, "Y3"]))
  c(p_missing_Y2 = pY2, p_missing_Y3 = pY3, p_missing_any = poverall, n = n)
}

get_cuts <- function(model) {
  if (model %in% c("M1","M2")) return(c(c2 = -2, c3 = -2))
  if (model == "M3") return(c(c2 = 2, c3 = 3))
  if (model == "M4") return(c(c2 = 2, c3 = 5))
  stop("unknown model")
}

## -------------------------
## TRUE PARAMETERS (Exact / 1D quadrature)
## -------------------------
pexgauss <- function(x, mu, sigma, tau) {
  if (tau <= 0) stop("tau must be > 0")
  z  <- (x - mu) / sigma
  z2 <- z - sigma / tau
  
  logPhi_z2 <- pnorm(z2, log.p = TRUE)
  a <- -(x - mu) / tau + (sigma^2) / (2 * tau^2)
  
  s <- a + logPhi_z2
  s[is.nan(s)] <- -Inf
  
  term2 <- exp(s)
  F <- pnorm(z) - term2
  F <- pmin(pmax(F, 0), 1)
  
  bad <- !is.finite(F)
  if (any(bad)) {
    F[bad] <- ifelse(x[bad] < mu, 0, 1)
  }
  F
}

true_params_exact <- function(model, rel.tol = 1e-10, abs.tol = 0) {
  model <- match.arg(model, c("M1","M2","M3","M4"))
  cuts <- get_cuts(model)
  c2 <- as.numeric(cuts["c2"]); c3 <- as.numeric(cuts["c3"])
  
  if (model == "M1") {
    alpha <- c(0.3, 0.3, 0.4)
    mu <- c(-3, 1, 5)
    theta <- sum(alpha * mu)
    P2 <- sum(alpha * pnorm(c2, mean = mu, sd = 1))
    P3 <- sum(alpha * pnorm(c3, mean = mu, sd = 1))
    return(c(theta2 = theta, theta3 = theta, P2 = P2, P3 = P3))
  }
  
  if (model == "M2") {
    alpha <- c(0.3, 0.3, 0.4)
    theta <- sum(alpha * c(-3, 1, 5))  # = 1.4
    
    P2 <- alpha[1] * pnorm(c2, mean = -3, sd = 1) +
      alpha[2] * pexp(c2, rate = 1) +
      alpha[3] * pnorm(c2, mean =  5, sd = 1)
    P3 <- alpha[1] * pnorm(c3, mean = -3, sd = 1) +
      alpha[2] * pexp(c3, rate = 1) +
      alpha[3] * pnorm(c3, mean =  5, sd = 1)
    
    return(c(theta2 = theta, theta3 = theta, P2 = P2, P3 = P3))
  }
  
  if (model == "M3") {
    mu <- 0.5
    sigma <- 0.5
    tau <- 1
    
    theta2 <- mu + tau
    theta3 <- theta2 + 1
    
    P2 <- pexgauss(c2, mu, sigma, tau)
    
    integrand <- function(t) {
      val <- pexgauss(c3 - t, mu, sigma, tau) * dchisq(t, df = 1)
      val[!is.finite(val)] <- 0
      val
    }
    P3 <- integrate(integrand, lower = 0, upper = Inf,
                    rel.tol = rel.tol, abs.tol = abs.tol)$value
    
    return(c(theta2 = theta2, theta3 = theta3, P2 = P2, P3 = P3))
  }
  
  if (model == "M4") {
    theta2 <- 2
    theta3 <- 1 + 2^2
    
    P2 <- pnorm(c2, mean = 2, sd = 1)
    
    integrand <- function(y) {
      val <- pnorm(c3 - y^2) * dnorm(y, mean = 2, sd = 1)
      val[!is.finite(val)] <- 0
      val
    }
    P3 <- integrate(integrand, lower = -Inf, upper = Inf,
                    rel.tol = rel.tol, abs.tol = abs.tol)$value
    
    return(c(theta2 = theta2, theta3 = theta3, P2 = P2, P3 = P3))
  }
  
  stop("unknown model")
}

## =========================================================
## PFI (Parametric Fractional Imputation) under SINGLE MVN
##  - EM for MVN with missing (observed likelihood)
##  - FI draws from conditional MVN, weights = 1/M
##  - Return (theta2, theta3, P2, P3)
## =========================================================
mvn_em_missing <- function(Y,
                           maxit = 200, tol = 1e-6,
                           ridge = 1e-6, verbose = FALSE) {
  n <- nrow(Y); p <- ncol(Y)
  stopifnot(p == 3)
  
  ## init: mean-impute
  Yimp <- Y
  for (j in 1:p) Yimp[is.na(Yimp[,j]), j] <- mean(Yimp[,j], na.rm = TRUE)
  
  mu <- colMeans(Yimp)
  Sigma <- cov(Yimp)
  Sigma <- (Sigma + t(Sigma)) / 2
  Sigma <- Sigma + diag(ridge, p)
  Sigma <- make_pd(Sigma)
  
  loglik_prev <- NA_real_
  
  for (it in 1:maxit) {
    Ey_sum  <- rep(0, p)
    Eyy_sum <- matrix(0, p, p)
    
    for (i in 1:n) {
      y <- Y[i,]
      obs <- which(!is.na(y))
      mis <- which(is.na(y))
      
      if (length(mis) == 0) {
        Ey_sum  <- Ey_sum  + y
        Eyy_sum <- Eyy_sum + tcrossprod(y)
        next
      }
      if (length(obs) == 0) stop("All entries missing: not expected here.")
      
      yobs  <- y[obs]
      muobs <- mu[obs]
      mumis <- mu[mis]
      
      S_oo <- Sigma[obs, obs, drop = FALSE]
      S_mo <- Sigma[mis, obs, drop = FALSE]
      S_om <- Sigma[obs, mis, drop = FALSE]
      S_mm <- Sigma[mis, mis, drop = FALSE]
      
      S_oo_inv <- safe_solve(S_oo)
      cond_mean_m <- as.vector(mumis + S_mo %*% S_oo_inv %*% (yobs - muobs))
      cond_cov_m  <- S_mm - S_mo %*% S_oo_inv %*% S_om
      cond_cov_m  <- make_pd(cond_cov_m)
      
      Ey <- rep(NA_real_, p)
      Ey[obs] <- yobs
      Ey[mis] <- cond_mean_m
      
      Eyy <- tcrossprod(Ey)
      Eyy[mis, mis] <- Eyy[mis, mis, drop = FALSE] + cond_cov_m
      
      Ey_sum  <- Ey_sum  + Ey
      Eyy_sum <- Eyy_sum + Eyy
    }
    
    mu_new <- Ey_sum / n
    Sigma_new <- Eyy_sum / n - tcrossprod(mu_new)
    Sigma_new <- (Sigma_new + t(Sigma_new)) / 2
    Sigma_new <- Sigma_new + diag(ridge, p)
    Sigma_new <- make_pd(Sigma_new)
    
    ## observed loglik for convergence
    loglik <- 0
    for (i in 1:n) loglik <- loglik + log_mvn_marginal(Y[i,], mu_new, Sigma_new)
    
    if (verbose) message(sprintf("MVN-EM iter=%d loglik=%.6f", it, loglik))
    
    converged <- FALSE
    if (is.finite(loglik_prev) && is.finite(loglik)) {
      relchg <- abs(loglik - loglik_prev) / (1 + abs(loglik_prev))
      converged <- is.finite(relchg) && (relchg < tol)
    }
    
    mu <- mu_new
    Sigma <- Sigma_new
    loglik_prev <- loglik
    if (converged) break
  }
  
  list(mu = mu, Sigma = Sigma, loglik = loglik_prev)
}

pfi_point_theta2_theta3_p2_p3 <- function(Ymis, c2, c3, M = 2000,
                                          em_args = list(maxit = 200, tol = 1e-6, ridge = 1e-6),
                                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(Ymis)
  
  fit <- do.call(mvn_em_missing, c(list(Y = Ymis), em_args))
  mu <- fit$mu
  Sigma <- fit$Sigma
  
  t2_sum <- 0
  t3_sum <- 0
  p2_sum <- 0
  p3_sum <- 0
  
  for (i in 1:n) {
    y <- Ymis[i,]
    obs <- which(!is.na(y))
    mis <- which(is.na(y))
    
    if (length(mis) == 0) {
      t2_sum <- t2_sum + unname(y[2])
      t3_sum <- t3_sum + unname(y[3])
      p2_sum <- p2_sum + as.numeric(unname(y[2]) < c2)
      p3_sum <- p3_sum + as.numeric(unname(y[3]) < c3)
      next
    }
    if (length(obs) == 0) stop("All entries missing: not expected here.")
    
    yobs  <- y[obs]
    muobs <- mu[obs]
    mumis <- mu[mis]
    
    S_oo <- Sigma[obs, obs, drop = FALSE]
    S_mo <- Sigma[mis, obs, drop = FALSE]
    S_om <- Sigma[obs, mis, drop = FALSE]
    S_mm <- Sigma[mis, mis, drop = FALSE]
    
    S_oo_inv <- safe_solve(S_oo)
    cond_mean_m <- as.vector(mumis + S_mo %*% S_oo_inv %*% (yobs - muobs))
    cond_cov_m  <- S_mm - S_mo %*% S_oo_inv %*% S_om
    cond_cov_m  <- make_pd(cond_cov_m)
    
    ymis_draw <- rmvnorm(M, mean = cond_mean_m, sigma = cond_cov_m)
    
    ## Y2 contribution
    if (is.na(y[2])) {
      pos2 <- which(mis == 2)
      y2_draw <- ymis_draw[, pos2]
      t2_sum <- t2_sum + mean(y2_draw)
      p2_sum <- p2_sum + mean(y2_draw < c2)
    } else {
      t2_sum <- t2_sum + unname(y[2])
      p2_sum <- p2_sum + as.numeric(unname(y[2]) < c2)
    }
    
    ## Y3 contribution
    if (is.na(y[3])) {
      pos3 <- which(mis == 3)
      y3_draw <- ymis_draw[, pos3]
      t3_sum <- t3_sum + mean(y3_draw)
      p3_sum <- p3_sum + mean(y3_draw < c3)
    } else {
      t3_sum <- t3_sum + unname(y[3])
      p3_sum <- p3_sum + as.numeric(unname(y[3]) < c3)
    }
  }
  
  c(theta2 = t2_sum / n,
    theta3 = t3_sum / n,
    P2     = p2_sum / n,
    P3     = p3_sum / n)
}

## -------------------------
## One replication (Full vs PFI; theta2, theta3, P2, P3)
## -------------------------
one_rep_full_pfi <- function(model, n = 500, M = 2000,
                             miss_frac = 0.25,
                             em_args = list(maxit = 200, tol = 1e-6, ridge = 1e-6),
                             seed = NULL,
                             check_missing = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  cuts <- get_cuts(model)
  c2 <- cuts["c2"]
  c3 <- cuts["c3"]
  
  Yfull <- gen_data(model, n)
  full_est <- c(theta2 = mean(Yfull[,2]),
                theta3 = mean(Yfull[,3]),
                P2     = mean(Yfull[,2] < c2),
                P3     = mean(Yfull[,3] < c3))
  
  Ymis <- impose_missing(Yfull, miss_frac = miss_frac)
  
  if (check_missing) {
    cat("Missing summary:", model, "\n")
    print(missing_summary(Ymis))
  }
  
  pfi_est <- pfi_point_theta2_theta3_p2_p3(Ymis, c2 = c2, c3 = c3,
                                           M = M, em_args = em_args, seed = NULL)
  
  list(full_est = unname(full_est),
       pfi_est  = unname(pfi_est))
}

## -------------------------
## Simulation driver (Full vs PFI)
## -------------------------
run_simulation_full_pfi <- function(models = c("M1","M2","M4"),
                                    B = 2000, n = 500, M = 2000,
                                    miss_frac = 0.25,
                                    ncores = max(1, detectCores() - 1),
                                    seed = 123,
                                    em_args = list(maxit = 200, tol = 1e-6, ridge = 1e-6),
                                    truth_rel_tol = 1e-10,
                                    truth_abs_tol = 0) {
  set.seed(seed)
  
  cat("Computing true parameters (Exact / 1D quadrature):\n")
  true_list <- vector("list", length(models))
  names(true_list) <- models
  for (m in models) {
    cat("  true_params_exact for", m, "...\n")
    tr <- true_params_exact(m, rel.tol = truth_rel_tol, abs.tol = truth_abs_tol)
    true_list[[m]] <- tr[c("theta2","theta3","P2","P3")]
  }
  cat("Done.\n\n")
  
  use_cluster <- (ncores >= 2)
  
  if (use_cluster) {
    cl <- makeCluster(ncores)
    on.exit(stopCluster(cl), add = TRUE)
    clusterEvalQ(cl, { library(mvtnorm) })
    
    to_export <- c(
      "one_rep_full_pfi",
      "gen_data","impose_missing","get_cuts","missing_summary",
      "pfi_point_theta2_theta3_p2_p3","mvn_em_missing",
      "expit","Sigma_AR","make_pd","safe_solve","log_mvn_marginal",
      "pexgauss","true_params_exact"
    )
    clusterExport(cl, varlist = to_export, envir = .GlobalEnv)
  }
  
  out <- list()
  std_names <- c("theta2","theta3","P2","P3")
  
  for (m in models) {
    cat("Running model:", m, "\n")
    seeds <- sample.int(1e9, B)
    
    if (use_cluster) {
      reps <- parLapply(cl, 1:B, function(b) {
        one_rep_full_pfi(m, n = n, M = M, miss_frac = miss_frac,
                         em_args = em_args, seed = seeds[b])
      })
    } else {
      pb <- txtProgressBar(min = 0, max = B, style = 3)
      reps <- vector("list", B)
      for (b in 1:B) {
        reps[[b]] <- one_rep_full_pfi(m, n = n, M = M, miss_frac = miss_frac,
                                      em_args = em_args, seed = seeds[b])
        setTxtProgressBar(pb, b)
      }
      close(pb)
      cat("\n")
    }
    
    full_mat <- do.call(rbind, lapply(reps, `[[`, "full_est"))
    pfi_mat  <- do.call(rbind, lapply(reps, `[[`, "pfi_est"))
    
    colnames(full_mat) <- std_names
    colnames(pfi_mat)  <- std_names
    
    out[[m]] <- list(true = true_list[[m]],
                     full = full_mat,
                     pfi  = pfi_mat)
  }
  
  out
}

## -------------------------
## Figure (paper-like): Full vs PFI (theta2, theta3, P2, P3)
##   - facet_grid(Model ~ Parameter): rows=Model, cols=Parameter
## -------------------------
make_figure_full_vs_pfi <- function(sim_out,
                                    models = c("M1","M2"),
                                    outfile = "Figure_like_full_vs_pfi_theta2_theta3_P2_P3.png",
                                    trim = 0.005) {
  models <- intersect(models, names(sim_out))
  if (length(models) == 0) stop("No requested models in sim_out.")
  
  ## ---- data
  df <- dplyr::bind_rows(lapply(models, function(m) {
    full <- as.data.frame(sim_out[[m]]$full); full$Method <- "Full"
    pfi  <- as.data.frame(sim_out[[m]]$pfi ); pfi$Method  <- "PFI"
    full$Model <- m; pfi$Model <- m
    dplyr::bind_rows(full, pfi)
  }))
  
  df$Method <- factor(df$Method, levels = c("Full","PFI"))
  df$Model  <- factor(df$Model, levels = models)
  
  df_long <- df |>
    tidyr::pivot_longer(cols = c(theta2, theta3, P2, P3),
                        names_to = "Parameter", values_to = "Estimate") |>
    dplyr::mutate(
      Parameter = factor(Parameter,
                         levels = c("theta2","theta3","P2","P3"),
                         labels = c("Theta 2","Theta 3","P2","P3"))
    )
  
  ## ---- (optional) plot-only winsorization within each panel (Model, Parameter)
  ##      trim=0 なら無効。外れ値で潰れる場合だけ 0.005〜0.01 推奨
  if (trim > 0) {
    df_long <- df_long |>
      dplyr::group_by(Model, Parameter) |>
      dplyr::mutate(
        lo = stats::quantile(Estimate, probs = trim, na.rm = TRUE, type = 7),
        hi = stats::quantile(Estimate, probs = 1 - trim, na.rm = TRUE, type = 7),
        Estimate_plot = pmin(pmax(Estimate, lo), hi)
      ) |>
      dplyr::ungroup()
  } else {
    df_long$Estimate_plot <- df_long$Estimate
  }
  
  ## ---- panel order: (M1: Theta2..P3) then (M2: Theta2..P3) ... (paper-like 2x4)
  param_levels <- levels(df_long$Parameter)
  panel_levels <- as.vector(outer(models, param_levels, paste, sep = "\n"))
  df_long$Panel <- factor(paste(df_long$Model, df_long$Parameter, sep = "\n"),
                          levels = panel_levels)
  
  true_df <- dplyr::bind_rows(lapply(models, function(m) {
    tr <- sim_out[[m]]$true
    out <- data.frame(
      Model = factor(m, levels = models),
      Parameter = factor(names(tr),
                         levels = c("theta2","theta3","P2","P3"),
                         labels = c("Theta 2","Theta 3","P2","P3")),
      True = as.numeric(tr)
    )
    out$Panel <- factor(paste(out$Model, out$Parameter, sep = "\n"),
                        levels = panel_levels)
    out
  }))
  
  ## ---- facet: panel-wise free_y (これが潰れ防止の本体)
  facet_fun <- if (requireNamespace("ggh4x", quietly = TRUE)) {
    ## ggh4x が入っていれば各パネルに y 軸を出せて論文っぽくなる
    function(...) ggh4x::facet_wrap2(..., axes = "y")
  } else {
    ggplot2::facet_wrap
  }
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Method, y = Estimate_plot)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.size = 0.5) +
    ggplot2::geom_hline(data = true_df, ggplot2::aes(yintercept = True),
                        linetype = "dashed", linewidth = 0.6, colour = "red") +
    facet_fun(~ Panel, nrow = length(models), ncol = 4, scales = "free_y") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
      strip.background = ggplot2::element_rect(fill = "grey90", color = "grey50"),
      strip.text = ggplot2::element_text(size = 11, lineheight = 0.95),
      panel.spacing = grid::unit(1.0, "lines"),
      plot.margin = grid::unit(c(8, 8, 8, 8), "pt")
    ) +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = paste0("Simulation 1 (", paste(models, collapse = ", "), "): Full vs PFI"),
      subtitle = "Parameters: Theta2, Theta3, P2, P3   |   Dashed lines = true values"
    )
  
  ggplot2::ggsave(outfile, p, width = 13.5, height = 6.8, dpi = 300)
  invisible(p)
}


## =========================================================
## RUN (FULL vs PFI)
## =========================================================

## 0) Missingness sanity check (optional)
# set.seed(1)
# Ytmp <- gen_data("M2", n = 500)
# Ytmp_mis <- impose_missing(Ytmp, miss_frac = 0.25)
# print(missing_summary(Ytmp_mis))

## 1) Main run: M1, M2, M3, M4
sim_all <- run_simulation_full_pfi(models = c("M1","M2","M3","M4"),
                                   B = 2000, n = 500, M = 2000,
                                   miss_frac = 0.25,
                                   ncores = 16, seed = 123)

## Figure1-like: M1, M2 (theta2, theta3, P2, P3)
make_figure_full_vs_pfi(sim_all, models = c("M1","M2"),
                        outfile = "Figure1_like_M1_M2_full_vs_pfi_theta2_theta3_P2_P3.png")

## Figure2-like: M3, M4 (theta2, theta3, P2, P3)
make_figure_full_vs_pfi(sim_all, models = c("M3","M4"),
                        outfile = "Figure2_like_M3_M4_full_vs_pfi_theta2_theta3_P2_P3.png")



saveRDS(sim_all, file = "sim_pfi.rds")




suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

make_figure_full_pfi_spfi <- function(sim_pfi, sim_spfi,
                                      models = c("M1","M2"),
                                      outfile = "Figure_combined_full_pfi_spfi.png",
                                      trim = 0.005) {
  models <- intersect(models, intersect(names(sim_pfi), names(sim_spfi)))
  if (length(models) == 0) stop("No requested models found in both sim objects.")
  
  ## Full は重複するので片方（ここでは SPFI 側）を採用（推定値自体は同一のはず）
  df <- dplyr::bind_rows(lapply(models, function(m) {
    full <- as.data.frame(sim_spfi[[m]]$full); full$Method <- "Full"
    pfi  <- as.data.frame(sim_pfi [[m]]$pfi ); pfi$Method  <- "PFI"
    spfi <- as.data.frame(sim_spfi[[m]]$spfi); spfi$Method <- "SPFI"
    
    full$Model <- m; pfi$Model <- m; spfi$Model <- m
    dplyr::bind_rows(full, pfi, spfi)
  }))
  
  df$Method <- factor(df$Method, levels = c("Full","PFI","SPFI"))
  df$Model  <- factor(df$Model,  levels = models)
  
  df_long <- df |>
    tidyr::pivot_longer(cols = c(theta2, theta3, P2, P3),
                        names_to = "Parameter", values_to = "Estimate") |>
    dplyr::mutate(
      Parameter = factor(Parameter,
                         levels = c("theta2","theta3","P2","P3"),
                         labels = c("Theta 2","Theta 3","P2","P3"))
    )
  
  ## plot-only winsorization（図の見やすさのため。元データは一切変更しない）
  if (trim > 0) {
    df_long <- df_long |>
      dplyr::group_by(Model, Parameter) |>
      dplyr::mutate(
        lo = stats::quantile(Estimate, probs = trim, na.rm = TRUE, type = 7),
        hi = stats::quantile(Estimate, probs = 1 - trim, na.rm = TRUE, type = 7),
        Estimate_plot = pmin(pmax(Estimate, lo), hi)
      ) |>
      dplyr::ungroup()
  } else {
    df_long$Estimate_plot <- df_long$Estimate
  }
  
  ## パネル順：モデル×パラメタを論文風に固定
  param_levels <- levels(df_long$Parameter)
  panel_levels <- as.vector(outer(models, param_levels, paste, sep = "\n"))
  df_long$Panel <- factor(paste(df_long$Model, df_long$Parameter, sep = "\n"),
                          levels = panel_levels)
  
  ## true 値（どちらから取っても同じ前提。ここでは SPFI 側）
  true_df <- dplyr::bind_rows(lapply(models, function(m) {
    tr <- sim_spfi[[m]]$true
    out <- data.frame(
      Model = factor(m, levels = models),
      Parameter = factor(names(tr),
                         levels = c("theta2","theta3","P2","P3"),
                         labels = c("Theta 2","Theta 3","P2","P3")),
      True = as.numeric(tr)
    )
    out$Panel <- factor(paste(out$Model, out$Parameter, sep = "\n"),
                        levels = panel_levels)
    out
  }))
  
  facet_fun <- if (requireNamespace("ggh4x", quietly = TRUE)) {
    function(...) ggh4x::facet_wrap2(..., axes = "y")
  } else {
    ggplot2::facet_wrap
  }
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Method, y = Estimate_plot)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.size = 0.5) +
    ggplot2::geom_hline(data = true_df, ggplot2::aes(yintercept = True),
                        linetype = "dashed", linewidth = 0.6, colour = "red") +
    facet_fun(~ Panel, nrow = length(models), ncol = 4, scales = "free_y") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
      strip.background = ggplot2::element_rect(fill = "grey90", color = "grey50"),
      strip.text = ggplot2::element_text(size = 11, lineheight = 0.95),
      panel.spacing = grid::unit(1.0, "lines"),
      plot.margin = grid::unit(c(8, 8, 8, 8), "pt")
    ) +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = paste0("Simulation 1 (", paste(models, collapse = ", "),
                     "): Full vs PFI vs SPFI"),
      subtitle = "Parameters: Theta2, Theta3, P2, P3   |   Dashed lines = true values"
    )
  
  ggplot2::ggsave(outfile, p, width = 13.5, height = 6.8, dpi = 300)
  invisible(p)
}

## ---- load outputs (別々に回した2本の結果)
sim_spfi <- readRDS("sim_spfi.rds")
sim_pfi  <- readRDS("sim_pfi.rds")

## ---- combined figures
make_figure_full_pfi_spfi(sim_pfi, sim_spfi, models = c("M1","M2"),
                          outfile = "Figure_combined_M1_M2_full_pfi_spfi.png")

make_figure_full_pfi_spfi(sim_pfi, sim_spfi, models = c("M3","M4"),
                          outfile = "Figure_combined_M3_M4_full_pfi_spfi.png")









make_figure_full_pfi_spfi_paper <- function(sim_pfi, sim_spfi,
                                            models = c("M1","M2"),
                                            outfile = "Figure_paper_like.png",
                                            trim = 0.005,
                                            show_x_all = FALSE) {
  suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(tidyr)
  })
  
  models <- intersect(models, intersect(names(sim_pfi), names(sim_spfi)))
  if (length(models) == 0) stop("No requested models found in both sim objects.")
  
  ## --- data (Full は sim_spfi 側を採用)
  df <- dplyr::bind_rows(lapply(models, function(m) {
    full <- as.data.frame(sim_spfi[[m]]$full); full$Method <- "Full"
    pfi  <- as.data.frame(sim_pfi [[m]]$pfi ); pfi$Method  <- "PFI"
    spfi <- as.data.frame(sim_spfi[[m]]$spfi); spfi$Method <- "SPFI"
    full$Model <- m; pfi$Model <- m; spfi$Model <- m
    dplyr::bind_rows(full, pfi, spfi)
  }))
  
  df$Method <- factor(df$Method, levels = c("Full","PFI","SPFI"))
  df$Model  <- factor(df$Model,  levels = models)
  
  df_long <- df |>
    tidyr::pivot_longer(cols = c(theta2, theta3, P2, P3),
                        names_to = "Parameter", values_to = "Estimate") |>
    dplyr::mutate(
      Parameter = factor(Parameter,
                         levels = c("theta2","theta3","P2","P3"),
                         labels = c("Theta 2","Theta 3","P2","P3"))
    )
  
  ## --- plot-only winsorization
  if (trim > 0) {
    df_long <- df_long |>
      dplyr::group_by(Model, Parameter) |>
      dplyr::mutate(
        lo = stats::quantile(Estimate, probs = trim, na.rm = TRUE, type = 7),
        hi = stats::quantile(Estimate, probs = 1 - trim, na.rm = TRUE, type = 7),
        Estimate_plot = pmin(pmax(Estimate, lo), hi)
      ) |>
      dplyr::ungroup()
  } else {
    df_long$Estimate_plot <- df_long$Estimate
  }
  
  ## --- Panel order: (M1: Theta2,Theta3,P2,P3) then (M2: ...)  ←論文順
  param_levels <- levels(df_long$Parameter)
  panel_levels <- as.vector(t(outer(models, param_levels, paste, sep = "\n")))
  df_long$Panel <- factor(paste(df_long$Model, df_long$Parameter, sep = "\n"),
                          levels = panel_levels)
  
  ## --- true values (ここでは sim_spfi 側)
  true_df <- dplyr::bind_rows(lapply(models, function(m) {
    tr <- sim_spfi[[m]]$true
    out <- data.frame(
      Model = factor(m, levels = models),
      Parameter = factor(names(tr),
                         levels = c("theta2","theta3","P2","P3"),
                         labels = c("Theta 2","Theta 3","P2","P3")),
      True = as.numeric(tr)
    )
    out$Panel <- factor(paste(out$Model, out$Parameter, sep = "\n"),
                        levels = panel_levels)
    out
  }))
  
  ## --- facet (Figure1: show_x_all=TRUE で上下両段に手法名)
  facet_layer <- if (requireNamespace("ggh4x", quietly = TRUE)) {
    axes_opt <- if (show_x_all) "all" else "y"
    ggh4x::facet_wrap2(~ Panel, nrow = length(models), ncol = 4,
                       scales = "free_y", axes = axes_opt)
  } else {
    if (show_x_all) {
      stop("Figure1 の上下両段に x 軸ラベルを出すには ggh4x が必要です。install.packages('ggh4x') を実行してください。")
    }
    ggplot2::facet_wrap(~ Panel, nrow = length(models), ncol = 4, scales = "free_y")
  }
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Method, y = Estimate_plot)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.size = 0.5) +
    ggplot2::geom_hline(data = true_df, ggplot2::aes(yintercept = True),
                        linetype = "dashed", linewidth = 0.6, colour = "red") +
    facet_layer +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
      strip.background = ggplot2::element_rect(fill = "grey90", color = "grey50"),
      strip.text = ggplot2::element_text(size = 11, lineheight = 0.95),
      panel.spacing = grid::unit(1.0, "lines"),
      plot.margin = grid::unit(c(8, 8, 8, 8), "pt")
    ) +
    ggplot2::labs(x = NULL, y = NULL,
                  subtitle = "Parameters: Theta2, Theta3, P2, P3   |   Dashed lines = true values")
  
  ## 縦サイズは行数に応じて少し調整
  h <- if (length(models) >= 2) 6.8 else 3.6
  ggplot2::ggsave(outfile, p, width = 13.5, height = h, dpi = 300)
  invisible(p)
}









sim_spfi <- readRDS("~/sim_spfi.rds")
sim_pfi  <- readRDS("~/sim_pfi.rds")

## Figure1: M1 & M2（論文の Figure1 と同じ並び + 上下両段に手法名）
make_figure_full_pfi_spfi_paper(sim_pfi, sim_spfi,
                                models = c("M1","M2"),
                                outfile = "Figure1_M1_M2_full_pfi_spfi_paper_order.png",
                                show_x_all = TRUE
)

## --- Figure2 用：M4 を “表示名 M3” に差し替えた一時オブジェクトを作る
sim_spfi_plot <- sim_spfi
sim_pfi_plot  <- sim_pfi

sim_spfi_plot[["M3"]] <- sim_spfi[["M4"]]
sim_pfi_plot [["M3"]] <- sim_pfi [["M4"]]

## Figure2: “M4 の結果”を “M3 と表示”して出力
make_figure_full_pfi_spfi_paper(sim_pfi_plot, sim_spfi_plot,
                                models  = c("M3"),
                                outfile = "Figure2_M4_full_pfi_spfi_paper_order.png",
                                show_x_all = FALSE)



