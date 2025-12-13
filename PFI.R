##############################
## Kim (2011) settings: MI vs PFI(FI) in ONE script
## Output: relative bias of variance estimator (%) for MI and FI
## Focus: eta1=E(y1), eta2=E(y2), eta4=Pr(y1<3)  (eta3 excluded)
##
## NOTE: With M=100 and B=2000 this can be computationally very heavy,
## especially FI (MH inside). For a quick test, set B <- 50 or 100 first.
##############################

set.seed(123)

## ===== 0. True parameters =====
beta0_true  <- 1
beta1_true  <- 0.7
sigma2_true <- 1
phi0_true   <- -3
phi1_true   <- 0.5
phi2_true   <- 0.7

## ===== 1. Data generation + missingness (Kim 2011) =====
generate_full_data <- function(n) {
  x  <- rnorm(n, mean = 2, sd = 1)
  y1 <- rnorm(n, mean = beta0_true + beta1_true * x, sd = sqrt(sigma2_true))
  lin <- phi0_true + phi1_true * x + phi2_true * y1
  p   <- plogis(lin)
  y2  <- rbinom(n, 1, p)
  data.frame(x=x, y1=y1, y2=y2)
}

impose_missing <- function(dat) {
  n  <- nrow(dat)
  pi1 <- plogis(0.5 * dat$x)   # logit(pi1)=0.5*x
  d1  <- rbinom(n, 1, pi1)     # y1 response indicator
  d2  <- rbinom(n, 1, 0.7)     # y2 response indicator (independent)
  y1_obs <- ifelse(d1==1, dat$y1, NA)
  y2_obs <- ifelse(d2==1, dat$y2, NA)
  cbind(dat, d1=d1, d2=d2, y1_obs=y1_obs, y2_obs=y2_obs)
}

## ===== 2. Target estimands (eta1, eta2, eta4 only) =====
complete_estimates_3 <- function(dat_full) {
  c(
    eta1 = mean(dat_full$y1),
    eta2 = mean(dat_full$y2),
    eta4 = mean(dat_full$y1 < 3)
  )
}

## ===== 3. MI (M = 100) =====
if (!requireNamespace("mvtnorm", quietly = TRUE)) {
  stop("Package 'mvtnorm' is required. Please install.packages('mvtnorm').")
}

draw_posterior_lm <- function(y, X) {
  n <- length(y)
  p <- ncol(X)
  XtX <- crossprod(X)
  XtY <- crossprod(X, y)
  beta_hat <- solve(XtX, XtY)
  res  <- y - X %*% beta_hat
  s2   <- sum(res^2) / (n - p)
  
  df_post <- n - p
  sigma2_draw <- df_post * s2 / rchisq(1, df_post)
  
  V_beta <- sigma2_draw * solve(XtX)
  beta_draw <- as.numeric(mvtnorm::rmvnorm(1, mean = beta_hat, sigma = V_beta))
  list(beta = beta_draw, sigma2 = sigma2_draw)
}

draw_posterior_logit <- function(y, X) {
  fit <- suppressWarnings(glm(y ~ X - 1, family = binomial))
  phi_hat <- coef(fit)
  I_inv   <- vcov(fit)
  as.numeric(mvtnorm::rmvnorm(1, mean = phi_hat, sigma = I_inv))
}

mi_fit_3 <- function(dat_miss, M = 100, burnin = 20, between = 1) {
  n <- nrow(dat_miss)
  
  idx_y1_obs <- which(!is.na(dat_miss$y1_obs))
  idx_y2_obs <- which(!is.na(dat_miss$y2_obs) & !is.na(dat_miss$y1_obs))
  
  ## initial values (respondents only)
  X1 <- cbind(1, dat_miss$x[idx_y1_obs])
  y1 <- dat_miss$y1_obs[idx_y1_obs]
  lm_init <- lm(y1 ~ X1 - 1)
  beta0   <- coef(lm_init)[1]
  beta1   <- coef(lm_init)[2]
  sigma2  <- summary(lm_init)$sigma^2
  
  X2 <- cbind(1, dat_miss$x[idx_y2_obs], dat_miss$y1_obs[idx_y2_obs])
  y2 <- dat_miss$y2_obs[idx_y2_obs]
  glm_init <- suppressWarnings(glm(y2 ~ X2 - 1, family = binomial))
  phi <- coef(glm_init)
  
  ## initialize missing values
  y1_imp <- dat_miss$y1_obs
  y2_imp <- dat_miss$y2_obs
  miss_y1 <- is.na(y1_imp)
  miss_y2 <- is.na(y2_imp)
  
  mu_y1_init <- beta0 + beta1 * dat_miss$x
  y1_imp[miss_y1] <- rnorm(sum(miss_y1), mu_y1_init[miss_y1], sqrt(sigma2))
  
  lin_y2_init <- phi[1] + phi[2] * dat_miss$x + phi[3] * y1_imp
  p_y2_init   <- plogis(lin_y2_init)
  y2_imp[miss_y2] <- rbinom(sum(miss_y2), 1, p_y2_init[miss_y2])
  
  completed_list <- vector("list", M)
  kept <- 0
  iter <- 0
  
  while (kept < M) {
    iter <- iter + 1
    
    ## P-step
    X1_all <- cbind(1, dat_miss$x[!is.na(y1_imp)])
    y1_all <- y1_imp[!is.na(y1_imp)]
    post1  <- draw_posterior_lm(y1_all, X1_all)
    beta0  <- post1$beta[1]
    beta1  <- post1$beta[2]
    sigma2 <- post1$sigma2
    
    X2_all <- cbind(1, dat_miss$x[!is.na(y2_imp)], y1_imp[!is.na(y2_imp)])
    y2_all <- y2_imp[!is.na(y2_imp)]
    phi    <- draw_posterior_logit(y2_all, X2_all)
    
    ## I-step
    mu_y1 <- beta0 + beta1 * dat_miss$x
    y1_imp[miss_y1] <- rnorm(sum(miss_y1), mu_y1[miss_y1], sqrt(sigma2))
    
    lin_y2 <- phi[1] + phi[2] * dat_miss$x + phi[3] * y1_imp
    p_y2   <- plogis(lin_y2)
    y2_imp[miss_y2] <- rbinom(sum(miss_y2), 1, p_y2[miss_y2])
    
    if (iter > burnin && ((iter - burnin) %% between == 0)) {
      kept <- kept + 1
      completed_list[[kept]] <- data.frame(x=dat_miss$x, y1=y1_imp, y2=y2_imp)
    }
  }
  
  ## Rubin combining rules (eta1, eta2, eta4)
  est_mat <- matrix(NA, M, 3)
  U_mat   <- matrix(NA, M, 3)
  colnames(est_mat) <- colnames(U_mat) <- c("eta1","eta2","eta4")
  
  for (m in 1:M) {
    dat_c <- completed_list[[m]]
    est_mat[m, ] <- complete_estimates_3(dat_c)
    
    n <- nrow(dat_c)
    U1 <- var(dat_c$y1) / n
    U2 <- var(dat_c$y2) / n
    U4 <- var(dat_c$y1 < 3) / n
    U_mat[m, ] <- c(U1, U2, U4)
  }
  
  est_bar <- colMeans(est_mat)
  U_bar   <- colMeans(U_mat)
  B_vec   <- apply(est_mat, 2, var)
  V_hat   <- U_bar + (1 + 1/M) * B_vec
  
  list(est=est_bar, var=V_hat)
}

## ===== 4. PFI(FI) (Kim 2011 Algorithm 1 + var eq(13) as in your FI code) =====
log_f1_norm <- function(y1, x, beta0, beta1, sigma2) {
  mu <- beta0 + beta1 * x
  dnorm(y1, mean=mu, sd=sqrt(sigma2), log=TRUE)
}
log_f2_bern <- function(y2, x, y1, phi0, phi1, phi2) {
  eta <- phi0 + phi1*x + phi2*y1
  y2*eta - log1p(exp(eta))
}

scores_all_6 <- function(y1_mat, y2_mat, x_vec, theta) {
  beta0 <- theta[1]; beta1 <- theta[2]; sigma2 <- theta[3]
  phi0  <- theta[4]; phi1  <- theta[5]; phi2  <- theta[6]
  
  n <- length(x_vec); M <- ncol(y1_mat)
  x_mat <- matrix(x_vec, n, M)
  
  mu_mat <- beta0 + beta1 * x_mat
  r_mat  <- y1_mat - mu_mat
  
  s_b0  <- r_mat / sigma2
  s_b1  <- x_mat * r_mat / sigma2
  s_s2  <- -0.5/sigma2 + 0.5*(r_mat^2)/(sigma2^2)
  
  eta_mat <- phi0 + phi1*x_mat + phi2*y1_mat
  p_mat   <- plogis(eta_mat)
  u_mat   <- (y2_mat - p_mat)
  s_p0 <- u_mat
  s_p1 <- x_mat * u_mat
  s_p2 <- y1_mat * u_mat
  
  list(s_b0=s_b0, s_b1=s_b1, s_s2=s_s2, s_p0=s_p0, s_p1=s_p1, s_p2=s_p2)
}

draw_y1_cond_mh <- function(x, y2_obs, theta1_prop, theta2_prop, M,
                            burn=200, thin=5) {
  beta0 <- theta1_prop[1]; beta1 <- theta1_prop[2]; sigma2 <- theta1_prop[3]
  phi0  <- theta2_prop[1]; phi1  <- theta2_prop[2]; phi2  <- theta2_prop[3]
  
  mu <- beta0 + beta1*x
  sd <- sqrt(sigma2)
  
  y_cur <- rnorm(1, mu, sd)
  loglik2 <- function(y) {
    eta <- phi0 + phi1*x + phi2*y
    y2_obs*eta - log1p(exp(eta))
  }
  
  keep <- numeric(M)
  t_keep <- 0
  T <- burn + M*thin
  
  for (t in 1:T) {
    y_prop <- rnorm(1, mu, sd)
    lr <- loglik2(y_prop) - loglik2(y_cur)
    if (log(runif(1)) < lr) y_cur <- y_prop
    
    if (t > burn && ((t - burn) %% thin == 0)) {
      t_keep <- t_keep + 1
      keep[t_keep] <- y_cur
    }
  }
  keep
}

init_theta2_imputed_score <- function(dat_miss, M=200, max_iter=50, tol=1e-9) {
  obs_y1 <- !is.na(dat_miss$y1_obs)
  obs_y2 <- !is.na(dat_miss$y2_obs)
  
  idx_y1 <- which(obs_y1)
  if (length(idx_y1) == 0L) stop("No y1 respondents (delta1=1).")
  
  idx_cc <- which(obs_y1 & obs_y2)
  if (length(idx_cc) < 5L) stop("Too few complete cases for initial logistic regression.")
  
  fit0 <- suppressWarnings(glm(y2_obs ~ x + y1_obs,
                               family=binomial, data=dat_miss, subset=idx_cc))
  theta_old <- as.numeric(coef(fit0))
  if (anyNA(theta_old)) stop("Initial logistic fit failed (NA coefficients).")
  
  for (it in 1:max_iter) {
    phi0 <- theta_old[1]; phi1 <- theta_old[2]; phi2 <- theta_old[3]
    
    x_i  <- dat_miss$x[idx_y1]
    y1_i <- dat_miss$y1_obs[idx_y1]
    d2_i <- dat_miss$d2[idx_y1]
    y2_i <- dat_miss$y2_obs[idx_y1]
    
    eta_i <- phi0 + phi1*x_i + phi2*y1_i
    p_i   <- plogis(eta_i)
    
    y2_rep <- numeric(0); x_rep <- numeric(0); y1_rep <- numeric(0); w_rep <- numeric(0)
    
    if (any(d2_i == 1)) {
      ii <- which(d2_i == 1)
      y2_rep <- c(y2_rep, y2_i[ii])
      x_rep  <- c(x_rep,  x_i[ii])
      y1_rep <- c(y1_rep, y1_i[ii])
      w_rep  <- c(w_rep,  rep(1, length(ii)))
    }
    
    if (any(d2_i == 0)) {
      ii <- which(d2_i == 0)
      y2_imp <- rbinom(length(ii)*M, size=1, prob=rep(p_i[ii], each=M))
      y2_rep <- c(y2_rep, y2_imp)
      x_rep  <- c(x_rep,  rep(x_i[ii], each=M))
      y1_rep <- c(y1_rep, rep(y1_i[ii], each=M))
      w_rep  <- c(w_rep,  rep(1/M, length(ii)*M))
    }
    
    fit <- suppressWarnings(glm(y2_rep ~ x_rep + y1_rep,
                                family=binomial, weights=w_rep))
    theta_new <- as.numeric(coef(fit))
    if (anyNA(theta_new)) break
    
    if (max(abs(theta_new - theta_old)) < tol) {
      theta_old <- theta_new
      break
    }
    theta_old <- theta_new
  }
  theta_old
}

init_theta0 <- function(dat_miss, M=200) {
  obs_y1 <- !is.na(dat_miss$y1_obs)
  
  lm0 <- lm(y1_obs ~ x, data=dat_miss, subset=obs_y1)
  beta0_0 <- coef(lm0)[1]
  beta1_0 <- coef(lm0)[2]
  
  ## sigma2_0 as MLE (no df correction)
  res0 <- residuals(lm0)
  n_obs1 <- sum(obs_y1)
  sigma2_0 <- sum(res0^2) / n_obs1
  
  theta1_0 <- c(beta0_0, beta1_0, sigma2_0)
  theta2_0 <- init_theta2_imputed_score(dat_miss, M=M, max_iter=50, tol=1e-9)
  
  list(theta1_0=theta1_0, theta2_0=theta2_0)
}

generate_imputations_conditional <- function(dat_miss, theta1_prop, theta2_prop, M,
                                             mh_burn=200, mh_thin=5) {
  n <- nrow(dat_miss)
  obs_y1 <- !is.na(dat_miss$y1_obs)
  obs_y2 <- !is.na(dat_miss$y2_obs)
  
  beta0 <- theta1_prop[1]; beta1 <- theta1_prop[2]; sigma2 <- theta1_prop[3]
  phi0  <- theta2_prop[1]; phi1  <- theta2_prop[2]; phi2  <- theta2_prop[3]
  
  y1_star <- matrix(NA, n, M)
  y2_star <- matrix(NA, n, M)
  
  for (i in 1:n) {
    x_i <- dat_miss$x[i]
    
    ## y1*
    if (obs_y1[i]) {
      y1_star[i,] <- dat_miss$y1_obs[i]
    } else {
      if (obs_y2[i]) {
        y1_star[i,] <- draw_y1_cond_mh(
          x=x_i, y2_obs=dat_miss$y2_obs[i],
          theta1_prop=theta1_prop, theta2_prop=theta2_prop,
          M=M, burn=mh_burn, thin=mh_thin
        )
      } else {
        mu_i <- beta0 + beta1*x_i
        y1_star[i,] <- rnorm(M, mu_i, sqrt(sigma2))
      }
    }
    
    ## y2*
    if (obs_y2[i]) {
      y2_star[i,] <- dat_miss$y2_obs[i]
    } else {
      eta_ij <- phi0 + phi1*x_i + phi2*y1_star[i,]
      p_ij   <- plogis(eta_ij)
      y2_star[i,] <- rbinom(M, 1, p_ij)
    }
  }
  
  list(y1_star=y1_star, y2_star=y2_star)
}

fi_fit_M <- function(dat_miss, M=100, C=5, max_iter=200, tol=1e-9,
                     mh_burn=200, mh_thin=5, use_step4=TRUE) {
  
  n <- nrow(dat_miss)
  obs_y1 <- !is.na(dat_miss$y1_obs)
  obs_y2 <- !is.na(dat_miss$y2_obs)
  
  ## SRS weights
  w_i <- rep(1/n, n)
  
  ## Step 1: theta(0)
  init <- init_theta0(dat_miss, M=max(M,200))
  theta1_0 <- init$theta1_0
  theta2_0 <- init$theta2_0
  
  ## Step 2: conditional proposal imputations
  imp <- generate_imputations_conditional(dat_miss, theta1_0, theta2_0, M,
                                          mh_burn=mh_burn, mh_thin=mh_thin)
  y1_star <- imp$y1_star
  y2_star <- imp$y2_star
  
  theta1 <- theta1_0
  theta2 <- theta2_0
  w_ij   <- matrix(1/M, n, M)
  
  for (iter in 1:max_iter) {
    beta0 <- theta1[1]; beta1 <- theta1[2]; sigma2 <- theta1[3]
    phi0  <- theta2[1]; phi1  <- theta2[2]; phi2  <- theta2[3]
    
    x_mat <- matrix(dat_miss$x, n, M)
    
    log_num <- log_f1_norm(y1_star, x_mat, beta0, beta1, sigma2) +
      log_f2_bern(y2_star, x_mat, y1_star, phi0, phi1, phi2)
    
    log_den <- log_f1_norm(y1_star, x_mat, theta1_0[1], theta1_0[2], theta1_0[3]) +
      log_f2_bern(y2_star, x_mat, y1_star, theta2_0[1], theta2_0[2], theta2_0[3])
    
    log_ratio <- log_num - log_den
    row_max <- apply(log_ratio, 1, max)
    ratio <- exp(log_ratio - row_max)
    row_sum <- rowSums(ratio)
    row_sum[row_sum == 0] <- 1
    w_ij <- ratio / row_sum
    
    ## Step 4 cap
    if (use_step4 && max(w_ij) > C/M) {
      theta1_0 <- theta1
      theta2_0 <- theta2
      imp <- generate_imputations_conditional(dat_miss, theta1_0, theta2_0, M,
                                              mh_burn=mh_burn, mh_thin=mh_thin)
      y1_star <- imp$y1_star
      y2_star <- imp$y2_star
      w_ij <- matrix(1/M, n, M)
      next
    }
    
    ## Step 5 M-step
    
    ## theta1 update
    x_rep1  <- c(dat_miss$x[obs_y1], rep(dat_miss$x[!obs_y1], each=M))
    y1_rep1 <- c(dat_miss$y1_obs[obs_y1], as.vector(t(y1_star[!obs_y1,,drop=FALSE])))
    
    w_rep1_obs <- w_i[obs_y1]
    w_rep1_mis <- as.vector(t(w_ij[!obs_y1,,drop=FALSE])) * rep(w_i[!obs_y1], each=M)
    w_rep1 <- c(w_rep1_obs, w_rep1_mis)
    
    lm_w <- lm(y1_rep1 ~ x_rep1, weights=w_rep1)
    beta0_new <- coef(lm_w)[1]
    beta1_new <- coef(lm_w)[2]
    res1 <- residuals(lm_w)
    sigma2_new <- sum(w_rep1 * res1^2) / sum(w_rep1)  # MLE form
    
    ## theta2 update
    idx_A <- which(obs_y2 & obs_y1)
    x_A  <- dat_miss$x[idx_A]
    y1_A <- dat_miss$y1_obs[idx_A]
    y2_A <- dat_miss$y2_obs[idx_A]
    w_A  <- w_i[idx_A]
    
    idx_B <- which(obs_y2 & !obs_y1)
    x_B  <- rep(dat_miss$x[idx_B], each=M)
    y1_B <- as.vector(t(y1_star[idx_B,,drop=FALSE]))
    y2_B <- rep(dat_miss$y2_obs[idx_B], each=M)
    w_B  <- as.vector(t(w_ij[idx_B,,drop=FALSE])) * rep(w_i[idx_B], each=M)
    
    idx_C <- which(!obs_y2)
    x_C  <- rep(dat_miss$x[idx_C], each=M)
    y1_C <- as.vector(t(y1_star[idx_C,,drop=FALSE]))
    y2_C <- as.vector(t(y2_star[idx_C,,drop=FALSE]))
    w_C  <- as.vector(t(w_ij[idx_C,,drop=FALSE])) * rep(w_i[idx_C], each=M)
    
    x_rep2  <- c(x_A,  x_B,  x_C)
    y1_rep2 <- c(y1_A, y1_B, y1_C)
    y2_rep2 <- c(y2_A, y2_B, y2_C)
    w_rep2  <- c(w_A,  w_B,  w_C)
    
    glm_w <- suppressWarnings(glm(y2_rep2 ~ x_rep2 + y1_rep2,
                                  family=binomial, weights=w_rep2))
    phi0_new <- coef(glm_w)[1]
    phi1_new <- coef(glm_w)[2]
    phi2_new <- coef(glm_w)[3]
    
    diff <- max(abs(c(beta0_new,beta1_new,sigma2_new) - theta1),
                abs(c(phi0_new,phi1_new,phi2_new) - theta2))
    
    theta1 <- c(beta0_new, beta1_new, sigma2_new)
    theta2 <- c(phi0_new, phi1_new, phi2_new)
    
    if (diff < tol) break
  }
  
  theta_hat <- c(theta1, theta2)
  
  ## point estimates (eta1, eta2, eta4)
  g1_mat <- y1_star
  g2_mat <- y2_star
  g4_mat <- (y1_star < 3) * 1
  
  gbar1 <- rowSums(w_ij * g1_mat)
  gbar2 <- rowSums(w_ij * g2_mat)
  gbar4 <- rowSums(w_ij * g4_mat)
  
  est <- c(
    eta1 = sum(w_i * gbar1),
    eta2 = sum(w_i * gbar2),
    eta4 = sum(w_i * gbar4)
  )
  
  ## variance (Kim 2011 eq. 13, SRS Omega)
  n <- nrow(dat_miss)
  Omega_diag <- 1/(n^2)
  Omega_off  <- -1/(n^2*(n-1))
  
  sc <- scores_all_6(y1_star, y2_star, dat_miss$x, theta_hat)
  scomp_list <- list(sc$s_b0, sc$s_b1, sc$s_s2, sc$s_p0, sc$s_p1, sc$s_p2)
  
  sbar <- sapply(scomp_list, function(S) rowSums(w_ij * S))
  colnames(sbar) <- c("b0","b1","s2","p0","p1","p2")
  
  varhat_eq13_one_g <- function(g_mat) {
    gbar <- rowSums(w_ij * g_mat)
    
    p <- ncol(sbar)
    M1 <- matrix(0, p, p)
    M2 <- rep(0, p)
    
    for (i in 1:n) {
      si_bar <- sbar[i,]
      M1 <- M1 + w_i[i] * (si_bar %*% t(si_bar))
      
      for (k in 1:p) {
        S_k <- scomp_list[[k]][i,]
        M2[k] <- M2[k] + w_i[i] * sum(w_ij[i,] * (S_k - si_bar[k]) * g_mat[i,])
      }
    }
    
    ## safer solve (still algebraically K1 = M1^{-1} M2)
    K1_hat <- as.numeric(tryCatch(qr.solve(M1, M2), error = function(e) {
      eps <- 1e-10
      solve(M1 + diag(eps, nrow(M1)), M2)
    }))
    
    e_i <- as.numeric(gbar + sbar %*% K1_hat)
    
    sum_e  <- sum(e_i)
    sum_e2 <- sum(e_i^2)
    sum_off <- (sum_e^2 - sum_e2)
    
    Omega_diag * sum_e2 + Omega_off * sum_off
  }
  
  varh <- c(
    eta1 = varhat_eq13_one_g(g1_mat),
    eta2 = varhat_eq13_one_g(g2_mat),
    eta4 = varhat_eq13_one_g(g4_mat)
  )
  
  list(est=est, var=varh)
}

## ===== 5. Monte Carlo: MI & FI together; output rel. bias of variance estimator =====
run_sim_MI_FI <- function(B=2000, n=200, M=100,
                          mi_burnin=20, mi_between=1,
                          fi_C=5, fi_tol=1e-9, fi_mh_burn=200, fi_mh_thin=5,
                          fi_use_step4=TRUE, verbose=TRUE) {
  
  est_mi <- matrix(NA, B, 3, dimnames=list(NULL, c("eta1","eta2","eta4")))
  var_mi <- matrix(NA, B, 3, dimnames=list(NULL, c("eta1","eta2","eta4")))
  est_fi <- matrix(NA, B, 3, dimnames=list(NULL, c("eta1","eta2","eta4")))
  var_fi <- matrix(NA, B, 3, dimnames=list(NULL, c("eta1","eta2","eta4")))
  
  for (b in 1:B) {
    if (verbose) cat("Replication", b, "/", B, "\r")
    
    dat_full <- generate_full_data(n)
    dat_miss <- impose_missing(dat_full)
    
    ## MI
    mi <- mi_fit_3(dat_miss, M=M, burnin=mi_burnin, between=mi_between)
    est_mi[b,] <- mi$est
    var_mi[b,] <- mi$var
    
    ## FI (PFI)
    fi <- fi_fit_M(dat_miss, M=M, C=fi_C, tol=fi_tol,
                   mh_burn=fi_mh_burn, mh_thin=fi_mh_thin,
                   use_step4=fi_use_step4)
    est_fi[b,] <- fi$est
    var_fi[b,] <- fi$var
  }
  if (verbose) cat("\n")
  
  ## Relative bias of variance estimator: 100*(mean(varhat)/MCVar - 1)
  mi_mc_var <- apply(est_mi, 2, var)
  fi_mc_var <- apply(est_fi, 2, var)
  
  mi_rel_bias_var <- 100 * (colMeans(var_mi) / mi_mc_var - 1)
  fi_rel_bias_var <- 100 * (colMeans(var_fi) / fi_mc_var - 1)
  
  list(
    mi_rel_bias_var = mi_rel_bias_var,
    fi_rel_bias_var = fi_rel_bias_var
  )
}

## ===== 6. Execute =====
B <- 2000
n <- 200
M <- 100

out <- run_sim_MI_FI(B=B, n=n, M=M, verbose=TRUE)

cat("\n=== Relative bias of variance estimator (%) ===\n")
cat("\nMI (M=100):\n")
print(round(out$mi_rel_bias_var, 2))

cat("\nPFI/FI (M=100):\n")
print(round(out$fi_rel_bias_var, 2))