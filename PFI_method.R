run_PFI <- function(dat_miss, M = 100, max_iter = 50, tol = 1e-6, verbose = FALSE) {
  n <- nrow(dat_miss)
  
  x      <- dat_miss$x
  y1_obs <- dat_miss$y1_obs
  y2_obs <- dat_miss$y2_obs
  
  ## -----------------------------
  ## Step 0: 初期パラメータ推定
  ## -----------------------------
  
  idx_y1_obs <- !is.na(y1_obs)
  fit_y1 <- lm(y1_obs ~ x, subset = idx_y1_obs)
  beta0_init <- coef(fit_y1)[1]
  beta1_init <- coef(fit_y1)[2]
  res1 <- resid(fit_y1)
  sigma_e_init <- sqrt(sum(res1^2) / length(res1))
  
  idx_y2_obs  <- !is.na(y2_obs)
  idx_both    <- idx_y1_obs & idx_y2_obs
  fit_y2 <- glm(y2_obs ~ x + y1_obs, family = binomial(), subset = idx_both)
  phi0_init <- coef(fit_y2)[1]
  phi1_init <- coef(fit_y2)[2]
  phi2_init <- coef(fit_y2)[3]
  
  theta1_0 <- c(beta0 = beta0_init, beta1 = beta1_init, sigma_e = sigma_e_init)
  theta2_0 <- c(phi0  = phi0_init,  phi1  = phi1_init,  phi2    = phi2_init)
  
  ## ------------------------------------------------
  ## 初期パラメータ θ^(0) の下で M 個の擬似値を固定生成
  ## ------------------------------------------------
  
  y1_imp <- matrix(NA_real_, nrow = n, ncol = M)
  y2_imp <- matrix(NA_real_, nrow = n, ncol = M)
  
  for (i in 1:n) {
    if (!is.na(y1_obs[i])) {
      y1_imp[i, ] <- rep(y1_obs[i], M)
    } else {
      mu_i <- beta0_init + beta1_init * x[i]
      y1_imp[i, ] <- rnorm(M, mean = mu_i, sd = sigma_e_init)
    }
    
    if (!is.na(y2_obs[i])) {
      y2_imp[i, ] <- rep(y2_obs[i], M)
    } else {
      linpred <- phi0_init + phi1_init * x[i] + phi2_init * y1_imp[i, ]
      p_ij <- 1 / (1 + exp(-linpred))
      y2_imp[i, ] <- rbinom(M, size = 1, prob = p_ij)
    }
  }
  
  ## -----------------------------
  ## E–M 反復
  ## -----------------------------
  
  theta1 <- theta1_0
  theta2 <- theta2_0
  
  diff <- Inf
  iter <- 0
  
  for (iter in 1:max_iter) {
    beta0 <- theta1["beta0"]
    beta1 <- theta1["beta1"]
    sigma_e <- theta1["sigma_e"]
    
    phi0 <- theta2["phi0"]
    phi1 <- theta2["phi1"]
    phi2 <- theta2["phi2"]
    
    mu_curr <- beta0 + beta1 * x
    mu_init <- theta1_0["beta0"] + theta1_0["beta1"] * x
    mu_curr_mat <- matrix(mu_curr, nrow = n, ncol = M)
    mu_init_mat <- matrix(mu_init, nrow = n, ncol = M)
    
    logf1_curr <- dnorm(y1_imp, mean = mu_curr_mat, sd = sigma_e, log = TRUE)
    logf1_init <- dnorm(y1_imp, mean = mu_init_mat, sd = theta1_0["sigma_e"], log = TRUE)
    
    eta_curr <- phi0 + phi1 * x + phi2 * y1_imp
    eta_init <- theta2_0["phi0"] + theta2_0["phi1"] * x + theta2_0["phi2"] * y1_imp
    
    p_curr <- 1 / (1 + exp(-eta_curr))
    p_init <- 1 / (1 + exp(-eta_init))
    
    logf2_curr <- dbinom(y2_imp, size = 1, prob = p_curr, log = TRUE)
    logf2_init <- dbinom(y2_imp, size = 1, prob = p_init, log = TRUE)
    
    logw <- (logf1_curr + logf2_curr) - (logf1_init + logf2_init)
    
    w <- matrix(NA_real_, nrow = n, ncol = M)
    for (i in 1:n) {
      li <- logw[i, ]
      finite_idx <- is.finite(li)
      
      if (!any(finite_idx)) {
        w[i, ] <- rep(1 / M, M)
      } else {
        maxlog <- max(li[finite_idx])
        wi <- rep(0, M)
        wi[finite_idx] <- exp(li[finite_idx] - maxlog)
        s <- sum(wi)
        if (s == 0 || !is.finite(s)) {
          w[i, ] <- rep(1 / M, M)
        } else {
          w[i, ] <- wi / s
        }
      }
    }
    
    x_rep  <- as.vector(matrix(x, nrow = n, ncol = M))
    y1_rep <- as.vector(y1_imp)
    y2_rep <- as.vector(y2_imp)
    w_rep  <- as.vector(w)
    
    df1 <- data.frame(y1 = y1_rep, x = x_rep)
    fit1 <- lm(y1 ~ x, data = df1, weights = w_rep)
    beta0_new <- coef(fit1)[1]
    beta1_new <- coef(fit1)[2]
    
    mu_hat <- beta0_new + beta1_new * x_rep
    resid1 <- y1_rep - mu_hat
    sigma_e_new <- sqrt(sum(w_rep * resid1^2) / sum(w_rep))
    
    theta1_new <- c(beta0 = beta0_new, beta1 = beta1_new, sigma_e = sigma_e_new)
    
    df2 <- data.frame(y2 = y2_rep, x = x_rep, y1 = y1_rep)
    fit2 <- glm(y2 ~ x + y1, data = df2,
                family = quasibinomial(), weights = w_rep)
    coefs2 <- coef(fit2)
    theta2_new <- c(phi0 = coefs2[1], phi1 = coefs2[2], phi2 = coefs2[3])
    
    diff <- max(abs(c(theta1_new - theta1, theta2_new - theta2)))
    theta1 <- theta1_new
    theta2 <- theta2_new
    
    if (verbose) {
      cat(sprintf("iter = %d, max param diff = %.3e\n", iter, diff))
    }
    if (diff < tol) break
  }
  
  EY1     <- mean(rowSums(w * y1_imp))
  EY2     <- mean(rowSums(w * y2_imp))
  PY1_lt3 <- mean(rowSums(w * (y1_imp < 3)))
  
  idx_beta1 <- grep("beta1", names(theta1))
  beta1_hat <- as.numeric(theta1[idx_beta1])
  
  list(
    theta1     = theta1,
    theta2     = theta2,
    EY1        = EY1,
    EY2        = EY2,
    PY1_lt3    = PY1_lt3,
    beta1_hat  = beta1_hat,
    w          = w,
    y1_imp     = y1_imp,
    y2_imp     = y2_imp,
    iterations = iter,
    converged  = (diff < tol)
  )
}


set.seed(123)

n <- 200

## ステップ1：完全データ生成
dat_full <- simulate_population_once(n)

## ステップ2：欠測付与
dat_miss <- apply_missingness(dat_full)

## ステップ3：PFI を適用
res_pfi <- run_PFI(dat_miss, M = 100, verbose = TRUE)

res_pfi$theta1      # (beta0, beta1, sigma_e)
res_pfi$theta2      # (phi0, phi1, phi2)
res_pfi$EY1         # E(Y1) の PFI 推定
res_pfi$EY2         # E(Y2) の PFI 推定
res_pfi$PY1_lt3     # P(Y1 < 3) の PFI 推定
res_pfi$beta1_hat   # β1 の PFI 推定
res_pfi$converged   # TRUE なら収束

