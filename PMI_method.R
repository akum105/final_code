library(MASS)
#--------------------------------------------
# 線形回帰の擬似事後から (beta, sigma_e) を1回サンプル
#   y: 応答ベクトル
#   X: デザイン行列（切片列を含める）
# 戻り値:
#   list(beta = ベクトル, sigma_e = スカラー)
#--------------------------------------------
draw_posterior_linear <- function(y, X) {
  n <- length(y)
  p <- ncol(X)
  
  # OLS 推定値
  XtX <- crossprod(X)            # X'X
  XtX_inv <- solve(XtX)
  beta_hat <- XtX_inv %*% crossprod(X, y)
  
  # 残差と残差平方和
  resid <- y - X %*% beta_hat
  s2_hat <- as.numeric(crossprod(resid) / (n - p))   # σ²の不偏推定
  
  # σ² | y ~ (n-p)*s2_hat / χ²_{n-p}
  sigma2 <- (n - p) * s2_hat / rchisq(1, df = n - p)
  sigma_e <- sqrt(sigma2)
  
  # β | σ², y ~ N(β_hat, σ² (X'X)^{-1})
  beta <- as.vector( MASS::mvrnorm(1, mu = as.vector(beta_hat),
                                   Sigma = sigma2 * XtX_inv) )
  
  list(beta = beta, sigma_e = sigma_e)
}

#--------------------------------------------
# Multiple Imputation (MI)
#
# dat_miss : apply_missingness() の出力
# M        : インプリケーション数
#
# 戻り値:
#   list(
#     EY1_MI, EY2_MI, PY1_lt3_MI, beta1_MI,   # MI 点推定
#     var_MI = named vector,                  # Rubin のルールによる分散推定
#     Q_mat  = M×4 行列 (各 m の推定値),
#     U_mat  = M×4 行列 (各 m の完全データ分散)
#   )
#--------------------------------------------
run_MI <- function(dat_miss, M = 20) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' が必要です。install.packages('MASS') してください。")
  }
  
  n <- nrow(dat_miss)
  x      <- dat_miss$x
  y1_obs <- dat_miss$y1_obs
  y2_obs <- dat_miss$y2_obs
  
  # 各インプリケーションごとの推定値と完全データ分散を保存する行列
  Q_mat <- matrix(NA_real_, nrow = M, ncol = 4)
  colnames(Q_mat) <- c("EY1", "EY2", "PY1_lt3", "beta1")
  
  U_mat <- matrix(NA_real_, nrow = M, ncol = 4)
  colnames(U_mat) <- c("EY1", "EY2", "PY1_lt3", "beta1")
  
  for (m in 1:M) {
    ## 1. Y1 | X の線形回帰（観測 Y1 のみ）
    idx_y1_obs <- !is.na(y1_obs)
    X1 <- cbind(1, x[idx_y1_obs])               # 切片＋x
    y1_vec <- y1_obs[idx_y1_obs]
    
    post1 <- draw_posterior_linear(y = y1_vec, X = X1)
    beta_lin <- post1$beta      # (beta0, beta1)
    sigma_e  <- post1$sigma_e
    
    beta0 <- beta_lin[1]
    beta1 <- beta_lin[2]
    
    ## 2. サンプルした θ1 で Y1 の欠測を代入
    y1_imp <- y1_obs
    idx_y1_miss <- is.na(y1_obs)
    mu_y1_miss <- beta0 + beta1 * x[idx_y1_miss]
    y1_imp[idx_y1_miss] <- rnorm(sum(idx_y1_miss),
                                 mean = mu_y1_miss, sd = sigma_e)
    
    ## 3. Y2 | (X, Y1) のロジスティック回帰（観測 Y2 のみ）
    idx_y2_obs <- !is.na(y2_obs)
    dat_y2 <- data.frame(
      y2 = y2_obs[idx_y2_obs],
      x  = x[idx_y2_obs],
      y1 = y1_imp[idx_y2_obs]
    )
    fit_logit <- glm(y2 ~ x + y1, data = dat_y2, family = binomial())
    
    phi_hat <- coef(fit_logit)
    V_phi   <- vcov(fit_logit)
    
    # 正規近似事後から φ を1回サンプル
    phi <- as.numeric(
      MASS::mvrnorm(1, mu = phi_hat, Sigma = V_phi)
    )
    names(phi) <- names(phi_hat)
    
    phi0 <- phi[1]
    phi1 <- phi[2]
    phi2 <- phi[3]
    
    ## 4. サンプルした θ2 で Y2 の欠測を代入
    y2_imp <- y2_obs
    idx_y2_miss <- is.na(y2_obs)
    linpred_miss <- phi0 + phi1 * x[idx_y2_miss] + phi2 * y1_imp[idx_y2_miss]
    p_miss <- 1 / (1 + exp(-linpred_miss))
    y2_imp[idx_y2_miss] <- rbinom(sum(idx_y2_miss), size = 1, prob = p_miss)
    
    ## 5. 完成データ上での推定量と完全データ分散
    
    # EY1
    EY1_hat <- mean(y1_imp)
    var_EY1 <- var(y1_imp) / n
    
    # EY2
    EY2_hat <- mean(y2_imp)
    var_EY2 <- EY2_hat * (1 - EY2_hat) / n   # 二値の平均の分散 (n が十分大きい想定)
    
    # P(Y1 < 3)
    I_lt3 <- as.numeric(y1_imp < 3)
    PY1_lt3_hat <- mean(I_lt3)
    var_PY1_lt3 <- PY1_lt3_hat * (1 - PY1_lt3_hat) / n
    
    # β1（Y1 ~ X の線形回帰）
    fit_beta1 <- lm(y1_imp ~ x)
    beta1_hat <- coef(fit_beta1)[2]
    var_beta1 <- vcov(fit_beta1)[2, 2]
    
    # 保存
    Q_mat[m, ] <- c(EY1_hat, EY2_hat, PY1_lt3_hat, beta1_hat)
    U_mat[m, ] <- c(var_EY1, var_EY2, var_PY1_lt3, var_beta1)
  }
  
  ## Rubin のルールによる集約
  
  Q_bar <- colMeans(Q_mat)           # 点推定の平均
  U_bar <- colMeans(U_mat)           # 完成データ分散の平均
  B     <- apply(Q_mat, 2, var)      # 完成データ間分散（分母 M-1）
  
  # MI 分散: T = U_bar + (1 + 1/M) * B
  T_var <- U_bar + (1 + 1 / M) * B
  names(T_var) <- colnames(Q_mat)
  
  list(
    EY1_MI     = Q_bar["EY1"],
    EY2_MI     = Q_bar["EY2"],
    PY1_lt3_MI = Q_bar["PY1_lt3"],
    beta1_MI   = Q_bar["beta1"],
    var_MI     = T_var,
    Q_mat      = Q_mat,
    U_mat      = U_mat
  )
}

#PFIとMIの比較
set.seed(123)
n <- 200

dat_full <- simulate_population_once(n)
dat_miss <- apply_missingness(dat_full)

# PFI
res_pfi <- run_PFI(dat_miss, M = 100, verbose = TRUE)

# MI
res_mi  <- run_MI(dat_miss, M = 20)

c(
  EY1_true   = mean(dat_full$y1),
  EY1_PFI    = res_pfi$EY1,
  EY1_MI     = res_mi$EY1_MI
)

c(
  EY2_true   = mean(dat_full$y2),
  EY2_PFI    = res_pfi$EY2,
  EY2_MI     = res_mi$EY2_MI
)

c(
  P_lt3_true = mean(dat_full$y1 < 3),
  P_lt3_PFI  = res_pfi$PY1_lt3,
  P_lt3_MI   = res_mi$PY1_lt3_MI
)

c(
  beta1_true = 0.7,
  beta1_PFI  = res_pfi$beta1_hat,
  beta1_MI   = res_mi$beta1_MI
)
