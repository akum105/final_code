#--------------------------------------------
# 1回レプリケーションをまとめて実行
#   n      : サンプルサイズ
#   M_pfi  : PFI の M
#   M_mi   : MI の M
#
# 戻り値: list(
#   truth = c(EY1, EY2, P_lt3, beta1),
#   PFI   = c(EY1, EY2, P_lt3, beta1),
#   MI    = c(EY1, EY2, P_lt3, beta1),
#   var_MI = c(EY1, EY2, PY1_lt3, beta1)  # Rubin の分散推定
# )
#--------------------------------------------
one_replication <- function(n = 200, M_pfi = 100, M_mi = 20) {
  ## 完全データ
  dat_full <- simulate_population_once(n)
  
  ## 欠測付与
  dat_miss <- apply_missingness(dat_full)
  
  ## PFI
  res_pfi <- run_PFI(dat_miss, M = M_pfi, verbose = FALSE)
  
  ## MI
  res_mi  <- run_MI(dat_miss, M = M_mi)
  
  ## 完全データからの「真値に対応する量」
  EY1_true   <- mean(dat_full$y1)
  EY2_true   <- mean(dat_full$y2)
  P_lt3_true <- mean(dat_full$y1 < 3)
  beta1_true <- coef(lm(y1 ~ x, data = dat_full))[2]
  
  truth_vec <- c(
    EY1   = EY1_true,
    EY2   = EY2_true,
    P_lt3 = P_lt3_true,
    beta1 = beta1_true
  )
  
  ## PFI 推定値
  PFI_vec <- c(
    EY1   = res_pfi$EY1,
    EY2   = res_pfi$EY2,
    P_lt3 = res_pfi$PY1_lt3,
    beta1 = res_pfi$beta1_hat
  )
  
  ## MI 推定値
  MI_vec <- c(
    EY1   = res_mi$EY1_MI,
    EY2   = res_mi$EY2_MI,
    P_lt3 = res_mi$PY1_lt3_MI,
    beta1 = res_mi$beta1_MI
  )
  
  ## MI の Rubin 分散推定（run_MI で返しているもの）
  var_MI_vec <- res_mi$var_MI[c("EY1", "EY2", "PY1_lt3", "beta1")]
  names(var_MI_vec) <- c("EY1", "EY2", "P_lt3", "beta1")
  
  list(
    truth  = truth_vec,
    PFI    = PFI_vec,
    MI     = MI_vec,
    var_MI = var_MI_vec
  )
}


#--------------------------------------------
# Monte Carlo ループ (PFI vs MI)
#
# B      : レプリケーション回数
# n      : サンプルサイズ
# M_pfi  : PFI の M
# M_mi   : MI の M
#
# 戻り値: list(
#   truth  = B×4 行列 (EY1, EY2, P_lt3, beta1),
#   PFI    = B×4 行列,
#   MI     = B×4 行列,
#   var_MI = B×4 行列 (Rubin 分散推定)
# )
#--------------------------------------------
run_MC <- function(B = 2000, n = 200, M_pfi = 100, M_mi = 20) {
  par_names <- c("EY1", "EY2", "P_lt3", "beta1")
  
  truth_mat  <- matrix(NA_real_, nrow = B, ncol = length(par_names),
                       dimnames = list(NULL, par_names))
  PFI_mat    <- matrix(NA_real_, nrow = B, ncol = length(par_names),
                       dimnames = list(NULL, par_names))
  MI_mat     <- matrix(NA_real_, nrow = B, ncol = length(par_names),
                       dimnames = list(NULL, par_names))
  var_MI_mat <- matrix(NA_real_, nrow = B, ncol = length(par_names),
                       dimnames = list(NULL, par_names))
  
  for (b in 1:B) {
    res <- one_replication(n = n, M_pfi = M_pfi, M_mi = M_mi)
    
    truth_mat[b, ]  <- res$truth
    PFI_mat[b, ]    <- res$PFI
    MI_mat[b, ]     <- res$MI
    var_MI_mat[b, ] <- res$var_MI
    
    if (b %% 50 == 0) {
      cat("b =", b, " / ", B, "\n")
    }
  }
  
  list(
    truth  = truth_mat,
    PFI    = PFI_mat,
    MI     = MI_mat,
    var_MI = var_MI_mat
  )
}


#--------------------------------------------
# Monte Carlo 結果の集計
#
# sim_res : run_MC() の戻り値
#
# 戻り値 : list(
#   summary_PFI = データフレーム,
#   summary_MI  = データフレーム
# )
#--------------------------------------------
summarize_MC <- function(sim_res) {
  par_names <- colnames(sim_res$truth)
  
  # 結果を入れる枠
  summary_PFI <- data.frame(
    param = par_names,
    bias  = NA_real_,
    mc_var = NA_real_,
    rmse   = NA_real_
  )
  
  summary_MI <- data.frame(
    param      = par_names,
    bias       = NA_real_,
    mc_var     = NA_real_,
    rmse       = NA_real_,
    avg_varhat = NA_real_,  # 推定分散の平均
    cov95      = NA_real_   # 95% CI 被覆率
  )
  
  for (j in seq_along(par_names)) {
    name <- par_names[j]
    
    theta_true <- sim_res$truth[, name]
    hat_PFI    <- sim_res$PFI[,   name]
    hat_MI     <- sim_res$MI[,    name]
    var_MI_hat <- sim_res$var_MI[, name]
    
    ## PFI
    bias_PFI <- mean(hat_PFI - theta_true)
    var_PFI  <- var(hat_PFI)
    rmse_PFI <- sqrt(bias_PFI^2 + var_PFI)
    
    summary_PFI$bias[j]  <- bias_PFI
    summary_PFI$mc_var[j] <- var_PFI
    summary_PFI$rmse[j]  <- rmse_PFI
    
    ## MI
    bias_MI <- mean(hat_MI - theta_true)
    var_MI  <- var(hat_MI)
    rmse_MI <- sqrt(bias_MI^2 + var_MI)
    
    # Rubin 分散推定の平均
    avg_varhat <- mean(var_MI_hat)
    
    # 95% CI の被覆率（正規近似）
    lower <- hat_MI - 1.96 * sqrt(var_MI_hat)
    upper <- hat_MI + 1.96 * sqrt(var_MI_hat)
    cov95 <- mean((theta_true >= lower) & (theta_true <= upper))
    
    summary_MI$bias[j]       <- bias_MI
    summary_MI$mc_var[j]     <- var_MI
    summary_MI$rmse[j]       <- rmse_MI
    summary_MI$avg_varhat[j] <- avg_varhat
    summary_MI$cov95[j]      <- cov95
  }
  
  list(
    summary_PFI = summary_PFI,
    summary_MI  = summary_MI
  )
}

set.seed(123)

## まずはテスト用に B=200 程度から
sim_res <- run_MC(B = 200, n = 200, M_pfi = 100, M_mi = 20)

mc_summary <- summarize_MC(sim_res)

mc_summary$summary_PFI
mc_summary$summary_MI

