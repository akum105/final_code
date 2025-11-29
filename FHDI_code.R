## =========================================================
##  Kim (2004) 第一シミュレーションのレプリケーション
##  - FI(M=10): fractional hot deck + fractional jackknife
##  - ABB(M=10): approximate Bayesian bootstrap + MI 分散
## =========================================================

## --- 真のパラメータ値（理論値） ----------------------------------------

true_params <- function() {
  mu1  <- 2.8
  mu2  <- 3.8
  sd1  <- sqrt(1.16)
  sd2  <- sqrt(1.735)
  p1   <- 0.5       # cell 1 の確率
  p2   <- 0.5       # cell 2 の確率
  
  # θ1 = E(Y)
  theta1 <- p1 * mu1 + p2 * mu2
  
  # z ⟂ Y なので θ2 = E(Y | z=1) = E(Y)
  theta2 <- theta1
  
  # θ3 = P(Y < 2)（2つの正規の混合）
  theta3 <- 0.5 * pnorm(2, mean = mu1, sd = sd1) +
    0.5 * pnorm(2, mean = mu2, sd = sd2)
  
  c(theta1 = theta1, theta2 = theta2, theta3 = theta3)
}

theta_true <- true_params()
theta_true


## --- サンプル生成（Kim (2004) の設定） -----------------------------------

generate_sample <- function(n = 120) {
  # imputation cell（g=1,2）
  g <- rbinom(n, size = 1, prob = 0.5) + 1L
  
  # Y の生成
  y <- ifelse(
    g == 1L,
    rnorm(n, mean = 2.8, sd = sqrt(1.16)),
    rnorm(n, mean = 3.8, sd = sqrt(1.735))
  )
  
  # domain 指標 z ~ Bernoulli(0.25)
  z <- rbinom(n, size = 1, prob = 0.25)
  
  # 応答確率
  resp_prob <- ifelse(g == 1L, 0.7, 0.6)
  r <- rbinom(n, size = 1, prob = resp_prob)
  
  y_obs <- y
  y_obs[r == 0L] <- NA
  
  list(y = y, y_obs = y_obs, z = z, g = g, r = r)
}


## --- θ1, θ2, θ3 の推定（完全データ用） -----------------------------------

estimate_theta <- function(y_comp, z) {
  idx_dom <- which(z == 1L)
  
  theta1_hat <- mean(y_comp)
  theta2_hat <- mean(y_comp[idx_dom])
  theta3_hat <- mean(as.numeric(y_comp < 2))
  
  c(theta1 = theta1_hat,
    theta2 = theta2_hat,
    theta3 = theta3_hat)
}


## --- 完全データに対する設計分散（MI 内の W_m 用） ------------------------

complete_var <- function(y_comp, z) {
  n  <- length(y_comp)
  idx_dom <- which(z == 1L)
  nz <- length(idx_dom)
  
  # θ1: 全体平均
  s2  <- var(y_comp)
  v1  <- s2 / n
  
  # θ2: domain mean
  v2  <- if (nz > 1L) var(y_comp[idx_dom]) / nz else NA_real_
  
  # θ3: proportion
  p_hat <- mean(as.numeric(y_comp < 2))
  v3    <- p_hat * (1 - p_hat) / n
  
  c(theta1 = v1, theta2 = v2, theta3 = v3)
}


## --- FI(M) : fractional version + donor 情報保持 --------------------------

impute_FI_M_frac <- function(y_obs, g, z, M = 10L) {
  n     <- length(y_obs)
  cells <- sort(unique(g))
  
  ## θ1, θ2 用の単一代入ベクトル
  y_imp <- y_obs
  
  ## θ3 用：各ユニットの “fractional indicator” ũ_i
  u_tilde <- numeric(n)
  
  ## 非欠測ユニットについてはそのままインジケータ
  for (i in seq_len(n)) {
    if (!is.na(y_obs[i])) {
      u_tilde[i] <- as.numeric(y_obs[i] < 2)
    }
  }
  
  ## 欠測ユニットの index と donor index の対応を保存
  recip_idx   <- integer(0)  # 欠測ユニットの元のインデックス i
  donors_list <- list()      # donors_list[[k]] が recip_idx[k] の donor ベクトル
  
  miss_counter <- 0L
  
  ## セルごとに FI(M) を実行
  for (cell in cells) {
    idx_cell <- which(g == cell)
    idx_resp <- idx_cell[!is.na(y_obs[idx_cell])]
    idx_miss <- idx_cell[is.na(y_obs[idx_cell])]
    
    r_g <- length(idx_resp)
    m_g <- length(idx_miss)
    if (m_g == 0L) next
    
    ## total donor uses = M * m_g = t_g * r_g + k_g
    total_uses <- M * m_g
    t_g <- total_uses %/% r_g
    k_g <- total_uses - t_g * r_g
    
    uses <- rep(t_g, r_g)
    if (k_g > 0L) {
      extra <- sample(seq_len(r_g), size = k_g, replace = FALSE)
      uses[extra] <- uses[extra] + 1L
    }
    
    donor_vec <- rep(idx_resp, times = uses)
    donor_vec <- sample(donor_vec, length(donor_vec), replace = FALSE)
    
    ## 欠測ユニットごとに M 個の donor を割り当て
    for (j in seq_len(m_g)) {
      miss_i   <- idx_miss[j]
      miss_counter <- miss_counter + 1L
      
      donors_j <- donor_vec[((j - 1L) * M + 1L):(j * M)]
      
      ## θ1, θ2 用：donor の平均で単一代入
      y_imp[miss_i] <- mean(y_obs[donors_j], na.rm = TRUE)
      
      ## θ3 用：fractional indicator ũ_i
      ind_donors       <- as.numeric(y_obs[donors_j] < 2)
      u_tilde[miss_i]  <- mean(ind_donors)
      
      ## donor 情報を保存
      recip_idx[miss_counter]      <- miss_i
      donors_list[[miss_counter]]  <- donors_j
    }
  }
  
  ## θ1, θ2：単一代入 y_imp から計算
  theta1_hat <- mean(y_imp)
  theta2_hat <- mean(y_imp[z == 1L])
  
  ## θ3：fractional indicator の平均
  theta3_hat <- mean(u_tilde)
  
  theta_hat <- c(theta1 = theta1_hat,
                 theta2 = theta2_hat,
                 theta3 = theta3_hat)
  
  list(
    y_imp     = y_imp,
    u_tilde   = u_tilde,
    theta_hat = theta_hat,
    recip_idx = recip_idx,
    donors    = donors_list
  )
}


## --- FI(M) に対する fractional jackknife ---------------------------------

var_jackknife_FI <- function(y_obs, g, z, M = 10L) {
  n <- length(y_obs)
  
  ## 1回だけ FI(M) を実行し，擬似完全データとドナー構造を作る
  full <- impute_FI_M_frac(y_obs, g, z, M = M)
  y_tilde   <- full$y_imp      # θ1, θ2 用の completed y
  u_tilde   <- full$u_tilde    # θ3 用の completed indicator
  theta_hat <- full$theta_hat  # full-sample estimator
  
  recip_idx <- full$recip_idx  # 欠測ユニットの index ベクトル
  donors    <- full$donors     # donors[[k]] は recip_idx[k] の donor index ベクトル
  
  ## leave-one-out 推定量を格納
  theta_jk <- matrix(NA_real_, nrow = n, ncol = 3L,
                     dimnames = list(NULL, names(theta_hat)))
  
  ## 各ユニットを1つずつ落とした jackknife replicate を作る
  for (k in seq_len(n)) {
    keep   <- (seq_len(n) != k)
    idx_ke <- which(keep)   # 残る元 index
    
    ## この replicate における completed y, u を作る
    y_rep <- numeric(length(idx_ke))
    u_rep <- numeric(length(idx_ke))
    
    for (pos in seq_along(idx_ke)) {
      i0 <- idx_ke[pos]  # 元の index
      
      if (!is.na(y_obs[i0])) {
        ## 元々応答していたユニット
        y_rep[pos] <- y_obs[i0]
        u_rep[pos] <- as.numeric(y_obs[i0] < 2)
      } else {
        ## 欠測ユニット → donors 情報から再計算
        kk <- which(recip_idx == i0)
        if (length(kk) != 1L) {
          stop("recipient index lookup failed in var_jackknife_FI")
        }
        donors_i <- donors[[kk]]
        
        ## この replicate で削除された k を donor から除く
        donors_i2 <- donors_i[donors_i != k]
        
        if (length(donors_i2) == 0L) {
          ## 全 donor が落ちた稀なケース：full の completed 値をそのまま使う
          y_rep[pos] <- y_tilde[i0]
          u_rep[pos] <- u_tilde[i0]
        } else {
          y_rep[pos] <- mean(y_obs[donors_i2], na.rm = TRUE)
          u_rep[pos] <- mean(as.numeric(y_obs[donors_i2] < 2), na.rm = TRUE)
        }
      }
    }
    
    ## z も delete-1 に対応させる
    z_rep <- z[keep]
    
    ## jackknife replicate の θ1, θ2, θ3
    theta1_loo <- mean(y_rep)
    theta2_loo <- mean(y_rep[z_rep == 1L])
    theta3_loo <- mean(u_rep)
    
    theta_jk[k, ] <- c(theta1 = theta1_loo,
                       theta2 = theta2_loo,
                       theta3 = theta3_loo)
  }
  
  ## delete-1 jackknife 分散
  theta_bar <- colMeans(theta_jk)
  diffs     <- sweep(theta_jk, 2L, theta_bar, FUN = "-")
  var_hat   <- (n - 1) / n * colSums(diffs^2)
  
  list(theta_hat = theta_hat, var_hat = var_hat)
}


## --- ABB(M) + MI variance --------------------------------------------------

impute_ABB_MI <- function(y_obs, g, z, M = 10L) {
  n     <- length(y_obs)
  cells <- sort(unique(g))
  
  theta_m <- matrix(NA_real_, nrow = M, ncol = 3L,
                    dimnames = list(NULL, c("theta1", "theta2", "theta3")))
  W_m <- matrix(NA_real_, nrow = M, ncol = 3L,
                dimnames = dimnames(theta_m))
  
  for (m in seq_len(M)) {
    y_imp <- y_obs
    
    for (cell in cells) {
      idx_cell <- which(g == cell)
      idx_resp <- idx_cell[!is.na(y_obs[idx_cell])]
      idx_miss <- idx_cell[is.na(y_obs[idx_cell])]
      
      r_g <- length(idx_resp)
      m_g <- length(idx_miss)
      if (m_g == 0L || r_g == 0L) next
      
      # ABB: 応答者からサイズ r_g のブートストラップサンプル
      boot_donors <- sample(idx_resp, size = r_g, replace = TRUE)
      
      # 欠測1つにつき1 donor をブートサンプルから選ぶ
      donors_for_missing <- sample(boot_donors, size = m_g, replace = TRUE)
      y_imp[idx_miss] <- y_obs[donors_for_missing]
    }
    
    theta_m[m, ] <- estimate_theta(y_imp, z)
    W_m[m, ]     <- complete_var(y_imp, z)
  }
  
  theta_bar <- colMeans(theta_m)
  W_bar     <- colMeans(W_m)
  B         <- apply(theta_m, 2L, var)
  
  V_MI <- W_bar + (1 + 1 / M) * B
  
  list(theta_hat = theta_bar, var_hat = V_MI)
}


## --- シミュレーション本体 --------------------------------------------------

run_simulation <- function(n = 120L,
                           R = 3000L,   # 本番は 30000L など
                           M = 10L,
                           seed = 123L) {
  set.seed(seed)
  
  methods   <- c("FI10", "ABB10")
  par_names <- c("theta1", "theta2", "theta3")
  
  est <- array(NA_real_,
               dim = c(R, length(par_names), length(methods)),
               dimnames = list(NULL, par_names, methods))
  
  varhat <- array(NA_real_,
                  dim = c(R, length(par_names), length(methods)),
                  dimnames = dimnames(est))
  
  for (r in seq_len(R)) {
    samp <- generate_sample(n)
    
    ## FI(M=10) + fractional jackknife
    fi <- var_jackknife_FI(y_obs = samp$y_obs,
                           g     = samp$g,
                           z     = samp$z,
                           M     = M)
    est[r, , "FI10"]    <- fi$theta_hat
    varhat[r, , "FI10"] <- fi$var_hat
    
    ## ABB(M=10) + MI variance
    abb <- impute_ABB_MI(y_obs = samp$y_obs,
                         g     = samp$g,
                         z     = samp$z,
                         M     = M)
    est[r, , "ABB10"]    <- abb$theta_hat
    varhat[r, , "ABB10"] <- abb$var_hat
    
    if (r %% 500L == 0L) {
      cat("rep", r, "of", R, "\n")
    }
  }
  
  list(est = est, varhat = varhat)
}


## --- Table 1, 2 の集計 -----------------------------------------------------

summarize_table1 <- function(sim_result) {
  est      <- sim_result$est
  methods  <- dimnames(est)[[3]]
  params   <- dimnames(est)[[2]]
  
  out <- data.frame()
  for (m in methods) {
    for (p in params) {
      vals    <- est[, p, m]
      mean_est <- mean(vals)
      var_est  <- var(vals)
      
      out <- rbind(out,
                   data.frame(
                     Method    = m,
                     Parameter = p,
                     Mean      = mean_est,
                     Variance  = var_est
                   ))
    }
  }
  out
}

summarize_table2 <- function(sim_result, table1) {
  est     <- sim_result$est
  varhat  <- sim_result$varhat
  methods <- dimnames(est)[[3]]
  params  <- dimnames(est)[[2]]
  
  out <- data.frame()
  for (m in methods) {
    for (p in params) {
      V_true <- table1$Variance[
        table1$Method == m & table1$Parameter == p
      ]
      vhat      <- varhat[, p, m]
      mean_vhat <- mean(vhat)
      bias      <- mean_vhat - V_true
      
      RB <- 100 * bias / mean_vhat    # 論文と同じ定義
      RV <- var(vhat) / (V_true^2)
      
      out <- rbind(out,
                   data.frame(
                     Method    = m,
                     Parameter = p,
                     Mean      = mean_vhat,
                     RB        = RB,
                     RV        = RV
                   ))
    }
  }
  out
}


## --- 実行例 ---------------------------------------------------------------

sim120  <- run_simulation(n = 120L, R = 30000L, M = 10L, seed = 1L)

tab1_120 <- summarize_table1(sim120)
tab2_120 <- summarize_table2(sim120, tab1_120)

tab1_120_sub <- subset(
  tab1_120,
  Method %in% c("FI10", "ABB10") &
    Parameter %in% c("theta1", "theta2", "theta3")
)

tab2_120_sub <- subset(
  tab2_120,
  Method %in% c("FI10", "ABB10") &
    Parameter %in% c("theta1", "theta2", "theta3")
)

tab1_120_sub
tab2_120_sub
