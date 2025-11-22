#--------------------------------------------------
# 母集団データを1回分生成する関数
#   n       : サンプルサイズ
#   params  : 真のパラメータを入れた list（省略可）
# 戻り値   : data.frame(x, y1, y2) （欠測なし）
#--------------------------------------------------

simulate_population_once <- function(
    n,
    params = list(
      beta0 = 1.0,
      beta1 = 0.7,
      sigma_e = 1.0,   # σ_ee
      phi0  = -3.0,
      phi1  = 0.5,
      phi2  = 0.7
    )
) {
  # 1. 共変量 x ~ N(2, 1)
  x <- rnorm(n, mean = 2, sd = 1)
  
  # 2. 連続応答 y1 | x ~ N(beta0 + beta1 * x, sigma_e^2)
  mu_y1 <- params$beta0 + params$beta1 * x
  y1 <- rnorm(n, mean = mu_y1, sd = params$sigma_e)
  
  # 3. 二値応答 y2 | (x, y1) ~ Bernoulli(p)
  #    logit(p) = phi0 + phi1 * x + phi2 * y1
  linpred <- params$phi0 + params$phi1 * x + params$phi2 * y1
  p <- 1 / (1 + exp(-linpred))  # inverse logit
  
  # rbinom で Bernoulli(p) を発生
  y2 <- rbinom(n, size = 1, prob = p)
  
  # 4. data.frame で返す
  dat <- data.frame(
    x  = x,
    y1 = y1,
    y2 = y2
  )
  
  return(dat)
}

set.seed(123)               # 再現性のため
n <- 200                    # 文書と同じ標本サイズ

dat_full <- simulate_population_once(n)

head(dat_full)
summary(dat_full)
colMeans(dat_full)          # 大まかに平均を見る（理論値と近いか確認）

#--------------------------------------------------
# 欠測を付与する関数
#   dat_full : data.frame(x, y1, y2) （完全データ）
# 戻り値    : 
#   data.frame(
#     x, y1, y2, 
#     delta1, delta2, 
#     y1_obs, y2_obs
#   )
#   - y1_obs, y2_obs は「観測された値（欠測は NA）」
#--------------------------------------------------

apply_missingness <- function(dat_full) {
  n <- nrow(dat_full)
  x  <- dat_full$x
  y1 <- dat_full$y1
  y2 <- dat_full$y2
  
  ## 1. y1 の応答指標 δ1 ~ Ber(π), logit(π) = 0.5 * x
  logit_pi <- 0.5 * x
  pi <- 1 / (1 + exp(-logit_pi))
  delta1 <- rbinom(n, size = 1, prob = pi)
  
  ## 2. y2 の応答指標 δ2 ~ Ber(0.7)
  delta2 <- rbinom(n, size = 1, prob = 0.7)
  
  ## 3. 観測値（欠測は NA にする）
  y1_obs <- ifelse(delta1 == 1, y1, NA)
  y2_obs <- ifelse(delta2 == 1, y2, NA)
  
  out <- data.frame(
    x      = x,
    y1     = y1,      # 真の y1（完全データ）
    y2     = y2,      # 真の y2（完全データ）
    delta1 = delta1,  # y1 の応答指標
    delta2 = delta2,  # y2 の応答指標
    y1_obs = y1_obs,  # 観測された y1（欠測は NA）
    y2_obs = y2_obs   # 観測された y2（欠測は NA）
  )
  
  return(out)
}

set.seed(123)

n <- 200

## ステップ1：完全データ生成
dat_full <- simulate_population_once(n)

## ステップ2：欠測付与
dat_miss <- apply_missingness(dat_full)

head(dat_miss)
summary(dat_miss[, c("x", "y1_obs", "y2_obs", "delta1", "delta2")])

## y1 の応答率（だいたい 0.72 付近になるはず）
mean(dat_miss$delta1)

## y2 の応答率（だいたい 0.7 付近）
mean(dat_miss$delta2)
