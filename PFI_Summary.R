set.seed(20241122)   # 好きな seed で良いです

B_big   <- 2000
n       <- 200
M_pfi   <- 100
M_mi    <- 20

system.time({
  sim_res_big <- run_MC(B = B_big, n = n,
                        M_pfi = M_pfi, M_mi = M_mi)
})

mc_summary_big <- summarize_MC(sim_res_big)

mc_summary_big$summary_PFI
mc_summary_big$summary_MI


## 1. summary から PFI / MI を取り出す
tab_PFI <- mc_summary_big$summary_PFI
tab_MI  <- mc_summary_big$summary_MI

## 2. PFI 側に MI と同じ列名を追加 (NA で埋める)
tab_PFI$avg_varhat <- NA_real_
tab_PFI$cov95      <- NA_real_
tab_PFI$method     <- "PFI"

## 3. MI 側にも method 列を追加
tab_MI$method <- "MI"

## 4. 列の順番を揃えて rbind
cols <- c("param", "method", "bias", "mc_var", "rmse", "avg_varhat", "cov95")

tab_all <- rbind(
  tab_PFI[, cols],
  tab_MI [, cols]
)

## 5. 数値を丸めて見やすくする
round_tab <- within(tab_all, {
  bias       <- round(bias,       4)
  mc_var     <- round(mc_var,     4)
  rmse       <- round(rmse,       4)
  avg_varhat <- round(avg_varhat, 4)
  cov95      <- round(cov95,      3)
})

## 6. param, method の順で並べ替え
round_tab <- round_tab[order(round_tab$param, round_tab$method), ]

round_tab

