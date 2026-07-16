---
title: "R14：機器學習的金融因子選擇"
output:
  github_document:
    toc: true
    toc_depth: 3
---

本附錄對應第 19 章。它使用 1967 年 1 月至 2021 年 11 月、659 個月乘 10 產業的真實凍結因子面板，以第 \(t\) 月的預測變數預測製造業第 \(t+1\) 月超額報酬。`ret`、Fama--French 與 Global-q 因子均為月小數報酬，例如 0.01 代表 1%；總體預測變數則沿用各自定義的比率、利差或對數尺度。重點是時間對齊、只依訓練資料進行前處理、展開視窗驗證、基準比較與選取穩定性。

資料來源方面，原課程 `fffqmacro.R` 從 Kenneth French Data Library 取得 FF3 與十產業投資組合、從 Global-q 取得 q5 因子、從 Welch--Goyal 工作簿取得總體預測變數，並以 FRED `CPIAUCNS` 輔助建檔；`ret` 是十產業總報酬減去當月 FF 無風險利率。各供應者版本與再散布條款必須分別保存。尤其 Fama--French 建構法後來由 CRSP FIZ 改為 CIZ，現在即時重抓不保證等於本 2021 工作快照。以下只評估固定時間切分下的預測表現；因子被選入不識別因果效果或風險價格，也不構成投資建議。

## 執行條件

- 主線手算版使用 base R 與 `knitr`；原課程捷徑另使用 `glmnet`。
- 執行時不安裝套件、不呼叫網路、不使用 `setwd()`。
- 固定資料：`data/processed/ff_qf_macro_industries_1967_2021.csv`。
- 資料建置與再散布注意事項：`data/DATA_SOURCES.md`。
- 公開界線：repo 隨附作者授權的 processed 合併 CSV、程式與執行結果，可離線自含重跑；若另由上游來源重建，必須固定各供應者版本與轉換規則。


``` r
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.width = 7, fig.height = 4.5
)
if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("本附錄的原課程套件捷徑需要 glmnet；請先在合法可重現環境中安裝。")
}
```


``` r
locate_project_file <- function(relative_path) {
  candidates <- c(
    relative_path,
    file.path("..", relative_path),
    file.path("../..", relative_path)
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) stop("找不到專案檔案：", relative_path)
  normalizePath(hit[1], mustWork = TRUE)
}

path <- locate_project_file(
  "data/processed/ff_qf_macro_industries_1967_2021.csv"
)
d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
d$month <- as.Date(d$month)
d <- d[order(d$industry, d$month), ]

stopifnot(nrow(d) == 6590L, ncol(d) == 24L)
table(d$industry)
```

```
## 
## Durbl Enrgy HiTec  Hlth Manuf NoDur Other Shops Telcm Utils 
##   659   659   659   659   659   659   659   659   659   659
```

``` r
range(d$month)
```

```
## [1] "1967-01-01" "2021-11-01"
```

``` r
data.frame(
  first_month = min(d$month),
  last_month = max(d$month),
  months = length(unique(d$month)),
  industries = length(unique(d$industry)),
  return_unit = "decimal per month"
)
```

```
##   first_month last_month months industries       return_unit
## 1  1967-01-01 2021-11-01    659         10 decimal per month
```

## 日期與面板完整性


``` r
key <- paste(d$industry, d$month)
stopifnot(!anyDuplicated(key))

months_by_industry <- tapply(d$month, d$industry, length)
stopifnot(length(unique(months_by_industry)) == 1L)

common_predictors <- grep("^(factor_|macro_)", names(d), value = TRUE)
stopifnot(length(common_predictors) == 21L)

# 同一月份的 common predictors 應在十個產業完全相同。
check_common <- aggregate(
  d[, common_predictors],
  by = list(month = d$month),
  FUN = function(x) max(x) - min(x)
)
stopifnot(max(as.matrix(check_common[, -1]), na.rm = TRUE) < 1e-12)
```

## 建立 (X_t\rightarrow Y_{t+1})

本示範先固定 `Manuf`。若改產業，必須保持同一測試起點與調校規則。


``` r
industry_name <- "Manuf"
one <- d[d$industry == industry_name, ]
one <- one[order(one$month), ]

factor_names <- grep("^factor_", names(one), value = TRUE)
macro_names <- grep("^macro_", names(one), value = TRUE)

X_main <- as.matrix(one[, c(factor_names, macro_names)])
storage.mode(X_main) <- "double"

interaction_list <- vector("list", length(factor_names) * length(macro_names))
interaction_names <- character(length(interaction_list))
k <- 1L
for (f in seq_along(factor_names)) {
  for (m in seq_along(macro_names)) {
    interaction_list[[k]] <- X_main[, factor_names[f]] * X_main[, macro_names[m]]
    interaction_names[k] <- paste(factor_names[f], macro_names[m], sep = ":")
    k <- k + 1L
  }
}
X_int <- do.call(cbind, interaction_list)
colnames(X_int) <- interaction_names
X <- cbind(X_main, X_int)

# 第 t 月 predictors 預測 t+1 月報酬；最後一列沒有 target。
y <- c(one$ret[-1], NA_real_)
keep <- complete.cases(X, y)
X <- X[keep, , drop = FALSE]
y <- y[keep]
dates <- one$month[keep]

stopifnot(ncol(X) == 125L)
c(observations = nrow(X), dictionary_columns = ncol(X))
```

```
##       observations dictionary_columns 
##                658                125
```

``` r
head(data.frame(predictor_month = dates, target_next = y), 3)
```

```
##   predictor_month target_next
## 1      1967-01-01      0.0085
## 2      1967-02-01      0.0541
## 3      1967-03-01      0.0342
```

## LASSO 座標下降與只依訓練資料進行的標準化


``` r
soft_threshold <- function(z, lambda) sign(z) * pmax(abs(z) - lambda, 0)

standardize_train <- function(X_train, X_new = NULL) {
  mu <- colMeans(X_train)
  s <- apply(X_train, 2, sd)
  keep <- is.finite(s) & s > 1e-10
  train <- sweep(sweep(X_train[, keep, drop = FALSE], 2, mu[keep]),
                 2, s[keep], "/")
  ans <- list(train = train, mean = mu[keep], sd = s[keep], keep = keep)
  if (!is.null(X_new)) {
    ans$new <- sweep(sweep(X_new[, keep, drop = FALSE], 2, mu[keep]),
                     2, s[keep], "/")
  }
  ans
}

lasso_cd <- function(X, y_centered, lambda, max_iter = 5000L, tol = 1e-7) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- numeric(p)
  residual <- y_centered
  x2 <- colSums(X^2) / n
  for (iter in seq_len(max_iter)) {
    old <- beta
    for (j in seq_len(p)) {
      residual <- residual + X[, j] * beta[j]
      z <- sum(X[, j] * residual) / n
      beta[j] <- soft_threshold(z, lambda) / x2[j]
      residual <- residual - X[, j] * beta[j]
    }
    if (max(abs(beta - old)) < tol) break
  }
  attr(beta, "iterations") <- iter
  beta
}

fit_predict_lasso <- function(X_train, y_train, X_new, lambda) {
  pp <- standardize_train(X_train, X_new)
  y_bar <- mean(y_train)
  beta <- lasso_cd(pp$train, y_train - y_bar, lambda)
  list(
    pred = as.numeric(y_bar + pp$new %*% beta),
    beta = beta,
    names = colnames(pp$train),
    prep = pp,
    intercept = y_bar
  )
}
```

## 保留測試集，再用展開視窗折調校


``` r
n <- length(y)
test_start <- floor(0.80 * n) + 1L
tv <- seq_len(test_start - 1L)
test <- test_start:n

validation_size <- 24L
train_ends <- unique(as.integer(seq(
  floor(0.50 * length(tv)),
  length(tv) - validation_size,
  length.out = 5
)))
folds <- lapply(train_ends, function(e) {
  list(train = seq_len(e), validation = (e + 1L):(e + validation_size))
})

data.frame(
  train_end = dates[vapply(folds, function(z) max(z$train), integer(1))],
  validation_start = dates[vapply(folds, function(z) min(z$validation), integer(1))],
  validation_end = dates[vapply(folds, function(z) max(z$validation), integer(1))]
)
```

```
##    train_end validation_start validation_end
## 1 1988-11-01       1988-12-01     1990-11-01
## 2 1993-10-01       1993-11-01     1995-10-01
## 3 1998-10-01       1998-11-01     2000-10-01
## 4 2003-10-01       2003-11-01     2005-10-01
## 5 2008-10-01       2008-11-01     2010-10-01
```


``` r
pp0 <- standardize_train(X[folds[[1]]$train, , drop = FALSE])
y0 <- y[folds[[1]]$train]
lambda_max <- max(abs(crossprod(pp0$train, y0 - mean(y0)))) / length(y0)
lambda_grid <- exp(seq(log(lambda_max), log(lambda_max * 0.015), length.out = 30))

fold_loss <- matrix(NA_real_, nrow = length(folds), ncol = length(lambda_grid))
for (v in seq_along(folds)) {
  tr <- folds[[v]]$train
  va <- folds[[v]]$validation
  for (j in seq_along(lambda_grid)) {
    fit <- fit_predict_lasso(
      X[tr, , drop = FALSE], y[tr],
      X[va, , drop = FALSE], lambda_grid[j]
    )
    fold_loss[v, j] <- mean((y[va] - fit$pred)^2)
  }
}

mean_loss <- colMeans(fold_loss)
best_lambda <- lambda_grid[which.min(mean_loss)]
data.frame(best_lambda = best_lambda, validation_mse = min(mean_loss))
```

```
##   best_lambda validation_mse
## 1 0.003798291    0.002486418
```

## 原課程套件捷徑：以 `glmnet` 跑相同的時序折

原課程 `slides/L11_Factor selection_via_ML/fffqmacro.R` 第 381--397 行直接以
`glmnet()` 估 LASSO／Ridge 路徑，第 428--497 行再以時間序列折調校；精簡版
`lasso_ff_fq_macro.R` 第 85--116 行則示範以統一工作流程呼叫同一個 `glmnet`
引擎。本節採最少依賴的第一種作法：保留上面完全相同的 expanding-window
`folds`、`lambda_grid`、預測期距與最終測試集，只把手寫標準化與座標下降換成
`glmnet::glmnet()`。不能直接使用隨機 K-fold 的 `cv.glmnet()`，因為那會讓較晚
月份進入較早月份的訓練資料。


``` r
glmnet_fold_loss <- matrix(
  NA_real_, nrow = length(folds), ncol = length(lambda_grid)
)
for (v in seq_along(folds)) {
  tr <- folds[[v]]$train
  va <- folds[[v]]$validation
  glmnet_fold_fit <- glmnet::glmnet(
    x = X[tr, , drop = FALSE],
    y = y[tr],
    alpha = 1,
    lambda = lambda_grid,
    standardize = TRUE,
    intercept = TRUE,
    thresh = 1e-10
  )
  glmnet_fold_prediction <- predict(
    glmnet_fold_fit,
    newx = X[va, , drop = FALSE],
    s = lambda_grid
  )
  glmnet_fold_loss[v, ] <- colMeans(
    (matrix(y[va], nrow = length(va), ncol = length(lambda_grid)) -
       glmnet_fold_prediction)^2
  )
}

glmnet_mean_loss <- colMeans(glmnet_fold_loss)
glmnet_best_lambda <- lambda_grid[which.min(glmnet_mean_loss)]

data.frame(
  implementation = c("manual coordinate descent", "glmnet package"),
  best_lambda = c(best_lambda, glmnet_best_lambda),
  validation_MSE = c(min(mean_loss), min(glmnet_mean_loss)),
  same_time_folds = TRUE
)
```

```
##              implementation best_lambda validation_MSE same_time_folds
## 1 manual coordinate descent 0.003798291    0.002486418            TRUE
## 2            glmnet package 0.003798291    0.002486430            TRUE
```

兩個版本使用相同目標函數族與資料切分，但停止準則、標準化細節及路徑計算的
數值實作不必逐位元相同。公平的對照是驗證／測試預測與非零集合；手算版係數
位於標準化尺度，`coef.glmnet()` 則已轉回原始變數尺度，不能直接逐格相減。


``` r
plot(log(lambda_grid), mean_loss, type = "l", lwd = 2,
     xlab = "log(lambda)", ylab = "Mean validation MSE",
     col = "#173B57")
lines(log(lambda_grid), glmnet_mean_loss, lwd = 1.5,
      lty = 3, col = "#1D6D73")
abline(v = log(best_lambda), lty = 2, col = "#A34045")
abline(v = log(glmnet_best_lambda), lty = 3, col = "#1D6D73")
legend(
  "topleft", c("Manual coordinate descent", "glmnet"),
  col = c("#173B57", "#1D6D73"), lty = c(1, 3),
  lwd = c(2, 1.5), bty = "n"
)
```

![製造業 expanding-window validation loss。](../R14_financial_factor_selection_files/figure-gfm/cv-plot-1.png)

## 最終測試評量


``` r
final <- fit_predict_lasso(
  X[tv, , drop = FALSE], y[tv],
  X[test, , drop = FALSE], best_lambda
)

glmnet_final_fit <- glmnet::glmnet(
  x = X[tv, , drop = FALSE],
  y = y[tv],
  alpha = 1,
  lambda = lambda_grid,
  standardize = TRUE,
  intercept = TRUE,
  thresh = 1e-10
)
glmnet_prediction <- as.numeric(predict(
  glmnet_final_fit,
  newx = X[test, , drop = FALSE],
  s = glmnet_best_lambda
))

# 固定在「手算版選到的同一 lambda」再比一次，隔離調校差異與估計器差異。
glmnet_prediction_at_manual_lambda <- as.numeric(predict(
  glmnet_final_fit,
  newx = X[test, , drop = FALSE],
  s = best_lambda
))

# 固定 forecast origin 可實現的歷史平均基準。
baseline <- rep(mean(y[tv]), length(test))
actual <- y[test]

score <- function(actual, forecast, baseline) {
  mse <- mean((actual - forecast)^2)
  c(
    MSE = mse,
    MAE = mean(abs(actual - forecast)),
    OOS_R2 = 1 - mse / mean((actual - baseline)^2)
  )
}

rbind(
  HistoricalMean = score(actual, baseline, baseline),
  LASSO_manual = score(actual, final$pred, baseline),
  LASSO_glmnet = score(actual, glmnet_prediction, baseline)
)
```

```
##                        MSE        MAE      OOS_R2
## HistoricalMean 0.002112639 0.03413601  0.00000000
## LASSO_manual   0.002136697 0.03393833 -0.01138775
## LASSO_glmnet   0.002136831 0.03393859 -0.01145105
```


``` r
glmnet_coef <- as.matrix(stats::coef(
  glmnet_final_fit, s = glmnet_best_lambda
))[, 1]
glmnet_selected <- setdiff(
  names(glmnet_coef)[abs(glmnet_coef) > 1e-8], "(Intercept)"
)
manual_selected <- final$names[abs(final$beta) > 1e-8]

data.frame(
  comparison = c(
    "chosen-lambda test prediction",
    "same-manual-lambda test prediction",
    "selected-variable overlap"
  ),
  value = c(
    max(abs(final$pred - glmnet_prediction)),
    max(abs(final$pred - glmnet_prediction_at_manual_lambda)),
    length(intersect(manual_selected, glmnet_selected)) /
      max(1L, length(union(manual_selected, glmnet_selected)))
  ),
  interpretation = c(
    "includes any difference in selected lambda",
    "isolates numerical implementation at a fixed lambda",
    "Jaccard share of the two nonzero sets"
  )
)
```

```
##                           comparison       value
## 1      chosen-lambda test prediction 4.77403e-05
## 2 same-manual-lambda test prediction 4.77403e-05
## 3          selected-variable overlap 1.00000e+00
##                                        interpretation
## 1          includes any difference in selected lambda
## 2 isolates numerical implementation at a fixed lambda
## 3               Jaccard share of the two nonzero sets
```

``` r
data.frame(
  implementation = c("manual", "glmnet"),
  nonzero_predictors = c(length(manual_selected), length(glmnet_selected)),
  best_lambda = c(best_lambda, glmnet_best_lambda)
)
```

```
##   implementation nonzero_predictors best_lambda
## 1         manual                  7 0.003798291
## 2         glmnet                  7 0.003798291
```

歷史平均相對自己的 `OOS_R2` 分母與分子相同，故為 0。兩個 LASSO 版本的負
`OOS_R2` 都必須原樣報告，不可因結果不理想而移動測試起點或重新選字典。


``` r
plot(dates[test], actual, type = "l", lwd = 1.8, col = "black",
     xlab = "Predictor month", ylab = "Next-month excess return")
lines(dates[test], final$pred, col = "#A34045", lwd = 1.5)
lines(dates[test], glmnet_prediction, col = "#1D6D73", lwd = 1.3, lty = 3)
lines(dates[test], baseline, col = "#62717E", lty = 2)
legend("topleft", c("Actual", "LASSO (manual)", "LASSO (glmnet)",
                     "Historical mean"),
       col = c("black", "#A34045", "#1D6D73", "#62717E"),
       lty = c(1, 1, 3, 2), lwd = c(1.8, 1.5, 1.3, 1), bty = "n")
```

![固定 test period 的製造業下一期報酬預測。](../R14_financial_factor_selection_files/figure-gfm/prediction-plot-1.png)

## 非零係數與跨折選取頻率


``` r
selected_final <- data.frame(
  predictor = final$names,
  coefficient = final$beta
)
selected_final <- selected_final[abs(selected_final$coefficient) > 1e-8, ]
selected_final <- selected_final[order(abs(selected_final$coefficient), decreasing = TRUE), ]
head(selected_final, 20)
```

```
##                   predictor   coefficient
## 122   factor_q_eg:macro_ltr  0.0038608918
## 104 factor_q_roe:macro_svar -0.0023348662
## 95    factor_q_ia:macro_lty -0.0018041462
## 91   factor_q_ia:macro_svar -0.0013755408
## 19                macro_tms  0.0009067502
## 97    factor_q_ia:macro_tms -0.0006899683
## 28  factor_ff_rf:macro_ntis -0.0003856468
```

``` r
selection <- matrix(0L, nrow = length(folds), ncol = ncol(X),
                    dimnames = list(NULL, colnames(X)))
for (v in seq_along(folds)) {
  tr <- folds[[v]]$train
  pp <- standardize_train(X[tr, , drop = FALSE])
  b <- lasso_cd(pp$train, y[tr] - mean(y[tr]), best_lambda)
  selection[v, pp$keep] <- as.integer(abs(b) > 1e-8)
}
frequency <- sort(colMeans(selection), decreasing = TRUE)
head(data.frame(predictor = names(frequency), selection_frequency = frequency), 20)
```

```
##                                       predictor selection_frequency
## factor_ff_rf:macro_ntis factor_ff_rf:macro_ntis                 1.0
## factor_q_roe:macro_svar factor_q_roe:macro_svar                 1.0
## factor_q_eg:macro_ltr     factor_q_eg:macro_ltr                 1.0
## macro_dfy                             macro_dfy                 0.8
## factor_q_ia:macro_lty     factor_q_ia:macro_lty                 0.8
## macro_tbl                             macro_tbl                 0.6
## macro_infl                           macro_infl                 0.6
## factor_ff_rf:macro_ep     factor_ff_rf:macro_ep                 0.6
## factor_q_ia:macro_tms     factor_q_ia:macro_tms                 0.6
## factor_q_roe:macro_tms   factor_q_roe:macro_tms                 0.6
## factor_q_eg:macro_infl   factor_q_eg:macro_infl                 0.6
## macro_dp                               macro_dp                 0.4
## macro_tms                             macro_tms                 0.4
## factor_q_ia:macro_ep       factor_q_ia:macro_ep                 0.4
## macro_ltr                             macro_ltr                 0.2
## factor_ff_rf:macro_ltr   factor_ff_rf:macro_ltr                 0.2
## factor_ff_rf:macro_tms   factor_ff_rf:macro_tms                 0.2
## factor_ff_smb:macro_dfy factor_ff_smb:macro_dfy                 0.2
## factor_q_me:macro_dfy     factor_q_me:macro_dfy                 0.2
## factor_q_ia:macro_dy       factor_q_ia:macro_dy                 0.2
```

選取頻率是這組歷史折、字典與固定 `best_lambda` 下的描述，不是 \(p\) 值或「被定價機率」。高度相關的預測變數可能互相替代。

## 累積損失差與結構穩定性


``` r
loss_difference <- (actual - final$pred)^2 - (actual - baseline)^2
plot(dates[test], cumsum(loss_difference), type = "l", lwd = 2,
     col = "#1D6D73", xlab = "Predictor month",
     ylab = "Cumulative loss difference")
abline(h = 0, lty = 2, col = "#62717E")
```

![LASSO 相對歷史平均的累積平方損失差；下降表示 LASSO 在該段累積較佳。](../R14_financial_factor_selection_files/figure-gfm/cumulative-loss-1.png)

## 延伸到十產業的安全規則

1. 預先固定預測期距、測試月份、字典、折與 lambda 網格。
2. 每個產業可在自己的驗證折選 lambda；不得看測試結果後另改網格。
3. 報告全部十產業，包括負的樣本外表現。
4. 公開版直接使用固定 processed CSV 與 manifest；若另作上游重建，應保存建置說明與版本，不可把即時下載默認為同一份歷史快照。


``` r
sessionInfo()
```

```
## R version 4.5.2 (2025-10-31)
## Platform: aarch64-apple-darwin20
## Running under: macOS Tahoe 26.5.1
## 
## Matrix products: default
## BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
## 
## locale:
## [1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8
## 
## time zone: Asia/Tokyo
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] tibble_3.3.0 dplyr_1.2.1 
## 
## loaded via a namespace (and not attached):
##  [1] shape_1.4.6.1       gtable_0.3.6        xfun_0.57          
##  [4] ggplot2_4.0.3       collapse_2.1.7      lattice_0.22-7     
##  [7] quadprog_1.5-8      vctrs_0.7.2         tools_4.5.2        
## [10] Rdpack_2.6.6        generics_0.1.4      curl_7.0.0         
## [13] parallel_4.5.2      sandwich_3.1-1      xts_0.14.2         
## [16] pkgconfig_2.0.3     gbutils_0.5.1       Matrix_1.7-4       
## [19] tidyverse_2.0.0     RColorBrewer_1.1-3  S7_0.2.1           
## [22] lifecycle_1.0.5     compiler_4.5.2      farver_2.1.2       
## [25] MatrixModels_0.5-4  maxLik_1.5-2.2      textshaping_1.0.5  
## [28] codetools_0.2-20    SparseM_1.84-2      quantreg_6.1       
## [31] htmltools_0.5.9     glmnet_4.1-10       Formula_1.2-5      
## [34] pillar_1.11.1       MASS_7.3-65         plm_2.6-7          
## [37] iterators_1.0.14    foreach_1.5.2       nlme_3.1-168       
## [40] fracdiff_1.5-4      pls_2.9-0           fBasics_4052.98    
## [43] tidyselect_1.2.1    bdsmatrix_1.3-7     digest_0.6.39      
## [46] labeling_0.4.3      splines_4.5.2       tseries_0.10-62    
## [49] miscTools_0.6-30    fastmap_1.2.0       grid_4.5.2         
## [52] colorspace_2.1-2    cli_3.6.5           magrittr_2.0.4     
## [55] utf8_1.2.6          survival_3.8-3      withr_3.0.2        
## [58] scales_1.4.0        forecast_9.0.2      TTR_0.24.4         
## [61] rmarkdown_2.31      quantmod_0.4.29     otel_0.2.0         
## [64] timeDate_4052.112   ragg_1.5.2          zoo_1.8-15         
## [67] timeSeries_4052.112 fGarch_4052.93      urca_1.3-4         
## [70] evaluate_1.0.5      knitr_1.51          rbibutils_2.4.1    
## [73] lmtest_0.9-40       rlang_1.1.7         spatial_7.3-18     
## [76] Rcpp_1.1.0          glue_1.8.0          R6_2.6.1           
## [79] cvar_0.6            systemfonts_1.3.2
```
