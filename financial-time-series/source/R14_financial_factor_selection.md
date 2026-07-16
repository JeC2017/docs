---
title: "R14：機器學習的金融因子選擇"
output:
  github_document:
    toc: true
    toc_depth: 3
---

本附錄對應第 19 章。它使用專案固定的十產業因子面板，以第 \(t\) 月的預測變數預測製造業第 \(t+1\) 月超額報酬。重點是時間對齊、只依訓練資料進行前處理、展開視窗驗證、基準比較與選取穩定性。

## 執行條件

- 只使用 base R 與 `knitr`。
- 不安裝套件、不呼叫網路、不使用 `setwd()`。
- 固定資料：`data/processed/ff_qf_macro_industries_1967_2021.csv`。
- 資料建置與再散布注意事項：`data/DATA_SOURCES.md`。


``` r
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.width = 7, fig.height = 4.5
)
```


``` r
path <- "data/processed/ff_qf_macro_industries_1967_2021.csv"
stopifnot(file.exists(path))
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


``` r
plot(log(lambda_grid), mean_loss, type = "l", lwd = 2,
     xlab = "log(lambda)", ylab = "Mean validation MSE",
     col = "#173B57")
abline(v = log(best_lambda), lty = 2, col = "#A34045")
```

![製造業 expanding-window validation loss。](./R14_financial_factor_selection_files/figure-gfm/cv-plot-1.png)

## 最終測試評量


``` r
final <- fit_predict_lasso(
  X[tv, , drop = FALSE], y[tv],
  X[test, , drop = FALSE], best_lambda
)

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
  LASSO = score(actual, final$pred, baseline)
)
```

```
##                        MSE        MAE      OOS_R2
## HistoricalMean 0.002112639 0.03413601  0.00000000
## LASSO          0.002136697 0.03393833 -0.01138775
```

歷史平均相對自己的 `OOS_R2` 分母與分子相同，故為 0。負的 LASSO `OOS_R2` 必須原樣報告，不可因結果不理想而移動測試起點或重新選字典。


``` r
plot(dates[test], actual, type = "l", lwd = 1.8, col = "black",
     xlab = "Predictor month", ylab = "Next-month excess return")
lines(dates[test], final$pred, col = "#A34045", lwd = 1.5)
lines(dates[test], baseline, col = "#62717E", lty = 2)
legend("topleft", c("Actual", "LASSO", "Historical mean"),
       col = c("black", "#A34045", "#62717E"),
       lty = c(1, 1, 2), lwd = c(1.8, 1.5, 1), bty = "n")
```

![固定 test period 的製造業下一期報酬預測。](./R14_financial_factor_selection_files/figure-gfm/prediction-plot-1.png)

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

![LASSO 相對歷史平均的累積平方損失差；下降表示 LASSO 在該段累積較佳。](./R14_financial_factor_selection_files/figure-gfm/cumulative-loss-1.png)

## 延伸到十產業的安全規則

1. 預先固定預測期距、測試月份、字典、折與 lambda 網格。
2. 每個產業可在自己的驗證折選 lambda；不得看測試結果後另改網格。
3. 報告全部十產業，包括負的樣本外表現。
4. 若資料供應條款在公開版要求讀者自行建置，應保留 manifest 與建置說明，不可改成執行時即時下載。


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
##  [1] vctrs_0.7.2        cli_3.6.5          knitr_1.51         rlang_1.1.7       
##  [5] xfun_0.57          otel_0.2.0         MatrixModels_0.5-4 generics_0.1.4    
##  [9] glue_1.8.0         grid_4.5.2         evaluate_1.0.5     SparseM_1.84-2    
## [13] MASS_7.3-65        lifecycle_1.0.5    compiler_4.5.2     pkgconfig_2.0.3   
## [17] quantreg_6.1       lattice_0.22-7     R6_2.6.1           tidyselect_1.2.1  
## [21] utf8_1.2.6         splines_4.5.2      pillar_1.11.1      magrittr_2.0.4    
## [25] Matrix_1.7-4       tools_4.5.2        withr_3.0.2        survival_3.8-3
```
