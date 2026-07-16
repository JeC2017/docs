---
title: "R13：Ridge、LASSO、Elastic Net 與 post-LASSO"
output:
  github_document:
    toc: true
    toc_depth: 3
---

本附錄對應第 17--18 章。目標是以 2007 年 10 月至 2018 年 10 月、共 133 個月的真實日本總體金融凍結快照，預測下一期日本股市報酬，並比較 OLS、Ridge、LASSO、Elastic Net 與 post-LASSO。資料來自原課程的 `data_t.csv` 與 `yield_10.csv`；原變數說明記載股價、匯率、工業生產、利率、美股、失業、CPI、M3、WTI、人口結構、外資、貿易與十年期殖利率的單位及季調狀態，但沒有保存原供應者、URL、下載日或 vintage。因此本附錄沿用來源欄位尺度，不把它改稱為可從公開來源完整重建的資料。這是固定時間切分下的預測比較，不識別任何預測變數的因果效果。

## 執行環境與資料

- R 4.1 以上；模型計算只使用 base R，`knitr` 負責轉檔，`ragg` 與 `systemfonts` 負責以 cwTeX 字型產生圖檔。
- 資料：`data/processed/japan_monthly_2007_2018.csv`。
- 建置紀錄與 MD5：`data/processed/manifest.csv`。
- 可從教科書專案根目錄或 `online_appendix/` knit；資料位置由相對路徑搜尋函數判定。
- 公開界線：repo 隨附作者授權的 processed CSV、程式與執行結果，可離線自含重跑；provider、URL、下載日與 vintage 的缺口仍限制從上游原始來源重建。


``` r
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.width = 8,
  fig.height = 5,
  dev = "ragg_png", dpi = 144,
  dev.args = list(background = "white")
)
set.seed(20260716)

root_candidates <- c(".", "..")
is_root <- vapply(root_candidates, function(x) {
  file.exists(file.path(x, "main.tex"))
}, logical(1))
stopifnot(any(is_root))
project_root <- root_candidates[which(is_root)[1]]
project_path <- function(...) file.path(project_root, ...)

stopifnot(
  requireNamespace("ragg", quietly = TRUE),
  requireNamespace("systemfonts", quietly = TRUE)
)
cwtex_file <- project_path("assets", "fonts", "cwTeXQKai-Medium.ttf")
stopifnot(file.exists(cwtex_file))
if (!"cwTeX Online" %in% systemfonts::registry_fonts()$family) {
  systemfonts::register_font("cwTeX Online", cwtex_file)
}
plot_family <- "cwTeX Online"
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

data_path <- locate_project_file(
  "data/processed/japan_monthly_2007_2018.csv"
)
manifest_path <- locate_project_file("data/processed/manifest.csv")

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
jp <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)
jp$date <- as.Date(jp$date)
jp <- jp[order(jp$date), ]

stopifnot(nrow(jp) == 133L, ncol(jp) == 30L)
manifest[grepl("japan_monthly", manifest$file), ]
```

```
##                                         file rows columns
## 8 data/processed/japan_monthly_2007_2018.csv  133      30
##                                md5
## 8 46b39f6fdde5d581ad31c83348d99933
##                                                          description
## 8 Japanese monthly macro-finance panel with 10-year government yield
##                  built_at
## 8 2026-07-16 09:57:40 UTC
```

``` r
data.frame(
  first_month = min(jp$date),
  last_month = max(jp$date),
  months = nrow(jp),
  observation_unit = "monthly Japan macro-finance observation",
  return_unit = "source-file scale; not independently relabelled"
)
```

```
##   first_month last_month months                        observation_unit
## 1  2007-10-01 2018-10-01    133 monthly Japan macro-finance observation
##                                       return_unit
## 1 source-file scale; not independently relabelled
```

## 建立可實現的下一期 target

第 \(t\) 月的預測變數只能預測第 \(t+1\) 月的 `return_j`。日期、當期目標與下一期才會知道的值不進入預測變數矩陣。缺值列在完成時間對齊後一次處理。


``` r
jp$target_next <- c(jp$return_j[-1], NA_real_)
predictor_names <- setdiff(
  names(jp),
  c("date", "return_j", "target_next")
)

model_df <- jp[, c("date", "target_next", predictor_names)]
model_df <- model_df[complete.cases(model_df), ]
X_raw <- as.matrix(model_df[, predictor_names])
storage.mode(X_raw) <- "double"
y <- model_df$target_next
dates <- model_df$date

c(observations = length(y), predictors = ncol(X_raw))
```

```
## observations   predictors 
##          131           28
```

``` r
head(data.frame(predictor_month = dates, target_next = y), 3)
```

```
##   predictor_month target_next
## 1      2007-11-01   0.1385835
## 2      2007-12-01  -4.5285968
## 3      2008-01-01  -1.0026484
```

## 只由訓練資料估計前處理


``` r
prep_x <- function(X_train, X_new = NULL) {
  mu <- colMeans(X_train)
  s <- apply(X_train, 2, sd)
  keep <- is.finite(s) & s > 1e-10
  X_train_s <- sweep(sweep(X_train[, keep, drop = FALSE], 2, mu[keep]),
                       2, s[keep], "/")
  ans <- list(
    train = X_train_s,
    mean = mu[keep],
    sd = s[keep],
    keep = keep
  )
  if (!is.null(X_new)) {
    ans$new <- sweep(sweep(X_new[, keep, drop = FALSE], 2, mu[keep]),
                     2, s[keep], "/")
  }
  ans
}

soft_threshold <- function(z, penalty) {
  sign(z) * pmax(abs(z) - penalty, 0)
}
```

## 三種正則化估計器

以下目標函數均為

\[
\frac{1}{2n}\lVert y-Xb\rVert_2^2+
\lambda\left\{\alpha\lVert b\rVert_1+
\frac{1-\alpha}{2}\lVert b\rVert_2^2\right\}.
\]


``` r
fit_ridge <- function(X, y_centered, lambda) {
  p <- ncol(X)
  solve(crossprod(X) / nrow(X) + lambda * diag(p),
        crossprod(X, y_centered) / nrow(X))[, 1]
}

fit_enet_cd <- function(X, y_centered, lambda, alpha = 1,
                        max_iter = 5000L, tol = 1e-8) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- numeric(p)
  residual <- y_centered
  x2 <- colSums(X^2) / n

  for (iter in seq_len(max_iter)) {
    beta_old <- beta
    for (j in seq_len(p)) {
      residual <- residual + X[, j] * beta[j]
      z <- sum(X[, j] * residual) / n
      beta[j] <- soft_threshold(z, lambda * alpha) /
        (x2[j] + lambda * (1 - alpha))
      residual <- residual - X[, j] * beta[j]
    }
    if (max(abs(beta - beta_old)) < tol) break
  }
  attr(beta, "iterations") <- iter
  beta
}

predict_scaled <- function(X_train, y_train, X_new,
                           lambda, method, alpha = 1) {
  pp <- prep_x(X_train, X_new)
  y_bar <- mean(y_train)
  yc <- y_train - y_bar
  beta <- switch(
    method,
    ridge = fit_ridge(pp$train, yc, lambda),
    enet = fit_enet_cd(pp$train, yc, lambda, alpha)
  )
  list(
    pred = as.numeric(y_bar + pp$new %*% beta),
    beta = beta,
    names = colnames(pp$train),
    y_bar = y_bar,
    prep = pp
  )
}
```

## 依時間排序的驗證折

最末 25% 留作最終測試集。前 75% 內使用擴展視窗；沒有任何月份被隨機洗牌。


``` r
n <- length(y)
test_start <- floor(0.75 * n) + 1L
idx_tv <- seq_len(test_start - 1L)
idx_test <- test_start:n

validation_size <- 8L
first_train <- max(45L, floor(0.55 * length(idx_tv)))
train_ends <- unique(as.integer(seq(
  first_train,
  length(idx_tv) - validation_size,
  length.out = 4
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
## 1 2012-03-01       2012-04-01     2012-11-01
## 2 2013-03-01       2013-04-01     2013-11-01
## 3 2014-03-01       2014-04-01     2014-11-01
## 4 2015-04-01       2015-05-01     2015-12-01
```

## 調校 \(\lambda\) 與 \(\alpha\)

LASSO 與 Elastic Net 的網格分別從「所有斜率恰為零」的門檻
`lambda_max / alpha` 向下延伸三個數量級；Ridge 的網格則涵蓋八個數量級。
每一折都重新估計中心與尺度。把零係數門檻明列為端點，可分辨資料真的偏好
只有截距的預測，或只是調校網格沒有包住最小值。


``` r
lambda_max_by_fold <- vapply(folds, function(fold) {
  pp <- prep_x(X_raw[fold$train, , drop = FALSE])
  yy <- y[fold$train]
  max(abs(crossprod(pp$train, yy - mean(yy)))) / length(yy)
}, numeric(1))
lambda_max <- max(lambda_max_by_fold)
lambda_lasso <- exp(seq(
  log(lambda_max), log(lambda_max * 0.001), length.out = 40
))
elastic_alpha <- 0.5
lambda_elastic <- exp(seq(
  log(lambda_max / elastic_alpha),
  log(lambda_max / elastic_alpha * 0.001),
  length.out = 40
))
lambda_ridge <- exp(seq(log(1e-4), log(1e4), length.out = 49))

cv_loss <- function(lambda_grid, method, alpha = 1) {
  out <- matrix(NA_real_, nrow = length(folds), ncol = length(lambda_grid))
  for (v in seq_along(folds)) {
    tr <- folds[[v]]$train
    va <- folds[[v]]$validation
    for (j in seq_along(lambda_grid)) {
      fit <- predict_scaled(
        X_raw[tr, , drop = FALSE], y[tr],
        X_raw[va, , drop = FALSE],
        lambda = lambda_grid[j], method = method, alpha = alpha
      )
      out[v, j] <- mean((y[va] - fit$pred)^2)
    }
  }
  colMeans(out)
}

loss_ridge <- cv_loss(lambda_ridge, "ridge")
loss_lasso <- cv_loss(lambda_lasso, "enet", alpha = 1)
loss_enet <- cv_loss(lambda_elastic, "enet", alpha = elastic_alpha)

best_ridge <- lambda_ridge[which.min(loss_ridge)]
best_lasso <- lambda_lasso[which.min(loss_lasso)]
best_enet <- lambda_elastic[which.min(loss_enet)]

data.frame(
  model = c("Ridge", "LASSO", "Elastic Net"),
  lambda = c(best_ridge, best_lasso, best_enet),
  validation_mse = c(min(loss_ridge), min(loss_lasso), min(loss_enet)),
  grid_location = c(
    "interior",
    if (which.min(loss_lasso) == 1L) "zero-slope threshold" else "interior",
    if (which.min(loss_enet) == 1L) "zero-slope threshold" else "interior"
  )
)
```

```
##         model    lambda validation_mse        grid_location
## 1       Ridge 31.622777       4.922462             interior
## 2       LASSO  1.138271       5.025106 zero-slope threshold
## 3 Elastic Net  2.276542       5.025106 zero-slope threshold
```

``` r
stopifnot(which.min(loss_ridge) %in% 2:(length(lambda_ridge) - 1L))
```

Ridge 的最小值落在擴大網格內部。LASSO 與 Elastic Net 則選到各自的
`zero-slope threshold`：這個上端點是依所有驗證折中最大的零係數門檻建立，
再增大 \(\lambda\) 也只會維持全部斜率為零。因此這兩個邊界解表示驗證資料
偏好只有截距的預測，不是網格仍未向上延伸的未解問題。


``` r
old_par <- par(family = plot_family)
plot(log(lambda_ridge), loss_ridge, type = "l", lwd = 2,
     xlab = expression(log(lambda)), ylab = "驗證期 MSE", col = "#173B57")
lines(log(lambda_lasso), loss_lasso, lwd = 2, col = "#A34045")
lines(log(lambda_elastic), loss_enet, lwd = 2, col = "#1D6D73")
legend("topleft", c("Ridge", "LASSO", "Elastic Net"),
       col = c("#173B57", "#A34045", "#1D6D73"), lwd = 2, bty = "n")
```

![擴展視窗驗證期的均方預測誤差。](../R13_regularization_files/figure-gfm/validation-plot-1.png)

``` r
par(old_par)
```

## 在完全未見的測試期間比較模型


``` r
X_tv <- X_raw[idx_tv, , drop = FALSE]
y_tv <- y[idx_tv]
X_test <- X_raw[idx_test, , drop = FALSE]
y_test <- y[idx_test]

ridge <- predict_scaled(X_tv, y_tv, X_test, best_ridge, "ridge")
lasso <- predict_scaled(X_tv, y_tv, X_test, best_lasso, "enet", alpha = 1)
enet <- predict_scaled(
  X_tv, y_tv, X_test, best_enet, "enet", alpha = elastic_alpha
)

# 以 QR pivot 選出線性獨立欄，再估計一個明確、可重現的降秩 OLS。
# 在秩虧設計下，完整係數向量不唯一；此基準不作結構性係數解讀。
rank_reduced_ols_predict <- function(X_train, y_train, X_new,
                                     tolerance = 1e-8) {
  X_train <- cbind(Intercept = 1, as.matrix(X_train))
  X_new <- cbind(Intercept = 1, as.matrix(X_new))
  qrx <- qr(X_train, tol = tolerance, LAPACK = FALSE)
  keep <- sort(qrx$pivot[seq_len(qrx$rank)])
  fit <- lm.fit(X_train[, keep, drop = FALSE], y_train)
  prediction <- as.numeric(
    X_new[, keep, drop = FALSE] %*% fit$coefficients
  )
  stopifnot(all(is.finite(prediction)))
  list(
    pred = prediction,
    rank = qrx$rank,
    columns = ncol(X_train),
    kept = colnames(X_train)[keep],
    dropped = setdiff(colnames(X_train), colnames(X_train)[keep])
  )
}

ols <- rank_reduced_ols_predict(X_tv, y_tv, X_test)
ols_pred <- ols$pred
data.frame(
  benchmark = "QR rank-reduced OLS",
  design_columns = ols$columns,
  numerical_rank = ols$rank,
  dropped_columns = length(ols$dropped)
)
```

```
##             benchmark design_columns numerical_rank dropped_columns
## 1 QR rank-reduced OLS             29             28               1
```

``` r
# post-LASSO：先選，再在相同 training+validation 上重新 OLS。
selected <- which(abs(lasso$beta) > 1e-8)
if (length(selected) == 0L) {
  post_pred <- rep(mean(y_tv), length(y_test))
  post_rank <- 1L
} else {
  pp_final <- lasso$prep
  post_fit <- rank_reduced_ols_predict(
    pp_final$train[, selected, drop = FALSE],
    y_tv,
    pp_final$new[, selected, drop = FALSE]
  )
  post_pred <- post_fit$pred
  post_rank <- post_fit$rank
}

pred <- data.frame(
  date = dates[idx_test],
  actual = y_test,
  HistoricalMean = rep(mean(y_tv), length(y_test)),
  OLS_QR = ols_pred,
  Ridge = ridge$pred,
  LASSO = lasso$pred,
  ElasticNet = enet$pred,
  PostLASSO = post_pred
)

stopifnot(all(vapply(
  pred[setdiff(names(pred), "date")],
  function(x) all(is.finite(x)),
  logical(1)
)))
```


``` r
score <- function(actual, forecast, baseline) {
  mse <- mean((actual - forecast)^2)
  c(
    MSE = mse,
    MAE = mean(abs(actual - forecast)),
    OOS_R2 = 1 - mse / mean((actual - baseline)^2)
  )
}

model_names <- setdiff(names(pred), c("date", "actual", "HistoricalMean"))
metrics <- t(vapply(
  model_names,
  function(nm) score(pred$actual, pred[[nm]], pred$HistoricalMean),
  numeric(3)
))
stopifnot(all(is.finite(metrics)))
round(metrics, 4)
```

```
##                MSE    MAE  OOS_R2
## OLS_QR     14.4723 3.3798 -6.8969
## Ridge       1.7429 1.0231  0.0490
## LASSO       1.8327 1.0738  0.0000
## ElasticNet  1.8327 1.0738  0.0000
## PostLASSO   1.8327 1.0738  0.0000
```


``` r
coef_table <- data.frame(
  predictor = lasso$names,
  lasso = lasso$beta,
  elastic_net = enet$beta
)
coef_table <- coef_table[order(abs(coef_table$lasso), decreasing = TRUE), ]
head(coef_table, 12)
```

```
##     predictor lasso elastic_net
## 1         spj     0           0
## 2         rer     0           0
## 3  rer_change     0           0
## 4         ipi     0           0
## 5  ipi_change     0           0
## 6         inr     0           0
## 7  inr_change     0           0
## 8         spf     0           0
## 9    return_f     0           0
## 10        unr     0           0
## 11 unr_change     0           0
## 12        cpi     0           0
```

``` r
cat("LASSO nonzero predictors:", sum(abs(lasso$beta) > 1e-8), "\n")
```

```
## LASSO nonzero predictors: 0
```


``` r
old_par <- par(family = plot_family)
matplot(
  pred$date,
  pred[, c("actual", "Ridge", "LASSO", "ElasticNet")],
  type = "l", lty = c(1, 2, 3, 4), lwd = c(2, 1.5, 1.5, 1.5),
  col = c("black", "#173B57", "#A34045", "#1D6D73"),
  xlab = "預測變數月份", ylab = "次月報酬率"
)
legend(
  "topleft", c("實現值", "Ridge", "LASSO", "Elastic Net"),
  lty = c(1, 2, 3, 4), lwd = c(2, 1.5, 1.5, 1.5),
  col = c("black", "#173B57", "#A34045", "#1D6D73"), bty = "n"
)
```

![最終 test period：實現值與正則化預測。](../R13_regularization_files/figure-gfm/test-plot-1.png)

``` r
par(old_par)
```

## 解讀限制

1. 這是小樣本預測示範，不是日本股市風險溢酬的結構估計。
2. `return_j` 的單位沿用凍結檔；不可與未核對單位的外部報酬直接合併。
3. 非零 LASSO 係數表示在這個字典、切割與損失下被使用，不等於顯著或因果。
4. 測試集沒有參與調校；若看完測試結果後改模型，必須另留新的測試期。


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
## loaded via a namespace (and not attached):
##  [1] compiler_4.5.2    cli_3.6.5         ragg_1.5.2        tools_4.5.2      
##  [5] otel_0.2.0        knitr_1.51        xfun_0.57         textshaping_1.0.5
##  [9] lifecycle_1.0.5   systemfonts_1.3.2 rlang_1.1.7       evaluate_1.0.5
```
