---
title: "R09：投資組合、CAPM 與財務因子"
output:
  github_document:
    toc: true
    toc_depth: 3
---

本附錄對應第 13--14 章，以固定種子的合成月資料示範投資組合、CAPM 與多因子迴歸。所有報酬都是教學模擬，不代表任何市場或基金；程式不安裝套件、不下載資料，也不使用 `setwd()`。


``` r
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.width = 7, fig.height = 4
)
stopifnot(getRversion() >= "4.3.0")
set.seed(909)
```

## 1. 軟體與資訊時點

只使用 R 內建套件。時間切分固定為：

- 前 480 月：訓練；
- 接續 120 月：驗證；
- 最後 120 月：測試。

任何平均數、共變異數、beta 與投資組合權重都只能由訓練期估計。驗證與測試期的同期因子若用於重建同期資產報酬，會明確稱為「條件重建」，不稱為事前預測。


``` r
data.frame(
  component = c("R", "stats", "graphics"),
  version = c(
    R.version.string,
    as.character(packageVersion("stats")),
    as.character(packageVersion("graphics"))
  )
)
```

```
##   component                      version
## 1         R R version 4.5.2 (2025-10-31)
## 2     stats                        4.5.2
## 3  graphics                        4.5.2
```

## 2. 建構方向正確的 SMB 與 HML

以六個規模 \(S/B\) 與帳面市值比 \(L/M/H\) 投資組合建立因子：

\[
\begin{aligned}
SMB&=(S_L+S_M+S_H)/3-(B_L+B_M+B_H)/3,\\
HML&=(S_H+B_H)/2-(S_L+B_L)/2.
\end{aligned}
\]


``` r
n_total <- 720
portfolio_names <- c("SL", "SM", "SH", "BL", "BM", "BH")

# 共同擾動讓六個基礎投資組合有合理的正相關。
common <- arima.sim(
  model = list(ar = 0.10),
  n = n_total, sd = 0.025
)
style_noise <- matrix(
  rnorm(n_total * length(portfolio_names), sd = 0.018),
  nrow = n_total,
  dimnames = list(NULL, portfolio_names)
)
base_portfolios <- sweep(style_noise, 1, common, "+")

# 加入小型平均規模與價值溢酬；數值純屬模擬設計。
base_portfolios[, c("SL", "SM", "SH")] <-
  base_portfolios[, c("SL", "SM", "SH")] + 0.0015
base_portfolios[, c("SH", "BH")] <-
  base_portfolios[, c("SH", "BH")] + 0.0020

SMB <- rowMeans(base_portfolios[, c("SL", "SM", "SH")]) -
  rowMeans(base_portfolios[, c("BL", "BM", "BH")])
HML <- rowMeans(base_portfolios[, c("SH", "BH")]) -
  rowMeans(base_portfolios[, c("SL", "BL")])

# 市場因子另行模擬；MKT 已是市場超額報酬。
MKT <- as.numeric(arima.sim(
  model = list(ar = 0.12),
  n = n_total, sd = 0.035
)) + 0.004

factors <- data.frame(
  month = seq_len(n_total),
  MKT = MKT,
  SMB = SMB,
  HML = HML
)
head(factors)
```

```
##   month          MKT          SMB          HML
## 1     1 -0.009539671 -0.007736603 -0.007997446
## 2     2 -0.029083750 -0.002171098 -0.036089712
## 3     3 -0.045304436  0.014443442  0.004514239
## 4     4 -0.003159349  0.001501598  0.013537953
## 5     5 -0.042868669  0.015459097 -0.017965287
## 6     6  0.007852383  0.018886008  0.014791961
```

### 2.1 HML 符號單元測試

若所有高帳面市值比投資組合各增加一個單位，HML 必須增加一個單位；若低比率投資組合各增加一個單位，HML 必須下降一個單位。


``` r
make_hml <- function(x) {
  rowMeans(x[, c("SH", "BH"), drop = FALSE]) -
    rowMeans(x[, c("SL", "BL"), drop = FALSE])
}

toy <- base_portfolios[1:5, , drop = FALSE]
high_up <- toy
high_up[, c("SH", "BH")] <- high_up[, c("SH", "BH")] + 1
low_up <- toy
low_up[, c("SL", "BL")] <- low_up[, c("SL", "BL")] + 1

stopifnot(
  isTRUE(all.equal(make_hml(high_up) - make_hml(toy), rep(1, 5))),
  isTRUE(all.equal(make_hml(low_up) - make_hml(toy), rep(-1, 5)))
)

data.frame(
  change = c("高 B/M 組合 +1", "低 B/M 組合 +1"),
  HML_change = c(
    mean(make_hml(high_up) - make_hml(toy)),
    mean(make_hml(low_up) - make_hml(toy))
  )
)
```

```
##           change HML_change
## 1 高 B/M 組合 +1          1
## 2 低 B/M 組合 +1         -1
```

## 3. 模擬六個資產報酬

令資產超額報酬由三因子曝險加上個別新衝擊產生。


``` r
beta_true <- rbind(
  Asset1 = c(MKT = 1.10, SMB = 0.40, HML = -0.20),
  Asset2 = c(MKT = 0.75, SMB = -0.15, HML = 0.60),
  Asset3 = c(MKT = 1.30, SMB = 0.70, HML = 0.20),
  Asset4 = c(MKT = 0.55, SMB = -0.30, HML = -0.40),
  Asset5 = c(MKT = 1.00, SMB = 0.10, HML = 0.80),
  Asset6 = c(MKT = 0.90, SMB = 0.55, HML = -0.55)
)
alpha_true <- c(0.0000, 0.0005, 0.0000, -0.0003, 0.0000, 0.0002)

factor_matrix <- as.matrix(factors[, c("MKT", "SMB", "HML")])
idiosyncratic <- matrix(
  rnorm(n_total * nrow(beta_true), sd = 0.025),
  nrow = n_total
)
excess_return <- factor_matrix %*% t(beta_true) +
  matrix(alpha_true, nrow = n_total, ncol = nrow(beta_true),
         byrow = TRUE) +
  idiosyncratic
colnames(excess_return) <- rownames(beta_true)

rf_monthly <- 0.001
total_return <- excess_return + rf_monthly

train_id <- 1:480
valid_id <- 481:600
test_id <- 601:720
stopifnot(max(train_id) < min(valid_id), max(valid_id) < min(test_id))
```

## 4. 投資組合代數

### 4.1 由訓練期估計平均數與共變異數


``` r
mu_train <- colMeans(total_return[train_id, , drop = FALSE])
Sigma_train <- cov(total_return[train_id, , drop = FALSE])

portfolio_moments <- function(w, mu, Sigma) {
  stopifnot(length(w) == length(mu), abs(sum(w) - 1) < 1e-8)
  c(
    expected_return = drop(crossprod(w, mu)),
    variance = drop(t(w) %*% Sigma %*% w),
    standard_deviation = sqrt(drop(t(w) %*% Sigma %*% w))
  )
}
```

### 4.2 等權重、GMV 與穩定化 GMV


``` r
n_asset <- ncol(total_return)
w_equal <- rep(1 / n_asset, n_asset)

gmv_weights <- function(Sigma) {
  one <- rep(1, nrow(Sigma))
  raw <- solve(Sigma, one)
  drop(raw / sum(raw))
}

w_gmv <- gmv_weights(Sigma_train)

# 向對角矩陣收縮，降低小特徵值造成的極端權重。
shrinkage <- 0.25
Sigma_shrunk <- (1 - shrinkage) * Sigma_train +
  shrinkage * diag(diag(Sigma_train))
w_gmv_shrunk <- gmv_weights(Sigma_shrunk)

weights <- rbind(
  Equal = w_equal,
  GMV = w_gmv,
  Shrunk_GMV = w_gmv_shrunk
)
colnames(weights) <- colnames(total_return)
round(weights, 4)
```

```
##             Asset1 Asset2  Asset3 Asset4 Asset5 Asset6
## Equal       0.1667 0.1667  0.1667 0.1667 0.1667 0.1667
## GMV        -0.0654 0.3320 -0.1576 0.6078 0.0745 0.2086
## Shrunk_GMV -0.0013 0.2767 -0.0543 0.5405 0.0719 0.1665
```


``` r
t(apply(weights, 1, portfolio_moments,
        mu = mu_train, Sigma = Sigma_train))
```

```
##            expected_return     variance standard_deviation
## Equal          0.006693737 0.0012141421         0.03484454
## GMV            0.003314328 0.0007166275         0.02676990
## Shrunk_GMV     0.003896037 0.0007451363         0.02729718
```

### 4.3 真正向前的測試期績效

權重在進入測試期前即鎖定。以下績效只使用測試期實現報酬；這是可行的靜態保留期評估。


``` r
portfolio_path <- function(return_matrix, w) {
  drop(return_matrix %*% w)
}

max_drawdown <- function(r) {
  wealth <- cumprod(1 + r)
  drawdown <- wealth / cummax(c(1, wealth))[-1] - 1
  min(drawdown)
}

performance <- t(apply(weights, 1, function(w) {
  r <- portfolio_path(total_return[test_id, , drop = FALSE], w)
  c(
    mean_monthly = mean(r),
    sd_monthly = sd(r),
    sharpe_monthly = mean(r - rf_monthly) / sd(r - rf_monthly),
    max_drawdown = max_drawdown(r)
  )
}))
performance
```

```
##            mean_monthly sd_monthly sharpe_monthly max_drawdown
## Equal       0.009111001 0.03694560      0.2195390   -0.1311641
## GMV         0.004536931 0.02942022      0.1202211   -0.1683296
## Shrunk_GMV  0.005601164 0.03019659      0.1523736   -0.1494338
```


``` r
wealth <- sapply(seq_len(nrow(weights)), function(j) {
  cumprod(1 + portfolio_path(
    total_return[test_id, , drop = FALSE],
    weights[j, ]
  ))
})
matplot(
  wealth, type = "l", lty = 1, lwd = 1.2,
  col = c("#173B57", "#A34045", "#3F7158"),
  xlab = "測試期月份", ylab = "累積財富",
  main = "只用訓練期權重的保留期績效"
)
legend(
  "topleft", legend = rownames(weights),
  col = c("#173B57", "#A34045", "#3F7158"),
  lty = 1, bty = "n"
)
```

![plot of chunk wealth-plot](./R09_portfolio_capm_factors_files/figure-gfm/wealth-plot-1.png)

模擬中的勝負不是一般結論；重點是權重估計與績效期間已分離。

## 5. CAPM 與三因子迴歸

### 5.1 只用訓練期估計


``` r
reg_data <- data.frame(
  y = excess_return[, "Asset1"],
  factors[, c("MKT", "SMB", "HML")]
)

fit_capm <- lm(y ~ MKT, data = reg_data, subset = train_id)
fit_three <- lm(y ~ MKT + SMB + HML,
                data = reg_data, subset = train_id)

rbind(
  truth = c(
    `(Intercept)` = alpha_true[1],
    MKT = beta_true["Asset1", "MKT"],
    SMB = beta_true["Asset1", "SMB"],
    HML = beta_true["Asset1", "HML"]
  ),
  CAPM = c(coef(fit_capm), SMB = NA, HML = NA),
  ThreeFactor = coef(fit_three)
)
```

```
##              (Intercept)      MKT      SMB        HML
## truth       0.0000000000 1.100000 0.400000 -0.2000000
## CAPM        0.0007645363 1.152039       NA         NA
## ThreeFactor 0.0007656593 1.155590 0.456667 -0.1836672
```

CAPM 遺漏 SMB 與 HML 時，市場 beta 與 alpha 可能吸收相關的共同變動。即使三因子模型接近資料生成式，有限樣本估計仍不會等於真值。

### 5.2 基礎 R 的 Newey--West 型 HAC 共變異數


``` r
newey_west_vcov <- function(model, lag = 6) {
  X <- model.matrix(model)
  e <- residuals(model)
  n <- nrow(X)
  stopifnot(lag >= 0, lag < n)

  Xe <- X * as.numeric(e)
  meat <- crossprod(Xe)
  if (lag > 0) {
    for (ell in seq_len(lag)) {
      weight <- 1 - ell / (lag + 1)
      Gamma <- crossprod(
        Xe[(ell + 1):n, , drop = FALSE],
        Xe[1:(n - ell), , drop = FALSE]
      )
      meat <- meat + weight * (Gamma + t(Gamma))
    }
  }
  bread <- solve(crossprod(X))
  bread %*% meat %*% bread
}

V_hac <- newey_west_vcov(fit_three, lag = 6)
data.frame(
  estimate = coef(fit_three),
  HAC_se = sqrt(diag(V_hac)),
  row.names = names(coef(fit_three))
)
```

```
##                  estimate      HAC_se
## (Intercept)  0.0007656593 0.001254279
## MKT          1.1555902572 0.035825092
## SMB          0.4566670111 0.077520367
## HML         -0.1836672223 0.066937074
```

這是教學實作；正式分析還要依資料頻率、重疊報酬與條件異質變異選擇落後期數及有限樣本修正。

## 6. HML 反向會發生什麼？

若錯把 HML 寫成 low minus high，整列因子乘以 \(-1\)。迴歸配適值可以完全相同，但 HML beta 反號，經濟解讀也必須反向。


``` r
reg_wrong <- transform(reg_data, HML_wrong = -HML)
fit_wrong <- lm(
  y ~ MKT + SMB + HML_wrong,
  data = reg_wrong, subset = train_id
)

c(
  correct_HML_beta = coef(fit_three)["HML"],
  reversed_HML_beta = coef(fit_wrong)["HML_wrong"],
  maximum_fitted_difference = max(abs(
    fitted(fit_three) - fitted(fit_wrong)
  ))
)
```

```
##        correct_HML_beta.HML reversed_HML_beta.HML_wrong 
##                  -0.1836672                   0.1836672 
##   maximum_fitted_difference 
##                   0.0000000
```

``` r
stopifnot(
  isTRUE(all.equal(
    unname(coef(fit_three)["HML"]),
    -unname(coef(fit_wrong)["HML_wrong"]),
    tolerance = 1e-10
  ))
)
```

## 7. 條件重建，不冒充事前預測

用訓練期係數與驗證／測試期「同期已實現」因子重建同期資產報酬。這可比較因子集合的條件解釋能力，但因子在月初尚未實現，所以不是可交易的月初預測。


``` r
reconstruction_error <- function(model, newdata, actual) {
  predicted <- predict(model, newdata = newdata)
  c(
    RMSE = sqrt(mean((actual - predicted)^2)),
    MAE = mean(abs(actual - predicted))
  )
}

reconstruction <- rbind(
  CAPM_validation = reconstruction_error(
    fit_capm, reg_data[valid_id, ],
    reg_data$y[valid_id]
  ),
  ThreeFactor_validation = reconstruction_error(
    fit_three, reg_data[valid_id, ],
    reg_data$y[valid_id]
  ),
  CAPM_test = reconstruction_error(
    fit_capm, reg_data[test_id, ],
    reg_data$y[test_id]
  ),
  ThreeFactor_test = reconstruction_error(
    fit_three, reg_data[test_id, ],
    reg_data$y[test_id]
  )
)
reconstruction
```

```
##                              RMSE        MAE
## CAPM_validation        0.02543935 0.02078985
## ThreeFactor_validation 0.02442712 0.01957688
## CAPM_test              0.02422435 0.01913302
## ThreeFactor_test       0.02442523 0.01974323
```

若要做真正事前預測，必須在每個預測起點另行預測 MKT、SMB 與 HML，或建立直接預測資產報酬的模型；該因子預測也只能使用當時資訊。

## 8. 可重現結論

本附錄建立了以下防線：

1. HML 由函數與單元測試固定為 high minus low；
2. 投資組合權重只由訓練期平均數與共變異數決定；
3. CAPM 與三因子 beta 只在訓練期估計；
4. 使用同期測試因子的結果只稱為條件重建；
5. 所有資料皆為明確標示的合成教學資料。


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
##  [1] utf8_1.2.6       R6_2.6.1         tidyselect_1.2.1 xfun_0.57       
##  [5] magrittr_2.0.4   glue_1.8.0       knitr_1.51       pkgconfig_2.0.3 
##  [9] generics_0.1.4   lifecycle_1.0.5  cli_3.6.5        vctrs_0.7.2     
## [13] withr_3.0.2      compiler_4.5.2   tools_4.5.2      evaluate_1.0.5  
## [17] pillar_1.11.1    otel_0.2.0       rlang_1.1.7
```
