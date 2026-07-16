---
title: "R17：共整合、Johansen 檢定與 VECM"
output:
  github_document:
    toc: true
    toc_depth: 3
---

本附錄對應第 22 章。第一部分以已知共整合向量的模擬核對利差、誤差修正與 Johansen 特徵值；第二部分用固定日本短長期利率作探索性示範。為避免偽造臨界值，本附錄計算共整合秩統計量，但不自行硬編套件表格或宣稱單一模型規格已證明共整合。

## 執行條件

- base R 與 `knitr`；不安裝套件、不下載資料。
- 固定資料：`data/processed/japan_monthly_2007_2018.csv`。
- Johansen 函數是教學版不受限截距（unrestricted constant）、無差分落後差分項的示範，不替代完整套件及其臨界值表。


``` r
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.width = 7, fig.height = 4.5
)
set.seed(20260716)
```

## 已知共整合秩為 1 的教學模擬

建立共同隨機趨勢 (q_t=q_{t-1}+\eta_t) 與定態 spread
(z_t=0.65z_{t-1}+e_t)，再令 (y_{2t}=q_t)、(y_{1t}=q_t+z_t)。因此真實共整合向量為 \(\beta=(1,-1)^\top\)。


``` r
Tn <- 500L
eta <- rnorm(Tn)
e <- rnorm(Tn, sd = 0.5)
q <- cumsum(eta)
z <- numeric(Tn)
for (t in 2:Tn) z[t] <- 0.65 * z[t - 1] + e[t]

Y <- cbind(y1 = q + z, y2 = q)
spread <- Y[, 1] - Y[, 2]

par(mfrow = c(2, 1))
matplot(Y, type = "l", lty = 1, col = c("#173B57", "#A34045"),
        ylab = "Level", xlab = "t")
legend("topleft", c("y1", "y2"), col = c("#173B57", "#A34045"),
       lty = 1, bty = "n")
plot(spread, type = "l", col = "#1D6D73", ylab = "y1 - y2", xlab = "t")
abline(h = mean(spread), lty = 2)
```

![plot of chunk simulate](./R17_cointegration_vecm_files/figure-gfm/simulate-1.png)

``` r
par(mfrow = c(1, 1))
```

## AR(1) persistence diagnostic

下列函數回報 (\Delta x_t=a+\rho x_{t-1}+\) lagged differences 的 OLS t-statistic。它不是完整 ADF 檢定器；臨界值依 deterministic terms 與是否為估計殘差而不同。


``` r
adf_regression <- function(x, diff_lags = 0L, include_trend = FALSE) {
  x <- as.numeric(x)
  dx <- diff(x)
  n <- length(dx)
  start <- diff_lags + 1L
  y <- dx[start:n]
  x_lag <- x[start:n]
  X <- cbind(Intercept = 1, LevelLag = x_lag)
  if (include_trend) X <- cbind(X, Trend = seq_along(y))
  if (diff_lags > 0L) {
    for (j in seq_len(diff_lags)) {
      X <- cbind(X, dx[(start - j):(n - j)])
      colnames(X)[ncol(X)] <- paste0("dLag", j)
    }
  }
  fit <- lm.fit(X, y)
  df <- length(y) - fit$rank
  sigma2 <- sum(fit$residuals^2) / df
  vcov <- sigma2 * chol2inv(qr.R(fit$qr)[seq_len(fit$rank),
                                      seq_len(fit$rank), drop = FALSE])
  # 此範例 X 滿秩；LevelLag 是第二欄。
  c(rho = fit$coefficients["LevelLag"],
    t_stat = fit$coefficients["LevelLag"] / sqrt(vcov[2, 2]))
}

rbind(
  y1_level = adf_regression(Y[, 1], diff_lags = 1),
  y2_level = adf_regression(Y[, 2], diff_lags = 1),
  spread = adf_regression(spread, diff_lags = 1)
)
```

```
##          rho.LevelLag t_stat.LevelLag
## y1_level -0.010508363       -1.690820
## y2_level -0.009143012       -1.614674
## spread   -0.345531963       -9.004964
```

解讀只比較相對 persistence：模擬 levels 應比 spread 更接近單根。不得把一般常態 1.96 當 ADF 臨界值。

## 兩步 ECM


``` r
# 第一步：y1 對 y2 的長期迴歸。
long_run <- lm(y1 ~ y2, data = as.data.frame(Y))
ec_error <- residuals(long_run)

# 第二步：兩個差分方程都含前一期估計均衡誤差。
dY <- diff(Y)
ect_lag <- ec_error[-length(ec_error)]
ecm_y1 <- lm(dY[, 1] ~ ect_lag)
ecm_y2 <- lm(dY[, 2] ~ ect_lag)

data.frame(
  equation = c("Delta y1", "Delta y2"),
  adjustment = c(coef(ecm_y1)["ect_lag"], coef(ecm_y2)["ect_lag"]),
  standard_error = c(summary(ecm_y1)$coef["ect_lag", "Std. Error"],
                     summary(ecm_y2)$coef["ect_lag", "Std. Error"])
)
```

```
##   equation  adjustment standard_error
## 1 Delta y1 -0.28585029     0.07841242
## 2 Delta y2  0.07958777     0.07250540
```

``` r
coef(long_run)
```

```
## (Intercept)          y2 
##   0.1055747   0.9904325
```

因 (y_1-y_2=z_t)，正 spread 應由兩個變數的相對變動縮小。正規化方向改變時，\(\alpha\) 符號也會一起改變。

## 教學版 Johansen eigenvalues

對 VECM(1)

\[
\Delta Y_t=\Pi Y_{t-1}+c+u_t,
\]

分別從 \(\Delta Y_t\) 與 \(Y_{t-1}\) 移除截距，再計算 generalized eigenvalues。完整 Johansen 程序需按差分 lag 與 deterministic case 改變 residualization 與臨界值。


``` r
johansen_eigen_demo <- function(Y) {
  Y <- as.matrix(Y)
  dY <- diff(Y)
  Ylag <- Y[-nrow(Y), , drop = FALSE]
  n <- nrow(dY)

  # 不受限截距的教學簡化：兩者都對常數殘差化（移除欄平均）。
  R0 <- sweep(dY, 2, colMeans(dY))
  R1 <- sweep(Ylag, 2, colMeans(Ylag))
  S00 <- crossprod(R0) / n
  S11 <- crossprod(R1) / n
  S01 <- crossprod(R0, R1) / n
  S10 <- t(S01)

  M <- solve(S11, S10 %*% solve(S00, S01))
  eig <- eigen(M)
  lambda <- sort(Re(eig$values), decreasing = TRUE)
  lambda <- pmin(pmax(lambda, 0), 1 - 1e-12)

  K <- ncol(Y)
  trace_stats <- vapply(0:(K - 1L), function(r) {
    -n * sum(log(1 - lambda[(r + 1L):K]))
  }, numeric(1))
  max_stats <- -n * log(1 - lambda)
  list(
    eigenvalues = lambda,
    trace = setNames(trace_stats, paste0("r<=", 0:(K - 1L))),
    max_eigen = setNames(max_stats, paste0(0:(K - 1L), " vs ", 1:K)),
    n = n
  )
}

j_sim <- johansen_eigen_demo(Y)
j_sim
```

```
## $eigenvalues
## [1] 0.181497909 0.005569738
## 
## $trace
##       r<=0       r<=1 
## 102.726452   2.787068 
## 
## $max_eigen
##    0 vs 1    1 vs 2 
## 99.939385  2.787068 
## 
## $n
## [1] 499
```

較大的第一特徵值與較小的第二特徵值符合共整合秩為 1 的資料生成結構，但正式的共整合秩判定必須使用與確定性規格、樣本數及落後階數相符的非標準臨界值。

## 已知 VECM 參數的矩陣核對


``` r
alpha <- matrix(c(-0.2, 0.1), ncol = 1)
beta <- matrix(c(1, -1), ncol = 1)
Pi <- alpha %*% t(beta)
eigen(Pi)$values
```

```
## [1] -0.3  0.0
```

``` r
Pi
```

```
##      [,1] [,2]
## [1,] -0.2  0.2
## [2,]  0.1 -0.1
```

``` r
previous_spread <- 0.5
expected_change <- as.numeric(alpha * previous_spread)
new_spread <- previous_spread + expected_change[1] - expected_change[2]
c(dy1 = expected_change[1], dy2 = expected_change[2], new_spread = new_spread)
```

```
##        dy1        dy2 new_spread 
##      -0.10       0.05       0.35
```

## 固定日本短長期利率：探索性示範


``` r
path <- "data/processed/japan_monthly_2007_2018.csv"
jp <- read.csv(path, stringsAsFactors = FALSE)
jp$date <- as.Date(jp$date)
jp <- jp[order(jp$date), ]
keep <- complete.cases(jp[, c("inr", "yield_10")])
Y_jp <- as.matrix(jp[keep, c("yield_10", "inr")])
dates <- jp$date[keep]
colnames(Y_jp) <- c("Yield10", "ShortRate")

matplot(dates, Y_jp, type = "l", lty = 1,
        col = c("#173B57", "#A34045"),
        xlab = "Date", ylab = "Percent")
legend("topright", colnames(Y_jp), col = c("#173B57", "#A34045"),
       lty = 1, bty = "n")
```

![plot of chunk japan-data](./R17_cointegration_vecm_files/figure-gfm/japan-data-1.png)


``` r
spread_jp <- Y_jp[, 1] - Y_jp[, 2]
rbind(
  Yield10 = adf_regression(Y_jp[, 1], diff_lags = 1),
  ShortRate = adf_regression(Y_jp[, 2], diff_lags = 1),
  Spread = adf_regression(spread_jp, diff_lags = 1)
)
```

```
##           rho.LevelLag t_stat.LevelLag
## Yield10    -0.01617041      -1.1870966
## ShortRate  -0.06474252      -3.4780026
## Spread     -0.01290881      -0.8331237
```

``` r
johansen_eigen_demo(Y_jp)
```

```
## $eigenvalues
## [1] 0.07896105 0.01295491
## 
## $trace
##      r<=0      r<=1 
## 12.578610  1.721221 
## 
## $max_eigen
##    0 vs 1    1 vs 2 
## 10.857389  1.721221 
## 
## $n
## [1] 132
```

這段輸出不應被寫成「已證明日本利率共整合」。樣本只有約十一年，零利率制度、結構變動、deterministic terms 與 lag choice 都可能改變判讀。

## 簡單 ECM 與殘差診斷


``` r
long_jp <- lm(Yield10 ~ ShortRate, data = as.data.frame(Y_jp))
ect_jp <- residuals(long_jp)
d_jp <- diff(Y_jp)
ect_lag_jp <- ect_jp[-length(ect_jp)]

ecm_long <- lm(d_jp[, "Yield10"] ~ ect_lag_jp)
ecm_short <- lm(d_jp[, "ShortRate"] ~ ect_lag_jp)

summary(ecm_long)$coef
```

```
##                Estimate Std. Error   t value  Pr(>|t|)
## (Intercept) -0.01117912 0.00714045 -1.565604 0.1198723
## ect_lag_jp  -0.02107412 0.01551867 -1.357985 0.1768207
```

``` r
summary(ecm_short)$coef
```

```
##                 Estimate  Std. Error    t value  Pr(>|t|)
## (Intercept) -0.003402197 0.002424581 -1.4032099 0.1629382
## ect_lag_jp  -0.001851880 0.005269455 -0.3514368 0.7258296
```

``` r
par(mfrow = c(1, 2))
acf(residuals(ecm_long), main = "Long-rate ECM residual ACF")
acf(residuals(ecm_short), main = "Short-rate ECM residual ACF")
```

![plot of chunk japan-ecm](./R17_cointegration_vecm_files/figure-gfm/japan-ecm-1.png)

``` r
par(mfrow = c(1, 1))
```

## 報告與退場規則

1. 先說明各 series 的整合階數證據與 deterministic terms。
2. VECM 的差分 lags 比水準 VAR lags 少一階。
3. 報告跡檢定與最大特徵值檢定所使用的確切規格與臨界值；本教學函數未內建它們。
4. 比較前先把共整合向量正規化為同一尺度與符號。
5. 共整合秩對落後階數、趨勢或子樣本翻轉時，不把單一結果包裝成長期法則。
6. 共整合不等於套利；VECM 創新項也不等於結構衝擊。


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
