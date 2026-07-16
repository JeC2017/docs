---
title: "R18：VAR、局部投影與進階 IRF"
output:
  github_document:
    toc: true
    toc_depth: 3
---

本附錄對應選讀第 23 章。它以**已知且可觀察的外生衝擊**建立教學模擬，比較 VAR 與局部投影（LP）IRF，並示範重疊預測期距的 HAC 標準誤。後段提供 DML/AIPW 評分函數與狀態相依性的小型核對例；不宣稱重現專案中三篇 2024--2026 年論文。

## 執行條件

- base R 與 `knitr`；不安裝套件、不下載資料。
- 所有結果由固定種子教學模擬產生。
- 衝擊識別由資料生成過程明示；若換成觀察性政策變數，這項條件不會自動保留。


``` r
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.width = 7, fig.height = 4.6
)
set.seed(20260716)
```

## 模擬結構 VAR 與可觀察衝擊


``` r
simulate_svarma <- function(Tn, A, B, C = matrix(0, nrow(B), ncol(B)),
                            burn = 300L) {
  K <- nrow(A)
  eps <- matrix(rnorm((Tn + burn) * K), ncol = K)
  Y <- matrix(0, Tn + burn, K)
  for (t in 2:nrow(Y)) {
    Y[t, ] <- A %*% Y[t - 1, ] + B %*% eps[t, ] + C %*% eps[t - 1, ]
  }
  keep <- (burn + 1L):(burn + Tn)
  list(Y = Y[keep, , drop = FALSE], shock = eps[keep, 1], eps = eps[keep, ])
}

true_irf <- function(A, B, C, horizon) {
  out <- vector("list", horizon + 1L)
  out[[1]] <- B
  if (horizon >= 1L) out[[2]] <- A %*% B + C
  if (horizon >= 2L) {
    for (h in 2:horizon) out[[h + 1L]] <- A %*% out[[h]]
  }
  out
}

extract_response <- function(irf_list, response = 1L, shock = 1L) {
  vapply(irf_list, function(M) M[response, shock], numeric(1))
}
```


``` r
A <- matrix(c(0.55, 0.12,
              -0.10, 0.35), 2, 2, byrow = TRUE)
B <- matrix(c(1.00, 0.00,
              0.40, 0.75), 2, 2, byrow = TRUE)
C0 <- matrix(0, 2, 2)
H <- 16L

sim <- simulate_svarma(900, A, B, C0)
Y <- sim$Y
shock <- sim$shock
colnames(Y) <- c("y1", "y2")
truth <- extract_response(true_irf(A, B, C0, H), 1, 1)
```

## LP 與 Newey--West covariance

對每個 horizon (h)，估

\[
y_{1,t+h}=a_h+\theta_hs_t+\Gamma_h'W_{t-1}+u_{t+h}^{(h)},
\]

其中控制變數是 \((Y_{t-1},\ldots,Y_{t-p})\)。衝擊在模擬中直接觀察且與創新項獨立，所以 \(\theta_h\) 有明確的因果 IRF 待估對象。


``` r
newey_west_vcov <- function(X, residuals, bandwidth) {
  score <- X * as.numeric(residuals)
  meat <- crossprod(score)
  n <- nrow(score)
  if (bandwidth > 0L) {
    for (lag in seq_len(min(bandwidth, n - 1L))) {
      weight <- 1 - lag / (bandwidth + 1)
      gamma <- crossprod(score[(lag + 1L):n, , drop = FALSE],
                         score[1:(n - lag), , drop = FALSE])
      meat <- meat + weight * (gamma + t(gamma))
    }
  }
  bread <- solve(crossprod(X))
  bread %*% meat %*% bread
}

estimate_lp <- function(Y, shock, horizon_max = 12L, p = 2L,
                        bandwidth_rule = function(h) max(h, 4L)) {
  Y <- as.matrix(Y)
  if (is.null(colnames(Y))) colnames(Y) <- paste0("y", seq_len(ncol(Y)))
  Tn <- nrow(Y)
  out <- data.frame(h = 0:horizon_max, theta = NA_real_, se_hac = NA_real_,
                    n = NA_integer_, bandwidth = NA_integer_)

  for (h in 0:horizon_max) {
    t_index <- (p + 1L):(Tn - h)
    y_h <- Y[t_index + h, 1]
    controls <- do.call(cbind, lapply(seq_len(p), function(lag) {
      Y[t_index - lag, , drop = FALSE]
    }))
    colnames(controls) <- unlist(lapply(seq_len(p), function(lag) {
      paste0(colnames(Y), "_L", lag)
    }))
    X <- cbind(Intercept = 1, Shock = shock[t_index], controls)
    fit <- lm.fit(X, y_h)
    L <- bandwidth_rule(h)
    V <- newey_west_vcov(X, fit$residuals, L)
    out[h + 1L, c("theta", "se_hac", "n", "bandwidth")] <-
      c(fit$coefficients["Shock"], sqrt(V["Shock", "Shock"]), length(y_h), L)
  }
  out
}

lp <- estimate_lp(Y, shock, H, p = 2)
head(lp)
```

```
##   h      theta       se_hac   n bandwidth
## 1 0 1.00000000 7.837684e-18 898         4
## 2 1 0.55460214 3.363868e-02 897         4
## 3 2 0.29132998 3.888439e-02 896         4
## 4 3 0.23219743 3.904792e-02 895         4
## 5 4 0.14518145 3.905648e-02 894         4
## 6 5 0.06446433 3.843893e-02 893         5
```

## VAR IRF


``` r
fit_var1 <- function(Y) {
  Y <- as.matrix(Y)
  X <- cbind(Intercept = 1, Y[-nrow(Y), , drop = FALSE])
  Ydep <- Y[-1, , drop = FALSE]
  coef <- qr.solve(X, Ydep)
  residuals <- Ydep - X %*% coef
  list(
    intercept = coef[1, ],
    A = t(coef[-1, , drop = FALSE]),
    sigma = crossprod(residuals) / nrow(residuals),
    residuals = residuals
  )
}

var1_irf <- function(fit, horizon) {
  B_hat <- t(chol(fit$sigma))
  out <- vector("list", horizon + 1L)
  out[[1]] <- B_hat
  if (horizon >= 1L) {
    for (h in seq_len(horizon)) out[[h + 1L]] <- fit$A %*% out[[h]]
  }
  out
}

fit_v <- fit_var1(Y)
var_estimate <- extract_response(var1_irf(fit_v, H), 1, 1)
```


``` r
ylim <- range(truth, var_estimate,
              lp$theta - 1.96 * lp$se_hac,
              lp$theta + 1.96 * lp$se_hac)
plot(0:H, truth, type = "l", lwd = 2, col = "black", ylim = ylim,
     xlab = "Horizon", ylab = "Response of y1 to shock 1")
lines(0:H, var_estimate, lwd = 2, lty = 2, col = "#173B57")
points(0:H, lp$theta, pch = 16, cex = 0.7, col = "#A34045")
segments(0:H, lp$theta - 1.96 * lp$se_hac,
         0:H, lp$theta + 1.96 * lp$se_hac,
         col = adjustcolor("#A34045", alpha.f = 0.45))
legend("topright", c("Truth", "VAR(1)", "LP with pointwise HAC interval"),
       col = c("black", "#173B57", "#A34045"),
       lty = c(1, 2, NA), pch = c(NA, NA, 16), lwd = c(2, 2, NA), bty = "n")
```

![正確低階 VAR DGP：真實、VAR 與 LP response。](./R18_local_projection_irf_files/figure-gfm/compare-correct-dgp-1.png)

在正確且簡約的 VAR 資料生成過程下，VAR 遞迴通常利用較多跨預測期距限制；LP 點估計可能較不平滑。單次圖不能比較涵蓋率，需重複模擬。

## 小型蒙地卡羅實驗：偏誤與標準差


``` r
monte_carlo <- function(repetitions = 100L, Tn = 350L, horizon = 8L,
                        C = matrix(0, 2, 2), seed = 20260716) {
  set.seed(seed)
  lp_draw <- var_draw <- matrix(NA_real_, repetitions, horizon + 1L)
  for (r in seq_len(repetitions)) {
    s <- simulate_svarma(Tn, A, B, C)
    lp_draw[r, ] <- estimate_lp(s$Y, s$shock, horizon, p = 2)$theta
    var_draw[r, ] <- extract_response(var1_irf(fit_var1(s$Y), horizon), 1, 1)
  }
  target <- extract_response(true_irf(A, B, C, horizon), 1, 1)
  data.frame(
    h = 0:horizon,
    truth = target,
    lp_bias = colMeans(lp_draw) - target,
    lp_sd = apply(lp_draw, 2, sd),
    var_bias = colMeans(var_draw) - target,
    var_sd = apply(var_draw, 2, sd)
  )
}

mc_correct <- monte_carlo(repetitions = 80, Tn = 350, horizon = 8, C = C0)
round(mc_correct, 3)
```

```
##   h truth lp_bias lp_sd var_bias var_sd
## 1 0 1.000   0.000 0.000    0.000  0.038
## 2 1 0.598  -0.015 0.051   -0.016  0.051
## 3 2 0.334  -0.012 0.065   -0.016  0.051
## 4 3 0.178  -0.019 0.067   -0.012  0.043
## 5 4 0.092  -0.023 0.070   -0.007  0.032
## 6 5 0.046  -0.016 0.073   -0.004  0.022
## 7 6 0.023  -0.020 0.078   -0.002  0.015
## 8 7 0.011  -0.007 0.060    0.000  0.009
## 9 8 0.005  -0.013 0.067    0.000  0.006
```

有限次蒙地卡羅實驗是方法的壓力測試，不是理論證明。增加重複次數可降低模擬雜訊，但會增加編譯時間。

## 局部 MA 動態錯設的教學壓力測試

加入小型 \(C\varepsilon_{t-1}\) 後，資料不再是精確 VAR(1)。真實當期反應為 \(B\)，一期反應為 \(AB+C\)，之後由 \(A\) 傳遞。


``` r
C_small <- 0.20 * B
sim_ma <- simulate_svarma(900, A, B, C_small)
truth_ma <- extract_response(true_irf(A, B, C_small, H), 1, 1)
lp_ma <- estimate_lp(sim_ma$Y, sim_ma$shock, H, p = 2)
var_ma <- extract_response(var1_irf(fit_var1(sim_ma$Y), H), 1, 1)

plot(0:H, truth_ma, type = "l", lwd = 2, col = "black",
     ylim = range(truth_ma, lp_ma$theta, var_ma),
     xlab = "Horizon", ylab = "Response")
lines(0:H, var_ma, lwd = 2, lty = 2, col = "#173B57")
lines(0:H, lp_ma$theta, lwd = 1.5, lty = 3, col = "#A34045")
legend("topright", c("Truth", "Misspecified VAR(1)", "LP"),
       col = c("black", "#173B57", "#A34045"),
       lwd = c(2, 2, 1.5), lty = c(1, 2, 3), bty = "n")
```

![plot of chunk misspecification](./R18_local_projection_irf_files/figure-gfm/misspecification-1.png)

``` r
mc_ma <- monte_carlo(repetitions = 80, Tn = 350, horizon = 8, C = C_small)
round(mc_ma, 3)
```

```
##   h truth lp_bias lp_sd var_bias var_sd
## 1 0 1.000   0.000 0.001    0.009  0.039
## 2 1 0.798  -0.016 0.051   -0.101  0.047
## 3 2 0.453  -0.015 0.071    0.002  0.053
## 4 3 0.245  -0.022 0.076    0.040  0.052
## 5 4 0.128  -0.027 0.079    0.045  0.046
## 6 5 0.065  -0.021 0.082    0.036  0.038
## 7 6 0.032  -0.023 0.088    0.026  0.030
## 8 7 0.016  -0.011 0.070    0.017  0.023
## 9 8 0.008  -0.014 0.073    0.010  0.017
```

這個資料生成過程只是可核對的小例，不能被描述成供應 PDF 的重現，也不能用來證明 LP 對所有錯誤設定都有正確涵蓋率。

## 頻寬敏感度


``` r
lp_Lh <- estimate_lp(Y, shock, H, p = 2,
                     bandwidth_rule = function(h) max(h, 1L))
lp_L4 <- estimate_lp(Y, shock, H, p = 2,
                     bandwidth_rule = function(h) 4L)
lp_L12 <- estimate_lp(Y, shock, H, p = 2,
                      bandwidth_rule = function(h) 12L)

data.frame(
  h = 0:H,
  se_Lh = lp_Lh$se_hac,
  se_L4 = lp_L4$se_hac,
  se_L12 = lp_L12$se_hac
)
```

```
##     h        se_Lh        se_L4       se_L12
## 1   0 8.662842e-18 7.837684e-18 7.628143e-18
## 2   1 3.415455e-02 3.363868e-02 3.150262e-02
## 3   2 3.790702e-02 3.888439e-02 3.928026e-02
## 4   3 3.889454e-02 3.904792e-02 3.872784e-02
## 5   4 3.905648e-02 3.905648e-02 3.758108e-02
## 6   5 3.843893e-02 3.898045e-02 3.804920e-02
## 7   6 3.607567e-02 3.617650e-02 3.613303e-02
## 8   7 4.080768e-02 3.948255e-02 4.045087e-02
## 9   8 4.020096e-02 3.930276e-02 4.008037e-02
## 10  9 4.102776e-02 4.072232e-02 4.109628e-02
## 11 10 3.931643e-02 3.995293e-02 3.851908e-02
## 12 11 3.799919e-02 3.988689e-02 3.753522e-02
## 13 12 3.504898e-02 3.903523e-02 3.504898e-02
## 14 13 3.919067e-02 4.001536e-02 3.928509e-02
## 15 14 3.678129e-02 3.934185e-02 3.713329e-02
## 16 15 3.976508e-02 4.144410e-02 4.049851e-02
## 17 16 4.244003e-02 4.264808e-02 4.244836e-02
```

HAC 頻寬是推論設定的一部分。重疊累積應變數至少要正視 \(h-1\) 階機械相依；本例應變數是未來水準，仍可能因動態控制變數與預測誤差而相關。

## DML/AIPW IRF 評分函數的可執行小例

以下為獨立同分配教學例，只核對專案所附 2024 年 PDF 使用的增廣評分函數形式。真實傾向分數與應變數迴歸由資料生成過程已知，因此不涉及機器學習調校。


``` r
set.seed(20260717)
n_aipw <- 3000L
X_aipw <- rnorm(n_aipw)
e_true <- plogis(-0.2 + 0.6 * X_aipw)
D_aipw <- rbinom(n_aipw, 1, e_true)
theta_aipw <- 0.7
mu0 <- 0.5 * X_aipw
mu1 <- theta_aipw + 0.5 * X_aipw
Y_aipw <- ifelse(D_aipw == 1, mu1, mu0) + rnorm(n_aipw)

aipw_score <- mu1 - mu0 +
  D_aipw / e_true * (Y_aipw - mu1) -
  (1 - D_aipw) / (1 - e_true) * (Y_aipw - mu0)

c(AIPW = mean(aipw_score), truth = theta_aipw,
  min_propensity = min(e_true), max_propensity = max(e_true))
```

```
##           AIPW          truth min_propensity max_propensity 
##      0.6216850      0.7000000      0.1182124      0.8925416
```

若 \(e(X)\) 接近 0 或 1，逆機率項會極端；若處置未滿足給定 \(X\) 的識別條件，評分函數也不能產生因果 IRF。時間序列版本還需要論文所列的相依與時間不變條件。

## 狀態相依性的教學模擬

建立 \(g(z)=0.4+0.3z^2\)，衝擊與預定狀態獨立。線性交互作用只能配適直線；已知的二次基底可以恢復彎曲。專案所附 2026 年 PDF 使用篩法與一致推論，本例不是其估計量的重現。


``` r
set.seed(20260718)
N <- 4000L
Z <- runif(N, -2, 2)
S <- rnorm(N)
g_true <- 0.4 + 0.3 * Z^2
Y_state <- g_true * S + 0.2 * Z + rnorm(N)

linear_fit <- lm(Y_state ~ S * Z + Z)
quadratic_fit <- lm(Y_state ~ S * Z + S:I(Z^2) + Z + I(Z^2))
grid <- seq(-2, 2, length.out = 101)

# 對 shock 係數取解析式；不以兩個任意 shock 值作差。
g_linear <- coef(linear_fit)["S"] + coef(linear_fit)["S:Z"] * grid
g_quad <- coef(quadratic_fit)["S"] + coef(quadratic_fit)["S:Z"] * grid +
  coef(quadratic_fit)["S:I(Z^2)"] * grid^2

plot(grid, 0.4 + 0.3 * grid^2, type = "l", lwd = 2, col = "black",
     xlab = "Predetermined state z", ylab = "Response g(z)")
lines(grid, g_linear, lwd = 2, lty = 2, col = "#A34045")
lines(grid, g_quad, lwd = 2, lty = 3, col = "#173B57")
legend("topleft", c("Truth", "Linear interaction", "Quadratic basis"),
       col = c("black", "#A34045", "#173B57"),
       lwd = 2, lty = c(1, 2, 3), bty = "n")
```

![plot of chunk state-dependent-demo](./R18_local_projection_irf_files/figure-gfm/state-dependent-demo-1.png)

## 解讀與退場規則

1. LP 衝擊必須先被識別；在本附錄由資料生成過程明示，不能移植到觀察性衝擊。
2. VAR/LP 比較必須使用同一待估對象、衝擊、樣本與正規化方式。
3. Pointwise 1.96 intervals 不是整條曲線的 simultaneous bands。
4. MA 壓力測試只說明特定資料生成過程，不概括專案所附 2026 年 LP/VAR 論文的定理。
5. AIPW/DML 仍需 consistency、無未觀察混淆、overlap 與相依資料推論。
6. 狀態相依的因果解讀需要對衝擊的線性條件平均數限制與衝擊外生性；基底選擇也屬估計程序。


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
