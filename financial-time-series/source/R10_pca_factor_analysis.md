---
title: "R10：PCA、低秩重建與因子分析"
output:
  github_document:
    toc: true
    toc_depth: 3
---

本附錄對應第 15 章，以固定種子的合成因子資料說明主成分分析（PCA）、低秩重建、因子個數選擇與最大概似因子分析。模擬的優點是共同因子與負荷量真值可核對；它不是任何真實市場的實證結果。


``` r
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.width = 7, fig.height = 4
)
stopifnot(getRversion() >= "4.3.0")
set.seed(1010)
```

## 1. 軟體與時間切分

本檔只使用 R 內建的 `stats` 與 `graphics`，不安裝套件、不下載資料、不使用 `setwd()`。700 期資料依時間固定切成：

- 1--400：估計期；
- 401--520：驗證期，只用來選主成分個數；
- 521--700：最終測試期。


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

## 2. 產生低秩共同結構

令

\[
\mathbf x_t=\mathbf B\mathbf f_t+\mathbf u_t,
\]

其中三個因子具有不同持續性，十個觀察變數有不同負荷量與個別變異。


``` r
T_total <- 700
N <- 10
r_true <- 3

simulate_ar1 <- function(n, phi, sd = 1) {
  stopifnot(abs(phi) < 1, sd > 0)
  x <- numeric(n)
  innovation_sd <- sd * sqrt(1 - phi^2)
  x[1] <- rnorm(1, sd = sd)
  for (t in 2:n) {
    x[t] <- phi * x[t - 1] + rnorm(1, sd = innovation_sd)
  }
  x
}

F_true <- cbind(
  F1 = simulate_ar1(T_total, phi = 0.55),
  F2 = simulate_ar1(T_total, phi = 0.25),
  F3 = simulate_ar1(T_total, phi = -0.20)
)

B_true <- rbind(
  c(0.90, 0.10, 0.00),
  c(0.85, 0.15, 0.05),
  c(0.75, 0.20, 0.10),
  c(0.10, 0.90, 0.05),
  c(0.15, 0.85, 0.10),
  c(0.20, 0.75, 0.15),
  c(0.05, 0.10, 0.90),
  c(0.10, 0.15, 0.85),
  c(0.20, 0.10, 0.75),
  c(0.45, 0.40, 0.35)
)
rownames(B_true) <- paste0("X", seq_len(N))
colnames(B_true) <- colnames(F_true)

idio_sd <- seq(0.35, 0.65, length.out = N)
U <- matrix(rnorm(T_total * N), nrow = T_total)
U <- sweep(U, 2, idio_sd, "*")
X <- F_true %*% t(B_true) + U
colnames(X) <- rownames(B_true)

estimate_id <- 1:400
validation_id <- 401:520
test_id <- 521:700
```


``` r
round(cov(X[estimate_id, ]), 2)
```

```
##       X1   X2   X3   X4   X5   X6   X7   X8   X9  X10
## X1  0.98 0.82 0.79 0.29 0.33 0.33 0.06 0.08 0.12 0.52
## X2  0.82 0.95 0.78 0.31 0.35 0.35 0.08 0.11 0.13 0.54
## X3  0.79 0.78 0.94 0.38 0.41 0.40 0.13 0.17 0.20 0.54
## X4  0.29 0.31 0.38 1.09 0.85 0.74 0.14 0.20 0.18 0.44
## X5  0.33 0.35 0.41 0.85 1.04 0.72 0.16 0.20 0.20 0.48
## X6  0.33 0.35 0.40 0.74 0.72 0.87 0.16 0.21 0.18 0.46
## X7  0.06 0.08 0.13 0.14 0.16 0.16 1.06 0.72 0.68 0.34
## X8  0.08 0.11 0.17 0.20 0.20 0.21 0.72 0.99 0.62 0.37
## X9  0.12 0.13 0.20 0.18 0.20 0.18 0.68 0.62 0.97 0.35
## X10 0.52 0.54 0.54 0.44 0.48 0.46 0.34 0.37 0.35 0.92
```

## 3. PCA 的特徵分解

`prcomp()` 以奇異值分解實作 PCA。這裡先標準化各欄，因此等價於對估計期相關矩陣做 PCA。


``` r
pca_estimate <- prcomp(
  X[estimate_id, , drop = FALSE],
  center = TRUE,
  scale. = TRUE
)

eigenvalues <- pca_estimate$sdev^2
PVE <- eigenvalues / sum(eigenvalues)
explained <- data.frame(
  component = seq_along(PVE),
  eigenvalue = eigenvalues,
  PVE = PVE,
  cumulative_PVE = cumsum(PVE)
)
explained
```

```
##    component eigenvalue        PVE cumulative_PVE
## 1          1  4.5348975 0.45348975      0.4534897
## 2          2  2.1167328 0.21167328      0.6651630
## 3          3  1.5119799 0.15119799      0.8163610
## 4          4  0.4102693 0.04102693      0.8573879
## 5          5  0.3653930 0.03653930      0.8939272
## 6          6  0.2930908 0.02930908      0.9232363
## 7          7  0.2507491 0.02507491      0.9483112
## 8          8  0.1965901 0.01965901      0.9679702
## 9          9  0.1687588 0.01687588      0.9848461
## 10        10  0.1515388 0.01515388      1.0000000
```


``` r
plot(
  explained$component, explained$eigenvalue,
  type = "b", pch = 19, col = "#173B57",
  xlab = "主成分", ylab = "特徵值",
  main = "只用估計期建立的陡坡圖"
)
abline(h = 1, lty = 2, col = "#A34045")
```

![plot of chunk scree-plot](./R10_pca_factor_analysis_files/figure-gfm/scree-plot-1.png)

### 3.1 以矩陣運算核對


``` r
Z_estimate <- scale(
  X[estimate_id, , drop = FALSE],
  center = pca_estimate$center,
  scale = pca_estimate$scale
)
S_estimate <- cov(Z_estimate)
eigen_direct <- eigen(S_estimate, symmetric = TRUE)

stopifnot(
  isTRUE(all.equal(
    unname(eigenvalues),
    unname(eigen_direct$values),
    tolerance = 1e-10
  ))
)

# 特徵向量可整欄反號，所以比較投影外積。
projection_prcomp <- tcrossprod(pca_estimate$rotation[, 1:3])
projection_eigen <- tcrossprod(eigen_direct$vectors[, 1:3])
max(abs(projection_prcomp - projection_eigen))
```

```
## [1] 6.661338e-16
```

## 4. 主成分與真因子的關係

PCA 不知道真因子名稱；只能由資料找共同子空間。下表取絕對相關，避免特徵向量任意反號造成誤解。


``` r
score_estimate <- pca_estimate$x[, 1:5, drop = FALSE]
round(abs(cor(score_estimate, F_true[estimate_id, ])), 3)
```

```
##        F1    F2    F3
## PC1 0.679 0.646 0.347
## PC2 0.415 0.129 0.834
## PC3 0.537 0.682 0.168
## PC4 0.096 0.023 0.024
## PC5 0.030 0.010 0.045
```

即使前三個主成分與三個真因子高度相關，欄次序也可能交換或混合。可識別的通常是共同子空間，不是未加限制的每一個因子名稱。

## 5. 無資料洩漏的低秩重建

### 5.1 固定估計期中心、尺度與負荷量


``` r
reconstruct_from_pca <- function(fit, newdata, r) {
  stopifnot(r >= 0, r <= ncol(fit$rotation))
  Z <- scale(newdata, center = fit$center, scale = fit$scale)
  if (r == 0) {
    Z_hat <- matrix(0, nrow = nrow(Z), ncol = ncol(Z))
  } else {
    V <- fit$rotation[, seq_len(r), drop = FALSE]
    score <- Z %*% V
    Z_hat <- score %*% t(V)
  }
  X_hat <- sweep(Z_hat, 2, fit$scale, "*")
  sweep(X_hat, 2, fit$center, "+")
}

reconstruction_mse <- function(actual, reconstructed) {
  mean((actual - reconstructed)^2)
}
```

### 5.2 只用驗證期選 \(r\)


``` r
r_grid <- 0:N
validation_mse <- vapply(r_grid, function(r) {
  X_hat <- reconstruct_from_pca(
    pca_estimate,
    X[validation_id, , drop = FALSE],
    r = r
  )
  reconstruction_mse(X[validation_id, ], X_hat)
}, numeric(1))

selection_table <- data.frame(
  r = r_grid,
  validation_MSE = validation_mse
)
selection_table
```

```
##     r validation_MSE
## 1   0   8.777857e-01
## 2   1   5.548663e-01
## 3   2   3.534909e-01
## 4   3   1.889625e-01
## 5   4   1.552875e-01
## 6   5   1.134014e-01
## 7   6   7.763963e-02
## 8   7   5.067196e-02
## 9   8   2.800476e-02
## 10  9   1.357753e-02
## 11 10   5.744278e-31
```

若只最小化同一組變數的重建誤差，保留全部 \(N\) 個主成分必然最好，因為沒有降維懲罰。為讓問題反映「用較少維度換取多少誤差」，本例使用一個事先指定的規則：選擇重建 MSE 不超過完整 PCA 驗證誤差加上總變異 10% 的最小 \(r\)。另一種方法是預先要求累積解釋比例達 80% 或 90%。


``` r
validation_total_variance <- mean(apply(
  X[validation_id, , drop = FALSE], 2, var
))
tolerance <- validation_mse[N + 1] + 0.10 * validation_total_variance
r_selected <- min(r_grid[validation_mse <= tolerance])
r_selected
```

```
## [1] 6
```

規則與容許誤差必須在查看測試期之前固定。

### 5.3 鎖定 \(r\) 後重新估計，最後一次評估測試期

可在選定 \(r\) 後合併估計與驗證期重新估計 PCA；測試期仍未參與。


``` r
development_id <- c(estimate_id, validation_id)
pca_development <- prcomp(
  X[development_id, , drop = FALSE],
  center = TRUE,
  scale. = TRUE
)
X_test_hat <- reconstruct_from_pca(
  pca_development,
  X[test_id, , drop = FALSE],
  r = r_selected
)
test_mse <- reconstruction_mse(X[test_id, ], X_test_hat)

data.frame(
  selected_r = r_selected,
  validation_tolerance = tolerance,
  final_test_MSE = test_mse
)
```

```
##   selected_r validation_tolerance final_test_MSE
## 1          6           0.08689857     0.07810104
```

測試期資料用來計算當期主成分分數與重建，這是同日降維；若要在前一期預測本期因子，還須對分數建立時間序列模型。

## 6. 因子分析與 varimax 旋轉

`factanal()` 是 R 內建的最大概似因子分析。它把每個變數變異拆成共同性與個別變異，目標不同於 PCA。


``` r
fa_three <- factanal(
  X[estimate_id, , drop = FALSE],
  factors = 3,
  rotation = "varimax",
  scores = "regression"
)
fa_three
```

```
## 
## Call:
## factanal(x = X[estimate_id, , drop = FALSE], factors = 3, scores = "regression",     rotation = "varimax")
## 
## Uniquenesses:
##    X1    X2    X3    X4    X5    X6    X7    X8    X9   X10 
## 0.149 0.152 0.185 0.192 0.207 0.271 0.267 0.327 0.388 0.434 
## 
## Loadings:
##     Factor1 Factor2 Factor3
## X1  0.910   0.149          
## X2  0.901   0.186          
## X3  0.859   0.257   0.109  
## X4  0.157   0.880          
## X5  0.216   0.857   0.107  
## X6  0.258   0.803   0.131  
## X7                  0.854  
## X8          0.122   0.809  
## X9  0.104   0.102   0.768  
## X10 0.533   0.386   0.366  
## 
##                Factor1 Factor2 Factor3
## SS loadings      2.814   2.455   2.158
## Proportion Var   0.281   0.246   0.216
## Cumulative Var   0.281   0.527   0.743
## 
## Test of the hypothesis that 3 factors are sufficient.
## The chi square statistic is 8.95 on 18 degrees of freedom.
## The p-value is 0.961
```


``` r
loadings_matrix <- unclass(fa_three$loadings)
round(loadings_matrix, 3)
```

```
##     Factor1 Factor2 Factor3
## X1    0.910   0.149   0.017
## X2    0.901   0.186   0.045
## X3    0.859   0.257   0.109
## X4    0.157   0.880   0.088
## X5    0.216   0.857   0.107
## X6    0.258   0.803   0.131
## X7    0.033   0.058   0.854
## X8    0.060   0.122   0.809
## X9    0.104   0.102   0.768
## X10   0.533   0.386   0.366
```

``` r
communality <- rowSums(loadings_matrix^2)
factor_summary <- data.frame(
  variable = rownames(loadings_matrix),
  communality = communality,
  uniqueness = fa_three$uniquenesses
)
factor_summary
```

```
##     variable communality uniqueness
## X1        X1   0.8511288  0.1488711
## X2        X2   0.8477483  0.1522515
## X3        X3   0.8153370  0.1846630
## X4        X4   0.8075700  0.1924301
## X5        X5   0.7930710  0.2069290
## X6        X6   0.7285846  0.2714154
## X7        X7   0.7331248  0.2668756
## X8        X8   0.6731115  0.3268863
## X9        X9   0.6117841  0.3882197
## X10      X10   0.5661422  0.4338530
```

因子旋轉可讓負荷量群組更清楚，但不會把統計因子自動變成已識別的經濟衝擊。若改變旋轉法，單一欄的名稱與符號可能改變，共同重建空間則可相近。

## 7. 共變異數 PCA 與相關矩陣 PCA

用相同估計期比較 `scale.=FALSE` 與 `scale.=TRUE`。由於本例各欄個別變異不同，兩組負荷量不會完全相同。


``` r
pca_cov <- prcomp(
  X[estimate_id, , drop = FALSE],
  center = TRUE, scale. = FALSE
)
pca_cor <- pca_estimate

comparison <- data.frame(
  variable = colnames(X),
  covariance_PC1 = pca_cov$rotation[, 1],
  correlation_PC1 = pca_cor$rotation[, 1]
)
comparison
```

```
##     variable covariance_PC1 correlation_PC1
## X1        X1     -0.3346065      -0.3385481
## X2        X2     -0.3423127      -0.3517636
## X3        X3     -0.3621640      -0.3732021
## X4        X4     -0.3520707      -0.3269951
## X5        X5     -0.3609405      -0.3445571
## X6        X6     -0.3345797      -0.3527932
## X7        X7     -0.2059879      -0.1915774
## X8        X8     -0.2218401      -0.2150694
## X9        X9     -0.2214534      -0.2179987
## X10      X10     -0.3644677      -0.3765357
```


``` r
matplot(
  seq_len(N),
  cbind(abs(comparison$covariance_PC1),
        abs(comparison$correlation_PC1)),
  type = "b", pch = c(16, 1), lty = 1,
  col = c("#173B57", "#A34045"),
  xaxt = "n", xlab = "變數", ylab = "|PC1 負荷量|",
  main = "尺度選擇會改變 PCA 問題"
)
axis(1, at = seq_len(N), labels = colnames(X))
legend(
  "topright",
  legend = c("共變異數 PCA", "相關矩陣 PCA"),
  col = c("#173B57", "#A34045"),
  pch = c(16, 1), lty = 1, bty = "n"
)
```

![plot of chunk loading-plot](./R10_pca_factor_analysis_files/figure-gfm/loading-plot-1.png)

## 8. 可重現結論

本附錄驗證：

1. `prcomp()` 的奇異值與相關矩陣特徵值一致；
2. 主成分符號可反轉，應比較投影空間而非單一正負號；
3. 驗證期只用來選 \(r\)，測試期直到最後才使用；
4. 測試期標準化沿用發展期中心與尺度；
5. PCA 與最大概似因子分析的目標、個別變異處理與旋轉解讀不同。


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
