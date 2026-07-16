---
title: "R11：S&P 500 平衡子樣本的 PCA 與跨股票 lag 修正"
output:
  github_document:
    toc: true
    toc_depth: 3
---

本附錄對應第 15 章，使用專案內凍結的日報酬面板，示範如何先修正跨股票使用 `lag()` 的錯誤，再以嚴格的時間切分進行主成分分析（PCA）。資料只包含 89 檔在共同日期均有觀察值的股票，**不是完整的 S&P 500 指數成分股歷史**；平衡化可能帶來成分股選擇與存活者偏誤。因此，以下結果只供方法教學，不用來衡量整體市場績效，也不構成交易策略證據。

依 `data/DATA_SOURCES.md` 的已決定政策，公開網站提供程式與執行結果，但不再散布凍結價格衍生面板；讀者須自行取得合法資料，再執行同一套清理程式。


``` r
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.width = 7, fig.height = 4
)
stopifnot(getRversion() >= "4.3.0")
set.seed(1111)
```

## 1. 找到凍結資料並核對檔案

本檔不使用 `setwd()`、不安裝套件，也不在執行時下載資料。下列函數容許讀者從專案根目錄或 `online_appendix/` 執行檔案。


``` r
locate_project_file <- function(relative_path) {
  candidates <- c(
    relative_path,
    file.path("..", relative_path),
    file.path("../..", relative_path)
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    stop("找不到專案檔案：", relative_path)
  }
  normalizePath(hit[1], mustWork = TRUE)
}

returns_file <- locate_project_file(
  "data/processed/sp500_returns_balanced_2013_2022.csv"
)
manifest_file <- locate_project_file("data/processed/manifest.csv")
```


``` r
return_df <- read.csv(
  returns_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)

dates <- as.Date(return_df$date)
R_all <- as.matrix(return_df[, -1, drop = FALSE])
storage.mode(R_all) <- "double"

manifest_key <- "data/processed/sp500_returns_balanced_2013_2022.csv"
manifest_row <- manifest[manifest$file == manifest_key, , drop = FALSE]
actual_md5 <- unname(tools::md5sum(returns_file))

stopifnot(
  nrow(manifest_row) == 1L,
  nrow(return_df) == manifest_row$rows,
  ncol(return_df) == manifest_row$columns,
  identical(actual_md5, manifest_row$md5),
  nrow(return_df) == 2384L,
  ncol(R_all) == 89L,
  !anyNA(dates),
  !anyDuplicated(dates),
  all(diff(dates) > 0),
  !anyNA(R_all),
  all(is.finite(R_all)),
  !anyDuplicated(colnames(R_all))
)

data.frame(
  first_date = min(dates),
  last_date = max(dates),
  trading_days = nrow(R_all),
  stocks = ncol(R_all),
  md5 = actual_md5
)
```

```
##   first_date  last_date trading_days stocks                              md5
## 1 2013-01-03 2022-06-22         2384     89 09c9690effb82b3fabdccaa982397e83
```

檢查碼的用途是確認本次分析讀到的確實是教材所凍結的版本；它不能替代資料授權或經濟意義的查核。

## 2. 為何不能在整張長表直接做 `lag()`

若長格式價格表依日期與股票代碼排序後，直接對整欄價格取落後值，某一列的「前一期」往往會是另一檔股票。下列微型資料刻意重現此錯誤。


``` r
toy_price <- data.frame(
  date = rep(as.Date("2024-01-01") + 0:2, each = 2),
  symbol = rep(c("AAA", "BBB"), times = 3),
  price = c(100, 50, 110, 45, 121, 49.5)
)
toy_price <- toy_price[order(toy_price$date, toy_price$symbol), ]

# 錯誤：完全沒有按股票分組。
toy_price$return_wrong <- c(
  NA_real_,
  diff(toy_price$price) / head(toy_price$price, -1)
)
toy_price
```

```
##         date symbol price return_wrong
## 1 2024-01-01    AAA 100.0           NA
## 2 2024-01-01    BBB  50.0   -0.5000000
## 3 2024-01-02    AAA 110.0    1.2000000
## 4 2024-01-02    BBB  45.0   -0.5909091
## 5 2024-01-03    AAA 121.0    1.6888889
## 6 2024-01-03    BBB  49.5   -0.5909091
```

例如 `BBB` 第一天的錯誤報酬拿 `AAA` 第一天的價格當分母；`AAA` 第二天又拿 `BBB` 第一天當分母。正確做法是先按股票分組、在每一組內按日期排序，然後才計算報酬。


``` r
within_symbol_return <- function(data) {
  required <- c("date", "symbol", "price")
  stopifnot(all(required %in% names(data)))

  groups <- split(data, data$symbol, drop = TRUE)
  corrected <- lapply(groups, function(one_stock) {
    one_stock <- one_stock[order(one_stock$date), , drop = FALSE]
    one_stock$return <- c(
      NA_real_,
      diff(one_stock$price) / head(one_stock$price, -1)
    )
    one_stock
  })
  corrected <- do.call(rbind, corrected)
  rownames(corrected) <- NULL
  corrected[order(corrected$date, corrected$symbol), ]
}

toy_correct <- within_symbol_return(toy_price[, c("date", "symbol", "price")])
toy_correct
```

```
##         date symbol price return
## 1 2024-01-01    AAA 100.0     NA
## 4 2024-01-01    BBB  50.0     NA
## 2 2024-01-02    AAA 110.0    0.1
## 5 2024-01-02    BBB  45.0   -0.1
## 3 2024-01-03    AAA 121.0    0.1
## 6 2024-01-03    BBB  49.5    0.1
```

``` r
# 可執行的方向與分組單元測試。
first_in_group <- !duplicated(toy_correct$symbol)
stopifnot(
  all(is.na(toy_correct$return[first_in_group])),
  isTRUE(all.equal(
    toy_correct$return[toy_correct$symbol == "AAA"][-1],
    c(0.10, 0.10), tolerance = 1e-12
  )),
  isTRUE(all.equal(
    toy_correct$return[toy_correct$symbol == "BBB"][-1],
    c(-0.10, 0.10), tolerance = 1e-12
  ))
)
```

本專案的凍結面板正是由「先按股票代碼分組，再計算報酬」的建檔程序產生；原始價格的第一期因無法計算報酬而不會硬填成 0。

## 3. 描述統計與時間切分

先查看每檔股票日報酬的平均數、標準差與極端值。這些摘要不等同於投資績效比較，尤其資料是經平衡化後的子樣本。


``` r
stock_summary <- data.frame(
  symbol = colnames(R_all),
  mean = colMeans(R_all),
  sd = apply(R_all, 2, sd),
  minimum = apply(R_all, 2, min),
  maximum = apply(R_all, 2, max)
)

stock_summary[order(stock_summary$sd, decreasing = TRUE)[1:10], ]
```

```
##      symbol         mean         sd    minimum   maximum
## AMD     AMD 0.0021343661 0.03689526 -0.2422907 0.5229008
## TSLA   TSLA 0.0025735523 0.03590150 -0.2106283 0.2439505
## MU       MU 0.0013139785 0.02875429 -0.1981856 0.1334149
## NVDA   NVDA 0.0020572043 0.02730032 -0.1875588 0.2980671
## PXD     PXD 0.0007199667 0.02655846 -0.3691970 0.2043432
## EOG     EOG 0.0006389926 0.02515343 -0.3200724 0.1657025
## BA       BA 0.0006138625 0.02438963 -0.2384841 0.2431861
## FTNT   FTNT 0.0013749644 0.02425128 -0.1926457 0.2191214
## COP     COP 0.0005927422 0.02348466 -0.2484006 0.2521384
## NEM     NEM 0.0004797220 0.02309381 -0.1222814 0.1401824
```

為避免未來資料影響中心化、尺度與主成分負荷量，本例固定切成連續三段：前 65% 為估計期、接著 15% 為驗證期、最後 20% 為測試期。不能先隨機打散日資料再切分。


``` r
T_total <- nrow(R_all)
train_end <- floor(0.65 * T_total)
validation_end <- floor(0.80 * T_total)

train_id <- seq_len(train_end)
validation_id <- seq.int(train_end + 1L, validation_end)
test_id <- seq.int(validation_end + 1L, T_total)

split_table <- data.frame(
  sample = c("估計期", "驗證期", "測試期"),
  first_date = dates[c(min(train_id), min(validation_id), min(test_id))],
  last_date = dates[c(max(train_id), max(validation_id), max(test_id))],
  observations = c(length(train_id), length(validation_id), length(test_id))
)
split_table
```

```
##   sample first_date  last_date observations
## 1 估計期 2013-01-03 2019-02-28         1549
## 2 驗證期 2019-03-01 2020-07-30          358
## 3 測試期 2020-07-31 2022-06-22          477
```

## 4. 只以估計期建立 PCA

個股波動尺度不同，因此本例先用**估計期**平均數與標準差把各欄標準化，再做 PCA。驗證期與測試期不得各自重新標準化。


``` r
pca_train <- prcomp(
  R_all[train_id, , drop = FALSE],
  center = TRUE,
  scale. = TRUE
)

eigenvalues <- pca_train$sdev^2
pve <- eigenvalues / sum(eigenvalues)
explained <- data.frame(
  component = seq_along(pve),
  eigenvalue = eigenvalues,
  PVE = pve,
  cumulative_PVE = cumsum(pve)
)
head(explained, 15)
```

```
##    component eigenvalue        PVE cumulative_PVE
## 1          1 30.5947665 0.34376142      0.3437614
## 2          2  4.5354152 0.05095972      0.3947211
## 3          3  3.0965127 0.03479228      0.4295134
## 4          4  2.5156583 0.02826582      0.4577792
## 5          5  1.9382370 0.02177794      0.4795572
## 6          6  1.6861223 0.01894519      0.4985024
## 7          7  1.4989376 0.01684200      0.5153444
## 8          8  1.3802086 0.01550796      0.5308523
## 9          9  1.2246514 0.01376013      0.5446125
## 10        10  1.1079264 0.01244861      0.5570611
## 11        11  1.0760845 0.01209084      0.5691519
## 12        12  1.0379729 0.01166262      0.5808145
## 13        13  0.9993702 0.01122888      0.5920434
## 14        14  0.9749119 0.01095407      0.6029975
## 15        15  0.8995359 0.01010715      0.6131046
```


``` r
Z_train <- scale(
  R_all[train_id, , drop = FALSE],
  center = pca_train$center,
  scale = pca_train$scale
)
direct_eigen <- eigen(cor(Z_train), symmetric = TRUE, only.values = TRUE)$values

stopifnot(isTRUE(all.equal(
  unname(eigenvalues), unname(direct_eigen), tolerance = 1e-9
)))
```


``` r
plot(
  explained$component[1:30],
  explained$eigenvalue[1:30],
  type = "b", pch = 19, col = "#173B57",
  xlab = "主成分", ylab = "特徵值",
  main = "估計期相關矩陣的前 30 個特徵值"
)
abline(h = 1, lty = 2, col = "#A34045")
```

![plot of chunk scree-plot](./R11_sp500_pca_corrected_files/figure-gfm/scree-plot-1.png)

### 4.1 負荷量與正負號不定性

特徵向量整欄乘以 (-1) 仍是同一個主成分。因此，不能把 PC1 負荷量的正負號當成可識別的經濟方向。下表依絕對值列出前兩個主成分的重要股票，並保留原符號以利重現。


``` r
largest_loadings <- function(fit, component, n = 12L) {
  loading <- fit$rotation[, component]
  keep <- order(abs(loading), decreasing = TRUE)[seq_len(n)]
  data.frame(
    component = paste0("PC", component),
    symbol = names(loading)[keep],
    loading = unname(loading[keep]),
    absolute_loading = abs(unname(loading[keep]))
  )
}

rbind(
  largest_loadings(pca_train, 1),
  largest_loadings(pca_train, 2)
)
```

```
##    component symbol    loading absolute_loading
## 1        PC1    HON  0.1440392        0.1440392
## 2        PC1     MA  0.1360299        0.1360299
## 3        PC1    MMC  0.1348022        0.1348022
## 4        PC1   FISV  0.1340674        0.1340674
## 5        PC1    TMO  0.1328993        0.1328993
## 6        PC1    ADP  0.1328473        0.1328473
## 7        PC1    DHR  0.1325424        0.1325424
## 8        PC1    ACN  0.1303539        0.1303539
## 9        PC1      V  0.1302504        0.1302504
## 10       PC1    JPM  0.1294404        0.1294404
## 11       PC1    EMR  0.1271807        0.1271807
## 12       PC1    ETN  0.1266206        0.1266206
## 13       PC2    DUK -0.3254492        0.3254492
## 14       PC2      D -0.2962972        0.2962972
## 15       PC2    NEE -0.2918643        0.2918643
## 16       PC2    SRE -0.2434229        0.2434229
## 17       PC2    PEP -0.2341949        0.2341949
## 18       PC2     KO -0.2203414        0.2203414
## 19       PC2     PG -0.2077691        0.2077691
## 20       PC2     VZ -0.1907459        0.1907459
## 21       PC2     PM -0.1756574        0.1756574
## 22       PC2   NVDA  0.1510703        0.1510703
## 23       PC2    AMT -0.1465163        0.1465163
## 24       PC2     MU  0.1441544        0.1441544
```

將估計期 PC1 分數與 89 檔股票的等權平均報酬相比，可作為「共同移動」的描述性核對；這個等權組合不是 S&P 500 指數。


``` r
equal_weight_train <- rowMeans(R_all[train_id, , drop = FALSE])
pc1_train <- pca_train$x[, 1]

data.frame(
  correlation = cor(pc1_train, equal_weight_train),
  absolute_correlation = abs(cor(pc1_train, equal_weight_train))
)
```

```
##   correlation absolute_correlation
## 1   0.9945361            0.9945361
```

## 5. 事先鎖定維度規則

若以同一組變數的重建誤差選維度，保留全部 89 個主成分必然達到零誤差。為避免看過測試期才改規則，本例在分析前指定：選擇估計期累積解釋比例首次達 80% 的最小主成分個數。


``` r
target_pve <- 0.80
r_selected <- which(explained$cumulative_PVE >= target_pve)[1]

data.frame(
  target_PVE = target_pve,
  selected_components = r_selected,
  achieved_training_PVE = explained$cumulative_PVE[r_selected]
)
```

```
##   target_PVE selected_components achieved_training_PVE
## 1        0.8                  39             0.8006484
```

80% 不是自然法則；它是透明、可事先登錄的教學規則。研究者也可用經濟損失函數或只使用驗證期的預測準則，但不能以最終測試期調參。

## 6. 驗證期診斷與最終測試

下列函數固定使用估計樣本所得到的中心、尺度與負荷量，把新資料投影到前 (r) 個主成分後再重建。


``` r
standardized_reconstruction <- function(fit, newdata, r) {
  stopifnot(r >= 0L, r <= ncol(fit$rotation))
  Z <- scale(newdata, center = fit$center, scale = fit$scale)
  if (r == 0L) {
    Z_hat <- matrix(0, nrow = nrow(Z), ncol = ncol(Z))
  } else {
    V <- fit$rotation[, seq_len(r), drop = FALSE]
    Z_hat <- (Z %*% V) %*% t(V)
  }
  list(actual = Z, reconstructed = Z_hat)
}

reconstruction_diagnostics <- function(object) {
  sse <- sum((object$actual - object$reconstructed)^2)
  sst <- sum(object$actual^2)
  c(
    standardized_MSE = mean((object$actual - object$reconstructed)^2),
    out_of_sample_fraction_reconstructed = 1 - sse / sst
  )
}
```

先用驗證期確認所選維度的表現，但不再改動 80% 規則。


``` r
r_grid <- sort(unique(c(1L, 2L, 5L, 10L, 20L, r_selected, 40L, 60L, 89L)))
r_grid <- r_grid[r_grid <= ncol(R_all)]

validation_table <- do.call(rbind, lapply(r_grid, function(r) {
  diagnostic <- reconstruction_diagnostics(
    standardized_reconstruction(
      pca_train,
      R_all[validation_id, , drop = FALSE],
      r
    )
  )
  data.frame(r = r, t(diagnostic), row.names = NULL)
}))
validation_table
```

```
##    r standardized_MSE out_of_sample_fraction_reconstructed
## 1  1     1.499017e+00                            0.5571620
## 2  2     1.347329e+00                            0.6019735
## 3  5     1.144295e+00                            0.6619536
## 4 10     9.659564e-01                            0.7146381
## 5 20     7.976416e-01                            0.7643615
## 6 39     5.120067e-01                            0.8487435
## 7 40     4.967072e-01                            0.8532633
## 8 60     2.357579e-01                            0.9303527
## 9 89     7.503907e-30                            1.0000000
```

鎖定 (r) 後，把估計期與驗證期合併重新估計中心、尺度與負荷量；最終測試期仍完全未參與估計或選擇。


``` r
development_id <- c(train_id, validation_id)
pca_development <- prcomp(
  R_all[development_id, , drop = FALSE],
  center = TRUE,
  scale. = TRUE
)

test_object <- standardized_reconstruction(
  pca_development,
  R_all[test_id, , drop = FALSE],
  r_selected
)
test_diagnostic <- reconstruction_diagnostics(test_object)

data.frame(
  selected_components = r_selected,
  test_standardized_MSE = unname(test_diagnostic[1]),
  test_fraction_reconstructed = unname(test_diagnostic[2])
)
```

```
##   selected_components test_standardized_MSE test_fraction_reconstructed
## 1                  39             0.2543771                   0.7877665
```

這是「同一期觀察值的降維重建」，不是 (t-1) 期預測 (t) 期報酬。若要宣稱預測能力，還須只用當時可取得的過去資料估計因子動態，並以真正的樣本外預測損失評估。

## 7. 測試期分數的時間圖

測試期分數用發展期的固定負荷量計算，不會重新估計 PCA。圖形只呈現共同波動的相對尺度，PC 的正負號仍可整欄翻轉。


``` r
Z_test <- scale(
  R_all[test_id, , drop = FALSE],
  center = pca_development$center,
  scale = pca_development$scale
)
test_scores <- Z_test %*% pca_development$rotation[, 1:2, drop = FALSE]

matplot(
  dates[test_id], test_scores,
  type = "l", lty = 1, col = c("#173B57", "#A34045"),
  xlab = "日期", ylab = "主成分分數",
  main = "以發展期負荷量計算的測試期分數"
)
legend(
  "topright", legend = c("PC1", "PC2"),
  col = c("#173B57", "#A34045"), lty = 1, bty = "n"
)
```

![plot of chunk test-score-plot](./R11_sp500_pca_corrected_files/figure-gfm/test-score-plot-1.png)

## 8. 可重現性與解讀清單

執行真實資料 PCA 時，至少應逐項確認：

1. 報酬是否確實在股票內計算，沒有跨股票落後值；
2. 股票宇宙與日期平衡規則是否清楚，是否有成分股或存活者偏誤；
3. 規格選擇時，中心、尺度與負荷量是否只由估計期決定；鎖定維度後的最終重估是否只用發展期（估計期＋驗證期），完全不使用測試期；
4. 維度規則是否在查看測試期前鎖定；
5. 是否把同日重建誤稱為未來報酬預測；
6. 是否承認 PCA 正負號與未旋轉共同子空間的可識別限制；
7. 原始資料是否允許再散布。


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
