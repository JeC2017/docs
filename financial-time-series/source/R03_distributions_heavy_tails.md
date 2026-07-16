---
title: "R03：分配、厚尾與實證特徵"
output:
  github_document:
    toc: true
    toc_depth: 3
---

本附錄對應第 3–4 章，要回答兩個彼此相關的問題：金融報酬的尾端是否比常態分配厚，以及大幅波動是否會在時間上成群出現？我們先用固定種子模擬常態與標準化 Student-$t$ 分配，確認「變異數相同，尾端仍可很不一樣」；再把同一套描述工具用在兩組真實日簡單報酬上。

第一組資料是 1986–2008 年的 Microsoft（MSFT），每一列代表一個交易日；第二組是由 89 檔股票組成的 2013–2022 年等權教學投資組合，每一列先記錄同一共同交易日的 89 個個股報酬，再取橫斷面平均。MSFT 資料源自 Ruey S. Tsay 教科書網站的 `d-msft8608.txt`，股票面板源自原課程 S&P 500 價格檔；固定版本與建置方式見 `data/DATA_SOURCES.md`。

兩組報酬都以小數表示，0.01 代表 1%。等權序列是教學用組合，不是官方 S&P 500 指數，也沒有校正成分股生存者偏誤。本頁使用完整歷史樣本做分配與相依性描述，不進行預測，因此不分訓練期、驗證期與測試期；結果也沒有識別任何市場事件對報酬的因果效果。


``` r
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.width = 8, fig.height = 5,
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

## 描述統計與 Jarque–Bera 函數

先寫兩個小函數，目的是讓每一個輸出都能與公式對上。`sample_moments()` 接受一條報酬序列，輸出有效觀察值數、平均數、標準差、偏態、峰度與尾端分位數；`jarque_bera()` 再用偏態與超額峰度組成檢定統計量。以下峰度採用「常態分配等於 3」的定義，超額峰度則是峰度減 3。有限樣本偏態與峰度有不同修正版本，因此比較手動版與套件版以前，必須先確認定義。


``` r
sample_moments <- function(x) {
  # 描述同一條序列時使用相同的有限值樣本，避免各統計量樣本數不同。
  x <- x[is.finite(x)]
  centered <- x - mean(x)
  m2 <- mean(centered^2)
  m3 <- mean(centered^3)
  m4 <- mean(centered^4)
  c(
    n = length(x),
    mean = mean(x),
    sd = sd(x),
    skewness = m3 / m2^(3 / 2),
    kurtosis = m4 / m2^2,
    excess_kurtosis = m4 / m2^2 - 3,
    q01 = unname(quantile(x, 0.01)),
    q05 = unname(quantile(x, 0.05)),
    median = median(x),
    q95 = unname(quantile(x, 0.95)),
    q99 = unname(quantile(x, 0.99))
  )
}

jarque_bera <- function(x) {
  m <- sample_moments(x)
  # 此處使用教科書中的 iid 漸近式；時間相依性稍後另做診斷。
  statistic <- m["n"] * (
    m["skewness"]^2 / 6 + m["excess_kurtosis"]^2 / 24
  )
  c(
    statistic = unname(statistic),
    p_value_iid_asymptotic = unname(
      pchisq(statistic, 2, lower.tail = FALSE)
    )
  )
}
```

## 先用已知分配看懂「厚尾」

模擬的好處是母體分配已知，可以先排除資料建置與模型選擇的干擾。自由度 5 的 Student-$t$ 變異數是 $5/(5-2)$；除以理論標準差後，它與標準常態具有相同變異數，但尾端機率仍較高。這一段要觀察的是峰度、尾端分位數與直方圖形狀，不是檢定金融市場。


``` r
n_sim <- 10000L
normal_draw <- rnorm(n_sim)
# 標準化後兩組理論變異數相同，尾端差異才不會只是尺度造成。
t5_draw <- rt(n_sim, df = 5) / sqrt(5 / 3)

simulation_summary <- rbind(
  Normal = sample_moments(normal_draw),
  Student_t5 = sample_moments(t5_draw)
)
knitr::kable(round(simulation_summary, 4))
```



|           |     n|   mean|     sd| skewness| kurtosis| excess_kurtosis|     q01|     q05| median|    q95|    q99|
|:----------|-----:|------:|------:|--------:|--------:|---------------:|-------:|-------:|------:|------:|------:|
|Normal     | 10000| 0.0101| 1.0136|  -0.0016|   3.0640|          0.0640| -2.4048| -1.6583| 0.0099| 1.6794| 2.3828|
|Student_t5 | 10000| 0.0007| 1.0129|  -0.1179|   8.1864|          5.1864| -2.6003| -1.5806| 0.0158| 1.5817| 2.6540|

樣本平均數與變異數會因這一次亂數抽樣略有誤差。標準化 Student-$t_5$ 的 1% 與 99% 分位數通常比常態分位數離零更遠，但 5% 與 95% 分位數反而更靠近零，不能籠統地說所有尾端分位數都更遠。直方圖呈現的是三個一起發生的特徵：中心尖峰較高、肩部較薄，而極端尾端較厚；峰度則把這種極端尾端的差異濃縮成一個數字。


``` r
old_par <- par(
  mfrow = c(1, 2), mar = c(4, 4, 2, 1),
  family = plot_family
)
hist(
  normal_draw, breaks = 80, probability = TRUE,
  xlim = c(-6, 6), col = "#9FC2D4", border = "white",
  main = "標準常態", xlab = "模擬值"
)
curve(dnorm(x), add = TRUE, lwd = 2, col = "#A34045")
hist(
  t5_draw, breaks = 100, probability = TRUE,
  xlim = c(-6, 6), col = "#D6B0A9", border = "white",
  main = "標準化 Student-t(5)", xlab = "模擬值"
)
curve(dnorm(x), add = TRUE, lwd = 2, col = "#173B57")
```

![相同理論變異數下，標準化 Student-t(5) 的尾端比常態分配厚。](../R03_distributions_heavy_tails_files/figure-gfm/simulation-plots-1.png)

``` r
par(old_par)
```

模擬讓我們看懂「相同變異數不代表相同尾端」，卻沒有提供 MSFT 或投資組合的實證結論。接下來才把這套觀察方式帶到固定的歷史資料。

## 讀取 MSFT 與股票面板

讀入兩份資料後，先依日期排序並確認沒有缺值。`R` 的每一列是一個共同交易日，每一欄是一檔股票；`rowMeans(R)` 因此是在同一天對 89 檔股票等權平均，而不是跨時間平均。這個方向若弄反，會得到完全不同的研究對象。


``` r
msft <- read.csv(project_path(
  "data", "processed", "msft_daily_returns_1986_2008.csv"
))
msft$date <- as.Date(msft$date)
# ACF 依賴時間次序，所以在任何診斷前先由早到晚排序。
msft <- msft[order(msft$date), ]

panel <- read.csv(
  project_path(
    "data", "processed", "sp500_returns_balanced_2013_2022.csv"
  ),
  check.names = FALSE
)
panel$date <- as.Date(panel$date)
R <- as.matrix(panel[, setdiff(names(panel), "date")])
storage.mode(R) <- "double"

stopifnot(
  !anyNA(msft$date), !anyNA(msft$simple_return),
  all(diff(msft$date) > 0),
  !anyNA(panel$date), !anyNA(R),
  all(diff(panel$date) > 0)
)

sp_equal <- rowMeans(R)

data_profile <- data.frame(
  序列 = c("MSFT", "89 檔股票等權教學投資組合"),
  起日 = c(min(msft$date), min(panel$date)),
  迄日 = c(max(msft$date), max(panel$date)),
  觀察值 = c(nrow(msft), nrow(panel)),
  資產數 = c(1L, ncol(R)),
  單位 = "日簡單報酬，小數",
  來源 = c(
    "Tsay 教科書網站 d-msft8608.txt",
    "原課程 S&P 500 價格檔的平衡面板"
  ),
  check.names = FALSE
)
knitr::kable(data_profile)
```



|序列                      |起日       |迄日       | 觀察值| 資產數|單位             |來源                            |
|:-------------------------|:----------|:----------|------:|------:|:----------------|:-------------------------------|
|MSFT                      |1986-03-14 |2008-12-31 |   5752|      1|日簡單報酬，小數 |Tsay 教科書網站 d-msft8608.txt  |
|89 檔股票等權教學投資組合 |2013-01-03 |2022-06-22 |   2384|     89|日簡單報酬，小數 |原課程 S&P 500 價格檔的平衡面板 |

`data_profile` 應顯示 MSFT 有 5,752 個交易日，等權組合有 2,384 個共同交易日。兩列的日期範圍不同，因此後面的統計量是在各自樣本內描述，不宜把數值差異全部歸因於資產本身；市場年代、樣本長度與成分股篩選也可能影響結果。

## 真實報酬的偏態、峰度與尾端

現在用同一組公式整理兩條真實報酬序列。分位數告訴我們左右尾的實際尺度，偏態描述兩側是否對稱，峰度則對極端觀察值特別敏感；三者要放在一起讀，不能只看單一 $p$ 值。


``` r
empirical_summary <- rbind(
  MSFT = sample_moments(msft$simple_return),
  SP_equal_weight = sample_moments(sp_equal)
)
knitr::kable(round(empirical_summary, 5))
```



|                |    n|    mean|      sd| skewness| kurtosis| excess_kurtosis|      q01|      q05|  median|     q95|     q99|
|:---------------|----:|-------:|-------:|--------:|--------:|---------------:|--------:|--------:|-------:|-------:|-------:|
|MSFT            | 5752| 0.00123| 0.02359| -0.13209| 12.92329|         9.92329| -0.06108| -0.03306| 0.00000| 0.03846| 0.06344|
|SP_equal_weight | 2384| 0.00076| 0.01075| -0.67951| 22.94458|        19.94458| -0.03040| -0.01555| 0.00106| 0.01469| 0.02507|

``` r
jb_table <- rbind(
  MSFT = jarque_bera(msft$simple_return),
  SP_equal_weight = jarque_bera(sp_equal)
)
knitr::kable(jb_table, digits = 6)
```



|                | statistic| p_value_iid_asymptotic|
|:---------------|---------:|----------------------:|
|MSFT            |  23617.11|                      0|
|SP_equal_weight |  39696.89|                      0|

MSFT 的樣本峰度約為 12.92，等權教學組合約為 22.94，都遠高於常態分配的 3；兩者的 Jarque–Bera 統計量也都很大。就這兩段固定樣本而言，常態分配無法妥善描述偏態與尾端。這項結果只指出常態模型的不足，還沒有在 Student-$t$、偏態 $t$ 或其他厚尾分配之間選出唯一答案。

Jarque–Bera 的卡方近似依賴獨立同分配與有限高階動差。金融報酬常有條件異質變異，因此極小的漸近 $p$ 值適合當作常態模型的警訊，後面仍要檢查報酬與平方報酬的時間相依性。

## 套件作法：用 `basicStats()` 與 `normalTest()` 整理分配特徵

原課程的實作程式在
`slides/L02_Return_properties/W1L2_R_scripts_Descriptive_stat_returns.R`
以 `fBasics::basicStats()` 整理描述統計，並以
`normalTest(..., method = "jb")` 執行 Jarque–Bera 檢定。下列程式沿用這個工作流程，但改讀本書的固定資料，不在執行時連線下載。套件會代為計算統計量並整理 `fHTEST` 輸出；資料期間、報酬尺度、峰度定義與時間相依性是否可忽略，仍要由使用者判斷。


``` r
stopifnot(requireNamespace("fBasics", quietly = TRUE))

empirical_series <- list(
  MSFT = msft$simple_return,
  `S&P 等權教學組合` = sp_equal
)

basic_row <- function(out, row_name) {
  out <- as.matrix(out)
  location <- match(tolower(row_name), tolower(row.names(out)))
  stopifnot(!is.na(location))
  as.numeric(out[location, 1])
}

fbasics_summary <- lapply(empirical_series, fBasics::basicStats)
manual_summary <- lapply(empirical_series, sample_moments)

# 按統計量名稱取值，避免套件輸出列次序改變時配錯數字。
moment_comparison <- do.call(rbind, lapply(names(empirical_series), function(nm) {
  package_result <- as.matrix(fbasics_summary[[nm]])
  manual_result <- manual_summary[[nm]]
  manual_values <- c(
    manual_result["mean"],
    manual_result["sd"],
    manual_result["skewness"],
    manual_result["excess_kurtosis"]
  )
  package_values <- c(
    basic_row(package_result, "Mean"),
    basic_row(package_result, "Stdev"),
    basic_row(package_result, "Skewness"),
    basic_row(package_result, "Kurtosis")
  )
  data.frame(
    序列 = nm,
    統計量 = c("平均數", "標準差", "偏態", "超額峰度"),
    手動版 = unname(manual_values),
    fBasics = package_values,
    套件減手動 = package_values - unname(manual_values),
    check.names = FALSE
  )
}))
row.names(moment_comparison) <- NULL
knitr::kable(moment_comparison, digits = 7)
```



|序列             |統計量   |     手動版|   fBasics| 套件減手動|
|:----------------|:--------|----------:|---------:|----------:|
|MSFT             |平均數   |  0.0012318|  0.001232|  0.0000002|
|MSFT             |標準差   |  0.0235944|  0.023594| -0.0000004|
|MSFT             |偏態     | -0.1320928| -0.132058|  0.0000348|
|MSFT             |超額峰度 |  9.9232911|  9.918798| -0.0044931|
|S&P 等權教學組合 |平均數   |  0.0007612|  0.000761| -0.0000002|
|S&P 等權教學組合 |標準差   |  0.0107466|  0.010747|  0.0000004|
|S&P 等權教學組合 |偏態     | -0.6795124| -0.679085|  0.0004274|
|S&P 等權教學組合 |超額峰度 | 19.9445772| 19.925332| -0.0192452|

``` r
# 對同一批固定觀察值執行套件版 JB，才能與手動公式比較。
fbasics_jb <- lapply(empirical_series, function(x) {
  fBasics::normalTest(x, method = "jb")
})

# 先看完整 fHTEST，認識統計量、自由度與 p 值在物件中的位置。
for (nm in names(fbasics_jb)) {
  cat("\n", nm, "\n", sep = "")
  print(fbasics_jb[[nm]])
}
```

```
## 
## MSFT
## 
## Title:
##  Jarque-Bera Normality Test
## 
## Test Results:
##   STATISTIC:
##     X-squared: 23617.1128
##   P VALUE:
##     Asymptotic p Value: < 2.2e-16 
## 
## 
## S&P 等權教學組合
## 
## Title:
##  Jarque-Bera Normality Test
## 
## Test Results:
##   STATISTIC:
##     X-squared: 39696.8887
##   P VALUE:
##     Asymptotic p Value: < 2.2e-16
```

``` r
extract_fhtest <- function(x) {
  test_result <- methods::slot(x, "test")
  c(
    statistic = unname(test_result$statistic),
    p_value = unname(test_result$p.value)
  )
}

jb_comparison <- do.call(rbind, lapply(names(empirical_series), function(nm) {
  manual_result <- jarque_bera(empirical_series[[nm]])
  package_result <- extract_fhtest(fbasics_jb[[nm]])
  data.frame(
    序列 = nm,
    手動JB = manual_result["statistic"],
    fBasics_JB = package_result["statistic"],
    手動p值 = manual_result["p_value_iid_asymptotic"],
    fBasics_p值 = package_result["p_value"],
    check.names = FALSE
  )
}))
row.names(jb_comparison) <- NULL
knitr::kable(jb_comparison, digits = 7)
```



|序列             |   手動JB| fBasics_JB| 手動p值| fBasics_p值|
|:----------------|--------:|----------:|-------:|-----------:|
|MSFT             | 23617.11|   23617.11|       0|           0|
|S&P 等權教學組合 | 39696.89|   39696.89|       0|           0|

平均數與標準差先幫我們確認資料與尺度是否一致；偏態、峰度與 Jarque–Bera 數值則可能因有限樣本修正與動差定義不同而略有差異。特別是 `basicStats()` 的 `Kurtosis` 採超額峰度，應與上文「常態等於 0」的欄位比較，不能直接對照「常態等於 3」的峰度。表中若出現小幅數值差，先查定義與修正方式；兩種作法在這份固定資料上的方向應一致，都指出常態模型對尾端的描述不足。


``` r
series_list <- list(
  MSFT = list(date = msft$date, return = msft$simple_return),
  `S&P 等權教學組合` = list(date = panel$date, return = sp_equal)
)

old_par <- par(
  mfrow = c(2, 3), mar = c(4.5, 3.7, 3, 1),
  family = plot_family
)
for (nm in names(series_list)) {
  z <- series_list[[nm]]
  plot(
    z$date, 100 * z$return, type = "l", col = "#173B57",
    xlab = "日期", ylab = "日報酬（%）", main = nm
  )
  hist(
    z$return, breaks = 70, probability = TRUE,
    col = "#9FC2D4", border = "white",
    xlab = "日簡單報酬", main = paste(nm, "分配")
  )
  curve(
    dnorm(x, mean(z$return), sd(z$return)),
    add = TRUE, lwd = 2, col = "#A34045"
  )
  qqnorm(
    z$return, pch = 16, cex = 0.4,
    col = "#173B57", main = paste(nm, "常態 Q–Q")
  )
  qqline(z$return, col = "#A34045", lwd = 2)
}
```

![MSFT 與等權教學投資組合的時間圖、經驗分配與常態 Q–Q 圖。](../R03_distributions_heavy_tails_files/figure-gfm/empirical-plots-1.png)

``` r
par(old_par)
```

時間圖先告訴我們極端波動落在哪些日期，直方圖比較整體形狀，Q–Q 圖則把左右尾偏離常態的方向分開呈現。若尾端點系統性離開參考線，應考慮厚尾創新或更合適的風險衡量方式；圖形本身還不能決定唯一的替代分配。

## 報酬與平方報酬的時間相依

厚尾描述的是無條件分配，波動群聚則問「大波動之後是否較常接著大波動」。因此同時畫原報酬與平方報酬的自相關函數（ACF）：前者主要看條件平均的線性相依，後者把正負號拿掉，用來觀察波動大小是否具有持續性。


``` r
old_par <- par(
  mfrow = c(2, 2), mar = c(4, 3.5, 4, 1),
  family = plot_family, cex.main = 0.85
)
for (nm in names(series_list)) {
  z <- series_list[[nm]]$return
  acf(z, lag.max = 30, main = paste(nm, "報酬 ACF"))
  acf(z^2, lag.max = 30, main = paste(nm, "平方報酬 ACF"))
}
```

![MSFT 與等權教學投資組合的報酬 ACF 與平方報酬 ACF。](../R03_distributions_heavy_tails_files/figure-gfm/acf-diagnostics-1.png)

``` r
par(old_par)
```


``` r
lb_row <- function(x, series_name) {
  q_return <- Box.test(x, lag = 20, type = "Ljung-Box")
  q_square <- Box.test(x^2, lag = 20, type = "Ljung-Box")
  data.frame(
    序列 = series_name,
    Q20_報酬 = unname(q_return$statistic),
    p_報酬 = q_return$p.value,
    Q20_平方報酬 = unname(q_square$statistic),
    p_平方報酬 = q_square$p.value,
    check.names = FALSE
  )
}

diagnostics <- rbind(
  lb_row(msft$simple_return, "MSFT"),
  lb_row(sp_equal, "S&P 等權教學組合")
)
knitr::kable(diagnostics, digits = 6)
```



|序列             | Q20_報酬|   p_報酬| Q20_平方報酬| p_平方報酬|
|:----------------|--------:|--------:|------------:|----------:|
|MSFT             |  41.3004| 0.003408|     1156.800|          0|
|S&P 等權教學組合 | 422.1970| 0.000000|     4037.729|          0|

在這份固定資料中，MSFT 報酬的 $Q(20)$ 約為 41.30，等權組合約為 422.20；平方報酬的統計量更分別達約 1,156.80 與 4,037.73。前 20 階同時為零的限制在兩條序列上都遭拒絕，平方報酬的相依尤其明顯。下一步可分開處理兩個問題：條件平均模型檢查殘差是否仍有線性相依，ARCH/GARCH 類模型則處理波動持續性。這些統計量沒有直接證明報酬可交易，也沒有提供造成波動群聚的因果機制。

## 定位極端日期

由於峰度對少數極端值很敏感，最後把絕對報酬最大的日期列出來。這一步是為了回查資料來源與理解哪些觀察值主導尾端，而不是根據結果大小自動刪除資料。


``` r
largest_moves <- function(date, x, label, number = 5L) {
  keep <- order(abs(x), decreasing = TRUE)[seq_len(number)]
  data.frame(
    序列 = label,
    日期 = date[keep],
    日簡單報酬 = x[keep],
    check.names = FALSE
  )
}

extreme_table <- rbind(
  largest_moves(msft$date, msft$simple_return, "MSFT"),
  largest_moves(panel$date, sp_equal, "S&P 等權教學組合")
)
knitr::kable(extreme_table, digits = 6)
```



|序列             |日期       | 日簡單報酬|
|:----------------|:----------|----------:|
|MSFT             |1987-10-19 |  -0.301158|
|MSFT             |2000-10-19 |   0.195652|
|MSFT             |1987-10-26 |  -0.186529|
|MSFT             |2008-10-13 |   0.186047|
|MSFT             |1987-10-21 |   0.179688|
|S&P 等權教學組合 |2020-03-16 |  -0.121581|
|S&P 等權教學組合 |2020-03-24 |   0.105691|
|S&P 等權教學組合 |2020-03-12 |  -0.096315|
|S&P 等權教學組合 |2020-03-13 |   0.087688|
|S&P 等權教學組合 |2020-03-09 |  -0.080871|

若回查後確認是建檔或價格調整錯誤，應修正資料並重新計算全部統計量；若極端值與原始來源一致，就應保留，因為它們正是厚尾與風險分析要解釋的現象。

## 這兩段歷史資料給我們的訊息

MSFT 與等權教學組合都呈現高峰度、常態 Q–Q 圖尾端偏離，以及顯著的平方報酬相依。這些證據支持在後續模型中正面處理厚尾與條件異質變異，卻還不能判定哪一個厚尾分配最好，也不能把歷史尾端頻率當成不變的未來機率。實務上的下一步，是在明確的訓練期內比較候選分配與波動模型，再用保留的測試期評估預測分配，而不是只依一次常態性檢定作決定。
