---
title: "R02：價格、日期與報酬率"
output:
  github_document:
    toc: true
    toc_depth: 3
---

本附錄對應第 2 章。使用完全固定的假想價格、股利與拆股調整例，示範日期排序、分組落後、簡單報酬、對數報酬、多期複利與投資組合報酬。所有數值都是教學資料，不代表真實資產。


``` r
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.width = 7, fig.height = 4.5
)
stopifnot(requireNamespace("dplyr", quietly = TRUE))
stopifnot(requireNamespace("tibble", quietly = TRUE))
library(dplyr)
library(tibble)
```

## 建立價格、股利與調整因子

price_raw 是觀察到的未調整收盤價；split_factor 表示當日一股舊股可換得多少新股。為了簡化教學，price_adjusted 已把拆股前價格換成可比較尺度。實務資料供應者的調整定義可能不同，必須查資料字典。


``` r
prices <- tribble(
  ~asset, ~date,       ~price_raw, ~dividend, ~split_factor, ~price_adjusted,
  "甲",    "2026-01-02", 100,        0,         1,             50,
  "甲",    "2026-01-05", 104,        0,         1,             52,
  "甲",    "2026-01-06",  53,        0,         2,             53,
  "甲",    "2026-01-07",  52,        1,         1,             52,
  "乙",    "2026-01-02",  80,        0,         1,             80,
  "乙",    "2026-01-05",  82,        0,         1,             82,
  "乙",    "2026-01-06",  81,        0,         1,             81,
  "乙",    "2026-01-07",  84,        0,         1,             84
) |>
  mutate(date = as.Date(date))

prices
```

```
## # A tibble: 8 × 6
##   asset date       price_raw dividend split_factor price_adjusted
##   <chr> <date>         <dbl>    <dbl>        <dbl>          <dbl>
## 1 甲    2026-01-02       100        0            1             50
## 2 甲    2026-01-05       104        0            1             52
## 3 甲    2026-01-06        53        0            2             53
## 4 甲    2026-01-07        52        1            1             52
## 5 乙    2026-01-02        80        0            1             80
## 6 乙    2026-01-05        82        0            1             82
## 7 乙    2026-01-06        81        0            1             81
## 8 乙    2026-01-07        84        0            1             84
```

## 入口檢查


``` r
validate_panel <- function(x) {
  needed <- c("asset", "date", "price_adjusted", "dividend")
  stopifnot(all(needed %in% names(x)))
  stopifnot(!anyNA(x[, needed]))
  stopifnot(all(is.finite(x$price_adjusted)), all(x$price_adjusted > 0))
  stopifnot(all(is.finite(x$dividend)), all(x$dividend >= 0))
  stopifnot(!anyDuplicated(paste(x$asset, x$date)))

  ordered <- x |> arrange(asset, date)
  by_asset <- split(ordered$date, ordered$asset)
  stopifnot(all(vapply(
    by_asset,
    function(z) length(z) < 2L || all(diff(z) > 0),
    logical(1)
  )))
  invisible(ordered)
}

prices <- validate_panel(prices)
```

## 務必先分組再取落後值

若直接對整欄取落後值，甲資產最後一天可能被錯接到乙資產第一天。正確做法是先依 asset 分組，再依日期排序。


``` r
returns <- prices |>
  group_by(asset) |>
  arrange(date, .by_group = TRUE) |>
  mutate(
    previous_price = lag(price_adjusted),
    gross_return = (price_adjusted + dividend) / previous_price,
    simple_return = gross_return - 1,
    log_return = log(gross_return),
    calendar_gap = as.integer(date - lag(date))
  ) |>
  ungroup()

returns
```

```
## # A tibble: 8 × 11
##   asset date       price_raw dividend split_factor price_adjusted previous_price
##   <chr> <date>         <dbl>    <dbl>        <dbl>          <dbl>          <dbl>
## 1 乙    2026-01-02        80        0            1             80             NA
## 2 乙    2026-01-05        82        0            1             82             80
## 3 乙    2026-01-06        81        0            1             81             82
## 4 乙    2026-01-07        84        0            1             84             81
## 5 甲    2026-01-02       100        0            1             50             NA
## 6 甲    2026-01-05       104        0            1             52             50
## 7 甲    2026-01-06        53        0            2             53             52
## 8 甲    2026-01-07        52        1            1             52             53
## # ℹ 4 more variables: gross_return <dbl>, simple_return <dbl>,
## #   log_return <dbl>, calendar_gap <int>
```

``` r
stopifnot(all.equal(
  returns$log_return[!is.na(returns$log_return)],
  log1p(returns$simple_return[!is.na(returns$simple_return)])
))
```

甲資產 2026-01-07 支付 1 元股利，總簡單報酬使用
\((P_t+D_t-P_{t-1})/P_{t-1}\)。若只算價格變化，會遺漏現金分配。

## 未調整價格會製造假報酬


``` r
split_comparison <- prices |>
  filter(asset == "甲") |>
  arrange(date) |>
  mutate(
    return_from_raw = price_raw / lag(price_raw) - 1,
    return_from_adjusted = price_adjusted / lag(price_adjusted) - 1
  ) |>
  select(date, price_raw, price_adjusted, return_from_raw, return_from_adjusted)

split_comparison
```

```
## # A tibble: 4 × 5
##   date       price_raw price_adjusted return_from_raw return_from_adjusted
##   <date>         <dbl>          <dbl>           <dbl>                <dbl>
## 1 2026-01-02       100             50         NA                   NA     
## 2 2026-01-05       104             52          0.0400               0.0400
## 3 2026-01-06        53             53         -0.490                0.0192
## 4 2026-01-07        52             52         -0.0189              -0.0189
```

拆股日的未調整價格約減半，不代表投資人損失一半。正式分析應優先使用明確記錄股利與公司行動的總報酬或調整價格。

## 多期複利與對數相加


``` r
asset_b <- returns |>
  filter(asset == "乙", !is.na(simple_return))

multi_simple <- prod(1 + asset_b$simple_return) - 1
multi_log <- sum(asset_b$log_return)

tibble(
  method = c("簡單報酬複利", "對數報酬相加後轉回"),
  total_simple_return = c(multi_simple, exp(multi_log) - 1)
)
```

```
## # A tibble: 2 × 2
##   method             total_simple_return
##   <chr>                            <dbl>
## 1 簡單報酬複利                    0.0500
## 2 對數報酬相加後轉回              0.0500
```

``` r
stopifnot(isTRUE(all.equal(multi_simple, exp(multi_log) - 1)))
```

簡單報酬不能直接跨期相加；對數報酬可以跨期相加，但轉回簡單報酬時仍要取指數減一。

## 簡單報酬與對數報酬的近似


``` r
grid <- seq(-0.5, 0.5, by = 0.01)
approximation <- tibble(
  simple = grid,
  log = log1p(grid),
  difference = log1p(grid) - grid
)

head(approximation, 4)
```

```
## # A tibble: 4 × 3
##   simple    log difference
##    <dbl>  <dbl>      <dbl>
## 1  -0.5  -0.693     -0.193
## 2  -0.49 -0.673     -0.183
## 3  -0.48 -0.654     -0.174
## 4  -0.47 -0.635     -0.165
```

``` r
tail(approximation, 4)
```

```
## # A tibble: 4 × 3
##   simple   log difference
##    <dbl> <dbl>      <dbl>
## 1   0.47 0.385    -0.0847
## 2   0.48 0.392    -0.0880
## 3   0.49 0.399    -0.0912
## 4   0.5  0.405    -0.0945
```


``` r
plot(
  approximation$simple, approximation$log,
  type = "l", lwd = 2, col = "#173B57",
  xlab = "簡單報酬", ylab = "對數報酬"
)
abline(0, 1, lty = 2, col = "#A34045")
```

![簡單報酬與對數報酬的差距在大幅負報酬時特別明顯。](./R02_prices_dates_returns_files/figure-gfm/approximation-plot-1.png)

## 投資組合報酬使用期初權重

在同一持有期間，簡單報酬的投資組合加總為
\(R_{p,t}=\sum_i w_{i,t-1}R_{i,t}\)。以下把兩資產報酬整理到共同日期，再使用固定期初權重 60\% 與 40\%。


``` r
common_dates <- Reduce(
  intersect,
  split(
    returns$date[!is.na(returns$simple_return)],
    returns$asset[!is.na(returns$simple_return)]
  )
)

ra <- returns |>
  filter(asset == "甲", date %in% common_dates) |>
  arrange(date) |>
  pull(simple_return)
rb <- returns |>
  filter(asset == "乙", date %in% common_dates) |>
  arrange(date) |>
  pull(simple_return)

portfolio <- tibble(
  date = sort(as.Date(common_dates, origin = "1970-01-01")),
  return_a = ra,
  return_b = rb,
  portfolio_return = 0.6 * ra + 0.4 * rb
)
portfolio
```

```
## # A tibble: 3 × 4
##   date       return_a return_b portfolio_return
##   <date>        <dbl>    <dbl>            <dbl>
## 1 2026-01-05   0.0400   0.0250          0.034  
## 2 2026-01-06   0.0192  -0.0122          0.00666
## 3 2026-01-07   0        0.0370          0.0148
```

如果期中不再平衡，權重會隨相對價格變動；不能每期都偷用固定權重卻聲稱是買進持有策略。

## 日期間隔不是固定交易頻率

calendar_gap 在週末前後會大於 1。若資料目標是每日交易報酬，這仍可能是相鄰交易日；若要年度化波動或合併不同市場，必須明確處理交易日曆、時區與休市差異。


``` r
returns |>
  filter(!is.na(calendar_gap)) |>
  count(calendar_gap)
```

```
## # A tibble: 2 × 2
##   calendar_gap     n
##          <int> <int>
## 1            1     4
## 2            3     2
```

## 建議輸出欄位

正式報酬資料至少保留 asset、date、原始價格、調整價格、股利／公司行動、簡單報酬、對數報酬、幣別與來源版本。不要只保存最後一欄報酬，否則很難追查異常值究竟來自市場、公司行動或程式錯置。
