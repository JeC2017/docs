---
title: "R00：環境、專案與可重現研究"
output:
  github_document:
    toc: true
    toc_depth: 3
---

本附錄對應第 1 章。目標是建立可審核的專案習慣：固定資料版本、使用相對路徑、記錄執行環境、讓亂數可重現，並讓錯誤在資料入口立即停止。全文不連線下載資料，也不在執行期間變更套件環境。

## 執行前提

- R 4.1 以上，並已具備 knitr 與 rmarkdown。
- 從教科書專案根目錄或 online_appendix 目錄執行皆可。
- 正式資料位於 data/processed，建置紀錄位於 manifest.csv。


``` r
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.width = 7,
  fig.height = 4.5
)

root_candidates <- c(".", "..")
is_root <- vapply(
  root_candidates,
  function(x) file.exists(file.path(x, "main.tex")) &&
    file.exists(file.path(x, "data", "processed", "manifest.csv")),
  logical(1)
)
stopifnot(any(is_root))
project_root <- root_candidates[which(is_root)[1]]
project_path <- function(...) file.path(project_root, ...)
```

## 專案結構是一項研究設計

建議把外部輸入、清理後資料、分析程式、圖表與正文分開。程式只讀固定輸入，產出寫到指定資料夾；不要讓手動修改的試算表成為無法追蹤的中間步驟。


``` r
required_paths <- c(
  "data/processed/manifest.csv",
  "data/DATA_SOURCES.md",
  "chapters/chapter01.tex",
  "online_appendix/README.md"
)

path_check <- data.frame(
  path = required_paths,
  exists = file.exists(vapply(required_paths, project_path, character(1))),
  row.names = NULL
)
path_check
```

```
##                          path exists
## 1 data/processed/manifest.csv   TRUE
## 2        data/DATA_SOURCES.md   TRUE
## 3      chapters/chapter01.tex   TRUE
## 4   online_appendix/README.md   TRUE
```

``` r
stopifnot(all(path_check$exists))
```

## 用 manifest 核對固定資料

manifest.csv 記錄檔名、列數、欄數、MD5 與建立時間。MD5 不能證明資料正確，但能確認讀者拿到的位元內容是否與教科書版本相同。


``` r
manifest <- read.csv(
  project_path("data", "processed", "manifest.csv"),
  stringsAsFactors = FALSE
)

actual_md5 <- unname(tools::md5sum(project_path(manifest$file)))
verification <- transform(
  manifest,
  actual_md5 = actual_md5,
  md5_match = tolower(md5) == tolower(actual_md5),
  file_exists = file.exists(project_path(file))
)
verification
```

```
##                                                  file  rows columns
## 1 data/processed/sp500_returns_balanced_2013_2022.csv  2384      90
## 2          data/processed/japan_monthly_2007_2018.csv   133      30
## 3 data/processed/ff_qf_macro_industries_1967_2021.csv  6590      24
## 4 data/processed/taiwan_icapm_second_stage_47x497.csv 23359       6
##                                md5                built_at
## 1 09c9690effb82b3fabdccaa982397e83 2026-07-16 07:35:23 UTC
## 2 3fd45a6a7a8d26e29d48f1c2f1497ad8 2026-07-16 07:35:23 UTC
## 3 4d9eea7ddeea063a7b635f238dd7ba24 2026-07-16 07:35:23 UTC
## 4 98fb791d16ee3b3e536ef0ce33381e93 2026-07-16 07:35:23 UTC
##                         actual_md5 md5_match file_exists
## 1 09c9690effb82b3fabdccaa982397e83      TRUE        TRUE
## 2 3fd45a6a7a8d26e29d48f1c2f1497ad8      TRUE        TRUE
## 3 4d9eea7ddeea063a7b635f238dd7ba24      TRUE        TRUE
## 4 98fb791d16ee3b3e536ef0ce33381e93      TRUE        TRUE
```

``` r
stopifnot(all(verification$file_exists), all(verification$md5_match))
```

再核對宣告的列數與欄數。這一步只讀資料，不修改正式檔。


``` r
dimensions <- lapply(manifest$file, function(f) {
  x <- read.csv(project_path(f), check.names = FALSE)
  c(rows = nrow(x), columns = ncol(x))
})
dimensions <- as.data.frame(do.call(rbind, dimensions))
dimensions$file <- manifest$file
dimensions$row_match <- dimensions$rows == manifest$rows
dimensions$column_match <- dimensions$columns == manifest$columns
dimensions[, c("file", "rows", "columns", "row_match", "column_match")]
```

```
##                                                  file  rows columns row_match
## 1 data/processed/sp500_returns_balanced_2013_2022.csv  2384      90      TRUE
## 2          data/processed/japan_monthly_2007_2018.csv   133      30      TRUE
## 3 data/processed/ff_qf_macro_industries_1967_2021.csv  6590      24      TRUE
## 4 data/processed/taiwan_icapm_second_stage_47x497.csv 23359       6      TRUE
##   column_match
## 1         TRUE
## 2         TRUE
## 3         TRUE
## 4         TRUE
```

``` r
stopifnot(all(dimensions$row_match), all(dimensions$column_match))
```

## 亂數種子與可重現模擬

同一個亂數種子、R 版本與亂數產生器設定會重現同一串模擬值。種子不是讓答案顯著的旋鈕；應在分析前固定，敏感度分析則另列多個事前指定種子或增加模擬次數。


``` r
set.seed(20260716)
draw_a <- rnorm(5)

set.seed(20260716)
draw_b <- rnorm(5)

data.frame(draw_a, draw_b, exactly_equal = draw_a == draw_b)
```

```
##       draw_a     draw_b exactly_equal
## 1 -0.3134385 -0.3134385          TRUE
## 2  0.5719575  0.5719575          TRUE
## 3 -0.1357020 -0.1357020          TRUE
## 4  2.1122746  2.1122746          TRUE
## 5  1.3737218  1.3737218          TRUE
```

``` r
stopifnot(identical(draw_a, draw_b))
```

## 讓錯誤提早發生

以下函數示範金融時間序列的最小資料契約：日期可轉換、嚴格遞增、沒有重複，價格為正且有限。


``` r
validate_price_data <- function(x) {
  needed <- c("date", "asset", "price")
  stopifnot(all(needed %in% names(x)))

  x$date <- as.Date(x$date)
  stopifnot(!anyNA(x$date), !anyNA(x$asset))
  stopifnot(all(is.finite(x$price)), all(x$price > 0))

  key <- paste(x$asset, x$date)
  stopifnot(!anyDuplicated(key))

  ordered <- x[order(x$asset, x$date), ]
  increasing <- vapply(
    split(ordered$date, ordered$asset),
    function(z) length(z) < 2L || all(diff(z) > 0),
    logical(1)
  )
  stopifnot(all(increasing))
  invisible(ordered)
}

toy_prices <- data.frame(
  date = as.Date("2026-01-01") + 0:3,
  asset = "教學資產",
  price = c(100, 102, 101, 104)
)
validate_price_data(toy_prices)
```

故意建立重複日期時，函數會停止。為了讓整份附錄仍能執行，這裡捕捉預期錯誤。


``` r
bad_prices <- rbind(toy_prices, toy_prices[2, ])
bad_result <- try(validate_price_data(bad_prices), silent = TRUE)
inherits(bad_result, "try-error")
```

```
## [1] TRUE
```

## 記錄環境

網站版輸出應保留 R 版本、作業系統與主要套件版本。若結果日後不同，這些資訊是追查起點。


``` r
environment_record <- list(
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  R = R.version.string,
  platform = R.version$platform,
  knitr = as.character(utils::packageVersion("knitr")),
  rmarkdown = as.character(utils::packageVersion("rmarkdown")),
  rng_kind = RNGkind()
)
environment_record
```

```
## $generated_at
## [1] "2026-07-16 08:35:58 UTC"
## 
## $R
## [1] "R version 4.5.2 (2025-10-31)"
## 
## $platform
## [1] "aarch64-apple-darwin20"
## 
## $knitr
## [1] "1.51"
## 
## $rmarkdown
## [1] "2.31"
## 
## $rng_kind
## [1] "Mersenne-Twister" "Inversion"        "Rejection"
```

## 可重現研究檢核表

1. 原始或凍結資料來源、授權與建立方式有書面說明。
2. 分析只讀固定檔，路徑由專案結構推導。
3. 日期、主鍵、單位、缺值與排序在入口檢查。
4. 模擬固定種子，正式結果記錄 R 與套件版本。
5. 預測的訓練、驗證、測試期間依時間分開。
6. 中間結果可由程式重建，不靠手動複製貼上。
7. 錯誤與警告被理解後處理，不以全域關閉訊息掩蓋。

本附錄只驗證「相同輸入與環境能否重算」，不等同驗證研究問題、模型假設或因果解讀正確。
