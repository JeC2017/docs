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
- 十個公開實證資料位於 data/processed，建置紀錄位於 manifest.csv；網站的「實證資料下載」頁提供相同版本。

## 兩條 R 學習路線

每份相關附錄保留兩種互補做法：

1. **手動建構版**把估計式、遞迴、損失函數與資訊邊界寫出來，適合用來理解方法與除錯；
2. **原課程套件捷徑**移植作者以前授課 R scripts 中實際使用的高階函數，例如 `forecast::auto.arima()`、`forecast::tsCV()`、`tseries::adf.test()`、`fGarch::garchFit()`、`stats::prcomp()` 與 `glmnet::glmnet()`。

套件版不重新從 Yahoo、FRED 或其他網站下載資料，而是讀取本書公開的固定 CSV。這項調整保留原課程的估計工作流，同時避免資料供應者日後修訂、交易日不同或網路中斷造成答案漂移。兩版結果若不同，附錄會交代模型集合、估計法、初始化或有限樣本修正的差異，不把「套件輸出」當作無條件正解。

### 一次安裝必要套件

安裝只在環境設定階段執行一次，不應放進正式分析 chunk。以下程式因此不會在製作網站時自動執行；學生可在 R console 執行。


``` r
install.packages(c(
  "knitr", "rmarkdown", "ragg", "systemfonts",
  "tidyverse", "fBasics", "forecast", "tseries",
  "urca", "fGarch", "plm", "glmnet", "pls"
))
```

`phtt` 與完整 `tidymodels` 工作流只在相應進階段落使用；若該段明列需要，再依附錄指示另行安裝。正式 Rmd 只檢查套件是否存在，不在執行中偷偷改變環境。


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
  "data/DATA_DICTIONARY.md",
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
## 3     data/DATA_DICTIONARY.md   TRUE
## 4      chapters/chapter01.tex   TRUE
## 5   online_appendix/README.md   TRUE
```

``` r
stopifnot(all(path_check$exists))
```

## 用 manifest 核對固定資料

manifest.csv 記錄十個公開實證檔的檔名、列數、欄數、MD5、內容說明與建立時間。MD5 不能證明資料正確，但能確認讀者拿到的位元內容是否與教科書版本相同。


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
##                                                            file  rows columns
## 1              data/processed/aapl_adjusted_daily_2019_2022.csv   875       7
## 2               data/processed/msft_daily_returns_1986_2008.csv  5752       2
## 3               data/processed/fred_jpy_twd_daily_2020_2022.csv   772       5
## 4  data/processed/tsay_five_stock_monthly_returns_1990_2008.csv   228       6
## 5       data/processed/tsay_barra_monthly_returns_1990_2003.csv   168      11
## 6                         data/processed/california_schools.csv  3932     110
## 7           data/processed/sp500_returns_balanced_2013_2022.csv  2384      90
## 8                    data/processed/japan_monthly_2007_2018.csv   133      30
## 9           data/processed/ff_qf_macro_industries_1967_2021.csv  6590      24
## 10          data/processed/taiwan_icapm_second_stage_47x497.csv 23359       6
##                                 md5
## 1  8205cf538da8ea57fd8c93264861c28d
## 2  d2f6b2d5d74be413385347a17c82b754
## 3  e50c18c906c6bbba9ddbfd10fd735080
## 4  8fff77a351970b335333a7fe88fba8b7
## 5  85eb8ec635cd19818ffa00d73acc5d39
## 6  28b3dff5db50448608925cad32feb18a
## 7  09c9690effb82b3fabdccaa982397e83
## 8  46b39f6fdde5d581ad31c83348d99933
## 9  69563611584d8a2dfd984ec6a53822a4
## 10 98fb791d16ee3b3e536ef0ce33381e93
##                                                                       description
## 1             AAPL adjusted prices and returns from the course S&P 500 price file
## 2            Microsoft daily simple returns used in the return-properties lecture
## 3         FRED JPY/USD and TWD/USD daily rates and the derived TWD/JPY cross rate
## 4                   Five-company monthly log returns used in Tsay PCA example 9.2
## 5            Ten-company monthly returns used in Tsay factor-analysis example 9.4
## 6    California school and district data used in the PCA/high-dimensional lecture
## 7  Balanced daily returns for 89 S&P 500 constituents; lag computed within symbol
## 8              Japanese monthly macro-finance panel with 10-year government yield
## 9                           Fama-French/global-q/macro/ten-industry monthly panel
## 10   Prepared Taiwan ICAPM second-stage panel; not a full first-stage replication
##                   built_at                       actual_md5 md5_match
## 1  2026-07-16 09:57:40 UTC 8205cf538da8ea57fd8c93264861c28d      TRUE
## 2  2026-07-16 09:57:40 UTC d2f6b2d5d74be413385347a17c82b754      TRUE
## 3  2026-07-16 09:57:40 UTC e50c18c906c6bbba9ddbfd10fd735080      TRUE
## 4  2026-07-16 09:57:40 UTC 8fff77a351970b335333a7fe88fba8b7      TRUE
## 5  2026-07-16 09:57:40 UTC 85eb8ec635cd19818ffa00d73acc5d39      TRUE
## 6  2026-07-16 09:57:40 UTC 28b3dff5db50448608925cad32feb18a      TRUE
## 7  2026-07-16 09:57:40 UTC 09c9690effb82b3fabdccaa982397e83      TRUE
## 8  2026-07-16 09:57:40 UTC 46b39f6fdde5d581ad31c83348d99933      TRUE
## 9  2026-07-16 09:57:40 UTC 69563611584d8a2dfd984ec6a53822a4      TRUE
## 10 2026-07-16 09:57:40 UTC 98fb791d16ee3b3e536ef0ce33381e93      TRUE
##    file_exists
## 1         TRUE
## 2         TRUE
## 3         TRUE
## 4         TRUE
## 5         TRUE
## 6         TRUE
## 7         TRUE
## 8         TRUE
## 9         TRUE
## 10        TRUE
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
##                                                            file  rows columns
## 1              data/processed/aapl_adjusted_daily_2019_2022.csv   875       7
## 2               data/processed/msft_daily_returns_1986_2008.csv  5752       2
## 3               data/processed/fred_jpy_twd_daily_2020_2022.csv   772       5
## 4  data/processed/tsay_five_stock_monthly_returns_1990_2008.csv   228       6
## 5       data/processed/tsay_barra_monthly_returns_1990_2003.csv   168      11
## 6                         data/processed/california_schools.csv  3932     110
## 7           data/processed/sp500_returns_balanced_2013_2022.csv  2384      90
## 8                    data/processed/japan_monthly_2007_2018.csv   133      30
## 9           data/processed/ff_qf_macro_industries_1967_2021.csv  6590      24
## 10          data/processed/taiwan_icapm_second_stage_47x497.csv 23359       6
##    row_match column_match
## 1       TRUE         TRUE
## 2       TRUE         TRUE
## 3       TRUE         TRUE
## 4       TRUE         TRUE
## 5       TRUE         TRUE
## 6       TRUE         TRUE
## 7       TRUE         TRUE
## 8       TRUE         TRUE
## 9       TRUE         TRUE
## 10      TRUE         TRUE
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
teaching_packages <- c(
  "knitr", "rmarkdown", "ragg", "systemfonts",
  "tidyverse", "fBasics", "forecast", "tseries",
  "urca", "fGarch", "plm", "glmnet", "pls"
)
package_record <- data.frame(
  package = teaching_packages,
  installed = vapply(
    teaching_packages, requireNamespace,
    logical(1), quietly = TRUE
  ),
  version = vapply(teaching_packages, function(package) {
    if (requireNamespace(package, quietly = TRUE)) {
      as.character(utils::packageVersion(package))
    } else {
      NA_character_
    }
  }, character(1)),
  row.names = NULL
)
package_record
```

```
##        package installed version
## 1        knitr      TRUE    1.51
## 2    rmarkdown      TRUE    2.31
## 3         ragg      TRUE   1.5.2
## 4  systemfonts      TRUE   1.3.2
## 5    tidyverse      TRUE   2.0.0
## 6      fBasics     FALSE    <NA>
## 7     forecast     FALSE    <NA>
## 8      tseries     FALSE    <NA>
## 9         urca     FALSE    <NA>
## 10      fGarch     FALSE    <NA>
## 11         plm     FALSE    <NA>
## 12      glmnet      TRUE  4.1.10
## 13         pls      TRUE   2.9.0
```

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
## [1] "2026-07-16 13:47:51 UTC"
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

1. 原始或固定資料的來源、公開範圍與建立方式有書面說明。
2. 分析只讀固定檔，路徑由專案結構推導。
3. 日期、主鍵、單位、缺值與排序在入口檢查。
4. 模擬固定種子，正式結果記錄 R 與套件版本。
5. 預測的訓練、驗證、測試期間依時間分開。
6. 中間結果可由程式重建，不靠手動複製貼上。
7. 錯誤與警告被理解後處理，不以全域關閉訊息掩蓋。

本附錄只驗證「相同輸入與環境能否重算」，不等同驗證研究問題、模型假設或因果解讀正確。
