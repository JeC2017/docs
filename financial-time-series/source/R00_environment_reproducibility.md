---
title: "R00：環境、專案與可重現研究"
output:
  github_document:
    toc: true
    toc_depth: 3
---

本附錄對應第 1 章。先想像一個很實際的情境：半年後重新打開這個專案時，我們能不能確定讀到的是同一份資料，並在另一部電腦上算出相同結果？要回答這個問題，光是保留 R 程式還不夠；資料版本、檔案路徑、亂數種子與套件版本都要一併留下。

這一頁不估計金融模型，而是先把後續實證分析的工作環境整理好。我們會確認固定資料是否齊全、使用相對路徑讀檔、記錄執行環境，並讓日期或價格問題在資料入口就顯示出來。全文不連線下載資料，也不在執行期間安裝或更新套件，因此不會因當天的網路狀況或套件版本變動而改寫分析環境。

## 執行前提

- R 4.1 以上，並已具備 knitr 與 rmarkdown。
- 從教科書專案根目錄或 online_appendix 目錄執行皆可。
- 十份固定實證資料位於 `data/processed`，建置紀錄位於 `manifest.csv`；網站的「實證資料下載」頁提供相同版本。

本頁的觀察單位是「一個資料檔」，不是一天或一檔股票；`manifest.csv` 的每一列描述一份檔案的大小與指紋。由於本頁尚未進行預測，沒有訓練期、驗證期與測試期之分。這三種期間的角色會在 R05 與 R06 中正式使用。

## 本書如何安排手動與套件作法

每份相關附錄保留兩種互補做法：

1. **手動建構版**把估計式、遞迴、損失函數與每個時點可使用的資料寫出來，適合用來理解方法與除錯；
2. **套件作法**沿用原課程 R 程式中實際使用的高階函數，例如 `forecast::auto.arima()`、`forecast::tsCV()`、`tseries::adf.test()`、`fGarch::garchFit()`、`stats::prcomp()`、`glmnet::glmnet()` 與 `quantreg::rq()`。

套件作法不重新從 Yahoo、FRED 或其他網站下載資料，而是讀取本書提供的固定 CSV。這樣仍可學到原課程的估計流程，也不會因供應者日後修訂歷史資料、交易日對齊方式改變或網路中斷而得到另一個答案。高階函數會替我們完成許多計算，卻不會替研究者決定資料版本、模型集合或資訊可得時點；兩種作法若出現差異，後文會從估計法、初始化與有限樣本修正逐項找原因。

### 一次安裝必要套件

安裝只需在環境設定階段做一次，不宜放進正式分析的程式區塊。以下程式因此不會隨網站頁面自動執行；第一次使用本專案時，可將它貼到 R 主控台執行。


``` r
install.packages(c(
  "knitr", "rmarkdown", "ragg", "systemfonts",
  "tidyverse", "fBasics", "forecast", "tseries",
  "urca", "fGarch", "plm", "glmnet", "quantreg", "pls"
))
```

`phtt` 與完整的 `tidymodels` 工作流程只出現在相應的進階段落；讀到該段時再依指示安裝即可。正式 Rmd 只檢查套件是否存在，不會在分析途中自動變更環境。

接下來的 `setup` 區塊先設定全頁共用的圖形選項，再從目前目錄與上一層目錄尋找專案根目錄。這一步讓同一份 Rmd 從兩個常見位置執行時，都會指向同一批資料；若找不到 `main.tex` 與 `manifest.csv`，程式會立刻停止，避免悄悄讀到其他同名檔案。


``` r
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.width = 7,
  fig.height = 4.5
)

root_candidates <- c(".", "..")
# 專案根目錄必須同時含主檔與資料清冊，單靠資料夾名稱不夠可靠。
is_root <- vapply(
  root_candidates,
  function(x) file.exists(file.path(x, "main.tex")) &&
    file.exists(file.path(x, "data", "processed", "manifest.csv")),
  logical(1)
)
stopifnot(any(is_root))
project_root <- root_candidates[which(is_root)[1]]
# 後續路徑都從同一個根目錄組成，搬動整個專案時不必逐行修改。
project_path <- function(...) file.path(project_root, ...)
```

## 先讓每一類檔案各就各位

一份實證專案通常同時包含外部輸入、整理後資料、分析程式、圖表與正文。把它們分開存放，日後看到一個檔案時，才知道它是原始輸入、可由程式重建的產出，還是需要人工撰寫的文字。下列程式先檢查後續附錄會用到的幾個關鍵路徑；只要缺少其中一項，就先補齊專案，而不要讓分析執行到一半才失敗。


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
# 缺少任何必要檔案時立即停止，避免後續錯把空白或舊檔當成輸入。
stopifnot(all(path_check$exists))
```

## 用 `manifest.csv` 確認拿到的是同一批資料

`manifest.csv` 記錄十份實證檔的檔名、列數、欄數、MD5、內容說明與建立時間。這裡要回答的是「目前電腦上的檔案，是否與本書製作時使用的版本完全相同」。MD5 相同表示檔案的位元內容一致；它不會告訴我們資料來源是否可靠，也不會替代經濟意義與建檔方法的檢查。


``` r
manifest <- read.csv(
  project_path("data", "processed", "manifest.csv"),
  stringsAsFactors = FALSE
)

actual_md5 <- unname(tools::md5sum(project_path(manifest$file)))
# 同時保留宣告值與實際值；不一致時才看得出問題出在哪一個檔案。
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

輸出的 `md5_match` 若全為 `TRUE`，表示目前的十份檔案與清冊版本相同。接著再比較列數與欄數，因為檔案若在試算表中被截短或多存了一欄，資料形狀會立刻改變。這一步只讀檔，不修改任何資料。


``` r
dimensions <- lapply(manifest$file, function(f) {
  # 每次重新讀取一份固定輸入，直接量測實際列數與欄數。
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

`row_match` 與 `column_match` 全為 `TRUE` 時，代表檔案形狀也符合清冊。若程式在這裡停止，合理的下一步是重新取得固定資料或查明哪個檔案曾被修改，而不是直接更新清冊來配合手上的檔案。

## 亂數種子與可重現模擬

模擬結果會隨亂數改變，所以在比較兩段程式以前，先固定亂數種子。相同種子配合相同的 R 版本與亂數產生器設定，應得到同一串模擬值。種子的用途是讓別人重做計算，不是反覆嘗試直到結果看起來理想；若擔心單一次抽樣太偶然，可以事前指定多個種子或提高模擬次數。


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

五列 `exactly_equal` 都應是 `TRUE`。這項結果只確認亂數流程可以重做，不表示模擬設計或統計結論已經正確。

## 讓錯誤提早發生

真正開始估計以前，先問一個更基本的問題：每一列資料是否代表清楚而且唯一的「資產—日期」觀察值？以下函數示範價格資料最基本的入口條件：日期可轉換、同一資產同一天不重複，價格為正且有限。函數最後依資產與日期排序，讓後續的落後值確實對應上一個觀察日。


``` r
validate_price_data <- function(x) {
  # 先確認欄位，否則後面的日期與價格檢查可能針對錯誤物件。
  needed <- c("date", "asset", "price")
  stopifnot(all(needed %in% names(x)))

  x$date <- as.Date(x$date)
  stopifnot(!anyNA(x$date), !anyNA(x$asset))
  stopifnot(all(is.finite(x$price)), all(x$price > 0))

  key <- paste(x$asset, x$date)
  # 同一資產同一天只能有一筆價格，否則報酬的上一期無法唯一決定。
  stopifnot(!anyDuplicated(key))

  # 先排序再檢查時間方向，避免原始列順序誤導 lag 或 diff。
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

`toy_prices` 的觀察單位是一個資產在一個日曆日的價格，價格尺度是任意的教學單位。函數沒有填補休市日，也沒有計算報酬；它只確定資料已具備進入下一步的條件。

接著故意放入一筆重複日期。正常分析遇到這種情況應停止並回查來源；為了讓整份附錄仍可繼續執行，這裡用 `try()` 留下錯誤訊號。


``` r
bad_prices <- rbind(toy_prices, toy_prices[2, ])
bad_result <- try(validate_price_data(bad_prices), silent = TRUE)
inherits(bad_result, "try-error")
```

```
## [1] TRUE
```

輸出為 `TRUE`，表示重複鍵確實被攔下。這類入口檢查應放在計算報酬以前；若先排序、取落後值再發現重複資料，錯誤可能已經傳到所有後續結果。

## 記錄環境

最後記下 R 版本、作業系統、亂數產生器與主要套件版本。若同一份資料與程式日後出現不同結果，這些資訊可以幫助我們分辨差異來自資料、程式，還是執行環境。


``` r
teaching_packages <- c(
  "knitr", "rmarkdown", "ragg", "systemfonts",
  "tidyverse", "fBasics", "forecast", "tseries",
  "urca", "fGarch", "plm", "glmnet", "quantreg", "pls"
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
## 6      fBasics      TRUE 4052.98
## 7     forecast      TRUE   9.0.2
## 8      tseries      TRUE 0.10.62
## 9         urca      TRUE   1.3.4
## 10      fGarch      TRUE 4052.93
## 11         plm      TRUE   2.6.7
## 12      glmnet      TRUE  4.1.10
## 13    quantreg      TRUE     6.1
## 14         pls      TRUE   2.9.0
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
## [1] "2026-07-16 21:52:07 UTC"
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

`package_record` 會把尚未安裝的套件標成 `FALSE`，並將版本留為缺值；這時只需安裝後續附錄真正會用到的套件。`environment_record` 則是本次執行的時間戳記與環境摘要，適合與分析結果一起保存。

## 把這些習慣帶到後續分析

到這裡，我們確認的是：相同的固定輸入、程式與執行環境，應能重算出相同結果。MD5、檔案形狀與亂數種子都通過，仍不代表研究問題選得恰當，也不代表模型假設或因果解讀成立；它們處理的是「能不能重做」，而不是「結論對不對」。

後續每次加入新資料時，先寫清楚觀察單位、期間、尺度、缺值與排序方式，再決定哪些資訊在預測起點已經可得。若要評估預測，訓練期、驗證期與測試期也要依時間分開。中間表格與圖形盡量由程式產生，遇到警告則先找出原因，再決定是否需要調整資料或模型。這樣即使分析結果不如預期，讀者仍能看懂它是如何得到的，以及下一步該從哪裡查起。
