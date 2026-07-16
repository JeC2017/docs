# 教科書實證資料：來源、版本與公開政策

本書的實證主線回到作者四個學期密集課程實際使用的資料與程式。`data/processed/` 內的十個 CSV 是可直接執行線上附錄的固定版本；`manifest.csv` 記錄列數、欄數與 MD5，避免下載、排序或清理差異被誤認為估計差異。

## 公開範圍

作者已明確同意公開本課程用於教學與研究重現的資料檔，因此 GitHub Pages 版會提供十個實證 CSV、資料字典、來源說明、R Markdown、執行後結果與圖表。這項同意不包含學生名冊、學生作業、評分檔、合約、行政資料或任何個人資料；這些檔案不會進入專案或網站。

部分資料源自 FRED、Kenneth French Data Library、global-q、Welch--Goyal、Ruey S. Tsay 教科書網站或課程整理檔。本站保留來源識別與轉換紀錄；下載者仍須遵守各原始供應者的條款。公開固定版本的目的，是讓讀者核對本書結果，不是替第三方授權作一般性保證。

## 原始來源鏈（查核日：2026-07-16）

- **FRED 匯率：**Board of Governors of the Federal Reserve System (US), H.10 Foreign Exchange Rates；[DEXJPUS](https://fred.stlouisfed.org/series/DEXJPUS) 與 [DEXTAUS](https://fred.stlouisfed.org/series/DEXTAUS)。官方頁定義兩者分別為一美元兌日圓與新臺幣數、日資料、未季調，並標示「Public Domain: Citation Requested」。本書保存 2020-01-01 至 2022-12-16 的 FRED CSV 回傳值。
- **Tsay 教科書資料：**Ruey S. Tsay, *Analysis of Financial Time Series*, 3rd ed. 的[官方 companion page](https://faculty.chicagobooth.edu/ruey-s-tsay/research/analysis-of-financial-time-series-3rd-edition)及原課程程式內的三個固定檔 URL：`d-msft8608.txt`、`m-5clog-9008.txt`、`m-barra-9003.txt`。
- **Fama--French 與十產業：**Kenneth R. French [Data Library](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html) 的 Fama/French 3 Factors 與 10 Industry Portfolios。現行網站說明 2025 年起 CRSP 來源格式由 FIZ 改為 CIZ，因此本書鎖定原課程截至 2021-11 的工作快照，不把今日重抓值當成同一版本。
- **global-q：**原課程 `fffqmacro.R` 指向 2021 版 q5 monthly CSV；因子定義與百分比單位可由 [global-q technical document](https://global-q.org/uploads/1/2/2/6/122679606/factorstd_2025feb.pdf) 查核。本書使用 ME、IA、ROE、EG 四欄，並於合併時由百分點除以 100。
- **Welch--Goyal 總體預測變數：**原課程依 [Tidy Finance 建檔流程](https://github.com/ramnathv/tidy-finance-website/blob/main/accessing-and-managing-financial-data.qmd) 讀取 Amit Goyal 公開工作簿的 Monthly sheet，再自行建立 `dp`、`dy`、`ep`、`de`、`tms` 與 `dfy` 等欄。本書鎖定原課程 2021 工作版，不宣稱是 2026 最新 vintage。
- **California schools：**原課程依 Stock and Watson, *Introduction to Econometrics*, 4th ed., Chapter 14 的學校／學區資料與變數表整理；原作者的[第四版資源頁](https://www.princeton.edu/~mwatson/Stock-Watson_4E/Stock-Watson-Resources-4e.html)提供課程資料下載。本書只用於 PCA 與高維控制示範，不把機構名稱或郵遞區號當成個人資料。
- **日本月資料、S&P 價格長表與臺灣 ICAPM：**保留自作者原課程資料夾；其中日本檔未附可核對的 provider/vintage，S&P 長表未附 vendor 授權文字，臺灣案例只保留匿名第二階段。這三組資料在正文均明示來源紀錄缺口與可支持的有限結論。

## 固定資料與用途

| 固定檔 | 樣本與用途 | 課程來源／處理重點 |
|---|---|---|
| `aapl_adjusted_daily_2019_2022.csv` | 875 個交易日；價格、報酬、ARMA、GARCH／GJR | 原 S&P 500 價格檔中的 AAPL，2019-01-02 至 2022-06-22；報酬依日期排序後計算 |
| `msft_daily_returns_1986_2008.csv` | 5,752 筆日簡單報酬；描述統計、厚尾與常態性 | Tsay 教科書網站 `d-msft8608.txt`；原課程第 2 講 |
| `fred_jpy_twd_daily_2020_2022.csv` | JPY/USD、TWD/USD、TWD/JPY 與交叉匯率對數報酬；ARMA、單根、預測評估 | FRED `DEXJPUS`、`DEXTAUS`，2020-01-02 至 2022-12-16；交叉匯率定義為 TWD per JPY |
| `tsay_five_stock_monthly_returns_1990_2008.csv` | IBM、HPQ、INTC、JPM、BAC 月對數報酬；PCA | Tsay 範例 9.2，`m-5clog-9008.txt` |
| `tsay_barra_monthly_returns_1990_2003.csv` | 十家公司月報酬；統計因子分析 | Tsay 範例 9.4，`m-barra-9003.txt` |
| `california_schools.csv` | 3,932 所學校、110 欄；高維 PCA、正則化、雙重選擇 | 原課程 `str_pca` hands-on；欄位定義見 `DATA_DICTIONARY.md` 與課程變數表 |
| `sp500_returns_balanced_2013_2022.csv` | 2,384 日、89 檔股票；PCA、共同因子與樣本外重建 | 原 S&P 500 長表；先在每一股票內取落後值，再取共同交易日，修正舊程式跨股票 `lag()` 問題 |
| `japan_monthly_2007_2018.csv` | 133 月、30 欄；LASSO、VAR、共整合與局部投影 | 原課程日本總體金融資料與十年期公債殖利率；日期轉為 ISO 並依時間排序 |
| `ff_qf_macro_industries_1967_2021.csv` | 659 月、十產業，合計 6,590 列；資產定價與金融因子選擇 | 原課程合併 Fama--French、global-q、Welch--Goyal 總體預測變數與十產業組合 |
| `taiwan_icapm_second_stage_47x497.csv` | 47 家公司、497 日；固定效果與條件共變數的第二階段案例 | 原課程 `datICAPM3.R`；只重現可核對的第二階段，不宣稱補出缺少的第一階段 DCC/PQR 原始資料與估計流程 |

## 可重建流程

從教科書專案根目錄執行：

```r
Rscript code/build_empirical_data.R
```

程式會直接讀取相鄰課程資料夾中的原課程檔，以及 `data/source_snapshots/` 中按原課程 URL 保存的 FRED／Tsay 快照，檢查必要欄位與預期維度，再重建十個 CSV 和 `manifest.csv`。公開網站則直接附上固定 CSV，因此一般讀者不必擁有作者的 Dropbox 路徑。

## 不能越過的實證界線

- 固定檔可重現「這個版本」的數值，不自動證明資料定義或識別假設正確。
- 臺灣 ICAPM 檔只有第二階段共變數與報酬；正文與 R appendix 均須明寫「部分重現」。
- FRED 交叉匯率以報價單位命名，避免把 `TWD per JPY` 誤寫成反方向。
- 資產報酬、總體變數與學校資料的迴歸，除非另有可信識別設計，不作因果宣稱。
