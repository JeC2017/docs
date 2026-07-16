# 教科書實證資料：來源與固定版本

本書的實證分析沿用「財務時間序列分析」課程實際使用的資料，並整理成十個可直接執行線上附錄的 CSV。檔案放在 `data/processed/`；同一目錄下的 `manifest.csv` 記錄每份檔案的列數、欄數與 MD5。如果讀者算出的結果與書中不同，可以先從這三項資訊確認是否讀到同一版本，再檢查排序、缺值處理與模型設定。

## 網站隨附哪些資料？

網站隨附十個實證 CSV、資料字典、來源說明、R Markdown、執行結果（Markdown）與圖表，讓讀者可以重跑書中的分析並比較輸出。下載區只包含本書使用的教材與固定資料，不含原課程的學生或行政紀錄。

部分欄位來自 FRED、Kenneth French Data Library、global-q、Welch–Goyal 與 Ruey S. Tsay 教科書網站，其餘則保留自原課程工作檔。下文整理目前可以確認的來源與轉換方式。固定版本說明了書中數值如何得到，原始供應者的變數定義與使用條款仍以來源網站的最新說明為準。若資料沒有留下完整的來源版本資訊，本書只解讀現存固定快照可以支持的結果，不主張能從目前的上游資料完整重建。

## 資料來源與版本（查閱日：2026-07-16）

- **FRED 匯率：**資料來自美國聯邦準備理事會（Board of Governors of the Federal Reserve System）的 H.10 Foreign Exchange Rates，包括 [DEXJPUS](https://fred.stlouisfed.org/series/DEXJPUS) 與 [DEXTAUS](https://fred.stlouisfed.org/series/DEXTAUS)。官方頁將兩者分別定義為一美元可兌換的日圓數與新臺幣數；資料頻率為日資料，未經季節調整，並標示「Public Domain: Citation Requested」。本書保存 2020-01-02 至 2022-12-16 的 FRED CSV 回傳值。
- **Tsay 教科書資料：**Ruey S. Tsay, *Analysis of Financial Time Series*, 3rd ed. 的 [官方配套資料頁](https://faculty.chicagobooth.edu/ruey-s-tsay/research/analysis-of-financial-time-series-3rd-edition)，以及原課程程式指定的三個檔案：`d-msft8608.txt`、`m-5clog-9008.txt` 與 `m-barra-9003.txt`。
- **Fama–French 與十產業：**資料取自 Kenneth R. French [Data Library](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html) 的 Fama/French 3 Factors 與 10 Industry Portfolios。該網站說明 CRSP 的月報酬資料自 2025 年起由 FIZ 改為 CIZ；這次變更牽涉資料格式與月報酬累積方式，因此新版序列可能和舊快照不同。本書固定使用原課程截至 2021-11 的工作快照，不把今日重新下載的值視為同一版本。現有快照沒有保留每月因子的發布日，也不足以用當時的成分股與權重重新建構投資組合。
- **global-q：**原課程的 `fffqmacro.R` 指向 2021 版 q5 monthly CSV；因子定義與百分比單位可參考 [global-q technical document](https://global-q.org/uploads/1/2/2/6/122679606/factorstd_2025feb.pdf)。本書使用 ME、IA、ROE 與 EG 四欄，合併資料時再將百分點除以 100，改成小數尺度。
- **Welch–Goyal 總體預測變數：**原課程依 [Tidy Finance 建檔流程](https://github.com/ramnathv/tidy-finance-website/blob/main/accessing-and-managing-financial-data.qmd) 讀取 Amit Goyal 公開工作簿的 Monthly 工作表，再建立 `dp`、`dy`、`ep`、`de`、`tms` 與 `dfy` 等欄。本書固定使用原課程 2021 年的工作版本；它不是 2026 年重新下載的最新資料版本。
- **加州學校資料：**原課程依 Stock and Watson, *Introduction to Econometrics*, 4th ed., Chapter 14 的學校與學區資料及變數表整理；原作者的 [第四版資源頁](https://www.princeton.edu/~mwatson/Stock-Watson_4E/Stock-Watson-Resources-4e.html) 提供課程資料下載。本書用這份資料示範 PCA、高維控制、雙重選擇與 DML；學校或學區識別碼不作為連續解釋變數。
- **日本月資料、S&P 價格長表與臺灣 ICAPM：**三者均保留自原課程資料夾。日本檔沒有留下可確認的原始供應者與資料版本；S&P 價格長表只保留課程工作檔，沒有完整的版本紀錄；臺灣案例則只保留匿名化的第二階段資料。相關附錄會分別說明這些缺口，並把結論限於固定資料版本實際能支持的範圍。

## 固定資料與用途

| 固定檔 | 樣本與用途 | 課程來源／處理重點 |
|---|---|---|
| `aapl_adjusted_daily_2019_2022.csv` | 875 個交易日；價格、報酬、ARMA、GARCH／GJR | 原 S&P 500 價格檔中的 AAPL，2019-01-02 至 2022-06-22；報酬依日期排序後計算 |
| `msft_daily_returns_1986_2008.csv` | 5,752 筆日簡單報酬；描述統計、厚尾與常態性 | Tsay 教科書網站 `d-msft8608.txt`；原課程第 2 講 |
| `fred_jpy_twd_daily_2020_2022.csv` | JPY/USD、TWD/USD、TWD/JPY 與交叉匯率對數報酬；ARMA、單根、預測評估 | FRED `DEXJPUS`、`DEXTAUS`，2020-01-02 至 2022-12-16；交叉匯率定義為 TWD per JPY |
| `tsay_five_stock_monthly_returns_1990_2008.csv` | IBM、HPQ、INTC、JPM、BAC 月對數報酬；PCA | Tsay 範例 9.2，`m-5clog-9008.txt` |
| `tsay_barra_monthly_returns_1990_2003.csv` | 十家公司月報酬；統計因子分析 | Tsay 範例 9.4，`m-barra-9003.txt` |
| `california_schools.csv` | 3,932 所學校、110 欄；高維 PCA、正則化、雙重選擇 | 原課程 `str_pca` 實作；欄位定義見 `DATA_DICTIONARY.md` 與課程變數表 |
| `sp500_returns_balanced_2013_2022.csv` | 2,384 日、89 檔股票；PCA、共同因子與樣本外重建 | 原 S&P 500 長表；先在每一股票內取落後值，再取共同交易日，修正舊程式跨股票 `lag()` 問題 |
| `japan_monthly_2007_2018.csv` | 133 月、30 欄；LASSO、VAR、共整合與局部投影 | 原課程日本總體金融資料與十年期公債殖利率；日期轉為 ISO 並依時間排序 |
| `ff_qf_macro_industries_1967_2021.csv` | 659 月、十產業，合計 6,590 列；資產定價與金融因子選擇 | 原課程合併 Fama–French、global-q、Welch–Goyal 總體預測變數與十產業組合 |
| `taiwan_icapm_second_stage_47x497.csv` | 47 家公司、497 日；固定效果與條件共變數的第二階段案例 | 原課程 `datICAPM3.R`；目前檔案可重算第二階段，第一階段 DCC/PQR 則缺少原始資料與估計流程 |

## 如何重建固定資料？

從教科書專案根目錄執行：

```r
Rscript code/build_empirical_data.R
```

這支程式會讀取原課程資料檔，以及 `data/source_snapshots/` 中保存的 FRED／Tsay 下載快照；確認必要欄位與資料維度後，再建立十個 CSV 和 `manifest.csv`。完整重建仍需要原課程工作目錄中的來源檔，一般讀者可以直接使用網站隨附的固定 CSV 執行附錄，不必另行取得原始工作路徑。

## 解讀結果前要先掌握的限制

- 固定檔讓讀者重算「這個版本」的數值；變數定義是否適合研究問題、識別假設是否成立，仍須另外判斷。
- 臺灣 ICAPM 檔只有第二階段共變數與報酬，因此相關附錄能重算的是第二階段分析。
- FRED 交叉匯率依報價單位命名；`TWD per JPY` 表示一日圓可兌換的新臺幣數，解讀時不可倒置方向。
- 資產報酬、總體變數與學校資料的迴歸，在沒有額外識別設計時，應解讀為條件關聯或預測關係。
- R14 可以確認因子欄的月份標籤與固定快照，卻無法得知每個歷史時點的實際發布日與當時可得版本。因此，該頁是在明示可得性假設下進行擬樣本外評估，不能直接視為即時交易回測。
