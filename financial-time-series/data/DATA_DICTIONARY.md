# 實證資料字典

本頁說明線上附錄實際使用的欄位、單位與時間索引。讀取一份新資料時，請先回答三個問題：每一列代表哪個觀察單位？報酬以小數還是百分點表示？日期指交易日、FRED 觀察日，還是月份？完整的列數、欄數與 MD5 記錄在 `processed/manifest.csv`，可以用來確認目前讀入的是哪一個固定版本。

## AAPL 與 MSFT

- `aapl_adjusted_daily_2019_2022.csv`
  - `date`：交易日。
  - `adjusted`：原價格檔的調整後收盤價，採美元／股的調整後價格尺度。
  - `simple_return`：$P_t/P_{t-1}-1$。
  - `log_return`：$\log P_t-\log P_{t-1}$。
  - `symbol`、`company`、`sector`：證券與產業識別欄。
- `msft_daily_returns_1986_2008.csv`
  - `date`：交易日。
  - `simple_return`：Microsoft 日簡單報酬，以小數表示；例如 `0.01` 代表 1%。

兩份檔案的每一列都代表一個交易日。AAPL 第一列沒有前一期價格，因此兩種報酬欄自然會是缺值；分析報酬時應保留這項時間關係，而不是把缺值填成 0。

## FRED 日匯率

`fred_jpy_twd_daily_2020_2022.csv`：

- `date`：FRED 工作日觀察日期；來源序列在該日沒有報價時保留為空值，R 讀入後記為 `NA`。
- `jpy_per_usd`：一美元可兌日圓數，FRED `DEXJPUS`。
- `twd_per_usd`：一美元可兌新臺幣數，FRED `DEXTAUS`。
- `twd_per_jpy`：`twd_per_usd / jpy_per_usd`。
- `log_return_twd_per_jpy`：交叉匯率的一期對數差分；遇缺值時不跨缺值補算。

每一列對應一個 FRED 工作日觀察日期。只有兩條來源匯率都有數值時，`twd_per_jpy` 才能形成；若其中一天缺值，`log_return_twd_per_jpy` 也不跨過該日計算。R06 篩選有限報酬後，預測的「下一期」是下一個有效報酬觀察日，不一定是下一個日曆日。

## Tsay 股票面板

- `tsay_five_stock_monthly_returns_1990_2008.csv`：每一列是一個月份，欄位為 `month` 與 IBM、HPQ、INTC、JPM、BAC 五家公司月對數報酬。
- `tsay_barra_monthly_returns_1990_2003.csv`：每一列是一個月份，欄位為 `month` 與 AGE、C、MWD、MER、DELL、HPQ、IBM、AA、CAT、PG 十家公司月報酬。原檔數值以百分點表示；程式在需要時除以 100。

## S&P 500 平衡報酬面板

`sp500_returns_balanced_2013_2022.csv` 的每一列是一個共同交易日。第一欄 `date` 記錄日期，其餘 89 欄以股票代碼命名，儲存小數尺度的日簡單報酬。建檔時先在每一檔股票內依日期排序並計算報酬，再保留 89 檔都有觀察值的日期。這個順序很重要；若直接對尚未分組的長表取落後值，前一期價格可能來自另一檔股票。

## 日本月總體金融資料

`japan_monthly_2007_2018.csv`：

每一列代表一個月份，以下水準值與變動率各自使用檔案中的原尺度：

- `date`：月份。
- `spj`、`return_j`：日本股價指數與報酬。
- `rer`、`rer_change`：日本實質匯率與變動率。
- `ipi`、`ipi_change`：工業生產指數與變動率。
- `inr`、`inr_change`：利率與變動率。
- `spf`、`return_f`：美國股價指數與報酬。
- `unr`、`unr_change`：失業率與變動率。
- `cpi`：消費者物價指數。
- `bmj`、`bmj_growth`：日本 M3 與成長率。
- `opi`、`opi_change`：WTI 油價與變動率。
- `yng`、`mid`、`old`：年輕、中年與老年人口比率。
- `dro`、`dro_growth`、`dry`、`dep`、`dep_growth`：老年、幼年與總扶養比及其變動。
- `fpi`：外國證券投資。
- `trd`：貿易差額。
- `yield_10`、`yield_10_change`：日本十年期公債殖利率與百分比變動。

## 因子、總體預測變數與十產業

`ff_qf_macro_industries_1967_2021.csv` 每一列代表某個月份的一個產業，同一月份因此會有十列：

- `month`、`industry`、`ret`：月份、十產業代碼、產業超額報酬。
- `factor_ff_rf`、`factor_ff_mkt_excess`、`factor_ff_smb`、`factor_ff_hml`：Fama–French 無風險利率與三因子。
- `factor_q_me`、`factor_q_ia`、`factor_q_roe`、`factor_q_eg`：global-q 規模、投資、ROE 與預期成長因子。
- `macro_dp`、`macro_dy`、`macro_ep`、`macro_de`、`macro_svar`、`macro_bm`、`macro_ntis`、`macro_tbl`、`macro_lty`、`macro_ltr`、`macro_tms`、`macro_dfy`、`macro_infl`：Welch–Goyal 月總體預測變數。

## 加州學校資料

`california_schools.csv` 的每一列是一所學校，共有 110 欄。R15 使用的主要應變數與重點解釋變數如下：

- `testscore`：五年級數學與英語／語文測驗分數之和。
- `elarts_score`、`math_score`：英文與數學成績。
- `str_s`、`str_d`：學校與學區師生比。
- `ell_frac_s`、`frpm_frac_s`、`freem_frac_s`：英語學習者、免費或減價午餐等學生背景比例。
- `te_*`：教師年資、人數與薪資。
- `exp_*`、`rev_*`：學區支出與收入。
- `re_*`：族群組成比例。
- `age_*`、`sex_*`、`ms_*`、`ed_*`、`hs_*`、`med_income_z`：郵遞區號人口、教育、住宅與所得特徵。
- `countycode`、`districtcode`、`schoolcode`：行政識別碼；正文迴歸不把識別碼當作連續控制變數。

原課程附有 Stock and Watson, *Introduction to Econometrics*, 4th ed., Chapter 14 的完整變數表；R15 只從具經濟意義的連續背景欄位建立控制字典。

## 臺灣 ICAPM 第二階段面板

- `day`、`firm`：共同識別「公司—日期」觀察值；資料為 47 家公司、497 日的平衡面板。
- `return`：公司報酬。
- `cov_X`、`cov_vix`、`cov_liq`：原課程第一階段預先計算的市場、VIX 與流動性條件共變數。

現有檔案沒有保存第一階段原始序列與完整估計紀錄，因此分析時將這三個共變數視為已給定的輸入。這份資料足以重算第二階段面板分析；DCC 或 PQR 第一階段則需要原始序列與估計設定，無法只靠目前檔案重新取得。
