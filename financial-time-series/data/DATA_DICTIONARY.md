# 實證資料字典

本頁列出線上附錄實際使用的欄位、單位與時間鍵。完整列數、欄數與 MD5 以 `processed/manifest.csv` 為準。

## AAPL 與 MSFT

- `aapl_adjusted_daily_2019_2022.csv`
  - `date`：交易日。
  - `adjusted`：原價格檔的調整後收盤價，美元。
  - `simple_return`：$P_t/P_{t-1}-1$。
  - `log_return`：$\log P_t-\log P_{t-1}$。
  - `symbol`、`company`、`sector`：證券與產業識別欄。
- `msft_daily_returns_1986_2008.csv`
  - `date`：交易日。
  - `simple_return`：Microsoft 日簡單報酬，小數表示。

## FRED 日匯率

`fred_jpy_twd_daily_2020_2022.csv`：

- `date`：日曆日；非交易日或原序列缺值保留為空值。
- `jpy_per_usd`：一美元可兌日圓數，FRED `DEXJPUS`。
- `twd_per_usd`：一美元可兌新臺幣數，FRED `DEXTAUS`。
- `twd_per_jpy`：`twd_per_usd / jpy_per_usd`。
- `log_return_twd_per_jpy`：交叉匯率的一期對數差分；遇缺值時不跨缺值補算。

## Tsay 股票面板

- `tsay_five_stock_monthly_returns_1990_2008.csv`：`month` 加 IBM、HPQ、INTC、JPM、BAC 五家公司月對數報酬。
- `tsay_barra_monthly_returns_1990_2003.csv`：`month` 加 AGE、C、MWD、MER、DELL、HPQ、IBM、AA、CAT、PG 十家公司月報酬。原檔數值以百分點表示；程式在需要時除以 100。

## S&P 500 平衡報酬面板

`sp500_returns_balanced_2013_2022.csv` 的第一欄 `date` 是共同交易日，其餘 89 欄是股票代碼，儲存日簡單報酬。建檔時先在 `symbol` 內依日期排序與計算報酬，再保留全部 89 檔都有觀察值的日期。

## 日本月總體金融資料

`japan_monthly_2007_2018.csv`：

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

`ff_qf_macro_industries_1967_2021.csv` 每月每產業一列：

- `month`、`industry`、`ret`：月份、十產業代碼、產業超額報酬。
- `factor_ff_rf`、`factor_ff_mkt_excess`、`factor_ff_smb`、`factor_ff_hml`：Fama--French 無風險利率與三因子。
- `factor_q_me`、`factor_q_ia`、`factor_q_roe`、`factor_q_eg`：global-q 規模、投資、ROE 與預期成長因子。
- `macro_dp`、`macro_dy`、`macro_ep`、`macro_de`、`macro_svar`、`macro_bm`、`macro_ntis`、`macro_tbl`、`macro_lty`、`macro_ltr`、`macro_tms`、`macro_dfy`、`macro_infl`：Welch--Goyal 月總體預測變數。

## California schools

`california_schools.csv` 有 110 欄。主要應變數與處置／解釋變數為：

- `testscore`：英文與數學平均測驗成績。
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

- `day`、`firm`：497 日與 47 家公司的平衡面板鍵。
- `return`：公司報酬。
- `cov_X`、`cov_vix`、`cov_liq`：原課程第一階段預先計算的市場、VIX 與流動性條件共變數。

因缺少第一階段原始序列與完整估計紀錄，這三個共變數只能視為給定輸入；不能由此檔獨立驗證 DCC 或 PQR 的第一階段。
