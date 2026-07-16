---
title: "R05：ARMA 估計、診斷與多步預測"
output:
  github_document:
    toc: true
    toc_depth: 3
---

本附錄對應第 6--7 章。以固定種子 ARMA(1,1) 資料示範：先保留時間排序的測試期，只在訓練樣本比較低階候選模型，再檢查根、殘差與平方殘差，最後形成一次多步預測與區間。因資料生成參數已知，可分辨估計誤差與程式錯誤。


``` r
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.width = 7, fig.height = 4.5
)
set.seed(20260716)
```

## 固定資料生成過程與時間切分


``` r
n <- 700L
true_phi <- 0.55
true_theta <- -0.35
innovation_sd <- 1.2

y <- as.numeric(arima.sim(
  model = list(ar = true_phi, ma = true_theta),
  n = n, sd = innovation_sd
))

train_end <- 550L
y_train <- y[1:train_end]
y_test <- y[(train_end + 1L):n]

c(
  training_observations = length(y_train),
  test_observations = length(y_test),
  first_test_index = train_end + 1L
)
```

```
## training_observations     test_observations      first_test_index 
##                   550                   150                   551
```

測試期在模型選擇時完全不使用。這個附錄比較的是「訓練期末一次形成的多步預測」，不是每期重新估計的一步策略；後者見 R06。

## 候選模型與共同樣本

候選集合事前限制為低階模型。stats::arima 對定態 ARMA 的 intercept 係數代表序列平均數，不是
\(Y_t=c+\phi Y_{t-1}+a_t\) 中的 \(c\)；判讀時必須注意參數化。


``` r
candidate_orders <- list(
  ARMA00 = c(0, 0, 0),
  AR10 = c(1, 0, 0),
  MA01 = c(0, 0, 1),
  ARMA11 = c(1, 0, 1),
  AR20 = c(2, 0, 0),
  MA02 = c(0, 0, 2)
)

fits <- lapply(candidate_orders, function(ord) {
  arima(
    y_train,
    order = ord,
    include.mean = TRUE,
    method = "ML"
  )
})

model_table <- data.frame(
  model = names(fits),
  p = vapply(candidate_orders, function(z) z[1], numeric(1)),
  q = vapply(candidate_orders, function(z) z[3], numeric(1)),
  log_likelihood = vapply(fits, function(z) as.numeric(logLik(z)), numeric(1)),
  AIC = vapply(fits, AIC, numeric(1)),
  BIC = vapply(fits, BIC, numeric(1)),
  row.names = NULL
)
model_table[order(model_table$AIC), ]
```

```
##    model p q log_likelihood      AIC      BIC
## 4 ARMA11 1 1      -910.9506 1829.901 1847.141
## 5   AR20 2 0      -911.6655 1831.331 1848.571
## 2   AR10 1 0      -913.3429 1832.686 1845.615
## 6   MA02 0 2      -912.8188 1833.638 1850.877
## 3   MA01 0 1      -915.3380 1836.676 1849.606
## 1 ARMA00 0 0      -926.4747 1856.949 1865.569
```

所有模型使用同一訓練序列；這避免 AIC 差異混入樣本不同。AIC/BIC 只負責縮小候選集合，不是測試期成績。

## 檢查 AR 與 MA 根


``` r
extract_roots <- function(fit) {
  b <- coef(fit)
  ar_coef <- b[grep("^ar[0-9]+$", names(b))]
  ma_coef <- b[grep("^ma[0-9]+$", names(b))]

  ar_roots <- if (length(ar_coef)) {
    polyroot(c(1, -ar_coef))
  } else {
    complex()
  }
  ma_roots <- if (length(ma_coef)) {
    polyroot(c(1, ma_coef))
  } else {
    complex()
  }

  rbind(
    if (length(ar_roots)) data.frame(
      part = "AR", root = ar_roots, modulus = Mod(ar_roots)
    ),
    if (length(ma_roots)) data.frame(
      part = "MA", root = ma_roots, modulus = Mod(ma_roots)
    )
  )
}

root_results <- lapply(names(fits), function(nm) {
  out <- extract_roots(fits[[nm]])
  if (!is.null(out) && nrow(out)) out$model <- nm
  out
})
root_results <- do.call(rbind, root_results)
root_results[, c("model", "part", "root", "modulus")]
```

```
##    model part                    root  modulus
## 1   AR10   AR  4.633781+0.000000e+00i 4.633781
## 2   MA01   MA -5.385256+0.000000e+00i 5.385256
## 3 ARMA11   AR  1.545550+0.000000e+00i 1.545550
## 4 ARMA11   MA  2.156871+0.000000e+00i 2.156871
## 5   AR20   AR  2.525975-6.286573e-23i 2.525975
## 6   AR20   AR -5.073672+6.286573e-23i 5.073672
## 7   MA02   MA -1.035531+3.120070e+00i 3.287425
## 8   MA02   MA -1.035531-3.120070e+00i 3.287425
```

``` r
stopifnot(all(root_results$modulus > 1))
```

AR 根大於 1 表示估計的定態表示；MA 根大於 1 表示可逆表示。根非常接近 1 時，即使形式上符合限制，有限樣本推論與遠期預測仍會不穩定。

## 殘差與平方殘差診斷


``` r
diagnose_fit <- function(fit, lag = 20L) {
  e <- residuals(fit)
  e <- e[is.finite(e)]
  fitted_parameters <- length(grep(
    "^(ar|ma)[0-9]+$", names(coef(fit))
  ))
  q_mean <- Box.test(
    e, lag = lag, type = "Ljung-Box",
    fitdf = fitted_parameters
  )
  q_square <- Box.test(
    e^2, lag = lag, type = "Ljung-Box"
  )
  c(
    residual_sd = sd(e),
    Q_mean = unname(q_mean$statistic),
    p_mean = q_mean$p.value,
    Q_square = unname(q_square$statistic),
    p_square = q_square$p.value
  )
}

diagnostic_table <- do.call(rbind, lapply(fits, diagnose_fit))
round(diagnostic_table, 4)
```

```
##        residual_sd  Q_mean p_mean Q_square p_square
## ARMA00      1.3053 67.6535 0.0000  14.5255   0.8029
## AR10        1.2745 23.1674 0.2301  13.5767   0.8513
## MA01        1.2791 30.1222 0.0503  13.9994   0.8305
## ARMA11      1.2689 15.7596 0.6093  13.8350   0.8388
## AR20        1.2706 18.0092 0.4550  13.1231   0.8720
## MA02        1.2733 21.9471 0.2343  12.6408   0.8923
```

殘差檢定應與圖形一起看。通過報酬殘差檢查只表示未留下明顯線性平均相依；平方殘差若相關，可能需要條件變異模型。


``` r
selected_name <- model_table$model[which.min(model_table$AIC)]
selected_fit <- fits[[selected_name]]
selected_residual <- residuals(selected_fit)

par(mfrow = c(1, 3))
plot(
  selected_residual, type = "l", col = "#173B57",
  xlab = "訓練期", ylab = "殘差", main = selected_name
)
acf(selected_residual, lag.max = 30, main = "殘差 ACF")
acf(selected_residual^2, lag.max = 30, main = "平方殘差 ACF")
```

![AIC 最小模型的殘差、殘差 ACF 與平方殘差 ACF。](./R05_arma_estimation_diagnostics_forecasting_files/figure-gfm/diagnostics-plots-1.png)

``` r
par(mfrow = c(1, 1))
```

## 多步預測與區間


``` r
h <- length(y_test)
forecast_object <- predict(selected_fit, n.ahead = h)
forecast_table <- data.frame(
  index = (train_end + 1L):n,
  actual = y_test,
  forecast = as.numeric(forecast_object$pred),
  standard_error = as.numeric(forecast_object$se)
)
forecast_table$lower95 <- forecast_table$forecast -
  1.96 * forecast_table$standard_error
forecast_table$upper95 <- forecast_table$forecast +
  1.96 * forecast_table$standard_error
forecast_table$error <- forecast_table$actual - forecast_table$forecast

head(forecast_table, 8)
```

```
##   index     actual   forecast standard_error   lower95  upper95      error
## 1   551 -1.2182263 0.04651117       1.267784 -2.438346 2.531368 -1.2647375
## 2   552  0.2705367 0.05047878       1.288925 -2.475815 2.576773  0.2200580
## 3   553 -0.8689805 0.05304589       1.297674 -2.490394 2.596486 -0.9220264
## 4   554 -1.1451046 0.05470686       1.301318 -2.495877 2.605291 -1.1998115
## 5   555  1.0903613 0.05578154       1.302841 -2.497787 2.609350  1.0345797
## 6   556 -0.4237374 0.05647688       1.303478 -2.498341 2.611294 -0.4802143
## 7   557 -0.3999396 0.05692678       1.303745 -2.498413 2.612267 -0.4568664
## 8   558  0.6426250 0.05721787       1.303856 -2.498341 2.612776  0.5854072
```

``` r
tail(forecast_table, 3)
```

```
##     index    actual   forecast standard_error   lower95  upper95     error
## 148   698 0.9368245 0.05775144       1.303937 -2.497965 2.613467 0.8790731
## 149   699 1.1217299 0.05775144       1.303937 -2.497965 2.613467 1.0639785
## 150   700 1.6704255 0.05775144       1.303937 -2.497965 2.613467 1.6126741
```


``` r
data.frame(
  selected_model = selected_name,
  test_RMSE = sqrt(mean(forecast_table$error^2)),
  test_MAE = mean(abs(forecast_table$error)),
  interval_coverage = mean(
    forecast_table$actual >= forecast_table$lower95 &
      forecast_table$actual <= forecast_table$upper95
  ),
  average_width = mean(
    forecast_table$upper95 - forecast_table$lower95
  )
)
```

```
##   selected_model test_RMSE test_MAE interval_coverage average_width
## 1         ARMA11  1.234145 1.012712              0.96      5.109814
```

這裡的常態區間主要反映未來創新，對參數與模型選擇不確定性的處理有限。涵蓋率只是這一個固定測試路徑的描述，不是母體涵蓋率證明。


``` r
plot(
  forecast_table$index, forecast_table$actual,
  type = "l", col = "gray35",
  xlab = "時間索引", ylab = "數值"
)
polygon(
  c(forecast_table$index, rev(forecast_table$index)),
  c(forecast_table$lower95, rev(forecast_table$upper95)),
  border = NA, col = adjustcolor("#9FC2D4", alpha.f = 0.45)
)
lines(
  forecast_table$index, forecast_table$forecast,
  col = "#A34045", lwd = 2
)
lines(forecast_table$index, forecast_table$actual, col = "gray35")
legend(
  "topright", c("實際值", "點預測", "95% 區間"),
  col = c("gray35", "#A34045", "#9FC2D4"),
  lty = c(1, 1, NA), pch = c(NA, NA, 15), bty = "n"
)
```

![訓練期末一次形成的多步預測與 95% 常態近似區間。](./R05_arma_estimation_diagnostics_forecasting_files/figure-gfm/forecast-plot-1.png)

## 手算核對 AR(1) 預測

另建一個已知平均數的 AR(1) 教學例，核對
\(\widehat Y_{T+h\mid T}=\mu+\phi^h(Y_T-\mu)\) 與
\(\sigma_h^2=\sigma_a^2(1-\phi^{2h})/(1-\phi^2)\)。


``` r
mu <- 2
phi_check <- 0.6
innovation_variance <- 4
y_T <- 5
h_check <- 1:5

manual_check <- data.frame(
  horizon = h_check,
  point_forecast = mu + phi_check^h_check * (y_T - mu),
  error_variance = innovation_variance *
    (1 - phi_check^(2 * h_check)) / (1 - phi_check^2)
)
manual_check$standard_error <- sqrt(manual_check$error_variance)
manual_check
```

```
##   horizon point_forecast error_variance standard_error
## 1       1        3.80000       4.000000       2.000000
## 2       2        3.08000       5.440000       2.332381
## 3       3        2.64800       5.958400       2.440983
## 4       4        2.38880       6.145024       2.478916
## 5       5        2.23328       6.212209       2.492430
```

## 建模紀錄

正式分析至少保存資料版本、訓練截止日、候選集合、估計方法、截距參數化、AIC/BIC、根、殘差診斷、預測期距與區間假設。若後來改看測試結果再換模型，應把該期間重新歸類為驗證期。
