# Forecasting GLD Using an ARMA-GJR-GARCH Model

We model the daily price
behavior of GLD (SPDR Gold Shares) using an ARMA-GJR-GARCH time-series model,
capturing volatility clustering and asymmetric shocks in daily returns.

**Authors:** Shreyes Balaji, Marwan Hegab, Mazin Hussein, Mikael Rotberg,
Roberto Rubio, Jordan Woda

## Summary

- Data: GLD daily OHLCV history from Yahoo Finance (via `quantmod`), January
  2014 through the essay's cutoff of April 29, 2026 — 3,099 price rows, no
  missing values or duplicate dates.
- Modeled daily log returns (`100 * ln(P_t / P_{t-1})`), not raw prices,
  since prices are non-stationary and returns are better suited to
  volatility modeling — 3,098 return observations. Tukey's IQR rule flagged
  144 extreme return days; these were kept (not deleted) since extreme
  moves are often real market events, with a winsorized version saved
  separately for sensitivity checks.
- Chronological 70/30 train/test split (time series can't be split
  randomly): 2,168 training observations (2014-01-03 to 2022-08-12), 930
  test observations (2022-08-15 to 2026-04-29).
- Diagnostics: ADF test supports stationarity of returns (p < 0.01);
  Ljung-Box on squared returns (325.4, p < 0.001) and the ARCH LM test
  (127.7, p ≈ 1.8e-21) both confirm strong ARCH effects, motivating a
  GARCH-type model.
- Model comparison by AIC/BIC across ARMA orders, GARCH/GJR-GARCH variance
  specs, and normal/Student-t/GED error distributions selected
  **ARMA(0,0)-GJR-GARCH(1,1) with Student-t errors** (GJR-GARCH chosen per
  the assignment's requirement; comparable BIC to the top sGARCH candidate).
- Final model: omega = 0.0127, alpha1 = 0.0633, beta1 = 0.9429 (highly
  persistent volatility), gamma1 = -0.0389 (significant asymmetric
  shock effect), shape (Student-t) = 5.02. Residual diagnostics (Ljung-Box
  on standardized and squared standardized residuals, ARCH LM) all pass
  (p > 0.05), meaning the model removed the volatility clustering it was
  built to capture.
- Test-set forecasting: return RMSE = 1.211, MAE = 0.850, direction accuracy
  = 54.6%. Price-level forecasts undershoot because GLD saw a large,
  sustained rally during the test period that a return/volatility model
  isn't built to predict — GARCH models conditional variance, not long-run
  price trend.

This is best understood as a **volatility model**, not a price predictor: it
successfully captures and forecasts *how much* GLD's returns are likely to
swing day to day, but it can't anticipate a multi-year directional rally.
See the essay for the full diagnostic plots (Q-Q, ACF/PACF pre- and
post-fit), the AIC/BIC candidate-model table, and discussion of limitations
(no macroeconomic drivers, single-asset scope).

## Files

- [`project_6_essay.pdf`](project_6_essay.pdf) — full write-up (introduction, data
  description, analysis, model evaluation, conclusion, references).
- [`analysis.R`](analysis.R) — R code for the full analysis: data download and
  cleaning, stationarity/ARCH diagnostics, ARMA order selection, GARCH
  candidate comparison, residual diagnostics, rolling test-set forecasting,
  and a 10-day future forecast.

## Reproducing the analysis

```r
install.packages(c("quantmod", "tidyverse", "lubridate", "tseries", "FinTS",
                    "forecast", "rugarch", "PerformanceAnalytics", "knitr", "zoo", "xts"))
```

No local data files are needed — `analysis.R` downloads GLD's full daily
price history directly from Yahoo Finance via `quantmod::getSymbols()`, so
just run it top to bottom. Because it pulls live market data, running it
today will include more history than the essay's April 2026 cutoff, so exact
figures (row counts, coefficients, forecasts) will differ slightly from the
essay above — the selected model and overall conclusions should still match.
Outputs (cleaned CSVs, diagnostic tables, and plots) are written to
`gld_project6_outputs/` in the working directory. The GARCH candidate grid
search and rolling forecast take a few minutes to run.

## References

1. Ghani, I. M. Md., & Rahim, H. A. (2019). Modeling and Forecasting of Volatility using ARMA-GARCH: Case Study on Malaysia Natural Rubber Prices. *IOP Conference Series: Materials Science and Engineering*, 548, 012023.
2. Hasan, M. (2021). Microsoft Stock Price Analysis using GARCH model in R. RPubs. https://rpubs.com/Mahmud_Hasan/778532
3. Yahoo Finance — GLD Historical Data: https://finance.yahoo.com/quote/GLD/history/
4. Ghalanos, A. — `rugarch`: Univariate GARCH Models (R package documentation).
