# Forecasting GLD Using an ARMA-GJR-GARCH Model
#
# How to run:
#   1. Install packages once:
#      install.packages(c("quantmod","tidyverse","lubridate","tseries","FinTS",
#                          "forecast","rugarch","PerformanceAnalytics","knitr","zoo","xts"))
#   2. Just run the script top to bottom (source it, or run section by section).
#      No local data files are needed -- it downloads GLD's full daily price
#      history from Yahoo Finance via quantmod::getSymbols().
#   3. Outputs (cleaned data, tables, and plots) are written to
#      gld_project6_outputs/ in the working directory.
#
# Note: because this pulls live market data, running it today will use more
# history than the essay's April 2026 cutoff, so exact figures (row counts,
# coefficients, forecasts) will differ slightly from the essay -- the model
# selected and the overall conclusions should still match.

# SECTION 0: setup ------------------------------------------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))
needed_packages <- c(
  "quantmod", "tidyverse", "lubridate", "tseries", "FinTS",
  "forecast", "rugarch", "PerformanceAnalytics", "knitr", "zoo", "xts"
)
installed <- rownames(installed.packages())
for (p in needed_packages) {
  if (!(p %in% installed)) install.packages(p)
}

library(quantmod)
library(tidyverse)
library(lubridate)
library(tseries)
library(FinTS)
library(forecast)
library(rugarch)
library(PerformanceAnalytics)
library(knitr)
library(zoo)
library(xts)

set.seed(456)

# SECTION 1: download and clean GLD daily data from Yahoo Finance -----------
symbol <- "GLD"
start_date <- as.Date("2014-01-01")
end_date <- Sys.Date()

out_dir <- "gld_project6_outputs"
dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)

gld_xts <- getSymbols(symbol, src = "yahoo", from = start_date, to = end_date,
                      auto.assign = FALSE)

gld_prices <- data.frame(date = index(gld_xts), coredata(gld_xts))
colnames(gld_prices) <- c("date", "open", "high", "low", "close", "volume", "adjusted")

# Basic raw data checks
raw_missing_counts <- sapply(gld_prices, function(x) sum(is.na(x)))
raw_duplicate_dates <- sum(duplicated(gld_prices$date))

# Clean prices: remove missing values, duplicate dates, and zero/negative volume days
gld_prices_clean <- gld_prices %>%
  arrange(date) %>%
  distinct(date, .keep_all = TRUE) %>%
  filter(if_all(c(open, high, low, close, volume, adjusted), ~ !is.na(.x))) %>%
  filter(volume > 0)

# Create daily log returns using adjusted close.
# Multiplying by 100 puts returns in percentage points.
gld_returns <- gld_prices_clean %>%
  mutate(
    log_return = 100 * (log(adjusted) - log(lag(adjusted))),
    simple_return = 100 * ((adjusted / lag(adjusted)) - 1)
  ) %>%
  drop_na(log_return, simple_return)

# Flag possible return outliers using Tukey's IQR rule.
Q1 <- quantile(gld_returns$log_return, 0.25, na.rm = TRUE)
Q3 <- quantile(gld_returns$log_return, 0.75, na.rm = TRUE)
IQR_val <- Q3 - Q1
lower_fence <- Q1 - 1.5 * IQR_val
upper_fence <- Q3 + 1.5 * IQR_val

gld_returns <- gld_returns %>%
  mutate(return_outlier_iqr = log_return < lower_fence | log_return > upper_fence)

# Optional sensitivity dataset: winsorized returns. Extreme return days in
# financial markets are often real market events rather than data errors, so
# they are flagged and kept in the main model rather than deleted -- a
# winsorized version is saved separately for anyone who wants to check
# sensitivity to that choice.
winsorize <- function(x, lower_prob = 0.01, upper_prob = 0.99) {
  lo <- quantile(x, lower_prob, na.rm = TRUE)
  hi <- quantile(x, upper_prob, na.rm = TRUE)
  pmin(pmax(x, lo), hi)
}
gld_returns_winsorized <- gld_returns %>%
  mutate(log_return_winsorized = winsorize(log_return))

# Save starting datasets
write.csv(gld_prices_clean, file.path(out_dir, "gld_daily_prices_clean.csv"), row.names = FALSE)
write.csv(gld_returns, file.path(out_dir, "gld_daily_returns.csv"), row.names = FALSE)
write.csv(gld_returns_winsorized, file.path(out_dir, "gld_daily_returns_winsorized_sensitivity.csv"), row.names = FALSE)

# Data summary table
summary_table <- tibble(
  ticker = symbol,
  start_date = min(gld_prices_clean$date),
  end_date = max(gld_prices_clean$date),
  price_rows = nrow(gld_prices_clean),
  return_rows = nrow(gld_returns),
  raw_duplicate_dates = raw_duplicate_dates,
  total_raw_missing_values = sum(raw_missing_counts),
  iqr_return_outliers_flagged = sum(gld_returns$return_outlier_iqr),
  mean_adjusted_close = mean(gld_prices_clean$adjusted),
  min_adjusted_close = min(gld_prices_clean$adjusted),
  max_adjusted_close = max(gld_prices_clean$adjusted),
  mean_daily_log_return = mean(gld_returns$log_return),
  sd_daily_log_return = sd(gld_returns$log_return)
)
write.csv(summary_table, file.path(out_dir, "data_summary_table.csv"), row.names = FALSE)
print(summary_table)
cat("\nMissing values by raw column:\n")
print(raw_missing_counts)

# SECTION 2: data visualizations ---------------------------------------------
p_price <- ggplot(gld_prices_clean, aes(x = date, y = adjusted)) +
  geom_line() +
  labs(title = "GLD Adjusted Closing Price Over Time",
       x = "Date", y = "Adjusted Close Price (USD)") +
  theme_minimal()
ggsave(file.path(out_dir, "plots", "01_gld_adjusted_price.png"), p_price, width = 9, height = 5)

p_returns <- ggplot(gld_returns, aes(x = date, y = log_return)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "GLD Daily Log Returns Over Time",
       x = "Date", y = "Daily Log Return (%)") +
  theme_minimal()
ggsave(file.path(out_dir, "plots", "02_gld_log_returns_time_series.png"), p_returns, width = 9, height = 5)

p_hist <- ggplot(gld_returns, aes(x = log_return)) +
  geom_histogram(bins = 60) +
  labs(title = "Histogram of GLD Daily Log Returns",
       x = "Daily Log Return (%)", y = "Count") +
  theme_minimal()
ggsave(file.path(out_dir, "plots", "03_gld_log_returns_histogram.png"), p_hist, width = 8, height = 5)

p_box <- ggplot(gld_returns, aes(y = log_return)) +
  geom_boxplot() +
  labs(title = "Boxplot of GLD Daily Log Returns", y = "Daily Log Return (%)") +
  theme_minimal()
ggsave(file.path(out_dir, "plots", "04_gld_log_returns_boxplot.png"), p_box, width = 6, height = 5)

# Base R diagnostic plots saved as PNG
png(file.path(out_dir, "plots", "05_acf_pacf_returns.png"), width = 1000, height = 700)
par(mfrow = c(1, 2))
acf(gld_returns$log_return, main = "ACF: GLD Log Returns")
pacf(gld_returns$log_return, main = "PACF: GLD Log Returns")
par(mfrow = c(1, 1))
dev.off()

png(file.path(out_dir, "plots", "06_acf_pacf_squared_returns.png"), width = 1000, height = 700)
par(mfrow = c(1, 2))
acf(gld_returns$log_return^2, main = "ACF: Squared GLD Returns")
pacf(gld_returns$log_return^2, main = "PACF: Squared GLD Returns")
par(mfrow = c(1, 1))
dev.off()

png(file.path(out_dir, "plots", "07_qq_plot_returns_diagonal_plot.png"), width = 800, height = 700)
qqnorm(gld_returns$log_return, main = "Normal Q-Q Plot of GLD Daily Log Returns")
qqline(gld_returns$log_return, lwd = 2)
dev.off()

# SECTION 3: time-series train/test split ------------------------------------
# Chronological split, not random, because time-series forecasting must train
# on the past and test on the future.
n <- nrow(gld_returns)
train_n <- floor(0.70 * n)
train_data <- gld_returns[1:train_n, ]
test_data <- gld_returns[(train_n + 1):n, ]

train_ret <- train_data$log_return
test_ret <- test_data$log_return

split_table <- tibble(
  dataset = c("Training set", "Test set"),
  rows = c(nrow(train_data), nrow(test_data)),
  start_date = c(min(train_data$date), min(test_data$date)),
  end_date = c(max(train_data$date), max(test_data$date))
)
write.csv(split_table, file.path(out_dir, "train_test_split_table.csv"), row.names = FALSE)
print(split_table)

# SECTION 4: assumption / preliminary tests ----------------------------------
# ADF test: checks stationarity of returns. For GARCH, returns should usually be stationary.
adf_result <- adf.test(train_ret)

# Ljung-Box test on returns: checks autocorrelation in the mean.
lb_returns <- Box.test(train_ret, lag = 20, type = "Ljung-Box")

# Ljung-Box test on squared returns and ARCH test: checks volatility clustering / ARCH effects.
lb_squared_returns <- Box.test(train_ret^2, lag = 20, type = "Ljung-Box")
arch_result <- ArchTest(train_ret, lags = 12)

prelim_tests <- tibble(
  test = c("ADF stationarity test", "Ljung-Box returns", "Ljung-Box squared returns", "ARCH LM test"),
  statistic = c(
    unname(adf_result$statistic),
    unname(lb_returns$statistic),
    unname(lb_squared_returns$statistic),
    unname(arch_result$statistic)
  ),
  p_value = c(
    adf_result$p.value,
    lb_returns$p.value,
    lb_squared_returns$p.value,
    arch_result$p.value
  )
)
write.csv(prelim_tests, file.path(out_dir, "preliminary_tests.csv"), row.names = FALSE)
print(prelim_tests)

# SECTION 5: ARMA mean model selection ---------------------------------------
# BIC is used because it favors a simpler model.
arma_fit_auto <- auto.arima(
  train_ret,
  max.p = 5, max.q = 5, d = 0,
  seasonal = FALSE, stationary = TRUE,
  ic = "bic", stepwise = FALSE, approximation = FALSE
)
print(arma_fit_auto)
arma_order <- arma_fit_auto$arma
arma_p <- arma_order[1]
arma_q <- arma_order[2]
cat("\nSelected ARMA order by BIC: ARMA(", arma_p, ",", arma_q, ")\n", sep = "")

# SECTION 6: candidate GARCH model comparison --------------------------------
# This is the model/subset selection step for this project: compare simple
# mean equations, variance models, and innovation distributions. Because the
# project specifically requires ARMA-GJR-GARCH, the final model is chosen
# from the converged GJR-GARCH candidates by BIC. Symmetric sGARCH is kept
# as a benchmark for comparison.
train_ret <- as.numeric(train_ret)
train_ret <- train_ret[is.finite(train_ret)]

aic_from_fit <- function(fit) {
  ic <- infocriteria(fit)
  as.numeric(ic[grep("Akaike", rownames(ic)), 1])
}
bic_from_fit <- function(fit) {
  ic <- infocriteria(fit)
  as.numeric(ic[grep("Bayes", rownames(ic)), 1])
}

fit_garch_candidate <- function(variance_model, ar_order, ma_order, dist_model) {
  spec <- ugarchspec(
    variance.model = list(model = variance_model, garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(ar_order, ma_order), include.mean = TRUE),
    distribution.model = dist_model
  )

  solvers_to_try <- c("hybrid", "nlminb", "solnp")
  last_error <- NA_character_

  for (solver_name in solvers_to_try) {
    fit <- tryCatch(
      suppressWarnings(
        ugarchfit(
          spec = spec,
          data = train_ret,
          solver = solver_name,
          solver.control = list(trace = 0),
          fit.control = list(scale = 1)
        )
      ),
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )

    if (!is.null(fit)) {
      aic_val <- tryCatch(aic_from_fit(fit), error = function(e) NA_real_)
      bic_val <- tryCatch(bic_from_fit(fit), error = function(e) NA_real_)
      return(tibble(
        variance_model = variance_model,
        ar = ar_order,
        ma = ma_order,
        distribution = dist_model,
        solver = solver_name,
        convergence = fit@fit$convergence,
        AIC = aic_val,
        BIC = bic_val,
        error_message = NA_character_,
        fit_object = list(fit)
      ))
    }
  }

  tibble(
    variance_model = variance_model,
    ar = ar_order,
    ma = ma_order,
    distribution = dist_model,
    solver = NA_character_,
    convergence = NA_integer_,
    AIC = NA_real_,
    BIC = NA_real_,
    error_message = last_error,
    fit_object = list(NULL)
  )
}

candidate_grid <- expand.grid(
  variance_model = c("sGARCH", "gjrGARCH"),
  ar_order = unique(c(0, arma_p, 1)),
  ma_order = unique(c(0, arma_q, 1)),
  dist_model = c("norm", "std", "ged"),
  stringsAsFactors = FALSE
)

candidate_results <- pmap_dfr(
  candidate_grid,
  function(variance_model, ar_order, ma_order, dist_model) {
    fit_garch_candidate(variance_model, ar_order, ma_order, dist_model)
  }
)

candidate_results_ranked <- candidate_results %>%
  filter(!is.na(BIC), !map_lgl(fit_object, is.null)) %>%
  arrange(convergence, BIC)

write.csv(
  candidate_results %>% select(-fit_object),
  file.path(out_dir, "garch_candidate_model_comparison_all_attempts.csv"),
  row.names = FALSE
)
write.csv(
  candidate_results_ranked %>% select(-fit_object),
  file.path(out_dir, "garch_candidate_model_comparison.csv"),
  row.names = FALSE
)

cat("\nTop GARCH candidate models ranked by convergence and BIC:\n")
print(candidate_results_ranked %>% select(-fit_object) %>% head(12))

if (nrow(candidate_results_ranked) == 0) {
  stop("No GARCH models could be fit. Check garch_candidate_model_comparison_all_attempts.csv for errors.")
}

# Choose the final model. We prefer the best converged GJR-GARCH candidate
# because the project requires an ARMA-GJR-GARCH model.
gjr_converged <- candidate_results_ranked %>%
  filter(variance_model == "gjrGARCH", convergence == 0)

if (nrow(gjr_converged) > 0) {
  best_row <- gjr_converged[1, ]
  final_selection_note <- "Selected best converged GJR-GARCH candidate by BIC because the project requires ARMA-GJR-GARCH."
} else {
  best_row <- candidate_results_ranked[1, ]
  final_selection_note <- "No GJR-GARCH candidate converged; selected best available GARCH candidate by BIC."
}

best_fit <- best_row$fit_object[[1]]
best_variance_model <- best_row$variance_model
best_ar <- as.integer(best_row$ar)
best_ma <- as.integer(best_row$ma)
best_distribution <- best_row$distribution

final_model_choice <- tibble(
  final_variance_model = best_variance_model,
  final_arma_order = paste0("ARMA(", best_ar, ",", best_ma, ")"),
  final_distribution = best_distribution,
  solver = best_row$solver,
  convergence = best_row$convergence,
  AIC = best_row$AIC,
  BIC = best_row$BIC,
  selection_note = final_selection_note
)
write.csv(final_model_choice, file.path(out_dir, "final_model_choice.csv"), row.names = FALSE)
print(final_model_choice)

final_spec <- ugarchspec(
  variance.model = list(model = best_variance_model, garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(best_ar, best_ma), include.mean = TRUE),
  distribution.model = best_distribution
)

# Reuse the chosen fitted object as the final fit.
final_fit <- best_fit
capture.output(show(final_fit), file = file.path(out_dir, "final_garch_model_summary.txt"))

# Standardized residual diagnostics
std_resid <- residuals(final_fit, standardize = TRUE)
std_resid <- as.numeric(std_resid)
std_resid <- std_resid[!is.na(std_resid)]

resid_diag <- tibble(
  test = c("Ljung-Box standardized residuals", "Ljung-Box squared standardized residuals", "ARCH LM standardized residuals"),
  statistic = c(
    unname(Box.test(std_resid, lag = 20, type = "Ljung-Box")$statistic),
    unname(Box.test(std_resid^2, lag = 20, type = "Ljung-Box")$statistic),
    unname(ArchTest(std_resid, lags = 12)$statistic)
  ),
  p_value = c(
    Box.test(std_resid, lag = 20, type = "Ljung-Box")$p.value,
    Box.test(std_resid^2, lag = 20, type = "Ljung-Box")$p.value,
    ArchTest(std_resid, lags = 12)$p.value
  )
)
write.csv(resid_diag, file.path(out_dir, "final_model_residual_diagnostics.csv"), row.names = FALSE)
print(resid_diag)

png(file.path(out_dir, "plots", "08_final_model_standardized_residuals_qq_plot.png"), width = 800, height = 700)
qqnorm(std_resid, main = "Q-Q Plot of Standardized GARCH Residuals")
qqline(std_resid, lwd = 2)
dev.off()

png(file.path(out_dir, "plots", "09_final_model_residual_acf.png"), width = 1000, height = 700)
par(mfrow = c(1, 2))
acf(std_resid, main = "ACF: Standardized Residuals")
acf(std_resid^2, main = "ACF: Squared Standardized Residuals")
par(mfrow = c(1, 1))
dev.off()

# SECTION 7: prediction / forecasting on the test set ------------------------
# First try rolling one-step-ahead forecasts. This is best for evaluation.
roll_success <- TRUE
roll <- tryCatch(
  suppressWarnings(
    ugarchroll(
      spec = final_spec,
      data = gld_returns$log_return,
      n.start = train_n,
      refit.every = 50,
      refit.window = "recursive",
      solver = "hybrid",
      solver.control = list(trace = 0),
      calculate.VaR = FALSE,
      keep.coef = TRUE
    )
  ),
  error = function(e) {
    roll_success <<- FALSE
    return(NULL)
  }
)

if (roll_success) {
  roll_df <- as.data.frame(roll)
  forecast_df <- tibble(
    date = test_data$date,
    actual_return = as.numeric(roll_df$Realized),
    predicted_return = as.numeric(roll_df$Mu),
    predicted_sigma = as.numeric(roll_df$Sigma),
    forecast_method = "rolling_one_step"
  )
} else {
  # Fallback: multi-step forecast from the end of the training set.
  fc_test <- ugarchforecast(final_fit, n.ahead = length(test_ret))
  forecast_df <- tibble(
    date = test_data$date,
    actual_return = test_ret,
    predicted_return = as.numeric(fitted(fc_test)),
    predicted_sigma = as.numeric(sigma(fc_test)),
    forecast_method = "multi_step_fallback"
  )
}

# Convert return forecasts into adjusted price forecasts.
last_train_price <- tail(train_data$adjusted, 1)
forecast_df <- forecast_df %>%
  mutate(
    predicted_adjusted_price = last_train_price * exp(cumsum(predicted_return / 100)),
    actual_adjusted_price = test_data$adjusted,
    return_error = actual_return - predicted_return,
    price_error = actual_adjusted_price - predicted_adjusted_price
  )

model_accuracy <- tibble(
  final_variance_model = best_variance_model,
  final_arma_order = paste0("ARMA(", best_ar, ",", best_ma, ")"),
  final_distribution = best_distribution,
  return_RMSE = sqrt(mean(forecast_df$return_error^2, na.rm = TRUE)),
  return_MAE = mean(abs(forecast_df$return_error), na.rm = TRUE),
  price_RMSE = sqrt(mean(forecast_df$price_error^2, na.rm = TRUE)),
  price_MAE = mean(abs(forecast_df$price_error), na.rm = TRUE),
  direction_accuracy = mean(sign(forecast_df$actual_return) == sign(forecast_df$predicted_return), na.rm = TRUE),
  volatility_proxy_RMSE_abs_return_vs_sigma = sqrt(mean((abs(forecast_df$actual_return) - forecast_df$predicted_sigma)^2, na.rm = TRUE))
)

write.csv(forecast_df, file.path(out_dir, "test_set_forecasts.csv"), row.names = FALSE)
write.csv(model_accuracy, file.path(out_dir, "model_accuracy_table.csv"), row.names = FALSE)
print(model_accuracy)

p_pred_returns <- ggplot(forecast_df, aes(x = date)) +
  geom_line(aes(y = actual_return), linewidth = 0.4) +
  geom_line(aes(y = predicted_return), linewidth = 0.4, linetype = "dashed") +
  labs(title = "Actual vs Predicted GLD Daily Log Returns on Test Set",
       x = "Date", y = "Daily Log Return (%)") +
  theme_minimal()
ggsave(file.path(out_dir, "plots", "10_actual_vs_predicted_returns.png"), p_pred_returns, width = 9, height = 5)

p_pred_price <- ggplot(forecast_df, aes(x = date)) +
  geom_line(aes(y = actual_adjusted_price), linewidth = 0.4) +
  geom_line(aes(y = predicted_adjusted_price), linewidth = 0.4, linetype = "dashed") +
  labs(title = "Actual vs Predicted GLD Adjusted Price on Test Set",
       x = "Date", y = "Adjusted Close Price (USD)") +
  theme_minimal()
ggsave(file.path(out_dir, "plots", "11_actual_vs_predicted_adjusted_price.png"), p_pred_price, width = 9, height = 5)

p_sigma <- ggplot(forecast_df, aes(x = date)) +
  geom_line(aes(y = abs(actual_return)), linewidth = 0.35) +
  geom_line(aes(y = predicted_sigma), linewidth = 0.45, linetype = "dashed") +
  labs(title = "Predicted Volatility vs Absolute Actual Return on Test Set",
       x = "Date", y = "Return Magnitude / Conditional Sigma") +
  theme_minimal()
ggsave(file.path(out_dir, "plots", "12_predicted_volatility_vs_abs_return.png"), p_sigma, width = 9, height = 5)

# SECTION 8: future forecast --------------------------------------------------
next_n <- 10
future_fc <- ugarchforecast(final_fit, n.ahead = next_n)
future_mu <- as.numeric(fitted(future_fc))
future_sigma <- as.numeric(sigma(future_fc))
last_price <- tail(gld_returns$adjusted, 1)

# Create approximate next trading dates by skipping weekends.
candidate_dates <- seq.Date(max(gld_returns$date) + 1, by = "day", length.out = 30)
future_dates <- candidate_dates[!weekdays(candidate_dates) %in% c("Saturday", "Sunday")][1:next_n]

future_forecast <- tibble(
  forecast_day = 1:next_n,
  approximate_trading_date = future_dates,
  predicted_return_percent = future_mu,
  predicted_sigma_percent = future_sigma,
  predicted_adjusted_price = last_price * exp(cumsum(future_mu / 100))
)

write.csv(future_forecast, file.path(out_dir, "future_10_day_forecast.csv"), row.names = FALSE)
print(future_forecast)
