# Simulation study for forecast-combination methods under model uncertainty
#
# This script implements the Chapter 3 methodology:
# - Sample sizes: 100, 200, 500, 1000
# - Monte Carlo replications: 50
# - Forecast horizons: 1, 3, 6, 12
# - Regimes: linear AR(2), nonlinear threshold AR, structural break, heteroskedastic AR-GARCH
# - Methods: Simple Average, Bates-Granger, Granger-Ramanathan,
#   Bayesian/Dynamic Model Averaging, Super Learner, XGBoost Stacking
#
# Output:
# - ranked mean metric tables printed directly in R, grouped by regime and sample size.
# - Diebold-Mariano p-value summaries printed and returned by sample size.
#
# The script uses base R. If the xgboost package is installed, it is used for
# XGBoost Stacking; otherwise, a small nonlinear boosted-stump fallback is used.

CONFIG <- list(
  sample_sizes = c(100, 300, 500,700, 1000),
  horizons = c(1, 3, 6, 12),
  regimes = c("linear_ar2", "nonlinear_threshold", "structural_break", "heteroskedastic"),
  replications = 50,
  seed = 20260629,
  train_fraction = 0.50,
  min_training_size = 50,
  min_meta_train = 25,
  bg_error_window = 50,
  dma_forgetting_factor = 0.97,
  evaluation_step = 1,
  xgb_nrounds = 40,
  xgb_eta = 0.05,
  ranking_metric = "mean_rmse",
  include_dm_tests = TRUE
)

METHODS <- c("SA", "BG", "GR", "BMA_DMA", "SL", "XGB_STACK")
BASE_MODELS <- c("naive", "mean", "ar1", "ar2", "ar4", "nonlinear_ar")
BASE_PARAM_COUNTS <- c(
  naive = 1,
  mean = 1,
  ar1 = 2,
  ar2 = 3,
  ar4 = 5,
  nonlinear_ar = 6
)

simulate_series <- function(regime, sample_size) {
  burn_in <- 200
  total <- sample_size + burn_in
  y <- numeric(total)
  
  if (regime == "linear_ar2") {
    eps <- rnorm(total, 0, 1)
    for (t in 3:total) {
      y[t] <- 0.60 * y[t - 1] + 0.20 * y[t - 2] + eps[t]
    }
  } else if (regime == "nonlinear_threshold") {
    eps <- rnorm(total, 0, 1)
    for (t in 3:total) {
      if (y[t - 1] <= 0) {
        y[t] <- 0.20 + 0.45 * y[t - 1] + eps[t]
      } else {
        y[t] <- -0.15 + 0.80 * y[t - 1] - 0.25 * y[t - 2] + eps[t]
      }
    }
  } else if (regime == "structural_break") {
    eps <- rnorm(total, 0, 1)
    tau <- burn_in + floor(sample_size / 2)
    for (t in 2:total) {
      if (t <= tau) {
        phi0 <- 0.00
        phi1 <- 0.65
      } else {
        phi0 <- 0.60
        phi1 <- -0.35
      }
      y[t] <- phi0 + phi1 * y[t - 1] + eps[t]
    }
  } else if (regime == "heteroskedastic") {
    omega <- 0.10
    alpha <- 0.12
    beta <- 0.82
    sigma2 <- rep(omega / (1 - alpha - beta), total)
    eps <- numeric(total)
    z <- rnorm(total, 0, 1)
    for (t in 2:total) {
      sigma2[t] <- omega + alpha * eps[t - 1]^2 + beta * sigma2[t - 1]
      eps[t] <- sqrt(max(sigma2[t], 1e-10)) * z[t]
      y[t] <- 0.55 * y[t - 1] + eps[t]
    }
  } else {
    stop("Unknown regime: ", regime)
  }
  
  y[(burn_in + 1):total]
}

feature_vector <- function(y, s, model_name) {
  if (model_name == "ar1") {
    c(1, y[s])
  } else if (model_name == "ar2") {
    c(1, y[s], y[s - 1])
  } else if (model_name == "ar4") {
    c(1, y[s], y[s - 1], y[s - 2], y[s - 3])
  } else if (model_name == "nonlinear_ar") {
    lag1 <- y[s]
    lag2 <- y[s - 1]
    c(1, lag1, lag2, lag1^2, lag1 * lag2, as.numeric(lag1 > 0))
  } else {
    stop("No feature vector for model: ", model_name)
  }
}

minimum_lag <- function(model_name) {
  if (model_name == "ar1") return(1)
  if (model_name == "ar2") return(2)
  if (model_name == "ar4") return(4)
  if (model_name == "nonlinear_ar") return(2)
  stop("Unknown model: ", model_name)
}

ridge_coef <- function(X, y, ridge = 1e-8) {
  penalty <- diag(ridge, ncol(X))
  penalty[1, 1] <- 0
  solve_result <- tryCatch(
    solve(crossprod(X) + penalty, crossprod(X, y)),
    error = function(e) qr.solve(crossprod(X) + penalty, crossprod(X, y))
  )
  as.numeric(solve_result)
}

expanding_ols_forecasts <- function(y, origins, h, model_name, min_obs = 20) {
  p <- minimum_lag(model_name)
  valid_s <- p:(length(y) - h)
  X_all <- t(vapply(valid_s, function(s) feature_vector(y, s, model_name), numeric(BASE_PARAM_COUNTS[[model_name]])))
  target <- y[valid_s + h]
  
  forecasts <- numeric(length(origins))
  for (i in seq_along(origins)) {
    origin <- origins[i]
    last_training_s <- origin - h
    n_train <- sum(valid_s <= last_training_s)
    min_needed <- max(min_obs, ncol(X_all) + 2)
    if (n_train < min_needed) {
      forecasts[i] <- y[origin]
    } else {
      X_train <- X_all[seq_len(n_train), , drop = FALSE]
      y_train <- target[seq_len(n_train)]
      beta <- ridge_coef(X_train, y_train)
      forecasts[i] <- sum(feature_vector(y, origin, model_name) * beta)
    }
  }
  forecasts
}

make_base_forecasts <- function(y, h, initial_train_size, evaluation_step) {
  origins <- seq(initial_train_size, length(y) - h, by = evaluation_step)
  actual <- y[origins + h]
  Fmat <- matrix(NA_real_, nrow = length(origins), ncol = length(BASE_MODELS))
  colnames(Fmat) <- BASE_MODELS
  
  Fmat[, "naive"] <- y[origins]
  cumulative_sum <- cumsum(y)
  Fmat[, "mean"] <- cumulative_sum[origins] / origins
  Fmat[, "ar1"] <- expanding_ols_forecasts(y, origins, h, "ar1")
  Fmat[, "ar2"] <- expanding_ols_forecasts(y, origins, h, "ar2")
  Fmat[, "ar4"] <- expanding_ols_forecasts(y, origins, h, "ar4")
  Fmat[, "nonlinear_ar"] <- expanding_ols_forecasts(y, origins, h, "nonlinear_ar")
  
  list(origins = origins, actual = actual, forecasts = Fmat)
}

project_to_simplex <- function(v) {
  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u)
  rho_candidates <- which(u * seq_along(u) > (cssv - 1))
  if (length(rho_candidates) == 0) return(rep(1 / length(v), length(v)))
  rho <- max(rho_candidates)
  theta <- (cssv[rho] - 1) / rho
  w <- pmax(v - theta, 0)
  if (sum(w) <= 0) rep(1 / length(v), length(v)) else w / sum(w)
}

super_learner_weights <- function(X, y, max_iter = 300) {
  k <- ncol(X)
  if (nrow(X) < k + 2) return(rep(1 / k, k))
  w <- rep(1 / k, k)
  gram <- crossprod(X) / nrow(X)
  lipschitz <- 2 * max(norm(gram, type = "2"), 1e-8)
  step <- 1 / lipschitz
  for (iter in seq_len(max_iter)) {
    grad <- -2 * as.numeric(crossprod(X, y - as.numeric(X %*% w))) / nrow(X)
    new_w <- project_to_simplex(w - step * grad)
    if (max(abs(new_w - w)) < 1e-8) {
      w <- new_w
      break
    }
    w <- new_w
  }
  w
}

fit_tiny_boosted_stumps <- function(X, y, nrounds = 40, eta = 0.05) {
  prediction <- rep(mean(y), length(y))
  stumps <- vector("list", nrounds)
  
  for (m in seq_len(nrounds)) {
    residual <- y - prediction
    best_sse <- Inf
    best <- NULL
    
    for (j in seq_len(ncol(X))) {
      thresholds <- unique(as.numeric(quantile(X[, j], probs = seq(0.1, 0.9, by = 0.1), names = FALSE)))
      for (threshold in thresholds) {
        left <- X[, j] <= threshold
        if (sum(left) == 0 || sum(!left) == 0) next
        left_value <- mean(residual[left])
        right_value <- mean(residual[!left])
        fitted <- ifelse(left, left_value, right_value)
        sse <- sum((residual - fitted)^2)
        if (sse < best_sse) {
          best_sse <- sse
          best <- list(feature = j, threshold = threshold, left = left_value, right = right_value)
        }
      }
    }
    
    if (is.null(best)) break
    stumps[[m]] <- best
    prediction <- prediction + eta * ifelse(X[, best$feature] <= best$threshold, best$left, best$right)
  }
  
  list(base = mean(y), stumps = stumps[!vapply(stumps, is.null, logical(1))], eta = eta)
}

predict_tiny_boosted_stumps <- function(model, row) {
  pred <- model$base
  for (stump in model$stumps) {
    pred <- pred + model$eta * ifelse(row[stump$feature] <= stump$threshold, stump$left, stump$right)
  }
  as.numeric(pred)
}

xgb_stack_predict <- function(X, y, row, config) {
  if (requireNamespace("xgboost", quietly = TRUE)) {
    dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X), label = y)
    model <- xgboost::xgb.train(
      data = dtrain,
      nrounds = config$xgb_nrounds,
      verbose = 0,
      params = list(
        objective = "reg:squarederror",
        max_depth = 2,
        eta = config$xgb_eta,
        subsample = 0.90,
        colsample_bytree = 1.0,
        nthread = 1
      )
    )
    pred <- predict(model, xgboost::xgb.DMatrix(data = matrix(row, nrow = 1)))
    list(prediction = as.numeric(pred), engine = "xgboost")
  } else {
    model <- fit_tiny_boosted_stumps(X, y, config$xgb_nrounds, config$xgb_eta)
    list(prediction = predict_tiny_boosted_stumps(model, row), engine = "tiny_boosted_stumps")
  }
}

bates_granger_weights <- function(errors) {
  mse <- colMeans(errors^2)
  inv <- 1 / pmax(mse, 1e-10)
  inv / sum(inv)
}

bma_weights <- function(errors) {
  n <- nrow(errors)
  if (n < 3) return(rep(1 / ncol(errors), ncol(errors)))
  bics <- numeric(ncol(errors))
  for (j in seq_len(ncol(errors))) {
    sse <- sum(errors[, j]^2)
    sigma2 <- max(sse / n, 1e-10)
    bics[j] <- n * log(sigma2) + BASE_PARAM_COUNTS[[colnames(errors)[j]]] * log(n)
  }
  weights <- exp(-0.5 * (bics - min(bics)))
  weights / sum(weights)
}

forecast_errors_matrix <- function(actual, forecasts) {
  matrix(actual, nrow = nrow(forecasts), ncol = ncol(forecasts)) - forecasts
}

granger_ramanathan_predict <- function(X, y, row) {
  design <- cbind(Intercept = 1, X)
  beta <- ridge_coef(design, y, ridge = 1e-6)
  as.numeric(sum(c(1, row) * beta))
}

combine_forecasts <- function(Fmat, actual, regime, config) {
  n <- nrow(Fmat)
  k <- ncol(Fmat)
  predictions <- setNames(vector("list", length(METHODS)), METHODS)
  elapsed <- setNames(rep(0, length(METHODS)), METHODS)
  for (method in METHODS) predictions[[method]] <- numeric(n)
  
  dma_posterior <- rep(1 / k, k)
  dma_sse <- rep(1, k)
  dma_count <- rep(1, k)
  xgb_engine <- "not_used"
  
  for (i in seq_len(n)) {
    row <- Fmat[i, ]
    
    t0 <- proc.time()[["elapsed"]]
    predictions$SA[i] <- mean(row)
    elapsed["SA"] <- elapsed["SA"] + proc.time()[["elapsed"]] - t0
    
    has_history <- i > config$min_meta_train
    if (has_history) {
      X_hist <- Fmat[seq_len(i - 1), , drop = FALSE]
      y_hist <- actual[seq_len(i - 1)]
    }
    
    t0 <- proc.time()[["elapsed"]]
    if (has_history) {
      begin <- max(1, i - config$bg_error_window)
      errors <- forecast_errors_matrix(
        actual[begin:(i - 1)],
        Fmat[begin:(i - 1), , drop = FALSE]
      )
      predictions$BG[i] <- sum(row * bates_granger_weights(errors))
    } else {
      predictions$BG[i] <- predictions$SA[i]
    }
    elapsed["BG"] <- elapsed["BG"] + proc.time()[["elapsed"]] - t0
    
    t0 <- proc.time()[["elapsed"]]
    if (has_history) {
      predictions$GR[i] <- granger_ramanathan_predict(X_hist, y_hist, row)
    } else {
      predictions$GR[i] <- predictions$SA[i]
    }
    elapsed["GR"] <- elapsed["GR"] + proc.time()[["elapsed"]] - t0
    
    t0 <- proc.time()[["elapsed"]]
    if (regime == "structural_break") {
      prior <- dma_posterior^config$dma_forgetting_factor
      prior <- prior / sum(prior)
      predictions$BMA_DMA[i] <- sum(row * prior)
    } else if (has_history) {
      errors <- forecast_errors_matrix(y_hist, X_hist)
      predictions$BMA_DMA[i] <- sum(row * bma_weights(errors))
    } else {
      predictions$BMA_DMA[i] <- predictions$SA[i]
    }
    elapsed["BMA_DMA"] <- elapsed["BMA_DMA"] + proc.time()[["elapsed"]] - t0
    
    t0 <- proc.time()[["elapsed"]]
    if (has_history) {
      w_sl <- super_learner_weights(X_hist, y_hist)
      predictions$SL[i] <- sum(row * w_sl)
    } else {
      predictions$SL[i] <- predictions$SA[i]
    }
    elapsed["SL"] <- elapsed["SL"] + proc.time()[["elapsed"]] - t0
    
    t0 <- proc.time()[["elapsed"]]
    if (has_history) {
      xgb_result <- xgb_stack_predict(X_hist, y_hist, row, config)
      predictions$XGB_STACK[i] <- xgb_result$prediction
      xgb_engine <- xgb_result$engine
    } else {
      predictions$XGB_STACK[i] <- predictions$SA[i]
    }
    elapsed["XGB_STACK"] <- elapsed["XGB_STACK"] + proc.time()[["elapsed"]] - t0
    
    if (regime == "structural_break") {
      error <- actual[i] - row
      sigma2 <- pmax(dma_sse / dma_count, 1e-8)
      likelihood <- exp(-0.5 * error^2 / sigma2) / sqrt(2 * pi * sigma2)
      posterior <- prior * pmax(likelihood, 1e-300)
      if (!is.finite(sum(posterior)) || sum(posterior) <= 0) {
        dma_posterior <- rep(1 / k, k)
      } else {
        dma_posterior <- posterior / sum(posterior)
      }
      dma_sse <- dma_sse + error^2
      dma_count <- dma_count + 1
    }
  }
  
  list(predictions = predictions, elapsed = elapsed, xgb_engine = xgb_engine)
}

rmse <- function(errors) sqrt(mean(errors^2))
mae <- function(errors) mean(abs(errors))
mase <- function(errors, train_y) mean(abs(errors)) / max(mean(abs(diff(train_y))), 1e-10)

diebold_mariano <- function(errors_a, errors_b, h, loss = "squared") {
  if (loss == "absolute") {
    loss_a <- abs(errors_a)
    loss_b <- abs(errors_b)
  } else {
    loss_a <- errors_a^2
    loss_b <- errors_b^2
  }
  
  d <- loss_a - loss_b
  n <- length(d)
  mean_d <- mean(d)
  centered <- d - mean_d
  max_lag <- max(h - 1, 0)
  long_run_var <- sum(centered^2) / n
  
  if (max_lag > 0) {
    for (lag in seq_len(max_lag)) {
      gamma <- sum(centered[(lag + 1):n] * centered[1:(n - lag)]) / n
      weight <- 1 - lag / (max_lag + 1)
      long_run_var <- long_run_var + 2 * weight * gamma
    }
  }
  
  if (!is.finite(long_run_var) || long_run_var <= 1e-14) {
    dm_stat <- 0
    p_value <- 1
  } else {
    dm_stat <- mean_d / sqrt(long_run_var / n)
    p_value <- 2 * (1 - pnorm(abs(dm_stat)))
  }
  
  favored <- if (mean_d < 0) "method_a" else if (mean_d > 0) "method_b" else "tie"
  data.frame(
    loss = loss,
    mean_loss_diff = mean_d,
    dm_stat = dm_stat,
    p_value = p_value,
    significant_5pct = p_value < 0.05,
    favored_method = favored,
    stringsAsFactors = FALSE
  )
}

evaluate_replication <- function(y, initial_train_size, h, predictions, actual, elapsed,
                                 ids, xgb_engine, include_dm = FALSE) {
  errors_by_method <- lapply(predictions, function(pred) actual - pred)
  rmse_sa <- rmse(errors_by_method$SA)
  train_y <- y[seq_len(initial_train_size)]
  
  metric_rows <- do.call(rbind, lapply(METHODS, function(method) {
    errors <- errors_by_method[[method]]
    data.frame(
      regime = ids$regime,
      sample_size = ids$sample_size,
      horizon = ids$horizon,
      replication = ids$replication,
      initial_train_size = ids$initial_train_size,
      evaluation_origins = ids$evaluation_origins,
      method = method,
      rmse = rmse(errors),
      mae = mae(errors),
      mase = mase(errors, train_y),
      relative_rmse = rmse(errors) / max(rmse_sa, 1e-10),
      computation_time_seconds = as.numeric(elapsed[[method]]),
      xgb_stack_engine = ifelse(method == "XGB_STACK", xgb_engine, ""),
      stringsAsFactors = FALSE
    )
  }))
  
  dm_rows <- NULL
  if (include_dm) {
    dm_rows <- list()
    idx <- 1
    pairs <- combn(METHODS, 2, simplify = FALSE)
    for (pair in pairs) {
      for (loss_name in c("squared", "absolute")) {
        test <- diebold_mariano(
          errors_by_method[[pair[1]]],
          errors_by_method[[pair[2]]],
          h,
          loss = loss_name
        )
        dm_rows[[idx]] <- cbind(
          data.frame(
            regime = ids$regime,
            sample_size = ids$sample_size,
            horizon = ids$horizon,
            replication = ids$replication,
            method_a = pair[1],
            method_b = pair[2],
            stringsAsFactors = FALSE
          ),
          test
        )
        idx <- idx + 1
      }
    }
    dm_rows <- do.call(rbind, dm_rows)
  }
  
  list(metrics = metric_rows, dm = dm_rows)
}

summarize_metrics <- function(metric_rows) {
  group_cols <- c("regime", "sample_size", "horizon", "method")
  
  win_keys <- paste(metric_rows$regime, metric_rows$sample_size, metric_rows$horizon,
                    metric_rows$replication, sep = "|")
  metric_rows$win <- FALSE
  for (key in unique(win_keys)) {
    rows <- which(win_keys == key)
    best <- min(metric_rows$rmse[rows])
    metric_rows$win[rows] <- abs(metric_rows$rmse[rows] - best) <= 1e-12
  }
  
  group_keys <- do.call(paste, c(metric_rows[group_cols], sep = "|"))
  groups <- split(metric_rows, group_keys, drop = TRUE)
  summary_rows <- lapply(groups, function(dat) {
    data.frame(
      regime = dat$regime[1],
      sample_size = dat$sample_size[1],
      horizon = dat$horizon[1],
      method = dat$method[1],
      replications = nrow(dat),
      mean_rmse = mean(dat$rmse),
      sd_rmse = sd(dat$rmse),
      mean_mae = mean(dat$mae),
      sd_mae = sd(dat$mae),
      mean_mase = mean(dat$mase),
      mean_relative_rmse = mean(dat$relative_rmse),
      mean_computation_time_seconds = mean(dat$computation_time_seconds),
      rmse_win_rate = mean(dat$win),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, summary_rows)
}

rank_mean_metrics_by_regime <- function(metric_rows) {
  group_cols <- c("regime", "sample_size", "method")
  
  win_keys <- paste(metric_rows$regime, metric_rows$sample_size, metric_rows$horizon,
                    metric_rows$replication, sep = "|")
  metric_rows$win <- FALSE
  for (key in unique(win_keys)) {
    rows <- which(win_keys == key)
    best <- min(metric_rows$rmse[rows])
    metric_rows$win[rows] <- abs(metric_rows$rmse[rows] - best) <= 1e-12
  }
  
  group_keys <- do.call(paste, c(metric_rows[group_cols], sep = "|"))
  groups <- split(metric_rows, group_keys, drop = TRUE)
  
  ranked <- lapply(groups, function(dat) {
    data.frame(
      regime = dat$regime[1],
      sample_size = dat$sample_size[1],
      method = dat$method[1],
      horizons_used = paste(sort(unique(dat$horizon)), collapse = ", "),
      replications_per_setting = length(unique(dat$replication)),
      mean_rmse = mean(dat$rmse),
      mean_mae = mean(dat$mae),
      mean_mase = mean(dat$mase),
      mean_relative_rmse = mean(dat$relative_rmse),
      mean_computation_time_seconds = mean(dat$computation_time_seconds),
      rmse_win_rate = mean(dat$win),
      stringsAsFactors = FALSE
    )
  })
  
  ranked <- do.call(rbind, ranked)
  ranked <- ranked[order(ranked$regime, ranked$sample_size, ranked$mean_rmse,
                         ranked$mean_mae), ]
  rank_group <- paste(ranked$regime, ranked$sample_size, sep = "|")
  ranked$rank <- ave(
    ranked$mean_rmse,
    rank_group,
    FUN = function(x) rank(x, ties.method = "first")
  )
  ranked <- ranked[order(ranked$regime, ranked$sample_size, ranked$rank), ]
  rownames(ranked) <- NULL
  ranked[, c(
    "regime",
    "sample_size",
    "rank",
    "method",
    "mean_rmse",
    "mean_mae",
    "mean_mase",
    "mean_relative_rmse",
    "mean_computation_time_seconds",
    "rmse_win_rate",
    "horizons_used",
    "replications_per_setting"
  )]
}

print_ranked_regime_tables <- function(ranked_metrics) {
  for (regime_name in unique(ranked_metrics$regime)) {
    for (n_value in sort(unique(ranked_metrics$sample_size[ranked_metrics$regime == regime_name]))) {
      cat("\n")
      cat("============================================================\n")
      cat("Regime:", regime_name, "| Sample size n =", n_value, "\n")
      cat("Ranked by lowest mean RMSE\n")
      cat("============================================================\n")
      regime_table <- ranked_metrics[
        ranked_metrics$regime == regime_name & ranked_metrics$sample_size == n_value,
      ]
      print(regime_table, row.names = FALSE)
    }
  }
}

print_best_models_by_sample_size <- function(ranked_metrics) {
  best_models <- ranked_metrics[ranked_metrics$rank == 1, ]
  best_models <- best_models[order(best_models$regime, best_models$sample_size), ]
  
  cat("\n")
  cat("============================================================\n")
  cat("Best Model Under Each Regime and Sample Size\n")
  cat("============================================================\n")
  print(best_models[, c(
    "regime",
    "sample_size",
    "method",
    "mean_rmse",
    "mean_mae",
    "mean_mase",
    "mean_relative_rmse",
    "rmse_win_rate"
  )], row.names = FALSE)
  
  invisible(best_models)
}

print_dm_benchmark_summary <- function(dm_summary, benchmark = "SA") {
  if (is.null(dm_summary) || nrow(dm_summary) == 0) return(invisible(NULL))
  
  benchmark_rows <- dm_summary[
    dm_summary$loss == "squared" &
      (dm_summary$method_a == benchmark | dm_summary$method_b == benchmark),
  ]
  
  if (nrow(benchmark_rows) == 0) return(invisible(NULL))
  
  benchmark_rows$comparison_method <- ifelse(
    benchmark_rows$method_a == benchmark,
    benchmark_rows$method_b,
    benchmark_rows$method_a
  )
  benchmark_rows$share_favoring_comparison <- ifelse(
    benchmark_rows$method_a == benchmark,
    benchmark_rows$share_favoring_method_b,
    benchmark_rows$share_favoring_method_a
  )
  benchmark_rows$share_favoring_benchmark <- ifelse(
    benchmark_rows$method_a == benchmark,
    benchmark_rows$share_favoring_method_a,
    benchmark_rows$share_favoring_method_b
  )
  
  display <- benchmark_rows[, c(
    "regime",
    "sample_size",
    "horizon",
    "comparison_method",
    "mean_dm_stat",
    "mean_p_value",
    "rejection_rate_5pct",
    "share_favoring_comparison",
    "share_favoring_benchmark"
  )]
  display <- display[order(display$regime, display$sample_size, display$horizon,
                           display$comparison_method), ]
  
  cat("\n")
  cat("============================================================\n")
  cat("Diebold-Mariano Test Summary\n")
  cat("Benchmark:", benchmark, "| Loss: squared error\n")
  cat("============================================================\n")
  print(display, row.names = FALSE)
  invisible(display)
}

prepare_dm_pvalues_by_sample_size <- function(dm_summary, loss = "squared") {
  if (is.null(dm_summary) || nrow(dm_summary) == 0) return(NULL)
  
  dm_pvalues <- dm_summary[dm_summary$loss == loss, c(
    "regime",
    "sample_size",
    "horizon",
    "method_a",
    "method_b",
    "mean_dm_stat",
    "mean_p_value",
    "rejection_rate_5pct",
    "share_favoring_method_a",
    "share_favoring_method_b"
  )]
  
  dm_pvalues <- dm_pvalues[order(
    dm_pvalues$sample_size,
    dm_pvalues$regime,
    dm_pvalues$horizon,
    dm_pvalues$method_a,
    dm_pvalues$method_b
  ), ]
  rownames(dm_pvalues) <- NULL
  dm_pvalues
}

print_dm_pvalues_by_sample_size <- function(dm_pvalues) {
  if (is.null(dm_pvalues) || nrow(dm_pvalues) == 0) return(invisible(NULL))
  
  for (n_value in sort(unique(dm_pvalues$sample_size))) {
    cat("\n")
    cat("============================================================\n")
    cat("Diebold-Mariano P-values | Sample size n =", n_value, "\n")
    cat("Loss: squared error | Pairwise method tests\n")
    cat("============================================================\n")
    print(dm_pvalues[dm_pvalues$sample_size == n_value, ], row.names = FALSE)
  }
  
  invisible(dm_pvalues)
}

summarize_dm <- function(dm_rows) {
  group_cols <- c("regime", "sample_size", "horizon", "method_a", "method_b", "loss")
  group_keys <- do.call(paste, c(dm_rows[group_cols], sep = "|"))
  groups <- split(dm_rows, group_keys, drop = TRUE)
  summary_rows <- lapply(groups, function(dat) {
    data.frame(
      regime = dat$regime[1],
      sample_size = dat$sample_size[1],
      horizon = dat$horizon[1],
      method_a = dat$method_a[1],
      method_b = dat$method_b[1],
      loss = dat$loss[1],
      replications = nrow(dat),
      mean_dm_stat = mean(dat$dm_stat),
      mean_p_value = mean(dat$p_value),
      rejection_rate_5pct = mean(dat$significant_5pct),
      share_favoring_method_a = mean(dat$favored_method == "method_a"),
      share_favoring_method_b = mean(dat$favored_method == "method_b"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, summary_rows)
}

run_study <- function(config) {
  all_metrics <- list()
  all_dm <- list()
  row_id <- 1
  dm_id <- 1
  
  set.seed(config$seed)
  
  for (regime in config$regimes) {
    for (sample_size in config$sample_sizes) {
      initial_train_size <- max(round(config$train_fraction * sample_size), config$min_training_size)
      if (initial_train_size >= sample_size - max(config$horizons) - 5) {
        initial_train_size <- floor(sample_size / 2)
      }
      
      for (replication in seq_len(config$replications)) {
        seed_value <- config$seed +
          100000 * match(regime, config$regimes) +
          1000 * sample_size +
          replication
        set.seed(seed_value)
        y <- simulate_series(regime, sample_size)
        
        for (h in config$horizons) {
          if (initial_train_size >= length(y) - h - 2) next
          
          base <- make_base_forecasts(y, h, initial_train_size, config$evaluation_step)
          combined <- combine_forecasts(base$forecasts, base$actual, regime, config)
          ids <- list(
            regime = regime,
            sample_size = sample_size,
            horizon = h,
            replication = replication,
            initial_train_size = initial_train_size,
            evaluation_origins = length(base$origins)
          )
          
          evaluated <- evaluate_replication(
            y = y,
            initial_train_size = initial_train_size,
            h = h,
            predictions = combined$predictions,
            actual = base$actual,
            elapsed = combined$elapsed,
            ids = ids,
            xgb_engine = combined$xgb_engine,
            include_dm = isTRUE(config$include_dm_tests)
          )
          
          all_metrics[[row_id]] <- evaluated$metrics
          if (isTRUE(config$include_dm_tests) && !is.null(evaluated$dm)) {
            all_dm[[dm_id]] <- evaluated$dm
            dm_id <- dm_id + 1
          }
          row_id <- row_id + 1
        }
        
        message(sprintf(
          "completed regime=%s, T=%s, replication=%s/%s",
          regime, sample_size, replication, config$replications
        ))
      }
    }
  }
  
  detailed_metrics <- do.call(rbind, all_metrics)
  ranked_metrics <- rank_mean_metrics_by_regime(detailed_metrics)
  print_ranked_regime_tables(ranked_metrics)
  best_models <- print_best_models_by_sample_size(ranked_metrics)
  
  dm_pvalues <- NULL
  if (isTRUE(config$include_dm_tests) && length(all_dm) > 0) {
    dm_tests <- do.call(rbind, all_dm)
    dm_summary <- summarize_dm(dm_tests)
    dm_pvalues <- prepare_dm_pvalues_by_sample_size(dm_summary, loss = "squared")
    print_dm_pvalues_by_sample_size(dm_pvalues)
  }
  
  invisible(list(
    best_models_by_sample_size = best_models,
    ranked_metrics_by_sample_size = ranked_metrics,
    dm_pvalues_by_sample_size = dm_pvalues
  ))
}

# To do a very small smoke test before the full run, uncomment these lines:
# CONFIG$sample_sizes <- c(200)
# CONFIG$horizons <- c(1, 3)
# CONFIG$replications <- 2

results <- run_study(CONFIG)

#Saving the output in excel
library(writexl)

write_xlsx(
  list(
    best_models_by_sample_size = results$best_models_by_sample_size,
    ranked_metrics_by_sample_size = results$ranked_metrics_by_sample_size,
    dm_pvalues_by_sample_size = results$dm_pvalues_by_sample_size
  ),
  path = "study_results.xlsx"
)
