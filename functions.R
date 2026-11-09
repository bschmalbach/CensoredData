#function definitions=========================================================

mice.impute.2l.tobit <- function(y, ry, x, type, lod, wy = NULL, log_scale = FALSE, ...) {
  if (is.null(wy)) wy <- !ry
  N <- length(y)
  if (sum(wy) == 0) return(numeric(0))
  
  # --- LOD per row -----------------------------------------------------------
  if (length(lod) == 1) {
    lod_all <- rep(lod, N)
  } else if (length(lod) == N) {
    lod_all <- lod
  } else {
    stop("lod must have length 1 or the same length as y.")
  }
  
  # --- cluster variable ------------------------------------------------------
  cluster_col <- which(type == -2)
  if (length(cluster_col) != 1) {
    stop("Exactly one cluster variable (type = -2) must be specified.")
  }
  cl_f <- factor(as.character(x[, cluster_col]))
  idx  <- as.integer(cl_f)
  J    <- nlevels(cl_f)
  n_i  <- tabulate(idx, nbins = J)
  
  # --- fixed-effect design matrix (intercept added, aliased columns dropped) --
  fixed_cols <- which(type == 1)
  if (length(fixed_cols) == 0) stop("No fixed-effect predictors (type = 1) specified.")
  x_fixed <- as.matrix(x[, fixed_cols, drop = FALSE])
  if (!is.null(colnames(x_fixed))) {
    x_fixed <- x_fixed[, !(colnames(x_fixed) %in% c("(Intercept)", "Intercept")),
                       drop = FALSE]
  }
  X <- cbind(`(Intercept)` = 1, x_fixed)
  if (anyNA(X)) stop("Fixed-effect predictors must be complete.")
  qrX  <- qr(X, tol = 1e-7)
  keep <- sort(qrX$pivot[seq_len(qrX$rank)])
  X    <- X[, keep, drop = FALSE]
  p    <- ncol(X)
  
  # --- completed outcome: y already contains the current imputations ---------
  na_y <- which(is.na(y))
  # safety net only (log_scale = TRUE: y and lod are on the log scale, so LOD/sqrt(2) becomes a subtraction)
  if (length(na_y)) y[na_y] <- if (log_scale) lod_all[na_y] - log(sqrt(2)) else lod_all[na_y] / sqrt(2)
  
  # --- 1. variance components from the completed data ------------------------
  dat <- as.data.frame(X[, -1, drop = FALSE])
  names(dat) <- make.names(names(dat), unique = TRUE)
  rhs <- if (ncol(dat) > 0) paste(sprintf("`%s`", names(dat)), collapse = " + ") else "1"
  dat$y       <- y
  dat$cluster <- cl_f
  fit <- suppressMessages(suppressWarnings(
    lme4::lmer(stats::as.formula(paste0("y ~ ", rhs, " + (1 | cluster)")),
               data = dat, REML = FALSE,
               control = lme4::lmerControl(optimizer = "bobyqa"))
  ))
  
  sigma2_hat <- stats::sigma(fit)^2
  tau2_hat   <- as.numeric(lme4::VarCorr(fit)$cluster[1, 1])   # already a variance
  if (!is.finite(sigma2_hat) || sigma2_hat <= 0) stop("Residual variance is invalid.")
  tau2_hat <- max(tau2_hat, 1e-8)
  
  # --- 2. draw variance components (scaled inverse chi-square, approximate) --
  df_sigma <- max(N - J - (p - 1), 1)
  df_tau   <- max(J - 1, 1)
  sigma2 <- sigma2_hat * df_sigma / stats::rchisq(1, df_sigma)
  tau2   <- tau2_hat   * df_tau   / stats::rchisq(1, df_tau)
  
  # --- 3. draw fixed effects from their GLS conditional posterior ------------
  #   beta | sigma2, tau2, y ~ N(beta_gls, (X' S^-1 X)^-1),
  #   S block-diagonal with sigma2*I + tau2*11'
  w_i <- tau2 / (sigma2 + n_i * tau2)
  Xs  <- rowsum(X, idx)                       # J x p cluster sums
  ys  <- as.numeric(rowsum(y, idx))           # J cluster sums
  A   <- (crossprod(X) - crossprod(Xs, w_i * Xs)) / sigma2
  bv  <- (crossprod(X, y) - crossprod(Xs, w_i * ys)) / sigma2
  V   <- solve(A)
  V   <- (V + t(V)) / 2
  beta_gls <- as.numeric(V %*% bv)
  L <- tryCatch(chol(V), error = function(e) NULL)
  beta <- if (is.null(L)) beta_gls else beta_gls + as.numeric(t(L) %*% stats::rnorm(p))
  
  # --- 4. draw random intercepts given the completed data --------------------
  xb <- as.numeric(X %*% beta)
  rs <- as.numeric(rowsum(y - xb, idx))
  v_i <- 1 / (1 / tau2 + n_i / sigma2)
  m_i <- v_i * rs / sigma2
  b   <- stats::rnorm(J, mean = m_i, sd = sqrt(v_i))
  
  # --- 5. draw censored cells from N(mu, sigma2) truncated to (0, LOD) ----
  mu_mis  <- (xb + b[idx])[wy]
  lod_mis <- lod_all[wy]
  sd_draw <- sqrt(sigma2)
  p_lo <- stats::pnorm(lod_mis, mean = mu_mis, sd = sd_draw)
  p_lo <- pmin(pmax(p_lo, 1e-12), 1 - 1e-12)
  u    <- pmax(stats::runif(sum(wy), min = 0, max = p_lo), 1e-300)
  imp  <- stats::qnorm(u, mean = mu_mis, sd = sd_draw)
  return(pmin(imp, lod_mis - 1e-8))
}

# Pool the R2 *difference*: compute Delta R2 within each imputation, then average
# (Rubin's point estimate). Same scale as the complete-data benchmark.
pool_delta_r2 <- function(r2_full, r2_null) {
  r2_full <- suppressWarnings(as.numeric(unlist(r2_full)))
  r2_null <- suppressWarnings(as.numeric(unlist(r2_null)))
  if (length(r2_full) != length(r2_null)) stop("pool_delta_r2: length mismatch")
  d <- r2_full - r2_null
  d <- d[is.finite(d)]
  if (length(d) == 0) stop("pool_delta_r2: no valid R2 values")
  mean(d)
}
# Wide (ID, group, X1..Xk) -> long analysis data. Takes the data frame as an
# explicit argument, so no reliance on how with.mids() scopes `data`.
to_long <- function(w) {
  tidyr::pivot_longer(w,
                      cols = tidyselect::starts_with("X"),
                      names_to = "time",
                      values_to = "outcome",
                      names_prefix = "X") %>%
    dplyr::mutate(time = factor(as.numeric(time)),
                  ID = factor(ID),
                  group = factor(group))
}

# Marginal R2 (first element of r.squaredGLMM) as a plain number; works whether
# MuMIn returns a 1x2 matrix or a named vector.
r2m_of <- function(f) as.numeric(MuMIn::r.squaredGLMM(f))[1]

# Fit one lmer per completed data set and return a mira object (for mice::D1).
fit_mira <- function(imp, formula) {
  mice::as.mira(lapply(seq_len(imp$m), function(i) {
    lme4::lmer(formula, data = to_long(mice::complete(imp, i)), REML = FALSE)
  }))
}

# Time x treatment interaction as a FACTOR interaction, i.e. one dummy per
# non-reference time level (trt * [time == level]). This mirrors the analysis
# model (outcome ~ time * group), where time is a factor. `time` itself enters
# the imputation model as a factor (mice expands it into dummies).
add_time_trt <- function(data) {
  lv <- if (is.factor(data$time)) levels(data$time) else sort(unique(as.character(data$time)))
  cols <- character(0)
  for (l in lv[-1]) {
    nm <- paste0("trt_", l)
    data[[nm]] <- as.numeric(as.character(data$time) == l) * data$trt
    cols <- c(cols, nm)
  }
  list(data = data, cols = cols)
}

impAnalysis_2L <- function(data, lod) {
  data <- as.data.frame(data) 
  
  # 1. Prepare treatment variable
  data$group <- factor(data$group)
  if (nlevels(data$group) != 2) {
    stop("This function currently assumes exactly two treatment groups.")
  }
  data$trt <- as.numeric(data$group == levels(data$group)[2])
  
  #2. Check ID
  data$ID <- as.numeric(data$ID)
  
  # 3. Create time x treatment interaction (factor interaction, one dummy per time level)
  tt <- add_time_trt(data)
  data <- tt$data
  int_cols <- tt$cols
  
  # 4. Predictor and type matrix
  pred <- mice::make.predictorMatrix(data)
  pred[, "ID"] <- -2
  pred["outcome", c("time", "trt", int_cols)] <- 1
  
  # 5. Imputation methods
  meth <- mice::make.method(data)
  meth["outcome"] <- "2l.tobit"
  meth[c("ID", "group", "trt", "time", int_cols)] <- ""
  
  # 6. Multiple imputation
  ret <- tryCatch({
    data_init <- data
    data_init$outcome[is.na(data_init$outcome)] <- lod / sqrt(2)
    
    imp <- mice(
      data, m = 10, method = meth, predictorMatrix = pred,
      maxit = 10, printFlag = FALSE, lod = lod, data.init = data_init
    )
    
    # 7. Fit analysis models on each completed (long) data set
    m <- imp$m
    fits_full <- vector("list", m)
    fits_null <- vector("list", m)
    
    for (i in seq_len(m)) {
      dat_i <- mice::complete(imp, i)
      dat_i$ID    <- factor(dat_i$ID)
      dat_i$group <- factor(dat_i$group)
      
      # Full model: time * group + (1 | ID)
      fits_full[[i]] <- lme4::lmer(outcome ~ time * group + (1 | ID),
                                   data = dat_i, REML = FALSE)
      # Null model: time + group + (1 | ID)
      fits_null[[i]] <- lme4::lmer(outcome ~ time + group + (1 | ID),
                                   data = dat_i, REML = FALSE)
    }
    
    # 8. D1 test
    test <- mitml::testModels(fits_full, fits_null, method = "D1")
    pval <- as.numeric(test$test[1, "P(>F)"])
    if (length(pval) != 1 || is.na(pval)) stop("Could not extract D1 p-value")
    
    # 9. Pooled R2 difference
    rsq_full <- vapply(fits_full, r2m_of, numeric(1))
    rsq_null <- vapply(fits_null, r2m_of, numeric(1))
    rsq_diff <- pool_delta_r2(rsq_full, rsq_null)
    
    return(c(pval = pval, rsq_diff = rsq_diff))
    
  }, error = function(e) {
    message("ERROR in impAnalysis_2L(): ", e$message)
    c(pval = NA_real_, rsq_diff = NA_real_)
  })
  return(ret)
}

# Lognormal multilevel Tobit MI: impute on the log scale (log(outcome), log(LOD)),
# back-transform every completed data set with exp(), then fit the SAME analysis
# models on the original scale as all other methods
impAnalysis_2L_log <- function(data, lod) {
  data <- as.data.frame(data)
  
  # 1. Prepare treatment variable
  data$group <- factor(data$group)
  if (nlevels(data$group) != 2) {
    stop("This function currently assumes exactly two treatment groups.")
  }
  data$trt <- as.numeric(data$group == levels(data$group)[2])
  
  # 2. Check ID
  data$ID <- as.numeric(data$ID)
  
  # 3. Create time x treatment interaction (factor interaction, one dummy per time level)
  tt <- add_time_trt(data)
  data <- tt$data
  int_cols <- tt$cols
  
  # 4. Move outcome and LOD to the log scale (NAs stay NA)
  lod_log <- log(lod)
  data$outcome <- log(data$outcome)
  
  # 5. Predictor and type matrix
  pred <- mice::make.predictorMatrix(data)
  pred[, "ID"] <- -2
  pred["outcome", c("time", "trt", int_cols)] <- 1
  
  # 6. Imputation methods
  meth <- mice::make.method(data)
  meth["outcome"] <- "2l.tobit"
  meth[c("ID", "group", "trt", "time", int_cols)] <- ""
  
  # 7. Multiple imputation on the log scale
  ret <- tryCatch({
    data_init <- data
    data_init$outcome[is.na(data_init$outcome)] <- lod_log - log(sqrt(2))
    
    imp <- mice(
      data, m = 10, method = meth, predictorMatrix = pred,
      maxit = 10, printFlag = FALSE, lod = lod_log, log_scale = TRUE,
      data.init = data_init
    )
    
    # 8. Back-transform and fit analysis models on the ORIGINAL scale
    m <- imp$m
    fits_full <- vector("list", m)
    fits_null <- vector("list", m)
    
    for (i in seq_len(m)) {
      dat_i <- mice::complete(imp, i)
      dat_i$ID      <- factor(dat_i$ID)
      dat_i$group   <- factor(dat_i$group)
      dat_i$outcome <- exp(dat_i$outcome)
      if (anyNA(dat_i$outcome) || any(!is.finite(dat_i$outcome))) {
        stop("Non-finite outcome after back-transformation")
      }
      
      fits_full[[i]] <- lme4::lmer(outcome ~ time * group + (1 | ID),
                                   data = dat_i, REML = FALSE)
      fits_null[[i]] <- lme4::lmer(outcome ~ time + group + (1 | ID),
                                   data = dat_i, REML = FALSE)
    }
    
    # 9. D1 test
    test <- mitml::testModels(fits_full, fits_null, method = "D1")
    pval <- as.numeric(test$test[1, "P(>F)"])
    if (length(pval) != 1 || is.na(pval)) stop("Could not extract D1 p-value")
    
    # 10. Pooled R2 difference
    rsq_full <- vapply(fits_full, r2m_of, numeric(1))
    rsq_null <- vapply(fits_null, r2m_of, numeric(1))
    rsq_diff <- pool_delta_r2(rsq_full, rsq_null)
    
    return(c(pval = pval, rsq_diff = rsq_diff))
    
  }, error = function(e) {
    message("ERROR in impAnalysis_2L_log(): ", e$message)
    c(pval = NA_real_, rsq_diff = NA_real_)
  })
  return(ret)
}

impAnalysis <- function(data, lod) {
  data <- tidyr::pivot_wider(data, id_cols = c("ID", "group"), 
                             names_from = "time", values_from = "outcome")
  
  pred <- mice::make.predictorMatrix(data)
  pred[, "ID"] <- 0
  meth <- mice::make.method(data)
  meth[c("ID", "group")] <- ""
  meth[-c(1, 2)] <- "tobit"
  
  ret <- tryCatch({
    data_init <- data
    for (cl in setdiff(names(data), c("ID", "group")))
      data_init[[cl]][is.na(data_init[[cl]])] <- lod / sqrt(2)
    imp <- suppressWarnings(mice(data, m = 10, method = meth, predictorMatrix = pred,
                                 lod = lod, maxit = 10, printFlag = FALSE, data.init=data_init))
    
    # Fit full model on all imputations (returns mira object)
    fit_full <- fit_mira(imp, outcome ~ time * group + (1 | ID))
    
    # Fit null model on all imputations
    fit_null <- fit_mira(imp, outcome ~ time + group + (1 | ID))
    
    # D1 test (approximate p-value at [1, 2])
    test <- mitml::testModels(fit_full$analyses, fit_null$analyses, method = "D1")
    pval <- as.numeric(test$test[1, "P(>F)"])
    if (length(pval) != 1 || is.na(pval)) stop("Could not extract D1 p-value")
    
    # Extract marginal R² from already-fitted models
    rsq_full <- vapply(fit_full$analyses, r2m_of, numeric(1))
    rsq_null <- vapply(fit_null$analyses, r2m_of, numeric(1))
    rsq_diff <- pool_delta_r2(rsq_full, rsq_null)
    
    c(pval = pval, rsq_diff = rsq_diff)
    
  }, error = function(e) {
    message("ERROR in impAnalysis(): ", e$message)
    c(NA_real_, NA_real_)
  })
  return(ret)
}

impAnalysis_pmm <- function(data, lod) {
  data <- tidyr::pivot_wider(data, id_cols = c("ID", "group"), 
                             names_from = "time", values_from = "outcome")
  
  pred <- mice::make.predictorMatrix(data)
  pred[, "ID"] <- 0
  meth <- mice::make.method(data)
  meth[c("ID", "group")] <- ""
  meth[-c(1, 2)] <- "pmm"
  
  ret <- tryCatch({
    data_init <- data
    for (cl in setdiff(names(data), c("ID", "group")))
      data_init[[cl]][is.na(data_init[[cl]])] <- lod / sqrt(2)
    imp <- mice(data, m = 10, method = meth, predictorMatrix = pred,
                maxit = 10, printFlag = FALSE, data.init=data_init)
    
    fit_full <- fit_mira(imp, outcome ~ time * group + (1 | ID))
    
    fit_null <- fit_mira(imp, outcome ~ time + group + (1 | ID))
    
    test <- mitml::testModels(fit_full$analyses, fit_null$analyses, method = "D1")
    pval <- as.numeric(test$test[1, "P(>F)"])
    if (length(pval) != 1 || is.na(pval)) stop("Could not extract D1 p-value")
    
    rsq_full <- vapply(fit_full$analyses, r2m_of, numeric(1))
    rsq_null <- vapply(fit_null$analyses, r2m_of, numeric(1))
    rsq_diff <- pool_delta_r2(rsq_full, rsq_null)
    
    c(pval = pval, rsq_diff = rsq_diff)
    
  }, error = function(e) {
    message("ERROR in impAnalysis_pmm(): ", e$message)
    c(NA_real_, NA_real_)
  })
  return(ret)
}

safe_rbind_results <- function(...) {
  do.call(rbind, lapply(list(...), function(x) {
    x <- suppressWarnings(as.numeric(unlist(x)))
    if (length(x) != 2) return(c(NA_real_, NA_real_))
    x
  }))
}

safe_imp <- function(f, data, lod) {
  tryCatch({
    res <- f(data, lod)
    if (length(res) == 2) res else c(NA_real_, NA_real_)
  }, error = function(e) {
    message("imputation failed: ", e$message)
    c(NA_real_, NA_real_)
  })
}

mckelvey_zavoina <- function(fit) {
  cf <- coef(fit, logSigma = FALSE)
  
  if ("sigmaMu" %in% names(cf)) {
    # panel random-effects
    beta   <- cf[seq_len(length(cf) - 2)]
    sigma2 <- cf["sigmaMu"]^2 + cf["sigmaNu"]^2
  } else {
    # cross-sectional Tobit
    beta   <- cf[seq_len(length(cf) - 1)]
    sigma2 <- cf["sigma"]^2
  }
  
  mf <- model.frame(fit)
  X  <- model.matrix(formula(fit), data = mf)
  X  <- X[, names(beta), drop = FALSE]
  yhat   <- as.numeric(X %*% beta)
  var_yh <- var(yhat)
  var_yh / (var_yh + sigma2)
}

rosImp <- function(dat, lod){
  dat$id <- seq_len(nrow(dat))
  dat$censored <- ifelse(is.na(dat$outcome), T, F)
  dat$outcome[dat$censored] <- lod
  r <- ros(outcome ~ censored, data = dat)
  ros_order <- order(dat$outcome, dat$censored)
  ros_results <- data.frame(id = dat$id[ros_order], ros_value = r$modeled)
  return(merge(dat, ros_results, by = "id", sort = FALSE))
}

# Tobit LR test + McKelvey-Zavoina R2 difference for one dataset.
# `response` is a column name in `data` (so ranks must be stored as a column).
safe_tobit <- function(data, response, left, label = "") {
  tryCatch({
    f_full <- stats::as.formula(paste(response, "~ group * time"))
    f_null <- stats::as.formula(paste(response, "~ group + time"))
    
    fit_full <- censReg::censReg(f_full, data = data, left = left, method = "BFGS")
    fit_null <- censReg::censReg(f_null, data = data, left = left, method = "BFGS")
    
    ok_codes <- c(0, 1, 2, 3, 8)
    rc_full <- maxLik::returnCode(fit_full); rc_null <- maxLik::returnCode(fit_null)
    if (!(rc_full %in% ok_codes) || !(rc_null %in% ok_codes)) {
      stop(sprintf("no convergence (codes %d/%d: %s | %s)", rc_full, rc_null,
                   maxLik::returnMessage(fit_full), maxLik::returnMessage(fit_null)))
    }
    
    p   <- as.numeric(lmtest::lrtest(fit_full, fit_null)[2, 5])
    rsq <- as.numeric(mckelvey_zavoina(fit_full) - mckelvey_zavoina(fit_null))
    
    if (length(p) != 1 || length(rsq) != 1) stop("unexpected output length")
    c(p = p, rsq = rsq)
  }, error = function(e) {
    message("safe_tobit failed (", label, "): ", conditionMessage(e))
    c(p = NA_real_, rsq = NA_real_)
  })
}

# ROS imputation + mixed model for one dataset.
safe_ros <- function(dat, lod, label = "") {
  tryCatch({
    r    <- rosImp(dat, lod = lod)
    full <- suppressMessages(lmer(ros_value ~ time * group + (1 | ID), r))
    null <- suppressMessages(lmer(ros_value ~ time + group + (1 | ID), r))
    
    p   <- as.numeric(suppressMessages(anova(full, null))[[8]][[2]])
    rsq <- as.numeric((MuMIn::r.squaredGLMM(full) - MuMIn::r.squaredGLMM(null))[[1]])
    
    if (length(p) != 1 || length(rsq) != 1) stop("unexpected output length")
    c(p = p, rsq = rsq)
  }, error = function(e) {
    message("safe_ros failed (", label, "): ", conditionMessage(e))
    c(p = NA_real_, rsq = NA_real_)
  })
}

log_progress <- function(label, current, total, file = "progress.txt") {
  
  mem_mb <- round(sum(gc()[, 2]), 1)
  
  line <- sprintf(
    "[%s] %s: %d / %d | Memory: %.1f MB\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    label,
    current,
    total,
    mem_mb
  )
  cat(line, file = file, append = TRUE)
}