#parameters---------------------------------------------------------------------
dir.create("out", showWarnings = FALSE)

pilot <- FALSE   # TRUE = one quick pilot cell, FALSE = full simulation

if (pilot) {
  reps      <- 1
  n_group   <- c(16)   # per group sample size
  prob_miss <- c(.12)  # threshold percentile (of the control group)
} else {
  reps      <- 1000
  n_group   <- c(16, 32, 64, 128)
  prob_miss <- c(.02, .04, .08, .16, .32, .64)
}

# Design factors (each combination is one "condition" within every replication)
T_levels  <- c(3, 9)                       # few vs many measurement points
cv_levels <- c(moderate = 0.32219, high = 0.81773) # CV of the outcome (control group, mean 20)
mu_base   <- 20                            # control mean (same at all time points)
target_r2 <- 0.013                         # target delta R2 of the group x time interaction
rho_lags  <- c(1, .5, .25, .125)           # correlation at lag 0, 1, 2, 3
rho_floor <- .1                            # correlation at lag >= 4

methods_names <- c("Real Data, No Missings" , "Pairwise Deletion", "Listwise Deletion", "Threshold", "Threshold/2", "Threshold/sqrt2", "PMM MICE", "Tobit MICE", "Multilevel Tobit", "Lognormal Multilevel Tobit", "ROS", "Tobit Reg", "Tobit on Ranks")

library(mvtnorm)
library(tidyverse)
library(lme4)
library(lmerTest)
library(MuMIn)
library(mice)
library(qgcomp)
library(future)
library(future.callr)
library(future.apply)
library(censReg)
library(plm)
library(NADA2)
library(lmtest)
library(mitml)
library(maxLik)

source("functions.R")

#data specification------------------------------------------------------------
# Effect size of the interaction is chosen so that the group x time interaction
# explains `target` of the total variance. v = total variance of the outcome
# (raw scale), so it must be passed in whenever the CV changes.
get_delta <- function(T_, target = 0.025, v = 400/9) {
  k <- (T_ + 1) / (48 * (T_ - 1))   # V_int = k * delta^2
  a <- 1 / (16 * k)                 # V_group = a * V_int
  x <- target * v / (1 - target * (a + 2))   # required V_int
  sqrt(x / k)
}

# Correlation matrix: decays with lag (1, .5, .25, .125) and then stays at a floor
make_cor <- function(n_time, rho = rho_lags, floor = rho_floor) {
  lag <- abs(outer(seq_len(n_time), seq_len(n_time), "-"))
  R <- matrix(floor, n_time, n_time)
  idx <- lag < length(rho)
  R[idx] <- rho[lag[idx] + 1]
  R
}

# Everything needed to simulate one condition (number of time points x skew)
make_design <- function(n_time, cv, skew_label, mu = mu_base, target = target_r2) {
  v     <- (cv * mu)^2                       # raw-scale variance
  delta <- get_delta(n_time, target, v)
  s2    <- log(1 + cv^2)                     # log-scale variance (control group)
  list(skew_label        = skew_label,
       n_time            = n_time,
       cv                = cv,
       v                 = v,
       delta             = delta,
       skewness_control  = (exp(s2) + 2) * sqrt(exp(s2) - 1),  # lognormal skewness
       mean_control      = rep(mu, n_time),
       mean_intervention = mu + (0:(n_time - 1)) * delta / (n_time - 1),
       cov               = v * make_cor(n_time))
}

designs <- list()
for (s in names(cv_levels)) {
  for (tt in T_levels) {
    designs[[paste0(s, "_T", tt)]] <- make_design(tt, cv_levels[[s]], s)
  }
}
saveRDS(designs, file = "out/designs.rds")

#data generation----------------------------------------------------------------
# Lognormal data with the specified raw-scale means and covariance (moment matching)
rlnorm_mv <- function(n, m, S) {
  rmvnorm(n,
          mean  = log(m) - 0.5 * log(1 + diag(S) / m^2),
          sigma = log(1 + S / outer(m, m)))
}

# One long-format data set (group, ID, time, outcome) for one design
simulate_long <- function(design, n_group) {
  tp <- paste0("X", seq_len(design$n_time))
  w <- exp(rbind(rlnorm_mv(n_group, design$mean_intervention, design$cov),
                 rlnorm_mv(n_group, design$mean_control,      design$cov)))
  w <- as.data.frame(w)
  names(w) <- tp
  w$group <- c(rep("Intervention", n_group), rep("Control", n_group))
  w$ID    <- 1:(n_group * 2)
  pivot_longer(w, cols = all_of(tp), names_to = "time", values_to = "outcome") %>%
    mutate(time = factor(time, levels = tp), ID = factor(ID), group = factor(group))
}

#analysis helpers---------------------------------------------------------------
as_pr <- function(x) {
  x <- suppressWarnings(as.numeric(unlist(x)))
  if (length(x) != 2) return(c(p = NA_real_, rsq = NA_real_))
  c(p = x[1], rsq = x[2])
}

# LRT p value and delta marginal R2 for time * group vs time + group (lmer)
lmer_pr <- function(dat, label = "") {
  if (all(is.na(dat$outcome))) return(c(p = NA_real_, rsq = NA_real_))
  tryCatch({
    full <- suppressMessages(lmer(outcome ~ time * group + (1 | ID), dat))
    null <- suppressMessages(lmer(outcome ~ time + group + (1 | ID), dat))
    c(p   = anova(full, null)[[8]][[2]],
      rsq = (r.squaredGLMM(full) - r.squaredGLMM(null))[[1]])
  }, error = function(e) {
    message("lmer failed (", label, "): ", conditionMessage(e))
    c(p = NA_real_, rsq = NA_real_)
  })
}

# Replace censored values by lod * f and analyse
sub_lod <- function(l, lod, f, label) {
  l$outcome[is.na(l$outcome)] <- lod * f
  lmer_pr(l, label)
}

#analyse one data set with all methods----------------------------------------
analyse_dataset <- function(l_full, prob_miss) {
  library(qgcomp)   # makes mice.impute.tobit visible to mice

  #0 all data present
  res0 <- lmer_pr(l_full, "complete")

  #introduce missings (LOD = percentile of the control group)
  lod <- quantile(l_full$outcome[l_full$group == "Control"], probs = prob_miss, na.rm = TRUE)
  l <- l_full
  l$outcome[l$outcome < lod] <- NA

  #1a pairwise deletion (lmer drops the NA rows)
  res1a <- lmer_pr(l, "pairwise")

  #1b listwise deletion (drop every ID with at least one censored value)
  r <- l
  r$outcome[r$ID %in% unique(r$ID[is.na(r$outcome)])] <- NA
  res1b <- lmer_pr(r, "listwise")

  #2-4 substitution
  res2 <- sub_lod(l, lod, 1,         "LOD")
  res3 <- sub_lod(l, lod, 1/2,       "LOD/2")
  res4 <- sub_lod(l, lod, 1/sqrt(2), "LOD/sqrt2")

  #5 multiple imputation
  res5a <- safe_imp(impAnalysis_pmm, l, lod)
  res5b <- safe_imp(impAnalysis,     l, lod)
  res5c <- safe_imp(impAnalysis_2L,  l, lod)
  res5d <- safe_imp(impAnalysis_2L_log, l, lod)   # imputes on log scale, analyses on original scale

  #6 regression on order statistics
  res6 <- safe_ros(l, lod, "ROS")

  #7 Tobit regression
  pdata <- pdata.frame(l, index = c("ID", "time"))
  pdata$outcome[is.na(pdata$outcome)] <- lod
  res7 <- safe_tobit(pdata, "outcome", lod, "Tobit")

  #8 Tobit on ranks
  pdata$y_rank <- rank(pdata$outcome, ties.method = "average")
  res8 <- safe_tobit(pdata, "y_rank", min(pdata$y_rank), "TobitRank")

  out <- list("Real Data, No Missings" = res0,
              "Pairwise Deletion"      = res1a,
              "Listwise Deletion"      = res1b,
              "Threshold"              = res2,
              "Threshold/2"            = res3,
              "Threshold/sqrt2"        = res4,
              "PMM MICE"               = res5a,
              "Tobit MICE"             = res5b,
              "Multilevel Tobit"       = res5c,
              "Lognormal Multilevel Tobit" = res5d,
              "ROS"                    = res6,
              "Tobit Reg"              = res7,
              "Tobit on Ranks"         = res8)
  stopifnot(identical(names(out), methods_names))
  out <- t(sapply(out, as_pr))

  data.frame(method          = methods_names,
             p               = unname(out[, "p"]),
             rsq             = unname(out[, "rsq"]),
             prop_cens       = mean(is.na(l$outcome)),
             prop_cens_ctrl  = mean(is.na(l$outcome[l$group == "Control"])),
             prop_cens_int   = mean(is.na(l$outcome[l$group == "Intervention"])),
             n_ids_complete  = length(unique(l$ID)) - length(unique(l$ID[is.na(l$outcome)])),
             stringsAsFactors = FALSE)
}

#one replication: all conditions (skew x time points)---------------------------
getData <- function(n_group, prob_miss){
  out <- lapply(designs, function(d) {
    l <- simulate_long(d, n_group)
    r <- analyse_dataset(l, prob_miss)
    r$skew  <- d$skew_label
    r$cv    <- d$cv
    r$n_time <- d$n_time
    r$delta <- d$delta
    r
  })
  return(list("results" = do.call(rbind, out)))
}


res <- vector(mode = "list", length = length(prob_miss))
names(res) <- paste0(100*prob_miss, "%")
for (l in 1:length(res)) {
  res[[l]] <- vector(mode = "list", length = length(n_group))
  names(res[[l]]) <- n_group
}

#run analysis-------------------------------------------------------------------

{plan(callr, workers = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "128")))
  for (probs in 1:length(prob_miss)) {
    for (ns in 1:length(n_group)) {
      res[[probs]][[ns]] <- future_replicate(
        n = reps,
        expr = {
          assign("mice.impute.2l.tobit", mice.impute.2l.tobit, envir = .GlobalEnv)
          getData(n_group = n_group[ns], prob_miss = prob_miss[probs])
        },
        simplify = FALSE,
        future.seed = TRUE,
        future.gc = TRUE,
        future.globals = TRUE,
        future.packages = c("mice", "mvtnorm", "lme4", "lmerTest", "MuMIn", "mitml", "qgcomp", "tidyverse", "censReg", "plm", "NADA2", "lmtest")
      )
      log_progress("n_group", ns, length(n_group))
    }
    log_progress("prob_miss", probs, length(prob_miss))
    
    # restart workers
    plan(sequential)
    plan(callr, workers = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "128")))
    gc()
  }
}

saveRDS(res, file = "out/res.rds")
