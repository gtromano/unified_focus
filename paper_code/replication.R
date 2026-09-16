## Replication script (R) for
## "focus and focus-cpt: Fast Online Changepoint Detection in R and Python"
## (Journal of Statistical Software)
##
## This file is generated from the code chunks of jss_paper.qmd by
## make_replication_scripts.py (make replication-scripts): do not edit by hand.
##
## This script reproduces all R results and figures of the manuscript, in their
## order of appearance: Section 3 (the interface), Section 4 (common use cases)
## and Section 5.1 (custom cost function on the NBA data). The Python results
## (Section 3, Sections 5.2 and 5.3, and Appendix A) are reproduced by the
## companion script replication.py.
##
## Requirements: R (>= 4.1.0) and the packages focus (>= 0.1.10), ggplot2,
## furrr, purrr and dplyr, all available from CRAN. The focus package can also
## be installed from the source package submitted with the manuscript:
##   R CMD INSTALL focus_0.1.10.tar.gz
##
## Run the script from the folder containing it, as the data are read from
## paper_data/, e.g. with
##   Rscript replication.R
## Printed results are written to the console, and figures are saved as PDF
## files in the folder figures/. Timings (system.time) are machine dependent.
## The Monte Carlo simulations of Sections 4.1 and 5.1 run in parallel on 4 and
## 8 workers, respectively; the whole script took about 6 minutes on a 20-core
## Linux workstation.

dir.create("figures", showWarnings = FALSE)


## ==========================================================================
## The Interface
## ==========================================================================

set.seed(42)
Y <- c(rnorm(100, mean = 0), rnorm(50, mean = 1))

library(focus)
det <- detector_create(type = "univariate")
for (i in seq_along(Y)) {
  detector_update(det, Y[i])
  result <- get_statistics(det, family = "gaussian")
  if (result$stat > 20)
    break
}
sprintf("Changepoint detected at time %d", i)

## ---- Changepoint Detectors -----------------------------------------------

## ---- Update the Detector -------------------------------------------------

det_uni <- detector_create(type = "univariate")
det_uni <- det_uni |> detector_update(0.5) 
det_uni <- det_uni |> detector_update(0.1)
print(det_uni)

## ---- Compute Statistics --------------------------------------------------

get_statistics(det_uni, family = "gaussian")$stat

get_statistics(det_uni, family = "gamma", shape = 2)$stat

## ---- Extracting Informations from the Detector Object --------------------

## ---- Accessing Changepoint Candidates ------------------------------------

candidates <- detector_candidates(det)
n_cand <- detector_cands_len(det)
print(n_cand)
print(head(candidates))

## ---- Number of Observations and Cumulative Statistics --------------------

n <- detector_info_n(det)
cumsum_stat <- detector_info_sn(det)
sprintf("%d", n)
sprintf("%.3f", cumsum_stat)

## ---- The Offline Interface -----------------------------------------------

result_offline <- focus_offline(Y, threshold = Inf,
                                type = "univariate",
                                family = "gaussian")

## Figure: fig-r-offline
pdf("figures/fig-r-offline.pdf", width = 7, height = 2)
library(ggplot2)
ggplot(data.frame(time = seq_along(Y), stat = result_offline$stat)) +
  aes(x = time, y = stat) +
  geom_line() +
  labs(x = "Time", y = "Statistic") +
  theme_minimal()
invisible(dev.off())


## ==========================================================================
## Common Use Cases of Focus
## ==========================================================================

## ---- Constrained Detection with One-sided Detectors ----------------------

library(furrr)
plan(multisession, workers = 4)
set.seed(42)
n_sim <- 500
runs <- future_map_dbl(
  1:n_sim, \(i) {
    Y_sim <- rnorm(1e5)
    res <- focus_offline(Y_sim, threshold = Inf,
      type = "univariate_one_sided", side = "right",
      family = "gaussian"
    )
    max(res$stat)
  },
  .options = furrr_options(seed = TRUE)
)
plan(sequential)
print(threshold_99 <- quantile(runs, 0.99))

set.seed(123)
Y <- c(rnorm(1e4, mean = 0), rnorm(1e4, mean = -1), rnorm(1e4, mean = 1))
det_one_sided <- detector_create(type = "univariate_one_sided",
                                 side = "right")
for (i in seq_along(Y)) {
  det_one_sided <- det_one_sided |> detector_update(Y[i])
  result <- det_one_sided |> get_statistics(family = "gaussian")
  if (result$stat > threshold_99)
    break
}
sprintf("Changepoint detected at time %d", i)

## ---- Multivariate Detectors ----------------------------------------------

set.seed(42)
n <- 1000
D <- 3
Y <- matrix(rnorm(n * D), nrow = n, ncol = D)
Y[501:1000, ] <- Y[501:1000, ] + 1

det_mv <- detector_create(type = "multivariate")
stat_trace <- vector("numeric", length = n)
for (i in seq_len(nrow(Y))) {
  det_mv <- det_mv |> detector_update(Y[i, ])
  result <- det_mv |> get_statistics(family = "gaussian")
  stat_trace[i] <- result$stat
  
}

## Figure: fig-multivariate
pdf("figures/fig-multivariate.pdf", width = 7, height = 2)
ggplot(data.frame(t = seq_len(n), trace = stat_trace)) +
  aes(x = t, y = trace) +
  geom_line() +
  geom_vline(xintercept = 500, color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Time", y = "Statistic") +
  theme_minimal()
invisible(dev.off())

dim_idx <- generate_projection_indexes(6, 2)
head(dim_idx, 2)
det_mv <- detector_create(type = "multivariate", dim_indexes = dim_idx)

set.seed(42)
n <- 1000
d <- 6

Y_multi <- rbind(
  matrix(rnorm(5000 * d, mean = -1, 1), ncol = d),
  matrix(rnorm(500 * d, mean = 1.2), ncol = d)
)

system.time(
  res_multi <- focus_offline(Y_multi, threshold = Inf,
                             type = "multivariate", family = "gaussian")
)
system.time(
  res_multi_approx <- focus_offline(Y_multi, threshold = Inf,
                                    type = "multivariate", family = "gaussian",
                                    dim_indexes = dim_idx)
)
all.equal(res_multi$stat, res_multi_approx$stat)

## ---- Anomaly Detection ---------------------------------------------------

set.seed(999)
n <- 1000
t <- seq_len(3*n/2 + 20)
mu_t <- 0.8 * sin(2 * pi * t / 400)
Y_anom <- c(
  rnorm(n/2, mean = mu_t[1:(n/2)]),  rnorm(10, mean = -3),
  rnorm(n/2, mean = mu_t[(n/2 + 11):(n + 10)]), rnorm(10, mean = 5),
  rnorm(n/2, mean = mu_t[(n + 21):(3*n/2 + 20)])
)

run_detector <- function(Y, anomaly_intensity = NULL) {
  det <- detector_create(type = "univariate",
                         anomaly_intensity = anomaly_intensity)
  stat_trace <- numeric(length(Y))
  for (i in seq_along(Y)) {
    det <- det |> detector_update(Y[i])
    res <- det |> get_statistics(family = "gaussian", theta0 = 0)
    stat_trace[i] <- res$stat
  }
  stat_trace
}
stat_no_thresh <- run_detector(Y_anom, anomaly_intensity = NULL)
stat_thresh <- run_detector(Y_anom, anomaly_intensity = 1.5)

## Figure: fig-anomaly_detection
pdf("figures/fig-anomaly_detection.pdf", width = 7, height = 4)
df <- data.frame(time = rep(seq_along(Y_anom), 3),
                 value = c(Y_anom, stat_no_thresh, stat_thresh),
                 type = rep(c("Data", "No threshold", "With threshold"),
                            each = length(Y_anom)))
ggplot(df) +
  aes(x = time, y = value) +
  geom_line() +
  facet_wrap(~type, ncol = 1, scales = "free_y") +
  labs(x = "Time", y = "") +
  theme_minimal()
invisible(dev.off())

## ---- Non-Parametric Changepoint Detection --------------------------------

set.seed(42)
Y1_pre <- rt(800, df = 2)
Y1_post <- rt(200, df = 2) + 0.5  # Location shift
Y1 <- c(Y1_pre, Y1_post)
Y2_pre <- rnorm(800)
Y2_post <- rnorm(200)
Y2_post[abs(Y2_post) > 2] <- Y2_post[abs(Y2_post) > 2] + 20  # Tail change
Y2 <- c(Y2_pre, Y2_post)

q1 <- qt(seq(0.01, 0.99, length.out = 5), df = 2)
det_np1 <- detector_create(type = "npfocus", quantiles = q1)
stat_trace1 <- matrix(nrow = length(Y1), ncol = 2)
for (i in seq_along(Y1)) {
  det_np1 <- det_np1 |> detector_update(Y1[i])
  res <- det_np1 |> get_statistics(family = "npfocus")
  stat_trace1[i, ] <- res$stat
}
q2 <- qnorm(seq(0.01, .99, length.out = 4))
det_np2 <- detector_create(type = "npfocus", quantiles = q2)
stat_trace2 <- matrix(nrow = length(Y2), ncol = 2)
for (i in seq_along(Y2)) {
  det_np2 <- det_np2 |> detector_update(Y2[i])
  res <- det_np2 |> get_statistics(family = "npfocus")
  stat_trace2[i, ] <- res$stat
}

## Figure: fig-npfocus
pdf("figures/fig-npfocus.pdf", width = 9, height = 5)
# Combine data for plotting
time1 <- seq_along(Y1)
time2 <- seq_along(Y2)
change_point1 <- 800
change_point2 <- 800

df_plot <- rbind(
  data.frame(scenario = "Location shift (t)", time = time1, value = Y1, statistic = "Data"),
  data.frame(scenario = "Location shift (t)", time = time1, value = stat_trace1[, 1], statistic = "Sum"),
  data.frame(scenario = "Location shift (t)", time = time1, value = stat_trace1[, 2], statistic = "Max"),
  data.frame(scenario = "Tail change (Normal)", time = time2, value = Y2, statistic = "Data"),
  data.frame(scenario = "Tail change (Normal)", time = time2, value = stat_trace2[, 1], statistic = "Sum"),
  data.frame(scenario = "Tail change (Normal)", time = time2, value = stat_trace2[, 2], statistic = "Max")
)

# Reorder factors for better visualization
df_plot$statistic <- factor(df_plot$statistic, levels = c("Data", "Sum", "Max"))

ggplot(df_plot) +
  aes(x = time, y = value) +
  geom_line(size = 0.5) +
  geom_vline(aes(xintercept = xint), data = data.frame(scenario = c("Location shift (t)", "Tail change (Normal)"), 
                                                         xint = c(change_point1, change_point2)), linetype = "dashed", linewidth = 1) +
  facet_grid(statistic ~ scenario, scales = "free_y") +
  labs(x = "Time", y = "") +
  theme_minimal()
invisible(dev.off())

## ---- Autoregressive Changepoint Detection --------------------------------

set.seed(123)
ar_coefs <- c(0.7, -0.3)
Y_pre <- arima.sim(n = 5000, model = list(ar = ar_coefs), sd = 1)
Y_post <- 2 + arima.sim(n = 100, model = list(ar = ar_coefs), sd = 1)
Y <- c(Y_pre, Y_post)

ar_coefs <- ar.yw(Y_pre[1:1000], order.max = 5)$ar
det_arp <- detector_create(type = "arp", rho = ar_coefs, mu0_arp = 0)
stat_trace <- numeric(length(Y))

for (i in 1001:length(Y)) {
  det_arp <- det_arp |> detector_update(Y[i])
  result <- det_arp |> get_statistics(family = "arp")
  stat_trace[i] <- result$stat
}

## Figure: fig-arp_detection
pdf("figures/fig-arp_detection.pdf", width = 7, height = 3)
# Create data frame for visualization
df <- data.frame(time = rep(seq_along(Y), 2),
                 value = c(Y, stat_trace),
                 type = rep(c("Data", "Statistic"), each = length(Y)))

ggplot(df) +
  aes(x = time, y = value) +
  geom_line() +
  geom_vline(xintercept = 5000, linetype = "dashed", linewidth = 1) +
  facet_wrap(~type, ncol = 1, scales = "free_y") +
  labs(x = "Time", y = "") +
  theme_minimal()
invisible(dev.off())


## ==========================================================================
## Some Examples of Real-Time Applications
## ==========================================================================

## ---- Implementing a Custom Cost Function: Change in Mean and Variance on NBA Plusminus Scores ---

get_seg_meanvar <- function(info_seg, min_var=1){
  if(info_seg[1] > 0){
    est_mse <-  info_seg[3] - info_seg[2]^2/info_seg[1]
    est_var <- est_mse/info_seg[1]
    if(est_var < min_var) est_var <- min_var
    return(info_seg[1]*log(est_var) + est_mse/est_var)
    } else { return(Inf) }
}

library(purrr)
get_stat_meanvar <- function(det_ptr, min_var=1) {
  mat <- detector_candidates(det_ptr)
  n <- detector_info_n(det_ptr)
  all <- detector_info_sn(det_ptr)
  stats <- map_dbl(seq_along(mat$tau)[-1], \(i) {
    info_seg_left <- c(mat$tau[i], mat$st[[i]])
    
    info_seg_right <- c(n, all) - info_seg_left
    return(
        - get_seg_meanvar(info_seg_left, min_var) -
        get_seg_meanvar(info_seg_right, min_var)  +
        get_seg_meanvar(info_seg_left + info_seg_right, min_var)
    )
  })
  max(stats)
}

library(dplyr)

dat <- readRDS("paper_data/cle_data.rds")

library(furrr)
plan(multisession, workers = 8)

set.seed(123)
dat_past <- dat %>% filter(yearSeason <= 2010, typeSeason == "Regular Season")

y <- future_map(1:10^3, \(i) {
  plusminus_team <- sample(dat_past$plusminusTeam, 1000, replace = TRUE)
  dat_for_focus <- rbind(plusminus_team, plusminus_team^2)
  
  det <- detector_create(type = "multivariate")
  stat_record <- vector("numeric", length = ncol(dat_for_focus))
  for (i in 1:ncol(dat_for_focus)) {
    detector_update(det, dat_for_focus[, i])
    stat_record[i] <- get_stat_meanvar(det)
  }
  max(stat_record)
}, .options = furrr_options(seed = TRUE))
results <- do.call(c, y)

threshold_99 <- quantile(results, 0.99)

print(threshold_99)

dat_recent <- dat %>% filter(yearSeason > 2010, typeSeason == "Regular Season")
Y_test <- rbind(dat_recent$plusminusTeam, dat_recent$plusminusTeam^2)
det <- detector_create(type = "multivariate")
stat_record <- vector("numeric", length = ncol(Y_test))
for (i in 1:ncol(Y_test)) {
  detector_update(det, Y_test[, i])
  stat_record[i] <- get_stat_meanvar(det, min_var=1)
}

## Figure: fig-nba-application
pdf("figures/fig-nba-application.pdf", width = 7, height = 3)
df <- data.frame(time = rep(seq_along(dat_recent$plusminusTeam), 2),
                 value = c(dat_recent$plusminusTeam, stat_record),
                 type = rep(c("Data", "Statistic"), each = length(dat_recent$plusminusTeam)))

ggplot(df) +
  aes(x = time, y = value) +
  geom_line() +
  geom_hline(yintercept = threshold_99, color = "gray", linetype = "dashed", linewidth = 1) +
  facet_wrap(~type, ncol = 1, scales = "free_y") +
  labs(x = "Time", y = "") +
  theme_minimal()
invisible(dev.off())
