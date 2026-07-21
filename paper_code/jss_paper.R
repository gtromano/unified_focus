#| include: false

# needed for the "R> " prompt in the code chunks
options(prompt = "R> ", continue = "+  ")


#| prompt: true

set.seed(42)
Y <- c(rnorm(100, mean = 0), rnorm(50, mean = 1))


#| prompt: true

library(focus)

det <- detector_create(type = "univariate")

for (i in seq_along(Y)) {
  detector_update(det, Y[i])
  result <- get_statistics(det, family = "gaussian")
  
  if (result$stat > 20) {
    print(
      paste(
        "Changepoint detected at time", i,
        "with estimated changepoint at", result$changepoint
      )
    )
    break
    
  }
  
}


#| prompt: true

## result <- det |>
##   detector_update(Y[i]) |>
##   get_statistics(family = "gaussian")


#| eval: false
## detector_create(type, dim_indexes = NULL, quantiles = NULL,
##   pruning_mult = 2L, pruning_offset = 1L, side = "right",
##   anomaly_intensity = NULL, rho = NULL, mu0_arp = NULL
## )


#| eval: false
## detector_update(info_ptr, y)


#| prompt: true

det_uni <- detector_create(type = "univariate")
det_uni <- detector_update(det_uni, 0.5) 
det_uni <- detector_update(det_uni, 0.1)
print(det_uni)


#| eval: false
## get_statistics(info_ptr, family, theta0 = NULL, shape = NULL)


#| prompt: true

get_statistics(det_uni, family = "gaussian")

get_statistics(det_uni, family = "gamma", shape = 2)


#| eval: false
## detector_info_n(info_ptr)
## detector_info_sn(info_ptr)


#| prompt: true

n <- detector_info_n(det)
cumsum_stat <- detector_info_sn(det)

sprintf("Processed %d observations", n)
sprintf("Current cumulative sum: %.3f", cumsum_stat)


#| eval: false
## detector_candidates(info_ptr)
## detector_cands_len(info_ptr)


#| prompt: true

# Get candidate information
candidates <- detector_candidates(det)
n_cand <- detector_cands_len(det)

sprintf("Number of candidate segments: %d", n_cand)
print(head(candidates))


#| eval: false

## focus_offline(Y, threshold, type = "univariate", family = "gaussian",
##   theta0 = NULL, dim_indexes = NULL, quantiles = NULL,
##   pruning_mult = 2L, pruning_offset = 1L, side = "right",
##   shape = NULL, anomaly_intensity = NULL, rho = NULL, mu0_arp = NULL
## )


#| label: fig-r-offline
#| fig-cap: "The trace of the statistics for one simulated time series using the R `focus_offline` interface."
#| fig-width: 7
#| fig-height: 3
#| prompt: true

result_offline <- focus_offline(Y,
                                threshold = Inf,
                                type = "univariate",
                                family = "gaussian")

library(ggplot2)
ggplot(data.frame(time = seq_along(Y), stat = result_offline$stat)) +
  aes(x = time, y = stat) +
  geom_line() +
  labs(x = "Time", y = "Statistic") +
  theme_minimal()


#| label: threshold_one_sided
#| message: false
#| warning: false
#| cache: true

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


#| prompt: true

set.seed(123)
Y <- c(rnorm(1e4, mean = 0), rnorm(1e4, mean = -1), rnorm(1e4, mean = 1))
det_one_sided <- detector_create(type = "univariate_one_sided",
                                 side = "right")

for (i in seq_along(Y)) {
  
  det_one_sided <- det_one_sided |> detector_update(Y[i])
  result <- det_one_sided |> get_statistics(family = "gaussian")
    
  if (result$stat > threshold_99) {
    print(
      paste(
        "Changepoint detected at time", i,
        "with estimated changepoint at", result$changepoint
      )
    )
    break
  }
}


#| prompt: true
#| label: fig-multivariate
#| fig-cap: "The trace of the statistics for one simulated 3-dimentional multivariate time series using the R online interface. Dashed red line illustrates the true change location."
#| fig-width: 7
#| fig-height: 2

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

ggplot(data.frame(t = seq_len(n), trace = stat_trace)) +
  aes(x = t, y = trace) +
  geom_line() +
  geom_vline(xintercept = 500, color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Time", y = "Statistic") +
  theme_minimal()



dim_idx <- generate_projection_indexes(D = 10, p = 2)
print(dim_idx[1:3])
det_mv <- detector_create(type = "multivariate", dim_indexes = dim_idx)


#| prompt: true
#| cache: true

set.seed(42)
n <- 1000
p <- 6

Y_multi <- rbind(
  matrix(rnorm(5000 * p, mean = -1, 1), ncol = p),
  matrix(rnorm(500 * p, mean = 1.2), ncol = p)
)

system.time(
  res_multi <- focus_offline(Y_multi, threshold = Inf,
                             type = "multivariate", family = "gaussian")
)

dim_indexes <- generate_projection_indexes(6, 2)
system.time(
  res_multi_approx <- focus_offline(Y_multi, threshold = Inf,
                                    type = "multivariate", family = "gaussian",
                                    dim_indexes = dim_indexes)
)

all.equal(res_multi$stat, res_multi_approx$stat)


#| prompt: true

set.seed(999)
n <- 1000
Y_anom <- c(
  rnorm(n/2, mean = 0), rnorm(10, mean = -3),
  rnorm(n/2, mean = 0), rnorm(10, mean = 5),
  rnorm(n/2, mean = 0)
)


#| prompt: true
#| label: fig-anomaly_detection
#| fig-cap: "The trace of the statistics for one simulated time series with transient anomalies using the R online interface."
#| fig-width: 7
#| fig-height: 3

det <- detector_create(type = "univariate", anomaly_intensity = 2)

stat_trace <- numeric(length(Y_anom))

for (i in seq_along(Y_anom)) {
  det <- det |> detector_update(Y_anom[i])
  res <- det |> get_statistics(family = "gaussian", theta0 = 0)
  stat <- res$stat
  stat_trace[i] <- stat
}

df <- data.frame(time = rep(seq_along(Y_anom), 2),
                 value = c(Y_anom, stat_trace), 
                 type = rep(c("Data", "Statistic"), each = length(Y_anom)))
ggplot(df) +
  aes(x = time, y = value) +
  geom_line() +
  facet_wrap(~type, ncol = 1, scales = "free_y") +
  labs(x = "Time", y = "") +
  theme_minimal()



#| prompt: true

quants <- qnorm(seq(0.01, .99, length.out = 5))


#| prompt: true
#| label: fig-npfocus
#| fig-cap: "The data and the resulting traces of the max and sum statistics for one simulated time series with a change in distribution using the R online interface. Dashed red line illustrates the true change location."
#| fig-width: 7
#| fig-height: 5

set.seed(42)
Y <- c(rnorm(1000), rcauchy(200))

det_np <- detector_create(type = "npfocus", quantiles = quants)

stat_trace <- matrix(nrow = length(Y), ncol = 2)

for (i in seq_along(Y)) {
  det_np <- det_np |> detector_update(Y[i])
  res <- det_np |> get_statistics(family = "npfocus")
  stat_trace[i, ] <- res$stat
}

# ggplot of the data and the two statistics, as in the snippet above 

df <- data.frame(time = rep(seq_along(Y), 3),
                 value = c(Y, stat_trace[,1], stat_trace[,2]), 
                 type = rep(c("Data", "Statistic (Max)", "Statistic (Sum)"),
                            each = length(Y)))
ggplot(df) +
  aes(x = time, y = value) +
  geom_line() +
  geom_vline(xintercept = 1000, color = "red", linetype = "dashed", linewidth = 1) +
  facet_wrap(~type, ncol = 1, scales = "free_y") +
  labs(x = "Time", y = "") +
  theme_minimal()


#| prompt: true

set.seed(123)
ar_coefs <- c(0.7, -0.3)  # AR(2) coefficients
n_pre <- 500
n_post <- 500

# Generate AR(2) process with change in mean
Y_pre <- arima.sim(n = n_pre, model = list(ar = ar_coefs), sd = 1)
Y_post <- 2 + arima.sim(n = n_post, model = list(ar = ar_coefs), sd = 1)
Y <- c(Y_pre, Y_post)


#| prompt: true
#| label: fig-arp_detection
#| fig-cap: "The AR(2) process with mean shift and the resulting detection statistic using the R online interface. Dashed red line illustrates the true change location."
#| fig-width: 7
#| fig-height: 3

det_arp <- detector_create(type = "arp", rho = ar_coefs)
stat_trace <- numeric(length(Y))

for (i in seq_along(Y)) {
  det_arp <- det_arp |> detector_update(Y[i])
  result <- det_arp |> get_statistics(family = "arp")
  stat_trace[i] <- result$stat
}

# Create data frame for visualization
df <- data.frame(time = rep(seq_along(Y), 2),
                 value = c(Y, stat_trace),
                 type = rep(c("Data", "Statistic"), each = length(Y)))

ggplot(df) +
  aes(x = time, y = value) +
  geom_line() +
  geom_vline(xintercept = n_pre, color = "red", linetype = "dashed", linewidth = 1) +
  facet_wrap(~type, ncol = 1, scales = "free_y") +
  labs(x = "Time", y = "") +
  theme_minimal()


#| include: false

# needed for the python ">>> " prompt in the code chunks
options(prompt = ">>> ", continue = "... ")


## from astropy.table import Table

## from astropy.io import fits

## import pandas as pd

## import numpy as np

## import matplotlib.pyplot as plt

## 

## import focus_cpt


## data_dict = {}

## step = 0.1

## 

## for hour in range(1, 2):

##     fp = 'glg_tte_n2_bn250814432_v00.fit'

## 

##     energy_bands = Table.read(fp, hdu=1)

##     table = Table.read(fp, hdu=2)

##     time = Table.read(fp, hdu=3)

##     #this is about ten minutes of data

##     df = table.to_pandas()

##     #df

##     grb_counts = df['TIME'] - df.loc[0, 'TIME']

## 

##     bins = np.arange(0, 360, step)

##     bin_map = pd.cut(grb_counts, bins=bins).apply(lambda I: I.left if pd.notna(I) else np.nan)

## 

## 

##     bin_series = pd.DataFrame(bin_map).dropna()

##     bin_series['count'] = 1

##     d = bin_series.groupby("TIME")["count"].count()

##     data_dict[hour]=d

## 

## data_full = pd.concat(list(data_dict.values())).reset_index(drop=True)

## count_times = data_full.index * pd.to_timedelta(step, "s") + pd.to_datetime("2025-08-14 10:00:00.000")

## data_full.index = count_times

## data_full[data_full==0]=np.nan


## plt.plot(data_full['2025-08-14 10:01':'2025-08-14 10:02'], linewidth=0.3, label="radiation count")

## plt.title("Gamma ray burst")

## plt.fill_between(data_full['2025-08-14 10:01':'2025-08-14 10:02'].index[768: 818], y1=70, y2=180, color="red", alpha=0.2, label="most significant interval")

## plt.tick_params("x", rotation=45)

## plt.legend()

## plt.show()


## data_snippet = data_full['2025-08-14 10:01':'2025-08-14 10:02']

## 

## fig, ax = plt.subplots()

## res = focus_cpt.focus_offline(data_snippet, theta0=data_full['2025-08-14 10:01':'2025-08-14 10:02'].mean(), threshold=np.inf, type="univariate_one_sided", family="poisson")

## plt.plot(data_snippet.index, np.sqrt(2*res["stat"]), color="C0", label="significance")

## plt.title("Using focus-cpt to detect a gamma-ray burst")

## 

## plt.axhline(5, color="red", linestyle="--", lw=2, label="5-sigma threshold")

## #plt.vlines(768, ymin=0, ymax=20, color="red", alpha=0.3)

## #plt.vlines(818, ymin=0, ymax=20, color="red", alpha=0.3)

## plt.fill_between(data_full['2025-08-14 10:01':'2025-08-14 10:02'].index[768: 818], y1=0, y2=20, color="red", alpha=0.2, label="most significant interval")

## plt.tick_params("x", rotation=45)

## plt.xlabel("time")

## plt.ylabel("significance")

## plt.legend()

## plt.show()


## data_window = data_full['2025-08-14 10:01':'2025-08-14 10:02']

## detector = focus_cpt.Detector(type="univariate")

## 

## theta0 = data_full['2025-08-14 10:01':'2025-08-14 10:02'].mean()

## 

## stats = []

## 

## for i, value in enumerate(data_window):

##     detector.update(value)

##     # Use running mean as theta0 at each iteration

##     # theta0 = data_window.iloc[:i+1].mean()

##     result = detector.get_statistics(family="poisson", theta0=theta0)

##     stats.append(result['stat'])

## 

## stat_array = np.array(stats)

## plt.figure()

## plt.plot(np.sqrt(2*stat_array), color="C0", label="significance")

## plt.title("Using focus-cpt to detect a gamma-ray burst")

## plt.hlines(5, xmin=0, xmax=1200, color="C1", label="5-sigma threshold")

## #plt.vlines(768, ymin=0, ymax=20, color="red", alpha=0.3)

## #plt.vlines(818, ymin=0, ymax=20, color="red", alpha=0.3)

## plt.fill_between(np.arange(768, 818), y1=0, y2=20, color="red", alpha=0.2, label="most significant interval")

## plt.legend()

## plt.show()


#| include: false

# needed for the python ">>> " prompt in the code chunks
options(prompt = ">>> ", continue = "... ")


## import numpy as np

## from focus_cpt import Detector

## 

## 

## np.random.seed(42)

## Y = np.concatenate([np.random.randn(100), np.random.randn(50) + 1])

## 

## det = Detector(type="univariate")

## 

## for i, y in enumerate(Y):

##     det.update(y)

##     result = det.get_statistics(family="gaussian")

## 

##     if result['stat'] > 20:

##         print(

##             f"Changepoint detected at time {i}: "

##             f"estimated changepoint at {result['changepoint']}"

##         )

##         break


## import pandas as pd

## from plotnine import *

## from focus_cpt import focus_offline

## 

## result_offline = focus_offline(Y,

##                                 threshold=np.inf,

##                                 type="univariate",

##                                 family="gaussian")

## 

## stat = result_offline['stat'].flatten()

## 

## df = pd.DataFrame({"time": range(1, len(Y) + 1), "stat": stat})

## (

##   ggplot(df) +

##     aes(x="time", y="stat") +

##     geom_line() +

##     labs(x="Time", y="Statistic") +

##     theme_minimal()

## )


#| include: false

# needed for the "R> " prompt in the code chunks
options(prompt = "R> ", continue = "+  ")

