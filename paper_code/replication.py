"""Replication script (Python) for
"focus and focus-cpt: Fast Online Changepoint Detection in R and Python"
(Journal of Statistical Software).

This file is generated from the code chunks of jss_paper.qmd by
make_replication_scripts.py (make replication-scripts): do not edit by hand.

This script reproduces all Python results and figures of the manuscript, in
their order of appearance: the Python example of Section 3 (the interface),
Section 5.2 (gamma-ray burst detection), Section 5.3 (constrained up-down model
for spike inference) and Appendix A (the Python interface). The R results are
reproduced by the companion script replication.R.

Requirements: Python (>= 3.8) and the packages focus-cpt (>= 0.1.10), numpy,
pandas, scipy, plotnine and astropy, all available from PyPI. The focus-cpt
package can also be installed from the source package submitted with the
manuscript (this requires a C++ compiler, CMake and the Qhull library, see its
README):
    pip install focus_cpt-0.1.10.tar.gz

Run the script from the folder containing it, as the data are read from
paper_data/, e.g. with
    python replication.py
Printed results are written to the console, and figures are saved as PDF files
in the folder figures/. The cross-validation of Section 5.3 takes about a
minute.
"""

import os

os.makedirs("figures", exist_ok=True)


# ==========================================================================
# The Interface
# ==========================================================================

import numpy as np
from focus_cpt import Detector

np.random.seed(42)
Y = np.concatenate([np.random.randn(100), np.random.randn(50) + 1])
det = Detector(type="univariate")
i = 0

det.update(Y[i])
result = det.get_statistics(family="gaussian")


# ==========================================================================
# Some Examples of Real-Time Applications
# ==========================================================================

# ---- Gamma-Ray Burst Detection -------------------------------------------

from astropy.table import Table
import pandas as pd
import numpy as np
from scipy import stats
from plotnine import (ggplot, aes, geom_line, geom_hline, geom_ribbon,
                      theme_minimal, xlim, labs)
import focus_cpt
import warnings
from astropy.units import UnitsWarning
warnings.filterwarnings('ignore', category=UnitsWarning) # needed to ignore warnings about time units in the Fermi-GBM

data_dict = {}
step = 0.01  # 10 milliseconds

for hour in range(1, 2):
    fp = 'paper_data/glg_tte_n2_bn250814432_v00.fit'
    table = Table.read(fp, hdu=2)
    df = table.to_pandas()
    grb_counts = df['TIME'] - df.loc[0, 'TIME']

    bins = np.arange(0, 360, step)
    bin_map = pd.cut(grb_counts, bins=bins).apply(lambda I: I.left if pd.notna(I) else np.nan)
    
    bin_series = pd.DataFrame(bin_map).dropna()
    bin_series['count'] = 1
    d = bin_series.groupby("TIME")["count"].count()
    data_dict[hour] = d

data_full = pd.concat(list(data_dict.values())).reset_index(drop=True)
count_times = data_full.index * pd.to_timedelta(step, "s") + pd.to_datetime("2025-08-14 10:00:00.000")
data_full.index = count_times
data_full[data_full == 0] = np.nan

data_window = data_full['2025-08-14 10:01':'2025-08-14 10:02']
theta0_est = np.mean(data_full['2025-08-14 10:00':'2025-08-14 10:00'])
detector = focus_cpt.Detector(type="univariate_one_sided",
                              side="right",
                              anomaly_intensity = 1.5)
stats_list = []
for i, value in enumerate(data_window):
    detector.update(value)
    result = detector.get_statistics(family="poisson", theta0=theta0_est)
    stats_list.append(result['stat'])
stat_array = np.array(stats_list)
significance = np.sqrt(2 * stat_array)

# Figure: fig-gamma-ray-burst
from plotnine import facet_wrap, geom_rect, geom_vline, theme

# Create dataframe for plotting
df_grb = pd.DataFrame({
    'index': range(len(data_window)),
    'count': data_window.values,
    'significance': significance
})

# Prepare data for faceting
df_plot = pd.DataFrame({
    'index': list(df_grb['index']) + list(df_grb['index']),
    'value': list(df_grb['count']) + list(df_grb['significance']),
    'type': ['Count'] * len(df_grb) + ['Significance'] * len(df_grb)
})

# Define the significant burst interval
grb_start = 7680
grb_end = 8180
threshold = 5

# Find the first point where significance exceeds the threshold
detection_point = np.where(significance > threshold)[0]
detection_point = detection_point[0] if len(detection_point) > 0 else None

# Create the faceted plot
p = (ggplot(df_plot, aes(x='index', y='value')) +
    geom_rect(aes(xmin=grb_start, xmax=grb_end, ymin=-np.inf, ymax=np.inf), 
              alpha=0.2, fill="lightgray", inherit_aes=False) +
    geom_line(size=0.7) +
    (geom_vline(xintercept=detection_point, linetype="dashed", size=1) if detection_point is not None else None) +
    facet_wrap('~type', ncol=1, scales='free_y') +
    xlim(0, 12000) +
    labs(x="Time (10 ms bins)", y="") +
    theme_minimal() +
    theme(figure_size=(7, 3))
)

p
p.save("figures/fig-gamma-ray-burst.pdf", verbose=False)

# ---- Constrained Up-Down Model for Spike Inference -----------------------

from plotnine import (ggplot, aes, geom_line, geom_vline, theme_minimal,
                      xlim, geom_segment, geom_point, geom_ribbon, labs)
import plotnine as pl
from focus_cpt import Detector, focus_offline
import pickle

with open("paper_data/example_trace.pkl", "rb") as f:
    neur_trace = pickle.load(f)


df_spikes = pd.DataFrame({
    't'   : neur_trace["spikes"],
    'y'   : [-0.5] * len(neur_trace["spikes"]),
    'yend': [-1.0] * len(neur_trace["spikes"])
})

df = pd.DataFrame({'t': neur_trace["time"], 'Y': neur_trace["trace"]})

def run_detector(Y, time, r_threshold, l_threshold):
    r_det = Detector(type='univariate_one_sided', side='right')
    l_det = Detector(type='univariate_one_sided', side='left')
    up_detections = []
    stp_times     = []
    stp_types     = []
    for t, y in enumerate(Y):
        r_det.update(y)
        l_det.update(y)
        r_s = r_det.get_statistics(family='gaussian')['stat']
        l_s = l_det.get_statistics(family='gaussian')['stat']
        cpt_found = False
        if r_s > r_threshold:
            cpt_found, cpt_type = True, 'right'
        elif l_s > l_threshold:
            cpt_found, cpt_type = True, 'left'
        if cpt_found:
            if cpt_type == 'right':
                up_detections.append(t)
            stp_times.append(time[t])
            stp_types.append(cpt_type)
            r_det = Detector(type='univariate_one_sided', side='right')
            l_det = Detector(type='univariate_one_sided', side='left')
    return {
        'up_detections': np.array(up_detections),
        'stp_times'    : stp_times,
        'stp_types'    : stp_types,
    }

def van_rossum_distance(detected, true_spikes, tau_samples):
    """Van Rossum distance between two spike trains (lower = better).

    Each train is convolved with a causal exponential kernel
    h(t) = exp(-t / tau) for t >= 0. The squared L2 distance between
    the resulting continuous functions is computed analytically:

        D = <det, det> - 2 <det, true> + <true, true>

    where the inner product of two trains A, B is

        <A, B> = (tau/2) * sum_{i in A, j in B} exp(-|t_i - t_j| / tau)

    This handles both failure modes better than toleranced F1:
      - A detection far from any true spike (e.g. during the calcium
        decay) is penalised more than one that is merely slightly off.
      - A detection during the decay is actually penalised MORE than
        a complete miss, because it adds a spurious kernel contribution
        far from any true spike, whereas a miss only incurs the cost of
        the unmatched true spike's self inner product.
      - Close detections near dense bursts are graded rather than
        all-or-nothing, so missing some burst spikes is not catastrophic.

    Parameters
    ----------
    detected    : array-like of sample indices of detected spikes
    true_spikes : array-like of sample indices of true spikes
    tau_samples : exponential time constant in samples
                  (e.g. 5 samples = 50 ms at 100 Hz)
    """
    def inner_product(a, b):
        if len(a) == 0 or len(b) == 0:
            return 0.0
        diff = np.abs(np.asarray(a, float)[:, None] - np.asarray(b, float)[None, :])
        return (tau_samples / 2.0) * np.sum(np.exp(-diff / tau_samples))

    return (inner_product(detected,    detected)
            - 2 * inner_product(detected, true_spikes)
            + inner_product(true_spikes, true_spikes))

time  = np.array(neur_trace["time"])
trace = np.array(neur_trace["trace"])

spike_indices = np.array(
    [np.argmin(np.abs(time - s)) for s in neur_trace["spikes"]]
)

# ── Train / test split (60 / 40 by time) ─────────────────────────────────────
SPLIT_FRAC = 0.3
split_idx  = int(len(trace) * SPLIT_FRAC)

Y_train    = trace[:split_idx]
Y_test     = trace[split_idx:]
time_train = time[:split_idx]
time_test  = time[split_idx:]

train_spikes = spike_indices[spike_indices <  split_idx]
test_spikes  = spike_indices[spike_indices >= split_idx] - split_idx

# ── 5-fold cross-validation over a 2D threshold grid ─────────────────────────
# Both r_threshold and l_threshold are tuned jointly.
# Each contiguous fold is treated as an independent segment so the streaming
# detector starts fresh — no warm-up leakage between folds.
# The scoring criterion is the van Rossum distance (minimise).

TAU_SAMPLES = 5     # exponential kernel time constant: 50 ms at 100 Hz
N_FOLDS     = 5

# Candidate grids
r_thresholds = np.concatenate([
    np.arange(0.005, 0.05, 0.005),
    np.arange(0.05,  0.30, 0.02),
])
l_thresholds = np.array([0.2, 0.5, 1.0, 2.0, 5.0])

n_train   = len(Y_train)
fold_size = n_train // N_FOLDS

cv_results = {}   # (r_thr, l_thr) -> mean CV van Rossum distance

for r_thr in r_thresholds:
    for l_thr in l_thresholds:
        fold_vr = []
        for fold in range(N_FOLDS):
            v_start    = fold * fold_size
            v_end      = v_start + fold_size if fold < N_FOLDS - 1 else n_train
            val_spikes = train_spikes[
                (train_spikes >= v_start) & (train_spikes < v_end)
            ] - v_start

            res = run_detector(Y_train[v_start:v_end],
                               time_train[v_start:v_end], r_thr, l_thr)
            fold_vr.append(
                van_rossum_distance(res['up_detections'], val_spikes, TAU_SAMPLES)
            )
        cv_results[(r_thr, l_thr)] = np.mean(fold_vr)

# ── Select best pair (minimum van Rossum distance) ────────────────────────────
(best_r_threshold, best_l_threshold) = min(cv_results, key=cv_results.get)
best_cv_vr = cv_results[(best_r_threshold, best_l_threshold)]

res_test = run_detector(Y_test, time_test,
                        best_r_threshold, best_l_threshold)

# Figure: fig-calcium_trace
test_spike_times = time_test[test_spikes]
df_test_changes = pd.DataFrame({
    't'   : res_test['stp_times'],
    'y'   : [-1.0] * len(res_test['stp_times']),
    'yend': [-1.5] * len(res_test['stp_times']),
    'type': res_test['stp_types'],
})
df_test = pd.DataFrame({'t': time_test, 'Y': Y_test})

p = (ggplot(df, aes(x='t', y='Y')) +
    geom_line() +
    geom_segment(df_spikes,
                 aes(x='t', xend='t', y='y', yend='yend'),
                 alpha=0.6, colour="orange") +
    geom_segment(df_test_changes.query('type == "right"'),
                 aes(x='t', xend='t', y='y', yend='yend'),
                 colour="blue") +
    geom_vline(xintercept=time_test[0], color="red", linetype="dashed") +             
    labs(title=f"Thresholds: right={best_r_threshold:.3f}, left={best_l_threshold:.3f}",
         x="Time (s)", y="Y") +
    theme_minimal() +
    theme(figure_size=(7, 3))
)
p.save("figures/fig-calcium_trace.pdf", verbose=False)


# ==========================================================================
# The Python Interface: Differences and Similarities
# ==========================================================================

# ---- Example Usage in Python ---------------------------------------------

import numpy as np
from focus_cpt import Detector


np.random.seed(42)
Y = np.concatenate([np.random.randn(100), np.random.randn(50) + 1])

det = Detector(type="univariate")

for i, y in enumerate(Y):
    det.update(y)
    result = det.get_statistics(family="gaussian")
    
    if result['stat'] > 20:
        print(
            f"Changepoint detected at time {i + 1}: "
            f"estimated changepoint at {result['changepoint']}"
        )
        break

# Figure: fig-python-offline
import pandas as pd
from plotnine import *
from focus_cpt import focus_offline

result_offline = focus_offline(Y,
                                threshold=np.inf,
                                type="univariate",
                                family="gaussian")

stat = result_offline['stat'].flatten()

df = pd.DataFrame({"time": range(1, len(Y) + 1), "stat": stat})
p = (
  ggplot(df) +
    aes(x="time", y="stat") +
    geom_line() +
    labs(x="Time", y="Statistic") +
    theme_minimal() +
    theme(figure_size=(7, 3))
)
p.save("figures/fig-python-offline.pdf", verbose=False)
