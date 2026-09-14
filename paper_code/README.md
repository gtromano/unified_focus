# Replication Materials for "focus and focus-cpt: Fast Online Changepoint Detection in R and Python"

This folder contains the replication materials for the article submitted to the
*Journal of Statistical Software*, together with the Quarto source of the
manuscript.

## Contents

| File | Description |
|------|-------------|
| `replication.R` | Plain R script reproducing all R results and figures (Sections 3, 4 and 5.1). |
| `replication.py` | Plain Python script reproducing all Python results and figures (Section 3, Sections 5.2 and 5.3, Appendix A). |
| `paper_data/` | Data sets used in Section 5 (see [Data](#data)). |
| `make_replication_scripts.py` | Generates the two replication scripts from the code chunks of `jss_paper.qmd` (not needed for replication). |
| `jss_paper.qmd`, `bibliography.bib`, `Makefile` | Quarto source of the manuscript and build tools (not needed for replication). |

The software is provided as source packages, submitted alongside the manuscript:

- `focus_0.1.10.tar.gz`: the R package **focus** (also on CRAN);
- `focus_cpt-0.1.10.tar.gz`: the Python package **focus-cpt** (also on PyPI).

## Replicating the Results

The scripts are generated from the manuscript with `make replication-scripts`
(`python make_replication_scripts.py`): every evaluated code chunk is copied
verbatim, in order of appearance, and the only added lines are those saving the
figures to files.

Both scripts must be run from this folder, as the data are read from
`paper_data/`. Results are printed to the console and figures are saved as PDF
files in `figures/`, named after the figure labels of the manuscript. All random
number generation is seeded. Timings reported by `system.time()` in Section 4.2
are machine dependent. On a 20-core Linux workstation, `replication.R` runs in
about 6 minutes and `replication.py` in about 15 seconds.

### R

Requirements: R (>= 4.1.0) and the CRAN packages ggplot2, furrr, purrr and
dplyr. The Monte Carlo simulations of Sections 4.1 and 5.1 run in parallel,
using 4 and 8 background R sessions respectively (via furrr).

```bash
R CMD INSTALL focus_0.1.10.tar.gz      # or: Rscript -e "install.packages('focus')"
Rscript -e "install.packages(c('ggplot2', 'furrr', 'purrr', 'dplyr'))"
Rscript replication.R
```

### Python

Requirements: Python (>= 3.8) and the PyPI packages numpy, pandas, scipy, plotnine
and astropy. Installing focus-cpt from the source package requires a C++17
compiler, CMake (>= 3.15) and the Qhull library (e.g., `libqhull-dev` on
Debian/Ubuntu, `qhull-devel` on Fedora, `brew install qhull` on macOS).

```bash
pip install focus_cpt-0.1.10.tar.gz    # or: pip install focus-cpt
pip install numpy pandas scipy plotnine astropy
python replication.py
```

## Data

The following data files are included in `paper_data/`:

- `cle_data.rds`: game logs of the 2089 regular season games of the Cleveland
  Cavaliers from the 1999-00 to the 2024-25 NBA season (Section 5.1), as an R
  data frame with columns `yearSeason`, `slugSeason`, `typeSeason`, `dateGame`,
  `nameTeam`, `slugTeam`, `isWin`, `ptsTeam`, `plusminusTeam` and `monthGame`.
  The game logs were retrieved with the `game_logs()` function of the
  **nbastatR** package (<https://github.com/abresler/nbastatR>), which is not
  required for replication.
- `glg_tte_n2_bn250814432_v00.fit`: time-tagged event data of the NaI detector
  n2 of the Fermi Gamma-ray Burst Monitor for the trigger bn250814432
  (Section 5.2), publicly available from
  <https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2025/bn250814432/current/>.
- `example_trace.pkl`: a pickled Python dictionary with the fluorescence trace
  (`trace`), its sampling times (`time`) and the true spike times (`spikes`) of
  one recording of the spikefinder challenge data (Section 5.3), available at
  <https://github.com/codeneuro/spikefinder>.

## Building the Manuscript (Optional)

The manuscript is written in Quarto and uses the
[JSS Quarto template](https://github.com/quarto-journals/jss). Building it
requires Quarto (>= 1.3.0), LaTeX, and the R and Python dependencies above (the
Python chunks are run through the R package reticulate).

```bash
quarto add quarto-journals/jss   # install the JSS template in this folder
make all                         # check dependencies and compile jss_paper.pdf
```

Other Make targets:

- `make check-deps`: check all system and package dependencies;
- `make install-deps`: install missing R and Python packages;
- `make compile`: compile the manuscript (`jss_paper.tex` and `jss_paper.pdf`);
- `make replicate`: run `replication.R` and `replication.py`;
- `make replication-scripts`: regenerate `replication.R` and `replication.py`
  from `jss_paper.qmd`;
- `make clean` and `make rebuild`: remove the Quarto cache and recompile.

## License

The content of this folder (including all code and replication materials) is
licensed under the Creative Commons Attribution 4.0 International License
(CC BY 4.0), see <https://creativecommons.org/licenses/by/4.0/>.
