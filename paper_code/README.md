# focus JSS Paper - Reproducible Build

This folder contains the reproducible source code for the JSS paper on `focus`: Fast Online Changepoint Detection in R and Python.

## Quick Start

The easiest way to build this paper is using the provided Makefile:

```bash
make all
```

This will:
1. Check for all system and package dependencies
2. Prompt you to install the JSS Quarto template (if not already installed)
3. Install any missing R and Python packages (with your confirmation)
4. Compile the paper to generate `jss_paper.tex` and `jss_paper.pdf`

## System Requirements

### Core Software

Before building, ensure you have installed:

- **Quarto** (≥ 1.3.0) - Document rendering engine
  - Installation: https://quarto.org/docs/get-started/
  - **Note:** The build process will check for and guide you through setting up the [JSS Quarto template](https://github.com/quarto-journals/jss) if needed.
  
- **R** (≥ 4.0) - Statistical computing language
  - Installation: https://www.r-project.org/
  
- **Python** (≥ 3.8) - Programming language
  - Installation: https://www.python.org/
  
- **LaTeX** - Required for PDF generation
  - On Ubuntu/Debian: `sudo apt-get install texlive-latex-extra`
  - On macOS: `brew install mactex` or install MacTeX
  - On Windows: Install MiKTeX or TeX Live

### R Package Dependencies

The following R packages are required:

- `focus` - The focus changepoint detection package
- `ggplot2` - Graphics system
- `furrr` - Parallel functional programming
- `purrr` - Functional programming tools
- `dplyr` - Data manipulation

These will be automatically installed by `make install-deps` if missing.

### Python Package Dependencies

The following Python packages are required:

- `focus-cpt` - The focus changepoint detection package (Python version)
- `pandas` - Data manipulation
- `numpy` - Numerical computing
- `scipy` - Scientific computing
- `plotnine` - Grammar of graphics plotting
- `astropy` - Astronomy/astrophysics library

These will be automatically installed by `make install-deps` if missing.

## Data Files

The following data files are included in this repository under paper_data/:

- `example_trace.pkl` - Example trace from the codeneuro challenge dataset. The original source can be found on GitHub at [this link](https://github.com/codeneuro/spikefinder}).
- `glg_tte_n2_bn250814432_v00.fit` - GRB detection data from Fermi-GBM detector. The exact data file can be found on [this page](https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2025/bn250814432//current/).
- `nba_data.rds` - NBA cavaliers data, obtained using the `nbastatR` package. Code to generate this file is included in the paper source (see `jss_paper.qmd`), but the pre-generated file is included here to speed up compilation. 

## Building the Paper

### Option 1: Using Make (Recommended)

```bash
# Check dependencies
make check-deps

# Install missing dependencies
make install-deps

# Build the paper
make compile

# Or do everything in one command
make all
```

### Option 2: Manual Build

If you don't have `make` available, you can build directly with Quarto:

```bash
# First, ensure all dependencies are installed
# Install R packages:
Rscript -e "install.packages(c('focus', 'ggplot2', 'furrr', 'purrr', 'dplyr'))"

# Install Python packages:
pip install focus-cpt pandas numpy scipy plotnine astropy

# Then render the paper
quarto render jss_paper.qmd --to jss-pdf
```

## Make Targets

- `make all` - Full pipeline: check deps, prompt to install missing, compile
- `make check-deps` - Check all dependencies (stops if JSS template missing)
- `make install-deps` - Install missing R and Python packages
- `make compile` - Compile the paper
- `make clean` - Remove generated cache and files directories
- `make rebuild` - Clean and recompile
- `make help` - Show help message

**Note on Installation:** The build system checks dependencies and provides clear instructions rather than silently installing packages. This ensures transparency and lets reviewers control what gets installed on their systems.

## Output Files

After successful compilation, you will have:

- **`jss_paper.tex`** - LaTeX source file (ready for submission to JSS)
- **`jss_paper.pdf`** - Compiled PDF version
- **`jss_paper_files/`** - Directory containing any embedded figures and assets
- **`jss_paper_cache/`** - Quarto cache directory (can be removed after build)

## Source Document

The paper is written in **Quarto** format (`.qmd`), a markup language that combines:

- **Markdown** for text and formatting
- **YAML** for metadata (title, authors, abstract, etc.)
- **R code blocks** for R examples and figures
- **Python code blocks** for Python examples and figures

The document is processed by Quarto, which:
1. Executes all code blocks (via `knitr`)
2. Generates figures and output
3. Converts to LaTeX format
4. Uses pdflatex to generate the final PDF

This paper uses the [JSS Quarto template](https://github.com/quarto-journals/jss).

## Troubleshooting

### JSS Quarto Extension Not Found

If the build fails with "Unable to read the extension 'jss'", you need to install the JSS template:

```bash
quarto add quarto-journals/jss
```

Then retry:
```bash
make compile
```

If the extension was installed correctly, in the same directory you should have the folder `_extensions/quarto-journals/jss` at the same level of `jss_paper.qmd`.

### LaTeX not found
If you get an error about pdflatex not found:
- **Ubuntu/Debian**: `sudo apt-get install texlive-latex-extra`
- **macOS**: `brew install mactex` (or download MacTeX from https://www.tug.org/mactex/)
- **Windows**: Install MiKTeX from https://miktex.org/

### R package installation fails
To install R packages manually:
```bash
Rscript -e "install.packages(c('focus', 'ggplot2', 'furrr', 'purrr', 'dplyr'), repos='https://cloud.r-project.org/')"
```

### Python package installation fails
To install Python packages manually:
```bash
pip install --upgrade pip
pip install focus-cpt pandas numpy scipy plotnine astropy
```

### Quarto not found
If Quarto is not in your PATH, download and install from: https://quarto.org/docs/get-started/


### Quitting from lines 1019-1060 [load-nba-data] (jss_paper.qmd)

If you see something like:

```bash
Error:
! invalid connection... Please make sure not to call closeAllConnections().
Backtrace:
```

This can happen in case the package `nbastatR` has downloaded all data, but did not close the connection. Simply restarting the `make compile` should fix this issue. 

## License

The paper in this folder (including all code and reproducibility materials) is licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0), in line with ArXiv and the Journal of Statistical Software. For full license terms, see https://creativecommons.org/licenses/by/4.0/

## Citation

If you use this code or paper, please cite:

```bibtex
@article{romano2026focus,
  title={focus and focus-cpt: Fast Online Changepoint Detection in R and Python},
  author={Romano, Gaetano and Ward, Kes and Fan, Yuntang and Rigaill, Guillem and Runge, Vincent and Eckley, Idris A and Fearnhead, Paul},
  journal={},
  year={2026}
}
```
