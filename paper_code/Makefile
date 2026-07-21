.PHONY: all check-deps install-deps compile clean help rebuild

# Colors for output
RED = \033[0;31m
GREEN = \033[0;32m
YELLOW = \033[1;33m
NC = \033[0m # No Color

help:
	@echo "$(GREEN)focus JSS Paper Build$(NC)"
	@echo ""
	@echo "Available targets:"
	@echo "  $(YELLOW)make all$(NC)          - Check dependencies, install if needed, and compile the paper"
	@echo "  $(YELLOW)make check-deps$(NC)    - Check all system and package dependencies"
	@echo "  $(YELLOW)make install-deps$(NC)  - Install missing dependencies"
	@echo "  $(YELLOW)make compile$(NC)       - Compile the paper (generate .tex and .pdf)"
	@echo "  $(YELLOW)make clean$(NC)         - Remove generated files (_cache and _files directories)"
	@echo "  $(YELLOW)make rebuild$(NC)       - Clean and recompile the paper"
	@echo "  $(YELLOW)make help$(NC)          - Show this help message"

all: check-deps install-deps compile
	@echo "$(GREEN)✓ Paper compiled successfully!$(NC)"

# Check all dependencies
check-deps: check-system-deps check-r-deps check-python-deps
	@echo "$(GREEN)✓ All dependency checks completed$(NC)"

# Check system dependencies (Quarto, R, Python)
check-system-deps:
	@echo "$(YELLOW)Checking system dependencies...$(NC)"
	@command -v quarto >/dev/null 2>&1 || { \
		echo "$(RED)✗ Quarto is not installed$(NC)"; \
		echo ""; \
		echo "$(YELLOW)To install Quarto:$(NC)"; \
		echo "  Ubuntu/Debian: Download the .deb from https://quarto.org/docs/download/"; \
		echo "                 sudo dpkg -i quarto-*.deb"; \
		echo "  Fedora/RHEL:   Download the .rpm from https://quarto.org/docs/download/"; \
		echo "                 sudo rpm -i quarto-*.rpm"; \
		echo "  macOS:         brew install --cask quarto"; \
		echo "                 (or download the .pkg from https://quarto.org/docs/download/)"; \
		echo "  Windows:       Download the installer from https://quarto.org/docs/download/"; \
		echo ""; \
		exit 1; \
	}
	@echo "$(GREEN)  ✓ Quarto found: $$(quarto --version)$(NC)"
	@command -v R >/dev/null 2>&1 || { echo "$(RED)✗ R is not installed$(NC)"; exit 1; }
	@echo "$(GREEN)  ✓ R found$(NC)"
	@command -v python3 >/dev/null 2>&1 || command -v python >/dev/null 2>&1 || { echo "$(RED)✗ Python is not installed$(NC)"; exit 1; }
	@echo "$(GREEN)  ✓ Python found$(NC)"
	@command -v pdflatex >/dev/null 2>&1 || { echo "$(RED)✗ LaTeX/pdflatex is not installed$(NC)"; exit 1; }
	@echo "$(GREEN)  ✓ LaTeX found$(NC)"
	@if [ -d "_extensions/quarto-journals/jss" ]; then \
		echo "$(GREEN)  ✓ JSS extension found$(NC)"; \
	else \
		echo "$(RED)✗ JSS extension not found$(NC)"; \
		echo ""; \
		echo "$(YELLOW)To set up the JSS template, run:$(NC)"; \
		echo "  quarto add quarto-journals/jss"; \
		echo ""; \
		exit 1; \
	fi
# Check R package dependencies
check-r-deps:
	@echo "$(YELLOW)Checking R packages...$(NC)"
	@Rscript --vanilla --slave -e "packages <- c('focus', 'ggplot2', 'furrr', 'purrr', 'dplyr'); missing <- packages[!sapply(packages, function(p) requireNamespace(p, quietly=TRUE))]; if(length(missing) > 0) { cat('Missing R packages:', paste(missing, collapse=', '), '\n'); quit(status=1) } else { cat('All R packages found\n') }" 2>/dev/null || { echo "$(RED)✗ Some R packages are missing. Running 'make install-deps'$(NC)"; }

# Check Python package dependencies
check-python-deps:
	@echo "$(YELLOW)Checking Python packages...$(NC)"
	@which python3 >/dev/null 2>&1 && PYTHON=python3 || PYTHON=python; \
	$$PYTHON -c "import focus_cpt, pandas, numpy, scipy, plotnine, astropy; print('All Python packages found')" 2>/dev/null || { echo "$(RED)✗ Some Python packages are missing. Running 'make install-deps'$(NC)"; }

# Install missing dependencies
install-deps: install-r-deps install-python-deps
	@echo "$(GREEN)✓ Dependency installation completed$(NC)"

install-r-deps:
	@echo "$(YELLOW)Installing/updating R packages...$(NC)"
	@Rscript --vanilla --slave -e "\
		packages <- c('focus', 'ggplot2', 'furrr', 'purrr', 'dplyr'); \
		failed <- c(); \
		for(pkg in packages) { \
			if(!requireNamespace(pkg, quietly=TRUE)) { \
				cat('  Installing', pkg, '...\n'); \
				tryCatch( \
					install.packages(pkg, repos='https://cloud.r-project.org/', quiet=TRUE), \
					error = function(e) { failed <<- c(failed, pkg) } \
				); \
				if(!requireNamespace(pkg, quietly=TRUE)) failed <- c(failed, pkg); \
			} \
		}; \
		if(length(failed) > 0) { \
			cat('Failed to install R packages:', paste(failed, collapse=', '), '\n'); \
			quit(status=1) \
		} else { cat('R packages ready\n') }" \
	|| { echo "$(RED)✗ Failed to install one or more R packages (see above)$(NC)"; exit 1; }

install-python-deps:
	@echo "$(YELLOW)Installing/updating Python packages...$(NC)"
	@which python3 >/dev/null 2>&1 && PYTHON=python3 || PYTHON=python; \
	INSTALL_OUT=$$($$PYTHON -m pip install focus-cpt pandas numpy scipy plotnine astropy 2>&1); \
	INSTALL_RC=$$?; \
	if [ $$INSTALL_RC -ne 0 ]; then \
		echo "$(RED)✗ pip install failed. Checking which packages are missing...$(NC)"; \
		for pkg in focus-cpt pandas numpy scipy plotnine astropy; do \
			import_name=$$(echo $$pkg | sed 's/-/_/g'); \
			$$PYTHON -c "import $$import_name" 2>/dev/null || \
				echo "$(RED)  ✗ Missing: $$pkg$(NC)"; \
		done; \
		echo ""; \
		echo "$(YELLOW)Note: the focus-cpt package requires Qhull dev libraries.$(NC)"; \
		echo "  Ubuntu/Debian: sudo apt-get install libqhull-dev"; \
		echo "  Fedora/RHEL:   sudo dnf install qhull-devel"; \
		echo "  macOS:         brew install qhull"; \
		echo ""; \
		exit 1; \
	fi; \
	MISSING=""; \
	for pkg in focus_cpt pandas numpy scipy plotnine astropy; do \
		$$PYTHON -c "import $$pkg" 2>/dev/null || MISSING="$$MISSING $$pkg"; \
	done; \
	if [ -n "$$MISSING" ]; then \
		echo "$(RED)✗ The following Python packages could not be imported after install:$$MISSING$(NC)"; \
		echo ""; \
		echo "$(YELLOW)Note: the focus-cpt package requires Qhull dev libraries.$(NC)"; \
		echo "  Ubuntu/Debian: sudo apt-get install libqhull-dev"; \
		echo "  Fedora/RHEL:   sudo dnf install qhull-devel"; \
		echo "  macOS:         brew install qhull"; \
		echo ""; \
		exit 1; \
	fi; \
	echo "Python packages ready"

# Compile the paper
compile: check-system-deps
	@echo "$(YELLOW)Compiling paper (jss_paper.qmd -> .tex and .pdf)...$(NC)"
	@quarto render jss_paper.qmd --to jss-pdf
	@if [ -f jss_paper.tex ]; then \
		echo "$(GREEN)✓ Generated jss_paper.tex$(NC)"; \
	else \
		echo "$(RED)✗ Failed to generate jss_paper.tex$(NC)"; \
		exit 1; \
	fi
	@if [ -f jss_paper.pdf ]; then \
		echo "$(GREEN)✓ Generated jss_paper.pdf$(NC)"; \
	else \
		echo "$(RED)✗ Failed to generate jss_paper.pdf$(NC)"; \
		exit 1; \
	fi

# Clean generated files
clean:
	@echo "$(YELLOW)Cleaning generated files...$(NC)"
	@rm -rf jss_paper_cache jss_paper_files
	@echo "$(GREEN)✓ Cleaned$(NC)"

# Force recompile even if files exist
rebuild: clean compile
	@echo "$(GREEN)✓ Paper rebuilt successfully!$(NC)"