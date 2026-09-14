"""Generate replication.R and replication.py from the code chunks of jss_paper.qmd.

Every evaluated R and Python chunk of the manuscript is copied verbatim, in
order of appearance, including the chunks whose code is hidden in the
manuscript. Chunks that are not evaluated (eval: false) and the R chunks that
only set the code prompt are skipped. The section headings of the manuscript
are added as comments.

The only lines added to the code are those needed to run it as a script and
save the figures: the creation of the folder figures/ and, for each figure chunk
(label starting with "fig-"), the lines saving the figure to
figures/<label>.pdf. For R figures, the chunk code is wrapped in pdf() and
dev.off(), with the fig-width and fig-height of the chunk. For Python figures,
p.save() is added after the chunk, and the plot is assigned to p when the chunk
ends with an unnamed plot expression.

Usage (from this folder):
    python make_replication_scripts.py    # or: make replication-scripts
"""

import pathlib
import re

HERE = pathlib.Path(__file__).resolve().parent
QMD = HERE / "jss_paper.qmd"

R_HEADER = """\
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
"""

PY_HEADER = '''\
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
'''


def parse_qmd(lines):
    """Return the headings and the evaluated code chunks, in order."""
    items = []
    i = 0
    while i < len(lines):
        chunk = re.match(r"^```\{(r|python)([^}]*)\}\s*$", lines[i])
        heading = re.match(r"^(#{2,4}) (.*?)\s*(\{.*\})?\s*$", lines[i])
        if chunk:
            engine, header = chunk.group(1), chunk.group(2)
            body = []
            i += 1
            while not lines[i].startswith("```"):
                body.append(lines[i])
                i += 1
            opts = {}
            for line in body:
                opt = re.match(r"#\|\s*([\w-]+):\s*(.*)$", line)
                if opt:
                    opts[opt.group(1)] = opt.group(2).strip()
            code = [line for line in body if not line.startswith("#|")]
            not_evaluated = "eval=FALSE" in header.replace(" ", "") or opts.get("eval") == "false"
            only_prompt = all("options(prompt" in line or line.strip() == "" or line.strip().startswith("#")
                              for line in code)
            if not not_evaluated and not (engine == "r" and only_prompt):
                label = opts.get("label") or (header.strip(" ,") or None)
                items.append(("chunk", engine, opts, label, code))
        elif heading:
            items.append(("heading", len(heading.group(1)), heading.group(2)))
        i += 1
    return items


def strip_blank_ends(code):
    while code and code[0].strip() == "":
        code = code[1:]
    while code and code[-1].strip() == "":
        code = code[:-1]
    return code


def build(items, engine):
    if engine == "r":
        comment = "##"
        out = [R_HEADER, 'dir.create("figures", showWarnings = FALSE)', ""]
    else:
        comment = "#"
        out = [PY_HEADER, "import os", "", 'os.makedirs("figures", exist_ok=True)', ""]

    pending = []  # headings not yet written, as (level, text)
    for item in items:
        if item[0] == "heading":
            level, text = item[1], item[2]
            pending = [h for h in pending if h[0] < level] + [(level, text)]
            continue
        _, chunk_engine, opts, label, code = item
        if chunk_engine != engine:
            continue
        for level, text in pending:
            if level == 2:
                rule = f"{comment} " + "=" * 74
                out += ["", rule, f"{comment} {text}", rule, ""]
            else:
                out += [f"{comment} ---- {text} " + "-" * max(3, 68 - len(text)), ""]
        pending = []

        code = strip_blank_ends(code)
        if label is not None and label.startswith("fig-"):
            if engine == "r":
                width, height = opts.get("fig-width", "7"), opts.get("fig-height", "5")
                out += [f"{comment} Figure: {label}",
                        f'pdf("figures/{label}.pdf", width = {width}, height = {height})']
                out += code + ["invisible(dev.off())", ""]
            else:
                if not any(re.match(r"^p\s*=", line) for line in code):
                    last = max(k for k, line in enumerate(code) if line.startswith("("))
                    code = code[:last] + ["p = " + code[last]] + code[last + 1:]
                out += [f"{comment} Figure: {label}"] + code
                out += [f'p.save("figures/{label}.pdf", verbose=False)', ""]
        else:
            out += code + [""]
    return "\n".join(out).rstrip() + "\n"


def main():
    items = parse_qmd(QMD.read_text().splitlines())
    for engine, name in (("r", "replication.R"), ("python", "replication.py")):
        (HERE / name).write_text(build(items, engine))
        print(f"Written {name}")


if __name__ == "__main__":
    main()
