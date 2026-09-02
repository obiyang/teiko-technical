# Teiko Technical: Cell-count analysis

This project turns the provided cell-count table into a normalized SQLite database,
reproducible analytical outputs, statistical comparisons, and a local interactive
dashboard. All paths are resolved from the repository itself; no data is uploaded or
sent to an external service.

## Reproduce in GitHub Codespaces or locally

From this directory, run:

```bash
make setup
make pipeline
make dashboard
```

The pinned dependency set supports Python 3.10–3.13 and was verified in a fresh
Python 3.11 virtual environment.

`make setup` creates `.venv` and installs the versions pinned in
`requirements.txt`. `make pipeline` rebuilds `cell_counts.db`, runs the full
analysis, and writes all results under `outputs/`. `make dashboard` starts the
Streamlit UI; use the local URL printed by Streamlit (normally
[`http://localhost:8501`](http://localhost:8501)). In Codespaces, open the forwarded port **8501** and use
the URL Codespaces provides.

The dashboard is intentionally local-only: this repository does not include a
public hosted deployment or a public dashboard URL.

## Data model

The SQLite schema separates reusable dimensions (`projects`, `conditions`,
`sexes`, `treatments`, `responses`, `sample_types`, and `cell_types`) from
`subjects`, `samples`, and the long-form `cell_counts` fact table. Foreign keys,
unique constraints, non-negative count checks, and indexes enforce integrity and
support common cohort filters.

The long-form count table is the key extensibility choice: adding a cell type means
adding one row to `cell_types` and count rows, rather than altering a wide table.
The subject/sample split avoids repeating demographics and keeps longitudinal
samples explicit. `dataset_loads` records the source filename, SHA-256 digest, row
count, and load time for provenance. The loader validates the source and atomically
replaces the database only after a successful load.

## Analysis choices

`outputs/sample_cell_frequencies.csv` contains one row per sample and cell type,
including the exact source population key (`population`: `b_cell`, `cd8_t_cell`,
`cd4_t_cell`, `nk_cell`, or `monocyte`), `count`, sample `total_count`, and
`percentage`.

The response analysis uses only **melanoma + miraclib + PBMC** records. It first
averages each cell-type percentage across samples for each subject, so repeated
longitudinal measurements are not incorrectly treated as independent observations.
Responders (`response=yes`) and non-responders (`response=no`) are compared with a
two-sided Mann-Whitney U test. This non-parametric test is suitable for bounded
relative-frequency data without a normality assumption. Benjamini-Hochberg false
discovery rate correction is applied across the five cell-type tests. The full
table—including group sizes, means, medians, median difference, U statistic,
rank-biserial effect size, raw p-value, adjusted p-value, and significance flag—is
written to `outputs/response_statistics.csv`. One boxplot per cell type is written
to `outputs/plots/`. A concise, generated statement of significant populations and
the baseline counts is available in `outputs/analysis_summary.md` after every run.

The requested baseline query is fixed to **melanoma + PBMC + time 0 + miraclib**.
Its outputs are:

- `outputs/baseline_matched_samples.csv`: all matching samples
- `outputs/baseline_samples_by_project.csv`: sample count by project
- `outputs/baseline_subjects_by_response.csv`: distinct subjects by response
- `outputs/baseline_subjects_by_sex.csv`: distinct subjects by sex

## Code structure

```text
.
├── data/cell-count.csv             # private, bundled input
├── load_data.py                    # no-argument validated SQLite loader
├── src/teiko_analysis/pipeline.py  # summaries, query, statistics, plots
├── dashboard/app.py                # interactive Streamlit dashboard
├── tests/                          # loader and analysis tests
├── outputs/                        # generated CSV tables and PNG plots
├── Makefile                        # setup / pipeline / dashboard entry points
└── requirements.txt                # pinned dependencies
```

Run tests after the pipeline with:

```bash
PYTHONPATH=src pytest -q
```

The pipeline performs its own row-count, percentage-sum, and SQLite integrity
checks in addition to the test suite.
