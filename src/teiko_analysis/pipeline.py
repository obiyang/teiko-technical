"""Reproducible summaries, cohort query, statistics, and plots."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parents[2]
DB_PATH = ROOT / "cell_counts.db"
OUTPUT_DIR = ROOT / "outputs"
PLOT_DIR = OUTPUT_DIR / "plots"


FREQUENCY_QUERY = """
SELECT
    p.name AS project,
    s.subject_code AS subject,
    c.name AS condition,
    sx.label AS sex,
    t.name AS treatment,
    r.label AS response,
    sm.sample_code AS sample,
    st.name AS sample_type,
    sm.time_from_treatment_start,
    ct.source_column AS population,
    ct.name AS cell_type,
    cc.count,
    SUM(cc.count) OVER (PARTITION BY sm.sample_id) AS total_count,
    100.0 * cc.count / NULLIF(SUM(cc.count) OVER (PARTITION BY sm.sample_id), 0)
        AS percentage
FROM samples sm
JOIN subjects s USING (subject_id)
JOIN projects p USING (project_id)
JOIN conditions c USING (condition_id)
JOIN sexes sx USING (sex_id)
JOIN treatments t USING (treatment_id)
LEFT JOIN responses r USING (response_id)
JOIN sample_types st USING (sample_type_id)
JOIN cell_counts cc USING (sample_id)
JOIN cell_types ct USING (cell_type_id)
ORDER BY sm.sample_code, ct.cell_type_id
"""


BASELINE_COHORT_QUERY = """
SELECT
    p.name AS project,
    s.subject_code AS subject,
    c.name AS condition,
    sx.label AS sex,
    t.name AS treatment,
    r.label AS response,
    sm.sample_code AS sample,
    st.name AS sample_type,
    sm.time_from_treatment_start
FROM samples sm
JOIN subjects s USING (subject_id)
JOIN projects p USING (project_id)
JOIN conditions c USING (condition_id)
JOIN sexes sx USING (sex_id)
JOIN treatments t USING (treatment_id)
LEFT JOIN responses r USING (response_id)
JOIN sample_types st USING (sample_type_id)
WHERE LOWER(c.name) = 'melanoma'
  AND LOWER(st.name) = 'pbmc'
  AND sm.time_from_treatment_start = 0
  AND LOWER(t.name) = 'miraclib'
ORDER BY p.name, s.subject_code, sm.sample_code
"""


def benjamini_hochberg(p_values: list[float] | np.ndarray) -> np.ndarray:
    """Return Benjamini-Hochberg adjusted p-values in original order."""
    values = np.asarray(p_values, dtype=float)
    if values.ndim != 1:
        raise ValueError("p_values must be one-dimensional")
    if np.any(~np.isfinite(values)) or np.any((values < 0) | (values > 1)):
        raise ValueError("p_values must be finite values between 0 and 1")
    if len(values) == 0:
        return values.copy()

    order = np.argsort(values)
    ranked = values[order]
    adjusted_ranked = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted_ranked = np.clip(adjusted_ranked, 0, 1)
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = adjusted_ranked
    return adjusted


def load_frequencies(connection: sqlite3.Connection) -> pd.DataFrame:
    data = pd.read_sql_query(FREQUENCY_QUERY, connection)
    expected_rows = pd.read_sql_query(
        "SELECT COUNT(*) AS n FROM samples", connection
    ).iloc[0, 0] * pd.read_sql_query(
        "SELECT COUNT(*) AS n FROM cell_types", connection
    ).iloc[0, 0]
    if len(data) != expected_rows:
        raise RuntimeError(
            f"Frequency summary has {len(data)} rows; expected {expected_rows}"
        )
    nonzero = data["total_count"] > 0
    if nonzero.any():
        percentage_sums = (
            data.loc[nonzero]
            .groupby("sample")["percentage"]
            .sum()
            .to_numpy(dtype=float)
        )
        if not np.allclose(percentage_sums, 100.0):
            raise RuntimeError("Nonzero sample percentages do not sum to 100")
    return data


def compute_response_statistics(frequencies: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare subject-level mean frequencies for responders and non-responders.

    Subject-level aggregation prevents repeated longitudinal samples from being
    treated as independent observations. A two-sided Mann-Whitney U test is used
    because relative frequencies are bounded and need not be normally distributed.
    """
    mask = (
        frequencies["condition"].str.casefold().eq("melanoma")
        & frequencies["treatment"].str.casefold().eq("miraclib")
        & frequencies["sample_type"].str.casefold().eq("pbmc")
        & frequencies["response"].str.casefold().isin(["yes", "no"])
    )
    subject_level = (
        frequencies.loc[mask]
        .groupby(
            ["project", "subject", "response", "population", "cell_type"],
            as_index=False,
        )["percentage"]
        .mean()
    )
    if subject_level.empty:
        raise RuntimeError("No samples match the response-comparison cohort")

    rows: list[dict[str, float | int | str]] = []
    for population in sorted(subject_level["population"].unique()):
        cell_data = subject_level[subject_level["population"] == population]
        cell_type = cell_data["cell_type"].iloc[0]
        responders = cell_data.loc[
            cell_data["response"].str.casefold() == "yes", "percentage"
        ].to_numpy()
        nonresponders = cell_data.loc[
            cell_data["response"].str.casefold() == "no", "percentage"
        ].to_numpy()
        if len(responders) == 0 or len(nonresponders) == 0:
            raise RuntimeError(f"Both response groups are required for {cell_type}")

        test = mannwhitneyu(
            responders, nonresponders, alternative="two-sided", method="auto"
        )
        denominator = len(responders) * len(nonresponders)
        rows.append(
            {
                "population": population,
                "cell_type": cell_type,
                "n_responder_subjects": len(responders),
                "n_nonresponder_subjects": len(nonresponders),
                "responder_mean_percentage": responders.mean(),
                "nonresponder_mean_percentage": nonresponders.mean(),
                "responder_median_percentage": np.median(responders),
                "nonresponder_median_percentage": np.median(nonresponders),
                "median_difference_percentage_points": (
                    np.median(responders) - np.median(nonresponders)
                ),
                "mann_whitney_u": float(test.statistic),
                "rank_biserial_correlation": 2 * float(test.statistic) / denominator - 1,
                "p_value": float(test.pvalue),
            }
        )

    statistics = pd.DataFrame(rows)
    statistics["p_value_bh"] = benjamini_hochberg(statistics["p_value"].to_numpy())
    statistics["significant_fdr_0_05"] = statistics["p_value_bh"] < 0.05
    statistics["test"] = "two-sided Mann-Whitney U on subject-level mean percentage"
    statistics["multiple_testing"] = "Benjamini-Hochberg across cell types"
    return statistics, subject_level


def write_response_plots(subject_level: pd.DataFrame, plot_dir: Path) -> None:
    plot_dir.mkdir(parents=True, exist_ok=True)
    colors = {"no": "#E07A5F", "yes": "#3D85C6"}
    for population in sorted(subject_level["population"].unique()):
        cell_data = subject_level[subject_level["population"] == population].copy()
        cell_type = cell_data["cell_type"].iloc[0]
        groups = []
        labels = []
        for response, label in [("no", "Non-responder"), ("yes", "Responder")]:
            values = cell_data.loc[
                cell_data["response"].str.casefold() == response, "percentage"
            ].to_numpy()
            groups.append(values)
            labels.append(label)

        fig, ax = plt.subplots(figsize=(6.4, 4.8))
        box = ax.boxplot(groups, tick_labels=labels, patch_artist=True, widths=0.55)
        for patch, response in zip(box["boxes"], ["no", "yes"]):
            patch.set_facecolor(colors[response])
            patch.set_alpha(0.75)
        rng = np.random.default_rng(42)
        for index, values in enumerate(groups, start=1):
            jitter = rng.normal(index, 0.035, size=len(values))
            ax.scatter(jitter, values, s=12, alpha=0.35, color="#263238")
        ax.set_title(f"{cell_type}: response groups")
        ax.set_ylabel("Subject-level mean relative frequency (%)")
        ax.grid(axis="y", alpha=0.2)
        fig.tight_layout()
        filename = population + ".png"
        fig.savefig(plot_dir / filename, dpi=160)
        plt.close(fig)


def write_baseline_query(
    connection: sqlite3.Connection, output_dir: Path
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cohort = pd.read_sql_query(BASELINE_COHORT_QUERY, connection)
    cohort.to_csv(output_dir / "baseline_matched_samples.csv", index=False)
    by_project = (
        cohort.groupby("project", dropna=False)
        .size()
        .rename("sample_count")
        .reset_index()
    )
    by_project.to_csv(output_dir / "baseline_samples_by_project.csv", index=False)
    subject_keys = cohort.assign(
        subject_key=cohort["project"].astype(str) + "::" + cohort["subject"].astype(str)
    )
    by_response = (
        subject_keys.groupby("response", dropna=False)["subject_key"]
        .nunique()
        .rename("unique_subject_count")
        .reset_index()
    )
    by_response.to_csv(output_dir / "baseline_subjects_by_response.csv", index=False)
    by_sex = (
        subject_keys.groupby("sex", dropna=False)["subject_key"]
        .nunique()
        .rename("unique_subject_count")
        .reset_index()
    )
    by_sex.to_csv(output_dir / "baseline_subjects_by_sex.csv", index=False)
    return cohort, by_project, by_response, by_sex


def _markdown_table(data: pd.DataFrame, columns: list[str]) -> list[str]:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in data[columns].itertuples(index=False, name=None):
        formatted = []
        for value in row:
            if isinstance(value, (float, np.floating)):
                formatted.append(f"{value:.6g}")
            else:
                formatted.append(str(value))
        lines.append("| " + " | ".join(formatted) + " |")
    return lines


def write_analysis_summary(
    frequencies: pd.DataFrame,
    statistics: pd.DataFrame,
    subject_level: pd.DataFrame,
    cohort: pd.DataFrame,
    by_project: pd.DataFrame,
    by_response: pd.DataFrame,
    by_sex: pd.DataFrame,
    output_path: Path,
) -> None:
    significant = statistics.loc[statistics["significant_fdr_0_05"], "population"].tolist()
    significance_text = ", ".join(significant) if significant else "None"
    response_subjects = subject_level[["project", "subject"]].drop_duplicates().shape[0]
    lines = [
        "# Cell-count analysis summary",
        "",
        "## Sample and population summary",
        "",
        f"- Samples: {frequencies['sample'].nunique():,}",
        f"- Populations: {frequencies['population'].nunique():,}",
        "- Each nonzero sample's population percentages sum to 100%.",
        "",
        "## Responder comparison",
        "",
        "Cohort: melanoma + miraclib + PBMC. Percentages were averaged within each "
        "subject before a two-sided Mann-Whitney U test; Benjamini-Hochberg correction "
        "was applied across populations.",
        "",
        f"- Independent subjects analyzed: {response_subjects:,}",
        f"- Populations significant at FDR < 0.05: **{significance_text}**",
        "",
    ]
    lines.extend(
        _markdown_table(
            statistics,
            [
                "population",
                "n_responder_subjects",
                "n_nonresponder_subjects",
                "p_value",
                "p_value_bh",
                "significant_fdr_0_05",
            ],
        )
    )
    lines.extend(
        [
            "",
            "## Baseline cohort query",
            "",
            "Filters: melanoma + PBMC + time_from_treatment_start = 0 + miraclib.",
            "",
            f"- Matching samples: {len(cohort):,}",
            "",
            "### Samples by project",
            "",
        ]
    )
    lines.extend(_markdown_table(by_project, ["project", "sample_count"]))
    lines.extend(["", "### Distinct subjects by response", ""])
    lines.extend(_markdown_table(by_response, ["response", "unique_subject_count"]))
    lines.extend(["", "### Distinct subjects by sex", ""])
    lines.extend(_markdown_table(by_sex, ["sex", "unique_subject_count"]))
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_pipeline(db_path: Path = DB_PATH, output_dir: Path = OUTPUT_DIR) -> None:
    if not db_path.is_file():
        raise FileNotFoundError(f"Database not found: {db_path}. Run load_data.py first.")
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = output_dir / "plots"
    with sqlite3.connect(db_path) as connection:
        frequencies = load_frequencies(connection)
        frequencies.to_csv(output_dir / "sample_cell_frequencies.csv", index=False)
        statistics, subject_level = compute_response_statistics(frequencies)
        statistics.to_csv(output_dir / "response_statistics.csv", index=False)
        subject_level.to_csv(output_dir / "response_subject_level_frequencies.csv", index=False)
        write_response_plots(subject_level, plot_dir)
        cohort, by_project, by_response, by_sex = write_baseline_query(
            connection, output_dir
        )
        write_analysis_summary(
            frequencies,
            statistics,
            subject_level,
            cohort,
            by_project,
            by_response,
            by_sex,
            output_dir / "analysis_summary.md",
        )

    print(f"Wrote analysis outputs to {output_dir}")


def main() -> None:
    run_pipeline()


if __name__ == "__main__":
    main()
