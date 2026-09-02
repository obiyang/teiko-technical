"""Interactive Streamlit dashboard for the normalized cell-count database."""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import pandas as pd
import streamlit as st


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from teiko_analysis.pipeline import FREQUENCY_QUERY  # noqa: E402


DB_PATH = ROOT / "cell_counts.db"
OUTPUT_DIR = ROOT / "outputs"


@st.cache_data(show_spinner=False)
def load_data(db_mtime_ns: int) -> pd.DataFrame:
    del db_mtime_ns
    with sqlite3.connect(DB_PATH) as connection:
        return pd.read_sql_query(FREQUENCY_QUERY, connection)


def multiselect_filter(data: pd.DataFrame, column: str, label: str) -> list:
    options = sorted(data[column].dropna().unique().tolist())
    return st.sidebar.multiselect(label, options, default=options)


st.set_page_config(page_title="Cell Count Explorer", page_icon="🧫", layout="wide")
st.title("Cell Count Explorer")
st.caption("Filter samples and inspect relative cell-type composition interactively.")

if not DB_PATH.is_file():
    st.error("Database not found. Run `make pipeline` from the project root first.")
    st.stop()

data = load_data(DB_PATH.stat().st_mtime_ns)
data["response"] = data["response"].fillna("not applicable")

st.sidebar.header("Cohort filters")
filters = {
    "project": multiselect_filter(data, "project", "Project"),
    "condition": multiselect_filter(data, "condition", "Condition"),
    "treatment": multiselect_filter(data, "treatment", "Treatment"),
    "response": multiselect_filter(data, "response", "Response"),
    "sample_type": multiselect_filter(data, "sample_type", "Sample type"),
    "time_from_treatment_start": multiselect_filter(
        data, "time_from_treatment_start", "Time from treatment start"
    ),
}

filtered = data.copy()
for column, values in filters.items():
    filtered = filtered[filtered[column].isin(values)]

tab_samples, tab_response, tab_baseline = st.tabs(
    ["Sample composition", "Responder analysis", "Baseline cohort"]
)

with tab_samples:
    st.caption("Part 2 — the sidebar filters apply to this tab.")
    if filtered.empty:
        st.warning("No samples match the selected filters.")
    else:
        summary = filtered.drop_duplicates("sample")
        subject_count = summary[["project", "subject"]].drop_duplicates().shape[0]
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Samples", f"{summary['sample'].nunique():,}")
        col2.metric("Subjects", f"{subject_count:,}")
        col3.metric("Projects", f"{summary['project'].nunique():,}")
        col4.metric("Median cells / sample", f"{summary['total_count'].median():,.0f}")

        st.subheader("Mean relative cell-type composition")
        composition = (
            filtered.groupby(["population", "cell_type"], as_index=False)["percentage"]
            .mean()
            .sort_values("percentage", ascending=False)
            .set_index("population")
        )
        st.bar_chart(composition[["percentage"]])

        st.subheader("Filtered sample-level percentages")
        wide = filtered.pivot_table(
            index=[
                "project",
                "subject",
                "condition",
                "sex",
                "treatment",
                "response",
                "sample",
                "sample_type",
                "time_from_treatment_start",
                "total_count",
            ],
            columns="population",
            values="percentage",
        ).reset_index()
        st.dataframe(wide, use_container_width=True, hide_index=True)

with tab_response:
    st.caption(
        "Part 3 — melanoma + miraclib + PBMC; subject-level averages, two-sided "
        "Mann-Whitney U tests, and Benjamini-Hochberg correction."
    )
    stats_path = OUTPUT_DIR / "response_statistics.csv"
    subject_path = OUTPUT_DIR / "response_subject_level_frequencies.csv"
    if not stats_path.is_file() or not subject_path.is_file():
        st.error("Analysis outputs are missing. Run `make pipeline` first.")
    else:
        statistics = pd.read_csv(stats_path)
        subject_level = pd.read_csv(subject_path)
        significance_mask = (
            statistics["significant_fdr_0_05"]
            .astype(str)
            .str.casefold()
            .eq("true")
        )
        significant = statistics.loc[significance_mask, "population"].tolist()
        if significant:
            st.success(
                "Populations significant at FDR < 0.05: " + ", ".join(significant)
            )
        else:
            st.info("No population is significant at FDR < 0.05.")
        st.dataframe(statistics, use_container_width=True, hide_index=True)

        labels = (
            subject_level[["population", "cell_type"]]
            .drop_duplicates()
            .set_index("population")["cell_type"]
            .to_dict()
        )
        selected = st.selectbox(
            "Population boxplot",
            sorted(subject_level["population"].unique()),
            format_func=lambda value: f"{labels[value]} ({value})",
        )
        plot_path = OUTPUT_DIR / "plots" / f"{selected}.png"
        if plot_path.is_file():
            st.image(str(plot_path), caption=f"Responder comparison: {selected}")
        else:
            st.warning(f"Plot not found: {plot_path.name}")

with tab_baseline:
    st.caption("Part 4 — melanoma + PBMC + baseline (time 0) + miraclib.")
    paths = {
        "matched": OUTPUT_DIR / "baseline_matched_samples.csv",
        "project": OUTPUT_DIR / "baseline_samples_by_project.csv",
        "response": OUTPUT_DIR / "baseline_subjects_by_response.csv",
        "sex": OUTPUT_DIR / "baseline_subjects_by_sex.csv",
    }
    if not all(path.is_file() for path in paths.values()):
        st.error("Baseline outputs are missing. Run `make pipeline` first.")
    else:
        matched = pd.read_csv(paths["matched"])
        by_project = pd.read_csv(paths["project"])
        by_response = pd.read_csv(paths["response"])
        by_sex = pd.read_csv(paths["sex"])
        st.metric("Matching samples", f"{len(matched):,}")
        left, middle, right = st.columns(3)
        with left:
            st.subheader("Samples by project")
            st.dataframe(by_project, use_container_width=True, hide_index=True)
        with middle:
            st.subheader("Subjects by response")
            st.dataframe(by_response, use_container_width=True, hide_index=True)
        with right:
            st.subheader("Subjects by sex")
            st.dataframe(by_sex, use_container_width=True, hide_index=True)
        st.subheader("Matching sample details")
        st.dataframe(matched, use_container_width=True, hide_index=True)
