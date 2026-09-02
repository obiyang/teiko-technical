from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from teiko_analysis.pipeline import (
    benjamini_hochberg,
    compute_response_statistics,
    load_frequencies,
)


def test_benjamini_hochberg_known_values() -> None:
    adjusted = benjamini_hochberg([0.01, 0.04, 0.03, 0.002])
    np.testing.assert_allclose(adjusted, [0.02, 0.04, 0.04, 0.008])


def test_real_database_frequencies_sum_to_100() -> None:
    db_path = ROOT / "cell_counts.db"
    if not db_path.exists():
        return
    with sqlite3.connect(db_path) as connection:
        frequencies = load_frequencies(connection)
    totals = frequencies.groupby("sample")["percentage"].sum()
    np.testing.assert_allclose(totals.to_numpy(), 100.0)


def test_response_analysis_aggregates_by_subject() -> None:
    rows = []
    for subject, response, values in [
        ("r1", "yes", [30.0, 40.0]),
        ("r2", "yes", [35.0, 45.0]),
        ("n1", "no", [10.0, 20.0]),
        ("n2", "no", [15.0, 25.0]),
    ]:
        for value in values:
            rows.append(
                {
                    "subject": subject,
                    "project": "prj1",
                    "response": response,
                    "population": "b_cell",
                    "cell_type": "B cell",
                    "percentage": value,
                    "condition": "melanoma",
                    "treatment": "miraclib",
                    "sample_type": "PBMC",
                }
            )
    statistics, subject_level = compute_response_statistics(pd.DataFrame(rows))
    assert len(subject_level) == 4
    assert statistics.loc[0, "n_responder_subjects"] == 2
    assert statistics.loc[0, "n_nonresponder_subjects"] == 2
    assert statistics.loc[0, "population"] == "b_cell"
    assert statistics.loc[0, "responder_mean_percentage"] == 37.5
    assert statistics.loc[0, "nonresponder_mean_percentage"] == 17.5
