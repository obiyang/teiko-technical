from __future__ import annotations

import csv
import sqlite3
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from load_data import REQUIRED_COLUMNS, build_database, read_and_validate

sys.path.insert(0, str(ROOT / "src"))
from teiko_analysis.pipeline import load_frequencies


def write_fixture(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=REQUIRED_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def sample_row(sample: str = "sample1") -> dict[str, str]:
    return {
        "project": "prj1",
        "subject": "subject1",
        "condition": "melanoma",
        "age": "50",
        "sex": "F",
        "treatment": "miraclib",
        "response": "yes",
        "sample": sample,
        "sample_type": "PBMC",
        "time_from_treatment_start": "0",
        "b_cell": "10",
        "cd8_t_cell": "20",
        "cd4_t_cell": "30",
        "nk_cell": "15",
        "monocyte": "25",
    }


def test_loader_creates_normalized_rows(tmp_path: Path) -> None:
    csv_path = tmp_path / "input.csv"
    db_path = tmp_path / "test.db"
    rows = [sample_row("sample1"), sample_row("sample2")]
    rows[1]["time_from_treatment_start"] = "7"
    write_fixture(csv_path, rows)

    assert build_database(csv_path, db_path) == 2

    with sqlite3.connect(db_path) as connection:
        assert connection.execute("SELECT COUNT(*) FROM subjects").fetchone()[0] == 1
        assert connection.execute("SELECT COUNT(*) FROM samples").fetchone()[0] == 2
        assert connection.execute("SELECT COUNT(*) FROM cell_types").fetchone()[0] == 5
        assert connection.execute("SELECT COUNT(*) FROM cell_counts").fetchone()[0] == 10
        assert connection.execute("PRAGMA integrity_check").fetchone()[0] == "ok"


def test_loader_rejects_duplicate_sample(tmp_path: Path) -> None:
    csv_path = tmp_path / "input.csv"
    write_fixture(csv_path, [sample_row(), sample_row()])
    with pytest.raises(ValueError, match="duplicate sample"):
        read_and_validate(csv_path)


def test_loader_rejects_negative_count(tmp_path: Path) -> None:
    csv_path = tmp_path / "input.csv"
    row = sample_row()
    row["b_cell"] = "-1"
    write_fixture(csv_path, [row])
    with pytest.raises(ValueError, match="must be non-negative"):
        read_and_validate(csv_path)


def test_zero_total_count_produces_missing_percentage(tmp_path: Path) -> None:
    csv_path = tmp_path / "input.csv"
    db_path = tmp_path / "test.db"
    row = sample_row()
    for population in ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]:
        row[population] = "0"
    write_fixture(csv_path, [row])
    build_database(csv_path, db_path)

    with sqlite3.connect(db_path) as connection:
        frequencies = load_frequencies(connection)
    assert (frequencies["total_count"] == 0).all()
    assert frequencies["percentage"].isna().all()
