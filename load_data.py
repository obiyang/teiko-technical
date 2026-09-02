"""Load the bundled cell-count CSV into a normalized SQLite database.

Run with no arguments from anywhere:

    python load_data.py

The loader always resolves paths relative to this file and atomically rebuilds
``cell_counts.db`` in the project root.
"""

from __future__ import annotations

import csv
import hashlib
import os
import sqlite3
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "data" / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"

METADATA_COLUMNS = [
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
]
CELL_COLUMNS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
REQUIRED_COLUMNS = METADATA_COLUMNS + CELL_COLUMNS


SCHEMA = """
PRAGMA foreign_keys = ON;

CREATE TABLE projects (
    project_id INTEGER PRIMARY KEY,
    name TEXT NOT NULL UNIQUE
);

CREATE TABLE conditions (
    condition_id INTEGER PRIMARY KEY,
    name TEXT NOT NULL UNIQUE
);

CREATE TABLE sexes (
    sex_id INTEGER PRIMARY KEY,
    label TEXT NOT NULL UNIQUE
);

CREATE TABLE treatments (
    treatment_id INTEGER PRIMARY KEY,
    name TEXT NOT NULL UNIQUE
);

CREATE TABLE responses (
    response_id INTEGER PRIMARY KEY,
    label TEXT NOT NULL UNIQUE
);

CREATE TABLE sample_types (
    sample_type_id INTEGER PRIMARY KEY,
    name TEXT NOT NULL UNIQUE
);

CREATE TABLE cell_types (
    cell_type_id INTEGER PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,
    source_column TEXT NOT NULL UNIQUE
);

CREATE TABLE subjects (
    subject_id INTEGER PRIMARY KEY,
    subject_code TEXT NOT NULL,
    project_id INTEGER NOT NULL REFERENCES projects(project_id),
    condition_id INTEGER NOT NULL REFERENCES conditions(condition_id),
    age INTEGER NOT NULL CHECK (age >= 0),
    sex_id INTEGER NOT NULL REFERENCES sexes(sex_id),
    UNIQUE (project_id, subject_code)
);

CREATE TABLE samples (
    sample_id INTEGER PRIMARY KEY,
    sample_code TEXT NOT NULL UNIQUE,
    subject_id INTEGER NOT NULL REFERENCES subjects(subject_id),
    treatment_id INTEGER NOT NULL REFERENCES treatments(treatment_id),
    response_id INTEGER REFERENCES responses(response_id),
    sample_type_id INTEGER NOT NULL REFERENCES sample_types(sample_type_id),
    time_from_treatment_start REAL NOT NULL
);

CREATE TABLE cell_counts (
    sample_id INTEGER NOT NULL REFERENCES samples(sample_id) ON DELETE CASCADE,
    cell_type_id INTEGER NOT NULL REFERENCES cell_types(cell_type_id),
    count INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, cell_type_id)
);

CREATE TABLE dataset_loads (
    load_id INTEGER PRIMARY KEY,
    source_file TEXT NOT NULL,
    source_sha256 TEXT NOT NULL,
    source_row_count INTEGER NOT NULL,
    loaded_at_utc TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_subjects_project_condition
    ON subjects(project_id, condition_id);
CREATE INDEX idx_samples_filters
    ON samples(treatment_id, sample_type_id, time_from_treatment_start, response_id);
CREATE INDEX idx_counts_cell_type
    ON cell_counts(cell_type_id);
"""


def _clean_text(value: str | None, column: str, row_number: int) -> str:
    cleaned = (value or "").strip()
    if not cleaned:
        raise ValueError(f"Row {row_number}: {column} must not be blank")
    return cleaned


def _nonnegative_int(value: str | None, column: str, row_number: int) -> int:
    try:
        parsed = int(_clean_text(value, column, row_number))
    except ValueError as exc:
        raise ValueError(f"Row {row_number}: {column} must be an integer") from exc
    if parsed < 0:
        raise ValueError(f"Row {row_number}: {column} must be non-negative")
    return parsed


def _number(value: str | None, column: str, row_number: int) -> float:
    try:
        return float(_clean_text(value, column, row_number))
    except ValueError as exc:
        raise ValueError(f"Row {row_number}: {column} must be numeric") from exc


def read_and_validate(csv_path: Path) -> list[dict[str, Any]]:
    """Read a source CSV, validate its contract, and return typed records."""
    if not csv_path.is_file():
        raise FileNotFoundError(f"Source data not found: {csv_path}")

    records: list[dict[str, Any]] = []
    sample_codes: set[str] = set()
    subject_attributes: dict[tuple[str, str], tuple[str, int, str]] = {}

    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        actual = reader.fieldnames or []
        missing = [column for column in REQUIRED_COLUMNS if column not in actual]
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(missing)}")

        for row_number, row in enumerate(reader, start=2):
            record: dict[str, Any] = {
                column: _clean_text(row.get(column), column, row_number)
                for column in METADATA_COLUMNS
                if column != "response"
            }
            record["response"] = (row.get("response") or "").strip() or None
            record["age"] = _nonnegative_int(row.get("age"), "age", row_number)
            record["time_from_treatment_start"] = _number(
                row.get("time_from_treatment_start"),
                "time_from_treatment_start",
                row_number,
            )
            for cell_column in CELL_COLUMNS:
                record[cell_column] = _nonnegative_int(
                    row.get(cell_column), cell_column, row_number
                )

            sample_code = record["sample"]
            if sample_code in sample_codes:
                raise ValueError(f"Row {row_number}: duplicate sample {sample_code!r}")
            sample_codes.add(sample_code)

            subject_key = (record["project"], record["subject"])
            attributes = (record["condition"], record["age"], record["sex"])
            previous = subject_attributes.setdefault(subject_key, attributes)
            if previous != attributes:
                raise ValueError(
                    f"Row {row_number}: inconsistent attributes for subject {subject_key}"
                )
            records.append(record)

    if not records:
        raise ValueError("Source CSV has no data rows")
    return records


def _insert_dimension(
    connection: sqlite3.Connection, table: str, value_column: str, values: set[str]
) -> dict[str, int]:
    connection.executemany(
        f"INSERT INTO {table} ({value_column}) VALUES (?)",
        [(value,) for value in sorted(values)],
    )
    id_column = {
        "projects": "project_id",
        "conditions": "condition_id",
        "sexes": "sex_id",
        "treatments": "treatment_id",
        "responses": "response_id",
        "sample_types": "sample_type_id",
    }[table]
    return {
        value: identifier
        for identifier, value in connection.execute(
            f"SELECT {id_column}, {value_column} FROM {table}"
        )
    }


def populate_database(
    connection: sqlite3.Connection, records: list[dict[str, Any]], source_path: Path
) -> None:
    """Create the schema and insert validated records in one transaction."""
    connection.executescript(SCHEMA)

    project_ids = _insert_dimension(
        connection, "projects", "name", {r["project"] for r in records}
    )
    condition_ids = _insert_dimension(
        connection, "conditions", "name", {r["condition"] for r in records}
    )
    sex_ids = _insert_dimension(
        connection, "sexes", "label", {r["sex"] for r in records}
    )
    treatment_ids = _insert_dimension(
        connection, "treatments", "name", {r["treatment"] for r in records}
    )
    response_ids = _insert_dimension(
        connection,
        "responses",
        "label",
        {r["response"] for r in records if r["response"] is not None},
    )
    sample_type_ids = _insert_dimension(
        connection, "sample_types", "name", {r["sample_type"] for r in records}
    )

    pretty_cell_names = {
        "b_cell": "B cell",
        "cd8_t_cell": "CD8 T cell",
        "cd4_t_cell": "CD4 T cell",
        "nk_cell": "NK cell",
        "monocyte": "Monocyte",
    }
    connection.executemany(
        "INSERT INTO cell_types (name, source_column) VALUES (?, ?)",
        [(pretty_cell_names[column], column) for column in CELL_COLUMNS],
    )
    cell_type_ids = {
        source_column: cell_type_id
        for cell_type_id, source_column in connection.execute(
            "SELECT cell_type_id, source_column FROM cell_types"
        )
    }

    subject_rows: dict[tuple[str, str], tuple[Any, ...]] = {}
    for record in records:
        key = (record["project"], record["subject"])
        subject_rows[key] = (
            record["subject"],
            project_ids[record["project"]],
            condition_ids[record["condition"]],
            record["age"],
            sex_ids[record["sex"]],
        )
    connection.executemany(
        """INSERT INTO subjects
           (subject_code, project_id, condition_id, age, sex_id)
           VALUES (?, ?, ?, ?, ?)""",
        subject_rows.values(),
    )
    subject_ids = {
        (project_name, subject_code): subject_id
        for subject_id, project_name, subject_code in connection.execute(
            """SELECT s.subject_id, p.name, s.subject_code
               FROM subjects s JOIN projects p USING (project_id)"""
        )
    }

    sample_rows = [
        (
            record["sample"],
            subject_ids[(record["project"], record["subject"])],
            treatment_ids[record["treatment"]],
            response_ids.get(record["response"]),
            sample_type_ids[record["sample_type"]],
            record["time_from_treatment_start"],
        )
        for record in records
    ]
    connection.executemany(
        """INSERT INTO samples
           (sample_code, subject_id, treatment_id, response_id,
            sample_type_id, time_from_treatment_start)
           VALUES (?, ?, ?, ?, ?, ?)""",
        sample_rows,
    )
    sample_ids = {
        sample_code: sample_id
        for sample_id, sample_code in connection.execute(
            "SELECT sample_id, sample_code FROM samples"
        )
    }

    count_rows = [
        (sample_ids[record["sample"]], cell_type_ids[cell_column], record[cell_column])
        for record in records
        for cell_column in CELL_COLUMNS
    ]
    connection.executemany(
        "INSERT INTO cell_counts (sample_id, cell_type_id, count) VALUES (?, ?, ?)",
        count_rows,
    )

    source_hash = hashlib.sha256(source_path.read_bytes()).hexdigest()
    try:
        displayed_source = str(source_path.relative_to(ROOT))
    except ValueError:
        # Tests and library callers may load an external fixture. Store only the
        # basename so provenance remains useful without persisting local paths.
        displayed_source = source_path.name
    connection.execute(
        """INSERT INTO dataset_loads
           (source_file, source_sha256, source_row_count)
           VALUES (?, ?, ?)""",
        (displayed_source, source_hash, len(records)),
    )

    integrity = connection.execute("PRAGMA integrity_check").fetchone()[0]
    if integrity != "ok":
        raise RuntimeError(f"SQLite integrity check failed: {integrity}")


def build_database(csv_path: Path = CSV_PATH, db_path: Path = DB_PATH) -> int:
    """Atomically rebuild ``db_path`` and return the number of loaded samples."""
    records = read_and_validate(csv_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = db_path.with_suffix(db_path.suffix + ".tmp")
    temporary_path.unlink(missing_ok=True)

    try:
        with sqlite3.connect(temporary_path) as connection:
            connection.execute("PRAGMA foreign_keys = ON")
            populate_database(connection, records, csv_path)
        os.replace(temporary_path, db_path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise
    return len(records)


def main() -> None:
    loaded = build_database()
    print(f"Loaded {loaded:,} samples into {DB_PATH.name}")


if __name__ == "__main__":
    main()
