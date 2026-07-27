#!/usr/bin/env python3
"""Load universal melt-curve tables (wide or long CSV) with validation."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

TEMP_ALIASES = {
    "temperature",
    "temp",
    "temp.",
    "t",
    "temperature (c)",
    "temperature (°c)",
    "temperature(c)",
}


class MeltInputError(ValueError):
    """Invalid melt-curve or sample-sheet input."""


def _find_temp_column(columns) -> str:
    for c in columns:
        if str(c).strip().lower() in TEMP_ALIASES:
            return c
    logger.warning(
        "No column named Temperature/Temp found; using first column %r as temperature",
        columns[0],
    )
    return columns[0]


def load_melt_csv(path: str | Path) -> pd.DataFrame:
    """
    Return DataFrame indexed by temperature (°C) with one column per well/sample.
    Accepts wide (Temperature, A1, A2, ...) or long (Temperature, Well, Fluorescence).
    """
    path = Path(path)
    if not path.exists():
        raise MeltInputError(f"Melt-curve file not found: {path}")
    if path.stat().st_size == 0:
        raise MeltInputError(f"Melt-curve file is empty: {path}")

    try:
        df = pd.read_csv(path)
    except pd.errors.EmptyDataError as e:
        raise MeltInputError(f"Could not parse CSV (empty or invalid): {path}") from e
    except pd.errors.ParserError as e:
        raise MeltInputError(f"CSV parse error in {path}: {e}") from e
    except OSError as e:
        raise MeltInputError(f"Cannot read {path}: {e}") from e

    if df.empty:
        raise MeltInputError(f"Melt-curve table has no rows: {path}")
    if df.shape[1] < 2:
        raise MeltInputError(
            f"Expected at least Temperature + one well column; got columns: {list(df.columns)}"
        )

    df.columns = [str(c).strip() for c in df.columns]
    if df.columns.duplicated().any():
        dups = df.columns[df.columns.duplicated()].tolist()
        raise MeltInputError(f"Duplicate column names in melt file: {dups}")

    well_col = None
    fluo_col = None
    for c in df.columns:
        cl = c.lower()
        if cl in {"well", "sample", "well_id", "id"} and well_col is None:
            well_col = c
        if cl in {"fluorescence", "fluo", "rfu", "signal", "fluor"} and fluo_col is None:
            fluo_col = c
    temp_col = _find_temp_column(df.columns)

    try:
        if well_col and fluo_col:
            if df[fluo_col].isna().all():
                raise MeltInputError("Fluorescence column is entirely missing/NaN")
            wide = df.pivot_table(
                index=temp_col, columns=well_col, values=fluo_col, aggfunc="mean"
            )
            wide = wide.sort_index()
            wide.index = pd.to_numeric(wide.index, errors="coerce")
            if wide.index.isna().all():
                raise MeltInputError("Temperature values could not be parsed as numbers")
            wide = wide.loc[wide.index.notna()].sort_index()
            wide.columns = [str(c) for c in wide.columns]
        else:
            wide = df.set_index(temp_col)
            wide.index = pd.to_numeric(wide.index, errors="coerce")
            if wide.index.isna().all():
                raise MeltInputError("Temperature values could not be parsed as numbers")
            wide = wide.loc[wide.index.notna()].sort_index()
            for c in list(wide.columns):
                wide[c] = pd.to_numeric(wide[c], errors="coerce")
            wide = wide.dropna(axis=1, how="all")
    except MeltInputError:
        raise
    except Exception as e:
        raise MeltInputError(f"Failed to reshape melt data from {path}: {e}") from e

    if wide.empty or wide.shape[1] == 0:
        raise MeltInputError(
            f"No numeric well/fluorescence columns found in {path}. "
            "Export raw RFU vs temperature, not amplification Ct tables."
        )
    if wide.shape[0] < 5:
        raise MeltInputError(
            f"Too few temperature points ({wide.shape[0]}); need at least 5 for Tm analysis"
        )

    # Drop wells that are all NaN after numeric conversion
    all_nan = [c for c in wide.columns if wide[c].isna().all()]
    if all_nan:
        logger.warning("Dropping wells with no numeric data: %s", all_nan)
        wide = wide.drop(columns=all_nan)
    if wide.shape[1] == 0:
        raise MeltInputError("All wells were empty/non-numeric after parsing")

    n_nan = int(wide.isna().sum().sum())
    if n_nan:
        logger.warning("Melt table contains %d missing RFU values (will be skipped per well)", n_nan)

    logger.info("Loaded melt curves: %d temperatures × %d wells from %s", *wide.shape, path)
    return wide


def load_samples(path: str | Path | None) -> pd.DataFrame | None:
    if path is None:
        return None
    path = Path(path)
    if not path.exists():
        logger.warning("Sample sheet not found (continuing without it): %s", path)
        return None
    try:
        s = pd.read_csv(path)
    except Exception as e:
        raise MeltInputError(f"Cannot read sample sheet {path}: {e}") from e

    if s.empty:
        logger.warning("Sample sheet is empty; ignoring")
        return None

    s.columns = [str(c).strip().lower() for c in s.columns]
    if "well" not in s.columns:
        raise MeltInputError(
            f"Sample sheet must contain a 'well' column; found: {list(s.columns)}"
        )
    s["well"] = s["well"].astype(str).str.strip()
    if s["well"].duplicated().any():
        dups = s.loc[s["well"].duplicated(), "well"].tolist()
        raise MeltInputError(f"Duplicate wells in sample sheet: {dups}")

    if "is_reference" in s.columns:
        s["is_reference"] = (
            s["is_reference"]
            .astype(str)
            .str.lower()
            .isin(["1", "true", "yes", "y", "ref", "reference"])
        )
    else:
        s["is_reference"] = False

    return s
