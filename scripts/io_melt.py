#!/usr/bin/env python3
"""Load universal melt-curve tables (wide or long CSV)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

TEMP_ALIASES = {
    "temperature",
    "temp",
    "temp.",
    "t",
    "temperature (c)",
    "temperature (°c)",
    "temperature(c)",
}


def _find_temp_column(columns) -> str:
    for c in columns:
        if str(c).strip().lower() in TEMP_ALIASES:
            return c
    # fallback: first column
    return columns[0]


def load_melt_csv(path: str | Path) -> pd.DataFrame:
    """
    Return DataFrame indexed by temperature (°C) with one column per well/sample.
    Accepts wide (Temperature, A1, A2, ...) or long (Temperature, Well, Fluorescence).
    """
    path = Path(path)
    df = pd.read_csv(path)
    df.columns = [str(c).strip() for c in df.columns]
    cols_lower = {c: c.lower() for c in df.columns}

    # Long format?
    well_col = None
    fluo_col = None
    for c in df.columns:
        cl = c.lower()
        if cl in {"well", "sample", "well_id", "id"}:
            well_col = c
        if cl in {"fluorescence", "fluo", "rfu", "signal", "fluor"}:
            fluo_col = c
    temp_col = _find_temp_column(df.columns)

    if well_col and fluo_col:
        wide = df.pivot_table(
            index=temp_col, columns=well_col, values=fluo_col, aggfunc="mean"
        )
        wide = wide.sort_index()
        wide.index = wide.index.astype(float)
        wide.columns = [str(c) for c in wide.columns]
        return wide

    # Wide format
    tcol = temp_col
    wide = df.set_index(tcol)
    wide.index = wide.index.astype(float)
    wide = wide.sort_index()
    # keep numeric columns only
    for c in list(wide.columns):
        wide[c] = pd.to_numeric(wide[c], errors="coerce")
    wide = wide.dropna(axis=1, how="all")
    return wide


def load_samples(path: str | Path | None) -> pd.DataFrame | None:
    if path is None or not Path(path).exists():
        return None
    s = pd.read_csv(path)
    s.columns = [str(c).strip().lower() for c in s.columns]
    if "well" not in s.columns:
        raise ValueError("samples CSV must contain a 'well' column")
    s["well"] = s["well"].astype(str)
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
