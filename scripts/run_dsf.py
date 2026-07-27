#!/usr/bin/env python3
"""Run universal DSF / Protein Thermal Shift analysis on melt-curve CSV."""

from __future__ import annotations

import argparse
import logging
import sys
import traceback
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from io_melt import MeltInputError, load_melt_csv, load_samples
from tm_calc import analyze_well, derivative_tm

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
)
logger = logging.getLogger("dsf")


def _load_config(path: str) -> dict:
    p = Path(path)
    if not p.exists():
        logger.warning("Config not found (%s); using defaults", path)
        return {}
    try:
        with open(p) as f:
            cfg = yaml.safe_load(f)
        if cfg is None:
            return {}
        if not isinstance(cfg, dict):
            raise ValueError("Config root must be a mapping/dict")
        return cfg
    except yaml.YAMLError as e:
        logger.error("Invalid YAML in %s: %s", path, e)
        raise SystemExit(2) from e
    except OSError as e:
        logger.error("Cannot read config %s: %s", path, e)
        raise SystemExit(2) from e


def _safe_plot_melt(melt: pd.DataFrame, path: Path) -> None:
    try:
        fig, ax = plt.subplots(figsize=(8, 5))
        for well in melt.columns:
            y = melt[well]
            ax.plot(melt.index, y, lw=1, alpha=0.8, label=str(well))
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel("Fluorescence (RFU)")
        ax.set_title("DSF melt curves")
        if len(melt.columns) <= 12:
            ax.legend(fontsize=7, frameon=False)
        fig.tight_layout()
        fig.savefig(path, dpi=150)
        plt.close(fig)
    except Exception as e:
        logger.warning("Could not write melt plot %s: %s", path, e)
        plt.close("all")


def _safe_plot_derivative(melt: pd.DataFrame, path: Path, sm: dict) -> None:
    try:
        fig, ax = plt.subplots(figsize=(8, 5))
        plotted = 0
        for well in melt.columns:
            T = melt.index.values.astype(float)
            F = melt[well].values.astype(float)
            try:
                d = derivative_tm(
                    T,
                    F,
                    smooth=sm.get("enabled", True),
                    window_length=sm.get("window_length", 11),
                    polyorder=sm.get("polyorder", 3),
                )
            except Exception as e:
                logger.debug("Derivative plot skip %s: %s", well, e)
                continue
            if d.get("T_d") is None or d.get("dF_dT") is None:
                continue
            ax.plot(d["T_d"], d["dF_dT"], lw=1, alpha=0.8, label=str(well))
            if np.isfinite(d.get("Tm_D", np.nan)):
                ax.axvline(d["Tm_D"], color="k", ls=":", lw=0.5, alpha=0.3)
            plotted += 1
        if plotted == 0:
            ax.text(0.5, 0.5, "No derivative curves", ha="center", transform=ax.transAxes)
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel("dF/dT")
        ax.set_title("First derivative (Tm_D = peak)")
        if plotted <= 12:
            ax.legend(fontsize=7, frameon=False)
        fig.tight_layout()
        fig.savefig(path, dpi=150)
        plt.close(fig)
    except Exception as e:
        logger.warning("Could not write derivative plot %s: %s", path, e)
        plt.close("all")


def _safe_plot_delta(tm: pd.DataFrame, path: Path) -> None:
    if "deltaTm_D" not in tm.columns or not tm["deltaTm_D"].notna().any():
        return
    try:
        fig, ax = plt.subplots(figsize=(7, 4))
        label = tm.apply(
            lambda r: f"{r.get('sample', r['well'])}|{r.get('condition', '')}",
            axis=1,
        )
        vals = tm["deltaTm_D"].fillna(0).values
        ax.bar(range(len(tm)), vals, tick_label=list(label))
        ax.set_ylabel("ΔTm_D (°C)")
        ax.set_title("Thermal shift vs reference")
        plt.xticks(rotation=45, ha="right", fontsize=8)
        fig.tight_layout()
        fig.savefig(path, dpi=150)
        plt.close(fig)
    except Exception as e:
        logger.warning("Could not write delta-Tm plot %s: %s", path, e)
        plt.close("all")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Universal DSF thermal shift analysis")
    parser.add_argument("--melt", required=True, help="Wide or long melt-curve CSV")
    parser.add_argument("--samples", default=None, help="Optional sample sheet CSV")
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument(
        "--method", default=None, choices=["derivative", "boltzmann", "both"]
    )
    parser.add_argument("--out", default="results")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        cfg = _load_config(args.config)
    except SystemExit:
        return 2

    method = args.method or cfg.get("method", "both")
    out_dir = Path(args.out)
    try:
        plot_dir = out_dir / "plots"
        plot_dir.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        logger.error("Cannot create output directory %s: %s", out_dir, e)
        return 1

    try:
        melt = load_melt_csv(args.melt)
    except MeltInputError as e:
        logger.error("%s", e)
        return 1
    except Exception as e:
        logger.error("Unexpected error loading melt file: %s", e)
        logger.debug(traceback.format_exc())
        return 1

    try:
        samples = load_samples(args.samples or cfg.get("samples_file"))
    except MeltInputError as e:
        logger.error("%s", e)
        return 1

    # Temperature window
    try:
        t_min = cfg.get("t_min")
        t_max = cfg.get("t_max")
        if t_min is not None:
            melt = melt.loc[melt.index.values >= float(t_min)]
        if t_max is not None:
            melt = melt.loc[melt.index.values <= float(t_max)]
        if melt.empty:
            logger.error("No data left after applying t_min/t_max window")
            return 1
    except (TypeError, ValueError) as e:
        logger.error("Invalid t_min/t_max in config: %s", e)
        return 1

    # Warn if sample wells missing from melt
    if samples is not None:
        missing = sorted(set(samples["well"]) - set(map(str, melt.columns)))
        extra = sorted(set(map(str, melt.columns)) - set(samples["well"]))
        if missing:
            logger.warning("Sample sheet wells not in melt file: %s", missing)
        if extra:
            logger.info("Melt wells not listed in sample sheet: %s", extra)

    sm = cfg.get("smoothing", {}) or {}
    der = cfg.get("derivative", {}) or {}
    bol = cfg.get("boltzmann", {}) or {}

    rows = []
    n_err = 0
    for well in melt.columns:
        T = melt.index.values.astype(float)
        F = melt[well].values.astype(float)
        if not np.any(np.isfinite(F)):
            logger.warning("Well %s has no finite RFU values; skipping analysis", well)
            rows.append(
                {
                    "well": well,
                    "amplitude": np.nan,
                    "Tm_D": np.nan,
                    "Tm_B": np.nan,
                    "derivative_flag": "no_data",
                    "boltzmann_flag": "no_data",
                    "derivative_peaks": "",
                }
            )
            n_err += 1
            continue
        try:
            res = analyze_well(
                T,
                F,
                method=method,
                multi_peak=der.get("multi_peak", False),
                min_prominence=der.get("min_prominence", 0.05),
                smooth=sm.get("enabled", True),
                window_length=sm.get("window_length", 11),
                polyorder=sm.get("polyorder", 3),
                roa_half_width=bol.get("roa_half_width", 12.0),
                roa_t_min=bol.get("roa_t_min"),
                roa_t_max=bol.get("roa_t_max"),
            )
        except Exception as e:
            logger.exception("Analysis failed for well %s: %s", well, e)
            res = {
                "amplitude": np.nan,
                "Tm_D": np.nan,
                "Tm_B": np.nan,
                "derivative_flag": "error",
                "boltzmann_flag": "error",
                "derivative_peaks": [],
            }
            n_err += 1
        row = {
            "well": well,
            **{k: res[k] for k in res if k not in ("derivative_peaks", "roa")},
        }
        peaks = res.get("derivative_peaks", []) or []
        row["derivative_peaks"] = ";".join(str(x) for x in peaks)
        rows.append(row)

    if not rows:
        logger.error("No wells could be analyzed")
        return 1

    tm = pd.DataFrame(rows)

    if samples is not None:
        try:
            tm = tm.merge(samples, on="well", how="left")
        except Exception as e:
            logger.error("Failed merging sample sheet: %s", e)
            return 1
    else:
        tm["sample"] = tm["well"]
        tm["condition"] = tm["well"]
        tm["is_reference"] = False

    # Reference Tm for delta
    ref_mask = tm.get("is_reference", False)
    if not isinstance(ref_mask, pd.Series):
        ref_mask = pd.Series([False] * len(tm))
    ref_wells = tm.loc[ref_mask == True, "well"]  # noqa: E712
    cfg_refs = (cfg.get("reference") or {}).get("wells") or []
    if len(ref_wells) == 0 and cfg_refs:
        ref_wells = pd.Series([str(w) for w in cfg_refs])

    for col_tm, col_d in [("Tm_D", "deltaTm_D"), ("Tm_B", "deltaTm_B")]:
        if col_tm not in tm.columns:
            continue
        if len(ref_wells):
            ref_vals = tm.loc[tm["well"].isin(ref_wells), col_tm]
            ref_mean = float(ref_vals.mean()) if ref_vals.notna().any() else np.nan
            if not np.isfinite(ref_mean):
                logger.warning("Reference wells have no valid %s; ΔTm not computed", col_tm)
        else:
            ref_mean = np.nan
            logger.info("No reference wells set; %s will be NaN", col_d)
        tm[col_d] = tm[col_tm] - ref_mean if np.isfinite(ref_mean) else np.nan

    # QC flags
    rel = (cfg.get("qc") or {}).get("min_rel_amplitude", 0.05)
    try:
        rel = float(rel)
    except (TypeError, ValueError):
        rel = 0.05
    amps = tm["amplitude"].values if "amplitude" in tm.columns else np.array([])
    finite_amps = amps[np.isfinite(amps)]
    amp_thr = rel * float(np.nanmax(finite_amps)) if len(finite_amps) else 0.0

    flags = []
    for _, r in tm.iterrows():
        f = []
        amp = r.get("amplitude", np.nan)
        if not np.isfinite(amp) or amp < amp_thr:
            f.append("low_signal")
        if r.get("boltzmann_flag") in {"fit_fail", "fit_out_of_range", "error"}:
            f.append(str(r.get("boltzmann_flag")))
        if r.get("derivative_flag") in {"multi_peak", "error", "no_data"}:
            f.append(str(r.get("derivative_flag")))
        if r.get("derivative_flag") == "too_few_points":
            f.append("too_few_points")
        flags.append(";".join(f) if f else "ok")
    tm["flags"] = flags

    try:
        tm.to_csv(out_dir / "tm_table.csv", index=False)
    except OSError as e:
        logger.error("Cannot write tm_table.csv: %s", e)
        return 1

    group_cols = [c for c in ["sample", "condition"] if c in tm.columns]
    if group_cols:
        try:
            agg = {}
            for c in ["Tm_D", "Tm_B", "deltaTm_D", "deltaTm_B"]:
                if c in tm.columns:
                    agg[c] = ["mean", "std", "count"]
            if agg:
                summary = tm.groupby(group_cols, dropna=False).agg(agg)
                summary.columns = ["_".join(map(str, col)).strip("_") for col in summary.columns.values]
                summary.to_csv(out_dir / "replicate_summary.csv")
        except Exception as e:
            logger.warning("Replicate summary failed: %s", e)

    _safe_plot_melt(melt, plot_dir / "melt_curves.png")
    _safe_plot_derivative(melt, plot_dir / "derivative.png", sm)
    _safe_plot_delta(tm, plot_dir / "delta_tm.png")

    # flags report
    try:
        tm[["well", "flags"]].to_csv(out_dir / "flags.csv", index=False)
    except OSError:
        pass

    logger.info("Wrote %s", out_dir / "tm_table.csv")
    logger.info("Plots in %s", plot_dir)
    if n_err:
        logger.warning("%d well(s) had analysis problems (see flags)", n_err)

    cols = [c for c in ["well", "Tm_D", "Tm_B", "flags"] if c in tm.columns]
    print(tm[cols].to_string(index=False))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        logger.error("Interrupted")
        sys.exit(130)
    except Exception as e:
        logger.error("Unhandled error: %s", e)
        logger.debug(traceback.format_exc())
        sys.exit(1)
