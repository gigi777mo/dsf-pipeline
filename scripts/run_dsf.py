#!/usr/bin/env python3
"""Run universal DSF / Protein Thermal Shift analysis on melt-curve CSV."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from io_melt import load_melt_csv, load_samples
from tm_calc import analyze_well, derivative_tm, smooth_signal


def main():
    parser = argparse.ArgumentParser(description="Universal DSF thermal shift analysis")
    parser.add_argument("--melt", required=True, help="Wide or long melt-curve CSV")
    parser.add_argument("--samples", default=None, help="Optional sample sheet CSV")
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--method", default=None, choices=["derivative", "boltzmann", "both"])
    parser.add_argument("--out", default="results")
    args = parser.parse_args()

    cfg = {}
    if Path(args.config).exists():
        with open(args.config) as f:
            cfg = yaml.safe_load(f) or {}

    method = args.method or cfg.get("method", "both")
    out_dir = Path(args.out)
    plot_dir = out_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    melt = load_melt_csv(args.melt)
    samples = load_samples(args.samples or cfg.get("samples_file"))

    t_min = cfg.get("t_min")
    t_max = cfg.get("t_max")
    T_all = melt.index.values.astype(float)
    if t_min is not None:
        melt = melt.loc[T_all >= t_min]
    if t_max is not None:
        melt = melt.loc[melt.index.values <= t_max]

    sm = cfg.get("smoothing", {})
    der = cfg.get("derivative", {})
    bol = cfg.get("boltzmann", {})

    rows = []
    for well in melt.columns:
        T = melt.index.values.astype(float)
        F = melt[well].values.astype(float)
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
        row = {"well": well, **{k: res[k] for k in res if k not in ("derivative_peaks", "roa")}}
        row["derivative_peaks"] = ";".join(str(x) for x in res.get("derivative_peaks", []))
        rows.append(row)

    tm = pd.DataFrame(rows)

    # Merge sample metadata
    if samples is not None:
        tm = tm.merge(samples, on="well", how="left")
    else:
        tm["sample"] = tm["well"]
        tm["condition"] = tm["well"]
        tm["is_reference"] = False

    # Reference Tm for delta
    ref_wells = tm.loc[tm.get("is_reference", False) == True, "well"]
    if len(ref_wells) == 0 and cfg.get("reference", {}).get("wells"):
        ref_wells = pd.Series(cfg["reference"]["wells"])

    for col_tm, col_d in [("Tm_D", "deltaTm_D"), ("Tm_B", "deltaTm_B")]:
        if col_tm not in tm.columns:
            continue
        if len(ref_wells):
            ref_mean = tm.loc[tm["well"].isin(ref_wells), col_tm].mean()
        else:
            ref_mean = np.nan
        tm[col_d] = tm[col_tm] - ref_mean if np.isfinite(ref_mean) else np.nan

    # QC flags
    rel = cfg.get("qc", {}).get("min_rel_amplitude", 0.05)
    amps = tm["amplitude"].values
    amp_thr = rel * np.nanmax(amps) if len(amps) else 0
    flags = []
    for _, r in tm.iterrows():
        f = []
        if r["amplitude"] < amp_thr:
            f.append("low_signal")
        if r.get("boltzmann_flag") == "fit_fail":
            f.append("fit_fail")
        if r.get("derivative_flag") == "multi_peak":
            f.append("multi_peak")
        flags.append(";".join(f) if f else "ok")
    tm["flags"] = flags

    out_dir.mkdir(parents=True, exist_ok=True)
    tm.to_csv(out_dir / "tm_table.csv", index=False)

    # Replicate summary
    group_cols = [c for c in ["sample", "condition"] if c in tm.columns]
    if group_cols:
        agg = {}
        for c in ["Tm_D", "Tm_B", "deltaTm_D", "deltaTm_B"]:
            if c in tm.columns:
                agg[c] = ["mean", "std", "count"]
        if agg:
            summary = tm.groupby(group_cols, dropna=False).agg(agg)
            summary.columns = ["_".join(col).strip("_") for col in summary.columns.values]
            summary.to_csv(out_dir / "replicate_summary.csv")

    # --- Plots ---
    fig, ax = plt.subplots(figsize=(8, 5))
    for well in melt.columns:
        ax.plot(melt.index, melt[well], lw=1, alpha=0.8, label=str(well))
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Fluorescence (RFU)")
    ax.set_title("DSF melt curves")
    if len(melt.columns) <= 12:
        ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(plot_dir / "melt_curves.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    for well in melt.columns:
        T = melt.index.values.astype(float)
        F = melt[well].values.astype(float)
        d = derivative_tm(
            T,
            F,
            smooth=sm.get("enabled", True),
            window_length=sm.get("window_length", 11),
            polyorder=sm.get("polyorder", 3),
        )
        if d["T_d"] is not None:
            ax.plot(d["T_d"], d["dF_dT"], lw=1, alpha=0.8, label=str(well))
            if np.isfinite(d["Tm_D"]):
                ax.axvline(d["Tm_D"], color="k", ls=":", lw=0.5, alpha=0.3)
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("dF/dT")
    ax.set_title("First derivative (Tm_D = peak)")
    if len(melt.columns) <= 12:
        ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(plot_dir / "derivative.png", dpi=150)
    plt.close(fig)

    # delta Tm plot
    if "deltaTm_D" in tm.columns and tm["deltaTm_D"].notna().any():
        fig, ax = plt.subplots(figsize=(7, 4))
        label = tm.apply(
            lambda r: f"{r.get('sample', r['well'])}|{r.get('condition', '')}", axis=1
        )
        ax.bar(range(len(tm)), tm["deltaTm_D"], tick_label=label)
        ax.set_ylabel("ΔTm_D (°C)")
        ax.set_title("Thermal shift vs reference")
        plt.xticks(rotation=45, ha="right", fontsize=8)
        fig.tight_layout()
        fig.savefig(plot_dir / "delta_tm.png", dpi=150)
        plt.close(fig)

    print(f"[+] Wrote {out_dir / 'tm_table.csv'}")
    print(f"[+] Plots in {plot_dir}")
    print(tm[["well", "Tm_D", "Tm_B", "flags"]].to_string(index=False))


if __name__ == "__main__":
    main()
