#!/usr/bin/env python3
"""
Generate example DSF curves with formulas annotated.

Ideal Boltzmann melt, first-derivative Tm, delta-Tm schematic,
and good vs problematic shapes for training/documentation.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def boltzmann(T, Fmin, Fmax, Tm, a):
    z = np.clip((Tm - T) / max(a, 1e-6), -60, 60)
    return Fmin + (Fmax - Fmin) / (1.0 + np.exp(z))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="docs/figures")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    T = np.linspace(25, 95, 400)
    Tm, a = 58.0, 2.2
    Fmin, Fmax = 1000.0, 3000.0
    F = boltzmann(T, Fmin, Fmax, Tm, a)
    # mild high-T drop
    F = F - 0.15 * (Fmax - Fmin) * np.clip((T - 75) / 20, 0, 1)
    dF = np.gradient(F, T)

    # --- Ideal + formula ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    ax.plot(T, F, color="#1f77b4", lw=2.5, label="Fluorescence F(T)")
    ax.axvline(Tm, color="#d62728", ls="--", lw=1.5, label=f"Tm = {Tm:.0f} °C")
    ax.annotate("native (folded)", xy=(32, Fmin + 80), fontsize=9)
    ax.annotate("unfolded", xy=(70, Fmax - 200), fontsize=9)
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Fluorescence (RFU)")
    ax.set_title("Ideal DSF melt curve (single transition)")
    ax.legend(frameon=False, loc="upper left")
    ax.set_xlim(25, 95)
    ax.text(
        0.98,
        0.05,
        r"$F(T)=F_{min}+\frac{F_{max}-F_{min}}{1+e^{(T_m-T)/a}}$"
        "\n\n"
        r"$T_m$: Boltzmann midpoint"
        "\n"
        r"$a$: steepness",
        transform=ax.transAxes,
        fontsize=10,
        va="bottom",
        ha="right",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#f7f7f7", edgecolor="#aaaaaa"),
    )

    ax = axes[1]
    ax.plot(T, dF, color="#2ca02c", lw=2.5, label="dF/dT")
    idx = int(np.argmax(dF))
    ax.axvline(T[idx], color="#d62728", ls="--", lw=1.5, label="Tm_D = peak")
    ax.scatter([T[idx]], [dF[idx]], color="#d62728", zorder=5, s=40)
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("dF/dT")
    ax.set_title("First derivative (Derivative Tm)")
    ax.legend(frameon=False, loc="upper left")
    ax.set_xlim(25, 95)
    ax.text(
        0.98,
        0.05,
        r"$T_{m,D}=\arg\max_T\,\frac{dF}{dT}$"
        "\n\nPeak of derivative =\ninflection of melt curve",
        transform=ax.transAxes,
        fontsize=10,
        va="bottom",
        ha="right",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#f7f7f7", edgecolor="#aaaaaa"),
    )

    fig.tight_layout()
    p1 = out / "ideal_dsf_curves.png"
    fig.savefig(p1, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"[+] {p1}")

    # --- Good vs bad ---
    F_hi = F + 800 - 0.5 * (T - 25)
    F2 = 1000 + 900 / (1 + np.exp((45 - T) / 2.0)) + 900 / (1 + np.exp((65 - T) / 2.5))
    F_flat = 1500 + 30 * np.sin((T - 25) / 10) + np.linspace(0, 40, len(T))

    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
    axes[0, 0].plot(T, F, color="#1f77b4", lw=2)
    axes[0, 0].set_title("A. Good — single sigmoidal melt")
    axes[0, 0].set_ylabel("RFU")
    axes[0, 0].text(
        0.05,
        0.9,
        "Use Tm_B or Tm_D",
        transform=axes[0, 0].transAxes,
        fontsize=9,
        bbox=dict(facecolor="#e8f5e9", edgecolor="#81c784"),
    )

    axes[0, 1].plot(T, F_hi, color="#ff7f0e", lw=2)
    axes[0, 1].set_title("B. High initial fluorescence")
    axes[0, 1].text(
        0.05,
        0.88,
        "Check dye/protein ratio;\npartial unfold or sticky dye",
        transform=axes[0, 1].transAxes,
        fontsize=8,
        bbox=dict(facecolor="#fff3e0", edgecolor="#ffb74d"),
    )

    axes[1, 0].plot(T, F2, color="#9467bd", lw=2)
    axes[1, 0].set_title("C. Multiple transitions")
    axes[1, 0].set_xlabel("Temperature (°C)")
    axes[1, 0].set_ylabel("RFU")
    axes[1, 0].text(
        0.05,
        0.88,
        "Use derivative multi-peak;\nnot single Boltzmann",
        transform=axes[1, 0].transAxes,
        fontsize=8,
        bbox=dict(facecolor="#f3e5f5", edgecolor="#ba68c8"),
    )

    axes[1, 1].plot(T, F_flat, color="#7f7f7f", lw=2)
    axes[1, 1].set_title("D. Flat / no clear melt")
    axes[1, 1].set_xlabel("Temperature (°C)")
    axes[1, 1].text(
        0.05,
        0.88,
        "Flag low_signal / flat;\ndo not report Tm",
        transform=axes[1, 1].transAxes,
        fontsize=8,
        bbox=dict(facecolor="#eeeeee", edgecolor="#9e9e9e"),
    )

    for ax in axes.ravel():
        ax.set_xlim(25, 95)
    fig.suptitle("How DSF curves should look (and common problems)", fontsize=13, y=1.01)
    fig.tight_layout()
    p2 = out / "good_vs_bad_dsf.png"
    fig.savefig(p2, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"[+] {p2}")

    # --- Delta Tm ---
    fig, ax = plt.subplots(figsize=(7, 4))
    for Tm_i, lab, col in [
        (55, "Reference (no ligand)", "#1f77b4"),
        (62, "Ligand-bound (+ΔTm)", "#d62728"),
    ]:
        Fi = boltzmann(T, Fmin, Fmax, Tm_i, a)
        ax.plot(T, Fi, lw=2.2, color=col, label=lab)
        ax.axvline(Tm_i, color=col, ls="--", lw=1.2)
    ax.annotate(
        "",
        xy=(62, 2200),
        xytext=(55, 2200),
        arrowprops=dict(arrowstyle="<->", color="black", lw=1.5),
    )
    ax.text(58.5, 2350, r"$\Delta T_m = T_m^{sample} - T_m^{ref}$", ha="center", fontsize=11)
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Fluorescence (RFU)")
    ax.set_title("Thermal shift (ΔTm) — stabilization by ligand/buffer")
    ax.legend(frameon=False)
    ax.set_xlim(25, 95)
    fig.tight_layout()
    p3 = out / "delta_tm_schematic.png"
    fig.savefig(p3, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"[+] {p3}")
    print("Done. See docs/example_curves.md for instructions.")


if __name__ == "__main__":
    main()
