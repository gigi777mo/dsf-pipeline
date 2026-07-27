#!/usr/bin/env python3
"""
Tm calculation for DSF / Protein Thermal Shift.

Methods (most published / Thermo PTS standard):
  1. Derivative Tm  — argmax of dF/dT
  2. Boltzmann Tm   — two-state sigmoidal midpoint fit

References: Niesen et al. Nat Protoc 2007; Pantoliano et al. J Biomol Screen 2001;
Thermo Fisher Protein Thermal Shift Software (Derivative Tm & Boltzmann Tm).
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, savgol_filter


def smooth_signal(y: np.ndarray, window_length: int = 11, polyorder: int = 3) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 5:
        return y
    wl = min(window_length, n if n % 2 == 1 else n - 1)
    if wl < polyorder + 2:
        wl = polyorder + 2 if (polyorder + 2) % 2 == 1 else polyorder + 3
    wl = min(wl, n if n % 2 == 1 else n - 1)
    if wl < 5:
        return y
    try:
        return savgol_filter(y, window_length=wl, polyorder=min(polyorder, wl - 1))
    except Exception:
        return uniform_filter1d(y, size=max(3, wl // 2))


def derivative_tm(
    temperature: np.ndarray,
    fluorescence: np.ndarray,
    multi_peak: bool = False,
    min_prominence: float = 0.05,
    smooth: bool = True,
    window_length: int = 11,
    polyorder: int = 3,
) -> dict:
    """
    First-derivative Tm: temperature at peak(s) of dF/dT.
    """
    T = np.asarray(temperature, dtype=float)
    F = np.asarray(fluorescence, dtype=float)
    mask = np.isfinite(T) & np.isfinite(F)
    T, F = T[mask], F[mask]
    if len(T) < 5:
        return {"Tm_D": np.nan, "peaks": [], "dF_dT": None, "T_d": None, "flag": "too_few_points"}

    order = np.argsort(T)
    T, F = T[order], F[order]
    if smooth:
        F_s = smooth_signal(F, window_length=window_length, polyorder=polyorder)
    else:
        F_s = F

    dF = np.gradient(F_s, T)
    # Prefer unfolding as positive peak (typical Sypro/PTS dye increase).
    # If overall trend is inverted, flip.
    if np.nanmax(dF) < -np.nanmin(dF):
        dF = -dF

    prom = min_prominence * (np.nanmax(np.abs(dF)) + 1e-12)
    peaks, props = find_peaks(dF, prominence=prom)
    if len(peaks) == 0:
        # fallback: global max
        idx = int(np.nanargmax(dF))
        peaks = np.array([idx])

    peak_T = T[peaks]
    peak_H = dF[peaks]
    order_p = np.argsort(peak_H)[::-1]
    peak_T = peak_T[order_p]
    peak_H = peak_H[order_p]

    if not multi_peak:
        peak_T = peak_T[:1]
        peak_H = peak_H[:1]

    return {
        "Tm_D": float(peak_T[0]),
        "peaks": [float(x) for x in peak_T],
        "peak_heights": [float(x) for x in peak_H],
        "dF_dT": dF,
        "T_d": T,
        "F_smooth": F_s,
        "flag": "multi_peak" if len(peak_T) > 1 else "ok",
    }


def _boltzmann(T, Fmin, Fmax, Tm, a):
    return Fmin + (Fmax - Fmin) / (1.0 + np.exp((Tm - T) / a))


def boltzmann_tm(
    temperature: np.ndarray,
    fluorescence: np.ndarray,
    roa_t_min: float | None = None,
    roa_t_max: float | None = None,
    center_hint: float | None = None,
    roa_half_width: float = 12.0,
    smooth: bool = True,
    window_length: int = 11,
    polyorder: int = 3,
) -> dict:
    """
    Boltzmann two-state Tm (midpoint of sigmoidal fit).
    """
    T = np.asarray(temperature, dtype=float)
    F = np.asarray(fluorescence, dtype=float)
    mask = np.isfinite(T) & np.isfinite(F)
    T, F = T[mask], F[mask]
    if len(T) < 8:
        return {"Tm_B": np.nan, "a": np.nan, "flag": "too_few_points"}

    order = np.argsort(T)
    T, F = T[order], F[order]
    if smooth:
        F_s = smooth_signal(F, window_length=window_length, polyorder=polyorder)
    else:
        F_s = F

    if roa_t_min is None or roa_t_max is None:
        if center_hint is None or not np.isfinite(center_hint):
            # use derivative peak as hint
            d = derivative_tm(T, F_s, smooth=False)
            center_hint = d["Tm_D"]
        if not np.isfinite(center_hint):
            center_hint = float(np.median(T))
        roa_t_min = float(center_hint - roa_half_width)
        roa_t_max = float(center_hint + roa_half_width)

    m = (T >= roa_t_min) & (T <= roa_t_max)
    if m.sum() < 6:
        m = np.ones(len(T), dtype=bool)

    Tt, Ft = T[m], F_s[m]
    Fmin0, Fmax0 = float(np.nanmin(Ft)), float(np.nanmax(Ft))
    if abs(Fmax0 - Fmin0) < 1e-9:
        return {"Tm_B": np.nan, "a": np.nan, "flag": "flat", "roa": (roa_t_min, roa_t_max)}

    Tm0 = float(center_hint) if center_hint is not None else float(np.median(Tt))
    a0 = max(1.0, (Tt.max() - Tt.min()) / 8.0)

    try:
        popt, _ = curve_fit(
            _boltzmann,
            Tt,
            Ft,
            p0=[Fmin0, Fmax0, Tm0, a0],
            bounds=(
                [-np.inf, -np.inf, Tt.min() - 5, 0.1],
                [np.inf, np.inf, Tt.max() + 5, 50.0],
            ),
            maxfev=10000,
        )
        Fmin, Fmax, Tm, a = [float(x) for x in popt]
        return {
            "Tm_B": Tm,
            "a": a,
            "Fmin": Fmin,
            "Fmax": Fmax,
            "roa": (float(roa_t_min), float(roa_t_max)),
            "flag": "ok",
            "fit_T": Tt,
            "fit_F": _boltzmann(Tt, *popt),
        }
    except Exception:
        return {
            "Tm_B": np.nan,
            "a": np.nan,
            "roa": (float(roa_t_min), float(roa_t_max)),
            "flag": "fit_fail",
        }


def analyze_well(
    temperature: np.ndarray,
    fluorescence: np.ndarray,
    method: str = "both",
    **kwargs,
) -> dict:
    out = {"amplitude": float(np.nanmax(fluorescence) - np.nanmin(fluorescence))}
    dres = derivative_tm(temperature, fluorescence, **{
        k: kwargs[k]
        for k in ("multi_peak", "min_prominence", "smooth", "window_length", "polyorder")
        if k in kwargs
    })
    out["Tm_D"] = dres["Tm_D"]
    out["derivative_peaks"] = dres.get("peaks", [])
    out["derivative_flag"] = dres.get("flag", "")

    if method in ("boltzmann", "both"):
        bres = boltzmann_tm(
            temperature,
            fluorescence,
            center_hint=dres.get("Tm_D"),
            **{
                k: kwargs[k]
                for k in (
                    "roa_t_min",
                    "roa_t_max",
                    "roa_half_width",
                    "smooth",
                    "window_length",
                    "polyorder",
                )
                if k in kwargs
            },
        )
        out["Tm_B"] = bres["Tm_B"]
        out["boltzmann_a"] = bres.get("a", np.nan)
        out["boltzmann_flag"] = bres.get("flag", "")
        out["roa"] = bres.get("roa")
    else:
        out["Tm_B"] = np.nan

    if method == "derivative":
        out["Tm_B"] = np.nan

    return out
