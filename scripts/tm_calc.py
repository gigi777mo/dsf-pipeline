#!/usr/bin/env python3
"""
Tm calculation for DSF / Protein Thermal Shift (with defensive error handling).

Methods (most published / Thermo PTS standard):
  1. Derivative Tm  — argmax of dF/dT
  2. Boltzmann Tm   — two-state sigmoidal midpoint fit

References: Niesen et al. Nat Protoc 2007; Pantoliano et al. J Biomol Screen 2001;
Thermo Fisher Protein Thermal Shift Software (Derivative Tm & Boltzmann Tm).
"""

from __future__ import annotations

import logging
import warnings

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.optimize import OptimizeWarning, curve_fit
from scipy.signal import find_peaks, savgol_filter

logger = logging.getLogger(__name__)


def smooth_signal(y: np.ndarray, window_length: int = 11, polyorder: int = 3) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 5:
        return y
    try:
        wl = int(window_length)
        po = int(polyorder)
    except (TypeError, ValueError):
        return y
    if wl < 5:
        wl = 5
    if wl % 2 == 0:
        wl += 1
    wl = min(wl, n if n % 2 == 1 else n - 1)
    if wl < 5:
        return y
    po = min(po, wl - 1)
    if po < 1:
        po = 1
    try:
        return savgol_filter(y, window_length=wl, polyorder=po)
    except Exception as e:
        logger.debug("Savitzky–Golay failed (%s); using uniform filter", e)
        return uniform_filter1d(y, size=max(3, wl // 2))


def _finite_series(temperature, fluorescence):
    T = np.asarray(temperature, dtype=float)
    F = np.asarray(fluorescence, dtype=float)
    if T.shape != F.shape:
        raise ValueError(
            f"Temperature and fluorescence length mismatch: {T.shape} vs {F.shape}"
        )
    mask = np.isfinite(T) & np.isfinite(F)
    return T[mask], F[mask]


def derivative_tm(
    temperature: np.ndarray,
    fluorescence: np.ndarray,
    multi_peak: bool = False,
    min_prominence: float = 0.05,
    smooth: bool = True,
    window_length: int = 11,
    polyorder: int = 3,
) -> dict:
    """First-derivative Tm: temperature at peak(s) of dF/dT."""
    empty = {
        "Tm_D": np.nan,
        "peaks": [],
        "peak_heights": [],
        "dF_dT": None,
        "T_d": None,
        "F_smooth": None,
        "flag": "too_few_points",
    }
    try:
        T, F = _finite_series(temperature, fluorescence)
    except ValueError as e:
        logger.warning("derivative_tm input error: %s", e)
        return {**empty, "flag": "input_error"}

    if len(T) < 5:
        return empty

    order = np.argsort(T)
    T, F = T[order], F[order]
    # require some temperature span
    if float(np.nanmax(T) - np.nanmin(T)) < 1.0:
        return {**empty, "flag": "temp_range_too_small"}

    if smooth:
        F_s = smooth_signal(F, window_length=window_length, polyorder=polyorder)
    else:
        F_s = F

    try:
        dF = np.gradient(F_s, T)
    except Exception as e:
        logger.warning("np.gradient failed: %s", e)
        return {**empty, "flag": "derivative_fail"}

    if not np.any(np.isfinite(dF)):
        return {**empty, "flag": "derivative_fail"}

    # Prefer unfolding as positive peak (typical PTS dye increase)
    try:
        if np.nanmax(dF) < -np.nanmin(dF):
            dF = -dF
    except ValueError:
        return {**empty, "flag": "derivative_fail"}

    absmax = float(np.nanmax(np.abs(dF))) + 1e-12
    prom = max(0.0, float(min_prominence)) * absmax
    try:
        peaks, _props = find_peaks(dF, prominence=prom)
    except Exception as e:
        logger.debug("find_peaks failed: %s", e)
        peaks = np.array([], dtype=int)

    if len(peaks) == 0:
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
    # clip exp argument to avoid overflow
    z = np.clip((Tm - T) / np.maximum(a, 1e-6), -60, 60)
    return Fmin + (Fmax - Fmin) / (1.0 + np.exp(z))


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
    """Boltzmann two-state Tm (midpoint of sigmoidal fit)."""
    fail = {"Tm_B": np.nan, "a": np.nan, "flag": "too_few_points"}
    try:
        T, F = _finite_series(temperature, fluorescence)
    except ValueError as e:
        logger.warning("boltzmann_tm input error: %s", e)
        return {**fail, "flag": "input_error"}

    if len(T) < 8:
        return fail

    order = np.argsort(T)
    T, F = T[order], F[order]
    if smooth:
        F_s = smooth_signal(F, window_length=window_length, polyorder=polyorder)
    else:
        F_s = F

    try:
        half_w = float(roa_half_width)
        if half_w <= 0:
            half_w = 12.0
    except (TypeError, ValueError):
        half_w = 12.0

    if roa_t_min is None or roa_t_max is None:
        if center_hint is None or not np.isfinite(center_hint):
            d = derivative_tm(T, F_s, smooth=False)
            center_hint = d["Tm_D"]
        if not np.isfinite(center_hint):
            center_hint = float(np.median(T))
        roa_t_min = float(center_hint - half_w)
        roa_t_max = float(center_hint + half_w)

    if roa_t_min >= roa_t_max:
        logger.warning("Invalid ROA (%s, %s); using full range", roa_t_min, roa_t_max)
        roa_t_min, roa_t_max = float(T.min()), float(T.max())

    m = (T >= roa_t_min) & (T <= roa_t_max)
    if int(m.sum()) < 6:
        m = np.ones(len(T), dtype=bool)

    Tt, Ft = T[m], F_s[m]
    Fmin0, Fmax0 = float(np.nanmin(Ft)), float(np.nanmax(Ft))
    if abs(Fmax0 - Fmin0) < 1e-9:
        return {
            "Tm_B": np.nan,
            "a": np.nan,
            "flag": "flat",
            "roa": (float(roa_t_min), float(roa_t_max)),
        }

    Tm0 = float(center_hint) if center_hint is not None and np.isfinite(center_hint) else float(
        np.median(Tt)
    )
    a0 = max(1.0, float(Tt.max() - Tt.min()) / 8.0)

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            popt, pcov = curve_fit(
                _boltzmann,
                Tt,
                Ft,
                p0=[Fmin0, Fmax0, Tm0, a0],
                bounds=(
                    [-np.inf, -np.inf, float(Tt.min()) - 5, 0.1],
                    [np.inf, np.inf, float(Tt.max()) + 5, 50.0],
                ),
                maxfev=10000,
            )
        Fmin, Fmax, Tm, a = [float(x) for x in popt]
        # sanity: Tm inside extended data range
        if not (float(T.min()) - 10 <= Tm <= float(T.max()) + 10):
            return {
                "Tm_B": np.nan,
                "a": a,
                "roa": (float(roa_t_min), float(roa_t_max)),
                "flag": "fit_out_of_range",
            }
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
    except Exception as e:
        logger.debug("Boltzmann fit failed: %s", e)
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
    method = (method or "both").lower()
    if method not in {"derivative", "boltzmann", "both"}:
        logger.warning("Unknown method %r; using 'both'", method)
        method = "both"

    try:
        F = np.asarray(fluorescence, dtype=float)
        amp = float(np.nanmax(F) - np.nanmin(F)) if np.any(np.isfinite(F)) else np.nan
    except Exception:
        amp = np.nan

    out: dict = {"amplitude": amp}

    try:
        dres = derivative_tm(
            temperature,
            fluorescence,
            **{
                k: kwargs[k]
                for k in ("multi_peak", "min_prominence", "smooth", "window_length", "polyorder")
                if k in kwargs
            },
        )
    except Exception as e:
        logger.exception("derivative_tm crashed: %s", e)
        dres = {
            "Tm_D": np.nan,
            "peaks": [],
            "flag": "error",
        }

    out["Tm_D"] = dres.get("Tm_D", np.nan)
    out["derivative_peaks"] = dres.get("peaks", [])
    out["derivative_flag"] = dres.get("flag", "")

    if method in ("boltzmann", "both"):
        try:
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
        except Exception as e:
            logger.exception("boltzmann_tm crashed: %s", e)
            bres = {"Tm_B": np.nan, "a": np.nan, "flag": "error"}
        out["Tm_B"] = bres.get("Tm_B", np.nan)
        out["boltzmann_a"] = bres.get("a", np.nan)
        out["boltzmann_flag"] = bres.get("flag", "")
        out["roa"] = bres.get("roa")
    else:
        out["Tm_B"] = np.nan
        out["boltzmann_flag"] = "skipped"

    if method == "derivative":
        out["Tm_B"] = np.nan

    return out
