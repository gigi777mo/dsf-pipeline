# Tm calculation methods

This pipeline implements the two most published Tm estimators used in Thermo Protein Thermal Shift Software and DSF literature.

## 1. Derivative Tm (Tm_D)

1. Optionally smooth fluorescence $F(T)$ (Savitzky–Golay).
2. Compute numerical first derivative $dF/dT$.
3. Find local maxima of $dF/dT$ within the analysis window.
4. **Primary Tm_D** = temperature of the strongest peak (or peaks if multi-Tm mode).

**Best for:** multi-domain proteins, multi-transition curves, noisy data where a global sigmoid fit is unreliable.

## 2. Boltzmann Tm (Tm_B)

Fit within a region of analysis (ROA):

$$F(T) = F_{min} + \frac{F_{max} - F_{min}}{1 + e^{(T_m - T)/a}}$$

- Auto ROA: window around the derivative peak, or rising portion of the curve.
- Manual ROA: set `t_min` / `t_max` in config.

**Best for:** clean single transitions; reports a midpoint of the two-state model.

## ΔTm

$$\Delta T_m = T_m(\text{sample}) - T_m(\text{reference})$$

Positive ΔTm → stabilization (typical ligand/buffer hit criterion in screening).

## Flags (typical)

| Flag | Meaning |
|------|--------|
| low_signal | Amplitude too small |
| flat | No clear transition |
| multi_peak | >1 derivative peak |
| fit_fail | Boltzmann fit did not converge |
| omit | User-marked exclude |
