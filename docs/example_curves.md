# Example DSF curves, formulas, and instructions

This page shows **how a good Protein Thermal Shift / DSF curve should look**, the **formulas** used in this pipeline, and how to interpret common problems.

Generate the figures locally:

```bash
conda activate dsf
python scripts/generate_example_curves.py --out docs/figures
```

---

## 1. Ideal single-transition melt curve

As the protein unfolds, the dye binds exposed hydrophobic regions and fluorescence rises. For a **single domain / two-state** unfold, the curve is sigmoidal.

### Boltzmann (two-state) formula — **Tm_B**

$$
F(T) = F_{\min} + \frac{F_{\max} - F_{\min}}{1 + e^{(T_m - T)/a}}
$$

| Symbol | Meaning |
|--------|--------|
| $F(T)$ | Fluorescence at temperature $T$ |
| $F_{\min}$, $F_{\max}$ | Baseline and upper fluorescence |
| $T_m$ | **Melting temperature** (midpoint of the transition) |
| $a$ | Steepness of the transition (larger $a$ → broader melt) |

This is the **Boltzmann Tm** used in Thermo Fisher Protein Thermal Shift Software and most published DSF analyses (Niesen et al., Nat Protoc 2007).

**What “good” looks like**

1. Relatively **flat** at low $T$ (folded / native)
2. **Steep sigmoidal rise** through the transition
3. Plateau (sometimes a mild drop at very high $T$ from aggregation or dye effects)
4. One clear transition (not many small wiggles)

---

## 2. First-derivative curve — **Tm_D**

### Derivative formula

$$
T_{m,D} = \arg\max_T \left(\frac{dF}{dT}\right)
$$

The **peak of dF/dT** is the inflection point of the melt curve — the **Derivative Tm** in Thermo PTS software.

| Situation | Recommendation |
|-----------|-----------------|
| Clean single peak | Tm_D ≈ Tm_B; either is fine if used consistently |
| Multiple peaks | Use **derivative** (multi-peak mode); do not force one Boltzmann fit |
| Noisy data | Smooth first (Savitzky–Golay), then take derivative |

---

## 3. Thermal shift (ΔTm)

$$
\Delta T_m = T_m(\mathrm{sample}) - T_m(\mathrm{reference})
$$

| Result | Typical interpretation |
|--------|------------------------|
| **ΔTm > 0** | Stabilization (ligand binding, better buffer, etc.) |
| **ΔTm ≈ 0** | No shift vs reference |
| **ΔTm < 0** | Destabilization |

Screening hits are often defined by a minimum positive ΔTm (e.g. +1 to +2 °C), with replicate consistency.

---

## 4. Visual guide: good vs problematic curves

| Panel | Shape | What to do in this pipeline |
|-------|--------|-----------------------------|
| **A. Good** | Flat → sharp sigmoid → plateau | Report **Tm_B** and/or **Tm_D**; use for ΔTm |
| **B. High initial RFU** | Starts high, weak rise | Check protein:dye ratio, buffer, aggregation; interpret with caution |
| **C. Multiple transitions** | Two or more steps | Enable multi-peak derivative; **do not** rely on single Boltzmann Tm |
| **D. Flat / no melt** | Almost no amplitude | Flag `low_signal` / `flat`; **do not report a Tm** |

---

## 5. Step-by-step instructions (your data)

1. **Export raw fluorescence vs temperature** from your qPCR (not Ct tables).
2. Save as wide CSV: `Temperature,A1,A2,...` (see [input_formats.md](input_formats.md)).
3. Optional: sample sheet with `well,sample,condition,replicate,is_reference`.
4. Run:

```bash
python scripts/run_dsf.py \
  --melt your_melt.csv \
  --samples your_samples.csv \
  --method both \
  --out results/
```

5. Open `results/plots/melt_curves.png` and `derivative.png`.
6. Compare shapes to the ideal panels above.
7. Use `results/tm_table.csv`:
   - Prefer wells with `flags = ok`
   - For multi-peak flags, re-run with `derivative.multi_peak: true` in config if needed
8. Report **one method consistently** (Tm_D or Tm_B) for all conditions in a study.

---

## 6. Controls to include (best practice)

| Control | Purpose |
|---------|--------|
| Protein + dye in reference buffer | Baseline Tm |
| Dye only (no protein) | Dye/buffer artifacts |
| Known stabilizer or lysozyme | Instrument / assay sanity |
| Replicates (n ≥ 3 when possible) | Trust ΔTm |

---

## 7. Citations for these formulas

- **Boltzmann / two-state Tm** — Thermo Protein Thermal Shift Software; standard in DSF literature  
- **Derivative Tm** — peak of dF/dT (Thermo PTS “Derivative Tm”)  
- **Assay & interpretation** — Niesen et al., *Nat Protoc* 2007; Pantoliano et al., *J Biomol Screen* 2001  

Full list: [CITATIONS.md](CITATIONS.md)
