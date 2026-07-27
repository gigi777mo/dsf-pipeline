# DSF / Protein Thermal Shift Pipeline — Instructions

Standalone user guide for running and interpreting the **dsf-pipeline**.
Formulas and literature are cited in [docs/CITATIONS.md](docs/CITATIONS.md).

---

## 1. What this pipeline is for

Differential scanning fluorimetry (**DSF**), also called a **protein thermal shift** assay: protein + environment-sensitive dye is heated in a qPCR machine; fluorescence vs temperature is used to estimate the melting temperature (**Tm**). Shifts in Tm (**ΔTm**) report relative stability (ligands, buffers, mutations).

This software analyzes **exported melt curves from any qPCR brand** (Thermo QuantStudio / Protein Thermal Shift chemistry, Bio-Rad, Roche, ABI, generic CSV). It does **not** control the instrument.

---

## 2. Required inputs

### 2.1 Melt-curve file (required)

**Wide CSV (preferred)**

```text
Temperature,A1,A2,A3,...
25.0,1200,1180,...
25.5,1210,1190,...
```

**Long CSV (also accepted)**

```text
Temperature,Well,Fluorescence
25.0,A1,1200
```

**Rules**

- Export **raw fluorescence (RFU) vs temperature**, not amplification Ct tables.
- At least **5 temperature points** (typically many more).
- Temperature in °C; well columns must be numeric RFU.

More detail: [docs/input_formats.md](docs/input_formats.md).

### 2.2 Sample sheet (optional but recommended)

```text
well,sample,condition,replicate,is_reference
A1,ProteinX,buffer,1,true
A2,ProteinX,buffer,2,true
B1,ProteinX,ligand,1,false
```

- `is_reference=true` defines the baseline for **ΔTm**.
- Replicates share the same `sample` + `condition` for summary statistics.

### 2.3 Config (optional)

Edit [config/config.yaml](config/config.yaml) for smoothing, multi-peak derivative, Boltzmann ROA, QC thresholds, and fixed temperature windows.

---

## 3. Installation

```bash
git clone https://github.com/gigi777mo/dsf-pipeline.git
cd dsf-pipeline
conda env create -f environment.yml
conda activate dsf
```

---

## 4. How to run

```bash
python scripts/run_dsf.py \
  --melt data/example_melt_wide.csv \
  --samples data/samples.csv \
  --method both \
  --out results/
```

| Argument | Meaning |
|----------|--------|
| `--melt` | Path to melt-curve CSV |
| `--samples` | Optional sample sheet |
| `--method` | `derivative` \| `boltzmann` \| `both` (default) |
| `--config` | YAML config (default `config/config.yaml`) |
| `--out` | Output directory |
| `-v` | Verbose logging |

**Exit codes:** `0` success · `1` data/analysis error · `2` config error · `130` interrupted.

---

## 5. Tm methods (use one consistently per study)

### 5.1 Derivative Tm (Tm_D)

$$
T_{m,D} = \arg\max_T \left(\frac{dF}{dT}\right)
$$

- Temperature at the **peak of the first derivative** of fluorescence.
- Standard **Derivative Tm** in Thermo Fisher Protein Thermal Shift Software.
- Prefer for **multi-transition** curves; enable multi-peak in config if needed.

### 5.2 Boltzmann Tm (Tm_B)

$$
F(T) = F_{\min} + \frac{F_{\max} - F_{\min}}{1 + e^{(T_m - T)/a}}
$$

- $T_m$ = midpoint of a two-state sigmoidal fit (region of analysis around the transition).
- Standard **Boltzmann Tm** in Thermo PTS software and most published DSF papers.
- Prefer for **clean single** transitions.

### 5.3 Thermal shift

$$
\Delta T_m = T_m(\mathrm{sample}) - T_m(\mathrm{reference})
$$

Positive ΔTm usually means relative **stabilization** vs the reference wells.

**Citations for these calculations:** [docs/CITATIONS.md](docs/CITATIONS.md) (Pantoliano et al. 2001; Niesen et al. 2007; Thermo PTS User Guide / Software).

---

## 6. How a good curve should look

Generate annotated example figures:

```bash
python scripts/generate_example_curves.py --out docs/figures
```

Full visual guide: [docs/example_curves.md](docs/example_curves.md).

| Quality | Shape | Action |
|---------|--------|--------|
| **Good** | Flat baseline → steep sigmoid → plateau | Report Tm_B and/or Tm_D |
| **High start RFU** | High fluorescence before the melt | Check dye/protein ratio; interpret cautiously |
| **Multiple steps** | Two or more transitions | Derivative multi-peak; avoid single Boltzmann |
| **Flat** | Almost no amplitude | Flag `low_signal` / `flat`; **do not report Tm** |

---

## 7. Outputs

```text
results/
├── tm_table.csv             # per-well Tm_D, Tm_B, ΔTm, flags
├── replicate_summary.csv    # mean ± SD by sample/condition
├── flags.csv
└── plots/
    ├── melt_curves.png
    ├── derivative.png
    └── delta_tm.png
```

**Flags (examples):** `ok` · `low_signal` · `multi_peak` · `fit_fail` · `fit_out_of_range` · `no_data` · `too_few_points`.

Use wells with `flags = ok` for primary conclusions.

---

## 8. Recommended experimental controls

| Control | Purpose |
|---------|--------|
| Protein + dye in reference buffer | Baseline Tm |
| Dye only (no protein) | Dye/buffer artifacts |
| Known stable protein (e.g. lysozyme) | Assay/instrument check |
| Replicates (ideally n ≥ 3) | Reliable ΔTm |

---

## 9. Troubleshooting

| Problem | What to try |
|---------|-------------|
| “Melt-curve file not found” | Check path to `--melt` |
| “No numeric well columns” | Ensure export is RFU vs T, not Ct |
| “Too few temperature points” | Export full melt table |
| All Tm NaN | Inspect raw curves; widen ROA or disable bad t_min/t_max |
| Boltzmann `fit_fail` | Rely on Tm_D or adjust ROA in config |
| No ΔTm | Set `is_reference=true` on reference wells |
| Sample wells missing | Well IDs in sample sheet must match melt column names |

---

## 10. Citation

If you use this pipeline, cite the underlying methods:

- Pantoliano MW, et al. *J Biomol Screen.* 2001;6:429–440. (ThermoFluor / thermal shift screening)
- Niesen FH, Berglund H, Vedadi M. *Nat Protoc.* 2007;2:2212–2221. (widely used DSF protocol)
- Thermo Fisher Scientific. Protein Thermal Shift Software / Studies User Guide (Derivative Tm and Boltzmann Tm).

Full reference list and suggested acknowledgment: **[docs/CITATIONS.md](docs/CITATIONS.md)**.

Optional:

> Analysis used the open DSF pipeline (https://github.com/gigi777mo/dsf-pipeline).

---

## 11. Related docs in this repo

| File | Content |
|------|--------|
| [README.md](README.md) | Overview and quick start |
| [docs/example_curves.md](docs/example_curves.md) | Ideal curves, formulas, good vs bad |
| [docs/tm_methods.md](docs/tm_methods.md) | Tm algorithm notes |
| [docs/input_formats.md](docs/input_formats.md) | Platform export tips |
| [docs/CITATIONS.md](docs/CITATIONS.md) | Literature |
| [config/config.yaml](config/config.yaml) | Analysis parameters |
