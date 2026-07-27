# DSF / Protein Thermal Shift Pipeline

**Universal differential scanning fluorimetry (DSF)** analysis for protein thermal shift assays.

Compatible with melt-curve exports from **any qPCR platform** (Thermo QuantStudio / Protein Thermal Shift chemistry, Bio-Rad CFX, Roche LightCycler, Applied Biosystems StepOne, Agilent, and generic CSV).

**Tm methods (most published approaches)**
1. **First-derivative Tm** — temperature at the peak of dF/dT (inflection of the melt)
2. **Boltzmann (two-state) Tm** — midpoint of a sigmoidal fit to the unfolding transition

These are the same two methods used in Thermo Fisher Protein Thermal Shift Software and throughout the DSF literature.

Full citations: **[docs/CITATIONS.md](docs/CITATIONS.md)**

---

## What this pipeline does

| Step | Description |
|------|-------------|
| **Import** | Wide or long CSV melt curves (Temperature + fluorescence per well) |
| **QC / smooth** | Optional Savitzky–Golay smoothing; flag low-signal or flat curves |
| **Tm (derivative)** | Numerical dF/dT; peak finding (single or multi-peak) |
| **Tm (Boltzmann)** | Two-state sigmoidal fit over an auto or user region of analysis (ROA) |
| **ΔTm** | vs reference condition / wells |
| **Replicates** | Mean, SD, SEM per sample group |
| **Plots** | Raw melt, derivative, overlay, ΔTm bar/dot plots |
| **Export** | Tm table (CSV), flags, figures |

---

## Input format (universal)

### Preferred: wide CSV

```text
Temperature,A1,A2,A3,...,H12
25.0,1200,1180,...
25.5,1210,1190,...
...
95.0,...
```

- First column: **Temperature** (°C)
- Other columns: well or sample IDs with fluorescence (RFU)

### Alternative: long CSV

```text
Temperature,Well,Fluorescence
25.0,A1,1200
25.0,A2,1180
...
```

### Sample sheet (optional)

```text
well,sample,condition,replicate,is_reference
A1,ProteinX,buffer,1,true
A2,ProteinX,buffer,2,true
A3,ProteinX,ligand,1,false
```

Export **raw fluorescence vs temperature** from your instrument software (not amplification Ct tables).

---

## Quick start

```bash
git clone https://github.com/gigi777mo/dsf-pipeline.git
cd dsf-pipeline

conda env create -f environment.yml
conda activate dsf

# Analyze a plate export
python scripts/run_dsf.py \
  --melt data/example_melt_wide.csv \
  --samples data/samples.csv \
  --method both \
  --out results/
```

Methods: `derivative` | `boltzmann` | `both` (default).

---

## Tm calculation (literature basis)

### 1. First-derivative Tm (Tm_D)

$$
T_{m,D} = \arg\max_T \left(\frac{dF}{dT}\right)
$$

- Widely used for single and **multiple** transitions
- Thermo Protein Thermal Shift Software: “Derivative Tm”
- Preferred when curves show more than one melt phase

### 2. Boltzmann two-state Tm (Tm_B)

$$
F(T) = F_{\min} + \frac{F_{\max} - F_{\min}}{1 + e^{(T_m - T)/a}}
$$

- $T_m$ = midpoint of the unfolding transition
- $a$ = slope / steepness parameter
- Standard in Thermo PTS software and most published ThermoFluor/DSF papers for **single-transition** curves

**Rule of thumb (Thermo + field practice):** use a **consistent method** for all samples in a study; for multi-peak curves use derivative Tm.

---

## Outputs

```
results/
├── tm_table.csv          # Tm_D, Tm_B, ΔTm, flags per well
├── replicate_summary.csv # mean ± SD by sample/condition
├── plots/
│   ├── melt_curves.png
│   ├── derivative.png
│   ├── delta_tm.png
│   └── overlays/
└── flags.csv
```

---

## Documentation

- **[Citations](docs/CITATIONS.md)** — primary literature for methods
- **[Input formats](docs/input_formats.md)** — Thermo, Bio-Rad, Roche, generic
- **[Tm methods](docs/tm_methods.md)** — derivative vs Boltzmann details

---

## Citation

If you use this pipeline, cite the underlying methods (see **docs/CITATIONS.md**):

- Protein Thermal Shift / DSF principle and Thermo workflow materials
- **Boltzmann Tm** and **derivative Tm** as implemented in standard PTS analysis
- Foundational ThermoFluor / DSF method papers (Pantoliano et al.; Niesen et al.)

---

## License

MIT
