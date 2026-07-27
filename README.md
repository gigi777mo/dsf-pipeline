# DSF / Protein Thermal Shift Pipeline

**Universal differential scanning fluorimetry (DSF)** analysis for protein thermal shift assays.

Works with melt-curve exports from **any qPCR platform** (Thermo QuantStudio / Protein Thermal Shift chemistry, Bio-Rad CFX, Roche LightCycler, ABI, generic CSV).

**Tm methods (most published)**
1. **First-derivative Tm** — peak of dF/dT  
2. **Boltzmann (two-state) Tm** — midpoint of sigmoidal fit  

Full citations: **[docs/CITATIONS.md](docs/CITATIONS.md)**  
**Example curves + formulas + instructions:** **[docs/example_curves.md](docs/example_curves.md)**

```bash
python scripts/generate_example_curves.py --out docs/figures
```

---

## Ideal curve (formulas)

**Boltzmann Tm**

$$
F(T) = F_{min} + \frac{F_{max} - F_{min}}{1 + e^{(T_m - T)/a}}
$$

**Derivative Tm**

$$
T_{m,D} = \arg\max_T \frac{dF}{dT}
$$

**Shift vs reference**

$$
\Delta T_m = T_m(\mathrm{sample}) - T_m(\mathrm{reference})
$$

Good data: flat baseline → sharp sigmoid → plateau. Multi-peak → use derivative. Flat → do not report Tm.

---

## Quick start

```bash
conda env create -f environment.yml && conda activate dsf

python scripts/run_dsf.py \
  --melt data/example_melt_wide.csv \
  --samples data/samples.csv \
  --method both \
  --out results/
```

Input: wide CSV `Temperature,A1,A2,...` (raw RFU). See [docs/input_formats.md](docs/input_formats.md).

---

## Documentation

- [Example curves, formulas, instructions](docs/example_curves.md)
- [Tm methods](docs/tm_methods.md)
- [Input formats](docs/input_formats.md)
- [Citations](docs/CITATIONS.md)

---

## License

MIT
