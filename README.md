# DSF / Protein Thermal Shift Pipeline

**Universal differential scanning fluorimetry (DSF)** analysis for protein thermal shift assays.

Works with melt-curve exports from **any qPCR platform** (Thermo QuantStudio / Protein Thermal Shift chemistry, Bio-Rad CFX, Roche LightCycler, ABI, generic CSV).

**Tm methods (most published)**
1. **First-derivative Tm** — peak of dF/dT  
2. **Boltzmann (two-state) Tm** — midpoint of sigmoidal fit  

---

## Start here

| Document | Purpose |
|----------|--------|
| **[INSTRUCTIONS.md](INSTRUCTIONS.md)** | **Full standalone user guide** (inputs, run, QC, troubleshooting) |
| [docs/example_curves.md](docs/example_curves.md) | Ideal curves, formulas, good vs bad shapes |
| [docs/CITATIONS.md](docs/CITATIONS.md) | Literature for assay and Tm calculations |
| [docs/input_formats.md](docs/input_formats.md) | CSV layouts by instrument brand |
| [docs/tm_methods.md](docs/tm_methods.md) | Algorithm notes |

```bash
conda env create -f environment.yml && conda activate dsf

python scripts/run_dsf.py \
  --melt data/example_melt_wide.csv \
  --samples data/samples.csv \
  --method both \
  --out results/

# Optional: example figures with formulas on the plots
python scripts/generate_example_curves.py --out docs/figures
```

---

## Formulas

$$
F(T) = F_{min} + \frac{F_{max}-F_{min}}{1+e^{(T_m-T)/a}}
\qquad
T_{m,D}=\arg\max_T\frac{dF}{dT}
\qquad
\Delta T_m = T_m^{\mathrm{sample}}-T_m^{\mathrm{ref}}
$$

---

## Citation

- Pantoliano et al., *J Biomol Screen* 2001  
- Niesen et al., *Nat Protoc* 2007  
- Thermo Fisher Protein Thermal Shift Software (Derivative + Boltzmann Tm)  

Full list and acknowledgment text: **[docs/CITATIONS.md](docs/CITATIONS.md)**  
How to run: **[INSTRUCTIONS.md](INSTRUCTIONS.md)**

---

## License

MIT
