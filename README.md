# DSF / Protein Thermal Shift Pipeline

---

> ## 🔴 USE MINICONDA
>
> **Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) before anything else.**  
> Then: `conda env create -f environment.yml` → `conda activate dsf`  
> This is the easiest, least error-prone path.

---

**Beginner path:** open **[START_HERE.md](START_HERE.md)** and follow the numbered steps.

Universal DSF analysis for protein thermal shift assays.  
Works with melt-curve CSVs from **any qPCR** (Thermo, Bio-Rad, Roche, ABI, …).

**Tm methods:** first-derivative peak + Boltzmann (two-state) fit — same ideas as Thermo Protein Thermal Shift software.

```bash
conda env create -f environment.yml
conda activate dsf
python scripts/run_dsf.py --melt data/example_melt_wide.csv --out results/
```

| Doc | Purpose |
|-----|--------|
| **[START_HERE.md](START_HERE.md)** | Simplest steps (start here) |
| [INSTRUCTIONS.md](INSTRUCTIONS.md) | Full user guide |
| [docs/example_curves.md](docs/example_curves.md) | What good curves look like |
| [docs/CITATIONS.md](docs/CITATIONS.md) | Papers to cite |
| [INSTALL_PIP.md](INSTALL_PIP.md) | Pip-only alternative (advanced) |

## License

MIT
