# DSF / Protein Thermal Shift Pipeline — Instructions

Standalone user guide. Citations: [docs/CITATIONS.md](docs/CITATIONS.md).

---

## 1. What this pipeline is for

**DSF / protein thermal shift**: fluorescence vs temperature → **Tm** and **ΔTm**. Works with melt CSVs from any qPCR brand.

---

## 2. Install (Python + pip from GitHub — recommended)

```bash
git clone https://github.com/gigi777mo/dsf-pipeline.git
cd dsf-pipeline

python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate

pip install -U pip
pip install -r requirements.txt
# optional editable install:
pip install -e .
```

See also [INSTALL_PIP.md](INSTALL_PIP.md). Conda is optional (`environment.yml`).

---

## 3. Inputs

**Melt CSV (required)** — wide preferred:

```text
Temperature,A1,A2,...
25.0,1200,1180,...
```

**Sample sheet (optional)** for ΔTm / replicates:

```text
well,sample,condition,replicate,is_reference
A1,ProteinX,buffer,1,true
```

---

## 4. Run

```bash
python scripts/run_dsf.py \
  --melt path/to/melt.csv \
  --samples path/to/samples.csv \
  --method both \
  --out results/
```

Methods: `derivative` | `boltzmann` | `both`.

Outputs: `tm_table.csv`, plots, flags. Prefer `flags = ok`. Use one Tm method consistently.

---

## 5. Formulas

**Boltzmann:** $F(T)=F_{min}+(F_{max}-F_{min})/(1+e^{(T_m-T)/a})$  
**Derivative:** $T_{m,D}=\arg\max_T dF/dT$  
**Shift:** $\Delta T_m = T_m^{sample}-T_m^{ref}$

Example figures: `python scripts/generate_example_curves.py --out docs/figures`  
Guide: [docs/example_curves.md](docs/example_curves.md)

---

## 6. Citation

Pantoliano et al., *J Biomol Screen* 2001; Niesen et al., *Nat Protoc* 2007; Thermo Protein Thermal Shift Software.  
Full list: [docs/CITATIONS.md](docs/CITATIONS.md).
