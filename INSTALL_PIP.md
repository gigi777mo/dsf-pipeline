# Install with Python + pip (from GitHub)

No Anaconda required.

## Recommended

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate

pip install -U pip
pip install git+https://github.com/gigi777mo/dsf-pipeline.git

# Still need the scripts tree for full CLI, so clone once:
git clone https://github.com/gigi777mo/dsf-pipeline.git
cd dsf-pipeline
pip install -e .

python scripts/run_dsf.py --melt data/example_melt_wide.csv --out results/
# or, after editable install from repo root:
# dsf-run --melt data/example_melt_wide.csv --out results/
```

## Minimal (ZIP / clone + requirements only)

```bash
git clone https://github.com/gigi777mo/dsf-pipeline.git
cd dsf-pipeline
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python scripts/run_dsf.py --melt your_melt.csv --out results/
```
