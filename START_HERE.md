# START HERE (no bioinformatics experience needed)

---

> ## 🔴 USE MINICONDA
>
> **Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) first.**  
> Do not start with random pip installs if you want the path of least resistance.  
> One tool → one environment → fewer errors.

---

## What this does (in plain English)

You heat a protein with a fluorescent dye in a PCR machine.  
This program reads the fluorescence table and finds the **melting temperature (Tm)**  
and how much it **shifted (ΔTm)** when you add a ligand or change the buffer.

You only need:
1. Miniconda installed  
2. A CSV of **Temperature + fluorescence** from your machine  
3. Three commands below  

---

## Step 1 — Install Miniconda (one time)

1. Download Miniconda for your computer:  
   https://docs.conda.io/en/latest/miniconda.html  
2. Install it (click Next / allow defaults).  
3. Open **Anaconda Prompt** (Windows) or **Terminal** (Mac/Linux).

---

## Step 2 — Download this pipeline

```bash
git clone https://github.com/gigi777mo/dsf-pipeline.git
cd dsf-pipeline
```

No `git`? Go to the GitHub page → green **Code** button → **Download ZIP** → unzip → open that folder in the terminal.

---

## Step 3 — Create the environment (one time)

```bash
conda env create -f environment.yml
conda activate dsf
```

Wait until it finishes. If it asks to proceed, type `y` and Enter.

---

## Step 4 — Put your data in

1. Export **raw fluorescence vs temperature** from your qPCR software (not Ct numbers).  
2. Save as CSV like this:

```text
Temperature,A1,A2,B1
25.0,1200,1180,990
30.0,1210,1190,1005
...
```

3. Put the file somewhere easy, e.g. `data/my_melt.csv`.

Optional sample sheet (for ΔTm vs a reference):

```text
well,sample,condition,replicate,is_reference
A1,MyProtein,buffer,1,true
A2,MyProtein,buffer,2,true
B1,MyProtein,ligand,1,false
```

---

## Step 5 — Run (every time)

```bash
conda activate dsf

python scripts/run_dsf.py \
  --melt data/my_melt.csv \
  --samples data/samples.csv \
  --method both \
  --out results/
```

No sample sheet? Omit `--samples`.

---

## Step 6 — Open the results

| File | What it is |
|------|------------|
| `results/tm_table.csv` | Tm for each well + flags |
| `results/plots/melt_curves.png` | Your curves |
| `results/plots/derivative.png` | Peaks = Tm |
| `results/plots/delta_tm.png` | Shifts vs reference |

**Use wells where `flags` says `ok`.**  
Positive ΔTm usually means the condition made the protein more stable.

---

## If something breaks

| Message / problem | Fix |
|-------------------|-----|
| `conda: command not found` | Install Miniconda; open a **new** terminal |
| `file not found` | Check the path after `--melt` |
| Empty / weird Tm | Confirm export is **RFU vs temperature**, not Ct |
| Environment errors | Run `conda activate dsf` again |

More detail: [INSTRUCTIONS.md](INSTRUCTIONS.md) · curves: [docs/example_curves.md](docs/example_curves.md)
