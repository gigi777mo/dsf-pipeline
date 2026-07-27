# Citations — DSF / Protein Thermal Shift Pipeline

Cite the methods appropriate to your use. Core references for the **assay** and **Tm calculations** implemented here:

---

## Assay principle (ThermoFluor / DSF)

- **Pantoliano et al. (foundational thermal shift / ThermoFluor screening)**  
  Pantoliano MW, Petrella EC, Kwasnoski JD, et al.  
  *High-density miniaturized thermal shift assays as a general strategy for drug discovery.*  
  Journal of Biomolecular Screening. 2001;6(6):429–440.  
  https://doi.org/10.1177/108705710100600609

- **Niesen et al. (most-cited practical DSF protocol)**  
  Niesen FH, Berglund H, Vedadi M.  
  *The use of differential scanning fluorimetry to detect ligand interactions that promote protein stability.*  
  Nature Protocols. 2007;2:2212–2221.  
  https://doi.org/10.1038/nprot.2007.321

- **Ericsson et al. (stability optimization for structural biology)**  
  Ericsson UB, Hallberg BM, Detitta GT, Dekker N, Nordlund P.  
  *Thermofluor-based high-throughput stability optimization of proteins for structural studies.*  
  Analytical Biochemistry. 2006;357(2):289–298.  
  https://doi.org/10.1016/j.ab.2006.07.027

---

## Tm calculation methods used in this pipeline

### First-derivative Tm (Tm_D)

The melting temperature is the temperature at the **maximum of the first derivative** of fluorescence with respect to temperature (peak of dF/dT), i.e. the inflection of the melt curve.

$$
T_{m,D} = \arg\max_T \left(\frac{dF}{dT}\right)
$$

This matches the **Derivative Tm** definition in Thermo Fisher Protein Thermal Shift Software and is the standard approach for single- and multi-transition DSF curves in the literature (including workflows based on Niesen et al.).

### Boltzmann (two-state) Tm (Tm_B)

Fluorescence in the transition region is fit to a two-state Boltzmann-type model:

$$
F(T) = F_{\min} + \frac{F_{\max} - F_{\min}}{1 + \exp((T_m - T)/a)}
$$

where $T_m$ is the midpoint of the unfolding transition and $a$ is a steepness parameter. This matches the **Boltzmann Tm** in Thermo Protein Thermal Shift Software and is the sigmoidal form most widely used in published DSF analysis tools and papers.

**Practice (Thermo + field):** for single transitions, derivative and Boltzmann Tm are usually similar but not identical. Use **one method consistently** for all samples in a study. For **multiple melt phases**, use derivative Tm.

### ΔTm

$$
\Delta T_m = T_m(\mathrm{sample}) - T_m(\mathrm{reference})
$$

As used throughout thermal-shift screening (Pantoliano et al.; Niesen et al.; Thermo PTS studies).

---

## Thermo Scientific Protein Thermal Shift

- Thermo Fisher Scientific. *Protein Thermal Shift Software* and *Protein Thermal Shift Studies User Guide* (e.g. MAN0025600 and related manuals).  
  Official descriptions of Boltzmann Tm, Derivative Tm, regions of analysis (ROA), multi-Tm mode, and ΔTm vs reference.  
  Product overview: https://www.thermofisher.com/us/en/home/life-science/pcr/real-time-pcr/real-time-pcr-applications/real-time-pcr-protein-analysis/protein-thermal-shift.html

---

## Analysis software & methods notes (secondary)

- SimpleDSFviewer and related tools describe first-derivative and half-maximal approaches for Tm extraction from DSF curves.
- Community and core-facility guides (e.g. Harvard CMI QuantStudio DSF guides) reiterate: Thermo PTS software uses derivative and Boltzmann models; use a consistent method when comparing samples.

---

## Suggested acknowledgment

> Differential scanning fluorimetry (protein thermal shift) data were analyzed using first-derivative and/or Boltzmann two-state Tm estimation consistent with Protein Thermal Shift workflows (Thermo Fisher Scientific) and established DSF protocols (Niesen et al., Nature Protocols 2007; Pantoliano et al., Journal of Biomolecular Screening 2001). Melting temperatures were compared as ΔTm relative to reference conditions.

## Cite this repository

> Analysis used the open DSF pipeline (https://github.com/gigi777mo/dsf-pipeline), which implements platform-agnostic melt-curve import and literature-standard derivative and Boltzmann Tm calculations.
