# Citations — DSF / Protein Thermal Shift Pipeline

Cite the methods appropriate to your use. Core references for **Tm calculation** and the assay:

## Assay principle & ThermoFluor / DSF

- **Pantoliano et al. (ThermoFluor)**  
  Pantoliano MW, Petrella EC, Kwasnoski JD, et al.  
  *High-density miniaturized thermal shift assays as a general strategy for drug discovery.*  
  Journal of Biomolecular Screening. 2001;6(6):429–440.  
  https://doi.org/10.1177/108705710100600609  
  Foundational high-throughput thermal shift (dye-based) screening paper.

- **Niesen et al. (widely cited DSF protocol)**  
  Niesen FH, Berglund H, Vedadi M.  
  *The use of differential scanning fluorimetry to detect ligand interactions that promote protein stability.*  
  Nature Protocols. 2007;2:2212–2221.  
  https://doi.org/10.1038/nprot.2007.321  
  Most-cited practical protocol for DSF / thermal shift ligand screening.

- **Ericsson et al.**  
  Ericsson UB, Hallberg BM, Detitta GT, Dekker N, Nordlund P.  
  *Thermofluor-based high-throughput stability optimization of proteins for structural studies.*  
  Analytical Biochemistry. 2006;357(2):289–298.

## Tm calculation methods (used in this pipeline)

### First-derivative Tm

The melting temperature is taken as the temperature at the **maximum of the first derivative** of fluorescence with respect to temperature (peak of dF/dT), i.e. the inflection of the melt curve. This is the **Derivative Tm** in Thermo Fisher Protein Thermal Shift Software and is the standard approach for single- and multi-transition curves in the DSF field.

- Documented in Thermo Fisher Protein Thermal Shift Software / User Guides (Tm-Derivative).
- Used throughout Niesen et al. and subsequent DSF protocols for reporting Tm from melt profiles.

### Boltzmann (two-state sigmoidal) Tm

Fluorescence in the transition region is fit to a two-state Boltzmann-type model:

$$F(T) = F_{min} + \frac{F_{max} - F_{min}}{1 + \exp((T_m - T)/a)}$$

where $T_m$ is the midpoint of the unfolding transition. This is the **Boltzmann Tm** in Thermo Protein Thermal Shift Software and the most common sigmoidal fit in published DSF analysis.

- Thermo Fisher Protein Thermal Shift Software (Boltzmann fit / ROA).
- Same functional form widely used in DSF analysis tools and papers (e.g. SimpleDSFviewer; isothermal/DSF fitting literature).

**Practice (Thermo + field):** for single transitions, derivative and Boltzmann Tm are usually similar but not identical; use **one method consistently** for all comparisons. For **multiple melt phases**, use derivative Tm (and/or multi-peak detection).

## Thermo Scientific Protein Thermal Shift

- Thermo Fisher Scientific. *Protein Thermal Shift Software* and *Protein Thermal Shift Studies User Guide*.  
  Assay chemistry (Protein Thermal Shift dye), instrument compatibility (QuantStudio, StepOne, 7500, ViiA 7, etc.), and official definitions of Boltzmann Tm vs Derivative Tm.

## Reviews & practical guides

- Bai N, et al. / fluorescence-based stability monitoring reviews (DSF vs nanoDSF, ICD).
- STAR Protocols and similar guides on performing and optimizing DSF (raw RFU export, controls, Tm extraction).

## Suggested acknowledgment

> Differential scanning fluorimetry (thermal shift) data were analyzed using first-derivative and/or Boltzmann two-state Tm estimation as implemented in standard Protein Thermal Shift workflows (Thermo Fisher) and consistent with established DSF protocols (Niesen et al., Nat Protoc 2007; Pantoliano et al., J Biomol Screen 2001). Melting temperatures were compared as ΔTm relative to reference conditions.

## Cite this pipeline

> Analysis used the open DSF pipeline (https://github.com/gigi777mo/dsf-pipeline), which implements platform-agnostic melt-curve import and literature-standard derivative and Boltzmann Tm calculations.
