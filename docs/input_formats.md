# Input formats (any qPCR brand)

## Universal wide CSV (recommended)

Export or convert instrument data to:

| Temperature | A1 | A2 | ... |
|-------------|----|----|-----|
| 25.0 | RFU | RFU | ... |

Column name for temperature may be `Temperature`, `Temp`, `T`, or `Temperature (°C)`.

## Long CSV

| Temperature | Well | Fluorescence |
|-------------|------|--------------|
| 25.0 | A1 | 1200 |

## Instrument notes

| Platform | What to export |
|----------|----------------|
| **Thermo QuantStudio / PTS** | Raw melt fluorescence (not only analyzed Tm). Convert EDS-derived tables to wide CSV if needed. |
| **Bio-Rad CFX** | Melt curve RFU export → CSV |
| **Roche LightCycler** | Absolute fluorescence vs temperature export |
| **ABI StepOne / 7500** | Multicomponent or raw data export → temperature + reporter |
| **Generic** | Any table with T and per-well fluorescence |

**Critical:** use **raw (or noise-reduced) fluorescence vs temperature**, not Ct / amplification results.

## Sample sheet

Optional but recommended for ΔTm and replicate stats:

```csv
well,sample,condition,replicate,is_reference
A1,Lysozyme,control,1,true
B1,Target,DMSO,1,false
B2,Target,ligand_10uM,1,false
```
