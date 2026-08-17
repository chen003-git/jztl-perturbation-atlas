# jztl-perturbation-atlas

Analysis scripts and processed data for the manuscript:

**"PPARG as a Dual-Axis Hub Linking Myeloid Neuroinflammation and Hepatic Lipid Metabolism: A Single-Cell Perturbation Atlas of the Jiangzhuo Tongluo Formula for Ischemic Stroke"**

*Submitted to Molecular Diversity (2026)*

---

## Overview

This repository contains the computational analysis scripts and processed outputs supporting a cell-type-resolved pharmacological atlas of the Jiangzhuo Tongluo (JZTL) formula, an eight-herb traditional Chinese medicine prescription for hyperlipidemia-associated ischemic stroke.

The study integrates:
- **Network pharmacology** (148 unique targets: Jun 140 via TCMSP, Chen–Zuo–Shi 13 via literature, 5 shared)
- **Dual-tissue single-cell transcriptomics** (brain GSE174574; liver Tabula Muris)
- **Drug-target Mendelian randomization** (PPARG Pro12Ala rs1801282; null in European-ancestry data)
- **Molecular docking** (smina; 5 compounds against PPARG, 1FM6)
- **Molecular dynamics** (GROMACS 2024.4, OPLS-AA, 100 ns, quercetin–PPARG)
- **Virtual knockout** (scTenifoldKnk; brain myeloid and hepatocyte Pparg)

---

## Repository structure

```
jztl-perturbation-atlas/
├── scripts/          # 19 analysis scripts
│   ├── extract_pparg2.py          # PPARG expression extraction
│   ├── m2_polarization.py         # M2 marker polarization (PPARG+ vs PPARG−)
│   ├── extract_myeloid2.py        # Brain myeloid subset extraction
│   ├── preprocess_brain.py        # Brain scRNA preprocessing
│   ├── preprocess_liver.py        # Liver scRNA preprocessing
│   ├── coverage_permutation3.py   # Jun-target permutation test (P = 0.005)
│   ├── czs_perm.py                # Chen–Zuo–Shi permutation test (P = 0.003)
│   ├── reverify_docking3.py       # Docking re-verification (alisol B −9.3)
│   ├── dock_ursolic.py            # Ursolic acid docking (−6.9)
│   └── ...                        # additional scripts
├── results/          # 5 processed result files
│   ├── coverage_perm_result.txt   # Jun permutation output
│   ├── enrich_myeloid_result.txt  # Myeloid DRG GO enrichment
│   ├── enrich_result.txt          # Enrichment results
│   ├── drgs_summary.txt           # vKO DRG summary
│   └── liver_go_bp.csv            # Liver vKO GO BP (adjusted P ≥ 0.14)
├── ligands/          # 12 ligand SDF files (alisol B, emodin, quercetin, ursolic acid, nuciferine)
└── data/             # Processed data (DRG tables, PPARG expression)
    ├── brain_DRGs.csv             # Brain vKO differentially regulated genes
    ├── brain_myeloid_DRGs.csv     # Brain myeloid vKO DRGs (75 genes)
    └── pparg_expression.json      # PPARG expression summary
```

---

## Data sources

| Data | Source |
|---|---|
| Herb targets | TCMSP (www.tcmsp-e.com), HERB 2.0, PubMed literature |
| Brain scRNA-seq | GSE174574 (GEO) |
| Liver scRNA-seq | Tabula Muris (GSE132042) |
| eQTL | GTEx v8 |
| Stroke GWAS | MEGASTROKE (GCST006908, GCST006910, GCST006907) |
| PPARG structure | RCSB PDB 1FM6 |

---

## Key results (honest reporting)

- 148 unique targets: Jun 140 (TCMSP) + Chen–Zuo–Shi 13 (literature) − 5 shared
- Jun-herb targets co-enriched across brain cell types (permutation P = 0.005)
- Chen–Zuo–Shi targets co-enriched across liver cell types (permutation P = 0.003)
- PPARG predominantly in hepatocytes (4.8–6.6%), relatively enriched in brain myeloid cells (2.1%, LAM-like)
- Drug-target MR: **null** (cardioembolic P = 0.28; any ischemic P = 0.87); PPARG is cis-eQTL-sparse (no significant liver/adipose eQTLs, no colocalisation with stroke GWAS)
- Docking: alisol B −9.3 kcal/mol > emodin −8.2 > quercetin −8.1 > ursolic acid −6.9 > nuciferine −5.8
- vKO: brain myeloid 75 DRGs (mitotic cell-cycle, P < 1e-16); hepatocyte 78 DRGs (metabolic, adjusted P ≥ 0.14, directional)

---

## Citation

If you use these scripts or data, please cite the associated manuscript (DOI to be assigned).

## License

MIT License — see LICENSE file (to be added).
