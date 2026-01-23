# Pilot Analysis: Contact-Weighted Mitochondrial Phylogenetics

Phylogenetic analysis of human mitochondrial DNA using structural contact information from Complex I to weight sequence distances.

## Overview

This project explores whether structural constraints in mitochondrial-encoded proteins can improve phylogenetic inference. The key idea: residues with many structural contacts are under stronger purifying selection, so mutations there should be downweighted when computing evolutionary distances.

**Data sources:**
- [Ancient mtDNA Database (amtDB)](https://amtdb.org/) — 2,022 ancient human mtDNA sequences
- [PDB 9I4I](https://www.rcsb.org/structure/9I4I) — Human respiratory Complex I structure (45 subunits)
- [rCRS (NC_012920.1)](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1) — Revised Cambridge Reference Sequence

## Notebooks

| Notebook | Description |
|----------|-------------|
| `contact_map.ipynb` | Parse Complex I structure, compute inter-chain contact maps, generate per-residue downweight vectors for mtDNA-encoded subunits (ND1-ND6, ND4L) |
| `mt_data.ipynb` | Extract ND gene sequences from ancient samples using GFF3 annotations, analyze sequence variability |
| `mt_phylo.ipynb` | Build distance matrices and UPGMA trees comparing unweighted vs contact-weighted approaches, QC filtering |
| `notebook.ipynb` | MSA preparation utilities, COX2 analysis (deprecated) |

## Methods

### Contact-Based Downweighting

For each mtDNA-encoded residue, compute total structural contacts:
- Intra-chain contacts (within same subunit)
- Inter-mtDNA contacts (between mtDNA-encoded subunits)
- Mito-nuclear contacts (between mtDNA and nuclear-encoded subunits)

Downweight = 1 / total_contacts (positions with more contacts contribute less to distance)

### Distance Metrics

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| Unweighted Hamming | Σ(mismatches) | All positions equal |
| Contact-weighted | Σ(weight × mismatch) | Buried residues downweighted |

### QC Filtering

Samples excluded based on:
- Frameshift detection (avg AA similarity < 50%)
- High missing data (any gene > 80% Ns)
- Zero similarity in any gene

## Data Structure

```
pilot_analysis/
├── data/
│   ├── amtdb_seqs/           # 2,022 ancient mtDNA FASTA files
│   ├── amtdb_metadata*.csv   # Sample metadata (dates, locations, haplogroups)
│   ├── mito_contact_summary.npz  # Per-residue contact counts and downweights
│   ├── NC_012920.1.gff3      # Gene annotations
│   └── rCRS.fa               # Reference sequence
├── 9I4I.cif                  # Complex I structure
├── contact_map.ipynb
├── mt_data.ipynb
├── mt_phylo.ipynb
└── notebook.ipynb
```

## Key Results

**Complex I subunit coverage:**
- 7 mtDNA-encoded subunits: MT-ND1, MT-ND2, MT-ND3, MT-ND4, MT-ND4L, MT-ND5, MT-ND6
- 38 nuclear-encoded subunits
- 2,109 structurally-resolved codons (6,327 nucleotide positions)

**Ancient samples:**
- 1,374 samples with reference-length sequences
- 1,365 pass QC filtering
- Date range: 200 - 13,980 years BP

**Distance statistics (n=1,365):**
| Metric | Mean | Median | Max |
|--------|------|--------|-----|
| Unweighted | 302.4 | 18.0 | 3,724 |
| Weighted | 55.3 | 3.5 | 708 |

## Requirements

```
numpy
pandas
scipy
matplotlib
biopython
```

## References

- Zhu et al. (2024) "Structure of human respiratory Complex I" [PDB 9I4I]
- Ehler et al. (2019) "AmtDB: a database of ancient human mitochondrial genomes"
