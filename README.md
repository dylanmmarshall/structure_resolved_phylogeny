# Contact-Weighted Mitochondrial Phylogenetics

Phylogenetic analysis of human mitochondrial DNA using structural contact information from Complex I to weight sequence distances.

## Overview

This project explores whether structural constraints in mitochondrial-encoded proteins can improve phylogenetic inference. The following is the key idea. Phylogenetic models, which take a multiple sequence alignment (MSA) as input, assume mutations occur independent of one another. In an adjacent field, protein structure prediction, which also takes MSAs as input, makes use of covariation between residues to predict structure. These two fields are in conflict with one another. This small repository aims to purify phylogenetic models with structural signal, which is indicative of such covariation. A structure-resolved phylogeny, in other words.

**Data sources:**
- [Ancient mtDNA Database (amtDB)](https://amtdb.org/) — 2,022 ancient human mtDNA sequences
- [PDB 9I4I](https://www.rcsb.org/structure/9I4I) — Human respiratory Complex I structure (45 subunits)
- [rCRS (NC_012920.1)](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1) — Revised Cambridge Reference Sequence

### Contact-Based Downweighting

For each mtDNA-encoded residue, compute total structural contacts:
- Intra-chain contacts (within same subunit)
- Inter-mtDNA contacts (between mtDNA-encoded subunits)
- Mito-nuclear contacts (between mtDNA and nuclear-encoded subunits)

Downweight = 1 / (total_contacts per residue)

## References

- Zhu et al. (2024) "Structure of human respiratory Complex I" [PDB 9I4I]
- Ehler et al. (2019) "AmtDB: a database of ancient human mitochondrial genomes"
