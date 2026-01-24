# Structure resolved phylogeny

## Overview

This miniproject explores whether confounding structural signal can be minimized prior to phylogenetic modeling. Phylogenetic models, which take a multiple sequence alignment (MSA) as input, assume mutations occur independent of one another. Protein structure prediction, an adjacent field which also takes MSAs as input, relies on distilling covariation between residues to predict structure. These two fields are in conflict with one another. This small repository is an attempt at purifying phylogenetic models of structural signal, which is indicative of such covariation. Phylogenetic trees generated from a UPGMA model are compared against a WPGMA model using a downweight coefficient vector that is inversely proportional to contact order. The tree topologies differ. Ancient human mitochondrial genomes and the human Complex I heteromer are used as data. Future work might investigate a more appropriate time-dependent phylogenetic model, such as a serial-sample WPGMA model. 

### Data sources:
- [Ancient mtDNA Database (amtDB)](https://amtdb.org/) — 2,022 ancient human mtDNA sequences
- [PDB 9I4I](https://www.rcsb.org/structure/9I4I) — Human respiratory Complex I structure (45 subunits)
- [rCRS (NC_012920.1)](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1) — Revised Cambridge Reference Sequence

### files:
- [contact_map.ipynb](https://github.com/dylanmmarshall/structure_resolved_phylogeny/blob/main/contact_map.ipynb): full contact map of 9I4I + downweight vector for each Complex I mito subunit
- [mt_data.ipynb](https://github.com/dylanmmarshall/structure_resolved_phylogeny/blob/main/mt_data.ipynb): extraction / subsampling of Complex I genes from ancient mtDNA sequences + QC
- [mt_phylo.ipynb](https://github.com/dylanmmarshall/structure_resolved_phylogeny/blob/main/mt_phylo.ipynb): collation of data generated from previous two notebooks + sequence comparison of ancient mito genes versus reference gene at nucleotide and amino acid level + alignment of downweight vector against ancient mito gene sequences + weighted versus unweighted phylo tree generation / visualization
- [notebook.ipynb](https://github.com/dylanmmarshall/structure_resolved_phylogeny/blob/main/mt_phylo.ipynb): deprecated R&D analysis
  

### Contact-Based Downweighting

For each mtDNA-encoded residue, compute total structural contacts:
- Intra-chain contacts (within same subunit)
- Inter-mtDNA contacts (between mtDNA-encoded subunits)
- Mito-nuclear contacts (between mtDNA and nuclear-encoded subunits)

Downweight = 1 / (total contacts per residue)

## References

- Zhu et al. (2024) "Structure of human respiratory Complex I" [PDB 9I4I]
- Ehler et al. (2019) "AmtDB: a database of ancient human mitochondrial genomes"

Note: not peer-reviewed work
