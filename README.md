# Structure Resolved Phylogeny

**_Overview_** This project explores whether confounding structural signal can be minimized in a phylogenetic model. Phylogenetic models, which take a multiple sequence alignment (MSA) as input, assume mutations occur independent of one another. Protein structure prediction, an adjacent field which also makes use of MSAs, functions by distilling covariation between residues to predict structure. These two fields are in conflict with one another. This repository's methodology is an attempt at purifying phylogenetic models of structural signal - which is indicative of such covariation. Phylogenetic trees from an unweighted hierarchical clustering model are compared against a weighted hierarchical clustering model that uses a downweight vector that is one over the number of contacts per residue. The tree topologies differ. Ancient mitochondrial genomes and the mitochondria originating subunits of the complex I heteromer are used as data.

  **_Downweight Justification_** Natural MSAs are a subsampling of a theoretically replete MSA, which would be constructed of all possible sequence variations composed of covariation that provide at least neutral fitness advantage. 


the more i think about how {coevo} \subset {contacts} the more i begin to wonder if that is a result of MSAs themselves inherently existing as subsampled observations over all viable sequence space. in other words, a theoretically perfect MSA would be constructed of all possible sequence variations composed of nth-order covariation that provide at least neutral fitness advantage. a MRF applied to such an MSA might result in {coevo} equating {contacts}. upon such a supposition does the 1/(contact order) WPGMA method (gamma) actually provide a structure-enabled advantage over pure sequence based phylo approaches, as opposed to a decomposed-MRF-enabled strategy.

Further solidification of this work, such as accounting for time and including more mitochondrial genome samples and including the rest of the mitochondrial genome sequence, might give rise to an improved definition of the human maternal haplotype.

### Data sources:
- [Ancient mtDNA Database (amtDB)](https://amtdb.org/) — 2,022 ancient human mtDNA sequences
- [PDB 9I4I](https://www.rcsb.org/structure/9I4I) — human respiratory complex I structure (7 mitochondria subunits + 38 nuclear subunits)
- [rCRS (NC_012920.1)](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1) — revised cambridge reference sequence

### Files:
- [contact_map.ipynb](https://github.com/dylanmmarshall/structure_resolved_phylogeny/blob/main/contact_map.ipynb)
  - full contact map of 9I4I + downweight vector for each complex I mito subunit
- [mt_data.ipynb](https://github.com/dylanmmarshall/structure_resolved_phylogeny/blob/main/mt_data.ipynb)
  - extraction / subsampling of complex I genes from ancient mtDNA sequences
  - QC
- [mt_phylo.ipynb](https://github.com/dylanmmarshall/structure_resolved_phylogeny/blob/main/mt_phylo.ipynb)
  - processing of data generated from previous two notebooks
  - sequence comparison of ancient mito genes versus reference gene at nucleotide and amino acid level
  - alignment of downweight vector against ancient mito gene sequences
  - weighted versus unweighted phylo tree generation / visualization
- [notebook.ipynb](https://github.com/dylanmmarshall/structure_resolved_phylogeny/blob/main/mt_phylo.ipynb)
  - deprecated R&D analysis
  

### Contact-Based Downweighting

For each mtDNA-encoded residue, compute total structural contacts:
- Intra-chain contacts (within same subunit)
- Inter-mtDNA contacts (between mtDNA-encoded subunits)
- Mito-nuclear contacts (between mtDNA and nuclear-encoded subunits)

downweight coefficient vector = 1 / (total contacts per residue)

### Ongoing work (20260123)

- A more appropriate time-dependent phylogenetic model, such as a serial-sample WPGMA model
- Additional analysis on those ancient mitochondrial genomes that were precluded due to sequence length mismatch with the reference

## References

- Zhu et al. (2024) "Structure of human respiratory Complex I" [PDB 9I4I]
- Ehler et al. (2019) "AmtDB: a database of ancient human mitochondrial genomes"

Note: not peer-reviewed work
