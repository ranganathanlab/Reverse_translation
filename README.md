# Reverse_translation
---
> 12.2022
>
> Copyright (C) 2022 Xinran Lian
>

## Introduction
Reverse translation to design gene oligos from protein sequences used in https://www.biorxiv.org/content/10.1101/2022.12.21.521443v1


Structure of the plasmid (*PRS316*) is depicted in

Zarrinpar, A., Park, S. H., & Lim, W. A. (2003). Optimization of specificity in a cellular protein interaction network by negative selection. Nature, 426(6967), 676-680.
https://www.nature.com/articles/nature02178

Plasmid sequence maps will be uploaded later.

## Getting Started

|            |                                                         |
| :---       | :---                                                    |
| Inputs/yeast_codon.xlsx      | Yeast codon usage table   |
| Inputs/Final_New_Proteins_nogap.fasta Inputs/test.fasta| Two demo fasta files |
| Inputs/twist_red_seqs_forblast.an | A local gene library as BLAST reference |
| Reverse_translation_300mer.ipynb    | Notebook for reverse translation, plase follow instruction inside  |
| local.m | See instructions in *Reverse_translation_300mer.ipynb*   |
| remove_gap_fasta.py | Script to remove gaps from input file |

## SH3 oligo specification
---
**Oligo sequence structure** (300 mer)
![oligo sequence structure](oligos300mer.png)

---
**check oligo workflow** (300 mer)
![check oligo workflow](seq_structure_check.png)