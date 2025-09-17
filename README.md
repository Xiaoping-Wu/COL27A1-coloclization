##--------------- Colocalization ---------------------------
## @Xiaoping WU
## 2025-08-20
Colocalization

Colocalization is used to verify whether two association signals - from GWAS and eQTL studies, in the same genomic region are likely to be driven by the same causal variant.

We selected the GWAS summary statistics from Moba data (DOI: 10.1101/2025.06.17.25329777) for COL27A1 region, chr9: 115 416 257 – 118 415 208 (hg19)

eQTL data was downloaded from eQTL Catalogue database (https://www.ebi.ac.uk/eqtl/).

R package coloc was used for colocalization. posterior probability PPH4 > 0.8 indicating same causal variant for two traits: GWAS and eQTL.




