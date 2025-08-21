# meerkatPaper

computational code accompanying 2025 Campbell et al:
Campbell, C. Ryan, Marta Manser, Mari Shiratori, Kelly Williams, Luis Barreiro, Tim Clutton‐Brock, and Jenny Tung. "A female‐biased gene expression signature of dominance in cooperatively breeding meerkats." Molecular Ecology 33, no. 21 (2024): e17467. https://doi.org/10.1111/mec.17467

## analysis.Rmd
base analyses of RNAseq data, starting with gene count matrix
gathers data, filters genes by logCPM, generates residual GE (to control for batch effects), run linear models (emmreml) - both basic and nested, permutation testing for FDR correction, gene set enrichment analysis (GSEA), TLR4 analysis, cibersort, and cross-species GSEA

## figures.Rmd, sFigures.Rmd, sTables.Rmd
figure and table generation from results of base analysis
using the output from the base analysis, these scripts generate the figures, supplemental figures, and supplemental tables found in the manuscript

## rScripts
.R files found here are functions used by the analysis.Rmd file
