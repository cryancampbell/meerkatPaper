# meerkatPaper

computational code accompanying 2025 Campbell et al:
Campbell, C. Ryan, Marta Manser, Mari Shiratori, Kelly Williams, Luis Barreiro, Tim Clutton‐Brock, and Jenny Tung. "A female‐biased gene expression signature of dominance in cooperatively breeding meerkats." Molecular Ecology 33, no. 21 (2024): e17467. https://doi.org/10.1111/mec.17467

Abstract:
Dominance is a primary determinant of social dynamics and resource access in social animals. Recent studies show that dominance is also reflected in the gene regulatory profiles of peripheral immune cells. However, the strength and direction of this relationship differs across the species and sex combinations investigated, potentially due to variation in the predictors and energetic consequences of dominance status. Here, we investigated the association between social status and gene expression in the blood of wild meerkats (Suricata suricatta; n=113 individuals), including in response to lipopolysaccharide, Gardiquimod (an agonist of TLR7, which detects single-stranded RNA in vivo) and glucocorticoid stimulation. Meerkats are cooperatively breeding social carnivores in which breeding females physically outcompete other females to suppress reproduction, resulting in high reproductive skew. They therefore present an opportunity to disentangle the effects of social dominance from those of sex per se. We identify a sex-specific signature of dominance, including 1045 differentially expressed genes in females but none in males. Dominant females exhibit elevated activity in innate immune pathways and a larger fold-change response to LPS challenge. Based on these results and a preliminary comparison to other mammals, we speculate that the gene regulatory signature of social status in the immune system depends on the determinants and energetic costs of social dominance, such that it is most pronounced in hierarchies where physical competition is important and reproductive skew is large. Such a pattern has the potential to mediate life history trade-offs between investment in reproduction versus somatic maintenance.


## analysis.Rmd
base analyses of RNAseq data, starting with gene count matrix
gathers data, filters genes by logCPM, generates residual GE (to control for batch effects), run linear models (emmreml) - both basic and nested, permutation testing for FDR correction, gene set enrichment analysis (GSEA), TLR4 analysis, cibersort, and cross-species GSEA

## figures.Rmd, sFigures.Rmd, sTables.Rmd
figure and table generation from results of base analysis
using the output from the base analysis, these scripts generate the figures, supplemental figures, and supplemental tables found in the manuscript

## rScripts
.R files found here are functions used by the analysis.Rmd file
