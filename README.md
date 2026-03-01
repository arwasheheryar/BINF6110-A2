# Yeast Biofilm RNA-seq Analysis
Arwa Sheheryar  
BINF6110 – Assignment 2
---
## Introduction

Biofilms are structured microbial communities encased in a self-produced extracellular matrix that confer coordinated resistance to environmental stressors, including antimicrobials and ethanol (Yang et al., 2018). Their formation has broad medical and industrial relevance. In clinical settings, biofilm-associated tolerance contributes to persistent infections, while in fermentation contexts, biofilm development influences community stability and metabolic output (Yang et al., 2018). In Saccharomyces cerevisiae, biofilm formation, including the velum that develops on the surface of aging wine, is largely governed by the FLO gene family, which differentially regulates cell–cell and cell–surface adhesion and makes yeast a strong genetic model for studying adhesion-mediated community behavior (Bester et al., 2006; Yang et al., 2018). Biofilm maturation represents a biologically distinct phase characterized by increased structural complexity and dynamic remodeling of gene expression, reflecting adaptive responses to changing chemical and nutritional conditions over time (Yang et al., 2018).

RNA-seq is particularly well suited for examining these transcriptional dynamics because it enables unbiased, genome-wide quantification of transcript abundance without relying on pre-designed probes (Nookaew et al., 2012). Compared to microarray platforms, RNA-seq offers a substantially broader dynamic range, allowing more accurate measurement of both low- and high-abundance transcripts and detecting a greater number of differentially expressed genes with larger fold changes in S. cerevisiae datasets (Nookaew et al., 2012; Zhao et al., 2014). Its sensitivity to rare and transient transcripts, including cryptic RNAs and alternative isoforms that often fall below array detection thresholds, strengthens its ability to capture the full transcriptional complexity of yeast biofilms (Zhao et al., 2014). Because RNA-seq avoids cross-hybridization artifacts, probe saturation, and fixed annotation constraints inherent to microarrays, it also provides improved statistical power for detecting differential expression among lowly expressed and downregulated genes when paired with count-based modeling frameworks (Nookaew et al., 2012).

While differential expression analysis identifies individual genes that change significantly between conditions, gene lists alone offer limited mechanistic insight without contextualization into coordinated biological functions (Reimand et al., 2019). Over-representation analysis (ORA) addresses this by statistically testing whether predefined biological categories from Gene Ontology (GO) or KEGG pathway databases are enriched among a threshold-defined set of significant genes relative to an appropriate background, typically using Fisher’s exact test (Reimand et al., 2019; Ziemann et al., 2024). ORA is appropriate when using a predefined significance threshold and aims to identify over-represented biological processes among significantly differentially expressed genes, whereas GSEA evaluates ranked gene lists without arbitrary cutoffs and detects coordinated but subtle shifts across all genes (Reimand et al., 2019; Yoon et al., 2016). ORA was selected here because this study focused on a rigorously defined DE gene set derived from stringent statistical modeling, enabling clear biological interpretation while maintaining transparency in background specification and multiple-testing correction (Ziemann et al., 2024).

This study reanalyzes publicly available RNA-seq data from Mardanov et al. (2020) to characterize transcriptional differences across three stages of yeast velum biofilm development: early, thin, and mature, sampled at 38, 83, and 109 days post-inoculation, respectively. Differential expression was performed using DESeq2, which applies a negative binomial model with regularized log fold change shrinkage to appropriately model overdispersion and limited replicate sizes typical of RNA-seq experiments (Love et al., 2014). Functional enrichment using GO biological process and KEGG pathway ORA was then applied to identify biological processes associated with biofilm progression.

---
## Methods
---
## Results
### Overall Data Structure
### Differential Expression
### Functional Enrichment (ORA)
---
## Discussion
---
## Reproducibility

## References
Bester, M. C., Pretorius, I. S., & Bauer, F. F. (2006). The regulation of Saccharomyces cerevisiae FLO gene expression and its role in yeast biofilm formation. Archives of Microbiology, 186(2), 99–111. 
Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550.
Nookaew, I., Papini, M., Pornputtapong, N., Scalcinati, G., & Nielsen, J. (2012). A comprehensive comparison of RNA-seq and microarray technologies for transcriptome analysis in Saccharomyces cerevisiae. BMC Genomics, 13, 621.

Reimand, J., Isserlin, R., Voisin, V., Kucera, M., Tannus-Lopes, C., Rostamianfar, A., Wadi, L., Meyer, M., Wong, J., Xu, C., Merico, D., & Bader, G. D. (2019). Pathway enrichment analysis and visualization of omics data using g:Profiler, GSEA, Cytoscape and EnrichmentMap. Nature Protocols, 14(2), 482–517. 

Yoon, S., Kim, S.-Y., & Nam, D. (2016). Improving gene-set enrichment analysis of RNA-seq data with small replicates. PLOS ONE, 11(11), e0165919. 

Zhao, S., Fung-Leung, W.-P., Bittner, A., Ngo, K., & Liu, X. (2014). Comparison of RNA-Seq and microarray in transcriptome profiling of activated T cells. PLoS ONE, 9(1), e78644.

Ziemann, M., Schroeter, B., & Bora, A. (2024). Two subtle problems with overrepresentation analysis. Bioinformatics Advances, 4(1), vbae159. 

Yang, Y., Guo, Z., & Wang, X. (2018). The regulatory network of Saccharomyces cerevisiae biofilm formation. Frontiers in Microbiology, 9, 1860.
