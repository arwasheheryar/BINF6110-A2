# ═══════════════════════════════════════════════════════════════════════
# Assignment 2 — Differential Expression Analysis
# Dataset:    Mardanov et al. 2020 (PRJNA592304)
# Organism:   Saccharomyces cerevisiae — yeast biofilm (velum) development
# Data type:  Single-end bulk RNA-seq
# Groups:     Early Biofilm (Day 38), Thin Biofilm (Day 83), Mature Biofilm (Day 109)
# ═══════════════════════════════════════════════════════════════════════

library(txdbmaker)
library(DESeq2)
library(tximport)
library(GenomicFeatures)
library(AnnotationDbi)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(here)

# here() anchors to the project root (where .here file is)
# Confirm it found the right root:
here() 


# STEP 1: Build tx2gene mapping from GTF ──────────────────────────

# ── STEP 1: Build tx2gene mapping from GTF ──────────────────────────
txdb <- txdbmaker::makeTxDbFromGFF(
            here("GCF_000146045.2_R64_genomic.gtf"),
            format = "gtf")

k       <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# Strip version numbers from TXNAME to match quant.sf IDs
tx2gene$TXNAME <- sub("\\..*", "", tx2gene$TXNAME)
head(tx2gene)

# STEP 2: Load metadata ────────────────────────────────────────────

sampleTable <- read.csv(here("scripts", "sample_table.csv"),
                         row.names = 1)

sampleTable$stage <- factor(sampleTable$stage,
                             levels = c("Early_Biofilm",
                                        "Thin_Biofilm",
                                        "Mature_Biofilm"))


# STEP 3: Build paths to quant.sf files ───────────────────────────
files <- file.path(here("quants"), rownames(sampleTable), "quant.sf")
names(files) <- rownames(sampleTable)
stopifnot(all(file.exists(files)))
cat("All 9 quant.sf files found.\n")


# STEP 4: Import counts with tximport ─────────────────────────────
txi <- tximport(files,
                type            = "salmon",
                tx2gene         = tx2gene,
                ignoreTxVersion = TRUE)

dim(txi$counts)
head(txi$counts)


# STEP 5: Create DESeq2 dataset ───────────────────────────────────
dds  <- DESeqDataSetFromTximport(txi, sampleTable, ~stage)
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat("Genes remaining after filtering:", nrow(dds), "\n")
dds  <- DESeq(dds)
resultsNames(dds)


# STEP 6: Extract pairwise results with LFC shrinkage ─────────────
res_thin_vs_early    <- results(dds,
                                 contrast = c("stage", "Thin_Biofilm", "Early_Biofilm"),
                                 alpha    = 0.05)
resLFC_thin_vs_early <- lfcShrink(dds,
                                   contrast = c("stage", "Thin_Biofilm", "Early_Biofilm"),
                                   type     = "ashr")

res_mature_vs_early    <- results(dds,
                                   contrast = c("stage", "Mature_Biofilm", "Early_Biofilm"),
                                   alpha    = 0.05)
resLFC_mature_vs_early <- lfcShrink(dds,
                                     contrast = c("stage", "Mature_Biofilm", "Early_Biofilm"),
                                     type     = "ashr")

res_mature_vs_thin    <- results(dds,
                                  contrast = c("stage", "Mature_Biofilm", "Thin_Biofilm"),
                                  alpha    = 0.05)
resLFC_mature_vs_thin <- lfcShrink(dds,
                                    contrast = c("stage", "Mature_Biofilm", "Thin_Biofilm"),
                                    type     = "ashr")

cat("\n── Thin vs Early ──\n");   summary(res_thin_vs_early)
cat("\n── Mature vs Early ──\n"); summary(res_mature_vs_early)
cat("\n── Mature vs Thin ──\n");  summary(res_mature_vs_thin)

# Export DE gene lists
write.csv(as.data.frame(resLFC_thin_vs_early),
          here("results", "DE_thin_vs_early.csv"))
write.csv(as.data.frame(resLFC_mature_vs_early),
          here("results", "DE_mature_vs_early.csv"))
write.csv(as.data.frame(resLFC_mature_vs_thin),
          here("results", "DE_mature_vs_thin.csv"))


# STEP 7 : Likelihood Ratio Test ────────────────────────
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~1)
res_lrt  <- results(dds_lrt)
summary(res_lrt)
write.csv(as.data.frame(res_lrt), here("results", "DE_LRT_timecourse.csv"))


# ══════════════════════════════════════════════════════════════════════
# FIGURES
# ═══════════════════════════════════════════
vsd <- vst(dds, blind = FALSE)


# ── Figure 1: PCA ────────────────────────────────────────────────────
pca_data   <- plotPCA(vsd, intgroup = "stage", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = stage)) +
    geom_point(size = 5, alpha = 0.9) +
    geom_text_repel(aes(label = name), size = 3, max.overlaps = 20) +
    scale_color_manual(values = c(
        "Early_Biofilm"  = "#2196F3",
        "Thin_Biofilm"   = "#FF9800",
        "Mature_Biofilm" = "#4CAF50"
    )) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    labs(title    = "Figure 1. PCA of Yeast Biofilm Development",
         subtitle = "VST-normalised counts; n=9 samples across three biofilm stages",
         color    = "Biofilm Stage") +
    theme_classic(base_size = 13)

ggsave(here("results", "figures", "Fig1_PCA.pdf"), p_pca, width = 7, height = 5)
ggsave(here("results", "figures", "Fig1_PCA.png"), p_pca, width = 7, height = 5, dpi = 300)


# ── Figure 2: Volcano Plot ───────────────────────────────────────────
res_df <- as.data.frame(resLFC_mature_vs_early)
res_df$gene <- rownames(res_df)
res_df <- na.omit(res_df)

res_df$category <- case_when(
    res_df$padj < 0.05 & res_df$log2FoldChange >  1 ~ "Upregulated",
    res_df$padj < 0.05 & res_df$log2FoldChange < -1 ~ "Downregulated",
    TRUE ~ "Not Significant"
)

top_labels <- res_df %>%
    filter(category != "Not Significant") %>%
    arrange(padj) %>%
    slice_head(n = 15)

p_volcano <- ggplot(res_df, aes(x = log2FoldChange,
                                 y = -log10(pvalue),
                                 color = category)) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_text_repel(data = top_labels,
                    aes(label = gene),
                    size = 2.8, color = "black",
                    max.overlaps = 20,
                    box.padding = 0.4) +
    scale_color_manual(values = c(
        "Upregulated"     = "#E53935",
        "Downregulated"   = "#1565C0",
        "Not Significant" = "grey70"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5) +
    labs(x     = "Log2 Fold Change (Mature vs Early Biofilm)",
         y     = "-Log10 p-value",
         title = "Figure 2. Volcano Plot: Mature vs Early Biofilm",
         color = "Expression") +
    theme_classic(base_size = 13)

ggsave(here("results", "figures", "Fig2_Volcano.pdf"), p_volcano, width = 8, height = 6)
ggsave(here("results", "figures", "Fig2_Volcano.png"), p_volcano, width = 8, height = 6, dpi = 300)


# ── Figure 3: Heatmap ────────────────────────────────────────────────
top40 <- resLFC_mature_vs_early %>%
    as.data.frame() %>%
    na.omit() %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    rownames() %>%
    head(40)

mat        <- assay(vsd)[top40, ]
mat_scaled <- t(scale(t(mat)))

annotation_col <- data.frame(
    Stage = sampleTable$stage,
    row.names = rownames(sampleTable)
)
ann_colors <- list(Stage = c(
    "Early_Biofilm"  = "#2196F3",
    "Thin_Biofilm"   = "#FF9800",
    "Mature_Biofilm" = "#4CAF50"
))

p_heat <- pheatmap(
    mat_scaled,
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    annotation_col    = annotation_col,
    annotation_colors = ann_colors,
    show_colnames     = FALSE,
    show_rownames     = TRUE,
    fontsize_row      = 7,
    color             = colorRampPalette(c("#1565C0", "white", "#E53935"))(100),
    main              = "Figure 3. Top 40 DE Genes: Mature vs Early Biofilm\n(row z-scores of VST-normalised counts)"
)

ggsave(here("results", "figures", "Fig3_Heatmap.pdf"), p_heat, width = 9, height = 11)
ggsave(here("results", "figures", "Fig3_Heatmap.png"), p_heat, width = 9, height = 11, dpi = 300)


# ── Figure 4: MA Plot ────────────────────────────────────────────────
pdf(here("results", "figures", "Fig4_MA_plots.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plotMA(res_mature_vs_early,
       ylim = c(-6, 6), main = "Before LFC Shrinkage", colSig = "red")
plotMA(resLFC_mature_vs_early,
       ylim = c(-6, 6), main = "After LFC Shrinkage (ashr)", colSig = "red")
dev.off()
