library(clusterProfiler)
library(org.Sc.sgd.db)
library(enrichplot)
library(tidyverse)
library(here)

# ── Load DE results ──────────────────────────────────────────────────
res_df <- read.csv(here("results", "DE_mature_vs_early.csv"),
                   row.names = 1)
res_df$gene <- rownames(res_df)
res_df <- na.omit(res_df)

# ── Check ID format ──────────────────────────────────────────────────
head(res_df$gene)
keytypes(org.Sc.sgd.db)

# ── Convert ORF IDs to Entrez IDs ───────────────────────────────────
# Gene IDs are yeast ORF names (e.g. YAL001C) — use fromType = "ORF"
# No version stripping needed for ORF IDs
gene_map <- bitr(res_df$gene,
                  fromType = "ORF",
                  toType   = c("ENTREZID", "GENENAME"),
                  OrgDb    = org.Sc.sgd.db)

res_df <- merge(res_df, gene_map,
                 by.x = "gene", by.y = "ORF",
                 all.x = TRUE)

cat("Genes mapped to Entrez IDs:", sum(!is.na(res_df$ENTREZID)),
    "/", nrow(res_df), "\n")

# ── Define gene sets ─────────────────────────────────────────────────
# Entrez IDs for GO enrichment
sig_genes <- res_df %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
    pull(ENTREZID) %>% na.omit() %>% unique()

all_genes <- res_df %>%
    pull(ENTREZID) %>% na.omit() %>% unique()

up_genes <- res_df %>%
    dplyr::filter(padj < 0.05 & log2FoldChange > 1) %>%
    pull(ENTREZID) %>% na.omit() %>% unique()

down_genes <- res_df %>%
    dplyr::filter(padj < 0.05 & log2FoldChange < -1) %>%
    pull(ENTREZID) %>% na.omit() %>% unique()

# ORF names for KEGG — yeast KEGG requires ORF format, not Entrez
sig_genes_orf <- res_df %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
    pull(gene) %>% na.omit() %>% unique()

cat("Significant:", length(sig_genes),
    "| Up:", length(up_genes),
    "| Down:", length(down_genes), "\n")

# ── GO Biological Process ORA ────────────────────────────────────────
ego_bp <- enrichGO(gene          = sig_genes,
                    universe      = all_genes,
                    OrgDb         = org.Sc.sgd.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.2,
                    readable      = FALSE)

# ── KEGG Pathway ORA ─────────────────────────────────────────────────
# Use ORF names directly — yeast KEGG requires this format
kegg_enrich <- enrichKEGG(gene         = sig_genes_orf,
                            organism     = "sce",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2)

# ── Up vs Down GO comparison ─────────────────────────────────────────
compare_obj <- compareCluster(
    geneCluster   = list(Upregulated   = up_genes,
                         Downregulated = down_genes),
    fun           = "enrichGO",
    OrgDb         = org.Sc.sgd.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
)

# ════════════════════════════════════════════════════════════════════
# FIGURE 5a — GO Biological Process dotplot
# ════════════════════════════════════════════════════════════════════
p_go <- dotplot(ego_bp, showCategory = 20,
                 title = "Figure 5a. GO Biological Process\nMature vs Early Biofilm") +
    theme(axis.text.y = element_text(size = 9))

ggsave(here("results", "figures", "Fig5a_GO_BP.pdf"), p_go, width = 8, height = 9)
ggsave(here("results", "figures", "Fig5a_GO_BP.png"), p_go, width = 8, height = 9, dpi = 300)
cat("Fig5a saved\n")

# ════════════════════════════════════════════════════════════════════
# FIGURE 5b — KEGG pathway dotplot (only if pathways found)
# ════════════════════════════════════════════════════════════════════
if (!is.null(kegg_enrich) && nrow(as.data.frame(kegg_enrich)) > 0) {
    p_kegg <- dotplot(kegg_enrich, showCategory = 15,
                       title = "Figure 5b. KEGG Pathway Enrichment\nMature vs Early Biofilm")
    ggsave(here("results", "figures", "Fig5b_KEGG.pdf"), p_kegg, width = 8, height = 7)
    ggsave(here("results", "figures", "Fig5b_KEGG.png"), p_kegg, width = 8, height = 7, dpi = 300)
    write.csv(as.data.frame(kegg_enrich), here("results", "KEGG_mature_vs_early.csv"))
    cat("KEGG enrichment: found", nrow(as.data.frame(kegg_enrich)), "pathways\n")
} else {
    cat("KEGG: no significantly enriched pathways found at current thresholds\n")
}

# ════════════════════════════════════════════════════════════════════
# FIGURE 5c — Up vs Downregulated GO terms compared
# ════════════════════════════════════════════════════════════════════
p_compare <- dotplot(compare_obj, showCategory = 15,
                      title = "Figure 5c. GO BP: Up vs Downregulated Genes\nMature vs Early Biofilm") +
  theme(
    axis.text.y = element_text(size = 8, lineheight = 0.9),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 60)
  )

ggsave(here("results", "figures", "Fig5c_GO_compare.pdf"), p_compare, width = 10, height = 9)
ggsave(here("results", "figures", "Fig5c_GO_compare.png"), p_compare, width = 10, height = 9, dpi = 300)
cat("Fig5c saved\n")

# ── Export results tables ────────────────────────────────────────────
write.csv(as.data.frame(ego_bp), here("results", "GO_BP_mature_vs_early.csv"))

cat("Done! All figures saved to results/figures/\n")

