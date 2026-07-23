# 03_dropout_and_wnn.R
# PBMC-CITE-seq-ADT-azimuth-mapping
#
# Purpose:
#   1. Dropout rescue: does ADT retain signal in cells where RNA fails to
#      detect a marker gene (dropout)?
#   2. WNN vs. RNA-only clustering: does incorporating the ADT dimension
#      change clustering, quantified via Adjusted Rand Index (ARI)?
#
# Both analyses load from the saved object at the end of notebook 02.

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

d2_pbmc_10x_CITE <- readRDS("data/processed/d2_pbmc_10x_CITE_full_analysis.rds")

# ---- 1. Dropout rescue: CD3 ADT vs. CD3D/E/G RNA dropout ----

# Ground truth for "is this a T cell" comes from Azimuth (predicted.celltype.l1),
# which uses the full transcriptome + protein context - not circular with the
# single-gene RNA dropout test below.
t_cells <- colnames(d2_pbmc_10x_CITE)[d2_pbmc_10x_CITE$predicted.celltype.l1 %in% c("CD4 T", "CD8 T")]
b_cells <- colnames(d2_pbmc_10x_CITE)[d2_pbmc_10x_CITE$predicted.celltype.l1 == "B"]

# raw RNA counts (not normalized) for CD3 complex genes, within T cells
cd3_genes <- c("CD3D", "CD3E", "CD3G")
rna_counts_t <- GetAssayData(d2_pbmc_10x_CITE, assay = "RNA", layer = "counts")[cd3_genes, t_cells, drop = FALSE]
cd3_rna_sum_t <- colSums(as.matrix(rna_counts_t))

dropout_cells   <- t_cells[cd3_rna_sum_t == 0]
detected_cells  <- t_cells[cd3_rna_sum_t > 0]

message("T cells with CD3 RNA dropout (all 3 genes = 0 counts): ", length(dropout_cells),
        " of ", length(t_cells), " T cells (",
        round(100 * length(dropout_cells) / length(t_cells), 1), "%)")

# CLR-normalized ADT CD3 signal for each group
adt_cd3 <- GetAssayData(d2_pbmc_10x_CITE, assay = "ADT", layer = "data")["CD3", ]

dropout_df <- bind_rows(
  data.frame(group = "T cell, CD3 RNA dropout", adt_cd3 = adt_cd3[dropout_cells]),
  data.frame(group = "T cell, CD3 RNA detected", adt_cd3 = adt_cd3[detected_cells]),
  data.frame(group = "B cell (true negative)", adt_cd3 = adt_cd3[b_cells])
)

dropout_df$group <- factor(dropout_df$group,
                            levels = c("T cell, CD3 RNA dropout",
                                       "T cell, CD3 RNA detected",
                                       "B cell (true negative)"))

fig_dropout <- ggplot(dropout_df, aes(x = group, y = adt_cd3, fill = group)) +
  geom_violin(color = NA, alpha = 0.8) +
  geom_jitter(size = 0.3, alpha = 0.3, width = 0.15, color = "black") +
  scale_fill_manual(values = c("T cell, CD3 RNA dropout" = "#F4A582",
                                "T cell, CD3 RNA detected" = "#92C5DE",
                                "B cell (true negative)" = "grey70")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8, angle = 20, hjust = 1),
        plot.title = element_text(size = 11, face = "plain")) +
  labs(title = "ADT CD3 signal rescues RNA dropout in T cells",
       subtitle = paste0(round(100 * length(dropout_cells) / length(t_cells), 1),
                          "% of Azimuth-defined T cells have zero RNA counts (raw) for CD3D/E/G"),
       x = NULL, y = "CD3 ADT (CLR-normalized)")

fig_dropout

ggsave("figures/03_fig_dropout_rescue.png", fig_dropout, width = 7, height = 5)

# ---- 2. WNN vs. RNA-only clustering: does ADT change clustering? ----

# ADT PCA (max 9 PCs since panel has 10 features)
d2_pbmc_10x_CITE <- ScaleData(d2_pbmc_10x_CITE, assay = "ADT")
d2_pbmc_10x_CITE <- RunPCA(d2_pbmc_10x_CITE, assay = "ADT", reduction.name = "apca",
                            npcs = 9, features = rownames(d2_pbmc_10x_CITE[["ADT"]]))

Reductions(d2_pbmc_10x_CITE)

ncol(Embeddings(d2_pbmc_10x_CITE, "pca"))
ncol(Embeddings(d2_pbmc_10x_CITE, "apca"))

# WNN: builds joint graph from RNA PCA (already computed as "pca" in notebook 02)
# and ADT PCA ("apca"), with per-cell modality weights
# there are only 8 (instead of 9) PCs for ADT, because the zero-variance isotype control was dropped
d2_pbmc_10x_CITE <- FindMultiModalNeighbors(d2_pbmc_10x_CITE,
                                            reduction.list = list("pca", "apca"),
                                            dims.list = list(1:30, 1:8))

d2_pbmc_10x_CITE <- RunUMAP(d2_pbmc_10x_CITE, nn.name = "weighted.nn",
                             reduction.name = "wnn.umap")

d2_pbmc_10x_CITE <- FindClusters(d2_pbmc_10x_CITE, graph.name = "RNA_snn",
                                 resolution = 0.8, cluster.name = "rna_only_clusters")

# manual ARI (Adjusted Rand Index) - avoids adding a new package dependency
adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  sum_comb_row <- sum(choose(rowSums(tab), 2))
  sum_comb_col <- sum(choose(colSums(tab), 2))
  sum_comb <- sum(choose(tab, 2))
  expected <- (sum_comb_row * sum_comb_col) / choose(n, 2)
  max_index <- 0.5 * (sum_comb_row + sum_comb_col)
  (sum_comb - expected) / (max_index - expected)
}

ari_rna_vs_wnn  <- adjusted_rand_index(d2_pbmc_10x_CITE$rna_only_clusters,
                                       d2_pbmc_10x_CITE$wnn_clusters)
ari_rna_vs_azimuth <- adjusted_rand_index(d2_pbmc_10x_CITE$rna_only_clusters,
                                          d2_pbmc_10x_CITE$predicted.celltype.l2)
ari_wnn_vs_azimuth <- adjusted_rand_index(d2_pbmc_10x_CITE$wnn_clusters,
                                          d2_pbmc_10x_CITE$predicted.celltype.l2)

message("ARI, RNA-only clusters vs. WNN clusters: ", round(ari_rna_vs_wnn, 3))
message("ARI, RNA-only clusters vs. Azimuth l2:   ", round(ari_rna_vs_azimuth, 3))
message("ARI, WNN clusters vs. Azimuth l2:        ", round(ari_wnn_vs_azimuth, 3))

Graphs(d2_pbmc_10x_CITE)

# side-by-side UMAPs, colored by predicted.celltype.l2 for visual comparison,
# reusing the labeled-UMAP helper style from notebook 02
plot_simple_umap <- function(embeddings, group_vals, title) {
  dim1_name <- colnames(embeddings)[1]
  dim2_name <- colnames(embeddings)[2]
  df <- as.data.frame(embeddings)
  df$group <- factor(group_vals)

  ggplot(df, aes(x = .data[[dim1_name]], y = .data[[dim2_name]], color = group)) +
    geom_point(size = 0.5, alpha = 0.6) +
    theme_classic() +
    theme(plot.title = element_text(size = 11, face = "plain"),
          legend.position = "none") +
    labs(title = title, x = dim1_name, y = dim2_name)
}

fig_rna_umap <- plot_simple_umap(Embeddings(d2_pbmc_10x_CITE, "umap"),
                                  d2_pbmc_10x_CITE$predicted.celltype.l2,
                                  "RNA-only clustering (colored by Azimuth l2)")

fig_wnn_umap <- plot_simple_umap(Embeddings(d2_pbmc_10x_CITE, "wnn.umap"),
                                  d2_pbmc_10x_CITE$predicted.celltype.l2,
                                  "WNN clustering (RNA + ADT, colored by Azimuth l2)")

fig_wnn_compare <- (fig_rna_umap + fig_wnn_umap) +
  plot_annotation(
    title = paste0("ARI (RNA-only vs. WNN): ", round(ari_rna_vs_wnn, 3),
                   "  |  ARI (RNA-only vs. Azimuth): ", round(ari_rna_vs_azimuth, 3),
                   "  |  ARI (WNN vs. Azimuth): ", round(ari_wnn_vs_azimuth, 3)),
    theme = theme(plot.title = element_text(size = 9))
  )

fig_wnn_compare

ggsave("figures/03_fig_wnn_vs_rna_simple_umap.png", fig_wnn_compare, width = 12, height = 5)



#######################
## re plot using notebook 02 more labeled umap fxn

fig_rna_umap2 <- plot_labeled_umap(Embeddings(d2_pbmc_10x_CITE, "umap"),
                                  d2_pbmc_10x_CITE$predicted.celltype.l2,
                                  "RNA-only clustering (colored by Azimuth l2)",
                                  show_legend = FALSE, label_size = 2.5)

fig_wnn_umap2 <- plot_labeled_umap(Embeddings(d2_pbmc_10x_CITE, "wnn.umap"),
                                  d2_pbmc_10x_CITE$predicted.celltype.l2,
                                  "WNN clustering (RNA + ADT, colored by Azimuth l2)",
                                  show_legend = FALSE, label_size = 2.5)

fig_wnn_compare2 <- (fig_rna_umap2 + fig_wnn_umap2) +
  plot_annotation(
    title = paste0("ARI (RNA-only vs. WNN): ", round(ari_rna_vs_wnn, 3),
                   "  |  ARI (RNA-only vs. Azimuth): ", round(ari_rna_vs_azimuth, 3),
                   "  |  ARI (WNN vs. Azimuth): ", round(ari_wnn_vs_azimuth, 3)),
    theme = theme(plot.title = element_text(size = 9))
  )

fig_wnn_compare2

ggsave("figures/03_fig_wnn_vs_rna_labeledByType.png", fig_wnn_compare2, width = 12, height = 5)


####### table with cell counts #############

library(gridExtra)
library(grid)

# helper: build count + % table for a cluster column, save as png
save_cluster_table <- function(obj, column, filename) {
  counts <- table(obj@meta.data[[column]])
  df <- data.frame(
    cluster = names(counts),
    n = as.integer(counts),
    pct = round(100 * as.integer(counts) / sum(counts), 1)
  )
  df <- df[order(-df$n), ]
  colnames(df) <- c("Cluster", "n", "%")
  
  tbl_grob <- tableGrob(df, rows = NULL,
                        theme = ttheme_minimal(base_size = 9))
  
  ggsave(filename, tbl_grob,
         width = 3.5, height = 0.3 * nrow(df) + 1, limitsize = FALSE)
}

save_cluster_table(d2_pbmc_10x_CITE, "rna_only_clusters",
                   "figures/table_rna_only_clusters.png")
save_cluster_table(d2_pbmc_10x_CITE, "wnn_clusters",
                   "figures/table_wnn_clusters.png")

############### 2 tables merged into one, after transforming to Azimuth labels l2 #########

library(dplyr)

# for each clustering method, find each de novo cluster's dominant l2 label,
# then tag every cell in that cluster with that dominant label
label_by_dominant_l2 <- function(obj, cluster_col) {
  obj@meta.data %>%
    group_by(.data[[cluster_col]]) %>%
    mutate(dominant_l2 = names(sort(table(predicted.celltype.l2), decreasing = TRUE))[1]) %>%
    ungroup() %>%
    pull(dominant_l2)
}

d2_pbmc_10x_CITE$rna_dominant_l2 <- label_by_dominant_l2(d2_pbmc_10x_CITE, "rna_only_clusters")
d2_pbmc_10x_CITE$wnn_dominant_l2 <- label_by_dominant_l2(d2_pbmc_10x_CITE, "wnn_clusters")

# combined table: one row per l2 identity, counts from both clusterings
combined_df <- full_join(
  as.data.frame(table(d2_pbmc_10x_CITE$rna_dominant_l2)) %>% rename(celltype = Var1, n_rna = Freq),
  as.data.frame(table(d2_pbmc_10x_CITE$wnn_dominant_l2)) %>% rename(celltype = Var1, n_wnn = Freq),
  by = "celltype"
) %>%
  mutate(n_rna = ifelse(is.na(n_rna), 0, n_rna),
         n_wnn = ifelse(is.na(n_wnn), 0, n_wnn),
         pct_rna = round(100 * n_rna / sum(n_rna), 1),
         pct_wnn = round(100 * n_wnn / sum(n_wnn), 1)) %>%
  arrange(desc(n_rna))

library(gridExtra)
tbl_grob <- tableGrob(combined_df, rows = NULL, theme = ttheme_minimal(base_size = 8))
ggsave("figures/table_rna_vs_wnn_by_celltype.png", tbl_grob,
       width = 6, height = 0.3 * nrow(combined_df) + 1, limitsize = FALSE)


# ---- Save updated object ----
d2_pbmc_10x_CITE$seurat_clusters <- d2_pbmc_10x_CITE$rna_only_clusters
saveRDS(d2_pbmc_10x_CITE, "data/processed/d2_pbmc_10x_CITE_full_analysis.rds")

# ---- Notes ----
# - Dropout analysis: T/B cell identity comes from Azimuth (predicted.celltype.l1),
#   which is not circular with the RNA-count-based dropout test (uses whole
#   transcriptome + protein context, not just CD3 counts).
# - ARI (Adjusted Rand Index): 0 = random agreement, 1 = identical clusterings.
#   Given this dataset's TBNK panel is lineage-only (no activation/memory
#   markers), a small ARI shift is the expected, honest result - this
#   demonstrates the WNN method and its effect size on this data, not a claim
#   that ADT resolves subtypes RNA can't (that requires markers like
#   CD45RA/CCR7, absent from this panel).
