# 02_adt_normalization_and_viz.R
# PBMC-CITE-seq-ADT-azimuth-mapping
#
# Purpose:
#   1. Load Dataset 2 (filtered, Azimuth-annotated) from notebook 01
#   2. CLR normalize ADT assay
#   3. Fig 1: de novo RNA-based UMAP (no Azimuth) - A) clusters, B) ADT overlay, C) marker dotplot
#   4. Fig 2: Azimuth refUMAP - A) l1, B) l2, C) l3, D) multi-marker ADT overlay on l3
#   5. Fig 3: ADT + corresponding gene(s) expression across predicted.celltype.l2 clusters (n >= 60)

library(Seurat)
library(Azimuth)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

# ---- 1. Load annotated object from notebook 01 ----

d2_pbmc_10x_CITE <- readRDS("data/processed/d2_pbmc_10x_CITE_filtered_annotated.rds")

d2_pbmc_10x_CITE

# ---- 2. CLR normalize ADT ----

# CLR (centered log-ratio), margin = 2 normalizes across antibodies within
# each cell - standard approach for ADT data (per-cell technical variation
# correction; RNA-style log-normalization assumes a shared library-size
# scaling that doesn't hold for antibody counts).
d2_pbmc_10x_CITE <- NormalizeData(d2_pbmc_10x_CITE, assay = "ADT",
                                   normalization.method = "CLR", margin = 2)

# ---- 3. Fig 1: de novo RNA-based UMAP (no Azimuth) ----

DefaultAssay(d2_pbmc_10x_CITE) <- "RNA"

d2_pbmc_10x_CITE <- NormalizeData(d2_pbmc_10x_CITE, assay = "RNA")
d2_pbmc_10x_CITE <- FindVariableFeatures(d2_pbmc_10x_CITE, assay = "RNA")
d2_pbmc_10x_CITE <- ScaleData(d2_pbmc_10x_CITE, assay = "RNA")
d2_pbmc_10x_CITE <- RunPCA(d2_pbmc_10x_CITE, assay = "RNA", reduction.name = "pca")
d2_pbmc_10x_CITE <- FindNeighbors(d2_pbmc_10x_CITE, reduction = "pca", dims = 1:30)
d2_pbmc_10x_CITE <- FindClusters(d2_pbmc_10x_CITE, resolution = 0.8)
d2_pbmc_10x_CITE <- RunUMAP(d2_pbmc_10x_CITE, reduction = "pca", dims = 1:30,
                             reduction.name = "umap")

# helper: scale a vector to [0,1] for RGB blending
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

# 1A - de novo clusters, unlabeled, square aspect
fig1a <- DimPlot(d2_pbmc_10x_CITE, reduction = "umap", group.by = "seurat_clusters",
                  label = TRUE) +
  labs(title = "De novo RNA clusters") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size = 10, face = "plain")
  )

# 1B - single UMAP, 3-marker ADT color blend overlay (RGB blend caps at 3
# channels for a clean, interpretable overlay - more markers would require
# either dropping the "single plot" constraint or accepting a muddier blend,
# so 3 is used here: CD14 (red), CD3 (blue), CD19 (green))
adt_clr_denovo <- GetAssayData(d2_pbmc_10x_CITE, assay = "ADT", layer = "data")

r_chan_1b <- scale01(adt_clr_denovo["CD14.1", ])
g_chan_1b <- scale01(adt_clr_denovo["CD19.1", ])
b_chan_1b <- scale01(adt_clr_denovo["CD3", ])

umap_coords_1b <- as.data.frame(Embeddings(d2_pbmc_10x_CITE, reduction = "umap"))
umap_coords_1b$blend <- rgb(r_chan_1b, g_chan_1b, b_chan_1b)

fig1b <- ggplot(umap_coords_1b, aes(x = umap_1, y = umap_2)) +
  geom_point(color = umap_coords_1b$blend, size = 0.6) +
  theme_classic() +
  labs(title = "ADT overlay: CD14 (red) / CD19 (green) / CD3 (blue)",
       x = "umap_1", y = "umap_2") +
  theme(plot.title = element_text(size = 11, face = "plain"))

# 1C - top 2 marker genes per de novo cluster, dotplot
# Idents(d2_pbmc_10x_CITE) <- "seurat_clusters"
# de_novo_markers <- FindAllMarkers(d2_pbmc_10x_CITE, assay = "RNA", only.pos = TRUE,
                                   # min.pct = 0.25, logfc.threshold = 0.5)

top2_markers <- de_novo_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 2) %>%
  pull(gene) %>%
  unique()

fig1c <- DotPlot(d2_pbmc_10x_CITE, features = top2_markers, group.by = "seurat_clusters",
                 dot.scale=4) +
  RotatedAxis() +
  labs(title = "Top 2 marker genes per de novo cluster",
       x = "Marker Gene",
       y = "Cluster") +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 4),
    plot.title = element_text(size = 11, face = "plain")
  )

# combine as one row: A (2.5) | B (2.5) | C (5), total width 10
fig1 <- fig1a + fig1b + fig1c + plot_layout(widths = c(2, 2, 6))
fig1

ggsave("figures/02_fig1_denovo_overview.png", fig1, width = 20, height = 4)


# ---- 4. Fig 2: Azimuth ref.umap (2x2 combined figure) ----

# helper: scale a vector to [0,1] for RGB blending
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

# helper: UMAP plot colored by group, with labels at each group's centroid,
# colored a darkened version of that group's own color
plot_labeled_umap <- function(embeddings, group_vals, title, show_legend = TRUE,
                              label_size = 3) {
  dim1_name <- colnames(embeddings)[1]
  dim2_name <- colnames(embeddings)[2]
  
  df <- as.data.frame(embeddings)
  df$group <- factor(group_vals)
  
  base_colors <- scales::hue_pal()(nlevels(df$group))
  names(base_colors) <- levels(df$group)
  dark_colors <- colorspace::darken(base_colors, amount = 0.35)
  names(dark_colors) <- levels(df$group)
  
  centroids <- df %>%
    group_by(group) %>%
    summarise(x = median(.data[[dim1_name]]), y = median(.data[[dim2_name]])) %>%
    mutate(label_color = dark_colors[as.character(group)],
           fill_color = colorspace::lighten(base_colors[as.character(group)], amount = 0.82))
  
  ggplot(df, aes(x = .data[[dim1_name]], y = .data[[dim2_name]], color = group)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(values = base_colors) +
    ggrepel::geom_label_repel(data = centroids, aes(x = x, y = y, label = group),
                              color = centroids$label_color, size = label_size,
                              fontface = "plain", inherit.aes = FALSE,
                              fill = centroids$fill_color, label.size = NA,
                              label.padding = unit(0.1, "lines")) +
    theme_classic() +
    theme(plot.title = element_text(size = 11, face = "plain"),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8),
          legend.position = if (show_legend) "right" else "none") +
    labs(title = title, x = dim1_name, y = dim2_name, color = NULL)
}

ref_umap_embed <- Embeddings(d2_pbmc_10x_CITE, reduction = "ref.umap")

fig2a <- plot_labeled_umap(ref_umap_embed, d2_pbmc_10x_CITE$predicted.celltype.l1,
                           "Azimuth ref.umap - predicted.celltype.l1",
                           show_legend = FALSE)

fig2b <- plot_labeled_umap(ref_umap_embed, d2_pbmc_10x_CITE$predicted.celltype.l2,
                           "Azimuth ref.umap - predicted.celltype.l2",
                           show_legend = FALSE)

fig2c <- plot_labeled_umap(ref_umap_embed, d2_pbmc_10x_CITE$predicted.celltype.l3,
                           "Azimuth ref.umap - predicted.celltype.l3",
                           show_legend = FALSE, label_size = 3)

# 2D - multi-marker ADT overlay (CD14 red, CD19 green, CD3 blue) on ref.umap,
# with light l3-cluster-colored background points underneath
adt_clr <- GetAssayData(d2_pbmc_10x_CITE, assay = "ADT", layer = "data")

r_chan <- scale01(adt_clr["CD14.1", ])
g_chan <- scale01(adt_clr["CD19.1", ])
b_chan <- scale01(adt_clr["CD3", ])

blend_colors <- rgb(r_chan, g_chan, b_chan)

umap_coords <- as.data.frame(ref_umap_embed)
dim1_name <- colnames(umap_coords)[1]
dim2_name <- colnames(umap_coords)[2]
umap_coords$blend <- blend_colors
umap_coords$l3 <- d2_pbmc_10x_CITE$predicted.celltype.l3

fig2d <- ggplot(umap_coords, aes(x = .data[[dim1_name]], y = .data[[dim2_name]])) +
  geom_point(aes(color = l3), size = 1.2, alpha = 0.15) +
  geom_point(color = umap_coords$blend, size = 0.6) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size = 11, face = "plain"),
        plot.subtitle = element_text(size = 8)) +
  labs(title = "ADT overlay: CD14 (red) / CD19 (green) / CD3 (blue)",
       subtitle = "background: light color per predicted.celltype.l3",
       x = dim1_name, y = dim2_name)

# combine as 2x2
fig2 <- (fig2a + fig2b) / (fig2c + fig2d)
fig2

ggsave("figures/02_fig2_refumap_2x2.png", fig2, width = 14, height = 12)

############## cell count tables to go with it

library(gridExtra)
library(grid)

# helper: build count + % table for any predicted.celltype.* column, save as png
save_celltype_table <- function(obj, column, filename, title) {
  counts <- table(obj@meta.data[[column]])
  df <- data.frame(
    celltype = names(counts),
    n = as.integer(counts),
    pct = round(100 * as.integer(counts) / sum(counts), 1)
  )
  df <- df[order(-df$n), ]  # sort descending by count
  colnames(df) <- c("Cell type", "n", "%")
  
  tbl_grob <- tableGrob(df, rows = NULL,
                        theme = ttheme_minimal(base_size = 9))
  
  ggsave(filename, tbl_grob,
         width = 4, height = 0.3 * nrow(df) + 1, limitsize = FALSE, bg = "white")
}

save_celltype_table(d2_pbmc_10x_CITE, "predicted.celltype.l1",
                    "figures/table_celltype_l1.png", "predicted.celltype.l1")
save_celltype_table(d2_pbmc_10x_CITE, "predicted.celltype.l2",
                    "figures/table_celltype_l2.png", "predicted.celltype.l2")
save_celltype_table(d2_pbmc_10x_CITE, "predicted.celltype.l3",
                    "figures/table_celltype_l3.png", "predicted.celltype.l3")


# combined cell count table

hierarchy_df <- d2_pbmc_10x_CITE@meta.data %>%
  count(predicted.celltype.l1, predicted.celltype.l2, predicted.celltype.l3) %>%
  arrange(predicted.celltype.l1, predicted.celltype.l2, desc(n)) %>%
  mutate(pct = round(100 * n / sum(n), 1))

# blank out repeated l1/l2 values so only the first row of each group shows the label
hierarchy_df <- hierarchy_df %>%
  mutate(
    l1_display = ifelse(duplicated(predicted.celltype.l1), "", as.character(predicted.celltype.l1)),
    l2_display = ifelse(duplicated(paste(predicted.celltype.l1, predicted.celltype.l2)), "",
                        as.character(predicted.celltype.l2))
  ) %>%
  select(l1_display, l2_display, predicted.celltype.l3, n, pct) %>%
  rename(l1 = l1_display, l2 = l2_display, l3 = predicted.celltype.l3, `%` = pct)

tbl_grob <- tableGrob(hierarchy_df, rows = NULL,
                      theme = ttheme_minimal(base_size = 8))

ggsave("figures/table_celltype_hierarchy.png", tbl_grob,
       width = 6, height = 0.25 * nrow(hierarchy_df) + 1, limitsize = FALSE, bg = "white")


# ---- 5. Fig 3: ADT vs. RNA agreement across predicted.celltype.l2 ----

# ADT marker -> actual ADT assay rowname (some got ".1" suffixes from Seurat's
# auto-renaming, since their gene symbol clashed with an RNA feature) and
# corresponding RNA gene(s). CD3/CD8 encode multi-gene protein complexes;
# CD16/CD56 have different gene symbols entirely. Z-scoring within modality
# lets us compare pattern/shape across clusters, not absolute magnitude -
# ADT (CLR) and RNA (log-normalized) are not on the same raw scale.
adt_gene_map <- list(
  CD3   = list(adt_row = "CD3",    genes = c("CD3D", "CD3E", "CD3G")),
  CD4   = list(adt_row = "CD4.1",  genes = c("CD4")),
  CD8   = list(adt_row = "CD8",    genes = c("CD8A", "CD8B")),
  CD11c = list(adt_row = "CD11c",  genes = c("ITGAX")),
  CD14  = list(adt_row = "CD14.1", genes = c("CD14")),
  CD16  = list(adt_row = "CD16",   genes = c("FCGR3A")),
  CD19  = list(adt_row = "CD19.1", genes = c("CD19")),
  CD56  = list(adt_row = "CD56",   genes = c("NCAM1")),
  CD45  = list(adt_row = "CD45",   genes = c("PTPRC"))
)

# filter l2 clusters to n >= 60 cells for this comparison
l2_counts <- table(d2_pbmc_10x_CITE$predicted.celltype.l2)
l2_keep <- names(l2_counts[l2_counts >= 60])
l2_dropped <- names(l2_counts[l2_counts < 60])

message("predicted.celltype.l2 clusters excluded (n < 60): ",
        paste(l2_dropped, collapse = ", "))

cells_keep <- colnames(d2_pbmc_10x_CITE)[d2_pbmc_10x_CITE$predicted.celltype.l2 %in% l2_keep]
cluster_labels <- d2_pbmc_10x_CITE$predicted.celltype.l2[cells_keep]

# helper: z-score features across kept cells, average multiple genes into
# one composite score, then mean per predicted.celltype.l2 cluster
mean_z_by_cluster <- function(obj, features, assay, layer, cells, clusters) {
  mat <- GetAssayData(obj, assay = assay, layer = layer)[features, cells, drop = FALSE]
  mat_z <- t(scale(t(as.matrix(mat))))
  composite <- if (length(features) > 1) colMeans(mat_z) else mat_z[1, ]
  tapply(composite, clusters, mean)
}

adt_matrix <- matrix(NA, nrow = length(adt_gene_map), ncol = length(l2_keep),
                     dimnames = list(names(adt_gene_map), l2_keep))
rna_matrix <- matrix(NA, nrow = length(adt_gene_map), ncol = length(l2_keep),
                     dimnames = list(names(adt_gene_map), l2_keep))

for (adt_name in names(adt_gene_map)) {
  adt_row <- adt_gene_map[[adt_name]]$adt_row
  genes <- adt_gene_map[[adt_name]]$genes
  genes_present <- intersect(genes, rownames(d2_pbmc_10x_CITE[["RNA"]]))
  
  adt_means <- mean_z_by_cluster(d2_pbmc_10x_CITE, adt_row, "ADT", "data",
                                 cells_keep, cluster_labels)
  adt_matrix[adt_name, names(adt_means)] <- adt_means
  
  if (length(genes_present) > 0) {
    rna_means <- mean_z_by_cluster(d2_pbmc_10x_CITE, genes_present, "RNA", "data",
                                   cells_keep, cluster_labels)
    rna_matrix[adt_name, names(rna_means)] <- rna_means
  }
}

# consistent ordering across all panels
cluster_order <- colnames(adt_matrix)[hclust(dist(t(adt_matrix)))$order]
marker_order  <- rownames(adt_matrix)[hclust(dist(adt_matrix))$order]

to_long <- function(mat, modality) {
  as.data.frame(mat) %>%
    tibble::rownames_to_column("marker") %>%
    pivot_longer(-marker, names_to = "cluster", values_to = "mean_z") %>%
    mutate(modality = modality,
           marker = factor(marker, levels = marker_order),
           cluster = factor(cluster, levels = cluster_order))
}

heat_long <- bind_rows(to_long(adt_matrix, "ADT (protein)"),
                       to_long(rna_matrix, "RNA (gene)"))

# ---- Panel A: heatmap ----
fig3a <- ggplot(heat_long, aes(x = cluster, y = marker, fill = mean_z)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, name = "mean z-score") +
  facet_wrap(~modality, ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 11, face = "plain"),
        strip.text = element_text(size = 9)) +
  labs(title = "ADT vs. RNA expression across predicted.celltype.l2",
       x = NULL, y = NULL)

# ---- Panel B: correlation across all marker/cluster pairs ----
scatter_df <- heat_long %>%
  select(marker, cluster, modality, mean_z) %>%
  pivot_wider(names_from = modality, values_from = mean_z) %>%
  rename(adt_z = `ADT (protein)`, rna_z = `RNA (gene)`) %>%
  filter(!is.na(adt_z), !is.na(rna_z))

cor_val <- cor(scatter_df$adt_z, scatter_df$rna_z, use = "complete.obs")

fig3b <- ggplot(scatter_df, aes(x = rna_z, y = adt_z, color = marker)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 1.8, alpha = 0.7) +
  annotate("text", x = min(scatter_df$rna_z), y = max(scatter_df$adt_z),
           label = paste0("r = ", round(cor_val, 2)), hjust = 0, size = 3.2) +
  theme_classic() +
  theme(plot.title = element_text(size = 11, face = "plain"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8)) +
  labs(title = "ADT vs. RNA agreement across all marker/cluster pairs",
       x = "mean RNA z-score", y = "mean ADT z-score", color = "Marker")

# ---- Panel C: CD8 detail, ADT + RNA combined (z-scored) on one axis ----
marker_focus <- "CD8"
adt_row_focus <- adt_gene_map[[marker_focus]]$adt_row
genes_focus <- intersect(adt_gene_map[[marker_focus]]$genes, rownames(d2_pbmc_10x_CITE[["RNA"]]))

adt_raw <- GetAssayData(d2_pbmc_10x_CITE, assay = "ADT", layer = "data")[adt_row_focus, cells_keep]
rna_raw <- colMeans(GetAssayData(d2_pbmc_10x_CITE, assay = "RNA", layer = "data")[genes_focus, cells_keep, drop = FALSE])

detail_df <- data.frame(
  cluster = factor(cluster_labels, levels = cluster_order),
  ADT = as.numeric(scale(adt_raw)),
  RNA = as.numeric(scale(rna_raw))
) %>%
  pivot_longer(cols = c(ADT, RNA), names_to = "modality", values_to = "z_expr")

modality_fill  <- c(ADT = "#D9A5A5", RNA = "lightblue")
modality_point <- c(ADT = "#7A1F2B", RNA = "darkblue")

fig3c <- ggplot(detail_df, aes(x = cluster, y = z_expr, fill = modality)) +
  geom_violin(position = position_dodge(width = 0.8), color = NA, alpha = 0.8) +
  geom_point(aes(color = modality), position = position_jitterdodge(jitter.width = 0.1,
                                                                    dodge.width = 0.8), size = 0.3, alpha = 0.3) +
  scale_fill_manual(values = modality_fill) +
  scale_color_manual(values = modality_point) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        plot.title = element_text(size = 11, face = "plain"),
        legend.title = element_blank()) +
  labs(title = paste0(marker_focus, ": ADT vs. RNA (z-scored) across predicted.celltype.l2"),
       x = NULL, y = "z-scored expression")

# ---- Combine A / B / C, tagged in top-left corner ----
fig3_final <- (fig3a + fig3b + fig3c) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold"),
        plot.tag.position = c(0, 1))

fig3_final

ggsave("figures/02_fig3_full.png", fig3_final, width = 15, height = 5)

# ---- Save updated object ----
saveRDS(d2_pbmc_10x_CITE, "data/processed/d2_pbmc_10x_CITE_full_analysis.rds")

# ---- Notes ----
# - Fig 3A averages z-scores across a marker's constituent RNA genes into one
#   composite score per marker; row/column order set by hierarchical
#   clustering of the ADT matrix, applied identically to both heatmap panels.
# - Fig 3B pools all marker/cluster pairs into one correlation plot; dashed
#   line marks perfect ADT/RNA agreement (y = x).
# - Fig 3C shows CD8 at single-cell resolution rather than cluster-mean
#   summary, both modalities z-scored independently onto one shared axis.
# - predicted.celltype.l2 clusters with n < 60 cells excluded from Fig 3 only
#   (still present in Fig 2 UMAPs) - see console message for excluded clusters.


