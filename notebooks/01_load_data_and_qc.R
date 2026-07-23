# 01_load_data_and_qc.R
# PBMC-CITE-seq-ADT-azimuth-mapping
#
# Purpose:
#   1. Load Dataset 2 (10x PBMC CITE-seq demo, query) from data/raw/
#   2. Basic QC on RNA + ADT
#   3. Map onto Dataset 1 (Hao et al. 2021 PBMC reference) via Azimuth
#   4. Save annotated object to data/processed/

library(Seurat)
library(Azimuth)
library(dplyr)
library(ggplot2)
library(patchwork)

# ---- 1. Load Dataset 2: 10x PBMC CITE-seq demo (query) ----

# Cell Ranger multi output contains two feature types:
# Gene Expression + Antibody Capture (ADT)
d2_pbmc_10x_CITE_raw <- Read10X_h5("data/raw/filtered_feature_bc_matrix.h5")

# confirm both Gene Expression (scRNA-seq) and Antibody Capture (ADT) are present
# in Cell Ranger multi output file
names(d2_pbmc_10x_CITE_raw)  # confirm both feature types present

# laod raw data into Seurat object (doesn't contain ADT yet)
d2_pbmc_10x_CITE <- CreateSeuratObject(counts = d2_pbmc_10x_CITE_raw$`Gene Expression`,
                                        project = "pbmc_10x_demo")

colnames(d2_pbmc_10x_CITE@meta.data)

# add the ADT data to the Seurat obj (as a new assay)
# ADT was seq'd and included in raw data: filtered_feature_bc_matrix.h5
# cell barcode order is shared between GEX and ADT matrices, guaranteed by Cell Ranger multi
d2_pbmc_10x_CITE[["ADT"]] <- CreateAssayObject(counts = d2_pbmc_10x_CITE_raw$`Antibody Capture`)

# Sanity checks
d2_pbmc_10x_CITE
rownames(d2_pbmc_10x_CITE[["ADT"]])  # confirm ADT panel names look right
dim(d2_pbmc_10x_CITE[["RNA"]])
dim(d2_pbmc_10x_CITE[["ADT"]])

# ---- 2. QC metrics ----

d2_pbmc_10x_CITE[["percent.mt"]] <- PercentageFeatureSet(d2_pbmc_10x_CITE, pattern = "^MT-")
d2_pbmc_10x_CITE[["nCount_ADT"]] <- colSums(GetAssayData(d2_pbmc_10x_CITE, assay = "ADT", layer = "counts"))

qc_long <- d2_pbmc_10x_CITE@meta.data %>%
  select(nFeature_RNA, nCount_RNA, percent.mt, nCount_ADT) %>%
  tidyr::pivot_longer(everything(), names_to = "metric", values_to = "value")

violin_colors <- c(
  nFeature_RNA = "lightblue",
  nCount_RNA   = "lightblue",
  percent.mt   = "lightblue",
  nCount_ADT   = "#D9A5A5"   # light maroon
)

point_colors <- c(
  nFeature_RNA = "darkblue",
  nCount_RNA   = "darkblue",
  percent.mt   = "darkblue",
  nCount_ADT   = "#7A1F2B"   # maroon
)

ggplot(qc_long, aes(x = metric, y = value)) +
  geom_violin(aes(fill = metric), color = NA) +
  geom_jitter(aes(color = metric), size = 0.3, alpha = 0.4, width = 0.15) +
  scale_fill_manual(values = violin_colors) +
  scale_color_manual(values = point_colors) +
  facet_wrap(~metric, scales = "free", ncol = 4) +
  labs(title = "QC Metrics — Dataset 2 (10x PBMC CITE-seq) unfiltered", x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

ggsave("figures/01_ds2_unfiltered_qc_violins.png", width = 12, height = 4)

# save the unfiltered dataset
saveRDS(d2_pbmc_10x_CITE, "data/processed/d2_pbmc_10x_CITE_unfiltered.rds")

# ---- 3. Filtering Dataset 2 ----

# get exact quantiles for where you start losing a lot of data
summary(d2_pbmc_10x_CITE$percent.mt)
quantile(d2_pbmc_10x_CITE$percent.mt, probs = c(0.5, 0.75, 0.9, 0.95, 0.99))
# median is 4.9, can afford to cut off around 15-20% to control for dying cells

summary(d2_pbmc_10x_CITE$nFeature_RNA)
quantile(d2_pbmc_10x_CITE$nFeature_RNA, probs = c(0.01, 0.05, 0.1, 0.25, 0.5))
# median ~2k is good but 1st percentile is only 193 --> a >200 cutoff removes almost nothing
# > 500 cuts closer to 7th percentile, a more standard "real cell" floor

summary(d2_pbmc_10x_CITE$nCount_RNA)
quantile(d2_pbmc_10x_CITE$nCount_RNA, probs = c(0.01, 0.05, 0.1, 0.25, 0.5))
# median 5k but 5th percentile only 1k.  >1000 cutoff floor removes bottom 5%.

sapply(c(10, 15, 20, 25), function(x) sum(d2_pbmc_10x_CITE$percent.mt < x))
sapply(c(200, 500, 1000), function(x) sum(d2_pbmc_10x_CITE$nFeature_RNA > x))
sapply(c(500, 1000, 1500), function(x) sum(d2_pbmc_10x_CITE$nCount_RNA > x))

# total cells before filtering, for reference: 11075
nrow(d2_pbmc_10x_CITE@meta.data)

# overwrites the unfiltered dataset (which was saved as RDS)
d2_pbmc_10x_CITE <- subset(d2_pbmc_10x_CITE,
                           subset = nFeature_RNA > 500 &
                             nCount_RNA > 1000 &
                             percent.mt < 20)

d2_pbmc_10x_CITE  # confirm cell count after filtering

# make new qc vlns fig
qc_long <- d2_pbmc_10x_CITE@meta.data %>%
  select(nFeature_RNA, nCount_RNA, percent.mt, nCount_ADT) %>%
  tidyr::pivot_longer(everything(), names_to = "metric", values_to = "value")

violin_colors <- c(
  nFeature_RNA = "lightblue",
  nCount_RNA   = "lightblue",
  percent.mt   = "lightblue",
  nCount_ADT   = "#D9A5A5"   # light maroon
)

point_colors <- c(
  nFeature_RNA = "darkblue",
  nCount_RNA   = "darkblue",
  percent.mt   = "darkblue",
  nCount_ADT   = "#7A1F2B"   # maroon
)

ggplot(qc_long, aes(x = metric, y = value)) +
  geom_violin(aes(fill = metric), color = NA) +
  geom_jitter(aes(color = metric), size = 0.3, alpha = 0.4, width = 0.15) +
  scale_fill_manual(values = violin_colors) +
  scale_color_manual(values = point_colors) +
  facet_wrap(~metric, scales = "free", ncol = 4) +
  labs(title = "QC Metrics — Dataset 2 (10x PBMC CITE-seq) filtered", x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

ggsave("figures/01_ds2_filtered_qc_violins.png", width = 12, height = 4)

# ---- 4. Azimuth - Map ds2 onto ds1 (Hao et al. reference) ----

# Downloads and caches the Hao et al. 2021 PBMC reference on first run
d2_pbmc_10x_CITE <- RunAzimuth(d2_pbmc_10x_CITE, reference = "pbmcref")

# Confirm predicted labels were added
colnames(d2_pbmc_10x_CITE@meta.data)

# Azimuth: predicted.celltype.l1/l2/l3
# l1 — broadest category: i.e. "T cell"
# l2 — intermediate res: i.e. "CD4 T", "CD8 T", "Treg", "CD14 Mono", "CD16 Mono".
# l3 — finest res: i.e. CD4 T --> "CD4 Naive", "CD4 TCM", "CD4 TEM"
table(d2_pbmc_10x_CITE$predicted.celltype.l1)

# Save annotated object
saveRDS(d2_pbmc_10x_CITE, "data/processed/d2_pbmc_10x_CITE_filtered_annotated.rds")

# ---- Notes ----
# - ADT normalization (CLR) not yet applied - handle in next notebook
# - nCount_ADT included in QC as a sanity check on ADT depth, not used in filtering here
