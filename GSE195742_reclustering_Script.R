setwd("~/Downloads/forMatt/")

set.seed(1234)
library(tidyverse)
library(Seurat)

dataSPF <- Read10X(data.dir = "SPF")
SPF <- CreateSeuratObject(counts = dataSPF, min.cells = 3, project = "SPF")

dataGF <- Read10X(data.dir = "GF")
GF <- CreateSeuratObject(counts = dataGF, min.cells = 3, project = "GF")

SPF[["percent.mt"]] <- PercentageFeatureSet(SPF, pattern = "^mt-")
GF[["percent.mt"]] <- PercentageFeatureSet(GF, pattern = "^mt-")

VlnPlot(SPF, features =  c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(GF, features =  c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

SPF <- subset(SPF, subset = nCount_RNA >= 1000 & nFeature_RNA >= 600 & percent.mt <= 25)
GF <- subset(GF, subset = nCount_RNA >= 1000 & nFeature_RNA >= 600 & percent.mt <= 25)

VlnPlot(SPF, features =  c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(GF, features =  c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

SPF <- NormalizeData(SPF, normalization.method = "LogNormalize", scale.factor = 10000)
GF <-  NormalizeData(GF, normalization.method = "LogNormalize", scale.factor = 10000)

SPF <- FindVariableFeatures(SPF, selection.method = "vst")
GF <- FindVariableFeatures(GF, selection.method = "vst")

VariableFeaturePlot(SPF)
VariableFeaturePlot(GF)

SPF <- ScaleData(SPF, features = rownames(SPF))
GF <- ScaleData(GF, features = rownames(SPF))

SPF <- RunPCA(SPF, features = VariableFeatures(object = SPF))
GF <- RunPCA(GF, features = VariableFeatures(object = GF))

feature_anchors <- SelectIntegrationFeatures(object.list = list(SPF, GF))
anchors <- FindIntegrationAnchors(object.list = list(SPF, GF),
                                  anchor.features = feature_anchors)

combined <- IntegrateData(anchorset = anchors)
DefaultAssay(combined) = "integrated"

combined <- ScaleData(combined)
combined <- RunPCA(combined)

ElbowPlot(combined)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:18)
combined <- FindClusters(combined, resolution = 0.2) # originally 0.3
combined <- RunUMAP(combined, reduction = "pca", dims = 1:18)
DimPlot(combined, reduction = "umap")
DimPlot(combined, reduction = "umap", group.by = "orig.ident")

saveRDS(combined, "Tsang_SPF_GF_11clusters_reduced.rds")
markers <- FindAllMarkers(combined, assay = "RNA", only.pos = TRUE)
top10 <- markers %>% group_by(cluster) %>% top_n(8, wt = avg_log2FC)

DefaultAssay(combined) <- "integrated"
DoHeatmap(combined, features = top10$gene) + NoLegend()

markers %>% group_by(cluster) %>% top_n(40, wt = avg_log2FC) %>%
  write_csv("Top50MarkerGenes.csv")

scTsang <- readRDS("Tsang_SPF_GF_11clusters.rds")
DimPlot(scTsang, reduction = "umap", label = TRUE)

ggsave("UMAP_Tsang_clusters.png", height = 3.5, width = 3.7)

DimPlot(scTsang, reduction = "umap", group.by = "orig.ident", cols = alpha(c("grey20",
                                                                       "seagreen3"), 0.4)) 
ggsave("UMAP_Tsang_condition.png", height = 3.5, width = 3.7)

for (i in 0:10L) {
  print(paste("Calculating SPF/GF for Cluster", i))
  FindMarkers(object = scTsang,
              assay = "RNA",
              ident.1 = "SPF",
              ident.2 = "GF",
              group.by = "orig.ident",
              subset.ident = i,
              test.use = "negbinom",
              min.pct = 0.2) %>%
    rownames_to_column(var = "gene") %>% mutate(Cluster = i) %>%
    write_csv(paste("negativeBinomial/SPFvGF_Cluster",
                    i, "FinalUMAP_NBinomial.csv", sep = "_"))
}

df_testing <- read_csv("NegativeBinomial/NegativeBinomial_SPFvGF.csv")

df_testing %>% mutate(signif = ifelse(p_val_adj < 0.05, "Yes", "No")) %>%
  ggplot(aes(
    x = avg_log2FC,
    color = signif
  )) +
  geom_freqpoly(binwidth = 0.01) +
  facet_wrap(~Cluster, ncol = 2) +
  theme_bw()
# ggsave("Distributions_NB_spfvgf.pdf", height = 12, width = 4)


df_testing %>% mutate(signif = ifelse(p_val_adj < 0.05, "Yes", "No")) %>%
  ggplot(aes(
    x = avg_log2FC,
    color = signif
  )) +
  geom_freqpoly(binwidth = 0.01) +
  theme_bw()