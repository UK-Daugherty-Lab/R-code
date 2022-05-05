#### This R code was used for the publication below
#### https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.121.058173

#Set library for scRNAseq data analysis
library(SeuratObject)
library(Seurat)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(cowplot)
library(zinbwave)
library(scRNAseq)
library(matrixStats)
library(edgeR)
library(stringr)
library(monocle)
library(Signac)
library(data.table)
library(devtools)
library(patchwork)
library(BiocParallel)
library(EnsDb.Mmusculus.v79)
library(org.Mm.eg.db)
library(clusterProfiler)
library(monocle)
library(igraph)
library(ensembldb)
library(MAST)
library(NMF)
library(rsvd)
library(GGally)

## creat function named set_up (used below)
set_up = function(UMI, project, stim){
  data = read.csv(UMI, header=T, row.names = "X", stringsAsFactors = F)
  each = CreateSeuratObject(counts = data, project = project, min.cells = 5)
  each$stim = stim
  each = NormalizeData(each, verbose = FALSE)
  each = FindVariableFeatures(each, selection.method = "vst", nfeatures = 2000)
}

# 1. Create a Seurat object for control --------------------------------------
## Read a CSV file for read counts of baseline-SHF cells
Ctrl = fread("SHF_baseline.csv", header = T, stringsAsFactors = F)

## Remove duplicate
Ctrl = Ctrl[!duplicated(Ctrl$V1),]
rown = Ctrl$V1
Ctrl = as.matrix(Ctrl[,-1])
rownames(Ctrl) = rown

## Create a Seurat object
Ctrl = Seurat::CreateSeuratObject(counts = Ctrl, project = "Ctrl", min.cells = 3, min.features = 200)
Ctrl$stim = "Ctrl"
Ctrl[["percent.mt"]] = PercentageFeatureSet(Ctrl, pattern = "^mt-")

## Quality check
VlnPlot(Ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 = FeatureScatter(Ctrl, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(Ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

## Exclude non-cell and aggregated cells
Ctrl = subset(Ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

## Normalization
Ctrl = NormalizeData(Ctrl, verbose = FALSE)
Ctrl = FindVariableFeatures(Ctrl, selection.method = "vst", nfeatures = 2000)

# 2. Create a Seurat object for AngII ----------------------------------------
## Read a CSV file for read counts of AngII-SHF cells
AngII = fread("SHF_AngII.csv", header = T, stringsAsFactors = F)

## Remove duplicate
AngII = AngII[!duplicated(AngII$V1),]
rown = AngII$V1
AngII = as.matrix(AngII[,-1])
rownames(AngII) = rown

## Create a Seurat object
AngII = CreateSeuratObject(counts = AngII, project = "AngII", min.cells = 3, min.features = 200)
AngII$stim = "AngII"
AngII[["percent.mt"]] = PercentageFeatureSet(AngII, pattern = "^mt-")

## Quality check
VlnPlot(AngII, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot3 = FeatureScatter(AngII, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 = FeatureScatter(AngII, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot3, plot4))

## Exclude non-cell and aggregated cells
AngII = subset(AngII, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

## Normalization
AngII = NormalizeData(AngII, verbose = FALSE)
AngII = FindVariableFeatures(AngII, selection.method = "vst", nfeatures = 2000)

# 3. Integrate two Seurat objects --------------------------------------------
features = SelectIntegrationFeatures(object.list = list(Ctrl, AngII))
anchors = FindIntegrationAnchors(object.list = list(Ctrl, AngII), dims = 1:20, anchor.features = features)
combined = IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(combined) = "integrated"

## set "Ctrl" as 0 and "AngII as 1
combined@meta.data$stim = factor(combined@meta.data$stim, c("Ctrl", "AngII")) 

# 4. Clustering --------------------------------------------------------------
combined = ScaleData(combined, verbose = FALSE)
combined = RunPCA(combined, npcs = 30, verbose = FALSE)
combined = RunUMAP(combined, reduction = "pca", dims = 1:20)
combined = FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined = FindClusters(combined, resolution = 0.5)

p1 = DimPlot(combined, reduction = "umap", group.by = "stim", pt.size = 0.4)
p2 = DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
ggsave("umap.png", p1 + p2, height=5, width=10)

DimPlot(combined, reduction = "umap", split.by = "stim", pt.size = 0.4)
ggsave("umap_stim.png", height=5, width=10)

# 5. Name clusters -----------------------------------------------------------
## Find markers in each cluster
markers = FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "marker genes - clusters.csv", quote = F)
top5 = markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 5, wt = avg_log2FC)
DoHeatmap(combined, features = top5$gene, size = 7) + theme(text = element_text(size = 18))
ggsave("heatmap - marker genes - clusters.png", width = 40, height = 15)

# --- Supplemental Figure 3A --- #
FeaturePlot(combined, features = c("Acta2", "Col1a1", "Cd68", "Pecam1", "Ly6a", "Myh11", "Col1a2", "Lyz2", "Vwf"),  min.cutoff = "q9", cols=c("lightgrey", "goldenrod","firebrick"), pt.size=0.2, ncol = 3)
ggsave("Sup_Fig_3A.png", width = 12, height = 11)

## cluster similarities (correlation)
combined = FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
varGene = VariableFeatures(combined)
avg = AverageExpression(combined)
avg = avg$RNA
data = avg[row.names(avg) %in% varGene,]
res = cor(data, method = "spearman")
png("heatmap - cluster similarity.png")
Heatmap_CorrelationClusters = pheatmap(res, main = "",
                                       cellwidth = 30,
                                       cellheight = 30,
                                       fontsize = 20,
                                       angle_col = 0)
ggsave("heatmap - cluster similarity.png", plot = Heatmap_CorrelationClusters, width=15, height=15)

## Conserved Markers in each cluster 
for(i in 1:15) {
  conserved_markers = FindConservedMarkers(combined, ident.1 = i, grouping.var = c("stim"), print.bar = FALSE)
  write.csv(conserved_markers, paste("comserved_markers_c", i, ".csv", sep = ""), quote = F)
}

## rename clusters
Idents(combined) = combined@meta.data$seurat_clusters

combined_sub = RenameIdents(combined, `0` = "SMC", `1` = "SMC", `2` = "FB2", `3` = "SMC", `4` = "FB2", `5` = "FB1", `6` = "FB1", `7` = "FB2", `8` = "FB2", `9` = "FB4", `10` = "FB3", `11` = "SMC", `12` = "SMC", `13` = "EC", `14` = "UI", `15` = "Mac")
combined = RenameIdents(combined, `0` = "SMC", `1` = "SMC", `2` = "FB", `3` = "SMC", `4` = "FB", `5` = "FB", `6` = "FB", `7` = "FB", `8` = "FB", `9` = "FB", `10` = "FB", `11` = "SMC", `12` = "SMC", `13` = "EC", `14` = "UI", `15` = "Mac")

# --- Figure 4C --- #
Idents(combined) = factor(Idents(combined), c("SMC", "FB", "EC", "Mac", "UI"))
u_col = c("#F8766D", "#00BFC4", "#7CAE00", "#619CFF", "#C77CFF")

DimPlot(combined_sub, reduction = "umap", pt.size = 0.4)
ggsave("umap_after renaming_sub.png", width = 4.5, height = 4)
DimPlot(combined_sub, reduction = "umap", pt.size = 0.4, split.by = "stim")
ggsave("umap_after renaming_sub stim.png", width = 8, height = 4)

DimPlot(combined, reduction = "umap", pt.size = 0.4, cols = u_col)
ggsave("umap_after renaming.png", width = 4.5, height = 4)
DimPlot(combined, reduction = "umap", pt.size = 0.4, split.by = "stim", cols = u_col)
ggsave("Fig_4C.png", width = 8, height = 4)

# --- Supplemental Figure 3B --- #
Idents(combined) = factor(Idents(combined), c("SMC", "FB", "EC", "UI", "Mac"))
CellMarkers = FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CellMarkers, "marker genes - clusters after renaming.csv", quote = F)

top5CellMarkers = CellMarkers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DotPlot(combined, features = unique(top5CellMarkers$gene), cols = c("gray95","royalblue4")) + RotatedAxis() 
ggsave("Sup_Fig_3B.png", width = 10, height = 2.5)

# --- Subclustering data --- #
Idents(combined_sub) = factor(Idents(combined_sub), c("SMC", "FB1", "FB2", "FB3", "FB4", "EC", "UI", "Mac"))
CellMarkers = FindAllMarkers(combined_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CellMarkers, "marker genes - clusters after renaming_sub.csv", quote = F)

# --- Figure 4D --- #
celltype = Idents(combined)
stim = combined@meta.data$stim
dat = data.frame(celltype, stim)

pie.Ctrl.dat = data.frame(with(dat %>% filter(stim == "Ctrl"), table(celltype, stim)))
pie.Ctrl.dat = filter (pie.Ctrl.dat, stim == "Ctrl")
pie.Ctrl.dat$celltype = factor(pie.Ctrl.dat$celltype,  c("SMC", "FB", "EC", "UI", "Mac"))

pie.AngII.dat = data.frame(with(dat %>% filter(stim == "AngII"), table(celltype, stim)))
pie.AngII.dat = filter (pie.AngII.dat, stim == "AngII")
pie.AngII.dat$celltype = factor(pie.AngII.dat$celltype,  c("SMC", "FB", "EC", "UI", "Mac"))

u_col = c("#7CAE00", "#7CAE00", "#7CAE00","#00BFC4", "#F8766D")
pie.Ctrl.dat$celltype = factor(pie.Ctrl.dat$celltype, rev(as.character(unique(pie.Ctrl.dat$celltype))))
p1 = ggplot(data = pie.Ctrl.dat, aes(x = "", y = Freq, fill = celltype)) +
  geom_bar(width = 1, size = 1, stat = "identity") +
  coord_polar("y") +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(name = "Ctrl", values = u_col) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

pie.AngII.dat$celltype = factor(pie.AngII.dat$celltype, rev(as.character(unique(pie.AngII.dat$celltype))))
p2 = ggplot(data = pie.AngII.dat, aes(x = "", y = Freq, fill = celltype)) +
  geom_bar(width = 1, size = 1, stat = "identity") +
  coord_polar("y") +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(name = "Ctrl", values = u_col) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
plot_grid(p1, p2)
ggsave("Fig_4D.png", width = 10, height = 5)

## chisq.test - Ctrl vs AngII in SMC/FB 
(tbl = with(dat, table(celltype, stim)))
chisq.test(tbl)
(tbl = with(dat, table(factor(celltype == "SMC"), stim)))
chisq.test(tbl)
(tbl = with(dat, table(factor(celltype == "FB"), stim)))
chisq.test(tbl)
(tbl = with(dat, table(factor(celltype == "EC"), stim)))
chisq.test(tbl)
(tbl = with(dat, table(factor(celltype == "Mac"), stim)))
chisq.test(tbl)

saveRDS(combined, "combined.rds")
saveRDS(combined, "combined_sub.rds")

# 6. SMC DEG analysis with edgeR -------------------------------------------------------------
## extract SMC cluster
SMC = subset(combined, ident = "SMC")
meta_SMC = SMC@meta.data

## output Ctrl and AngII in SMC
for(i in unique(meta_SMC$orig.ident)){
  each = meta_SMC[meta_SMC$orig.ident==i,]
  UMI = SMC@assays$RNA@counts
  UMI_each = UMI[,colnames(UMI) %in% row.names(each)]
  write.csv(UMI_each, paste0(i,"_SMC_fromIntegrative.csv"), quote = F, row.names = T)
}

## create Seurat objects for Ctrl and AngII in SMC serpartely
Ctrl_SMC  = set_up("Ctrl_SMC_fromIntegrative.csv", "Ctrl_SMC", "Ctrl_SMC")
AngII_SMC  = set_up("AngII_SMC_fromIntegrative.csv", "AngII_SMC", "AngII_SMC")

## Integrate two Seurat objects
anchors = FindIntegrationAnchors(object.list = list(Ctrl_SMC, AngII_SMC), dims = 1:20)
SMC = IntegrateData(anchorset = anchors, dims = 1:20)
SMC = ScaleData(SMC, verbose = FALSE)

## Extract UMI for DEG analysis
counts = as.matrix(SMC@assays$RNA@counts)
meta = SMC@meta.data
sce = SingleCellExperiment(assays = list(counts = counts))

## filtering
filter = rowSums(assay(sce)>0) > dim(counts)[2]/4
sce = sce[filter,]

## fit a ZINB regression model
zinb = zinbFit(sce, K = 2, epsilon = 1000, BPPARAM = BiocParallel::bpparam())
saveRDS(zinb, "zinb_SMC.rds")
sce_zinb = zinbwave(sce, fitted_model = zinb, K = 2, epsilon = 1000, observationalWeights = TRUE)

## DEGs by edgeR
weights = assay(sce_zinb, "weights")
dge = DGEList(assay(sce_zinb))
dge = calcNormFactors(dge)
meta$stim = factor(meta$stim, c("Ctrl_SMC", "AngII_SMC"))
design = model.matrix(~meta$stim)
dge$weights = weights
dge = estimateDisp(dge, design)
fit = glmFit(dge, design)
lrt = glmWeightedF(fit, coef = 2)
lrt_SMC = topTags(lrt, n = nrow(lrt))
write.csv(lrt_SMC, "DEG_SMC.csv", quote = F, row.names = T)

# --- Supplemental Figure 4A --- #
lrt_SMC = data.frame(lrt_SMC)
lrt_SMC$color = with(lrt_SMC, ifelse(FDR > 0.05, 0, ifelse(logFC > 0, 1, 2)))

ggplot(data = lrt_SMC, aes(x = logFC, y = -log10(FDR)))+
  geom_point(size = 1, aes(color=factor(color))) +
  scale_color_manual(values = c("lightgray", "tomato", "dodgerblue", "black"))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  scale_x_continuous(limits = c(-10, 10), breaks = seq(-10, 10, 5))+
  theme_bw()+
  theme(axis.text = element_text(size=16, color = "black"),
        axis.title = element_text(size=16, color = "black"),
        legend.position="none")+
  labs(y = "-Log10(FDR)", x = "Log2FC(AngII/Ctrl)")
ggsave("Sup_Fig_4A.png", width=3.1, height=3)

# --- Supplemental Figure 4B --- #
lrt_SMC_top10 = lrt_SMC
lrt_SMC_top10 = lrt_SMC_top10 %>% filter(FDR < 0.05)
lrt_SMC_top10_up = lrt_SMC_top10 %>% arrange(-logFC) %>% head(10)
lrt_SMC_top10_down = lrt_SMC_top10 %>% arrange(logFC) %>% head(10)
lrt_SMC_top10 = rbind(lrt_SMC_top10_up, lrt_SMC_top10_down) %>% arrange(-logFC)
lrt_SMC_top10$color = with(lrt_SMC_top10, ifelse(logFC > 0, "tomato", "dodgerblue"))
lrt_SMC_top10$X = row.names(lrt_SMC_top10)
lrt_SMC_top10$X = factor(lrt_SMC_top10$X, lrt_SMC_top10$X)

ggplot(lrt_SMC_top10, aes(x = X, y = logFC)) +
  geom_bar(stat = "identity",  width=0.7, fill = lrt_SMC_top10$color, color = lrt_SMC_top10$color) +
  theme_bw()+
  scale_y_continuous(limits = c(-8, 8), breaks = seq(-8, 8, 4)) +
    theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") +
  labs(y = "Log2FC(AngII/Ctrl)", x = "")
ggsave("Sup_Fig_4B.png", width=7.5, height=3.5)

# --- Figure 4E --- #
## Contraction
genes = read.table("input_contraction.txt", header=F, stringsAsFactors = F)
lrt_SMC_Contraction = lrt_SMC
lrt_SMC_Contraction$color = with(lrt_SMC_Contraction, ifelse(logFC > 0, "tomato", "dodgerblue"))
lrt_SMC_Contraction = lrt_SMC_Contraction %>% filter(FDR < 0.05)
lrt_SMC_Contraction$X = row.names(lrt_SMC_Contraction)
lrt_SMC_Contraction = lrt_SMC_Contraction %>% dplyr::filter(X %in% genes$V1) %>% arrange(-logFC)
lrt_SMC_Contraction$X = factor(lrt_SMC_Contraction$X, lrt_SMC_Contraction$X)

ggplot(lrt_SMC_Contraction, aes(x = X, y = logFC)) +
  geom_bar(stat = "identity",  width=0.7, fill = lrt_SMC_Contraction$color, color = lrt_SMC_Contraction$color) +
  scale_y_continuous(limits = c(-1.0, 1.0), breaks = seq(-2, 2, 0.5)) +
  theme_bw()+
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") +
  labs(y = "Log2FC(AngII/Ctrl)", x = "")
ggsave("Fig_4E_Contraction.png", width = 2.3, height = 4.0)

## Proliferation
genes = read.table("input_proliferation.txt", header=F, stringsAsFactors = F)
lrt_SMC_Proliferation = lrt_SMC
lrt_SMC_Proliferation$color = with(lrt_SMC_Proliferation, ifelse(logFC > 0, "tomato", "dodgerblue"))
lrt_SMC_Proliferation = lrt_SMC_Proliferation %>% filter(FDR < 0.05)
lrt_SMC_Proliferation$X = row.names(lrt_SMC_Proliferation)
lrt_SMC_Proliferation = lrt_SMC_Proliferation %>% dplyr::filter(X %in% genes$V1) %>% arrange(-logFC)
lrt_SMC_Proliferation$X = factor(lrt_SMC_Proliferation$X, lrt_SMC_Proliferation$X)

ggplot(lrt_SMC_Proliferation, aes(x = X, y = logFC)) +
  geom_bar(stat = "identity",  width=0.7, fill = lrt_SMC_Proliferation$color, color = lrt_SMC_Proliferation$color) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  theme_bw()+
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") +
  labs(y = "Log2FC(AngII/Ctrl)", x = "")
ggsave("Fig_4E_Proliferation.png", width = 2.5, height = 4.0)

## LRP1 ligands
genes = read.table("input_lrp1_ligands.txt", header=F, stringsAsFactors = F)
lrt_SMC_LRP1 = lrt_SMC
lrt_SMC_LRP1$color = with(lrt_SMC_LRP1, ifelse(logFC > 0, "tomato", "dodgerblue"))
lrt_SMC_LRP1 = lrt_SMC_LRP1 %>% filter(FDR < 0.05)
lrt_SMC_LRP1$X = row.names(lrt_SMC_LRP1)
lrt_SMC_LRP1 = lrt_SMC_LRP1 %>% dplyr::filter(X %in% genes$V1) %>% arrange(-logFC)
lrt_SMC_LRP1$X = factor(lrt_SMC_LRP1$X, lrt_SMC_LRP1$X)

ggplot(lrt_SMC_LRP1, aes(x = X, y = logFC)) +
  geom_bar(stat = "identity",  width=0.7, fill = lrt_SMC_LRP1$color, color = lrt_SMC_LRP1$color) +
  scale_y_continuous(limits = c(-8, 4), breaks = seq(-8, 4, 4)) +
  theme_bw()+
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") +
  labs(y = "Log2FC(AngII/Ctrl)", x = "")
ggsave("Fig_4E_LRP1.png", width = 2.6, height = 4.3)

## ECM components
genes = read.table("input_ECM_comp.txt", header=F, stringsAsFactors = F)
lrt_SMC_ECM = lrt_SMC
lrt_SMC_ECM$color = with(lrt_SMC_ECM, ifelse(logFC > 0, "tomato", "dodgerblue"))
lrt_SMC_ECM = lrt_SMC_ECM %>% filter(FDR < 0.05)
lrt_SMC_ECM$X = row.names(lrt_SMC_ECM)
lrt_SMC_ECM = lrt_SMC_ECM %>% dplyr::filter(X %in% genes$V1) %>% arrange(-logFC)
lrt_SMC_ECM$X = factor(lrt_SMC_ECM$X, lrt_SMC_ECM$X)

ggplot(lrt_SMC_ECM, aes(x = X, y = logFC)) +
  geom_bar(stat = "identity",  width=0.7, fill = lrt_SMC_ECM$color, color = lrt_SMC_ECM$color) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-2, 2, 0.5)) +
  theme_bw()+
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") +
  labs(y = "Log2FC(AngII/Ctrl)", x = "")
ggsave("Fig_4E_ECM.png", width = 3.8, height = 4)

## TGF beta
genes = read.table("input_TGFb.txt", header=F, stringsAsFactors = F)
lrt_SMC_TGFb = lrt_SMC
lrt_SMC_TGFb$color = with(lrt_SMC_TGFb, ifelse(logFC > 0, "tomato", "dodgerblue"))
lrt_SMC_TGFb = lrt_SMC_TGFb %>% filter(FDR < 0.05)
lrt_SMC_TGFb$X = row.names(lrt_SMC_TGFb)
lrt_SMC_TGFb = lrt_SMC_TGFb %>% dplyr::filter(X %in% genes$V1) %>% arrange(-logFC)
lrt_SMC_TGFb$X = factor(lrt_SMC_TGFb$X, lrt_SMC_TGFb$X)

ggplot(lrt_SMC_TGFb, aes(x = X, y = logFC)) +
  geom_bar(stat = "identity",  width=0.7, fill = lrt_SMC_TGFb$color, color = lrt_SMC_TGFb$color) +
  scale_y_continuous(limits = c(-1.0, 1.0), breaks = seq(-2, 2, 0.5)) +
  theme_bw()+
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") +
  labs(y = "Log2FC(AngII/Ctrl)", x = "")
ggsave("Fig_4E_TGFB.png", width = 3.6, height = 4)

# --- Supplemental Figure 4C --- #
lrt_SMC_go = lrt_SMC %>% filter(FDR < 0.05)
lrt_SMC_go = lrt_SMC_go %>% filter(logFC > 0.5 | logFC < -0.5)
lrt_SMC_go$X = row.names(lrt_SMC_go)
geneIDs = ensembldb::select(EnsDb.Mmusculus.v79, keys = as.character(lrt_SMC_go$X), keytype = "SYMBOL", columns = c("ENTREZID","SYMBOL"))
lrt_SMC_go = lrt_SMC_go %>% dplyr::left_join(geneIDs, c("X" = "SYMBOL")) 
lrt_SMC_go = na.omit(lrt_SMC_go)
gene_Enrich = as.character(lrt_SMC_go$ENTREZID)

dat_background = toTable(org.Mm.egSYMBOL)
go_bp = enrichGO(gene = gene_Enrich,
                 universe = dat_background$gene_id,
                 OrgDb = org.Mm.eg.db,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable = TRUE)
write.csv(go_bp, "go_bp_SMC.csv", row.names = F, quote = F)

go_bp_5 = go_bp %>% arrange(p.adjust) %>% head(5)
go_bp_5$Description = factor(go_bp_5$Description, rev(go_bp_5$Description))
ggplot(data = go_bp_5, aes(x = Description, y = -log10(p.adjust))) +
  geom_point(aes(size = Count)) +
  scale_size(range = c(5, 10)) +
  scale_y_continuous(limits = c(10, 18), breaks = seq(10, 18, 4)) +
  theme_bw()+
  theme(axis.text = element_text(size=20, colour = "black"), 
        legend.position="right", text = element_text(size = 20)) + 
  labs(y = "-Log10FDR", x = "") +
  coord_flip()
ggsave("Sup_Fig_4C.png", width = 9.5, height = 3)

# 7. FB DEG analysis with edgeR -------------------------------------------------------------
## extract FB cluster
FB = subset(combined, ident = "FB")
meta_FB = FB@meta.data

## output Ctrl and AngII in FB
for(i in unique(meta_FB$orig.ident)){
  each = meta_FB[meta_FB$orig.ident==i,]
  UMI = FB@assays$RNA@counts
  UMI_each = UMI[,colnames(UMI) %in% row.names(each)]
  write.csv(UMI_each, paste0(i,"_FB_fromIntegrative.csv"), quote = F, row.names = T)
}

## create Seurat objects for Ctrl and AngII in FB serpartely
Ctrl_FB  = set_up("Ctrl_FB_fromIntegrative.csv", "Ctrl_FB", "Ctrl_FB")
AngII_FB  = set_up("AngII_FB_fromIntegrative.csv", "AngII_FB", "AngII_FB")

## Integrate two Seurat objects
anchors = FindIntegrationAnchors(object.list = list(Ctrl_FB, AngII_FB), dims = 1:20)
FB = IntegrateData(anchorset = anchors, dims = 1:20)
FB = ScaleData(FB, verbose = FALSE)
FB@meta.data$stim = gsub("_FB","", FB@meta.data$orig.ident)

## Extract UMI for DEG analysis
counts = as.matrix(FB@assays$RNA@counts)
meta = FB@meta.data
sce = SingleCellExperiment(assays = list(counts = counts))

## filtering
filter = rowSums(assay(sce)>0) > dim(counts)[2]/4
sce = sce[filter,]

## fit a ZINB regression model
zinb = zinbFit(sce, K = 2, epsilon = 1000, BPPARAM = BiocParallel::bpparam())
saveRDS(zinb, "zinb_FB.rds")
sce_zinb = zinbwave(sce, fitted_model = zinb, K = 2, epsilon = 1000, observationalWeights = TRUE)

## DEGs by edgeR
weights = assay(sce_zinb, "weights")
dge = DGEList(assay(sce_zinb))
dge = calcNormFactors(dge)
meta$stim = factor(meta$stim, c("Ctrl", "AngII"))
design = model.matrix(~meta$stim)
dge$weights = weights
dge = estimateDisp(dge, design)
fit = glmFit(dge, design)
lrt = glmWeightedF(fit, coef = 2)
lrt_FB = topTags(lrt, n = nrow(lrt))
write.csv(lrt_FB, "DEG_FB.csv", quote = F, row.names = T)

# --- Supplemental Figure 5A --- #
lrt_FB = data.frame(lrt_FB)
lrt_FB$color = with(lrt_FB, ifelse(FDR > 0.05, 0, ifelse(logFC > 0, 1, 2)))

ggplot(data = lrt_FB, aes(x = logFC, y = -log10(FDR)))+
  geom_point(size = 1, aes(color=factor(color))) +
  scale_color_manual(values = c("lightgray", "tomato", "dodgerblue", "black"))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  scale_x_continuous(limits = c(-10, 10), breaks = seq(-10, 10, 5))+
  theme_bw()+
  theme(axis.text = element_text(size=16, color = "black"),
        axis.title = element_text(size=16, color = "black"),
        legend.position="none")+
  labs(y = "-Log10(FDR)", x = "Log2FC(AngII/Ctrl)")
ggsave("Sup_Fig_5A.png", width=3.1, height=3)

# --- Supplemental Figure 5B --- #
lrt_FB_top10 = lrt_FB
lrt_FB_top10 = lrt_FB_top10 %>% filter(FDR < 0.05)
lrt_FB_top10_up = lrt_FB_top10 %>% arrange(-logFC) %>% head(10)
lrt_FB_top10_down = lrt_FB_top10 %>% arrange(logFC) %>% head(10)
lrt_FB_top10 = rbind(lrt_FB_top10_up, lrt_FB_top10_down) %>% arrange(-logFC)
lrt_FB_top10$color = with(lrt_FB_top10, ifelse(logFC > 0, "tomato", "dodgerblue"))
lrt_FB_top10$X = row.names(lrt_FB_top10)
lrt_FB_top10$X = factor(lrt_FB_top10$X, lrt_FB_top10$X)

ggplot(lrt_FB_top10, aes(x = X, y = logFC)) +
  geom_bar(stat = "identity",  width=0.7, fill = lrt_FB_top10$color, color = lrt_FB_top10$color) +
  theme_bw()+
  scale_y_continuous(limits = c(-8, 8), breaks = seq(-8, 8, 4)) +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") +
  labs(y = "Log2FC(AngII/Ctrl)", x = "")
ggsave("Sup_Fig_5B.png", width=7.5, height=4)

# --- Figure 4F --- #
## Proliferation
genes = read.table("input_proliferation.txt", header=F, stringsAsFactors = F)
lrt_FB_Proliferation = lrt_FB
lrt_FB_Proliferation$color = with(lrt_FB_Proliferation, ifelse(logFC > 0, "tomato", "dodgerblue"))
lrt_FB_Proliferation = lrt_FB_Proliferation %>% filter(FDR < 0.05)
lrt_FB_Proliferation$X = row.names(lrt_FB_Proliferation)
lrt_FB_Proliferation = lrt_FB_Proliferation %>% dplyr::filter(X %in% genes$V1) %>% arrange(-logFC)
lrt_FB_Proliferation$X = factor(lrt_FB_Proliferation$X, lrt_FB_Proliferation$X)

ggplot(lrt_FB_Proliferation, aes(x = X, y = logFC)) +
  geom_bar(stat = "identity",  width=0.7, fill = lrt_FB_Proliferation$color, color = lrt_FB_Proliferation$color) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  theme_bw()+
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") +
  labs(y = "Log2FC(AngII/Ctrl)", x = "")
ggsave("Fig_4F_Proliferation.png", width = 2.6, height = 4.0)

## LRP1 ligands
genes = read.table("input_lrp1_ligands.txt", header=F, stringsAsFactors = F)
lrt_FB_LRP1 = lrt_FB
lrt_FB_LRP1$color = with(lrt_FB_LRP1, ifelse(logFC > 0, "tomato", "dodgerblue"))
lrt_FB_LRP1 = lrt_FB_LRP1 %>% filter(FDR < 0.05)
lrt_FB_LRP1$X = row.names(lrt_FB_LRP1)
lrt_FB_LRP1 = lrt_FB_LRP1 %>% dplyr::filter(X %in% genes$V1) %>% arrange(-logFC)
lrt_FB_LRP1$X = factor(lrt_FB_LRP1$X, lrt_FB_LRP1$X)

ggplot(lrt_FB_LRP1, aes(x = X, y = logFC)) +
  geom_bar(stat = "identity",  width=0.7, fill = lrt_FB_LRP1$color, color = lrt_FB_LRP1$color) +
  scale_y_continuous(limits = c(-1.02, 1.02), breaks = seq(-1, 1, 0.5)) +
  theme_bw()+
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") +
  labs(y = "Log2FC(AngII/Ctrl)", x = "")
ggsave("Fig_4F_LRP1.png", width = 2.2, height = 4.2)

## ECM components
genes = read.table("input_ECM_comp.txt", header=F, stringsAsFactors = F)
lrt_FB_ECM = lrt_FB
lrt_FB_ECM$color = with(lrt_FB_ECM, ifelse(logFC > 0, "tomato", "dodgerblue"))
lrt_FB_ECM = lrt_FB_ECM %>% filter(FDR < 0.05)
lrt_FB_ECM$X = row.names(lrt_FB_ECM)
lrt_FB_ECM = lrt_FB_ECM %>% dplyr::filter(X %in% genes$V1) %>% arrange(-logFC)
lrt_FB_ECM$X = factor(lrt_FB_ECM$X, lrt_FB_ECM$X)

ggplot(lrt_FB_ECM, aes(x = X, y = logFC)) +
  geom_bar(stat = "identity",  width=0.7, fill = lrt_FB_ECM$color, color = lrt_FB_ECM$color) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  theme_bw()+
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") +
  labs(y = "Log2FC(AngII/Ctrl)", x = "")
ggsave("Fig_4F_ECM.png", width = 4.2, height = 4)

## TGF beta
genes = read.table("input_TGFb.txt", header=F, stringsAsFactors = F)
lrt_FB_TGFb = lrt_FB
lrt_FB_TGFb$color = with(lrt_FB_TGFb, ifelse(logFC > 0, "tomato", "dodgerblue"))
lrt_FB_TGFb = lrt_FB_TGFb %>% filter(FDR < 0.05)
lrt_FB_TGFb$X = row.names(lrt_FB_TGFb)
lrt_FB_TGFb = lrt_FB_TGFb %>% dplyr::filter(X %in% genes$V1) %>% arrange(-logFC)
lrt_FB_TGFb$X = factor(lrt_FB_TGFb$X, lrt_FB_TGFb$X)

ggplot(lrt_FB_TGFb, aes(x = X, y = logFC)) +
  geom_bar(stat = "identity",  width=0.7, fill = lrt_FB_TGFb$color, color = lrt_FB_TGFb$color) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  theme_bw()+
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") +
  labs(y = "Log2FC(AngII/Ctrl)", x = "")
ggsave("Fig_4F_TGFB.png", width = 3, height = 4)

# --- Supplemental Figure 5C --- #
lrt_FB_go = lrt_FB %>% filter(FDR < 0.05)
lrt_FB_go = lrt_FB_go %>% filter(logFC > 0.5 | logFC < -0.5)
lrt_FB_go$X = row.names(lrt_FB_go)
geneIDs = ensembldb::select(EnsDb.Mmusculus.v79, keys= as.character(lrt_FB_go$X), keytype = "SYMBOL", columns = c("ENTREZID","SYMBOL"))
lrt_FB_go = lrt_FB_go %>% dplyr::left_join(geneIDs, c("X" = "SYMBOL")) 
lrt_FB_go = na.omit(lrt_FB_go)
gene_Enrich = as.character(lrt_FB_go$ENTREZID)

dat_background = toTable(org.Mm.egSYMBOL)
go_bp = enrichGO(gene = gene_Enrich,
                 universe = dat_background$gene_id,
                 OrgDb = org.Mm.eg.db,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable = TRUE)

go_bp_5 = go_bp %>% arrange(p.adjust) %>% head(5)
go_bp_5$Description = factor(go_bp_5$Description, rev(go_bp_5$Description))
ggplot(data = go_bp_5, aes(x = Description, y = -log10(p.adjust))) +
  geom_point(aes(size = Count)) +
  scale_size(range = c(5, 10)) +
  scale_y_continuous(limits = c(15, 18), breaks = seq(15, 18, 1)) +
  theme_bw()+
  theme(axis.text = element_text(size=20, colour = "black"), 
        legend.position="right", text = element_text(size = 20)) + 
  labs(y = "-Log10FDR", x = "") +
  coord_flip()
ggsave("Sup_Fig_5C.png", width = 9.5, height = 3)

# 8. sub-clustering in SMC ----------------------------------------------------
Ctrl_SMC  = set_up("Ctrl_SMC_fromIntegrative.csv", "Ctrl", "Ctrl")
AngII_SMC  = set_up("AngII_SMC_fromIntegrative.csv", "AngII", "AngII")

## Change the order Ctrl -> AngII
SMC@meta.data$stim = factor(SMC@meta.data$stim, c("Ctrl", "AngII")) 

## Combined Ctrl and AngII SMC data sets
anchors = FindIntegrationAnchors(object.list = list(Ctrl_SMC, AngII_SMC), dims = 1:20)
SMC = IntegrateData(anchorset = anchors, dims = 1:20)
SMC = ScaleData(SMC, verbose = FALSE)

## Clustering and UMAP
SMC = RunPCA(SMC, npcs = 30, verbose = FALSE)
SMC = RunUMAP(SMC, reduction = "pca", dims = 1:20, seed.use = 100)
SMC = FindNeighbors(SMC, reduction = "pca", dims = 1:20)
SMC = FindClusters(SMC, resolution = 0.1)
SMC@meta.data$stim = factor(SMC@meta.data$stim, c("Ctrl", "AngII"))

p1 = DimPlot(SMC, reduction = "umap", group.by = "stim", pt.size = 0.8)
p2 = DimPlot(SMC, reduction = "umap", label = TRUE, repel = TRUE)
ggsave("umap_SMC.png", p1 + p2, height=4, width=8)

DimPlot(SMC, reduction = "umap", split.by = "stim", pt.size = 0.8)
ggsave("Sup_Fig_7.png", height=4, width=8)

saveRDS(SMC, "SMC.rds")

# 9. sub-clustering in FB ----------------------------------------------------
FB = RunPCA(FB, npcs = 30, verbose = FALSE)
FB = RunUMAP(FB, reduction = "pca", dims = 1:25, seed.use = 10000)
FB = FindNeighbors(FB, reduction = "pca", dims = 1:25)
FB = FindClusters(FB, resolution = 0.1)

# --- Figure 5A --- #
FB@meta.data$stim = factor(FB@meta.data$stim, c("Ctrl", "AngII"))
DimPlot(FB, reduction = "umap", split.by = "stim")
Idents(FB) = FB@meta.data$seurat_clusters
FB = RenameIdents(FB, `0` = "FB1", `1` = "FB2", `2` = "FB3", `3` = "FB4")
FB@meta.data$celltype = FB@active.ident
Idents(FB) = factor(Idents(FB), paste0("FB", 1:4))
FB@meta.data$stim = factor(FB@meta.data$stim, c("Ctrl", "AngII"))

DimPlot(FB, reduction = "umap", split.by = "stim", pt.size = 0.8, dims = c(2,1))
ggsave("Fig_5A.png", height=4, width=7)

# --- Figure 5B --- #
markers_FB = FindAllMarkers(FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers_FB, "marker genes - FB sub clusters.csv", quote = F)

top10_markers = markers_FB %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_markers_FB4 = top10_markers %>% filter(cluster == "FB4") %>% arrange(-avg_log2FC)
top10_markers_FB4$gene = factor(top10_markers_FB4$gene, top10_markers_FB4$gene)
ggplot(top10_markers_FB4, aes(x = gene, y = avg_log2FC)) +
  geom_bar(stat = "identity", width=0.7) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 4, 1)) +
  theme_bw()+
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.ticks.x = element_blank(),
        #axis.line.x = element_blank(),
        legend.position="none") +
  labs(y = "Log2FC", x = "")
ggsave("Fig_5B.png", width = 4.2, height = 3)

# --- Figure 5C --- #
FeaturePlot(FB, features = c("Col1a1", "Ly6a"), min.cutoff = "q9", cols = c("lightgrey", "goldenrod","firebrick"), pt.size = 0.2, ncol = 2, dims = c(2,1))
ggsave("Fig_5C.png", width = 6.5, height = 3)

# --- Figure 5D --- #
FB = FindVariableFeatures(FB, selection.method = "vst", nfeatures = 2000)
varGene = VariableFeatures(FB)
avg = AverageExpression(FB, group.by = "celltype")
avg = avg$RNA
data = avg[row.names(avg) %in% varGene,]
res = cor(data, method = "spearman")
png("Fig_5D.png")
pheatmap(res, main = "", cellwidth = 30, cellheight = 30, fontsize = 20, angle_col = 90)
dev.off()

# --- Figure 5E --- #
## ligand and receptor
genes = read.table("input_ligand_receptor.txt", sep = "\t", h = T)
ligand = genes$ligand
receptor = genes$receptor
combined@meta.data$celltype = combined@active.ident
celltype = data.frame(name = names(combined@active.ident), celltype = combined@active.ident)
FB = SetIdent(FB, value = "celltype")
celltype_FB = data.frame(name = names(FB@active.ident), celltype = FB@active.ident)
celltype_FB$name = gsub(".1_1", "-1", celltype_FB$name)
celltype_FB$name = gsub(".1_2", ".1", celltype_FB$name)

celltype = celltype %>% left_join(celltype_FB, c("name" = "name"))
celltype$celltype.x = as.character(celltype$celltype.x)
celltype$celltype.y = as.character(celltype$celltype.y)
celltype$celltype = with(celltype, ifelse(is.na(celltype.y), celltype.x, celltype.y))

combined@meta.data$celltype = (celltype$celltype)
combined = SetIdent(combined, value = "celltype")
celltype = c("SMC", "FB1", "FB2", "FB3", "FB4", "EC", "UI", "Mac")

for (g in c("AngII", "Ctrl")){
  eval(call("<-", as.name(g), SetIdent(combined, value = "stim")))
  eval(call("<-", as.name(g), subset(eval(parse(text=g)), ident = g)))
  eval(call("<-", as.name(g), SetIdent(eval(parse(text=g)), value = "celltype")))
  for (cell in celltype){
    objname = paste0(g, "_", cell)
    #eval(call("<-", as.name(objname), SetIdent(eval(parse(text=g)), value = cell)))
    eval(call("<-", as.name(objname), subset(eval(parse(text=g)), ident = cell)))
    eval(call("<-", as.name(objname), eval(parse(text=objname))@assays$RNA@counts))
  }
}

re_AngII = matrix(nrow = length(celltype), nco=length(celltype))
re_Ctrl = matrix(nrow = length(celltype), nco=length(celltype))

g="Ctrl"
for (i in 1:length(celltype)) {
  for (j in 1:length(celltype)) {
    objname_ligand = paste0(g, "_", celltype[i])
    counts_ligand = eval(parse(text=objname_ligand))
    genename_ligand = rownames(counts_ligand)
    counts_ligand = as.matrix(counts_ligand)
    m_ligand = ncol(counts_ligand)
    
    objname_receptor = paste0(g, "_", celltype[j])
    counts_receptor = eval(parse(text=objname_receptor))
    genename_receptor = rownames(counts_receptor)
    counts_receptor = as.matrix(counts_receptor)
    m_receptor = ncol(counts_receptor)
    
    re = c()
    for (k in 1:nrow(genes)) {
      flag_ligand = (genename_ligand %in% genes$ligand[k])
      flag_receptor = (genename_receptor %in%  genes$receptor[k])
      
      if (sum(flag_ligand) == 1) {
        tmp_ligand =  counts_ligand[flag_ligand,]
        re_ligand = as.numeric(sum(tmp_ligand == 0) < m_ligand*0.2)
      } else {
        re_ligand = 0
      }
      if (sum(flag_receptor) == 1) {
        tmp_receptor =  counts_receptor[flag_receptor,]
        re_receptor = as.numeric(sum(tmp_receptor == 0) < m_receptor*0.2)
      } else {
        re_receptor = 0
      }
      
      re = c(re, re_ligand*re_receptor)
    }
    re_sum = sum(re)
    re_Ctrl[i,j] = re_sum
  }
}

g="AngII"
for (i in 1:length(celltype)) {
  for (j in 1:length(celltype)) {
    objname_ligand = paste0(g, "_", celltype[i])
    counts_ligand = eval(parse(text=objname_ligand))
    genename_ligand = rownames(counts_ligand)
    counts_ligand = as.matrix(counts_ligand)
    m_ligand = ncol(counts_ligand)
    
    objname_receptor = paste0(g, "_", celltype[j])
    counts_receptor = eval(parse(text=objname_receptor))
    genename_receptor = rownames(counts_receptor)
    counts_receptor = as.matrix(counts_receptor)
    m_receptor = ncol(counts_receptor)
    
    re = c()
    for (k in 1:nrow(genes)) {
      flag_ligand = (genename_ligand %in% genes$ligand[k])
      flag_receptor = (genename_receptor %in%  genes$receptor[k])
      
      if (sum(flag_ligand) == 1) {
        tmp_ligand =  counts_ligand[flag_ligand,]
        re_ligand = as.numeric(sum(tmp_ligand == 0) < m_ligand*0.2)
      } else {
        re_ligand = 0
      }
      if (sum(flag_receptor) == 1) {
        tmp_receptor =  counts_receptor[flag_receptor,]
        re_receptor = as.numeric(sum(tmp_receptor == 0) < m_receptor*0.2)
      } else {
        re_receptor = 0
      }
      
      re = c(re, re_ligand*re_receptor)
    }
    re_sum = sum(re)
    re_AngII[i,j] = re_sum
  }
}

rownames(re_Ctrl) = celltype
colnames(re_Ctrl) = celltype

rownames(re_AngII) = celltype
colnames(re_AngII) = celltype

d = re_Ctrl
rownames(d) = colnames(d) = c("SMC", "FB2", "FB1", "FB3", "FB4", "EC", "UI", "Mac")
cell_order = c("FB4", "FB3", "FB2", "FB1", "SMC", "EC", "Mac", "UI")
d = d[cell_order, cell_order]
d = d[1:5, 1:5]
d = d*0.4
color = c("#C77CFF", "#00BFC4", "#7CAE00", "#F8766D", "#CD9600")

g = graph.adjacency(d, weighted = T, mode = "directed")
E(g)$weight = E(g)$weight*5

l = t(matrix(c(0,1,-cos(pi/10),sin(pi/10),-cos(3*pi/10),-sin(3*pi/10),cos(3*pi/10),-sin(3*pi/10),cos(pi/10),sin(pi/10)), ncol = 5, nrow = 2))

png("Fig_5E_Ctrl.png", width = 3500, height = 3000)
plot(g, edge.width = E(g)$weight,
     ylim = c(-1.4,1.4),
     
     vertex.size = c(5, 15, 40, 21, 33),
     vertex.shape = "circle",
     vertex.label = V(g)$name,
     vertex.color = color,
     vertex.label.color = "black",
     vertex.label.family = "sans",
     vertex.label.font = 1,
     vertex.frame.color = "darkgray",
     vertex.label.cex = 20,
     edge.loop.angle = c(c(rep(0, 4)),
                         c(rep(0, 1), -pi*3/4, rep(0, 3)),
                         c(rep(0, 1), -pi*5/4, rep(0, 2)),
                         c(rep(0, 3), -pi*7/4, rep(0, 1)),
                         c(rep(0, 4), -pi*9/4)), 
     edge.curved = 0.2,
     edge.arrow.size = 0,
     edge.width = E(g)$weight,
     edge.color = c(rep("#C77CFF", 4),
                    rep("#00BFC4", 5),
                    rep("#7CAE00", 4),
                    rep("#F8766D", 5),
                    rep("#CD9600", 5)),
     layout = l)
dev.off()

d = re_AngII
rownames(d) = colnames(d) = c("SMC", "FB2", "FB1", "FB3", "FB4", "EC", "UI", "Mac")
cell_order = c("FB4", "FB3", "FB2", "FB1", "SMC", "EC", "Mac", "UI")
d = d[cell_order, cell_order]
d = d[1:5, 1:5]
d = d*0.4
color = c("#C77CFF", "#00BFC4", "#7CAE00", "#F8766D", "#CD9600")

g = graph.adjacency(d, weighted = T, mode = "directed")
E(g)$weight = E(g)$weight*5

l = t(matrix(c(0,1,-cos(pi/10),sin(pi/10),-cos(3*pi/10),-sin(3*pi/10),cos(3*pi/10),-sin(3*pi/10),cos(pi/10),sin(pi/10)), ncol = 5, nrow = 2))

png("Fig_5E_AngII.png", width = 3500, height = 3000)
plot(g, edge.width = E(g)$weight,
     ylim = c(-1.4,1.4),
     
     vertex.size = c(5, 15, 40, 21, 33),
     vertex.shape = "circle",
     vertex.label = V(g)$name,
     vertex.color = color,
     vertex.label.color = "black",
     vertex.label.family = "sans",
     vertex.label.font = 1,
     vertex.frame.color = "darkgray",
     vertex.label.cex = 20,
     edge.loop.angle = c(-pi*2/4, c(rep(0, 4)),
                       c(rep(0, 1), -pi*3/4, rep(0, 3)),
                       c(rep(0, 2), -pi*5/4, rep(0, 2)),
                       c(rep(0, 3), -pi*7/4, rep(0, 1)),
                       c(rep(0, 4), -pi*9/4)), 
     edge.curved = 0.2,
     edge.arrow.size = 0, 
     edge.width = E(g)$weight,
     edge.color = c(rep("#C77CFF", 5),
                    rep("#00BFC4", 5),
                    rep("#7CAE00", 5),
                    rep("#F8766D", 5),
                    rep("#CD9600", 5)),
     layout = l)
dev.off()

# --- Figure 5F --- #
combined@meta.data$celltype = Idents(combined)
SMC_FB1_4 = subset(combined, ident = c("SMC", "FB1", "FB2", "FB3", "FB4"))
SMC_FB1_4@meta.data$celltype = as.factor(as.character(SMC_FB1_4@meta.data$celltype))
SMC_FB1_4 = SetIdent(SMC_FB1_4, value = "orig.ident")
SMC_FB1_4_AngII = subset(SMC_FB1_4, ident = "AngII")
UMI = SMC_FB1_4_AngII@assays$RNA@counts
meta = SMC_FB1_4_AngII@meta.data
sample_sheet = data.frame(meta)
gene_anno = data.frame(row.names(UMI))
row.names(gene_anno) = gene_anno[,1]
colnames(gene_anno) = "gene_short_name"
pd = new("AnnotatedDataFrame", data = sample_sheet)
fd = new("AnnotatedDataFrame", data = gene_anno)
Immu = newCellDataSet(as(as.matrix(UMI), "sparseMatrix"), phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())
Immu = estimateSizeFactors(Immu)
Immu = estimateDispersions(Immu)
Immu = detectGenes(Immu, min_expr = 0.1)
fData(Immu)$use_for_ordering = fData(Immu)$num_cells_expressed > 0.05 * ncol(Immu)
Immu = reduceDimension(Immu, max_components = 2, norm_method = "log", num_dim = 3, reduction_method = "tSNE", verbose = T)
Immu = clusterCells(Immu, verbose = F)
expressed_genes = row.names(subset(fData(Immu), num_cells_expressed >= 10))
diff_test_res = differentialGeneTest(Immu[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
ordering_genes = row.names (subset(diff_test_res, qval < 0.01))
Immu = setOrderingFilter(Immu, ordering_genes)
Immu = reduceDimension(Immu, max_components = 2,   method = 'DDRTree')
Immu = orderCells(Immu)  ###, reverse = T
plot_cell_trajectory(Immu, color_by = "celltype", show_branch_points = FALSE, cell_size = 0.8)
ggsave("Fig_5F_cells.png", width=3, height=3.5)
plot_cell_trajectory(Immu, color_by = "Pseudotime", show_branch_points = FALSE, cell_size = 0.8)
ggsave("Fig_5F_time.png", width=3, height=3.7)

# --- Figure 5G --- #
FB1_FB4 = subset(combined_sub, ident = c("FB1", "FB4"))
FB1_FB4@meta.data$cell_cluster = FB1_FB4@active.ident
FB1_FB4 = SetIdent(FB1_FB4, value = "stim")
FB1_FB4_AngII = subset(FB1_FB4, ident = c("AngII"))

# Extract UMI for DEG analysis
UMI = FB1_FB4_AngII@assays$RNA@counts
meta = FB1_FB4_AngII@meta.data
counts = as.matrix(UMI)
sce = SingleCellExperiment(assays = list(counts = counts))

## filtering
filter = rowSums(assay(sce)>0) > dim(counts)[2]*0.25
sce_tmp = sce[filter,] 
zinb = zinbFit(sce_tmp, K = 2, epsilon = 1000, BPPARAM = BiocParallel::bpparam())
saveRDS(zinb, "zinb_FB4FB1.rds")
sce_zinb = zinbwave(sce_tmp, fitted_model = zinb, K = 2, epsilon = 1000, observationalWeights = TRUE)

## DEGs by edgeR
weights = assay(sce_zinb, "weights")
dge = DGEList(assay(sce_zinb))
dge = calcNormFactors(dge)

design = model.matrix(~meta$cell_cluster)
dge$weights = weights
dge = estimateDisp(dge, design)
fit = glmFit(dge, design)
lrt = glmWeightedF(fit, coef = 2)
lrt_FB1_FB4_AngII = topTags(lrt, n = nrow(lrt))
write.csv(lrt_FB1_FB4_AngII, "DEG_FB4_FB1.csv", quote = F, row.names = T)

## Set DEGs for FB4 DEG analysis
lrt_FB1_FB4_AngII = read.csv("DEG_FB4_FB1.csv", header = T, stringsAsFactors = F)
lrt_FB1_FB4_AngII$color = with(lrt_FB1_FB4_AngII, ifelse(logFC > 0, 1, 2))

## LRP1 ligands
input_lrp1 = read.table("input_lrp1_ligands.txt", header=F, stringsAsFactors = F)
lrt_FB1_FB4_AngII_lrp1 = lrt_FB1_FB4_AngII %>% filter(FDR < 0.05)
lrt_FB1_FB4_AngII_lrp1 = lrt_FB1_FB4_AngII_lrp1 %>% dplyr::filter(X %in% input_lrp1$V1)
lrt_FB1_FB4_AngII_lrp1 = lrt_FB1_FB4_AngII_lrp1[order(lrt_FB1_FB4_AngII_lrp1$logFC, decreasing = T), ]
lrt_FB1_FB4_AngII_lrp1$X = factor(lrt_FB1_FB4_AngII_lrp1$X, lrt_FB1_FB4_AngII_lrp1$X)

ggplot(lrt_FB1_FB4_AngII_lrp1, aes(x = X, y = logFC)) +
  geom_bar(stat = "identity", width=0.7, aes(fill = factor(color), color = factor(color))) +
  scale_fill_manual(values = "dodgerblue") +
  scale_color_manual(values = "dodgerblue") +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  theme_bw()+
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") +
  labs(y = "Log2FC(FB4/FB1)", x = "")
ggsave("Fig_5G_Lrp1.png", width = 2.1, height = 4.2)

## ECM components
input_ecm = read.table("input_ECM_comp.txt", header=F, stringsAsFactors = F)
lrt_FB1_FB4_AngII_ecm = lrt_FB1_FB4_AngII %>% filter(FDR < 0.05)
lrt_FB1_FB4_AngII_ecm = lrt_FB1_FB4_AngII_ecm %>% filter(X %in% input_ecm$V1)
lrt_FB1_FB4_AngII_ecm = lrt_FB1_FB4_AngII_ecm[order(lrt_FB1_FB4_AngII_ecm$logFC, decreasing = T), ]
lrt_FB1_FB4_AngII_ecm$X = factor(lrt_FB1_FB4_AngII_ecm$X, lrt_FB1_FB4_AngII_ecm$X)

ggplot(lrt_FB1_FB4_AngII_ecm, aes(x = X, y = logFC)) +
  geom_bar(stat = "identity",  width=0.7, aes(fill = factor(color), color = factor(color))) +
  scale_fill_manual(values = "dodgerblue") +
  scale_color_manual(values = "dodgerblue") +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  theme_bw()+
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") +
  labs(y = "Log2FC(FB4/FB1)", x = "")
ggsave("Fig_5G_ECM_hs.png", width = 4.4, height = 4)

## TGF
input_tgf = read.table("input_TGFb.txt", header=F, stringsAsFactors = F)
lrt_FB1_FB4_AngII_tgf = lrt_FB1_FB4_AngII %>% filter(FDR < 0.05)
lrt_FB1_FB4_AngII_tgf = lrt_FB1_FB4_AngII_tgf %>% filter(X %in% input_tgf$V1)
lrt_FB1_FB4_AngII_tgf = lrt_FB1_FB4_AngII_tgf[order(lrt_FB1_FB4_AngII_tgf$logFC, decreasing = T), ]
lrt_FB1_FB4_AngII_tgf$X = factor(lrt_FB1_FB4_AngII_tgf$X, lrt_FB1_FB4_AngII_tgf$X)

ggplot(lrt_FB1_FB4_AngII_tgf, aes(x = X, y = logFC)) +
  geom_bar(stat = "identity",  width=0.7, aes(fill = factor(color), color = factor(color))) +
  scale_fill_manual(values = "dodgerblue") +
  scale_color_manual(values = "dodgerblue") +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
  theme_bw()+
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") +
  labs(y = "Log2FC(FB4/FB1)", x = "")
ggsave("Fig_5G_TGFB_hs.png", width = 2.2, height = 4)

# 10. FB4 DEG analysis with edgeR -------------------------------------------------------------
# --- Supplemental Figure 8A --- #
tbl = table(combined_sub@active.ident, combined_sub@meta.data$stim)
write.csv(tbl, "Cell# of FB subclusters.csv", quote = F, row.names = T)

## extract SMC cluster
FB4 = subset(combined_sub, ident = "FB4")

# Extract UMI for DEG analysis
UMI = FB4@assays$RNA@counts
meta = FB4@meta.data
counts = as.matrix(UMI)
sce = SingleCellExperiment(assays = list(counts = counts))

## filtering
filter = rowSums(assay(sce)>0) > dim(counts)[2]*0.25
sce_tmp = sce[filter,] 
zinb = zinbFit(sce_tmp, K = 2, epsilon = 1000, BPPARAM = BiocParallel::bpparam())
saveRDS(zinb, "zinb_FB4.rds")
sce_zinb = zinbwave(sce_tmp, fitted_model = zinb, K = 2, epsilon = 1000, observationalWeights = TRUE)

## DEGs by edgeR
weights = assay(sce_zinb, "weights")
dge = DGEList(assay(sce_zinb))
dge = calcNormFactors(dge)

design = model.matrix(~meta$stim)
dge$weights = weights
dge = estimateDisp(dge, design)
fit = glmFit(dge, design)
lrt = glmWeightedF(fit, coef = 2)
lrt_FB1_FB4_AngII = topTags(lrt, n = nrow(lrt))
write.csv(lrt_FB1_FB4_AngII, "DEG_FB4.csv", quote = F, row.names = T)

# --- Supplemental Figure 8B --- #
deg = read.csv("DEG_FB4.csv", header = T, stringsAsFactors = F)
deg$color = with(deg, ifelse(FDR>0.05, 0, ifelse(logFC>0, 1, 2)))

volplot = ggplot(data = deg, aes(x = logFC, y = -log10(FDR)))+
  geom_point(size = 1, aes(color=factor(color))) +
  scale_color_manual(values = c("lightgray", "tomato", "dodgerblue", "black"))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  scale_x_continuous(limits = c(-4, 8), breaks = seq(-4, 8, 4)) +
  scale_y_continuous(limits = c(0, 30))+
  theme_bw()+
  theme(axis.text = element_text(size=16, color = "black"),
        axis.title = element_text(size=16, color = "black"),
        legend.position="none")+
  labs(y = "-Log10(FDR)", x = "Log2FC(AngII/Ctrl)")
ggsave("Sup_Fig_8B.png", plot = volplot, width=3, height=3)

# --- Supplemental Figure 5C --- #
lrt_FB4_go = read.csv("DEG_FB4.csv", header = T, stringsAsFactors = F) 
lrt_FB4_go = lrt_FB4_go %>% filter(FDR < 0.05)
lrt_FB4_go = lrt_FB4_go %>% filter(logFC > 0.5 | logFC < -0.5)
geneIDs = ensembldb::select(EnsDb.Mmusculus.v79, keys= as.character(lrt_FB4_go$X), keytype = "SYMBOL", columns = c("ENTREZID","SYMBOL"))
lrt_FB4_go = lrt_FB4_go %>% dplyr::left_join(geneIDs, c("X" = "SYMBOL")) 
lrt_FB4_go = na.omit(lrt_FB4_go)
gene_Enrich = as.character(lrt_FB4_go$ENTREZID)

dat_background = toTable(org.Mm.egSYMBOL)
go_bp = enrichGO(gene = gene_Enrich,
                 universe = dat_background$gene_id,
                 OrgDb = org.Mm.eg.db,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable = TRUE)

go_bp_5 = go_bp %>% arrange(p.adjust) %>% head(5)
go_bp_5$Description = factor(go_bp_5$Description, rev(go_bp_5$Description))
ggplot(data = go_bp_5, aes(x = Description, y = -log10(p.adjust))) +
  geom_point(aes(size = Count)) +
  scale_size(range = c(5, 10)) +
  scale_y_continuous(limits = c(8, 9), breaks = seq(8, 9, 0.5)) +
  theme_bw()+
  theme(axis.text = element_text(size=20, colour = "black"), 
        legend.position="right", text = element_text(size = 20)) + 
  labs(y = "-Log10FDR", x = "") +
  coord_flip()
ggsave("Sup_Fig_8B.png", width = 9.5, height = 3)

##### Protemics vs SMC/FB
proteomics_t_test = read.csv("Proteomics_t_test.csv", header = T, stringsAsFactors = F)
proteomics_t_test = proteomics_t_test %>% filter(FDR < 0.05)

scRNAseq_smc = read.csv("DEG_SMC.csv", header = T, stringsAsFactors = F)
scRNAseq_smc = scRNAseq_smc %>% filter(scRNAseq_smc$FDR < 0.05)

scRNAseq_fb = read.csv("DEG_FB.csv", header = T, stringsAsFactors = F)
scRNAseq_fb = scRNAseq_fb %>% filter(scRNAseq_fb$FDR < 0.05)

# --- Figure 6A --- #
##### Proteomics vs SMC
deg = data.frame(gene = proteomics_t_test$X, proteomics_t_test) %>% left_join (scRNAseq_smc, c("gene" = "X"))
deg = subset(deg, !(is.na(deg$logFC)))

input_lrp1 = c("Serpine1", "Thbs1", "Cd44", "Ctgf", "Ctsb", "Tgfb2", "Fn1")
deg$ligand = ifelse(deg$gene %in% input_lrp1, deg$gene, NA) 

ggplot(data = deg, aes(x = LogeFC, y = logFC))+
  geom_point(size = 2.5, shape = 21, color = "gray60", fill = "lightgray") +
  geom_point(data = deg %>% filter(!is.na(ligand)), size = 2.5, shape = 21, color = "gray60", fill = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1))+
  scale_y_continuous(limits = c(-8, 4), breaks = seq(-8, 4, 4))+
  theme_bw()+
  theme(axis.text = element_text(size=16, color = "black"),
        axis.title = element_text(size=16, color = "black"),
        legend.position="none")+
  labs(y = "scRNAseq \n Log2FC(AngII/Ctrl)", x = "Proteomics \n LogFC(AngII/Saline)")
ggsave("Fig_6A.png", width=3.1, height=3)

# --- Figure 6B --- #
##### Proteomics vs FB
deg = data.frame(gene = proteomics_t_test$X, proteomics_t_test) %>% left_join (scRNAseq_fb, c("gene" = "X"))
deg = subset(deg, !(is.na(deg$logFC)))
deg$ligand = ifelse(deg$gene %in% input_lrp1, deg$gene, NA) 

ggplot(data = deg, aes(x = LogeFC, y = logFC))+
  geom_point(size = 2.5, shape = 21, color = "gray60", fill = "lightgray") +
  geom_point(data = deg %>% filter(!is.na(ligand)), size = 2.5, shape = 21, color = "gray60", fill = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1))+
  scale_y_continuous(limits = c(-8, 4), breaks = seq(-8, 4, 4))+
  theme_bw()+
  theme(axis.text = element_text(size=16, color = "black"),
        axis.title = element_text(size=16, color = "black"),
        legend.position="none")+
  labs(y = "scRNAseq \n Log2FC(AngII/Ctrl)", x = "Proteomics \n LogFC(AngII/Saline)")
ggsave("Fig_6B.png", width=3.1, height=3)

# Human scRNA seq analysis ------------------------------------------------
## Read a CSV file for read counts
filename = c("GSM4704931_Con4", "GSM4704932_Con6", "GSM4704933_Con9", "GSM4704934_TAA1", "GSM4704935_TAA2", "GSM4704936_TAA3", "GSM4704937_TAA4", "GSM4704938_TAA5", "GSM4704939_TAA6", "GSM4704940_TAA7", "GSM4704941_TAA8")
for (f in filename){
  dat = read.table(paste0(f,".txt"), header = T, stringsAsFactors = F)
  dat = dat[!duplicated(rownames(dat)),]
  dat = Seurat::CreateSeuratObject(counts = dat, project = f, min.cells = 3, min.features = 200)
  dat$stim = substr(f, 12,14)
  dat[["percent.mt"]] = PercentageFeatureSet(dat, pattern = "^mt-")
  dat = subset(dat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
  dat = NormalizeData(dat, verbose = FALSE)
  dat = FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
  eval(call("<-", as.name(f), dat))
}

## Integrate
objectlist = list(GSM4704931_Con4, GSM4704932_Con6, GSM4704933_Con9, GSM4704934_TAA1, GSM4704935_TAA2, GSM4704936_TAA3, GSM4704937_TAA4, GSM4704938_TAA5, GSM4704939_TAA6, GSM4704940_TAA7, GSM4704941_TAA8)
features = SelectIntegrationFeatures(object.list = objectlist)
anchors = FindIntegrationAnchors(object.list = objectlist, dims = 1:20, anchor.features = features)
combined_hm = IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(combined_hm) = "integrated"

## Clustering and UMAP
combined_hm = ScaleData(combined_hm, verbose = FALSE)
combined_hm = RunPCA(combined_hm, npcs = 30, verbose = FALSE)
combined_hm = RunUMAP(combined_hm, reduction = "pca", dims = 1:20, seed.use = 1)
combined_hm = FindNeighbors(combined_hm, reduction = "pca", dims = 1:20)
combined_hm = FindClusters(combined_hm, resolution = 0.5)

p1 = DimPlot(combined_hm, reduction = "umap", group.by = "stim", pt.size = 0.4)
p2 = DimPlot(combined_hm, reduction = "umap", label = TRUE, repel = TRUE)
ggsave("umap_human.png", p1 + p2, height=5, width=10)

DimPlot(combined_hm, reduction = "umap", split.by = "stim", pt.size = 0.4)
ggsave("umap__human_stim.png", height=5, width=10)

saveRDS(combined_hm, "combined_hm.rds")

# --- Figure 6E --- #
SMC_FB = subset(combined_hm, ident = c("5", "6", "8", "12", "15"))
SMC_FB@meta.data$celltype = ifelse(SMC_FB@active.ident == "6" |  SMC_FB@active.ident == "15", "FB", "SMC")
SMC_FB = ScaleData(SMC_FB, verbose = FALSE)
SMC_FB = RunPCA(SMC_FB, npcs = 30, verbose = FALSE)
SMC_FB = RunUMAP(SMC_FB, reduction = "pca", dims = 1:20, seed.use = 2)
SMC_FB = FindNeighbors(SMC_FB, reduction = "pca", dims = 1:20)
SMC_FB = FindClusters(SMC_FB, resolution = 0.08)
SMC_FB@meta.data$celltype = factor(SMC_FB@meta.data$celltype, c("SMC", "FB"))
DimPlot(SMC_FB, reduction = "umap", split.by = "stim", label = T, group.by = "celltype")
FeaturePlot(SMC_FB, features = c("SERPINE1"), min.cutoff = "q9",cols=c("lightgrey", "goldenrod","firebrick"),pt.size=1, ncol = 2, split.by = "stim")

# --- Figure 6F --- #
counts = as.matrix(SMC_FB@assays$RNA@counts)
meta = SMC_FB@meta.data
sce = SingleCellExperiment(assays = list(counts = counts))
filter = rowSums(assay(sce)>0) > dim(counts)[2]*0.25
sce_tmp = sce[filter,] 
zinb = zinbFit(sce_tmp, K = 2, epsilon = 1000, BPPARAM = BiocParallel::bpparam())
saveRDS(zinb, "zinb_hm.rds")
sce_zinb = zinbwave(sce_tmp, fitted_model = zinb, K = 2, epsilon = 1000, observationalWeights = TRUE)
weights = assay(sce_zinb, "weights")
dge = DGEList(assay(sce_zinb))
dge = calcNormFactors(dge)
tmm = edgeR::cpm(dge, log = T)
tmm_gene = tmm[rownames(tmm) == "SERPINE1",]
SMC_FB = SetIdent(SMC_FB, value = "celltype")
SMC_FB@meta.data$tmm_gene = tmm_gene
VlnPlot(SMC_FB, features = "tmm_gene", split.by = "stim", split.plot = F, pt = 0.1)+
  theme(axis.text.x = element_text(size = 24, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22))

# --- Supplemental Figure 6A --- #
SMC = subset(SMC_FB, ident = "SMC")
counts = as.matrix(SMC@assays$RNA@counts)
meta = SMC@meta.data
sce = SingleCellExperiment(assays = list(counts = counts))
filter = rowSums(assay(sce)>0) > dim(counts)[2]*0.25
gene.name = c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "LRP1", "SERPINE1")
filter[which (names(filter) %in% gene.name)] = TRUE
sce_tmp = sce[filter, ]
zinb = zinbFit(sce_tmp, K = 2, epsilon = 1000, BPPARAM = BiocParallel::bpparam())
saveRDS(zinb, "zinb_hm_Tgfb.rds")
sce_zinb = zinbwave(sce_tmp, fitted_model = zinb, K = 2, epsilon = 1000, observationalWeights = TRUE)
weights = assay(sce_zinb, "weights")
dge = DGEList(assay(sce_zinb))
dge = calcNormFactors(dge)
meta = meta %>% filter(rownames(meta) %in% colnames(dge))
design = model.matrix(~meta$stim)
dge$weights = weights
dge = estimateDisp(dge, design)

tmm = edgeR::cpm(dge, log = T)
tmm = tmm[rownames(tmm) %in% gene.name,]
tmm = round(tmm, 10)
tmm = tmm - min(tmm)
SMC = subset(SMC_FB, ident = "SMC")

zlmfit = zlm(~stim + age + (1|stim:id), scaRaw, method='glmer', ebayes=FALSE)
zlmfit = zlm(~stim + age, scaRaw)

sumzlmfit = summary(zlm(~stim + age, scaRaw), doLRT = 'stimTAA')
print(sumzlmfit, n = length(gene.name))
summaryDt <- sumzlmfit$datatable
data.frame(summaryDt)

fcHurdle <- merge(summaryDt[contrast=='stimTAA' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='stimTAA' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
print(fcHurdle[fcHurdle$primerid %in% gene.name,])

### violin plot
tmm = edgeR::cpm(dge, log = T)
tmm = tmm[rownames(tmm) %in% gene.name,]
dat_tmm = rbind(data.frame(cell = colnames(tmm), tmm = tmm[rownames(tmm) == "TGFB1",], gene = "TGFB1"),
                data.frame(cell = colnames(tmm), tmm = tmm[rownames(tmm) == "TGFB2",], gene = "TGFB2"),
                data.frame(cell = colnames(tmm), tmm = tmm[rownames(tmm) == "TGFB3",], gene = "TGFB3"),
                data.frame(cell = colnames(tmm), tmm = tmm[rownames(tmm) == "TGFBR1",], gene = "TGFBR1"),
                data.frame(cell = colnames(tmm), tmm = tmm[rownames(tmm) == "TGFBR2",], gene = "TGFBR2"),
                data.frame(cell = colnames(tmm), tmm = tmm[rownames(tmm) == "LRP1",], gene = "LRP1"))

meta$cell =rownames(meta) 
dat_tmm = dat_tmm %>% inner_join(meta, by = "cell")
dat_tmm$gene = factor(dat_tmm$gene, c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "LRP1"))

ggplot(data = dat_tmm, aes(x = stim,  y = tmm))+
  geom_violin(aes(color = stim, fill = stim))+
  geom_jitter(width = 0.05, size = 0.1, color = "grey60")+
  facet_wrap(~ gene)+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank())

# --- Supplemental Figure 6B --- #
FB = subset(SMC_FB, ident = "FB")
counts = as.matrix(FB@assays$RNA@counts)
meta = FB@meta.data
sce = SingleCellExperiment(assays = list(counts = counts))
filter = rowSums(assay(sce)>0) > dim(counts)[2]*0.25
gene.name = c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "LRP1", "SERPINE1")
filter[which (names(filter) %in% gene.name)] = TRUE
sce_tmp = sce[filter, ]
zinb = zinbFit(sce_tmp, K = 2, epsilon = 1000, BPPARAM = BiocParallel::bpparam())
sce_zinb = zinbwave(sce_tmp, fitted_model = zinb, K = 2, epsilon = 1000, observationalWeights = TRUE)
weights = assay(sce_zinb, "weights")
dge = DGEList(assay(sce_zinb))
dge = calcNormFactors(dge)
meta = meta %>% filter(rownames(meta) %in% colnames(dge))
design = model.matrix(~meta$stim)
dge$weights = weights
dge = estimateDisp(dge, design)

tmm = edgeR::cpm(dge, log = T)
tmm = tmm[rownames(tmm) %in% gene.name,]
tmm = round(tmm, 10)
tmm = tmm - min(tmm)

sumzlmfit = summary(zlm(~stim + age, scaRaw), doLRT = 'stimTAA')
print(sumzlmfit, n = length(gene.name))
summaryDt <- sumzlmfit$datatable
data.frame(summaryDt)

fcHurdle <- merge(summaryDt[contrast=='stimTAA' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='stimTAA' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
print(fcHurdle)
print(fcHurdle[fcHurdle$primerid %in% gene.name,])

### violin plot
tmm = edgeR::cpm(dge, log = T)
tmm = tmm[rownames(tmm) %in% gene.name,]
dat_tmm = rbind(data.frame(cell = colnames(tmm), tmm = tmm[rownames(tmm) == "TGFB1",], gene = "TGFB1"),
                data.frame(cell = colnames(tmm), tmm = tmm[rownames(tmm) == "TGFB2",], gene = "TGFB2"),
                data.frame(cell = colnames(tmm), tmm = tmm[rownames(tmm) == "TGFB3",], gene = "TGFB3"),
                data.frame(cell = colnames(tmm), tmm = tmm[rownames(tmm) == "TGFBR1",], gene = "TGFBR1"),
                data.frame(cell = colnames(tmm), tmm = tmm[rownames(tmm) == "TGFBR2",], gene = "TGFBR2"),
                data.frame(cell = colnames(tmm), tmm = tmm[rownames(tmm) == "LRP1",], gene = "LRP1"))

meta$cell =rownames(meta) 
dat_tmm = dat_tmm %>% inner_join(meta, by = "cell")
dat_tmm$gene = factor(dat_tmm$gene, c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "LRP1"))

ggplot(data = dat_tmm, aes(x = stim,  y = tmm))+
  geom_violin(aes(color = stim, fill = stim))+
  geom_jitter(width = 0.05, size = 0.1, color = "grey60")+
  facet_wrap(~ gene)+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank())
