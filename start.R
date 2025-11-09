library(Seurat)
library(dplyr)


no2 <- no1

if(any(duplicated(rownames(no2)))){
  rownames(no2) <- make.unique(rownames(no2))
}
no2[["nCount_RNA_log10"]] <- log10(no2$nCount_RNA + 1)
no2[["nFeature_RNA_log10"]] <- log10(no2$nFeature_RNA + 1)

mt_genes <- rownames(no2)[grep("^MT-", rownames(no2), ignore.case = TRUE)]
if(length(mt_genes) == 0){
  warning(" 未检测到 MT 基因，将跳过 percent.mt 过滤")
  no2[["percent.mt"]] <- 0
} else {
  no2[["percent.mt"]] <- PercentageFeatureSet(no2, features = mt_genes)
}
if(!"batch" %in% colnames(no2@meta.data)){
  if("orig.ident" %in% colnames(no2@meta.data)){
    no2$batch <- no2$orig.ident
  } else {
    no2$batch <- "batch1"
  }
}

qc_cells <- WhichCells(no2, expression = nCount_RNA <= 1.2e5 &
                         nFeature_RNA >= 150 & nFeature_RNA <= 10000 &
                         percent.mt <= 20)
cat("放宽 QC 后通过的细胞数: ", length(qc_cells), "\n")
no2 <- subset(no2, cells = qc_cells)
no2 <- NormalizeData(no2, normalization.method = "RC", scale.factor = 1e4)
cat("已完成文库大小归一化\n")
cat(" 放宽 QC 后 Seurat 对象信息")
print(no2)
cat("细胞数: ", ncol(no2), " 基因数: ", nrow(no2), "\n")
library(Seurat)
library(dplyr)
no4 <- FindVariableFeatures(no4, selection.method = "vst", nfeatures = 3000)
no4 <- ScaleData(no4, features = VariableFeatures(no4))
no4 <- RunPCA(no4, features = VariableFeatures(no4), npcs = 50)
no4 <- FindNeighbors(no4, dims = 1:50, k.param = 10)
no4 <- FindClusters(no4, resolution = 0.5)
markers <- FindAllMarkers(no4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
markers_sig <- markers[markers$p_val_adj < 0.05 & markers$avg_log2FC > 1, ]
top10_markers <- markers_sig %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(no4, features = top10_markers$gene) + NoLegend()
assign_CAF_subtypes <- function(caf_cells, mat){
 pos_myCAF <- c("ACTA2","TAGLN","MYH11","CNN1")
  pos_iCAF  <- c("IL6","CXCL12","PDGFRA","CXCL1")
  pos_apCAF <- c("CD74","HLA-DRA","HLA-DRB1")
  neg_immune <- c("CD3D","MS4A1","CD1C")
  avg_expr <- function(genes) colMeans(mat[rownames(mat) %in% genes, , drop=FALSE])
  myCAF_score <- avg_expr(pos_myCAF)
  iCAF_score  <- avg_expr(pos_iCAF)
  apCAF_score <- avg_expr(pos_apCAF) - avg_expr(neg_immune)
  subtype <- rep("CAF_general", length(caf_cells))
  names(subtype) <- caf_cells
  subtype[myCAF_score > iCAF_score & myCAF_score > apCAF_score] <- "myCAF"
  subtype[iCAF_score > myCAF_score & iCAF_score > apCAF_score] <- "iCAF"
  subtype[apCAF_score > myCAF_score & apCAF_score > iCAF_score] <- "apCAF"
  return(subtype)
}

cluster2celltype <- c(
  "0" = "T_cells",
  "1" = "CD8_T_cells",
  "2" = "B_cells",
  "3" = "CAF_general",
  "4" = "iCAF",
  "5" = "apCAF",
  "6" = "Plasma_cells",
  "7" = "NK_cells",
  "8" = "CAF_general",
  "9" = "CD4_T_cells",
  "10"= "Tregs")
no4$predicted_celltype <- cluster2celltype[as.character(no4$seurat_clusters)]
caf_cells <- names(no4$predicted_celltype[no4$predicted_celltype %in% c("CAF_general","iCAF","apCAF","myCAF")])
caf_mat <- GetAssayData(no4, slot="data")[, caf_cells]
caf_subtypes <- assign_CAF_subtypes(caf_cells, caf_mat)
no4$predicted_celltype[caf_cells] <- caf_subtypes
table(no4$predicted_celltype)[order(-table(no4$predicted_celltype))]
no4 <- RunUMAP(no4, dims = 1:50)
DimPlot(no4, group.by = "predicted_celltype", label = TRUE, repel = TRUE)
markers_sig <- markers[markers$p_val_adj < 0.05 & markers$avg_log2FC > 1, ]
top50_markers <- markers_sig %>% 
  group_by(cluster) %>% 
  top_n(50, avg_log2FC)
no5withmarkers <- no5
no5withmarkers@misc$markers <- markers
no5withmarkers@misc$markers_sig <- markers_sig
no5withmarkers@misc$top50_markers <- top50_markers
saveRDS(no5withmarkers, "no5withmarkers.rds")
cluster_gene_lists <- split(top50_markers$gene, top50_markers$cluster)
marker_db <- list(
  Hepatocyte = c("ALB", "APOA1", "APOC3", "APOB", "KRT8", "KRT18", "GPC3", "EPCAM"),
  Malignant = c("KRT19", "EPCAM", "MDK", "LGALS3", "AFP", "SPP1"),
  T_NK = c("CD3D", "CD3E", "CD2", "CD7", "TRAC", "NKG7", "GNLY", "GZMB","KLRD1"),
  Bcell = c("MS4A1", "CD79A", "CD79B", "MZB1", "IGHM", "IGKC", "TNFRSF13B"),
  Plasma = c("MZB1", "JCHAIN", "IGKC", "DERL3", "SDC1"),
  Monocyte_Macro = c("LYZ","FCGR3A","CSF1R","S100A8","S100A9","CD68"),
  DC = c("LILRA4","IRF8","CLEC9A","XCR1","ITGAX","ITGAM"),
  Endothelial = c("PECAM1","VWF","ENG","KDR","ESAM"),
  Fibroblast_CAF = c("COL1A1","COL1A2","COL3A1","ACTA2","TAGLN","PDPN","FAP")
)
cluster_annotation <- sapply(cluster_gene_lists, function(genes){
  scores <- sapply(marker_db, function(ref) sum(ref %in% genes))
  names(which.max(scores))
})
cluster_annotation
myCAF_genes <- c("ACTA2","TAGLN","MYH11","CNN1")
iCAF_genes  <- c("IL6","CXCL12","PDGFRA","CXCL1")
apCAF_genes <- c("CD74","HLA-DRA","HLA-DRB1")
no5withmarkers <- AddModuleScore(no5withmarkers, features = list(myCAF_genes), name = "score_myCAF")
no5withmarkers <- AddModuleScore(no5withmarkers, features = list(iCAF_genes), name = "score_iCAF")
no5withmarkers <- AddModuleScore(no5withmarkers, features = list(apCAF_genes), name = "score_apCAF")
library(dplyr)
cluster_scores <- no5withmarkers@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(myCAF = mean(score_myCAF1, na.rm=TRUE),
            iCAF  = mean(score_iCAF1, na.rm=TRUE),
            apCAF = mean(score_apCAF1, na.rm=TRUE),
            n = n())
print(cluster_scores)

library(Seurat)
library(dplyr)
marker_db <- list(
  Hepatocyte = c("ALB", "APOA1", "APOC3", "APOB", "KRT8", "KRT18", "GPC3", "EPCAM"),
  Malignant = c("KRT19", "EPCAM", "MDK", "LGALS3", "AFP", "SPP1"),
  T_NK = c("CD3D", "CD3E", "CD2", "CD7", "TRAC", "NKG7", "GNLY", "GZMB","KLRD1"),
  Bcell = c("MS4A1", "CD79A", "CD79B", "MZB1", "IGHM", "IGKC", "TNFRSF13B"),
  Plasma = c("MZB1", "JCHAIN", "IGKC", "DERL3", "SDC1"),
  Monocyte_Macro = c("LYZ","FCGR3A","CSF1R","S100A8","S100A9","CD68"),
  DC = c("LILRA4","IRF8","CLEC9A","XCR1","ITGAX","ITGAM"),
  Endothelial = c("PECAM1","VWF","ENG","KDR","ESAM"),
  Fibroblast_CAF = c("COL1A1","COL1A2","COL3A1","ACTA2","TAGLN","PDPN","FAP")
)
myCAF_genes <- c("ACTA2","TAGLN","MYH11","CNN1")
iCAF_genes  <- c("IL6","CXCL12","PDGFRA","CXCL1")
apCAF_genes <- c("CD74","HLA-DRA","HLA-DRB1")

no5withmarkers <- AddModuleScore(no5withmarkers, list(myCAF_genes), name="score_myCAF")
no5withmarkers <- AddModuleScore(no5withmarkers, list(iCAF_genes), name="score_iCAF")
no5withmarkers <- AddModuleScore(no5withmarkers, list(apCAF_genes), name="score_apCAF")


CD4_genes <- c("CD4","IL7R")
CD8_genes <- c("CD8A","CD8B","GZMB","PRF1")
Treg_genes <- c("FOXP3","IL2RA","CTLA4")
Exhaustion_genes <- c("PDCD1","CTLA4","LAG3","HAVCR2")
no5withmarkers <- AddModuleScore(no5withmarkers, list(CD4_genes), name="score_CD4")
no5withmarkers <- AddModuleScore(no5withmarkers, list(CD8_genes), name="score_CD8")
no5withmarkers <- AddModuleScore(no5withmarkers, list(Treg_genes), name="score_Treg")
no5withmarkers <- AddModuleScore(no5withmarkers, list(Exhaustion_genes), name="score_exhaustion")
Bcell_genes <- c("MS4A1","CD79A","CD79B","MZB1","IGHM","IGKC","TNFRSF13B")
Plasma_genes <- c("MZB1","JCHAIN","IGKC","DERL3","SDC1")
no5withmarkers <- AddModuleScore(no5withmarkers, list(Bcell_genes), name="score_Bcell")
no5withmarkers <- AddModuleScore(no5withmarkers, list(Plasma_genes), name="score_Plasma")
cluster_annotation <- sapply(cluster_gene_lists, function(genes){
  scores <- sapply(marker_db, function(ref) sum(ref %in% genes))
  names(which.max(scores))
})
no5withmarkers@meta.data$celltype_big <- cluster_annotation[as.character(no5withmarkers@meta.data$seurat_clusters)]
no5withmarkers@meta.data <- no5withmarkers@meta.data %>%
  mutate(
     Tsubtype = case_when(
      celltype_big == "T_NK" & score_CD41 > score_CD81 & score_CD41 > score_Treg1 & score_CD41 > score_exhaustion1 ~ "CD4",
      celltype_big == "T_NK" & score_CD81 > score_CD41 & score_CD81 > score_Treg1 & score_CD81 > score_exhaustion1 ~ "CD8",
      celltype_big == "T_NK" & score_Treg1 > score_CD41 & score_Treg1 > score_CD81 & score_Treg1 > score_exhaustion1 ~ "Treg",
      celltype_big == "T_NK" & score_exhaustion1 > score_CD41 & score_exhaustion1 > score_CD81 & score_exhaustion1 > score_Treg1 ~ "Exhaustion",
      TRUE ~ NA_character_
    ),
    Bsubtype = case_when(
      celltype_big %in% c("Bcell","Plasma") & score_Bcell1 > score_Plasma1 ~ "Bcell",
      celltype_big %in% c("Bcell","Plasma") & score_Plasma1 > score_Bcell1 ~ "Plasma",
      TRUE ~ NA_character_
    ),
    CAF_subtype_final = case_when(
      celltype_big == "Fibroblast_CAF" & score_myCAF1 > score_iCAF1 & score_myCAF1 > score_apCAF1 ~ "myCAF",
      celltype_big == "Fibroblast_CAF" & score_iCAF1 > score_myCAF1 & score_iCAF1 > score_apCAF1 ~ "iCAF",
      celltype_big == "Fibroblast_CAF" & score_apCAF1 > score_myCAF1 & score_apCAF1 > score_iCAF1 ~ "apCAF",
      celltype_big == "Fibroblast_CAF" ~ "CAF_other",
      TRUE ~ NA_character_
    )
  )
no5withmarkers@meta.data <- no5withmarkers@meta.data %>%
  mutate(
    celltype_final = case_when(
      celltype_big == "Fibroblast_CAF" ~ "CAF",
      celltype_big == "T_NK" ~ ifelse(!is.na(Tsubtype), Tsubtype, "T_other"),
      celltype_big %in% c("Bcell","Plasma") ~ ifelse(!is.na(Bsubtype), Bsubtype, "B_other"),
      TRUE ~ celltype_big
    )
  )
no6 <- no5withmarkers
save(no6, file="no6.RData")
write.csv(no6@meta.data, file="no6_meta.csv", row.names=TRUE)
no6@meta.data %>%
  group_by(celltype_final) %>%
  summarise(n = n()) %>%
  arrange(desc(n))
no6@meta.data %>%
  filter(celltype_big == "Fibroblast_CAF") %>%
  group_by(CAF_subtype_final) %>%
  summarise(n = n()) %>%
  arrange(desc(n))


meta <- no5withmarkers@meta.data

meta$CAF_subtype <- NA
meta$Tsubtype   <- NA
meta$Bsubtype   <- NA
CAF_clusters <- names(cluster_annotation)[cluster_annotation=="Fibroblast_CAF"]
meta$CAF_subtype[meta$seurat_clusters %in% CAF_clusters] <- apply(
  meta[meta$seurat_clusters %in% CAF_clusters, c("score_myCAF1","score_iCAF1","score_apCAF1")],
  1,
  function(x) names(which.max(x))
)
TNK_clusters <- names(cluster_annotation)[cluster_annotation=="T_NK"]
meta$Tsubtype[meta$seurat_clusters %in% TNK_clusters] <- apply(
  meta[meta$seurat_clusters %in% TNK_clusters, c("score_CD41","score_CD81","score_Treg1","score_exhaustion1")],
  1,
  function(x) names(which.max(x))
  B_clusters <- names(cluster_annotation)[cluster_annotation %in% c("Bcell","Plasma")]
  meta$Bsubtype[meta$seurat_clusters %in% B_clusters] <- apply(
  meta[meta$seurat_clusters %in% B_clusters, c("score_Bcell1","score_Plasma1")],
  1,
  function(x) names(which.max(x))
))

no5withmarkers@meta.data <- meta
library(dplyr)
CAF_stats <- meta %>% filter(!is.na(CAF_subtype)) %>%
  group_by(seurat_clusters, CAF_subtype) %>%
  summarise(n=n()) %>%
  mutate(freq=n/sum(n))
T_stats <- meta %>% filter(!is.na(Tsubtype)) %>%
  group_by(seurat_clusters, Tsubtype) %>%
  summarise(n=n()) %>%
  mutate(freq=n/sum(n))
B_stats <- meta %>% filter(!is.na(Bsubtype)) %>%
  group_by(seurat_clusters, Bsubtype) %>%
  summarise(n=n()) %>%
 mutate(freq=n/sum(n))
 CAF_stats
 T_stats
 B_stats

library(dplyr)
no5withmarkers@meta.data %>%
  group_by(celltype_big) %>%
  summarise(n = n()) %>%
  arrange(desc(n))
no5withmarkers@meta.data %>%
  group_by(celltype_final) %>%
  summarise(n = n()) %>%
  arrange(desc(n))
no5withmarkers@meta.data %>%
  filter(celltype_big == "Fibroblast_CAF") %>%
  group_by(CAF_subtype_final) %>%
  summarise(n = n()) %>%
  arrange(desc(n))
no5withmarkers@meta.data %>%
  filter(celltype_big == "T_NK") %>%
  group_by(celltype_final) %>%
  summarise(n = n()) %>%
  arrange(desc(n))
no5withmarkers@meta.data %>%
  filter(celltype_big %in% c("Bcell","Plasma")) %>%
  group_by(celltype_final) %>%
  summarise(n = n()) %>%
  arrange(desc(n))


library(dplyr)
library(ggplot2)
library(scales)


lymph_cells <- no6_with_celltype@meta.data %>%
  filter(celltype_big %in% c("T_NK", "Bcell", "Plasma")) %>%
  group_by(celltype_final) %>%
  summarise(n = n()) %>%
  ungroup()
lymph_colors <- c(
  "CD41" = "#1f78b4",        # T CD4
  "CD81" = "#33a02c",        # T CD8
  "Treg1" = "#e31a1c",       # Treg
  "exhaustion1" = "#ff7f00", # T Exhaustion
  "Bcell1" = "#6a3d9a",      # Bcell
  "Plasma1" = "#b15928"      # Plasma
)
ggplot(lymph_cells, aes(x = celltype_final, y = n, fill = celltype_final)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = lymph_colors) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(x = "Lymphocyte Subtype", y = "Cell Number")
caf_cells <- no6_with_celltype@meta.data %>%
  filter(celltype_big == "Fibroblast_CAF") %>%
  group_by(CAF_subtype_final) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(prop = n / sum(n),
         label = paste0(CAF_subtype_final, " (", n, ")"))
ggplot(caf_cells, aes(x = "", y = prop, fill = CAF_subtype_final)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("myCAF"="#1f78b4", "iCAF"="#33a02c", "apCAF"="#e31a1c")) +
  theme_void(base_size = 14) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4)

library(dplyr)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
t_cells <- subset(obj, subset = celltype_big == "T_NK")

t_subtype_counts <- t_cells@meta.data %>%
  filter(!is.na(celltype_final)) %>%
  group_by(celltype_final) %>%
  summarise(count = n())
p_T <- ggplot(t_subtype_counts, aes(x = "", y = count, fill = celltype_final)) +
  geom_bar(stat = "identity", width = 1, color="white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = paste0(celltype_final, "\n", count)),
            position = position_stack(vjust = 0.5), size = 3) +
  theme_void() +
  scale_fill_brewer(palette = "Set3")
b_cells <- subset(obj, subset = celltype_big == "Bcell")

b_subtype_counts <- b_cells@meta.data %>%
  filter(!is.na(Bsubtype)) %>%
  group_by(Bsubtype) %>%
  summarise(count = n())

p_B <- ggplot(b_subtype_counts, aes(x = "", y = count, fill = Bsubtype)) +
  geom_bar(stat = "identity", width = 1, color="white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = paste0(Bsubtype, "\n", count)),
            position = position_stack(vjust = 0.5), size = 3) +
  theme_void() +
  scale_fill_brewer(palette = "Set3")
caf_cells <- subset(obj, subset = celltype_big == "Fibroblast_CAF")

caf_subtype_counts <- caf_cells@meta.data %>%
  filter(!is.na(CAF_subtype_final)) %>%
  group_by(CAF_subtype_final) %>%
  summarise(count = n())

p_CAF <- ggplot(caf_subtype_counts, aes(x = "", y = count, fill = CAF_subtype_final)) +
  geom_bar(stat = "identity", width = 1, color="white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = paste0(CAF_subtype_final, "\n", count)),
            position = position_stack(vjust = 0.5), size = 3) +
  theme_void() +
  scale_fill_brewer(palette = "Set3")
pdf("T_B_CAF_subtypes_piecharts.pdf", width = 12, height = 4)
gridExtra::grid.arrange(p_T, p_B, p_CAF, ncol = 3)
dev.off()

library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)  


table(obj@meta.data$celltype_big)

t_cells <- subset(obj, subset = celltype_big == "Tcell")
df_t <- data.frame(
  UMAP_1 = Embeddings(t_cells, "umap")[,1],
  UMAP_2 = Embeddings(t_cells, "umap")[,2],
  Subtype = t_cells@meta.data$Tsubtype  
)

p_t <- ggplot(df_t, aes(x=UMAP_1, y=UMAP_2, color=Subtype)) +
  geom_point(size=0.5, alpha=0.8) +
  theme_classic(base_size = 12) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right") +
  guides(color=guide_legend(title="T cell subtype"))
b_cells <- subset(obj, subset = celltype_big == "Bcell")
df_b <- data.frame(
  UMAP_1 = Embeddings(b_cells, "umap")[,1],
  UMAP_2 = Embeddings(b_cells, "umap")[,2],
  Subtype = b_cells@meta.data$Bsubtype)  

p_b <- ggplot(df_b, aes(x=UMAP_1, y=UMAP_2, color=Subtype)) +
  geom_point(size=0.5, alpha=0.8) +
  theme_classic(base_size = 12) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right") +
  guides(color=guide_legend(title="B cell subtype"))

caf_cells <- subset(obj, subset = celltype_big == "Fibroblast_CAF")
df_caf <- data.frame(
  UMAP_1 = Embeddings(caf_cells, "umap")[,1],
  UMAP_2 = Embeddings(caf_cells, "umap")[,2],
  Subtype = caf_cells@meta.data$CAF_subtype_final ) 

p_caf <- ggplot(df_caf, aes(x=UMAP_1, y=UMAP_2, color=Subtype)) +
  geom_point(size=0.5, alpha=0.8) +
  theme_classic(base_size = 12) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right") +
  guides(color=guide_legend(title="CAF subtype"))
p_all <- plot_grid(p_t, p_b, p_caf, ncol=3, align="hv")
pdf("UMAP_T_B_CAF_subtypes.pdf", width=15, height=5)
print(p_all)
dev.off()


library(ggplot2)

library(dplyr)

library(Seurat)
t_cells <- subset(obj, subset = celltype_big == "Tcell")
t_counts <- table(t_cells@meta.data$subtype) %>% as.data.frame()
colnames(t_counts) <- c("subtype","count")
b_cells <- subset(obj, subset = celltype_big == "Bcell")
b_counts <- table(b_cells@meta.data$subtype) %>% as.data.frame()
colnames(b_counts) <- c("subtype","count")
caf_cells <- subset(obj, subset = celltype_big == "Fibroblast_CAF")
caf_counts <- table(caf_cells@meta.data$CAF_subtype_final) %>% as.data.frame()
colnames(caf_counts) <- c("subtype","count")
plot_pie <- function(df){
  df <- df %>% 
    arrange(desc(count)) %>% 
    mutate(prop = count / sum(count),
           ypos = cumsum(prop) - 0.5*prop)
  
  ggplot(df, aes(x="", y=prop, fill=subtype)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    geom_text(aes(y=ypos, label=paste0(subtype,"\n",count)), color="black", size=3) +
    theme_void() +
    scale_fill_brewer(palette="Set3")  
}
p_t <- plot_pie(t_counts)
p_b <- plot_pie(b_counts)
p_caf <- plot_pie(caf_counts)
pdf("T_B_CAF_subcluster_pie.pdf", width=12, height=4)
print(p_t)
print(p_b)
print(p_caf)
dev.off()

library(dplyr)
library(ggplot2)
library(stringr)  

all_cells <- no6_with_celltype@meta.data %>%
  filter(celltype_final != "Hepatocyte") %>%
  mutate(celltype_clean = str_remove_all(celltype_final, "\\d+")) %>%  
  group_by(celltype_clean) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  arrange(n)

nature_colors <- c(
  "CD81" = "#33a02c",
  "Monocyte_Macro" = "#e31a1c",
  "Malignant" = "#ff7f0e",
  "Endothelial" = "#9467bd",
  "Plasma" = "#8c564b",
  "exhaustion" = "#d62728",
  "Treg" = "#2ca02c",
  "CD41" = "#1f78b4",
  "Bcell" = "#6a3d9a",
  "CAF" = "#17becf",
  "DC" = "#bcbd22"
)
ggplot(all_cells, aes(x = n, y = reorder(celltype_clean, n), fill = celltype_clean)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = nature_colors) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(x = "Number of Cells")

