
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

obj <- readRDS("no6_with_celltype.rds")


caf_cells <- WhichCells(obj, expression = celltype_big == "Fibroblast_CAF")
obj_caf <- subset(obj, cells = caf_cells)


Idents(obj_caf) <- obj_caf@meta.data$CAF_subtype_final
caf_subtypes <- unique(obj_caf@meta.data$CAF_subtype_final)


deg_list <- list()
for(sub in caf_subtypes){
  deg <- FindMarkers(
    obj_caf,
    ident.1 = sub,
    ident.2 = NULL,   
    assay = "RNA",
    slot = "data",
    logfc.threshold = 0.25,
    min.pct = 0.1
  )
  deg <- deg %>% rownames_to_column("gene") %>% mutate(subtype = sub)
  deg_list[[sub]] <- deg
}
caf_deg_all <- bind_rows(deg_list)

gene_list <- caf_deg_all %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
  pull(gene)

gene_entrez <- bitr(gene_list, fromType="SYMBOL",
                    toType="ENTREZID", OrgDb="org.Hs.eg.db")

go_res <- enrichGO(
  gene         = gene_entrez$ENTREZID,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  pAdjustMethod= "BH",
  qvalueCutoff = 0.05,
  readable     = TRUE
)

tcell_go <- go_res@result %>%
  filter(
    grepl("T cell|T-cell|lymphocyte|CD4|CD8|Treg|B cell|B-cell", Description),
    p.adjust < 0.05
  ) %>%
  arrange(p.adjust) %>%
  head(20)   # 前20个 GO term

ecm_go <- go_res@result %>%
  filter(
    grepl("ECM|extracellular matrix|collagen|fibronectin|integrin|basement membrane|matrix organization|cell adhesion", Description, ignore.case = TRUE),
    p.adjust < 0.05
  ) %>%
  arrange(p.adjust) %>%
  head(20)

plot_go <- function(df, title_text, output_file){
  p <- ggplot(df, aes(x=reorder(Description, -Count), y=Count,
                      size=FoldEnrichment, color=p.adjust)) +
    geom_point() +
    coord_flip() +
    scale_color_gradient(low="red", high="blue") 
    theme_bw(base_size = 14) +
    labs(
      title = title_text,
      x = NULL, y = "Gene Count",
      color = "Adjusted p-value", size = "Fold Enrichment"
    ) +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title  = element_text(size = 16, face="bold", hjust=0.5)
    )
  
  ggsave(output_file, p, width = 15, height = 15)
  message("已保存图 ", output_file)
}

plot_go(tcell_go, "CAF DEG enriched T/B cell GO terms (top 20, p<0.05)", "CAF_GO_TBcell_top20_p05.pdf")
plot_go(ecm_go, "CAF DEG enriched ECM-related GO terms (top 20, p<0.05)", "CAF_GO_ECM_top20_p05.pdf")
library(Seurat)
library(dplyr)
library(pheatmap)
library(tidyr)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
obj <- readRDS("no6_with_celltype.rds")
caf_cells <- WhichCells(obj, expression = celltype_big == "Fibroblast_CAF")
obj_caf <- subset(obj, cells = caf_cells)
Idents(obj_caf) <- obj_caf@meta.data$CAF_subtype_final
caf_subtypes <- unique(obj_caf@meta.data$CAF_subtype_final)
deg_list <- list()
for(sub in caf_subtypes){
  deg <- FindMarkers(
    obj_caf,
    ident.1 = sub,
    ident.2 = NULL,   
    assay = "RNA",
    slot = "data",
    logfc.threshold = 0.25,
    min.pct = 0.1
  )
  deg <- deg %>% rownames_to_column("gene") %>% mutate(subtype = sub)
  deg_list[[sub]] <- deg
}
caf_deg_all <- bind_rows(deg_list)
gene_list <- caf_deg_all %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
  pull(gene)

gene_entrez <- bitr(gene_list, fromType="SYMBOL",
                    toType="ENTREZID", OrgDb="org.Hs.eg.db")
go_res <- enrichGO(
  gene         = gene_entrez$ENTREZID,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  pAdjustMethod= "BH",
  qvalueCutoff = 0.05,
  readable     = TRUE
)
tcell_go <- go_res@result %>%
  filter(
    grepl("T cell|T-cell|lymphocyte|CD4|CD8|Treg|B cell|B-cell", Description),
    p.adjust < 0.05
  ) %>%
  arrange(p.adjust)   

if(nrow(tcell_go) == 0) stop("放宽条件。")
ggplot(tcell_go, aes(x=reorder(Description, -Count), y=Count,
                     size=FoldEnrichment, color=p.adjust)) +
  geom_point() +
  coord_flip() +
  scale_color_gradient(low="red", high="blue") +
  theme_bw(base_size = 14) +
  labs(
    title="CAF DEG enriched T/B cell GO terms",
    x=NULL, y="Gene Count",
    color="Adjusted p-value", size="Fold Enrichment"
  ) +
  theme(
    axis.text.y = element_text(size=10),
    plot.title = element_text(size=15, face="bold", hjust=0.5)
  )
ggsave("CAF_GO_TBcell_p05_redLow_blueHigh.pdf",
       width=10, height=8)


library(Seurat)
library(dplyr)
for (celltype in names(deg_results_list)) {
  deg <- deg_results_list[[celltype]]
  
  deg_filtered <- deg %>% 
    filter(adj.P.Val < 0.05, abs(logFC) > 0.25) %>%
    arrange(desc(logFC))
  
  top_genes_raw <- head(rownames(deg_filtered), 20)
  top_genes <- sub("^[^_]+_", "", top_genes_raw)
  top_genes <- top_genes[top_genes %in% rownames(filtered_seurat_obj)]
  
  if(length(top_genes) == 0) next
  
  cells_use <- WhichCells(filtered_seurat_obj, idents = celltype)
  
  p <- DoHeatmap(filtered_seurat_obj, features = top_genes, cells = cells_use) + 
    ggtitle(paste("DEG Heatmap -", celltype))
  
  print(p)
  ggsave(filename = paste0("heatmap_", celltype, ".pdf"), plot = p, width = 8, height = 6)
}

library(Seurat)
library(dplyr)
library(limma)
celltypes <- unique(filtered_seurat_obj$celltype_annotated)
expr_mat <- GetAssayData(filtered_seurat_obj, slot = "data")
treatment <- filtered_seurat_obj$treatment_status
deg_results_list <- list()
for (celltype in celltypes) 
  message("Processing cell type: ", celltype)
  cells_use <- WhichCells(filtered_seurat_obj, ident = NULL, expression = celltype_annotated == celltype)
  expr_subset <- expr_mat[, cells_use]
  group <- treatment[cells_use]
  if (length(unique(group)) < 2) {
    message("  Skipping ", celltype, ": only one treatment group present.")
    next
  }
  
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(factor(group))
  fit <- lmFit(expr_subset, design)
  contrast <- makeContrasts(Chemotherapy - Naive, levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ReactomePA)
  library(ggsci)
  library(cowplot)
  library(data.table)
  obj <- readRDS("~/no6_with_celltype.rds")
  dir.create("sc_results_CAF", showWarnings = FALSE)
  caf_subclusters <- na.omit(unique(obj@meta.data$CAF_subtype_final))
  run_caf_analysis <- function(subclusters){
    all_top50_genes <- c()
    
    combs <- combn(subclusters, 2, simplify = FALSE)
    
    for(cmp in combs){
      ident1 <- cmp[1]
      ident2 <- cmp[2]
      message("🔹 Comparing ", ident1, " vs ", ident2, " in CAF")
      
      deg <- FindMarkers(obj,
                         ident.1 = ident1,
                         ident.2 = ident2,
                         group.by = "CAF_subtype_final",
                         min.pct = 0.1,
                         logfc.threshold = 0.25)
      if(nrow(deg) == 0) next
      
      deg <- deg %>% rownames_to_column("gene")
      top50 <- deg %>% arrange(desc(avg_log2FC)) %>% head(50)
      
    
      write.csv(top50,
                paste0("sc_results_CAF/CAF_", ident1, "_vs_", ident2, "_Top50_DEG.csv"),
                row.names = FALSE)
      
  
      p <- ggplot(top50, aes(x=reorder(gene, avg_log2FC), y=avg_log2FC)) +
        geom_bar(stat="identity", fill = pal_npg("nrc")(10)[2]) +
        coord_flip() +
        theme_classic(base_size=12) +
        labs(title=paste0("CAF ", ident1, " vs ", ident2, " Top50 DEG"),
             x="Gene", y="Avg log2FC") +
        theme(plot.title = element_text(hjust=0.5, size=14))
      
      ggsave(paste0("sc_results_CAF/CAF_", ident1, "_vs_", ident2, "_Top50_DEG.png"),
             p, width=6, height=5)
      
      all_top50_genes <- c(all_top50_genes, top50$gene)
    }
    
    return(unique(all_top50_genes))
  }
  caf_top50_genes <- run_caf_analysis(caf_subclusters)
  caf_entrez <- bitr(caf_top50_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  go_res <- enrichGO(caf_entrez$ENTREZID, OrgDb=org.Hs.eg.db,
                     ont="BP", pvalueCutoff=0.05)
  write.csv(go_res, "sc_results_CAF/CAF_GO_BP_Top50.csv", row.names = FALSE)
  png("sc_results_CAF/CAF_GO_BP_dotplot.png", width=6, height=5, units="in", res=300)
  print(dotplot(go_res, showCategory=15) + theme_classic())
  dev.off()
  kegg_res <- enrichKEGG(caf_entrez$ENTREZID, organism="hsa", pvalueCutoff=0.05)
  write.csv(kegg_res, "sc_results_CAF/CAF_KEGG_Top50.csv", row.names = FALSE)
  png("sc_results_CAF/CAF_KEGG_dotplot.png", width=6, height=5, units="in", res=300)
  print(dotplot(kegg_res, showCategory=15) + theme_classic())
  dev.off()
  ecm_genes <- c("COL1A1","COL1A2","COL3A1","COL4A1","COL5A1","FN1","LAMA1","LAMC1")
  ecm_deg <- caf_top50_genes[caf_top50_genes %in% ecm_genes]
  
  message(" CAF Top50 中 ECM/Collagen 相关基因: ", paste(ecm_deg, collapse=", "))
  
  
  
  library(Seurat)
  library(dplyr)
  library(pheatmap)
  library(tidyr)
  library(readr)
  obj <- readRDS("no6_with_celltype.rds")
  ecm_genes <- c("COL1A1","COL1A2","COL3A1","COL4A1","COL5A1","FN1","LAMA1","LAMC1")
  caf_cells <- WhichCells(obj, expression = celltype_big == "Fibroblast_CAF")
  caf_subtypes <- obj@meta.data[caf_cells, "CAF_subtype_final"]
  expr_mat <- GetAssayData(obj, slot="data")[ecm_genes, caf_cells, drop=FALSE]
  avg_expr <- t(apply(expr_mat, 1, function(g) tapply(g, caf_subtypes, mean)))
  avg_expr_scaled <- t(scale(t(avg_expr)))
  breaksList <- seq(-2, 2, by=0.05)
  colors <- colorRampPalette(c("white","firebrick3"))(length(breaksList))
  pheatmap(
    avg_expr_scaled,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_colnames = TRUE,
    show_rownames = TRUE,
    main = "CAF Subtypes vs ECM/Collagen (Expr)",
    color = colors,
    breaks = breaksList,
    cellwidth = 60,
    cellheight = 30,
    fontsize_row = 12,
    fontsize_col = 14
  )
  caf_deg_list <- list(
    myCAF_vs_iCAF = read_csv("sc_results_CAF/myCAF_vs_iCAF_DEG.csv"),
    myCAF_vs_apCAF = read_csv("sc_results_CAF/myCAF_vs_apCAF_DEG.csv"),
    iCAF_vs_apCAF = read_csv("sc_results_CAF/iCAF_vs_apCAF_DEG.csv")
  )
  caf_ecm_summary <- lapply(names(caf_deg_list), function(comp){
    deg <- caf_deg_list[[comp]]
    deg %>% filter(gene %in% ecm_genes) %>%
      mutate(comparison = comp)
  }) %>% bind_rows()
  heatmat <- caf_ecm_summary %>%
    select(gene, comparison, avg_log2FC) %>%
    tidyr::pivot_wider(names_from = comparison, values_from = avg_log2FC) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  heatmat[is.na(heatmat)] <- 0
  breaksList_deg <- seq(-3, 3, by=0.05)
  colors_deg <- colorRampPalette(c("white","firebrick3"))(length(breaksList_deg))
  pheatmap(
    heatmat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_colnames = TRUE,
    show_rownames = TRUE,
    main = "CAF Subtypes vs ECM/Collagen (DEG avg_log2FC)",
    color = colors_deg,
    breaks = breaksList_deg,
    cellwidth = 60,
    cellheight = 30,
    fontsize_row = 12,
    fontsize_col = 14
  )

  library(dplyr)
  library(ggplot2)
  library(readr)
  caf_deg_list <- list(
    myCAF_vs_iCAF = read_csv("sc_results_CAF/myCAF_vs_iCAF_DEG.csv"),
    myCAF_vs_apCAF = read_csv("sc_results_CAF/myCAF_vs_apCAF_DEG.csv"),
    iCAF_vs_apCAF = read_csv("sc_results_CAF/iCAF_vs_apCAF_DEG.csv")
  )
  ecm_genes <- c("COL1A1","COL1A2","COL3A1","COL4A1","COL5A1","FN1","LAMA1","LAMC1")
  caf_ecm_summary <- lapply(names(caf_deg_list), function(comp){
    deg <- caf_deg_list[[comp]]
    deg %>% filter(gene %in% ecm_genes) %>%
      mutate(comparison = comp)
  }) %>% bind_rows()
  caf_effect <- data.frame(
    CAF = c("myCAF", "iCAF", "apCAF"),
    ECM_Strength = c(
      mean(caf_ecm_summary$avg_log2FC[caf_ecm_summary$comparison == "myCAF_vs_iCAF" &
                                        caf_ecm_summary$gene %in% ecm_genes], na.rm=TRUE),
      mean(caf_ecm_summary$avg_log2FC[caf_ecm_summary$comparison == "iCAF_vs_apCAF" &
                                        caf_ecm_summary$gene %in% ecm_genes], na.rm=TRUE),
      mean(caf_ecm_summary$avg_log2FC[caf_ecm_summary$comparison == "myCAF_vs_apCAF" &
                                        caf_ecm_summary$gene %in% ecm_genes], na.rm=TRUE)
    )
  )
  
 
  ggplot(caf_effect, aes(x = CAF, y = ECM_Strength, group=1)) +
    geom_line(color="firebrick3", size=1.2) +
    geom_point(color="firebrick3", size=4) +
    ylim(min(caf_effect$ECM_Strength)-0.5, max(caf_effect$ECM_Strength)+0.5) +
    labs(title="CAF Subtypes Promote ECM/Collagen Strength",
         x="CAF Subtype", y="Average log2FC of ECM Genes") +
    theme_minimal(base_size = 14)
  library(Seurat)
  library(dplyr)
  library(pheatmap)
  library(tidyr)
  library(readr)
  obj <- readRDS("no6_with_celltype.rds")
  
 
  caf_markers <- c(
    "IL6", "CXCL12", "CCL5", "CCL19", "CCL21", 
    "TGFB1", "PDPN", "FAP", "ACTA2", "VCAM1",
    "ICAM1", "SPP1", "MMP14", "PTGS2", "CX3CL1",
    "TNFSF4", "TNFSF9", "SELP", "IL7", "IL15"
  )
  caf_cells <- WhichCells(obj, expression = celltype_big == "Fibroblast_CAF")
  caf_subtypes <- obj@meta.data[caf_cells, "CAF_subtype_final"]
  expr_mat <- GetAssayData(obj, slot="data")[caf_markers, caf_cells, drop=FALSE]
  
  avg_expr <- t(apply(expr_mat, 1, function(g) tapply(g, caf_subtypes, mean)))
  
  avg_expr_scaled <- t(scale(t(avg_expr)))
 
  breaksList <- seq(-2, 2, by=0.05)
  colors <- colorRampPalette(c("white","firebrick3"))(length(breaksList))
  

  pheatmap(
    avg_expr_scaled,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_colnames = TRUE,
    show_rownames = TRUE,
    main = "CAF Subtypes vs CAF Marker Genes (Expr)",
    color = colors,
    breaks = breaksList,
    cellwidth = 34,
    cellheight = 18,
    fontsize_row = 10,
    fontsize_col = 8
  )
  caf_deg_list <- list(
    myCAF_vs_iCAF = read_csv("sc_results_CAF/myCAF_vs_iCAF_DEG.csv"),
    myCAF_vs_apCAF = read_csv("sc_results_CAF/myCAF_vs_apCAF_DEG.csv"),
    iCAF_vs_apCAF = read_csv("sc_results_CAF/iCAF_vs_apCAF_DEG.csv")
  )
  
  caf_marker_deg <- lapply(names(caf_deg_list), function(comp){
    deg <- caf_deg_list[[comp]]
    deg %>% filter(gene %in% caf_markers) %>%
      mutate(comparison = comp)
  }) %>% bind_rows()
  heatmat <- caf_marker_deg %>%
    select(gene, comparison, avg_log2FC) %>%
    tidyr::pivot_wider(names_from = comparison, values_from = avg_log2FC) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  heatmat[is.na(heatmat)] <- 0
  breaksList_deg <- seq(-3, 3, by=0.05)
  colors_deg <- colorRampPalette(c("white","firebrick3"))(length(breaksList_deg))
  pheatmap(
    heatmat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_colnames = TRUE,
    show_rownames = TRUE,
    main = "CAF Subtypes vs CAF Marker Genes (DEG avg_log2FC)",
    color = colors_deg,
    breaks = breaksList_deg,
    cellwidth = 34,
    cellheight = 18,
    fontsize_row = 10,
    fontsize_col = 8
  )
  

  library(Seurat)
  library(dplyr)
  library(pheatmap)
  library(tidyr)
  library(readr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  obj <- readRDS("no6_with_celltype.rds")
  
 
  caf_cells <- WhichCells(obj, expression = celltype_big == "Fibroblast_CAF")
  obj_caf <- subset(obj, cells = caf_cells)
  
 
  Idents(obj_caf) <- obj_caf@meta.data$CAF_subtype_final
  caf_subtypes <- unique(obj_caf@meta.data$CAF_subtype_final)
  deg_list <- list()
  
  for(sub in caf_subtypes){
    deg <- FindMarkers(
      obj_caf,
      ident.1 = sub,
      ident.2 = NULL,  
      assay = "RNA",
      slot = "data",
      logfc.threshold = 0.25,
      min.pct = 0.1
    )
    deg <- deg %>% rownames_to_column("gene") %>% mutate(subtype = sub)
    deg_list[[sub]] <- deg
  }
  
  caf_deg_all <- bind_rows(deg_list)
  
  gene_list <- caf_deg_all %>% filter(p_val_adj < 0.05 & avg_log2FC > 0) %>% pull(gene)
  
  gene_entrez <- bitr(gene_list, fromType="SYMBOL",
                      toType="ENTREZID", OrgDb="org.Hs.eg.db")
  go_res <- enrichGO(
    gene         = gene_entrez$ENTREZID,
    OrgDb        = org.Hs.eg.db,
    ont          = "BP",
    pAdjustMethod= "BH",
    qvalueCutoff = 0.05,
    readable     = TRUE
  )
  tcell_go <- go_res@result %>%
    filter(grepl("T cell|T-cell|lymphocyte|CD4|CD8|Treg|B cell|B-cell", Description)) %>%
    arrange(p.adjust) %>%
    head(20)
  ggplot(tcell_go, aes(x=reorder(Description, -Count), y=Count, 
                       size=FoldEnrichment, color=p.adjust)) +
    geom_point() +
    coord_flip() +
    scale_color_gradient(low="blue", high="red") +
    theme_bw() +
    labs(title="CAF DEG enriched T/B cell GO terms", x="", y="Gene Count") +
    theme(axis.text.y = element_text(size=10))

  ecm_go <- go_res@result %>%
    filter(grepl("extracellular matrix|collagen", Description)) %>%
    arrange(p.adjust) %>%
    head(20)
  ggplot(ecm_go, aes(x=reorder(Description, -Count), y=Count, 
                     size=FoldEnrichment, color=p.adjust)) +
    geom_point() +
    coord_flip() +
    scale_color_gradient(low="blue", high="red") +
    theme_bw() +
    labs(title="CAF DEG enriched ECM/Collagen GO terms", x="", y="Gene Count") +
    theme(axis.text.y = element_text(size=10))  
  
  
  write_csv(myCAF_deg, "sc_results_CAF/myCAF_vs_All_DEG.csv")
  write_csv(iCAF_deg, "sc_results_CAF/iCAF_vs_All_DEG.csv")
  write_csv(apCAF_deg, "sc_results_CAF/apCAF_vs_All_DEG.csv")
  
  pheatmap(
    avg_expr_scaled,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_colnames = TRUE,
    show_rownames = TRUE,
    main = "CAF Subtypes vs CAF Marker Genes (Expr)",
    color = colors,
    breaks = breaksList,
    cellwidth = 34,
    cellheight = 18,
    fontsize_row = 10,
    fontsize_col = 8,
    filename = "figures/CAF_marker_expr_heatmap.pdf"  
  )
  

  pheatmap(
    heatmat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_colnames = TRUE,
    show_rownames = TRUE,
    main = "CAF Subtypes vs CAF Marker Genes (DEG avg_log2FC)",
    color = colors_deg,
    breaks = breaksList_deg,
    cellwidth = 34,
    cellheight = 18,
    fontsize_row = 10,
    fontsize_col = 8,
    filename = "figures/CAF_marker_deg_heatmap.pdf"  
  )
  
  library(dplyr)
  library(pheatmap)
  library(readr)
  library(tidyr)
  library(RColorBrewer)
  iCAF_deg <- read_csv("sc_results_CAF/iCAF_vs_All_DEG.csv")
  myCAF_deg <- read_csv("sc_results_CAF/myCAF_vs_All_DEG.csv")
  apCAF_deg <- read_csv("sc_results_CAF/apCAF_vs_All_DEG.csv")
  ecm_genes <- c("TGFB1","TGFB3","FN1","COL1A1","COL1A2",
                 "COL4A1","COL4A2","COL6A1","COL6A2","COL6A3",
                 "THBS1","THBS2","LAMA2","LAMA4","LAMA5",
                 "LAMB1","LAMB2","LAMC1","LAMC3","APP")
  iCAF_deg <- iCAF_deg %>% mutate(CAF = "iCAF") %>% select(gene, avg_log2FC, CAF)
  myCAF_deg <- myCAF_deg %>% mutate(CAF = "myCAF") %>% select(gene, avg_log2FC, CAF)
  apCAF_deg <- apCAF_deg %>% mutate(CAF = "apCAF") %>% select(gene, avg_log2FC, CAF)
  
  deg_all <- bind_rows(iCAF_deg, myCAF_deg, apCAF_deg) %>%
    filter(gene %in% ecm_genes)
  expr_mat <- deg_all %>%
    pivot_wider(names_from = CAF, values_from = avg_log2FC, values_fill = 0) %>%
    column_to_rownames("gene")
  
  
  expr_mat_scaled <- t(scale(t(expr_mat)))
  
  colors <- colorRampPalette(c("white", "yellow", "green"))(50)  
  pdf("CAF_ECM_Collagen_DEG_heatmap.pdf", width=6, height=8)
  pheatmap(expr_mat_scaled,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           color = colors,
           main = "CAF Subtypes ECM/Collagen DEGs",
           fontsize_row = 10,
           fontsize_col = 12)
  dev.off()
  
  
  
  
  library(Seurat)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  obj <- readRDS("no6_with_celltype.rds")
  caf_cells <- WhichCells(obj, expression = celltype_big == "Fibroblast_CAF")
  obj_caf <- subset(obj, cells = caf_cells)
  Idents(obj_caf) <- obj_caf@meta.data$CAF_subtype_final
  caf_subtypes <- unique(na.omit(obj_caf@meta.data$CAF_subtype_final))
  deg_list <- list()
  for(sub in caf_subtypes){
    deg <- FindMarkers(
      obj_caf,
      ident.1 = sub,
      ident.2 = NULL,
      assay = "RNA",
      slot = "data",
      logfc.threshold = 0.5,      
      min.pct = 0.15
    )
    deg <- deg %>% rownames_to_column("gene") %>% mutate(subtype = sub)
    deg_list[[sub]] <- deg
  }
  caf_deg_all <- bind_rows(deg_list)
  
  gene_list <- caf_deg_all %>%
    filter(p_val_adj < 0.01 & avg_log2FC > 0.5) %>%
    pull(gene)
  
  
  gene_entrez <- bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
 
  go_res <- enrichGO(
    gene         = gene_entrez$ENTREZID,
    OrgDb        = org.Hs.eg.db,
    ont          = "BP",
    pAdjustMethod= "BH",
    qvalueCutoff = 0.2,          
    readable     = TRUE
  )ecm_go <- go_res@result %>%
    filter(grepl("extracellular matrix|collagen|basement membrane|ECM organization|fibroblast", Description, ignore.case = TRUE)) %>%
    arrange(p.adjust) %>%
    head(20)
  pdf("CAF_ECM_GO_dotplot_improved.pdf", width=7, height=6)
  
  ggplot(ecm_go, aes(x=reorder(Description, -Count), y=Count, 
                     size=FoldEnrichment, color=-log10(p.adjust))) +
    geom_point() +
    coord_flip() +
    scale_color_gradient(low="#2166AC", high="#B2182B") +  
    theme_bw(base_size=13) +
    labs(title="CAF DEGs enriched in ECM/Collagen processes",
         x=NULL, y="Gene Count", color="-log10(adj.p)", size="Fold Enrichment") +
    theme(axis.text.y = element_text(size=9),
          plot.title = element_text(size=14, face="bold", hjust=0.5))
  
  dev.off()
  
  cat("one: CAF_ECM_GO_dotplot_improved.pdf generated\n")
  
  
  
  
  
  
  
  