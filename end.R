
library(Seurat)
library(CellChat)
library(dplyr)
library(ggplot2)
obj <- readRDS("no6_with_celltype.rds")
caf_subtypes <- c("myCAF","iCAF","apCAF")
tb_subtypes  <- c("CD81","CD41","exhaustion1","Treg1","Bcell1","Plasma1")
cells_keep <- WhichCells(obj, expression = CAF_subtype_final %in% caf_subtypes |
                           celltype_final %in% tb_subtypes)
obj_sub <- subset(obj, cells = cells_keep)

obj_sub$cellgroup <- ifelse(obj_sub$CAF_subtype_final %in% caf_subtypes,
                            obj_sub$CAF_subtype_final,
                            obj_sub$celltype_final)
cellchat <- createCellChat(object = obj_sub, group.by = "cellgroup")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
interaction_tb <- subsetCommunication(cellchat,
                                      sources = caf_subtypes,
                                      targets = tb_subtypes)
if(nrow(interaction_tb) == 0){
  warning("降低阈值或检查子簇名")
}
if(nrow(interaction_tb) > 0){
  netVisual_bubble(cellchat,
                   sources.use = caf_subtypes,
                   targets.use = tb_subtypes,
                   remove.isolate = FALSE)}

df_lr <- subsetCommunication(cellchat)
write.csv(df_lr, "CellChat_LR_interactions_all.csv", row.names = FALSE)
cat("CellChat_LR_interactions_all.csv\n")
interaction_all <- read.csv("CAF_to_TB_interactions_all.csv", stringsAsFactors = FALSE)
ecm_kw <- c("COLLAGEN", "ECM", "Collagen", "extracellular matrix", "ECM-Receptor", "FN1", "COL", "MMP", "LOX", "CD44", "ITGA", "ITGB")

interaction_all$ecm_related <- apply(interaction_all[, c("pathway_name","interaction_name","interaction_name_2","ligand","receptor","annotation")],
                                     1,
                                     function(row){
                                       any(sapply(ecm_kw, function(k) any(grepl(k, row, ignore.case = TRUE))))
                                     })

ecm_interactions <- subset(interaction_all, ecm_related == TRUE)
nrow(ecm_interactions)
head(ecm_interactions[, c("source","target","ligand","receptor","pathway_name","prob","annotation","evidence")])
table(ecm_interactions$pathway_name)

library(tidyverse)
library(igraph)
library(ggraph)
lr_data <- read.csv("CellChat_LR_interactions_all.csv", stringsAsFactors = FALSE)
head(lr_data)
edges <- lr_data %>%
  group_by(source, target) %>%
  summarise(weight = sum(prob, na.rm = TRUE)) %>%
  ungroup()

graph <- graph_from_data_frame(edges, directed = TRUE)
ggraph(graph, layout = "circle") +
  geom_edge_link(aes(width = weight), alpha = 0.7, color = "steelblue") +
  geom_node_point(size = 5, color = "tomato") +
  geom_node_text(aes(label = name), vjust = -1) +
  theme_void() +
  ggtitle("CellChat 通讯网络")



library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
CellChat_LR_interactions_all <- read_csv("CellChat_LR_interactions_all.csv")
caf_types <- c("CAF", "iCAF", "apCAF", "myCAF")
tb_types <- c("CD4", "CD8", "Treg", "Exhaustion", "Bcell", "Plasm")
ecm_targets <- c("TGFB1","TGFB3","FN1","COL1A1","COL1A2",
                 "COL4A1","COL4A2","COL6A1","COL6A2","COL6A3",
                 "THBS1","THBS2","LAMA2","LAMA4","LAMA5",
                 "LAMB1","LAMB2","LAMC1","LAMC3","APP")

caf_tb <- CellChat_LR_interactions_all %>%
  filter(source %in% caf_types & target %in% tb_types)
caf_tb_promote <- caf_tb %>%
  filter(str_detect(annotation, regex("promot|activat", ignore_case = TRUE)))
caf_tb_promote <- caf_tb_promote %>%
  mutate(ECM_related = ifelse(ligand %in% ecm_targets | receptor %in% ecm_targets, TRUE, FALSE))
ecm_stats <- caf_tb_promote %>%
  group_by(ECM_related) %>%
  summarise(
    count = n(),
    mean_prob = mean(prob, na.rm = TRUE),
    total_prob = sum(prob, na.rm = TRUE)
  ) %>%
  mutate(prop_count = count / sum(count),
         prop_prob = total_prob / sum(total_prob))

print(ecm_stats)
ggplot(ecm_stats, aes(x = ECM_related, y = prop_prob, fill = ECM_related)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#999999", "#1f78b4")) +
  labs(x = "是否为ECM相关通路", y = "影响占比（按prob加权）",
       title = "CAF→TB促进作用中ECM相关通路影响占比") +
  theme_minimal(base_size = 14)
ecm_gene_effect <- caf_tb_promote %>%
  filter(ECM_related) %>%
  group_by(ligand, target) %>%
  summarise(mean_prob = mean(prob, na.rm = TRUE)) %>%
  ungroup()

ggplot(ecm_gene_effect, aes(x = ligand, y = target, fill = mean_prob)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  labs(x = "CAF ECM配体", y = "T/B亚群", fill = "平均影响值 (prob)",
       title = "CAF ECM相关通路对T/B亚群的促进影响热图") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(forcats)
df <- read_csv("CellChat_LR_interactions_all.csv")
caf_subtypes <- c("myCAF","iCAF","apCAF")
tb_subtypes  <- c("CD4","CD8","Treg","exhaustion","Bcell","Plasma")

df <- df %>% filter(source %in% caf_subtypes & target %in% tb_subtypes)
promoting_keywords <- c("MIF","CD40","IL6","CXCL","CCL","ICOS")    
suppressive_keywords <- c("TGFB","PDL","PD1","IL10","VEGF")        
bidirectional_keywords <- c("CX3C","CXCL12")                        # 示例，可自定义
classify_pathway <- function(pathway){
  if(any(str_detect(pathway, promoting_keywords))) return("Promoting")
  if(any(str_detect(pathway, suppressive_keywords))) return("Suppressive")
  if(any(str_detect(pathway, bidirectional_keywords))) return("Bidirectional")
  return("Neutral")
}

df <- df %>% mutate(effect = sapply(pathway_name, classify_pathway))
pathway_stats <- df %>%
  group_by(effect, pathway_name) %>%
  summarise(mean_prob = mean(prob, na.rm = TRUE), .groups="drop") %>%
  group_by(effect) %>%
  mutate(prop = mean_prob / sum(mean_prob)) %>%
  ungroup()

write_csv(pathway_stats, "CAF_TBChat_pathway_stats.csv")
ggplot(pathway_stats, aes(x=fct_reorder(pathway_name, mean_prob), y=mean_prob, fill=effect)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_minimal() +
  labs(title="CAF → TBChat: 通路影响值", y="平均影响值", x="通路") +
  scale_fill_brewer(palette="Set2")
ggsave("CAF_TBChat_pathway_effects.png", width=10, height=6)
ecm_targets <- c("TGFB1","TGFB3","FN1","COL1A1","COL1A2","COL4A1","COL4A2","COL6A1",
                 "COL6A2","COL6A3","THBS1","THBS2","LAMA2","LAMA4","LAMA5",
                 "LAMB1","LAMB2","LAMC1","LAMC2")

df <- df %>% mutate(ECM_related = ifelse(ligand %in% ecm_targets | receptor %in% ecm_targets, TRUE, FALSE))

ecm_stats <- df %>%
  group_by(ECM_related) %>%
  summarise(count=n(),
            mean_prob=mean(prob, na.rm=TRUE),
            total_prob=sum(prob, na.rm=TRUE)) %>%
  mutate(prop_count=count/sum(count),
         prop_prob=total_prob/sum(total_prob))

write_csv(ecm_stats, "CAF_TBChat_ECM_stats.csv")
ggplot(ecm_stats, aes(x=factor(ECM_related), y=prop_prob, fill=ECM_related)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(title="CAF → TBChat ECM markers 占比", x="ECM_related", y="影响值比例") +
  scale_fill_brewer(palette="Set2")
ggsave("CAF_TBChat_ECM_prop.png", width=6, height=6)
promoting_df <- df %>% filter(effect=="Promoting")

promoting_ecm <- promoting_df %>%
  group_by(ECM_related) %>%
  summarise(count=n(),
            mean_prob=mean(prob, na.rm=TRUE),
            total_prob=sum(prob, na.rm=TRUE)) %>%
  mutate(prop_count=count/sum(count),
         prop_prob=total_prob/sum(total_prob))
write_csv(promoting_ecm, "CAF_TBChat_promoting_ECM_stats.csv")
ggplot(promoting_ecm, aes(x=factor(ECM_related), y=prop_prob, fill=ECM_related)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(title="CAF → TBChat 促进通路中 ECM markers 占比", x="ECM_related", y="影响值比例") +
  scale_fill_brewer(palette="Set2")
ggsave("CAF_TBChat_promoting_ECM_prop.png", width=6, height=6)



library(dplyr)
library(readr)
library(igraph)
library(ggraph)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(scales)
output_dir <- "猫猫的生信文件夹"
if (!dir.exists(output_dir)) dir.create(output_dir)
df_all <- read_csv("CellChat_LR_interactions_all.csv")
pathway_info <- read_csv("CAF_TB_Pathway_Impact_All.csv")

caf_subtypes <- c("myCAF","iCAF","apCAF")
tb_subtypes  <- c("CD4","CD8","Treg","exhaustion","Bcell","Plasma")

df <- df_all %>% filter(source %in% caf_subtypes & target %in% tb_subtypes)

ecm_genes <- c("TGFB1","TGFB3","FN1","COL1A1","COL1A2",
               "COL4A1","COL4A2","COL6A1","COL6A2","COL6A3",
               "THBS1","THBS2","LAMA2","LAMA4","LAMA5",
               "LAMB1","LAMB2","LAMC1","LAMC3","APP")
df2 <- left_join(df, pathway_info, by="pathway_name")

g <- graph_from_data_frame(df2, directed=TRUE)
p_net <- ggraph(g, layout="fr") +
  geom_edge_link(aes(width=prob, color=effect), alpha=0.8) +
  geom_node_point(size=6, color="steelblue") +
  geom_node_text(aes(label=name), repel=TRUE, size=4) +
  scale_edge_color_manual(values=c("Promoting"="#E69F00","Inhibiting"="#56B4E9")) +
  theme_void() +
  theme(legend.position="bottom",
        legend.title=element_text(size=12),
        legend.text=element_text(size=10))
promote_df <- df2 %>% filter(effect=="Promoting") %>%
  mutate(ECM_related = ifelse(ligand %in% ecm_genes | receptor %in% ecm_genes, "ECM", "Other"))

ecm_stats <- promote_df %>%
  group_by(ECM_related) %>%
  summarise(n = n(), mean_prob = mean(prob, na.rm = TRUE)) %>%
  mutate(ratio = n / sum(n),
         prob_ratio = mean_prob / sum(mean_prob),
         ratio_label = paste0(round(ratio*100,1), "%"))
p_ecm_pie <- ggplot(ecm_stats, aes(x="", y=ratio, fill=ECM_related)) +
  geom_col(color="white") +
  coord_polar(theta="y") +
  geom_text(aes(label=ratio_label), position=position_stack(vjust=0.5), size=5, fontface="bold") +
  scale_fill_brewer(palette="Dark2") +
  theme_void(base_size=14) +
  labs(title="Promoting通路中ECM数量占比") +
  theme(plot.title = element_text(face="bold", hjust=0.5),
        legend.position="none")
y_max <- max(ecm_stats$mean_prob) * 1.2

p_ecm_bar <- ggplot(ecm_stats, aes(x=ECM_related, y=mean_prob, fill=ECM_related)) +
  geom_col(width=0.6, color="black") +
  geom_text(aes(label=paste0(round(prob_ratio*100,1), "%")), vjust=-0.5, size=5, fontface="bold") +
  scale_fill_brewer(palette="Dark2") +
  scale_y_continuous(limits=c(0, y_max), labels = scales::percent_format(accuracy=1)) +
  theme_minimal(base_size=14) +
  labs(title="Promoting通路中ECM平均通信强度贡献", x=NULL, y="Mean probability (%)") +
  theme(legend.position="none",
        plot.title = element_text(face="bold", hjust=0.5))
combined_plot <- p_net / (p_ecm_pie | p_ecm_bar) + 
  plot_layout(heights = c(2, 1))
print(combined_plot)
ggsave(file.path(output_dir, "CAF_TB_ECM_Combined_Nature.png"), combined_plot, width=12, height=10, dpi=300)
ecm_stats <- promote_df %>%
  group_by(ECM_related) %>%
  summarise(n = n(),
            sum_prob = sum(prob, na.rm = TRUE)) %>%
  mutate(prob_ratio = sum_prob / sum(sum_prob),  # 真实贡献占比
         ratio_label = paste0(round(prob_ratio*100,1), "%"))
p_ecm_bar <- ggplot(ecm_stats, aes(x=ECM_related, y=prob_ratio, fill=ECM_related)) +
  geom_col(width=0.6, color="black") +
  geom_text(aes(label=ratio_label), vjust=-0.5, size=5, fontface="bold") +
  scale_fill_brewer(palette="Dark2") +
  scale_y_continuous(labels = scales::percent_format(accuracy=1), limits=c(0,1)) +
  theme_minimal(base_size=14) +
  labs(title="Promoting通路中ECM真实贡献占比", x=NULL, y="Contribution (%)") +
  theme(legend.position="none",
        plot.title = element_text(face="bold", hjust=0.5))
print(p_ecm_bar)
ggsave(file.path(output_dir, "CAF_TB_ECM_bar_RealContribution.png"), 
       p_ecm_bar, width=6, height=5, dpi=300)




library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
data <- read_csv("CAF_TB_interactions.csv")r
genes_of_interest <- c(
  "TGFB1","TGFB3","FN1","COL1A1","COL1A2",
  "COL4A1","COL4A2","COL6A1","COL6A2","COL6A3",
  "THBS1","THBS2","LAMA2","LAMA4","LAMA5",
  "LAMB1","LAMB2","LAMC1","LAMC3","APP"
)
plot_tb <- data %>%
  filter(ligand %in% genes_of_interest)
targets <- unique(plot_tb$target)
complete_tb <- expand.grid(
  ligand = genes_of_interest,
  target = targets,
  stringsAsFactors = FALSE
)
plot_tb_full <- complete_tb %>%
  left_join(plot_tb, by = c("ligand", "target"))
plot_tb_full$prob[is.na(plot_tb_full$prob)] <- 0

# 颜色方案：
#   0   → 极淡蓝 (无表达)
#   低值 → 天蓝 (#9ecae1)
#   中值 → 蓝 (#3182bd)
#   高值 → 藏蓝 (#08306b)
colors <- c("#f0f7ff", "#9ecae1", "#3182bd", "#08306b")
pdf("CAF_to_TB_collagen_maker_heatmap_enhanced_blue.pdf", width = 10, height = 6)

ggplot(plot_tb_full, aes(x = ligand, y = target, fill = prob)) +
  geom_tile(color = NA) +
  scale_fill_gradientn(
    colours = colors,
    limits = c(0, max(plot_tb_full$prob, na.rm = TRUE)),
    na.value = "#f7fbff"
  ) +
  theme_classic(base_size = 12) +
  labs(
    title = "CAF → TB (collagen/ECM maker interactions)",
    x = "CAF ligand/maker",
    y = "Target TB subtype",
    fill = "Sum prob"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )

dev.off()
pkgs <- c("data.table","dplyr","ggplot2","pheatmap","RColorBrewer","scales","igraph", "stringr")
for(p in pkgs) if(!requireNamespace(p, quietly=TRUE)) install.packages(p)
library(data.table); library(dplyr); library(ggplot2); library(pheatmap)
library(RColorBrewer); library(scales); library(igraph); library(stringr)
infile <- "CellChat_LR_interactions_all.csv"
outdir <- "CAF_TB_promote_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
if(!file.exists(infile)) stop("找不到文件: ", infile)
dt <- fread(infile)
message("Loaded interactions: ", nrow(dt), " rows")
dt <- dt %>% mutate(source = as.character(source), target = as.character(target))
dt <- dt %>% mutate(source = str_replace(source, "1$", ""),
                    target = str_replace(target, "1$", ""))
caf_subtypes <- c("myCAF","iCAF","apCAF")
tb_subtypes  <- c("CD8","CD4","exhaustion","Treg","Bcell","Plasma")
dt_caf_tb <- dt %>% filter(source %in% caf_subtypes & target %in% tb_subtypes)
message("CAF->TB interactions kept: ", nrow(dt_caf_tb))
fwrite(dt_caf_tb, file.path(outdir, "CAF_to_TB_all_interactions.csv"))
promote_keywords <- c("activ", "stimulat", "promot", "recruit", "chemo", "migra", "prolif", "growth",
                      "surviv", "costim", "co-stim", "TNF", "IL", "CXCL","CCL","MIF","EGF","HGF","VEGF",
                      "PDGF","IFN","ICOS","CD80","CD86","OX40","4-1BB","TNFSF","IL6","IL7","IL15","SPP1","APP")
promote_annotation_kw <- c("Secreted", "Cytokine", "Chemokine", "Growth factor", "TNF", "IL")
has_kw <- function(txt, kws){
  if(is.na(txt) || txt == "") return(FALSE)
  t <- tolower(as.character(txt))
  any(sapply(tolower(kws), function(k) grepl(k, t, fixed = FALSE)))
}
dt_caf_tb <- dt_caf_tb %>%
  mutate(
    is_promote = mapply(function(lig, rec, path, annot, iname){
      any(c(
        has_kw(lig, promote_keywords),
        has_kw(rec, promote_keywords),
        has_kw(path, promote_keywords),
        has_kw(annot, promote_annotation_kw),
        has_kw(iname, promote_keywords)
      ))
    }, ligand, receptor, pathway_name, annotation, interaction_name)
  )
dt_promote <- dt_caf_tb %>% filter(is_promote)
dt_other   <- dt_caf_tb %>% filter(!is_promote)

message("Promoting interactions (heuristic): ", nrow(dt_promote),
        " ; Non-promoting: ", nrow(dt_other))
fwrite(dt_promote, file.path(outdir, "CAF_to_TB_promoting_interactions.csv"))
fwrite(dt_other,   file.path(outdir, "CAF_to_TB_nonpromoting_interactions.csv"))
dt_promote <- dt_promote %>% mutate(prob = as.numeric(prob))
summary_tbl <- dt_promote %>%
  group_by(source, target) %>%
  summarise(n_interactions = n(),
            sum_prob = sum(ifelse(is.na(prob), 0, prob)),
            mean_prob = mean(ifelse(is.na(prob), 0, prob)),
            top_ligands = paste(head(unique(ligand[!is.na(ligand)]), 6), collapse="; "),
            top_pathways = paste(head(unique(pathway_name[!is.na(pathway_name)]),6), collapse="; "),
            .groups="drop") %>%
  arrange(desc(sum_prob))

fwrite(summary_tbl, file.path(outdir, "CAF_to_TB_promote_summary_by_pair.csv"))
plot_tb <- dt_promote %>%
  mutate(pathway = ifelse(is.na(pathway_name) | pathway_name=="", "Unknown", pathway_name)) %>%
  group_by(source, target, interaction_name, ligand, receptor, pathway) %>%
  summarise(prob = max(as.numeric(prob), na.rm=TRUE), .groups="drop")
plot_tb$prob[is.na(plot_tb$prob)] <- 0
nature_pal <- c(RColorBrewer::brewer.pal(8,"Dark2"), RColorBrewer::brewer.pal(8,"Set2"))
names(nature_pal) <- unique(plot_tb$pathway)[1:length(unique(plot_tb$pathway))]
p_bubble <- ggplot(plot_tb, aes(x = target, y = source, size = prob, color = pathway)) +
  geom_point(alpha=0.95) +
  scale_size_continuous(range=c(2,9), breaks = pretty(plot_tb$prob, n=4)) +
  scale_color_manual(values = rep(nature_pal, length.out=length(unique(plot_tb$pathway)))) +
  theme_classic(base_size = 12) +
  labs(title = "CAF -> T/B promoting interactions (bubble)",
       x = "Target (T/B subtypes)", y = "Source (CAF subtypes)",
       size = "prob (weight)", color = "pathway") +
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggsave(filename = file.path(outdir, "CAF_to_TB_promote_bubble.pdf"), plot = p_bubble, width=8, height=5)
ggsave(filename = file.path(outdir, "CAF_to_TB_promote_bubble.png"), plot = p_bubble, width=8, height=5, dpi=300)
mat <- summary_tbl %>% select(source, target, sum_prob) %>%
  tidyr::pivot_wider(names_from = target, values_from = sum_prob, values_fill = 0) %>%
  as.data.frame()
rownames(mat) <- mat$source; mat$source <- NULL
mat_m <- as.matrix(mat)

if(all(mat_m == 0)) {
  warning("heatmap matrix all zero; will plot but interpret carefully")
}
ph_colors <- colorRampPalette(c("white","#fde0dd","#fa9fb5","#c51b8a"))(50)
pdf(file.path(outdir, "CAF_to_TB_promote_heatmap.pdf"), width=6 + ncol(mat_m)*0.4, height=4 + nrow(mat_m)*0.4)
pheatmap(mat_m, cluster_rows = FALSE, cluster_cols = FALSE, color = ph_colors,
         main = "CAF -> T/B (sum prob of promoting interactions)",
         fontsize_row = 10, fontsize_col = 10)
dev.off()

png(file.path(outdir, "CAF_to_TB_promote_heatmap.png"), width=1000, height=600, res=150)
pheatmap(mat_m, cluster_rows = FALSE, cluster_cols = FALSE, color = ph_colors,
         main = "CAF -> T/B (sum prob of promoting interactions)",
         fontsize_row = 10, fontsize_col = 10)
dev.off()
edges <- summary_tbl %>% select(source, target, sum_prob) %>% filter(sum_prob > 0)
if(nrow(edges) == 0){
  warning("No weighted edges to plot in network")
} else {
  g <- graph_from_data_frame(edges, directed = TRUE)
  E(g)$width <- scales::rescale(edges$sum_prob, to = c(1,8))
  V(g)$color <- ifelse(V(g)$name %in% caf_subtypes, "#8dd3c7", "#fb8072") # CAF vs TB simple coloring
  pdf(file.path(outdir, "CAF_to_TB_promote_network.pdf"), width=8, height=6)
  plot(g, edge.arrow.size=0.6, vertex.size=30, vertex.label.cex=0.9,
       main="CAF -> T/B promoting network (edge width ~ sum_prob)")
  dev.off()
  png(file.path(outdir, "CAF_to_TB_promote_network.png"), width=1200, height=900, res=150)
  plot(g, edge.arrow.size=0.6, vertex.size=30, vertex.label.cex=0.9,
       main="CAF -> T/B promoting network (edge width ~ sum_prob)")
  dev.off()
}
top_lig_by_target <- dt_promote %>%
  group_by(target, ligand) %>%
  summarise(total_prob = sum(as.numeric(prob), na.rm=TRUE), count = n(), .groups="drop") %>%
  arrange(target, desc(total_prob)) %>%
  group_by(target) %>%
  slice_head(n=10)

fwrite(top_lig_by_target, file.path(outdir, "CAF_promote_top_ligands_per_target.csv"))
report <- list(
  generated = Sys.time(),
  input_rows = nrow(dt),
  caf_tb_rows = nrow(dt_caf_tb),
  promote_rows = nrow(dt_promote),
  summary_pairs = nrow(summary_tbl)
)
fwrite(as.data.table(report), file.path(outdir, "CAF_to_TB_promote_report_summary.csv"))
message("All outputs saved to: ", normalizePath(outdir))


library(tidyverse)
df <- read_csv("CellChat_LR_interactions_all.csv")
df$source <- gsub("1$", "", df$source)
df$target <- gsub("1$", "", df$target)

caf_subtypes <- c("myCAF","iCAF","apCAF")
tb_subtypes  <- c("CD4","CD8","Treg","exhaustion","Bcell","Plasma")
df <- df %>%
  filter(source %in% caf_subtypes & target %in% tb_subtypes)
all_pathways <- unique(df$pathway_name)
promoting_keywords <- c("MIF","CD40","CXCL","CCL","ICAM","TNF","IFN","IL2","IL7","CD28","CD86","CD48",
                        "APP","FN1","ANXA","SPP1","EGF","FGF","VEGF","IGF","ANGPT","LAMININ",
                        "COLLAGEN","VCAM","CD44","SEMA","GAS6","WNT","DLL","NOTCH")

suppressive_keywords <- c("TGFB","PDL","PDCD","CD47","LGALS","SIRP","HAVCR","LAG3","IL10",
                          "IDO","ADORA","NT5E","PTGER","SIGLEC","CEACAM","CD200")

bidirectional_keywords <- c("CX3C","CXCL","SEMA","NOTCH","EPHA","Ephrins","FGF","IGF")
classify_pathway <- function(name){
  up <- toupper(name)
  if(any(str_detect(up, promoting_keywords))) return("Promoting")
  if(any(str_detect(up, suppressive_keywords))) return("Suppressive")
  if(any(str_detect(up, bidirectional_keywords))) return("Bidirectional")
  return("Unclassified")
}
pathway_annot <- tibble(
  pathway_name = all_pathways,
  effect = sapply(all_pathways, classify_pathway)
)
pathway_annot <- pathway_annot %>%
  mutate(note = case_when(
    effect == "Promoting" ~ "促进免疫活化、T/B细胞增殖或迁移",
    effect == "Suppressive" ~ "抑制免疫反应或诱导Treg",
    effect == "Bidirectional" ~ "具有环境依赖的双向调节作用",
    TRUE ~ "功能未明确或非典型免疫信号"
  ))

write_csv(pathway_annot, "CAF_TB_Pathway_Annotation.csv")
df2 <- left_join(df, pathway_annot, by="pathway_name") %>%
  group_by(pathway_name, effect) %>%
  summarise(mean_prob = mean(prob, na.rm=TRUE),
            n_interactions = n(),
            .groups = "drop")

write_csv(df2, "CAF_TB_Pathway_Impact_All.csv")
ggplot(df2, aes(x=reorder(pathway_name, mean_prob), y=mean_prob, fill=effect)) +
  geom_col(color="black", width=0.7) +
  coord_flip() +
  scale_fill_manual(values=c("Promoting"="#e74c3c",
                             "Suppressive"="#3498db",
                             "Bidirectional"="#f1c40f",
                             "Unclassified"="grey70")) +
  theme_classic(base_size=14) +
  labs(x="Pathway", y="Mean Communication Probability",
       title="CAF→T/B 通信通路促进/抑制影响力比较",
       fill="Effect") +
  theme(plot.title = element_text(face="bold", size=16, hjust=0.5))

cat("CAF_TB_Pathway_Annotation.csv 和 CAF_TB_Pathway_Impact_All.csv\n")
pathway_summary <- df2 %>%
  group_by(effect) %>%
  summarise(
    n_pathways = n(),
    total_prob = sum(mean_prob, na.rm=TRUE),
    mean_prob = mean(mean_prob, na.rm=TRUE)
  ) %>%
  mutate(
    pathway_ratio = n_pathways / sum(n_pathways),
    prob_ratio = total_prob / sum(total_prob)
  )

print(pathway_summary)
write_csv(pathway_summary, "CAF_TB_Pathway_Effect_Summary.csv")

ggplot(pathway_summary, aes(x="", y=prob_ratio, fill=effect)) +
  geom_bar(width=1, stat="identity", color="white") +
  coord_polar(theta="y") +
  scale_fill_manual(values=c("Promoting"="#e74c3c",
                             "Suppressive"="#3498db",
                             "Bidirectional"="#f1c40f",
                             "Unclassified"="grey70")) +
  theme_void() +
  labs(title="CAF→T/B 通信通路影响占比") +
  geom_text(aes(label=scales::percent(prob_ratio, accuracy=0.1)),
            position=position_stack(vjust=0.5), size=5, color="white") +
  theme(plot.title = element_text(face="bold", size=16, hjust=0.5))

ggplot(df2 %>% filter(effect == "Promoting"),
       aes(x=reorder(pathway_name, mean_prob), y=mean_prob)) +
  geom_col(fill="#e74c3c", color="black", width=0.7) +
  coord_flip() +
  theme_classic(base_size=14) +
  labs(x="Promoting Pathway", y="Mean Communication Probability",
       title="CAF→T/B 促进性通路影响强度") +
  theme(plot.title = element_text(face="bold", size=16, hjust=0.5))

cat(" CAF_TB_Pathway_Effect_Summary.csv \n")








