############################################################################
# ANALYSIS OF NON-EC DATA
############################################################################

#load libraries
library(knitr)
library(dplyr)
library(Seurat)
library(SeuratData)
library(patchwork)
library(harmony)
library(stringr)
library(base)
library(Matrix)
library(ggplot2)
library(clusterProfiler)
library(AnnotationDbi)
library(dittoSeq)
library(enrichplot)
library(stats)
library(org.Mm.eg.db)
library(msigdbr)
library(fgsea)
library(gridExtra)



# Non-EC data
#Load and read .rds files
file3 <- choose.files()
Non_EC <- readRDS(file3)
Non_EC

# 2.Quality control and filtering
View(Non_EC@meta.data)

NonEC <- Non_EC

NonEC@meta.data$orig.ident

# Change orig.ident column name to Genotype
colnames(NonEC@meta.data)[which(colnames(NonEC@meta.data) == "orig.ident")] <- "Genotype"
unique(NonEC@meta.data$Genotype)

# change Dll4LOF to Dll4-iDEC in genotype column
library(forcats)
NonEC@meta.data$Genotype <- fct_recode(NonEC@meta.data$Genotype,
                                       "Dll4-iDEC" = "Dll4LOF")



# Visualize QC metrics as a violin plot
NonEC <- SetIdent(NonEC, value = NonEC@meta.data$Genotype)
pdf("QC_NonEC_VlnPlot.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
VlnPlot(NonEC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        cols = c('green', 'red'), ncol = 3)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("QC_NonEC_ScatterPlot.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
plot1 <- FeatureScatter(NonEC, feature1 = "nCount_RNA", 
                        cols = c('green', 'red'), feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(NonEC, feature1 = "nCount_RNA", 
                        cols = c('green', 'red'), feature2 = "percent.mt")
plot1 + plot2
dev.off()

# Filtering the data
# default of mitochondrial genes is usually 5% but the heart is a special tissue, where cardiomyocytes
#can have up to 30% in mitochondrial genes
NonEC_filtered <- subset(NonEC,
                         subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 15 & DF.classifications_0.25_0.09_441 == 'Singlet')



#-------------
# 3. Normalize merged dataset
NonEC_n <- NormalizeData(NonEC_filtered,
                         normalization.method = "LogNormalize",
                         scale.factor = 10000)
# store names of all genes
AllGenes <- rownames(NonEC_n)

#-------------
# 4. Identify features that are outliers on a 'mean variability plot'
NonEC_n <- FindVariableFeatures(NonEC_n, selection.method = "vst",
                                nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(NonEC_n), 10)

# plot variable features with and without labels
pdf("VariableFeatures_NonEC.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
plot1 <- VariableFeaturePlot(NonEC_n)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()


#-------------
# 5. Scale Data
NonEC_n <- ScaleData(NonEC_n, features = AllGenes)


#-------------
# 6. PCA analysis of dataset
NonEC_n <- RunPCA(NonEC_n)

pdf("PCA_NonEC.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(NonEC_n, reduction = "pca", cols = c('green', 'red')) + NoLegend()
DimHeatmap(NonEC_n, dims = 1:30, cells = 500, balanced = TRUE)
ElbowPlot(NonEC_n, ndims = 30, reduction = "pca")
dev.off()

# 7. Check Dll4 deletion and Notch ligands for possible compensation
pdf("Expression of Notch pathway components_NonECs.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

VlnPlot(NonEC_n, features = c("Dll4", "Dll1", "Dll3", "Jag1", "Jag2",
                              "Notch1", "Notch2", "Notch3", 'Notch4', "Hes1",
                              "Hey1", "Hey2"),
        cols = c('green', 'red'), slot = "data", log = TRUE)

DotPlot(NonEC_n, features = c("Dll4", "Dll1", "Dll3", "Jag1", "Jag2",
                              "Notch1", "Notch2", "Notch3", 'Notch4', "Hes1",
                              "Hey1", "Hey2")) + coord_flip() +
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06")
dev.off()


#-------------
# 8. Find neighbours and clusters
NonEC_n <- FindNeighbors(NonEC_n, dims = 1:25, reduction = "pca")
NonEC_n <- FindClusters(NonEC_n,
                        resolution = c(0.1, 0.2, 0.4, 0.6, 0.8))

#-------------
# 9. Run non-linear dimensional reduction
# prepare UMAP
Idents(NonEC_n) <- 'RNA_snn_res.0.6'
NonEC_n <- RunUMAP(NonEC_n, dims = 1:25, metric = "euclidean",
                   min.dist = 0.3)

DimPlot(NonEC_n, reduction = "umap")


#-------------
# 10. Check sample distribution in UMAP
NonEC_n_samples <- SetIdent(NonEC_n, value = NonEC_n@meta.data$Genotype)
pdf("UMAP_NonEC_Samples dist.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(NonEC_n_samples, reduction = "umap", cols = c("green", "red"))
dev.off()


#-------------
# 11. Find cluster markers to annotate clusters (report only the positive ones)



NonEC_n.markers <- FindAllMarkers(NonEC_n, only.pos = TRUE)
NonEC_n.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05)

NonEC_n.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

pdf("HeatmapCluster_NonEC_res0.6.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DoHeatmap(NonEC_n, features = top10$gene, size = 3)  + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(text = element_text(size = 5))
dev.off()



#cluster0 fibroblasts
#cluster1 = fibroblasts
#cluster2 = macrophage
# cluster3 = SMC
# cluster 4 =) b cells
# cluster5 
# cluster 6 pericytes
# cluster 7 t cells
# cluster 8 leukocyte
# cluster9 
#cluster10 monocyte/macrophage
# cluster11
#cluster 12= macrophage
#cluster 13 monocyte  NK cells
#cluster 14 glial cells
#cluster15 


#-------------
# 12. Attribute new cluster names based on top markers
new.cluster.ids_NonEC0.6 <- c("Fibroblast1", "Fibroblast2", "Macrophage1", "SMC",
                              "B-cells", "?1", "Pericytes", "T-cells", "Leukocytes",
                              "?2", "Monocyte/macrophage", "?3", "Macrophage2",
                              "Monocyte NK cells", "Glial cells", "?4")
names(new.cluster.ids_NonEC0.6) <- levels(NonEC_n)
NonEC_n <- RenameIdents(NonEC_n, new.cluster.ids_NonEC0.6)
DimPlot(NonEC_n, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

pdf("UMAPCluster_NonECs_res0.6.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(NonEC_n, reduction = "umap", label = FALSE, pt.size = 0.5)
dev.off()


pdf("ViolinPlot_Top nonEC genes_res0.6.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
VlnPlot(NonEC_n, features <- c('Smoc2', 'Col15a1', 'Fbln1', 'Gsn',
                                    'C1qa', 'Mrc1', 'Myh11', 'Tagln',
                                    'Cd79a', 'Ly6d', 'Rspo3', 'Ecrg4',
                                    'Kcnj8', 'Vtn', 'Cd3e', 'Ms4a4b',
                                    'Cd300e', 'Ear2', 'Duoxa1', 'Efhd1',
                                    'Cxcr2', 'Csf3r', 'Npnt', 'Ccn3',
                                    'Ccr2', 'Ly6c2','Flt3', 'Cd209a',
                                    'Kcna1', 'Sox10', 'Msln', 'Upk3b'),
        group.by = 'RNA_snn_res.0.6', stack = TRUE, flip = TRUE)
dev.off()


pdf("DotPlot_Top nonEC genes_res0.6.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DotPlot(NonEC_n, features = c('Smoc2', 'Col15a1', 'Fbln1', 'Gsn',
                              'C1qa', 'Mrc1', 'Myh11', 'Tagln',
                              'Cd79a', 'Ly6d', 'Rspo3', 'Ecrg4',
                              'Kcnj8', 'Vtn', 'Cd3e', 'Ms4a4b',
                              'Cd300e', 'Ear2', 'Duoxa1', 'Efhd1',
                              'Cxcr2', 'Csf3r', 'Npnt', 'Ccn3',
                              'Ccr2', 'Ly6c2','Flt3', 'Cd209a',
                              'Kcna1', 'Sox10', 'Msln', 'Upk3b')) + 
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()
dev.off()


#-------------------
#-------------------
#repeat with resolution 0.4
# 9. Run non-linear dimensional reduction
# prepare UMAP
Idents(NonEC_n) <- 'RNA_snn_res.0.4'
NonEC_n <- RunUMAP(NonEC_n, dims = 1:25, metric = "euclidean",
                   min.dist = 0.3)

DimPlot(NonEC_n, reduction = "umap")

#UMAP with personalized colors
pdf("UMAP_NonEC_res0.4.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

DimPlot(object = NonEC_n, reduction = 'umap',
        cols = c('0' = 'red', '1' = 'grey', '2' = 'blue', '3' = '#57b000',
                 '4' = 'magenta', '5' = '#e6cc00', '6' = '#009dff',
                 '7' = '#8b4513', '8' = '#e47200', '9' = 'lightgreen', 
                 '10' = 'black', '11' = 'purple', '12' = 'darkblue',
                 '13' = 'darkred', '14' = 'pink'))
dev.off()


#-------------
# 10. Check sample distribution in UMAP
NonEC_n_samples <- SetIdent(NonEC_n, value = NonEC_n@meta.data$orig.ident)
pdf("UMAP_NonEC_Samples dist.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(NonEC_n_samples, reduction = "umap", cols = c("green", "red"))
dev.off()


#-------------
# 11. Find cluster markers to annotate clusters (report only the positive ones)
NonEC_n.markers <- FindAllMarkers(NonEC_n, only.pos = TRUE)
NonEC_n.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05)

NonEC_n.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

pdf("HeatmapCluster_NonEC_res0.4.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DoHeatmap(NonEC_n, features = top10$gene, size = 3)  + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(text = element_text(size = 5))
dev.off()

#NOT UPDATED THE ANNOTATION BELOW

#cluster0 fibroblasts
#cluster1 = fibroblasts
#cluster2 = macrophage
# cluster3 = SMC
# cluster 4 =) b cells
# cluster5 
# cluster 6 pericytes
# cluster 7 t cells
# cluster 8 leukocyte
# cluster9 
#cluster10 monocyte/macrophage
# cluster11
#cluster 12= macrophage
#cluster 13 monocyte  NK cells
#cluster 14 glial cells
#cluster15 


#-------------
# 12. Attribute new cluster names based on top markers
new.cluster.ids_NonEC0.4 <- c("Fibroblast1", "Fibroblast2", "Macrophage1", "SMC",
                              "B-cells", "?1", "Pericytes", "T-cells", "Leukocytes",
                              "?2", "Monocyte/macrophage", "?3", "Macrophage2",
                              "Monocyte NK cells", "Glial cells")
names(new.cluster.ids_NonEC0.4) <- levels(NonEC_n)
NonEC_n <- RenameIdents(NonEC_n, new.cluster.ids_NonEC0.4)
DimPlot(NonEC_n, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

pdf("HeatmapCluster_AllGenes_res0.4.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(NonEC_n, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

#-------------------
#-------------------
#repeat with resolution 0.1
# 9. Run non-linear dimensional reduction
# prepare UMAP
Idents(NonEC_n) <- 'RNA_snn_res.0.1'
NonEC_n <- RunUMAP(NonEC_n, dims = 1:25, metric = "euclidean",
                   min.dist = 0.3)

DimPlot(NonEC_n, reduction = "umap")

#-------------
# 11. Find cluster markers to annotate clusters (report only the positive ones)
NonEC_n.markers <- FindAllMarkers(NonEC_n, only.pos = TRUE)
NonEC_n.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05)

NonEC_n.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

pdf("HeatmapCluster_NonEC_res0.1.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DoHeatmap(NonEC_n, features = top10$gene, size = 3)  + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(text = element_text(size = 5))
dev.off()

#cluster0 Fibroblasts
#cluster1 = SMC/pericytes
#cluster2 = Macrophage
# cluster3 = B Cells
# cluster 4 = Leukocyte
# cluster5 = T Cells
# cluster 6 = Monocyte/macrophage?
# cluster 7 = Myofibroblast/fibroblast?
# cluster 8 = Schwann cells
# cluster9 = Epicardial cells



#-------------
# 12. Attribute new cluster names based on top markers
new.cluster.ids_NonEC0.1 <- c("Fibroblasts", "SMC/pericytes", "Macrophage", "B Cells",
                              "Leukocyte", "T Cells",
                              "Monocyte/macrophage", "Myofibroblast/fibroblast",
                              "Schwann cells", "Epicardial cells")
names(new.cluster.ids_NonEC0.1) <- levels(NonEC_n)
NonEC_n <- RenameIdents(NonEC_n, new.cluster.ids_NonEC0.1)
DimPlot(NonEC_n, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

pdf("HeatmapCluster_AllGenes_res0.1.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(NonEC_n, reduction = "umap", label = FALSE, pt.size = 0.7)
dev.off()


pdf("ViolinPlot_Top nonEC genes_res0.1.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
VlnPlot(NonEC_n, features <- c('Dcn', 'Mgp', 'Rgs5', 'Myh11',
                               'C1qa', 'Mrc1', 'Cd79a', 'Ly6d',
                               'Cd300e', 'Ear2', 'Cd3e', 'Cd3d',
                               'Cxcr2', 'Csf3r', 'Npnt', 'Ccn3',
                               'Kcna1', 'Sox10', 'Msln', 'Upk3b'),
        group.by = 'RNA_snn_res.0.1', stack = TRUE, flip = TRUE)
dev.off()


pdf("DotPlot_Top nonEC genes_res0.1.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DotPlot(NonEC_n, features = c('Dcn', 'Mgp', 'Rgs5', 'Myh11',
                              'C1qa', 'Mrc1', 'Cd79a', 'Ly6d',
                              'Cd300e', 'Ear2', 'Cd3e', 'Cd3d',
                              'Cxcr2', 'Csf3r', 'Npnt', 'Ccn3',
                              'Kcna1', 'Sox10', 'Msln', 'Upk3b')) + 
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()
dev.off()





# create metadata column with cluster ids assigned to each cell
NonEC_n@meta.data[["SR_res0.1_eu"]] <- NonEC_n@active.ident

# create new metadata column with celltype and genotype to do DEG analysis
NonEC_n@meta.data[["Celltype_genotype"]] <- paste0(NonEC_n@meta.data$SR_res0.1_eu,'_',NonEC_n@meta.data$Genotype)

#######################

pdf("Barplot with non EC_cells_res0.1_euc.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

dittoBarPlot(NonEC_n, var = "SR_res0.1_eu", group.by = "Genotype")
dev.off()



table(NonEC_n@meta.data$Celltype_genotype)

Control_Numb <- c(186, 27, 601, 62, 215, 11, 9, 24, 233, 64)
Dll4_Numb <- c(95, 5, 1160, 162, 282, 84, 74, 23, 444, 103)
Cell_freq <- cbind(Control_Numb, Dll4_Numb)
rownames(Cell_freq) <- c('B-cells', 'Epicardial', 'Fibroblasts', 'Leukocyte',
                         'Macrophage', 'Monocyte/Macrophage', 'Myofibroblast/fibroblast',
                         'Schwann cell', 'SMC/Pericytes', 'T-cells')
colnames(Cell_freq) <- c('Control', 'Dll4-iDEC')

# Export the table to a PDF file
Cell_freq_g <- tableGrob(Cell_freq)
pdf("Cell prop in nonEC clusters table_1.pdf")
grid.arrange(Cell_freq_g)
dev.off()

#-------------
# make celltype and genotype the active idents
Idents(NonEC_n) <- "Celltype_genotype"

#------------------------------------------------------------------------------
#DEG ANALYSIS OF FIBROBLASTS
#------------------------------------------------------------------------------
Fibroblast_DEG <- NonEC_n.markers <- FindMarkers(NonEC_n, ident.1 = 'Fibroblasts_Dll4-iDEC',
                               ident.2 = 'Fibroblasts_Control')
head(Fibroblast_DEG)
Fibroblast_DEG_genes <- rownames(Fibroblast_DEG)

# store data of number total DEG, Up significant and down significant
FibTotal <- nrow(Fibroblast_DEG)
nrow(subset(Fibroblast_DEG, Fibroblast_DEG$p_val_adj < 0.05))
Fib_UpSig <- nrow(subset(Fibroblast_DEG, Fibroblast_DEG$p_val_adj < 0.05 & Fibroblast_DEG$avg_log2FC > 0))
Fib_DwSig <- nrow(subset(Fibroblast_DEG, Fibroblast_DEG$p_val_adj < 0.05 & Fibroblast_DEG$avg_log2FC < 0))

Cluster <- c()
Total_DEG <- c()
UpSig_DEG <- c()
DwSig_DEG <- c()

Cluster <- c(Cluster, 'Fibroblast')
Total_DEG <- c(Total_DEG, FibTotal)
UpSig_DEG <- c(UpSig_DEG, Fib_UpSig)
DwSig_DEG <- c(DwSig_DEG, Fib_DwSig)








# ORA - Over-representation analysis
pdf("ORA_Fibroblasts.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

Fibro_GO_resultsCC <- enrichGO(gene = Fibroblast_DEG_genes,
                              universe = AllGenes, OrgDb = org.Mm.eg.db,
                              keyType = "SYMBOL", ont = "CC",
                              pAdjustMethod = 'BH',
                              pvalueCutoff = 0.05, qvalueCutoff = 0.05)


par(mar=c(8,20,3,1))
barplot(Fibro_GO_resultsCC, showCategory = 30, cex.names = 0.3,
        title = 'Fibroblast cellular component')
dotplot(Fibro_GO_resultsCC, showCategory=30) +
  ggtitle("dotplot for Fibroblast ORA CC")


Fibro_GO_resultsBP <- enrichGO(gene = Fibroblast_DEG_genes,
                              universe = AllGenes, OrgDb = org.Mm.eg.db,
                              keyType = "SYMBOL", ont = "BP",
                              pAdjustMethod = 'BH',
                              pvalueCutoff = 0.05, qvalueCutoff = 0.05)

par(mar=c(8,20,3,1))
barplot(Fibro_GO_resultsBP, showCategory = 30, cex.names = 0.3,
        title = 'Fibroblast biological processes')
dotplot(Fibro_GO_resultsBP, showCategory=30) +
  ggtitle("dotplot for Fibroblast ORA BP")

Fibro_GO_resultsMF <- enrichGO(gene = Fibroblast_DEG_genes,
                               universe = AllGenes, OrgDb = org.Mm.eg.db,
                              keyType = "SYMBOL", ont = "MF",
                              pAdjustMethod = 'BH',
                              pvalueCutoff = 0.05, qvalueCutoff = 0.05)

par(mar=c(8,30,3,1))
barplot(Fibro_GO_resultsMF, showCategory = 30, cex.names = 0.3,
        title = 'Fibroblast molecular function')
dotplot(Fibro_GO_resultsMF, showCategory=30) +
  ggtitle("dotplot for Fibroblast ORA MF")


dev.off()

# GO Gene Set Enrichment Analysis

## feature 1: numeric vector with Log2 FC
geneList = Fibroblast_DEG[,2]

## feature 2: named vector with gene names converted to ENTEZID
# Convert gene symbols to Entrez Gene IDs
entrez_ids <- mapIds(org.Mm.eg.db, keys = Fibroblast_DEG_genes,
                     keytype = "SYMBOL", column = "ENTREZID")

# Print the results
names(geneList) = as.character(Fibroblast_DEG_genes)

## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)



# get table wiht names and genes foldchange
names <- rownames(Fibroblast_DEG)
logfc <- Fibroblast_DEG$avg_log2FC
names(logfc) <- names


#get gene set to find hallmarks GO

m_df<- msigdbr(species = "Mus musculus", category = 'H')


fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#find enriched Hallmarks for our geneset
fgseaRes_Fibro <- fgsea(fgsea_sets, stats=logfc, minSize  = 2 )


fgseaResTidy_Fibro <- fgseaRes_Fibro %>%  as_tibble() %>%  arrange(desc(NES))
fgseaResTidy_Fibro %>%   dplyr::select(-leadingEdge, -ES) %>%   arrange(padj) %>%   head()


# Select top 10 upregulated and top 10 downregulated
a <- head(subset(fgseaResTidy_Fibro, fgseaResTidy_Fibro$NES > 0), 10)
b <- tail(subset(fgseaResTidy_Fibro, fgseaResTidy_Fibro$NES < 0), 10)
fgseaResTidy_FibroNEW <- rbind(a,b)

pdf("GSEA_Fibroblasts_Hallmarks.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

ggplot(fgseaResTidy_FibroNEW, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES)) +
  coord_flip() +
  scale_fill_gradient2(low='red3', mid='yellow2', high='green3')+
  labs(x="Hallmark", y="Normalized Enrichment Score",
       title="Hallmark NES from GSEA_Fibroblasts") + 
  theme_minimal()

dev.off()


#------------------------------------------------------------------------------
#DEG ANALYSIS OF SMC_PERICYTES
#------------------------------------------------------------------------------
SMC_pericytes_DEG <- NonEC_n.markers <- FindMarkers(NonEC_n, ident.1 = 'SMC/pericytes_Dll4-iDEC',
                                                 ident.2 = 'SMC/pericytes_Control')

head(SMC_pericytes_DEG)
SMC_pericytes_DEG_genes <- rownames(SMC_pericytes_DEG)

# store data of number total DEG, Up significant and down significant
SMATotal <- nrow(SMC_pericytes_DEG)
nrow(subset(SMC_pericytes_DEG, SMC_pericytes_DEG$p_val_adj < 0.05))
SMA_UpSig <- nrow(subset(SMC_pericytes_DEG, SMC_pericytes_DEG$p_val_adj < 0.05 & SMC_pericytes_DEG$avg_log2FC > 0))
SMA_DwSig <- nrow(subset(SMC_pericytes_DEG, SMC_pericytes_DEG$p_val_adj < 0.05 & SMC_pericytes_DEG$avg_log2FC < 0))

Cluster <- c(Cluster, 'SMA/Pericytes')
Total_DEG <- c(Total_DEG, SMATotal)
UpSig_DEG <- c(UpSig_DEG, SMA_UpSig)
DwSig_DEG <- c(DwSig_DEG, SMA_DwSig)








# ORA - Over-representation analysis
pdf("ORA_Fibroblasts.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

SMC_pericytes_GO_resultsCC <- enrichGO(gene = SMC_pericytes_DEG_genes,
                               universe = AllGenes, OrgDb = org.Mm.eg.db,
                               keyType = "SYMBOL", ont = "CC",
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05)


par(mar=c(8,20,3,1))
barplot(SMC_pericytes_GO_resultsCC, showCategory = 30, cex.names = 0.3,
        title = 'Fibroblast cellular component')
dotplot(SMC_pericytes_GO_resultsCC, showCategory=30) +
  ggtitle("dotplot for Fibroblast ORA CC")


SMC_pericytes_GO_resultsBP <- enrichGO(gene = SMC_pericytes_DEG_genes,
                               universe = AllGenes, OrgDb = org.Mm.eg.db,
                               keyType = "SYMBOL", ont = "BP",
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05)

par(mar=c(8,20,3,1))
barplot(SMC_pericytes_GO_resultsBP, showCategory = 30, cex.names = 0.3,
        title = 'Fibroblast biological processes')
dotplot(SMC_pericytes_GO_resultsBP, showCategory=30) +
  ggtitle("dotplot for Fibroblast ORA BP")

SMC_pericytes_GO_resultsMF <- enrichGO(gene = SMC_pericytes_DEG_genes,
                               universe = AllGenes, OrgDb = org.Mm.eg.db,
                               keyType = "SYMBOL", ont = "MF",
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05)

par(mar=c(8,30,3,1))
barplot(SMC_pericytes_GO_resultsMF, showCategory = 30, cex.names = 0.3,
        title = 'Fibroblast molecular function')
dotplot(SMC_pericytes_GO_resultsMF, showCategory=30) +
  ggtitle("dotplot for Fibroblast ORA MF")


dev.off()

# GO Gene Set Enrichment Analysis

## feature 1: numeric vector with Log2 FC
geneList = SMC_pericytes_DEG[,2]

## feature 2: named vector with gene names converted to ENTEZID
# Convert gene symbols to Entrez Gene IDs
entrez_ids <- mapIds(org.Mm.eg.db, keys = SMC_pericytes_DEG_genes,
                     keytype = "SYMBOL", column = "ENTREZID")

# Print the results
names(geneList) = as.character(SMC_pericytes_DEG_genes)

## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)



# get table wiht names and genes foldchange
names <- rownames(SMC_pericytes_DEG)
logfc <- SMC_pericytes_DEG$avg_log2FC
names(logfc) <- names


#get gene set to find hallmarks GO

m_df<- msigdbr(species = "Mus musculus", category = 'H')


fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#find enriched Hallmarks for our geneset
fgseaRes_SMC_pericytes <- fgsea(fgsea_sets, stats=logfc, minSize  = 2 )


fgseaResTidy_SMC_pericytes <- fgseaRes_SMC_pericytes %>%  as_tibble() %>%  arrange(desc(NES))
fgseaResTidy_SMC_pericytes %>%   dplyr::select(-leadingEdge, -ES) %>%   arrange(padj) %>%   head()


# Select top 10 upregulated and top 10 downregulated
a <- head(subset(fgseaResTidy_SMC_pericytes, fgseaResTidy_SMC_pericytes$NES > 0), 10)
b <- tail(subset(fgseaResTidy_SMC_pericytes, fgseaResTidy_SMC_pericytes$NES < 0), 10)
fgseaResTidy_SMC_pericytesNEW <- rbind(a,b)

pdf("GSEA_SMC_pericytes_Hallmarks.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

ggplot(fgseaResTidy_SMC_pericytesNEW, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES)) +
  coord_flip() +
  scale_fill_gradient2(low='red3', mid='yellow2', high='green3')+
  labs(x="Hallmark", y="Normalized Enrichment Score",
       title="Hallmark NES from GSEA_SMC_pericytes") + 
  theme_minimal()

dev.off()




#------------------------------------------------------------------------------
#DEG ANALYSIS OF MACROPHAGES
#------------------------------------------------------------------------------
Macrophage_DEG <- NonEC_n.markers <- FindMarkers(NonEC_n, ident.1 = 'Macrophage_Dll4-iDEC',
                                                    ident.2 = 'Macrophage_Control')

head(Macrophage_DEG)
Macrophage_DEG_genes <- rownames(Macrophage_DEG)

# store data of number total DEG, Up significant and down significant
MacTotal <- nrow(Macrophage_DEG)
nrow(subset(Macrophage_DEG, Macrophage_DEG$p_val_adj < 0.05))
Mac_UpSig <- nrow(subset(Macrophage_DEG, Macrophage_DEG$p_val_adj < 0.05 & Macrophage_DEG$avg_log2FC > 0))
Mac_DwSig <- nrow(subset(Macrophage_DEG, Macrophage_DEG$p_val_adj < 0.05 & Macrophage_DEG$avg_log2FC < 0))

Cluster <- c(Cluster, 'Macrophage')
Total_DEG <- c(Total_DEG, MacTotal)
UpSig_DEG <- c(UpSig_DEG, Mac_UpSig)
DwSig_DEG <- c(DwSig_DEG, Mac_DwSig)

# ORA - Over-representation analysis
pdf("ORA_Macrophage.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

Macrophage_GO_resultsCC <- enrichGO(gene = Macrophage_DEG_genes,
                               universe = AllGenes, OrgDb = org.Mm.eg.db,
                               keyType = "SYMBOL", ont = "CC",
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05)


par(mar=c(8,20,3,1))
barplot(Macrophage_GO_resultsCC, showCategory = 30, cex.names = 0.3,
        title = 'Macrophage cellular component')
dotplot(Macrophage_GO_resultsCC, showCategory=30) +
  ggtitle("dotplot for Macrophage ORA CC")


Macrophage_GO_resultsBP <- enrichGO(gene = Macrophage_DEG_genes,
                               universe = AllGenes, OrgDb = org.Mm.eg.db,
                               keyType = "SYMBOL", ont = "BP",
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05)

par(mar=c(8,20,3,1))
barplot(Macrophage_GO_resultsBP, showCategory = 30, cex.names = 0.3,
        title = 'Macrophage biological processes')
dotplot(Macrophage_GO_resultsBP, showCategory=30) +
  ggtitle("dotplot for Macrophage ORA BP")

Macrophage_GO_resultsMF <- enrichGO(gene = Macrophage_DEG_genes,
                               universe = AllGenes, OrgDb = org.Mm.eg.db,
                               keyType = "SYMBOL", ont = "MF",
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05)

par(mar=c(8,30,3,1))
barplot(Macrophage_GO_resultsMF, showCategory = 30, cex.names = 0.3,
        title = 'Macrophage molecular function')
dotplot(Macrophage_GO_resultsMF, showCategory=30) +
  ggtitle("dotplot for Macrophage ORA MF")


dev.off()

# GO Gene Set Enrichment Analysis

## feature 1: numeric vector with Log2 FC
geneList = Macrophage_DEG[,2]

## feature 2: named vector with gene names converted to ENTEZID
# Convert gene symbols to Entrez Gene IDs
entrez_ids <- mapIds(org.Mm.eg.db, keys = Macrophage_DEG_genes,
                     keytype = "SYMBOL", column = "ENTREZID")

# Print the results
names(geneList) = as.character(Macrophage_DEG_genes)

## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)



# get table wiht names and genes foldchange
names <- rownames(Macrophage_DEG)
logfc <- Macrophage_DEG$avg_log2FC
names(logfc) <- names


#get gene set to find hallmarks GO

m_df<- msigdbr(species = "Mus musculus", category = 'H')


fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#find enriched Hallmarks for our geneset
fgseaRes_Macrophage <- fgsea(fgsea_sets, stats=logfc, minSize  = 2 )


fgseaResTidy_Macrophage <- fgseaRes_Macrophage %>%  as_tibble() %>%  arrange(desc(NES))
fgseaResTidy_Macrophage %>%   dplyr::select(-leadingEdge, -ES) %>%   arrange(padj) %>%   head()


# Select top 10 upregulated and top 10 downregulated
a <- head(subset(fgseaResTidy_Macrophage, fgseaResTidy_Macrophage$NES > 0), 10)
b <- tail(subset(fgseaResTidy_Macrophage, fgseaResTidy_Macrophage$NES < 0), 10)
fgseaResTidy_MacrophageNEW <- rbind(a,b)

pdf("GSEA_Macrophage_Hallmarks.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

ggplot(fgseaResTidy_MacrophageNEW, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES)) +
  coord_flip() +
  scale_fill_gradient2(low='red3', mid='yellow2', high='green3')+
  labs(x="Hallmark", y="Normalized Enrichment Score",
       title="Hallmark NES from GSEA_Macrophage") + 
  theme_minimal()

dev.off()



#------------------------------------------------------------------------------
#DEG ANALYSIS OF LEUKOCYTES
#------------------------------------------------------------------------------
Leukocyte_DEG <- NonEC_n.markers <- FindMarkers(NonEC_n, ident.1 = 'Leukocyte_Dll4-iDEC',
                                                 ident.2 = 'Leukocyte_Control')
head(Leukocyte_DEG)
Leukocyte_DEG_genes <- rownames(Leukocyte_DEG)

# store data of number total DEG, Up significant and down significant
LeukTotal <- nrow(Leukocyte_DEG)
nrow(subset(Leukocyte_DEG, Leukocyte_DEG$p_val_adj < 0.05))
Leuk_UpSig <- nrow(subset(Leukocyte_DEG, Leukocyte_DEG$p_val_adj < 0.05 & Leukocyte_DEG$avg_log2FC > 0))
Leuk_DwSig <- nrow(subset(Leukocyte_DEG, Leukocyte_DEG$p_val_adj < 0.05 & Leukocyte_DEG$avg_log2FC < 0))

Cluster <- c(Cluster, 'Leukocytes')
Total_DEG <- c(Total_DEG, LeukTotal)
UpSig_DEG <- c(UpSig_DEG, Leuk_UpSig)
DwSig_DEG <- c(DwSig_DEG, Leuk_DwSig)

# ORA - Over-representation analysis
pdf("ORA_Leukocyte.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

Leukocyte_GO_resultsCC <- enrichGO(gene = Leukocyte_DEG_genes,
                               universe = AllGenes, OrgDb = org.Mm.eg.db,
                               keyType = "SYMBOL", ont = "CC",
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05)


par(mar=c(8,20,3,1))
barplot(Leukocyte_GO_resultsCC, showCategory = 30, cex.names = 0.3,
        title = 'Leukocyte cellular component')
dotplot(Leukocyte_GO_resultsCC, showCategory=30) +
  ggtitle("dotplot for Leukocyte ORA CC")


Leukocyte_GO_resultsBP <- enrichGO(gene = Leukocyte_DEG_genes,
                               universe = AllGenes, OrgDb = org.Mm.eg.db,
                               keyType = "SYMBOL", ont = "BP",
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05)

par(mar=c(8,20,3,1))
barplot(Leukocyte_GO_resultsBP, showCategory = 30, cex.names = 0.3,
        title = 'Leukocyte biological processes')
dotplot(Leukocyte_GO_resultsBP, showCategory=30) +
  ggtitle("dotplot for Leukocyte ORA BP")

Leukocyte_GO_resultsMF <- enrichGO(gene = Leukocyte_DEG_genes,
                               universe = AllGenes, OrgDb = org.Mm.eg.db,
                               keyType = "SYMBOL", ont = "MF",
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05)

par(mar=c(8,30,3,1))
barplot(Leukocyte_GO_resultsMF, showCategory = 30, cex.names = 0.3,
        title = 'Leukocyte molecular function')
dotplot(Leukocyte_GO_resultsMF, showCategory=30) +
  ggtitle("dotplot for Leukocyte ORA MF")


dev.off()

# GO Gene Set Enrichment Analysis

## feature 1: numeric vector with Log2 FC
geneList = Leukocyte_DEG[,2]

## feature 2: named vector with gene names converted to ENTEZID
# Convert gene symbols to Entrez Gene IDs
entrez_ids <- mapIds(org.Mm.eg.db, keys = Leukocyte_DEG_genes,
                     keytype = "SYMBOL", column = "ENTREZID")

# Print the results
names(geneList) = as.character(Leukocyte_DEG_genes)

## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)



# get table wiht names and genes foldchange
names <- rownames(Leukocyte_DEG)
logfc <- Leukocyte_DEG$avg_log2FC
names(logfc) <- names


#get gene set to find hallmarks GO

m_df<- msigdbr(species = "Mus musculus", category = 'H')


fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#find enriched Hallmarks for our geneset
fgseaRes_Leukocyte <- fgsea(fgsea_sets, stats=logfc, minSize  = 2 )


fgseaResTidy_Leukocyte <- fgseaRes_Leukocyte %>%  as_tibble() %>%  arrange(desc(NES))
fgseaResTidy_Leukocyte %>%   dplyr::select(-leadingEdge, -ES) %>%   arrange(padj) %>%   head()


# Select top 10 upregulated and top 10 downregulated
a <- head(subset(fgseaResTidy_Leukocyte, fgseaResTidy_Leukocyte$NES > 0), 10)
b <- tail(subset(fgseaResTidy_Leukocyte, fgseaResTidy_Leukocyte$NES < 0), 10)
fgseaResTidy_LeukocyteNEW <- rbind(a,b)

pdf("GSEA_Leukocyte_Hallmarks.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

ggplot(fgseaResTidy_LeukocyteNEW, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES)) +
  coord_flip() +
  scale_fill_gradient2(low='red3', mid='yellow2', high='green3')+
  labs(x="Hallmark", y="Normalized Enrichment Score",
       title="Hallmark NES from GSEA_Leukocyte") + 
  theme_minimal()

dev.off()


#------------------------------------------------------------------------------
#DEG ANALYSIS OF Monocyte_macrophage
#------------------------------------------------------------------------------
Monocyte_macrophage_DEG <- NonEC_n.markers <- FindMarkers(NonEC_n, ident.1 = 'Monocyte/macrophage_Dll4-iDEC',
                                                ident.2 = 'Monocyte/macrophage_Control')

head(Monocyte_macrophage_DEG)
Monocyte_macrophage_DEG_genes <- rownames(Monocyte_macrophage_DEG)

# store data of number total DEG, Up significant and down significant
MMTotal <- nrow(Monocyte_macrophage_DEG)
nrow(subset(Monocyte_macrophage_DEG, Monocyte_macrophage_DEG$p_val_adj < 0.05))
MM_UpSig <- nrow(subset(Monocyte_macrophage_DEG, Monocyte_macrophage_DEG$p_val_adj < 0.05 & Monocyte_macrophage_DEG$avg_log2FC > 0))
MM_DwSig <- nrow(subset(Monocyte_macrophage_DEG, Monocyte_macrophage_DEG$p_val_adj < 0.05 & Monocyte_macrophage_DEG$avg_log2FC < 0))

Cluster <- c(Cluster, 'Monocyte/Macrophage')
Total_DEG <- c(Total_DEG, MMTotal)
UpSig_DEG <- c(UpSig_DEG, MM_UpSig)
DwSig_DEG <- c(DwSig_DEG, MM_DwSig)



#------------------------------------------------------------------------------
#DEG ANALYSIS OF SCHWANN
#------------------------------------------------------------------------------
Schwann_DEG <- NonEC_n.markers <- FindMarkers(NonEC_n, ident.1 = 'Schwann cells_Dll4-iDEC',
                                                          ident.2 = 'Schwann cells_Control')

head(Schwann_DEG)

# store data of number total DEG, Up significant and down significant
SchwannTotal <- nrow(Schwann_DEG)
nrow(subset(Schwann_DEG, Schwann_DEG$p_val_adj < 0.05))
Schwann_UpSig <- nrow(subset(Schwann_DEG, Schwann_DEG$p_val_adj < 0.05 & Schwann_DEG$avg_log2FC > 0))
Schwann_DwSig <- nrow(subset(Schwann_DEG, Schwann_DEG$p_val_adj < 0.05 & Schwann_DEG$avg_log2FC < 0))

Cluster <- c(Cluster, 'Schwann')
Total_DEG <- c(Total_DEG, SchwannTotal)
UpSig_DEG <- c(UpSig_DEG, Schwann_UpSig)
DwSig_DEG <- c(DwSig_DEG, Schwann_DwSig)


#------------------------------------------------------------------------------
#DEG ANALYSIS OF B-CELLS
#------------------------------------------------------------------------------
Bcells_DEG <- NonEC_n.markers <- FindMarkers(NonEC_n, ident.1 = 'B Cells_Dll4-iDEC',
                                              ident.2 = 'B Cells_Control')

head(Bcells_DEG)
Bcells_DEG_genes <- rownames(Bcells_DEG)

# store data of number total DEG, Up significant and down significant
BTotal <- nrow(Bcells_DEG)
nrow(subset(Bcells_DEG, Bcells_DEG$p_val_adj < 0.05))
B_UpSig <- nrow(subset(Bcells_DEG, Bcells_DEG$p_val_adj < 0.05 & Bcells_DEG$avg_log2FC > 0))
B_DwSig <- nrow(subset(Bcells_DEG, Bcells_DEG$p_val_adj < 0.05 & Bcells_DEG$avg_log2FC < 0))

Cluster <- c(Cluster, 'B cells')
Total_DEG <- c(Total_DEG, BTotal)
UpSig_DEG <- c(UpSig_DEG, B_UpSig)
DwSig_DEG <- c(DwSig_DEG, B_DwSig)

# ORA - Over-representation analysis
pdf("ORA_Bcells.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

Bcells_GO_resultsCC <- enrichGO(gene = Bcells_DEG_genes,
                               universe = AllGenes, OrgDb = org.Mm.eg.db,
                               keyType = "SYMBOL", ont = "CC",
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05)


par(mar=c(8,20,3,1))
barplot(Bcells_GO_resultsCC, showCategory = 30, cex.names = 0.3,
        title = 'Bcells cellular component')
dotplot(Bcells_GO_resultsCC, showCategory=30) +
  ggtitle("dotplot for Bcells ORA CC")


Bcells_GO_resultsBP <- enrichGO(gene = Bcells_DEG_genes,
                               universe = AllGenes, OrgDb = org.Mm.eg.db,
                               keyType = "SYMBOL", ont = "BP",
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05)

par(mar=c(8,20,3,1))
barplot(Bcells_GO_resultsBP, showCategory = 30, cex.names = 0.3,
        title = 'Bcells biological processes')
dotplot(Bcells_GO_resultsBP, showCategory=30) +
  ggtitle("dotplot for Bcells ORA BP")

Bcells_GO_resultsMF <- enrichGO(gene = Bcells_DEG_genes,
                               universe = AllGenes, OrgDb = org.Mm.eg.db,
                               keyType = "SYMBOL", ont = "MF",
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05)

par(mar=c(8,30,3,1))
barplot(Bcells_GO_resultsMF, showCategory = 30, cex.names = 0.3,
        title = 'Bcells molecular function')
dotplot(Bcells_GO_resultsMF, showCategory=30) +
  ggtitle("dotplot for Bcells ORA MF")


dev.off()

# GO Gene Set Enrichment Analysis

## feature 1: numeric vector with Log2 FC
geneList = Bcells_DEG[,2]

## feature 2: named vector with gene names converted to ENTEZID
# Convert gene symbols to Entrez Gene IDs
entrez_ids <- mapIds(org.Mm.eg.db, keys = Bcells_DEG_genes,
                     keytype = "SYMBOL", column = "ENTREZID")

# Print the results
names(geneList) = as.character(Bcells_DEG_genes)

## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)



# get table wiht names and genes foldchange
names <- rownames(Bcells_DEG)
logfc <- Bcells_DEG$avg_log2FC
names(logfc) <- names


#get gene set to find hallmarks GO

m_df<- msigdbr(species = "Mus musculus", category = 'H')


fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#find enriched Hallmarks for our geneset
fgseaRes_Bcells <- fgsea(fgsea_sets, stats=logfc, minSize  = 2 )


fgseaResTidy_Bcells <- fgseaRes_Bcells %>%  as_tibble() %>%  arrange(desc(NES))
fgseaResTidy_Bcells %>%   dplyr::select(-leadingEdge, -ES) %>%   arrange(padj) %>%   head()


# Select top 10 upregulated and top 10 downregulated
a <- head(subset(fgseaResTidy_Bcells, fgseaResTidy_Bcells$NES > 0), 10)
b <- tail(subset(fgseaResTidy_Bcells, fgseaResTidy_Bcells$NES < 0), 10)
fgseaResTidy_BcellsNEW <- rbind(a,b)

pdf("GSEA_Bcells_Hallmarks.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

ggplot(fgseaResTidy_BcellsNEW, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES)) +
  coord_flip() +
  scale_fill_gradient2(low='red3', mid='yellow2', high='green3')+
  labs(x="Hallmark", y="Normalized Enrichment Score",
       title="Hallmark NES from GSEA_Bcells") + 
  theme_minimal()

dev.off()



#------------------------------------------------------------------------------
#DEG ANALYSIS OF T CELLS
#------------------------------------------------------------------------------
Tcells_DEG <- NonEC_n.markers <- FindMarkers(NonEC_n, ident.1 = 'T Cells_Dll4-iDEC',
                                             ident.2 = 'T Cells_Control')

head(Tcells_DEG)

# store data of number total DEG, Up significant and down significant
TTotal <- nrow(Tcells_DEG)
nrow(subset(Tcells_DEG, Tcells_DEG$p_val_adj < 0.05))
T_UpSig <- nrow(subset(Tcells_DEG, Tcells_DEG$p_val_adj < 0.05 & Tcells_DEG$avg_log2FC > 0))
T_DwSig <- nrow(subset(Tcells_DEG, Tcells_DEG$p_val_adj < 0.05 & Tcells_DEG$avg_log2FC < 0))

Cluster <- c(Cluster, 'T cells')
Total_DEG <- c(Total_DEG, TTotal)
UpSig_DEG <- c(UpSig_DEG, T_UpSig)
DwSig_DEG <- c(DwSig_DEG, T_DwSig)



#------------------------------------------------------------------------------
#DEG ANALYSIS OF epicardial
#------------------------------------------------------------------------------
Epicardial_DEG <- NonEC_n.markers <- FindMarkers(NonEC_n, ident.1 = 'Epicardial cells_Dll4-iDEC',
                                             ident.2 = 'Epicardial cells_Control')

head(Epicardial_DEG)

# store data of number total DEG, Up significant and down significant
EpiTotal <- nrow(Epicardial_DEG)
nrow(subset(Epicardial_DEG, Epicardial_DEG$p_val_adj < 0.05))
Epi_UpSig <- nrow(subset(Epicardial_DEG, Epicardial_DEG$p_val_adj < 0.05 & Epicardial_DEG$avg_log2FC > 0))
Epi_DwSig <- nrow(subset(Epicardial_DEG, Epicardial_DEG$p_val_adj < 0.05 & Epicardial_DEG$avg_log2FC < 0))

Cluster <- c(Cluster, 'Epicardial cells')
Total_DEG <- c(Total_DEG, EpiTotal)
UpSig_DEG <- c(UpSig_DEG, Epi_UpSig)
DwSig_DEG <- c(DwSig_DEG, Epi_DwSig)

#------------------------------------------------------------------------------
#DEG ANALYSIS OF Myofibroblast_fibroblas
#------------------------------------------------------------------------------
Myofibroblast_fibroblast_DEG <- NonEC_n.markers <- FindMarkers(NonEC_n, ident.1 = 'Myofibroblast/fibroblast_Dll4-iDEC',
                                                 ident.2 = 'Myofibroblast/fibroblast_Control')

head(Myofibroblast_fibroblast_DEG)


# store data of number total DEG, Up significant and down significant
MFTotal <- nrow(Myofibroblast_fibroblast_DEG)
nrow(subset(Myofibroblast_fibroblast_DEG, Myofibroblast_fibroblast_DEG$p_val_adj < 0.05))
MF_UpSig <- nrow(subset(Myofibroblast_fibroblast_DEG, Myofibroblast_fibroblast_DEG$p_val_adj < 0.05 & Myofibroblast_fibroblast_DEG$avg_log2FC > 0))
MF_DwSig <- nrow(subset(Myofibroblast_fibroblast_DEG, Myofibroblast_fibroblast_DEG$p_val_adj < 0.05 & Myofibroblast_fibroblast_DEG$avg_log2FC < 0))

Cluster <- c(Cluster, 'Myofibroblast/Fibroblast')
Total_DEG <- c(Total_DEG, MFTotal)
UpSig_DEG <- c(UpSig_DEG, MF_UpSig)
DwSig_DEG <- c(DwSig_DEG, MF_DwSig)


# MAKE DF WITH deg FREQUENCIES

DEG_Freq <- cbind(Total_DEG, UpSig_DEG, DwSig_DEG)
rownames(DEG_Freq) <- Cluster
colnames(DEG_Freq) <- c('Total', 'Up (p_adj<0.05)', 'Down (p_adj<0.05)')
# Order the data frame based on row names
DEG_Freq <- DEG_Freq[order(rownames(DEG_Freq)), ]

# Export the table to a PDF file
DEG_Freq_g <- tableGrob(DEG_Freq)
pdf("DEG frequency.pdf")
grid.arrange(DEG_Freq_g)
dev.off()



saveRDS(NonEC_n, file='C:\\Users\\scferreira\\Desktop\\Non-EC analysis\\Non_ECs_nichenet.rds')

