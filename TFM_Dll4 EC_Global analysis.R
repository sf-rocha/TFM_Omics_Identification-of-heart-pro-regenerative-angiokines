############################################################################
## ANALYSIS OF CONTROL AND DLL4-IDEC IN ALL HEART ECS
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

#-----------------------------------------------------------------------------
# 1.Load and read .rds files
#-----------------------------------------------------------------------------
file1 <- choose.files()
ContCUB <- readRDS(file1) # Control seurat object
ContCUB

file2 <- choose.files()
Dll4CUB <- readRDS(file2) # Dll4-iDEC seurat object
Dll4CUB

# 1. Merge the two datasets
Cont_Dll4_merge <- merge(ContCUB, y = Dll4CUB,
                    add.cell.ids = c("Control", "Dll4_iDEC"))
View(Cont_Dll4_merge@meta.data)



#-----------------------------------------------------------------------------
# 2.Quality control and filtering
#-----------------------------------------------------------------------------
# Visualize QC metrics as a violin plot
pdf("QC_All cells_VlnPlot.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
VlnPlot(Cont_Dll4_merge, cols = c("green", "red"),
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


# FeatureScatter is typically used to visualize feature-feature relationships,
# but can be used for anything calculated by the object,
# i.e. columns in object metadata, PC scores etc.

pdf("QC_All cells_ScatterPlot.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
plot1 <- FeatureScatter(Cont_Dll4_merge, cols = c("green", "red"),
                        feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(Cont_Dll4_merge, cols = c("green", "red"),
                        feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1 + plot2
dev.off()

# Filtering the data
# default of mitochondrial genes is usually 5% but the heart is a special 
# tissue, where cardiomyocytes can have up to 30% in mitochondrial genes
Cont_Dll4_filtered <- subset(Cont_Dll4_merge,
                             subset = nFeature_RNA > 200 & nFeature_RNA < 4500 &
                               percent.mt < 15 & 
                               DF.classifications_0.25_0.09_408 == 'Singlet')



#-----------------------------------------------------------------------------
# 3. Normalize merged dataset
#-----------------------------------------------------------------------------
Cont_Dll4_new <- NormalizeData(Cont_Dll4_filtered,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000)
# store names of all genes
AllGenes <- rownames(Cont_Dll4_new)



#-----------------------------------------------------------------------------
# 4. Identify features that are outliers on a 'mean variability plot'
#-----------------------------------------------------------------------------
Cont_Dll4_new <- FindVariableFeatures(Cont_Dll4_new, selection.method = "vst",
                                      nfeatures = 2000)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Cont_Dll4_new), 10)


# plot variable features all cells with and without labels
pdf("VariableFeatures_All cells.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
plot1 <- VariableFeaturePlot(Cont_Dll4_new)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,
                     xnudge = 0, ynudge = 0)
plot1 + plot2
dev.off()


#-----------------------------------------------------------------------------
# 5. Scale Data
#-----------------------------------------------------------------------------
Cont_Dll4_new <- ScaleData(Cont_Dll4_new, features = AllGenes)




#-----------------------------------------------------------------------------
# 6. PCA analysis of dataset
#-----------------------------------------------------------------------------
Cont_Dll4_new <- RunPCA(Cont_Dll4_new)


# generate plots/heatmaps to decide on PC for dimension reduction
pdf("PCA_All cells.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(Cont_Dll4_new, reduction = "pca", cols = c('green', 'red')) + NoLegend() 
DimHeatmap(Cont_Dll4_new, dims = 1:25, cells = 500, balanced = TRUE)
ElbowPlot(Cont_Dll4_new, ndims = 30, reduction = "pca")
dev.off()


#-----------------------------------------------------------------------------
# 7. Before continuing, check/validate that Dll4 is effectively deleted
# check also downstream targets and other Notch ligands/receptors for
# possible compensation
#-----------------------------------------------------------------------------
Cont_Dll4_new <- SetIdent(Cont_Dll4_new,
                          value = Cont_Dll4_new@meta.data$Genotype)

pdf("Expression of Notch pathway components_All cells.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

VlnPlot(Cont_Dll4_new, features = c("Dll4", "Dll1", "Dll3", "Jag1", "Jag2",
                                    "Notch1", "Notch2", "Notch3", "Notch4",
                                    "Hes1", "Hey1", "Hey2"),slot = "data",
        log = TRUE,  cols = c('green', 'red'))

DotPlot(Cont_Dll4_new, features = c("Dll4", "Dll1", "Dll3", "Jag1", "Jag2",
                                   "Notch1", "Notch2", "Notch3", "Notch4",
                                   "Hes1", "Hey1", "Hey2")) + 
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()
dev.off()


#-----------------------------------------------------------------------------
# 8. Find neighbours and clusters
#-----------------------------------------------------------------------------
Cont_Dll4_new <- FindNeighbors(Cont_Dll4_new, dims = 1:20, reduction = "pca")
Cont_Dll4_new <- FindClusters(Cont_Dll4_new,
                          resolution = c(0.4, 0.6, 0.8, 1.0, 1.2))




################## CHECK UMAP WITH RESOLUTION 0.6 ##################

#-----------------------------------------------------------------------------
# 10. Run non-linear dimensional reduction using cosine metric
#-----------------------------------------------------------------------------
# prepare UMAP
Idents(Cont_Dll4_new) <- 'RNA_snn_res.0.6'
Cont_Dll4_new <- RunUMAP(Cont_Dll4_new, dims = 1:20)

DimPlot(Cont_Dll4_new, reduction = "umap")


#-----------------------------------------------------------------------------
# Check sample distribution in UMAP by genotype
#-----------------------------------------------------------------------------
Cont_Dll4_new_samples <- SetIdent(Cont_Dll4_new,
                                  value = Cont_Dll4_new@meta.data$Genotype)
pdf("UMAP_All cells_Samples dist_n.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(Cont_Dll4_new_samples, reduction = "umap", cols = c("green", "red"))
dev.off()



#-----------------------------------------------------------------------------
# 11. Find cluster markers to annotate clusters (report only the positive ones)
# and with a adjusted p-value < 0.05
#-----------------------------------------------------------------------------
Cont_Dll4_new.markers <- FindAllMarkers(Cont_Dll4_new, only.pos = TRUE)
Cont_Dll4_new.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05)

Cont_Dll4_new.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10


# plot heatmap of top 10 genes of each cluster
pdf("HeatmapCluster_AllGenes_res0.6.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DoHeatmap(Cont_Dll4_new, features = top10$gene, size = 4)  + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(text = element_text(size = 7))
dev.off()


#-----------------------------------------------------------------------------
# 12. Attribute new cluster names based on top markers
#-----------------------------------------------------------------------------
new.cluster.ids0.6 <- c("Capillaries 1", "Capillaries 2", "Endocardium",
                        "Lymphatic", "Artery", "Vein", "Pericytes")
names(new.cluster.ids0.6) <- levels(Cont_Dll4_new)
Cont_Dll4_new <- RenameIdents(Cont_Dll4_new, new.cluster.ids0.6)

# Plot UMAP of all ECs with resolution 0.6 and with cluster labeling
pdf("UMAP_AllCells_res0.6.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")
# Paper size
DimPlot(Cont_Dll4_new, reduction = "umap", label = FALSE, pt.size = 1,
        cols = c('Capillaries 1' = '#e49400', 'Capillaries 2' = '#C1B80C',
                 'Endocardium' = 'purple', 'Lymphatic' = '#009dff',
                 'Artery' = '#FC0808', 'Vein' = '#0E47D8',
                 'Pericytes' = '#A80519'))
dev.off()


# Plot stacked violin plot of cluster from resolution 0.6, for 2 selected 
# top10 markers
pdf("ViolinPlot_Top genes_res0.6.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

VlnPlot(Cont_Dll4_new, features <- c('Cd36', 'Mgll', 'Cd34', 'Sparc',
                                     'Npr3','Tmem108', 'Prox1', 'Foxp2',
                                     'Sox17', 'Efhd1', 'Zbtb7c', 'Sema3d',
                                     'Rgs4', 'Notch3'),
        group.by = 'RNA_snn_res.0.6', stack = TRUE, flip = TRUE)
dev.off()


# Plot dotplot of cluster from resolution 0.6, for 2 selected top10 markers
pdf("DotPlot_Top genes_res0.6.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DotPlot(Cont_Dll4_new, features = c('Cd36', 'Mgll', 'Cd34', 'Sparc',
                                    'Npr3','Tmem108', 'Prox1', 'Foxp2',
                                    'Sox17', 'Efhd1', 'Zbtb7c', 'Sema3d',
                                    'Rgs4', 'Notch3')) + 
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()
dev.off()


################## CHECK UMAP WITH RESOLUTION 0.4 ##################

#-----------------------------------------------------------------------------
# 10. Run non-linear dimensional reduction
#-----------------------------------------------------------------------------
# prepare UMAP
Idents(Cont_Dll4_new) <- 'RNA_snn_res.0.4'
Cont_Dll4_new <- RunUMAP(Cont_Dll4_new, dims = 1:20)

DimPlot(Cont_Dll4_new, reduction = "umap")


#-----------------------------------------------------------------------------
# 11. Find cluster markers to annotate clusters (report only the positive ones)
#-----------------------------------------------------------------------------
Cont_Dll4_new.markers <- FindAllMarkers(Cont_Dll4_new, only.pos = TRUE)
Cont_Dll4_new.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05)

Cont_Dll4_new.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

# plot heatmap of top 10 genes of each cluster
pdf("HeatmapCluster_AllGenes_res0.4.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DoHeatmap(Cont_Dll4_new, features = top10$gene, size = 4)  + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(text = element_text(size = 7))
dev.off()


#-----------------------------------------------------------------------------
# 12. Attribute new cluster names based on top markers
#-----------------------------------------------------------------------------
new.cluster.ids0.4 <- c("Capillaries 1", "Activated Capillaries 2",
                        "Endocardium", "Lymphatic")
names(new.cluster.ids0.4) <- levels(Cont_Dll4_new)
Cont_Dll4_new <- RenameIdents(Cont_Dll4_new, new.cluster.ids0.4)
DimPlot(Cont_Dll4_new, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend()

# Plot UMAP of all ECs with resolution 0.4 and with cluster labeling
pdf("UMAP_AllCells_res0.4.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(Cont_Dll4_new, reduction = "umap", label = FALSE, pt.size = 1,
        cols = c('Capillaries 1' = '#e49400', 'Capillaries 2' = '#C1B80C',
                 'Endocardium' = 'purple', 'Lymphatic' = '#009dff'))
dev.off()


# Plot stacked violin plot of cluster from resolution 0.4, for 2 selected 
# top10 markers
pdf("ViolinPlot_Top genes_res0.4.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
VlnPlot(Cont_Dll4_new, features <- c('Cd36', 'Mgll', 'Cd34', 'Sparc',
                                     'Npr3','Tmem108', 'Prox1', 'Foxp2'),
        group.by = 'RNA_snn_res.0.4', stack = TRUE, flip = TRUE)
dev.off()


# Plot dotplot of cluster from resolution 0.4, for 2 selected top10 markers
pdf("DotPlot_Top genes_res0.4.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DotPlot(Cont_Dll4_new, features = c('Cd36', 'Mgll', 'Cd34', 'Sparc',
                                    'Npr3','Tmem108', 'Prox1', 'Foxp2')) + 
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()
dev.off()



################## CHECK UMAP WITH RESOLUTION 0.8 ##################

#-----------------------------------------------------------------------------
# 10. Run non-linear dimensional reduction
#-----------------------------------------------------------------------------
# prepare UMAP
Idents(Cont_Dll4_new) <- 'RNA_snn_res.0.8'
Cont_Dll4_new <- RunUMAP(Cont_Dll4_new, dims = 1:20)

DimPlot(Cont_Dll4_new, reduction = "umap")


#-----------------------------------------------------------------------------
# 11. Find cluster markers to annotate clusters (report only the positive ones)
#-----------------------------------------------------------------------------
Cont_Dll4_new.markers <- FindAllMarkers(Cont_Dll4_new, only.pos = TRUE)
Cont_Dll4_new.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05)

Cont_Dll4_new.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

# plot heatmap of top 10 genes of each cluster
pdf("HeatmapCluster_AllGenes_res0.8.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DoHeatmap(Cont_Dll4_new, features = top10$gene, size = 4)  + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(text = element_text(size = 7))
dev.off()


#-----------------------------------------------------------------------------
# 12. Attribute new cluster names based on top markers
#-----------------------------------------------------------------------------
new.cluster.ids0.8 <- c("Endocardium", "Capillaries1", "Capillaries3",
                        "Capillaries2", "Lymphatics", "Arteries",
                         "Veins", "Pericytes")
names(new.cluster.ids0.8) <- levels(Cont_Dll4_new)
Cont_Dll4_new <- RenameIdents(Cont_Dll4_new, new.cluster.ids0.8)
DimPlot(Cont_Dll4_new, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


# Plot UMAP of all ECs with resolution 0.8 and with cluster labeling
pdf("UMAP_AllGenes_res0.8.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(Cont_Dll4_new, reduction = "umap", label = FALSE, pt.size = 1,
        cols = c('Endocardium' = 'purple', 'Capillaries1' = '#C1B80C',
                 'Capillaries3' = '#e49400', "Capillaries2" = '#2AD868',
                 'Lymphatics' = '#009dff', 'Arteries' = '#FC0808',
                 'Veins' = '#0E47D8', 'Pericytes' = '#A80519'))
dev.off()


# Plot stacked violin plot of cluster from resolution 0.8, for 2 selected 
# top10 markers
pdf("ViolinPlot_Top genes_res0.8.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
VlnPlot(Cont_Dll4_new, features <- c('Cd36', 'Mgll', 'Timp4', 'Cd34', 'Sparc',
                                     'Npr3','Tmem108', 'Prox1', 'Foxp2',
                                     'Sox17', 'Glul', 'Zbtb7c', 'Sema3d',
                                     'Rgs4', 'Notch3'),
        group.by = 'RNA_snn_res.0.8', stack = TRUE, flip = TRUE)
dev.off()


# Plot dotplot of cluster from resolution 0.8, for 2 selected top10 markers
pdf("DotPlot_Top genes_res0.8.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DotPlot(Cont_Dll4_new, features = c('Cd36', 'Mgll', 'Timp4', 'Cd34', 'Sparc',
                                    'Npr3','Tmem108', 'Prox1', 'Foxp2')) + 
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()
dev.off()


#-----------------------------------------------------------------------------
# 13. Create metadata column with cluster ids assigned to each cell
#-----------------------------------------------------------------------------
Cont_Dll4_new@meta.data[["SR_res0.8_euALL"]] <- Cont_Dll4_new@active.ident




#-----------------------------------------------------------------------------
# 14. Confirm expression of known cell type-specific markers with
# identified/annotated clusters
#-----------------------------------------------------------------------------

pdf("Dot_Feature_Cluster markers_All Cells.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

# Endocardial markers
DotPlot(Cont_Dll4_new, features = c("Npr3", "Cytl1")) +
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06")
FeaturePlot(Cont_Dll4_new, features = c("Npr3", "Cytl1"),
            cols = c('#e1e1e1', '#fee391', '#fe9929', '#d25706', '#7e2c06'))
                                                                            
# Endothelial markers
DotPlot(Cont_Dll4_new, features = c("Fabp4", "Cd36")) +
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06")
FeaturePlot(Cont_Dll4_new, features = c("Fabp4", "Cd36"),
            cols = c('#e1e1e1', '#fee391', '#fe9929', '#d25706', '#7e2c06'))
                                                                            
# Lymphatics markers
DotPlot(Cont_Dll4_new, features = c("Lyve1", "Flt4", "Prox1")) +
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06")
FeaturePlot(Cont_Dll4_new,
            features = c("Lyve1", "Flt4", "Prox1"),
            cols = c('#e1e1e1', '#fee391', '#fe9929', '#d25706', '#7e2c06'))
                                                                                                     
# Pericyte markers
DotPlot(Cont_Dll4_new, features = c("Rgs4", "Cspg4", "Pdgfrb")) +
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06")
FeaturePlot(Cont_Dll4_new, features = c("Rgs4", "Cspg4", "Pdgfrb"),
            cols = c('#e1e1e1', '#fee391', '#fe9929', '#d25706', '#7e2c06'))
dev.off()




#-----------------------------------------------------------------------------
# 15. Confirm proportion of control and Dll4-iDEC cells in each cluster
#-----------------------------------------------------------------------------

pdf("Barplot with all cells_res0.8_euc.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

dittoBarPlot(Cont_Dll4_new, var = "SR_res0.8_euALL", group.by = "Genotype",
             color.panel = c('purple', '#C1B80C', '#2AD868', '#e49400', '#009dff',
                                '#FC0808', '#0E47D8', '#A80519'),
                                var.labels.reorder = c(5,2,3,4,6,1,8,7))
dev.off()



#-----------------------------------------------------------------------------
# 16. Subset Seurat object on cell types
#-----------------------------------------------------------------------------

# set annotations of 0.8res euclidean metric as current Identity
Cont_Dll4_new <- SetIdent(Cont_Dll4_new,
                          value = Cont_Dll4_new@meta.data$SR_res0.8_euALL)

# remove non-coronary ECs
Cont_Dll4_EC <- subset(x = Cont_Dll4_new, idents = c("Endocardium", "Pericytes",
                                                 "Lymphatics"), invert = TRUE)

# Keep only lymphatics
Cont_Dll4_Lymph <- subset(x = Cont_Dll4_new, idents = "Lymphatics") 

# Keep only Endocardium
Cont_Dll4_Endo <- subset(x = Cont_Dll4_new, idents = "Endocardium") 

# Keep only Pericytes
Cont_Dll4_Peri <- subset(x = Cont_Dll4_new, idents = "Pericytes")



#-----------------------------------------------------------------------------
# 17. Check Dll4 and downstream targets of Notch signalling in each of the
# cell type subsets created
#-----------------------------------------------------------------------------
# analysis of Dll4 and Notch component pathway members in ECs
Idents(Cont_Dll4_EC) <- 'Genotype'

# Violin plot and dotplot of Dll4, Hes1, Hey1/2 in ECs
pdf("Endothelial only pop_notch pathway.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
VlnPlot(Cont_Dll4_EC, features <- c('Dll4','Hes1', 'Hey1', 'Hey2'),
        group.by = 'Genotype', col = c('green', 'red'))

DotPlot(Cont_Dll4_EC, features = c('Dll4','Hes1', 'Hey1', 'Hey2')) + 
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()
dev.off()

# absolute number of control and mutant cells in EC subset
table(Cont_Dll4_EC@meta.data$Genotype)



# analysis of Dll4 and Notch component pathway members in Lymphatics
Idents(Cont_Dll4_Lymph) <- 'Genotype'

# Violin plot and dotplot of Dll4, Hes1, Hey1/2 in Lymphatics
pdf("Lymphatic only pop_notch pathway.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
VlnPlot(Cont_Dll4_Lymph, features <- c('Dll4','Hes1', 'Hey1', 'Hey2'),
        group.by = 'Genotype', col = c('green', 'red'))

DotPlot(Cont_Dll4_Lymph, features = c('Dll4','Hes1', 'Hey1', 'Hey2')) + 
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()
dev.off()

# absolute number of control and mutant cells in Lymphatics subset
table(Cont_Dll4_Lymph@meta.data$Genotype)



# analysis of Dll4 and Notch component pathway members in Endocardium
Idents(Cont_Dll4_Endo) <- 'Genotype'

# Violin plot and dotplot of Dll4, Hes1, Hey1/2 in Endocardium
pdf("Endocardial only pop_notch pathway_4.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
VlnPlot(Cont_Dll4_Endo, features <- c('Dll4','Hes1', 'Hey1', 'Hey2'),
        group.by = 'Genotype', col = c('green', 'red'))

DotPlot(Cont_Dll4_Endo, features = c('Dll4','Hes1', 'Hey1', 'Hey2')) + 
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()
dev.off

# absolute number of control and mutant cells in Endocardium subset
table(Cont_Dll4_Endo@meta.data$Genotype)



# analysis of Dll4 and Notch component pathway members in Pericytes
Idents(Cont_Dll4_Peri) <- 'Genotype'

# Violin plot and dotplot of Dll4, Hes1, Hey1/2 in Pericytes
pdf("Pericyte only pop_notch pathway.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
VlnPlot(Cont_Dll4_Peri, features <- c('Dll4','Hes1', 'Hey1', 'Hey2'),
        group.by = 'Genotype', col = c('green', 'red'))

DotPlot(Cont_Dll4_Peri, features = c('Dll4','Hes1', 'Hey1', 'Hey2')) + 
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()
dev.off()

# absolute number of control and mutant cells in Pericytes subset
table(Cont_Dll4_Peri@meta.data$Genotype)


#-----------------------------------------------------------------------------
# 18. Save endocardial population and EC population in separate .rds files
#-----------------------------------------------------------------------------
saveRDS(Cont_Dll4_EC,
        file='S:/LAB_RB/LAB/Susana/TFM/Heart scRNAseq objects - to train/ContDll4_EC_1/Cont_Dll4_EC.rds')
saveRDS(Cont_Dll4_Endo,
        file='S:/LAB_RB/LAB/Susana/TFM/Heart scRNAseq objects - to train/ContDll4_EC_1/Cont_Dll4_Endo.rds')













