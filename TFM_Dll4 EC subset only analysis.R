############################################################################
## ANALYSIS OF CORONARY ECS ONLY SUBSET
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
library(org.Hs.eg.db)
library(Orthology.eg.db)
library(nichenetr)


#-----------------------------------------------------------------------------
# 1.Load and read .rds files
#-----------------------------------------------------------------------------
file3 <- choose.files()
Cont_Dll4_EC <- readRDS(file3)
Cont_Dll4_EC



#-----------------------------------------------------------------------------
# 2. Find variable features again for just the EC population without
# without 'contaminant' lymphatics, endocardium and pericytes
#-----------------------------------------------------------------------------
# store gene names of genes expressed in ECs
AllGenes_EC <- rownames(Cont_Dll4_EC)
Cont_Dll4_EC <- FindVariableFeatures(Cont_Dll4_EC)

# Identify the 10 most highly variable genes
top10_EC <- head(VariableFeatures(Cont_Dll4_EC), 10)

# plot variable features with and without labels
pdf("VariableFeatures_ECs.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
plot1 <- VariableFeaturePlot(Cont_Dll4_EC)
plot2 <- LabelPoints(plot = plot1, points = top10_EC, repel = TRUE)
plot1 + plot2
dev.off()


#-----------------------------------------------------------------------------
# 2. Scale data
#-----------------------------------------------------------------------------
Cont_Dll4_EC <- ScaleData(Cont_Dll4_EC)


#-----------------------------------------------------------------------------
# 3. PCA analysis of dataset
#-----------------------------------------------------------------------------
Cont_Dll4_EC <- RunPCA(Cont_Dll4_EC)

# generate plots/heatmaps to decide on PC for dimension reduction
pdf("PCA_ECs.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(Cont_Dll4_EC, reduction = "pca") + NoLegend()
DimHeatmap(Cont_Dll4_EC, dims = 1:24, cells = 500, balanced = TRUE)
ElbowPlot(Cont_Dll4_EC, ndims = 30, reduction = "pca")
dev.off()


#----------------------------------------------------------------------------
# 4. Find neighbours and clusters
#----------------------------------------------------------------------------
Cont_Dll4_EC <- FindNeighbors(Cont_Dll4_EC, dims = 1:20, reduction = "pca")
Cont_Dll4_EC <- FindClusters(Cont_Dll4_EC,
                             resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))


################## CHECK UMAP WITH RESOLUTION 0.4 ##################

#----------------------------------------------------------------------------
# 5. Run non-linear dimensional reduction
#----------------------------------------------------------------------------
# prepare UMAP
Idents(Cont_Dll4_EC) <- 'RNA_snn_res.0.4'
Cont_Dll4_EC <- RunUMAP(Cont_Dll4_EC, dims = 1:20)

DimPlot(Cont_Dll4_EC, reduction = "umap")



#----------------------------------------------------------------------------
# Check sample distribution in UMAP by genotype
#----------------------------------------------------------------------------
Cont_Dll4_EC_samples <- SetIdent(Cont_Dll4_EC,
                                 value = Cont_Dll4_EC@meta.data$Genotype)
pdf("UMAP_ECs_Samples dist_final.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(Cont_Dll4_EC_samples, reduction = "umap", cols = c("green", "red"))
dev.off()


#----------------------------------------------------------------------------
# 6. Find cluster markers to annotate clusters (report only the positive ones)
#----------------------------------------------------------------------------
Cont_Dll4_EC.markers <- FindAllMarkers(Cont_Dll4_EC, only.pos = TRUE)
Cont_Dll4_EC.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05)

Cont_Dll4_EC.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10


# plot heatmap of top 10 genes of each cluster
pdf("HeatmapCluster_ECs_res0.4def.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DoHeatmap(Cont_Dll4_EC, features = top10$gene, size = 4)  + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(text = element_text(size = 7))
dev.off()





#----------------------------------------------------------------------------
# 7. Attribute new cluster names based on top markers
#----------------------------------------------------------------------------

new.cluster.ids_EC0.4 <- c("Capillaries1", "Capillaries2", "Artery")
names(new.cluster.ids_EC0.4) <- levels(Cont_Dll4_EC)
Cont_Dll4_EC <- RenameIdents(Cont_Dll4_EC, new.cluster.ids_EC0.4)
DimPlot(Cont_Dll4_EC, reduction = "umap", label = TRUE, pt.size = 0.5) +
  NoLegend()


# Plot UMAP of all ECs with resolution 0.4 and with cluster labeling
pdf("UMAP_ECs_Res0.4def.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(Cont_Dll4_EC, reduction = "umap", pt.size = 0.7,
        cols = c('Capillaries1' = '#45b05c', 'Capillaries2' = '#009dff',
                 'Artery' = '#bf2630'))
dev.off()


# Plot stacked violin plot of cluster from resolution 0.4, for 2 selected 
# top10 markers
pdf("ViolinPlot_Top EC genes_res0.4.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
VlnPlot(Cont_Dll4_EC, features <- c('Ptprb', 'Flt1', 'Sparc', 'Cd34',
                                    'Glul', 'Efhd1'),
        group.by = 'RNA_snn_res.0.4', stack = TRUE, flip = TRUE)
dev.off()


# Plot dotplot of cluster from resolution 0.4, for 2 selected top10 markers
pdf("DotPlot_Top EC genes_res1.2.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DotPlot(Cont_Dll4_EC, features = c('Ptprb', 'Flt1', 'Sparc', 'Cd34',
                                   'Glul', 'Efhd1')) + 
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()
dev.off()




#----------------------------------------------------------------------------
# 8. Do feature and dotplots with known markers of EC subtypes to 
# confirm annotation
#----------------------------------------------------------------------------

pdf("UMAP_ECs_Markers_ALL.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
FeaturePlot(Cont_Dll4_EC, features = c('Bmx','Gja5',"Fbln5", "Gja4", "Hey1", 
                                       "Sox17", 'Efnb2', 'Dach1',
                                       "Vwf", "Vcam1", "Nr2f2",
                                       "Mki67", "Cdca8", "Ccnb2",
                                       "Apln", "Kcne3", "Esm1", "Odc1"),
            cols = c('#e1e1e1', '#fee391', '#fe9929', '#d25706','#7e2c06'))
                                                       
dev.off()                              


pdf("DotPlot_ECs_Markers_ALL.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DotPlot(Cont_Dll4_EC, features = c('Bmx','Gja5',"Fbln5", "Gja4", "Hey1", 
                                   "Sox17", 'Efnb2', 'Dach1',
                                   "Vwf", "Vcam1", "Nr2f2",
                                   "Mki67", "Cdca8", "Ccnb2",
                                   "Apln", "Kcne3", "Esm1", "Odc1"),
        group.by = 'Genotype') + 
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()
dev.off()

                          
################## CHECK UMAP WITH RESOLUTION 0.8 ##################                              
#---------------------------------------------------------------------------
# 5. prepare UMAP
#---------------------------------------------------------------------------
Idents(Cont_Dll4_EC) <- 'RNA_snn_res.0.8'
Cont_Dll4_EC <- RunUMAP(Cont_Dll4_EC, dims = 1:20)
                              
DimPlot(Cont_Dll4_EC, reduction = "umap")


#---------------------------------------------------------------------------
# 6. Find cluster markers to annotate clusters (report only the positive ones)
#---------------------------------------------------------------------------
Cont_Dll4_EC.markers <- FindAllMarkers(Cont_Dll4_EC, only.pos = TRUE)
Cont_Dll4_EC.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05)

Cont_Dll4_EC.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

pdf("HeatmapCluster_ECs_res0.8_def.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DoHeatmap(Cont_Dll4_EC, features = top10$gene, size = 4)  + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(text = element_text(size = 7))
dev.off()


#---------------------------------------------------------------------------
# 7. Attribute new cluster names based on top markers
#---------------------------------------------------------------------------

new.cluster.ids_EC0.8 <- c("Capillaries1", "Capillaries2", "Angiogenic",
                           "Activated Cap", "Artery1", "Vein",
                           "Artery2", "Interferon")
names(new.cluster.ids_EC0.8) <- levels(Cont_Dll4_EC)
Cont_Dll4_EC <- RenameIdents(Cont_Dll4_EC, new.cluster.ids_EC0.8)
DimPlot(Cont_Dll4_EC, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend()


# Plot UMAP of all ECs with resolution 0.8 and with cluster labeling
pdf("UMAP_ECs_Res0.8 def.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(Cont_Dll4_EC, reduction = "umap", pt.size = 0.7,
        cols = c('Capillaries1' = '#e6cc00', 'Capillaries2' = '#009dff','Angiogenic' = '#f21bcb',
                 'Activated Cap' = '#f2981b', 'Artery1' = '#bf2630',
                 'Vein' = '#1b5bf2', 'Artery2' = '#fa4854', 'Interferon' = '#45b05c'))
dev.off()


# Plot stacked violin plot of cluster from resolution 0.8, for 2 selected 
# top10 markers
pdf("ViolinPlot_Top EC genes_res0.8.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
VlnPlot(Cont_Dll4_EC, features <- c('Ptprb', 'Flt1', 'Apln', 'Mmp14',
                                    'Sparc', 'Cd34', 'Glul', 'Efhd1',
                                    'Sema3d', 'Nr2f2', 'Gja5', 'Fbln5',
                                    'Ifi206', 'Ifit2'),
        group.by = 'RNA_snn_res.0.8', stack = TRUE, flip = TRUE)
dev.off()


# Plot dotplot of cluster from resolution 0.4, for 2 selected top10 markers
pdf("DotPlot_Top EC genes_res0.8.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DotPlot(Cont_Dll4_EC, features = c('Ptprb', 'Flt1', 'Apln', 'Mmp14',
                                   'Sparc', 'Cd34', 'Glul', 'Efhd1',
                                   'Sema3d', 'Nr2f2', 'Gja5', 'Fbln5',
                                   'Ifi206', 'Ifit2')) + 
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()
dev.off()




################## CHECK UMAP WITH RESOLUTION 0.6 ##################                              
#---------------------------------------------------------------------------
# 5. prepare UMAP
#---------------------------------------------------------------------------                     
# prepare UMAP
Idents(Cont_Dll4_EC) <- 'RNA_snn_res.0.6'
Cont_Dll4_EC <- RunUMAP(Cont_Dll4_EC, dims = 1:20)
                              
DimPlot(Cont_Dll4_EC, reduction = "umap")
                              


#---------------------------------------------------------------------------
# 6. Find cluster markers to annotate clusters (report only the positive ones)
#---------------------------------------------------------------------------
Cont_Dll4_EC.markers <- FindAllMarkers(Cont_Dll4_EC, only.pos = TRUE)
Cont_Dll4_EC.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05)

Cont_Dll4_EC.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10


pdf("HeatmapCluster_ECs_res0.6def.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DoHeatmap(Cont_Dll4_EC, features = top10$gene, size = 4)  + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(text = element_text(size = 7))
dev.off()


#---------------------------------------------------------------------------
# 7. Attribute new cluster names based on top markers
#---------------------------------------------------------------------------
new.cluster.ids_EC0.6 <- c("Capillaries1", 'Angiogenic', "Capillaries2",
                           "Artery 1", "Artery 2", "Vein", "Interferon")
names(new.cluster.ids_EC0.6) <- levels(Cont_Dll4_EC)
Cont_Dll4_EC <- RenameIdents(Cont_Dll4_EC, new.cluster.ids_EC0.6)
DimPlot(Cont_Dll4_EC, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend()


pdf("UMAP_ECs_Res0.6 def.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(Cont_Dll4_EC, reduction = "umap", pt.size = 0.7,
        cols = c('Capillaries1' = '#e6cc00', 'Angiogenic' = '#f21bcb',
                 'Capillaries2' = '#f2981b', 'Artery 1' = '#bf2630',
                 'Artery 2' = '#fa4854', 'Vein' = '#1b5bf2',
                 'Interferon' = '#45b05c'))

dev.off()





########## CHECK UMAP WITH RESOLUTION 0.6, spread=0.1, distmin=0.3 ############

#---------------------------------------------------------------------------
# 5. Run non-linear dimensional reduction
#---------------------------------------------------------------------------
# prepare UMAP
Idents(Cont_Dll4_EC) <- 'RNA_snn_res.0.6'
Cont_Dll4_EC <- RunUMAP(Cont_Dll4_EC, dims = 1:20, min.dist = 0.3, spread = 0.1)

DimPlot(Cont_Dll4_EC, reduction = "umap")

#---------------------------------------------------------------------------
# Visualize new UMAP with annotation from res=0.6
#---------------------------------------------------------------------------
new.cluster.ids_EC0.6 <- c("Capillaries1", 'Angiogenic', "Capillaries2",
                           "Artery 1", "Artery 2", "Vein", "Interferon")
names(new.cluster.ids_EC0.6) <- levels(Cont_Dll4_EC)
Cont_Dll4_EC <- RenameIdents(Cont_Dll4_EC, new.cluster.ids_EC0.6)
DimPlot(Cont_Dll4_EC, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend()

pdf("UMAP_ECs_Res0.6_03_01.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(Cont_Dll4_EC, reduction = "umap", pt.size = 0.7,
        cols = c('Capillaries1' = '#e6cc00', 'Angiogenic' = '#f21bcb',
                 'Capillaries2' = '#f2981b', 'Artery 1' = '#bf2630',
                 'Artery 2' = '#fa4854', 'Vein' = '#1b5bf2',
                 'Interferon' = '#45b05c'))
dev.off()



########## CHECK UMAP WITH RESOLUTION 0.6, spread=0.3, distmin=0.5 ############

#---------------------------------------------------------------------------
# 5. Run non-linear dimensional reduction
#---------------------------------------------------------------------------
# prepare UMAP
Idents(Cont_Dll4_EC) <- 'RNA_snn_res.0.6'
Cont_Dll4_EC <- RunUMAP(Cont_Dll4_EC, dims = 1:20, min.dist = 0.3, spread = 0.5)

DimPlot(Cont_Dll4_EC, reduction = "umap")

#---------------------------------------------------------------------------
# Visualize new UMAP with annotation from res=0.6
#---------------------------------------------------------------------------
new.cluster.ids_EC0.6 <- c("Capillaries1", 'Angiogenic', "Capillaries2",
                           "Artery 1", "Artery 2", "Vein", "Interferon")
names(new.cluster.ids_EC0.6) <- levels(Cont_Dll4_EC)
Cont_Dll4_EC <- RenameIdents(Cont_Dll4_EC, new.cluster.ids_EC0.6)
DimPlot(Cont_Dll4_EC, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

pdf("UMAP_ECs_Res0.6_03_05.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(Cont_Dll4_EC, reduction = "umap", pt.size = 0.7,
        cols = c('Capillaries1' = '#e6cc00', 'Angiogenic' = '#f21bcb',
                 'Capillaries2' = '#f2981b', 'Artery 1' = '#bf2630',
                 'Artery 2' = '#fa4854', 'Vein' = '#1b5bf2',
                 'Interferon' = '#45b05c'))

dev.off()


########## CHECK UMAP WITH RESOLUTION 0.6, spread=0.1, distmin=0.1 ############

#---------------------------------------------------------------------------
# 5. Run non-linear dimensional reduction
#--------------------------------------------------------------------------
# prepare UMAP
Idents(Cont_Dll4_EC) <- 'RNA_snn_res.0.6'
Cont_Dll4_EC <- RunUMAP(Cont_Dll4_EC, dims = 1:20, min.dist = 0.1, spread = 0.1)

DimPlot(Cont_Dll4_EC, reduction = "umap")

#---------------------------------------------------------------------------
# Visualize new UMAP with annotation from res=0.6
#---------------------------------------------------------------------------
new.cluster.ids_EC0.6 <- c("Capillaries1", 'Angiogenic', "Capillaries2",
                           "Artery 1", "Artery 2", "Vein", "Interferon")
names(new.cluster.ids_EC0.6) <- levels(Cont_Dll4_EC)
Cont_Dll4_EC <- RenameIdents(Cont_Dll4_EC, new.cluster.ids_EC0.6)
DimPlot(Cont_Dll4_EC, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

pdf("UMAP_ECs_Res0.6_03_05.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(Cont_Dll4_EC, reduction = "umap", pt.size = 0.7,
        cols = c('Capillaries1' = '#e6cc00', 'Angiogenic' = '#f21bcb',
                 'Capillaries2' = '#f2981b', 'Artery 1' = '#bf2630',
                 'Artery 2' = '#fa4854', 'Vein' = '#1b5bf2',
                 'Interferon' = '#45b05c'))

dev.off()




########## CHECK UMAP WITH RESOLUTION 0.6, spread=0.4, distmin=0.2 ############

#---------------------------------------------------------------------------
# 5. Run non-linear dimensional reduction
#--------------------------------------------------------------------------
# prepare UMAP
Idents(Cont_Dll4_EC) <- 'RNA_snn_res.0.6'
Cont_Dll4_EC <- RunUMAP(Cont_Dll4_EC, dims = 1:20, min.dist = 0.4, spread = 0.2)

DimPlot(Cont_Dll4_EC, reduction = "umap")

#---------------------------------------------------------------------------
# Visualize new UMAP with annotation from res=0.6
#---------------------------------------------------------------------------
new.cluster.ids_EC0.6 <- c("Capillaries1", 'Angiogenic', "Capillaries2",
                           "Artery 1", "Artery 2", "Vein", "Interferon")
names(new.cluster.ids_EC0.6) <- levels(Cont_Dll4_EC)
Cont_Dll4_EC <- RenameIdents(Cont_Dll4_EC, new.cluster.ids_EC0.6)
DimPlot(Cont_Dll4_EC, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

pdf("UMAP_ECs_Res0.6_04_02.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(Cont_Dll4_EC, reduction = "umap", pt.size = 0.7,
        cols = c('Capillaries1' = '#e6cc00', 'Angiogenic' = '#f21bcb',
                 'Capillaries2' = '#f2981b', 'Artery 1' = '#bf2630',
                 'Artery 2' = '#fa4854', 'Vein' = '#1b5bf2',
                 'Interferon' = '#45b05c'))

dev.off()


#---------------------------------------------------------------------------
# 8. Violin and dot plots of 2 top10 markers of each EC cluster
#---------------------------------------------------------------------------
pdf("ViolinPlot_Top EC genes_res0.6_new.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
VlnPlot(Cont_Dll4_EC, features <- c('Ptprb', 'Flt1', 'Apln', 'Lgals1', 'Cd34',
                                    'Glul', 'Efhd1', 'Gja5', 'Fbln5',
                                    'Sema3d', 'Vcam1', 'Ifi206', 'Ifit3b'),
        group.by = 'RNA_snn_res.0.6', stack = TRUE, flip = TRUE)
dev.off()


pdf("DotPlot_Top EC genes_res0.6_new.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DotPlot(Cont_Dll4_EC, features = c('Ptprb', 'Flt1', 'Apln', 'Lgals1', 'Cd34',
                                    'Glul', 'Efhd1', 'Gja5', 'Fbln5',
                                    'Sema3d', 'Vcam1', 'Ifi206', 'Ifit3b')) + 
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()
dev.off()


#---------------------------------------------------------------------------
# 9. Feature and dotplots of known markers of EC subtypes
#---------------------------------------------------------------------------
pdf("Feature_ECs_Markers_ALL final.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
FeaturePlot(Cont_Dll4_EC, features = c('Bmx','Gja5',"Fbln5", "Gja4", "Hey1", 
                                       "Sox17", 'Efnb2', 'Dach1',
                                       "Vwf", "Vcam1", "Nr2f2",
                                       "Mki67", "Cdca8", "Ccnb2",
                                       "Apln", "Kcne3", "Esm1", "Odc1"),
            cols = c('#e1e1e1', '#fee391', '#fe9929', '#d25706',
                              '#7e2c06'))
dev.off()                              
      

                        
pdf("DotPlot_ECs_Markers_ALL final.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size                              
DotPlot(Cont_Dll4_EC, features = c('Bmx','Gja5',"Fbln5", "Gja4", "Hey1",
                                         "Sox17", 'Efnb2', 'Dach1',
                                         "Vwf", "Vcam1", "Nr2f2",
                                         "Mki67", "Cdca8", "Ccnb2",
                                         "Apln", "Kcne3", "Esm1", "Odc1"),
              group.by = 'Genotype') +
        scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
        coord_flip()
dev.off()




#-----------------------------------------------------------------------------
# 10. Create metadata column with cluster ids assigned to each cell
#-----------------------------------------------------------------------------
Cont_Dll4_EC@meta.data[["SR_res_0.6 EC_cos"]] <- Cont_Dll4_EC@active.ident



#-----------------------------------------------------------------------------
# 11. find arterial DEGs without filtering for logfc threshold and cell minimums  
# to have all genes
#-----------------------------------------------------------------------------
Idents(Cont_Dll4_EC) <- "SR_res_0.6 EC_cos"
ContvsDll4_DEG <- FindMarkers(Cont_Dll4_EC, ident.1 = 'Artery 1',
                              ident.2 = 'Artery 2', logfc.threshold = 0,
                              min.cells.group = 1, min.pct = 0)

# find classical arterial markers
arterial_index <- which(rownames(ContvsDll4_DEG) %in% c('Bmx', 'Gja5', 'Fbln5',
                                                        'Gja4', 'Hey1', 'Sox17', 
                                                        'Efnb2', 'Dach1'))

# extract DEG data of specified arterial markers
arterial_FCs <- ContvsDll4_DEG[arterial_index,]

# Specify the desired order of columns
arterial_FCn <- as.data.frame(arterial_FCs)
arterial_FCn$pct.1 <- NULL
arterial_FCn$pct.2 <- NULL
arterial_FCn$p_val <- NULL

# order according to log2FC values
arterial_FCn <- arterial_FCn[order(arterial_FCn$avg_log2FC), ]

# save the df above as a table in a PDF
# Create a table
arterial_FCn$avg_log2FC <- round(arterial_FCn$avg_log2FC, 3)
arterial_FCn$p_val_adj <- format(arterial_FCn$p_val_adj, scientific = TRUE)
table_arterial_FCn <- tableGrob(arterial_FCn)

# Export the table to a PDF file
pdf("Arterial log2Fc table_1.pdf")
grid.arrange(table_arterial_FCn)
dev.off()


#-----------------------------------------------------------------------------
# 12. get number of control and mutant cells per cluster  
#-----------------------------------------------------------------------------
# get number of cells per cluster
n <- as.character(unique(Cont_Dll4_EC@meta.data$`SR_res_0.6 EC_cos`))
Clust <- c()
C_num <- c() # number of control cells
D_num <- c() # number of Dll4-iDEC cells

for(i in 1:length(n)){
  Clust <- c(Clust, n[i])
  C_num <- c(C_num, nrow(subset(Cont_Dll4_EC@meta.data, 
                                Cont_Dll4_EC@meta.data$`SR_res_0.6 EC_cos`==n[i] &
                                  Cont_Dll4_EC@meta.data$Genotype == 'Control')))
  D_num <- c(D_num, nrow(subset(Cont_Dll4_EC@meta.data,
                                Cont_Dll4_EC@meta.data$`SR_res_0.6 EC_cos`==n[i] &
                                  Cont_Dll4_EC@meta.data$Genotype == 'Dll4-iDEC')))
}

Cell_clust_freq <- cbind(C_num, D_num)
rownames(Cell_clust_freq) <- Clust
colnames(Cell_clust_freq) <- c('Control', 'Dll4-iDEC')

# Order the data frame based on row names
Cell_clust_freq <- Cell_clust_freq[order(rownames(Cell_clust_freq)), ]

# Export the table to a PDF file
Cell_clust_freq_g <- tableGrob(Cell_clust_freq)
pdf("Cell prop in EC original clusters table.pdf")
grid.arrange(Cell_clust_freq_g)
dev.off()


#-----------------------------------------------------------------------------
# 13. Do UMAP separated by genotype  
#-----------------------------------------------------------------------------
Idents(Cont_Dll4_EC) <- 'SR_res_0.6 EC_cos'
Cont_Dll4_EC <- RunUMAP(Cont_Dll4_EC, dims = 1:20, min.dist = 0.4, spread = 0.2)

DimPlot(Cont_Dll4_EC, reduction = "umap")


pdf("UMAP_ECs_Res0.6_04_02_by genotype.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
DimPlot(Cont_Dll4_EC, reduction = "umap", pt.size = 0.7,
        cols = c('Capillaries1' = '#e6cc00', 'Angiogenic' = '#f21bcb',
                 'Capillaries2' = '#f2981b', 'Artery 1' = '#bf2630',
                 'Artery 2' = '#fa4854', 'Vein' = '#1b5bf2',
                 'Interferon' = '#45b05c'), split.by = 'Genotype')

dev.off()


#-----------------------------------------------------------------------------
# 14. Do barplot separated by genotype with frequency of subtype
#-----------------------------------------------------------------------------
pdf("Barplot with all ECs_1.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

dittoBarPlot(Cont_Dll4_EC, var = "SR_res_0.6 EC_cos", group.by = "Genotype",
             color.panel = c('#e6cc00', '#f21bcb', '#f2981b', '#bf2630',
                                      '#fa4854', '#1b5bf2', '#45b05c'),
                                      var.labels.reorder = c(4,1,5,2,3,7,6))
dev.off()




#-----------------------------------------------------------------------------
# 15. create new column with manual annotated celltypes 
#-----------------------------------------------------------------------------
Cont_Dll4_EC@meta.data$SR_Manual_Clusters <- Cont_Dll4_EC@meta.data$'SR_res_0.6 EC_cos'


# JOIN Artery1 and Artery2 into a single 'Arterial' cluster

# Replace 'Artery 1' and 'Artery 2' in the 'metadata_column' with a 'Arterial'
# Convert the factor column to character
Cont_Dll4_EC@meta.data$SR_Manual_Clusters <- as.character(Cont_Dll4_EC@meta.data$SR_Manual_Clusters)

# Replace 'Artery 1/2 with 'Arterial'
A1 <- which(Cont_Dll4_EC@meta.data$SR_Manual_Clusters == 'Artery 1')
Cont_Dll4_EC@meta.data$SR_Manual_Clusters[A1] <- 'Arterial'

A2 <- which(Cont_Dll4_EC@meta.data$SR_Manual_Clusters == 'Artery 2')
Cont_Dll4_EC@meta.data$SR_Manual_Clusters[A2] <- 'Arterial'



# JOIN Capillaries1/2 and Angiogenic into a single 'Capillary' cluster

# Replace 'Capillaries 1/2 with 'Capillary'
C1 <- which(Cont_Dll4_EC@meta.data$SR_Manual_Clusters == 'Capillaries1')
Cont_Dll4_EC@meta.data$SR_Manual_Clusters[C1] <- 'Capillary'

C2 <- which(Cont_Dll4_EC@meta.data$SR_Manual_Clusters == 'Capillaries2')
Cont_Dll4_EC@meta.data$SR_Manual_Clusters[C2] <- 'Capillary'

A <- which(Cont_Dll4_EC@meta.data$SR_Manual_Clusters == 'Angiogenic')
Cont_Dll4_EC@meta.data$SR_Manual_Clusters[A] <- 'Capillary'

# Convert back to factor column if necessary
Cont_Dll4_EC@meta.data$SR_Manual_Clusters <- as.factor(Cont_Dll4_EC@meta.data$SR_Manual_Clusters)


#-----------------------------------------------------------------------------
# 16. create new column with manual annotated celltypes 
#-----------------------------------------------------------------------------
# create new metadata column with celltype and genotype to do DEG analysis
Cont_Dll4_EC@meta.data[["Celltype_genotype"]] <- paste0(Cont_Dll4_EC@meta.data$SR_Manual_Clusters,'_',Cont_Dll4_EC@meta.data$Genotype)




#-----------------------------------------------------------------------------
# 17. get number of control and mutant cells per cluster  
#-----------------------------------------------------------------------------
# make celltype and genotype the active idents
Idents(Cont_Dll4_EC) <- "Celltype_genotype"


table(Cont_Dll4_EC@meta.data$Celltype_genotype)
Control_Numb <- c(31, 284, 3, 6)
Dll4_Numb <- c(88, 274, 21, 20)
Cell_freq <- cbind(Control_Numb, Dll4_Numb)
rownames(Cell_freq) <- c('Artery', 'Capillary', 'Interferon', 'Vein')
colnames(Cell_freq) <- c('Control', 'Dll4-iDEC')

# Export the table to a PDF file
Cell_freq_g <- tableGrob(Cell_freq)
pdf("Cell prop in EC clusters table_1.pdf")
grid.arrange(Cell_freq_g)
dev.off()





#------------------------------------------------------------------------------
# 18. DEG ANALYSIS OF ARTERIAL CELLS
#------------------------------------------------------------------------------
Arterial_DEG <- FindMarkers(Cont_Dll4_EC, ident.1 = 'Arterial_Dll4-iDEC',
                                                 ident.2 = 'Arterial_Control')
head(Arterial_DEG, 20)
# store data of number total DEG, Up significant and down significant
ATotal <- nrow(Arterial_DEG)
nrow(subset(Arterial_DEG, Arterial_DEG$p_val_adj < 0.05))
A_UpSig <- nrow(subset(Arterial_DEG, Arterial_DEG$p_val_adj < 0.05 & Arterial_DEG$avg_log2FC > 0))
A_DwSig <- nrow(subset(Arterial_DEG, Arterial_DEG$p_val_adj < 0.05 & Arterial_DEG$avg_log2FC < 0))

Cluster <- 'Arterial'
Total_DEG <- ATotal
UpSig_DEG <- A_UpSig
DwSig_DEG <- A_DwSig

Arterial_DEG_genes <- rownames(Arterial_DEG)

write.csv(Arterial_DEG, file='Arterial DEGs.csv')

# Export the table to a PDF file
Art_DEG <- tableGrob(Arterial_DEG)
pdf("Arterial DEGs table.pdf")
grid.arrange(Art_DEG)
dev.off()



#------------------------------------------------------------------------------
# 19. ORA ANALYSIS OF ARTERIAL CELLS
#------------------------------------------------------------------------------
# ORA - Over-representation analysis
pdf("ORA_Arterial.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

Arterial_GO_resultsCC <- enrichGO(gene = Arterial_DEG_genes,
                               universe = AllGenes_EC, OrgDb = org.Mm.eg.db,
                               keyType = "SYMBOL", ont = "CC",
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05)

Art_GO_CC_table <- Arterial_GO_resultsCC@result

# Export the table to a PDF file
Art_GO_CC <- tableGrob(Art_GO_CC_table)
pdf("Arterial ORA CC table.pdf")
grid.arrange(Art_GO_CC)
dev.off()


par(mar=c(8,20,3,1))
barplot(Arterial_GO_resultsCC, showCategory = 30, cex.names = 0.3,
        title = 'Arterial cellular component')
dotplot(Arterial_GO_resultsCC, showCategory=30) +
  ggtitle("dotplot for Arterial ORA CC")


Arterial_GO_resultsBP <- enrichGO(gene = Arterial_DEG_genes,
                               universe = AllGenes_EC, OrgDb = org.Mm.eg.db,
                               keyType = "SYMBOL", ont = "BP",
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05)

Art_GO_BP_table <- Arterial_GO_resultsBP@result

# Export the table to a PDF file
Art_GO_BP <- tableGrob(Art_GO_BP_table)
pdf("Arterial ORA BP table.pdf")
grid.arrange(Art_GO_BP)
dev.off()

# Find the rows in the BP analysis where arterial DEGs are present
rows_with_gene_Col15a1 <- grep('Col15a1', Art_GO_BP_table$geneID)
Art_GO_BP_table[rows_with_gene_Col15a1, 2]

rows_with_gene_Cd34 <- grep('Cd34', Art_GO_BP_table$geneID)
Art_GO_BP_table[rows_with_gene_Cd34, 2]

rows_with_gene_S100a10 <- grep('S100a10', Art_GO_BP_table$geneID)
Art_GO_BP_table[rows_with_gene_S100a10, 2]



# Print the indices of the rows
print(rows_with_gene_Y)



par(mar=c(8,20,3,1))
barplot(Arterial_GO_resultsBP, showCategory = 30, cex.names = 0.3,
        title = 'Arterial biological processes')
dotplot(Arterial_GO_resultsBP, showCategory=30) +
  ggtitle("dotplot for Arterial ORA BP")

Arterial_GO_resultsMF <- enrichGO(gene = Arterial_DEG_genes,
                               universe = AllGenes_EC, OrgDb = org.Mm.eg.db,
                               keyType = "SYMBOL", ont = "MF",
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05)

par(mar=c(8,30,3,1))
barplot(Arterial_GO_resultsMF, showCategory = 30, cex.names = 0.3,
        title = 'Arterial molecular function')
dotplot(Arterial_GO_resultsMF, showCategory=30) +
  ggtitle("dotplot for Arterial ORA MF")


dev.off()



#------------------------------------------------------------------------------
# 19. GSEA ANALYSIS OF ARTERIAL CELLS
#------------------------------------------------------------------------------
## feature 1: numeric vector with Log2 FC
geneList = Arterial_DEG[,2]

## feature 2: named vector with gene names converted to ENTEZID
# Convert gene symbols to Entrez Gene IDs
entrez_ids <- mapIds(org.Mm.eg.db, keys = Arterial_DEG_genes,
                     keytype = "SYMBOL", column = "ENTREZID")

# Print the results
names(geneList) = as.character(Arterial_DEG_genes)

## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)



# get table with names and genes foldchange
names <- rownames(Arterial_DEG)
logfc <- Arterial_DEG$avg_log2FC
names(logfc) <- names


#get gene set to find hallmarks GO

m_df<- msigdbr(species = "Mus musculus", category = 'H')


fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#find enriched Hallmarks for our geneset
fgseaRes_Arterial <- fgsea(fgsea_sets, stats=logfc, minSize  = 2 )


fgseaResTidy_Arterial <- fgseaRes_Arterial %>%  as_tibble() %>%  arrange(desc(NES))
fgseaResTidy_Arterial %>%   dplyr::select(-leadingEdge, -ES) %>%   arrange(padj) %>%   head()


# Select top 10 upregulated and top 10 downregulated
a <- head(subset(fgseaResTidy_Arterial, fgseaResTidy_Arterial$NES > 0), 10)
b <- tail(subset(fgseaResTidy_Arterial, fgseaResTidy_Arterial$NES < 0), 10)
fgseaResTidy_ArterialNEW <- rbind(a,b)

pdf("GSEA_Arterial_Hallmarks.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

ggplot(fgseaResTidy_ArterialNEW, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES)) +
  coord_flip() +
  scale_fill_gradient2(low='red3', mid='yellow2', high='green3')+
  labs(x="Hallmark", y="Normalized Enrichment Score",
       title="Hallmark NES from GSEA_Arterial") + 
  theme_minimal()

dev.off()


#------------------------------------------------------------------------------
# 20. DEG ANALYSIS OF VENOUS CELLS
#------------------------------------------------------------------------------
Vein_DEG <- FindMarkers(Cont_Dll4_EC, ident.1 = 'Vein_Dll4-iDEC',
                            ident.2 = 'Vein_Control')
head(Vein_DEG, 20)
VTotal <- nrow(Vein_DEG)
# store data of number total DEG, Up significant and down significant
nrow(subset(Vein_DEG, Vein_DEG$p_val_adj < 0.05))
V_UpSig <- nrow(subset(Vein_DEG, Vein_DEG$p_val_adj < 0.05 & Vein_DEG$avg_log2FC > 0))
V_DwSig <- nrow(subset(Vein_DEG, Vein_DEG$p_val_adj < 0.05 & Vein_DEG$avg_log2FC < 0))

Cluster <- c(Cluster, 'Vein')
Total_DEG <- c(Total_DEG, VTotal)
UpSig_DEG <- c(UpSig_DEG, V_UpSig)
DwSig_DEG <- c(DwSig_DEG, V_DwSig)



#------------------------------------------------------------------------------
# 21. DEG ANALYSIS OF CAPILLARY CELLS
#------------------------------------------------------------------------------
Cap_DEG <- FindMarkers(Cont_Dll4_EC, ident.1 = 'Capillary_Dll4-iDEC',
                            ident.2 = 'Capillary_Control')
head(Cap_DEG, 20)
CTotal <- nrow(Cap_DEG)
# store data of number total DEG, Up significant and down significant
nrow(subset(Cap_DEG, Cap_DEG$p_val_adj < 0.05))
C_UpSig <- nrow(subset(Cap_DEG, Cap_DEG$p_val_adj < 0.05 & Cap_DEG$avg_log2FC > 0))
C_DwSig <- nrow(subset(Cap_DEG, Cap_DEG$p_val_adj < 0.05 & Cap_DEG$avg_log2FC < 0))

Cluster <- c(Cluster, 'Capillary')
Total_DEG <- c(Total_DEG, CTotal)
UpSig_DEG <- c(UpSig_DEG, C_UpSig)
DwSig_DEG <- c(DwSig_DEG, C_DwSig)

Cap_DEG_genes <- rownames(Cap_DEG)

write.csv(Cap_DEG, file='Capillary DEGs.csv')

#------------------------------------------------------------------------------
# 22. ORA ANALYSIS OF CAPILLARY CELLS
#------------------------------------------------------------------------------
pdf("ORA_Cap.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

Cap_GO_resultsCC <- enrichGO(gene = Cap_DEG_genes,
                                  universe = AllGenes_EC, OrgDb = org.Mm.eg.db,
                                  keyType = "SYMBOL", ont = "CC",
                                  pAdjustMethod = 'BH',
                                  pvalueCutoff = 0.05, qvalueCutoff = 0.05)

Cap_GO_CC_table <- Cap_GO_resultsCC@result

# Export the table to a PDF file
Cap_GO_CC <- tableGrob(Cap_GO_CC_table)
pdf("Capillary ORA CC table.pdf")
grid.arrange(Cap_GO_CC)
dev.off()

par(mar=c(8,20,3,1))
barplot(Angio_GO_resultsCC, showCategory = 30, cex.names = 0.3,
        title = 'Capillary cellular component')
dotplot(Angio_GO_resultsCC, showCategory=30) +
  ggtitle("dotplot for Capillary ORA CC")


Cap_GO_resultsBP <- enrichGO(gene = Cap_DEG_genes,
                                  universe = AllGenes_EC, OrgDb = org.Mm.eg.db,
                                  keyType = "SYMBOL", ont = "BP",
                                  pAdjustMethod = 'BH',
                                  pvalueCutoff = 0.05, qvalueCutoff = 0.05)

par(mar=c(8,20,3,1))
barplot(Cap_GO_resultsBP, showCategory = 30, cex.names = 0.3,
        title = 'Capillary biological processes')
dotplot(Cap_GO_resultsBP, showCategory=30) +
  ggtitle("dotplot for Capillary ORA BP")

Cap_GO_resultsMF <- enrichGO(gene = Cap_DEG_genes,
                                  universe = AllGenes_EC, OrgDb = org.Mm.eg.db,
                                  keyType = "SYMBOL", ont = "MF",
                                  pAdjustMethod = 'BH',
                                  pvalueCutoff = 0.05, qvalueCutoff = 0.05)

par(mar=c(8,30,3,1))
barplot(Cap_GO_resultsMF, showCategory = 30, cex.names = 0.3,
        title = 'Capillary molecular function')
dotplot(Cap_GO_resultsMF, showCategory=30) +
  ggtitle("dotplot for Capillary ORA MF")


dev.off()

#------------------------------------------------------------------------------
# 23. GSEA ANALYSIS OF CAPILLARY CELLS
#------------------------------------------------------------------------------

## feature 1: numeric vector with Log2 FC
geneList = Cap_DEG[,2]

## feature 2: named vector with gene names converted to ENTEZID
# Convert gene symbols to Entrez Gene IDs
entrez_ids <- mapIds(org.Mm.eg.db, keys = Cap_DEG_genes,
                     keytype = "SYMBOL", column = "ENTREZID")

# Print the results
names(geneList) = as.character(Cap_DEG_genes)

## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)



# get table wiht names and genes foldchange
names <- rownames(Cap_DEG)
logfc <- Cap_DEG$avg_log2FC
names(logfc) <- names


#get gene set to find hallmarks GO

m_df<- msigdbr(species = "Mus musculus", category = 'H')


fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#find enriched Hallmarks for our geneset
fgseaRes_Cap <- fgsea(fgsea_sets, stats=logfc, minSize  = 2 )


fgseaResTidy_Cap <- fgseaRes_Cap %>%  as_tibble() %>%  arrange(desc(NES))
fgseaResTidy_Cap %>%   dplyr::select(-leadingEdge, -ES) %>%   arrange(padj) %>%   head()


# Select top 10 upregulated and top 10 downregulated
a <- head(subset(fgseaResTidy_Cap, fgseaResTidy_Cap$NES > 0), 10)
b <- tail(subset(fgseaResTidy_Cap, fgseaResTidy_Cap$NES < 0), 10)
fgseaResTidy_CapNEW <- rbind(a,b)

pdf("GSEA_Capillary_Hallmarks.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

ggplot(fgseaResTidy_CapNEW, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES)) +
  coord_flip() +
  scale_fill_gradient2(low='red3', mid='yellow2', high='green3')+
  labs(x="Hallmark", y="Normalized Enrichment Score",
       title="Hallmark NES from GSEA_Capillary") + 
  theme_minimal()

dev.off()


#------------------------------------------------------------------------------
# 24. DEG ANALYSIS OF INTERFERON CELLS
#------------------------------------------------------------------------------
Interferon_DEG <- FindMarkers(Cont_Dll4_EC, ident.1 = 'Interferon_Dll4-iDEC',
                            ident.2 = 'Interferon_Control')
head(Interferon_DEG, 20)
nrow(Interferon_DEG)

# store data of number total DEG, Up significant and down significant
ITotal <- nrow(Interferon_DEG)
nrow(subset(Interferon_DEG, Interferon_DEG$p_val_adj < 0.05))
I_UpSig <- nrow(subset(Interferon_DEG, Interferon_DEG$p_val_adj < 0.05 & Interferon_DEG$avg_log2FC > 0))
I_DwSig <- nrow(subset(Interferon_DEG, Interferon_DEG$p_val_adj < 0.05 & Interferon_DEG$avg_log2FC < 0))

Cluster <- c(Cluster, 'Interferon')
Total_DEG <- c(Total_DEG, ITotal)
UpSig_DEG <- c(UpSig_DEG, I_UpSig)
DwSig_DEG <- c(DwSig_DEG, I_DwSig)




#------------------------------------------------------------------------------
# 25. MAKE DF WITH DEG FREQUENCIES PER CLUSTER
#------------------------------------------------------------------------------
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


saveRDS(Cont_Dll4_EC, file='C:\\Users\\scferreira\\Desktop\\EC only analysis\\EC_nichenet.rds')



#------------------------------------------------------------------------------
# 26. Identify arterial DEGS that encode for ligands using Nichenet's databases
#------------------------------------------------------------------------------
# load databases/networks
organism <- "mouse"

lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
  
# identify ligands that are DE in arteries by comparing with nichenets lr_network
DE_arterial_ligands <- intersect(Arterial_DEG_genes, lr_network$from)
Arterial_DEG_genes_sig <- rownames(subset(Arterial_DEG, Arterial_DEG$p_val_adj < 0.05))
DE_arterial_ligands_sig <- intersect(Arterial_DEG_genes_sig, lr_network$from)

# extract table of significantly DEG that are ligands
Arterial_DEG_ligands_index <- which(rownames(Arterial_DEG) %in% DE_arterial_ligands_sig)
Arterial_DEG_ligands <- Arterial_DEG[Arterial_DEG_ligands_index,]
# order according to log2FC values
Arterial_DEG_ligands <- Arterial_DEG_ligands[rev(order(Arterial_DEG_ligands$avg_log2FC)), ]

# Export the table to a PDF file
Art_l <- Arterial_DEG_ligands[,c(2,5)]
colnames(Art_l) <- c('Log2FC', 'p_val adj.')
Art_l$Log2FC <- round(Art_l$Log2FC, 3)
Art_l$'p_val adj.'<- format(Art_l$'p_val adj.', scientific = TRUE)
table_Arterial_DEG_ligands <- tableGrob(Art_l)

pdf("Arterial log2Fc table_new.pdf",
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4")
grid.arrange(table_Arterial_DEG_ligands)
dev.off()


#------------------------------------------------------------------------------
# 27. Identify capillary DEGS that encode for ligands using Nichenet's databases
#------------------------------------------------------------------------------
# identify ligands that are DE in capillaries by comparing with nichenets lr_network
DE_cap_ligands <- intersect(Cap_DEG_genes, lr_network$from)
Cap_DEG_genes_sig <- rownames(subset(Cap_DEG, Cap_DEG$p_val_adj < 0.05))
DE_cap_ligands_sig <- intersect(Cap_DEG_genes_sig, lr_network$from)

# extract table of significantly DEG that are ligands
Cap_DEG_ligands_index <- which(rownames(Cap_DEG) %in% DE_cap_ligands_sig)
Cap_DEG_ligands <- Cap_DEG[Cap_DEG_ligands_index,]
# order according to log2FC values
Cap_DEG_ligands <- Cap_DEG_ligands[rev(order(Cap_DEG_ligands$avg_log2FC)), ]

# Export the table to a PDF file
Cap_l <- Cap_DEG_ligands[,c(2,5)]
colnames(Cap_l) <- c('Log2FC', 'p_val adj.')
Cap_l$Log2FC <- round(Cap_l$Log2FC, 3)
Cap_l$'p_val adj.'<- format(Cap_l$'p_val adj.', scientific = TRUE)
table_Cap_DEG_ligands <- tableGrob(Cap_l)

pdf("Capillary log2Fc table_new.pdf",
    width = 11.69, height = 30)   # Width and height in inches
    
grid.arrange(table_Cap_DEG_ligands)
dev.off()



#------------------------------------------------------------------------------
# 28. dotplot of DE ligands according to clustering
#------------------------------------------------------------------------------
# draw dotplot of ligands according to clustering
Idents(Cont_Dll4_EC) <- "SR_res_0.6 EC_cos"

pdf("DotPlot_arterial ligands.pdf",         # File name
    width = 11.69, height = 16,   # Width and height in inches
    paper = "A4")          # Paper size                              
DotPlot(Cont_Dll4_EC, features = rownames(Art_l),group.by = 'Celltype_genotype') +
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()
dev.off()


pdf("DotPlot_capillary ligands.pdf",         # File name
    width = 11.69, height = 16,   # Width and height in inches
    paper = "A4")          # Paper size                              
DotPlot(Cont_Dll4_EC, features = rownames(Cap_l),group.by = 'Celltype_genotype') +
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()
dev.off()



#------------------------------------------------------------------------------
# 29. Extract corresponding receptors of DE ligands
#------------------------------------------------------------------------------
# create a vector with all unique artery and capillary DEG ligands
ligands <- c(DE_cap_ligands_sig, DE_arterial_ligands_sig)
ligands <- unique(ligands)

# create list of receptors of the DE ligands
receptors <- c()
for (i in 1:length(ligands)){
    ind <- which(lr_network$from == ligands[i])
    receptors <- c(receptors, lr_network$to[ind])
}
receptors <- unique(receptors)


#create table with mouse ligands and receptors
receptors_per_ligand <- c()
for (i in 1:length(ligands)){
  ind <- which(lr_network$from == ligands[i])
  if (length(ind > 1)){
    receptors_per_ligand[i] <- paste(c(lr_network$to[ind]), collapse = ",")
  } else {
    receptors_per_ligand[i] <- lr_network$to[ind]
  }
}

DEG_ligand_receptors <- cbind(ligands, receptors_per_ligand)
colnames(DEG_ligand_receptors) <- c('Mouse ligands', 'Mouse receptors')
DEG_ligand_receptors_df <- as.data.frame(DEG_ligand_receptors)


DEG_LR <- tableGrob(DEG_ligand_receptors_df)

pdf("Table with mouse DEG ligands and corresponding receptors.pdf",
    width = 42, height = 60,   # Width and height in inches
    paper = "A4")
grid.arrange(DEG_LR)
dev.off()

write.csv(DEG_ligand_receptors_df,
          file='Table with mouse DEG ligands and corresponding receptors.csv')



#------------------------------------------------------------------------------
# 30. Convert mouse gene symbols to human gene symbols to have human gene names
# of DE ligand receptors for downstream analysis of Kuppe et al 2020 data
#------------------------------------------------------------------------------
Species_conversion <- function(mouseids, horg, morg, orth){
  mouseg <- mapIds(morg, mouseids, "ENTREZID", "SYMBOL")
  mapped <- select(orth, mouseg, "Homo_sapiens","Mus_musculus")
  names(mapped) <- c("Mus_egid","Homo_egid")
  husymb <- select(horg, as.character(mapped[,2]), "SYMBOL","ENTREZID")
  return(data.frame(Mus_symbol = mouseids,
                    mapped,
                    Homo_symbol = husymb[,2]))
}

human_mouse_ortho <- Species_conversion(receptors, org.Hs.eg.db, org.Mm.eg.db,
                                        Orthology.eg.db)
receptors_human <- human_mouse_ortho$Homo_symbol

