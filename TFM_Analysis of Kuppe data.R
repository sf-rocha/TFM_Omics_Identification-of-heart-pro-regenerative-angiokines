############################################################################
## ANALYSIS OF KUPEE at al DATA TO CHECK FOR EXPRESSION OF HUMAN
## ORTHOLOGUES OF THE RECEPTORS OF THE DE LIGANDS IN DLL4-IDEC SAMPLES
############################################################################

#-----------------------------------------------------------------------------
#Load and read .rds files
#-----------------------------------------------------------------------------
file4 <- choose.files()
Kuppe <- readRDS(file4)
Kuppe


#-----------------------------------------------------------------------------
# Kuppe data has genes with Ensembl ID and no gene symbols
#-----------------------------------------------------------------------------
library(org.Hs.eg.db)
library(AnnotationDbi)
library(Seurat)
library(biomaRt)
library(RColorBrewer)
library(dittoSeq)
library(ggplot2)
library(pheatmap)


#-----------------------------------------------------------------------------
# get Ensemlb IDs for human receptor list (receptors_human)
#-----------------------------------------------------------------------------
receptors_ensembl <- mapIds(org.Hs.eg.db, keys = receptors_human,
                            keytype = "SYMBOL", column="ENSEMBL")


# convert the dictionary int  a normal list/vector
receptors_ensembl_n <- unname(receptors_ensembl)
receptors_ensembl_n <- unlist(receptors_ensembl_n)

Human_recep <- cbind(receptors_human, receptors_ensembl)
Human_recep <- as.data.frame(Human_recep)
rownames(Human_recep) <- NULL



#-----------------------------------------------------------------------------
# plot all identified receptor according to cardiac region
#-----------------------------------------------------------------------------

Idents(Kuppe) <- 'major_labl'   # cardiac region

DotPlot(Kuppe, features = c(receptors_ensembl_n[1:50])) + 
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()

DotPlot(Kuppe, features = c(receptors_ensembl_n[51:100])) + 
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()

DotPlot(Kuppe, features = c(receptors_ensembl_n[101:153])) + 
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()



#-----------------------------------------------------------------------------
# Calculate Log2FC of receptors to filetr for those that are more expressed in 
# ischaemic and boundary zones than in healthy patienets or remote regions
#-----------------------------------------------------------------------------
object<-Kuppe
object<-SetIdent(object, value="major_labl")
top_markers <- receptors_ensembl_n

#Create a matrix with the average expression of each gene and each condition
assay<-"RNA"
Avexpre<-AverageExpression(object, features = top_markers, group.by = "ident",
                           assays = assay, return.seurat = T)
x<-Avexpre@assays[[assay]]@data

# Reorder the x columns as you want the conditions to be plotted
desired_order <- c("CTRL", "RZ", "BZ", "IZ", "FZ")
x <- x[, desired_order]  # Rearrange columns according to desired order


#Transform the matrix to log2FC values
# Assuming 'x' is your expression matrix and 'Control' is the condition
control_column <- "CTRL"  # Replace with the actual column name if needed

# Extract the 'Control' column and replicate it to match the dimensions of
# the matrix
control_values <- x[, control_column]
control_matrix <- matrix(rep(control_values, ncol(x)),
                         ncol = ncol(x), byrow = FALSE)

# Divide each value by the 'Control' column and take log2 transformation
x_log2fc <- log2(x / control_matrix)

#Remove rows with a value=NaN or =+-Inf
x_log2fc <- na.omit(x_log2fc)

rows_to_remove <- which(apply(x_log2fc, 1, function(row) any(row == -Inf)))

if (length(rows_to_remove) > 0) {
  # Remove rows with -Inf values
  x_log2fc <- x_log2fc[-rows_to_remove, ]
}


# filter matrix for genes that are up in BZ or IZ in relation to control but
# not in BZ
x_log2fc_up <- x_log2fc[x_log2fc[, "RZ"] <= 0, ]
x_log2fc_up <- x_log2fc_up[x_log2fc_up[, "BZ"] > 0 | x_log2fc_up[, "IZ"] > 0, ]


# convert rownames from ensembl to symbol
ensembl_up <- row.names(x_log2fc_up)
gene_id_up <- c()
for (i in 1:nrow(x_log2fc_up)){
  ind <- which(Human_recep$receptors_ensembl == ensembl_up[i])
  gene_id_up <- c(gene_id_up, Human_recep$receptors_human[ind])
}

x_log2fc_up_ENS <- x_log2fc_up
row.names(x_log2fc_up) <- gene_id_up

#-----------------------------------------------------------------------------
#Do the heatmap with log2FC of receptors by region
#-----------------------------------------------------------------------------

# Define the breaks corresponding to the desired color transitions
breaks <- seq(-3, 3, length.out = 101)  # 100 breaks from -1.5 to 1.5 inclusive
breaks[51] <- 0  # Set the 51st break (at 0) to ensure it's completely white

# Create a custom color palette using colorRamp and the breaks
my_palette <- colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1)

pdf("Heatmap of upregulated recpetor in BZ or IZ_5.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size  
# Create the heatmap using the custom color scale and breaks
pheatmap(x_log2fc_up,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = FALSE,
         color = my_palette,
         breaks = breaks,
         main = "Log2-Fold Changes Compared to Control",
         fontsize = 6,
         cellnote_col = "red", 
         annotation_legend = TRUE,
         annotation_names_col = TRUE
)
dev.off()

Idents(Kuppe) <- 'major_labl'



#-----------------------------------------------------------------------------
# Do dotplot with log2FC of receptors by region
#-----------------------------------------------------------------------------
# Change order of levels in major_labl metadata column
Kuppe@meta.data$major_labl <- factor(Kuppe@meta.data$major_labl,
                                     levels = desired_order)
Idents(Kuppe) <- 'major_labl'


pdf("DotPlot_Human receptors",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size
Idents(Kuppe) <- 'major_labl'
DotPlot(Kuppe, features = rownames(x_log2fc_up_ENS)) + 
  scale_colour_gradient2(low="#e1e1e1", mid="#fe9929", high="#7e2c06") +
  coord_flip()
dev.off()




#-----------------------------------------------------------------------------
# Do violin plot with expression of receptors by cell type
# (normalized cell count)
#-----------------------------------------------------------------------------
#check min cell number by cell type
dwnsamp <- min(table(Kuppe@meta.data$cell_type_original))
Idents(Kuppe) <- 'cell_type_original'
Kuppe_dwns_celltype <- subset(Kuppe, downsample=dwnsamp)
# confirm downsampling
table(Kuppe_dwns_celltype@meta.data$cell_type_original)

pdf("ViolinPlot_Human receptors",         # File name
    width = 21, height = 30,   # Width and height in inches
    paper = "A4r")          # Paper size
Idents(Kuppe_dwns_celltype) <- 'cell_type_original'
VlnPlot(Kuppe_dwns_celltype, features <- rownames(x_log2fc_up_ENS),
        stack = TRUE, flip = TRUE)
dev.off()



#-----------------------------------------------------------------------------
# Do violin plot with expression of receptors by zone
# (normalized cell count)
#-----------------------------------------------------------------------------
#check min cell number by zone
dwnsampz <- min(table(Kuppe@meta.data$major_labl))
Idents(Kuppe) <- 'major_labl'
Kuppe_dwns_zone <- subset(Kuppe, downsample=dwnsampz)
# confirm downsampling
table(Kuppe_dwns_zone@meta.data$major_labl)

pdf("ViolinPlot_Human receptors by zone.pdf",         # File name
    width = 21, height = 30,   # Width and height in inches
    paper = "A4r")          # Paper size
Idents(Kuppe_dwns_zone) <- 'major_labl'
VlnPlot(Kuppe_dwns_zone, features <- rownames(x_log2fc_up_ENS),
        stack = TRUE, flip = TRUE)
dev.off()




#-----------------------------------------------------------------------------
# Do barplot with frequency of cel types by zone
# (normalized cell count)
#-----------------------------------------------------------------------------
pdf("Barplot with human cell types.pdf",         # File name
    width = 11.69, height = 8.27,   # Width and height in inches
    paper = "A4r")          # Paper size

dittoBarPlot(Kuppe, var = "cell_type_original", group.by = "major_labl",
             color.panel = c('#e6cc00', '#CE299C', '#f2981b', '#1FC116',
                                      'red', '#1b5bf2', '#45b05c', 'purple',
                                      '#0ACDED', 'yellow', 'grey'))
dev.off()




zone_count <- table(Kuppe@meta.data$major_labl)
