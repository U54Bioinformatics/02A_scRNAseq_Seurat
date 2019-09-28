
#Load required libraries
library(Seurat)
library(data.table)

#Counts
counts.df <- fread("X965X9.noaggr.counts.txt")
counts.df <- subset(counts.df, counts.df$'Gene Symbol' != "")
counts.df <- unique(counts.df, by='Gene Symbol')
counts.df <- as.data.frame(counts.df)
rownames(counts.df) <- counts.df$`Gene Symbol`
counts <- counts.df[,4:ncol(counts.df)]
rm(counts.df)

save(counts, file="counts.RData")

#Annotation
Annot.df <- fread("X965X9.metadata/cell_metadata.txt")
Annot.df <- Annot.df[match(colnames(counts), Annot.df$Cell), ]
rownames(Annot.df) <- Annot.df$Cell

save(Annot.df, file="Annot.df.RData")

table(Annot.df$Sample)

ggplot(as.data.frame(table(Annot.df$Sample)), aes(x=Var1, y=Freq)) + geom_bar(stat="identity", fill = "#1f1d66") + 
  theme_classic(base_size = 18) + xlab("Sample") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Seurat pipeline 
## Create Seurat object, perform QC analysis and generate plots 
seu.coh051 <- CreateSeuratObject(counts,  project = "COH051",  meta.data = Annot.df,  min.cells = 5, min.features = 100)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seu.coh051[["percent.mt"]] <- PercentageFeatureSet(seu.coh051, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(seu.coh051, features = c("nFeature_RNA"), group.by="Sample")
VlnPlot(seu.coh051, features = c("nCount_RNA"), group.by="Sample")
VlnPlot(seu.coh051, features = c("percent.mt"), group.by="Sample")

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
FeatureScatter(seu.coh051, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by="Sample")
FeatureScatter(seu.coh051, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="Sample")

## Filter based on features, count and percent MT
seu.coh051 <- subset(seu.coh051, subset = nFeature_RNA < 5000 & nFeature_RNA > 100 & nCount_RNA < 2e4 & percent.mt < 25)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
FeatureScatter(seu.coh051, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by="Sample")
FeatureScatter(seu.coh051, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="Sample")

## Mark cell cycle genes
#Get cell cycle genes
cc.genes <- readLines("regev_lab_cell_cycle_genes.txt")

# Segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
#Create scores
seu.coh051 <- CellCycleScoring(seu.coh051, s.genes, g2m.genes, set.ident = TRUE)

# Visualize QC metrics as a violin plot
VlnPlot(seu.coh051, features = c("nFeature_RNA"), group.by="Sample")
VlnPlot(seu.coh051, features = c("nCount_RNA"), group.by="Sample")
VlnPlot(seu.coh051, features = c("percent.mt"), group.by="Sample")

## NormalizeData and ScaleData
#normalize
seu.coh051 <- NormalizeData(seu.coh051, normalization.method = "LogNormalize", scale.factor = 10000)

#identify variable features
seu.coh051 <- FindVariableFeatures(seu.coh051, selection.method = "vst", nfeatures=1000)
top10 <- head(VariableFeatures(seu.coh051), 25) #top 10 variable features

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seu.coh051)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#scaling
seu.coh051 <- ScaleData(seu.coh051)

## PCA analysis and plots
seu.coh051 <- RunPCA(seu.coh051, features = VariableFeatures(object = seu.coh051))
VizDimLoadings(seu.coh051, dims = 1:2, reduction = "pca")
DimPlot(seu.coh051, reduction = "pca", group.by = "Sample")
DimHeatmap(seu.coh051, dims = 1:10, cells = 500, balanced = TRUE)
ElbowPlot(seu.coh051)

## Cluster data, run UMAP and TSNE
seu.coh051 <- FindNeighbors(seu.coh051, dims = 1:15)
seu.coh051 <- FindClusters(seu.coh051, resolution = 0.7)

seu.coh051 <- RunUMAP(seu.coh051, dims = 1:15)
DimPlot(seu.coh051, reduction = "umap", group.by = "Sample", label=FALSE, pt.size = 0.5) + 
  scale_color_brewer(palette="Set2") +
  theme(legend.text = element_text(size=32)) +
  theme(text = element_text(size=24, face="bold")) +
  theme(axis.line = element_line(colour = 'black', size = 2)) +
  theme(axis.text = element_text(size=24))

## Plot cell types
x <- match(colnames(seu.coh051), Annot.df.singler$Cell)

seu.coh051[["HPCA"]] <- Annot.df.singler$SingleR_HPCA[x]
seu.coh051[["ENCODE"]] <- Annot.df.singler$SingleR_ENCODE[x]

DimPlot(seu.coh051, reduction = "umap", group.by = "ENCODE", shape.by = "Sample", label=TRUE, label.size=6, repel=TRUE, pt.size = 0.75) + 
  theme(legend.text = element_text(size=18)) +
  theme(text = element_text(size=12, face="bold")) +
  theme(axis.line = element_line(colour = 'black', size = 1.5)) +
  theme(axis.text = element_text(size=18))

#plot features on top of cluster and cells
FeaturePlot(object = seu.coh051, features = c("ReadsByFeatures"), pt.size = 0.25)
FeaturePlot(object = seu.coh051, features = c("nCount_RNA"), pt.size = 0.25)

seu.coh051[["ReadsByFeatures"]] <- seu.coh051$`Total Reads` / seu.coh051$`Expressed Features`

############## Individual markers #############  
#Macrophages
FeaturePlot(seu.coh051, features = c("CD14", "AIF1", "FCER1G", "FCGR3A", "TYROBP", "CSF1R"), order=TRUE, pt.size=0.25)

#T-cells
FeaturePlot(seu.coh051, features = c("CD2", "CD3D", "CD3E", "CD3G"), order=TRUE, pt.size=0.25)

#Oligodendrocytes
FeaturePlot(seu.coh051, features = c("MBP", "TF", "PLP1", "MAG", "MOG", "CLDN11"), pt.size=0.25)

#Fibroblasts 
FeaturePlot(seu.coh051, features = c("PDPN", "THY1", "CD34", "CDH11"), order=TRUE, pt.size=0.25)

#Epithelial cell
FeaturePlot(seu.coh051, features = c("EPCAM", "ESR1", "KRT7", "KRT8", "KRT18", "KRT19"), order=TRUE, pt.size=0.25)

#EMT
FeaturePlot(seu.coh051, features = c("CDH1", "DSP", "CTNNB1", "ZEB1", "ZEB2", "SNAI1", "TWIST1"), order=TRUE, pt.size=0.25)

#Oncogenes
FeaturePlot(seu.coh051, features = c("MET", "EGFR", "BRAF", "SMO", "CDK6"), order=TRUE, pt.size=0.25)

#Normal vs Cancer
FeaturePlot(seu.coh051, features = c("PTPRC", "MKI67", "PCNA"), order=TRUE, pt.size=0.25)

#COSMIC amplified
FeaturePlot(seu.coh051, features = c("AKT2", "ALK", "CCNE1", "DROSHA", "EGFR", "ERBB2"), order=TRUE, pt.size=1)
FeaturePlot(seu.coh051, features = c("ERG", "FLT4", "JUN", "KAT6A", "LMO1", "MDM2"), order=TRUE, pt.size=1)
FeaturePlot(seu.coh051, features = c("MDM4", "MITF", "MYC", "MYCL", "MYCN", "NKX2-1"), order=TRUE, pt.size=1)
FeaturePlot(seu.coh051, features = c("NSD3", "NTRK1", "PPM1D", "RAF1", "REL", "SOX2"), order=TRUE, pt.size=1)
