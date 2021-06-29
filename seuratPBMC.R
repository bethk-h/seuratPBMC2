# 10X Genomics PBMC data set. 2700 cells sequenced by Illumina nextseq500
install.packages('Seurat')
yes
library(dplyr)
library(Seurat)
library(patchwork)
#Load PBMC dataset
pbmc.data <- Read10X(data.dir = "~/filtered_gene_bc_matrices/hg19")
# Initialise the seurat object with raw data
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells=3, min.features=200)                     
pbmc  
#Exploring count matrix genes in first 30 cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
#The number of unique genes, total molecules and percentage mitochondrial counts can allow cell filtering.
#QC metrics for the first 5 cells
head (pbmc@meta.data, 5)
#Note mitochondrial data is missing
#Visualising QC metrics as a violin plot
VlnPlot(pbmc, features =c("nFeature_RNA", "nCount_RNA", "percent.mt", ncol=3))
#To visualise feature feature relationships use FeatureScatter
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2
#Filtering cells with unique feature counts over 2,500 (>2500 genes suggests doublets) or <200
pbmc <- subset(pbmc, subset = nFeature_RNA >200 & nFeature_RNA <2500)
# Normalising data - normalises the feature expression measurements for each cell by total expression, multiplies by a scale factor and log-transforms
#pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10,000)
pbmc <- NormalizeData(pbmc)
#Selecting highly variable features - high cell to cell variation
pbmc <- FindVariableFeatures(pbmc, selection.method ="vst", nfeatures=2000)
#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
#Plot variable features without labels
plot1 <- VariableFeaturePlot(pbmc)
#Plot variable features with labels
plot2 <- LabelPoints(plot=plot1, points = top10, repel=TRUE)
plot1
plot2
#Scaling the data before dimensionality reduction
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#Performing a PCA
pbmc <- RunPCA(pbmc, features=VariableFeatures(object = pbmc))
DimPlot(pbmc, reduction = "pca")
#To determine which PCs to use down the line can use heatmap. Setting 'cells' to a number plots extreme cells on both ends.
DimHeatmap(pbmc, dims=1, cells=500, balanced=TRUE)
DimHeatmap(pbmc, dims=1:15, cells=500, balanced=TRUE)
#To determine dimensionality you can use approx. methods such as elbow plot, or JackStraw method.Note: JackStraw is very slow.
pbmc <- JackStraw(pbmc, num.replicate=100)
pbmc <- ScoreJackStraw(pbmc, dims=1:20)
ElbowPlot(pbmc)
#Use K-nearest neighbours based on PCA distance in FindNeighbours() function based on previously definited dimensions (10) from elbow plot
pbmc <- FindNeighbors(pbmc, dims =1:10)
#group cells using findclusters function
pbmc <- FindClusters(pbmc, resolution = 0.5)
#Look at cluster ID's of first 5 cells
head(Idents(pbmc), 5)
#Running a UMAP
pbmc <- RunUMAP(pbmc, dims=1:10)
DimPlot(pbmc, reduction = "umap")
#Finding DEGs/cluster biomarkers
#Find all markers of cluster 2. min.pct requires a feature to be detected at a minimum percentage in either 2 groups of cells
cluster2.markers <- FindMarkers(pbmc, ident.1 =2, min.pct =0.25)
head(cluster2.markers, n=5)
#Find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1=5, ident.2 =c(0,3), min.pct=0.25)
head(cluster5.markers, n=5)
#Find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
pbmc.markers %>%
  group_by(cluster)%>%
  top_n(n=2, wt = avg_log2FC)
#VlnPLot shows expression probability distributions across clusters
VlnPlot(pbmc, features =c("MS4A1", "CD79A"))
#This can also be plotted as raw counts
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot="counts", log = TRUE)
#FeaturePLot visualises feature expression on tSNE or PCA
FeaturePlot(pbmc, features=c("MS4A1", "GNLY", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
#Assigning cell type identity to clusters using canonical markers
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction ="umap", label = TRUE, pt.size=0.5) + NoLegend()
