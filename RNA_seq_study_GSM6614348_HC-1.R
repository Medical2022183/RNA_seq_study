library(dplyr)
library(Seurat)

setwd('/Users/yanyucheng/Desktop/RNA_seq/')

## load data
# method one
dir.10x = './data/filtered_gene_bc_matrices/GSM6614348_HC-1/'
genes = read.table(paste0(dir.10x,'features.tsv'), stringsAsFactors = F, header = F)$V2
genes = make.unique(genes, sep = '.')
barcodes = readLines(paste0(dir.10x,'barcodes.tsv'))
mtx = Matrix::readMM(paste0(dir.10x,'matrix.mtx'))
mtx = as(mtx,"CsparseMatrix")
colnames(mtx) = barcodes
rownames(mtx) = genes
pbmc = CreateSeuratObject(counts = mtx, project = "pbmc3k",
                          min.cells = 3, min.features = 200)

# method two
dir.10x = './data/filtered_gene_bc_matrices/GSM6614348_HC-1/'
pbmc.data <- Read10X(data.dir = dir.10x)
pbmc = CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",
                          min.cells = 3, min.features = 200)

# QC
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Feature selection
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc),10)
top10
#plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T, xnudge = 0, ynudge = 0)
plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimPlot(pbmc, reduction = "pca")
ElbowPlot(pbmc,ndims = 50)

pbmc <- RunUMAP(pbmc, dims = 1:10) # umap tsne
FeaturePlot(pbmc, features = c('FCGR3A', 'CD14'), reduction = 'umap')

# Cluster
pbmc <- FindNeighbors(pbmc, dims = 1:10) # louvain cluster, graph based
pbmc <- FindClusters(pbmc, resolution = 0.5) # resolution 越大群越多
DimPlot(pbmc, reduction = "umap", group.by = 'seurat_clusters', label = T)
#FeaturePlot(pbmc, features = c("MS4A1", "TYROBP", "CD14", 'FCGRA', "FCER1A",
#                                "CCR7", "IL7R", "PPBP", "CD8А"))
FeaturePlot(pbmc, features = c("IGLL5","IGLC2","IGLC3","ADAMDEC1","CXCL14","IGHM",    
                               "IGHG2","IGHG1","CFD","CCL13"))
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = T)

# DE analysis
pbmc.markers <- FindAllMarkers(pbmc, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)

CD4.mem.DEGS <- FindMarkers(pbmc, ident.1 = 'Memory CD4 T', ident.2 = 'Naive CD4 T', min.pct = 0.25)

# gene signature analysis
exhaustion_genes = list(c('PDCD1', 'CD160', 'FASLG', 'CD244', 'LAG3', 'TNFRSF1B', 'CCR5', 'CCL3',
                          'MX1', 'IRF4', 'EOMES', 'PBX3', 'NFATC1'))
pbmc = Seurat::AddModuleScore(pbmc, features = exhaustion_genes, name='exhaustion.score')
FeaturePlot(pbmc, features = 'exhaustion.score1', reduction = 'umap')
VlnPlot(pbmc, features = 'exhaustion.score1')
