library(reshape2)
library(Seurat)
library(SeuratData)
library(cowplot)
library(patchwork)
library(ggplot2)
library(dplyr)
library("BiocFileCache")
require(gridExtra)
library(EnhancedVolcano)
library(viridis)
library(RColorBrewer)
library("writexl")

IVD3=read.csv("IVD3.csv", sep=",", header=T)
IVD3m= as.matrix(IVD3[, 3:2495])
row.names(IVD3m)=IVD3$HGNC.symbol
IVD3m=as.data.frame(IVD3m)
IVD3m=as.matrix(IVD3m)
nIVD1 <- CreateSeuratObject(counts = IVD3m, project = "IVD3m", min.cells = 3, min.features = 200)

IVD4=read.csv("IVD4.csv", sep=",", header=T)
IVD4m= as.matrix(IVD4[, 3:2059])
row.names(IVD4m)=IVD4$HGNC.symbol
IVD4m=as.data.frame(IVD4m)
IVD4m=as.matrix(IVD4m)
nIVD2 <- CreateSeuratObject(counts = IVD4m, project = "IVD4m", min.cells = 3, min.features = 200)

IVD5=read.csv("IVD5.csv", sep=",", header=T)
IVD5m= as.matrix(IVD5[, 3:3197])
row.names(IVD5m)=IVD5$HGNC.symbol
IVD5m=as.data.frame(IVD5m)
IVD5m=as.matrix(IVD5m)
nIVD3 <- CreateSeuratObject(counts = IVD5m, project = "IVD5m", min.cells = 3, min.features = 200)

aIVD=read.csv("NPCX.csv", sep=",", header=T)
aIVDm= as.matrix(aIVD[, 3:2702])
row.names(aIVDm)=aIVD$HGNC.symbol
aIVDm=as.data.frame(aIVDm)
aIVDm=as.matrix(aIVDm)
aIVD <- CreateSeuratObject(counts = aIVDm, project = "aIVDm", min.cells = 3, min.features = 200)

DataAll <- c(nIVD1, nIVD2, nIVD3, aIVD)

for (i in 1:length(DataAll)) {
  DataAll[[i]] <- NormalizeData(DataAll[[i]], verbose = FALSE)
  DataAll[[i]] <- FindVariableFeatures(DataAll[[i]], selection.method = "vst", 
                                       nfeatures = 6000, verbose = FALSE)
}

IVD.anchors <- FindIntegrationAnchors(object.list = DataAll, dims = 1:30)
IVD.integrated <- IntegrateData(anchorset = IVD.anchors, dims = 1:30)
DefaultAssay(IVD.integrated) <- "integrated"

IVD.integrated <- ScaleData(IVD.integrated, verbose = FALSE)
IVD.integrated <- RunPCA(IVD.integrated, npcs = 50, verbose = FALSE)
IVD.integrated <- RunUMAP(IVD.integrated, reduction = "pca", dims = 1:10)
DimPlot(IVD.integrated, reduction = "umap", pt.size = 2, label.size = 20) + theme(aspect.ratio = "1")

Idents(object = IVD.integrated)
levels(x = IVD.integrated)
IVD.integrated2 <- RenameIdents(object = IVD.integrated, "IVD3m" = "nIVD", "IVD4m" = "nIVD", 'IVD5m' = 'nIVD', 'aIVDm' = 'aIVD')
levels(IVD.integrated2)

Plot1a <- DimPlot(IVD.integrated2, reduction = "umap", pt.size = 2, label = FALSE, cols = c("brown", "dodgerblue4")) + theme(aspect.ratio = "1")
Plot1a

IVD.integrated3 <- FindNeighbors(IVD.integrated2, dims = 1:10)
IVD.integrated3 <- FindClusters(IVD.integrated3, resolution = 0.5)
IVD.integrated3 <- RunUMAP(IVD.integrated3, reduction = "pca", dims = 1:10)

IVD.integrated3$ClusterNumber <- Idents(IVD.integrated3)
Idents(IVD.integrated3) <- "orig.ident"
IVD.integrated3 <- RenameIdents(object = IVD.integrated3, "IVD3m" = "nIVD", "IVD4m" = "nIVD", 'IVD5m' = 'nIVD', 'aIVDm' = 'aIVD')
IVD.integrated3$MajorType <- Idents(IVD.integrated3)
table(IVD.integrated3$MajorType)

IVD.nIVD <- subset(x = IVD.integrated3, idents = "nIVD")
IVD.aIVD <- subset(x = IVD.integrated3, idents = "aIVD")

Idents(IVD.nIVD) <- "ClusterNumber"
Idents(IVD.aIVD) <- "ClusterNumber"
Idents(IVD.integrated3) <- "ClusterNumber"

IVD.integrated3 <- RenameIdents(object = IVD.integrated3, "0" = "NP2", "1" = "non-MC", "2" = "OAF2", "3" = "NP3", "4" = "OAF1", "5" = "NP4", "6" =  "NP1", "7" = "IAF1", "8" = "IAF2", "9" = "NC1", "10" = "OAF3", "11" = "NC2", "12" = "RBC", "13" =  "IC")
IVD.nIVD <- RenameIdents(object = IVD.nIVD, "0" = "NP2", "1" = "non-MC", "2" = "OAF2", "3" = "NP3", "4" = "OAF1", "5" = "NP4", "6" =  "NP1", "7" = "IAF1", "8" = "IAF2", "9" = "NC1", "10" = "OAF3", "11" = "NC2", "12" = "RBC", "13" =  "IC")
IVD.aIVD <- RenameIdents(object = IVD.aIVD, "0" = "NP2", "1" = "non-MC", "2" = "OAF2", "3" = "NP3", "4" = "OAF1", "5" = "NP4", "6" =  "NP1", "7" = "IAF1", "8" = "IAF2", "9" = "NC1", "10" = "OAF3", "11" = "NC2", "12" = "RBC", "13" =  "IC")

IVD.integrated3$Subtype <- Idents(IVD.integrated3)
IVD.nIVD$Subtype <- Idents(IVD.nIVD)
IVD.aIVD$Subtype <- Idents(IVD.aIVD)

Idents(IVD.integrated3) <- IVD.integrated3$Subtype
Idents(IVD.nIVD) <- IVD.nIVD$Subtype
Idents(IVD.aIVD) <- IVD.aIVD$Subtype

levels(x = IVD.integrated3) <- c('RBC', 'IC', "non-MC", "OAF3", 'OAF2', "OAF1", "IAF1", "IAF2", "NP4", "NP3", "NP2", "NP1", "NC2", "NC1")
levels(x = IVD.aIVD) <- c('RBC', 'IC', "non-MC", "OAF3", 'OAF2', "OAF1", "IAF1", "IAF2", "NP4", "NP3", "NP2", "NP1", "NC2", "NC1")
levels(x = IVD.nIVD) <- c('RBC', 'IC', "non-MC", "OAF3", 'OAF2', "OAF1", "IAF1", "IAF2", "NP4", "NP3", "NP2", "NP1", "NC2", "NC1")

DimPlot(IVD.integrated3, reduction = "umap", pt.size = 0.3, label.size = FALSE) + theme(aspect.ratio = "1")

Idents(IVD.integrated3) <- "ClusterNumber"
DimPlot(IVD.integrated3, reduction = "umap", pt.size = 0.3, label.size = FALSE, split.by = "MajorType") + theme(aspect.ratio = "1")

Idents(IVD.nIVD) <- "orig.ident"
DimPlot(IVD.nIVD, reduction = "umap", pt.size = 0.6, label.size = FALSE, cols = c("red", "limegreen", "lightslateblue")) + theme(aspect.ratio = "1")
Idents(IVD.nIVD) <- "Subtype"

IVD.All.RBC.markers <- FindMarkers(IVD.integrated3, ident.1 = "RBC", min.pct = 0.25)
IVD.All.nonMC.markers <- FindMarkers(IVD.integrated3, ident.1 = "non-MC", min.pct = 0.25)
IVD.All.OAF3.markers <- FindMarkers(IVD.integrated3, ident.1 = "OAF3", min.pct = 0.25)
IVD.All.OAF2.markers <- FindMarkers(IVD.integrated3, ident.1 = "OAF2", min.pct = 0.25)
IVD.All.OAF1.markers <- FindMarkers(IVD.integrated3, ident.1 = "OAF1", min.pct = 0.25)
IVD.All.IAF2.markers <- FindMarkers(IVD.integrated3, ident.1 = "IAF2", min.pct = 0.25)
IVD.All.IAF1.markers <- FindMarkers(IVD.integrated3, ident.1 = "IAF1", min.pct = 0.25)
IVD.All.NP4.markers <- FindMarkers(IVD.integrated3, ident.1 = "NP4", min.pct = 0.25)
IVD.All.NP3.markers <- FindMarkers(IVD.integrated3, ident.1 = "NP3", min.pct = 0.25)
IVD.All.NP2.markers <- FindMarkers(IVD.integrated3, ident.1 = "NP2", min.pct = 0.25)
IVD.All.NP1.markers <- FindMarkers(IVD.integrated3, ident.1 = "NP1", min.pct = 0.25)
IVD.All.IC.markers <- FindMarkers(IVD.integrated3, ident.1 = "IC", min.pct = 0.25)
IVD.All.NC2.markers <- FindMarkers(IVD.integrated3, ident.1 = "NC2", min.pct = 0.25)
IVD.All.NC1.markers <- FindMarkers(IVD.integrated3, ident.1 = "NC1", min.pct = 0.25)

IVD.All.RBC.top10 <- rownames(head(IVD.All.RBC.markers, n = 10))
IVD.All.nonMC.top10 <- rownames(head(IVD.All.nonMC.markers, n = 10))
IVD.All.OAF3.top10 <- rownames(head(IVD.All.OAF3.markers, n = 10))
IVD.All.OAF2.top10 <- rownames(head(IVD.All.OAF2.markers, n = 10))
IVD.All.OAF1.top10 <- rownames(head(IVD.All.OAF1.markers, n = 10))
IVD.All.IAF2.top10 <- rownames(head(IVD.All.IAF2.markers, n = 10))
IVD.All.IAF1.top10 <- rownames(head(IVD.All.IAF1.markers, n = 10))
IVD.All.NP4.top10 <- rownames(head(IVD.All.NP4.markers, n = 10))
IVD.All.NP3.top10 <- rownames(head(IVD.All.NP3.markers, n = 10))
IVD.All.NP2.top10 <- rownames(head(IVD.All.NP2.markers, n = 10))
IVD.All.NP1.top10 <- rownames(head(IVD.All.NP1.markers, n = 10))
IVD.All.IC.top10 <- rownames(head(IVD.All.IC.markers, n = 10))
IVD.All.NC2.top10 <- rownames(head(IVD.All.NC2.markers, n = 10))
IVD.All.NC1.top10 <- rownames(head(IVD.All.NC1.markers, n = 10))

nIVD.RBC.markers <- FindMarkers(IVD.nIVD, ident.1 = "RBC", min.pct = 0.25)
nIVD.IC.markers <- FindMarkers(IVD.nIVD, ident.1 = "IC", min.pct = 0.25)
nIVD.OAF3.markers <- FindMarkers(IVD.nIVD, ident.1 = "OAF3", min.pct = 0.25)
nIVD.OAF2.markers <- FindMarkers(IVD.nIVD, ident.1 = "OAF2", min.pct = 0.25)
nIVD.OAF1.markers <- FindMarkers(IVD.nIVD, ident.1 = "OAF1", min.pct = 0.25)
nIVD.IAF2.markers <- FindMarkers(IVD.nIVD, ident.1 = "IAF2", min.pct = 0.25)
nIVD.IAF1.markers <- FindMarkers(IVD.nIVD, ident.1 = "IAF1", min.pct = 0.25)
nIVD.NP4.markers <- FindMarkers(IVD.nIVD, ident.1 = "NP4", min.pct = 0.25)
nIVD.NP3.markers <- FindMarkers(IVD.nIVD, ident.1 = "NP3", min.pct = 0.25)
nIVD.NP2.markers <- FindMarkers(IVD.nIVD, ident.1 = "NP2", min.pct = 0.25)
nIVD.NP1.markers <- FindMarkers(IVD.nIVD, ident.1 = "NP1", min.pct = 0.25)
nIVD.DDSC.markers <- FindMarkers(IVD.nIVD, ident.1 = "DDSC", min.pct = 0.25)
nIVD.NC2.markers <- FindMarkers(IVD.nIVD, ident.1 = "NC2", min.pct = 0.25)
nIVD.NC1.markers <- FindMarkers(IVD.nIVD, ident.1 = "NC1", min.pct = 0.25)

aIVD.RBC.markers <- FindMarkers(IVD.aIVD, ident.1 = "RBC", min.pct = 0.25)
aIVD.IC.markers <- FindMarkers(IVD.aIVD, ident.1 = "IC", min.pct = 0.25)
aIVD.OAF3.markers <- FindMarkers(IVD.aIVD, ident.1 = "OAF3", min.pct = 0.25)
aIVD.OAF2.markers <- FindMarkers(IVD.aIVD, ident.1 = "OAF2", min.pct = 0.25)
aIVD.OAF1.markers <- FindMarkers(IVD.aIVD, ident.1 = "OAF1", min.pct = 0.25)
aIVD.IAF2.markers <- FindMarkers(IVD.aIVD, ident.1 = "IAF2", min.pct = 0.25)
aIVD.IAF1.markers <- FindMarkers(IVD.aIVD, ident.1 = "IAF1", min.pct = 0.25)
aIVD.NP4.markers <- FindMarkers(IVD.aIVD, ident.1 = "NP4", min.pct = 0.25)
aIVD.NP3.markers <- FindMarkers(IVD.aIVD, ident.1 = "NP3", min.pct = 0.25)
aIVD.NP2.markers <- FindMarkers(IVD.aIVD, ident.1 = "NP2", min.pct = 0.25)
aIVD.NP1.markers <- FindMarkers(IVD.aIVD, ident.1 = "NP1", min.pct = 0.25)
aIVD.DDSC.markers <- FindMarkers(IVD.aIVD, ident.1 = "DDSC", min.pct = 0.25)
aIVD.NC2.markers <- FindMarkers(IVD.aIVD, ident.1 = "NC2", min.pct = 0.25)
aIVD.NC1.markers <- FindMarkers(IVD.aIVD, ident.1 = "NC1", min.pct = 0.25)

levels(x = IVD.integrated3) <- c('RBC', 'IC', "non-MC", "OAF3", 'OAF2', "OAF1", "IAF1", "IAF2", "NP4", "NP3", "NP2", "NP1", "NC2", "NC1")
levels(x = IVD.aIVD) <- c('RBC', 'IC', "non-MC", "OAF3", 'OAF2', "OAF1", "IAF1", "IAF2", "NP4", "NP3", "NP2", "NP1", "NC2", "NC1")
levels(x = IVD.nIVD) <- c('RBC', 'IC', "non-MC", "OAF3", 'OAF2', "OAF1", "IAF1", "IAF2", "NP4", "NP3", "NP2", "NP1", "NC2", "NC1")

DimPlot(IVD.integrated3, reduction = "umap", pt.size = 0.2, label.size = FALSE, split.by = "MajorType",
                  cols = c("azure4","darkolivegreen","azure4","azure4","azure4","azure4","azure4","azure4","azure4","azure4","azure4","azure4", "springgreen3","forestgreen")
) + theme(aspect.ratio = "1")

DimPlot(IVD.integrated3, reduction = "umap", pt.size = 0.2, label.size = FALSE, split.by = "MajorType",
                  cols = c("azure4","azure4","azure4","azure4","azure4","azure4","azure4","azure4","dodgerblue3","dodgerblue2","dodgerblue1","dodgerblue4", "azure4","azure4")
) + theme(aspect.ratio = "1")


DimPlot(IVD.integrated3, reduction = "umap", pt.size = 0.2, label.size = FALSE, split.by = "MajorType",
                  cols = c("azure4","azure4","azure4","darkorange3","darkorange1","orange","orchid3","plum3","azure4","azure4","azure4", "azure4","azure4", "azure4")
) + theme(aspect.ratio = "1")

Plot9b1 <- FeaturePlot(IVD.integrated3, features = "SOX17", ncol = 1, label = FALSE, split.by = "MajorType",  order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
Plot9b1

Plot9b2 <- FeaturePlot(IVD.integrated3, features = "BASP1", ncol = 1, label = FALSE, split.by = "MajorType", order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
Plot9b2

Plot9b3 <- FeaturePlot(IVD.integrated3, features = "CD44", ncol = 1, label = FALSE, split.by = "MajorType", order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
Plot9b3

Plot9b4 <- FeaturePlot(IVD.integrated3, features = "HSPA1A", ncol = 1, label = FALSE, split.by = "MajorType", order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
Plot9b4

Plot9b5 <- FeaturePlot(IVD.integrated3, features = "OGN", ncol = 1, label = FALSE, split.by = "MajorType", order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
Plot9b5

Plot9b6 <- FeaturePlot(IVD.integrated3, features = "FGFBP2", ncol = 1, label = FALSE, split.by = "MajorType", order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
Plot9b6

Plot9b7 <- FeaturePlot(IVD.integrated3, features = "ANXA1", ncol = 1, label = FALSE, split.by = "MajorType", order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
Plot9b7

Plot9b8 <- FeaturePlot(IVD.integrated3, features = "GADD45G", ncol = 1, label = FALSE, split.by = "MajorType", order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
Plot9b8

Plot9b9 <- FeaturePlot(IVD.integrated3, features = "S100A6", ncol = 1, label = FALSE, split.by = "MajorType", order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
Plot9b9

Plot9b10 <- FeaturePlot(IVD.integrated3, features = "LGALS1", ncol = 1, label = FALSE, split.by = "MajorType", order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
Plot9b10

Plot9b11 <- FeaturePlot(IVD.integrated3, features = "CCN2", ncol = 1, label = FALSE, split.by = "MajorType", order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
Plot9b11

Plot9b12 <- FeaturePlot(IVD.integrated3, features = "MRGPRF", ncol = 1, label = FALSE, split.by = "MajorType", order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
Plot9b12

Plot9b13 <- FeaturePlot(IVD.integrated3, features = "FTH1", ncol = 1, label = FALSE, split.by = "MajorType", order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
Plot9b13

Plot9b14 <- FeaturePlot(IVD.integrated3, features = "ZNF385D", ncol = 1, label = FALSE, split.by = "MajorType", order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
Plot9b14

Plot9b15 <- FeaturePlot(IVD.integrated3, features = "HSPB1", ncol = 1, label = FALSE, split.by = "MajorType", order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
Plot9b15

Plot9b16 <- FeaturePlot(IVD.integrated3, features = "HSPA6", ncol = 1, label = FALSE, split.by = "MajorType", order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
Plot9b16

Plot9b17 <- FeaturePlot(IVD.integrated3, features = "HES1", ncol = 1, label = FALSE, split.by = "MajorType", order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
Plot9b17

Plot9b18 <- FeaturePlot(IVD.integrated3, features = "FABP5", ncol = 1, label = FALSE, split.by = "MajorType", order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
Plot9b18

levels(x = IVD.nIVD) <- c('NC1', 'NC2',  'NP1', "NP2", "NP3", "NP4", "IAF2", "IAF1", "OAF1", "OAF2", "OAF3", "non-MC", "IC", "RBC")
levels(x = IVD.aIVD) <- c('NC1', 'NC2',  'NP1', "NP2", "NP3", "NP4", "IAF2", "IAF1", "OAF1", "OAF2", "OAF3", "non-MC", "IC", "RBC")

Plot12b1 <- VlnPlot(IVD.nIVD, features = "LGALS1", idents = c("IAF2", "IAF1", "OAF1", "OAF2", "OAF3"),
                    slot = "counts", log = TRUE, cols = c("plum2", "orchid3", "darkorange", "orange2", "chocolate3") , pt.size = 0.01, sort = FALSE, ncol=1)

Plot12b1

Plot12b2 <- VlnPlot(IVD.nIVD, features = "OGN", idents = c("IAF2", "IAF1", "OAF1", "OAF2", "OAF3"),
                    slot = "counts", log = TRUE, cols = c("plum2", "orchid3", "darkorange", "orange2", "chocolate3") , pt.size = 0.01, sort = FALSE,  ncol=1)
Plot12b2

Plot12b3 <- VlnPlot(IVD.nIVD, features = "HES1", idents = c("IAF2", "IAF1", "OAF1", "OAF2", "OAF3"),
                    slot = "counts", log = TRUE, cols = c("plum2", "orchid3", "darkorange", "orange2", "chocolate3") , pt.size = 0.01, sort = FALSE,  ncol=1)
Plot12b3

Plot12b4 <- VlnPlot(IVD.nIVD, features = "ZNF385D", idents = c("IAF2", "IAF1", "OAF1", "OAF2", "OAF3"),
                    slot = "counts", log = TRUE, cols = c("plum2", "orchid3", "darkorange", "orange2", "chocolate3") , pt.size = 0.01, sort = FALSE,  ncol=1)
Plot12b4

Plot12b5 <- VlnPlot(IVD.nIVD, features = "HSPA1A", idents = c("IAF2", "IAF1", "OAF1", "OAF2", "OAF3"),
                    slot = "counts", log = TRUE, cols = c("plum2", "orchid3", "darkorange", "orange2", "chocolate3") , pt.size = 0.01, sort = FALSE,  ncol=1)
Plot12b5


Plot12b6 <- VlnPlot(IVD.nIVD, features = c("OGN", "LGALS1", "HES1", "ZNF385D", "HSPA1A"), idents = c("IAF2", "IAF1", "OAF1", "OAF2", "OAF3"),
                    slot = "counts", log = TRUE, cols = c("plum2", "orchid3", "darkorange", "orange2", "chocolate3") , pt.size = 0.01, sort = FALSE,  ncol=5)
Plot12b6


Plot12c1 <- VlnPlot(IVD.aIVD, features = "LGALS1", idents = c("IAF2", "IAF1", "OAF1", "OAF2", "OAF3"),
                    slot = "counts", log = TRUE, cols = c("plum2", "orchid3", "darkorange", "orange2", "chocolate3") , pt.size = 0.01, sort = FALSE, ncol=1)

Plot12c1

Plot12c2 <- VlnPlot(IVD.aIVD, features = "OGN", idents = c("IAF2", "IAF1", "OAF1", "OAF2", "OAF3"),
                    slot = "counts", log = TRUE, cols = c("plum2", "orchid3", "darkorange", "orange2", "chocolate3") , pt.size = 0.01, sort = FALSE,  ncol=1)
Plot12c2

Plot12c3 <- VlnPlot(IVD.aIVD, features = "HES1", idents = c("IAF2", "IAF1", "OAF1", "OAF2", "OAF3"),
                    slot = "counts", log = TRUE, cols = c("plum2", "orchid3", "darkorange", "orange2", "chocolate3") , pt.size = 0.01, sort = FALSE,  ncol=1)
Plot12c3

Plot12c4 <- VlnPlot(IVD.aIVD, features = "ZNF385D", idents = c("IAF2", "IAF1", "OAF1", "OAF2", "OAF3"),
                    slot = "counts", log = TRUE, cols = c("plum2", "orchid3", "darkorange", "orange2", "chocolate3") , pt.size = 0.01, sort = FALSE,  ncol=1)
Plot12c4

Plot12c5 <- VlnPlot(IVD.aIVD, features = "HSPA1A", idents = c("IAF2", "IAF1", "OAF1", "OAF2", "OAF3"),
                    slot = "counts", log = TRUE, cols = c("plum2", "orchid3", "darkorange", "orange2", "chocolate3") , pt.size = 0.01, sort = FALSE,  ncol=1)
Plot12c5



Plot12c6 <- VlnPlot(IVD.aIVD, features = c("OGN", "LGALS1", "HES1", "ZNF385D", "HSPA1A"), idents = c("IAF2", "IAF1", "OAF1", "OAF2", "OAF3"),
                    slot = "counts", log = TRUE, cols = c("plum2", "orchid3", "darkorange", "orange2", "chocolate3") , pt.size = 0.01, sort = FALSE,  ncol=5)
Plot12c6



Plot14b1 <- VlnPlot(IVD.nIVD, idents = c("IAF1", "IAF2", "OAF1", "OAF2", "OAF3"), same.y.lims = TRUE,
                    features = c("COL12A1", "COL6A3", 'COL5A2', "COL5A1", "COL3A1", "COL1A1"),
                    cols = c("plum2", "orchid3", "darkorange", "orange2", "chocolate3"), stack = TRUE,
                    slot = "counts", log = TRUE,  pt.size = 0, sort = FALSE, fill.by = "ident") 


Plot14b2 <- VlnPlot(IVD.aIVD, idents = c("IAF1", "IAF2", "OAF1", "OAF2", "OAF3"), same.y.lims = TRUE,
                    features = c("COL12A1", "COL6A3", 'COL5A2', "COL5A1", "COL3A1", "COL1A1"),
                    cols = c("plum2", "orchid3", "darkorange", "orange2", "chocolate3"), stack = TRUE,
                    slot = "counts", log = TRUE,  pt.size = 0, sort = FALSE, fill.by = "ident") 

Plot14b1
Plot14b2

Plot14c1 <- VlnPlot(IVD.nIVD, idents = c("IAF1", "IAF2", "OAF1", "OAF2", "OAF3"), y.max = 1000, same.y.lims = TRUE,
                    features = c("COL12A1", "COL6A3", 'COL5A2', "COL5A1", "COL3A1", "COL1A1"),
                    cols = c("plum2", "orchid3", "darkorange", "orange2", "chocolate3"), stack = TRUE,
                    slot = "counts", log = TRUE,  pt.size = 0, sort = FALSE, fill.by = "ident")
                    

Plot14c2 <- VlnPlot(IVD.aIVD, idents = c("IAF1", "IAF2", "OAF1", "OAF2", "OAF3"), same.y.lims = TRUE, y.max = 1000,
                    features = c("COL12A1", "COL6A3", 'COL5A2', "COL5A1", "COL3A1", "COL1A1"),
                    cols = c("plum2", "orchid3", "darkorange", "orange2", "chocolate3"), stack = TRUE,
                    slot = "counts", log = TRUE,  pt.size = 0, sort = FALSE, fill.by = "ident") 

Plot14c1
Plot14c2


nIVD.NC1.markers.volcano.highlight <- c('ZEB2','TFPI', "IGFBP7", "AKR1C1", "FABP4", "TMSB4X", "COL4A1")
nIVD.NC2.markers.volcano.highlight <- c('BMX','PCAT19', "MYCT1", "TMSB10", "CAVIN2", "SOX17", "GNG11", "CXCL8", "FABP4", "TMSB4X")

aIVD.NC1.markers.volcano.highlight <- c('FABP4','COL4A1', "APOE", "AKR1C1")
aIVD.NC2.markers.volcano.highlight <- c("None")



Plot17a1 <- EnhancedVolcano(nIVD.NC1.markers,
                lab = rownames(nIVD.NC1.markers),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "NC1 (Neonatal IVD)",
                pCutoff = 0.0001,
                FCcutoff = 0.5,
                labSize = 4,
                pointSize = 2.5,
                xlim = c(-7,9),
                col = c("grey30", "grey30", "grey30", "red2"),
                encircle = nIVD.NC1.markers.volcano.highlight,
                encircleCol = 'black',
                encircleSize = 0.5,
                encircleFill = 'pink',
                encircleAlpha = 1/3,
                legendPosition = 'none',
                drawConnectors = TRUE,
                widthConnectors = 0.35,
                selectLab = nIVD.NC1.markers.volcano.highlight) + theme(aspect.ratio = 1)

Plot17a2 <- EnhancedVolcano(nIVD.NC2.markers,
                lab = rownames(nIVD.NC2.markers),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "NC2 (Neonatal IVD)",
                pCutoff = 0.0001,
                FCcutoff = 0.5,
                labSize = 4,
                pointSize = 2.5,
                xlim = c(-7.5,13),
                col = c("grey30", "grey30", "grey30", "red2"),
                encircle = nIVD.NC2.markers.volcano.highlight,
                encircleCol = 'black',
                encircleSize = 0.5,
                encircleFill = 'pink',
                encircleAlpha = 1/3,
                legendPosition = 'none',
                drawConnectors = TRUE,
                widthConnectors = 0.35,
                selectLab = nIVD.NC2.markers.volcano.highlight) + theme(aspect.ratio = 1)



Plot17b1 <- EnhancedVolcano(aIVD.NC1.markers,
                            lab = rownames(aIVD.NC1.markers),
                            x = 'avg_log2FC',
                            y = 'p_val_adj',
                            title = "NC1 (Adult IVD)",
                            pCutoff = 0.0001,
                            FCcutoff = 0.5,
                            labSize = 4,
                            pointSize = 2.5,
                            xlim = c(-7.5,10),
                            ylim = c(0, 11.5),
                            col = c("grey30", "grey30", "grey30", "red2"),
                            encircle = aIVD.NC1.markers.volcano.highlight,
                            encircleCol = 'black',
                            encircleSize = 0.5,
                            encircleFill = 'pink',
                            encircleAlpha = 1/3,
                            legendPosition = 'none',
                            drawConnectors = TRUE,
                            widthConnectors = 0.35,
                            selectLab = aIVD.NC1.markers.volcano.highlight)+ theme(aspect.ratio = 1)

Plot17b2 <- EnhancedVolcano(aIVD.NC2.markers,
                            lab = rownames(aIVD.NC2.markers),
                            x = 'avg_log2FC',
                            y = 'p_val_adj',
                            title = "NC2 (Adult IVD)",
                            pCutoff = 0.0001,
                            FCcutoff = 0.5,
                            labSize = 4,
                            pointSize = 2.5,
                            xlim = c(-8,8),
                            ylim = c(0, 6),
                            col = c("grey30", "grey30", "grey30", "red2"),

                            legendPosition = 'none',
                            drawConnectors = TRUE,
                            widthConnectors = 0.35
                     )+ theme(aspect.ratio = 1)




Plot17a1
Plot17a2

Plot17b1
Plot17b2

Idents(IVD.nIVD) <- "Subtype"
Idents(IVD.aIVD) <- "Subtype"

levels(x = IVD.nIVD) <- c('RBC', 'IC', "non-MC", "OAF3", 'OAF2', "OAF1", "IAF1", "IAF2", "NP4", "NP3", "NP2", "NP1", "NC2", "NC1")
levels(x = IVD.aIVD) <- c('RBC', 'IC', "non-MC", "OAF3", 'OAF2', "OAF1", "IAF1", "IAF2", "NP4", "NP3", "NP2", "NP1", "NC2", "NC1")

IVD.nIVD <- RenameIdents(object = IVD.nIVD, 
                         "RBC" = "Other Clusters", "IC" = "Other Clusters", "OAF3" = "Other Clusters", "OAF2" = "Other Clusters",
                         "OAF1" = "Other Clusters","IAF1" = "Other Clusters","IAF2" = "Other Clusters","NP4" = "Other Clusters",
                         "NP3" = "Other Clusters","NP2" = "Other Clusters", "NP1" = "Other Clusters", "non-MC" = "Other Clusters")
IVD.nIVD@active.ident
Plot23a1 <- VlnPlot(IVD.nIVD, features = c("MAP1B", "SOX4","SOX17", "ITGA6", "BASP1", "CD44", "CD14"
                                           ), 
                    slot = "counts", log = TRUE, same.y.lims = TRUE,
                    cols = c("azure4", "darkolivegreen", "forestgreen", "springgreen3"),
                    pt.size = 0.2, sort = FALSE, stack = TRUE, fill.by = "ident") + theme(aspect.ratio = 1)

IVD.aIVD <- RenameIdents(object = IVD.aIVD, 
                         "RBC" = "Other Clusters", "IC" = "Other Clusters", "OAF3" = "Other Clusters", "OAF2" = "Other Clusters",
                         "OAF1" = "Other Clusters","IAF1" = "Other Clusters","IAF2" = "Other Clusters","NP4" = "Other Clusters",
                         "NP3" = "Other Clusters","NP2" = "Other Clusters", "NP1" = "Other Clusters", "non-MC" = "Other Clusters")
IVD.aIVD@active.ident
Plot23a2 <- VlnPlot(IVD.aIVD, features = c("MAP1B", "SOX4","SOX17", "ITGA6", "BASP1", "CD44", "CD14"
), 
slot = "counts", log = TRUE, same.y.lims = TRUE,
cols = c("azure4", "darkolivegreen", "forestgreen", "springgreen3"),
pt.size = 0.2, sort = FALSE, stack = TRUE, fill.by = "ident") + theme(aspect.ratio = 1)

Plot23a1
Plot23a2

Idents(IVD.nIVD) <- "Subtype"
Idents(IVD.aIVD) <- "Subtype"

levels(x = IVD.nIVD) <- c('RBC', 'IC', "non-MC", "OAF3", 'OAF2', "OAF1", "IAF1", "IAF2", "NP4", "NP3", "NP2", "NP1", "NC2", "NC1")
levels(x = IVD.aIVD) <- c('RBC', 'IC', "non-MC", "OAF3", 'OAF2', "OAF1", "IAF1", "IAF2", "NP4", "NP3", "NP2", "NP1", "NC2", "NC1")

IVD.nIVD <- RenameIdents(object = IVD.nIVD, 
                         "RBC" = "Others", "IC" = "Others", "OAF3" = "Others", "OAF2" = "Others",
                         "OAF1" = "Others","IAF1" = "Others","IAF2" = "Others",
                        "NC1" = "Others", "NC2" = "Others", "non-MC" = "Others")
IVD.nIVD@active.ident
Plot23b1 <- VlnPlot(IVD.nIVD, features = c("COL2A1", "ACAN",  "FABP5", "OGN", "FGFBP2", "ANXA1", "COL1A1"), 
                    slot = "counts", log = TRUE, same.y.lims = TRUE,
                    cols = c("azure4", "dodgerblue4", "dodgerblue1", "dodgerblue3", "dodgerblue2"),
                    pt.size = 0.2, sort = FALSE, stack = TRUE, fill.by = "ident") + theme(aspect.ratio = 1)


IVD.aIVD <- RenameIdents(object = IVD.aIVD, 
                         "RBC" = "Others", "IC" = "Others", "OAF3" = "Others", "OAF2" = "Others",
                         "OAF1" = "Others","IAF1" = "Others","IAF2" = "Others",
                         "NC1" = "Others", "NC2" = "Others", "non-MC" = "Others")
IVD.aIVD@active.ident
Plot23b2 <- VlnPlot(IVD.aIVD, features = c("COL2A1", "ACAN",  "FABP5", "OGN", "FGFBP2", "ANXA1", "COL1A1"), 
                    slot = "counts", log = TRUE, same.y.lims = TRUE,
                    cols = c("azure4", "dodgerblue4", "dodgerblue1", "dodgerblue3", "dodgerblue2"),
                    pt.size = 0.2, sort = FALSE, stack = TRUE, fill.by = "ident") + theme(aspect.ratio = 1)

Plot23b1
Plot23b2

Idents(IVD.nIVD) <- "Subtype"
Idents(IVD.aIVD) <- "Subtype"

levels(x = IVD.nIVD) <- c('RBC', 'IC', "non-MC", "OAF3", 'OAF2', "OAF1", "IAF1", "IAF2", "NP4", "NP3", "NP2", "NP1", "NC2", "NC1")
levels(x = IVD.aIVD) <- c('RBC', 'IC', "non-MC", "OAF3", 'OAF2', "OAF1", "IAF1", "IAF2", "NP4", "NP3", "NP2", "NP1", "NC2", "NC1")

IVD.nIVD <- RenameIdents(object = IVD.nIVD, 
                         "RBC" = "Others", "IC" = "Others", "NP4" = "Others",
                         "NP3" = "Others","NP2" = "Others", "NP1" = "Others",
                         "NC1" = "Others", "NC2" = "Others", "non-MC" = "Others")
IVD.nIVD@active.ident
Plot23c1 <- VlnPlot(IVD.nIVD, features = c("COL1A1","CALR", "HSPA6","COL2A1", "ACAN", "MAP1B", "CD44"), 
                    slot = "counts", log = TRUE, same.y.lims = TRUE,
                    cols = c("Azure4",  "chocolate3", "orange2","darkorange","orchid3",  "plum2"),
                    pt.size = 0.2, sort = FALSE, stack = TRUE, fill.by = "ident") + theme(aspect.ratio = 1)


IVD.aIVD <- RenameIdents(object = IVD.aIVD, 
                         "RBC" = "Others", "IC" = "Others", "NP4" = "Others",
                         "NP3" = "Others","NP2" = "Others", "NP1" = "Others",
                         "NC1" = "Others", "NC2" = "Others", "non-MC" = "Others")
IVD.aIVD@active.ident
Plot23c2 <- VlnPlot(IVD.aIVD, features = c("COL1A1","CALR", "HSPA6","COL2A1", "ACAN", "MAP1B", "CD44"), 
                    slot = "counts", log = TRUE, same.y.lims = TRUE,
                    cols = c("Azure4",  "chocolate3", "orange2","darkorange","orchid3",  "plum2"),
                    pt.size = 0.2, sort = FALSE, stack = TRUE, fill.by = "ident") + theme(aspect.ratio = 1)

Plot23c1
Plot23c2

Idents(IVD.nIVD) <- "Subtype"
Idents(IVD.aIVD) <- "Subtype"

levels(x = IVD.nIVD) <- c('RBC', 'IC', "non-MC", "OAF3", 'OAF2', "OAF1", "IAF1", "IAF2", "NP4", "NP3", "NP2", "NP1", "NC2", "NC1")
levels(x = IVD.aIVD) <- c('RBC', 'IC', "non-MC", "OAF3", 'OAF2', "OAF1", "IAF1", "IAF2", "NP4", "NP3", "NP2", "NP1", "NC2", "NC1")

Plot24a <- VlnPlot(IVD.integrated3, features = c("AKR1C1", "FABP4", "APOE"), 
                    slot = "counts", log = TRUE, same.y.lims = TRUE, split.by = "MajorType", 
                    cols = c("lightpink3", "lightblue4", "springgreen3", "azure4", "azure4", "azure4", "azure4", "azure4", "azure4", "azure4", "azure4", "darkolivegreen", "forestgreen", "springgreen3"),
                    pt.size = 0.2, sort = FALSE, stack = TRUE, fill.by = "feature") + theme(aspect.ratio = 5)



Plot24a


Idents(IVD.nIVD) <- "Subtype"
Idents(IVD.aIVD) <- "Subtype"

levels(x = IVD.nIVD) <- c('RBC', 'IC', "non-MC", "OAF3", 'OAF2', "OAF1", "IAF1", "IAF2", "NP4", "NP3", "NP2", "NP1", "NC2", "NC1")
levels(x = IVD.aIVD) <- c('RBC', 'IC', "non-MC", "OAF3", 'OAF2', "OAF1", "IAF1", "IAF2", "NP4", "NP3", "NP2", "NP1", "NC2", "NC1")

Plot25b1 <- DotPlot(IVD.nIVD, features = c(
  "MAP1B", "SOX4","SOX17", "ITGA6", "BASP1", 
  "ACAN", "COL2A1", "SOX9",
  "COL1A1", "CALR", "HSPA6",
  "MEG3", "CD44", "CD14", "HBB", "HBA1"), cols = c("blue4",  "sienna2"),
  scale.by = "radius", dot.scale = 9, col.min = -1, col.max = 1) + 
  theme(axis.text.x = element_text(angle = 90)) 

Plot25b2 <- DotPlot(IVD.aIVD, features = c(
  "MAP1B", "SOX4","SOX17", "ITGA6", "BASP1", 
  "ACAN", "COL2A1", "SOX9",
  "COL1A1", "CALR", "HSPA6",
  "MEG3", "CD44", "CD14", "HBB", "HBA1"), cols = c("blue4",  "sienna2"),
  scale.by = "radius", dot.scale = 9, col.min = -1, col.max = 1) + 
  theme(axis.text.x = element_text(angle = 90)) 

Plot25c1 <- DotPlot(IVD.nIVD, features = c("MAP1B","SOX4","SOX17", "ITGA6", "BASP1", "AKR1C1",  "FABP4", "APOE", "ACAN", "COL2A1", "SOX9",
                                           "COL1A1", "CALR", "HSPA6", "OGN", "LGALS1", "HES1",  "ZNF385D", "HSPA1A", "MEG3", "CD44", "CD14", "HBB", "HBA1"), cols = c("blue4",  "sienna2"),
  scale.by = "radius", dot.scale = 9, col.min = -1, col.max = 1) + 
  theme(axis.text.x = element_text(angle = 90)) 

Plot25c2 <- DotPlot(IVD.aIVD, features = c("MAP1B", "SOX4","SOX17", "ITGA6", "BASP1", "AKR1C1",  "FABP4", "APOE", "ACAN", "COL2A1", "SOX9",
                                           "COL1A1", "CALR", "HSPA6", "OGN", "LGALS1", "HES1",  "ZNF385D", "HSPA1A", "MEG3", "CD44", "CD14", "HBB", "HBA1"), cols = c("blue4",  "sienna2"),
                    scale.by = "radius", dot.scale = 9, col.min = -1, col.max = 1) + 
  theme(axis.text.x = element_text(angle = 90)) 

Plot25b1
Plot25b2
Plot25c1
Plot25c2

#---------The End-----------------