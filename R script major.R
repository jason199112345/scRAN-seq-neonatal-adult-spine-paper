#===============Reading Data and Create Seurat Objects ==================

IVD3=read.csv("/Users/jiangw2/Dropbox/Research/Dmitriy Sheyn lab/03b-scRNAseq for notochordal/00_Original Data/IVD3.csv", sep=",", header=T)
IVD3m= as.matrix(IVD3[, 3:2495])
row.names(IVD3m)=IVD3$HGNC.symbol
IVD3m=as.data.frame(IVD3m)
IVD3m=as.matrix(IVD3m)
nIVD1 <- CreateSeuratObject(counts = IVD3m, project = "IVD3m", min.cells = 3, min.features = 200)


IVD4=read.csv("/Users/jiangw2/Dropbox/Research/Dmitriy Sheyn lab/03b-scRNAseq for notochordal/00_Original Data/IVD4.csv", sep=",", header=T)
IVD4m= as.matrix(IVD4[, 3:2059])
row.names(IVD4m)=IVD4$HGNC.symbol
IVD4m=as.data.frame(IVD4m)
IVD4m=as.matrix(IVD4m)
nIVD2 <- CreateSeuratObject(counts = IVD4m, project = "IVD4m", min.cells = 3, min.features = 200)


IVD5=read.csv("/Users/jiangw2/Dropbox/Research/Dmitriy Sheyn lab/03b-scRNAseq for notochordal/00_Original Data/IVD5.csv", sep=",", header=T)
IVD5m= as.matrix(IVD5[, 3:3197])
row.names(IVD5m)=IVD5$HGNC.symbol
IVD5m=as.data.frame(IVD5m)
IVD5m=as.matrix(IVD5m)
nIVD3 <- CreateSeuratObject(counts = IVD5m, project = "IVD5m", min.cells = 3, min.features = 200)


aIVD=read.csv("/Users/jiangw2/Dropbox/Research/Dmitriy Sheyn lab/03b-scRNAseq for notochordal/00_Original Data/NPCX.csv", sep=",", header=T)
aIVDm= as.matrix(aIVD[, 3:2702])
row.names(aIVDm)=aIVD$HGNC.symbol
aIVDm=as.data.frame(aIVDm)
aIVDm=as.matrix(aIVDm)
aIVD1 <- CreateSeuratObject(counts = aIVDm, project = "aIVDm", min.cells = 3, min.features = 200)

Data_A <- Read10X(
  data.dir = "/Users/jiangw2/Dropbox/Research/Dmitriy Sheyn lab/03b-scRNAseq for notochordal/00_Original Data/Cadaver-#A2-011.2/",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

aIVD2 = CreateSeuratObject(counts = Data_A, project = "aIVD2")

Data_B <- Read10X(
  data.dir = "/Users/jiangw2/Dropbox/Research/Dmitriy Sheyn lab/03b-scRNAseq for notochordal/00_Original Data/NPC-#119.1/",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

aIVD3 = CreateSeuratObject(counts = Data_B, project = "aIVD3")

DataAll <- c(nIVD1, nIVD2, nIVD3, aIVD1, aIVD2, aIVD3)


for (i in 1:length(DataAll)) {
  DataAll[[i]] <- NormalizeData(DataAll[[i]], verbose = FALSE)
  DataAll[[i]] <- FindVariableFeatures(DataAll[[i]], selection.method = "vst", 
                                       nfeatures = 6000, verbose = FALSE)
}

IVD.anchors <- FindIntegrationAnchors(object.list = DataAll, anchor.features = features)
IVD.integrated <- IntegrateData(anchorset = IVD.anchors)
DefaultAssay(IVD.integrated) <- "integrated"

IVD.integrated <- ScaleData(IVD.integrated, verbose = FALSE)
IVD.integrated <- RunPCA(IVD.integrated, npcs = 30, verbose = FALSE)
IVD.integrated <- RunUMAP(IVD.integrated, reduction = "pca", dims = 1:30)
DimPlot(IVD.integrated, reduction = "umap", pt.size = 2, label.size = 20) + theme(aspect.ratio = "1")

IVD.integrated3 <- FindNeighbors(IVD.integrated, reduction = "pca", dims = 1:30)
IVD.integrated3 <- FindClusters(IVD.integrated3, resolution = 0.5)
DefaultAssay(IVD.integrated3) <- "RNA"
DimPlot(IVD.integrated3, reduction = "umap", pt.size = 0.3, label = TRUE, label.size = 10) + theme(aspect.ratio = "1")

IVD.integrated3$ClusterNumber <- Idents(IVD.integrated3)
Idents(IVD.integrated3) <- "orig.ident"
IVD.integrated3 <- RenameIdents(object = IVD.integrated3, "IVD3m" = "nIVD", "IVD4m" = "nIVD", 'IVD5m' = 'nIVD', 'aIVDm' = 'aIVD', 'aIVD2' = 'aIVD', 'aIVD3' = 'aIVD')
IVD.integrated3$MajorType <- Idents(IVD.integrated3)
table(IVD.integrated3$MajorType)
Idents(IVD.integrated3) <- "ClusterNumber"

Idents(IVD.integrated3) <- "MajorType"
IVD.nIVD <- subset(x = IVD.integrated3, idents = "nIVD")
IVD.aIVD <- subset(x = IVD.integrated3, idents = "aIVD")
Idents(IVD.integrated3) <- "ClusterNumber"
Idents(IVD.nIVD) <- "ClusterNumber"
Idents(IVD.aIVD) <- "ClusterNumber"
DimPlot(IVD.nIVD,  reduction = "umap", pt.size = 0.1, split.by = "MajorType", label = TRUE, label.size = 5,
        cols = c("red", "red3", "gray35", "gray70", "darkorange3", "darkorange2", "darkorange1", "plum", "dodgerblue4", "dodgerblue3", "dodgerblue2", "dodgerblue1", "cyan3", "forestgreen", "darkolivegreen")) + theme(aspect.ratio = "1") 


#Assigning cluster identity & Figure 1

IVD.integrated3 <- RenameIdents(object = IVD.integrated3, "0" = "NP1", "1" = "NP2", "2" = "NC2", "3" = "OAF3", "4" = "OAF2", "5" = "non-MC", "6" =  "OAF1", "7" = "IAF1", "8" = "NP3", "9" = "NP4", "10" = "RBC2", "11" = "NC1", "12" = "RBC1", "13" = "IC", "14" = "NCxNP")
IVD.nIVD <- RenameIdents(object = IVD.nIVD, "0" = "NP1", "1" = "NP2", "2" = "NC2", "3" = "OAF3", "4" = "OAF2", "5" = "non-MC", "6" =  "OAF1", "7" = "IAF1", "8" = "NP3", "9" = "NP4", "10" = "RBC2", "11" = "NC1", "12" = "RBC1", "13" = "IC", "14" = "NCxNP")
IVD.aIVD <- RenameIdents(object = IVD.aIVD, "0" = "NP1", "1" = "NP2", "2" = "NC2", "3" = "OAF3", "4" = "OAF2", "5" = "non-MC", "6" =  "OAF1", "7" = "IAF1", "8" = "NP3", "9" = "NP4", "10" = "RBC2", "11" = "NC1", "12" = "RBC1", "13" = "IC", "14" = "NCxNP")
IVD.integrated3$Subtype <- Idents(IVD.integrated3)
IVD.nIVD$Subtype <- Idents(IVD.nIVD)
IVD.aIVD$Subtype <- Idents(IVD.aIVD)

levels(x = IVD.integrated3) <- c('RBC2', "RBC1", "IC", "non-MC", 'OAF3', "OAF2", "OAF1", "IAF1", "NP4", "NP3", "NP2", "NP1", "NCxNP", "NC2", "NC1")

p00a <- DimPlot(IVD.integrated3,  reduction = "umap", pt.size = 0.1, split.by = "MajorType", label = FALSE, label.size = 5,
                cols = c("red", "red3", "gray35", "gray70", "darkorange3", "darkorange2", "darkorange1", "plum", "dodgerblue4", "dodgerblue3", "dodgerblue2", "dodgerblue1", "cyan3", "forestgreen", "darkolivegreen")) + theme(aspect.ratio = "1") 

FeaturePlot(IVD.integrated3, features = c("MAP1B", "ACAN", "COL1A1", "CALR", "HSPA6", "MEG3", "HBB", "CD24", "TBXT"), ncol = 3, label = FALSE,  order = TRUE, cols = c("gray75", "royalblue4"),  pt.size = 0.1) 
VlnPlot(IVD.integrated3, split.by = "MajorType", features = c("MAP1B", "ACAN", "COL1A1", "CALR", "HSPA6", "MEG3", "HBB",  "CD24", "TBXT"), slot = "counts", log = TRUE, same.y.lims = TRUE, pt.size = 0.5, sort = FALSE, stack = TRUE, fill.by = "ident") + theme(aspect.ratio = 3)
Idents(IVD.integrated3) <- "ClusterNumber"
Idents(IVD.nIVD) <- "ClusterNumber"
Idents(IVD.aIVD) <- "ClusterNumber"

levels(x = IVD.integrated3) <- c("10", "12", '13', "5", "3", "4", "6", "7", "9", "8", "1", "0", "14", "2", "11")
levels(x = IVD.nIVD) <- c("10", "12", '13', "5", "3", "4", "6", "7", "9", "8", "1", "0", "14", "2", "11")
levels(x = IVD.aIVD) <- c("10", "12", '13', "5", "3", "4", "6", "7", "9", "8", "1", "0", "14", "2", "11")


p00b <- DotPlot(IVD.nIVD, features = c("MAP1B", "SOX4","SOX17", "ITGA6", "BASP1", 
                                       "ACAN", "COL2A1", "SOX9",
                                       "COL1A1", "CALR", "HSPA6",
                                       "MEG3", "CD44", "CD14", "HBB", "HBA1"), cols = c("royalblue4", "rosybrown2"), scale.by = "radius", dot.scale = 6, col.min = -2, col.max = 2) + 
  theme(axis.text.x = element_text(angle = 90)) 
p00c <- DotPlot(IVD.aIVD, features = c("MAP1B", "SOX4","SOX17", "ITGA6", "BASP1", 
                                       "ACAN", "COL2A1", "SOX9",
                                       "COL1A1", "CALR", "HSPA6",
                                       "MEG3", "CD44", "CD14", "HBB", "HBA1"), cols = c("royalblue4", "rosybrown2"), scale.by = "radius", dot.scale = 6, col.min = -2, col.max = 2) + 
  theme(axis.text.x = element_text(angle = 90)) 

p00b + p00c

Idents(IVD.integrated3) <- "Subtype"
Idents(IVD.nIVD) <- "Subtype"
Idents(IVD.aIVD) <- "Subtype"

levels(x = IVD.integrated3) <- c('RBC2', "RBC1", "IC", "non-MC", 'OAF3', "OAF2", "OAF1", "IAF1", "NP4", "NP3", "NP2", "NP1", "NCxNP", "NC2", "NC1")
levels(x = IVD.nIVD) <- c('RBC2', "RBC1", "IC", "non-MC", 'OAF3', "OAF2", "OAF1", "IAF1", "NP4", "NP3", "NP2", "NP1", "NCxNP", "NC2", "NC1")
levels(x = IVD.aIVD) <- c('RBC2', "RBC1", "IC", "non-MC", 'OAF3', "OAF2", "OAF1", "IAF1", "NP4", "NP3", "NP2", "NP1", "NCxNP", "NC2", "NC1")

P00 <- DimPlot(IVD.integrated3,  reduction = "umap", pt.size = 0.1, label.size = FALSE,
               cols = c("red", "red3", "gray35", "gray70", "darkorange3", "darkorange2", "darkorange1", "plum", "dodgerblue4", "dodgerblue3", "dodgerblue2", "dodgerblue1", "cyan3", "forestgreen", "darkolivegreen")) + theme(aspect.ratio = "1")


p01 <- FeaturePlot(IVD.integrated3, features = c("MAP1B", "ACAN", "COL1A1", "CALR", "HSPA6", "MEG3", "HBB"), 
            ncol = 3, label = FALSE,  order = TRUE, cols = c("gray75", "royalblue4"),  pt.size = 0.1) 
table(IVD.integrated3$MajorType)
table(IVD.integrated3$MajorType, IVD.integrated3$Subtype)


#Figure 2
IVD.nIVD <- RenameIdents(object = IVD.nIVD, "NC1" = "NC1", "NC2" = "NC2", "NP1" = "Other Cells", "NP2" = "Other Cells", "NP3" = "Other Cells", "NP4" = "Other Cells", "IAF1" = "Other Cells", "OAF1" = "Other Cells", "OAF2" = "Other Cells", "OAF3" = "Other Cells", "NCxNP" = "Other Cells", "RBC2" = "Other Cells", "RBC1" = "Other Cells", "IC" = "Other Cells", "non-MC" = "Other Cells")
IVD.aIVD <- RenameIdents(object = IVD.aIVD, "NC1" = "NC1", "NC2" = "NC2", "NP1" = "Other Cells", "NP2" = "Other Cells", "NP3" = "Other Cells", "NP4" = "Other Cells", "IAF1" = "Other Cells", "OAF1" = "Other Cells", "OAF2" = "Other Cells", "OAF3" = "Other Cells", "NCxNP" = "Other Cells", "RBC2" = "Other Cells", "RBC1" = "Other Cells", "IC" = "Other Cells", "non-MC" = "Other Cells")

levels(x = IVD.nIVD) <- c('Other Cells', 'NC2', "NC1")
levels(x = IVD.aIVD) <- c('Other Cells', 'NC2', "NC1")

p03a <- VlnPlot(IVD.nIVD, features = c("MAP1B", "SOX4", "SOX17", "ITGA6",  "CD44", "CD14"),
        slot = "counts", log = TRUE, same.y.lims = TRUE,
cols = c("azure4",  "mediumseagreen", "springgreen4"), idents = c("Other Cells", "NC2", "NC1"), pt.size = 0.5, sort = FALSE, stack = TRUE, fill.by = "ident") + theme(aspect.ratio = 1.5)

p03b <- VlnPlot(IVD.aIVD, features =  c("MAP1B", "SOX4", "SOX17", "ITGA6",  "CD44", "CD14"),
                slot = "counts", log = TRUE, same.y.lims = TRUE,
                cols = c("azure4",  "mediumseagreen", "springgreen4"), idents = c("Other Cells", "NC2", "NC1"), pt.size = 0.5, sort = FALSE, stack = TRUE, fill.by = "ident") + theme(aspect.ratio = 1.5)


P04a <- DimPlot(IVD.nIVD, reduction = "umap", pt.size = 0.1, label.size = FALSE, 
               cols = c("azure4",  "mediumseagreen", "springgreen4")) + theme(aspect.ratio = "1")

P04b <- DimPlot(IVD.aIVD, reduction = "umap", pt.size = 0.1, label.size = FALSE, 
                cols = c("azure4",  "mediumseagreen", "springgreen4")) + theme(aspect.ratio = "1")

Idents(IVD.nIVD) <- "Subtype"
Idents(IVD.aIVD) <- "Subtype"

p05a <- FeaturePlot(IVD.nIVD, features = c("MAP1B", "SOX4", "CD24", "FABP4", "AKR1C1", "APOE"), keep.scale = "all",
            ncol = 3, label = FALSE,  order = TRUE, cols = c("gray75", "royalblue4"),  pt.size = 0.8) 

p05b <- FeaturePlot(IVD.aIVD, features = c("MAP1B", "SOX4", "CD24", "FABP4", "AKR1C1", "APOE"), keep.scale = "all",
                    ncol = 3, label = FALSE,  order = TRUE, cols = c("gray75", "royalblue4"),  pt.size = 0.8) 

IVD.All.RBC2.markers <- FindMarkers(IVD.integrated3, ident.1 = "RBC2", min.pct = 0.25)
IVD.All.RBC1.markers <- FindMarkers(IVD.integrated3, ident.1 = "RBC1", min.pct = 0.25)
IVD.All.IC.markers <- FindMarkers(IVD.integrated3, ident.1 = "IC", min.pct = 0.25)
IVD.All.nonMC.markers <- FindMarkers(IVD.integrated3, ident.1 = "non-MC", min.pct = 0.25)
IVD.All.OAF3.markers <- FindMarkers(IVD.integrated3, ident.1 = "OAF3", min.pct = 0.25)
IVD.All.OAF2.markers <- FindMarkers(IVD.integrated3, ident.1 = "OAF2", min.pct = 0.25)
IVD.All.OAF1.markers <- FindMarkers(IVD.integrated3, ident.1 = "OAF1", min.pct = 0.25)
IVD.All.IAF1.markers <- FindMarkers(IVD.integrated3, ident.1 = "IAF1", min.pct = 0.25)
IVD.All.NP4.markers <- FindMarkers(IVD.integrated3, ident.1 = "NP4", min.pct = 0.25)
IVD.All.NP3.markers <- FindMarkers(IVD.integrated3, ident.1 = "NP3", min.pct = 0.25)
IVD.All.NP2.markers <- FindMarkers(IVD.integrated3, ident.1 = "NP2", min.pct = 0.25)
IVD.All.NP1.markers <- FindMarkers(IVD.integrated3, ident.1 = "NP1", min.pct = 0.25)
IVD.All.NCxNP.markers <- FindMarkers(IVD.integrated3, ident.1 = "NCxNP", min.pct = 0.25)
IVD.All.NC2.markers <- FindMarkers(IVD.integrated3, ident.1 = "NC2", min.pct = 0.25)
IVD.All.NC1.markers <- FindMarkers(IVD.integrated3, ident.1 = "NC1", min.pct = 0.25)


IVD.All.RBC2.top10 <- rownames(head(IVD.All.RBC2.markers, n = 10))
IVD.All.RBC1.top10 <- rownames(head(IVD.All.RBC1.markers, n = 10))
IVD.All.IC.top10 <- rownames(head(IVD.All.IC.markers, n = 10))
IVD.All.nonMC.top10 <- rownames(head(IVD.All.nonMC.markers, n = 10))
IVD.All.OAF3.top10 <- rownames(head(IVD.All.OAF3.markers, n = 10))
IVD.All.OAF2.top10 <- rownames(head(IVD.All.OAF2.markers, n = 10))
IVD.All.OAF1.top10 <- rownames(head(IVD.All.OAF1.markers, n = 10))
IVD.All.IAF1.top10 <- rownames(head(IVD.All.IAF1.markers, n = 10))
IVD.All.NP4.top10 <- rownames(head(IVD.All.NP4.markers, n = 10))
IVD.All.NP3.top10 <- rownames(head(IVD.All.NP3.markers, n = 10))
IVD.All.NP2.top10 <- rownames(head(IVD.All.NP2.markers, n = 10))
IVD.All.NP1.top10 <- rownames(head(IVD.All.NC1.markers, n = 10))
IVD.All.NCxNP.top10 <- rownames(head(IVD.All.NCxNP.markers, n = 10))
IVD.All.NC2.top10 <- rownames(head(IVD.All.NC2.markers, n = 10))
IVD.All.NC1.top10 <- rownames(head(IVD.All.NC1.markers, n = 10))

IVD.nIVD.RBC2.markers <- FindMarkers(IVD.nIVD, ident.1 = "RBC2", min.pct = 0.25)
IVD.nIVD.RBC1.markers <- FindMarkers(IVD.nIVD, ident.1 = "RBC1", min.pct = 0.25)
IVD.nIVD.IC.markers <- FindMarkers(IVD.nIVD, ident.1 = "IC", min.pct = 0.25)
IVD.nIVD.nonMC.markers <- FindMarkers(IVD.nIVD, ident.1 = "non-MC", min.pct = 0.25)
IVD.nIVD.OAF3.markers <- FindMarkers(IVD.nIVD, ident.1 = "OAF3", min.pct = 0.25)
IVD.nIVD.OAF2.markers <- FindMarkers(IVD.nIVD, ident.1 = "OAF2", min.pct = 0.25)
IVD.nIVD.OAF1.markers <- FindMarkers(IVD.nIVD, ident.1 = "OAF1", min.pct = 0.25)
IVD.nIVD.IAF1.markers <- FindMarkers(IVD.nIVD, ident.1 = "IAF1", min.pct = 0.25)
IVD.nIVD.NP4.markers <- FindMarkers(IVD.nIVD, ident.1 = "NP4", min.pct = 0.25)
IVD.nIVD.NP3.markers <- FindMarkers(IVD.nIVD, ident.1 = "NP3", min.pct = 0.25)
IVD.nIVD.NP2.markers <- FindMarkers(IVD.nIVD, ident.1 = "NP2", min.pct = 0.25)
IVD.nIVD.NP1.markers <- FindMarkers(IVD.nIVD, ident.1 = "NP1", min.pct = 0.25)
IVD.nIVD.NCxNP.markers <- FindMarkers(IVD.nIVD, ident.1 = "NCxNP", min.pct = 0.25)
IVD.nIVD.NC2.markers <- FindMarkers(IVD.nIVD, ident.1 = "NC2", min.pct = 0.25)
IVD.nIVD.NC1.markers <- FindMarkers(IVD.nIVD, ident.1 = "NC1", min.pct = 0.25)

IVD.aIVD.RBC2.markers <- FindMarkers(IVD.aIVD, ident.1 = "RBC2", min.pct = 0.25)
IVD.aIVD.RBC1.markers <- FindMarkers(IVD.aIVD, ident.1 = "RBC1", min.pct = 0.25)
IVD.aIVD.IC.markers <- FindMarkers(IVD.aIVD, ident.1 = "IC", min.pct = 0.25)
IVD.aIVD.nonMC.markers <- FindMarkers(IVD.aIVD, ident.1 = "non-MC", min.pct = 0.25)
IVD.aIVD.OAF3.markers <- FindMarkers(IVD.aIVD, ident.1 = "OAF3", min.pct = 0.25)
IVD.aIVD.OAF2.markers <- FindMarkers(IVD.aIVD, ident.1 = "OAF2", min.pct = 0.25)
IVD.aIVD.OAF1.markers <- FindMarkers(IVD.aIVD, ident.1 = "OAF1", min.pct = 0.25)
IVD.aIVD.IAF1.markers <- FindMarkers(IVD.aIVD, ident.1 = "IAF1", min.pct = 0.25)
IVD.aIVD.NP4.markers <- FindMarkers(IVD.aIVD, ident.1 = "NP4", min.pct = 0.25)
IVD.aIVD.NP3.markers <- FindMarkers(IVD.aIVD, ident.1 = "NP3", min.pct = 0.25)
IVD.aIVD.NP2.markers <- FindMarkers(IVD.aIVD, ident.1 = "NP2", min.pct = 0.25)
IVD.aIVD.NP1.markers <- FindMarkers(IVD.aIVD, ident.1 = "NP1", min.pct = 0.25)
IVD.aIVD.NCxNP.markers <- FindMarkers(IVD.aIVD, ident.1 = "NCxNP", min.pct = 0.25)
IVD.aIVD.NC2.markers <- FindMarkers(IVD.aIVD, ident.1 = "NC2", min.pct = 0.25)
IVD.aIVD.NC1.markers <- FindMarkers(IVD.aIVD, ident.1 = "NC1", min.pct = 0.25)

IVD.nIVD <- RenameIdents(object = IVD.nIVD, "NC1" = "NC", "NC2" = "NC", "NP1" = "Other Cells", "NP2" = "Other Cells", "NP3" = "Other Cells", "NP4" = "Other Cells", "IAF1" = "Other Cells", "OAF1" = "Other Cells", "OAF2" = "Other Cells", "OAF3" = "Other Cells", "NCxNP" = "Other Cells", "RBC2" = "Other Cells", "RBC1" = "Other Cells", "IC" = "Other Cells", "non-MC" = "Other Cells")
IVD.aIVD <- RenameIdents(object = IVD.aIVD, "NC1" = "NC", "NC2" = "NC", "NP1" = "Other Cells", "NP2" = "Other Cells", "NP3" = "Other Cells", "NP4" = "Other Cells", "IAF1" = "Other Cells", "OAF1" = "Other Cells", "OAF2" = "Other Cells", "OAF3" = "Other Cells", "NCxNP" = "Other Cells", "RBC2" = "Other Cells", "RBC1" = "Other Cells", "IC" = "Other Cells", "non-MC" = "Other Cells")
IVD.nIVD.NC.markers <- FindMarkers(IVD.nIVD, ident.1 = "NC", min.pct = 0.05)
IVD.aIVD.NC.markers <- FindMarkers(IVD.aIVD, ident.1 = "NC", min.pct = 0.05)

p06a1 <-  EnhancedVolcano(IVD.nIVD.NC.markers,
                lab = rownames(IVD.nIVD.NC.markers),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "NC (Neonatal IVD)",
                xlim = c(0,6),
                ylim = c(0,300),
                pCutoff = 10e-24,
                FCcutoff = 1,
                labSize = 4,
                pointSize = 2.5,
                col = c("grey30", "grey30", "grey30", "red2"),
                colAlpha = 3/5,
                legendPosition = 'none',
                drawConnectors = TRUE,
                widthConnectors = 1,    
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                encircle = c("FABP4", "COL4A1", "EDNRB", "TFPI"
                             ,"CD36", "MAP1B", "FTL", "IGFBP7",
                             "TMSB4X", "TMSB10", "AKR1C1", "CD24", "APOE"),
                encircleCol = 'black',
                encircleSize = 0.5,
                encircleFill = 'pink',
                encircleAlpha = 1/3,
                selectLab = c("FABP4", "COL4A1", "EDNRB", "TFPI"
                              ,"CD36", "MAP1B", "FTL", "IGFBP7",
                              "TMSB4X", "TMSB10", "AKR1C1", "CD24", "APOE")
) + theme(aspect.ratio = 1) + coord_flip()

p06a2 <-  EnhancedVolcano(IVD.aIVD.NC.markers,
                          lab = rownames(IVD.aIVD.NC.markers),
                          x = 'avg_log2FC',
                          y = 'p_val_adj',
                          title = "NC (Adult IVD)",
                          xlim = c(0,6),
                          pCutoff = 10e-24,
                          FCcutoff = 1,
                          labSize = 4,
                          pointSize = 2.5,
                          col = c("grey30", "grey30", "grey30", "red2"),
                          colAlpha = 3/5,
                          legendPosition = 'none',
                          drawConnectors = TRUE,
                          widthConnectors = 1,    
                          colConnectors = 'black',
                          gridlines.major = FALSE,
                          gridlines.minor = FALSE,
                          encircle = c("AKR1C1"),
                          encircleCol = 'black',
                          encircleSize = 0.5,
                          encircleFill = 'pink',
                          encircleAlpha = 1/3,
                          selectLab = c("FABP4",  "MAP1B",  "AKR1C1", "CD24", "APOE")
) + theme(aspect.ratio = 1) + coord_flip()

p06b1 <- EnhancedVolcano(IVD.nIVD.NC1.markers,
                         lab = rownames(IVD.nIVD.NC1.markers),
                         x = 'avg_log2FC',
                         y = 'p_val_adj',
                         title = "NC1 (Neonatal IVD)",
                         xlim = c(0,6),
                         ylim = c(0,350),
                         pCutoff = 10e-24,
                         FCcutoff = 1,
                         labSize = 4,
                         pointSize = 2.5,
                         col = c("grey30", "grey30", "grey30", "red2"),
                         colAlpha = 3/5,
                         legendPosition = 'none',
                         drawConnectors = TRUE,
                         widthConnectors = 1,    
                         colConnectors = 'black',
                         gridlines.major = FALSE,
                         gridlines.minor = FALSE,
                         encircle = c("FABP4", "STC1", "MYCT1", "CXCL8"
                                      ,"CALCRL", "EMCN", "EDNRB", "CD36", "GNG11",
                                      "TMSB4X", "TMSB10", "TM4SF1", "AKR1C1", "AKR1C1", "APOE"),
                         encircleCol = 'black',
                         encircleSize = 0.5,
                         encircleFill = 'pink',
                         encircleAlpha = 1/3,
                         selectLab = c("FABP4", "STC1", "MYCT1", "CXCL8"
                                       ,"CALCRL", "EMCN", "EDNRB", "CD36", "GNG11",
                                       "TMSB4X", "TMSB10", "TM4SF1", "MAP1B", "AKR1C1", "SOX17", "IGTA6", "AKR1C1", "APOE")
                           ) + theme(aspect.ratio = 1) + coord_flip()

p06b2 <- EnhancedVolcano(IVD.nIVD.NC2.markers,
                         lab = rownames(IVD.nIVD.NC2.markers),
                         x = 'avg_log2FC',
                         y = 'p_val_adj',
                         title = "NC2 (Neonatal IVD)",
                         xlim = c(0, 3.5),
                         pCutoff = 10e-24,
                         FCcutoff = 1,
                         labSize = 4,
                         pointSize = 2.5,
                         col = c("grey30", "grey30", "grey30", "red2"),
                         colAlpha = 3/5,
                         legendPosition = 'none',
                         drawConnectors = TRUE,
                         widthConnectors = 1,    
                         colConnectors = 'black',
                         gridlines.major = FALSE,
                         gridlines.minor = FALSE,
                         
                         encircle = c("AKR1C1", "COL4A1", "IGFBP7", "FTL"
                                      ,"TFPI", "KCNE4", "COL4A2", "OLFM2", "UCHL1",
                                      "MAP1B", "APOE", "FABP4", "CD24"),
                         encircleCol = 'black',
                         encircleSize = 0.5,
                         encircleFill = 'pink',
                         encircleAlpha = 1/3,
                         selectLab = c("AKR1C1", "COL4A1", "IGFBP7", "FTL"
                                       ,"TFPI", "KCNE4", "COL4A2", "OLFM2", "UCHL1",
                                       "MAP1B", "APOE", "FABP4", "CD24")) + theme(aspect.ratio = 1) + coord_flip()

p06c1 <- EnhancedVolcano(IVD.aIVD.NC1.markers,
                         lab = rownames(IVD.aIVD.NC1.markers),
                         x = 'avg_log2FC',
                         y = 'p_val_adj',
                         title = "NC1 (Adult IVD)",
                         xlim = c(0, 6),
                         ylim = c(0, 350),
                         pCutoff = 10e-24,
                         FCcutoff = 1,
                         labSize = 4,
                         pointSize = 2.5,
                         col = c("grey30", "grey30", "grey30", "red2"),
                         colAlpha = 3/5,
                         legendPosition = 'none',
                         drawConnectors = TRUE,
                         widthConnectors = 1,    
                         colConnectors = 'black',
                         gridlines.major = FALSE,
                         gridlines.minor = FALSE,
                         encircle = c( "BMX", "MYCT1", "CLDN14", "FABP4"
                                      ),
                         encircleCol = 'black',
                         encircleSize = 0.5,
                         encircleFill = 'pink',
                         encircleAlpha = 1/3,
                         selectLab = c("AKR1C1", "BMX", "MYCT1", "CLDN14",
                                       "MAP1B", "APOE", "FABP4", "CD24")
) + theme(aspect.ratio = 1) + coord_flip()

p06c2 <- EnhancedVolcano(IVD.aIVD.NC2.markers,
                         lab = rownames(IVD.aIVD.NC2.markers),
                         x = 'avg_log2FC',
                         y = 'p_val_adj',
                         title = "NC2 (Adult IVD)",
                         xlim = c(0, 2),
                         pCutoff = 10e-24,
                         FCcutoff = 1,
                         labSize = 4,
                         pointSize = 2.5,
                         col = c("grey30", "grey30", "grey30", "red2"),
                         colAlpha = 3/5,
                         legendPosition = 'none',
                         drawConnectors = TRUE,
                         widthConnectors = 1,    
                         colConnectors = 'black',
                         gridlines.major = FALSE,
                         gridlines.minor = FALSE,
                         encircle = c("AKR1C1"),
                         encircleCol = 'black',
                         encircleSize = 0.5,
                         encircleFill = 'pink',
                         encircleAlpha = 1/3,
                         selectLab = c("AKR1C1", 
                                       "MAP1B", "APOE", "FABP4", "CD24")) + theme(aspect.ratio = 1) + coord_flip()


#Figure 3

levels(x = IVD.integrated3) <- c('RBC2', "RBC1", "IC", "non-MC", 'OAF3', "OAF2", "OAF1", "IAF1", "NP4", "NP3", "NP2", "NP1", "NCxNP", "NC2", "NC1")
levels(x = IVD.nIVD) <- c('RBC2', "RBC1", "IC", "non-MC", 'OAF3', "OAF2", "OAF1", "IAF1", "NP4", "NP3", "NP2", "NP1", "NCxNP", "NC2", "NC1")
levels(x = IVD.aIVD) <- c('RBC2', "RBC1", "IC", "non-MC", 'OAF3', "OAF2", "OAF1", "IAF1", "NP4", "NP3", "NP2", "NP1", "NCxNP", "NC2", "NC1")

Idents(IVD.integrated3) <- "Subtype"
Idents(IVD.nIVD) <- "Subtype"
Idents(IVD.aIVD) <- "Subtype"

levels(x = IVD.integrated3) <- c('RBC2', "RBC1", "IC", "non-MC", 'OAF3', "OAF2", "OAF1", "IAF1", "NP4", "NP3", "NP2", "NP1", "NCxNP", "NC2", "NC1")
levels(x = IVD.nIVD) <- c('RBC2', "RBC1", "IC", "non-MC", 'OAF3', "OAF2", "OAF1", "IAF1", "NP4", "NP3", "NP2", "NP1", "NCxNP", "NC2", "NC1")
levels(x = IVD.aIVD) <- c('RBC2', "RBC1", "IC", "non-MC", 'OAF3', "OAF2", "OAF1", "IAF1", "NP4", "NP3", "NP2", "NP1", "NCxNP", "NC2", "NC1")

p07 <- VlnPlot(IVD.integrated3, features = c("FABP4", "AKR1C1", "APOE"), 
        slot = "counts", log = TRUE, same.y.lims = TRUE, split.by = "MajorType", 
        cols = c("lightpink3", "lightcyan4"), flip = FALSE,
        pt.size = 0.1, sort = FALSE, stack = TRUE, fill.by = "feature") + theme(aspect.ratio = 5)

# Figure 4

Idents(IVD.integrated3) <- "Subtype"
Idents(IVD.nIVD) <- "Subtype"
Idents(IVD.aIVD) <- "Subtype"

IVD.nIVD <- RenameIdents(object = IVD.nIVD, "NC1" = "Other Cells", "NC2" = "Other Cells", "NP1" = "NPC1", "NP2" = "NPC2", "NP3" = "NPC3", "NP4" = "NPC4", "IAF1" = "Other Cells", "OAF1" = "Other Cells", "OAF2" = "Other Cells", "OAF3" = "Other Cells", "NCxNP" = "Other Cells", "RBC2" = "Other Cells", "RBC1" = "Other Cells", "IC" = "Other Cells", "non-MC" = "Other Cells")
IVD.aIVD <- RenameIdents(object = IVD.aIVD, "NC1" = "Other Cells", "NC2" = "Other Cells", "NP1" = "NPC1", "NP2" = "NPC2", "NP3" = "NPC3", "NP4" = "NPC4", "IAF1" = "Other Cells", "OAF1" = "Other Cells", "OAF2" = "Other Cells", "OAF3" = "Other Cells", "NCxNP" = "Other Cells", "RBC2" = "Other Cells", "RBC1" = "Other Cells", "IC" = "Other Cells", "non-MC" = "Other Cells")
levels(x = IVD.nIVD) <- c('Other Cells', "NPC4", "NPC3", "NPC2", 'NPC1')
levels(x = IVD.aIVD) <- c('Other Cells', "NPC4", "NPC3", "NPC2", 'NPC1')

p08a <- VlnPlot(IVD.nIVD, features = c("ACAN","COL2A1", "SOX9", "OGN", "MMP13", "COL9A1", "DPT", "COL1A1"), 
        slot = "counts", log = TRUE, same.y.lims = TRUE,  cols = c("azure4", "dodgerblue4", "dodgerblue1", "dodgerblue3", "dodgerblue2"),
        pt.size = 0.2, sort = FALSE, stack = TRUE, fill.by = "ident") + theme(aspect.ratio = 2)

p08b <- VlnPlot(IVD.aIVD, features = c("ACAN","COL2A1", "SOX9", "S100A2", "TIMP3", "APOD", "S100A4", "COL1A1"), 
        slot = "counts", log = TRUE, same.y.lims = TRUE,  cols = c("azure4", "dodgerblue4", "dodgerblue1", "dodgerblue3", "dodgerblue2"),
        pt.size = 0.2, sort = FALSE, stack = TRUE, fill.by = "ident") + theme(aspect.ratio = 2)
p08c1 <- FeaturePlot(IVD.nIVD, features = "OGN", ncol = 1, label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 1, coord.fixed = TRUE) 
p08c2 <- FeaturePlot(IVD.nIVD, features = "MMP13", ncol = 1, label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 1, coord.fixed = TRUE) 
p08c3 <- FeaturePlot(IVD.nIVD, features = "COL9A1", ncol = 1, label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 1, coord.fixed = TRUE) 
p08c4 <- FeaturePlot(IVD.nIVD, features = "DPT", ncol = 1, label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 1, coord.fixed = TRUE) 
p08d1 <- FeaturePlot(IVD.aIVD, features = "S100A2", ncol = 1, label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 1, coord.fixed = TRUE) 
p08d2 <- FeaturePlot(IVD.aIVD, features = "TIMP3", ncol = 1, label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 1, coord.fixed = TRUE) 
p08d3 <- FeaturePlot(IVD.aIVD, features = "APOD", ncol = 1, label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 1, coord.fixed = TRUE) 
p08d4 <- FeaturePlot(IVD.aIVD, features = "S100A4", ncol = 1, label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 1, coord.fixed = TRUE) 

p09a <- DimPlot(IVD.nIVD, reduction = "umap", pt.size = 0.2, label.size = FALSE, 
        cols = c("azure4","dodgerblue3","dodgerblue2","dodgerblue1","dodgerblue4")
) + theme(aspect.ratio = "1")

p09b <- DimPlot(IVD.aIVD, reduction = "umap", pt.size = 0.2, label.size = FALSE, 
        cols = c("azure4","dodgerblue3","dodgerblue2","dodgerblue1","dodgerblue4")
) + theme(aspect.ratio = "1")


# Figure 5

IVD.nIVD <- RenameIdents(object = IVD.nIVD, "NC1" = "Other Cells", "NC2" = "Other Cells", "NP1" = "Other Cells", "NP2" = "Other Cells", "NP3" = "Other Cells", "NP4" = "Other Cells", "IAF1" = "iAFC1", "OAF1" = "oAFC1", "OAF2" = "oAFC2", "OAF3" = "oAFC3", "NCxNP" = "Other Cells", "RBC2" = "Other Cells", "RBC1" = "Other Cells", "IC" = "Other Cells", "non-MC" = "Other Cells")
IVD.aIVD <- RenameIdents(object = IVD.aIVD, "NC1" = "Other Cells", "NC2" = "Other Cells", "NP1" = "Other Cells", "NP2" = "Other Cells", "NP3" = "Other Cells", "NP4" = "Other Cells", "IAF1" = "iAFC1", "OAF1" = "oAFC1", "OAF2" = "oAFC2", "OAF3" = "oAFC3", "NCxNP" = "Other Cells", "RBC2" = "Other Cells", "RBC1" = "Other Cells", "IC" = "Other Cells", "non-MC" = "Other Cells")
levels(x = IVD.nIVD) <- c('Other Cells', "iAFC1", "oAFC1", "oAFC2", 'oAFC3')
levels(x = IVD.aIVD) <- c('Other Cells', "iAFC1", "oAFC1", "oAFC2", 'oAFC3')

p10a <-  VlnPlot(IVD.nIVD, features = c("COL1A1","CALR", "HSPA6", "ACAN", "COL2A1", "SOX9", "MAP1B"), 
                 slot = "counts", log = TRUE, same.y.lims = TRUE,  cols = c("azure4", "plum", "darkorange1", "darkorange2", "darkorange3"),
                 pt.size = 0.2, sort = FALSE, stack = TRUE, fill.by = "ident") + theme(aspect.ratio = 1.6)

p10b <-   VlnPlot(IVD.aIVD, features = c("COL1A1","CALR", "HSPA6", "ACAN", "COL2A1", "SOX9", "MAP1B"), 
                  slot = "counts", log = TRUE, same.y.lims = TRUE,  cols = c("azure4", "plum", "darkorange1", "darkorange2", "darkorange3"),
                  pt.size = 0.2, sort = FALSE, stack = TRUE, fill.by = "ident") + theme(aspect.ratio = 1.6)


p11a1 <- DimPlot(IVD.nIVD, reduction = "umap", pt.size = 0.2, label.size = FALSE, 
        cols = c("azure4", "plum", "darkorange1", "darkorange2", "darkorange3")) + theme(aspect.ratio = "1")

p11a2 <- DimPlot(IVD.aIVD, reduction = "umap", pt.size = 0.2, label.size = FALSE, 
                 cols = c("azure4", "plum", "darkorange1", "darkorange2", "darkorange3")) + theme(aspect.ratio = "1")

p11c1 <- FeaturePlot(IVD.nIVD, features = "LGALS1", ncol = 1, label = FALSE, order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
p11c2 <- FeaturePlot(IVD.nIVD, features = "HES1", ncol = 1, label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
p11c3 <- FeaturePlot(IVD.nIVD, features = "HERPUD1", ncol = 1, label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
p11c4 <- FeaturePlot(IVD.nIVD, features = "DNAJB1", ncol = 1, label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 

p11d1 <- FeaturePlot(IVD.aIVD, features = "ASPN", ncol = 1, label = FALSE, order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
p11d2 <- FeaturePlot(IVD.aIVD, features = "MANF", ncol = 1, label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
p11d3 <- FeaturePlot(IVD.aIVD, features = "HERPUD1", ncol = 1, label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 
p11d4 <- FeaturePlot(IVD.aIVD, features = "DNAJB1", ncol = 1, label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"),  pt.size = 0.4, coord.fixed = TRUE) 

p11e1 <- VlnPlot(IVD.nIVD, features = c("LGALS1", "HES1", "HERPUD1", "DNAJB1"), idents = c( "iAFC1", "oAFC1", "oAFC2", "oAFC3"), slot = "counts", log = TRUE, cols = c("plum2","darkorange1", "darkorange2", "darkorange3") , pt.size = 0, sort = FALSE, ncol=2)
p11e2 <-VlnPlot(IVD.aIVD, features = c("ASPN", "MANF", "HERPUD1", "DNAJB1"), idents = c( "iAFC1", "oAFC1", "oAFC2", "oAFC3"), slot = "counts", log = TRUE, cols = c("plum2","darkorange1", "darkorange2", "darkorange3") , pt.size = 0, sort = FALSE, ncol=2)

p12a <-   VlnPlot(IVD.nIVD, idents = c("iAFC1", "oAFC1", "oAFC2", "oAFC3"), 
                  features =c("COL12A1", "COL6A3", "COL5A2", "COL5A1", "COL3A1", "COL1A1"),
                  cols = c("plum", "darkorange1", "darkorange2", "darkorange3"), fill.by = "ident", same.y.lims = TRUE,
                  slot = "counts", log = TRUE,  pt.size = 0.01, sort = FALSE, stack = TRUE, ncol = 1)

p12b <-   VlnPlot(IVD.aIVD, idents = c("iAFC1", "oAFC1", "oAFC2", "oAFC3"), 
                  features =c("COL12A1", "COL6A3", "COL5A2", "COL5A1", "COL3A1", "COL1A1"),
                  cols = c("plum", "darkorange1", "darkorange2", "darkorange3"), fill.by = "ident", same.y.lims = TRUE,
                  slot = "counts", log = TRUE,  pt.size = 0.01, sort = FALSE, stack = TRUE, ncol = 1)


#Figure 6

p13a <- DotPlot(IVD.nIVD, features = c("MAP1B", "SOX4","SOX17", "ITGA6", "BASP1", "FABP4", "AKR1C1", "APOE", 
                                       "ACAN", "COL2A1", "SOX9",
                                       "COL1A1", "CALR", "HSPA6", "LGALS1", "HES1", "HERPUD1", "DNAJB1",  
                                       "MEG3", "CD44", "CD14", "HBB", "HBA1"), cols = c("royalblue4", "rosybrown2"), scale.by = "radius", dot.scale = 6, col.min = -2, col.max = 2) + 
  theme(axis.text.x = element_text(angle = 90)) 

p13b <- DotPlot(IVD.aIVD, features = c("MAP1B", "SOX4","SOX17", "ITGA6", "BASP1", "FABP4", "AKR1C1", "APOE", 
                                       "ACAN", "COL2A1", "SOX9",
                                       "COL1A1", "CALR", "HSPA6", "ASPN", "MANF", "HERPUD1", "DNAJB1",  
                                       "MEG3", "CD44", "CD14", "HBB", "HBA1"), cols = c("royalblue4", "rosybrown2"), scale.by = "radius", dot.scale = 6, col.min = -2, col.max = 2) + 
  theme(axis.text.x = element_text(angle = 90)) 

IVD.nIVD@active.ident

IVD.nIVD.3DUMAP <- IVD.nIVD
IVD.nIVD.3DUMAP <- RunUMAP(IVD.nIVD.3DUMAP,
                           dims = 1:30,
                           n.components = 3L)
umap_1 <- IVD.nIVD.3DUMAP[["umap"]]@cell.embeddings[,1]
umap_2 <- IVD.nIVD.3DUMAP[["umap"]]@cell.embeddings[,2]
umap_3 <- IVD.nIVD.3DUMAP[["umap"]]@cell.embeddings[,3]
Embeddings(object = IVD.nIVD.3DUMAP, reduction = "umap")
plot.data <- FetchData(object = IVD.nIVD.3DUMAP, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters"))
plot.data$label <- paste(rownames(plot.data))
IVD.nIVD.3DUMAP$ClusterNumber
IVD.nIVD.3DUMAP@active.ident
IVD.nIVD.3DUMAP$seurat_clusters
plot.data@active.ident


p14a <- plot_ly(data = plot.data, 
                   x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                   color = ~seurat_clusters, 
                   colors = c( "dodgerblue1",
                               "dodgerblue2",
                               "darkolivegreen" ,
                               "darkorange3" ,
                               "darkorange2" ,
                               "gray35",
                                "darkorange1",
                                "plum",
                                "dodgerblue3",
                                "dodgerblue4",
                                "red3",
                               "forestgreen",
                                "red",
                               "gray70",
                               "cyan3"),
                   type = "scatter3d", 
                   mode = "markers", 
                   marker = list(size = 0.8, width=0.8), 
                   text=~label, 
                   hoverinfo="seurat_clusters") 


IVD.nIVD.3DUMAP <- IVD.aIVD
IVD.nIVD.3DUMAP <- RunUMAP(IVD.nIVD.3DUMAP,
                           dims = 1:30,
                           n.components = 3L)
umap_1 <- IVD.nIVD.3DUMAP[["umap"]]@cell.embeddings[,1]
umap_2 <- IVD.nIVD.3DUMAP[["umap"]]@cell.embeddings[,2]
umap_3 <- IVD.nIVD.3DUMAP[["umap"]]@cell.embeddings[,3]
Embeddings(object = IVD.nIVD.3DUMAP, reduction = "umap")
plot.data <- FetchData(object = IVD.nIVD.3DUMAP, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters"))
plot.data$label <- paste(rownames(plot.data))
IVD.nIVD.3DUMAP$ClusterNumber
IVD.nIVD.3DUMAP@active.ident
IVD.nIVD.3DUMAP$seurat_clusters
plot.data@active.ident


p14b <- plot_ly(data = plot.data, 
                x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                color = ~seurat_clusters, 
                colors = c( "dodgerblue1",
                            "dodgerblue2",
                            "darkolivegreen" ,
                            "darkorange3" ,
                            "darkorange2" ,
                            "gray35",
                            "darkorange1",
                            "plum",
                            "dodgerblue3",
                            "dodgerblue4",
                            "red3",
                            "forestgreen",
                            "red",
                            "gray70",
                            "cyan3"),
                type = "scatter3d", 
                mode = "markers", 
                marker = list(size = 0.8, width=0.8), 
                text=~label, 
                hoverinfo="seurat_clusters") 


#Export for pathway analysis in Qiagen IPA

Idents(IVD.integrated3) <- "Subtype"
Idents(IVD.nIVD) <- "Subtype"
Idents(IVD.aIVD) <- "Subtype"

MarkersToExcel.nIVD.NC1 <- cbind("Gene ID"=rownames(IVD.nIVD.NC1.markers), IVD.nIVD.NC1.markers)
MarkersToExcel.nIVD.NC2 <- cbind("Gene ID"=rownames(IVD.nIVD.NC2.markers), IVD.nIVD.NC2.markers)
MarkersToExcel.nIVD.NP1 <- cbind("Gene ID"=rownames(IVD.nIVD.NP1.markers), IVD.nIVD.NP1.markers)
MarkersToExcel.nIVD.NP2 <- cbind("Gene ID"=rownames(IVD.nIVD.NP2.markers), IVD.nIVD.NP2.markers)
MarkersToExcel.nIVD.NP3 <- cbind("Gene ID"=rownames(IVD.nIVD.NP3.markers), IVD.nIVD.NP3.markers)
MarkersToExcel.nIVD.NP4 <- cbind("Gene ID"=rownames(IVD.nIVD.NP4.markers), IVD.nIVD.NP4.markers)
MarkersToExcel.nIVD.IAF1 <- cbind("Gene ID"=rownames(IVD.nIVD.IAF1.markers), IVD.nIVD.IAF1.markers)
MarkersToExcel.nIVD.OAF1 <- cbind("Gene ID"=rownames(IVD.nIVD.OAF1.markers), IVD.nIVD.OAF1.markers)
MarkersToExcel.nIVD.OAF2 <- cbind("Gene ID"=rownames(IVD.nIVD.OAF2.markers), IVD.nIVD.OAF2.markers)
MarkersToExcel.nIVD.OAF3 <- cbind("Gene ID"=rownames(IVD.nIVD.OAF3.markers), IVD.nIVD.OAF3.markers)

MarkersToExcel.aIVD.NC1 <- cbind("Gene ID"=rownames(IVD.aIVD.NC1.markers), IVD.aIVD.NC1.markers)
MarkersToExcel.aIVD.NC2 <- cbind("Gene ID"=rownames(IVD.aIVD.NC2.markers), IVD.aIVD.NC2.markers)
MarkersToExcel.aIVD.NP1 <- cbind("Gene ID"=rownames(IVD.aIVD.NP1.markers), IVD.aIVD.NP1.markers)
MarkersToExcel.aIVD.NP2 <- cbind("Gene ID"=rownames(IVD.aIVD.NP2.markers), IVD.aIVD.NP2.markers)
MarkersToExcel.aIVD.NP3 <- cbind("Gene ID"=rownames(IVD.aIVD.NP3.markers), IVD.aIVD.NP3.markers)
MarkersToExcel.aIVD.NP4 <- cbind("Gene ID"=rownames(IVD.aIVD.NP4.markers), IVD.aIVD.NP4.markers)
MarkersToExcel.aIVD.IAF1 <- cbind("Gene ID"=rownames(IVD.aIVD.IAF1.markers), IVD.aIVD.IAF1.markers)
MarkersToExcel.aIVD.OAF1 <- cbind("Gene ID"=rownames(IVD.aIVD.OAF1.markers), IVD.aIVD.OAF1.markers)
MarkersToExcel.aIVD.OAF2 <- cbind("Gene ID"=rownames(IVD.aIVD.OAF2.markers), IVD.aIVD.OAF2.markers)
MarkersToExcel.aIVD.OAF3 <- cbind("Gene ID"=rownames(IVD.aIVD.OAF3.markers), IVD.aIVD.OAF3.markers)

write_xlsx(MarkersToExcel.nIVD.NC1,  path = "Outputs_v2/For Pathway Analysis/Marker list for nIVD NC1.xlsx")
write_xlsx(MarkersToExcel.nIVD.NC2,  path = "Outputs_v2/For Pathway Analysis/Marker list for nIVD NC2.xlsx")
write_xlsx(MarkersToExcel.nIVD.NP1,  path = "Outputs_v2/For Pathway Analysis/Marker list for nIVD NPC1.xlsx")
write_xlsx(MarkersToExcel.nIVD.NP2,  path = "Outputs_v2/For Pathway Analysis/Marker list for nIVD NPC2.xlsx")
write_xlsx(MarkersToExcel.nIVD.NP3,  path = "Outputs_v2/For Pathway Analysis/Marker list for nIVD NPC3.xlsx")
write_xlsx(MarkersToExcel.nIVD.NP4,  path = "Outputs_v2/For Pathway Analysis/Marker list for nIVD NPC4.xlsx")
write_xlsx(MarkersToExcel.nIVD.IAF1,  path = "Outputs_v2/For Pathway Analysis/Marker list for nIVD IAFC1.xlsx")
write_xlsx(MarkersToExcel.nIVD.OAF1,  path = "Outputs_v2/For Pathway Analysis/Marker list for nIVD OAFC1.xlsx")
write_xlsx(MarkersToExcel.nIVD.OAF2,  path = "Outputs_v2/For Pathway Analysis/Marker list for nIVD OAFC2.xlsx")
write_xlsx(MarkersToExcel.nIVD.OAF3,  path = "Outputs_v2/For Pathway Analysis/Marker list for nIVD OAFC3.xlsx")

write_xlsx(MarkersToExcel.aIVD.NC1,  path = "Outputs_v2/For Pathway Analysis/Marker list for aIVD NC1.xlsx")
write_xlsx(MarkersToExcel.aIVD.NC2,  path = "Outputs_v2/For Pathway Analysis/Marker list for aIVD NC2.xlsx")
write_xlsx(MarkersToExcel.aIVD.NP1,  path = "Outputs_v2/For Pathway Analysis/Marker list for aIVD NPC1.xlsx")
write_xlsx(MarkersToExcel.aIVD.NP2,  path = "Outputs_v2/For Pathway Analysis/Marker list for aIVD NPC2.xlsx")
write_xlsx(MarkersToExcel.aIVD.NP3,  path = "Outputs_v2/For Pathway Analysis/Marker list for aIVD NPC3.xlsx")
write_xlsx(MarkersToExcel.aIVD.NP4,  path = "Outputs_v2/For Pathway Analysis/Marker list for aIVD NPC4.xlsx")
write_xlsx(MarkersToExcel.aIVD.IAF1,  path = "Outputs_v2/For Pathway Analysis/Marker list for aIVD IAFC1.xlsx")
write_xlsx(MarkersToExcel.aIVD.OAF1,  path = "Outputs_v2/For Pathway Analysis/Marker list for aIVD OAFC1.xlsx")
write_xlsx(MarkersToExcel.aIVD.OAF2,  path = "Outputs_v2/For Pathway Analysis/Marker list for aIVD OAFC2.xlsx")
write_xlsx(MarkersToExcel.aIVD.OAF3,  path = "Outputs_v2/For Pathway Analysis/Marker list for aIVD OAFC3.xlsx")

IVD.nIVD <- RenameIdents(object = IVD.nIVD, "NC1" = "NC", "NC2" = "NC", "NP1" = "NPC", "NP2" = "NPC", "NP3" = "NPC", "NP4" = "NPC", "IAF1" = "Other Cells", "OAF1" = "Other Cells", "OAF2" = "Other Cells", "OAF3" = "Other Cells", "NCxNP" = "Other Cells", "RBC2" = "Other Cells", "RBC1" = "Other Cells", "IC" = "Other Cells", "non-MC" = "Other Cells")
IVD.aIVD <- RenameIdents(object = IVD.aIVD, "NC1" = "NC", "NC2" = "NC", "NP1" = "NPC", "NP2" = "NPC", "NP3" = "NPC", "NP4" = "NPC", "IAF1" = "Other Cells", "OAF1" = "Other Cells", "OAF2" = "Other Cells", "OAF3" = "Other Cells", "NCxNP" = "Other Cells", "RBC2" = "Other Cells", "RBC1" = "Other Cells", "IC" = "Other Cells", "non-MC" = "Other Cells")

IVD.nIVD.NC.markers <- FindMarkers(IVD.nIVD, ident.1 = "NC", min.pct = 0.05)
IVD.aIVD.NC.markers <- FindMarkers(IVD.aIVD, ident.1 = "NC", min.pct = 0.05)
IVD.nIVD.NPC.markers <- FindMarkers(IVD.nIVD, ident.1 = "NPC", min.pct = 0.05)
IVD.aIVD.NPC.markers <- FindMarkers(IVD.aIVD, ident.1 = "NPC", min.pct = 0.05)

MarkersToExcel.nIVD.NC <- cbind("Gene ID"=rownames(IVD.nIVD.NC.markers), IVD.nIVD.NC.markers)
MarkersToExcel.aIVD.NC <- cbind("Gene ID"=rownames(IVD.aIVD.NC.markers), IVD.aIVD.NC.markers)
MarkersToExcel.nIVD.NPC <- cbind("Gene ID"=rownames(IVD.nIVD.NPC.markers), IVD.nIVD.NPC.markers)
MarkersToExcel.aIVD.NPC <- cbind("Gene ID"=rownames(IVD.aIVD.NPC.markers), IVD.aIVD.NPC.markers)

write_xlsx(MarkersToExcel.nIVD.NC,  path = "Outputs_v2/For Pathway Analysis/Marker list for nIVD NC.xlsx")
write_xlsx(MarkersToExcel.aIVD.NC,  path = "Outputs_v2/For Pathway Analysis/Marker list for aIVD NC.xlsx")
write_xlsx(MarkersToExcel.nIVD.NPC,  path = "Outputs_v2/For Pathway Analysis/Marker list for nIVD NPC.xlsx")
write_xlsx(MarkersToExcel.aIVD.NPC,  path = "Outputs_v2/For Pathway Analysis/Marker list for aIVD NPC.xlsx")

#SI Heatmap
DefaultAssay(IVD.integrated3) <- "integrated"
levels(x = IVD.integrated3) <- c('NC1', "NC2", "NCxNP", "NP1", 'NP2', "NP3", "NP4", "IAF1", "OAF1", "OAF2", "OAF3", "non-MC", "IC", "RBC2", "RBC1")

p15 <- DoHeatmap(IVD.integrated3, features = c(IVD.All.NC1.top10, IVD.All.NC2.top10 , IVD.All.NCxNP.top10, 
                                        IVD.All.NP1.top10,  IVD.All.NP2.top10, IVD.All.NP3.top10, IVD.All.NP4.top10,
                                        IVD.All.IAF1.top10, IVD.All.OAF1.top10, IVD.All.OAF2.top10,  IVD.All.OAF3.top10,
                                        IVD.All.nonMC.top10,  IVD.All.IC.top10, IVD.All.RBC2.top10 , IVD.All.RBC1.top10) , group.by = "ident", label = FALSE) 

DefaultAssay(IVD.integrated3) <- "RNA"

#----------End Here------------
