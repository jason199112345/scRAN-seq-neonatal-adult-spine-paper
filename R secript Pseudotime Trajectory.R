#Six samples were created as Seurat objects by code shown in "R script major.R".
#The six seurat objects are: nIVD1, nIVD2, nIVD3, aIVD1, aIVD2, aIVD3.

DataAll <- c(nIVD1, nIVD2, nIVD3, aIVD1, aIVD2, aIVD3)


for (i in 1:length(DataAll)) {
  DataAll[[i]] <- NormalizeData(DataAll[[i]], verbose = FALSE)
  DataAll[[i]] <- FindVariableFeatures(DataAll[[i]], selection.method = "vst", 
                                       nfeatures = 6000, verbose = FALSE)
}

IVD.anchors <- FindIntegrationAnchors(object.list = DataAll, anchor.features = 6000)
IVD.integrated <- IntegrateData(anchorset = IVD.anchors)
DefaultAssay(IVD.integrated) <- "integrated"

IVD.integrated <- ScaleData(IVD.integrated, verbose = FALSE)
IVD.integrated <- RunPCA(IVD.integrated, npcs = 30, verbose = FALSE)
IVD.integrated <- RunUMAP(IVD.integrated, reduction = "pca", dims = 1:30)
DimPlot(IVD.integrated, reduction = "umap", pt.size = 2, label.size = 20) 

IVD.integrated3 <- FindNeighbors(IVD.integrated, reduction = "pca", dims = 1:30)
IVD.integrated3 <- FindClusters(IVD.integrated3, resolution = 0.5)
DefaultAssay(IVD.integrated3) <- "RNA"
DimPlot(IVD.integrated3, reduction = "umap", pt.size = 0.3, label = TRUE, label.size = 10) 


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
        cols = c("red", "red3", "gray35", "gray70", "darkorange3", "darkorange2", "darkorange1", "plum", "dodgerblue4", "dodgerblue3", "dodgerblue2", "dodgerblue1", "cyan3", "forestgreen", "darkolivegreen"))


#Assigning cluster identity & Figure 1

IVD.integrated3 <- RenameIdents(object = IVD.integrated3, "0" = "NP1", "1" = "NP2", "2" = "NC2", "3" = "OAF3", "4" = "OAF2", "5" = "non-MC", "6" =  "OAF1", "7" = "IAF1", "8" = "NP3", "9" = "NP4", "10" = "RBC2", "11" = "NC1", "12" = "RBC1", "13" = "IC", "14" = "NCxNP")
IVD.nIVD <- RenameIdents(object = IVD.nIVD, "0" = "NP1", "1" = "NP2", "2" = "NC2", "3" = "OAF3", "4" = "OAF2", "5" = "non-MC", "6" =  "OAF1", "7" = "IAF1", "8" = "NP3", "9" = "NP4", "10" = "RBC2", "11" = "NC1", "12" = "RBC1", "13" = "IC", "14" = "NCxNP")
IVD.aIVD <- RenameIdents(object = IVD.aIVD, "0" = "NP1", "1" = "NP2", "2" = "NC2", "3" = "OAF3", "4" = "OAF2", "5" = "non-MC", "6" =  "OAF1", "7" = "IAF1", "8" = "NP3", "9" = "NP4", "10" = "RBC2", "11" = "NC1", "12" = "RBC1", "13" = "IC", "14" = "NCxNP")
IVD.integrated3$Subtype <- Idents(IVD.integrated3)
IVD.nIVD$Subtype <- Idents(IVD.nIVD)
IVD.aIVD$Subtype <- Idents(IVD.aIVD)

levels(x = IVD.integrated3) <- c('RBC2', "RBC1", "IC", "non-MC", 'OAF3', "OAF2", "OAF1", "IAF1", "NP4", "NP3", "NP2", "NP1", "NCxNP", "NC2", "NC1")

#===============Monocle Trajectory==============

#Starting HERE
library(monocle)
library(reshape2)

#Normal start package
library(Seurat)
library(SeuratData)
library(cowplot)
library(patchwork)
library(ggplot2)
library(dplyr)
library("BiocFileCache")
require(gridExtra)
library(viridis)
library("writexl")



#=======================Pre-processing and loading data================================



IVD.integrated3 <- RenameIdents(object = IVD.integrated3,  
                                "NP1" = "NP", "NP2" = "NP", 'NP3' = 'NP', 'NP4' = 'NP', "IAF1" = "IAF", "OAF1" = "OAF", "OAF2" = "OAF", "OAF3" = "OAF", "NCxNP" = "NC/NP", "NC1" = "NC", "NC2" = "NC"
)

IVD.integrated3$Subtype <- Idents(IVD.integrated3)

#Loading data

IVD.all.newtraj <- subset(x = IVD.integrated3, idents = c("NC", "NP", "IAF", "OAF", "NC/NP"))

data <- as(as.matrix(IVD.all.newtraj@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = IVD.all.newtraj@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monocle_cds <- newCellDataSet(data, 
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

#Filter low-quality cells
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
print(head(fData(monocle_cds)))
expressed_genes <- row.names(subset(fData(monocle_cds),
                                    num_cells_expressed >= 10))
print(head(pData(monocle_cds)))



#====Universal Gene Marker============ (can be used after creating monocle_cds)
COL2A1_id <- row.names(subset(fData(monocle_cds), gene_short_name == "COL2A1"))
COL1A1_id <- row.names(subset(fData(monocle_cds), gene_short_name == "COL1A1"))
ACAN_id <- row.names(subset(fData(monocle_cds), gene_short_name == "ACAN"))
CD24_id <- row.names(subset(fData(monocle_cds), gene_short_name == "CD24"))
SCX_id <- row.names(subset(fData(monocle_cds), gene_short_name == "SCX"))
SOX9_id <- row.names(subset(fData(monocle_cds), gene_short_name == "SOX9"))

CXCL2_id <- row.names(subset(fData(monocle_cds), gene_short_name == "CXCL2"))
MAP1B_id <- row.names(subset(fData(monocle_cds), gene_short_name == "MAP1B"))
KRT19_id <- row.names(subset(fData(monocle_cds), gene_short_name == "KRT19"))
ITGA6_id <- row.names(subset(fData(monocle_cds), gene_short_name == "ITGA6"))
SOX17_id <- row.names(subset(fData(monocle_cds), gene_short_name == "SOX17"))
CYTL1_id <- row.names(subset(fData(monocle_cds), gene_short_name == "CYTL1"))
COL4A1_id <- row.names(subset(fData(monocle_cds), gene_short_name == "COL4A1"))
MYCT1_id <- row.names(subset(fData(monocle_cds), gene_short_name == "MYCT1"))
IL17B_id <- row.names(subset(fData(monocle_cds), gene_short_name == "IL17B"))
BIRC3_id <- row.names(subset(fData(monocle_cds), gene_short_name == "BIRC3"))

CXCL14_id <- row.names(subset(fData(monocle_cds), gene_short_name == "CXCL14"))
S100A2_id <- row.names(subset(fData(monocle_cds), gene_short_name == "S100A2"))
FABP4_id <- row.names(subset(fData(monocle_cds), gene_short_name == "FABP4"))
FABP5_id <- row.names(subset(fData(monocle_cds), gene_short_name == "FABP5"))
S100A4_id <- row.names(subset(fData(monocle_cds), gene_short_name == "S100A4"))

MAP1B_id <- row.names(subset(fData(monocle_cds), gene_short_name == "MAP1B"))
CD44_id <- row.names(subset(fData(monocle_cds), gene_short_name == "CD44"))
CALR_id <- row.names(subset(fData(monocle_cds), gene_short_name == "CALR"))
HSPB1_id <- row.names(subset(fData(monocle_cds), gene_short_name == "HSPB1"))
HSPA6_id <- row.names(subset(fData(monocle_cds), gene_short_name == "HSPA6"))
MEG3_id <- row.names(subset(fData(monocle_cds), gene_short_name == "MEG3"))
HBB_id <- row.names(subset(fData(monocle_cds), gene_short_name == "HBB"))
HBA1_id <- row.names(subset(fData(monocle_cds), gene_short_name == "HBA1"))
TNFRSF11B_id <- row.names(subset(fData(monocle_cds), gene_short_name == "TNFRSF11B"))
BASP1_id <- row.names(subset(fData(monocle_cds), gene_short_name == "BASP1"))
HSPA1A_id <- row.names(subset(fData(monocle_cds), gene_short_name == "HSPA1A"))
OGN_id <- row.names(subset(fData(monocle_cds), gene_short_name == "OGN"))
FGFBP2_id <- row.names(subset(fData(monocle_cds), gene_short_name == "FGFBP2"))
ANXA1_id <- row.names(subset(fData(monocle_cds), gene_short_name == "ANXA1"))
GADD45G_id <- row.names(subset(fData(monocle_cds), gene_short_name == "GADD45G"))
S100A6_id <- row.names(subset(fData(monocle_cds), gene_short_name == "S100A6"))
FTH1_id <- row.names(subset(fData(monocle_cds), gene_short_name == "FTH1"))
ZNF385D_id <- row.names(subset(fData(monocle_cds), gene_short_name == "ZNF385D"))
HSP90AA1_id <- row.names(subset(fData(monocle_cds), gene_short_name == "HSP90AA1"))
MEG3_id <- row.names(subset(fData(monocle_cds), gene_short_name == "MEG3"))
NEAT1_id <- row.names(subset(fData(monocle_cds), gene_short_name == "NEAT1"))
COL12A1_id <- row.names(subset(fData(monocle_cds), gene_short_name == "COL12A1"))
TFPI_id <- row.names(subset(fData(monocle_cds), gene_short_name == "TFPI"))
LGALS1_id <- row.names(subset(fData(monocle_cds), gene_short_name == "LGALS1"))
MRGPRF_id <- row.names(subset(fData(monocle_cds), gene_short_name == "MRGPRF"))
COL3A1_id <- row.names(subset(fData(monocle_cds), gene_short_name == "COL3A1"))
EVI2B_id <- row.names(subset(fData(monocle_cds), gene_short_name == "EVI2B"))
IGFBP5_id <- row.names(subset(fData(monocle_cds), gene_short_name == "IGFBP5"))
HES1_id <- row.names(subset(fData(monocle_cds), gene_short_name == "HES1"))
COL6A3_id <- row.names(subset(fData(monocle_cds), gene_short_name == "COL6A3"))
AKR1C1_id <- row.names(subset(fData(monocle_cds), gene_short_name == "AKR1C1"))
APOE_id <- row.names(subset(fData(monocle_cds), gene_short_name == "APOE"))
SOX4_id <- row.names(subset(fData(monocle_cds), gene_short_name == "SOX4"))
FABP4_id <- row.names(subset(fData(monocle_cds), gene_short_name == "FABP4"))



#==================Analysis ====================



pData(monocle_cds)$Total_mRNAs <- Matrix::colSums(exprs(monocle_cds))

monocle_cds <- monocle_cds[,pData(monocle_cds)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(monocle_cds)$Total_mRNAs)) +
                     2*sd(log10(pData(monocle_cds)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(monocle_cds)$Total_mRNAs)) -
                     2*sd(log10(pData(monocle_cds)$Total_mRNAs)))

qplot(Total_mRNAs, data = pData(monocle_cds), color = "CellTypes", geom =
        "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)


monocle_cds <- monocle_cds[,pData(monocle_cds)$Total_mRNAs > lower_bound &
                             pData(monocle_cds)$Total_mRNAs < upper_bound]
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)


# Log-transform each value in the expression matrix.
L <- log(exprs(monocle_cds[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily


library(reshape2)
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values. (this not working, what is Value?)
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")





#------------------Cluster cells------------------


cth2 <- newCellTypeHierarchy()
cth2 <- addCellType(cth2, "NC", classify_func = function(x) { 
  
  x[MAP1B_id,] >= 1 
    
    
}
)



cth2 <- addCellType(cth2, "NP", classify_func = function(x) {
  
  x[ACAN_id,] >= 1 |    x[COL2A1_id,] >= 1 
  
}
)



cth2 <- addCellType(cth2, "IAF", classify_func = function(x) {
  
  x[COL1A1_id,] >= 1
  
  
}
)

cth2 <- addCellType(cth2, "OAF", classify_func = function(x) {
  
  x[CALR_id,] >= 1 | x[HSPA6_id,] >= 1
  
  
}
)



cth2 <- addCellType(cth2, "NC/NP", classify_func = function(x) {
  
  x[SOX4_id,] >= 1 & x[COL2A1_id,] >= 1
  
  
}
)






monocle_cds2 <- classifyCells(monocle_cds, cth2, 0.1)

table(pData(monocle_cds2)$CellType)

#This pie chart shows cell type
pie <- ggplot(pData(monocle_cds2),
              aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

#clustering using marker genes 

expressed_genes <- row.names(subset(fData(monocle_cds2),
                                    num_cells_expressed >= 1))
marker_diff <- markerDiffTable(monocle_cds2[expressed_genes,],
                               cth2,
                               remove_ambig = TRUE,
                               remove_unknown = TRUE,
                               cores = 1)


candidate_clustering_genes <- row.names(subset(marker_diff, qval < 0.01))
marker_spec <- calculateMarkerSpecificity(monocle_cds2[candidate_clustering_genes,], cth2)
head(selectTopMarkers(marker_spec, 10))

semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
monocle_cds2 <- setOrderingFilter(monocle_cds2, semisup_clustering_genes)
plot_ordering_genes(monocle_cds2)



# irlba package should be re-installed if this line takes too long 
plot_pc_variance_explained(monocle_cds2, return_all = F)
#We don't include cds3
monocle_cds4 <- reduceDimension(monocle_cds2, max_components = 2, norm_method = "log", num_dim = 3,
                                reduction_method = 'tSNE', verbose = T)

monocle_cds4 <- clusterCells(monocle_cds4, verbose = F)
plot_cell_clusters(monocle_cds4, color_by = "Subtype") 
monocle_cds4 <- clusterCells(monocle_cds4,
                             num_clusters = 3,
                             frequency_thresh = 0.1,
                             cell_type_hierarchy = cth2)
plot_cell_clusters(monocle_cds4, 1, 2, 
                   markers = c("MAP1B", "CD44"))

pie <- ggplot(pData(monocle_cds4), aes(x = factor(1), fill =
                                         factor(Subtype))) +
  geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())


#------------Constructing Single Cell Trajectory---------
#--------choosing gene that defines cell progress-----------
colnames(pData(monocle_cds4))


diff_test_res <- differentialGeneTest(monocle_cds4[expressed_genes,],
                                      fullModelFormulaStr = "~Subtype")

ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

monocle_cds5 <- setOrderingFilter(monocle_cds4, ordering_genes)
plot_ordering_genes(monocle_cds5)

#----------Reduce data dimensionality---------------

#Trajectory
monocle_cds5 <- reduceDimension(monocle_cds5, max_components = 2,
                                method = 'DDRTree')

monocle_cds6 <- orderCells(monocle_cds5, reverse = TRUE)
Plot1a <- plot_cell_trajectory(monocle_cds6, color_by = "Subtype", cell_size = 1,  show_branch_points = FALSE, show_backbone=FALSE) + theme(aspect.ratio = 1) + scale_color_manual(values=c("dodgerblue1", "plum","darkorange2","cyan2", "forestgreen"))

Plot1c <- plot_cell_trajectory(monocle_cds6, color_by = "MajorType", cell_size = 1,  show_branch_points = FALSE, show_backbone=FALSE) + theme(aspect.ratio = 1) + scale_color_manual(values=c("lightcoral","gray30"))

Plot1b <- plot_cell_trajectory(monocle_cds6, color_by = "State", cell_size = 1,  show_branch_points = FALSE, show_backbone=FALSE) + theme(aspect.ratio = 1) + scale_color_manual(values=c( "gray30","gray40","gray60","mistyrose3", "lightcoral", "lightslategrey", "mistyrose4"))

Plot1d <- plot_cell_trajectory(monocle_cds6, color_by = "Pseudotime", cell_size = 1,  show_branch_points = FALSE, show_backbone=FALSE) + theme(aspect.ratio = 1) 

Plot1d1 <- plot_cell_trajectory(monocle_cds6,  cell_size = 3,  show_branch_points = FALSE, color_by = "MajorType",show_backbone=FALSE,
                                markers= c("MAP1B", "SOX4", "SOX17", "ITGA6", "BASP1", "AKR1C1", "FABP4", "APOE"),
                                use_color_gradient= TRUE ) + theme(aspect.ratio = 1)

Plot1d2 <- plot_cell_trajectory(monocle_cds6,  cell_size = 3, show_branch_points = FALSE, color_by = "MajorType",show_backbone=FALSE,
                                markers= c("ACAN", "COL2A1", "SOX9", "COL1A1", "CALR", "HSPA6", "CD44", "CD14"),
                                use_color_gradient= TRUE ) + theme(aspect.ratio = 1)

#Jitter plot
to_be_tested <- row.names(subset(fData(monocle_cds6),
                                 gene_short_name %in% c("ACAN", "COL2A1", "COL1A1", "SOX4", "SOX9", "MAP1B", "FABP4", "AKR1C1", "APOE")))
cds6_subset <- monocle_cds6[to_be_tested,]
Plot2b1 <- plot_genes_jitter(cds6_subset,
                             grouping = "MajorType",
                             color_by = "Subtype",
                             ncol = NULL,
                             plot_trend = TRUE, cell_size = 1) + scale_color_manual(values=c("palegreen3", "dodgerblue","orchid3","sienna1", "black"))

Plot2b2 <- plot_genes_jitter(cds6_subset,
                             grouping = "Subtype",
                             color_by = "MajorType",
                             ncol = NULL,
                             plot_trend = TRUE, cell_size = 2) + scale_color_manual(values=c("tomato3","black"))
Plot2b3 <- plot_genes_jitter(cds6_subset,
                             grouping = "State",
                             color_by = "MajorType",
                             ncol = NULL,
                             plot_trend = TRUE, cell_size = 2) + scale_color_manual(values=c("tomato3","black"))
Plot2b4 <- plot_genes_jitter(cds6_subset,
                             grouping = "MajorType",
                             color_by = "State",
                             ncol = NULL,
                             plot_trend = TRUE, cell_size = 2) + scale_color_manual(values=c( "red2","lightcoral","darkgoldenrod1","gray60", "black"))


Plot2b5 <- plot_genes_jitter(cds6_subset,
                             grouping = "State",
                             color_by = "Subtype",
                             ncol = NULL,
                             plot_trend = TRUE, cell_size = 2) + scale_color_manual(values=c("palegreen3", "dodgerblue","orchid3","sienna1"))

Plot2b6 <- plot_genes_jitter(cds6_subset,
                             grouping = "Subtype",
                             color_by = "State",
                             ncol = NULL,
                             plot_trend = TRUE, cell_size = 2) + scale_color_manual(values=c( "red2","lightcoral","darkgoldenrod1","gray60", "black"))

Plot2b5 + Plot2b6
IVD.integrated3@meta.data

#Pseudotime

to_be_tested <- row.names(subset(fData(monocle_cds6),
                                 gene_short_name %in% c("MAP1B", "AKR1C1", "FABP4", "APOE", "SOX4")))
to_be_tested2 <- row.names(subset(fData(monocle_cds6),
                                  gene_short_name %in% c("ACAN", "COL2A1", "SOX9", "COL1A1", "SOX4")))

cds6_pseudotime_subset <- monocle_cds6[to_be_tested2,]

diff_test_res <- differentialGeneTest(cds6_pseudotime_subset,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res[,c("gene_short_name", "pval", "qval")]

NC_cells <- row.names(subset(pData(monocle_cds6), Subtype == "NC"))
monocle_cds6_NC <- monocle_cds6[,NC_cells]
cds6_NC_pseudotime_subset <- monocle_cds6_NC[to_be_tested,]


NP_cells <- row.names(subset(pData(monocle_cds6), Subtype == "NP"))
monocle_cds6_NP <- monocle_cds6[,NP_cells]
cds6_NP_pseudotime_subset <- monocle_cds6_NP[to_be_tested2,]



Plot2a1 <- plot_genes_in_pseudotime(cds6_pseudotime_subset, color_by = "Subtype", ncol = 6, cell_size = 0.5, label_by_short_name = FALSE) + scale_color_manual(values=c("palegreen3", "dodgerblue","orchid3","sienna1"))

Plot2a1_NC <- plot_genes_in_pseudotime(cds6_NC_pseudotime_subset, min_expr = 0.01, color_by = "Subtype",  ncol = 5, cell_size = 0.1, label_by_short_name = FALSE) + scale_color_manual(values=c("palegreen3", "dodgerblue","orchid3","sienna1")) + scale_x_reverse()

Plot2a1_NP <- plot_genes_in_pseudotime(cds6_NP_pseudotime_subset, min_expr = 0.1, color_by = "Subtype",  ncol = 5, cell_size = 0.1, label_by_short_name = FALSE) + scale_color_manual(values=c("dodgerblue", "dodgerblue","orchid3","sienna1")) + scale_x_reverse()


Plot2a2 <- plot_genes_in_pseudotime(cds6_pseudotime_subset, ncol = 5, color_by = "MajorType",  cell_size = 0.1, label_by_short_name = FALSE) + scale_color_manual(values=c("lightcoral","gray30")) + scale_x_reverse()


Plot2a3 <- plot_genes_in_pseudotime(cds6_pseudotime_subset, color_by = "State", ncol = 6, cell_size = 0.1, label_by_short_name = FALSE) + scale_color_manual(values=c("gray30","gray40","gray60","mistyrose3", "lightcoral", "lightslategrey", "mistyrose4"))+ scale_x_reverse()


Plot2a1 + Plot2a2 + Plot2a3


#Pseudotime branched
to_be_tested <- row.names(subset(fData(monocle_cds6),
                                 gene_short_name %in% c("SOX4", "SOX9", "MAP1B", "ACAN")))
plot_genes_branched_pseudotime(monocle_cds6[to_be_tested,],
                               branch_states = NULL, branch_point = 1,
                               branch_labels = NULL, method = "fitting", min_expr = NULL,
                               cell_size = 2, nrow = NULL, ncol = 1, panel_order = NULL,
                               color_by = "MajorType", expression_curve_linetype_by = "Branch",
                               trend_formula = "~ sm.ns(Pseudotime, df=3) * Branch",
                               reducedModelFormulaStr = NULL, label_by_short_name = TRUE,
                               relative_expr = TRUE)



plot_cell_trajectory(monocle_cds6, color_by = "MajorType") + theme(aspect.ratio = 1)
table(pData(monocle_cds6)$Subtype)
table(pData(monocle_cds6)$CellType)
plot_cell_trajectory(monocle_cds6, color_by = "State") + theme(aspect.ratio = 1)
plot_cell_trajectory(monocle_cds6, color_by = "CellType") + theme(aspect.ratio = 1)
monocle_cds6$seurat_clusters
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

plot_cell_trajectory(monocle_cds6, color_by = "Pseudotime")



#Stop Here-----------------
