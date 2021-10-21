library(monocle)
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
library("writexl")

IVD.integrated3 <- RenameIdents(object = IVD.integrated3,  
                                "NP1" = "NP", "NP2" = "NP", 'NP3' = 'NP', 'NP4' = 'NP', "IAF1" = "IAF", "IAF2" = "IAF", "OAF1" = "OAF", "OAF2" = "OAF", "OAF3" = "OAF"
                                )
IVD.integrated3 <- RenameIdents(object = IVD.integrated3,  "NC1" = "NC", "NC2" = "NC")
IVD.integrated3$Subtype <- Idents(IVD.integrated3)

IVD.all.newtraj <- subset(x = IVD.integrated3, idents = c("NC", "NP", "IAF", "OAF"))
data <- as(as.matrix(IVD.all.newtraj@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = IVD.all.newtraj@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

monocle_cds <- newCellDataSet(data, 
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
print(head(fData(monocle_cds)))
expressed_genes <- row.names(subset(fData(monocle_cds),
                                    num_cells_expressed >= 10))
print(head(pData(monocle_cds)))

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

L <- log(exprs(monocle_cds[expressed_genes,]))

melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")

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

cth2 <- newCellTypeHierarchy()

monocle_cds2 <- classifyCells(monocle_cds, cth2, 0.1)

pie <- ggplot(pData(monocle_cds2),
              aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

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

plot_pc_variance_explained(monocle_cds2, return_all = F)

monocle_cds4 <- reduceDimension(monocle_cds2, max_components = 2, norm_method = "log", num_dim = 3,
                                reduction_method = 'tSNE', verbose = T)

monocle_cds4 <- clusterCells(monocle_cds4, verbose = F)
plot_cell_clusters(monocle_cds4, color_by = "Subtype") + theme(aspect.ratio = "1")

monocle_cds4 <- clusterCells(monocle_cds4,
                             num_clusters = 3,
                             frequency_thresh = 0.1,
                             cell_type_hierarchy = cth)
plot_cell_clusters(monocle_cds4, 1, 2, 
                   markers = c("MAP1B", "CD44"))

pie <- ggplot(pData(monocle_cds4), aes(x = factor(1), fill =
                                         factor(Subtype))) +
  geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

colnames(pData(monocle_cds4))
diff_test_res <- differentialGeneTest(monocle_cds4[expressed_genes,],
                                      fullModelFormulaStr = "~Subtype")

ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

monocle_cds5 <- setOrderingFilter(monocle_cds4, ordering_genes)
plot_ordering_genes(monocle_cds5)

monocle_cds5 <- reduceDimension(monocle_cds5, max_components = 2,
                                method = 'DDRTree')

monocle_cds6 <- orderCells(monocle_cds5)
Plot1a <- plot_cell_trajectory(monocle_cds6, color_by = "Subtype", cell_size = 1,  show_branch_points = FALSE, show_backbone=FALSE) + theme(aspect.ratio = 1) + scale_color_manual(values=c("palegreen3", "dodgerblue","orchid3","sienna1"))

Plot1c <- plot_cell_trajectory(monocle_cds6, color_by = "MajorType", cell_size = 1,  show_branch_points = FALSE, show_backbone=FALSE) + theme(aspect.ratio = 1) + scale_color_manual(values=c("tomato3","black"))

Plot1b <- plot_cell_trajectory(monocle_cds6, color_by = "State", cell_size = 1,  show_branch_points = FALSE, show_backbone=FALSE) + theme(aspect.ratio = 1) + scale_color_manual(values=c( "red2","lightcoral","darkgoldenrod1","gray60", "black"))

Plot1d1 <- plot_cell_trajectory(monocle_cds6,  cell_size = 3,  show_branch_points = FALSE, color_by = "MajorType",show_backbone=FALSE,
                                markers= c("MAP1B", "SOX4", "SOX17", "ITGA6", "BASP1", "AKR1C1", "FABP4", "APOE"),
                                use_color_gradient= TRUE ) + theme(aspect.ratio = 1)

Plot1d2 <- plot_cell_trajectory(monocle_cds6,  cell_size = 3, show_branch_points = FALSE, color_by = "MajorType",show_backbone=FALSE,
                                markers= c("ACAN", "COL2A1", "SOX9", "COL1A1", "CALR", "HSPA6", "CD44", "CD14"),
                                use_color_gradient= TRUE ) + theme(aspect.ratio = 1)

#Jitter plot
to_be_tested <- row.names(subset(fData(monocle_cds6),
                                 gene_short_name %in% c("MAP1B", "SOX4", "SOX17", "ITGA6", "BASP1", "AKR1C1", "FABP4", "APOE","ACAN", "COL2A1", "SOX9", "COL1A1", "CALR", "HSPA6", "CD44", "CD14")))
cds6_subset <- monocle_cds6[to_be_tested,]
Plot2b1 <- plot_genes_jitter(cds6_subset,
                  grouping = "MajorType",
                  color_by = "Subtype",
                  ncol = NULL,
                  plot_trend = TRUE, cell_size = 1) + scale_color_manual(values=c("palegreen3", "dodgerblue","orchid3","sienna1"))
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
                                 gene_short_name %in% c("MAP1B", "AKR1C1", "FABP4", "APOE")))
to_be_tested2 <- row.names(subset(fData(monocle_cds6),
                                 gene_short_name %in% c("ACAN", "COL2A1", "SOX9")))
cds6_pseudotime_subset <- monocle_cds6[to_be_tested,]
diff_test_res <- differentialGeneTest(cds6_pseudotime_subset,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_res[,c("gene_short_name", "pval", "qval")]

NC_cells <- row.names(subset(pData(monocle_cds6), Subtype == "NC"))
monocle_cds6_NC <- monocle_cds6[,NC_cells]
cds6_NC_pseudotime_subset <- monocle_cds6_NC[to_be_tested,]


NP_cells <- row.names(subset(pData(monocle_cds6), Subtype == "NP"))
monocle_cds6_NP <- monocle_cds6[,NP_cells]
cds6_NP_pseudotime_subset <- monocle_cds6_NP[to_be_tested2,]

Plot2a1 <- plot_genes_in_pseudotime(cds6_pseudotime_subset, color_by = "Subtype",  ncol = 6, cell_size = 0.5, label_by_short_name = FALSE) + scale_color_manual(values=c("palegreen3", "dodgerblue","orchid3","sienna1"))
Plot2a1_NC <- plot_genes_in_pseudotime(cds6_NC_pseudotime_subset, color_by = "Subtype",  ncol = 6, cell_size = 0.5, label_by_short_name = FALSE) + scale_color_manual(values=c("palegreen3", "dodgerblue","orchid3","sienna1"))
Plot2a1_NP <- plot_genes_in_pseudotime(cds6_NP_pseudotime_subset, color_by = "Subtype",  ncol = 6, cell_size = 0.5, label_by_short_name = FALSE) + scale_color_manual(values=c("dodgerblue", "dodgerblue","orchid3","sienna1"))


Plot2a2 <- plot_genes_in_pseudotime(cds6_pseudotime_subset, color_by = "MajorType", ncol = 6, cell_size = 0.5, label_by_short_name = FALSE) + scale_color_manual(values=c("tomato3","black"))

Plot2a3 <- plot_genes_in_pseudotime(cds6_pseudotime_subset, color_by = "State", ncol = 6, cell_size = 0.5, label_by_short_name = FALSE) + scale_color_manual(values=c("red2","lightcoral","darkgoldenrod1","gray60", "black"))


Plot2a1 + Plot2a2 + Plot2a3


#----------------Stop Here-----------------



