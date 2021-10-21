library(plotly)
levels(x = IVD.aIVD) <- c("NC1","NC2","NP1", "NP2", "NP3","NP4","IAF2","IAF1", "OAF1",'OAF2',"OAF3", "non-MC", 'IC', 'RBC')
levels(x = IVD.nIVD) <- c("NC1","NC2","NP1", "NP2", "NP3","NP4","IAF2","IAF1", "OAF1",'OAF2',"OAF3", "non-MC", 'IC', 'RBC')

IVD.nIVD.3DUMAP <- IVD.nIVD
IVD.nIVD.3DUMAP <- RunUMAP(IVD.nIVD.3DUMAP,
                            dims = 1:10,
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
# Seurat ckyster number and thier assignment to sub-populations:  0-NP1 1-NP2 2-OAF2 3-IC 4-OAF1 5-NP3 6-IAF1 7-NP4 8-IAF2 9-OAF3 10-NC1 11-NC2 12-RBC 13-DDSC
Plot1d1 <- plot_ly(data = plot.data, 
               x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
               color = ~seurat_clusters, 
               colors = c("dodgerblue1",
                          "gray70",
                          "darkorange1",
                          "dodgerblue2",
                          "orange",
                          "dodgerblue3",
                          "dodgerblue4",
                          "orchid3",
                          "plum3",
                          "forestgreen",
                          "darkorange3",
                          "springgreen3",
                          "darkolivegreen",
                          "red"),
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 0.8, width=0.8), 
               text=~label, 
               hoverinfo="seurat_clusters") 


IVD.aIVD.3DUMAP <- IVD.aIVD
IVD.aIVD.3DUMAP <- RunUMAP(IVD.aIVD.3DUMAP,
                           dims = 1:10,
                           n.components = 3L)
umap_1 <- IVD.aIVD.3DUMAP[["umap"]]@cell.embeddings[,1]
umap_2 <- IVD.aIVD.3DUMAP[["umap"]]@cell.embeddings[,2]
umap_3 <- IVD.aIVD.3DUMAP[["umap"]]@cell.embeddings[,3]
Embeddings(object = IVD.aIVD.3DUMAP, reduction = "umap")
plot.data <- FetchData(object = IVD.aIVD.3DUMAP, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters"))
plot.data$label <- paste(rownames(plot.data))
IVD.aIVD.3DUMAP$ClusterNumber
IVD.aIVD.3DUMAP@active.ident
# Seurat ckyster number and thier assignment to sub-populations:  0-NP1 1-NP2 2-OAF2 3-IC 4-OAF1 5-NP3 6-IAF1 7-NP4 8-IAF2 9-OAF3 10-NC1 11-NC2 12-RBC 13-DDSC
Plot1d2 <- plot_ly(data = plot.data, 
                   x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                   color = ~seurat_clusters, 
                   colors = c("dodgerblue1",
                              "gray70",
                              "darkorange1",
                              "dodgerblue2",
                              "orange",
                              "dodgerblue3",
                              "dodgerblue4",
                              "orchid3",
                              "plum3",
                              "forestgreen",
                              "darkorange3",
                              "springgreen3",
                              "darkolivegreen",
                              "red"),
                   type = "scatter3d", 
                   mode = "markers", 
                   marker = list(size = 0.8, width=0.8), 
                   text=~label, 
                   hoverinfo="seurat_clusters") 


Plot1d2.1 <- plot_ly(data = plot.data, 
                   x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                   color = ~seurat_clusters, 
                   colors = c("dodgerblue1",
                              "gray70",
                              "darkorange1",
                              "dodgerblue2",
                              "orange",
                              "dodgerblue3",
                              "dodgerblue4",
                              "orchid3",
                              "plum3",
                              "forestgreen",
                              "darkorange3",
                              "springgreen3",
                              "darkolivegreen",
                              "red"),
                   type = "scatter3d", 
                   mode = "markers", 
                   marker = list(size = 3, width=3), 
                   text=~label, 
                   hoverinfo="seurat_clusters") 


