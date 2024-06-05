# spatial RNAseq clustering
# Interfaculty Bioinformatics Unit (IBU) (https://www.bioinformatics.unibe.ch

library(Seurat)     #v4.1
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggpubr)
library(clusterProfiler)
library(clustree)
library(future)
library(openxlsx)	
library(RColorBrewer)	




#library(STutility)
#library(topGO) 
#library(limma)
#library(biomaRt)
#library(msigdbr)
#library(scales)


# setup --------------------------

# Specify input directories and files
workPath <- "/data/projects/p907_SIBdays_tutorial_-_spatial_transcriptomics/"            ### ADJUST

analysisFolder <- paste0(workPath, "scripts/analysis")

name <- "Posterior"            ### ADJUST
#---------------------------------


seu.norm <- readRDS(paste0(analysisFolder, "/normalized.rds"))

#parallelization
plan("multisession", workers = availableCores())
options(future.globals.maxSize = 8000 * 1024^2)


seu.norm <- RunPCA(seu.norm, assay = "SCT", npcs = 50, verbose = FALSE)
seu.norm <- RunUMAP(seu.norm, reduction = "pca", dims = 1:50)

DimPlot(seu.norm, reduction = "pca", group.by = "orig.ident") + NoLegend()
DimPlot(seu.norm, reduction = "umap", group.by = "orig.ident")+ ggtitle("") + NoLegend()

#seu.norm <- ScaleData(seu.norm, features=rownames(seu.norm), verbose = FALSE)


# 1. Dimensionality reduction  #########################
  # Tool: [Seurat](https://satijalab.org/seurat/articles/spatial_vignette.html)

# Selecting which PCs to use: ---------
# To overcome the extensive technical noise in any single gene, Seurat clusters spots based 
# on their PCA scores, with each PC essentially representing a metagene that combines information 
# across a correlated group of genes. 

# ElbowPlot ranks the principal components based on the variance explained by each. This plot typically 
# shows an "elbow", which can be used to assess how many PCs are needed to capture most of the signal in the data.

use.pcs  <- 50
ElbowPlot(seu.norm, ndims=use.pcs)

# The sctransform workflow performs more effective normalization, strongly removing technical effects from the data,
# this means that higher PCs are more likely to represent subtle, but biologically relevant, sources of heterogeneity
# - so including them may improve downstream analysis. Therefore, higher number of PC can be used.
# By default we are using the first 50 PCs



# 2. Identifying clusters  #########################
# Seurat implements a graph-based clustering approach. Distances between the spots are calculated based on
# previously identified PCs. Briefly, Seurat identifies clusters of spots by a shared nearest neighbor (SNN)
# modularity optimization based clustering algorithm. First, it identifies k-nearest neighbors (KNN) and constructs
# the SNN graph. Then it optimizes the modularity function to determine clusters. For a full description of the
# algorithms, see Waltman and van Eck (2013) The European Physical Journal B.

# The FindClusters function implements the procedure, and contains a resolution parameter that sets the granularity
# of the downstream clustering, with increased values leading to a greater number of clusters.

# Selecting which resolution to use: -------
resolution_vector <- seq(0.2,2,0.2)
seu.norm <- FindNeighbors(seu.norm, reduction="pca", dims=1:use.pcs)
seu.norm <- FindClusters(object=seu.norm, resolution=resolution_vector, verbose=FALSE)

resTable <- as.data.frame(sapply(grep("res",colnames(seu.norm@meta.data),value = TRUE), function(x) length(unique(seu.norm@meta.data[,x]))))
colnames(resTable) <- "number of clusters"

# how many clusters in each resolution
resTable

for(i in seq(resolution_vector)){
  print(DimPlot(seu.norm, reduction = "umap", label=TRUE, group.by=paste0("SCT_snn_res.", resolution_vector[i])) + labs(title=paste0("res_", resolution_vector[i])))
  print(SpatialDimPlot(seu.norm, label = TRUE, combine=FALSE, label.size = 3, group.by=paste0("SCT_snn_res.", resolution_vector[i])))
}

# Plotting clustering trees
# To build a clustering tree we need to look at how spots move as the clustering resolution is increased. 
# Each cluster forms a node in the tree and edges are constructed by considering the spots in a cluster 
# at a lower resolution that end up in a cluster at the next higher resolution. By connecting clusters 
# in this way we can see how clusters are related to each other, which are clearly distinct and which are unstable.

# If there are nodes with multiple incoming edges, it is a good indication that we have over-clustered the data.

clustree(seu.norm, prefis="res.")



# -> Set the resolution to 0.4 ############
resolution <- "0.4" 

# Final clustering:
Idents(seu.norm) <- paste0("SCT_snn_res.", resolution)
seu.norm$finalCluster <- Idents(seu.norm)

DimPlot(seu.norm, reduction = "umap", label=TRUE, group.by = "ident") 
SpatialDimPlot(seu.norm, label = TRUE, combine=FALSE, label.size = 3)


# As there are many colors, it can be challenging to visualize which spot belongs to which cluster. 
# Therefore, we highlight spots of each cluster seperately:
SpatialDimPlot(seu.norm, combine=TRUE, images=name, cells.highlight = CellsByIdentities(object = seu.norm), cols.highlight = c("red", "NA"), facet.highlight = TRUE, ncol = 4)


# Total number of cells per cluster per slice
seu.norm@meta.data %>% dplyr::select(orig.ident, finalCluster) %>% table() 



# 3. Identify marker genes for each cluster #######################
# Identify unique upregulated marker genes for each clusters (based on normalized data)

## gvg: clustering done, set default assay back to SCT for DGE analyses
DefaultAssay(seu.norm) <- "SCT"

#prepare object to run differential expression on SCT assay with multiple models
integrated <- PrepSCTFindMarkers(seu.norm)

saveRDS(seu.norm, paste0(analysisFolder, "/clustered_res", resolution, ".rds"))

# find marker genes of all clusters (only overexpressed genes -> only.pos=TRUE)
markers_all <- FindAllMarkers(object=seu.norm, only.pos=TRUE, min.pct = 0.25)
markers_all <- markers_all[markers_all$p_val_adj < 0.05,]

#average log2 fold change is calculated:
# data1 = log2(mean(expm1(ExpressionData1)) + 1)
# data2 = log2(mean(expm1(ExpressionData2)) + 1)
# avg_log2FC = data1-data2

#average log2 fold change is calculated (for "negbinom", "poisson", or "DESeq2" DE tests):
# data1 = log2(mean(ExpressionData1) + 1)
# data2 = log2(mean(ExpressionData2) + 1)
# avg_log2FC = data1-data2


# only genes unique between clusters
markers_all_single <- markers_all[markers_all$gene %in% names(table(markers_all$gene))[table(markers_all$gene) == 1],]

# nb of unique marker genes per cluster
summary <- as.data.frame(table(markers_all_single$cluster))
colnames(summary) <- c("cluster", "Nb unique marker genes")
summary 

write.xlsx(markers_all_single %>% dplyr::select(cluster, everything()), paste0(analysisFolder, "/uniqueUpregulatedGenes_perCluster_res", resolution, ".xlsx"), sheetName="Sheet1", colNames=TRUE, rowNames=TRUE, append=FALSE)
# *Note: p_val: p-value, avg_log2FC: log2 fold-change in the average expression between the two groups.
# Positive values indicate that the gene is more highly expressed within the cluster, pct.1: percentage 
# of spots where the gene is detected within the cluster, pct.2: percentage of spots where the gene is 
# detected outside the cluster, p_val_adj: adjusted p-value based on Bonferroni correction using all 
# genes in the dataset, cluster: cluster ID, gene: gene name
  
# Heatmap of top 8 unique marker genes between clusters:
# heatmap of genes by cluster for the top 8 marker genes per cluster
top <- markers_all_single %>% group_by(cluster) %>% top_n(8, avg_log2FC)
DoHeatmap(object = seu.norm, features = top$gene)


# Dotplot of top 8 **unique marker genes** between clusters:	
DotPlot(seu.norm, features = top$gene, dot.scale = 8) + coord_flip() +  RotatedAxis()	


# Identify all upregulated genes for each cluster ----------------
# For each cluster, we identify all genes that are upregulated in spots in the cluster compared to spots outside the cluster. 
# In contrast to the analysis of unique upregulated marker genes above, we do not require that a gene is specific to a single 
# cluster, i.e. a given gene can be detected for more than one cluster.

for(cluster in sort(unique(markers_all$cluster))){
  
  markers_all_cluster <- markers_all[markers_all$cluster == cluster, ]
  rownames(markers_all_cluster) <- markers_all_cluster$gene
  
  print(paste0("### Cluster: ",cluster))
  print(paste0("Markers genes for cluster ",cluster, ":"))
  write.xlsx(markers_all_cluster, paste0(analysisFolder, "/upregulatedGenes_res", resolution, "_cluster", cluster, ".xlsx"), sheetName="Sheet1", colNames=TRUE, rowNames=TRUE, append=FALSE)
  # *Note: p_val: p-value, avg_log2FC: log2 fold-change in the average expression between the two groups. Positive values 
  # indicate that the gene is more highly expressed within the cluster, pct.1: percentage of spots where the gene is detected 
  # within the cluster, pct.2: percentage of spots where the gene is detected outside the cluster, p_val_adj: adjusted p-value 
  # based on Bonferroni correction using all genes in the dataset, cluster: cluster ID, gene: gene name
  
  print(FeaturePlot(seu.norm, head(markers_all_cluster$gene, n=4), cols = c("lightgrey", "blue"), ncol = 2))
  print(SpatialFeaturePlot(seu.norm, features = head(markers_all_cluster$gene, n=4), alpha = c(0.1, 1)))
}
  


  
# 4. Identify spatial variable features ##########
# Search for features exhibiting spatial patterning in the absence of pre-annotation. The method, is inspired 
# by the Trendsceek, which models spatial transcriptomics data as a mark point process and computes a variogram, 
# which identifies genes whose expression level is dependent on their spatial location. More specifically, this 
# process calculates gamma(r) values measuring the dependence between two spots a certain "r" distance apart. 
# By default, we use an r-value of '5' in these analyses, and only compute these values for variable genes 
# (where variation is calculated independently of spatial location) to save time.\n\nTop 9 of spatial variable genes:

seu.norm <- FindSpatiallyVariableFeatures(seu.norm, assay = "SCT", features = rownames(GetAssayData(seu.norm, slot = "scale.data")), selection.method = "markvariogram")
  
spatialFeatures <- SVFInfo(seu.norm, selection.method="markvariogram", status=TRUE)
spatialFeatures <- na.omit(spatialFeatures)
spatialFeatures <- spatialFeatures %>% filter(variable == TRUE) %>% arrange(rank)
spatialFeatures <- spatialFeatures[,c("r.metric.5"), drop=FALSE]
 
write.xlsx(spatialFeatures, paste0(analysisFolder, "/spatialVariableFeatures_res", resolution, "_", names, ".xlsx"), sheetName="Sheet1", colNames=TRUE, rowNames=TRUE, append=FALSE) 
top.features <- head(SpatiallyVariableFeatures(integrated, selection.method = "markvariogram"), 9)
SpatialFeaturePlot(seu.norm, features=top.features, ncol=3, alpha = c(0.1, 1))







############ Heidi: until here #####################




# 5. Cell type identification (only human or mouse)
# ScType a computational method for automated selection of marker genes based merely on scRNA-seq data.
# Reference from CellMarker 2.0: One of the largest available databases for cell-type marker. Cell markers 
# of different cell types from different tissues with more than one evidence were used.

# Be aware this is an automatic annotation and requires manual curation!
  

# DB file
doAnnotation <- TRUE
species <- ""
if(msigdb_species == "Homo sapiens"){
  species <- "Human"
} else if(msigdb_species == "Mus musculus"){
  species <- "Mouse"
} else {
  cat("**--> nothing done, as species is neither human nor mouse!")
  doAnnotation <- FALSE
}
# ```

# ```{r databaseSetup , eval=doAnnotation, echo=FALSE, warning=FALSE, message=FALSE, results="hide"}
# load gene set preparation function
library("HGNChelper")
source("utils/sc-type/gene_sets_prepare.R")
# load cell type annotation function
source("utils/sc-type/sctype_score_.R")

# create own db (downloaded from CellMarker 2.0)
tryCatch({
  own_db_file <- read.xlsx("http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download_files/file/Cell_marker_All2.xlsx")       
},error =function(e){
  own_db_file <<- read.xlsx("utils/Cell_marker_All.xlsx")
}
)
own_db_file <- own_db_file[own_db_file$species == species, c("tissue_class", "cell_name", "Symbol")] 

# remove empty symbol lines
own_db_file <- own_db_file[!is.na(own_db_file$Symbol), ]

# merge marker names into one line and remove entries with single support
own_db_file$tissue_cell_name <- paste0(own_db_file$tissue_class, own_db_file$cell_name)
own_db <- own_db_file %>% group_by(tissue_cell_name) %>% 
  mutate(Symbol = paste(unique(Symbol[duplicated(Symbol)]), collapse=",")) %>%
  distinct() %>% arrange(tissue_class)

# remove empty symbol lines
own_db <- own_db[own_db$Symbol != "", ]


own_db <- own_db[,c("tissue_class", "cell_name", "Symbol")]
colnames(own_db) <- c("tissueType", "cellName", "geneSymbolmore1")
own_db$geneSymbolmore2 <- ""
own_db$shortName <- ""

#remove entries with < 3 cellName
own_db <- own_db[own_db$tissueType %in% names(table(own_db$tissueType)[table(own_db$tissueType) > 2]),]

write.xlsx(own_db, paste0(analysisFolder, "/referenceDB.xlsx"), sheetName="Sheet1", colNames=TRUE, rowNames=FALSE, append=FALSE)

# automatically detect tissue type of the dataset if not set
if(tissue == "" || length(grep(tissue, own_db$tissueType, ignore.case=TRUE, value=TRUE)) == 0){
  source("utils/sc-type/auto_detect_tissue_type.R")
  # guess a tissue type
  # if saled = TRUE, make sure the data is scaled, as seuratObject[[assay]]@scale.data is used. 
  # If you just created a Seurat object, without any scaling and normalization, set scaled = FALSE, seuratObject[[assay]]@counts will be used         
  
  tissue_guess <- auto_detect_tissue_type(path_to_db_file = paste0(analysisFolder, "/referenceDB.xlsx"), seuratObject = integrated, scaled = TRUE, assay = "SCT") 
  write.table(tissue_guess, paste0(analysisFolder, "/tissue_guess.txt"), quote=FALSE, row.names=FALSE)
  
  tissue <-  tissue_guess$tissue[1]
} else {
  #make sure its written the correct way
  tissue <- grep(tissue, own_db$tissueType, ignore.case=TRUE, value=TRUE)[1]
}
# ```

# ```{r autodetectedTissue , eval=doAnnotation,echo=FALSE, warning=FALSE, message=FALSE,results="asis"}
cat("\n")
cat(paste0("**-> Tissue autodetection: ", tissue, "**"))
# 
# ```

# `r if(doAnnotation){"<br>"}`

# `r if(doAnnotation){"### Top 5 cell type annotation for each cluster"}`

# ```{r cellTypeIdentification , eval=doAnnotation, echo=FALSE, warning=FALSE, message=FALSE}
# prepare gene sets
gs_list <- gene_sets_prepare(paste0(analysisFolder, "/referenceDB.xlsx"), tissue)

# get cell-type by cell matrix
es.max <- sctype_score(scRNAseqData = integrated[["SCT"]]@scale.data, scaled = TRUE, 
                       gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(integrated@meta.data$finalCluster), function(cl){
  es.max.cl <- sort(rowSums(es.max[ ,rownames(integrated@meta.data[integrated@meta.data$finalCluster==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(integrated@meta.data$finalCluster==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
sctype_scores_top5 <- cL_resutls %>% group_by(cluster) %>% top_n(n = 5, wt = scores)

sctype_scores_top5 %>%  kable() %>% kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
# ```
# `r if(doAnnotation){"*Note: cluster: cluster id, type: cell type annotation, scores: Cell-type specificity score provides a quantitative measure of how uniquely markers identify a specific cell-type of the given tissue. A cell-type with the highest ScType score is used for assignment to the cluster. Low ScType score (less than quarter the number of cells in a cluster), or a negative ScType score indicates a low-confidence cell-type annotation, which will be assigned as “unknown” cell type, ncells: number of cells in this cluster*"}`
# 
# `r if(doAnnotation){"<br>"}`
# 
# ```{r cellTypeIdentification_dimPlot , eval=doAnnotation, echo=FALSE, warning=FALSE, message=FALSE, results="hide", fig.height=5, fig.width = 9}
# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

# UMAP
integrated@meta.data$clusterAnnotation = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  integrated@meta.data$clusterAnnotation[integrated@meta.data$finalCluster == j] = as.character(cl_type$type[1])
}

DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'clusterAnnotation')
```


# `r if(doAnnotation){"<br>"}`

# `r if(doAnnotation){"Bubble plot showing all the cell types that were considered by ScType for cluster annotation. The outter (grey) bubbles correspond to each cluster (the bigger bubble, the more cells in the cluster), while the inner bubbles correspond to considered cell types for each cluster, with the biggest bubble corresponding to assigned cell type."}`

# ```{r cellTypeIdentification_bubblePlot , eval=doAnnotation, echo=FALSE, warning=FALSE, message=FALSE, results="hide", fig.height=7, fig.width = 9}
# load libraries
library(ggraph)
library(igraph)

# prepare edges
cL_resutls <- cL_resutls[order(cL_resutls$cluster),]
edges <- cL_resutls
edges$type <- paste0(edges$type,"_",edges$cluster)
edges$cluster <- paste0("cluster ", edges$cluster)
edges <- edges[,c("cluster", "type")]
colnames(edges) <- c("from", "to")
rownames(edges) <- NULL

# prepare nodes
nodes_lvl1 <- sctype_scores[,c("cluster", "ncells")]
nodes_lvl1$cluster <- paste0("cluster ", nodes_lvl1$cluster)
nodes_lvl1$Colour <- "#f1f1ef"
nodes_lvl1$ord <- 1
nodes_lvl1$realname <- nodes_lvl1$cluster
nodes_lvl1 <- as.data.frame(nodes_lvl1)
nodes_lvl2 <- c() 
ccolss <- hue_pal()(length(unique(cL_resutls$cluster)))
for (i in 1:length(unique(cL_resutls$cluster))){
  dt_tmp <- cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]
  nodes_lvl2 <- rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes <- rbind(nodes_lvl1, nodes_lvl2)
nodes$ncells[nodes$ncells<1] <- 1
files_db <- own_db[,c("cellName","shortName")]
#files_db <- openxlsx::read.xlsx(db_)[,c("cellName","shortName")]
files_db <- unique(files_db)
nodes <- merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName) | nodes$shortName == ""] <- nodes$realname[is.na(nodes$shortName) | nodes$shortName == ""]
nodes <- nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

mygraph <- graph_from_data_frame(edges, vertices=nodes)

# Make the graph
plotList <- list()
plotList[[1]] <- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + 
  geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + 
  geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#000000"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.2))) + 
  geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")

plotList[[2]] <- DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss)
ggarrange(plotlist=plotList, ncol=1, nrow=1, hjust=-1)



# reset finalCluster
integrated$finalCluster_number <- integrated$finalCluster

Idents(integrated) <- integrated@meta.data$clusterAnnotation
integrated$finalCluster <- Idents(integrated)
saveRDS(integrated, paste0(analysisFolder, "/clustered_cellTypeAnnotated_res", resolution, "_", comparison, ".rds"))
# ```


# <br>
#   
#   `r if(length(uniqueConditions) <= 1){"# 8.  Cell-type deconvolution"}`
# `r if(length(uniqueConditions) > 1){"# 10.  Cell-type deconvolution"}`
Tools:
  
  # * [STdeconvolve](https://bioconductor.org/packages/release/bioc/html/STdeconvolve.html)
# * [CellMarker](http://bio-bigdata.hrbmu.edu.cn/CellMarker/index.html)

# STdeconvolve is an unsupervised, reference-free approach to infer latent cell-type proportions and transcriptional profiles within multi-cellular spatially-resolved pixels from spatial transcriptomics (ST) datasets. STdeconvolve builds on latent Dirichlet allocation (LDA), a generative statistical model commonly used in natural language processing. In the context of ST data, given a count matrix of gene expression in multi-cellular ST pixels, STdeconvolve applies LDA to infer the putative transcriptional profile for each cell-type and the proportional representation of each cell-type in each multi-cellular ST pixel. 

# <br>
  
  # ```{r STdeconvolve_selectK , eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE, results="hide"}

library(STdeconvolve)

# this is the genes x barcode sparse count matrix
cd <- integrated[["SCT"]]@counts

# xy coordinates
addition <- 0
pos <- data.frame()
posList <- list()
for(name in names){
  posImage <- GetTissueCoordinates(integrated, image = name)
  colnames(posImage) <- c("y", "x")
  #reverse y axis
  posImage$y <- 2*tail(posImage$y,1)-posImage$y
  posList[[name]] <- posImage
  posImage$y <- posImage$y + addition
  pos <- rbind(pos, posImage)
  addition <- addition + 1000
}


## feature select for genes
corpus <- restrictCorpus(cd,
                         removeAbove = 1.0,
                         removeBelow = 0.05,
                         alpha = 0.05,
                         plot = FALSE,
                         verbose = TRUE)

#choose optimal K         
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2, 15, by = 1),
               perc.rare.thresh = 0.05,
               plot=TRUE,
               verbose=TRUE)

## select model 
optLDA <- optimalModel(models = ldas, opt = "kneed")
optLDA2 <- optimalModel(models = ldas, opt = "min")

if(optLDA2@k > optLDA@k){
  optLDA <- optimalModel(models = ldas, opt = round(mean(c(optLDA2@k, optLDA@k))))
}
# ```

# **-> K of `r optLDA@k` chosen**
  
  # <br>
  
  #### Scatterpies of proportions of all deconvolved cell-types:
  # Cell-types from spots that contribute less than 5% of the spot proportions are removed and the deconvolved transcriptional profiles are scaled by 1000 

# ```{r STdeconvolve_proportionPlot , eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE, results="hide", fig.height=9, fig.width = 9}               
## Extract pixel cell-type proportions (theta) and cell-type gene expression
## profiles (beta) for the given dataset.
## We can also remove cell-types from pixels that contribute less than 5% of the
## pixel proportion and scale the deconvolved transcriptional profiles by 1000 
results <- getBetaTheta(optLDA,
                        perc.filt = 0.05,
                        betaScale = 1000)

deconProp <- results$theta
deconGexp <- results$beta

plotList <- list()
for(name in names){
  cells <- rownames(posList[[name]])
  #visualize the proportion of each deconvolved cell-type across the original spatially resolved pixels. (long runtime!)
  plotList[[name]] <- vizAllTopics(deconProp[cells,], posList[[name]], r=1.8, lwd = 0, showLegend = TRUE, plotTitle = name) +
    ggplot2::guides(colour = "none")
}
ggarrange(plotlist=plotList, ncol=1, nrow=1, hjust=-1)

#plot seperately for each celltype
#vizTopic(theta = deconProp, pos = pos, topic = "5", plotTitle = "X5",
#         size = 5, stroke = 1, alpha = 0.5,
#         low = "white",
#         high = "red")
# ```


#### Compare deconvolved cell-types to clustering
Correlation matrix of the correlations between the proportions of each cell-type and the clustering

# ```{r STdeconvolve_correlation , eval=TRUE,echo=FALSE, warning=FALSE, message=FALSE, results="asis"}
cluster_proxyTheta <- model.matrix(~ 0 + integrated@meta.data$finalCluster_number)
rownames(cluster_proxyTheta) <- rownames(integrated@meta.data)

corMat_prop <- getCorrMtx(m1 = as.matrix(cluster_proxyTheta),
                          m2 = deconProp,
                          type = "t")

rownames(corMat_prop) <- paste0("cluster_", levels(integrated@meta.data$finalCluster_number))
colnames(corMat_prop) <- paste0("decon_", seq(ncol(corMat_prop)))

correlationPlot(mat = corMat_prop, colLabs = "Transcriptional clusters", rowLabs = "STdeconvolve") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)
  )
# ```

# <br>
  
  ## Marker genes for each deconvoluted cell-type
  Marker genes are defined here as genes highly expressed in the deconvolved cell-type (count > 2) that also have the highest log2(fold change) (>1) when comparing the deconvolved cell-type’s expression profile to the average of all other deconvolved cell-types expression profiles.

# ```{r STdeconvolve_markerGenes , eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE, results="hide"} 

#Now, let’s get the differentially expressed genes for each deconvolved cell-type transcriptional profile and label the top expressed genes for each cell-type:

markerList <- lapply(colnames(deconProp), function(celltype) {
  celltype <- as.numeric(celltype)
  ## highly expressed in cell-type of interest
  highgexp <- names(which(deconGexp[celltype,] > 2))
  
  ## high log2(fold-change) compared to other deconvolved cell-types
  markerTable <- data.frame(celltype=celltype, gene=highgexp, log2FC=log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])),
                            exp.profile.celltype=deconGexp[celltype,highgexp],
                            mean.exp.profile.other.celltypes=colMeans(deconGexp[-celltype,highgexp]))
  
  markerTable <- markerTable[markerTable$log2FC > 1,]
  
  markerTable <- markerTable %>% arrange(desc(log2FC))
  return(markerTable)
})


# ```


# ```{r STdeconvolve_markerGenesPlot , eval=TRUE,echo=FALSE, warning=FALSE, message=FALSE, fig.height=10, fig.width=13, results="asis"}
celltype <- 0
for(markerTable in markerList){
  celltype <- celltype + 1;
  cat("\n")
  cat("\n")
  cat(paste0("### Cell-tpye: ",celltype))
  cat("\n")
  cat(paste0("Markers genes for cell-type ",celltype, ":"))
  cat("\n")
  cat("\n")
  write.xlsx(markerTable, paste0(analysisFolder, "/STdeconvolve_markerGenes_cellType", i, "_", comparison, ".xlsx"), sheetName="Sheet1", colNames=TRUE, rowNames=TRUE, append=FALSE)
  
  cat(markerTable %>%
        kable() %>%
        kable_styling() %>%
        scroll_box(width = "100%", height = "300px"))
  cat("\n")
  cat("*Note: celltype: deconvolved cell-type, gene: gene name, log2FC: log2 fold-change in comparing deconvolved cell-type expression profile to the average of all other deconvolved cell-types expression profiles. Positive values indicate that the gene is more highly expressed within the deconvolved cell-type, exp.profile.celltype: expression profile of deconvolved cell-type, mean.exp.profile.other.celltypes: mean expression profile of all other deconvolved cell-types*\n")
  cat("\n")
  cat("\n")
  
  if(nrow(markerTable) > 0){
    print(FeaturePlot(integrated, head(markerTable$gene, n=4), cols = c("lightgrey", "blue"), ncol = 2))
  }
  
  if(nrow(markerTable) > 1){
    for(i in c(1:2)){
      if(i <= length(markerTable$gene)){
        #split into images of 2
        imageNames <- names(integrated@images)
        j <- 1
        while(j <= length(imageNames)){
          if(j+1 <= length(imageNames)){
            print(SpatialFeaturePlot(integrated, features = markerTable$gene[i], pt.size.factor=point_size, alpha = c(0.1, 1), images=imageNames[j:(j+1)]))
          } else {
            print(SpatialFeaturePlot(integrated, features = markerTable$gene[i], pt.size.factor=point_size, alpha = c(0.1, 1), images=imageNames[j]) + ggtitle(imageNames[j]) + theme(plot.title = element_text(hjust = 0.5)))
          }
          j <- j + 2
        }
        
      }
    }
    cat("\n")
    cat("\n")
  }
  
  
  # GO term enrichment of marker genes in a cell-type
  all.genes <- colnames(deconGexp)
  
  marker.genes <- markerTable$gene
  
  # define geneList as 1 if gene is in marker.genes, 0 otherwise
  geneList <- ifelse(all.genes %in% marker.genes, 1, 0)
  names(geneList) <- all.genes
  
  onts <- c("BP", "MF", "CC")	
  gotable <- data.frame()	
  for(ont in onts){	
    tryCatch({
      # Create topGOdata object	
      GOdata <- new("topGOdata", ontology = ont, # use biological process ontology	
                    allGenes= geneList, geneSelectionFun = function(x)(x == 1),	
                    annot=annFUN.org, mapping=topgo_db, ID="symbol")	
      
      # Test for enrichment using Fisher's Exact Test	
      resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher", cutOff=0.05)	
      gotable.ont <- GenTable(GOdata, Fisher = resultFisher, numChar = 60, topNodes=20)	
      gotable.ont <- gotable.ont[as.numeric(gotable.ont$Fisher) <= 0.05,]	
      if(nrow(gotable.ont) > 0){	
        gotable <- rbind(gotable, cbind(ont, gotable.ont))	
      }	
    },
    error = function(e) {
      cat(paste0("**-> GO enrichment with ", onts, " failed!**"))
      cat("\n")
      cat("\n")
    }
    )
  }
  
  if(nrow(gotable) > 0){
    write.xlsx(gotable, paste0(analysisFolder, "/STdeconvolve_markerGenes_cellType", i,"_enrichedGOterms_", comparison, ".xlsx"), sheetName="Sheet1", colNames=TRUE, rowNames=TRUE, append=FALSE)
    
    cat("Enriched GO-terms (top 20 of each ontology domain):")
    cat("\n")
    cat("\n")
    cat(gotable %>%
          kable() %>%
          kable_styling() %>%
          scroll_box(width = "100%", height = "300px"))
    cat("\n")
    cat("*Note: ont: gene ontology domain, GO.ID: Gene Ontology ID, Term: GO term, Annotated: Number of genes annotated within this GO term, Significant: Number of upregulated genes within this GO term, Expected: Expected number of genes (under null hypothesis) within this GO term, Fisher: p-value of Fisher's exact test*\n")
    cat("\n")
    cat("\n")
  }  else{
    cat("\n")
    cat("\n")
    cat("**-> no GO-terms enriched**")
    cat("\n")
    cat("<br>")
    cat("<br>")
  }
}
# ```

# <br>
#   
#   
#   `r if(doAnnotation){"## Cell-type annotation"}`
# `r if(doAnnotation){"Given a list of reference gene sets for different cell types, we can performed gene set enrichment analysis (GSEA) on the deconvolved transcriptional profiles to test for significant enrichment of any known ground truth cell-types."}`
# 
# `r if(doAnnotation){"Reference from CellMarker 2.0: One of the largest available databases for cell-type marker. Human cell markers of different cell types from different tissues with more than one evidence were used."}`
# 
# `r if(doAnnotation){"**Be aware this is an automatic annotation and requires manual curation!**"}`
# 
# `r if(doAnnotation){"<br>"}`
# ```{r STdeconvolve_annotation , eval=doAnnotation, echo=FALSE, warning=FALSE, message=FALSE, results="hide"}
#GSEA annotation 
# prepare gene sets
gs_list <- gene_sets_prepare(paste0(analysisFolder, "/referenceDB.xlsx"), tissue)

celltype_annotations <- annotateCellTypesGSEA(beta = deconGexp , gset = gs_list$gs_positive, qval = 0.05)

celltype_annotation_table <- data.frame()
i <- 0
for(result in celltype_annotations$results){
  i <- i + 1
  result <- cbind(data.frame(celltype=i, annotation=rownames(result)), result)
  celltype_annotation_table <- rbind(celltype_annotation_table, result)
}
```

# ```{r STdeconvolve_annotation_table , eval=doAnnotation, echo=FALSE, warning=FALSE, message=FALSE, results="asis"}
celltype_annotation_table %>% kable(row.names=FALSE) %>% kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
# ```
# `r if(doAnnotation){"*Note: celltype: deconvolved cell-type, annotation: annotation, p.val: p-value of the enrichment score, q.val: q-value estimation for false discovery rate. This is the significance threshold that should be considered, sscore: enrichment score represent the degree to which a set is over-represented at the top or bottom of the ranked list, edge: edge score*"}`


# ```{r STdeconvolve_proportionPlot2 , eval=doAnnotation, echo=FALSE, warning=FALSE, message=FALSE, results="hide", fig.height=9, fig.width = 9}               

annotationNames <- as.character(celltype_annotations$predictions)
annotationNames[is.na(annotationNames)] <- "unknown"
colnames(deconProp) <- annotationNames

plotList <- list()
for(name in names){
  cells <- rownames(posList[[name]])
  #visualize the proportion of each deconvolved cell-type across the original spatially resolved pixels. (long runtime!)
  plotList[[name]] <- vizAllTopics(deconProp[cells,], posList[[name]], r=1.8, lwd = 0, showLegend = TRUE, plotTitle = name) +
    ggplot2::guides(colour = "none")
}
ggarrange(plotlist=plotList, ncol=1, nrow=1, hjust=-1)
# ```


# `r if(doAnnotation){"#### Compare deconvolved cell-types annotation to clustering annotation"}`
# `r if(doAnnotation){"Correlation matrix of the correlations between the proportions of each cell-type annotation and the clustering annotation"}`

# ```{r STdeconvolve_correlationAnnotation , eval=doAnnotation, echo=FALSE, warning=FALSE, message=FALSE, results="asis"}

deconPropMerged <- t(apply(deconProp, 1, function(x) tapply(x,colnames(deconProp),sum)))

# correlation with clustering
cluster_proxyTheta <- model.matrix(~ 0 + integrated@meta.data$finalCluster)
rownames(cluster_proxyTheta) <- rownames(integrated@meta.data)

corMat_prop <- getCorrMtx(m1 = as.matrix(cluster_proxyTheta),
                          m2 = deconPropMerged,
                          type = "t")

rownames(corMat_prop) <- levels(integrated@meta.data$finalCluster)

correlationPlot(mat = corMat_prop,  colLabs = "Transcriptional clusters", rowLabs = "STdeconvolve") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)
  )
# ```


# <!--
#   # Integration with single-cell data
#   -> see https://satijalab.org/seurat/articles/spatial_vignette.html
# -->
  
  
  


sessionInfo()

