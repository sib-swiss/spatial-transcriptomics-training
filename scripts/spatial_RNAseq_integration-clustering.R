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

name1 <- "Anterior" 
name2 <- "Posterior"            ### ADJUST
#---------------------------------


seu.norm1 <- readRDS(paste0(analysisFolder, "/normalized_", name1, ".rds"))
seu.norm2 <- readRDS(paste0(analysisFolder, "/normalized_", name2, ".rds"))

seu.norm1@meta.data$orig.ident <- name1
seu.norm2@meta.data$orig.ident <- name2

#parallelization
plan("multisession", workers = availableCores())
options(future.globals.maxSize = 8000 * 1024^2)


# 1. Integration analysis  #########################
# Quite often there are strong batch effects between different ST sections, so it may be a good
# idea to integrate the data across sections.
# In Seurat, the methods aim to identify shared cell states that are present across different datasets,
# even if they were collected from different individuals or experimental conditions. The method aims to
# first identify 'anchors' between pairs of datasets utilizing reciprocal PCA (‘RPCA’). When determining
# anchors between any two datasets using RPCA, each dataset is projected into the others PCA space and 
# the anchors are constrained by the same mutual neighborhood requirement. These 'anchors' are then used
# to harmonize the datasets, or transfer information from one dataset to another

# UMAP before integration
allData <- merge(seu.norm1, seu.norm2, merge.data = TRUE)

varFeatures <- VariableFeatures(seu.norm1)
varFeatures <- c(varFeatures, VariableFeatures(seu.norm2))

#remove duplicates
varFeatures <- varFeatures[!duplicated(varFeatures)]
VariableFeatures(allData) <- varFeatures

allData <- RunPCA(allData, npcs = 50, verbose = FALSE)
allData <- RunUMAP(allData, reduction = "pca", dims = 1:50)

DimPlot(allData, reduction = "pca", group.by = "orig.ident") + ggtitle("Before integration")
DimPlot(allData, reduction = "umap", group.by = "orig.ident") + ggtitle("Before integration") 

integrated <- IntegrateLayers(object = allData, method = RPCAIntegration, normalization.method = "SCT", new.reduction = "integrated.rpca", verbose = TRUE)
integrated <- RunUMAP(integrated, dims = 1:50, reduction = "integrated.rpca")
  
DimPlot(integrated, reduction = "umap", group.by = "orig.ident") + ggtitle("After integration")



# 2. Dimensionality reduction  #########################

# Selecting which PCs to use: ---------
# To overcome the extensive technical noise in any single gene, Seurat clusters spots based 
# on their PCA scores, with each PC essentially representing a metagene that combines information 
# across a correlated group of genes. 

# ElbowPlot ranks the principal components based on the variance explained by each. This plot typically 
# shows an "elbow", which can be used to assess how many PCs are needed to capture most of the signal in the data.

use.pcs  <- 50
ElbowPlot(integrated, ndims=use.pcs)

# The sctransform workflow performs more effective normalization, strongly removing technical effects from the data, 
# this means that higher PCs are more likely to represent subtle, but biologically relevant, sources of heterogeneity 
# - so including them may improve downstream analysis. Therefore, higher number of PC can be used. 
# By default we are using the first 50 PCs



# 3. Identifying clusters  #########################
# Seurat implements a graph-based clustering approach. Distances between the spots are calculated based on 
# previously identified PCs. Briefly, Seurat identifies clusters of spots by a shared nearest neighbor (SNN) 
# modularity optimization based clustering algorithm. First, it identifies k-nearest neighbors (KNN) and constructs 
# the SNN graph. Then it optimizes the modularity function to determine clusters. For a full description of the 
# algorithms, see Waltman and van Eck (2013) The European Physical Journal B.

# The FindClusters function implements the procedure, and contains a resolution parameter that sets the granularity 
# of the downstream clustering, with increased values leading to a greater number of clusters.

# Selecting which resolution to use: -------
resolution_vector <- seq(0.2,2,0.2)
integrated <- FindNeighbors(integrated, reduction="pca", dims=1:use.pcs)
integrated <- FindClusters(object=integrated, resolution=resolution_vector, verbose=FALSE)

resTable <- as.data.frame(sapply(grep("res",colnames(integrated@meta.data),value = TRUE), function(x) length(unique(integrated@meta.data[,x]))))
colnames(resTable) <- "number of clusters"

# how many clusters in each resolution
resTable

for(i in seq(resolution_vector)){
  print(DimPlot(integrated, reduction = "umap", label=TRUE, group.by=paste0("SCT_snn_res.", resolution_vector[i])) + labs(title=paste0("res_", resolution_vector[i])))
  print(SpatialDimPlot(integrated, label = TRUE, combine=FALSE, label.size = 3, group.by=paste0("SCT_snn_res.", resolution_vector[i])))
}

# Plotting clustering trees
# To build a clustering tree we need to look at how spots move as the clustering resolution is increased. 
# Each cluster forms a node in the tree and edges are constructed by considering the spots in a cluster 
# at a lower resolution that end up in a cluster at the next higher resolution. By connecting clusters 
# in this way we can see how clusters are related to each other, which are clearly distinct and which are unstable.

# If there are nodes with multiple incoming edges, it is a good indication that we have over-clustered the data.

clustree(integrated, prefis="res.")



# -> Set the resolution to 0.4 ############
resolution <- "0.4" 

# Final clustering:
Idents(integrated) <- paste0("SCT_snn_res.", resolution)
integrated$finalCluster <- Idents(integrated)

DimPlot(integrated, reduction = "umap", label=TRUE, group.by = "ident") 
SpatialDimPlot(integrated, label = TRUE, combine=FALSE, label.size = 3)


# As there are many colors, it can be challenging to visualize which spot belongs to which cluster. 
# Therefore, we highlight spots of each cluster seperately:
SpatialDimPlot(integrated, combine=TRUE, images=name1, cells.highlight = CellsByIdentities(object = integrated), cols.highlight = c("red", "NA"), facet.highlight = TRUE, ncol = 4)
SpatialDimPlot(integrated, combine=TRUE, images=name2, cells.highlight = CellsByIdentities(object = integrated), cols.highlight = c("red", "NA"), facet.highlight = TRUE, ncol = 4)

# Total number of cells per cluster per slice
integrated@meta.data %>% dplyr::select(orig.ident, finalCluster) %>% table() 


# 4. Identify marker genes for each cluster #######################
# Identify unique upregulated marker genes for each clusters (based on normalized data)

## gvg: clustering done, set default assay back to SCT for DGE analyses
DefaultAssay(integrated) <- "SCT"

#prepare object to run differential expression on SCT assay with multiple models
integrated <- PrepSCTFindMarkers(integrated)

saveRDS(integrated, paste0(analysisFolder, "/clustered_res", resolution, ".rds"))

# find marker genes of all clusters (only overexpressed genes -> only.pos=TRUE)
markers_all <- FindAllMarkers(object=integrated, only.pos=TRUE, min.pct = 0.25)
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
DoHeatmap(object = integrated, features = top$gene)


# Dotplot of top 8 **unique marker genes** between clusters:	
DotPlot(integrated, features = top$gene, dot.scale = 8) + coord_flip() +  RotatedAxis()	


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
  
  print(FeaturePlot(integrated, head(markers_all_cluster$gene, n=4), cols = c("lightgrey", "blue"), ncol = 2))
  print(SpatialFeaturePlot(integrated, features = head(markers_all_cluster$gene, n=4), alpha = c(0.1, 1)))
}
  


  
# 5. Identify spatial variable features ##########
# Search for features exhibiting spatial patterning in the absence of pre-annotation. The method, is inspired 
# by the Trendsceek, which models spatial transcriptomics data as a mark point process and computes a variogram, 
# which identifies genes whose expression level is dependent on their spatial location. More specifically, this 
# process calculates gamma(r) values measuring the dependence between two spots a certain "r" distance apart. 
# By default, we use an r-value of '5' in these analyses, and only compute these values for variable genes 
# (where variation is calculated independently of spatial location) to save time.\n\nTop 9 of spatial variable genes:

seu.norm1 <- FindSpatiallyVariableFeatures(seu.norm1, assay = "SCT", features = rownames(GetAssayData(seu.norm1, slot = "scale.data")), selection.method = "markvariogram")
  
spatialFeatures <- SVFInfo(seu.norm1, selection.method="markvariogram", status=TRUE)
spatialFeatures <- na.omit(spatialFeatures)
spatialFeatures <- spatialFeatures %>% filter(variable == TRUE) %>% arrange(rank)
spatialFeatures <- spatialFeatures[,c("r.metric.5"), drop=FALSE]
 
write.xlsx(spatialFeatures, paste0(analysisFolder, "/spatialVariableFeatures_res", resolution, "_", names, ".xlsx"), sheetName="Sheet1", colNames=TRUE, rowNames=TRUE, append=FALSE) 
top.features <- head(SpatiallyVariableFeatures(seu.norm1, selection.method = "markvariogram"), 9)
SpatialFeaturePlot(seu.norm1, features=top.features, ncol=3, alpha = c(0.1, 1))




# 6. Cell type identification (only human or mouse) ##########
# ScType a computational method for automated selection of marker genes based merely on scRNA-seq data.
# Reference from CellMarker 2.0: One of the largest available databases for cell-type marker. Cell markers 
# of different cell types from different tissues with more than one evidence were used.

# Be aware this is an automatic annotation and requires manual curation!
 
species <- "Mouse" 

# load gene set preparation function
library("HGNChelper")
source(paste0(workPath, "scripts/sc-type/gene_sets_prepare.R"))
# load cell type annotation function
source(paste0(workPath, "scripts/sc-type/sctype_score_.R"))
# load tissue type detection function
source(paste0(workPath, "scripts/sc-type/auto_detect_tissue_type.R"))


# create own db (downloaded from CellMarker 2.0) ------
tryCatch({
  own_db_file <- read.xlsx("http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download_files/file/Cell_marker_All2.xlsx")       
},error =function(e){
  own_db_file <<- read.xlsx(paste0(workPath, "scripts/Cell_marker_All.xlsx"))
})
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


# automatically detect tissue type of the dataset ------
# guess a tissue type
# if saled = TRUE, make sure the data is scaled, as seuratObject[[assay]]@scale.data is used. 
# If you just created a Seurat object, without any scaling and normalization, set scaled = FALSE, seuratObject[[assay]]@counts will be used         
tissue_guess <- auto_detect_tissue_type(path_to_db_file = paste0(analysisFolder, "/referenceDB.xlsx"), seuratObject = integrated, scaled = TRUE, assay = "SCT") 
head(tissue_guess)
tissue <-  tissue_guess$tissue[1]



# Top 5 cell type annotation for each cluster
# prepare gene sets
gs_list <- gene_sets_prepare(paste0(analysisFolder, "/referenceDB.xlsx"), tissue)

# get cell-type by cell matrix
es.max <- sctype_score(scRNAseqData = integrated[["SCT"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(integrated@meta.data$finalCluster), function(cl){
  es.max.cl <- sort(rowSums(es.max[ ,rownames(integrated@meta.data[integrated@meta.data$finalCluster==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(integrated@meta.data$finalCluster==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
sctype_scores_top5 <- cL_resutls %>% group_by(cluster) %>% top_n(n = 5, wt = scores)

sctype_scores_top5
# Note: cluster: cluster id, type: cell type annotation, scores: Cell-type specificity score provides a quantitative measure of 
# how uniquely markers identify a specific cell-type of the given tissue. A cell-type with the highest ScType score is used for 
# assignment to the cluster. Low ScType score (less than quarter the number of cells in a cluster), or a negative ScType score 
# indicates a low-confidence cell-type annotation, which will be assigned as “unknown” cell type, ncells: number of cells in this cluster

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

# UMAP
integrated@meta.data$clusterAnnotation = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  integrated@meta.data$clusterAnnotation[integrated@meta.data$finalCluster == j] = as.character(cl_type$type[1])
}

DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'clusterAnnotation')
SpatialDimPlot(integrated, label = TRUE, combine=FALSE, group.by = 'clusterAnnotation', label.size = 3)

# reset finalCluster
integrated$finalCluster_number <- integrated$finalCluster

Idents(integrated) <- integrated@meta.data$clusterAnnotation
integrated$finalCluster <- Idents(integrated)
saveRDS(integrated, paste0(analysisFolder, "/clustered_cellTypeAnnotated_res", resolution, ".rds"))




# 7.  Cell-type deconvolution ##############
# STdeconvolve is an unsupervised, reference-free approach to infer latent cell-type proportions and transcriptional profiles 
# within multi-cellular spatially-resolved pixels from spatial transcriptomics (ST) datasets. STdeconvolve builds on latent 
# Dirichlet allocation (LDA), a generative statistical model commonly used in natural language processing. In the context of 
# ST data, given a count matrix of gene expression in multi-cellular ST pixels, STdeconvolve applies LDA to infer the putative 
# transcriptional profile for each cell-type and the proportional representation of each cell-type in each multi-cellular ST pixel. 

library(STdeconvolve)

# this is the genes x barcode sparse count matrix
cd <- integrated[["SCT"]]@counts

# xy coordinates
addition <- 0
pos <- data.frame()
posList <- list()
for(name in c(name1, name2)){
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
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(6, 12, by = 1),
               perc.rare.thresh = 0.05,
               plot=TRUE,
               verbose=TRUE)

## select model 
optLDA <- optimalModel(models = ldas, opt = "kneed")
optLDA2 <- optimalModel(models = ldas, opt = "min")

if(optLDA2@k > optLDA@k){
  optLDA <- optimalModel(models = ldas, opt = round(mean(c(optLDA2@k, optLDA@k))))
}

# Scatterpies of proportions of all deconvolved cell-types: ------
# Cell-types from spots that contribute less than 5% of the spot proportions are removed and the deconvolved transcriptional profiles are scaled by 1000 

# Extract pixel cell-type proportions (theta) and cell-type gene expression
# profiles (beta) for the given dataset.
# We can also remove cell-types from pixels that contribute less than 5% of the
# pixel proportion and scale the deconvolved transcriptional profiles by 1000 
results <- getBetaTheta(optLDA,
                        perc.filt = 0.05,
                        betaScale = 1000)

deconProp <- results$theta
deconGexp <- results$beta

vizAllTopics(deconProp[rownames(posList[[name1]]),], posList[[name1]], r=3, lwd = 0, showLegend = TRUE, plotTitle = name1) + ggplot2::guides(colour = "none")
vizAllTopics(deconProp[rownames(posList[[name2]]),], posList[[name2]], r=3, lwd = 0, showLegend = TRUE, plotTitle = name2) + ggplot2::guides(colour = "none")



# Compare deconvolved cell-types to clustering ----
# Correlation matrix of the correlations between the proportions of each cell-type and the clustering

cluster_proxyTheta <- model.matrix(~ 0 + integrated@meta.data$finalCluster_number)
rownames(cluster_proxyTheta) <- rownames(integrated@meta.data)

corMat_prop <- getCorrMtx(m1 = as.matrix(cluster_proxyTheta),
                          m2 = deconProp,
                          type = "t")

rownames(corMat_prop) <- paste0("cluster_", levels(integrated@meta.data$finalCluster_number))
colnames(corMat_prop) <- paste0("decon_", seq(ncol(corMat_prop)))

correlationPlot(mat = corMat_prop, colLabs = "Transcriptional clusters", rowLabs = "STdeconvolve") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))


# Marker genes for each deconvoluted cell-type ----
# Marker genes are defined here as genes highly expressed in the deconvolved cell-type (count > 2) that also have the highest 
# log2(fold change) (>1) when comparing the deconvolved cell-type’s expression profile to the average of all other deconvolved 
# cell-types expression profiles.

# Now, let’s get the differentially expressed genes for each deconvolved cell-type transcriptional profile and label the top expressed genes for each cell-type:

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

celltype <- 0
for(markerTable in markerList){
  write.xlsx(markerTable, paste0(analysisFolder, "/STdeconvolve_markerGenes_cellType", i, ".xlsx"), sheetName="Sheet1", colNames=TRUE, rowNames=TRUE, append=FALSE)
  # Note: celltype: deconvolved cell-type, gene: gene name, log2FC: log2 fold-change in comparing deconvolved cell-type expression profile to 
  # the average of all other deconvolved cell-types expression profiles. Positive values indicate that the gene is more highly expressed within 
  # the deconvolved cell-type, exp.profile.celltype: expression profile of deconvolved cell-type, mean.exp.profile.other.celltypes: mean expression 
  # profile of all other deconvolved cell-types

  if(nrow(markerTable) > 0){
    print(FeaturePlot(integrated, head(markerTable$gene, n=4), cols = c("lightgrey", "blue"), ncol = 2))
    print(SpatialFeaturePlot(integrated, features = markerTable$gene[1], alpha = c(0.1, 1)))
  }
}  
  
  
# Cell-type annotation ----
# Given a list of reference gene sets for different cell types, we can performed gene set enrichment analysis (GSEA) 
# on the deconvolved transcriptional profiles to test for significant enrichment of any known ground truth cell-types.
# Reference from CellMarker 2.0: One of the largest available databases for cell-type marker. Human cell markers of 
# different cell types from different tissues with more than one evidence were used.

# Be aware this is an automatic annotation and requires manual curation!

#GSEA annotation 
# prepare gene sets
gs_list <- gene_sets_prepare(paste0(analysisFolder, "/referenceDB.xlsx"), tissue)
colnames(deconGexp) <- toupper(colnames(deconGexp))
celltype_annotations <- annotateCellTypesGSEA(beta = deconGexp , gset = gs_list$gs_positive, qval = 0.05)
celltype_annotation_table <- data.frame()
i <- 0
for(result in celltype_annotations$results){
  i <- i + 1
  if(nrow(result) > 0){
    result2 <- cbind(data.frame(celltype=i, annotation=rownames(result)), result)
    celltype_annotation_table <- rbind(celltype_annotation_table, result2)
  }
}
celltype_annotation_table
# Note: celltype: deconvolved cell-type, annotation: annotation, p.val: p-value of the enrichment score, 
# q.val: q-value estimation for false discovery rate. This is the significance threshold that should be considered,
# sscore: enrichment score represent the degree to which a set is over-represented at the top or bottom of the ranked list, 
# edge: edge score


annotationNames <- as.character(celltype_annotations$predictions)
annotationNames[is.na(annotationNames)] <- "unknown"
colnames(deconProp) <- annotationNames

vizAllTopics(deconProp[rownames(posList[[name1]]),], posList[[name1]], r=2.5, lwd = 0, showLegend = TRUE, plotTitle = name1) + ggplot2::guides(colour = "none")
vizAllTopics(deconProp[rownames(posList[[name2]]),], posList[[name2]], r=2.5, lwd = 0, showLegend = TRUE, plotTitle = name2) + ggplot2::guides(colour = "none")


# Compare deconvolved cell-types annotation to clustering annotation ----
# Correlation matrix of the correlations between the proportions of each cell-type annotation and the clustering annotation
deconPropMerged <- t(apply(deconProp, 1, function(x) tapply(x,colnames(deconProp),sum)))

# correlation with clustering
cluster_proxyTheta <- model.matrix(~ 0 + integrated@meta.data$finalCluster)
rownames(cluster_proxyTheta) <- rownames(integrated@meta.data)

corMat_prop <- getCorrMtx(m1 = as.matrix(cluster_proxyTheta),
                          m2 = deconPropMerged,
                          type = "t")

rownames(corMat_prop) <- levels(integrated@meta.data$finalCluster)

correlationPlot(mat = corMat_prop,  colLabs = "Transcriptional clusters", rowLabs = "STdeconvolve") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

sessionInfo()

