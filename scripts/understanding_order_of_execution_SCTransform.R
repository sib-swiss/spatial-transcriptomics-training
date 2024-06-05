# loading > merging > SCTransform > PrepSCTFindMarkers

seu_list <- lapply(c("Anterior", "Posterior"), function(x) {
  seu <- Load10X_Spatial(data.dir = file.path("raw_data", x),
                  slice = x)
  seu$orig.ident <- x
  return(seu)
})

seu <- merge(seu_list[[1]], seu_list[[2]])
seu <- SetIdent(seu, value = "orig.ident")

## does it on the two counts layers separately:
seu <- SCTransform(seu, assay = "Spatial")

VlnPlot(seu, features = "nCount_SCT")

data <- GetAssayData(seu, assay = "SCT", layer = "data")

# Loading > merging > SplitObject > lapply SCTransform > merge

seu_list_v2 <- lapply(c("Anterior", "Posterior"), function(slice) {
  seu <- Load10X_Spatial(data.dir = file.path("raw_data", slice),
                         slice = slice)
  seu$orig.ident <- slice
  return(seu)
})

seu_v2 <- merge(seu_list_v2[[1]], seu_list_v2[[2]])

seu_v2 <- SetIdent(seu_v2, value = "orig.ident")

seu_list_v2 <- SplitObject(seu_v2, split.by = "orig.ident")

# images aren't split with SplitObject. Resetting the images. 
for(slice in names(seu_list_v2)) {
  seu_list_v2[[slice]]@images <- setNames(
    list(seu_list_v2[[slice]]@images[[slice]]),
    slice)
}

seu_list_v2 <- lapply(X = seu_list_v2, FUN = SCTransform, assay = "Spatial")

seu_v2 <- merge(seu_list_v2[[1]], seu_list_v2[[2]])

seu_v2 <- PrepSCTFindMarkers(seu_v2)

data_v2 <- GetAssayData(seu_v2, assay = "SCT", layer = "data")

# loading > SCTransform > merging > PrepSCTFindMarkers 

seu_list_v3 <- lapply(c("Anterior", "Posterior"), function(x) {
  seu <- Load10X_Spatial(data.dir = file.path("raw_data", x),
                         slice = x)
  seu$orig.ident <- x
  
  seu <- SCTransform(seu, assay = "Spatial")
  
  return(seu)
})

seu_v3 <- merge(seu_list_v3[[1]], seu_list_v3[[2]])

seu_v3 <- PrepSCTFindMarkers(seu_v3)

seu_v3 <- SetIdent(seu_v3, value = "orig.ident")

VlnPlot(seu_v3, features = "nCount_SCT")

data_v3 <- GetAssayData(seu_v3, assay = "SCT", layer = "data")

identical(data_v2, data_v3)

plot(colSums(data_v3), colSums(data_v2), col = as.factor(seu$orig.ident))
abline(0,1)
plot(rowSums(data_v3), rowSums(data_v2), col = as.factor(seu$orig.ident))
abline(0,1)

plot(seu$nCount_SCT, seu_v3$nCount_SCT, col = as.factor(seu$orig.ident))

seu_v3@assays$SCT@SCTModel.list

seu@assays$SCT@SCTModel.list
