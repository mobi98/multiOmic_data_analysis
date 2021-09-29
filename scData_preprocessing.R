## scRNA-seq and scATAC-seq integration ##

library(Signac)
library(Seurat)
library(ggplot2)
library(EnsDb.Mmusculus.v79) # load genome to map to, can be human, mouse etc
set.seed(1234)

# load data 

base <- "/Users/marthaobrien/Documents/imperial2020_21/data_analysis_project"
rna <- Read10X(data.dir = paste0(base, '/data/GSE126074_AdBrainCortex/cDNA'), gene.column = 1) # default Read10X selects gene.column 2  
atac <- Read10X(data.dir = paste0(base, "/data/GSE126074_AdBrainCortex/chromatin"), gene.column = 1)
fragments <- "/Users/marthaobrien/Documents/imperial2020_21/data_analysis_project/data/GSE126074_AdBrainCortex/chromatin/fragments.sort.bed.gz"   
# need 2 files: 'fragments.sort.bed.gz' and 'fragments.sort.bed.gz.tbi'; they need to be in same directory  
# if fragment files are not available for SNARE-seq, files can be created using following code: https://github.com/timoast/SNARE-seq

# create seurat object and add the assays

snare <- CreateSeuratObject(counts = rna)
snare[['ATAC']] <- CreateChromatinAssay(
  counts = atac,
  sep = c(':','-'),
  genome = 'mm10',
  fragments = fragments 
)

# extract the annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(snare[["ATAC"]]) <- annotations

# QC and filtering of the ATAC assay within the Seurat object

DefaultAssay(snare) <- "ATAC" 
snare <- TSSEnrichment(snare) # Transcriptional start site (TSS) enrichment score
snare <- NucleosomeSignal(snare) # Nucleosome banding pattern
snare$blacklist_fraction <- FractionCountsInRegion(
  object = snare,
  assay = 'ATAC',
  regions = blacklist_mm10
) #blacklist represents reads that are associated with aretfactual signals


Idents(snare) <- "all"  # group all cells together, rather than by replicate # idents gets the cell identity classses
VlnPlot( # creating a violin plot of the data
  snare,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment",
               "nucleosome_signal", "blacklist_fraction"),
  pt.size = 0.1,
  ncol = 5
)

# filtering the data so that only features with certain conditions are retained 
snare <- subset(
  x = snare,
  subset = blacklist_fraction < 0.03 &
    TSS.enrichment < 20 &
    nCount_ATAC > 500
  )

snare <- RunTFIDF(snare) # Run term frequency inverse document frequency (TF-IDF) normalization on the matrix.
snare <- FindTopFeatures(snare, assay = "ATAC", min.cutoff = 10) # min num cells must have certain feature 
snare <- RunSVD(snare, assay = "ATAC", n = 50) # Run partial singular value decomposition using irlba
# snare <- RunPCA(snare, assay = 'ATAC', npcs= 50, reduction.name = 'atac.pca')



# QC and filtering of RNA data within Seurat object

# can do additional cell and gene level filtering if required
DefaultAssay(snare) <- "RNA"
snare$log10GenesPerUmi <- log10(rna_seurat$nFeature_RNA) / log10(rna_seurat$nCount_RNA)
snare <- subset(x = snare, 
                         subset= (nCount_RNA >= 500) & 
                           (nFeature_RNA >= 250) & 
                           (log10GenesPerUmi > 0.80))


snare <- FindVariableFeatures(snare, nfeatures = 3000)
snare <- NormalizeData(snare)
snare <- ScaleData(snare)
snare <- RunPCA(snare, npcs = 50, reduction.name = "rna.pca")


# plotting of features 

#Plotting RNA with t-SNE and UMAP
snare <- RunTSNE(snare, reduction.name = 'rna.tsne', reduction = "rna.pca", tsne.method = 'Rtsne')
snare <- FindNeighbors(snare, dims = 1:30)
snare <- FindClusters(snare, resolution = 0.5, algorithm = 3)
DimPlot(snare, reduction = 'umap.rna') + ggtitle('UMAP of Adult Mouse Brain RNA Expression') + theme(title = element_text(hjust = 0.5))

DimPlot(snare, reduction = 'rna.tsne') + ggtitle('t-SNE of scRNA-seq data from Adult Mouse Brain')

#snare <- RunUMAP(snare, dims = 1:30, reduction.name = "umap.rna")
#DimPlot(snare, reduction = 'umap.rna', label = T) + ggtitle("UMAP of scRNA-seq data from Adult Mouse Brain")

# plotting ATAC with t-SNE and UMAP
snare <- RunTSNE(object = snare, assay = 'ATAC', reduction = 'lsi', dims = 2:30, reduction.name = 'atac.tsne', tsne.method = 'Rtsne')
DimPlot(snare, reduction = 'atac.tsne', label = T) + ggtitle('t-SNE of scATAC-seq data from Adult Mouse Brain') + labs(x = 't_SNE1', y = 't_SNE2')
#snare <- RunUMAP(snare, assay = "ATAC", reduction.name = 'atac.umap', dims = 2:30)
#DimPlot(snare, reduction = 'atac.umap', label = T) + ggtitle('UMAP of scATAc-seq data from Adult Mouse Brain') + labs(x = 'UMAP_1', y = "UMAP_2")


# Extracting data to run in Multi-SNE

# for RNA can select all variable features (8055 x 3000)
rna_assay <- GetAssay(object = snare, assay = "RNA")
rna_var_genes <- VariableFeatures(rna_assay)
rna_df <- as.data.frame(rna_assay@scale.data[rna_var_genes,])
rna_df <- t(rna_df)
rna_df <- unname(rna_df, force = T)

# OR select PCA embeddings (8055 x 50)

rna_pca <- as.matrix(unname(snare@reductions$rna.pca@cell.embeddings, force = TRUE))


# For ATAC can do similar but 'Variable Features' are often still too high

atac_assay <- GetAssay(object = snare2, assay = "ATAC")
atac_var.features <- VariableFeatures(atac_assay)
atac_df <- as.data.frame(atac_assay@data[atac_var.features,]) # no scale.data slot for ATAC as this Signac vignette didn't run ScaleData on ATAC (but I think you can scale it)
atac_df <- t(atac_df)
atac_df <- unname(atac_df, force = T)

# running PCA
atac_pca <- prcomp(atac_df)
atac_pca <- atac_pca$x[,1:50]

# better to use CisTopics for ATAC

library(cisTopic)

# inbuilt cisTopic function didnt work so redefined it a bit
createcisTopicObject_multiOme <- function (count.matrix, project.name = "cisTopicProject", getRegionRange=FALSE,
min.cells = 1, min.regions = 1, is.acc = 1, keepCountsMatrix = TRUE,
...) {
cisTopic.version <- packageVersion("cisTopic")
object <- new(Class = "cisTopic", is.acc = is.acc,
project.name = project.name, version = cisTopic.version)
object.binary.count.matrix <- Matrix::Matrix(1 * (count.matrix >=
is.acc), sparse = TRUE)
num.acc.cells <- Matrix::rowSums(object.binary.count.matrix)
num.acc.regions <- Matrix::colSums(object.binary.count.matrix)
cells.use <- which(num.acc.regions >= min.regions)
regions.use <- which(num.acc.cells >= min.cells)
count.matrix <- count.matrix[regions.use, cells.use]
object@cell.names <- colnames(count.matrix)
object@region.names <- rownames(count.matrix)
seqnames <- sapply(strsplit(rownames(count.matrix), split = ":"),
"[", 1)
coord <- sapply(strsplit(rownames(count.matrix), split = ":"),
"[", 2)
start <- sapply(strsplit(coord, split = "-"), "[",
1)
end <- sapply(strsplit(coord, split = "-"), "[",
2)
bed_coord <- cbind(seqnames, start, end)
rownames(bed_coord) <- rownames(count.matrix)
if (getRegionRange == TRUE) {
bed_coord[,2] <- as.factor(bed_coord[,2])
bed_coord[,3] <- as.factor(bed_coord[,3])
object@region.ranges <- makeGRangesFromDataFrame(as.data.frame(bed_coord))
}
nCounts_celldata <- Matrix::colSums(count.matrix)
nCounts_regiondata <- Matrix::rowSums(count.matrix)
if (keepCountsMatrix == TRUE) {
object@count.matrix <- count.matrix
}
rm(count.matrix)
object.binary.count.matrix <- object.binary.count.matrix[regions.use,
cells.use]
nCounts <- nCounts_celldata
rm(nCounts_celldata)
nAcc <- Matrix::colSums(object.binary.count.matrix)
object.cell.data <- cbind(nCounts, nAcc)
object.cell.data <- apply(object.cell.data, 2, function(x) as.numeric(as.character(x)))
rownames(object.cell.data) <- object@cell.names
object@cell.data <- as.data.frame(object.cell.data)
nCounts <- nCounts_regiondata
rm(nCounts_regiondata)
nCells <- as.numeric(Matrix::rowSums(object.binary.count.matrix))
width <- abs(as.numeric(end) - as.numeric(start))
object.region.data <- cbind(seqnames, start, end, width,
nCounts, nCells)
rownames(object.region.data) <- object@region.names
object@region.data <- as.data.frame(object.region.data)
object@region.data[, 2:ncol(object@region.data)] <- apply(object@region.data[,
2:ncol(object@region.data)], 2, function(x) as.numeric(as.character(x)))
object@binary.count.matrix <- object.binary.count.matrix
rm(object.binary.count.matrix)
object@calc.params[["createcisTopicObject"]] <- list(min.cells = min.cells,
min.regions = min.regions, is.acc = is.acc)
return(object)
}

atac_counts <- snare@assays$ATAC@counts # extract the ATAC counts from the snare assay

cisTopicObject <- createcisTopicObject_multiOme(atac_counts, project.name = "cisTopicRevised") # create object with counts

cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(10:15, 20, 25, 35, 40, 45, 50), seed=123, nCores=20, addModels=FALSE)

# cisObject <- runCGSModels(cisObject, topic=c(2,10,15,20,25,30,35,40), seed=987, nCores=9, burnin = 120, iterations = 150, addModels=FALSE)

cisTopicObject <- selectModel(cisTopicObject, type='perplexity') # selecting model with best log likelihood if don't specify

atac_data <- cisTopicObject@selected.model$document_expects

# Running multi-SNE

snareList <- list(rna_pca, atac_data)
adBrain_multiSNE <- multiSNE::multiSNE(snareList)


## scRNA-seq and surface protein expression integration ##


setwd("/Users/marthaobrien/Documents/imperial2020_21/data_analysis_project/data/cite_seq_data")

cbmc.rna <- as.sparse(read.csv(file = "GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", sep = ",",
    header = TRUE, row.names = 1))

cbmc.adt <-  as.sparse(read.csv(file = "GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", sep = ",",
    header = TRUE, row.names = 1))

cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)

# initialise Seurat object from RNA assay
cbmc <- CreateSeuratObject(counts = cbmc.rna)

adt_assay <- CreateAssayObject(counts = cbmc.adt) # add the ADT surface protein expression assay 
cbmc[['ADT']] <- adt_assay
Assays(cbmc)

# RNA gene level filtering
cbmc <- subset(x = cbmc, subset = nFeature_RNA > 250 &
                 nFeature_RNA < 2500)

# standard workflow: RNA
DefaultAssay(cbmc) <- "RNA"
cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc, nfeatures = 1000)
cbmc <- ScaleData(cbmc)
cbmc <- Seurat::RunPCA(cbmc, npcs = 50, reduction.name = 'rna.pca')

cbmc <- Seurat::RunTSNE(object = cbmc, assay = 'RNA', reduction = 'rna.pca', dims = 1:30, reduction.name = 'rna.tsne', tsne.method = 'Rtsne')
cbmc <- Seurat::FindNeighbors(cbmc, dims = 1:30, reduction = 'rna.pca', assay = 'RNA')
cbmc <- Seurat::FindClusters(cbmc)

# standard workflow: ADT
DefaultAssay(cbmc) <- "ADT"
cbmc <- NormalizeData(cbmc, normalization.method = "CLR", margin = 2, assay = 'ADT')
cbmc <- ScaleData(cbmc)
cbmc <- FindVariableFeatures(cbmc, nfeatures = 10)
cbmc <- RunPCA(object = cbmc, npcs = 10, reduction.name = 'adt.pca')
cbmc <- RunTSNE(object = cbmc, assay = "ADT", reduction = 'adt.pca', dims = 1:10, reduction.name = 'adt.tsne', tsne.method = 'Rtsne')
adt.dist <- as.matrix(t(cbmc@assays$ADT@data))
cbmc <- FindNeighbors(cbmc, assay = "ADT", reduction = 'adt.pca', dims = 1:9)
cbmc <- FindClusters(object = cbmc, assay = "ADT")

# running multiSNE

# with raw features 
rna_assay <- Seurat::GetAssay(cbmc, assay = 'RNA')
adt_assay <- Seurat::GetAssay(cbmc, assay = 'ADT')

rna.var.genes <- VariableFeatures(rna_assay) # extracting top 1000 variable genes
adt.var.feat <- VariableFeatures(adt_assay) # extracting 10 variable proteins 

cbmc_rna <- t(as.matrix(unname(rna_assay@scale.data[rna.var.genes,], force = T)))
cbmc_adt <- t(as.matrix(unname(adt_assay@scale.data[adt.var.feat,], force = T)))
# cbmc_adt <- t(as.matrix(unname(adt_assay@scale.data, force = T))) using all 13 proteins instead of variable 10 

# using PCA features
rna_pca <- as.matrix(unname(cbmc@reductions$rna.pca@cell.embeddings, force = TRUE))
adt_pca <- as.matrix(unname(cbmc@reductions$adt.pca@cell.embeddings, force =  TRUE))

cbmcList <- list(rna_pca, adt_pca)
# cbmcList <- list(cbmc_rna, cbmc_adt)

cbmc_multiSNE <- multiSNE::multiSNE(cbmcList)
