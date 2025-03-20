.libPaths(.libPaths()[c(2,1)])
library(Seurat)     # single cell RNA-Seq analysis
library(kableExtra) # format tables
library(ggplot2)   # create graphics
library(viridis)   # accessible color palettes

## Set up proper project path
#data.dir <- "/share/bioinfo/Workshops/Spatial_SingleCell_March2025/00-RawData/Xenium/Xenium_Prime_Human_Lung_Cancer_FFPE_outs"
data.dir <- "/share/bioinfo/projects/Internal_Jessie_UCD/Workshops/scRNA_Spatial/Mar2025/00-RawData/Xenium.Multi"
project.dir <- "/share/bioinfo/projects/Internal_Jessie_UCD/Workshops/scRNA_Spatial/Mar2025/02-Seurat-Xenium"

## samples
samples <- c("TgCRND8", "wildtype")
sample.dir <- paste0("Xenium_V1_FFPE_", samples, "_5_7_months_outs")

## Source a modified version of Seurat's LoadXenium function to accommodate the new Xenium output files
source("Seurat_functions.R")

## load Xenium data
### old version Xenium use Seurat::LoadXenium with uncompressed cell_feature_matrix folder;
### new version Xenium use myLoadXenium because of the missing transcripts.csv.gz, but still require uncompressed cell_feature_matrix folder
#xenium <- myLoadXenium(file.path(data.dir, sample.dir[1], segmentations = "cell", cell.centroids = T)
#xenium <- Seurat::LoadXenium(data.dir = file.path(data.dir, sample.dir[1]), segmentations = "cell", cell.centroids = T)

options(future.globals.maxSize = 2000 * 1024 ^ 2)
## load multiple slices
merge.data <- lapply(seq_along(samples), function(i){
	## specify the fov slot names for easier access later
	tmp <- Seurat::LoadXenium(data.dir = file.path(data.dir, sample.dir[i]), segmentations = "cell", cell.centroids = T, fov = paste0("fov-", samples[i]))
	tmp <- RenameCells(tmp, new.names = paste0(samples[i], "-", Cells(tmp)))
	## change the active.idents slot
	Idents(tmp) <- samples[i]
	return(tmp)
})

for (i in seq_along(merge.data)){
	if (i == 1){
		experiment.merged <- merge.data[[i]]
	}else{
		experiment.merged <- merge(experiment.merged, merge.data[[i]])
	}
}


## explore the Seurat object

### 10X documentation on codewords: https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/algorithms-overview/xoa-algorithms#qvs
### genes: designed probes for the gene panel
VlnPlot(experiment.merged, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0, log = T)

### blank, unassigned codewords: are unused codewords. No probe in the corresponding gene panel that will generate the codeword
VlnPlot(experiment.merged, features = c("nFeature_BlankCodeword", "nCount_BlankCodeword"), ncol = 2, pt.size = 0.01)

### negative control codewords: codewords that do not have any probes matching that code. They are used to assess the specificity of the decoding algoritm
VlnPlot(experiment.merged, features = c("nFeature_ControlCodeword", "nCount_ControlCodeword"), ncol = 2, pt.size = 0.01)

### control, negative control probe: probes that exist in the panel but do not target any biological sequences. They are used to assess the specificity of the assay
VlnPlot(experiment.merged, features = c("nFeature_ControlProbe", "nCount_ControlProbe"), ncol = 2, pt.size = 0.01)

### genomic control: designed to bind to intergenic genomic DNA but not to any transcript sequence present in the tissue. They are present in the Xenium Prime assay, but not in earlier Xenium assays. The human lung data has it.
#VlnPlot(experiment.merged, features = c("nFeature_GenomicControl", "nCount_GenomicControl"), ncol = 2, pt.size = 0.01)

## filtering cells based on above features
experiment.merged <- subset(experiment.merged, subset = nCount_Xenium > 20)

## explore data: FOV metadata
## FOV-class documentation from Seurat: https://search.r-project.org/CRAN/refmans/SeuratObject/html/FOV-class.html
str(experiment.merged)
## plot some gene expressions
ImageDimPlot(experiment.merged, fov = "fov", molecules = c(
## zoom in to specific FOV
fov.coords <- Crop(experiment.merged[["fov"]], x = c(290, 480), y = c(11451, 7190), coords = "plot")
experiment.merged[["fov.zoom"]] <- fov.coords
DefaultBoundary(experiment.merged[["fov.zoom"]]) <- "segmentation"
ImageDimPlot(experiment.merged, fov = "fov.zoom", axes = T, border.color = "white", border.size = 0.1, cols = "polychrome",
	coord.fixed = F, molecules = c(

## Normalization
experiment.merged <- SCTransform(experiment.merged, assay = "Xenium")

## Dimentionality reduction
experiment.merged <- RunPCA(experiment.merged, npcs = 30, features = rownames(experiment.merged))
experiment.merged <- RunUMAP(experiment.merged, dims = 1:30)
experiment.merged <- FindNeighbors(experiment.merged, reduction = "pca", dims = 1:30)
## resolution 0.8:34, 0.9:36
## resolution 1.0:35, 1.5:47, 2.0:55 too high
## resolution 0.3:20, 0.5:27
experiment.merged <- FindClusters(experiment.merged, resolution = c(0.3, 0.5))
saveRDS(experiment.merged, file="merged.rds")

## Visualization of the clustering results
Idents(experiment.merged) <- "SCT_snn_res.0.5"
DimPlot(experiment.merged)
ImageDimPlot(experiment.merged, fov = "fov", cols = DiscretePalette(n = 27), group.by = "SCT_snn_res.0.5", size = 0.5, alpha = 0.7)
ImageDimPlot(experiment.merged, fov = "fov", cols = "polychrome", group.by = "SCT_snn_res.0.5", size = 0.5, alpha = 0.7)
Idents(experiment.merged) <- "orig.ident"
ImageDimPlot(experiment.merged, fov = "fov", molecules = c("Satb2", "Lamp5", "Car4"), group.by = NULL, size = 0.5, alpha = 0.5, axes = T)
ImageFeaturePlot(experiment.merged, fov = "fov", features = c("Satb2", "Lamp5", "Car4"), size = 0.5, alpha = 0.7, axes = T)
### Seurat by default flip x and y axes, using flip_xy to restore original layout
ImageDimPlot(experiment.merged, fov = "fov", molecules = c("Satb2", "Lamp5", "Car4"), group.by = NULL, size = 0.5, alpha = 0.5, axes = T, flip_xy = F, coord.fixed = T)
### my modified ImageFeaturePlot to flip the x and y axes back to original
myImageFeaturePlot(experiment.merged, fov = "fov", features = c("Satb2", "Lamp5", "Car4"), size = 0.5, alpha = 0.7, axes = T, flip_xy = F, coord.fixed = T, combine = F, blend = F)

markers.2 <- FindMarkers(experiment.merged, ident.1 = 2)

##> head(markers.2)
##        p_val avg_log2FC pct.1 pct.2 p_val_adj
##Satb2       0   2.791416 0.844 0.160         0
##Lamp5       0   2.986377 0.849 0.191         0
##Pdzrn3      0   3.549878 0.699 0.084         0
##Car4        0   2.941589 0.763 0.165         0
##Slc17a7     0   1.610459 0.963 0.391         0
##Cux2        0   2.874990 0.845 0.280         0


## crop to a specific fov
crop <- Crop(experiment.merged[["fov"]], x = 

## save expression data for cell type annotation using Allen Brain atlas
## save dimreductions, graphs and images for scNiche, not working because the spatial information is not saved
#sxenium <- DietSeurat(object = experiment.merged, counts = T, data = T, scale.data = F, assays = c("SCT"))
#saveRDS(sxenium, file="mergeddata4scNiche.rds")
#as.anndata(x = sxenium, file_path="./", file_name="mergeddata4scNiche.h5ad", assay = "SCT", main_layer = "counts", other_layer = "data", transer_dimreduc = T)
## save smaller data for cell type annotation on web
sxenium <- DietSeurat(object = experiment.merged, counts = T, data = T, scale.data = F, assays = c("SCT"), dimreducs = NULL, graphs = NULL)
slot(sxenium, "images") <- list(fov.TgCRND8 = NULL, fov.wildtype = NULL)
slot(sxenium, "commands") <- list(NULL)
saveRDS(sxenium, file="mergeddata4celltype.rds")

## convert seurat to anndata for Allen Brain atlas online use
##library(scCustomize)
use_python("/share/bioinfo/projects/Internal_Jessie_UCD/software/anaconda3.2430.anndata/bin/python")
seurat <- readRDS("mergeddata4celltype.rds")
## change feature names to ensembl ids
## This section is not necessary because of --map_to_ensembl parameter
jsonlite::read_json(file.path(sample.dir[1], "gene_panel.json")
out <- do.call(rbind, lapply(geneid2ensembl[["payload"]][["targets"]], function(x){
	return(data.frame(Ensembl = ifelse(length(x$type$data$id) < 1, NA, x$type$data$id[1]), GeneName = ifelse(length(x$type$data$name) < 1, NA, x$type$data$name[1])))
}))
sct <- seurat[["SCT"]]
newnms <- out$Ensembl[match(seurat[["SCT"]]@counts@Dimnames[[1]], out$GeneName)]
sct@counts@Dimnames[[1]] <- newnms
sct@data@Dimnames[[1]] <- newnms
rownames(sct@scale.data) <- newnms
sct@var.features <- out$Ensembl[match(sct@var.features, out$GeneName)]
rownames(sct@meta.features) <- newnms
rownames(sct[["SCT"]]@SCTModel.list$counts@feature.attributes) <- sct[["SCT"]]@counts@Dimnames[[1]]
seurat[["SCT"]] <- sct
saveRDS(seurat, file="data4celltype_ensembl.rds")
seurat <- readRDS("data4celltype_ensembl.rds")
## anndata requires python anndata module /share/bioinfo/projects/Internal_Jessie_UCD/software/anaconda3.2430.celltypemapper has it
as.anndata(x = seurat, file_path="./", file_name="mergeddata4celltype.h5ad", assay = "SCT", main_layer = "counts", other_layer = "data", transer_dimreduc = F)

## cell type annotation using Allen Brain atlas 10X whole mouse brain dataset
## API: https://knowledge.brain-map.org/mapmycells/process
## CLI:
## module load anaconda3./24.3.0
## conda create -p /share/bioinfo/projects/Internal_Jessie_UCD/software/anaconda3.2430.celltypemapper python=3.10	the lowest working python version
## conda install anndata/matplotlib -c conda-forge
## git clone https://github.com/AllenInstitute/cell_type_mapper.git
## cd cell_type_mapper
## pip install .	"pip install -e ." suggested in github doesn't work because of "-e" editable installation, which is not necessary
## export PYTHONPATH=/share/bioinfo/projects/Internal_Jessie_UCD/software/anaconda3.2430.celltypemapper/lib/python3.10/site-packages/:$PYTHONPATH
## the precomputed and markers were downloaded using the commands from section 8 at /share/bioinfo/projects/Internal_Jessie_UCD/Workshops/scRNA_Spatial/Mar2025/02-Seurat-Xenium
## command: python -m cell_type_mapper.cli.from_specified_markers --query_path /share/bioinfo/projects/Internal_Jessie_UCD/Workshops/scRNA_Spatial/Mar2025/02-Seurat-Xenium/data4celltype.h5ad --extended_result_path /share/bioinfo/projects/Internal_Jessie_UCD/Workshops/scRNA_Spatial/Mar2025/02-Seurat-Xenium/celltypes_from_10X_whole_brain.json --csv_result_path /share/bioinfo/projects/Internal_Jessie_UCD/Workshops/scRNA_Spatial/Mar2025/02-Seurat-Xenium/celltypes_from_10X_whole_brain.csv --drop_level CCN20230722_SUPT --cloud_safe False --query_markers.serialized_lookup /share/bioinfo/projects/Internal_Jessie_UCD/Workshops/scRNA_Spatial/Mar2025/00-RawData/Allen_atlas/mouse_markers_230821.json --precomputed_stats.path /share/bioinfo/projects/Internal_Jessie_UCD/Workshops/scRNA_Spatial/Mar2025/00-RawData/Allen_atlas/precomputed_stats_ABC_revision_230821.h5 --type_assignment.normalization raw --type_assignment.n_processors 12
## --map_to_ensembl provides mapping from gene symbol to ensembl id
## --type_assignment.normalization raw will convert data to log2(CPM + 1) before mapping, log2CPM will use the data as it is
## the web version is 1.4.0 and the github version is 1.5.0, the results are mostly identical with 53015 shared celltype annotations out of 57983 cells.

celltypes <- read.csv("mergedcelltypes.web.csv", header=T, quote="")


## multiple slices
wt <- Seurat::LoadXenium(data.dir = file.path(data.dir, sample.dir[2]), segmentations = "cell", cell.centroids = T)


## scNiche
### module load anaconda3/24.3.0
### conda activate /share/bioinfo/projects/Internal_Jessie_UCD/software/anaconda-install/envs/scniche
### export PYTHONPATH=/share/bioinfo/projects/Internal_Jessie_UCD/software/anaconda-install/envs/scniche/lib/python3.9/site-packages

