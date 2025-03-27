import spatialdata as sd
from spatialdata_io import xenium

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq

# define Xenium output folder
xenium_path = "/share/workshop/Spatial_scRNA_workshop/Data/Xenium.Multi/Xenium_V1_FFPE_TgCRND8_5_7_months_outs"

# read in Xenium data
sdata = xenium(xenium_path)
adata = sdata.tables["table"]

# filter cells using minimum transcripts / cell > 20
sc.pp.filter_cells(adata, min_counts = 20)

# normalization, log-transform, dimensionality reduction and clustering
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, inplace = True)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)

# visualize the the total number of transcripts / cell, the clustering results on UMAP and spatial visualization of the clusters
from matplotlib.backends.backend_pdf import PdfPages

with PdfPages('output.pdf') as pdf:
  pdf.savefig(sc.pl.umap(adata, color = ["total_counts", "leiden",], wspace = 0.4,))
  pdf.savefig(sq.pl.spatial_scatter(adata, library_id="spatial", shape = None, color=["leiden",], wspace = 0.4,))

pdf.close()

import pandas as pd

## read in celltypes
celltypes = pd.read_csv("TgCRND8.celltypes.csv")
keys = celltypes['cell_id'].tolist()
values = celltypes['class_name'].tolist()
d = dict([(k, v) for k, v in zip(keys, values)])
adata.obs["celltype_annotation"] = adata.obs['cell_id'].map(d)

# save anndata object to file
import anndata as anndata

anndata.io.write_h5ad("TgCRND8.h5ad", adata)

