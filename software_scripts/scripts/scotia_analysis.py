import scotia
import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
adata = sc.read_h5ad("TgCRND8.h5ad")

lr_pairs = pd.read_csv("lr.txt")
adata.obs['sample'] = "TgCRND8"   
scotia.run_scotia.lr_score(adata = adata, lr_list = lr_pairs, sample_col = 'sample', fov_col = 'region', celltype_col = "celltype_annotation", output_path = "./scotia")


