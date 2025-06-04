import scanpy as sc
import pandas as pd
import anndata as ad
import harmonypy as hm
import numpy as np
import matplotlib.pyplot as plt
import pertpy as pt
from scipy import sparse

import warnings
warnings.filterwarnings('ignore')

metadata = "replace with metadata"
adata = "prepped adata object"
required_metadata_obs = []

# Keep raw
adata.layers["raw"] = adata.X.copy()

# Convert to dense matrix if needed
raw_counts = adata.X
if sparse.issparse(raw_counts):
    raw_counts = raw_counts.toarray()

# DataFrame of raw counts
counts_df = pd.DataFrame(raw_counts, index=adata.obs_names, columns=adata.var_names)
counts_df['sample'] = adata.obs['sample'].values

# Sum counts per gene per sample (pseudobulk)
pseudobulk_counts = counts_df.groupby('sample').sum().T  # Genes x Samples
pseudobulk_counts.to_csv("pseudobulk_counts.csv")

# Prepare metadata
sample_metadata = adata.obs[required_metadata_obs].drop_duplicates().set_index('sample')
sample_metadata.to_csv("pseudobulk_metadata.csv")