#!/usr/bin/env python
"""Quick sanity test of py-scmetabolism package."""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import anndata as ad
from py_scmetabolism import sc_metabolism_anndata

print("Testing py-scmetabolism import and basic functionality...")

# Create tiny dataset
n_cells, n_genes = 10, 50
counts = np.random.poisson(5.0, (n_cells, n_genes))
adata = ad.AnnData(
    counts,
    obs=pd.DataFrame(index=[f"Cell{i}" for i in range(n_cells)]),
    var=pd.DataFrame(index=[f"Gene{i}" for i in range(n_genes)])
)

# Test VISION method
adata_vision = sc_metabolism_anndata(
    adata.copy(),
    method="VISION",
    imputation=False,
    ncores=1,
    metabolism_type="KEGG",
    key_added="metabolism_vision"
)

print(f"  VISION scores shape: {adata_vision.obsm['X_metabolism_vision'].shape}")
print(f"  Number of pathways: {len(adata_vision.uns['metabolism_vision_pathways'])}")
print("  Test passed!")

# Test AUCell method
adata_aucell = sc_metabolism_anndata(
    adata.copy(),
    method="AUCell",
    imputation=False,
    ncores=1,
    metabolism_type="KEGG",
    key_added="metabolism_aucell"
)
print(f"  AUCell scores shape: {adata_aucell.obsm['X_metabolism_aucell'].shape}")
print("  All tests completed successfully.")