from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import anndata as ad
import pathlib


def _load_gmt(gmt_path):
    pathways = {}
    with open(gmt_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                pathways[parts[0]] = parts[2:]
    return pathways


@pytest.fixture(scope="session")
def kegg_gmt_path():
    data_dir = pathlib.Path(__file__).parent.parent / "py_scmetabolism" / "data"
    gmt_path = data_dir / "KEGG_metabolism_nc.gmt"
    if not gmt_path.exists():
        pytest.skip(f"KEGG GMT file not found at {gmt_path}")
    return gmt_path


@pytest.fixture(scope="session")
def reactome_gmt_path():
    data_dir = pathlib.Path(__file__).parent.parent / "py_scmetabolism" / "data"
    gmt_path = data_dir / "REACTOME_metabolism.gmt"
    if not gmt_path.exists():
        pytest.skip(f"REACTOME GMT file not found at {gmt_path}")
    return gmt_path


@pytest.fixture(scope="session")
def kegg_pathways(kegg_gmt_path):
    return _load_gmt(kegg_gmt_path)


@pytest.fixture(scope="session")
def kegg_gene_names(kegg_pathways):
    all_genes = set()
    for genes in kegg_pathways.values():
        all_genes.update(genes)
    return sorted(all_genes)


@pytest.fixture(scope="session")
def kegg_adata(kegg_gene_names):
    rng = np.random.default_rng(42)
    n_cells = 30
    n_genes = len(kegg_gene_names)
    counts = rng.poisson(5.0, (n_cells, n_genes)).astype(np.float64)
    obs = pd.DataFrame(index=[f"Cell{i+1}" for i in range(n_cells)])
    var = pd.DataFrame(index=kegg_gene_names)
    return ad.AnnData(counts, obs=obs, var=var)


@pytest.fixture(scope="session")
def kegg_expr(kegg_adata):
    X = kegg_adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    return X.T.copy()
