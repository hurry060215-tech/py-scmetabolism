"""Shared fixtures for py-scmetabolism tests.

We build small synthetic datasets with known pathway gene expression patterns
so that tests can verify scoring algorithms produce expected trends.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import anndata as ad
import pathlib


@pytest.fixture(scope="session")
def small_random_adata():
    """Small synthetic 200 cells × 500 genes with random negative binomial counts.

    Used for basic functionality tests of scoring methods.
    """
    rng = np.random.default_rng(42)
    n_cells, n_genes = 200, 500
    # Negative binomial counts with some zero inflation
    counts = rng.negative_binomial(n=10, p=0.1, size=(n_cells, n_genes))
    # Add sparsity (~70% zeros)
    counts = counts * (rng.random(size=(n_cells, n_genes)) > 0.7)

    obs = pd.DataFrame(index=[f"Cell{i+1}" for i in range(n_cells)])
    var = pd.DataFrame(index=[f"Gene{i+1}" for i in range(n_genes)])
    return ad.AnnData(counts, obs=obs, var=var)


@pytest.fixture(scope="session")
def adata_with_pathway_genes():
    """Synthetic data where specific gene sets (mimicking pathways) are upregulated.

    Creates three artificial 'pathways':
      - Pathway_A: genes 0-9 upregulated in cells 0-49
      - Pathway_B: genes 10-19 upregulated in cells 50-99
      - Pathway_C: genes 20-29 upregulated in cells 100-149
    """
    rng = np.random.default_rng(123)
    n_cells, n_genes = 150, 100
    X = rng.poisson(5.0, (n_cells, n_genes)).astype(np.float64)

    # Pathway A markers
    X[:50, :10] += rng.poisson(15.0, (50, 10)).astype(np.float64)
    # Pathway B markers
    X[50:100, 10:20] += rng.poisson(15.0, (50, 10)).astype(np.float64)
    # Pathway C markers
    X[100:, 20:30] += rng.poisson(15.0, (50, 10)).astype(np.float64)

    obs = pd.DataFrame(
        {"group": ["A"] * 50 + ["B"] * 50 + ["C"] * 50},
        index=[f"c{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"g{i}" for i in range(n_genes)])
    return ad.AnnData(X=X, obs=obs, var=var)


@pytest.fixture(scope="session")
def kegg_gmt_path():
    """Return path to KEGG metabolism GMT file included with the package."""
    data_dir = pathlib.Path(__file__).parent.parent / "py_scmetabolism" / "data"
    gmt_path = data_dir / "KEGG_metabolism_nc.gmt"
    if not gmt_path.exists():
        pytest.skip(f"KEGG GMT file not found at {gmt_path}")
    return gmt_path


@pytest.fixture(scope="session")
def reactome_gmt_path():
    """Return path to REACTOME metabolism GMT file included with the package."""
    data_dir = pathlib.Path(__file__).parent.parent / "py_scmetabolism" / "data"
    gmt_path = data_dir / "REACTOME_metabolism_nc.gmt"
    if not gmt_path.exists():
        pytest.skip(f"REACTOME GMT file not found at {gmt_path}")
    return gmt_path