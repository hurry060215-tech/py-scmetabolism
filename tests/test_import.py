import numpy as np
from py_scmetabolism import sc_metabolism, sc_metabolism_anndata


def test_import_sc_metabolism():
    assert callable(sc_metabolism)


def test_import_sc_metabolism_anndata():
    assert callable(sc_metabolism_anndata)


def test_import_methods():
    from py_scmetabolism.methods import aucell_score, gsva_score, ssgsea_score, vision_score

    assert callable(aucell_score)
    assert callable(ssgsea_score)
    assert callable(gsva_score)
    assert callable(vision_score)


def test_import_visualize():
    from py_scmetabolism import boxplot_metabolism, dimplot_metabolism, dotplot_metabolism

    assert callable(dimplot_metabolism)
    assert callable(dotplot_metabolism)
    assert callable(boxplot_metabolism)


def test_load_kegg_gmt():
    from py_scmetabolism.compute import load_metabolism_gmt

    pathways = load_metabolism_gmt("KEGG")
    assert isinstance(pathways, dict)
    assert len(pathways) > 0
    for name, genes in pathways.items():
        assert isinstance(name, str)
        assert isinstance(genes, list)
        assert len(genes) > 0


def test_vision_via_anndata(kegg_adata):
    result = sc_metabolism_anndata(
        kegg_adata.copy(),
        method="VISION",
        imputation=False,
        ncores=1,
        metabolism_type="KEGG",
    )
    assert "X_metabolism" in result.obsm
    assert "metabolism_pathways" in result.uns
    assert result.obsm["X_metabolism"].shape[0] == kegg_adata.n_obs
    assert len(result.uns["metabolism_pathways"]) > 0
    assert np.all(np.isfinite(result.obsm["X_metabolism"]))


def test_aucell_via_anndata(kegg_adata):
    result = sc_metabolism_anndata(
        kegg_adata.copy(),
        method="AUCell",
        imputation=False,
        ncores=1,
        metabolism_type="KEGG",
    )
    assert "X_metabolism" in result.obsm
    assert result.obsm["X_metabolism"].shape[0] == kegg_adata.n_obs
    assert np.all(np.isfinite(result.obsm["X_metabolism"]))
