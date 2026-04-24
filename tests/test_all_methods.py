import numpy as np
import pandas as pd
import pytest
from py_scmetabolism.methods import aucell_score, gsva_score, ssgsea_score, vision_score
from py_scmetabolism import sc_metabolism, sc_metabolism_anndata


class TestVISION:
    def test_output_shape(self, kegg_expr, kegg_pathways, kegg_gene_names):
        result = vision_score(kegg_expr, kegg_pathways, kegg_gene_names)
        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == len(kegg_pathways)
        assert result.shape[1] == kegg_expr.shape[1]

    def test_scores_finite(self, kegg_expr, kegg_pathways, kegg_gene_names):
        result = vision_score(kegg_expr, kegg_pathways, kegg_gene_names)
        assert np.all(np.isfinite(result.values))

    def test_pathway_names(self, kegg_expr, kegg_pathways, kegg_gene_names):
        result = vision_score(kegg_expr, kegg_pathways, kegg_gene_names)
        assert list(result.index) == list(kegg_pathways.keys())

    def test_nonzero_variance(self, kegg_expr, kegg_pathways, kegg_gene_names):
        result = vision_score(kegg_expr, kegg_pathways, kegg_gene_names)
        row_vars = result.var(axis=1)
        assert (row_vars > 0).sum() > 0


class TestAUCell:
    def test_output_shape(self, kegg_expr, kegg_pathways, kegg_gene_names):
        result = aucell_score(kegg_expr, kegg_pathways, kegg_gene_names)
        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == len(kegg_pathways)
        assert result.shape[1] == kegg_expr.shape[1]

    def test_scores_range(self, kegg_expr, kegg_pathways, kegg_gene_names):
        result = aucell_score(kegg_expr, kegg_pathways, kegg_gene_names)
        assert np.all(result.values >= -1e-10)
        assert np.all(result.values <= 1.0 + 1e-10)

    def test_scores_finite(self, kegg_expr, kegg_pathways, kegg_gene_names):
        result = aucell_score(kegg_expr, kegg_pathways, kegg_gene_names)
        assert np.all(np.isfinite(result.values))

    def test_pathway_names(self, kegg_expr, kegg_pathways, kegg_gene_names):
        result = aucell_score(kegg_expr, kegg_pathways, kegg_gene_names)
        assert list(result.index) == list(kegg_pathways.keys())


class TestSSGSEA:
    def test_output_shape(self, kegg_expr, kegg_pathways, kegg_gene_names):
        result = ssgsea_score(kegg_expr, kegg_pathways, kegg_gene_names, normalize=True)
        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == len(kegg_pathways)
        assert result.shape[1] == kegg_expr.shape[1]

    def test_scores_finite_normalized(self, kegg_expr, kegg_pathways, kegg_gene_names):
        result = ssgsea_score(kegg_expr, kegg_pathways, kegg_gene_names, normalize=True)
        assert np.all(np.isfinite(result.values))

    def test_scores_finite_raw(self, kegg_expr, kegg_pathways, kegg_gene_names):
        result = ssgsea_score(kegg_expr, kegg_pathways, kegg_gene_names, normalize=False)
        assert np.all(np.isfinite(result.values))

    def test_normalize_reduces_range(self, kegg_expr, kegg_pathways, kegg_gene_names):
        result_norm = ssgsea_score(kegg_expr, kegg_pathways, kegg_gene_names, normalize=True)
        result_raw = ssgsea_score(kegg_expr, kegg_pathways, kegg_gene_names, normalize=False)
        assert result_norm.shape == result_raw.shape
        norm_range = result_norm.values.max() - result_norm.values.min()
        raw_range = result_raw.values.max() - result_raw.values.min()
        if raw_range > 0:
            assert norm_range <= raw_range + 1e-10


class TestGSVA:
    def test_output_shape(self, kegg_expr, kegg_pathways, kegg_gene_names):
        result = gsva_score(kegg_expr, kegg_pathways, kegg_gene_names)
        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == len(kegg_pathways)
        assert result.shape[1] == kegg_expr.shape[1]

    def test_scores_finite(self, kegg_expr, kegg_pathways, kegg_gene_names):
        result = gsva_score(kegg_expr, kegg_pathways, kegg_gene_names)
        assert np.all(np.isfinite(result.values))

    def test_pathway_names(self, kegg_expr, kegg_pathways, kegg_gene_names):
        result = gsva_score(kegg_expr, kegg_pathways, kegg_gene_names)
        assert list(result.index) == list(kegg_pathways.keys())


class TestScMetabolism:
    def test_dataframe_input(self, kegg_adata):
        count_df = pd.DataFrame(
            kegg_adata.X,
            index=kegg_adata.obs_names,
            columns=kegg_adata.var_names,
        )
        result = sc_metabolism(count_df, method="VISION", metabolism_type="KEGG")
        assert isinstance(result, pd.DataFrame)
        assert result.shape[1] == kegg_adata.n_obs

    def test_anndata_input(self, kegg_adata):
        result = sc_metabolism(kegg_adata, method="VISION", metabolism_type="KEGG")
        assert isinstance(result, pd.DataFrame)

    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_ndarray_input(self, kegg_adata):
        result = sc_metabolism(kegg_adata.X, method="VISION", metabolism_type="KEGG")
        assert isinstance(result, pd.DataFrame)

    def test_invalid_method(self, kegg_adata):
        with pytest.raises(ValueError, match="Unknown method"):
            sc_metabolism(kegg_adata, method="INVALID")

    def test_invalid_metabolism_type(self, kegg_adata):
        with pytest.raises(ValueError):
            sc_metabolism(kegg_adata, method="VISION", metabolism_type="INVALID")


class TestScMetabolismAnnData:
    def test_key_added(self, kegg_adata):
        result = sc_metabolism_anndata(
            kegg_adata.copy(),
            method="VISION",
            metabolism_type="KEGG",
            key_added="metabolism_test",
        )
        assert "X_metabolism_test" in result.obsm
        assert "metabolism_test_pathways" in result.uns

    def test_vision(self, kegg_adata):
        result = sc_metabolism_anndata(
            kegg_adata.copy(),
            method="VISION",
            metabolism_type="KEGG",
        )
        assert "X_metabolism" in result.obsm
        assert "metabolism_pathways" in result.uns
        assert result.obsm["X_metabolism"].shape[0] == kegg_adata.n_obs
        assert np.all(np.isfinite(result.obsm["X_metabolism"]))

    def test_aucell(self, kegg_adata):
        result = sc_metabolism_anndata(
            kegg_adata.copy(),
            method="AUCell",
            metabolism_type="KEGG",
        )
        assert "X_metabolism" in result.obsm
        assert result.obsm["X_metabolism"].shape[0] == kegg_adata.n_obs
        assert np.all(np.isfinite(result.obsm["X_metabolism"]))

    def test_ssgsea(self, kegg_adata):
        result = sc_metabolism_anndata(
            kegg_adata.copy(),
            method="ssGSEA",
            metabolism_type="KEGG",
        )
        assert "X_metabolism" in result.obsm
        assert result.obsm["X_metabolism"].shape[0] == kegg_adata.n_obs
        assert np.all(np.isfinite(result.obsm["X_metabolism"]))

    def test_gsva(self, kegg_adata):
        result = sc_metabolism_anndata(
            kegg_adata.copy(),
            method="GSVA",
            metabolism_type="KEGG",
        )
        assert "X_metabolism" in result.obsm
        assert result.obsm["X_metabolism"].shape[0] == kegg_adata.n_obs
        assert np.all(np.isfinite(result.obsm["X_metabolism"]))
