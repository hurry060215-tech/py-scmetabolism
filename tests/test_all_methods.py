"""Comprehensive test: compare Python implementations vs R results."""
import numpy as np
import pandas as pd
from scipy import stats
import sys
sys.path.insert(0, '..')

from py_scmetabolism.methods import aucell_score, ssgsea_score, gsva_score, vision_score

test_dir = "tests/test_data_real_genes"
counts = pd.read_csv(f"{test_dir}/test_counts_real.csv", index_col=0)
expr = counts.values.T  # genes x cells
gene_names = list(counts.columns)
n_genes, n_cells = expr.shape
print(f"Expression: {n_genes} genes x {n_cells} cells")

def load_gmt(gmt_file):
    pathways = {}
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                pathways[parts[0]] = parts[2:]
    return pathways

pathways = load_gmt("py_scmetabolism/data/KEGG_metabolism_nc.gmt")

def compare(py_scores, r_scores_df, method_name):
    common = [p for p in r_scores_df.index if p in py_scores.index]
    if not common:
        print(f"  {method_name}: No common pathways!")
        return

    per_pathway_cors = []
    worst_pathway = None
    worst_corr = 1.0
    for p_name in common:
        r_vals = r_scores_df.loc[p_name].values.astype(float)
        py_vals = py_scores.loc[p_name].values.astype(float)
        if np.std(r_vals) > 1e-10 and np.std(py_vals) > 1e-10:
            c = np.corrcoef(r_vals, py_vals)[0, 1]
            per_pathway_cors.append(c)
            if c < worst_corr:
                worst_corr = c
                worst_pathway = p_name
        elif np.std(r_vals) < 1e-10 and np.std(py_vals) < 1e-10:
            per_pathway_cors.append(1.0)

    if per_pathway_cors:
        cors = np.array(per_pathway_cors)
        print(f"  {method_name}:")
        print(f"    Per-pathway corr: mean={np.mean(cors):.6f}, min={np.min(cors):.6f}, "
              f"median={np.median(cors):.6f}")
        print(f"    Pathways >=0.99: {np.sum(cors >= 0.99)}/{len(cors)}")
        print(f"    Pathways >=0.95: {np.sum(cors >= 0.95)}/{len(cors)}")
        if worst_pathway:
            print(f"    Worst: {worst_pathway} ({worst_corr:.6f})")

        # Overall correlation
        all_r = np.concatenate([r_scores_df.loc[p].values.astype(float) for p in common])
        all_py = np.concatenate([py_scores.loc[p].values.astype(float) for p in common])
        valid = np.isfinite(all_r) & np.isfinite(all_py) & (np.abs(all_r) > 1e-15)
        if np.sum(valid) > 2:
            overall = np.corrcoef(all_r[valid], all_py[valid])[0, 1]
            print(f"    Overall corr: {overall:.6f}")

# AUCell
print("\n=== AUCell ===")
r_aucell = pd.read_csv(f"{test_dir}/R_results/R_AUCell_scores.csv", index_col=0)
py_aucell = aucell_score(expr, pathways, gene_names)
compare(py_aucell, r_aucell, "AUCell")

# ssGSEA (no normalization)
print("\n=== ssGSEA (no normalization) ===")
r_ssgsea = pd.read_csv(f"{test_dir}/R_results/R_ssGSEA_nonorm_scores.csv", index_col=0)
py_ssgsea = ssgsea_score(expr, pathways, gene_names)
compare(py_ssgsea, r_ssgsea, "ssGSEA")

# ssGSEA (normalized)
print("\n=== ssGSEA (normalized) ===")
r_ssgsea_norm = pd.read_csv(f"{test_dir}/R_results/R_ssGSEA_norm_scores.csv", index_col=0)
py_ssgsea_norm = ssgsea_score(expr, pathways, gene_names, normalize=True)
compare(py_ssgsea_norm, r_ssgsea_norm, "ssGSEA (norm)")

# GSVA
print("\n=== GSVA ===")
r_gsva = pd.read_csv(f"{test_dir}/R_results/R_GSVA_scores.csv", index_col=0)
py_gsva = gsva_score(expr, pathways, gene_names)
compare(py_gsva, r_gsva, "GSVA")

print("\nDone!")
