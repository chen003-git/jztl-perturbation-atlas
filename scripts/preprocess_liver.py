import scanpy as sc
import numpy as np
import pandas as pd
from scipy.io import mmread

def read10x_liver(sample_dir, inner):
    base = sample_dir + '/' + inner + '/'
    mtx = mmread(base + 'matrix.mtx').tocsr().T  # cell x gene
    genes = pd.read_csv(base + 'genes.tsv', sep='\t', header=None)
    bc = pd.read_csv(base + 'barcodes.tsv', sep='\t', header=None)
    ad = sc.AnnData(X=mtx,
                    obs=pd.DataFrame(index=bc[0].values),
                    var=pd.DataFrame(index=genes[1].values.astype(str)))
    ad.var_names_make_unique()
    return ad

samples = [
    ('GSM3040892_Liver-10X_P4_2', 'Liver-10X_P4_2'),
    ('GSM3040898_Liver-10X_P7_0', 'Liver-10X_P7_0'),
    ('GSM3040899_Liver-10X_P7_1', 'Liver-10X_P7_1'),
]
ads = []
for d, inner in samples:
    ad = read10x_liver(d, inner)
    ad.obs['sample'] = d
    ads.append(ad)
liver = sc.concat(ads, join='outer')
liver.obs_names_make_unique()
print('raw cells:', liver.n_obs, 'genes:', liver.n_vars)

liver.var['mt'] = liver.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(liver, qc_vars=['mt'], inplace=True)
liver = liver[(liver.obs.n_genes_by_counts >= 200) &
              (liver.obs.n_genes_by_counts < 5000) &
              (liver.obs.pct_counts_mt < 20)].copy()
print('after QC:', liver.n_obs)

sc.pp.normalize_total(liver, target_sum=1e4)
sc.pp.log1p(liver)
sc.pp.highly_variable_genes(liver, n_top_genes=2000, flavor='seurat')
sc.pp.pca(liver, n_comps=50)
sc.pp.neighbors(liver, n_neighbors=15)
sc.tl.umap(liver)
sc.tl.leiden(liver, resolution=0.5)
print('clusters:', liver.obs['leiden'].nunique())

markers = {
    'Hepatocyte': ['Alb', 'Ttr', 'Apoc3'],
    'Stellate': ['Pdgfrb', 'Lrat', 'Des'],
    'Kupffer': ['Clec4f', 'Vsig4', 'Cd163'],
    'LSEC': ['Stab2', 'Cd32b', 'Fcgr2b'],
}
for ct, gs in markers.items():
    sc.tl.score_genes(liver, gs, score_name=ct)

print('\n=== cluster annotation ===')
for cl in sorted(liver.obs['leiden'].unique(), key=int):
    sub = liver.obs[liver.obs['leiden'] == cl]
    scs = {ct: round(float(sub[ct].mean()), 3) for ct in markers}
    top = max(scs, key=scs.get)
    print(f"cluster {cl} (n={len(sub)}): {scs} -> {top}")

pp = 'Pparg' if 'Pparg' in liver.var_names else ('PPARG' if 'PPARG' in liver.var_names else None)
print('\nPPARG gene:', pp)
if pp:
    print('=== PPARG expression per cluster ===')
    for cl in sorted(liver.obs['leiden'].unique(), key=int):
        sub = liver[liver.obs['leiden'] == cl]
        me = float(np.asarray(sub[:, pp].X.mean()))
        pct = float((np.asarray(sub[:, pp].X.toarray()) > 0).mean() * 100)
        print(f"  cl{cl}: Pparg mean={me:.4f} pct={pct:.1f}%")

liver.write('liver_processed.h5ad')
print('\nsaved liver_processed.h5ad')
