import scanpy as sc
import numpy as np
import pandas as pd
from scipy.io import mmread

def read10x(prefix):
    mtx = mmread(prefix + '_matrix.mtx.gz').tocsr().T  # cell x gene
    genes = pd.read_csv(prefix + '_genes.tsv.gz', sep='\t', header=None)
    bc = pd.read_csv(prefix + '_barcodes.tsv.gz', sep='\t', header=None)
    ad = sc.AnnData(X=mtx,
                    obs=pd.DataFrame(index=bc[0].values),
                    var=pd.DataFrame(index=genes[1].values.astype(str)))
    ad.var_names_make_unique()
    return ad

samples = ['GSM5319987_sham1', 'GSM5319988_sham2', 'GSM5319989_sham3',
           'GSM5319990_MCAO1', 'GSM5319991_MCAO2', 'GSM5319992_MCAO3']
ads = [read10x(s) for s in samples]
for ad, s in zip(ads, samples):
    ad.obs['sample'] = s
    ad.obs['condition'] = 'MCAO' if 'MCAO' in s else 'sham'
brain = sc.concat(ads, join='outer')
brain.obs_names_make_unique()
print('raw cells:', brain.n_obs, 'genes:', brain.n_vars)

# QC
brain.var['mt'] = brain.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(brain, qc_vars=['mt'], inplace=True)
brain = brain[(brain.obs.n_genes_by_counts >= 200) &
              (brain.obs.n_genes_by_counts < 5000) &
              (brain.obs.pct_counts_mt < 20)].copy()
print('after QC:', brain.n_obs)

# normalize + cluster
sc.pp.normalize_total(brain, target_sum=1e4)
sc.pp.log1p(brain)
sc.pp.highly_variable_genes(brain, n_top_genes=2000, flavor='seurat')
sc.pp.pca(brain, n_comps=50)
sc.pp.neighbors(brain, n_neighbors=15)
sc.tl.umap(brain)
sc.tl.leiden(brain, resolution=0.5)
print('clusters:', brain.obs['leiden'].nunique())

# marker score annotation
markers = {
    'Endo': ['Cldn5', 'Pecam1', 'Cdh5'],
    'Micro': ['Cx3cr1', 'Tmem119', 'P2ry12'],
    'Astro': ['Gfap', 'Aqp4', 'Aldh1l1'],
    'Peri': ['Pdgfrb', 'Cspg4', 'Rgs5'],
    'Oligo': ['Mog', 'Mbp', 'Plp1'],
}
for ct, gs in markers.items():
    sc.tl.score_genes(brain, gs, score_name=ct)

print('\n=== cluster annotation ===')
for cl in sorted(brain.obs['leiden'].unique(), key=int):
    sub = brain.obs[brain.obs['leiden'] == cl]
    scs = {ct: round(float(sub[ct].mean()), 3) for ct in markers}
    top = max(scs, key=scs.get)
    print(f"cluster {cl} (n={len(sub)}): {scs} -> {top}")

# PPARG expression
pp = 'Pparg' if 'Pparg' in brain.var_names else ('PPARG' if 'PPARG' in brain.var_names else None)
print('\nPPARG gene:', pp)
if pp:
    print('=== PPARG expression per cluster ===')
    for cl in sorted(brain.obs['leiden'].unique(), key=int):
        sub = brain[brain.obs['leiden'] == cl]
        me = float(np.asarray(sub[:, pp].X.mean()))
        pct = float((np.asarray(sub[:, pp].X.toarray()) > 0).mean() * 100)
        print(f"  cl{cl}: Pparg mean={me:.4f} pct={pct:.1f}%")

brain.write('brain_processed.h5ad')
print('\nsaved brain_processed.h5ad')
