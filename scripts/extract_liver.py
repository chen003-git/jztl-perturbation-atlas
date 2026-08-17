import scanpy as sc
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix

def read10x_liver(base_dir):
    mtx = mmread(base_dir + '/matrix.mtx').tocsr().T
    genes = pd.read_csv(base_dir + '/genes.tsv', sep='\t', header=None)
    bc = pd.read_csv(base_dir + '/barcodes.tsv', sep='\t', header=None)
    ad = sc.AnnData(X=mtx, obs=pd.DataFrame(index=bc[0].values),
                    var=pd.DataFrame(index=genes[1].values.astype(str)))
    ad.var_names_make_unique()
    return ad

dirs = ['GSM3040892_Liver-10X_P4_2/Liver-10X_P4_2',
        'GSM3040898_Liver-10X_P7_0/Liver-10X_P7_0',
        'GSM3040899_Liver-10X_P7_1/Liver-10X_P7_1']
ads = [read10x_liver(d) for d in dirs]
for ad, d in zip(ads, dirs):
    ad.obs['sample'] = d
liver = sc.concat(ads, join='outer')
liver.obs_names_make_unique()

liver.var['mt'] = liver.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(liver, qc_vars=['mt'], inplace=True)
liver = liver[(liver.obs.n_genes_by_counts >= 200) &
              (liver.obs.n_genes_by_counts < 5000) &
              (liver.obs.pct_counts_mt < 20)].copy()
print('肝 QC后:', liver.n_obs)

liver.raw = liver.copy()
sc.pp.normalize_total(liver, target_sum=1e4)
sc.pp.log1p(liver)
sc.pp.highly_variable_genes(liver, n_top_genes=2000, flavor='seurat')
sc.pp.pca(liver, n_comps=50, random_state=42)
sc.pp.neighbors(liver, n_neighbors=15, random_state=42)
sc.tl.umap(liver, random_state=42)
sc.tl.leiden(liver, resolution=0.5, random_state=42)

markers = {'Hepatocyte': ['Alb', 'Ttr', 'Apoc3'], 'Stellate': ['Pdgfrb', 'Lrat', 'Des'],
           'Kupffer': ['Clec4f', 'Vsig4', 'Cd163'], 'LSEC': ['Stab2', 'Fcgr2b']}
for ct, gs in markers.items():
    sc.tl.score_genes(liver, gs, score_name=ct)

cluster_scores = {}
for cl in liver.obs['leiden'].unique():
    sub = liver.obs[liver.obs['leiden'] == cl]
    scs = {ct: float(sub[ct].mean()) for ct in markers}
    top = max(scs, key=scs.get)
    cluster_scores[cl] = (top, scs)
    print(f"  cl{cl} (n={len(sub)}): {scs} -> {top}")

target_clusters = [cl for cl, (top, scs) in cluster_scores.items() if top == 'Hepatocyte']
mask = liver.obs['leiden'].isin(target_clusters)
sub = liver[mask].copy()
print('肝细胞:', sub.n_obs, 'clusters', target_clusters)

raw_counts = sub.raw[:, sub.var_names].X
mat = csr_matrix(raw_counts)
df = pd.DataFrame(mat.toarray().T, index=sub.var_names, columns=sub.obs_names)
df.to_csv('liver_counts.csv')
print('liver_counts.csv:', df.shape)
